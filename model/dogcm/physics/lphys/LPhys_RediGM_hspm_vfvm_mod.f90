!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief This module to calculate tendecies of potential temperature and salinity due to isopycnal diffusion and GM advection.  
!! 
!! @author Yuta Kawai
!!
!!
module LPhys_RediGM_hspm_vfvm_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  !* Dennou-OGCM

  use ProfUtil_mod
  
  use DOGCM_Admin_Constants_mod, only: &
       & PI,                           &
       & RPlanet, Omega, Grav,         &
       & RefDens

  use SpmlUtil_mod, only: &
       & w_AlphaOptr_xy,  &
       & xy_w, w_xy,      &
       & xya_wa, wa_xya,  &
       & xy_GradLon_w,    &
       & xy_GradLat_w,    &
       & xyz_DSig_xyz,    &
       & xy_CosLat

  use EOSDriver_mod, only: &
       & EOSDriver_Eval,      &
       & EOSDriver_alpha_beta

  use DOGCM_Admin_Grid_mod, only: &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,   &
       & IBLOCK, JBLOCK, KBLOCK,               &
       & iMax, jMax, kMax, lMax, tMax,         &
       & xyz_Lat, xyz_Lon,                     &
       & x_CI, y_CJ, z_CK,                     &
       & z_CDK, z_RCDK, z_FDK, z_RFDK
#include "../../admin/DOGCM_Admin_GaussSpmGridIndexDef.h"
  
  use DOGCM_Admin_TInteg_mod, only: &
       & CurrentTime
  
  use LPhys_RediGMHelper_mod, only: &
       & LPhys_RediGMHelper_Init, LPhys_RediGMHelper_Final, &
       & prepare_SlopeTapering,                               &       
       & xyz_Dz_xyz,                                          &
       & DFM08Info, prepare_DFM08Info, TaperingDFM08_GM, TaperingDFM08_IDIFF
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: LPhys_RediGM_hspm_vfvm_Init, LPhys_RediGM_hspm_vfvm_Final
  public :: LPhys_RediGM_hspm_vfvm_GetParameters
  public :: LPhys_RediGM_hspm_vfvm_Output, LPhys_RediGM_hspm_vfvm_PrepareOutput

  public :: LPhys_RediGM_hspm_vfvm_AddMixingTerm

  interface calc_IsoNeutralSlope
     module procedure calc_IsoNeutralSlope_new
  end interface
  public :: calc_IsoNeutralSlope
  
  public :: calc_IsopycDiffFlux
  public :: calc_BolusVelocity
  
  ! 公開変数
  ! Public variable
  !

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'LPhys_RediGM_hspm_vfvm_mod' !< Module Name

  real(DP) :: SGSEddyMixType
  real(DP) :: Kappa_Redi              !< Isopycnal diffusivity with Redi scheme [m^2/s]
  real(DP) :: Kappa_GM                !< Diffusivity with GM scheme             [m^2/s]

  type(gt_history), save :: hst_SGSEddyMix
  logical :: OutputFlag


  logical :: DFM08Flag

  integer, parameter :: LON = 1
  integer, parameter :: LAT = 2
  integer, parameter :: PTEMP = 1
  integer, parameter :: SALT  = 2

  
contains
  
  !> Initialize this module.
  !!
  !! Choose a scheme to parametrize sub-grid scale eddy mixing.
  !! If Redi scheme is used, set SGSEddyMixParamType to \link #LPhys_RediGM_hspm_vfvm_Redi \endlink. 
  !! If GM scheme is used, set SGSEddyMixParamType to \link #LPhys_RediGM_hspm_vfvm_GM \endlink. 
  !! KappaRedi and KappaGM are the parameters associated with the magnitude of diffusivity tensor in Redi or GM scheme,
  !! respectively(see a document of Dennou-OGCM for details).
  !! If you want, you can also set the parameters and output flag for analysis through namelist. In that case, specify
  !! the confignmlFileName.
  !!
  subroutine LPhys_RediGM_hspm_vfvm_Init( &
       & KappaRedi, KappaGM, isVarsOutput, confignmlFileName & ! (in)
       & )

    use DOGCM_Admin_TInteg_mod, only: &
         & OriginTime => RestartTime, &
         & EndTime
    
    use DOGCM_IO_History_mod, only: &
         & FilePrefix, OutputIntrvalSec

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), optional  :: KappaRedi           !< Isopycnal diffusivity, parameter associated with the magnitude of diffusivity tensor with Redi scheme [m^2/s]
    real(DP), intent(in), optional  :: KappaGM             !< Parameter associated with the magnitude of diffusivity tensor with GM scheme [m^2/s]
    logical, intent(in), optional :: isVarsOutput
    character(*), intent(in), optional :: confignmlFileName   !< Namelist name

    ! 実行文; Executable statements
    !

    ! Set default values.
    Kappa_Redi = 1000d0
    Kappa_GM   = 1000d0
    OutputFlag = .false.

    ! Set some parameters from arguments
    if(present(KappaRedi)) Kappa_Redi = KappaRedi
    if(present(KappaGM)) Kappa_GM = KappaGM
    if(present(isVarsOutput)) OutputFlag = isVarsOutput

    ! If configNmlFileName is specfied, we read namelist file to set the values of parameters.
    if(present(configNmlFileName)) then
       call read_nmlData(configNmlFileName)
    end if

    ! Initialize a module to help this module
    call LPhys_RediGMHelper_Init()

    if (OutputFlag) then
       call LPhys_RediGM_hspm_vfvm_PrepareOutput(   &
            & OriginTime, EndTime, OutputIntrvalSec, FilePrefix & ! (in)
            & )
       call LPhys_RediGM_hspm_vfvm_Output()
    end if
    
  end subroutine LPhys_RediGM_hspm_vfvm_Init

  !-------------------------------------------------------------------
  
  !>
  !!
  !!
  subroutine LPhys_RediGM_hspm_vfvm_Final()

    ! 実行文; Executable statements
    !

    call LPhys_RediGMHelper_Final()
    
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)
    
  end subroutine LPhys_RediGM_hspm_vfvm_Final

  !-------------------------------------------------------------------
  
  !>
  !!
  !!
  subroutine LPhys_RediGM_hspm_vfvm_GetParameters( &
       & KappaRedi, KappaGM                   & ! (out)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out), optional :: KappaRedi
    real(DP), intent(out), optional :: KappaGM

    ! 実行文; Executable statements
    !
    
    if(present(KappaRedi)) KappaRedi = Kappa_Redi
    if(present(KappaGM)) KappaGM = Kappa_GM

  end subroutine LPhys_RediGM_hspm_vfvm_GetParameters
  
  !> @brief 
  !!
  !!
  subroutine LPhys_RediGM_hspm_vfvm_AddMixingTerm( &
       & xyz_PTemp_RHS, xyz_Salt_RHS, xyz_VDiffCoef,           &  ! (inout)
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo            &  ! (in)
       & )

    ! モジュール引用; Use statements
    !    
    use DOGCM_Admin_Constants_mod

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_PTemp_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xyz_Salt_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DensPot(IA,JA,KA)
    real(DP) :: xy_RefPress(IA,JA)

    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(IA,JA,KA)
    real(DP) :: xyz_SLat(IA,JA,KA)
    real(DP) :: xyr_SLon(IA,JA,KA)
    real(DP) :: xyr_SLat(IA,JA,KA)

    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(IA,JA,KA)

    !
    real(DP) :: xyzaa_HGradTRC(IA,JA,KA,2,2)
    real(DP) :: xyra_DSigTRC(IA,JA,KA,2)
    
    integer :: k
    
    ! 実行文; Executable statement
    !
    
    ! Calculate the potential density.
    !

    !$omp parallel do private(xy_RefPress)
    do k=KS,KE
       xy_RefPress(:,:) = 0d0
       call EOSDriver_Eval( xyz_DensPot(IS:IE,JS:JE,k),                                    & ! (out)
            & xyz_PTemp(IS:IE,JS:JE,k), xyz_Salt(IS:IE,JS:JE,k), xy_RefPress(IS:IE,JS:JE) )  ! (in)
    end do

    ! Calculate the components of the isoneutral slope. 
    !

    call calc_GradTRC( &
       & xyzaa_HGradTRC, xyra_DSigTRC,            & ! (out)
       & xyz_PTemp, xyz_Salt                      & ! (in)
       & )    

    call calc_IsoNeutralSlope_new( &
         & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,            & ! (out)
         & xyz_PTemp, xyz_Salt, xyzaa_HGradTRC, xyra_DSigTRC, & ! (in)
         & xyz_H, xyz_Z )                                       !(in)

    call prepare_SlopeTapering( xyz_T, xyz_SLon, xyz_SLat, & ! (inout)
         & xyz_DensPot, xyz_Z                              & ! (in)
         & )
!!$    !$omp parallel do
!!$    do k = KS, KE
!!$       xyz_VDiffCoef(:,:,k) = xyz_VDiffCoef(:,:,k) +                             &
!!$            & Kappa_Redi*xyz_T(:,:,k)**2 *(xyz_SLon(:,:,k)**2 + xyz_SLat(:,:,k)**2)
!!$    end do

    !    call check_StaticStability(xyz_T, xyz_DensPot)

    
    ! Calculate the tendency due to Gent-McWilliams and Redi flux
    !

    call append_Redi_GM_RHS( xyz_PTemp_RHS,    & ! (inout)
         & xyz_PTemp, PTEMP )                    ! (in)    

    call append_Redi_GM_RHS( xyz_Salt_RHS,     & ! (inout)
         & xyz_Salt, SALT )                      ! (in)    
    
  contains

    subroutine append_Redi_GM_RHS( xyz_RHS, & ! (inout)
         xyz_TRC, TRCID )                     ! (in)

      ! 宣言文; Declaration statement
      !
      real(DP), intent(inout) :: xyz_RHS(IA,JA,KA)
      real(DP), intent(in) :: xyz_TRC(IA,JA,KA)
      integer, intent(in) :: TRCID
      
      ! 局所変数
      ! Local variables
      !
      real(DP) :: xyz_FLonRedi(IA,JA,KA)
      real(DP) :: xyz_FLatRedi(IA,JA,KA)
      real(DP) :: xyr_FSigRedi(IA,JA,KA)

      real(DP) :: xyz_FLonGM(IA,JA,KA)
      real(DP) :: xyz_FLatGM(IA,JA,KA)
      real(DP) :: xyr_FSigGM(IA,JA,KA)

      integer :: i
      integer :: j
      integer :: k
      
      ! 実行文; Executable statement
      !

      call calc_IsopycDiffFlux( &
           & xyz_FLonRedi, xyz_FLatRedi, xyr_FSigRedi,                     &  ! (out)      
           & xyz_TRC, xyz_H, xyz_Z,                                        &  ! (in)
           & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat, xyz_T,                &  ! (in)
           & xyzaa_HGradTRC(:,:,:,:,TRCID), xyra_DSigTRC(:,:,:,TRCID)      &  ! (in)
           & )

      call calc_SkewFlux( &
           & xyz_FLonGM, xyz_FLatGM, xyr_FSigGM,                           &  ! (out)      
           & xyz_TRC, xyz_H, xyz_Z,                                        &  ! (in)
           & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat, xyz_T,                &  ! (in)
           & xyzaa_HGradTRC(:,:,:,:,TRCID), xyra_DSigTRC(:,:,:,TRCID)      &  ! (in)
           & )
      
      !$omp parallel private(i,j,k)
      !$omp do
      do k=KS, KE
         xyz_RHS(IS:IE,JS:JE,k) = xyz_RHS(IS:IE,JS:JE,k) &
              & + xy_w( w_AlphaOptr_xy(                                                  &
              &     (xyz_FLonRedi(IS:IE,JS:JE,k) + xyz_FLonGM(IS:IE,JS:JE,k))*xy_CosLat, &
              &     (xyz_FLatRedi(IS:IE,JS:JE,k) + xyz_FLatGM(IS:IE,JS:JE,k))*xy_CosLat  &
              &    ) )                                                                   

         !$omp simd
         do j=JS, JE
         do i=IS, IS+_IM_-1
            xyz_RHS(i,j,k) = xyz_RHS(i,j,k)  &
                 & + (   (xyr_FSigRedi(i,j,k-1) + xyr_FSigGM(i,j,k-1))            &
                 &     - (xyr_FSigRedi(i,j,k  ) + xyr_FSigGM(i,j,k  ))            & 
                 &   )*z_RCDK(k)/xyz_H(i,j,k)
         end do
         end do
      end do
      !$omp end parallel
      
    end subroutine append_Redi_GM_RHS

  end subroutine LPhys_RediGM_hspm_vfvm_AddMixingTerm

  !-------------------------------------------------------------------
  
  subroutine calc_GradTRC( &
       & xyzaa_HGradTRC, xyra_DSigTRC,            & ! (out)
       & xyz_PTemp, xyz_Salt                      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyzaa_HGradTRC(IA,JA,KA,2,2)
    real(DP), intent(out) :: xyra_DSigTRC(IA,JA,KA,2)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: n
    integer :: l
    
    real(DP) :: m1
    real(DP) :: m2

    real(DP) :: w_TRC(lMax)
    
    ! 実行文; Executable statement
    !
    
    !$omp parallel private(i,j,jj,k,n)

    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyra_DSigTRC(i,j,k,PTEMP) = (xyz_PTemp(i,j,k) - xyz_PTemp(i,j,k+1))*z_RFDK(k)
       xyra_DSigTRC(i,j,k,SALT ) = (xyz_Salt (i,j,k) - xyz_Salt (i,j,k+1))*z_RFDK(k)       
    end do
    end do
    end do

    !$omp do collapse(2)
    do n=1, 2
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyra_DSigTRC(i,j,KS-1,n) = xyra_DSigTRC(i,j,KS  ,n)
       xyra_DSigTRC(i,j,KE  ,n) = xyra_DSigTRC(i,j,KE-1,n)
    end do
    end do
    end do
    !$omp end parallel

    !$omp parallel private(w_TRC,i,j,k,l,n)
    !$omp do
    do k=KS, KE
       w_TRC(:) = w_xy(xyz_PTemp(IS:IE,JS:JE,k))
       xyzaa_HGradTRC(IS:IE,JS:JE,k,LON,PTEMP) = xy_GradLon_w(w_TRC)/RPlanet
       xyzaa_HGradTRC(IS:IE,JS:JE,k,LAT,PTEMP) = xy_GradLat_w(w_TRC)/RPlanet
    end do
    !$omp do
    do k=KS, KE
       w_TRC(:) = w_xy(xyz_Salt(IS:IE,JS:JE,k))
       xyzaa_HGradTRC(IS:IE,JS:JE,k,LON,SALT) = xy_GradLon_w(w_TRC)/RPlanet
       xyzaa_HGradTRC(IS:IE,JS:JE,k,LAT,SALT) = xy_GradLat_w(w_TRC)/RPlanet
    end do

    !$omp do collapse(3)
    do n=1, 2
    do l=1, 2
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyzaa_HGradTRC(i,j,KS-1,l,n) = xyzaa_HGradTRC(i,j,KS,l,n)
       xyzaa_HGradTRC(i,j,KE+1,l,n) = xyzaa_HGradTRC(i,j,KE,l,n)
    end do
    end do
    end do
    end do
    !$omp end parallel
    
  end subroutine calc_GradTRC

  !-------------------------------------------------------------------
  
  subroutine calc_IsopycDiffFlux( &
       & xyz_FLon, xyz_FLat, xyr_FSig,            &  ! (out)
       & xyz_TRC, xyz_H, xyz_Z,                   &  ! (in)
       & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (in)
       & xyz_T,                                   &  ! (in)
       & xyza_HGradTRC, xyr_DSigTRC               &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(IA,JA,KA)
    real(DP), intent(out) :: xyz_FLat(IA,JA,KA)
    real(DP), intent(out) :: xyr_FSig(IA,JA,KA)
    real(DP), intent(in)  :: xyz_TRC(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)    
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyr_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyr_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_T(IA,JA,KA)
    real(DP), intent(in) :: xyza_HGradTRC(IA,JA,KA,2)
    real(DP), intent(in) :: xyr_DSigTRC(IA,JA,KA)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_DzTRC(IA,JA)
    real(DP) :: DzTRC
    real(DP) :: r_m1(KA)
    real(DP) :: r_m2(KA)
    real(DP) :: SLon
    real(DP) :: SLat

    integer :: i
    integer :: j
    integer :: jj
    integer :: k

    
    ! 実行文; Executable statement
    !

    
    ! Calculate the components of diffusive flux along isopycnal surface.
    !

    do k=KS, KE-1
       r_m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       r_m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
    end do

    !$omp parallel private(DzTRC,SLon,SLat,i,j,jj,k)

    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       DzTRC = 0.5d0*(xyr_DSigTRC(i,j,k-1) + xyr_DSigTRC(i,j,k))/xyz_H(i,j,k)

       xyz_FLon(i,j,k) = Kappa_Redi*(                                                      &
            &   xyza_HGradTRC(i,j,k,LON) + xyz_T(i,j,k)*xyz_SLon(i,j,k)*DzTRC )  

       xyz_FLat(i,j,k) = Kappa_Redi*(                                                      &
            &   xyza_HGradTRC(i,j,k,LAT) + xyz_T(i,j,k)*xyz_SLat(i,j,k)*DzTRC )  
    end do
    end do
    end do

    !$omp do collapse(2)
    do k=KS, KE-1
    do j=JS, JE
    do i=IS, IS+_IM_-1
       SLon = r_m1(k)*xyz_T(i,j,k)*xyz_SLon(i,j,k) + r_m2(k)*xyz_T(i,j,k+1)*xyz_SLon(i,j,k+1)
       SLat = r_m1(k)*xyz_T(i,j,k)*xyz_SLat(i,j,k) + r_m2(k)*xyz_T(i,j,k+1)*xyz_SLat(i,j,k+1)

       xyr_FSig(i,j,k) = Kappa_Redi*( &
            &     SLon*(r_m1(k)*xyza_HGradTRC(i,j,k,LON) + r_m2(k)*xyza_HGradTRC(i,j,k+1,LON))                     &
            &   + SLat*(r_m1(k)*xyza_HGradTRC(i,j,k,LAT) + r_m2(k)*xyza_HGradTRC(i,j,k+1,LAT))                     &
            &   + (SLon**2 + SLat**2)*xyr_DSigTRC(i,j,k)/(r_m1(k)*xyz_H(i,j,k) + r_m2(k)*xyz_H(i,j,k+1))           &
            & ) 
    end do
    end do
    end do

    !$omp do
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyr_FSig(i,j,KS-1) = 0d0
       xyr_FSig(i,j,KE  ) = 0d0
    end do
    end do

    !$omp end parallel
    
  end subroutine calc_IsopycDiffFlux

  !-------------------------------------------------------------------
  
  subroutine calc_SkewFlux( &
       & xyz_FLon, xyz_FLat, xyr_FSig,            &  ! (out)
       & xyz_TRC, xyz_H, xyz_Z,                   &  ! (in)
       & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (in)
       & xyz_T,                                   &  ! (in)
       & xyza_HGradTRC, xyr_DSigTRC               &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(IA,JA,KA)
    real(DP), intent(out) :: xyz_FLat(IA,JA,KA)
    real(DP), intent(out) :: xyr_FSig(IA,JA,KA)
    real(DP), intent(in)  :: xyz_TRC(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)    
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyr_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyr_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_T(IA,JA,KA)
    real(DP), intent(in) :: xyza_HGradTRC(IA,JA,KA,2)
    real(DP), intent(in) :: xyr_DSigTRC(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xyza_Psi(IA,JA,KA,2)
    real(DP) :: xy_DzTRC(IA,JA)
    real(DP) :: xya_Psi(IA,JA,2)

    real(DP) :: DzTRC
    real(DP) :: PsiLat
    real(DP) :: PsiLon    
    real(DP) :: r_m1(KA)
    real(DP) :: r_m2(KA)

    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: l
    
    ! 実行文; Executable statement
    !

    do k=KS, KE-1
       r_m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       r_m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
    end do

    !$omp parallel private(PsiLon,PsiLat,DzTRC,i,j,jj,k,l)
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyza_Psi(i,j,k,LON) = - Kappa_GM*xyz_T(i,j,k)*xyz_SLat(i,j,k)
       xyza_Psi(i,j,k,LAT) =   Kappa_GM*xyz_T(i,j,k)*xyz_SLon(i,j,k)       
    end do
    end do
    end do
    !$omp do collapse(2)
    do l=1, 2
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyza_Psi(i,j,KS-1,l) = - xyza_Psi(i,j,KS,l)
       xyza_Psi(i,j,KE+1,l) = - xyza_Psi(i,j,KE,l)
    end do
    end do
    end do
    
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       DzTRC = 0.5d0*(xyr_DSigTRC(i,j,k-1) + xyr_DSigTRC(i,j,k))/xyz_H(i,j,k)
       xyz_FLon(i,j,k) = - xyza_Psi(i,j,k,LAT)*DzTRC
       xyz_FLat(i,j,k) =   xyza_Psi(i,j,k,LON)*DzTRC
    end do
    end do
    end do

    !$omp do collapse(2)
    do k=KS, KE-1
    do j=JS, JE
    do i=IS, IS+_IM_-1
       PsiLon = r_m1(k)*xyza_Psi(i,j,k,LON) + r_m2(k)*xyza_Psi(i,j,k+1,LON)
       PsiLat = r_m1(k)*xyza_Psi(i,j,k,LAT) + r_m2(k)*xyza_Psi(i,j,k+1,LAT)
       xyr_FSig(i,j,k) = ( &
            &   PsiLat*(r_m1(k)*xyza_HGradTRC(i,j,k,LON) + r_m2(k)*xyza_HGradTRC(i,j,k+1,LON)) &
            & - PsiLon*(r_m1(k)*xyza_HGradTRC(i,j,k,LAT) + r_m2(k)*xyza_HGradTRC(i,j,k+1,LAT)) &
            & )
    end do
    end do
    end do
    !$omp do
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyr_FSig(i,j,KS-1) = 0d0
       xyr_FSig(i,j,KE  ) = 0d0
    end do
    end do
    !$omp end parallel

  end subroutine calc_SkewFlux

  !-------------------------------------------------------------------
  
  subroutine calc_IsoNeutralSlope_new( xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (out)
       & xyz_PTemp, xyz_Salt, xyzaa_HGradTRC, xyra_DSigTRC, xyz_H, xyz_Z        &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(out) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(out) :: xyr_SLon(IA,JA,KA)
    real(DP), intent(out) :: xyr_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyzaa_HGradTRC(IA,JA,KA,2,2)
    real(DP), intent(in) :: xyra_DSigTRC(IA,JA,KA,2)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: fac
    real(DP) :: xyz_alpha(IA,JA,KA)
    real(DP) :: xyz_beta(IA,JA,KA)

    integer :: i
    integer :: j
    integer :: jj
    integer :: k

    real(DP) :: r_m1(KA)
    real(DP) :: r_m2(KA)

    real(DP), parameter :: EPS = 1d-10
    
    ! 実行文; Executable statement
    !

    
    !* Preparation
    !

    do k=KS, KE-1
       r_m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       r_m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
    end do

    !$omp parallel private(fac,i,j,jj,k)

    !$omp do
    do k=KS, KE
       call EOSDriver_alpha_beta( alpha=xyz_alpha(:,:,k), beta=xyz_beta(:,:,k),   & ! (out)
            & theta=xyz_PTemp(:,:,k), S=xyz_Salt(:,:,k),                          & ! (in)
            & p=-RefDens*Grav*xyz_Z(:,:,k) )                                        ! (in)
    end do
    
    !* Calculate the iso-netural slopes
    !

    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       fac = -  xyz_H(i,j,k) / ( EPS +                                                            &
            &   xyz_alpha(i,j,k) *0.5d0*(xyra_DSigTRC(i,j,k-1,PTEMP) + xyra_DSigTRC(i,j,k,PTEMP))  &
            & - xyz_beta (i,j,k) *0.5d0*(xyra_DSigTRC(i,j,k-1,SALT ) + xyra_DSigTRC(i,j,k,SALT ))  &
            & )
!!$
!!$       xy_Tmp(:,:) = -  xyz_H(:,:,k) / ( EPS +                                                     &
!!$            &   xyz_alpha(:,:,k) *0.5d0*(xyra_DSigTRC(:,:,k-1,PTEMP) + xyra_DSigTRC(:,:,k,PTEMP))  &
!!$            & - xyz_beta (:,:,k) *0.5d0*(xyra_DSigTRC(:,:,k-1,SALT ) + xyra_DSigTRC(:,:,k,SALT ))  &
!!$            & )
       
       xyz_SLon(i,j,k) =   fac*( &
            &   xyz_alpha(i,j,k)*xyzaa_HGradTRC(i,j,k,LON,PTEMP)   &
            & - xyz_beta (i,j,k)*xyzaa_HGradTRC(i,j,k,LON,SALT ) )

       xyz_SLat(i,j,k) =   fac*( &
            &   xyz_alpha(i,j,k)*xyzaa_HGradTRC(i,j,k,LAT,PTEMP)   &
            & - xyz_beta (i,j,k)*xyzaa_HGradTRC(i,j,k,LAT,SALT ) )
    end do
    end do
    end do
    
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1

!!$       xy_alpha(:,:) = m1*xyz_alpha(:,:,k) + m2*xyz_alpha(:,:,k+1)
!!$       xy_beta (:,:) = m1*xyz_beta (:,:,k) + m2*xyz_beta (:,:,k+1)
!!$       
!!$       !----------------------
!!$       
!!$       xy_Tmp(:,:) = - (m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1)) / (EPS +                                   &
!!$            &   xy_alpha(:,:)*xyra_DSigTRC(:,:,k,PTEMP) - xy_beta(:,:)*xyra_DSigTRC(:,:,k,SALT ) )

       xyr_SLon(i,j,k) = r_m1(k)*xyz_SLon(i,j,k) + r_m2(k)*xyz_SLon(i,j,k+1)
       xyr_SLat(i,j,k) = r_m1(k)*xyz_SLat(i,j,k) + r_m2(k)*xyz_SLat(i,j,k+1)
!!$       
!!$       xyr_SLon(:,:,k) = xy_Tmp(:,:)*( &
!!$            &   xy_alpha(:,:)*(m1*xyzaa_HGradTRC(:,:,k,LON,PTEMP) + m2*xyzaa_HGradTRC(:,:,k+1,LON,PTEMP))   &
!!$            & - xy_beta (:,:)*(m1*xyzaa_HGradTRC(:,:,k,LON,SALT ) + m2*xyzaa_HGradTRC(:,:,k+1,LON,SALT )) )
!!$
!!$       xyr_SLat(:,:,k) = xy_Tmp(:,:)*( &
!!$            &   xy_alpha(:,:)*(m1*xyzaa_HGradTRC(:,:,k,LAT,PTEMP) + m2*xyzaa_HGradTRC(:,:,k+1,LAT,PTEMP))   &
!!$            & - xy_beta (:,:)*(m1*xyzaa_HGradTRC(:,:,k,LAT,SALT ) + m2*xyzaa_HGradTRC(:,:,k+1,LAT,SALT )) )
    end do
    end do
    end do

    !$omp do
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyr_SLon(i,j,KS-1) = xyr_SLon(i,j,KS  )
       xyr_SLon(i,j,KE  ) = xyr_SLon(i,j,KE-1)
       xyr_SLat(i,j,KS-1) = xyr_SLat(i,j,KS  )
       xyr_SLat(i,j,KE  ) = xyr_SLat(i,j,KE-1)
    end do
    end do

    !$omp end parallel

  end subroutine calc_IsoNeutralSlope_new

  !-----------------------------------------------------------
  
  subroutine calc_BolusVelocity( &
       & xyz_BolusU, xyz_BolusV, xyz_BolusW,                & ! (out)
       & xyz_H, xyz_Z, xyz_SLon, xyz_SLat, xyz_T            & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_BolusU(IA,JA,KA)
    real(DP), intent(out) :: xyz_BolusV(IA,JA,KA)
    real(DP), intent(out) :: xyz_BolusW(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_T(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_PsiLon(IA,JA,KA)
    real(DP) :: xyz_PsiLat(IA,JA,KA)
    integer :: k

    ! 実行文; Executable statement
    !
!!$
!!$    !$omp parallel
!!$    !$omp do
!!$    do k=1, kMax-1
!!$       xyz_PsiLon(:,:,k) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
!!$       xyz_PsiLat(:,:,k) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
!!$    end do
!!$
!!$    !$omp workshare
!!$    xyz_PsiLon(:,:,KS) = 0d0; xyz_PsiLat(:,:,0) = 0d0    
!!$    xyz_PsiLon(:,:,k) = 0d0; xyz_PsiLat(:,:,kMax) = 0d0
!!$    !$omp end workshare
!!$    !$omp end parallel
!!$    
!!$    if(DFM08Flag) then
!!$       call TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
!!$            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Z )
!!$    end if
!!$
!!$    xyz_BolusU(:,:,:) = - xyz_DSig_xyz(xyz_PsiLat)/xyz_H !-xyz_Dz_xyz(xyz_PsiLat, xyz_H)
!!$    xyz_BolusV(:,:,:) =   xyz_DSig_xyz(xyz_PsiLon)/xyz_H ! xyz_Dz_xyz(xyz_PsiLon, xyz_H)
!!$                                                         ! xyz_Dz_xyz(xyz_PsiLon*xyz_CosLat, xyz_H)/xyz_CosLat
!!$                                                         ! = xyz_Dz_xyz(xyz_PsiLon, xyz_H)
!!$    !$omp parallel do
!!$    do k = 0, kMax
!!$       xyz_BolusW(:,:,k) = xy_w( &
!!$            & w_AlphaOptr_xy(xyz_PsiLat(:,:,k)*xy_CosLat, -xyz_PsiLon(:,:,k)*xy_CosLat) &
!!$            & )
!!$    end do

  end subroutine calc_BolusVelocity

  !--------------------------------------------------------------------------------------------------

  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !
    use LPhys_RediGMHelper_mod, only: &
         & init_taperingScheme, &
         & TAPERINGTYPE_DM95_NAME, PBLTAPERINGTYPE_LDD95_NAME
    
    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    use dc_calendar, only: DCCalConvertByUnit

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 

    ! IOSTAT of NAMELIST read
    character(TOKEN) :: EddyMixTypeName
    real(DP) :: DiffCoefRedi
    real(DP) :: DiffCoefGM
    character(TOKEN) :: InteriorTaperingName
    character(TOKEN) :: PBLTaperingName 
    real(DP) :: SlopeMax
    real(DP) :: Sd
    logical :: isUsedDFM08
    logical :: DiagOutputFlag
    character(STRING) :: msg

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /RediGM_nml/ &
         & DiffCoefRedi, DiffCoefGM, &
         & InteriorTaperingName, PBLTaperingName, &
         & SlopeMax, Sd, &
         & isUsedDFM08, &
         & DiagOutputFlag

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    DiffCoefRedi = Kappa_Redi
    DiffCoefGM = Kappa_GM

    InteriorTaperingName = TAPERINGTYPE_DM95_NAME
    PBLTaperingName = PBLTAPERINGTYPE_LDD95_NAME
    SlopeMax = 4d-3; Sd = 1d-3;
    isUsedDFM08 = .false.
    
    DiagOutputFlag = OutputFlag
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                        ! (in)
            & nml = RediGM_nml,              &  ! (out)
            & iostat = iostat_nml, iomsg=msg )  ! (out)
       close( unit_nml )
    end if


    ! Initialize tapering scheme
    call init_taperingScheme( &
         & interiorTaperingName, PBLTaperingName, &
         & SlopeMax, Sd, isUsedDFM08 )
         

    ! Set diffusivity and some flags   
    Kappa_Redi = DiffCoefRedi
    Kappa_GM = DiffCoefGM
    OutputFlag = DiagOutputFlag
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, ' DiffCoefRedi      = %f', d=(/ DiffCoefRedi  /) )
    call MessageNotify( 'M', module_name, ' DiffCoefGM        = %f', d=(/ DiffCoefGM  /) )
    call MessageNotify( 'M', module_name, ' InteriorTaperngName = %c', c1=InteriorTaperingName )
    call MessageNotify( 'M', module_name, ' PBLTaperngName      = %c', c1=PBLTaperingName )
    call MessageNotify( 'M', module_name, ' SlopeMax          = %f', d=(/ SlopeMax  /) )
    call MessageNotify( 'M', module_name, ' Sd                = %f', d=(/ Sd  /) )    
    call MessageNotify( 'M', module_name, ' isUsedDFM08       = %b', L=(/ isUsedDFM08 /) )    
    call MessageNotify( 'M', module_name, ' DiagOutputFlag    = %b', L=(/ DiagOutputFlag /) )
    
  end subroutine read_nmlData

  !-----------------------------------------------------------------------

  subroutine LPhys_RediGM_hspm_vfvm_Output()

    ! モジュール引用; Use statement
    !

    use DOGCM_Admin_TInteg_mod, only: &
         & TL_N => TIMELV_ID_N

    use DOGCM_Admin_Grid_mod, only: &
         & xyz_Z, xy_Topo
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC, xyza_H,     &
         & TRCID_SALT, TRCID_PTEMP

    use DOGCM_IO_History_mod, only:      &
       & DOGCM_IO_History_IsOutputTiming
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_PTemp(IA,JA,KA)
    real(DP) :: xyz_Salt(IA,JA,KA)
    real(DP) :: xyz_H(IA,JA,KA)

    real(DP) :: xyz_RefPress(IA,JA,KA)

    real(DP) :: xyz_DensPot(IA,JA,KA)
    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(IA,JA,KA)
    real(DP) :: xyz_SLat(IA,JA,KA)
    real(DP) :: xyr_SLon(IA,JA,KA)
    real(DP) :: xyr_SLat(IA,JA,KA)

    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(IA,JA,KA)

    !
    real(DP) :: xyzaa_HGradTRC(IA,JA,KA,2,2)
    real(DP) :: xyra_DSigTRC(IA,JA,KA,2)
    
    real(DP) :: xyz_FLon(IA,JA,KA)
    real(DP) :: xyz_FLat(IA,JA,KA)
    real(DP) :: xyz_FSig(IA,JA,KA)

    real(DP) :: xyz_BolusU(IA,JA,KA)
    real(DP) :: xyz_BolusV(IA,JA,KA)
    real(DP) :: xyz_BolusW(IA,JA,KA)

    real(DP) :: xyz_G(IA,JA,KA)
    real(DP) :: xyz_C(IA,JA,KA)
    real(DP) :: xy_BLD(IA,JA)
    
    integer :: k

    ! 実行文; Executable statement
    !

    if( .not. DOGCM_IO_History_isOutputTiming(CurrentTime) ) return
    if( .not. OutputFlag) return 

    xyz_PTemp(:,:,:) = xyzaa_TRC(:,:,:, TRCID_PTEMP,TL_N)
    xyz_Salt(:,:,:) = xyzaa_TRC(:,:,:, TRCID_SALT,TL_N)
    xyz_H(:,:,:) = xyza_H(:,:,:, TL_N)
    
    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot,            &  !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress )     !(in)

    ! Calculate the components of the isoneutral slope. 

    call calc_GradTRC( &
       & xyzaa_HGradTRC, xyra_DSigTRC,            & ! (out)
       & xyz_PTemp, xyz_Salt                      & ! (in)
       & )    
    
    call calc_IsoNeutralSlope_new( &
         & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,            & ! (out)
         & xyz_PTemp, xyz_Salt, xyzaa_HGradTRC, xyra_DSigTRC, & ! (in)
         & xyz_H, xyz_Z )                                       !(in)

    call prepare_SlopeTapering( xyz_T(:,:,KS:KE), xyz_SLon(:,:,KS:KE), xyz_SLat(:,:,KS:KE), & ! (inout)
         & xyz_DensPot(:,:,KS:KE), xyz_Z(:,:,KS:KE)                                         & ! (in)
         & )
!!$
!!$    call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig, &  ! (out)      
!!$         & xyz_PTemp, wa_xya(xyz_PTemp), xyz_H, xyz_Z,     &  ! (in)
!!$         & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
!!$         & )
!!$
!!$    call calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig,       &  ! (out)      
!!$         & xyz_PTemp, wa_xya(xyz_PTemp), xyz_H, xyz_Z,     &  ! (in)
!!$         & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
!!$         & )
    
    call HistoryPut('DensPot', xyz_DensPot(IS:IE,JS:JE,KS:KE), hst_SGSEddyMix)
!!$    call HistoryPut('MixLyrDepth', xy_BLD, hst_SGSEddyMix)

    call HistoryPut('SlopeLon', xyz_SLon(IS:IE,JS:JE,KS:KE), hst_SGSEddyMix)
    call HistoryPut('SlopeLat', xyz_T(IS:IE,JS:JE,KS:KE)*xyz_SLat(IS:IE,JS:JE,KS:KE), hst_SGSEddyMix)
    call HistoryPut('diffCoefEff', xyz_T(IS:IE,JS:JE,KS:KE), hst_SGSEddyMix)

!!$    !
!!$    call calc_BolusVelocity( &
!!$         & xyz_BolusU, xyz_BolusV, xyz_BolusW,                & ! (out)
!!$         & xyz_H, xyz_Z, xyz_SLon, xyz_SLat, xyz_T            & ! (in)
!!$         & )
!!$       
!!$    call HistoryPut('BolusU', xyz_BolusU(:,:,KS:KE), hst_SGSEddyMix)
!!$    call HistoryPut('BolusV', xyz_BolusV(:,:,KS:KE), hst_SGSEddyMix)
!!$    call HistoryPut('BolusW', xyz_BolusW(:,:,KS:KE), hst_SGSEddyMix)

  end subroutine LPhys_RediGM_hspm_vfvm_Output

  subroutine LPhys_RediGM_hspm_vfvm_PrepareOutput(   &
       & OriginTime, EndTime, Intrv, FilePrefix & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: OriginTime
    real(DP), intent(in) :: EndTime
    real(DP), intent(in) :: Intrv
    character(*), intent(in) :: FilePrefix

    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: lonName
    character(TOKEN) :: latName
    character(TOKEN) :: sigName
    character(TOKEN) :: timeName

    ! 実行文; Executable statement
    !

    if(.not. OutputFlag) return 

    
    !
    lonName = 'lon'; latName='lat'; sigName='sig'; timeName='t'

    call HistoryCreate( &                            ! ヒストリー作成
         & file=trim(FilePrefix) // 'SGSEddyMixOutput.nc', title='OGCM Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/lonName,latName,sigName,timeName/), dimsizes=(/iMax,jMax,KE-KS+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', 'sec         '/),  &
         & origin=real(OriginTime), interval=real(Intrv),  &        
         & history=hst_SGSEddyMix  )    

    call HistoryPut(lonName, x_CI(IS:IE), hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'topology', 'circular', hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'modulo', 360.0, hst_SGSEddyMix)
    call HistoryPut(latName, y_CJ(JS:JE), hst_SGSEddyMix)
    call HistoryPut(sigName, z_CK(KS:KE), hst_SGSEddyMix)

!!$    call regist_XYZTVariable('Etc1', 'tendency term1', 's-1')
!!$    call regist_XYZTVariable('Etc2', 'tendency term2', 's-1')
!!$    call regist_XYZTVariable('Etc3', 'tendency term3', 's-1')
!!$    call regist_XYZTVariable('Etc4', 'tendency term1', 's-1')
!!$    call regist_XYZTVariable('Etc5', 'tendency term2', 's-1')
!!$    call regist_XYZTVariable('Etc6', 'tendency term3', 's-1')

    call regist_XYZTVariable('DensPot', 'potential density', 'kg/m3')
    call regist_XYZTVariable('SlopeLon', 'the longitude component of slope of isoneutral', '1')
    call regist_XYZTVariable('SlopeLat', 'the latitude component of slope of isoneutral', '1')
    call regist_XYZTVariable('diffCoefEff', 'rescale factor of diffusivity', '1')
    call regist_XYTVariable('MixLyrDepth', 'depth of mixed layer used in linear tapering near the surface', 'm')
    
    call regist_XYZTVariable('BolusU', 'bolus velocity(longitude)', 'm/s')
    call regist_XYZTVariable('BolusV', 'bolus velocity(meridional)', 'm/s')
    call regist_XYZTVariable('BolusW', 'bolus velocity(vertical)', 'm/s')

  contains
    subroutine regist_XYZTVariable(varName, long_name, units)

      ! 宣言文; Declaration statement
      !      
      character(*), intent(in) :: varName, long_name, units

      ! 局所変数
      ! Local variables
      !      
      character(TOKEN) :: dims_XYZT(4)

      ! 実行文; Executable statement
      !
      dims_XYZT = (/ lonName, latName, sigName, timeName /)      
      call HistoryAddVariable(varName, dims_XYZT, &
           & long_name, units, xtype='float', history=hst_SGSEddyMix)

    end subroutine regist_XYZTVariable

    subroutine regist_XYTVariable(varName, long_name, units)

      ! 宣言文; Declaration statement
      !      
      character(*), intent(in) :: varName, long_name, units

      ! 局所変数
      ! Local variables
      !
      character(TOKEN) :: dims_XYT(3)

      ! 実行文; Executable statement
      !      
      dims_XYT = (/ lonName, latName, timeName /)      
      call HistoryAddVariable(varName, dims_XYT, &
           & long_name, units, xtype='float', history=hst_SGSEddyMix)

    end subroutine regist_XYTVariable

  end subroutine LPhys_RediGM_hspm_vfvm_PrepareOutput

end module LPhys_RediGM_hspm_vfvm_mod

