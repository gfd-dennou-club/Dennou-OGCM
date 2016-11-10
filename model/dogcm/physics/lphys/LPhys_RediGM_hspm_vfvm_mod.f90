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
       & IS, IE, JS, JE, KS, KE, KA,   &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon,             &
       & x_CI, y_CJ, z_CK,             &
       & z_CDK, z_RCDK, z_FDK, z_RFDK

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
    real(DP), intent(inout) :: xyz_PTemp_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyz_Salt_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyz_VDiffCoef(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xy_Topo(0:iMax-1,jMax)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DensPot(0:iMax-1,jMax,KA)
    real(DP) :: xyz_RefPress(0:iMax-1,jMax,KA)

    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP) :: xyr_SLon(0:iMax-1,jMax,KA)
    real(DP) :: xyr_SLat(0:iMax-1,jMax,KA)

    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(0:iMax-1,jMax,KA)

    !
    real(DP) :: xyzaa_HGradTRC(0:iMax-1,jMax,KA,2,2)
    real(DP) :: xyra_DSigTRC(0:iMax-1,jMax,KA,2)
    
    integer :: k
    
    ! 実行文; Executable statement
    !
    
    ! Calculate the potential density.
    !
    xyz_RefPress = 0d0
    call EOSDriver_Eval( xyz_DensPot, &                ! (out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress )         ! (in)

    
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

    call prepare_SlopeTapering( xyz_T(:,:,KS:KE), xyz_SLon(:,:,KS:KE), xyz_SLat(:,:,KS:KE), & ! (inout)
         & xyz_DensPot(:,:,KS:KE), xyz_Z(:,:,KS:KE)                                         & ! (in)
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

    subroutine append_Redi_GM_RHS( xyz_RHS, xyz_TRC, TRCID )

      ! 宣言文; Declaration statement
      !
      real(DP), intent(inout) :: xyz_RHS(0:iMax-1,jMax,KA)
      real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,KA)
      integer, intent(in) :: TRCID
      
      ! 局所変数
      ! Local variables
      !
      real(DP) :: xyz_FLonRedi(0:iMax-1,jMax,KA)
      real(DP) :: xyz_FLatRedi(0:iMax-1,jMax,KA)
      real(DP) :: xyr_FSigRedi(0:iMax-1,jMax,KA)

      real(DP) :: xyz_FLonGM(0:iMax-1,jMax,KA)
      real(DP) :: xyz_FLatGM(0:iMax-1,jMax,KA)
      real(DP) :: xyr_FSigGM(0:iMax-1,jMax,KA)

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

!      !$omp parallel do
      do k = KS, KE
         xyz_RHS(:,:,k) = xyz_RHS(:,:,k)  &
              & + xy_w( w_AlphaOptr_xy(                                        &
              &     (xyz_FLonRedi(:,:,k) + xyz_FLonGM(:,:,k))*xy_CosLat,       &
              &     (xyz_FLatRedi(:,:,k) + xyz_FLatGM(:,:,k))*xy_CosLat        &
              &    ) )                                                         &
              & + (   (xyr_FSigRedi(:,:,k-1) + xyr_FSigGM(:,:,k-1))            &
              &     - (xyr_FSigRedi(:,:,k  ) + xyr_FSigGM(:,:,k  ))            & 
              &   )*z_RCDK(k)/xyz_H(:,:,k)
      end do
      
    end subroutine append_Redi_GM_RHS

  end subroutine LPhys_RediGM_hspm_vfvm_AddMixingTerm

  subroutine calc_GradTRC( &
       & xyzaa_HGradTRC, xyra_DSigTRC,            & ! (out)
       & xyz_PTemp, xyz_Salt                      & ! (in)
       & )

    real(DP), intent(out) :: xyzaa_HGradTRC(0:iMax-1,jMax,KA,2,2)
    real(DP), intent(out) :: xyra_DSigTRC(0:iMax-1,jMax,KA,2)
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,KA)
    
    integer :: k

    real(DP) :: m1
    real(DP) :: m2

    real(DP) :: w_TRC(lMax)
    
 !   !$omp parallel
 !   !$omp do
    do k=KS, KE-1
       xyra_DSigTRC(:,:,k,PTEMP) = (xyz_PTemp(:,:,k) - xyz_PTemp(:,:,k+1))*z_RFDK(k)
       xyra_DSigTRC(:,:,k,SALT ) = (xyz_Salt (:,:,k) - xyz_Salt (:,:,k+1))*z_RFDK(k)       
    end do
 !   !$omp workshare
    xyra_DSigTRC(:,:,KS-1,:) = xyra_DSigTRC(:,:,KS  ,:)
    xyra_DSigTRC(:,:,KE  ,:) = xyra_DSigTRC(:,:,KE-1,:)
 !   !$omp end workshare
 !   !$omp end parallel

 !   !$omp parallel private(w_TRC)
!    !$omp do
    do k=KS, KE
       w_TRC(:) = w_xy(xyz_PTemp(:,:,k))
       xyzaa_HGradTRC(:,:,k,LON,PTEMP) = xy_GradLon_w(w_TRC)/RPlanet
       xyzaa_HGradTRC(:,:,k,LAT,PTEMP) = xy_GradLat_w(w_TRC)/RPlanet

       w_TRC(:) = w_xy(xyz_Salt(:,:,k))
       xyzaa_HGradTRC(:,:,k,LON,SALT) = xy_GradLon_w(w_TRC)/RPlanet
       xyzaa_HGradTRC(:,:,k,LAT,SALT) = xy_GradLat_w(w_TRC)/RPlanet
    end do
!    !$omp workshare
    xyzaa_HGradTRC(:,:,KS-1,:,:) = xyzaa_HGradTRC(:,:,KS,:,:)
    xyzaa_HGradTRC(:,:,KE+1,:,:) = xyzaa_HGradTRC(:,:,KE,:,:)
!    !$omp end workshare
!    !$omp end parallel
    
  end subroutine calc_GradTRC

  
  subroutine calc_IsopycDiffFlux( &
       & xyz_FLon, xyz_FLat, xyr_FSig,            &  ! (out)
       & xyz_TRC, xyz_H, xyz_Z,                   &  ! (in)
       & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (in)
       & xyz_T,                                   &  ! (in)
       & xyza_HGradTRC, xyr_DSigTRC               &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_FLat(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyr_FSig(0:iMax-1,jMax,KA)
    real(DP), intent(in)  :: xyz_TRC(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyza_HGradTRC(0:iMax-1,jMax,KA,2)
    real(DP), intent(in) :: xyr_DSigTRC(0:iMax-1,jMax,KA)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_C(0:iMax-1,jMax,KA)
    real(DP) :: xy_DzTRC(0:iMax-1,jMax)
    real(DP) :: m1
    real(DP) :: m2

    real(DP) :: xy_SLon(0:iMax-1,jMax)
    real(DP) :: xy_SLat(0:iMax-1,jMax)
    
    integer :: k

    ! 実行文; Executable statement
    !

    ! Calculate the components of diffusive flux along isopycnal surface.
    !

    !$omp parallel private(xy_DzTRC, m1, m2, xy_SLon, xy_SLat)

    !$omp do
    do k = KS, KE
       xy_DzTRC(:,:) = 0.5d0*(xyr_DSigTRC(:,:,k-1) + xyr_DSigTRC(:,:,k))/xyz_H(:,:,k)

       xyz_FLon(:,:,k) = Kappa_Redi*(                                                      &
            &   xyza_HGradTRC(:,:,k,LON) + xyz_T(:,:,k)*xyz_SLon(:,:,k)*xy_DzTRC(:,:) )  

       xyz_FLat(:,:,k) = Kappa_Redi*(                                                      &
            &   xyza_HGradTRC(:,:,k,LAT) + xyz_T(:,:,k)*xyz_SLat(:,:,k)*xy_DzTRC(:,:) )  
    end do

    !$omp do
    do k=KS, KE-1
       m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
       
!!$       xyr_FSig(:,:,k) = Kappa_Redi*(m1*xyz_T(:,:,k) + m2*xyz_T(:,:,k+1)) * ( &
!!$            &    xyr_SLon(:,:,k)*(m1*xyza_HGradTRC(:,:,k,LON) + m2*xyza_HGradTRC(:,:,k+1,LON))                        &
!!$            &  + xyr_SLat(:,:,k)*(m1*xyza_HGradTRC(:,:,k,LAT) + m2*xyza_HGradTRC(:,:,k+1,LAT))                        &
!!$            &  + (xyr_SLon(:,:,k)**2 + xyr_SLat(:,:,k)**2)*xyr_DSigTRC(:,:,k)/(m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1))   &
!!$            &  )

       xy_SLon(:,:) = m1*xyz_T(:,:,k)*xyz_SLon(:,:,k) + m2*xyz_T(:,:,k)*xyz_SLon(:,:,k)
       xy_SLat(:,:) = m1*xyz_T(:,:,k)*xyz_SLat(:,:,k) + m2*xyz_T(:,:,k)*xyz_SLat(:,:,k)

       xyr_FSig(:,:,k) = Kappa_Redi*( &
            &     xy_SLon(:,:)*(m1*xyza_HGradTRC(:,:,k,LON) + m2*xyza_HGradTRC(:,:,k+1,LON))                     &
            &   + xy_SLat(:,:)*(m1*xyza_HGradTRC(:,:,k,LAT) + m2*xyza_HGradTRC(:,:,k+1,LAT))                     &
            &   + (xy_SLon(:,:)**2 + xy_SLat(:,:)**2)*xyr_DSigTRC(:,:,k)/(m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1))   &
            & )       
    end do
    !$omp workshare
    xyr_FSig(:,:,KS-1) = 0d0
    xyr_FSig(:,:,KE  ) = 0d0
    !$omp end workshare
    !$omp end parallel
    
!!$    if(DFM08Flag) then
!!$       call TaperingDFM08_IDIFF(xyz_C, &
!!$            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, xyz_Z )
!!$
!!$       !$omp parallel do private(xy_GradLonTRC, xy_GradLatTRC)
!!$       do k = 0, kMax
!!$          xy_GradLonTRC(:,:) = xy_GradLon_w(wz_TRC(:,k))/RPlanet
!!$          xy_GradLatTRC(:,:) = xy_GradLat_w(wz_TRC(:,k))/RPlanet
!!$          
!!$          xyz_FLon(:,:,k) = xyz_C(:,:,k)*Kappa_Redi*xyz_T(:,:,k)*xy_GradLonTRC &
!!$               & + (1d0 - xyz_C(:,:,k))*xyz_FLon(:,:,k)
!!$          xyz_FLat(:,:,k) = xyz_C(:,:,k)*Kappa_Redi*xyz_T(:,:,k)*xy_GradLatTRC &
!!$               & + (1d0 - xyz_C(:,:,k))*xyz_FLat(:,:,k)
!!$       end do
!!$    end if
    
  end subroutine calc_IsopycDiffFlux

  subroutine calc_SkewFlux( &
       & xyz_FLon, xyz_FLat, xyr_FSig,            &  ! (out)
       & xyz_TRC, xyz_H, xyz_Z,                   &  ! (in)
       & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (in)
       & xyz_T,                                   &  ! (in)
       & xyza_HGradTRC, xyr_DSigTRC               &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_FLat(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyr_FSig(0:iMax-1,jMax,KA)
    real(DP), intent(in)  :: xyz_TRC(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyza_HGradTRC(0:iMax-1,jMax,KA,2)
    real(DP), intent(in) :: xyr_DSigTRC(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xyza_Psi(0:iMax-1,jMax,KA,2)
    real(DP) :: xy_DzTRC(0:iMax-1,jMax)
    real(DP) :: xya_Psi(0:iMax-1,jMax,2)
    real(DP) :: m1
    real(DP) :: m2
    
    integer :: k

    ! 実行文; Executable statement
    !

    !$omp parallel private(m1 ,m2, xya_Psi, xy_DzTRC)

    !$omp do
    do k=KS, KE
       xyza_Psi(:,:,k,LON) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
       xyza_Psi(:,:,k,LAT) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
    end do
    !$omp workshare
    xyza_Psi(:,:,KS-1,:) = - xyza_Psi(:,:,KS,:)
    xyza_Psi(:,:,KE+1,:) = - xyza_Psi(:,:,KE,:)
    !$omp end workshare


!!$    if(DFM08Flag) then
!!$       call TaperingDFM08_GM(xyza_Psi(:,:,:,LON), xyza_Psi(:,:,:,LAT), &
!!$            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Z )
!!$    end if

    !$omp do
    do k = KS, KE
       xy_DzTRC(:,:) = 0.5d0*(xyr_DSigTRC(:,:,k-1) + xyr_DSigTRC(:,:,k+1))/xyz_H(:,:,k)
       xyz_FLon(:,:,k) = - xyza_Psi(:,:,k,LAT)*xy_DzTRC(:,:)
       xyz_FLat(:,:,k) =   xyza_Psi(:,:,k,LON)*xy_DzTRC(:,:)       
    end do

    !$omp do
    do k = KS, KE-1
       m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))

       xya_Psi(:,:,:) = m1*xyza_Psi(:,:,k,:) + m2*xyza_Psi(:,:,k+1,:)
       xyr_FSig(:,:,k) = ( &
            &   xya_Psi(:,:,LAT)*(m1*xyza_HGradTRC(:,:,k,LON) + m2*xyza_HGradTRC(:,:,k,LON)) &
            & - xya_Psi(:,:,LON)*(m1*xyza_HGradTRC(:,:,k,LAT) + m2*xyza_HGradTRC(:,:,k,LAT)) &
            & )
    end do
    !$omp workshare
    xyr_FSig(:,:,KS-1) = 0d0
    xyr_FSig(:,:,KE  ) = 0d0
    !$omp end workshare

    !$omp end parallel

  end subroutine calc_SkewFlux

  subroutine calc_IsoNeutralSlope_new( xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,  &  ! (out)
       & xyz_PTemp, xyz_Salt, xyzaa_HGradTRC, xyra_DSigTRC, xyz_H, xyz_Z        &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyr_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyr_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyzaa_HGradTRC(0:iMax-1,jMax,KA,2,2)
    real(DP), intent(in) :: xyra_DSigTRC(0:iMax-1,jMax,KA,2)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_Tmp(0:iMax-1,jMax)    
    real(DP) :: xyz_alpha(0:iMax-1,jMax,KA)
    real(DP) :: xyz_beta(0:iMax-1,jMax,KA)
    real(DP) :: xy_alpha(0:iMax-1,jMax)
    real(DP) :: xy_beta(0:iMax-1,jMax)
    
    integer :: k

    real(DP) :: m1
    real(DP) :: m2

    real(DP), parameter :: EPS = 1d-10
    
    ! 実行文; Executable statement
    !

!!$    xyr_SLon = 1d100
!!$    xyr_SLat = 1d100
!!$    xyz_SLon = 1d100
!!$    xyz_SLat = 1d100
    
    !$omp parallel do private(xy_Tmp)
    do k = KS, KE
       call EOSDriver_alpha_beta( alpha=xyz_alpha(:,:,k), beta=xyz_beta(:,:,k),   & ! (out)
            & theta=xyz_PTemp(:,:,k), S=xyz_Salt(:,:,k),                          & ! (in)
            & p=-RefDens*Grav*xyz_Z(:,:,k) )                                        ! (in)

       !--------------------
       
       xy_Tmp(:,:) = -  xyz_H(:,:,k) / ( EPS +                                                     &
            &   xyz_alpha(:,:,k) *0.5d0*(xyra_DSigTRC(:,:,k-1,PTEMP) + xyra_DSigTRC(:,:,k,PTEMP))  &
            & - xyz_beta (:,:,k) *0.5d0*(xyra_DSigTRC(:,:,k-1,SALT ) + xyra_DSigTRC(:,:,k,SALT ))  &
            & )
       
       xyz_SLon(:,:,k) =   xy_Tmp(:,:)*( &
            &   xyz_alpha(:,:,k)*xyzaa_HGradTRC(:,:,k,LON,PTEMP)   &
            & - xyz_beta (:,:,k)*xyzaa_HGradTRC(:,:,k,LON,SALT ) )

       xyz_SLat(:,:,k) =   xy_Tmp(:,:)*( &
            &   xyz_alpha(:,:,k)*xyzaa_HGradTRC(:,:,k,LAT,PTEMP)   &
            & - xyz_beta (:,:,k)*xyzaa_HGradTRC(:,:,k,LAT,SALT ) )
    end do

    
    !$omp parallel private(xy_Tmp, m1, m2, xy_alpha, xy_beta)
    !$omp do
    do k = KS, KE-1
       m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))

!!$       xy_alpha(:,:) = m1*xyz_alpha(:,:,k) + m2*xyz_alpha(:,:,k+1)
!!$       xy_beta (:,:) = m1*xyz_beta (:,:,k) + m2*xyz_beta (:,:,k+1)
!!$       
!!$       !----------------------
!!$       
!!$       xy_Tmp(:,:) = - (m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1)) / (EPS +                                   &
!!$            &   xy_alpha(:,:)*xyra_DSigTRC(:,:,k,PTEMP) - xy_beta(:,:)*xyra_DSigTRC(:,:,k,SALT ) )

       xyr_SLon(:,:,k) = m1*xyz_SLon(:,:,k) + m2*xyz_SLon(:,:,k+1)
       xyr_SLat(:,:,k) = m1*xyz_SLat(:,:,k) + m2*xyz_SLat(:,:,k+1)
!!$       
!!$       xyr_SLon(:,:,k) = xy_Tmp(:,:)*( &
!!$            &   xy_alpha(:,:)*(m1*xyzaa_HGradTRC(:,:,k,LON,PTEMP) + m2*xyzaa_HGradTRC(:,:,k+1,LON,PTEMP))   &
!!$            & - xy_beta (:,:)*(m1*xyzaa_HGradTRC(:,:,k,LON,SALT ) + m2*xyzaa_HGradTRC(:,:,k+1,LON,SALT )) )
!!$
!!$       xyr_SLat(:,:,k) = xy_Tmp(:,:)*( &
!!$            &   xy_alpha(:,:)*(m1*xyzaa_HGradTRC(:,:,k,LAT,PTEMP) + m2*xyzaa_HGradTRC(:,:,k+1,LAT,PTEMP))   &
!!$            & - xy_beta (:,:)*(m1*xyzaa_HGradTRC(:,:,k,LAT,SALT ) + m2*xyzaa_HGradTRC(:,:,k+1,LAT,SALT )) )
    end do
    !$omp workshare
    xyr_SLon(:,:,KS-1) = xyr_SLon(:,:,KS  )
    xyr_SLon(:,:,KE  ) = xyr_SLon(:,:,KE-1)
    xyr_SLat(:,:,KS-1) = xyr_SLat(:,:,KS  )
    xyr_SLat(:,:,KE  ) = xyr_SLat(:,:,KE-1)
    !$omp end workshare
    !$omp end parallel

  end subroutine calc_IsoNeutralSlope_new

  !-----------------------------------------------------------
  
  subroutine calc_BolusVelocity( &
       & xyz_BolusU, xyz_BolusV, xyz_BolusW,                & ! (out)
       & xyz_H, xyz_Z, xyz_SLon, xyz_SLat, xyz_T            & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_BolusU(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_BolusV(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_BolusW(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_PsiLon(0:iMax-1,jMax,KA)
    real(DP) :: xyz_PsiLat(0:iMax-1,jMax,KA)
    integer :: k

    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp do
    do k=1, kMax-1
       xyz_PsiLon(:,:,k) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
       xyz_PsiLat(:,:,k) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
    end do

    !$omp workshare
    xyz_PsiLon(:,:,0) = 0d0; xyz_PsiLat(:,:,0) = 0d0    
    xyz_PsiLon(:,:,kMax) = 0d0; xyz_PsiLat(:,:,kMax) = 0d0
    !$omp end workshare
    !$omp end parallel
    
    if(DFM08Flag) then
       call TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Z )
    end if

    xyz_BolusU(:,:,:) = - xyz_DSig_xyz(xyz_PsiLat)/xyz_H !-xyz_Dz_xyz(xyz_PsiLat, xyz_H)
    xyz_BolusV(:,:,:) =   xyz_DSig_xyz(xyz_PsiLon)/xyz_H ! xyz_Dz_xyz(xyz_PsiLon, xyz_H)
                                                         ! xyz_Dz_xyz(xyz_PsiLon*xyz_CosLat, xyz_H)/xyz_CosLat
                                                         ! = xyz_Dz_xyz(xyz_PsiLon, xyz_H)
    !$omp parallel do
    do k = 0, kMax
       xyz_BolusW(:,:,k) = xy_w( &
            & w_AlphaOptr_xy(xyz_PsiLat(:,:,k)*xy_CosLat, -xyz_PsiLon(:,:,k)*xy_CosLat) &
            & )
    end do

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
    real(DP) :: xyz_PTemp(0:iMax-1,jMax,KA)
    real(DP) :: xyz_Salt(0:iMax-1,jMax,KA)
    real(DP) :: xyz_H(0:iMax-1,jMax,KA)

    real(DP) :: xyz_RefPress(0:iMax-1,jMax,KA)

    real(DP) :: xyz_DensPot(0:iMax-1,jMax,KA)
    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(0:iMax-1,jMax,KA)
    real(DP) :: xyz_SLat(0:iMax-1,jMax,KA)
    real(DP) :: xyr_SLon(0:iMax-1,jMax,KA)
    real(DP) :: xyr_SLat(0:iMax-1,jMax,KA)

    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(0:iMax-1,jMax,KA)

    !
    real(DP) :: xyzaa_HGradTRC(0:iMax-1,jMax,KA,2,2)
    real(DP) :: xyra_DSigTRC(0:iMax-1,jMax,KA,2)
    
    real(DP) :: xyz_FLon(0:iMax-1,jMax,KA)
    real(DP) :: xyz_FLat(0:iMax-1,jMax,KA)
    real(DP) :: xyz_FSig(0:iMax-1,jMax,KA)

    real(DP) :: xyz_BolusU(0:iMax-1,jMax,KA)
    real(DP) :: xyz_BolusV(0:iMax-1,jMax,KA)
    real(DP) :: xyz_BolusW(0:iMax-1,jMax,KA)

    real(DP) :: xyz_G(0:iMax-1,jMax,KA)
    real(DP) :: xyz_C(0:iMax-1,jMax,KA)
    real(DP) :: xy_BLD(0:iMax-1,jMax)
    
    integer :: k

    ! 実行文; Executable statement
    !

    if( .not. DOGCM_IO_History_isOutputTiming(CurrentTime) ) return
    if( .not. OutputFlag) return 

    xyz_PTemp(:,:,:) = xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TL_N)
    xyz_Salt(:,:,:) = xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TL_N)
    xyz_H(:,:,:) = xyza_H(IS:IE,JS:JE,KS:KE, TL_N)
    
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
    
    call HistoryPut('DensPot', xyz_DensPot(:,:,KS:KE), hst_SGSEddyMix)
!!$    call HistoryPut('MixLyrDepth', xy_BLD, hst_SGSEddyMix)

    call HistoryPut('SlopeLon', xyz_SLon(:,:,KS:KE), hst_SGSEddyMix)
    call HistoryPut('SlopeLat', xyz_T(:,:,KS:KE)*xyz_SLat(:,:,KS:KE), hst_SGSEddyMix)
    call HistoryPut('diffCoefEff', xyz_T(:,:,KS:KE), hst_SGSEddyMix)

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

