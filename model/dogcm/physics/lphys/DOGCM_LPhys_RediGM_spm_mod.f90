!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief This module to calculate tendecies of potential temperature and salinity due to isopycnal diffusion and GM advection.  
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_LPhys_RediGM_spm_mod 

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
       & PI, RPlanet, Omega

  use SpmlUtil_mod, only: &
       & w_AlphaOptr_xy,  &
       & xy_w, w_xy,      &
       & xya_wa, wa_xya,  &
       & xy_GradLon_w,    &
       & xy_GradLat_w,    &
       & xy_CosLat

  use EOSDriver_mod, only: &
       & EOSDriver_Eval

  use DOGCM_Admin_Grid_mod, only: &
       & IS, IE, JS, JE, KS, kE,       &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon,             &
       & x_CI, y_CJ, z_CK

  use DOGCM_LPhys_RediGMHelper_mod, only: &
       & DOGCM_LPhys_RediGMHelper_Init, DOGCM_LPhys_RediGMHelper_Final, &
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
  public :: DOGCM_LPhys_RediGM_spm_Init, DOGCM_LPhys_RediGM_spm_Final
  public :: DOGCM_LPhys_RediGM_spm_GetParameters
  public :: DOGCM_LPhys_RediGM_spm_Output, DOGCM_LPhys_RediGM_spm_PrepareOutput

  public :: DOGCM_LPhys_RediGM_spm_AddMixingTerm

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
  character(*), parameter:: module_name = 'DOGCM_LPhys_RediGM_spm_mod' !< Module Name

  real(DP) :: SGSEddyMixType
  real(DP) :: Kappa_Redi              !< Isopycnal diffusivity with Redi scheme [m^2/s]
  real(DP) :: Kappa_GM                !< Diffusivity with GM scheme             [m^2/s]

  type(gt_history), save :: hst_SGSEddyMix
  logical :: OutputFlag


  logical :: DFM08Flag
  
contains
  
  !> Initialize this module.
  !!
  !! Choose a scheme to parametrize sub-grid scale eddy mixing.
  !! If Redi scheme is used, set SGSEddyMixParamType to \link #DOGCM_LPhys_RediGM_spm_Redi \endlink. 
  !! If GM scheme is used, set SGSEddyMixParamType to \link #DOGCM_LPhys_RediGM_spm_GM \endlink. 
  !! KappaRedi and KappaGM are the parameters associated with the magnitude of diffusivity tensor in Redi or GM scheme,
  !! respectively(see a document of Dennou-OGCM for details).
  !! If you want, you can also set the parameters and output flag for analysis through namelist. In that case, specify
  !! the confignmlFileName.
  !!
  subroutine DOGCM_LPhys_RediGM_spm_Init( &
       & KappaRedi, KappaGM, isVarsOutput, confignmlFileName & ! (in)
       & )

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
    call DOGCM_LPhys_RediGMHelper_Init()
    
  end subroutine DOGCM_LPhys_RediGM_spm_Init

  
  !>
  !!
  !!
  subroutine DOGCM_LPhys_RediGM_spm_Final()

    ! 実行文; Executable statements
    !

    call DOGCM_LPhys_RediGMHelper_Final()
    
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)
    
  end subroutine DOGCM_LPhys_RediGM_spm_Final

  !-------------------------------------------------------------------
  
  !>
  !!
  !!
  subroutine DOGCM_LPhys_RediGM_spm_GetParameters( &
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

  end subroutine DOGCM_LPhys_RediGM_spm_GetParameters
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_LPhys_RediGM_spm_AddMixingTerm( &
       & xyz_PTemp_RHS, xyz_Salt_RHS,                                &  ! (inout)
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo                  &  ! (in)
       & )

    ! モジュール引用; Use statements
    !    
    use DOGCM_Admin_Constants_mod
    use DOGCM_Admin_TInteg_mod, only: &
         & CurrentTime

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_PTemp_RHS(0:iMax-1,jMax, 0:kMax)
    real(DP), intent(inout) :: xyz_Salt_RHS(0:iMax-1,jMax, 0:kMax)
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_Topo(0:iMax-1,jMax)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DensPot(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_RefPress(0:iMax-1,jMax,0:kMax)

    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(0:iMax-1,jMax,0:kMax)

    integer :: k

    ! 実行文; Executable statement
    !
    
    ! Calculate the potential density.
    !
    xyz_RefPress = 0d0
    call EOSDriver_Eval( xyz_DensPot, &                !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress )         !(in)

    
    ! Calculate the components of the isoneutral slope. 
    !
    
    call calc_IsoNeutralSlope( xyz_SLon, xyz_SLat, &   !(out)
         & xyz_DensPot, xyz_H )                        !(in)

    call prepare_SlopeTapering( xyz_T, xyz_SLon, xyz_SLat, & ! (inout)
         & xyz_DensPot, xyz_Z                              & ! (in)
         & )

!    call check_StaticStability(xyz_T, xyz_DensPot)

    
    ! Calculate the tendency due to Gent-McWilliams and Redi flux
    !
    
    call append_Redi_GM_RHS( xyz_PTemp_RHS,    & ! (inout)
         & xyz_PTemp )                           ! (in)    

    call append_Redi_GM_RHS( xyz_Salt_RHS,     & ! (inout)
         & xyz_Salt )                            ! (in)    

  contains

    subroutine append_Redi_GM_RHS( xyz_RHS, xyz_TRC )

      ! 宣言文; Declaration statement
      !
      real(DP), intent(inout) :: xyz_RHS(0:iMax-1,jMax,0:kMax)
      real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP) :: wz_TRC(lMax,0:kMax)
      
      real(DP) :: xyz_FLonRedi(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_FLatRedi(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_FSigRedi(0:iMax-1,jMax,0:kMax)

      real(DP) :: xyz_FLonGM(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_FLatGM(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_FSigGM(0:iMax-1,jMax,0:kMax)

      integer :: k

      real(DP) :: xyz_DzFSig(0:iMax-1,jMax,0:kMax)
      
      ! 実行文; Executable statement
      !

      wz_TRC(:,:) = wa_xya(xyz_TRC)
      
      call calc_IsopycDiffFlux( &
           & xyz_FLonRedi, xyz_FLatRedi, xyz_FSigRedi,       &  ! (out)      
           & xyz_TRC, wz_TRC, xyz_H, xyz_Z,                  &  ! (in)
           & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
           & )

      call calc_SkewFlux( &
           & xyz_FLonGM, xyz_FLatGM, xyz_FSigGM,             &  ! (out)      
           & xyz_TRC, wz_TRC, xyz_H, xyz_Z,                  &  ! (in)
           & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
           & )
      
      xyz_DzFSig(:,:,:) = xyz_Dz_xyz(xyz_FSigRedi + xyz_FSigGM, xyz_H)
      !$omp parallel do
      do k = 0, kMax
         xyz_RHS(:,:,k) = xyz_RHS(:,:,k)  &
              & + xy_w( w_AlphaOptr_xy(                                  &
              &     (xyz_FLonRedi(:,:,k) + xyz_FLonGM(:,:,k))*xy_CosLat, &
              &     (xyz_FLatRedi(:,:,k) + xyz_FLatGM(:,:,k))*xy_CosLat  &
              &    ) )                                                   &
              & + xyz_DzFSig(:,:,k)
      end do

    end subroutine append_Redi_GM_RHS

  end subroutine DOGCM_LPhys_RediGM_spm_AddMixingTerm
  
  subroutine calc_IsopycDiffFlux( &
       & xyz_FLon, xyz_FLat, xyz_FSig,            &  ! (out)
       & xyz_TRC, wz_TRC, xyz_H, xyz_Z,           &  ! (in)
       & xyz_SLon, xyz_SLat, xyz_T                &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_FLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_FSig(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in)  :: xyz_TRC(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: wz_TRC(lMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,0:kMax)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_GradLonTRC(0:iMax-1,jMax)
    real(DP) :: xy_GradLatTRC(0:iMax-1,jMax)
    real(DP) :: xyz_DzTRC(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_C(0:iMax-1,jMax,0:kMax)

    integer :: k

    ! 実行文; Executable statement
    !

    ! Calculate the components of diffusive flux along isopycnal surface.
    !
    
    xyz_DzTRC(:,:,:) = xyz_Dz_xyz(xyz_TRC, xyz_H)
    !$omp parallel do private(xy_GradLonTRC, xy_GradLatTRC)
    do k = 1, kMax-1
       xy_GradLonTRC(:,:) = xy_GradLon_w(wz_TRC(:,k))/RPlanet
       xy_GradLatTRC(:,:) = xy_GradLat_w(wz_TRC(:,k))/RPlanet
       
       xyz_FLon(:,:,k) = Kappa_Redi*( &
            & xy_GradLonTRC + xyz_T(:,:,k)*xyz_SLon(:,:,k)*xyz_DzTRC(:,:,k) )

       xyz_FLat(:,:,k) = Kappa_Redi*( &
            & xy_GradLatTRC + xyz_T(:,:,k)*xyz_SLat(:,:,k)*xyz_DzTRC(:,:,k) )

       xyz_FSig(:,:,k) = Kappa_Redi*xyz_T(:,:,k)*( &
            &    xyz_SLon(:,:,k)*xy_GradLonTRC + xyz_SLat(:,:,k)*xy_GradLatTRC  &
            & + (xyz_SLon(:,:,k)**2 + xyz_SLat(:,:,k)**2)*xyz_DzTRC(:,:,k)      &
            & )
    end do

    ! Set zero for the vertical component of isopycnal diffusive flux at sea surface and bottom.
    !
    
    !$omp parallel
    !$omp workshare
    xyz_FSig(:,:,0) = 0d0
    xyz_FSig(:,:,kMax) = 0d0
    !$omp end workshare
    !$omp end parallel
    
    if(DFM08Flag) then
       call TaperingDFM08_IDIFF(xyz_C, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, xyz_Z )

       !$omp parallel do private(xy_GradLonTRC, xy_GradLatTRC)
       do k = 0, kMax
          xy_GradLonTRC(:,:) = xy_GradLon_w(wz_TRC(:,k))/RPlanet
          xy_GradLatTRC(:,:) = xy_GradLat_w(wz_TRC(:,k))/RPlanet
          
          xyz_FLon(:,:,k) = xyz_C(:,:,k)*Kappa_Redi*xyz_T(:,:,k)*xy_GradLonTRC &
               & + (1d0 - xyz_C(:,:,k))*xyz_FLon(:,:,k)
          xyz_FLat(:,:,k) = xyz_C(:,:,k)*Kappa_Redi*xyz_T(:,:,k)*xy_GradLatTRC &
               & + (1d0 - xyz_C(:,:,k))*xyz_FLat(:,:,k)
       end do
    end if
    
  end subroutine calc_IsopycDiffFlux

  subroutine calc_SkewFlux( &
       & xyz_FLon, xyz_FLat, xyz_FSig,            &  ! (out)
       & xyz_TRC, wz_TRC, xyz_H, xyz_Z,           &  ! (in)
       & xyz_SLon, xyz_SLat, xyz_T                &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_FLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_FLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_FSig(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in)  :: xyz_TRC(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: wz_TRC(lMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DzTRC(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PsiLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PsiLat(0:iMax-1,jMax,0:kMax)

    integer :: k

    ! 実行文; Executable statement
    !

    xyz_DzTRC(:,:,:) = xyz_Dz_xyz(xyz_TRC, xyz_H)

    !$omp parallel 
    !$omp do
    do k=1, kMax-1
       xyz_PsiLon(:,:,k) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
       xyz_PsiLat(:,:,k) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
    end do

    !$omp workshare
    xyz_PsiLat(:,:,0) = 0d0; xyz_PsiLon(:,:,0) = 0d0
    xyz_PsiLat(:,:,kMax) = 0d0; xyz_PsiLon(:,:,kMax) = 0d0
    !$omp end workshare
    !$omp end parallel

    if(DFM08Flag) then
       call TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Z )
    end if
    
    !$omp parallel do    
    do k=0, kMax
       xyz_FLon(:,:,k) = - xyz_PsiLat(:,:,k)*xyz_DzTRC(:,:,k)
       xyz_FLat(:,:,k) =   xyz_PsiLon(:,:,k)*xyz_DzTRC(:,:,k)
       xyz_FSig(:,:,k) = ( &
            &   xyz_PsiLat(:,:,k)*xy_GradLon_w(wz_TRC(:,k)) &
            & - xyz_PsiLon(:,:,k)*xy_GradLat_w(wz_TRC(:,k)) &
            & )/RPlanet
    end do

  end subroutine calc_SkewFlux
  
  subroutine calc_IsoNeutralSlope( xyz_SLon, xyz_SLat,          &  ! (out)
       & xyz_DensPot, xyz_H                                     &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_DensPot(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DzDensPot(0:iMax-1,jMax,0:kMax)
    real(DP) :: w_DensPot(lMax)
    integer :: i, j, k
    real(DP), parameter :: EPS = 1d-10

    ! 実行文; Executable statement
    !

    xyz_DzDensPot(:,:,:) = xyz_Dz_xyz(xyz_DensPot, xyz_H)

    !$omp parallel do private(w_DensPot)
    do k=0, kMax
       w_DensPot(:) = w_xy(xyz_DensPot(:,:,k))
       xyz_SLon(:,:,k) = - xy_GradLon_w(w_DensPot)/(xyz_DzDensPot(:,:,k) - EPS)/RPlanet
       xyz_SLat(:,:,k) = - xy_GradLat_w(w_DensPot)/(xyz_DzDensPot(:,:,k) - EPS)/RPlanet
    end do
    
  end subroutine calc_IsoNeutralSlope


  
!!$  subroutine check_StaticStability(xyz_KappRho, xyz_DensPot)
!!$
!!$    ! 宣言文; Declaration statement
!!$    !
!!$    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: xyz_KappRho
!!$    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot
!!$
!!$    ! 局所変数
!!$    ! Local variables
!!$    !
!!$    integer :: k
!!$
!!$    ! 実行文; Executable statement
!!$    !
!!$
!!$    return
!!$    
!!$!    !$omp parallel do
!!$    do k=0, kMax-1
!!$       where(xyz_DensPot(:,:,k) > xyz_DensPot(:,:,k+1))
!!$          xyz_KappRho(:,:,k) = 0d0
!!$          xyz_KappRho(:,:,k+1) = 0d0
!!$       end where
!!$    end do
!!$
!!$  end subroutine check_StaticStability
  
  subroutine calc_BolusVelocity( &
       & xyz_BolusU, xyz_BolusV, xyz_BolusW,                & ! (out)
       & xyz_H, xyz_Z, xyz_SLon, xyz_SLat, xyz_T            & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_BolusU(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_BolusV(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_BolusW(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_T(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_PsiLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PsiLat(0:iMax-1,jMax,0:kMax)
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

    xyz_BolusU(:,:,:) = - xyz_Dz_xyz(xyz_PsiLat, xyz_H)
    xyz_BolusV(:,:,:) =   xyz_Dz_xyz(xyz_PsiLon, xyz_H) ! xyz_Dz_xyz(xyz_PsiLon*xyz_CosLat, xyz_H)/xyz_CosLat
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
    use DOGCM_LPhys_RediGMHelper_mod, only: &
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

  subroutine DOGCM_LPhys_RediGM_spm_Output()

    ! モジュール引用; Use statement
    !

    use DOGCM_Admin_TInteg_mod, only: &
         & TL_N => TIMELV_ID_N

    use DOGCM_Admin_Grid_mod, only: &
         & xyz_Z, xy_Topo
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC, xyza_H,     &
         & TRCID_SALT, TRCID_PTEMP

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_H(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_RefPress(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_DensPot(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SLat(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_FLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_FLat(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_FSig(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_BolusU(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_BolusV(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_BolusW(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_T(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_G(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_C(0:iMax-1,jMax,0:kMax)
    real(DP) :: xy_BLD(0:iMax-1,jMax)
    
    integer :: k

    ! 実行文; Executable statement
    !

    if(.not. OutputFlag) return 

    xyz_PTemp(:,:,:) = xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TL_N)
    xyz_Salt(:,:,:) = xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TL_N)
    xyz_H(:,:,:) = xyza_H(IS:IE,JS:JE,KS:KE, TL_N)
    
    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot,            &  !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress )     !(in)

    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_H                )       !(in)

    !
    !
    call prepare_SlopeTapering(xyz_T, xyz_SLon, xyz_SLat,  & ! (inout)
         & xyz_DensPot, xyz_Z,                             & ! (in)
         & xy_BLD                                          & ! (out)
         & )
    
!    call check_StaticStability(xyz_T, xyz_DensPot)    

    call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig, &  ! (out)      
         & xyz_PTemp, wa_xya(xyz_PTemp), xyz_H, xyz_Z,     &  ! (in)
         & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
         & )

    call calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig,       &  ! (out)      
         & xyz_PTemp, wa_xya(xyz_PTemp), xyz_H, xyz_Z,     &  ! (in)
         & xyz_SLon, xyz_SLat, xyz_T                       &  ! (in)
         & )
    
!!$    call HistoryPut('DensPot', xyz_DensPot, hst_SGSEddyMix)
!!$    call HistoryPut('MixLyrDepth', xy_BLD, hst_SGSEddyMix)

    call HistoryPut('SlopeLon', xyz_T*xyz_SLon, hst_SGSEddyMix)
    call HistoryPut('SlopeLat', xyz_T*xyz_SLat, hst_SGSEddyMix)
    call HistoryPut('diffCoefEff', xyz_T, hst_SGSEddyMix)

    !
    call calc_BolusVelocity( &
         & xyz_BolusU, xyz_BolusV, xyz_BolusW,                & ! (out)
         & xyz_H, xyz_Z, xyz_SLon, xyz_SLat, xyz_T            & ! (in)
         & )
       
    call HistoryPut('BolusU', xyz_BolusU, hst_SGSEddyMix)
    call HistoryPut('BolusV', xyz_BolusV, hst_SGSEddyMix)
    call HistoryPut('BolusW', xyz_BolusW, hst_SGSEddyMix)

  end subroutine DOGCM_LPhys_RediGM_spm_Output

  subroutine DOGCM_LPhys_RediGM_spm_PrepareOutput(   &
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
         & dims=(/lonName,latName,sigName,timeName/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', 'sec         '/),  &
         & origin=real(OriginTime), interval=real(Intrv),  &        
         & history=hst_SGSEddyMix  )    

    call HistoryPut(lonName, x_CI(IS:IE)*180d0/PI, hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'topology', 'circular', hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'modulo', 360.0, hst_SGSEddyMix)
    call HistoryPut(latName, y_CJ(JS:JE)*180d0/PI, hst_SGSEddyMix)
    call HistoryPut(sigName, z_CK(KS:KE), hst_SGSEddyMix)

!!$    call regist_XYZTVariable('Etc1', 'tendency term1', 's-1')
!!$    call regist_XYZTVariable('Etc2', 'tendency term2', 's-1')
!!$    call regist_XYZTVariable('Etc3', 'tendency term3', 's-1')
!!$    call regist_XYZTVariable('Etc4', 'tendency term1', 's-1')
!!$    call regist_XYZTVariable('Etc5', 'tendency term2', 's-1')
!!$    call regist_XYZTVariable('Etc6', 'tendency term3', 's-1')

!!$    call regist_XYZTVariable('DensPot', 'potential density', 'kg/m3')
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

  end subroutine DOGCM_LPhys_RediGM_spm_PrepareOutput

end module DOGCM_LPhys_RediGM_spm_mod

