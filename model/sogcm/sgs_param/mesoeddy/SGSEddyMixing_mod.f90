!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief This module to calculate tendecies of potential temperature and salinity due to isopycnal diffusion and GM advection.  
!! 
!! @author Yuta Kawai
!!
!!
module SGSEddyMixing_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  !* Dennou-OGCM

  use Constants_mod, only: &
       & PI, RPlanet, Omega

  use SpmlUtil_mod

  use DiagnoseUtil_mod

  use EOSDriver_mod, only: &
       & EOSDriver_Eval

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon

  use VariableSet_mod

  use SGSEddyMixingHelper_mod, only: &
       & SGSEddyMixingHelper_Init, SGSEddyMixingHelper_Final, &
       & xyz_Dz_xyz, &
       & prepare_SlopeTapering, &
       & DFM08Info, prepare_DFM08Info, TaperingDFM08_GM, TaperingDFM08_IDIFF
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SGSEddyMixing_Init, SGSEddyMixing_Final
  public :: SGSEddyMixing_GetParameters
  public :: SGSEddyMixing_Output, SGSEddyMixing_PrepareOutput
  public :: SGSEddyMixing_AddMixingTerm, SGSEddyMixing_AddMixingTerm_wz
  public :: calc_IsoNeutralSlope
  public :: calc_IsopycDiffFlux
  public :: calc_BolusVelocity
  
  ! 公開変数
  ! Public variable
  !
  character(*), parameter, public :: SGSEddyMixing_Redi_NAME = 'Redi'
  integer, parameter, public :: SGSEddyMixing_Redi = 1
  character(*), parameter, public :: SGSEddyMixing_GM_NAME = 'GM90'
  integer, parameter, public :: SGSEddyMixing_GM90 = 2

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSEddyMixing_mod' !< Module Name

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
  !! If Redi scheme is used, set SGSEddyMixParamType to \link #SGSEddyMixing_Redi \endlink. 
  !! If GM scheme is used, set SGSEddyMixParamType to \link #SGSEddyMixing_GM \endlink. 
  !! KappaRedi and KappaGM are the parameters associated with the magnitude of diffusivity tensor in Redi or GM scheme,
  !! respectively(see a document of Dennou-OGCM for details).
  !! If you want, you can also set the parameters and output flag for analysis through namelist. In that case, specify
  !! the confignmlFileName.
  !!
  subroutine SGSEddyMixing_Init(SGSEddyMixParamType, KappaRedi, KappaGM, isVarsOutput, confignmlFileName)

    ! 宣言文; Declaration statement
    !
    integer, intent(in), optional :: SGSEddyMixParamType   !< Type of parameterization for sub-grid scale eddy mixing
    real(DP), intent(in), optional  :: KappaRedi           !< Isopycnal diffusivity, parameter associated with the magnitude of diffusivity tensor with Redi scheme [m^2/s]
    real(DP), intent(in), optional  :: KappaGM             !< Parameter associated with the magnitude of diffusivity tensor with GM scheme [m^2/s]
    logical, intent(in), optional :: isVarsOutput
    character(*), intent(in), optional :: confignmlFileName   !< Namelist name

    ! 実行文; Executable statements
    !

    ! Set default values.
    SGSEddyMixType = SGSEddyMixing_Redi
    Kappa_Redi = 1000d0
    Kappa_GM   = 1000d0
    OutputFlag = .false.

    !

    if(present(SGSEddyMixParamType)) SGSEddyMixType = SGSEddyMixParamType
    if(present(KappaRedi)) Kappa_Redi = KappaRedi
    if(present(KappaGM)) Kappa_GM = KappaGM
    if(present(isVarsOutput)) OutputFlag = isVarsOutput

    if(present(configNmlFileName)) call read_nmlData(configNmlFileName)

    call SGSEddyMixingHelper_Init()
    
  end subroutine SGSEddyMixing_Init

  
  !>
  !!
  !!
  subroutine SGSEddyMixing_GetParameters( &
       & KappaRedi, KappaGM )
    real(DP), intent(out), optional :: KappaRedi, KappaGM

    if(present(KappaRedi)) KappaRedi = Kappa_Redi
    if(present(KappaGM)) KappaGM = Kappa_GM

  end subroutine SGSEddyMixing_GetParameters
  
  !>
  !!
  !!
  subroutine SGSEddyMixing_Final()

    ! 実行文; Executable statements
    !

    call SGSEddyMixingHelper_Final()
    
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)

    
  end subroutine SGSEddyMixing_Final

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixing_AddMixingTerm_wz( wz_PTempRHS, wz_SaltRHS, &  ! (inout)
       & wz_PTemp, wz_Salt, xy_totDepth                            &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_PTempRHS, wz_SaltRHS
    real(DP), dimension(lMax,0:kMax), intent(in)  :: wz_PTemp, wz_Salt
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_totDepth

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTempRHS, xyz_SaltRHS    

    xyz_PTempRHS = 0d0; xyz_SaltRHS = 0d0
    call SGSEddyMixing_AddMixingTerm(xyz_PTempRHS, xyz_SaltRHS, &
         & xyz_wz(wz_PTemp), xyz_wz(wz_Salt), xy_totDepth)

    wz_PTempRHS(:,:) = wz_PTempRHS + wz_xyz(xyz_PTempRHS)
    wz_SaltRHS(:,:) = wz_SaltRHS + wz_xyz(xyz_SaltRHS)
    
  end subroutine SGSEddyMixing_AddMixingTerm_wz

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixing_AddMixingTerm( xyz_PTempRHS, xyz_SaltRHS,   &  ! (inout)
       & xyz_PTemp, xyz_Salt, xy_totDepth                              &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: xyz_PTempRHS, xyz_SaltRHS
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in)  :: xyz_PTemp, xyz_Salt
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_totDepth

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DensPot, xyz_RefPress, xyz_CosLat, xyz_totDepth

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, &  !< The longitude componet of th slope of isopycnal surface
         & xyz_SLat     !< The latitude componet of th slope of isopycnal surface

    real(DP), dimension(lMax,0:kMax) :: wz_PTemp, wz_Salt
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_KappRho, xyz_Depth
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_T    
    integer :: k

    ! 実行文; Executable statement
    !

    xyz_CosLat(:,:,:) = spread(cos(xyz_Lat(:,:,0)),3,kMax+1)

    !$omp parallel do
    do k=0, kMax
       xyz_totDepth(:,:,k) = xy_totDepth
       xyz_Depth(:,:,k) = xy_totDepth*g_Sig(k)
       wz_PTemp(:,k) = w_xy(xyz_PTemp(:,:,k))
       wz_Salt(:,k) = w_xy(xyz_Salt(:,:,k))
    end do
    
    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot, &         !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress ) !(in)

    
    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_totDepth )               !(in)

    !
    call prepare_SlopeTapering(xyz_T, xyz_SLon, xyz_SLat, & ! (inout)
         & xyz_DensPot, xyz_Depth                         & ! (in)
         & )

!    call check_StaticStability(xyz_T, xyz_DensPot)

    !
    xyz_PTempRHS(:,:,:) = xyz_PTempRHS + calc_IsopycDiffTerm(wz_PTemp, xyz_PTemp)
    xyz_SaltRHS(:,:,:)  = xyz_SaltRHS  + calc_IsopycDiffTerm(wz_Salt, xyz_Salt)

    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       xyz_PTempRHS(:,:,:) = xyz_PTempRHS + calc_EddyInducedVelAdvTerm(wz_PTemp, xyz_PTemp)
       xyz_SaltRHS(:,:,:)  = xyz_SaltRHS + calc_EddyInducedVelAdvTerm(wz_Salt, xyz_Salt)
    end if

  contains
    function calc_IsopycDiffTerm(wz_Tracer, xyz_Tracer) result(xyz_DiffTerm)

      ! 宣言文; Declaration statement
      !
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
      real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_DiffTerm(0:iMax-1,jMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_FLon, xyz_FLat, xyz_FSig

      ! 実行文; Executable statement
      !

      call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,                &  !(out)      
           & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth &  !(in)
           & )

      xyz_DiffTerm(:,:,:) = &
           &    xyz_wz(wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat)) &
           & +  xyz_Dz_xyz(xyz_FSig)

    end function calc_IsopycDiffTerm

    function calc_EddyInducedVelAdvTerm(wz_Tracer, xyz_Tracer) result(xyz_AdvTerm)

      ! 宣言文; Declaration statement
      !
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
      real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_AdvTerm(0:iMax-1,jMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
           & xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, &
           & xyz_DzT

      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_FLon, xyz_FLat, xyz_FSig

      
      ! 実行文; Executable statement
      !

      call calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig,                      & ! (out)
           & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth & ! (in)
           & )

!!$      call calc_BolusVelocity( &
!!$           & xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, & ! (out)
!!$           & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T      & ! (in)
!!$           & )
!!$      xyz_EddInducedW(:,:,0) = 0d0
!!$      xyz_EddInducedW(:,:,kMax) = 0d0
!!$      xyz_FLon = -xyz_EddInducedU*xyz_Tracer
!!$      xyz_FLat = -xyz_EddInducedV*xyz_Tracer
!!$      xyz_FSig = -xyz_EddInducedW*xyz_Tracer
!!$
!!$      wz_AdvTerm(:,:) = &
!!$           &   wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat) &
!!$           & + wz_xyz( &
!!$           &    xyz_Tracer*xyz_wz(wz_AlphaOptr_xyz(xyz_EddInducedU*xyz_CosLat, xyz_EddInducedV*xyz_CosLat)) &
!!$           &  - xyz_EddInducedW*xyz_Dz_xyz(xyz_Tracer) &
!!$           &   )
      xyz_AdvTerm(:,:,:) = &
           &   xyz_wz(wz_AlphaOptr_xyz( xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat )) &
           & + xyz_Dz_xyz(xyz_FSig)

    end function calc_EddyInducedVelAdvTerm

  end subroutine SGSEddyMixing_AddMixingTerm

  subroutine calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,          &  ! (out)
       & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth    &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: &
         & xyz_FLon, xyz_FLat, xyz_FSig
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth
    real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_GradLonT, xyz_GradLatT, xyz_DzT, xyz_C, xyz_AI
    integer :: k

    ! 実行文; Executable statement
    !

    xyz_GradLonT(:,:,:) = xyz_GradLon_wz(wz_Tracer)
    xyz_GradLatT(:,:,:) = xyz_GradLat_wz(wz_Tracer)
    xyz_DzT(:,:,:) = xyz_Dz_xyz(xyz_Tracer)
    
    xyz_AI(:,:,:) = Kappa_Redi*xyz_T
    !$omp parallel workshare
    xyz_FLon(:,:,:) = Kappa_Redi*(xyz_GradLonT + xyz_T*xyz_SLon*xyz_DzT)
    xyz_FLat(:,:,:) = Kappa_Redi*(xyz_GradLatT + xyz_T*xyz_SLat*xyz_DzT)
!!$    xyz_FLon(:,:,:) = xyz_AI*(xyz_GradLonT + xyz_SLon*xyz_DzT)
!!$    xyz_FLat(:,:,:) = xyz_AI*(xyz_GradLatT + xyz_SLat*xyz_DzT)    
    xyz_FSig(:,:,:) = xyz_AI*( &
         &             xyz_SLon*xyz_GradLonT + xyz_SLat*xyz_GradLatT &
         &          + (xyz_SLon**2 + xyz_SLat**2)*xyz_DzT            &
         &         ) 
    !$omp end parallel workshare

    if(DFM08Flag) then
       call TaperingDFM08_IDIFF(xyz_C, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, xyz_Depth )
       xyz_FLon(:,:,:) = xyz_C*xyz_AI*xyz_GradLonT + (1d0 - xyz_C)*xyz_FLon
       xyz_FLat(:,:,:) = xyz_C*xyz_AI*xyz_GradLatT + (1d0 - xyz_C)*xyz_FLat
    end if

    xyz_FSig(:,:,0) = 0d0; xyz_FSig(:,:,kMax) = 0d0

  end subroutine calc_IsopycDiffFlux

  subroutine calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig,                  &  ! (out)
       & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth      &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: &
         & xyz_FLon, xyz_FLat, xyz_FSig
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth
    real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DzT, xyz_PsiLon, xyz_PsiLat
    integer :: k

    ! 実行文; Executable statement
    !

    xyz_DzT(:,:,:) = xyz_Dz_xyz(xyz_Tracer)

    !$omp parallel do
    do k=1, kMax-1
       xyz_PsiLon(:,:,k) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
       xyz_PsiLat(:,:,k) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
    end do

    xyz_PsiLat(:,:,0) = 0d0; xyz_PsiLon(:,:,0) = 0d0
    xyz_PsiLat(:,:,kMax) = 0d0; xyz_PsiLon(:,:,kMax) = 0d0
    
    if(DFM08Flag) then
       call TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Depth )
    end if
    
    !$omp parallel do    
    do k=0, kMax
       xyz_FLon(:,:,k) = - xyz_PsiLat(:,:,k)*xyz_DzT(:,:,k)
       xyz_FLat(:,:,k) =   xyz_PsiLon(:,:,k)*xyz_DzT(:,:,k)
       xyz_FSig(:,:,k) = ( &
            &   xyz_PsiLat(:,:,k)*xy_GradLon_w(wz_Tracer(:,k)) &
            & - xyz_PsiLon(:,:,k)*xy_GradLat_w(wz_Tracer(:,k)) &
            & )/RPlanet
    end do

  end subroutine calc_SkewFlux

  
  subroutine calc_IsoNeutralSlope( xyz_SLon, xyz_SLat,          &  ! (out)
       & xyz_DensPot, xyz_totDepth                              &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_SLon, xyz_SLat
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot, xyz_totDepth

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DzDensPot, xyz_Slimt
    real(DP) :: w_DensPot(lMax)
    integer :: i, j, k
    real(DP), parameter :: EPS = 1d-10

    ! 実行文; Executable statement
    !

    xyz_DzDensPot(:,:,:) = xyz_Dz_xyz(xyz_DensPot, isUsedDF=.true.)

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
       & xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, & ! (out)
       & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T

    ! 局所変数
    ! Local variables
    !
    integer :: k
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PsiLon, xyz_PsiLat, xyz_CosLat

    ! 実行文; Executable statement
    !

    xyz_CosLat(:,:,:) = cos(xyz_Lat)


    !$omp parallel do
    do k=1, kMax-1
       xyz_PsiLon(:,:,k) = - Kappa_GM*xyz_T(:,:,k)*xyz_SLat(:,:,k)
       xyz_PsiLat(:,:,k) =   Kappa_GM*xyz_T(:,:,k)*xyz_SLon(:,:,k)       
    end do

    xyz_PsiLon(:,:,0) = 0d0; xyz_PsiLat(:,:,0) = 0d0    
    xyz_PsiLon(:,:,kMax) = 0d0; xyz_PsiLat(:,:,kMax) = 0d0

    if(DFM08Flag) then
       call TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
            & DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb, xyz_Depth )
    end if
    
    xyz_EddInducedU(:,:,:) = - xyz_Dz_xyz(xyz_PsiLat)
    xyz_EddInducedV(:,:,:) =   xyz_Dz_xyz(xyz_PsiLon*xyz_CosLat)/xyz_CosLat 
    xyz_EddInducedW(:,:,:) = xyz_wz( &
         & wz_AlphaOptr_xyz(xyz_PsiLat*xyz_CosLat, -xyz_PsiLon*xyz_CosLat) &
         & )

  end subroutine calc_BolusVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !
    use SGSEddyMixingHelper_mod, only: &
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
    real(DP) :: DiffCoefRedi, DiffCoefGM
    character(TOKEN) :: InteriorTaperingName, PBLTaperingName 
    real(DP) :: SlopeMax, Sd
    logical :: isUsedDFM08
    logical :: DiagOutputFlag
    character(STRING) :: msg

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /SGSEddyMix_nml/ &
         & EddyMixTypeName, &
         & DiffCoefRedi, DiffCoefGM, &
         & InteriorTaperingName, PBLTaperingName, &
         & SlopeMax, Sd, &
         & isUsedDFM08, &
         & DiagOutputFlag

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    EddyMixTypeName = SGSEddyMixing_Redi_NAME
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
            & nml = SGSEddyMix_nml,          &  ! (out)
            & iostat = iostat_nml, iomsg=msg )  ! (out)
       close( unit_nml )
    end if

    !
    select case(EddyMixTypeName)
    case(SGSEddyMixing_Redi_NAME)
       SGSEddyMixType = SGSEddyMixing_Redi
    case(SGSEddyMixing_GM_NAME)
       SGSEddyMixType = SGSEddyMixing_GM90
    case default
       call MessageNotify('E', module_name, &
            & ' EddyMixTypeName ''%c'' is invalid.', c1=trim(EddyMixTypeName) )
    end select

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
    call MessageNotify( 'M', module_name, ' EddyMixTypeName   = %c', c1=EddyMixTypeName )
    call MessageNotify( 'M', module_name, ' DiffCoefRedi      = %f', d=(/ DiffCoefRedi  /) )
    call MessageNotify( 'M', module_name, ' DiffCoefGM        = %f', d=(/ DiffCoefGM  /) )
    call MessageNotify( 'M', module_name, ' InteriorTaperngName = %c', c1=InteriorTaperingName )
    call MessageNotify( 'M', module_name, ' PBLTaperngName      = %c', c1=PBLTaperingName )
    call MessageNotify( 'M', module_name, ' SlopeMax          = %f', d=(/ SlopeMax  /) )
    call MessageNotify( 'M', module_name, ' Sd                = %f', d=(/ Sd  /) )    
    call MessageNotify( 'M', module_name, ' isUsedDFM08       = %b', L=(/ isUsedDFM08 /) )    
    call MessageNotify( 'M', module_name, ' DiagOutputFlag    = %b', L=(/ DiagOutputFlag /) )
    
  end subroutine read_nmlData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine SGSEddyMixing_Output()

    use VariableSet_mod, only: &
         & z_PTempBasic, xyz_PTempEddN, xyz_SaltN, xy_SurfHeightN, &
         & xy_totDepthBasic

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
         & xyz_DensPot, xyz_PTemp, xyz_RefPress, &
         & xyz_SLon, xyz_SLat, xyz_totDepth, xyz_Depth, &
         & xyz_FLon, xyz_FLat, xyz_FSig, xyz_Tracer

    real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
         & xyz_BolusU, xyz_BolusV, xyz_BolusW, &
         & xyz_T, xyz_G, xyz_C

    real(DP), dimension(0:iMax-1,jMax)  :: xy_BLD

    integer :: k

    ! 実行文; Executable statement
    !

    if(.not. OutputFlag) return 

    do k=0, kMax
       xyz_PTemp(:,:,k) = z_PTempBasic(k) + xyz_PTempEddN(:,:,k)
       xyz_totDepth(:,:,k) = xy_totDepthBasic(:,:) + xy_SurfHeightN(:,:)
       xyz_Depth(:,:,k) = xyz_totDepth(:,:,k)*g_Sig(k)
    end do

    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot, &  !(out)
         & xyz_PTemp, xyz_SaltN, xyz_RefPress ) !(in)

    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_totDepth ) !(in)

    !
    !
    call prepare_SlopeTapering(xyz_T, xyz_SLon, xyz_SLat,  & ! (inout)
         & xyz_DensPot, xyz_Depth,                         & ! (in)
         & xy_BLD                                          & ! (out)
         & )
    
!    call check_StaticStability(xyz_T, xyz_DensPot)

    
    xyz_Tracer = xyz_PTemp
    call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,                &  !(out)      
         & xyz_Tracer, wz_xyz(xyz_Tracer), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth &  !(in)
         & )
!!$    call HistoryPut('Etc1', xyz_FLat, hst_SGSEddyMix)
!!$    call HistoryPut('Etc2', xyz_FSig, hst_SGSEddyMix)
!!$    call HistoryPut('Etc3', &
!!$         & +   xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*cos(xyz_Lat), xyz_FLat*cos(xyz_Lat)) ) & 
!!$         & +   xyz_Dz_xyz(xyz_FSig)  &
!!$         &    , hst_SGSEddyMix )
    
    call calc_SkewFlux( &
         & xyz_FLon, xyz_FLat, xyz_FSig, &
         & xyz_Tracer, wz_xyz(xyz_Tracer), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )
!!$    call HistoryPut('Etc4', xyz_FLat, hst_SGSEddyMix)
!!$    call HistoryPut('Etc5', xyz_FSig, hst_SGSEddyMix)
!!$    call HistoryPut('Etc6', &
!!$         & +   xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*cos(xyz_Lat), xyz_FLat*cos(xyz_Lat)) ) & 
!!$         & +   xyz_Dz_xyz(xyz_FSig)  &
!!$         &    , hst_SGSEddyMix )

!!$    call HistoryPut('DensPot', xyz_DensPot, hst_SGSEddyMix)
    call HistoryPut('SlopeLon', xyz_T*xyz_SLon, hst_SGSEddyMix)
    call HistoryPut('SlopeLat', xyz_T*xyz_SLat, hst_SGSEddyMix)
    call HistoryPut('diffCoefEff', xyz_T, hst_SGSEddyMix)
!    call HistoryPut('MixLyrDepth', xy_BLD, hst_SGSEddyMix)

    !
    if(SGSEddyMixType == SGSEddyMixing_GM90) then

       call calc_BolusVelocity(xyz_BolusU, xyz_BolusV, xyz_BolusW, & !(out)
            & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T )

       call HistoryPut('BolusU', xyz_BolusU, hst_SGSEddyMix)
       call HistoryPut('BolusV', xyz_BolusV, hst_SGSEddyMix)
       call HistoryPut('BolusW', xyz_BolusW, hst_SGSEddyMix)
    end if

  end subroutine SGSEddyMixing_Output

  subroutine SGSEddyMixing_PrepareOutput(OriginTime, EndTime, Intrv, FilePrefix)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: OriginTime, EndTime, Intrv
    character(*), intent(in) :: FilePrefix

    ! 局所変数
    ! Local variables
    !
    Character(TOKEN) :: lonName, latName, sigName, timeName


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

    call HistoryPut(lonName, xyz_Lon(:,1,0)*180d0/PI, hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'topology', 'circular', hst_SGSEddyMix)
    call HistoryAddAttr(lonName, 'modulo', 360.0, hst_SGSEddyMix)
    call HistoryPut(latName, xyz_Lat(0,:,0)*180d0/PI, hst_SGSEddyMix)
    call HistoryPut(sigName, g_Sig, hst_SGSEddyMix)

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
    
    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       call regist_XYZTVariable('BolusU', 'bolus velocity(longitude)', 'm/s')
       call regist_XYZTVariable('BolusV', 'bolus velocity(meridional)', 'm/s')
       call regist_XYZTVariable('BolusW', 'bolus velocity(vertical)', 'm/s')
    end if

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

  end subroutine SGSEddyMixing_PrepareOutput

end module SGSEddyMixing_mod

