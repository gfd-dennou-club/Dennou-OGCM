!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
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

  ! Dennou-OGCM


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
  public :: SGSEddyMixing_AddMixingTerm
  public :: calc_IsoNeutralSlope, calc_BolusVelocity
  public :: TaperingFunc_DM95, TaperingFunc_LDD95
  public :: prepare_DFM08Info
  
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
  real(DP) :: SlopeMaxVal
  real(DP) :: Sd

  type(gt_history), save :: hst_SGSEddyMix
  logical :: OutputFlag

  type :: DiabaticLyrInfo_DFM08
     real(DP), pointer :: xy_BLD(:,:)
     real(DP), pointer :: xy_TLT(:,:)
     real(DP), pointer :: xy_Lamb(:,:)
  end type DiabaticLyrInfo_DFM08

  type(DiabaticLyrInfo_DFM08), save :: DFM08Info
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
    SlopeMaxVal = 4d-3
    Sd = 1d-3
    OutputFlag = .false.

    !

    if(present(SGSEddyMixParamType)) SGSEddyMixType = SGSEddyMixParamType
    if(present(KappaRedi)) Kappa_Redi = KappaRedi
    if(present(KappaGM)) Kappa_GM = KappaGM
    if(present(isVarsOutput)) OutputFlag = isVarsOutput

    if(present(configNmlFileName)) call read_nmlData(configNmlFileName)

    DFM08Flag = .false.
    allocate( &
         & DFM08Info%xy_BLD(0:iMax-1,jMax), DFM08Info%xy_TLT(0:iMax-1,jMax), &
         & DFM08Info%xy_Lamb(0:iMax-1,jMax) )
    
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
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)

    if(DFM08Flag) then
       deallocate(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
       nullify(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
    end if
    
  end subroutine SGSEddyMixing_Final

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixing_AddMixingTerm( wz_PTempRHS, wz_SaltRHS, &  ! (inout)
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
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DensPot, xyz_RefPress, xyz_CosLat, xyz_totDepth

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, &  !< The longitude componet of th slope of isopycnal surface
         & xyz_SLat     !< The latitude componet of th slope of isopycnal surface

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTemp, xyz_Salt, xyz_KappRho, xyz_Depth

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_T
    
    integer :: k

    ! 実行文; Executable statement
    !

    xyz_CosLat(:,:,:) = spread(cos(xyz_Lat(:,:,0)),3,kMax+1)

    !$omp parallel do
    do k=0, kMax
       xyz_totDepth(:,:,k) = xy_totDepth
       xyz_Depth(:,:,k) = xy_totDepth*g_Sig(k)
       xyz_PTemp(:,:,k) = xy_w(wz_PTemp(:,k))
       xyz_Salt(:,:,k) = xy_w(wz_Salt(:,k))
    end do
    
    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot, &         !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress ) !(in)

    
    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_totDepth )               !(in)

    !
    
    !$omp parallel do
    do k=0, kMax
       xyz_T(:,:,k) = TaperingFunc_DM95( &
            & xyz_SLon(:,:,k), xyz_SLat(:,:,k), xyz_Depth(:,:,k), xyz_Depth(:,:,kMax), &
            & 2d0/(2d0*Omega*sqrt(1d0 - xyz_CosLat(:,:,k)**2)) &
            & )
    end do
    xyz_T(:,:,:) = xyz_T(:,:,:)*TaperingFunc_LDD95(xyz_SLon, xyz_SLat, xyz_Depth, &
         & 2d0/(2d0*Omega*sqrt(1d0 - xyz_CosLat(:,:,0)**2)) )

    if(DFM08Flag) then
       call prepare_DFM08Info( &
            & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth )
    end if
    
    call check_StaticStability(xyz_T, xyz_DensPot)

    !
    wz_PTempRHS(:,:) = wz_PTempRHS + calc_IsopycDiffTerm(wz_PTemp, xyz_PTemp)
    wz_SaltRHS(:,:)  = wz_SaltRHS + calc_IsopycDiffTerm(wz_Salt, xyz_Salt)

    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       wz_PTempRHS(:,:) = wz_PTempRHS + calc_EddyInducedVelAdvTerm(wz_PTemp, xyz_PTemp)
       wz_SaltRHS(:,:)  = wz_SaltRHS + calc_EddyInducedVelAdvTerm(wz_Salt, xyz_Salt)
    end if

  contains
    function calc_IsopycDiffTerm(wz_Tracer, xyz_Tracer) result(wz_DiffTerm)

      ! 宣言文; Declaration statement
      !
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
      real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
      real(DP) :: wz_DiffTerm(lMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_FLon, xyz_FLat, xyz_FSig

      ! 実行文; Executable statement
      !

      call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,                &  !(out)      
           & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth &  !(in)
           & )

      wz_DiffTerm(:,:) = &
           &    wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat) & 
           & +  wz_xyz(xyz_Dz_xyz(xyz_FSig))

    end function calc_IsopycDiffTerm

    function calc_EddyInducedVelAdvTerm(wz_Tracer, xyz_Tracer) result(wz_AdvTerm)

      ! 宣言文; Declaration statement
      !
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
      real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
      real(DP) :: wz_AdvTerm(lMax,0:kMax)

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
      wz_AdvTerm(:,:) = &
           &   wz_AlphaOptr_xyz( xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat ) &
           & + wz_xyz( xyz_Dz_xyz(xyz_FSig) )

    end function calc_EddyInducedVelAdvTerm

  end subroutine SGSEddyMixing_AddMixingTerm

  subroutine calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,          &  ! (out)
       & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth &  ! (in)
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
    xyz_FLon(:,:,:) = xyz_AI*(xyz_GradLonT + xyz_SLon*xyz_DzT)
    xyz_FLat(:,:,:) = xyz_AI*(xyz_GradLatT + xyz_SLat*xyz_DzT)
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
    
    xyz_FSig(:,:,0) = 0d0
    xyz_FSig(:,:,kMax) = 0d0

  end subroutine calc_IsopycDiffFlux

  subroutine calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig,                  &  ! (out)
       & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth   &  ! (in)
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

  subroutine check_StaticStability(xyz_KappRho, xyz_DensPot)

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: xyz_KappRho
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot

    ! 局所変数
    ! Local variables
    !
    integer :: k

    ! 実行文; Executable statement
    !

    return
    
!    !$omp parallel do
    do k=0, kMax-1
       where(xyz_DensPot(:,:,k) > xyz_DensPot(:,:,k+1))
          xyz_KappRho(:,:,k) = 0d0
          xyz_KappRho(:,:,k+1) = 0d0
       end where
    end do

  end subroutine check_StaticStability

  elemental function TaperingFunc_DM95(SLon, SLat, Depth, totDepth, BarocEddDispH) result(f)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: SLon, SLat, Depth, totDepth, BarocEddDispH
    real(DP) :: f

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: r, BarocEddyDispZ, SlopeABS

    ! 実行文; Executable statement
    !

    SlopeABS = sqrt(SLon**2 + SLat**2)
    f = 0.5d0*(1d0 + tanh((SlopeMaxVal - SlopeABS)/Sd))
    
    !
!!$    BarocEddyDispZ = min(max(15d3, BarocEddDispH), 100d3)*SlopeABS
!!$    r =  min(0d0 - Depth, Depth - totDepth)/BarocEddyDispZ
!!$    if(r < 1d0) then
!!$       f = f* 0.5d0*(1d0 + sin(PI*(r - 0.5d0)))
!!$    end if

  end function TaperingFunc_DM95

  function TaperingFunc_LDD95(xyz_SLon, xyz_SLat, xyz_Depth, xy_BarocEddDispH) result(xyz_f)

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_SLon, xyz_SLat, xyz_Depth
    real(DP), intent(in) :: xy_BarocEddDispH(0:iMax-1,jMax)
    real(DP) :: xyz_f(0:iMax-1, jMax, 0:kMax)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: r, xy_BarocEddyDispZ(0:iMax-1,jMax), tmpBarocEddyDispZ, SlopeABS
    real(DP) :: xyz_r(0:iMax-1,jMax,0:kMax)
    integer :: i, j, k, kStart(1)
    
    ! 実行文; Executable statement
    !

    
    do j=1,jMax
       do i=0, iMax-1
          xy_BarocEddyDispZ(i,j) = SlopeMaxVal*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
       end do
    end do

    xyz_r(:,:,:) = 1d0
    do j=1,jMax
       do i=0, iMax-1

          !
          kStart = minloc(abs(xy_BarocEddyDispZ(i,j)-(0d0-xyz_Depth(i,j,:))))
          do k=kStart(1), 0, -1
             SlopeABS = sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
             tmpBarocEddyDispZ = SlopeABS*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
             if((0d0 - xyz_Depth(i,j,k))/tmpBarocEddyDispZ < 1d0) then
                xyz_r(i,j,0:k) = (0d0 - xyz_Depth(i,j,0:k))/tmpBarocEddyDispZ; exit
             end if
          end do

          kStart = minloc(abs(xy_BarocEddyDispZ(i,j)-(xyz_Depth(i,j,:)-xyz_Depth(i,j,kMax))))
          do k=kStart(1), kMax
             SlopeABS = sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
             tmpBarocEddyDispZ = SlopeABS*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
             if((xyz_Depth(i,j,k)-xyz_Depth(i,j,kMax))/tmpBarocEddyDispZ < 1d0) then
                xyz_r(i,j,k:kMax) = (xyz_Depth(i,j,k:kMax) - xyz_Depth(i,j,kMax))/tmpBarocEddyDispZ; exit
             end if
          end do


       end do
    end do

    !
    xyz_f(:,:,:) = 1d0
    where (xyz_r < 1d0)
       xyz_f = 0.5d0*(1d0 + sin(PI*(xyz_r - 0.5d0)))
    end where

  end function TaperingFunc_LDD95

  
  subroutine prepare_DFM08Info( &
       & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth, &
       & xyz_G )

    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth
    real(DP), intent(out), optional :: xyz_G(0:iMax-1,jMax,0:kMax)
    
    integer :: i, j
    integer :: MLD_k, DLD_k
    real(DP) :: DLD, BLD, TLT

    real(DP), parameter :: c = 2d0
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_DzDensPot, xyz_DzzDensPot
    real(DP), dimension(0:kMax)  :: z_DzDensPot, z_DzzDensPot
    real(DP) :: wt(2)
    real(DP) :: lamb, z, SlopeABS, ddzDensPot, ddzDensPotMax
    integer :: k, DzMax_k, minLocId(1)
    
    xyz_DzDensPot = xyz_Dz_xyz(xyz_DensPot)
    xyz_DzzDensPot = xyz_Dz_xyz(xyz_DzDensPot)

    do j=1, jMax
       do i=0, iMax-1

          z_DzDensPot(:) = xyz_DzDensPot(i,j,:)
          z_DzzDensPot(:) = xyz_DzzDensPot(i,j,:)

          ddzDensPotMax = 0d0
          do k=1, kMax
             ddzDensPot = -(xyz_DensPot(i,j,0) - xyz_DensPot(i,j,k))/(xyz_Depth(i,j,0) - xyz_Depth(i,j,k))
             if(ddzDensPotMax < ddzDensPot) then
                ddzDensPotMax = ddzDensPot
             else if(k > 1) then
                DzMax_k = k
                exit
             end if
          end do
          
          !
          do k=1, kMax
             if((z_DzDensPot(k) + ddzDensPotMax)*(z_DzDensPot(k-1) + ddzDensPotMax) <= 0d0) then
                MLD_k = k; exit
             end if
             if(k==kMax) then
                minLocId(:) = minloc(abs(z_DzzDensPot(:) + ddzDensPotMax))
                MLD_k = minLocId(1)
!!$                write(*,*) "j=", j, ": ddzDensPotMax=", ddzDensPotMax
!!$                write(*,*) z_DzDensPot
             end if
          end do
          wt(:) = abs(z_DzDensPot(MLD_k:MLD_k-1:-1) + ddzDensPotMax)/sum(abs(z_DzDensPot(MLD_k-1:MLD_k) + ddzDensPotMax)) 

          !
          BLD = abs(sum(wt(:)*xyz_Depth(i,j,MLD_k-1:MLD_k)))

          SlopeABS = sqrt(xyz_SLon(i,j,MLD_k)**2 + xyz_SLat(i,j,MLD_k)**2)
          TLT =  min(SlopeABS,SlopeMaxVal)*min(max(15d3,c/(2d0*Omega*abs(sin(xyz_Lat(i,j,0))))), 100d3) + 1d0

          DLD = min(BLD + TLT, abs(xyz_Depth(i,j,kMax)))

          DLD_k = MLD_k
          do k=MLD_k, kMax
             if( xyz_Depth(i,j,k) <= -DLD ) then
                DLD_k = k; exit
             end if
             if(k==kMax) DLD_k = k
          end do
          
          wt(:) = abs(xyz_Depth(i,j,DLD_k:DLD_k-1:-1) + DLD)/sum(abs(xyz_Depth(i,j,DLD_k-1:DLD_k) + DLD))

          !
          lamb = abs( - sum(wt(:)*z_DzDensPot(DLD_k-1:DLD_k))/sum(wt(:)*z_DzzDensPot(DLD_k-1:DLD_k)) )

          !
          DFM08Info%xy_BLD(i,j) = BLD; DFM08Info%xy_TLT(i,j) = TLT; 
          DFM08Info%xy_Lamb(i,j) = lamb

!!$          if(DLD_k == kMax) then
!!$             write(*,*) "DensPot", xyz_DensPot(i,j,:)
!!$             write(*,*) "Dz", z_DzDensPot
!!$             write(*,*) "Max", DzMax_k, ddzDensPotMax
!!$             stop
!!$          end if
       end do
    end do


    if(present(xyz_G)) then
       xyz_G = 1d0
       do k=0, kMax
          where (-DFM08Info%xy_BLD(:,:) < xyz_Depth(:,:,k))
             xyz_G(:,:,k) = -xyz_Depth(:,:,k)/(2d0*DFM08Info%xy_BLD + DFM08Info%xy_TLT) &
                  &        *(2d0 + DFM08Info%xy_TLT/DFM08Info%xy_Lamb)
          end where
       end do
    end if
!!$
!!$    write(*,*) "-- DFM08INFO ---------"
!!$    do j=1, jMax/2+1
!!$    write(*,*) "j=", j, "BLD=", DFM08Info%xy_BLD(0,j), "TLT=", DFM08Info%xy_TLT(0,j), "lamb=", DFM08Info%xy_Lamb(0,j)
!!$    end do

  end subroutine prepare_DFM08Info

  subroutine TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
       & xy_BLD, xy_TLT, xy_Lamb, xyz_Depth )

    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PsiLon, xyz_PsiLat
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_BLD, xy_TLT, xy_Lamb
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Depth

    integer :: k, i, j
    real(DP), dimension(0:iMax-1,jMax) :: xy_DLD, xy_G, xy_PsiILon, xy_PsiILat
    real(DP) :: G
    logical :: terminate
    
    xy_DLD = xy_BLD + xy_TLT

    xy_PsiILon = 0d0; xy_PsiILat = 0d0
    do k=1, kMax
       where(xyz_Depth(:,:,k) < -xy_DLD .and. xyz_Depth(:,:,k-1) >= -xy_DLD)
          xy_PsiILon(:,:) = xyz_PsiLon(:,:,k);  xy_PsiILat(:,:) = xyz_PsiLat(:,:,k)
       end where
    end do

    terminate = .false.
    do k=0, kMax
       do j=1, jMax
          do i=0, iMax-1
             if(-xy_BLD(i,j) < xyz_Depth(i,j,k)) then
                G = - xyz_Depth(i,j,k)/(2d0*xy_BLD(i,j) + xy_TLT(i,j))*(2d0 + xy_TLT(i,j)/xy_Lamb(i,j))
                xyz_PsiLon(i,j,k) = G*xy_PsiILon(i,j); xyz_PsiLat(i,j,k) = G*xy_PsiILat(i,j);
!!$                if(j== 5) then
!!$                   write(*,*) "*", j, k, ":", xy_BLD(i,j), xy_DLD(i,j), G
!!$                end if

             else if(-xy_DLD(i,j) < xyz_Depth(i,j,k)) then
                G =   - (xyz_Depth(i,j,k) + xy_BLD(i,j))**2/(xy_DLD(i,j)**2 - xy_BLD(i,j)**2)*(1d0 + xy_DLD(i,j)/xy_Lamb(i,j)) &
                     & - xyz_Depth(i,j,k)/(2d0*xy_BLD(i,j) + xy_TLT(i,j))*(2d0 + xy_TLT(i,j)/xy_Lamb(i,j))

!                G = 1d0
                xyz_PsiLon(i,j,k) = G*xy_PsiILon(i,j); xyz_PsiLat(i,j,k) = G*xy_PsiILat(i,j)
!!$                if(j<= 32) then
!!$                   write(*,*) j, k, ":", xy_BLD(i,j), xy_DLD(i,j), G
!!$                end if
             else if(-xy_DLD(i,j) <  xyz_Depth(i,j,kMax))then
!!$                write(*,*) "Error", j, xy_DLD(i,j)
                terminate = .true.
             end if
          end do
       end do
       if(terminate) stop
!!$       where(-xy_BLD < xyz_Depth(:,:,k))
!!$          xy_G = -xyz_Depth(:,:,k)/(2d0*xy_BLD + xy_TLT)*(2d0 + xy_TLT/xy_Lamb)
!!$          xyz_PsiLon(:,:,k) = xy_G*xy_PsiILon; xyz_PsiLat(:,:,k) = xy_G*xy_PsiILat;
!!$       end where
!!$       where(-xy_DLD < xyz_Depth(:,:,k) .and. -xy_BLD > xyz_Depth(:,:,k) )
!!$          xy_G = 1d0!(xyz_Depth(:,:,k) + xy_BLD)**2/(xy_DLD**2 - xy_BLD**2)*(1d0 + xy_DLD/xy_Lamb) &
!!$!               & - xyz_Depth(:,:,k)/(2d0*xy_BLD + xy_TLT)*(2d0 + xy_TLT/xy_Lamb)
!!$          xyz_PsiLon(:,:,k) = xy_G*xy_PsiILon; xyz_PsiLat(:,:,k) = xy_G*xy_PsiILat;          
!!$       end where
    end do

  end subroutine TaperingDFM08_GM

  subroutine TaperingDFM08_IDIFF(xyz_C, &
       & xy_BLD, xy_TLT, xyz_Depth )

    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: xyz_C
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_BLD, xy_TLT
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Depth

    integer :: k
    real(DP), dimension(0:iMax-1,jMax) :: xy_DLD, xy_C

    
    xy_DLD(:,:) = xy_BLD + xy_TLT

    xyz_C(:,:,:) = 0d0
    
    do k=0, kMax
       where(-xy_BLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = 1d0
       elsewhere(-xy_DLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = (xyz_Depth(:,:,k) + xy_DLD)/xy_TLT
       end where
    end do
    
  end subroutine TaperingDFM08_IDIFF
  
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

  function xyz_Dz_xyz(xyz, isUSedDF) 

    use VariableSet_mod

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
    logical, intent(in), optional :: isUSedDF
    real(DP) :: xyz_Dz_xyz(0:iMax-1,jMax,0:kMax)
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: s(0:kMax), t(0:kMax), xyt(0:iMax-1,jMax,0:tMax)
    integer :: k

    ! 実行文; Executable statement
    !

!!$    if(present(isUsedDF) .and. isUsedDF) then
    
    
    t(1:kMax-1) = g_Sig(0:kMax-2) - g_Sig(1:kMax-1)
    s(1:kMax-1) = g_Sig(1:kMax-1) - g_Sig(2:kMax)

    !$omp parallel do
    do k=1,kMax-1
       xyz_Dz_xyz(:,:,k) = &
            & (s(k)**2*xyz(:,:,k-1) - (s(k)**2-t(k)**2)*xyz(:,:,k) - t(k)**2*xyz(:,:,k+1)) &
            & /(s(k)*t(k)*(s(k) + t(k)))/xy_totDepthBasic
    end do
    xyz_Dz_xyz(:,:,0) = &
         & (xyz(:,:,0) - xyz(:,:,1))/(g_Sig(0) - g_Sig(1))/xy_totDepthBasic
    xyz_Dz_xyz(:,:,kMax) = &
         & (xyz(:,:,kMax-1) - xyz(:,:,kMax))/(g_Sig(kMax-1) - g_Sig(kMax))/xy_totDepthBasic

!!$ else
!!$        xyz_Dz_xyz = xyz_DSig_xyz(xyz)/spread(xy_totDepthBasic,3,kMax+1)
!!$     end if
  end function xyz_Dz_xyz

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
    logical :: DiagOutputFlag
    character(STRING) :: msg

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /SGSEddyMix_nml/ &
         & EddyMixTypeName, &
         & DiffCoefRedi, DiffCoefGM, DiagOutputFlag

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    EddyMixTypeName = SGSEddyMixing_Redi_NAME
    DiffCoefRedi = Kappa_Redi
    DiffCoefGM = Kappa_GM
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

    select case(EddyMixTypeName)
    case(SGSEddyMixing_Redi_NAME)
       SGSEddyMixType = SGSEddyMixing_Redi
    case(SGSEddyMixing_GM_NAME)
       SGSEddyMixType = SGSEddyMixing_GM90
    case default
       call MessageNotify('E', module_name, &
            & ' EddyMixTypeName ''%c'' is invalid.', c1=trim(EddyMixTypeName) )
    end select

    Kappa_Redi = DiffCoefRedi
    Kappa_GM = DiffCoefGM
    OutputFlag = DiagOutputFlag
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, ' EddyMixTypeName   = %c', c1=EddyMixTypeName )
    call MessageNotify( 'M', module_name, ' DiffCoefRedi      = %f', d=(/ DiffCoefRedi  /) )
    call MessageNotify( 'M', module_name, ' DiffCoefGM        = %f', d=(/ DiffCoefGM  /) )
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
    do k=0, kMax
       xyz_T(:,:,k) = TaperingFunc_DM95( &
            & xyz_SLon(:,:,k), xyz_SLat(:,:,k), xyz_Depth(:,:,k), xyz_Depth(:,:,kMax), &
            & 2d0/(2d0*Omega*abs(sin(xyz_Lat(:,:,k)))) &
            & )
    end do
    xyz_T(:,:,:) = xyz_T(:,:,:)*TaperingFunc_LDD95( xyz_SLon, xyz_SLat, xyz_Depth, &
         & 2d0/(2d0*Omega*abs(sin(xyz_Lat(:,:,0)))) )
    xyz_G = 1d0
    if(DFM08Flag) then
       call prepare_DFM08Info( &
            & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth, &
            & xyz_G )
    end if
    
    call check_StaticStability(xyz_T, xyz_DensPot)

    
    xyz_Tracer = xyz_PTemp
    call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,                &  !(out)      
         & xyz_Tracer, wz_xyz(xyz_Tracer), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth &  !(in)
         & )
    
    call HistoryPut('Etc1', xyz_FLat, hst_SGSEddyMix)
    call HistoryPut('Etc2', xyz_FSig, hst_SGSEddyMix)
    call HistoryPut('Etc3', &
         & +   xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*cos(xyz_Lat), xyz_FLat*cos(xyz_Lat)) ) & 
         & +   xyz_Dz_xyz(xyz_FSig)  &
         &    , hst_SGSEddyMix )

    call calc_SkewFlux( &
         & xyz_FLon, xyz_FLat, xyz_FSig, &
         & xyz_Tracer, wz_xyz(xyz_Tracer), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )
    call HistoryPut('Etc4', xyz_FLat, hst_SGSEddyMix)
    call HistoryPut('Etc5', xyz_FSig, hst_SGSEddyMix)
    call HistoryPut('Etc6', &
         & +   xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*cos(xyz_Lat), xyz_FLat*cos(xyz_Lat)) ) & 
         & +   xyz_Dz_xyz(xyz_FSig)  &
         &    , hst_SGSEddyMix )
      
    call HistoryPut('DensPot', xyz_DensPot, hst_SGSEddyMix)
    call HistoryPut('SlopeLon', xyz_T*xyz_SLon, hst_SGSEddyMix)
    call HistoryPut('SlopeLat', xyz_T*xyz_SLat, hst_SGSEddyMix)
    call HistoryPut('diffCoefEff', xyz_T*xyz_G, hst_SGSEddyMix)

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

    !
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: lonName, latName, sigName, timeName
    character(TOKEN) :: dims_XYZT(4)

    ! 実行文; Executable statement
    !

    if(.not. OutputFlag) return 

    
    !
    lonName = 'lon'; latName='lat'; sigName='sig'; timeName='t'
    dims_XYZT = (/ lonName, latName, sigName, timeName /)

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

    call regist_Variable('Etc1', 'tendency term1', 's-1')
    call regist_Variable('Etc2', 'tendency term2', 's-1')
    call regist_Variable('Etc3', 'tendency term3', 's-1')
    call regist_Variable('Etc4', 'tendency term1', 's-1')
    call regist_Variable('Etc5', 'tendency term2', 's-1')
    call regist_Variable('Etc6', 'tendency term3', 's-1')

    call regist_Variable('DensPot', 'potential density', 'kg/m3')
    call regist_Variable('SlopeLon', 'the longitude component of slope of isoneutral', '1')
    call regist_Variable('SlopeLat', 'the latitude component of slope of isoneutral', '1')
    call regist_Variable('diffCoefEff', 'rescale factor of diffusivity', '1')

    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       call regist_Variable('BolusU', 'bolus velocity(longitude)', 'm/s')
       call regist_Variable('BolusV', 'bolus velocity(meridional)', 'm/s')
       call regist_Variable('BolusW', 'bolus velocity(vertical)', 'm/s')
    end if

  contains
    subroutine regist_Variable(varName, long_name, units)
      character(*), intent(in) :: varName, long_name, units

      call HistoryAddVariable(varName, dims_XYZT, &
           & long_name, units, xtype='float', history=hst_SGSEddyMix)

    end subroutine regist_Variable
  end subroutine SGSEddyMixing_PrepareOutput

end module SGSEddyMixing_mod

