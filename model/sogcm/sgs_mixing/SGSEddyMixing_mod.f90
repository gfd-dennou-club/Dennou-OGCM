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
  use dc_types, only: &
       & DP, TOKEN, STRING

  use gtool_history

  !
  use Constants_mod, only: &
       & RPlanet

  use SpmlUtil_mod

  use DiagnoseUtil_mod

  use EOSDriver_mod

  use GridSet_mod

  

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SGSEddyMixing_Init, SGSEddyMixing_Final
  public :: SGSEddyMixing_PrepareOutput, SGSEddyMixing_Output
  public :: SGSEddyMixing_AddMixingTerm

  public :: calc_IsoNeutralSlope, calc_BolusVelocity

  ! 公開変数
  ! Public variable
  !
  integer, parameter, public :: SGSEddyMixing_Redi = 1
  integer, parameter, public :: SGSEddyMixing_GM90 = 2

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSEddyMixing_mod' !< Module Name

  real(DP) :: SGSEddyMixType
  real(DP) :: Kapper_rho  !< Isopycnal diffusivity
  real(DP) :: SlopeMaxVal, Sd

  type(gt_history), save :: hst_SGSEddyMix
  logical :: OutputFlag
  
contains

  !>
  !!
  !!
  subroutine SGSEddyMixing_Init(SGSEddyMixParamType, KapperRho)

    ! 実行文; Executable statements
    !
    integer, intent(in) :: SGSEddyMixParamType
    real(DP), intent(in) :: KapperRho

    SGSEddyMixType = SGSEddyMixParamType
    Kapper_rho = 1000d0!KapperRho
    SlopeMaxVal = 4d-3
    Sd = 1d-3!1d-3 !SlopeMaxVal/4d0
    OutputFlag = .false.

  end subroutine SGSEddyMixing_Init

  !>
  !!
  !!
  subroutine SGSEddyMixing_Final()

    ! 実行文; Executable statements
    !
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)

  end subroutine SGSEddyMixing_Final

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixing_AddMixingTerm(wz_PTempRHS, wz_SaltRHS, &
       & wz_PTemp, wz_Salt, xy_totDepth )
    
    use VariableSet_mod

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_PTempRHS, wz_SaltRHS
    real(DP), dimension(lMax,0:kMax), intent(in)  :: wz_PTemp, wz_Salt
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_totDepth

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DensPot, xyz_RefPress, xyz_CosLat, xyz_totDepth, xyz_T

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, &  !< The longitude componet of th slope of isopycnal surface
         & xyz_SLat     !< The latitude componet of th slope of isopycnal surface

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTemp, xyz_Salt, xyz_KappRho, xyz_Depth
    
    integer :: k

    ! 実行文; Executable statement
    !

    xyz_CosLat = cos(xyz_Lat)
    do k=0, kMax
       xyz_totDepth(:,:,k) = xy_totDepth
       xyz_Depth(:,:,k) = xy_totDepth*g_Sig(k)
    end do
    xyz_PTemp = xyz_wz(wz_PTemp)
    xyz_Salt = xyz_wz(wz_Salt)

    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot, &  !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress ) !(in)

    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_totDepth ) !(in)

    !
    xyz_T = eval_TaperingFunc(xyz_SLon, xyz_SLat, xyz_Depth)

!!$    xyz_KappRho = Kapper_rho*min(1d0, SlopeMaxVal/sqrt(xyz_SLon**2 + xyz_SLat**2))
    call check_StaticStability(xyz_T, xyz_DensPot)

    ! 
    wz_PTempRHS = wz_PTempRHS + calc_IsoNeutralDiffTerm(wz_PTemp, xyz_PTemp)
    wz_SaltRHS = wz_SaltRHS + calc_IsoNeutralDiffTerm(wz_Salt, xyz_Salt)

    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       wz_PTempRHS = wz_PTempRHS + calc_EddyInducedVelAdvTerm(wz_PTemp, xyz_PTemp)
       wz_SaltRHS = wz_SaltRHS + calc_EddyInducedVelAdvTerm(wz_Salt, xyz_Salt)
    end if

!!$    write(*,*) "MeanSlope", AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(sqrt(xyz_SLon**2+xyz_Slat**2)))
!!$    write(*,*) "Max Slope", maxval(sqrt(xyz_SLon**2+xyz_Slat**2)), maxloc(sqrt(xyz_SLon**2+xyz_Slat**2))

    contains


      function calc_IsoNeutralDiffTerm(wz_Tracer, xyz_Tracer) result(wz_DiffTerm)

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

        call calc_NeutralDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig,     &  !(out)      
             & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth &  !(in)
             & )

!!$        wz_DiffTerm = &
!!$             &    wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat) & 
!!$             & +  wz_xyz( &
!!$             &        xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_FSig)))/xyz_totDepth &
!!$             &    )
        wz_DiffTerm = &
             &    wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat) & 
             & +  wz_xyz(xyz_Dz_xyz(xyz_FSig))
        
      end function calc_IsoNeutralDiffTerm
    
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

        xyz_DzT = xyz_Dz_xyz(xyz_Tracer) !xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))/xyz_totDepth

        xyz_FLon = (Kapper_rho*xyz_T)*xyz_SLon*xyz_DzT
        xyz_FLat = (Kapper_rho*xyz_T)*xyz_SLat*xyz_DzT
        xyz_FSig = (Kapper_rho*xyz_T)*( &
             &   xyz_SLon*xyz_GradLon_wz(wz_Tracer) + xyz_SLat*xyz_GradLat_wz(wz_Tracer)  &
             &   )

        xyz_FSig(:,:,0) = 0d0
        xyz_FSig(:,:,kMax) = 0d0


        wz_AdvTerm = &
             !             & +  wz_AlphaOptr_xyz(xyz_EddInducedU*xyz_Tracer*xyz_CosLat, xyz_EddInducedV*xyz_Tracer*xyz_CosLat) &
             !             & +  wz_wt(wt_DSig_wt(wt_xyz(xyz_EddInducedW*xyz_Tracer/xyz_totDepth)))
             & - wz_AlphaOptr_xyz( xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat ) &
             & + wz_xyz( xyz_Dz_xyz(xyz_FSig) ) !xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_FSig)))/xyz_totDepth )

      end function calc_EddyInducedVelAdvTerm

    end subroutine SGSEddyMixing_AddMixingTerm

    subroutine calc_NeutralDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig, &
         & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth )

      ! 宣言文; Declaration statement
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: &
           & xyz_FLon, xyz_FLat, xyz_FSig
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
           & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
           & xyz_GradLonT, xyz_GradLatT, xyz_DzT

        
      ! 実行文; Executable statement
      !

      xyz_GradLonT = xyz_GradLon_wz(wz_Tracer)
      xyz_GradLatT = xyz_GradLat_wz(wz_Tracer)
!!$      xyz_DzT = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))/xyz_totDepth
      xyz_DzT = xyz_Dz_xyz(xyz_Tracer)

      xyz_FLon = Kapper_rho*(xyz_GradLonT + xyz_T*xyz_SLon*xyz_DzT)
      xyz_FLat = Kapper_rho*(xyz_GradLatT + xyz_T*xyz_SLat*xyz_DzT)
      xyz_FSig = Kapper_rho*xyz_T*( &
           &             xyz_SLon*xyz_GradLonT + xyz_SLat*xyz_GradLatT &
           &          + (xyz_SLon**2 + xyz_SLat**2)*xyz_DzT            &
           &         ) 

      xyz_FSig(:,:,0) = 0d0
      xyz_FSig(:,:,kMax) = 0d0

    end subroutine calc_NeutralDiffFlux

    subroutine calc_SkewFlux(xyz_FLon, xyz_FLat, xyz_FSig, &
         & xyz_Tracer, wz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth )

      ! 宣言文; Declaration statement
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: &
           & xyz_FLon, xyz_FLat, xyz_FSig
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
           & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth
      real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
           & xyz_GradLonT, xyz_GradLatT, xyz_DzT

        
      ! 実行文; Executable statement
      !
      xyz_GradLonT = xyz_GradLon_wz(wz_Tracer)
      xyz_GradLatT = xyz_GradLat_wz(wz_Tracer)
!!$      xyz_DzT = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))/xyz_totDepth
      xyz_DzT = xyz_Dz_xyz(xyz_Tracer)

      xyz_FLon = Kapper_rho*xyz_T*xyz_SLon*xyz_DzT
      xyz_FLat = Kapper_rho*xyz_T*xyz_SLat*xyz_DzT
      xyz_FSig = Kapper_rho*xyz_T*( &
           &             xyz_SLon*xyz_GradLonT + xyz_SLat*xyz_GradLatT &
           &         ) 

      xyz_FSig(:,:,0) = 0d0
      xyz_FSig(:,:,kMax) = 0d0

    end subroutine calc_SkewFlux

    subroutine check_StaticStability(xyz_KappRho, xyz_DensPot)
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: xyz_KappRho
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot

      integer :: k

      do k=0, kMax-1
         where(xyz_DensPot(:,:,k) >= xyz_DensPot(:,:,k+1))
            xyz_KappRho(:,:,k+1) = 0d0
         end where
      end do

    end subroutine check_StaticStability

    function eval_TaperingFunc(xyz_SLon, xyz_SLat, xyz_Depth) result(xyz_f)

      use Constants_mod, only: Omega, PI

      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_SLon, xyz_SLat, xyz_Depth
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_f

      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_r, xyz_BarocEddyDisplZ, xyz_SlopeABS
      real(DP), parameter :: c = 2d0
      integer :: k

      xyz_SlopeABS = sqrt(xyz_SLon**2 + xyz_SLat**2)
      xyz_f = 0.5d0*(1d0 + tanh((SlopeMaxVal - xyz_SlopeABS)/Sd ))

      xyz_BarocEddyDisplZ = min(max(15d3, c/abs(2d0*Omega*sin(xyz_Lat))), 100d3)

      xyz_r =  min(-xyz_Depth, xyz_Depth-spread(xyz_Depth(:,:,kMax),3,kMax))&
           &  / (xyz_BarocEddyDisplZ*xyz_SlopeABS)
      where(xyz_r < 1d0)
         xyz_f = xyz_f* 0.5d0*(1d0 + sin(PI*(xyz_r - 0.5d0)))
      end where

    end function eval_TaperingFunc


    subroutine calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &
         & xyz_DensPot, xyz_totDepth )

      ! 宣言文; Declaration statement
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_SLon, xyz_SLat
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot, xyz_totDepth

      ! 局所変数
      ! Local variables
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
           & xyz_GradDensPotLon, xyz_GradDensPotLat, &
           & xyz_DzDensPot, xyz_Slimt

      integer :: i, j, k
      real(DP), dimension(0:iMax-1,jMax,0:tMax) :: xyt_DensPot
      real(DP), parameter :: EPS = 1d-10


      ! 実行文; Executable statement
      !

      xyz_GradDensPotLon = xyz_GradLon_wz(wz_xyz(xyz_DensPot))
      xyz_GradDensPotLat = xyz_GradLat_wz(wz_xyz(xyz_DensPot))


!!$      xyt_DensPot = xyt_xyz(xyz_DensPot)
!!$      xyz_DzDensPot = xyz_xyt(xyt_DSig_xyt(xyt_DensPot))/xyz_totDepth
      xyz_DzDensPot = xyz_Dz_xyz(xyz_DensPot)

!!$      xyz_Slimt = -sqrt(xyz_GradDensPotLon**2 + xyz_GradDensPotLat**2)/SlopeMaxVal
!!$      xyz_DzDensPot = min(xyz_Slimt, xyz_DzDensPot)

      xyz_SLon = - xyz_GradDensPotLon/(xyz_DzDensPot - EPS)
      xyz_SLat = - xyz_GradDensPotLat/(xyz_DzDensPot - EPS)
!!$      xyz_SLon(:,:,0) = 0d0
!!$      xyz_SLat(:,:,0) = 0d0
!!$      xyz_SLon(:,:,kMax) = 0d0
!!$      xyz_SLat(:,:,kMax) = 0d0

    end subroutine calc_IsoNeutralSlope

    function xyz_Dz_xyz(xyz) 
      use VariableSet_mod

      real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
      real(DP) :: xyz_Dz_xyz(0:iMax-1,jMax,0:kMax)

      real(DP) :: s(0:kMax), t(0:kMax), xyt(0:iMax-1,jMax,0:tMax)
      integer :: k
     
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
      
    end function xyz_Dz_xyz

    subroutine calc_BolusVelocity(xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, &
         & xyz_SLon, xyz_SLat, xyz_totDepth, xyz_KappRho )

      ! 宣言文; Declaration statement
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_SLon, xyz_SLat, xyz_totDepth, xyz_KappRho

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_CosLat

      ! 実行文; Executable statement
      !

      xyz_CosLat = cos(xyz_Lat)
      xyz_EddInducedU = - xyz_Dz_xyz(xyz_KappRho*xyz_SLon) !xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_KappRho*xyz_SLon)))/xyz_totDepth
      xyz_EddInducedV = - xyz_Dz_xyz(xyz_KappRho*xyz_SLat*xyz_CosLat)/xyz_CosLat !xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_KappRho*xyz_SLat*xyz_CosLat)))/(xyz_totDepth*xyz_CosLat)
      xyz_EddInducedW = xyz_wz( &
           & wz_AlphaOptr_xyz(xyz_KappRho*xyz_SLon*xyz_CosLat, xyz_KappRho*xyz_SLat*xyz_CosLat) &
           & )

    end subroutine calc_BolusVelocity

    !!!!!!!!!!!!!!!!!!!!!

    subroutine SGSEddyMixing_Output()
      
      use VariableSet_mod, only: &
           & z_PTempBasic, xyz_PTempEddN, xyz_SaltN, xy_SurfHeightN, &
           & xy_totDepthBasic

      real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
           & xyz_DensPot, xyz_PTemp, xyz_RefPress, &
           & xyz_SLon, xyz_SLat, xyz_totDepth, xyz_Depth, &
           & xyz_FLon, xyz_FLat, xyz_FSig

      real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
           & xyz_BolusU, xyz_BolusV, xyz_BolusW, xyz_T

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
      xyz_T = eval_TaperingFunc(xyz_SLon, xyz_SLat, xyz_Depth)
      call check_StaticStability(xyz_T, xyz_DensPot)

!      call calc_NeutralDiffFlux( &
      call calc_SkewFlux( &
           & xyz_FLon, xyz_FLat, xyz_FSig, &
           & xyz_SaltN, wz_xyz(xyz_SaltN), xyz_SLon, xyz_SLat, xyz_T, xyz_totDepth )

      !
      call HistoryPut('Etc1', -xyz_FLat, hst_SGSEddyMix)
      call HistoryPut('Etc2', xyz_FSig, hst_SGSEddyMix)
      call HistoryPut('Etc3', &
             & -   xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*cos(xyz_Lat), xyz_FLat*cos(xyz_Lat)) ) & 
             & +  xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_FSig)))/xyz_totDepth &
             &    , hst_SGSEddyMix )
!!$      call HistoryPut('Etc1', xya_GradLat_wa(wz_xyz(xyz_DensPot))/RPlanet, hst_SGSEddyMix)
!!$      call HistoryPut('Etc2', xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_DensPot)))/xyz_totDepth, hst_SGSEddyMix)

      call HistoryPut('DensPot', xyz_DensPot, hst_SGSEddyMix)
      call HistoryPut('SlopeLon', xyz_SLon, hst_SGSEddyMix)
      call HistoryPut('SlopeLat', xyz_SLat, hst_SGSEddyMix)
      call HistoryPut('diffCoefEff', xyz_T*Kapper_rho, hst_SGSEddyMix)
      
      !
      if(SGSEddyMixType == SGSEddyMixing_GM90) then


         call calc_BolusVelocity(xyz_BolusU, xyz_BolusV, xyz_BolusW, & !(out)
              & xyz_SLon, xyz_SLat, xyz_totDepth, Kapper_rho*xyz_T )

         call HistoryPut('BolusU', xyz_BolusU, hst_SGSEddyMix)
         call HistoryPut('BolusV', xyz_BolusV, hst_SGSEddyMix)
         call HistoryPut('BolusW', xyz_BolusW, hst_SGSEddyMix)
      end if

    end subroutine SGSEddyMixing_Output

    subroutine SGSEddyMixing_PrepareOutput(OriginTime, EndTime, Intrv, FilePrefix)

      use Constants_mod, only: PI

      use GridSet_mod, only: &
           & xyz_Lon, xyz_Lat

      use SpmlUtil_mod, only: g_Sig

      ! 局所変数
      ! Local variables
      !

      !
      !
      real(DP), intent(in) :: OriginTime, EndTime, Intrv
      character(*), intent(in) :: FilePrefix

      !
      !
      character(TOKEN) :: lonName, latName, sigName, timeName
      character(TOKEN) :: dims_XYZT(4)

      ! 実行文; Executable statement
      !

      !
      OutputFlag = .true.

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

      call HistoryAddVariable( varname='Etc1', &
         & dims=dims_XYZT, longname='potential density', units='kg/m3', &
         & history=hst_SGSEddyMix, xtype='double' )
      call HistoryAddVariable( varname='Etc2', &
         & dims=dims_XYZT, longname='potential density', units='kg/m3', &
         & history=hst_SGSEddyMix, xtype='double' )
      call HistoryAddVariable( varname='Etc3', &
         & dims=dims_XYZT, longname='potential density', units='kg/m3', &
         & history=hst_SGSEddyMix, xtype='double' )


      call HistoryAddVariable( varname='DensPot', &
         & dims=dims_XYZT, longname='potential density', units='kg/m3', &
         & history=hst_SGSEddyMix, xtype='double' )

      call HistoryAddVariable( varname='SlopeLon', &
         & dims=dims_XYZT, longname='the longitude component of slope of isoneutral ', units='1', &
         & history=hst_SGSEddyMix, xtype='double' )

      call HistoryAddVariable( varname='SlopeLat', &
         & dims=dims_XYZT, longname='the latitude component of slope of isoneutral ', units='1', &
         & history=hst_SGSEddyMix, xtype='double' )

      call HistoryAddVariable( varname='diffCoefEff', &
           & dims=dims_XYZT, longname='effective diffusivity with isoneutral mixing', units='m2/s', &
           & history=hst_SGSEddyMix, xtype='double' )

      if(SGSEddyMixType == SGSEddyMixing_GM90) then

         call HistoryAddVariable( varname='BolusU', &
              & dims=dims_XYZT, longname='bolus velocity(longitude) ', units='m/s', &
              & history=hst_SGSEddyMix, xtype='double' )

         call HistoryAddVariable( varname='BolusV', &
              & dims=dims_XYZT, longname='bolus velocity(meridional) ', units='m/s', &
              & history=hst_SGSEddyMix, xtype='double' )

         call HistoryAddVariable( varname='BolusW', &
              & dims=dims_XYZT, longname='bolus velocity(vertical) ', units='m/s', &
              & history=hst_SGSEddyMix, xtype='double' )
      end if
    end subroutine SGSEddyMixing_PrepareOutput


end module SGSEddyMixing_mod

