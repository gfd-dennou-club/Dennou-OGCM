!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolver_mod 

  ! モジュール引用; Use statements
  !

  use dc_types, only: &
       & DP, STRING 

  use Constants_mod, only: &
       & Omega, Grav, RPlanet

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use TemporalIntegSet_mod, only: &
       & DelTime, SubCycleNum, &
       & nShortTimeLevel

  use EqState_JM95_mod, only: &
       & EqState_JM95_Eval


  use VariableSet_mod, only: &
       & SaltTracerID, PTempTracerID, TracerNum, &       
       & refDens

  use HydroBouEqSolverRHS_mod

!!$
!!$  use BarotModeTimeFilter_mod, only: &
!!$       & BarotModeTimeFilter_Init, BarotModeTimeFilter_Final

  use SpmlUtil_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolver_Init, HydroBouEqSolver_Final
  public :: HydroBouEqSolver_AdvanceTStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  real(DP), allocatable :: xy_Cori(:,:)

  character(*), parameter:: module_name = 'HydroBouEqSolver_mod' !< Module Name

  real(DP), parameter :: LFAM3_GAM = 1d0/12d0
  real(DP), parameter :: LFAM3_BETA = 0d0
  real(DP), parameter :: LFAM3_EPS = 0d0
  real(DP), parameter :: GFB_AM3BETA = 0.281105d0
  real(DP), parameter :: GFB_DELTA = 0.614d0
  real(DP), parameter :: GFB_EPS = 0.013d0
  real(DP), parameter :: GFB_GAM = 0.088d0
  

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Init()

    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    call HydroBouEqSolverRHS_Init()

  end subroutine HydroBouEqSolver_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Final()

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    call HydroBouEqSolverRHS_Final()

    !

  end subroutine HydroBouEqSolver_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTStep()
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & Al, Nl, Bl, DelTime

    use VariableSet_mod, only: &
         & xyz_UB, xyz_UN, xyz_UA, &
         & xyz_VB, xyz_VN, xyz_VA, &
         & xyz_PTempB, xyz_PTempN, xyz_PTempA, &
         & xyz_SaltB, xyz_SaltN, xyz_SaltA, &
         & xy_SurfHeightB, xy_SurfHeightN, xy_SurfHeightA, &
         & xy_SurfPress, xy_totDepthBasic, &
         & xyz_SigDot

    use at_module, only: at_BoundariesGrid_NN
    
    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_UrfN(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_VrfN(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_VorN(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DivN(0:iMax-1, jMax, 0:kMax)
    real(DP) :: wt_VorA(lMax, 0:tMax)
    real(DP) :: wt_DivA(lMax, 0:tMax)
    real(DP) :: wt_PTempA(lMax, 0:tMax)
    real(DP) :: wz_Psi(lMax, 0:kMax)
    real(DP) :: wz_Chi(lMax, 0:kMax)
    real(DP) :: xyz_GeoPotN(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_PressEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: wz_VorRHSN(lMax, 0:kMax)
    real(DP) :: wz_DivRHSN(lMax, 0:kMax)
    real(DP) :: wz_PTempRHSN(lMax, 0:kMax)
    real(DP) :: w_SurfHeightRHSN(lMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)


    ! 実行文; Executable statement
    !
    
    xyz_UrfN = xyz_UN*cos(xyz_Lat)
    xyz_VrfN = xyz_VN*cos(xyz_Lat)

    xyz_VorN = xyz_wz( wz_AlphaOptr_xyz(xyz_VrfN, -xyz_UrfN) )
    xyz_DivN = xyz_wz( wz_AlphaOptr_xyz(xyz_UrfN,  xyz_VrfN) )

    xy_totDepth = xy_totDepthBasic + xy_SurfHeightN
    xyz_SigDot  = diagnose_SigDot( xy_totDepth, xyz_UrfN, xyz_VrfN, xyz_DivN )
    xyz_GeoPotN  = diagnose_GeoPot( xy_totDepth, xy_SurfHeightN )

    xyz_DensEdd = 0d0!
!!$    = EqState_JM95_Eval( &
!!$         & theta=xyz_PTempN, s=xyz_SaltN, &
!!$         & p=spread(xy_SurfPress,3,kMax+1) - RefDens*xyz_GeoPotN) &
!!$         & - RefDens

    xyz_PressEdd = diagnose_PressEdd( xy_totDepth, xyz_DensEdd )

    call calc_VorEqDivEqRHS(wz_VorRHSN, wz_DivRHSN, &
         & xyz_VorN, xyz_UrfN, xyz_VrfN, xy_SurfHeightN, xyz_DensEdd, xyz_PressEdd, xyz_GeoPotN, xyz_SigDot)

    call calc_TracerEqRHS(wz_PTempRHSN, &
         & xyz_PTempN, xyz_UrfN, xyz_VrfN, xyz_DivN, xyz_SigDot )
    
    call calc_SurfHeightRHS(w_SurfHeightRHSN, &
         & xyz_UrfN, xyz_VrfN, xy_totDepth ) 

    !
    call wt_timeIntEuler_wt(wt_VorA, wt_xyz(xyz_VorN), wt_wz(wz_VorRHSN), DelTime)
    call wt_timeIntEuler_wt(wt_DivA, wt_xyz(xyz_DivN), wt_wz(wz_DivRHSN), DelTime)
    call wt_timeIntEuler_wt(wt_PTempA, wt_xyz(xyz_PTempN), wt_wz(wz_PTempRHSN), DelTime)
    call timeInt_Euler3(xy_SurfHeightA, xy_SurfHeightN, xy_w(w_SurfHeightRHSN), DelTime)

    !
!!$    call at_BoundariesGrid_NN(wt_VorA)
!!$    call at_BoundariesGrid_NN(wt_DivA)
    call at_BoundariesGrid_NN(wt_PTempA)

    wz_Psi = wz_InvLapla2D_wz( wz_wt(wt_VorA) )
    wz_Chi = wz_InvLapla2D_wz( wz_wt(wt_DivA) )

    xyz_UA = cos(xyz_Lat) * xyz_AlphaOptr_wz(wz_Chi, -wz_Psi)
    xyz_VA = cos(xyz_Lat) * xyz_AlphaOptr_wz(wz_Psi,  wz_Chi)
    xyz_PTempA = xyz_wt(wt_PTempA)

    contains
      subroutine timeInt_Euler1(xyzA, xyzN, xyzRHSN, dt)
        real(DP), intent(out) :: xyzA(0:iMax-1,jMax,0:kMax)
        real(DP), intent(in) :: xyzN(0:iMax-1,jMax,0:kMax)
        real(DP), intent(in) :: xyzRHSN(0:iMax-1,jMax,0:kMax)
        real(DP), intent(in) :: dt

        xyzA = xyzN + xyzRHSN*dt

      end subroutine timeInt_Euler1


      subroutine wt_timeIntEuler_wt(wtA, wtN, wtRHSN, dt)
        real(DP), intent(out) :: wtA(lMax,0:tMax)
        real(DP), intent(in) :: wtN(lMax,0:tMax)
        real(DP), intent(in) :: wtRHSN(lMax,0:tMax)
        real(DP), intent(in) :: dt

        wtA = wtN + wtRHSN*dt
      end subroutine wt_timeIntEuler_wt

      subroutine timeInt_Euler3(xyA, xyN, xyRHSN, dt)
        real(DP), intent(out) :: xyA(0:iMax-1,jMax)
        real(DP), intent(in) :: xyN(0:iMax-1,jMax)
        real(DP), intent(in) :: xyRHSN(0:iMax-1,jMax)
        real(DP), intent(in) :: dt

        xyA = xyN + xyRHSN*dt

      end subroutine timeInt_Euler3


  end subroutine HydroBouEqSolver_AdvanceTStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div ) result(xyz_SigDot)


    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    integer :: k
    real(DP) :: xyz_UrfHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_VrfHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DivHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xy_UrfHat(0:iMax-1,jMax)
    real(DP) :: xy_VrfHat(0:iMax-1,jMax)
    real(DP) :: xy_DivHat(0:iMax-1,jMax)
    real(DP) :: sigWeight
    real(DP) :: xy_DtotDepthDmu(0:iMax-1,jMax)
    real(DP) :: xy_DtotDepthDLambda(0:iMax-1,jMax)
    real(DP) :: sig

    xyz_UrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Urf)
    xyz_VrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Vrf)
    xyz_DivHatSig = xyz_IntSig_SigToTop_xyz(xyz_Div)
    xy_UrfHat = xy_IntSig_BtmToTop_xyz(xyz_Urf)
    xy_VrfHat = xy_IntSig_BtmToTop_xyz(xyz_Vrf)
    xy_DivHat = xy_IntSig_BtmToTop_xyz(xyz_Div)

    xy_DtotDepthDLambda = xy_w(w_DivLambda_xy(xy_totDepth))
    xy_DtotDepthDmu = xy_w(w_DivMu_xy(xy_totDepth))

    do k=1, kMax
       sig = g_Sig(k)
       xyz_SigDot(:,:,k) = &
            &   (1d0 + sig)*( &
            &           (xy_UrfHat*xy_DtotDepthDLambda + xy_VrfHat*xy_DtotDepthDmu)/xy_totDepth &
            &        +  xy_DivHat )  &
            & - xyz_DivHatSig(:,:,k) &
            & - (xyz_UrfHatSig(:,:,k)*xy_DtotDepthDLambda + xyz_VrfHatSig(:,:,k)*xy_DtotDepthDmu)/xy_totDepth 
    end do
    
  end function diagnose_SigDot

  function diagnose_PressEdd( xy_totDepth, xyz_DensEdd) result(xyz_PressEdd)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PressEdd(0:iMax-1,jMax,0:kMax)

    xyz_PressEdd = - Grav*spread(xy_totDepth,3,kMax+1)*xyz_IntSig_SigToTop_xyz(xyz_DensEdd)

  end function diagnose_PressEdd

  function diagnose_GeoPot( xy_totDepth, xy_SurfHeight ) result(xyz_GeoPot)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SurfHeight(0:iMax-1,jMax)
    real(DP) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)

    integer :: k

    do k=1, kMax
       xyz_GeoPot(:,:,k) = g_Sig(k)*xy_totDepth + xy_SurfHeight
    end do

  end function diagnose_GeoPot

end module HydroBouEqSolver_mod

