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

  use dc_message, only: &
       & MessageNotify

  use Constants_mod, only: &
       & Omega, Grav, RPlanet, &
       & hDiffCoef, vDiffCoef, &
       & RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use BoundCondSet_mod, only: &
       & KinBC_Surface, DynBC_Surface, DynBC_Bottom

  use TemporalIntegSet_mod, only: &
       & DelTime, SubCycleNum, &
       & nShortTimeLevel, &
       & timeIntMode_Euler, timeIntMode_LFTR, timeIntMode_LFAM3, &
       & timeIntMode_RK4

  use GovernEqSet_mod, only: &
       & EOSType, &
       & EOSTYPE_LINEAR, EOSTYPE_JM95

  use EOSDriver_mod, only: &
       & EOSDriver_Eval


  use VariableSet_mod, only: &
       & SaltTracerID, PTempTracerID, TracerNum

  use HydroBouEqSolverRHS_mod
  use HydroBouEqSolverVDiffProc_mod

!!$
!!$  use BarotModeTimeFilter_mod, only: &
!!$       & BarotModeTimeFilter_Init, BarotModeTimeFilter_Final

  use SpmlUtil_mod

  use TemporalIntegUtil_mod

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
    call TemporalIntegUtil_Init(iMax, jMax, kMax, lMax, tMax, DelTime)
    call HydroBouEqSolverVDiffProc_Init()

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
    call TemporalIntegUtil_Final()
    call HydroBouEqSolverVDiffProc_Final()

    !

  end subroutine HydroBouEqSolver_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTStep(timeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt)
    
    ! モジュール引用; Use statements
    !

    use VariableSet_mod, only: &
         & xyz_UB, xyz_UN, xyz_UA, &
         & xyz_VB, xyz_VN, xyz_VA, &
         & xyz_PTempEddB, xyz_PTempEddN, xyz_PTempEddA, &
         & xyz_SaltB, xyz_SaltN, xyz_SaltA, &
         & xy_SurfHeightB, xy_SurfHeightN, xy_SurfHeightA, &
         & xy_SurfPress, xy_WindStressU, xy_WindStressV, xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic
    !

    use at_module, only: &
         & at_BoundariesGrid_NN, &
         & at_BoundariesGrid_DD, &
         & at_BoundariesGrid_ND
    
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: timeIntMode
    integer, intent(in) :: nStage_BarocTimeInt
    logical, intent(in) :: isVarBUsed_BarocTimeInt

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_CosLat(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_SurfHeight(0:iMax-1, jMax)

    real(DP), dimension(0:iMax-1, jMax, 0:kMax) :: xyz_Urf, xyz_Vrf
    real(DP), dimension(lMax, 0:tMax) :: wt_Vor, wt_VorN, wt_VorB, wt_VorRKTmp
    real(DP), dimension(lMax, 0:tMax) :: wt_Div, wt_DivN, wt_DivB, wt_DivRKTmp
    real(DP), dimension(lMax, 0:tMax) :: wt_PTempEdd, wt_PTempEddN, wt_PTempEddB, wt_PTempEddRKTmp
    real(DP), dimension(lMax, 0:tMax) :: wt_Salt, wt_SaltN, wt_SaltB, wt_SaltRKTmp

    real(DP) :: wz_Psi(lMax, 0:kMax)
    real(DP) :: wz_Chi(lMax, 0:kMax)
    real(DP) :: wz_VorRHS(lMax, 0:kMax)
    real(DP) :: wz_DivRHS(lMax, 0:kMax)
    real(DP) :: wz_PTempRHS(lMax, 0:kMax)
    real(DP) :: w_SurfHeightRHS(lMax), xy_SurfHeightRKTmp(0:iMax-1,jMax)

    integer :: Stage

    ! 実行文; Executable statement
    !

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_CosLat = cos(xyz_Lat)
    xyz_Urf = xyz_UN*xyz_CosLat; xyz_Vrf = xyz_VN*xyz_CosLat
    wt_VorN = wt_wz( wz_AlphaOptr_xyz(xyz_Vrf, -xyz_Urf) )
    wt_DivN = wt_wz( wz_AlphaOptr_xyz(xyz_Urf,  xyz_Vrf) )
    wt_PTempEddN = wt_xyz(xyz_PTempEddN)
    wt_SaltN = wt_xyz(xyz_SaltN)

    if( isVarBUsed_BarocTimeInt ) then
       wt_VorB = wt_wz( wz_AlphaOptr_xyz(xyz_VB*xyz_CosLat, -xyz_UB*xyz_CosLat) )
       wt_DivB = wt_wz( wz_AlphaOptr_xyz(xyz_UB*xyz_CosLat,  xyz_VB*xyz_CosLat) )
       wt_PTempEddB = wt_xyz(xyz_PTempEddB)
       wt_SaltB = wt_xyz(xyz_SaltB)
    end if
    
    wt_Div = wt_DivN;  wt_Vor = wt_VorN
    wt_PTempEdd = wt_PTempEddN
    wt_Salt = wt_SaltN
    xy_SurfHeight = xy_SurfHeightN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do Stage=1, nStage_BarocTimeInt

       call calc_GovernEqInvisRHS(wz_VorRHS, wz_DivRHS, wz_PTempRHS, w_SurfHeightRHS, &
            & xyz_Urf, xyz_Vrf, xyz_wt(wt_Vor), xyz_wt(wt_Div), xyz_wt(wt_PTempEdd), xyz_wt(wt_Salt), &
            & xy_SurfHeight )


       if(timeIntMode == timeIntMode_LFAM3 ) then
          if(Stage==2) then
             call calc_GovernEqDiffRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
                  wz_wt(wt_VorN), wz_wt(wt_DivN), wz_wt(wt_PTempEddN) )

             !write(*,*) "Stage=", Stage
             call correct_DivEqRHSUnderRigidLid(wz_DivRHS, &
                  & xy_SurfPress, wz_wt(wt_DivN), xy_totDepthBasic, vDiffCoef, DelTime)

          else
             !write(*,*) "Stage=", Stage
             call correct_DivEqRHSUnderRigidLid(wz_DivRHS, &
                  & xy_SurfPress, wz_wt(wt_DivB), xy_totDepthBasic, vDiffCoef, DelTime*2)
          end if

       else 
          call calc_GovernEqDiffRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
               & wz_wt(wt_Vor), wz_wt(wt_Div), wz_wt(wt_PTempEdd) )
!          write(*,*) "Stage=", Stage
          call correct_DivEqRHSUnderRigidLid(wz_DivRHS, &
               & xy_SurfPress, wz_wt(wt_Div), xy_totDepthBasic, vDiffCoef, DelTime)

       end if



       select case(timeIntMode)
       case(timeIntMode_Euler)
          call timeInt_Euler()
       case(timeIntMode_RK4)
          call timeInt_RK4(Stage)
       case(timeIntMode_LFTR)
          call timeInt_LFTR(Stage)
       case(timeIntMode_LFAM3)
          call timeInt_LFAM3(Stage)
       end select

       !
       call apply_boundaryConditions(wt_Vor, wt_Div, wt_PTempEdd)

       if( Stage == nStage_BarocTimeInt ) then
          call Advance_VDiffProc( wt_Vor, wt_Div, wt_PTempEdd, &
               & xy_WindStressU, xy_WindStressV, xy_totDepthBasic+xy_SurfHeight, vDiffCoef, DelTime, &
               & DynBC_Surface, DynBC_Bottom )
       end if

       wz_Psi = wz_InvLapla2D_wz( wz_wt(wt_Vor) )
       wz_Chi = wz_InvLapla2D_wz( wz_wt(wt_Div) )
       xyz_Urf = xyz_CosLat**2 * xyz_AlphaOptr_wz(wz_Chi, -wz_Psi)
       xyz_Vrf = xyz_CosLat**2 * xyz_AlphaOptr_wz(wz_Psi,  wz_Chi)
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_UA = xyz_Urf / xyz_CosLat;  xyz_VA = xyz_Vrf / xyz_CosLat 
    xyz_PTempEddA = xyz_wt(wt_PTempEdd)
    xy_SurfHeightA = xy_SurfHeight

    contains

      subroutine timeInt_Euler()
        wt_Vor = wt_timeIntEuler( wt_VorN, wt_wz(wz_VorRHS) )
        wt_Div = wt_timeIntEuler( wt_DivN, wt_wz(wz_DivRHS) )
        wt_PTempEdd = wt_timeIntEuler( wt_PTempEddN, wt_wz(wz_PTempRHS) )
        xy_SurfHeight = xy_timeIntEuler( xy_SurfHeightN, xy_w(w_SurfHeightRHS) )
      end subroutine timeInt_Euler

      subroutine timeInt_RK4(RKStage)
        integer, intent(in) :: RKStage
        wt_Vor = wt_timeIntRK4( wt_VorN, wt_wz(wz_VorRHS), RKStage, wt_VorRKTmp )
        wt_Div = wt_timeIntRK4( wt_DivN, wt_wz(wz_DivRHS), RKStage, wt_DivRKTmp )
        wt_PTempEdd = wt_timeIntRK4( wt_PTempEddN, wt_wz(wz_PTempRHS), RKStage, wt_PTempEddRKTmp )
        xy_SurfHeight = xy_timeIntRK4( xy_SurfHeightN, xy_w(w_SurfHeightRHS), RKStage, xy_SurfHeightRKTmp )
      end subroutine timeInt_RK4

      subroutine timeInt_LFTR(Stage)
        integer, intent(in) :: Stage
        wt_Vor = wt_timeIntLFTR( wt_VorN, wt_VorB, wt_wz(wz_VorRHS), Stage )
        wt_Div = wt_timeIntLFTR( wt_DivN, wt_DivB, wt_wz(wz_DivRHS), Stage )
        wt_PTempEdd = wt_timeIntLFTR( wt_PTempEddN, wt_PTempEddB, wt_wz(wz_PTempRHS), Stage )
        xy_SurfHeight = xy_timeIntLFTR( xy_SurfHeightN, xy_SurfHeightB, xy_w(w_SurfHeightRHS), Stage )
      end subroutine timeInt_LFTR

      subroutine timeInt_LFAM3(Stage)
        integer, intent(in) :: Stage
        wt_Vor = wt_timeIntLFAM3( wt_VorN, wt_VorB, wt_wz(wz_VorRHS), Stage )
        wt_Div = wt_timeIntLFAM3( wt_DivN, wt_DivB, wt_wz(wz_DivRHS), Stage )
        wt_PTempEdd = wt_timeIntLFAM3( wt_PTempEddN, wt_PTempEddB, wt_wz(wz_PTempRHS), Stage )
        xy_SurfHeight = xy_timeIntLFAM3( xy_SurfHeightN, xy_SurfHeightB, xy_w(w_SurfHeightRHS), Stage )
      end subroutine timeInt_LFAM3

  end subroutine HydroBouEqSolver_AdvanceTStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine calc_GovernEqInvisRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, w_SurfHeightRHS, &
       & xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt, &
       & xy_SurfHeight )
    
    !
    !
    use DiagnoseUtil_mod, only: &
         & Diagnose_SigDot, Diagnose_PressBaroc, Diagnose_GeoPot

    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic, xy_SurfPress, xyz_SigDot

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP), intent(out), dimension(lMax) :: w_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_GeoPot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_PressEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)    
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)

    ! 実行文; Executable statement
    !

    xy_totDepth = xy_totDepthBasic! + xy_SurfHeightN
    xyz_SigDot  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )
    xyz_GeoPot  = Diagnose_GeoPot( xy_totDepth )

    xyz_PTemp = xyz_PTempEdd + spread(spread(z_PTempBasic,1,jMax), 1, iMax)

    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)

    xyz_PressEdd = Diagnose_PressBaroc( xy_totDepth, xyz_DensEdd )

    call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_PressEdd, xyz_GeoPot, xyz_SigDot)

    call calc_TracerEqInvisRHS(wz_PTempRHS, &
         & xyz_PTemp, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )

!!$    call calc_TracerEqInvisRHS(wz_SaltRHS, &
!!$         & xyz_Salt, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    
    call calc_SurfHeightRHS(w_SurfHeightRHS, &
         & xyz_Urf, xyz_Vrf, xy_totDepth ) 
    
  end subroutine calc_GovernEqInvisRHS

  !> @brief 
  !!
  !!
  subroutine calc_GovernEqDiffRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
       & wz_Vor, wz_Div, wz_PTempEdd )
    
    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic, xy_SurfPress, xyz_SigDot

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTempEdd


    ! 局所変数
    ! Local variables
    !
    real(DP) :: Ah
    real(DP) :: wz_totPTempEdd(lMax,0:kMax)
    integer :: k

    ! 実行文; Executable statement
    !

    call calc_VorEqDivEqDiffRHS(wz_VorRHS, wz_DivRHS, &
         & wz_Vor, wz_Div, hDiffCoef )

!!$    do k=0,kMax
!!$       wz_totPTempEdd(:,k) = w_xy( xy_w(wz_PTempEdd(:,k)) + z_PTempBasic(k) )
!!$    end do
    call calc_TracerEqDiffRHS(wz_PTempRHS, &
         & wz_PTempEdd, hDiffCoef)

  end subroutine calc_GovernEqDiffRHS

  !> @brief 
  !!
  !!
  subroutine apply_boundaryConditions( &
       & wt_Vor, wt_Div, wt_PTempEdd )
    
    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & xy_WindStressU, xy_WindStressV, xy_totDepthBasic

    use BoundCondSet_mod, only: &
         & DynBC_Surface, DynBC_Bottom, &
         & DynBCTYPE_NoSlip, DynBCTYPE_Slip

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wt_Vor(lMax, 0:tMax)
    real(DP), intent(inout) :: wt_Div(lMax, 0:tMax)
    real(DP), intent(inout) :: wt_PTempEdd(lMax, 0:tMax)

    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_CosLat(0:iMax-1,jMax)
    character :: BtmBCType_HVel

    ! 実行文; Executable statement
    !
    
    if( DynBC_Bottom == DynBCTYPE_NoSlip ) then
       BtmBCType_HVel = 'D'
    else
       BtmBCType_HVel = 'N'
    end if

    if ( DynBC_Surface == DynBCTYPE_NoSlip ) then
       xy_CosLat = cos(xyz_Lat(:,:,1))
       call apply_ZBoundaryCond( wt_Vor, &
            & 'N', BtmBCType_HVel, &
            & w_SurfBCWork = &
            &   w_AlphaOptr_xy(xy_WindStressV*xy_CosLat, -xy_WindStressU*xy_CosLat)/(RefDens*vDiffCoef) &
            & )

       call apply_ZBoundaryCond( wt_Div, &
            & 'N', BtmBCType_HVel, &
            & w_SurfBCWork = &
            &   w_AlphaOptr_xy(xy_WindStressU*xy_CosLat, xy_WindStressV*xy_CosLat)/(RefDens*vDiffCoef) &
            & )

    else
       call apply_ZBoundaryCond(wt_Vor, &
            & 'N', BtmBCType_HVel )

       call apply_ZBoundaryCond(wt_Div, &
            & 'N', BtmBCType_HVel )
    end if

    call apply_ZBoundaryCond(wt_PTempEdd, &
!!$         & 'D', 'D' )
         & 'N', 'N' )

  end subroutine apply_boundaryConditions

end module HydroBouEqSolver_mod

