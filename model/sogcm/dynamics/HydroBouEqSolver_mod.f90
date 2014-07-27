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
       & DP, TOKEN, STRING 

  use dc_message, only: &
       & MessageNotify

  use Constants_mod, only: &
       & Omega, Grav, RPlanet, &
       & hViscCoef, vViscCoef, &
       & hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef, &
       & RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use BoundCondSet_mod, only: &
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom 

  use TemporalIntegSet_mod, only: &
       & CurrentTimeStep, SubCycleNum, &
       & nShortTimeLevel, &
       & timeIntMode_Euler, timeIntMode_LFTR, timeIntMode_LFAM3, &
       & timeIntMode_RK2, timeIntMode_RK4

  use GovernEqSet_mod, only: &
       & EOSType

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
  public :: apply_boundaryConditions

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
    call TemporalIntegUtil_Init( iMax, jMax, kMax, lMax, tMax, 0d0 )
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
  subroutine HydroBouEqSolver_AdvanceTStep(DelTime, timeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt)
    
    ! モジュール引用; Use statements
    !

    use VariableSet_mod, only: &
         & xyz_UB, xyz_UN, xyz_UA, &
         & xyz_VB, xyz_VN, xyz_VA, &
         & xyz_PTempEddB, xyz_PTempEddN, xyz_PTempEddA, &
         & xyz_SaltB, xyz_SaltN, xyz_SaltA, &
         & xy_SurfHeightB, xy_SurfHeightN, xy_SurfHeightA, &
         & xy_SurfPressB, xy_SurfPressN, xy_SurfPressA, &
         & xy_WindStressU, xy_WindStressV, xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic
    !


    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTime
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

    real(DP), dimension(lMax, 0:kMax)  :: wz_Psi, wz_Chi
    real(DP), dimension(lMax, 0:kMax) ::  wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP) :: w_SurfHeightRHS(lMax), xy_SurfHeightRKTmp(0:iMax-1,jMax)

    integer :: Stage
    character(TOKEN) :: TIntType_SurfPressTerm

real(DP) :: wt_Tmp(lMax,0:tMax), xy_SurfPress(0:iMax-1,jMax)


    ! 実行文; Executable statement
    !

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_CosLat = cos(xyz_Lat)
    xyz_Urf = xyz_UN*xyz_CosLat; xyz_Vrf = xyz_VN*xyz_CosLat
    wt_VorN = wt_wz( wz_AlphaOptr_xyz(xyz_Vrf, -xyz_Urf) )
    wt_DivN = wt_wz( wz_AlphaOptr_xyz(xyz_Urf,  xyz_Vrf) )
    wt_PTempEddN = wt_xyz(xyz_PTempEddN)
    wt_SaltN = wt_xyz(xyz_SaltN)
    xy_SurfPress = xy_SurfPressN

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

    call TemporalIntegUtil_SetDelTime(DelTime)

    do Stage=1, nStage_BarocTimeInt

       call calc_GovernEqInvisRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, w_SurfHeightRHS, &                 ! (inout)
            & xyz_Urf, xyz_Vrf, xyz_wt(wt_Vor), xyz_wt(wt_Div), xyz_wt(wt_PTempEdd), xyz_wt(wt_Salt), &  ! (in)
            & xy_SurfHeight, xy_SurfPress                                                             & ! (in)
            & )


       call calc_GovernEqHViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS,      & ! (inout)  
            & wz_wt(wt_Vor), wz_wt(wt_Div), wz_wt(wt_PTempEdd),            & ! (in)
            & hViscCoef, hHyperViscCoef, hDiffCoef                         & ! (in)
            & )

       TIntType_SurfPressTerm = "CRANKNIC"

       if(timeIntMode == timeIntMode_LFAM3 ) then

          if(Stage==2) then

             call calc_GovernEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS,         & ! (inout)  
                  & wz_wt(wt_VorN), wz_wt(wt_DivN), wz_wt(wt_PTempEddN),            & ! (in)
                  & 0.5d0*vViscCoef, vHyperViscCoef, 0.5d0*vDiffCoef                & ! (in)
                  & )

          else

             call calc_GovernEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS,         & ! (inout)  
                  & wz_wt(wt_VorN), wz_wt(wt_DivN), wz_wt(wt_PTempEddN),            & ! (in)
                  & vViscCoef, vHyperViscCoef, vDiffCoef                            & ! (in)
                  & )

             TIntType_SurfPressTerm = "CRANKNIC_WithLF"
             call correct_DivEqRHSUnderRigidLid( wz_DivRHS, xy_SurfPressA, &
                  & wz_wt(wt_DivB), xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
                  & 1d0 / (TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime), &
                  & TIntType_SurfPressTerm )
             xy_SurfPress = 0.5d0*xy_SurfPressN

          end if

       else 

          call calc_GovernEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS,         & ! (inout)  
               & wz_wt(wt_Vor), wz_wt(wt_Div), wz_wt(wt_PTempEdd),               & ! (in)
               & vViscCoef, vHyperViscCoef, vDiffCoef                            & ! (in)
               & )

          call correct_DivEqRHSUnderRigidLid( wz_DivRHS, xy_SurfPressA, &
               & wz_wt(wt_DivN), xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
               & 1d0 / (TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime), &
               & TIntType_SurfPressTerm )
          xy_SurfPress = xy_SurfPressN

       end if


       select case(timeIntMode)
       case(timeIntMode_Euler)
          call timeInt_Euler()
       case(timeIntMode_RK2)
          call timeInt_RK2(Stage)
       case(timeIntMode_RK4)
          call timeInt_RK4(Stage)
       case(timeIntMode_LFTR)
          call timeInt_LFTR(Stage)
       case(timeIntMode_LFAM3)
          call timeInt_LFAM3(Stage)
       end select

       if( timeIntMode == timeIntMode_LFAM3 ) then
          if ( Stage /= nStage_BarocTimeInt ) then
          else

             xy_SurfPressA = &
                  & 2d0*RefDens/DelTime * xy_w( w_InvLapla2D_w(w_IntSig_BtmToTop_wz(wz_wt(wt_Div))) )

             call Advance_VDiffProc( wt_Vor, wt_Div, wt_PTempEdd, &  !(inout)
                  & xy_WindStressU, xy_WindStressV, xy_totDepthBasic+xy_SurfHeight, & 
                  & 0.5d0*vViscCoef, 0.5d0*vDiffCoef, & 
                  & vDiffCoef, DelTime, &
                  & DynBC_Surface, DynBC_Bottom )

!!$             call correct_DivEqRHSUnderRigidLid2( wt_Div, xy_SurfPressA, &
!!$                  & wz_wt(wt_Div), xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
!!$                  & 1d0 / (TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime), &
!!$                  & 'CRANKNIC' )

!!$
!!$             xy_SurfPressA = xy_SurfPressA &
!!$             if( mod(CurrentTimeStep, 10) == 0) then
!!$write(*,*) "-------------"
                  xy_SurfPressA = xy_SurfPressA + &
                  &  2d0*RefDens * xy_w( w_InvLapla2D_w( 0.5d0*vViscCoef * &
                  &     w_IntSig_BtmToTop_wz(wz_wt(wt_DSigDSig_wt(wt_Div))) &
                  & ) )/xy_totDepthBasic**2
!!$               end if
          end if
       end if

       call apply_boundaryConditions(wt_Vor, wt_Div, wt_PTempEdd)


       wz_Psi = wz_InvLapla2D_wz( wz_wt(wt_Vor) )
       wz_Chi = wz_InvLapla2D_wz( wz_wt(wt_Div) )
       xyz_Urf = xyz_CosLat**2 * xyz_AlphaOptr_wz(wz_Chi, -wz_Psi)
       xyz_Vrf = xyz_CosLat**2 * xyz_AlphaOptr_wz(wz_Psi,  wz_Chi)

    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_UA = xyz_Urf / xyz_CosLat;  xyz_VA = xyz_Vrf / xyz_CosLat 
    xyz_PTempEddA = xyz_wt(wt_PTempEdd)
!    xy_SurfPressA = xy_SurfPress
    xy_SurfHeightA = xy_SurfHeight

    contains

      subroutine timeInt_Euler()
        wt_Vor = wt_timeIntEuler( wt_VorN, wt_wz(wz_VorRHS) )
        wt_Div = wt_timeIntEuler( wt_DivN, wt_wz(wz_DivRHS) )
        wt_PTempEdd = wt_timeIntEuler( wt_PTempEddN, wt_wz(wz_PTempRHS) )
        xy_SurfHeight = xy_timeIntEuler( xy_SurfHeightN, xy_w(w_SurfHeightRHS) )
      end subroutine timeInt_Euler

      subroutine timeInt_RK2(RKStage)
        integer, intent(in) :: RKStage
        wt_Vor = wt_timeIntRK2( wt_VorN, wt_wz(wz_VorRHS), RKStage, wt_VorRKTmp )
        wt_Div = wt_timeIntRK2( wt_DivN, wt_wz(wz_DivRHS), RKStage, wt_DivRKTmp )
        wt_PTempEdd = wt_timeIntRK2( wt_PTempEddN, wt_wz(wz_PTempRHS), RKStage, wt_PTempEddRKTmp )
        xy_SurfHeight = xy_timeIntRK2( xy_SurfHeightN, xy_w(w_SurfHeightRHS), RKStage, xy_SurfHeightRKTmp )
      end subroutine timeInt_RK2

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
       & xy_SurfHeight, xy_SurfPress )

    ! モジュール引用; Use statements
    !
    use DiagnoseUtil_mod, only: &
         & Diagnose_SigDot, Diagnose_PressBaroc, Diagnose_GeoPot

    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic, xyz_SigDot

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP), intent(out), dimension(lMax) :: w_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight, xy_SurfPress

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_GeoPot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_Press(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)    
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)
    integer :: i, j

    ! 実行文; Executable statement
    !

    xy_totDepth = xy_totDepthBasic! + xy_SurfHeightN
    xyz_SigDot  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )
    xyz_GeoPot  = Diagnose_GeoPot( xy_totDepth )

    forAll(i=0:iMax-1,j=1:jMax) &
         & xyz_PTemp(i,j,:) = xyz_PTempEdd(i,j,:) + z_PTempBasic(:)

    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)

    ! Calculate the pressure which is deviation from -RefDens*Grav*z).
    !
    xyz_Press(:,:,:) =   spread(xy_SurfPress, 3, kMax+1)                 &    ! barotropic component
                &      + Diagnose_PressBaroc( xy_totDepth, xyz_DensEdd )      ! baroclinic component

    !
    !
    call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot)

    Call calc_TracerEqInvisRHS(wz_PTempRHS, &
         & xyz_PTemp, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )

!!$    call calc_TracerEqInvisRHS(wz_SaltRHS, &
!!$         & xyz_Salt, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    
    call calc_SurfHeightRHS(w_SurfHeightRHS, &
         & xyz_Urf, xyz_Vrf, xy_totDepth ) 
    
  end subroutine calc_GovernEqInvisRHS

  !> @brief 
  !!
  !!
  subroutine calc_GovernEqHViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
       & wz_Vor, wz_Div, wz_PTempEdd, &
       & hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef  )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTempEdd
    real(DP), intent(in) :: hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef

    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !

    call calc_HDiffRHS(wz_VorRHS,                         &  !(inout)
         & wz_Vor, hViscTermCoef, hHyperViscTermCoef, 'V' &  !(in)
         & )

    call calc_HDiffRHS(wz_DivRHS,                               &  !(inout)
         & wz_Div, hViscTermCoef, hHyperViscTermCoef, 'V', 2d0  &  !(in)
         & )

    call calc_HDiffRHS(wz_PTempRHS,                    &  !(inout)
         & wz_PTempEdd, hDiffTermCoef, 0d0,  'S'       &  !(in)
         & )


  end subroutine calc_GovernEqHViscRHS


  !> @brief 
  !!
  !!
  subroutine calc_GovernEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, &
       & wz_Vor, wz_Div, wz_PTempEdd, &
       & vViscTermCoef, vHyperViscTermCoef, vDiffTermCoef  )
    
    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTempEdd
    real(DP), intent(in) :: vViscTermCoef, vHyperViscTermCoef, vDiffTermCoef

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !


    xyz_totDepth = spread(xy_totDepthBasic, 3, kMax+1)

    call calc_VDiffRHS(wz_VorRHS,                                   &  !(inout)
         & wz_Vor, vViscTermCoef, vHyperViscTermCoef, xyz_totDepth  &  !(in)
         & )

    call calc_VDiffRHS(wz_DivRHS,                                   &  !(inout)
         & wz_Div, vViscTermCoef, vHyperViscTermCoef, xyz_totDepth  &  !(in)
         & )

    call calc_VDiffRHS(wz_PTempRHS,                              &  !(inout)
         & wz_PTempEdd, vDiffTermCoef, 0d0, xyz_totDepth         &  !(in)
         & )

  end subroutine calc_GovernEqVViscRHS


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
         & DynBCTYPE_NoSlip, DynBCTYPE_Slip, &
         & ThermBCTYPE_Adiabat

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wt_Vor(lMax, 0:tMax)
    real(DP), intent(inout) :: wt_Div(lMax, 0:tMax)
    real(DP), intent(inout) :: wt_PTempEdd(lMax, 0:tMax)

    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_Coef(0:iMax-1,jMax)
    character :: BtmBCType_HVel
    character :: BtmBCType_Therm, SurfBCType_Therm

    ! 実行文; Executable statement
    !
    
    if( DynBC_Bottom == DynBCTYPE_NoSlip ) then
       BtmBCType_HVel = 'D'
    else
       BtmBCType_HVel = 'N'
    end if

    if ( DynBC_Surface == DynBCTYPE_NoSlip ) then
       xy_Coef = xy_totDepthBasic*cos(xyz_Lat(:,:,1))
       call apply_ZBoundaryCond( wt_Vor, &
            & 'N', BtmBCType_HVel, &
            & w_SurfBCWork = &
            &   w_AlphaOptr_xy(xy_WindStressV*xy_Coef, -xy_WindStressU*xy_Coef)/(RefDens*vDiffCoef) &
            & )

       call apply_ZBoundaryCond( wt_Div, &
            & 'N', BtmBCType_HVel, &
            & w_SurfBCWork = &
            &   w_AlphaOptr_xy(xy_WindStressU*xy_Coef, xy_WindStressV*xy_Coef)/(RefDens*vDiffCoef) &
            & )

    else
       call apply_ZBoundaryCond(wt_Vor, &
            & 'N', BtmBCType_HVel )

       call apply_ZBoundaryCond(wt_Div, &
            & 'N', BtmBCType_HVel )
    end if

    if ( ThermBC_Surface == ThermBCTYPE_Adiabat ) then
       SurfBCType_Therm = 'N'
    end if

    if ( ThermBC_Bottom == ThermBCTYPE_Adiabat ) then
       BtmBCType_Therm = 'N'
    end if

    call apply_ZBoundaryCond(wt_PTempEdd, &
            & SurfBCType_Therm, BtmBCType_Therm )

  end subroutine apply_boundaryConditions

end module HydroBouEqSolver_mod

