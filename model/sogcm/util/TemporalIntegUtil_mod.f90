!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module TemporalIntegUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP
  use dc_message, only: MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: TemporalIntegUtil_Init, TemporalIntegUtil_Final
  public :: TemporalIntegUtil_SetDelTime
  public :: TemporalIntegUtil_GetDDtCoef
  public :: TemporalIntegUtil_getInfo
  
  ! * For Euler
  public :: xyz_timeIntEuler, xy_timeIntEuler, wt_timeIntEuler

  ! * For RK2
  public :: xyz_timeIntRK2, xy_timeIntRK2, wt_timeIntRK2

  ! * For RK4
  public :: xyz_timeIntRK4, xy_timeIntRK4, wt_timeIntRK4

  ! * For LFTR
  public :: xyz_timeIntLFTR, xy_timeIntLFTR, wt_timeIntLFTR

  ! * For LFAM3
  interface xy_timeIntLFAM3
     module procedure xy_timeIntLFAM3_EX
     module procedure xy_timeIntLFAM3_IMEX
  end interface xy_timeIntLFAM3

  interface wt_timeIntLFAM3
     module procedure wt_timeIntLFAM3_EX
     module procedure wt_timeIntLFAM3_IMEX
  end interface wt_timeIntLFAM3

  public :: xyz_timeIntLFAM3, xy_timeIntLFAM3, wt_timeIntLFAM3

  ! 公開変数
  ! Public variable
  !
  integer, parameter, public :: TimeIntMode_Euler = 1
  character(*), parameter, public :: TimeIntModeLBL_Euler = 'TimeIntMode_Eluer' 

  integer, parameter, public :: TimeIntMode_RK4 = 2
  character(*), parameter, public :: TimeIntModeLBL_RK4 = 'TimeIntMode_RK4' 

  integer, parameter, public :: TimeIntMode_LFTR = 3
  character(*), parameter, public :: TimeIntModeLBL_LFTR = 'TimeIntMode_LFTR' 

  integer, parameter, public :: TimeIntMode_LFAM3 = 4
  character(*), parameter, public :: TimeIntModeLBL_LFAM3 = 'TimeIntMode_LFAM3' 

  integer, parameter, public :: TimeIntMode_RK2 = 5
  character(*), parameter, public :: TimeIntModeLBL_RK2 = 'TimeIntMode_RK2' 

  integer, parameter, public :: TimeIntMode_RK3 = 6
  character(*), parameter, public :: TimeIntModeLBL_RK3 = 'TimeIntMode_RK3' 

  integer, parameter, public :: TimeIntMode_PC23_AB2AM3CR = 7
  character(*), parameter, public :: TimeIntModeLBL_PC23_AB2AM3CR = 'TimeIntMode_PC23_AB2AM3CR' 


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'TemporalIntegUtil_mod' !< Module Name
  integer :: iMax, jMax, kMax, lMax, tMax

  real(DP) :: dt
  real(DP), parameter  :: Euler_coef(1) = (/ 1d0 /)
  real(DP), parameter :: RK2_coef(2) = (/ 0.5d0, 0.5d0 /)
  real(DP), parameter :: RK3_coef(3) = (/ 0.5d0, 0.5d0, 0.5d0 /)
  real(DP), parameter  :: RK4_coef(4) = (/ 1d0, 0.5d0, 0.5d0, 1d0 /)
  real(DP), parameter  :: LFTR_coef(2) = (/ 2d0, 1d0 /)
  real(DP), parameter  :: LFAM3_coef(2) = (/ 2d0, 1d0 /)
  real(DP), parameter  :: PC23_AB2AM3CR_coef(2) = (/ 1d0, 1d0 /)

contains

  !>
  !!
  !!
  subroutine TemporalIntegUtil_Init( & 
       & im, jm, km, lm, tm, delTime )
    
    ! 
    !
    integer, intent(in) :: im, jm, km, lm, tm
    real(DP), intent(in) :: DelTime

    ! 実行文; Executable statements
    !
    iMax = im; jMax = jm; kMax = km
    lMax = lm; tMax = tm
    dt = DelTime

  end subroutine TemporalIntegUtil_Init

  !>
  !!
  !!
  subroutine TemporalIntegUtil_Final()

    ! 実行文; Executable statements
    !

  end subroutine TemporalIntegUtil_Final

  !> @brief 
  !!
  !!
  function TemporalIntegUtil_getInfo( tIntModeLabel, & !(in)
       & tIntModeID, is_VarB_Used, nStage ) result(is_Registered)
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: tIntModeLabel
    integer, intent(out) :: tIntModeID
    logical, intent(out) :: is_VarB_Used
    integer, intent(out) :: nStage
    logical :: is_Registered

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    is_Registered = .true.
    select case(tIntModeLabel)
       case(TimeIntModeLBL_Euler)
          tIntModeID = TimeIntMode_Euler
          is_VarB_Used = .false.
          nStage = 1
       case(TimeIntModeLBL_LFTR)
          tIntModeID = timeIntMode_LFTR
          is_VarB_Used = .true.
          nStage = 2
       case(TimeIntModeLBL_LFAM3)
          tIntModeID = timeIntMode_LFAM3
          is_VarB_Used = .true.
          nStage = 2
       case(TimeIntModeLBL_RK2)
          tIntModeID = timeIntMode_RK2
          is_VarB_Used = .false.
          nStage = 2
       case(TimeIntModeLBL_RK4)
          tIntModeID = timeIntMode_RK4
          is_VarB_Used = .false.
          nStage = 4
       case(TimeIntModeLBL_PC23_AB2AM3CR)
          tIntModeID = TimeIntMode_PC23_AB2AM3CR
          is_VarB_Used = .true.
          nStage = 2
       case default
          call MessageNotify( "W", module_name, &
               & "Specified name of temporal integration method '%a' is invalid", ca=(/tIntModeLabel/) )  
          is_Registered = .false.
    end select
    
  end function TemporalIntegUtil_getInfo


  !> @brief 
  !!
  !!
  subroutine TemporalIntegUtil_SetDelTime(newDelTime)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: newDelTime
    
    
    ! 実行文; Executable statement
    !
    dt = newDelTime

  end subroutine TemporalIntegUtil_SetDelTime

  !> @brief 
  !!
  !!
  function TemporalIntegUtil_GetDDtCoef(timeIntMode, stage) result(coef)
    
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: timeIntMode
    integer, intent(in) :: stage
    real(DP) :: coef
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    select case(timeIntMode)
    case(timeIntMode_Euler)
       coef = Euler_coef(stage)
    case(timeIntMode_RK2)
       coef = RK2_coef(stage)
    case(timeIntMode_RK4)
       coef = RK4_coef(stage)
    case(timeIntMode_LFTR)
       coef = LFTR_coef(stage)
    case(timeIntMode_LFAM3)
       coef = LFAM3_coef(stage)
    case(TimeIntMode_PC23_AB2AM3CR)
       coef = PC23_AB2AM3CR_coef(stage)
    end select

  end function TemporalIntegUtil_GetDDtCoef


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  !*** Euler method

  function xyz_timeIntEuler(xyzN, xyzRHSN) result(xyzA)
    real(DP), intent(in) :: xyzN(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyzRHSN(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyzA(0:iMax-1,jMax,0:kMax)
    
    xyzA = xyzN + xyzRHSN*dt

  end function xyz_timeIntEuler
  
  function wt_timeIntEuler(wtN, wtRHSN) result(wtA) 
    real(DP), intent(in) :: wtN(lMax,0:tMax)
    real(DP), intent(in) :: wtRHSN(lMax,0:tMax)
    real(DP) :: wtA(lMax,0:tMax)

    wtA = wtN + wtRHSN*dt
  end function wt_timeIntEuler

  function xy_timeIntEuler(xyN, xyRHSN) result(xyA)
    real(DP), intent(in) :: xyN(0:iMax-1,jMax)
    real(DP), intent(in) :: xyRHSN(0:iMax-1,jMax)
    real(DP) :: xyA(0:iMax-1, jMax)

    xyA = xyN + xyRHSN*dt

  end function xy_timeIntEuler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  !*** Runge=Kutta 2nd order method

  function xyz_timeIntRK2(xyzN0, xyzRHS, rkStage, rkTmp) result(xyzA)
    real(DP), intent(in) :: xyzN0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyzRHS(0:iMax-1,jMax,0:kMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyzA(0:iMax-1,jMax,0:kMax)
    
    if(rkStage==1) rkTmp = xyzN0

    if(rkStage==2) then
       xyzA = rkTmp + dt*0.5d0*xyzRHS
       return
    end if

    rkTmp = rkTmp + dt*0.5d0*xyzRHS
    xyzA = xyzN0 + dt*xyzRHS

  end function xyz_timeIntRK2

  function wt_timeIntRK2(wtN0, wtRHS, rkStage, rkTmp) result(wtA)
    real(DP), intent(in) :: wtN0(lMax,0:tMax)
    real(DP), intent(in) :: wtRHS(lMax,0:tMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(lMax,0:tMax)
    real(DP) :: wtA(lMax,0:tMax)

    if(rkStage==1) rkTmp = wtN0

    if(rkStage==2) then
       wtA = rkTmp + dt*0.5d0*wtRHS
       return
    end if

    rkTmp = rkTmp + dt*0.5d0*wtRHS
    wtA = wtN0 + dt*wtRHS

  end function wt_timeIntRK2

  function xy_timeIntRK2(xyN0, xyRHS, rkStage, rkTmp) result(xyA)
    real(DP), intent(in) :: xyN0(0:iMax-1,jMax)
    real(DP), intent(in) :: xyRHS(0:iMax-1,jMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(0:iMax-1,jMax)
    real(DP) :: xyA(0:iMax-1,jMax)
    
    if(rkStage==1) rkTmp = xyN0

    if(rkStage==2) then
       xyA = rkTmp + dt*0.5d0*xyRHS
       return
    end if

    rkTmp = rkTmp + dt*0.5d0*xyRHS
    xyA = xyN0 + dt*xyRHS

  end function xy_timeIntRK2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  !*** Runge=Kutta 4th order method

  function xyz_timeIntRK4(xyzN0, xyzRHS, rkStage, rkTmp) result(xyzA)
    real(DP), intent(in) :: xyzN0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyzRHS(0:iMax-1,jMax,0:kMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyzA(0:iMax-1,jMax,0:kMax)
    
    if(rkStage==1) rkTmp = xyzN0

    if(rkStage==4) then
       xyzA = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyzRHS
       return
    end if

    rkTmp = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyzRHS
    xyzA = xyzN0 + dt*RK4_coef(rkStage+1)*xyzRHS

  end function xyz_timeIntRK4

  function wt_timeIntRK4(wtN0, wtRHS, rkStage, rkTmp) result(wtA)
    real(DP), intent(in) :: wtN0(lMax,0:tMax)
    real(DP), intent(in) :: wtRHS(lMax,0:tMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(lMax,0:tMax)
    real(DP) :: wtA(lMax,0:tMax)
    
    if(rkStage==1) rkTmp = wtN0

    if(rkStage==4) then
       !$omp parallel workshare
       wtA = rkTmp + dt/(6d0*RK4_coef(rkStage))*wtRHS
       !$omp end parallel workshare
       return
    end if

    !$omp parallel workshare
    rkTmp = rkTmp + dt/(6d0*RK4_coef(rkStage))*wtRHS
    wtA = wtN0 + dt*RK4_coef(rkStage+1)*wtRHS
    !$omp end parallel workshare

  end function wt_timeIntRK4

  function xy_timeIntRK4(xyN0, xyRHS, rkStage, rkTmp) result(xyA)
    real(DP), intent(in) :: xyN0(0:iMax-1,jMax)
    real(DP), intent(in) :: xyRHS(0:iMax-1,jMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(0:iMax-1,jMax)
    real(DP) :: xyA(0:iMax-1,jMax)
    
    if(rkStage==1) rkTmp = xyN0

    if(rkStage==4) then
       !$omp parallel workshare
       xyA = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyRHS
       !$omp end parallel workshare
       return
    end if

    !$omp parallel workshare
    rkTmp = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyRHS
    xyA = xyN0 + dt*RK4_coef(rkStage+1)*xyRHS
    !$omp end parallel workshare

  end function xy_timeIntRK4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  !*** Leap Frog Trapezoidal Rule(LF-TR)

  function xyz_timeIntLFTR( valN, valB, RHS, stage) result(val)
    real(DP), intent(in) :: valN(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: valB(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: RHS(0:iMax-1,jMax,0:kMax)
    integer, intent(in) :: stage
    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFTR", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       val = valB + 2d0*dt*RHS
       val = 0.5d0*(val + valN)
    case(2)
       val = valN + dt*RHS
    end select

  end function xyz_timeIntLFTR

  function wt_timeIntLFTR( valN, valB, RHS, stage) result(val)
    real(DP), intent(in) :: valN(lMax,0:tMax)
    real(DP), intent(in) :: valB(lMax,0:tMax)
    real(DP), intent(in) :: RHS(lMax,0:tMax)
    integer, intent(in) :: stage
    real(DP) :: val(lMax,0:tMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFTR", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       val = valB + 2d0*dt*RHS
       val = 0.5d0*(val + valN)
    case(2)
       val = valN + dt*RHS
    end select

  end function wt_timeIntLFTR


  function xy_timeIntLFTR( valN, valB, RHS, stage) result(val)
    real(DP), intent(in) :: valN(0:iMax-1,jMax)
    real(DP), intent(in) :: valB(0:iMax-1,jMax)
    real(DP), intent(in) :: RHS(0:iMax-1,jMax)
    integer, intent(in) :: stage
    real(DP) :: val(0:iMax-1,jMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFTR", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       val = valB + 2d0*dt*RHS
       val = 0.5d0*(val + valN)
    case(2)
       val = valN + dt*RHS
    end select

  end function xy_timeIntLFTR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  !*** Leap Frog with 3rd order Adams-Bashforth extrapolation (LF-AM3)

  function xyz_timeIntLFAM3( valN, valB, RHS, stage ) result(val)
    real(DP), intent(in) :: valN(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: valB(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: RHS(0:iMax-1,jMax,0:kMax)
    integer, intent(in) :: stage
    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFAM3", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       !$omp parallel workshare
       val = ( 5d0*(valB + 2d0*dt*RHS) + 8d0*valN - valB )/12d0
       !$omp end parallel workshare
    case(2)
       !$omp parallel workshare
       val = valN + dt*RHS
       !$omp end parallel workshare
    end select

  end function xyz_timeIntLFAM3


  function wt_timeIntLFAM3_EX( valN, valB, RHS, stage) result(val)
    real(DP), intent(in) :: valN(lMax,0:tMax)
    real(DP), intent(in) :: valB(lMax,0:tMax)
    real(DP), intent(in) :: RHS(lMax,0:tMax)
    integer, intent(in) :: stage
    real(DP) :: val(lMax,0:tMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFAM3", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       !$omp parallel workshare
       val = ( 5d0*(valB + 2d0*dt*RHS) + 8d0*valN - valB )/12d0
       !$omp end parallel workshare
    case(2)
       !$omp parallel workshare
       val = valN + dt*RHS
       !$omp end parallel workshare
    end select

  end function wt_timeIntLFAM3_EX


  function xy_timeIntLFAM3_EX( valN, valB, RHS, stage) result(val)
    real(DP), intent(in) :: valN(0:iMax-1,jMax)
    real(DP), intent(in) :: valB(0:iMax-1,jMax)
    real(DP), intent(in) :: RHS(0:iMax-1,jMax)
    integer, intent(in) :: stage
    real(DP) :: val(0:iMax-1,jMax)
    
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFAM3", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       !$omp parallel workshare
       val = ( 5d0*(valB + 2d0*dt*RHS) + 8d0*valN - valB )/12d0
       !$omp end parallel workshare
    case(2)
       !$omp parallel workshare
       val = valN + dt*RHS
       !$omp end parallel workshare
    end select

  end function xy_timeIntLFAM3_EX

  function wt_timeIntLFAM3_IMEX( valN, valB, RHS, stage, implicitSolver) result(val)
    real(DP), intent(in) :: valN(lMax,0:tMax)
    real(DP), intent(in) :: valB(lMax,0:tMax)
    real(DP), intent(in) :: RHS(lMax,0:tMax)
    integer, intent(in) :: stage
    real(DP) :: val(lMax,0:tMax)
    
    interface
       ! Call a solver for the homogeneous problem, A (val**) = (val*). 
       function implicitSolver(wt_val) result(wt_ret)
         use dc_types, only: DP
         use GridSet_mod, only: lMax, tMax
         real(DP), intent(in) :: wt_val(lMax,0:tMax)
         real(DP) :: wt_ret(lMax,0:tMax)
       end function implicitSolver
    end interface
    
    real(DP) :: wt_valImpl(lMax,0:tMax)

#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFAM3", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       wt_valImpl(:,:) = implicitSolver(valB + 2d0*dt*RHS)
       !$omp parallel workshare
       val = ( 5d0*wt_valImpl  + 8d0*valN - valB )/12d0
       !$omp end parallel workshare
    case(2)
       val = implicitSolver(valN + dt*RHS)
    end select

  end function wt_timeIntLFAM3_IMEX


  function xy_timeIntLFAM3_IMEX( valN, valB, RHS, stage, implicitSolver) result(val)
    real(DP), intent(in) :: valN(0:iMax-1,jMax)
    real(DP), intent(in) :: valB(0:iMax-1,jMax)
    real(DP), intent(in) :: RHS(0:iMax-1,jMax)
    integer, intent(in) :: stage
    real(DP) :: val(0:iMax-1,jMax)

    interface
       ! Call a solver for the homogeneous problem, A (val**) = (val*). 
       function implicitSolver(xy_val) result(xy_ret)
         use dc_types, only: DP
         use GridSet_mod, only: iMax, jMax
         real(DP), intent(in) :: xy_val(0:iMax-1,jMax)
         real(DP) :: xy_ret(0:iMax-1,jMax)
       end function implicitSolver
    end interface
    
    real(DP) :: xy_valImpl(0:iMax-1,jMax)
#ifdef DEBUG
    if(stage > 2 .or. stage < 1) call MessageNotify("E", module_name//"::LFAM3", "The number of stage is invalid") 
#endif

    select case(stage)
    case(1)
       xy_valImpl = implicitSolver(valB + 2d0*dt*RHS)
       !$omp parallel workshare
       val = ( 5d0*xy_valImpl + 8d0*valN - valB )/12d0
       !$omp end parallel workshare
    case(2)
       val = implicitSolver(valN + dt*RHS)
    end select

  end function xy_timeIntLFAM3_IMEX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module TemporalIntegUtil_mod
