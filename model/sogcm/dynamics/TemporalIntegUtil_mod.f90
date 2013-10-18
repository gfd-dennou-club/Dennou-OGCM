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

  public :: xyz_timeIntEuler, xy_timeIntEuler, wt_timeIntEuler
  public :: xyz_timeIntRK4, xy_timeIntRK4, wt_timeIntRK4
  public :: xyz_timeIntLFTR, xy_timeIntLFTR, wt_timeIntLFTR
  public :: xyz_timeIntLFAM3, xy_timeIntLFAM3, wt_timeIntLFAM3

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'TemporalIntegUtil_mod' !< Module Name
  integer :: iMax, jMax, kMax, lMax, tMax

  real(DP) :: dt
  real(DP) :: RK4_coef(4)

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

    RK4_coef(:) = (/ 1d0, 0.5d0, 0.5d0, 1d0 /)

  end subroutine TemporalIntegUtil_Init

  !>
  !!
  !!
  subroutine TemporalIntegUtil_Final()

    ! 実行文; Executable statements
    !

  end subroutine TemporalIntegUtil_Final

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
       wtA = rkTmp + dt/(6d0*RK4_coef(rkStage))*wtRHS
       return
    end if

    rkTmp = rkTmp + dt/(6d0*RK4_coef(rkStage))*wtRHS
    wtA = wtN0 + dt*RK4_coef(rkStage+1)*wtRHS

  end function wt_timeIntRK4

  function xy_timeIntRK4(xyN0, xyRHS, rkStage, rkTmp) result(xyA)
    real(DP), intent(in) :: xyN0(0:iMax-1,jMax)
    real(DP), intent(in) :: xyRHS(0:iMax-1,jMax)
    integer, intent(in) :: rkStage
    real(DP), intent(inout) :: rkTmp(0:iMax-1,jMax)
    real(DP) :: xyA(0:iMax-1,jMax)
    
    if(rkStage==1) rkTmp = xyN0

    if(rkStage==4) then
       xyA = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyRHS
       return
    end if

    rkTmp = rkTmp + dt/(6d0*RK4_coef(rkStage))*xyRHS
    xyA = xyN0 + dt*RK4_coef(rkStage+1)*xyRHS

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

  function xyz_timeIntLFAM3( valN, valB, RHS, stage) result(val)
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
       val = ( 5d0*val + 8d0*valN - valB )/12d0
    case(2)
       val = valN + dt*RHS
    end select

  end function xyz_timeIntLFAM3

  function wt_timeIntLFAM3( valN, valB, RHS, stage) result(val)
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
       val = ( 5d0*val + 8d0*valN - valB )/12d0
    case(2)
       val = valN + dt*RHS
    end select

  end function wt_timeIntLFAM3


  function xy_timeIntLFAM3( valN, valB, RHS, stage) result(val)
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
       val = ( 5d0*val + 8d0*valN - valB )/12d0
    case(2)
       val = valN + dt*RHS
    end select

  end function xy_timeIntLFAM3


end module TemporalIntegUtil_mod
