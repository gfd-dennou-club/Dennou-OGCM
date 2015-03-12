!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module TemporalIntegUtil_mod2

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

  interface timeIntEuler
     module procedure timeIntEuler_ary1D
     module procedure timeIntEuler_ary2D
     module procedure timeIntEuler_ary3D
  end interface timeIntEuler

  interface timeIntRK
     module procedure timeIntRK_ary1D
     module procedure timeIntRK_ary2D
     module procedure timeIntRK_ary3D
  end interface timeIntRK

  interface timeIntLF_IMEX
     module procedure timeIntLF_IMEX_ary1D
     module procedure timeIntLF_IMEX_ary2D
     module procedure timeIntLF_IMEX_ary3D
  end interface timeIntLF_IMEX

  interface timeIntLFAM3
     module procedure timeIntLFAM3_ary1D
     module procedure timeIntLFAM3_ary2D
     module procedure timeIntLFAM3_ary3D
  end interface timeIntLFAM3

  interface timeIntLFAM3_IMEX
     module procedure timeIntLFAM3_IMEX_ary1D
     module procedure timeIntLFAM3_IMEX_ary2D
     module procedure timeIntLFAM3_IMEX_ary3D
  end interface timeIntLFAM3_IMEX

  public :: timeIntEuler, timeIntRK, timeIntLFAM3, timeIntLFAM3_IMEX, timeIntLF_IMEX
  

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

  integer, parameter, public :: TimeIntMode_LF = 8
  character(*), parameter, public :: TimeIntModeLBL_LF = 'TimeIntMode_LF' 

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'TemporalIntegUtil_mod' !< Module Name

  real(DP) :: dt

  !
  real(DP), parameter :: RK2_coefA(2) = (/ 0d0, 1d0 /)
  real(DP), parameter :: RK2_coefB(2)   = (/ 0.5d0, 0.5d0 /)

  real(DP), parameter :: RK3_coefA(3) = (/ 0d0, 1d0/3d0, 2d0/3d0 /)
  real(DP), parameter :: RK3_coefB(3)   = (/ 1d0, 0d0, 3d0/4d0 /)

  real(DP), parameter :: RK4_coefA(4) = (/ 0d0, 0.5d0, 0.5d0, 1d0 /)  
  real(DP), parameter :: RK4_coefB(4) = (/ 1d0, 2d0, 2d0, 1d0 /)/6d0

  real(DP) :: RKCoefA(4,4)
  real(DP) :: RKCoefB(4,4)

  !
  real(DP), parameter :: DtCoef_Euler(1) = (/ 1d0 /)
  real(DP), parameter :: DtCoef_RK2(2) = (/ 1d0, 1d0 /)
  real(DP), parameter :: DtCoef_RK3(3) = (/ 1d0/3d0, 2d0/3d0, 1d0 /)
  real(DP), parameter :: DtCoef_RK4(4) = (/ 0.5d0, 0.5d0, 1d0, 1d0  /)
  real(DP), parameter :: DtCoef_LF(1) = (/ 2d0 /)
  real(DP), parameter :: DtCoef_LFTR(2) = (/ 2d0, 1d0 /)
  real(DP), parameter :: DtCoef_LFAM3(2) = (/ 2d0, 1d0 /)
  real(DP), parameter :: DtCoef_PC23_AB3AM3CR(2) = (/ 1d0, 1d0 /)

contains

!!!!!!!!!!!!!!
  
  !>
  !!
  !!
  subroutine TemporalIntegUtil_Init( delTime )

    ! 
    !
    real(DP), intent(in) :: DelTime

    ! 実行文; Executable statements
    !
    dt = DelTime

    RKCoefA = 0d0; RKCoefB = 0d0
    RKCoefA(1:2,2) = RK2_coefA; RKCoefB(1:2,2) = RK2_coefB
    RKCoefA(1:3,3) = RK3_coefA; RKCoefB(1:3,3) = RK3_coefB
    RKCoefA(1:4,4) = RK4_coefA; RKCoefB(1:4,4) = RK4_coefB

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
    case(TimeIntModeLBL_LF)
       tIntModeID = timeIntMode_LF
       is_VarB_Used = .true.
       nStage = 1
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
       coef = DtCoef_Euler(Stage)
    case(timeIntMode_RK2)
       coef = DtCoef_RK2(Stage)
    case(timeIntMode_RK3)
       coef = DtCoef_RK3(Stage)
    case(timeIntMode_RK4)
       coef = DtCoef_RK4(Stage)
    case(timeIntMode_LFTR)
       coef = DtCoef_LFTR(Stage)
    case(timeIntMode_LF)
       coef = DtCoef_LF(Stage)
    case(timeIntMode_LFAM3)
       coef = DtCoef_LFAM3(Stage)
    case(TimeIntMode_PC23_AB2AM3CR)
       coef = DtCoef_PC23_AB3AM3CR(Stage)
    end select

  end function TemporalIntegUtil_GetDDtCoef

!!!!!!!!!* Euler scheme *******************************

  function timeIntEuler_ary1D(N, RHSN) result(A)
    real(DP), intent(in) :: N(:), RHSN(:)
    real(DP) :: A(size(N,1))
    A(:) = N + dt*RHSN
  end function timeIntEuler_ary1D

  function timeIntEuler_ary2D(N, RHSN) result(A)
    real(DP), intent(in) :: N(:,:), RHSN(:,:)
    real(DP) :: A(size(N,1),size(N,2))
    A(:,:) = N + dt*RHSN
  end function timeIntEuler_ary2D

  function timeIntEuler_ary3D(N, RHSN) result(A)
    real(DP), intent(in) :: N(:,:,:), RHSN(:,:,:)
    real(DP) :: A(size(N,1),size(N,2),size(N,3))
    A(:,:,:) = N + dt*RHSN
  end function timeIntEuler_ary3D

!!!!!!!!!!* Runge=Kutta scheme **********************************

  function timeIntRK_ary1D(N0, RHS, RKOrder, Stage, RKTmp) result(A)
    real(DP), intent(in) :: N0(:), RHS(:)
    integer, intent(in) :: RKOrder, Stage
    real(DP), intent(inout) :: RKTmp(:)
    real(DP) :: A(size(N0,1))

    if(Stage==1) RKTmp(:) = N0
    if(Stage==RKOrder) then
       A(:) = RKTmp + dt*RKCoefB(Stage,RKOrder)*RHS; return
    end if

    RKTmp(:) = RKTmp(:) + dt*RKCoefB(Stage,RKOrder)*RHS
    A(:) = N0 + dt*RKCoefA(Stage+1,RKOrder)*RHS
  end function timeIntRK_ary1D

  function timeIntRK_ary2D(N0, RHS, RKOrder, Stage, RKTmp) result(A)
    real(DP), intent(in) :: N0(:,:), RHS(:,:)
    integer, intent(in) :: RKOrder, Stage
    real(DP), intent(inout) :: RKTmp(:,:)
    real(DP) :: A(size(N0,1),size(N0,2))

    if(Stage==1) RKTmp(:,:) = N0
    if(Stage==RKOrder) then
       A(:,:) = RKTmp + dt*RKCoefB(Stage,RKOrder)*RHS; return
    end if

    RKTmp(:,:) = RKTmp + dt*RKCoefB(Stage,RKOrder)*RHS
    A(:,:) = N0 + dt*RKCoefA(Stage+1,RKOrder)*RHS
  end function timeIntRK_ary2D

  function timeIntRK_ary3D(N0, RHS, RKOrder, Stage, RKTmp) result(A)
    real(DP), intent(in) :: N0(:,:,:), RHS(:,:,:)
    integer, intent(in) :: RKOrder, Stage
    real(DP), intent(inout) :: RKTmp(:,:,:)
    real(DP) :: A(size(N0,1),size(N0,2),size(N0,3))
    
    if(Stage==1) RKTmp(:,:,:) = N0
    if(Stage==RKOrder) then
       A(:,:,:) = RKTmp + dt*RKCoefB(Stage,RKOrder)*RHS; return
    end if

    RKTmp(:,:,:) = RKTmp + dt*RKCoefB(Stage,RKOrder)*RHS
    A(:,:,:) = N0 + dt*RKCoefA(Stage+1,RKOrder)*RHS
  end function timeIntRK_ary3D

  !*** Semi Implicit LF Scheme **************

  function timeIntLF_IMEX_ary1D(B, DImpl, Stage) result(A)
    real(DP), intent(in) :: B(:), DImpl(:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(B,1))

    A(:) = B + DImpl
  end function timeIntLF_IMEX_ary1D

  function timeIntLF_IMEX_ary2D(B, DImpl, Stage) result(A)
    real(DP), intent(in) :: B(:,:), DImpl(:,:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(B,1),size(B,2))

    A(:,:) = B + DImpl
  end function timeIntLF_IMEX_ary2D

  function timeIntLF_IMEX_ary3D(B, DImpl, Stage) result(A)
    real(DP), intent(in) :: B(:,:,:), DImpl(:,:,:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(B,1),size(B,2),size(B,3))

    A(:,:,:) = B + DImpl
  end function timeIntLF_IMEX_ary3D

  !*** LFAM3 Scheme **************

  function timeIntLFAM3_ary1D(N, B, RHS, Stage) result(A)
    real(DP), intent(in) :: N(:), B(:), RHS(:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(N,1))

    select case(Stage)
    case(1)
       A(:) = (5d0*(B + 2d0*dt*RHS) + 8d0*N - B)/12d0 
    case(2)
       A(:) = N + dt*RHS
    end select
  end function timeIntLFAM3_ary1D

  function timeIntLFAM3_ary2D(N, B, RHS, Stage) result(A)
    real(DP), intent(in) :: N(:,:), B(:,:), RHS(:,:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(N,1),size(N,2))

    select case(Stage)
    case(1)
       A(:,:) = (5d0*(B + 2d0*dt*RHS) + 8d0*N - B)/12d0 
    case(2)
       A(:,:) = N + dt*RHS
    end select
  end function timeIntLFAM3_ary2D

  function timeIntLFAM3_ary3D(N, B, RHS, Stage) result(A)
    real(DP), intent(in) :: N(:,:,:), B(:,:,:), RHS(:,:,:)
    integer, intent(in) :: Stage
    real(DP) :: A(size(N,1),size(N,2),size(N,3))

    select case(Stage)
    case(1)
       A(:,:,:) = (5d0*(B + 2d0*dt*RHS) + 8d0*N - B)/12d0 
    case(2)
       A(:,:,:) = N + dt*RHS
    end select
  end function timeIntLFAM3_ary3D

  !*** Semi Implicit LFAM3 Scheme **************
  ! When AM3 is .false., DTmp means the increment u(t+dt) - u(t-N*dt),
  ! where N equals 1 and 0 in the first stage and second stage respectively.
  ! On the other hand, when AM3 is .true., DTmp means the variable at time
  ! level t+dt obtained in the first stage.
  !
  function timeIntLFAM3_IMEX_ary1D(N, B, DTmp, Stage, AM3) result(ret)
    real(DP), intent(in) :: N(:), B(:), DTmp(:)
    integer, intent(in) :: Stage
    logical, intent(in) :: AM3
    real(DP) :: ret(size(N,1))

    select case(Stage)
    case(1)
       if(.not. AM3) then
          ret(:) = B + DTmp
       else
          ret(:) = (5d0*DTmp + 8d0*N - B)/12d0
       end if
    case(2)
       ret(:) = N + DTmp
    end select
  end function timeIntLFAM3_IMEX_ary1D

  function timeIntLFAM3_IMEX_ary2D(N, B, DTmp, Stage, AM3) result(ret)
    real(DP), intent(in) :: N(:,:), B(:,:), DTmp(:,:)
    integer, intent(in) :: Stage
    logical, intent(in) :: AM3
    real(DP) :: ret(size(N,1),size(N,2))

    select case(Stage)
    case(1)
       if(.not. AM3) then
          ret(:,:) = B + DTmp
       else
          ret(:,:) = (5d0*DTmp + 8d0*N - B)/12d0
       end if
    case(2)
       ret(:,:) = N + DTmp
    end select
  end function timeIntLFAM3_IMEX_ary2D

  function timeIntLFAM3_IMEX_ary3D(N, B, DTmp, Stage, AM3) result(ret)
    real(DP), intent(in) :: N(:,:,:), B(:,:,:), DTmp(:,:,:)
    integer, intent(in) :: Stage
    logical, intent(in) :: AM3
    real(DP) :: ret(size(N,1),size(N,2),size(N,3))

    select case(Stage)
    case(1)
       if(.not. AM3) then
          ret(:,:,:) = B + DTmp
       else
          ret(:,:,:) = (5d0*DTmp + 8d0*N - B)/12d0
       end if
    case(2)
       ret(:,:,:) = N + DTmp
    end select
  end function timeIntLFAM3_IMEX_ary3D
  
end module TemporalIntegUtil_mod2

