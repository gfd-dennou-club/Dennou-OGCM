!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module BarotModeTimeFilter_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BarotModeTimeFilter_Init, BarotModeTimeFilter_Final

  real(DP), public, save, allocatable :: Am(:), Bm(:)
  integer, public, save :: nTotBarotStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BarotModeTimeFilter_mod' !< Module Name

  integer, parameter :: p = 2
  integer, parameter :: q = 4
  real(DP), parameter :: r = 2.846158d-01
  real(DP), parameter :: A0 = 11.61885325553501857598348578903824090958
  real(DP), parameter :: Tau0=1.445999918256672067684576177271082997322
  real(DP), parameter :: TauMax = 1.316661018018372653060055199603084474802

contains

  !>
  !!
  !!
  subroutine BarotModeTimeFilter_Init(DelTimeLong, DelTimeShort)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTimeLong
    real(DP), intent(in) :: DelTimeShort
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: tau
    real(DP) :: dTau
    integer :: m
    integer :: nOriBarotStep

    ! 実行文; Executable statements
    !

    nOriBarotStep = int( DelTimeLong/DelTimeShort )
    nTotBarotStep = int( DelTimeLong*TauMax/DelTimeShort )
    allocate( Am(nTotBarotStep), Bm(nTotBarotStep) )

    dTau = 1d0/dble(nOriBarotStep)
    tau = 0d0

    do m=1, nTotBarotStep
       tau = tau + dTau
       Am(m) = ShapeFunc(tau)
    end do

    ! Normalization for Am ( sum_m A_m = 1 )
    Am = Am / sum(Am)
    
    do m=1, nTotBarotStep
       Bm(m) = sum( Am(m:nTotBarotStep) )/dble(nOriBarotStep)
    end do

#ifdef DEBUG
    call MessageNotify("M", module_name, &
         & "Total number of barotropic step=%d", i=(/ nTotBarotStep /))
#endif

  end subroutine BarotModeTimeFilter_Init

  !>
  !!
  !!
  subroutine BarotModeTimeFilter_Final()

    ! 実行文; Executable statements
    !

    deallocate( Am, Bm )

  end subroutine BarotModeTimeFilter_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !! @return 
  !!
  pure function ShapeFunc(tau) result(val)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: tau
    real(DP) :: val

    ! 局所変数
    ! Local variables
    !
    real(DP) :: x
    
    ! 実行文; Executable statement
    !
    
    x = tau/Tau0
    val = A0*( x**p - x**(p+q) - r*x )

  end function ShapeFunc

end module BarotModeTimeFilter_mod
