!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module EOSDriver_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP
       

  use EqState_JM95_mod, only: &
       & EqState_JM95_Init, EqState_JM95_Final, &
       & EqState_JM95_Eval

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  interface EOSDriver_Eval
     module procedure EOSDriver_Eval_element
     module procedure EOSDriver_Eval_volfield
  end interface EOSDriver_Eval

  public :: EOSDriver_Init, EOSDriver_Final
  public :: EOSDriver_Eval

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EOSDriver_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine EOSDriver_Init()

    ! 実行文; Executable statements
    !

    call EqState_JM95_Init()


!    write(*,*) "Check EOS val:", EqState_JM95_Eval(S=35.5d0, p=300d0, theta=3d0)

  end subroutine EOSDriver_Init

  !>
  !!
  !!
  subroutine EOSDriver_Final()

    ! 実行文; Executable statements
    !

    call EqState_JM95_Final()

  end subroutine EOSDriver_Final

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_element(rho, theta, S, p)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: rho, theta, S, p
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    rho = EqState_JM95_Eval(theta, S, p)

  end subroutine EOSDriver_Eval_element

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_volfield(rho, theta, S, p)
    
    ! モジュール引用; Use statement
    !
    use GeometricField_mod, only: &
         volScalarField

    use GridSet_mod, only: &
         & GridSet_getLocalMeshInfo, nVzLyr

    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: rho, theta, S, p
    
    ! 局所変数
    ! Local variables
    !
    integer :: i, nCell
    
    ! 実行文; Executable statement
    !
    
    call GridSet_getLocalMeshInfo(rho%mesh, nCell=nCell)

    !$omp parallel do
    do i=1, nCell
       rho%data%v_(1:nVzLyr,i) = EqState_JM95_Eval( &
            & theta%data%v_(1:nVzLyr,i), S%data%v_(1:nVzLyr,i), p%data%v_(1:nVzLyr,i))
    end do

  end subroutine EOSDriver_Eval_volfield

end module EOSDriver_mod

