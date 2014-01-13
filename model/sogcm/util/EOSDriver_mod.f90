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

  use EqState_Linear_mod, only: &
       & EqState_Linear_Init, EqState_Linear_Final, &
       & EqState_Linear_Eval

  use EqState_JM95_mod, only: &
       & EqState_JM95_Init, EqState_JM95_Final, &
       & EqState_JM95_Eval

  use Constants_mod 

  use GridSet_mod, only: &
       & iMax, jMax, kMax

  use GovernEqSet_mod, only: &
       & EOSTYPE_LINEAR, EOSTYPE_JM95

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  interface EOSDriver_Eval
     module procedure EOSDriver_Eval_element
     module procedure EOSDriver_Eval_array3d
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
  integer :: EOSType

contains

  !>
  !!
  !!
  subroutine EOSDriver_Init(EqOfStateType)

    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: EqOfStateType

    ! 実行文; Executable statements
    !

    EOSType = EqOfStateType

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call EqState_Linear_Init( & 
            & refDens_=RefDens, refTemp_=RefTemp, BetaT_=ThermalExpanCoef, Cp0_=Cp0, Cs0_=RefSoundSpeed  )
    case(EOSTYPE_JM95)
       call EqState_JM95_Init()
    end select


!    write(*,*) "Check EOS val:", EqState_JM95_Eval(S=35.5d0, p=300d0, theta=3d0)

  end subroutine EOSDriver_Init

  !>
  !!
  !!
  subroutine EOSDriver_Final()

    ! 実行文; Executable statements
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call EqState_Linear_Final()
    case(EOSTYPE_JM95)
       call EqState_JM95_Final()
    end select

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

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       rho = EqState_Linear_Eval(theta, S, p)
    case(EOSTYPE_JM95)
       rho = EqState_JM95_Eval(theta, S, p)
    end select



  end subroutine EOSDriver_Eval_element

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_array3d(rhoEdd, theta, S, p)
    
    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: rhoEdd
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) ::  theta, S, p

    
    ! 局所変数
    ! Local variables
    !

    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       rhoEdd = EqState_Linear_Eval(theta, S, p)
    case(EOSTYPE_JM95)
       rhoEdd = EqState_JM95_Eval(theta, S, p) - RefDens
    end select

  end subroutine EOSDriver_Eval_array3d

end module EOSDriver_mod

