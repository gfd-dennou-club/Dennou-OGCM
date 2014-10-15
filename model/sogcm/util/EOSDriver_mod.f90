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

  use EOS_Linear_mod, only: &
       & EOSTYPE_LINEAR, EOS_Linear_Init, EOS_Linear_Final, &
       & EOS_Linear_Eval, EOS_Linear_PTemp2Temp

  use EOS_SimpleNonLinear_mod, only: &
       & EOSTYPE_SIMPLENONLINEAR, EOS_SimpleNonLinear_Init, EOS_SimpleNonLinear_Final, &
       & EOS_SimpleNonLinear_Eval, EOS_SimpleNonLinear_PTemp2Temp

  use EOS_JM95_mod, only: &
       & EOSTYPE_JM95, EOS_JM95_Init, EOS_JM95_Final, &
       & EOS_JM95_Eval, EOS_JM95_PTemp2Temp

  use Constants_mod 

  use GridSet_mod, only: &
       & iMax, jMax, kMax


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

  interface EOSDriver_PTemp2Temp
     module procedure EOSDriver_PTemp2Temp_array3d
  end interface EOSDriver_PTemp2Temp

  public :: EOSDriver_Init, EOSDriver_Final
  public :: EOSDriver_Eval
  public :: EOSDriver_PTemp2Temp

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
       call EOS_Linear_Init( & 
            & refDens_=RefDens, refTemp_=RefTemp, BetaT_=ThermalExpanCoef, Cp0_=Cp0, Cs0_=RefSoundSpeed  )
    case(EOSTYPE_JM95)
       call EOS_JM95_Init()
    end select


!    write(*,*) "Check EOS val:", EOS_JM95_Eval(S=35.5d0, p=300d0, theta=3d0)

  end subroutine EOSDriver_Init

  !>
  !!
  !!
  subroutine EOSDriver_Final()

    ! 実行文; Executable statements
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call EOS_Linear_Final()
    case(EOSTYPE_JM95)
       call EOS_JM95_Final()
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
       rho = EOS_Linear_Eval(theta, S, p)
    case(EOSTYPE_JM95)
       rho = EOS_JM95_Eval(theta, S, p)
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

       !$omp parallel workshare
       rhoEdd = EOS_Linear_Eval(theta, S, p)
       !$omp end parallel workshare

    case(EOSTYPE_JM95)

       !$omp parallel workshare
       rhoEdd = EOS_JM95_Eval(theta, S, p) - RefDens
       !$omp end parallel workshare

    end select

  end subroutine EOSDriver_Eval_array3d

  !> @brief 
  !!
  !!
  subroutine EOSDriver_PTemp2Temp_array3d(InSituTemp, theta, S, p)
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: InSituTemp
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) ::  theta, S, p
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)

       !$omp parallel workshare
       InSituTemp = EOS_Linear_PTemp2Temp(p, S, theta, 0d0)
       !$omp end parallel workshare

    case(EOSTYPE_JM95)

       !$omp parallel workshare
       InSituTemp = EOS_JM95_PTemp2Temp(p, S, theta, 0d0)
       !$omp end parallel workshare

    end select
    
  end subroutine EOSDriver_PTemp2Temp_array3d

end module EOSDriver_mod

