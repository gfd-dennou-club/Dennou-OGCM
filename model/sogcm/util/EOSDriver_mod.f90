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

  !* gtool5
  !
  use dc_types, only: DP

  !* Dennou-OGCM
  !

  use UnitConversion_mod, only: &
       & Pa2bar, K2degC

  use EOS_Linear_mod, only: &
       & EOSTYPE_LINEAR, EOS_Linear_Init, EOS_Linear_Final, &
       & EOS_Linear_Eval, EOS_Linear_PTemp2Temp

  use EOS_SimpleNonLinear_mod, only: &
       & EOSTYPE_SIMPLENONLINEAR, EOS_SimpleNonLinear_Init, EOS_SimpleNonLinear_Final, &
       & EOS_SimpleNonLinear_Eval, EOS_SimpleNonLinear_PTemp2Temp

  use EOS_JM95_mod, only: &
       & EOSTYPE_JM95, EOS_JM95_Init, EOS_JM95_Final, &
       & EOS_JM95_Eval, EOS_JM95_PTemp2Temp

  use Constants_mod, only: &
       & RefDens, RefTemp, ThermalExpanCoef, Cp0, RefSoundSpeed

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
     module procedure EOSDriver_Eval_array1d
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
    case (EOSTYPE_SIMPLENONLINEAR)
       call EOS_SimpleNonLinear_Init()
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
    case(EOSTYPE_SIMPLENONLINEAR)
       call EOS_SimpleNonLinear_Final()
    case(EOSTYPE_JM95)
       call EOS_JM95_Final()
    end select

  end subroutine EOSDriver_Final

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_element(rhoEdd, theta, S, p)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: rhoEdd
    real(DP), intent(in) :: theta, S, p
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       rhoEdd = EOS_Linear_Eval(theta, S, p)
    case(EOSTYPE_SIMPLENONLINEAR)
       rhoEdd = EOS_SimpleNonLinear_Eval(theta, S, p)
    case(EOSTYPE_JM95)
       rhoEdd = EOS_JM95_Eval(K2degC(theta), S, Pa2bar(p)) - RefDens
    end select

  end subroutine EOSDriver_Eval_element

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_array1d(rhoEdd, theta, S, p)
    
    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(:), intent(in) ::  theta, S, p
    real(DP), dimension(size(theta)), intent(out) :: rhoEdd
    
    ! 局所変数
    ! Local variables
    !

    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       rhoEdd = EOS_Linear_Eval(theta, S, p)
    case(EOSTYPE_SIMPLENONLINEAR)
       rhoEdd = EOS_SimpleNonLinear_Eval(theta, S, p)
    case(EOSTYPE_JM95)
       rhoEdd = EOS_JM95_Eval(K2degC(theta), S, Pa2bar(p)) - RefDens
    end select

  end subroutine EOSDriver_Eval_array1d

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
    integer :: k
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)

       !$omp parallel do
       do k=0, kMax
          rhoEdd(:,:,k) = EOS_Linear_Eval(theta(:,:,k), S(:,:,k), p(:,:,k))
       end do

    case(EOSTYPE_SIMPLENONLINEAR)

       !$omp parallel do
       do k=0, kMax
          rhoEdd(:,:,k) = EOS_SimpleNonLinear_Eval(theta(:,:,k), S(:,:,k), p(:,:,k))
       end do

    case(EOSTYPE_JM95)

       !$omp parallel do
       do k=0, kMax
          rhoEdd(:,:,k) = EOS_JM95_Eval(K2degC(theta(:,:,k)), S(:,:,k), Pa2bar(p(:,:,k))) - RefDens
       end do

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

    case(EOSTYPE_SIMPLENONLINEAR)

       !$omp parallel workshare
       InSituTemp = EOS_SimpleNonLinear_PTemp2Temp(p, S, theta, 0d0)
       !$omp end parallel workshare

    case(EOSTYPE_JM95)

       !$omp parallel workshare
       InSituTemp = EOS_JM95_PTemp2Temp(p, S, theta, 0d0)
       !$omp end parallel workshare

    end select
    
  end subroutine EOSDriver_PTemp2Temp_array3d

end module EOSDriver_mod

