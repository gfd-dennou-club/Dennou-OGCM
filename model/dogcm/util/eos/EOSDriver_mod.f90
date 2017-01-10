!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief A driver module to access some utility modules of sea-water EOS. 
!! 
!! @author Kawai Yuta
!!
!!
module EOSDriver_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify
  
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
       & EOS_JM95_Eval, EOS_JM95_PTemp2Temp,          &
       & EOS_JM95_alpha_beta
  

  use DOGCM_Admin_Constants_mod, only: &
       & RefDens, RefTemp, ThermalExpanCoef, Cp0, RefSoundSpeed



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
     module procedure EOSDriver_Eval_array2d
     module procedure EOSDriver_Eval_array3d
  end interface EOSDriver_Eval

  interface EOSDriver_PTemp2Temp
     module procedure EOSDriver_PTemp2Temp_array3d
  end interface EOSDriver_PTemp2Temp

  interface EOSDriver_alpha_beta
     module procedure EOSDriver_Eval_alpha_beta_element
     module procedure EOSDriver_Eval_alpha_beta_array2d     
     module procedure EOSDriver_Eval_alpha_beta_array3d
  end interface EOSDriver_alpha_beta
  
  public :: EOSDriver_Init, EOSDriver_Final
  public :: EOSDriver_Eval
  public :: EOSDriver_alpha_beta
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

  !----------------------------------------------
  
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

  !----------------------------------------------
  
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
    real(DP), intent(in) ::  theta(:)
    real(DP), intent(in) ::  S(:)
    real(DP), intent(in) ::  p(:)
    real(DP), intent(out) :: rhoEdd(size(theta))
    
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

  !----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_array2d(rhoEdd, theta, S, p)
    
    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: rhoEdd(:,:)
    real(DP), intent(in) ::  theta(:,:)
    real(DP), intent(in) ::  S(:,:)
    real(DP), intent(in) ::  p(:,:)

    
    ! 局所変数
    ! Local variables
    !
    integer :: k
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       rhoEdd(:,:) = EOS_Linear_Eval(theta(:,:), S(:,:), p(:,:))
    case(EOSTYPE_SIMPLENONLINEAR)
       rhoEdd(:,:) = EOS_SimpleNonLinear_Eval(theta(:,:), S(:,:), p(:,:))
    case(EOSTYPE_JM95)
       rhoEdd(:,:) = EOS_JM95_Eval(K2degC(theta(:,:)), S(:,:), Pa2bar(p(:,:))) - RefDens
    end select

  end subroutine EOSDriver_Eval_array2d

  !----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_array3d(rhoEdd, theta, S, p)
    
    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: rhoEdd(:,:,:)
    real(DP), intent(in) ::  theta(:,:,:)
    real(DP), intent(in) ::  S(:,:,:)
    real(DP), intent(in) ::  p(:,:,:)

    
    ! 局所変数
    ! Local variables
    !
    integer :: k
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)

       !$omp parallel do
       do k=1, size(theta,3)
          rhoEdd(:,:,k) = EOS_Linear_Eval(theta(:,:,k), S(:,:,k), p(:,:,k))
       end do

    case(EOSTYPE_SIMPLENONLINEAR)

       !$omp parallel do
       do k=1, size(theta,3)
          rhoEdd(:,:,k) = EOS_SimpleNonLinear_Eval(theta(:,:,k), S(:,:,k), p(:,:,k))
       end do

    case(EOSTYPE_JM95)

       !$omp parallel do
       do k=1, size(theta,3)
          rhoEdd(:,:,k) = EOS_JM95_Eval(K2degC(theta(:,:,k)), S(:,:,k), Pa2bar(p(:,:,k))) - RefDens
       end do

    end select

  end subroutine EOSDriver_Eval_array3d

  !----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine EOSDriver_PTemp2Temp_array3d(InSituTemp, theta, S, p)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: InSituTemp(:,:,:)
    real(DP), intent(in) ::  theta(:,:,:)
    real(DP), intent(in) ::  S(:,:,:)
    real(DP), intent(in) ::  p(:,:,:)
    
    
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

  !----------------------------------------------

  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_alpha_beta_element( alpha, beta, & ! (out)
       & theta, S, p )                                 ! (in)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: alpha
    real(DP), intent(out) :: beta
    real(DP), intent(in) :: theta
    real(DP), intent(in) :: s
    real(DP), intent(in) :: p
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call MessageNotify( 'E', module_name, "EOS_Linear_alpha_beta has not been implemented.")
    case(EOSTYPE_SIMPLENONLINEAR)
       call MessageNotify( 'E', module_name, "EOS_SimpleNonLinear_alpha_beta has not been implemented.")
    case(EOSTYPE_JM95)
       call EOS_JM95_alpha_beta( alpha, beta, &
            &  K2degC(theta), S, Pa2bar(p)*10d0 )
    end select

  end subroutine EOSDriver_Eval_alpha_beta_element

  !----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_alpha_beta_array2d( alpha, beta, & ! (out)
       & theta, S, p )                                         ! (in)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: alpha(:,:)
    real(DP), intent(out) :: beta(:,:)
    real(DP), intent(in) :: theta(:,:)
    real(DP), intent(in) :: s(:,:)
    real(DP), intent(in) :: p(:,:)
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call MessageNotify( 'E', module_name, "EOS_Linear_alpha_beta has not been implemented.")
    case(EOSTYPE_SIMPLENONLINEAR)
       call MessageNotify( 'E', module_name, "EOS_SimpleNonLinear_alpha_beta has not been implemented.")
    case(EOSTYPE_JM95)
       call EOS_JM95_alpha_beta( alpha, beta, &
            &  K2degC(theta), S, Pa2bar(p)*10d0 )
    end select

  end subroutine EOSDriver_Eval_alpha_beta_array2d
  
  !----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine EOSDriver_Eval_alpha_beta_array3d( alpha, beta, & ! (out)
       & theta, S, p )                                 ! (in)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: alpha(:,:,:)
    real(DP), intent(out) :: beta(:,:,:)
    real(DP), intent(in) :: theta(:,:,:)
    real(DP), intent(in) :: s(:,:,:)
    real(DP), intent(in) :: p(:,:,:)

    ! 作業変数
    ! Work variables
    !
    integer :: k
    
    ! 実行文; Executable statement
    !

    select case(EOSType)
    case(EOSTYPE_LINEAR)
       call MessageNotify( 'E', module_name, "EOS_Linear_alpha_beta has not been implemented.")
    case(EOSTYPE_SIMPLENONLINEAR)
       call MessageNotify( 'E', module_name, "EOS_SimpleNonLinear_alpha_beta has not been implemented.")
    case(EOSTYPE_JM95)
       !$omp parallel do
       do k = 1, size(theta,3)
          call EOS_JM95_alpha_beta( alpha(:,:,k), beta(:,:,k),             & ! (out)
               &  K2degC(theta(:,:,k)), S(:,:,k), Pa2bar(p(:,:,k))*10d0   )  ! (in)
       end do
    end select
    
  end subroutine EOSDriver_Eval_alpha_beta_array3d
  
  
end module EOSDriver_mod

