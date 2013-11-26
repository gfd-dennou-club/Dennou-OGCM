!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_W94_Case2_mod 

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message

  use Constants_mod, only: &
       & Grav, PI, RPlanet, Omega

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_W94_Case2_Init, Exp_W94_Case2_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_W94_Case2_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine Exp_W94_Case2_Init()

    ! 実行文; Executable statements
    !

  end subroutine Exp_W94_Case2_Init

  !> @brief 
  !!
  !!
  subroutine setInitCondition()
    
    !
    !
    use GridSet_mod, only: xyz_Lat

    use VariableSet_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0
    real(DP) :: h0, u0
    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 2.94d04 / Grav
    u0 = 2d0*PI*RPlanet / (3600d0*24d0*12d0)

    xy_totDepthBasic = h0
    xyz_UN = u0*cos(xyz_Lat)
    xy_SurfHeightN = - (RPlanet*Omega*u0 + 0.5d0*u0**2) * sin(xyz_Lat(:,:,1))**2 / Grav

  end subroutine setInitCondition

  !>
  !!
  !!
  subroutine Exp_W94_Case2_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_W94_Case2_Final

end module Exp_W94_Case2_mod

