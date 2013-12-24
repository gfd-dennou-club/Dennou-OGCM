!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_InternalGravWave_mod 

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message

  use Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefTemp

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_InternalGravWave_Init, Exp_InternalGravWave_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_InternalGravWave_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine Exp_InternalGravWave_Init()

    ! 実行文; Executable statements
    !

  end subroutine Exp_InternalGravWave_Init

  !>
  !!
  !!
  subroutine Exp_InternalGravWave_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_InternalGravWave_Final

  !> @brief 
  !!
  !!
  subroutine setInitCondition()
    
    !
    !
    use GridSet_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & xyz_Lat, xyz_Lon

    use VariableSet_mod

    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0
    real(DP) :: h0, u0, z
    integer :: k
    real(DP), parameter :: N2 = 1.5d-04 

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 8d03
    xy_totDepthBasic = h0
    xy_SurfHeightN = 0d0

    U0 = 0.000001d0
    xyz_UN = - U0*sin(xyz_Lon) 
    xyz_VN = - U0*sin(xyz_Lat)*cos(xyz_Lon)
    
    do k=0, kMax
       z = h0*g_Sig(k)
       xyz_UN(:,:,k) = xyz_UN(:,:,k)*cos(PI*z/h0)
       xyz_VN(:,:,k) = xyz_VN(:,:,k)*cos(PI*z/h0)
       z_PTempBasic(k) = refTemp + refTemp*N2/Grav*z
    end do

  end subroutine setInitCondition

end module Exp_InternalGravWave_mod
