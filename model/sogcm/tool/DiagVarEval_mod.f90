!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DiagVarEval_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  

  use Constants_mod, only: &
       & RPlanet, Grav

  use GridSet_mod, only: &
       & iMax, jMax, kMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DiagVarEval_Init, DiagVarEval_Final
  public :: eval_Vor, eval_Div, eval_totPress

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DiagVarEval_mod' !< Module Name

real(DP), parameter :: Cp = 3986d0
real(DP), parameter :: BetaT = 1.67d-04
real(DP), parameter :: T0 = 283d0
real(DP), parameter :: RefDens = 1d03
real(DP) :: H_T

contains

  !>
  !!
  !!
  subroutine DiagVarEval_Init()

    ! 実行文; Executable statements
    !

  end subroutine DiagVarEval_Init

  !>
  !!
  !!
  subroutine DiagVarEval_Final()

    ! 実行文; Executable statements
    !

  end subroutine DiagVarEval_Final

  !> @brief 
  !!
  !! @return 
  !!
  function eval_Div(xyz_u, xyz_v) result(xyz_Div)
    
    ! 宣言文; Declaration statement
    !
    real(DP) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_v(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Div(0:iMax-1,jMax,0:kMax)

    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    xyz_Div = xyz_wz( &
         & wz_AlphaOptr_xyz(xyz_u*cos(xyz_Lat), xyz_v*cos(xyz_Lat)) &
         & )

  end function eval_Div

  !> @brief 
  !!
  !! @return 
  !!
  function eval_Vor(xyz_u, xyz_v) result(xyz_Vor)
    
    ! 宣言文; Declaration statement
    !
    real(DP) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_v(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Vor(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !
    xyz_Vor = xyz_wz( &
         & wz_AlphaOptr_xyz(xyz_v*cos(xyz_Lat), -xyz_u*cos(xyz_Lat)) &
         & )
    
  end function eval_Vor

  !> @brief 
  !!
  !! @return 
  !!
  function eval_totPress(xy_surfPress, xyz_barocPress) result(xyz_totPress)
    
    ! 宣言文; Declaration statement
    !
    real(DP) :: xy_surfPress(0:iMax-1,jMax)
    real(DP) :: xyz_barocPress(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_totPress(0:iMax-1,jMax,0:kMax)

    !
    !
    !
    integer :: k

    ! 実行文; Executable statement
    !

    do k=0, kMax
       xyz_totPress(:,:,k) = xy_surfPress(:,:) + xyz_barocPress(:,:,k)
    end do
  end function eval_totPress

end module DiagVarEval_mod
