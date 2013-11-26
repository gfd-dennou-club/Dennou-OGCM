!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_WindDrivenCirculation_mod 

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
  public :: Exp_WindDrivenCirculation_Init, Exp_WindDrivenCirculation_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_WindDrivenCirculation_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirculation_Init()

    ! 実行文; Executable statements
    !

  end subroutine Exp_WindDrivenCirculation_Init

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirculation_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_WindDrivenCirculation_Final

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
    real(DP), parameter :: Tau0 = 0.1d0

    integer :: m
    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 5.2d03
    xy_totDepthBasic = h0
    xy_SurfHeightN = 0d0

    
    
    refDens = 1000.0d0
    refPTemp = 288.543d0

    

    xy_WindStressU = construct_WindStressU()
    xy_WindStressV = 0d0

    do k=0, kMax
       z_PTempBasic(k) = eval_PTempBasic(g_Sig(k))
    end do

!!$write(*,*) "-- WindStressU --"
!!$write(*,*) xy_WindStressU(1,:)
!!$stop
!!$write(*,*) "-- PTempBasic --"
!!$write(*,*) z_PTempBasic
!!$stop
  contains
    function construct_WindStressU() result(xy)
      real(DP) :: xy(0:iMax-1,jMax)

      real(DP), parameter :: coef(8) = &
         & (/ 0.0265682, -0.0784899, -0.00880389, 0.0343205, 0.0233334, &
         & 0.000641955, -0.00387676, -0.00150998 /)
      integer :: m

      xy = 0d0
      do m=1, size(coef)
         xy = xy + coef(m)*cos((2*m-1)*xyz_Lat(:,:,1))
      end do

    end function construct_WindStressU

    function eval_PTempBasic(sig) result(z)
      real(DP), intent(in) :: sig
      real(DP) :: z

      real(DP), parameter :: coef(13) = &
         & (/ 277.121, 4.73219, 2.93132, 1.67006, 0.945594, 0.566825, &
         &    0.382828, 0.295956, 0.197681, 0.128343, 0.0627121, 0.0400944, -0.0106413 /)
      integer :: m

      z = 0d0
      do m=1, size(coef)
         z = z + coef(m)*cos((m-1)*acos(1d0+2*sig))
      end do

    end function eval_PTempBasic

  end subroutine setInitCondition

end module Exp_WindDrivenCirculation_mod

