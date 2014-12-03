!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!

#define POLYORDER_Nk3

module LagrangePolyn_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  interface LagrangePolyn_Init
#ifdef POLYORDER_Nk3
     module procedure LagrangePolyn_Init_TriNk3
#elif defined  POLYORDER_Nk2
     module procedure LagrangePolyn_Init_TriNk2
#else
     module procedure LagrangePolyn_Init_TriNk1
#endif
  end interface LagrangePolyn_Init

  interface Nk1_interpolate
     module procedure Nk1_interpolate_1
     module procedure Nk1_interpolate_2
  end interface Nk1_interpolate

  interface Nk2_interpolate
     module procedure Nk2_interpolate_1
     module procedure Nk2_interpolate_2
  end interface Nk2_interpolate

  interface Nk3_interpolate
     module procedure Nk3_interpolate_1
     module procedure Nk3_interpolate_2
  end interface Nk3_interpolate

  interface TriNk_basis
#ifdef POLYORDER_Nk3
     module procedure TriNk3_basis
#elif defined  POLYORDER_Nk2
     module procedure TriNk2_basis
#else
     module procedure TriNk1_basis
#endif
  end interface TriNk_basis

  interface TriNk_basis_dy1
#ifdef POLYORDER_Nk3
     module procedure TriNk3_basis_dy1
#elif defined  POLYORDER_Nk2
     module procedure TriNk2_basis_dy1
#else
     module procedure TriNk1_basis_dy1
#endif
  end interface TriNk_basis_dy1

  interface TriNk_basis_dy2
#ifdef POLYORDER_Nk3
     module procedure TriNk3_basis_dy2
#elif defined POLYORDER_Nk2
     module procedure TriNk2_basis_dy2
#else
     module procedure TriNk1_basis_dy2
#endif
  end interface TriNk_basis_dy2
  
  interface TriNk_interpolate
#ifdef POLYORDER_Nk3
     module procedure TriNk3_interpolate_1
     module procedure TriNk3_interpolate_2
#elif defined POLYORDER_Nk2
     module procedure TriNk2_interpolate_1
     module procedure TriNk2_interpolate_2
#else
     module procedure TriNk1_interpolate_1
     module procedure TriNk1_interpolate_2
#endif
  end interface TriNk_interpolate

  interface TriNk_sinteg
#ifdef POLYORDER_Nk3
     module procedure TriNk3_sinteg_1
     module procedure TriNk3_sinteg_2
     module procedure TriNk3_sinteg_3
#elif defined POLYORDER_Nk2
     module procedure TriNk2_sinteg_1
     module procedure TriNk2_sinteg_2
     module procedure TriNk2_sinteg_3
#else
     module procedure TriNk1_sinteg_1
     module procedure TriNk1_sinteg_2
     module procedure TriNk1_sinteg_3
#endif
  end interface TriNk_sinteg
  

  interface TriNk_sinteg_dotProdWt
#ifdef POLYORDER_Nk3
     module procedure TriNk3_sinteg_dotProdWt
#elif defined POLYORDER_Nk2
     module procedure TriNk2_sinteg_dotProdWt
#else
     module procedure TriNk1_sinteg_dotProdWt
#endif
  end interface TriNk_sinteg_dotProdWt



  interface TriNk_linteg
#ifdef POLYORDER_Nk3
     module procedure TriNk3_linteg_1
     module procedure TriNk3_linteg_2
     module procedure TriNk3_linteg_3
#elif defined POLYORDER_Nk2
     module procedure TriNk2_linteg_1
     module procedure TriNk2_linteg_2
     module procedure TriNk2_linteg_3
#else
     module procedure TriNk1_linteg_1
     module procedure TriNk1_linteg_2
     module procedure TriNk1_linteg_3
#endif
  end interface TriNk_linteg

  public :: LagrangePolyn_Init
  public :: TriNk_basis, TriNk_basis_dy1, TriNk_basis_dy2
  public :: TriNk_interpolate, TriNk_sinteg, TriNk_sinteg_dotProdWt, TriNk_linteg

  ! 非公開手続き
  ! Private procedure
  !

  ! 公開変数
  ! Public variable
  !


  !
  real(DP), public, parameter :: GLQd_order3_y(2) = (/ -1d0/sqrt(3d0), 1d0/sqrt(3d0) /)
  real(DP), public, parameter :: GLQd_order3_Wt(2) = (/ 1d0, 1d0 /)

  real(DP), public, parameter :: TriNk1_y1(3) = (/ 0d0, 1d0, 0d0 /)
  real(DP), public, parameter :: TriNk1_y2(3) = (/ 0d0, 0d0, 1d0 /)

  real(DP), public, parameter :: TriNk1_SInt_y1(3) = (/ 1d0/6d0, 2d0/3d0, 1d0/6d0 /)
  real(DP), public, parameter :: TriNk1_SInt_y2(3) = (/ 1d0/6d0, 1d0/6d0, 2d0/3d0 /)
  real(DP), public, parameter :: TriNk1_SInt_Wt(3) = (/ 1d0/3d0, 1d0/3d0, 1d0/3d0 /)


  !
  real(DP), public, parameter :: GLQd_order5_y(3) = (/ -sqrt(3d0/5d0), 0d0, sqrt(3d0/5d0) /)
  real(DP), public, parameter :: GLQd_order5_Wt(3) = (/ 5d0/9d0, 8d0/9d0, 5d0/9d0 /)

  real(DP), public, parameter :: TriNk2_y1(6) = (/ 0d0, 0.5d0, 1d0, 0.5d0, 0d0,   0d0 /)
  real(DP), public, parameter :: TriNk2_y2(6) = (/ 0d0,   0d0, 0d0, 0.5d0, 1d0, 0.5d0 /)

  real(DP), public, parameter :: TriNk2_SInt_y1(6) = &
         & (/ 0.44594849091597d0, 0.44594849091597d0, 0.10810301816807d0, &
         &    0.09157621350977d0, 0.09157621350977d0, 0.81684757298046d0 /)
  real(DP), public, parameter :: TriNk2_SInt_y2(6) = &
         & (/ 0.44594849091597d0, 0.10810301816807d0, 0.44594849091597d0, &
         &    0.09157621350977d0, 0.81684757298046d0, 0.09157621350977d0 /)
  real(DP), public, parameter :: TriNk2_SInt_Wt(6) = &
         & (/ 0.22338158967801d0, 0.22338158967801d0, 0.22338158967801d0, &
         &    0.10995174365532d0, 0.10995174365532d0, 0.10995174365532d0 /)



  !
  real(DP), public, parameter :: GLQd_order7_y(4) = &
       & (/ -sqrt((3d0+2d0*sqrt(6d0/5d0))/7d0), -sqrt((3d0-2d0*sqrt(6d0/5d0))/7d0), &
       &     sqrt((3d0-2d0*sqrt(6d0/5d0))/7d0),  sqrt((3d0+2d0*sqrt(6d0/5d0))/7d0) /)
  real(DP), public, parameter :: GLQd_order7_Wt(4) = &
       & (/ (18d0 - sqrt(30d0))/36d0, (18d0 + sqrt(30d0))/36d0, &
       &    (18d0 + sqrt(30d0))/36d0, (18d0 - sqrt(30d0))/36d0 /)

  real(DP), parameter :: Phi = 0.5d0*(-1d0+sqrt(5d0))   ! golden ratio conjugate
  real(DP), public, parameter :: TriNk3_y1(10) = (/ 0d0,      Phi/sqrt(5d0), (Phi+1d0)/sqrt(5d0),        1d0,  &
       &                                                (Phi+1d0)/sqrt(5d0),       Phi/sqrt(5d0),        0d0,  &
       &                                                                                0d0,               0d0,  &
       &                                                                                               1d0/3d0 /)
  real(DP), public, parameter :: TriNk3_y2(10) = (/ 0d0,            0d0,                   0d0,         0d0,  &
       &                                                  Phi/sqrt(5d0),  (Phi+1d0)/sqrt(5d0),          1d0,  &
       &                                                                  (Phi+1d0)/sqrt(5d0), Phi/sqrt(5d0),  &
       &                                                                                             1d0/3d0 /)

  real(DP) :: TriNk3_LagrangeBasisCoef(10,10)

  real(DP), public, parameter :: TriNk3_SInt_y1(12) = &
         & (/ 0.24928674517091d0, 0.24928674517091d0, 0.50142650965818d0, &
         &    0.06308901449150d0, 0.06308901449150d0, 0.87382197101700d0, &
         &    0.31035245103378d0, 0.63650249912140d0, 0.05314504984482d0, &
         &    0.63650249912140d0, 0.31035245103378d0, 0.05314504984482d0 /)
  real(DP), public, parameter :: TriNk3_SInt_y2(12) = &
         & (/ 0.24928674517091d0, 0.50142650965818d0, 0.24928674517091d0, &
         &    0.06308901449150d0, 0.87382197101700d0, 0.06308901449150d0, &
         &    0.63650249912140d0, 0.05314504984482d0, 0.31035245103378d0, &
         &    0.31035245103378d0, 0.05314504984482d0, 0.63650249912140d0  /)
  real(DP), public, parameter :: TriNk3_SInt_Wt(12) = &
         & (/ 0.116786275726378d0, 0.116786275726378d0, 0.116786275726378d0, &
         &    0.050844906370206d0, 0.050844906370206d0, 0.050844906370206d0, &
         &    0.082851075618374d0, 0.082851075618374d0, 0.082851075618374d0,    &
         &    0.082851075618374d0, 0.082851075618374d0, 0.082851075618374d0  /)

  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'LagrangePolyn_mod' !< Module Name

contains

  !>
  !!
  !!

  !>
  !!
  !!
  subroutine LagrangePolyn_Final()

    ! 実行文; Executable statements
    !

  end subroutine LagrangePolyn_Final


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LagrangePolyn_Init_TriNk1()

    ! 実行文; Executable statements
    !

  end subroutine LagrangePolyn_Init_TriNk1

  function Nk1_interpolate_1(y, f) result(val)
    real(DP), intent(in) :: y, f(2)
    real(DP) :: val

    val = (1d0 - y)*f(1) + y*f(2)

  end function Nk1_interpolate_1

  function Nk1_interpolate_2(y, f) result(val)
    real(DP), intent(in) :: y(:), f(2)
    real(DP) :: val(size(y))

    val(:) = (1d0 - y(:))*f(1) + y(:)*f(2)

  end function Nk1_interpolate_2

  function TriNk1_interpolate_1(y1, y2, f) result(val)
    real(DP), intent(in) :: y1, y2, f(3)
    real(DP) :: val

    val = (1d0 - y1 - y2)*f(1) + y1*f(2) + y2*f(3)

  end function TriNk1_interpolate_1

  function TriNk1_interpolate_2(y1, y2, f) result(val)
    real(DP), intent(in) :: y1(:), y2(:), f(3)
    real(DP) :: val(size(y1))

    val = (1d0 - y1 - y2)*f(1) + y1*f(2) + y2*f(3)

  end function TriNk1_interpolate_2

  function TriNk1_basis(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(3)
    
    val = (/ 1d0 - y1 - y2, y1, y2 /)

  end function TriNk1_basis

  function TriNk1_basis_dy1(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(3)
    
    val = (/ -1d0, 1d0, 0d0 /)

  end function TriNk1_basis_dy1

  function TriNk1_basis_dy2(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(3)
    
    val = (/ -1d0, 0d0, 1d0 /)

  end function TriNk1_basis_dy2

  function TriNk1_sinteg_1(f) result(val)
    real(DP), intent(in) :: f(3)
    real(DP) :: val
    
    integer :: i

    val = 0d0
    do i=1, 3
       val = val + TriNk1_SInt_Wt(i)*&
            & TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), f)
    end do
    val = 0.5d0*val
         
  end function TriNk1_sinteg_1

  function TriNk1_sinteg_2(f, g) result(val)
    real(DP), intent(in) :: f(3), g(3)
    real(DP) :: val

    integer :: i

    val = 0d0
    do i=1, 3
       val = val + TriNk1_SInt_Wt(i) &
            & *TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), f) &
            & *TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), g)
    end do
    val = 0.5d0*val

  end function TriNk1_sinteg_2

  function TriNk1_sinteg_3(f, g, h) result(val)
    real(DP), intent(in) :: f(3), g(3), h(3)
    real(DP) :: val

    integer :: i

    val = 0d0
    do i=1, 3
       val = val + TriNk1_SInt_Wt(i) &
            & *TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), f) &
            & *TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), g) &
            & *TriNk1_interpolate_1(TriNk1_SInt_y1(i), TriNk1_SInt_y2(i), h)
    end do
    val = 0.5d0*val

  end function TriNk1_sinteg_3

  function TriNk1_sinteg_dotProdWt(f) result(val)
    real(DP) :: f(3)
    real(DP) :: val
    
    val = 0.5d0*sum( TriNk1_SInt_Wt(:)*f(:) )
         
  end function TriNk1_sinteg_dotProdWt

  function TriNk1_linteg_1(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(2),f2(2),f3(2)
    real(DP) :: val, val_
    
    real(DP), parameter, dimension(2) :: gaussNodes = &
         & 0.5d0*(1d0 + GLQd_order3_y(:))

!!$    val = 0.5d0*sum( GLQd_order3_Wt(:)*( &
!!$          &   Nk1_interpolate(gaussNodes, f1)           &
!!$          & + sqrt(2d0)*Nk1_interpolate(gaussNodes, f2) &
!!$          & + Nk1_interpolate(gaussNodes, f3)           &
!!$          & ))

    val = 0.5d0*(sum(f1+f3) + sqrt(2d0)*sum(f2))

  end function TriNk1_linteg_1

  function TriNk1_linteg_2(f1,g1, f2,g2, f3,g3) result(val)
    real(DP), dimension(2), intent(in) :: f1,g1, f2,g2, f3,g3
    real(DP) :: val
    
    real(DP), parameter, dimension(2) :: GLNodes = &
         & 0.5d0*(1d0 + GLQd_order3_y(:))

    val = 0.5d0*sum( GLQd_order3_Wt(:)*( &
          &              Nk1_interpolate(GLNodes, f1)*Nk1_interpolate(GLNodes, g1)  &
          & +  sqrt(2d0)*Nk1_interpolate(GLNodes, f2)*Nk1_interpolate(GLNodes, g2)  &
          & +            Nk1_interpolate(GLNodes, f3)*Nk1_interpolate(GLNodes, g3)  &
          & ))

  end function TriNk1_linteg_2

  function TriNk1_linteg_3(f1,g1,h1, f2,g2,h2, f3,g3,h3) result(val)
    real(DP), dimension(2), intent(in) :: f1,g1,h1, f2,g2,h2, f3,g3,h3
    real(DP) :: val
    
    real(DP), parameter, dimension(2) :: GLNodes = &
         & 0.5d0*(1d0 + GLQd_order3_y(:))

    val = 0.5d0*sum( GLQd_order3_Wt(:)*( &
          &              Nk1_interpolate(GLNodes, f1)*Nk1_interpolate(GLNodes, g1)*Nk1_interpolate(GLNodes, h1)  &
          & +  sqrt(2d0)*Nk1_interpolate(GLNodes, f2)*Nk1_interpolate(GLNodes, g2)*Nk1_interpolate(GLNodes, h2)  &
          & +            Nk1_interpolate(GLNodes, f3)*Nk1_interpolate(GLNodes, g3)*Nk1_interpolate(GLNodes, h3)  &
          & ))

  end function TriNk1_linteg_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LagrangePolyn_Init_TriNk2()

    ! 実行文; Executable statements
    !

  end subroutine LagrangePolyn_Init_TriNk2

  function Nk2_interpolate_1(y, f) result(val)
    real(DP), intent(in) :: y, f(3)
    real(DP) :: val

    real(DP) :: l1, l2

    l1 = 1d0 - y; l2 = y;
    val = l1*(2d0*l1 - 1d0)*f(1) + 4d0*l1*l2*f(2) + l2*(2d0*l2 - 1d0)*f(3)

  end function Nk2_interpolate_1

  function Nk2_interpolate_2(y, f) result(val)
    real(DP), intent(in) :: y(:), f(3)
    real(DP) :: val(size(y))

    real(DP), dimension(size(y)) :: l1, l2

    l1(:) = 1d0 - y(:); l2(:) = y(:);
    val(:) = l1(:)*(2d0*l1(:) - 1d0)*f(1) + 4d0*l1(:)*l2(:)*f(2) + l2(:)*(2d0*l2(:) - 1d0)*f(3)

  end function Nk2_interpolate_2


  function TriNk2_interpolate_1(y1, y2, f) result(val)
    real(DP), intent(in) :: y1, y2, f(6)
    real(DP) :: val

    real(DP) :: l1, l2, l3

    l1 = 1d0 - y1 -y2
    l2 = y1; l3 = y2;

    val =    l1*(2d0*l1-1d0)*f(1) + 4d0*l1*l2*f(2) &
         & + l2*(2d0*l2-1d0)*f(3) + 4d0*l2*l3*f(4) &
         & + l3*(2d0*l3-1d0)*f(5) + 4d0*l3*l1*f(6)

  end function TriNk2_interpolate_1

  function TriNk2_interpolate_2(y1, y2, f) result(val)
    real(DP), intent(in) :: y1(:), y2(:), f(6)
    real(DP) :: val(size(y1))

    real(DP), dimension(size(y1)) :: l1, l2, l3

    l1 = 1d0 - y1 -y2
    l2 = y1; l3 = y2;

    val =    l1*(2d0*l1-1d0)*f(1) + 4d0*l1*l2*f(2) &
         & + l2*(2d0*l2-1d0)*f(3) + 4d0*l2*l3*f(4) &
         & + l3*(2d0*l3-1d0)*f(5) + 4d0*l3*l1*f(6)

  end function TriNk2_interpolate_2
  
  function TriNk2_basis(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(6)

    real(DP) :: l1, l2, l3

    l1 = 1d0 - y1 -y2
    l2 = y1; l3 = y2;

    val = (/ l1*(2d0*l1-1d0), 4d0*l1*l2, &
         &   l2*(2d0*l2-1d0), 4d0*l2*l3, &
         &   l3*(2d0*l3-1d0), 4d0*l3*l1 /)

  end function TriNk2_basis

  function TriNk2_basis_dy1(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(6)

    real(DP) :: l1, l2, l3

    l1 = 1d0 - y1 -y2
    l2 = y1; l3 = y2;

    val = (/ -4d0*l1+1d0,    4d0*(l1-l2), &
         &    4d0*l2-1d0,         4d0*l3, &
         &           0d0,        -4d0*l3 /)

  end function TriNk2_basis_dy1

  function TriNk2_basis_dy2(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(6)

    real(DP) :: l1, l2, l3

    l1 = 1d0 - y1 -y2
    l2 = y1; l3 = y2;

    val = (/ -4d0*l1+1d0,        -4d0*l2, &
         &           0d0,         4d0*l2, &
         &    4d0*l3-1d0,    4d0*(l1-l3) /)
    
  end function TriNk2_basis_dy2

  function TriNk2_sinteg_1(f) result(val)
    real(DP), intent(in) :: f(6)
    real(DP) :: val


    integer :: i

    val = 0d0
    do i=1, 6
       val = val + TriNk2_SInt_Wt(i) &
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), f)
    end do
    val = 0.5d0*val

  end function TriNk2_sinteg_1

  function TriNk2_sinteg_2(f, g) result(val)
    real(DP), intent(in) :: f(6), g(6)
    real(DP) :: val

    integer :: i

    val = 0d0
    do i=1, 6
       val = val + TriNk2_SInt_Wt(i)&
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), f) &
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), g)
    end do
    val = 0.5d0*val

  end function TriNk2_sinteg_2

  function TriNk2_sinteg_3(f, g, h) result(val)
    real(DP), intent(in) :: f(6), g(6), h(6)
    real(DP) :: val

    integer :: i

    val = 0d0    
    do i=1, 6
       val = val + TriNk2_SInt_Wt(i)&
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), f) &
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), g) &
            & *TriNk2_interpolate_1(TriNk2_SInt_y1(i), TriNk2_SInt_y2(i), h)
    end do
    val = 0.5d0*val

  end function TriNk2_sinteg_3

  function TriNk2_sinteg_dotProdWt(f) result(val)
    real(DP) :: f(6)
    real(DP) :: val
    
    val = 0.5d0*sum( TriNk2_SInt_Wt(:)*f(:) )
         
  end function TriNk2_sinteg_dotProdWt

  function TriNk2_linteg_1(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(3),f2(3),f3(3)
    real(DP) :: val
    
    real(DP), parameter, dimension(3) :: GLNodes = &
         & 0.5d0*(1d0 + GLQd_order5_y(:))

    val = 0.5d0*sum( GLQd_order5_Wt(:)*( &
          &              Nk2_interpolate(GLNodes, f1) &
          & +  sqrt(2d0)*Nk2_interpolate(GLNodes, f2) &
          & +            Nk2_interpolate(GLNodes, f3) &
          & ))

  end function TriNk2_linteg_1

  function TriNk2_linteg_2(f1,g1, f2,g2, f3,g3) result(val)
    real(DP), dimension(3), intent(in) :: f1,g1, f2,g2, f3,g3
    real(DP) :: val
    
    real(DP), parameter, dimension(3) :: GLNodes = &
         & 0.5d0*(1d0 + GLQd_order5_y(:))

    val = 0.5d0*sum( GLQd_order5_Wt(:)*( &
          &              Nk2_interpolate(GLNodes, f1)*Nk2_interpolate(GLNodes, g1)  &
          & +  sqrt(2d0)*Nk2_interpolate(GLNodes, f2)*Nk2_interpolate(GLNodes, g2)  &
          & +            Nk2_interpolate(GLNodes, f3)*Nk2_interpolate(GLNodes, g3)  &
          & ))

  end function TriNk2_linteg_2

  function TriNk2_linteg_3(f1,g1,h1, f2,g2,h2, f3,g3,h3) result(val)
    real(DP), dimension(3), intent(in) :: f1,g1,h1, f2,g2,h2, f3,g3,h3
    real(DP) :: val
    
    real(DP), parameter, dimension(3) :: GLNodes = &
         & 0.5d0*(1d0 + GLQd_order5_y(:))

    val = 0.5d0*sum( GLQd_order5_Wt(:)*( &
          &              Nk2_interpolate(GLNodes, f1)*Nk2_interpolate(GLNodes, g1)*Nk2_interpolate(GLNodes, h1)  &
          & +  sqrt(2d0)*Nk2_interpolate(GLNodes, f2)*Nk2_interpolate(GLNodes, g2)*Nk2_interpolate(GLNodes, h2)  &
          & +            Nk2_interpolate(GLNodes, f3)*Nk2_interpolate(GLNodes, g3)*Nk2_interpolate(GLNodes, h3)  &
          & ))

  end function TriNk2_linteg_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LagrangePolyn_Init_TriNk3()

    integer :: i

    ! 実行文; Executable statements
    !
!!$    
    open(10, file='LagrangePolyn_TriNk3.dat', status='old')
    do i=1, 10
       read(10, *) TriNk3_LagrangeBasisCoef(i,:)
    end do
    close(10)

  end subroutine LagrangePolyn_Init_TriNk3

  function Nk3_interpolate_1(y, f) result(val)
    real(DP), intent(in) :: y, f(4)
    real(DP) :: val

    real(DP) :: phi(4)
    real(DP), parameter :: glPt(4) = (/ -1d0, -1d0/sqrt(5d0), 1d0/sqrt(5d0), 1d0 /)

    val = (y**2 - 1d0)/12d0 * dx_L3(y) * sum( f(:)/(L3(glPt(:))*(y - glPt(:))) )

    contains 
      real(DP) elemental function L3(x) 
        real(DP), intent(in) :: x
        L3 = 0.5d0*x*(5d0*x**2 - 3d0)
      end function L3
      real(DP) elemental function dx_L3(x) 
        real(DP), intent(in) :: x
        dx_L3 = 1.5d0*(5d0*x**2 - 1d0)
      end function dx_L3
  end function Nk3_interpolate_1

  function Nk3_interpolate_2(y, f) result(val)
    real(DP), intent(in) :: y(:), f(4)
    real(DP) :: val(size(y))

    real(DP) :: phi(4)
    real(DP), parameter :: glPt(4) = (/ -1d0, -1d0/sqrt(5d0), 1d0/sqrt(5d0), 1d0 /)
    integer :: i

    do i=1, size(y)
       val(i) = (y(i)**2 - 1d0)/12d0 * dx_L3(y(i)) * sum( f(:)/(L3(glPt(:))*(y(i) - glPt(:))) )
    end do

  contains 
    real(DP) elemental function L3(x) 
      real(DP), intent(in) :: x
      L3 = 0.5d0*x*(5d0*x**2 - 3d0)
    end function L3
    real(DP) elemental function dx_L3(x) 
      real(DP), intent(in) :: x
      dx_L3 = 1.5d0*(5d0*x**2 - 1d0)
    end function dx_L3
  end function Nk3_interpolate_2

  function TriNk3_interpolate_1(y1, y2, f) result(val)
    real(DP), intent(in) :: y1, y2, f(10)
    real(DP) :: val

    real(DP) :: y(10)

    y = (/ 1d0, y1, y2, y1**2, y1*y2, y2**2, y1**3, y1**2*y2, y1*y2**2, y2**3 /)
    val = &
         &   f(1)*(y(1)-6*(y(2)+y(3))+10*(y(4)+y(6))+21*y(5)-5*(y(7)+y(10))-16*(y(8)+y(9))) &
         & + f(2)*5d0*((Phi+1)*y(2)-(3*Phi+2)*y(4)-(2*Phi+3)*y(5)+(2*Phi+1)*y(7)+3*(Phi+1)*y(8)+(Phi+2)*y(9)) &
         & + f(3)*5d0*(-Phi*y(2)+(3*Phi+1)*y(4)+(2*Phi-1)*y(5)-(2*Phi+1)*y(7)-3*Phi*y(8)+(1-Phi)*y(9)) &
         & + f(4)*(y(2)-5*(y(4)-y(7))+(y(5)-y(8)-y(9))) &
         & + f(5)*5d0*(-y(5)+(Phi+2)*y(8)+(1-Phi)*y(9)) &
         & + f(6)*5d0*(-y(5)+(Phi+2)*y(9)+(1-Phi)*y(8)) &
         & + f(7)*(y(3)-5*(y(6)-y(10))+(y(5)-y(8)-y(9))) &
         & + f(8)*5d0*(-Phi*y(3)+(3*Phi+1)*y(6)+(2*Phi-1)*y(5)-(2*Phi+1)*y(10)-3*Phi*y(9)+(1-Phi)*y(8)) &
         & + f(9)*5d0*((Phi+1)*y(3)-(3*Phi+2)*y(6)-(2*Phi+3)*y(5)+(2*Phi+1)*y(10)+3*(Phi+1)*y(9)+(Phi+2)*y(8)) &
         & + f(10)*27d0*(y(5)-y(8)-y(9))

!!$    val = dot_product(matmul(TriNk3_LagrangeBasisCoef, y), f)

  end function TriNk3_interpolate_1


  function TriNk3_interpolate_2(y1, y2, f) result(val)
    real(DP), intent(in) :: y1(:), y2(:), f(10)
    real(DP) :: val(size(y1))

    real(DP) :: y(size(y1),10)
    integer :: i 

    do i=1, size(y1)
       y(i,:) = (/ 1d0, y1(i), y2(i), y1(i)**2, y1(i)*y2(i), y2(i)**2, y1(i)**3, y1(i)**2*y2(i), y1(i)*y2(i)**2, y2(i)**3 /)
    end do

    val = &
         &   f(1)*(y(:,1)-6*(y(:,2)+y(:,3))+10*(y(:,4)+y(:,6))+21*y(:,5)-5*(y(:,7)+y(:,10))-16*(y(:,8)+y(:,9))) &
         & + f(2)*5*((Phi+1)*y(:,2)-(3*Phi+2)*y(:,4)-(2*Phi+3)*y(:,5)+(2*Phi+1)*y(:,7)+3*(Phi+1)*y(:,8)+(Phi+2)*y(:,9)) &
         & + f(3)*5*(-Phi*y(:,2)+(3*Phi+1)*y(:,4)+(2*Phi-1)*y(:,5)-(2*Phi+1)*y(:,7)-3*Phi*y(:,8)+(1-Phi)*y(:,9)) &
         & + f(4)*(y(:,2)-5*(y(:,4)-y(:,7))+(y(:,5)-y(:,8)-y(:,9))) &
         & + f(5)*5*(-y(:,5)+(Phi+2)*y(:,8)+(1-Phi)*y(:,9)) &
         & + f(6)*5*(-y(:,5)+(Phi+2)*y(:,9)+(1-Phi)*y(:,8)) &
         & + f(7)*(y(:,3)-5*(y(:,6)-y(:,10))+(y(:,5)-y(:,8)-y(:,9))) &
         & + f(8)*5*(-Phi*y(:,3)+(3*Phi+1)*y(:,6)+(2*Phi-1)*y(:,5)-(2*Phi+1)*y(:,10)-3*Phi*y(:,9)+(1-Phi)*y(:,8)) &
         & + f(9)*5*((Phi+1)*y(:,3)-(3*Phi+2)*y(:,6)-(2*Phi+3)*y(:,5)+(2*Phi+1)*y(:,10)+3*(Phi+1)*y(:,9)+(Phi+2)*y(:,8)) &
         & + f(10)*27*(y(:,5)-y(:,8)-y(:,9))

  end function TriNk3_interpolate_2
  
  function TriNk3_basis(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(10)

    integer :: i
    real(DP) :: y(10)

    y = (/ 1d0, y1, y2, y1**2, y1*y2, y2**2, y1**3, y1**2*y2, y1*y2**2, y2**3 /)
    val = (/ y(1)-6*(y(2)+y(3))+10*(y(4)+y(6))+21*y(5)-5*(y(7)+y(10))-16*(y(8)+y(9)), &
         &   5*((Phi+1)*y(2)-(3*Phi+2)*y(4)-(2*Phi+3)*y(5)+(2*Phi+1)*y(7)+3*(Phi+1)*y(8)+(Phi+2)*y(9)), &
         &   5*(-Phi*y(2)+(3*Phi+1)*y(4)+(2*Phi-1)*y(5)-(2*Phi+1)*y(7)-3*Phi*y(8)+(1-Phi)*y(9)), &
         &   y(2)-5*(y(4)-y(7))+(y(5)-y(8)-y(9)), &
         &   5*(-y(5)+(Phi+2)*y(8)+(1-Phi)*y(9)), &
         &   5*(-y(5)+(Phi+2)*y(9)+(1-Phi)*y(8)), &
         &   y(3)-5*(y(6)-y(10))+(y(5)-y(8)-y(9)), &
         &   5*(-Phi*y(3)+(3*Phi+1)*y(6)+(2*Phi-1)*y(5)-(2*Phi+1)*y(10)-3*Phi*y(9)+(1-Phi)*y(8)), &
         &   5*((Phi+1)*y(3)-(3*Phi+2)*y(6)-(2*Phi+3)*y(5)+(2*Phi+1)*y(10)+3*(Phi+1)*y(9)+(Phi+2)*y(8)), &
         &   27*(y(5)-y(8)-y(9)) /)

  end function TriNk3_basis

  function TriNk3_basis_dy1(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(10)

    integer :: i
    real(DP) :: y(10)

    y = (/ 0d0, 1d0, 0d0, 2d0*y1, y2, 0d0, 3d0*y1**2, 2d0*y1*y2, y2**2, 0d0 /)
    val = (/ -6*y(2)+10*(y(4)+y(6))+21*y(5)-5*(y(7)+y(10))-16*(y(8)+y(9)), &
         &   5*((Phi+1)*y(2)-(3*Phi+2)*y(4)-(2*Phi+3)*y(5)+(2*Phi+1)*y(7)+3*(Phi+1)*y(8)+(Phi+2)*y(9)), &
         &   5*(-Phi*y(2)+(3*Phi+1)*y(4)+(2*Phi-1)*y(5)-(2*Phi+1)*y(7)-3*Phi*y(8)+(1-Phi)*y(9)), &
         &   y(2)-5*(y(4)-y(7))+(y(5)-y(8)-y(9)), &
         &   5*(-y(5)+(Phi+2)*y(8)+(1-Phi)*y(9)), &
         &   5*(-y(5)+(Phi+2)*y(9)+(1-Phi)*y(8)), &
         &   -5*(-y(10))+(y(5)-y(8)-y(9)), &
         &   5*((2*Phi-1)*y(5)-(2*Phi+1)*y(10)-3*Phi*y(9)+(1-Phi)*y(8)), &
         &   5*(-(2*Phi+3)*y(5)+(2*Phi+1)*y(10)+3*(Phi+1)*y(9)+(Phi+2)*y(8)), &
         &   27*(y(5)-y(8)-y(9)) /)

  end function TriNk3_basis_dy1

  function TriNk3_basis_dy2(y1, y2) result(val)
    real(DP), intent(in) :: y1, y2
    real(DP) :: val(10)

    integer :: i
    real(DP) :: y(10)

    y = (/ 0d0, 0d0, 1d0, 0d0, y1, 2d0*y2, 0d0, y1**2, 2d0*y1*y2, 3d0*y2**2 /)
    val = (/ -6*y(3)+10*y(6)+21*y(5)-5*y(10)-16*(y(8)+y(9)), &
         &   5*(-(2*Phi+3)*y(5)+3*(Phi+1)*y(8)+(Phi+2)*y(9)), &
         &   5*((2*Phi-1)*y(5)-3*Phi*y(8)+(1-Phi)*y(9)), &
         &   -5*(-y(7))+(y(5)-y(8)-y(9)), &
         &   5*(-y(5)+(Phi+2)*y(8)+(1-Phi)*y(9)), &
         &   5*(-y(5)+(Phi+2)*y(9)+(1-Phi)*y(8)), &
         &   y(3)-5*(y(6)-y(10))+(y(5)-y(8)-y(9)), &
         &   5*(-Phi*y(3)+(3*Phi+1)*y(6)+(2*Phi-1)*y(5)-(2*Phi+1)*y(10)-3*Phi*y(9)+(1-Phi)*y(8)), &
         &   5*((Phi+1)*y(3)-(3*Phi+2)*y(6)-(2*Phi+3)*y(5)+(2*Phi+1)*y(10)+3*(Phi+1)*y(9)+(Phi+2)*y(8)), &
         &   27*(y(5)-y(8)-y(9)) /)
   
  end function TriNk3_basis_dy2

  function TriNk3_sinteg_1(f) result(val)
    real(DP), intent(in) :: f(10)
    real(DP) :: val

    integer :: i

    val = 0d0
    do i=1, 12
       val = val + TriNk3_SInt_Wt(i) &
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), f)
    end do
    val = 0.5d0*val

  end function TriNk3_sinteg_1

  function TriNk3_sinteg_2(f, g) result(val)
    real(DP), intent(in) :: f(10), g(10)
    real(DP) :: val

    integer :: i

    val = 0d0
    do i=1, 12
       val = val + TriNk3_SInt_Wt(i)&
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), f) &
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), g)
    end do
    val = 0.5d0*val

  end function TriNk3_sinteg_2

  function TriNk3_sinteg_3(f, g, h) result(val)
    real(DP), intent(in) :: f(10), g(10), h(10)
    real(DP) :: val

    integer :: i

    val = 0d0    
    do i=1, 12
       val = val + TriNk3_SInt_Wt(i)&
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), f) &
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), g) &
            & *TriNk3_interpolate_1(TriNk3_SInt_y1(i), TriNk3_SInt_y2(i), h)
    end do
    val = 0.5d0*val

  end function TriNk3_sinteg_3

  function TriNk3_sinteg_dotProdWt(f) result(val)
    real(DP) :: f(12)
    real(DP) :: val
    
    val = 0.5d0*sum( TriNk3_SInt_Wt(:)*f(:) )
         
  end function TriNk3_sinteg_dotProdWt

  function TriNk3_linteg_1(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(4),f2(4),f3(4)
    real(DP) :: val
    
    real(DP), parameter, dimension(4) :: GLNodes = GLQd_order7_y(:)

    val = 0.5d0*sum( GLQd_order7_Wt(:)*( &
          &              Nk3_interpolate(GLNodes, f1) &
          & +  sqrt(2d0)*Nk3_interpolate(GLNodes, f2) &
          & +            Nk3_interpolate(GLNodes, f3) &
          & ))

  end function TriNk3_linteg_1

  function TriNk3_linteg_2(f1,g1, f2,g2, f3,g3) result(val)
    real(DP), dimension(4), intent(in) :: f1,g1, f2,g2, f3,g3
    real(DP) :: val
    
    real(DP), parameter, dimension(4) :: GLNodes = GLQd_order7_y(:)

    val = 0.5d0*sum( GLQd_order7_Wt(:)*( &
          &              Nk3_interpolate(GLNodes, f1)*Nk3_interpolate(GLNodes, g1)  &
          & +  sqrt(2d0)*Nk3_interpolate(GLNodes, f2)*Nk3_interpolate(GLNodes, g2)  &
          & +            Nk3_interpolate(GLNodes, f3)*Nk3_interpolate(GLNodes, g3)  &
          & ))

  end function TriNk3_linteg_2

  function TriNk3_linteg_3(f1,g1,h1, f2,g2,h2, f3,g3,h3) result(val)
    real(DP), dimension(4), intent(in) :: f1,g1,h1, f2,g2,h2, f3,g3,h3
    real(DP) :: val
    
    real(DP), parameter, dimension(4) :: GLNodes = GLQd_order7_y(:)

    val = 0.5d0*sum( GLQd_order7_Wt(:)*( &
          &              Nk3_interpolate(GLNodes, f1)*Nk3_interpolate(GLNodes, g1)*Nk3_interpolate(GLNodes, h1)  &
          & +  sqrt(2d0)*Nk3_interpolate(GLNodes, f2)*Nk3_interpolate(GLNodes, g2)*Nk3_interpolate(GLNodes, h2)  &
          & +            Nk3_interpolate(GLNodes, f3)*Nk3_interpolate(GLNodes, g3)*Nk3_interpolate(GLNodes, h3)  &
          & ))

  end function TriNk3_linteg_3

end module LagrangePolyn_mod


