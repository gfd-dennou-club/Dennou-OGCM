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
     module procedure TriNk3_linteg
#elif defined POLYORDER_Nk2
     module procedure TriNk2_linteg
#else
     module procedure TriNk1_linteg
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
  real(DP), public, parameter :: TriNk1_y1(3) = (/ 0d0, 1d0, 0d0 /)
  real(DP), public, parameter :: TriNk1_y2(3) = (/ 0d0, 0d0, 1d0 /)

  real(DP), public, parameter :: TriNk2_y1(6) = (/ 0d0, 0.5d0, 1d0, 0.5d0, 0d0,   0d0 /)
  real(DP), public, parameter :: TriNk2_y2(6) = (/ 0d0,   0d0, 0d0, 0.5d0, 1d0, 0.5d0 /)

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


  real(DP), public, parameter :: TriNk1_SInt_y1(3) = (/ 1d0/6d0, 2d0/3d0, 1d0/6d0 /)
  real(DP), public, parameter :: TriNk1_SInt_y2(3) = (/ 1d0/6d0, 1d0/6d0, 2d0/3d0 /)
  real(DP), public, parameter :: TriNk1_SInt_Wt(3) = (/ 1d0/3d0, 1d0/3d0, 1d0/3d0 /)


  real(DP), public, parameter :: TriNk2_SInt_y1(6) = &
         & (/ 0.44594849091597, 0.44594849091597, 0.10810301816807, &
         &    0.09157621350977, 0.09157621350977, 0.81684757298046 /)
  real(DP), public, parameter :: TriNk2_SInt_y2(6) = &
         & (/ 0.44594849091597, 0.10810301816807, 0.44594849091597, &
         &    0.09157621350977, 0.81684757298046, 0.09157621350977 /)
  real(DP), public, parameter :: TriNk2_SInt_Wt(6) = &
         & (/ 0.22338158967801, 0.22338158967801, 0.22338158967801, &
         &    0.10995174365532, 0.10995174365532, 0.10995174365532 /)

  real(DP), public, parameter :: TriNk3_SInt_y1(12) = &
         & (/ 0.24928674517091, 0.24928674517091, 0.50142650965818, 0.06308901449150, 0.06308901449150, 0.87382197101700, &
         &    0.31035245103378, 0.63650249912140, 0.05314504984482, 0.63650249912140, 0.31035245103378, 0.05314504984482 /)
  real(DP), public, parameter :: TriNk3_SInt_y2(12) = &
         & (/ 0.24928674517091, 0.50142650965818, 0.24928674517091, 0.06308901449150, 0.87382197101700, 0.06308901449150,  &
         &    0.63650249912140, 0.05314504984482, 0.31035245103378, 0.31035245103378, 0.05314504984482, 0.63650249912140  /)
  real(DP), public, parameter :: TriNk3_SInt_Wt(12) = &
         & (/ 0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021, 0.05084490637021, 0.05084490637021, &
         &    0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837  /)

  
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

  function TriNk1_linteg(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(2),f2(2),f3(2)
    real(DP) :: val
    
    val = 0.5d0*(sum(f1+f3) + sqrt(2d0)*sum(f2))

  end function TriNk1_linteg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LagrangePolyn_Init_TriNk2()

    ! 実行文; Executable statements
    !

  end subroutine LagrangePolyn_Init_TriNk2


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

  function TriNk2_linteg(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(3),f2(3),f3(3)
    real(DP) :: val
    
    real(DP), parameter :: w1 = 1d0/3d0
    real(DP), parameter :: w0 = 4d0/3d0
    real(DP) :: f2_(3)

    f2_ = sqrt(2d0)*f2
    val = 0.5d0*( &
         &   w1*(f1(1) + f2_(1) + f3(1) + f1(3) + f2_(3) + f3(3)) &
         & + w0*(f1(2) + f2_(2) + f3(2)) &
         & )

  end function TriNk2_linteg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LagrangePolyn_Init_TriNk3()

    integer :: i

    ! 実行文; Executable statements
    !
!!$    
!!$    open(10, file='LagrangePolyn_TriNk3.dat', status='old')
!!$    do i=1, 10
!!$       read(10, *) TriNk3_LagrangeBasisCoef(:,i)
!!$    end do
!!$    close(10)

  end subroutine LagrangePolyn_Init_TriNk3


  function TriNk3_interpolate_1(y1, y2, f) result(val)
    real(DP), intent(in) :: y1, y2, f(10)
    real(DP) :: val

    real(DP) :: y(10)

    y = (/ 1d0, y1, y2, y1**2, y1*y2, y2**2, y1**3, y1**2*y2, y1*y2**2, y2**3 /)
    val = &
         &   f(1)*(y(1)-6*(y(2)+y(3))+10*(y(4)+y(6))+21*y(5)-5*(y(7)+y(10))-16*(y(8)+y(9))) &
         & + f(2)*5*((Phi+1)*y(2)-(3*Phi+2)*y(4)-(2*Phi+3)*y(5)+(2*Phi+1)*y(7)+3*(Phi+1)*y(8)+(Phi+2)*y(9)) &
         & + f(3)*5*(-Phi*y(2)+(3*Phi+1)*y(4)+(2*Phi-1)*y(5)-(2*Phi+1)*y(7)-3*Phi*y(8)+(1-Phi)*y(9)) &
         & + f(4)*(y(2)-5*(y(4)-y(7))+(y(5)-y(8)-y(9))) &
         & + f(5)*5*(-y(5)+(Phi+2)*y(8)+(1-Phi)*y(9)) &
         & + f(6)*5*(-y(5)+(Phi+2)*y(9)+(1-Phi)*y(8)) &
         & + f(7)*(y(3)-5*(y(6)-y(10))+(y(5)-y(8)-y(9))) &
         & + f(8)*5*(-Phi*y(3)+(3*Phi+1)*y(6)+(2*Phi-1)*y(5)-(2*Phi+1)*y(10)-3*Phi*y(9)+(1-Phi)*y(8)) &
         & + f(9)*5*((Phi+1)*y(3)-(3*Phi+2)*y(6)-(2*Phi+3)*y(5)+(2*Phi+1)*y(10)+3*(Phi+1)*y(9)+(Phi+2)*y(8)) &
         & + f(10)*27*(y(5)-y(8)-y(9))

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
    real(DP), intent(in) :: f(12)
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
    real(DP), intent(in) :: f(12), g(12)
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
    real(DP), intent(in) :: f(12), g(12), h(12)
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

  function TriNk3_linteg(f1,f2,f3) result(val)
    real(DP), intent(in) :: f1(4),f2(4),f3(4)
    real(DP) :: val
    
    real(DP), parameter :: w1 = 1d0/6d0
    real(DP), parameter :: w2 = 5d0/6d0
    real(DP) :: f2_(4)

    f2_ = sqrt(2d0)*f2
    val = 0.5d0*( &
         &   w1*(f1(1) + f2_(1) + f3(1) + f1(4) + f2_(4) + f3(4)) &
         & + w2*(f1(2) + f2_(2) + f3(2) + f1(3) + f2_(3) + f3(3)) &
         & )

  end function TriNk3_linteg

end module LagrangePolyn_mod


