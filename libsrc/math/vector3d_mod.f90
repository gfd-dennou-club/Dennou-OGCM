module vector3d_mod
  use dc_types
  implicit none
  private
  type, public :: vector3d
    real(DP) :: v_(3)
  end type vector3d
  interface operator(+)
    module procedure optr_add1
    module procedure optr_add2
    module procedure optr_add3
  end interface operator(+)
  interface operator(-)
    module procedure optr_sub1
    module procedure optr_sub2
    module procedure optr_sub3
  end interface operator(-)
  interface operator(*)
    module procedure optr_scalarMul
  end interface operator(*)
  interface operator(/)
    module procedure optr_scalarDiv
  end interface operator(/)
  interface assignment(=)
    module procedure optr_assign
    module procedure optr_assign_scalar
    module procedure optr_assign_vecElem
    module procedure optr_assign_objArray
  end interface assignment(=)
  interface toArray
     module procedure toElemTypeArray
  end interface toArray
  interface print
     module procedure printVec
  end interface print
  public :: operator(+), operator(-), operator(*), operator(/), assignment(=)
  public :: toArray
  public :: print
contains
subroutine printVec(v)
  type(vector3d), intent(in) :: v
  write(*,*) "[", v%v_, "]"
end subroutine printVec
!
! Provide the oparation to add two elements of V.
!
pure function optr_add1(v1, v2) result(ret)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: ret
  ret%v_ = v1%v_ + v2%v_
end function optr_add1
pure function optr_add2(val, v2) result(ret)
  real(DP), intent(in) :: val
  type(vector3d), intent(in) :: v2
  type(vector3d) :: ret
  ret%v_ = val + v2%v_
end function optr_add2
pure function optr_add3(v1, val) result(ret)
  type(vector3d), intent(in) :: v1
  real(DP), intent(in) :: val
  type(vector3d) :: ret
  ret%v_ = v1%v_ + val
end function optr_add3
!
! Provide the oparation to substract two elements of V.
!
pure function optr_sub1(v1, v2) result(ret)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: ret
  ret%v_ = v1%v_ - v2%v_
end function optr_sub1
pure function optr_sub2(val, v2) result(ret)
  real(DP), intent(in) :: val
  type(vector3d), intent(in) :: v2
  type(vector3d) :: ret
  ret%v_ = val - v2%v_
end function optr_sub2
pure function optr_sub3(v1, val) result(ret)
  type(vector3d), intent(in) :: v1
  real(DP), intent(in) :: val
  type(vector3d) :: ret
  ret%v_ = v1%v_ - val
end function optr_sub3
!
!
!
pure function optr_scalarMul(scalar, v) result(ret)
  real(DP), intent(in) :: scalar
  type(vector3d), intent(in) :: v
  type(vector3d) :: ret
  ret%v_ = scalar*v%v_
end function optr_scalarMul
pure function optr_scalarDiv(v, scalar) result(ret)
  type(vector3d), intent(in) :: v
  real(DP), intent(in) :: scalar
  type(vector3d) :: ret
  ret%v_ = v%v_/scalar
end function optr_scalarDiv
pure subroutine optr_assign(v, rhs)
  type(vector3d), intent(out) :: v
  type(vector3d), intent(in) :: rhs
  v%v_ = rhs%v_
end subroutine optr_assign
pure subroutine optr_assign_scalar(v, rhs)
  type(vector3d), intent(out) :: v
  real(DP), intent(in) :: rhs
  v%v_ = rhs
end subroutine optr_assign_scalar
pure subroutine optr_assign_vecElem(v, rhs)
  type(vector3d), intent(out) :: v
  real(DP), intent(in) :: rhs(3)
  v%v_ = rhs
end subroutine optr_assign_vecElem
pure subroutine optr_assign_objArray(v, rhs)
  type(vector3d), intent(out) :: v(:)
  type(vector3d), intent(in) :: rhs(:)
  integer :: i
  do i=lbound(v,1), ubound(v,1)
     v(i) = rhs(i)
  end do
end subroutine optr_assign_objArray
function toElemTypeArray(v) result(ary)
  type(vector3d), intent(in) :: v
  real(DP) :: ary(3)
  ary(:) = v%v_(:)
end function toElemTypeArray
end module vector3d_mod
