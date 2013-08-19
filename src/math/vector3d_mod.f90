module vector3d_mod
  use dc_types
  implicit none
  private
  type, public :: vector3d
    real(DP) :: v_(3)
  end type vector3d
  interface operator(+)
    module procedure optr_add
  end interface operator(+)
  interface operator(-)
    module procedure optr_sub
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
  end interface assignment(=)
  interface print
     module procedure printVec
  end interface print
  public :: operator(+), operator(-), operator(*), operator(/), assignment(=)
  public :: print
contains
subroutine printVec(v)
  type(vector3d), intent(in) :: v
  write(*,*) "[", v%v_, "]"
end subroutine printVec
function optr_add(v1, v2) result(ret)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: ret
  ret%v_ = v1%v_ + v2%v_
end function optr_add
!
! Provide the oparation to substract two elements of V.
!
function optr_sub(v1, v2) result(ret)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: ret
  ret%v_ = v1%v_ - v2%v_
end function optr_sub
function optr_scalarMul(scalar, v) result(ret)
  real(DP), intent(in) :: scalar
  type(vector3d), intent(in) :: v
  type(vector3d) :: ret
  ret%v_ = scalar*v%v_
end function optr_scalarMul
function optr_scalarDiv(v, scalar) result(ret)
  type(vector3d), intent(in) :: v
  real(DP), intent(in) :: scalar
  type(vector3d) :: ret
  ret%v_ = v%v_/scalar
end function optr_scalarDiv
subroutine optr_assign(v, rhs)
  type(vector3d), intent(out) :: v
  type(vector3d), intent(in) :: rhs
  v%v_ = rhs%v_
end subroutine optr_assign
subroutine optr_assign_scalar(v, rhs)
  type(vector3d), intent(out) :: v
  real(DP), intent(in) :: rhs
  v%v_ = rhs
end subroutine optr_assign_scalar
end module vector3d_mod
