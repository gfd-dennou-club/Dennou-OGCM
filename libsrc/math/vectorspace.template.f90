#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _vectorspace_typename(sizesig, typesig) vector##sizesig##typesig
#define vectorspaceTypeName _eval2(_vectorspace_typename, vecspace_elem_size, vecspace_elem_type_sig)
#define _moduleName(s) s ## _mod
#define moduleName _eval1(_moduleName, vectorspaceTypeName)

module moduleName
  use dc_types

  implicit none
  private

  type, public :: vectorspaceTypeName
    vecspace_elem_type :: v_(vecspace_elem_size)
  end type vectorspaceTypeName

  interface operator(+)
    module procedure optr_add1
    module procedure optr_add2
    module procedure optr_add3
  end interface	operator(+)

  interface operator(-)
    module procedure optr_sub1
    module procedure optr_sub2
    module procedure optr_sub3
  end interface	operator(-)

  interface operator(*)
    module procedure optr_scalarMul
  end interface	operator(*)


  interface operator(/)
    module procedure optr_scalarDiv
  end interface	operator(/)

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
  type(vectorspaceTypeName), intent(in) :: v
  
  write(*,*) "[", v%v_, "]"
  
end subroutine printVec

!
! Provide the oparation to add two elements of V. 
!
pure function optr_add1(v1, v2) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1, v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ + v2%v_
  
end function optr_add1

pure function optr_add2(val, v2) result(ret)
  real(DP), intent(in) :: val
  type(vectorspaceTypeName), intent(in) :: v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = val + v2%v_
  
end function optr_add2

pure function optr_add3(v1, val) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1
  real(DP), intent(in) :: val
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ + val
  
end function optr_add3


!
! Provide the oparation to substract two elements of V. 
!
pure function optr_sub1(v1, v2) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1, v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ - v2%v_
  
end function optr_sub1

pure function optr_sub2(val, v2) result(ret)
  real(DP), intent(in) :: val
  type(vectorspaceTypeName), intent(in) :: v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = val - v2%v_
  
end function optr_sub2

pure function optr_sub3(v1, val) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1
  real(DP), intent(in) :: val
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ - val
  
end function optr_sub3

!
!
!
pure function optr_scalarMul(scalar, v) result(ret)
  real(DP), intent(in) :: scalar
  type(vectorspaceTypeName), intent(in) :: v
  type(vectorspaceTypeName) :: ret

  ret%v_ = scalar*v%v_
  
end function optr_scalarMul

pure function optr_scalarDiv(v, scalar) result(ret)
  type(vectorspaceTypeName), intent(in) :: v
  real(DP), intent(in) :: scalar
  type(vectorspaceTypeName) :: ret

  ret%v_ = v%v_/scalar
  
end function optr_scalarDiv

pure subroutine optr_assign(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v
  type(vectorspaceTypeName), intent(in) :: rhs
  
  v%v_ = rhs%v_

end subroutine optr_assign

pure subroutine optr_assign_scalar(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v
  real(DP), intent(in) :: rhs
  
  v%v_ = rhs

end subroutine optr_assign_scalar


pure subroutine optr_assign_vecElem(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v
  vecspace_elem_type, intent(in) :: rhs(vecspace_elem_size)
  
  v%v_ = rhs

end subroutine optr_assign_vecElem

pure subroutine optr_assign_objArray(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v(:)
  type(vectorspaceTypeName), intent(in) :: rhs(:)
  
  integer :: i
  
  do i=lbound(v,1), ubound(v,1)
     v(i) = rhs(i)
  end do

end subroutine optr_assign_objArray

function toElemTypeArray(v) result(ary)
  type(vectorspaceTypeName), intent(in) :: v  
  vecspace_elem_type :: ary(vecspace_elem_size)

  ary(:) = v%v_(:)
end function toElemTypeArray

end module moduleName
