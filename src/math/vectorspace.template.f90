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
    module procedure optr_add
  end interface	operator(+)

  interface operator(-)
    module procedure optr_sub
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
  end interface assignment(=)

  interface print
     module procedure printVec
  end interface print
  
  public :: operator(+), operator(-), operator(*), operator(/), assignment(=)
  public :: print

contains
subroutine printVec(v)
  type(vectorspaceTypeName), intent(in) :: v
  
  write(*,*) "[", v%v_, "]"
  
end subroutine printVec

function optr_add(v1, v2) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1, v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ + v2%v_
  
end function optr_add

!
! Provide the oparation to substract two elements of V. 
!
function optr_sub(v1, v2) result(ret)
  type(vectorspaceTypeName), intent(in) :: v1, v2
  type(vectorspaceTypeName) :: ret

  ret%v_ = v1%v_ - v2%v_
  
end function optr_sub

function optr_scalarMul(scalar, v) result(ret)
  real(DP), intent(in) :: scalar
  type(vectorspaceTypeName), intent(in) :: v
  type(vectorspaceTypeName) :: ret

  ret%v_ = scalar*v%v_
  
end function optr_scalarMul

function optr_scalarDiv(v, scalar) result(ret)
  type(vectorspaceTypeName), intent(in) :: v
  real(DP), intent(in) :: scalar
  type(vectorspaceTypeName) :: ret

  ret%v_ = v%v_/scalar
  
end function optr_scalarDiv

subroutine optr_assign(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v
  type(vectorspaceTypeName), intent(in) :: rhs
  
  v%v_ = rhs%v_

end subroutine optr_assign

subroutine optr_assign_scalar(v, rhs)

  type(vectorspaceTypeName), intent(out) :: v
  real(DP), intent(in) :: rhs
  
  v%v_ = rhs

end subroutine optr_assign_scalar

end module moduleName
