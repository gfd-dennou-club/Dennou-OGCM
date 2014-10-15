program DGBaseMod_test

  ! モジュール引用; Use statement
  !
  use dc_test
  use dc_message
  use dc_string
  use dc_types

  use VectorSpace_mod
  use DGElement_mod
  use LagrangePolyn_mod

  implicit none

  type(DGTriElement) :: DGElemInfo
  integer, parameter :: PolyDegree = 3
  integer :: nNode

  call LagrangePolyn_Init()
  call DGTriElement_Init(DGElemInfo, PolyDegree)
  nNode = DGElemInfo%nNode

  call check_basis1

contains

function func_basis1(y) result(val)
  type(vector2d), intent(in) :: y
  real(DP) :: val
  val = y%v_(1)**2 + y%v_(2)**2 + (0.5d0 - y%v_(1))*(0.25d0 - y%v_(2))
end function func_basis1

function func_basis1_dy1(y) result(val)
  type(vector2d), intent(in) :: y
  real(DP) :: val
  val = 2*y%v_(1) - (0.25d0 - y%v_(2))
end function func_basis1_dy1

function func_basis1_dy2(y) result(val)
  type(vector2d), intent(in) :: y
  real(DP) :: val
  val = 2*y%v_(2) - (0.5d0 - y%v_(1))
end function func_basis1_dy2

subroutine check_basis1
  integer :: nk
  real(DP) :: y1, y2, f(nNode), g(nNode), h(nNode), I(nNode)
  type(vector2d) :: pos

  real(DP) :: f_ip

  do nk=1,nNode
!     call print(DGElemInfo%node(nk))
     y1 = DGElemInfo%node(nk)%v_(1)
     y2 = DGElemInfo%node(nk)%v_(2)
     f(nk) = func_basis1(DGElemInfo%node(nk))
     g(nk) = y1
     h(nk) = 1-y2
     I(nk) = 1d0
!     write(*,*) "nk=",nk,":", TriNk_basis(y1,y2)
  end do

pos =(/0d0,0d0/)
f_ip = TriNk_interpolate(0d0, 0d0, f)
write(*,*) func_basis1(pos), f_ip

  pos = (/ 1d0/3d0, 1d0/3d0 /)
  call AssertLessThan(message="check_basis1",  &
       & answer=1d-10, check=abs(func_basis1(pos)-TriNk_interpolate(1d0/3d0, 1d0/3d0, f)) &
       & )

  call AssertLessThan(message="check_basis1_dy1",  &
       & answer=1d-10, check=abs(func_basis1_dy1(pos)-dot_product(TriNk_basis_dy1(1d0/3d0, 1d0/3d0), f)) &
       & )

  call AssertLessThan(message="check_basis1_dy2",  &
       & answer=1d-10, check=abs(func_basis1_dy2(pos)-dot_product(TriNk_basis_dy2(1d0/3d0, 1d0/3d0), f)) &
       & )

!write(*,*) TriNk_sinteg(g,h,I)
  call AssertLessThan(message="check_basis1_sinteg",  &
       & answer=1d-9, check=abs(TriNk_sinteg(g,h,I)-1d0/8d0) &
       & )
  call AssertLessThan(message="check_basis1_sinteg",  &
       & answer=1d-9, check=abs(TriNk_sinteg(g,I,h)-1d0/8d0) &
       & )
  call AssertLessThan(message="check_basis1_sinteg",  &
       & answer=1d-9, check=abs(TriNk_sinteg(I,h,g)-1d0/8d0) &
       & )

  
end subroutine check_basis1

end program DGBaseMod_test
