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
  integer :: nk, k
  real(DP) :: y1, y2, f(nNode), g(nNode), h(nNode), I(nNode)
  type(vector2d) :: check_pos(DGElemInfo%nSIntNode)
  real(DP) :: pos_y1, pos_y2, check_val

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


  
  check_pos(1) = (/ 1d0/4d0, 3d0/4d0 /)
  check_pos(2) = (/ 3d0/4d0, 1d0/4d0 /)
  check_pos(3) = (/ 1d0/2d0, 1d0/4d0 /)
  do k=1, DGElemInfo%nSIntNode
     check_pos(k) = DGElemInfo%sIntNode(k)
  end do

  do k=1, size(check_pos)
     pos_y1 = check_pos(k)%v_(1)
     pos_y2 = check_pos(k)%v_(2)

     check_val = abs(func_basis1(check_pos(k))-TriNk_interpolate(pos_y1, pos_y2, f))
     call MessageNotify('M', 'check_interpolate', 'check value=%f', d=(/ check_val /))
     call AssertLessThan(message="check_interpolate",  answer=1d-10, check=check_val)

     check_val = abs(func_basis1_dy1(check_pos(k))-dot_product(TriNk_basis_dy1(pos_y1, pos_y2), f))
     call MessageNotify('M', 'check_basis_dy1', 'check value=%f', d=(/ check_val /))
     call AssertLessThan(message="check_basis1_dy1", answer=1d-10, check=check_val) 

     check_val = abs(func_basis1_dy2(check_pos(k))-dot_product(TriNk_basis_dy2(pos_y1, pos_y2), f))
     call MessageNotify('M', 'check_basis_dy2', 'check value=%f', d=(/ check_val /))
     call AssertLessThan(message="check_basis1_dy2", answer=1d-10, check=check_val) 

  end do

!write(*,*) TriNk_sinteg(g,h,I)
  check_val = abs(TriNk_sinteg(g,h,I)-1d0/8d0)
  call MessageNotify('M', 'check_sinteg', 'check value=%f', d=(/ check_val /))
  call AssertLessThan(message="check_basis1_sinteg",  answer=1d-9,  check=check_val )

  check_val = abs(TriNk_sinteg(h,g,I)-1d0/8d0)
  call MessageNotify('M', 'check_sinteg', 'check value=%f', d=(/ check_val /))
  call AssertLessThan(message="check_basis1_sinteg",  answer=1d-9,  check=check_val )

  check_val = abs(TriNk_sinteg(h,I,g)-1d0/8d0)
  call MessageNotify('M', 'check_sinteg', 'check value=%f', d=(/ check_val /))
  call AssertLessThan(message="check_basis1_sinteg",  answer=1d-9,  check=check_val )

  
end subroutine check_basis1

end program DGBaseMod_test
