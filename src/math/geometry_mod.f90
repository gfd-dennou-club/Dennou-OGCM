module geometry_mod
  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod

  implicit none
  private

  public :: icosahedron_vertex

contains
function icosahedron_vertex() result(orth_icvertex)
 
  
  ! 宣言文 ; Declaration statements
  !
  type(Vector3d) :: orth_icvertex(12)

  ! 作業変数
  ! Work variable
  !
  real(DP), parameter :: PI = acos(-1d0)
  real(DP) :: unit_radius = 1.0d0
  real(DP) :: dLambda = 2.0d0 * PI / 5.0d0
  integer :: i

  ! 実行文 ; Executable statements
  ! 
  
  orth_icvertex(1)%v_(:) = (/ 0.0d0, 0.0d0, unit_radius /)
  orth_icvertex(12)%v_(:) = (/ 0.0d0, 0.0d0, - unit_radius /)

  do i=0,4
     orth_icvertex(i+2) = SphToCartPos(dlambda*i, PI/6d0)
     orth_icvertex(i+7) = SphToCartPos(dlambda*(i+0.5d0), -PI/6d0)
  end do

end function icosahedron_vertex

end module geometry_mod
