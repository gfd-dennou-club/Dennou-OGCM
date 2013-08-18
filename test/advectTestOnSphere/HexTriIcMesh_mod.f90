module HexTriIcMesh_mod

  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod

  implicit none
  private

  type, public :: HexTriIcMesh
     integer :: glevel
     type(PolyMesh) :: mesh
  end type HexTriIcMesh

  type :: VorRegion
     integer :: siteID
     integer :: nodeIDs(6)
     integer :: edgesID(6)
  end type VorRegion
  
  public :: HexTriIcMesh_Init, HexTriIcMesh_Final
  public :: HexTriIcMesh_generate

contains
subroutine HexTriIcMesh_Init(mesh, glevel_)
  type(HexTriIcMesh), intent(inout) :: mesh
  integer, intent(in) :: glevel_

  type(Vector3d), allocatable :: points(:)
  type(Face), allocatable :: faces(:)
  type(Cell), allocatable :: cells(:)

  mesh%glevel = glevel_
  !call PolyMesh_Init(mesh%mesh, points, faces, cells)

end subroutine HexTriIcMesh_Init

subroutine HexTriIcMesh_Final(mesh)
  type(HexTriIcMesh), intent(inout) :: mesh

  call PolyMesh_Final(mesh%mesh)

end subroutine HexTriIcMesh_Final

subroutine HexTriIcMesh_generate(mesh)

  use SVoronoiGen_mod

  type(HexTriIcMesh), intent(inout) :: mesh

  type(Vector3d), allocatable :: pts(:)

  call construct_icosahedralGrid(mesh%glevel, pts)

  call SVoronoiGen_Init(size(pts))
  call SVoronoiDiagramGen(pts, (/ 1, 7, 8, 10 /))
  call SVoronoiGen_Final()

end subroutine HexTriIcMesh_generate

subroutine construct_icosahedralGrid(glevel, pts)
  integer, intent(in) :: glevel
  type(Vector3d), intent(inout), allocatable :: pts(:)

  integer :: icgridNum
  type(Vector3d) :: icLv0Vx(12)

  icgridNum = 10*4**glevel + 2
  allocate(pts(icgridNum))

  icLv0Vx(:) = icosahedron_vertex()
  
  !
  pts(1:12) = icLv0Vx(:)

end subroutine construct_icosahedralGrid

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

end module HexTriIcMesh_mod
