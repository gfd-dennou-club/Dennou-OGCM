module PolyMesh_mod

  use dc_types
  use VectorSpace_mod

  implicit none
  private

  integer, parameter :: MAX_FACE_VERTEX_NUM = 6

  type, public :: Face
     integer :: vertNum
     integer :: vertIdList(MAX_FACE_VERTEX_NUM)

     integer :: ownCellId
     integer :: neighCellId
  end type Face

  integer, parameter :: MAX_CELL_FACE_NUM = 6

  type, public :: Cell
     integer :: faceNum
     integer :: faceIdList(MAX_CELL_FACE_NUM)

  end type Cell

  type, public :: PolyMesh
     type(Face), pointer :: faceList(:)
     type(Cell), pointer :: cellList(:)
     type(Vector3d), pointer :: pointList(:)
  end type PolyMesh

  public :: Face_Init, Cell_Init
  public :: PolyMesh_Init, PolyMesh_Final

contains

subroutine Face_Init(face_, vertIds)
  type(Face), intent(inout) :: face_
  integer, intent(in) :: vertIds(:)
  
  face_%vertNum = size(vertIds)
  face_%vertIdList(1:face_%vertNum) = vertIds

end subroutine Face_Init

subroutine Cell_Init(cell_, faceIds)
  type(Cell), intent(inout) :: cell_
  integer, intent(in) :: faceIds(:)
  
  cell_%faceNum = size(faceIds)
  cell_%faceIdList(1:cell_%faceNum) = faceIds

end subroutine Cell_Init


subroutine PolyMesh_Init(mesh, points, faces, cells)
  type(PolyMesh), intent(inout) :: mesh
  type(Vector3d), intent(in) :: points(:)
  type(Face), intent(in) :: faces(:)
  type(Cell), intent(in) :: cells(:)

  if(associated(mesh%pointList)) deallocate(mesh%pointList)
  if(associated(mesh%faceList)) deallocate(mesh%faceList)
  if(associated(mesh%cellList)) deallocate(mesh%cellList)

  allocate( mesh%pointList(size(points)) )
  allocate( mesh%faceList(size(faces)) )
  allocate( mesh%cellList(size(cells)) )

  mesh%cellList(:) = cells(:)
  mesh%faceList(:) = faces(:)  
  mesh%pointList(:) = points(:)

end subroutine PolyMesh_Init

subroutine PolyMesh_Final(mesh)
  type(PolyMesh), intent(inout) :: mesh

  if(associated(mesh%pointList)) deallocate(mesh%pointList)
  if(associated(mesh%faceList)) deallocate(mesh%faceList)
  if(associated(mesh%cellList)) deallocate(mesh%cellList)

end subroutine PolyMesh_Final

end module PolyMesh_mod
