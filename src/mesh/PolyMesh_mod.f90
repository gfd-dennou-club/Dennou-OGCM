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
     type(Face), pointer :: faceList(:) => null()
     type(Cell), pointer :: cellList(:) => null()
     type(Vector3d), pointer :: pointList(:) => null()
  end type PolyMesh

  interface PolyMesh_Init
    module procedure PolyMesh_Init1
    module procedure PolyMesh_Init2
  end interface PolyMesh_Init
  
  interface getFaceVertex
    module procedure getFaceVertex1
    module procedure getFaceVertex2
  end interface getFaceVertex
  
  public :: Face_Init, Cell_Init
  public :: PolyMesh_Init, PolyMesh_Final
  public :: PolyMesh_SetPoint, PolyMesh_SetFace, PolyMesh_SetCell
  public :: getPointListSize, getFaceListSize, getCellListSize
  public :: getFaceVertex

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

subroutine PolyMesh_Init1(mesh, ptNum, faceNum, cellNum)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: ptNum, faceNum, cellNum

  call PolyMesh_dataAlloc(mesh, ptNum, faceNum, cellNum)

end subroutine PolyMesh_Init1

subroutine PolyMesh_Init2(mesh, points, faces, cells)
  type(PolyMesh), intent(inout) :: mesh
  type(Vector3d), intent(in) :: points(:)
  type(Face), intent(in) :: faces(:)
  type(Cell), intent(in) :: cells(:)

  call PolyMesh_dataAlloc(mesh, size(points), size(faces), size(cells))

  mesh%cellList(:) = cells(:)
  mesh%faceList(:) = faces(:)  
  mesh%pointList(:) = points(:)

end subroutine PolyMesh_Init2

function getFaceListSize(mesh) result(listSize)
  type(PolyMesh), intent(in) :: mesh
  integer :: listSize

  listSize = size(mesh%faceList)
end function getFaceListSize

function getPointListSize(mesh) result(listSize)
  type(PolyMesh), intent(in) :: mesh
  integer :: listSize

  listSize = size(mesh%PointList)
end function getPointListSize

function getCellListSize(mesh) result(listSize)
  type(PolyMesh), intent(in) :: mesh
  integer :: listSize

  listSize = size(mesh%cellList)
end function getCellListSize

function getFaceVertex1(mesh, faceGId) result(vxs)
  type(PolyMesh), intent(in) :: mesh
  integer, intent(in) :: faceGId
  type(Vector3d) :: vxs(mesh%FaceList(faceGId)%vertNum)

  type(Face) :: face_

  face_ = mesh%FaceList(faceGId)
  vxs(:) = mesh%PointList( face_%vertIdList(1:face_%vertNum) )  

end function getFaceVertex1

function getFaceVertex2(mesh, cellGId, faceLId) result(vxs)
  type(PolyMesh), intent(in) :: mesh
  integer, intent(in) :: cellGId
  integer, intent(in) :: faceLId
  type(Vector3d) :: vxs(mesh%FaceList(mesh%cellList(cellGId)%faceIdList(faceLId))%vertNum)

  integer :: faceGId
  type(Vector3d) :: tmp
  type(Face) :: face_


  faceGId = mesh%cellList(cellGId)%faceIdList(faceLId)
  face_ = mesh%FaceList(faceGId)
  vxs(:) = mesh%PointList( face_%vertIdList(1:face_%vertNum) )  

  if( face_%neighCellId == cellGId ) then
    ! Swap
    tmp = vxs(1); vxs(1) = vxs(2); vxs(2) = tmp 
  end if

end function getFaceVertex2

subroutine PolyMesh_dataAlloc(mesh, ptNum, faceNum, cellNum)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: ptNum, faceNum, cellNum

  if(associated(mesh%pointList)) deallocate(mesh%pointList)
  if(associated(mesh%faceList)) deallocate(mesh%faceList)
  if(associated(mesh%cellList)) deallocate(mesh%cellList)

  allocate( mesh%pointList(ptNum) )
  allocate( mesh%faceList(faceNum) )
  allocate( mesh%cellList(cellNum) )

end subroutine PolyMesh_dataAlloc

subroutine PolyMesh_Final(mesh)
  type(PolyMesh), intent(inout) :: mesh

  if(associated(mesh%pointList)) deallocate(mesh%pointList)
  if(associated(mesh%faceList)) deallocate(mesh%faceList)
  if(associated(mesh%cellList)) deallocate(mesh%cellList)

end subroutine PolyMesh_Final

subroutine PolyMesh_setPoint(mesh, ptId, pt)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: ptId
  type(Vector3d), intent(in) :: pt

  mesh%pointList(ptId) = pt

end subroutine PolyMesh_setPoint

subroutine PolyMesh_setFace(mesh, faceId, face_)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: faceId
  type(Face), intent(in) :: face_

  mesh%faceList(faceId) = face_

end subroutine PolyMesh_setFace


subroutine PolyMesh_setCell(mesh, cellId, cell_)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: cellId
  type(Cell), intent(in) :: cell_

  mesh%cellList(cellId) = cell_

end subroutine PolyMesh_setCell

end module PolyMesh_mod
