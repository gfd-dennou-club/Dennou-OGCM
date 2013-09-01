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
     type(Vector3d), pointer :: cellPosList(:) => null()
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
  public :: PolyMesh_dataAlloc
  public :: PolyMesh_SetPoint, PolyMesh_SetFace, PolyMesh_SetCell
  public :: PolyMesh_getConnectivity, PolyMesh_setConnectivity
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
  if(associated(mesh%cellPosList)) deallocate(mesh%cellPosList)

  allocate( mesh%pointList(ptNum) )
  allocate( mesh%faceList(faceNum) )
  allocate( mesh%cellList(cellNum) )
  allocate( mesh%cellPosList(cellNum) )

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


subroutine PolyMesh_setCell(mesh, cellId, cell_, cellPos)
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: cellId
  type(Cell), intent(in) :: cell_
  type(Vector3d), intent(in) :: cellPos

  mesh%cellList(cellId) = cell_
  mesh%cellPosList(cellId) = cellPos

end subroutine PolyMesh_setCell

subroutine PolyMesh_getConnectivity( &
     & mesh, cell_points, cell_faces, face_points, face_links )
  type(PolyMesh), intent(in) :: mesh
  integer, optional, allocatable :: cell_points(:,:)
  integer, optional, allocatable :: cell_faces(:,:)
  integer, optional, allocatable :: face_points(:,:)
  integer, optional, allocatable :: face_links(:,:)

  integer :: cellId, faceId, i
  
  if( present(cell_points) ) then
     allocate( cell_points(MAX_CELL_FACE_NUM, getCellListSize(mesh)) )
     
     cell_points = -1
     do cellId=1, getCellListSize(mesh)
        do faceId=1, mesh%cellList(cellId)%faceNum
           cell_points(faceId, cellId) = mesh%faceList( mesh%cellList(cellId)%faceIdList(faceId) )%vertIdList(1)
        end do
     end do
  end if

  if( present(cell_faces) ) then
     allocate( cell_faces(MAX_CELL_FACE_NUM, getCellListSize(mesh)) )
     
     cell_faces = -1
     do cellId=1, getCellListSize(mesh)
        do faceId=1, mesh%cellList(cellId)%faceNum
           cell_faces(faceId, cellId) = mesh%cellList(cellId)%faceIdList(faceId)
        end do
     end do
  end if

  if( present(face_points) ) then
     allocate( face_points(MAX_FACE_VERTEX_NUM, getFaceListSize(mesh)) )
     
     face_points = -1
     do faceId=1, getFaceListSize(mesh)
        do i=1, mesh%faceList(faceId)%vertNum
           face_points(i, faceId) = mesh%faceList(faceId)%vertIdList(i)
        end do
     end do
  end if

  if( present(face_links) ) then
     allocate( face_links(2, getFaceListSize(mesh)) )
     
     face_links = -1
     do faceId=1, getFaceListSize(mesh)
        face_links(1, faceId) = mesh%faceList(faceId)%ownCellId
        face_links(2, faceId) = mesh%faceList(faceId)%neighCellId
     end do
  end if
  
end subroutine PolyMesh_getConnectivity

subroutine PolyMesh_setConnectivity( &
     & mesh, cell_points, cell_faces, face_points, face_links )
  type(PolyMesh), intent(inout) :: mesh
  integer, intent(in) :: cell_points(:,:)
  integer, intent(in) :: cell_faces(:,:)
  integer, intent(in) :: face_points(:,:)
  integer, intent(in) :: face_links(:,:)

  integer :: cellId, faceId
  integer :: faceNum, fvertNum

  do cellId=1, getCellListSize(mesh)
     do faceNum=1, MAX_CELL_FACE_NUM-1
        if(cell_points(faceNum+1,cellId) == -1) exit
     end do
     mesh%cellList(cellId)%faceNum = faceNum
     mesh%cellList(cellId)%faceIdList(1:faceNum) = cell_points(1:faceNum, cellId)
  end do

  do faceId=1, getfaceListSize(mesh)
     mesh%faceList(faceId)%vertNum = size(face_points, 1)
     mesh%faceList(faceId)%vertIdList(:) = face_points(:, faceId)
  end do

  do faceId=1, getfaceListSize(mesh)
     mesh%faceList(faceId)%ownCellId = face_links(1,faceId)
     mesh%faceList(faceId)%neighCellId = face_links(2,faceId)
  end do

end subroutine PolyMesh_setConnectivity

end module PolyMesh_mod
