!
!
module HexTriIcMesh_mod

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod
  use SphericalCoord_mod
  use geometry_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use GeometricField_mod

  implicit none
  private

  type, public :: HexTriIcLocalMesh
     integer :: localMeshID
     type(PolyMesh) :: mesh
  end type HexTriIcLocalMesh

  type, public :: HexTriIcMesh
     type(PolyMesh), pointer :: globalMesh => null()
     type(HexTriIcLocalMesh), pointer :: localMeshs(:)
     logical :: internalGenFlag = .true.
     real(DP) :: radius = 1d0
  end type HexTriIcMesh

  
  public :: HexTriIcMesh_Init, HexTriIcMesh_Final
  public :: HexTriIcMesh_getLocalMesh
  public :: HexTriIcMesh_generate, HexTriIcMesh_ConfigFvMeshInfo

  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'HexTriMesh_mod' !< Module Name

contains
subroutine HexTriIcMesh_Init(htiMesh, globalMesh, radius, nLocalPMesh)
  type(HexTriIcMesh), intent(inout) :: htiMesh
  type(PolyMesh), target, optional :: globalMesh
  real(DP), intent(in), optional :: radius
  integer, intent(in), optional :: nLocalPMesh

  integer :: lcPMeshNum

  if(present(globalMesh)) then
     htiMesh%globalMesh => globalMesh
     htiMesh%internalGenFlag = .false.
  end if

  if( present(radius) ) then
     htiMesh%radius = radius
     if(.not. htiMesh%internalGenFlag) call projectPosVecIntoSphere(htiMesh)
  end if

  !
  lcPMeshNum = 0
  if( present(nLocalPMesh) ) lcPMeshNum = nLocalPMesh  
  allocate(htiMesh%localMeshs(lcPMeshNum))

end subroutine HexTriIcMesh_Init

subroutine HexTriIcMesh_Final(htiMesh, radius)
  type(HexTriIcMesh), intent(inout) :: htiMesh
  real(DP), intent(in), optional :: radius

  if( htiMesh%internalGenFlag ) then
     call PolyMesh_Final(htiMesh%globalMesh)
     deallocate(htiMesh%globalMesh)
  end if

  deallocate(htiMesh%localMeshs)

end subroutine HexTriIcMesh_Final

subroutine HexTriIcMesh_generate(htiMesh, glevel, scvMaxItrNum)

  use SCVoronoiGen_mod
!  use SVoronoiGen2_mod

  type(HexTriIcMesh), intent(inout) :: htiMesh
  integer, intent(in) :: glevel
  integer, intent(in), optional :: scvMaxItrNum

  integer :: i
  integer :: siteNum
  integer :: vertNum
  integer :: maxItrNum

  siteNum = 10*4**glevel + 2
  vertNum = 2*siteNum - 4

  !
  if( htiMesh%internalGenFlag ) then
     allocate(htiMesh%globalMesh)
     call PolyMesh_Init(htiMesh%globalMesh, vertNum, siteNum+vertNum-2, siteNum)
     call construct_icosahedralGrid(glevel, htiMesh%globalMesh%cellPosList)
  end if

  !
  !

  call SCVoronoiGen_Init()

  maxItrNum = 1
  if(present(scvMaxItrNum)) maxItrNum = scvMaxItrNum
  call SCVoroniDiagram_Generate(htiMesh%globalMesh%cellPosList, getInit4IcGridPtId(glevel), maxItrNum)
  call SCVoronoi_SetTopology(htiMesh%globalMesh)

  call SCVoronoiGen_Final()

  !
  call projectPosVecIntoSphere(htiMesh)

end subroutine HexTriIcMesh_generate

function HexTriIcMesh_getLocalMesh(htiMesh, localMeshID) result(localMesh)
  type(HexTriIcMesh), intent(in), target :: htiMesh
  integer, intent(in) :: localMeshID
  type(HexTriIcLocalMesh), pointer :: localMesh

  if( localMeshID <= 0 .or. localMeshID > size(htiMesh%localMeshs) ) &
       & call MessageNotify("E", module_name, "Specified ID for local mesh is invalid.")

  localMesh => htiMesh%localMeshs(localMeshID)

end function HexTriIcMesh_getLocalMesh

subroutine HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo, localMeshID)
  type(HexTriIcMesh), intent(in), target :: htiMesh
  type(fvMeshInfo), intent(inout) :: fvInfo
  integer, intent(in), optional :: localMeshID

  type(volScalarField) :: v_cellVolume
  type(surfaceVectorField) :: s_faceAreaVec
  type(surfaceVectorField) :: s_faceCenter
  type(surfaceScalarField) :: s_dualMeshFaceArea
  type(PointScalarField) :: p_dualMeshCellVol

  integer :: cellGId, faceGId, pointGId, faceLId, faceNum
  type(HexTriIcLocalMesh), pointer :: localMesh
  type(PolyMesh), pointer :: mesh
  real(DP) :: areas(6), faceArea
  type(Vector3d) :: faceNormal, edgeVxs(2)
  type(Face), pointer :: face_
  integer :: cellIds(2), ptIds(2)

  if( present(localMeshID) ) then
     localMesh => HexTriIcMesh_getLocalMesh(htiMesh, localMeshID)
     mesh => localMesh%mesh
  else
     mesh => htiMesh%globalMesh
  end if

  !
  call GeometricField_Init( v_cellVolume, mesh, "v_cellVolume")
  call GeometricField_Init( s_faceAreaVec, mesh, "s_faceAreaVec") 
  call GeometricField_Init( s_faceCenter, mesh, "s_faceCenter") 
  call GeometricField_Init( s_dualMeshFaceArea, mesh, "s_dualMeshFaceArea")
  call GeometricField_Init( p_dualMeshCellVol, mesh, "p_dualMeshCellVol")

  !

  !
  do cellGId=1, getCellListSize(mesh)
    faceNum = mesh%CellList(cellGId)%faceNum
    do faceLId=1, faceNum
      edgeVxs(:) = getFaceVertex(mesh, cellGId, faceLId)
      areas(faceLId) = sphericalTriArea( mesh%cellPosList(cellGId), edgeVxs(1), edgeVxs(2) ) 
    end do
    
    v_CellVolume%data%v_(1,cellGId) = sum(areas(1:faceNum))

  end do


  do faceGId=1, getFaceListSize(mesh)
     face_ => mesh%faceList(faceGId)
     edgeVxs(:) = getFaceVertex(mesh, faceGId)

!!     s_faceCenter%data%v_(faceGId) = htimesh%radius * normalizedVec( 0.5d0*(edgeVxs(1) + edgeVxs(2)) )
     cellIds(:) = fvInfo%Face_CellId(1:2, faceGId)
     s_faceCenter%data%v_(1,faceGId) = htimesh%radius * &
          & normalizedVec( 0.5d0*(mesh%CellPosList(cellIds(1)) + mesh%CellPosList(cellIds(2))) )

     ptIds(:) = fvInfo%Face_PointId(1:2,faceGId)
     faceArea = geodesicArcLength( mesh%PointPosList(ptIds(1)), mesh%PointPosList(ptIds(2)) )
!!     faceNormal = normalizedVec( (edgeVxs(2) - edgeVxs(1)) .cross. (mesh%cellPosList(face_%ownCellId)) )
     faceNormal = normalizedVec( mesh%CellPosList(cellIds(2)) - mesh%CellPosList(cellIds(1)) )
     s_faceAreaVec%data%v_(1,faceGId) =  faceArea * faceNormal



     s_dualMeshFaceArea%data%v_(1,faceGId) = &
          & geodesicArcLength(mesh%cellPosList(face_%ownCellId), mesh%cellPosList(face_%neighCellId))
  end do

  do pointGId=1, getPointListSize(mesh)
     do faceLId=1, 3
        cellIds(:) = fvInfo%Face_CellId(1:2, fvInfo%Point_FaceId(faceLId,pointGId))
        areas(faceLId) = sphericalTriArea( mesh%pointPosList(pointGId), mesh%cellPosList(cellIds(1)), mesh%cellPosList(cellIds(2)) )
     end do
     p_dualMeshCellVol%data%v_(1,pointGId) = sum(areas(1:3))
  end do

  !
  call fvMeshInfo_prepair(fvInfo, &
       & v_cellVolume, s_faceAreaVec, s_faceCenter, &
       & s_dualMeshFaceArea, p_dualMeshCellVol)

  !
  call GeometricField_Final(v_cellVolume)
  call GeometricField_Final(s_faceAreaVec)
  call GeometricField_Final(s_faceCenter)
  call GeometricField_Final(s_dualMeshFaceArea)
  call GeometricField_Final(p_dualMeshCellVol)

end subroutine Hextriicmesh_configfvMeshInfo

!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine projectPosVecIntoSphere(htimesh)

  type(HexTriIcMesh), target, intent(inout) :: htimesh

  integer :: cellNum, pointNum
  integer :: i
  type(PolyMesh), pointer :: mesh

  mesh => htimesh%globalMesh
  cellNum = getCellListSize(mesh)
  pointNum = getPointListSize(mesh)

  do i=1, cellNum
     mesh%cellPosList(i) = htimesh%radius * normalizedVec(mesh%cellPosList(i))
  end do

  do i=1, pointNum
     mesh%pointPosList(i) = htimesh%radius * normalizedVec(mesh%pointPosList(i))
  end do

end subroutine projectPosVecIntoSphere

function getInit4IcGridPtId(glevel) result(iniPtsId4)
  integer, intent(in) :: glevel
  integer :: iniPtsId4(4)

  integer :: lcEnd

  lcEnd = 2**glevel+1

  iniPtsId4(:) = (/ 1, 2 + (lcEnd-1)*(6*(lcEnd-1) - 1), &
    & 2 + (lcEnd-1)*(7*(lcEnd-1) - 1), 2 + (lcEnd-1)*(9*(lcEnd-1) - 1) /)  

end function getInit4IcGridPtId

subroutine construct_icosahedralGrid(glevel, pts)
  integer, intent(in) :: glevel
  type(Vector3d), intent(inout), pointer :: pts(:)

  integer :: icgridNum, lcEnd
  type(Vector3d) :: icLv0Vx(12)
  type(Vector3d) :: icPts(2**glevel+1, 2**glevel+1,10)
  integer :: divId, recId, recEnd
  integer :: lcRCId, i, j
  real(DP), parameter :: PI = acos(-1d0)
  integer :: ptId

  icgridNum = 10*4**glevel + 2
  allocate(pts(icgridNum))

  icLv0Vx(:) = icosahedron_vertex()
  
  !
  lcEnd = 2**glevel+1
  icpts(1,1, 1) = icLv0Vx(1)
  icpts(lcEnd,1, 1) = icLv0Vx(2)
  icpts(1,lcEnd, 1) = icLv0Vx(3)
  icpts(lcEnd,lcEnd, 1) = icLv0Vx(7)
  icpts(1,1, 6) = icLv0Vx(3)
  icpts(lcEnd,1, 6) = icLv0Vx(7)
  icpts(1,lcEnd, 6) = icLv0Vx(8)
  icpts(lcEnd,lcEnd, 6) = icLv0Vx(12)

   call split_LCRCSubLCRC(1,lCEnd, 1, lcEnd, icpts(:,:,1))
  call split_LCRCSubLCRC(1,lCEnd, 1, lcEnd, icpts(:,:,6))
  do j=1, lcEnd
    do i=1, lcEnd
      icpts(i,j,1) = normalizedVec(icpts(i,j,1))
      icpts(i,j,6) = normalizedVec(icpts(i,j,6))
    end do
  end do

  do lcRCId=2, 5
    do j=1, lcEnd
      do i=1, lcEnd
        icpts(i,j,lcRCId) = rotateZ(icPts(i,j,1), (lcRCId-1)*2d0*PI/5d0)
        icpts(i,j,5+lcRCId) = rotateZ(icPts(i,j,6), (lcRCId-1)*2d0*PI/5d0)
      end do
    end do
  end do

  pts(1) = icLv0Vx(1)
  pts(icgridNum) = icLv0Vx(12)

  do ptId=2, icgridNum-1
    lcRCId = 1 + (ptId - 2)/(lcENd - 1)**2
    i = 2 + (ptId-2 - (lcRCID-1)*(lcEnd-1)**2)/(lcENd - 1)
    j = ptId - 1 - (((lcRcID - 1)*(lcEnd - 1) + (i - 2))*(lcEnd-1))

    if(.not.(lcRCID<=5 .and. i==1 .and. j==1) .and. &
      & .not.(lcRCID>5 .and. i==lcEnd .and. j==lcEnd) ) then
      pts(ptId) = icpts(i,j,lcRCId)
      !write(*,*) "ptId=, ", ptId, "lcRCId=", lcRCId, ", i=", i, ", j=", j, &
      !  & radToDegUnit(cartToSphPos(pts(ptId)))
      
    end if

  end do

contains
function rotateZ(vec, angle) result(rotVec)
  type(Vector3d), intent(in) :: vec
  real(DP), intent(in) :: angle
  type(Vector3d) :: rotVec

  rotVec%v_(1) = cos(angle)*vec%v_(1) - sin(angle)*vec%v_(2)
  rotVec%v_(2) = sin(angle)*vec%v_(1) + cos(angle)*vec%v_(2)
  rotVec%v_(3) = vec%v_(3)
end function rotateZ

recursive subroutine split_LCRCSubLCRC(i1, i2, j1, j2, icpts)
  integer, intent(in) :: i1, i2, j1, j2
  type(Vector3d), intent(inout) :: icpts(:,:)
  integer :: miId, mjId
  
  if( i1 + 1 == i2 ) return

  miId = (i1 + i2)/2
  mjId = (j1 + j2)/2

  icPts(miId, j1) = 0.5d0*(icPts(i1,j1) + icPts(i2,j1))
  icPts(miId, j2) = 0.5d0*(icPts(i1,j2) + icPts(i2,j2))
  icPts(i1, mjId) = 0.5d0*(icPts(i1,j1) + icPts(i1,j2))
  icPts(i2, mjId) = 0.5d0*(icPts(i2,j1) + icPts(i2,j2))
  icPts(miId, mjId) = 0.5d0*(icPts(i2,j1) + icPts(i1,j2))

  call split_LCRCSubLCRC(i1, miId, j1, mjId, icpts)
  call split_LCRCSubLCRC(i1, miId, mjId, j2, icpts)
  call split_LCRCSubLCRC(miId, i2, j1, mjId, icpts)
  call split_LCRCSubLCRC(miId, i2, mjId, j2, icpts)

end subroutine split_LCRCSubLCRC

end subroutine construct_icosahedralGrid

end module HexTriIcMesh_mod
