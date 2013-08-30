!
!
module HexTriIcMesh_mod

  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod
  use geometry_mod
  use PolyMesh_mod
  use fvMeshInfo_mod

  use GeometricField_mod

  implicit none
  private

  type, public :: HexTriIcMesh
     integer :: glevel
     type(PolyMesh) :: mesh
     type(Vector3d), pointer :: pts(:) => null()
  end type HexTriIcMesh

  
  public :: HexTriIcMesh_Init, HexTriIcMesh_Final
  public :: HexTriIcMesh_generate, HexTriIcMesh_ConfigFvMeshInfo

contains
subroutine HexTriIcMesh_Init(mesh, glevel_)
  type(HexTriIcMesh), intent(inout) :: mesh
  integer, intent(in) :: glevel_

  integer :: siteNum
  integer :: vertNum

  mesh%glevel = glevel_
  siteNum = 10*4**glevel_ + 2
  vertNum = 2*siteNum - 4
  call PolyMesh_Init(mesh%mesh, vertNum, siteNum+vertNum-2, siteNum)
  
  allocate(mesh%pts(siteNum))

end subroutine HexTriIcMesh_Init

subroutine HexTriIcMesh_Final(mesh)
  type(HexTriIcMesh), intent(inout) :: mesh

  if(associated(mesh%pts)) deallocate(mesh%pts)

  call PolyMesh_Final(mesh%mesh)
  
end subroutine HexTriIcMesh_Final

subroutine HexTriIcMesh_generate(mesh)

  use SCVoronoiGen_mod
!  use SVoronoiGen2_mod

  type(HexTriIcMesh), intent(inout) :: mesh

  integer :: iniPtsId4(4)
  integer :: i

  call construct_icosahedralGrid(mesh%glevel, mesh%pts, iniPtsId4)

  call SCVoronoiGen_Init()
  call SCVoroniDiagram_Generate(mesh%pts, iniPtsId4)
  call SCVoronoi_SetTopology(mesh%mesh)

  call SCVoronoiGen_Final()

end subroutine HexTriIcMesh_generate

subroutine HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)
  type(HexTriIcMesh), intent(in), target :: htiMesh
  type(fvMeshInfo), intent(inout) :: fvInfo

  type(volVectorField) :: v_cellPos
  type(volScalarField) :: v_cellVolume
  type(surfaceVectorField) :: s_surfAreaVec
  
  integer :: cellGId, faceGId, faceLId, faceNum
  type(PolyMesh), pointer :: mesh
  real(DP) :: areas(6), faceArea
  type(Vector3d) :: faceNormal, edgeVxs(2)
  
  mesh => htiMesh%mesh

  !
  call GeometricField_Init( v_cellVolume, mesh, &
      & name="v_CellVol", long_name="volume of Cell", units="m3")

  call GeometricField_Init( v_CellPos, mesh, &
    & name="v_CellVec", long_name="coordinate of cell center", units="m")

  call GeometricField_Init( s_surfAreaVec, mesh, & 
    & name="s_surfAreaVec", long_name="area vector of surface", units="m3")


  !
  v_CellPos = htiMesh%pts  
  
  !
  do cellGId=1, getCellListSize(mesh)
    faceNum = mesh%CellList(cellGId)%faceNum
    do faceLId=1, faceNum
      edgeVxs(:) = getFaceVertex(mesh, cellGId, faceLId)
      areas(faceLId) = sphericalTriArea( v_CellPos.At.(cellGId), edgeVxs(1), edgeVxs(2) ) 
    end do
    
    v_CellVolume%data%v_(cellGId) = sum(areas(1:faceNum))
  end do

  do faceGId=1, getFaceListSize(mesh)
    edgeVxs(:) = getFaceVertex(mesh, faceGId)
    faceArea = geodesicArcLength( edgeVxs(1), edgeVxs(2) )
    faceNormal = normalizedVec( (edgeVxs(2) - edgeVxs(1)) .cross. (v_CellPos.At.(mesh%faceList(faceGId)%ownCellId)) )
    s_surfAreaVec%data%v_(faceGId) =  faceArea * faceNormal
  end do

  !
  call fvMeshInfo_prepair(fvInfo, v_cellPos, v_cellVolume, s_surfAreaVec)

  !
  call GeometricField_Final(v_cellPos)
  call GeometricField_Final(v_cellVolume)
  call GeometricField_Final(s_surfAreaVec)

end subroutine Hextriicmesh_configfvMeshInfo

!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine construct_icosahedralGrid(glevel, pts, iniPtsId4)
  integer, intent(in) :: glevel
  type(Vector3d), intent(inout), pointer :: pts(:)
  integer, intent(out) :: iniPtsId4(4)

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

  iniPtsId4(:) = (/ 1, 2 + (lcEnd-1)*(6*(lcEnd-1) - 1), &
    & 2 + (lcEnd-1)*(7*(lcEnd-1) - 1), 2 + (lcEnd-1)*(9*(lcEnd-1) - 1) /)  
  
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
