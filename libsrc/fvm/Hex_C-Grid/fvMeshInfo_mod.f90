module fvMeshInfo_mod

  use dc_types
  use dc_message
  use VectorSpace_mod
  use SphericalCoord_mod
  use List_mod
  use GeometricField_mod
  use PolyMesh_mod
  
  implicit none
  private

  type, public :: fvMeshInfo
     type(PolyMesh), pointer :: mesh => null()

     type(volScalarField) :: v_CellVol
     type(surfaceVectorField) :: s_faceAreaVec
     type(surfaceVectorField) :: s_faceCenter
     type(surfaceScalarField) :: s_dualMeshFaceArea
     type(PointScalarField) :: p_dualMeshCellVol

     integer, pointer :: Cell_FaceId(:,:) => null()
     integer, pointer :: Point_FaceId(:,:) => null()
     integer, pointer :: Face_CellId(:,:) => null()
     integer, pointer :: Point_CellId(:,:) => null()
     integer, pointer :: Cell_PointId(:,:) => null()
     integer, pointer :: Face_PointId(:,:) => null()
     integer, pointer :: Face_PairCellFaceId(:,:) => null()
     integer, pointer :: CellPoint_PairFaceId(:,:,:) => null()

     integer, pointer :: n_fv(:,:)
     integer, pointer :: t_fv(:,:)
    
     logical :: dualMeshFlag = .false.

  end type fvMeshInfo

  public :: fvMeshInfo_Init, fvMeshInfo_Final
  public :: fvMeshInfo_prepair

contains
subroutine fvMeshInfo_Init(fvMesh, mesh, dualMeshFlag)
  type(fvMeshInfo), intent(inout) :: fvMesh
  type(PolyMesh), intent(in), target :: mesh
  logical, optional, intent(in) :: dualMeshFlag

  fvMesh%mesh => mesh
  if(present(dualMeshFlag)) fvMesh%dualMeshFlag = dualMeshFlag

  call GeometricField_Init( &
    & fvMesh%v_CellVol, fvMesh%mesh, name="v_CellVol", &
    & long_name="volume of Cell", units="m3", vLayerNum=1)

  call GeometricField_Init( &
    & fvMesh%s_faceAreaVec, fvMesh%mesh, name="s_faceAreaVec", &
    & long_name="area vector of face", units="m2", vLayerNum=1)

  call GeometricField_Init( &
    & fvMesh%s_faceCenter, fvMesh%mesh, name="s_faceCenter", &
    & long_name="center of face", units="m", vLayerNum=1)

  if( fvMesh%dualMeshFlag ) then
     call GeometricField_Init( &
          & fvMesh%s_dualMeshFaceArea, fvMesh%mesh, name="s_dualMeshFaceArea", &
          & long_name="area of face in dual mesh", units="m2", vLayerNum=1)

     call GeometricField_Init( &
          & fvMesh%p_dualMeshCellVol, fvMesh%mesh, name="p_dualMeshCellVol", &
          & long_name="area of face in dual mesh", units="m3", vLayerNum=1)

  end if

  call set_ConectivityInfo(fvMesh)
  call set_edgeNormalDirSign(fvMesh)

end subroutine fvMeshInfo_Init

subroutine fvMeshInfo_prepair( fvMesh, & 
  & v_cellVolume, s_faceAreaVec, s_faceCenter, &
  & s_dualMeshFaceArea, p_dualMeshCellVol &
  & )
  type(fvMeshInfo), intent(inout) :: fvMesh
  type(volScalarField), intent(in) :: v_cellVolume
  type(surfaceVectorField), intent(in) :: s_faceAreaVec
  type(surfaceVectorField), intent(in) :: s_faceCenter
  type(surfaceScalarField), intent(in), optional :: s_dualMeshFaceArea
  type(PointScalarField), intent(in), optional :: p_dualMeshCellVol

  fvMesh%v_CellVol = v_cellVolume
  fvMesh%s_faceAreaVec = s_faceAreaVec
  fvMesh%s_faceCenter = s_faceCenter

  if ( fvMesh%dualMeshFlag ) then
     if( present(s_dualMeshFaceArea) .and. present(p_dualMeshCellVol) ) then
        fvMesh%s_dualMeshFaceArea = s_dualMeshFaceArea
        fvMesh%p_dualMeshCellVol = p_dualMeshCellVol
     else
        call MessageNotify('E', 'DQGModel::fvMeshInfo_mod::fvMeshInfo_prepair', &
             & "Both p_dualMeshCellVol and s_dualMeshFaceArea must be specified if the flag of dual mesh is true.")
        stop
     end if
  end if
end subroutine fvMeshInfo_prepair

subroutine fvMeshInfo_Final(fvMesh)
  type(fvMeshInfo), intent(inout) :: fvMesh

  call GeometricField_Final(fvMesh%v_CellVol)
  call GeometricField_Final(fvMesh%s_faceAreaVec)
  call GeometricField_Final(fvMesh%s_faceCenter)

  if( fvMesh%dualMeshFlag ) then
     call GeometricField_Final(fvMesh%s_dualMeshFaceArea)
     call GeometricField_Final(fvMesh%p_dualMeshCellVol)
  end if

  deallocate(fvMesh%Cell_FaceId)
  deallocate(fvMesh%Point_FaceId)
  deallocate(fvMesh%Face_CellId)
  deallocate(fvMesh%Point_CellId)
  deallocate(fvMesh%Face_PointId)
  deallocate(fvMesh%Cell_PointId)
  deallocate(fvMesh%Face_PairCellFaceId)
  deallocate(fvMesh%CellPoint_PairFaceId)

  deallocate(fvMesh%n_fv, fvMesh%t_fv)
end subroutine fvMeshInfo_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine set_ConectivityInfo(fvm)
  type(fvMeshInfo), intent(inout), target :: fvm

  type(PolyMesh), pointer :: mesh => null()
  integer :: cellNum, faceNum, pointNum
  integer :: lfaceNum
  integer :: i, j, k

  integer :: faceId, ptId, nowCellId, prevFaceId, nowFaceId
  logical, allocatable :: ptFlag(:)
  type(Face), pointer :: nowFace
  integer :: walkCounter, pairCellId(2), lfaceNum1, lfaceNum2

  mesh => fvm%mesh
  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  pointNum = getPointListSize(mesh)
  allocate(ptFlag(pointNum))

  allocate( fvm%Cell_FaceId(MAX_CELL_FACE_NUM, cellNum) )
  allocate( fvm%Point_FaceId(3, pointNum) )
  allocate( fvm%Face_CellId(2, faceNum) )
  allocate( fvm%Point_CellId(3, pointNum) )
  allocate( fvm%Face_PointId(2, faceNum) )
  allocate( fvm%Cell_PointId(MAX_CELL_FACE_NUM, cellNum) )
  allocate( fvm%Face_PairCellFaceId(2*MAX_CELL_FACE_NUM, faceNum) )
  allocate( fvm%CellPoint_PairFaceId(2, MAX_CELL_FACE_NUM, cellNum) )

  fvm%Cell_FaceId = -1
  do i=1, cellNum
     lfaceNum = mesh%cellList(i)%faceNum
     fvm%Cell_FaceId(1:lfaceNum, i) = mesh%cellList(i)%faceIdList(1:lfaceNum)
     do j=1, lfaceNum
        nowFace => mesh%faceList(fvm%Cell_FaceId(j,i))
        if( nowFace%ownCellId==i) then
           fvm%Cell_PointId(j, i) = nowFace%vertIdList(1)
        else
           fvm%Cell_PointId(j, i) = nowFace%vertIdList(2)
        end if
     end do
  end do

  fvm%Face_CellId = -1
  fvm%Face_PointId = -1
  do i=1, faceNum
     fvm%Face_CellId(1:2,i) = (/ mesh%faceList(i)%ownCellId, mesh%faceList(i)%neighCellId /)
     fvm%Face_PointId(1:2,i) = mesh%faceList(i)%vertIdList(1:2)
  end do

  !
  ptFlag = .false.
  do i=1, CellNum
     lfaceNum = mesh%cellList(i)%faceNum
     do j=1, lfaceNum
        faceId = fvm%Cell_FaceId(j,i)
        ptId = mesh%faceList(faceId)%vertIdList(2)

        if( .not. ptFlag(ptId) ) then
           walkCounter = 1
           nowCellId = i
           nowFaceId = faceId

           ! Walk around a vertex to Collect a information about the cells and face around it.
           do while(.true.)
              fvm%Point_CellId(walkCounter, ptId) = nowCellId
              fvm%Point_FaceId(walkCounter, ptId) = nowfaceId
              
!              write(*,*) "cell pair", fvm%Face_CellId(1:2, nowFaceId), "(nowcell=", nowCellId, ")"
              if( fvm%Face_CellId(1, nowFaceId) == nowCellId ) then
                 nowCellId = fvm%Face_CellId(2, nowFaceId)
              else
                 nowCellId = fvm%Face_CellId(1, nowFaceId)
              end if

              if( nowCellId == i) exit; ! Finish exproring the cells and face around a vertex(ptId)

              prevFaceId = nowFaceId
              do k=1, mesh%CellList(nowCellId)%faceNum
                 nowFaceId = fvm%Cell_FaceId(k,nowCellId)
                 nowFace => mesh%FaceList(nowFaceID)
                 if(  nowFaceId /= prevFaceId  .and. & 
                      & (nowFace%vertIdList(1) == ptId  .or. nowFace%vertIdList(2) == ptId) ) exit;
              end do

              walkCounter = walkCounter + 1
           end do

           ptFlag(ptId) = .true.
        end if

     end do
  end do


  !
  fvm%Face_PairCellFaceId = -1
  do i=1, faceNum
     pairCellId(:) = fvm%Face_CellId(:,i)
     lfaceNum1 = mesh%CellList(pairCellId(1))%faceNum
     lfaceNum2 = mesh%CellList(pairCellId(2))%faceNum

     do j=1, lfaceNum1
        if( i == fvm%Cell_FaceId(j,pairCellId(1)) ) then
           fvm%Face_PairCellFaceId(1:lfaceNum1, i) =  cshift(fvm%Cell_FaceId(1:lfaceNum1, pairCellId(1)), j-1) 
           exit
        end if
     end do

     do j=1, lfaceNum2
        if( i == fvm%Cell_FaceId(j,pairCellId(2)) ) then
           fvm%Face_PairCellFaceId(lfaceNum1+1:lfaceNum1+lfaceNum2, i) =  cshift(fvm%Cell_FaceId(1:lfaceNum2, pairCellId(2)), j) 
           exit
        end if
     end do
     fvm%Face_PairCellFaceId(lfaceNum1+lfaceNum2, i) = -1

!write(*,*) "faceID=", faceId, fvm%Face_PairCellFaceId(:,i)
  end do

  fvm%CellPoint_PairFaceId = -1
  do i=1, cellNum
     lfaceNum = mesh%CellList(i)%faceNum
     do j=1, lfaceNum
        faceId = fvm%Cell_FaceId(j,i)
        if( mesh%FaceList(faceId)%ownCellId == i) then
           ptId = fvm%Face_PointId(2,faceId) 
        else
           ptId = fvm%Face_PointId(1,faceId) 
        end if

        do k=1, lfaceNum
           if( fvm%Cell_PointId(k,i) == ptId) then
              fvm%CellPoint_PairFaceId(1:2, k, i) = (/ faceId, fvm%Cell_FaceId(mod(j,lfaceNum)+1, i) /)
              exit
           end if
           if( lfaceNum==k ) then
              write(*,*) "error", ptId, fvm%Cell_PointId(:,i)
              write(*,*) "**", fvm%Face_PointId(:,faceId), "faceNum=", lfaceNum
              stop
           end if
        end do

     end do
  end do

end subroutine set_ConectivityInfo

subroutine set_edgeNormalDirSign(fvm)
  type(fvMeshInfo), intent(inout), target :: fvm

  type(PolyMesh), pointer :: mesh
  integer :: cellNum, faceNum, lfaceNum, pointNum
  integer :: cellId, ptId, faceId, lfaceId

  mesh => fvm%mesh
  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  pointNum = getPointListSize(mesh)
  
  allocate( fvm%n_fv(MAX_CELL_FACE_NUM, cellNum) )
  allocate( fvm%t_fv(MAX_CELL_FACE_NUM, pointNum) )

  fvm%n_fv = 0
  fvm%t_fv = 0
  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     
     do lfaceId=1, lfaceNum
        faceId = mesh%cellList(cellId)%faceIdList(lfaceId)
        if( mesh%faceList(faceId)%ownCellId==cellId) then
           fvm%n_fv(lfaceId, cellId) = 1
        else
           fvm%n_fv(lfaceId, cellId) = -1
        end if
     end do

  end do

  do ptId=1, pointNum
     do lfaceId=1, 3!size(fvInfo%Point_FaceId,1)
        faceID = fvm%Point_FaceId(lfaceId, ptId)
        if(fvm%Face_PointId(2,faceId)==ptId) then
           fvm%t_fv(lfaceId, ptId) = 1
        else
           fvm%t_fv(lfaceId, ptId) = -1
        end if
     end do
  end do

end subroutine set_edgeNormalDirSign

end module fvMeshInfo_mod
