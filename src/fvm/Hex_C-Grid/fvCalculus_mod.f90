module fvCalculus_mod

  use dc_types
  
  use VectorSpace_mod
  use fvMeshInfo_mod
  use PolyMesh_mod
  use GeometricField_mod

  implicit none
  private

  type(fvMeshInfo), pointer :: fvmInfo => null()
  
  interface div
     module procedure vdiv_volScalar_surfNormVec
     module procedure vdiv_surfNormVec
  end interface div

  interface grad
     module procedure sgrad_volScalar
  end interface grad

  interface curl
     module procedure pcurl_surfNormVec
  end interface curl

  interface interpolate
 !    module procedure interpolate_VolSVal_SurfSVal
  end interface interpolate

  public :: fvCalculus_Init, fvCalculus_Final
  public :: div, grad, curl
  public :: interpolate

contains
subroutine fvCalculus_Init(fvmInfo_)
  type(fvMeshInfo), target, intent(in) :: fvmInfo_

  fvmInfo => fvmInfo_
  
end subroutine fvCalculus_Init

subroutine fvCalculus_Final()

  fvmInfo => null()

end subroutine fvCalculus_Final

function sgrad_volScalar( v_Scalar ) result(s_grad)

  type(volScalarField), intent(in) :: v_scalar
  type(surfaceScalarField) :: s_grad

  integer :: faceNum, faceId, cellIds(2), layerNum
  type(PolyMesh), pointer :: mesh => null()

  mesh => fvmInfo%mesh
  faceNum = getFaceListSize(mesh)
  layerNum = v_Scalar%vLayerNum

  call GeometricField_Init(s_grad, mesh, "s_gradOptrTmpData", vLayerNum=layerNum)
  s_grad%tempDataFlag = .true.

  !$omp parallel do private(cellIds)
  do faceId=1, faceNum
     cellIds(:) = fvmInfo%Face_CellId(:,faceId)
     s_grad%data%v_(1:layerNum,faceId) = ( v_Scalar%data%v_(1:layerNum, CellIds(2)) - v_Scalar%data%v_(1:layerNum, CellIds(1)) ) &
          &                / fvmInfo%s_dualMeshFaceArea%data%v_(1, faceId)  
  end do

  if( v_Scalar%tempDataFlag ) call Release(v_Scalar)
  
end function sgrad_volScalar

function vdiv_surfNormVec( surfNormalVec ) result(v_div)

  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(volScalarField) :: v_div
  
  integer :: cellNum, faceNum, lfaceNum, layerNum
  integer :: cellId, faceId, lfaceId, lyrId, k
  type(PolyMesh), pointer :: mesh
  real(DP), allocatable :: surfaceflux(:,:)
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  layerNum = surfNormalVec%vLayerNum

  allocate(surfaceflux(layerNum, faceNum))

  !$omp parallel do
  do faceId=1, faceNum
     surfaceflux(1:layerNum,faceId) = surfNormalVec%data%v_(1:layerNum,faceId) * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId) )                                 !
  end do


  call GeometricField_Init(v_div, mesh, "v_divOptrTmpData")
  v_div%tempDataFlag = .true.

  !$omp parallel do private(lfaceNum,k)
  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum

     forall(k=1:layerNum) &
          & v_div%data%v_(k,cellId) = &
          & sum( fvmInfo%n_fv(1:lfaceNum,cellId)*surfaceflux(k,fvmInfo%Cell_FaceId(1:lfaceNum,cellId))) &
          & / fvmInfo%v_CellVol%data%v_(k, cellId)
  end do

  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function vdiv_surfNormVec

function vdiv_volScalar_surfNormVec( volScalar, surfNormalVec) result(v_div)
  type(volScalarField), intent(in) :: volScalar
  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(volScalarField) :: v_div
  
  integer :: cellNum, faceNum, lfaceNum, layerNum
  integer :: cellId, faceId, lfaceId, lyrId, k
  type(PolyMesh), pointer :: mesh
  real(DP), allocatable :: surfaceflux(:,:)
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  layerNum = volScalar%vLayerNum

  allocate(surfaceflux(layerNum, faceNum))
  call GeometricField_Init(v_div, mesh, "v_divOptrTmpData")
  v_div%tempDataFlag = .true.

  !$omp parallel do 
  do faceId=1, faceNum
     surfaceflux(1:layerNum,faceId) =   &
          & 0.5d0*sum( volScalar%data%v_(1:layerNum, fvmInfo%Face_CellId(1:2,faceId)), 2)  & ! value interpolated at a face of cell 
          &                * surfNormalVec%data%v_(1:layerNum, faceId)                     & ! 
          &                * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId) )               !
  end do

  !$omp parallel do private(lfaceNum,k)
  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     forall(k=1:layerNum) &
          & v_div%data%v_(k,cellId) = &
          & sum( fvmInfo%n_fv(1:lfaceNum,cellId)*surfaceflux(k,fvmInfo%Cell_FaceId(1:lfaceNum,cellId))) &
          & / fvmInfo%v_CellVol%data%v_(k, cellId)
  end do

  if( volScalar%tempDataFlag ) call Release(volScalar)
  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function vdiv_volScalar_surfNormVec

function pcurl_surfNormVec(surfNormalVec) result(p_curl)
  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(pointScalarField) :: p_curl

  type(PolyMesh), pointer :: mesh
  integer :: ptNum, ptId, faceIds(3), k
  integer :: layerNum
  
  mesh => fvmInfo%mesh
  ptNum = getPointListSize(mesh)
  layerNum = surfNormalVec%vLayerNum

  call GeometricField_Init(p_curl, mesh, "s_curlOptrTmpDat", vLayerNum=layerNum)
  p_curl%TempDataFlag = .true.

  !$omp parallel do private(faceIds, k)
  do ptId=1, ptNum
     faceIds(:) = fvmInfo%Point_FaceId(1:3,ptId)
     forall(k=1:layerNum) &
          & p_curl%data%v_(k,ptId) = &
          & sum( fvmInfo%t_fv(1:3,ptId)*fvmInfo%s_dualMeshFaceArea%data%v_(1,faceIds) * surfNormalVec%data%v_(k,faceIds) ) &
          & / fvmInfo%p_dualMeshCellVol%data%v_(k,ptId)

  end do

  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function pcurl_surfNormVec

!!$subroutine interpolate_VolSVal_SurfSVal(v_scalar, s_scalar)
!!$  type(volScalarField), intent(in) :: v_scalar
!!$  type(surfaceScalarField), intent(out) :: s_scalar
!!$
!!$  integer :: faceId, faceNum
!!$  type(PolyMesh), pointer :: mesh
!!$  type(Face), pointer :: face_
!!$  
!!$  mesh => fvmInfo%mesh
!!$  faceNum = getFaceListSize(mesh)
!!$
!!$  call GeometricField_Init(s_scalar, mesh, "interpSurfVal", "value interpolated at surface")
!!$
!!$  do faceId=1, faceNum
!!$     face_ => mesh%faceList(faceId)
!!$     s_scalar%data%v_(faceId) =   0.5*((v_Scalar.vSlice.face_%ownCellId) + (v_scalar.vSlice.face_%neighCellId))
!!$  end do
!!$
!!$end subroutine interpolate_VolSVal_SurfSVal

end module fvCalculus_mod
