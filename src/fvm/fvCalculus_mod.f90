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
     module procedure interpolate_VolSVal_SurfSVal
  end interface interpolate

  public :: fvCalculus_Init, fvCalculus_Final
  public :: div, grad, curl
  public :: interpolate

  integer, allocatable :: n_fv(:,:)
  integer, allocatable :: t_fv(:,:)

  public :: n_fv, t_fv

contains
subroutine fvCalculus_Init(fvmInfo_)
  type(fvMeshInfo), target, intent(in) :: fvmInfo_

  integer :: cellNum, pointNum, faceNum, lfaceNum
  integer :: cellId, faceId, lfaceId, ptId 
  type(PolyMesh), pointer :: mesh

  fvmInfo => fvmInfo_
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  pointNum = getPointListSize(mesh)

  allocate( n_fv(6,cellNum) )
  allocate( t_fv(6, pointNum) )

  n_fv = 0
  t_fv = 0
  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     
     do lfaceId=1, lfaceNum
        faceId = mesh%cellList(cellId)%faceIdList(lfaceId)
        if( mesh%faceList(faceId)%ownCellId==cellId) then
           n_fv(lfaceId, cellId) = 1
        else
           n_fv(lfaceId, cellId) = -1
        end if
     end do

  end do

  do ptId=1, pointNum
     do lfaceId=1, 3!size(fvInfo%Point_FaceId,1)
        faceID = fvmInfo%Point_FaceId(lfaceId, ptId)
        if(fvmInfo%Face_PointId(2,faceId)==ptId) then
           t_fv(lfaceId, ptId) = 1
        else
           t_fv(lfaceId, ptId) = -1
        end if
     end do
  end do

end subroutine fvCalculus_Init

subroutine fvCalculus_Final()

  fvmInfo => null()
  if(allocated(n_fv)) deallocate(n_fv)

end subroutine fvCalculus_Final

function sgrad_volScalar( v_Scalar ) result(s_grad)
  use SphericalCoord_mod

  type(volScalarField), intent(in) :: v_scalar
  type(surfaceScalarField) :: s_grad

  integer :: faceNum, faceId, cellIds(2)
  type(PolyMesh), pointer :: mesh
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh

  faceNum = getFaceListSize(mesh)

  call GeometricField_Init(s_grad, mesh, "s_gradOptrTmpData")
  s_grad%tempDataFlag = .true.

  do faceId=1, faceNum
     face_ => mesh%faceList(faceId)
     cellIds(:) = fvmInfo%Face_CellId(:,faceId)
     s_grad%data%v_(faceId) = ( (v_Scalar.At.CellIds(2)) - (v_Scalar.At.CellIds(1)) ) &
          &                / (fvmInfo%s_dualMeshFaceArea.At.faceId)  
  end do

  if( v_Scalar%tempDataFlag ) call Release(v_Scalar)
  
end function sgrad_volScalar

function vdiv_surfNormVec( surfNormalVec) result(v_div)

  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(volScalarField) :: v_div
  
  integer :: cellNum, faceNum, lfaceNum
  integer :: cellId, faceId, lfaceId 
  type(PolyMesh), pointer :: mesh
  real(DP), allocatable :: surfaceflux(:)
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  allocate(surfaceflux(faceNum))

  do faceId=1, faceNum
     face_ => mesh%faceList(faceId)
     surfaceflux(faceId) =  (surfNormalVec.At.faceId) * l2norm( fvmInfo%s_faceAreaVec.At.faceId )                                 !
  end do


  call GeometricField_Init(v_div, mesh, "v_divOptrTmpData")
  v_div%tempDataFlag = .true.

  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     v_div%data%v_(cellId) = &
          &   dot_product( dble(n_fv(1:lfaceNum,cellId)), surfaceflux(fvmInfo%Cell_FaceId(1:lfaceNum,cellId)) ) &
          & / (fvmInfo%v_CellVol.At.cellId)
  end do

  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function vdiv_surfNormVec

function vdiv_volScalar_surfNormVec( volScalar, surfNormalVec) result(v_div)
  type(volScalarField), intent(in) :: volScalar
  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(volScalarField) :: v_div
  
  integer :: cellNum, faceNum, lfaceNum
  integer :: cellId, faceId, lfaceId 
  type(PolyMesh), pointer :: mesh
  real(DP), allocatable :: surfaceflux(:)
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  allocate(surfaceflux(faceNum))

  do faceId=1, faceNum
     face_ => mesh%faceList(faceId)
     surfaceflux(faceId) =   0.5d0*((volScalar.At.face_%ownCellId) + (volScalar.At.face_%neighCellId)) & ! value interpolated at a face of cell 
          &                * (surfNormalVec.At.faceId)                                                 & ! 
          &                * l2norm( fvmInfo%s_faceAreaVec.At.faceId )                                 !

!     write(*,*) faceID,":", geodesicArcLength(

  end do


  call GeometricField_Init(v_div, mesh, "v_divOptrTmpData")
  v_div%tempDataFlag = .true.

  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     v_div%data%v_(cellId) = &
          &   dot_product( n_fv(1:lfaceNum,cellId), surfaceflux(mesh%cellList(cellId)%faceIdList(1:lfaceNum)) ) &
          & / (fvmInfo%v_CellVol.At.cellId)
    
  end do


  if( volScalar%tempDataFlag ) call Release(volScalar)
  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function vdiv_volScalar_surfNormVec

function pcurl_surfNormVec(surfNormalVec) result(p_curl)
  type(surfaceScalarField), intent(in) :: surfNormalVec
  type(pointScalarField) :: p_curl

  type(PolyMesh), pointer :: mesh
  integer :: ptNum, ptId, faceIds(3)
  
  mesh => fvmInfo%mesh
  ptNum = getPointListSize(mesh)

  call GeometricField_Init(p_curl, mesh, "s_curlOptrTmpDat")
  p_curl%TempDataFlag = .true.

  do ptId=1, ptNum
     faceIds(:) = fvmInfo%Point_FaceId(1:3,ptId)
     p_curl%data%v_(ptId) = &
          & sum( t_fv(1:3,ptId) * (fvmInfo%s_dualMeshFaceArea.At.faceIds) * (surfNormalVec.At.faceIds) ) &
          & / (fvmInfo%p_dualMeshCellVol.At.ptId)
  end do

  if( surfNormalVec%tempDataFlag ) call Release(surfNormalVec)

end function pcurl_surfNormVec

subroutine interpolate_VolSVal_SurfSVal(v_scalar, s_scalar)
  type(volScalarField), intent(in) :: v_scalar
  type(surfaceScalarField), intent(out) :: s_scalar

  integer :: faceId, faceNum
  type(PolyMesh), pointer :: mesh
  type(Face), pointer :: face_
  
  mesh => fvmInfo%mesh
  faceNum = getFaceListSize(mesh)

  call GeometricField_Init(s_scalar, mesh, "interpSurfVal", "value interpolated at surface")

  do faceId=1, faceNum
     face_ => mesh%faceList(faceId)
     s_scalar%data%v_(faceId) =   0.5*((v_Scalar.At.face_%ownCellId) + (v_scalar.At.face_%neighCellId))
  end do

end subroutine interpolate_VolSVal_SurfSVal

end module fvCalculus_mod
