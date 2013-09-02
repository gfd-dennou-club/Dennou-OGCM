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
     module procedure div_volScalar_surfNormVec
  end interface div

  interface interpolate
     module procedure interpolate_VolSVal_SurfSVal
  end interface interpolate

  public :: fvCalculus_Init, fvCalculus_Final
  public :: div
  public :: interpolate

  real(DP), allocatable :: n_fv(:,:)

contains
subroutine fvCalculus_Init(fvmInfo_)
  type(fvMeshInfo), target, intent(in) :: fvmInfo_

  integer :: cellNum, faceNum, lfaceNum
  integer :: cellId, faceId, lfaceId 
  type(PolyMesh), pointer :: mesh

  fvmInfo => fvmInfo_
  mesh => fvmInfo%mesh

  cellNum = getCellListSize(mesh)
  faceNum = getFaceListSize(mesh)
  
  allocate( n_fv(6,cellNum) )
  
  n_fv = 0d0
  do cellId=1, cellNum
     lfaceNum = mesh%cellList(cellId)%faceNum
     
     do lfaceId=1, lfaceNum
        faceId = mesh%cellList(cellId)%faceIdList(lfaceId)
        if( mesh%faceList(faceId)%ownCellId==cellId) then
           n_fv(lfaceId, cellId) = 1d0
        else
           n_fv(lfaceId, cellId) = -1d0
        end if
     end do

  end do

end subroutine fvCalculus_Init

subroutine fvCalculus_Final()

  fvmInfo => null()
  if(allocated(n_fv)) deallocate(n_fv)

end subroutine fvCalculus_Final

function div_volScalar_surfNormVec( volScalar, surfNormalVec) result(v_div)
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
     surfaceflux(faceId) =   0.5*((volScalar.At.face_%ownCellId) + (volScalar.At.face_%neighCellId)) & ! value interpolated at a face of cell 
          &                * (surfNormalVec.At.faceId)                                               & ! 
          &                * l2norm( fvmInfo%s_faceAreaVec.At.faceId )                                 !
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

end function div_volScalar_surfNormVec

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
