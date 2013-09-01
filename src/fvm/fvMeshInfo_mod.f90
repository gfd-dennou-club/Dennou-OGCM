module fvMeshInfo_mod

  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod
  use List_mod
  use GeometricField_mod
  use PolyMesh_mod
  
  implicit none
  private

  type, public :: fvMeshInfo
    type(volScalarField) :: v_CellVol
    type(surfaceVectorField) :: s_surfAreaVec
    type(PolyMesh), pointer :: mesh => null()
  end type fvMeshInfo

  public :: fvMeshInfo_Init, fvMeshInfo_Final
  public :: fvMeshInfo_prepair

contains
subroutine fvMeshInfo_Init(fvMesh, mesh)
  type(fvMeshInfo), intent(inout) :: fvMesh
  type(PolyMesh), intent(in), target :: mesh

  fvMesh%mesh => mesh
  
  call GeometricField_Init( &
    & fvMesh%v_CellVol, fvMesh%mesh, name="v_CellVol", &
    & long_name="volume of Cell", units="m3")

  call GeometricField_Init( &
    & fvMesh%s_surfAreaVec, fvMesh%mesh, name="s_surfAreaVec", &
    & long_name="area vector of surface", units="m3")

end subroutine fvMeshInfo_Init

subroutine fvMeshInfo_prepair( fvMesh, & 
  & v_cellVolume, s_surfAreaVec &
  & )
  type(fvMeshInfo), intent(inout) :: fvMesh
  type(volScalarField), intent(in) :: v_cellVolume
  type(surfaceVectorField), intent(in) :: s_surfAreaVec

  fvMesh%v_CellVol = v_cellVolume
  fvMesh%s_surfAreaVec = s_surfAreaVec

end subroutine fvMeshInfo_prepair

subroutine fvMeshInfo_Final(fvMesh)
  type(fvMeshInfo), intent(inout) :: fvMesh

  call GeometricField_Final(fvMesh%v_CellVol)
  call GeometricField_Final(fvMesh%s_surfAreaVec)

end subroutine fvMeshInfo_Final

end module fvMeshInfo_mod
