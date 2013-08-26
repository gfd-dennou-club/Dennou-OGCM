module fvCalculus_mod

  use dc_types
  
  use VectorSpace_mod
  use fvMeshInfo_mod
  use PolyMesh_mod

  implicit none
  private

  type(fvMeshInfo), pointer :: fvmInfo => null()
  
  public :: fvCalculus_Init, fvCalculus_Final

contains
subroutine fvCalculus_Init(fvmInfo_)
  type(fvMeshInfo), target, intent(in) :: fvmInfo_

  fvmInfo => fvmInfo_

end subroutine fvCalculus_Init

subroutine fvCalculus_Final()

  fvmInfo => null()

end subroutine fvCalculus_Final

end module fvCalculus_mod
