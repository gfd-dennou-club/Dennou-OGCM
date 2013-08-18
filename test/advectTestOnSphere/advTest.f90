program advTest

  use dc_types
  use VectorSpace_mod
  use HexTriIcMesh_mod

  implicit none

  type(HexTriIcMesh) :: mesh
  integer, parameter :: glevel = 0

  write(*,*) 'glevel:', glevel
  call HexTriIcMesh_Init(mesh, glevel)
  call HexTriIcMesh_generate(mesh)

contains
end program advTest
