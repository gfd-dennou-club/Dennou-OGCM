program test_geometricField_mod

  use dc_types
!  use List_d_mod
  use geometricField_mod
!  use volScalarField_mod
  use PolyMesh_mod
use tmp

  implicit none
  
  integer, parameter :: listSize = 100
  type(PolyMesh) :: mesh
  type(volScalarField) :: field1, field2
!  type(List_d) :: list


  call PolyMesh_Init(mesh, listSize, listSize, listSize)

  call GeometricField_Init(field1, mesh, name="sfield1")
  call GeometricField_Init(field2, mesh, name="sfield2")

  field1 = field1 + field2

  call GeometricField_Final(field1)
  call GeometricField_Final(field2)
  call PolyMesh_Final(mesh)

end program test_geometricField_mod
