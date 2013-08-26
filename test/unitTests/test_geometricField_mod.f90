program test_geometricField_mod

  use dc_types
  use geometricField_mod
  use PolyMesh_mod

  implicit none
  
  integer, parameter :: listSize = 100
  type(PolyMesh) :: mesh
  type(volScalarField) :: field1, field2, field3


  call PolyMesh_Init(mesh, listSize, listSize, listSize)

  call GeometricField_Init(field1, mesh, name="sfield1")
  call GeometricField_Init(field2, mesh, name="sfield2")
  call GeometricField_Init(field3, mesh, name="sfield3")

  field1 = 1d0
  field2 = 0d0

  field3 = field1 + field2 + field1
  write(*,*) field3%data%v_(1:10)

  field3 = field1 + (field2 - field1) 
  write(*,*) field3%data%v_(1:10)

  call GeometricField_Final(field1)
  call GeometricField_Final(field2)
  call GeometricField_Final(field3)
  call PolyMesh_Final(mesh)

end program test_geometricField_mod
