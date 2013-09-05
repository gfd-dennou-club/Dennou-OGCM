module VariableSet_mod

  use GeometricField_mod


  use GridSet_mod, only: plMesh
  implicit none


  type(volScalarField) :: v_height, v_div
  type(surfaceScalarField) :: s_normalVel

contains

subroutine VariableSet_Init()

  !
  call GeometricField_Init(v_height, plMesh, "height", "height", "m")
  call GeometricField_Init(v_div, plMesh, "div", "div", "1")
  call GeometricField_Init(s_normalVel, plMesh, "normalVel", "velocity normal to surface", "ms-1")

end subroutine VariableSet_Init

subroutine VariableSet_Final()

  call GeometricField_Final(v_height)
  call GeometricField_Final(v_div)
  call GeometricField_Final(s_normalVel)

end subroutine VariableSet_Final

end module VariableSet_mod
