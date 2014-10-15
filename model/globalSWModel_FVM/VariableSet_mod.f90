module VariableSet_mod

  use GeometricField_mod


  use GridSet_mod, only: plMesh
  implicit none


  type(volScalarField) :: v_h, v_div, v_hb
  type(surfaceScalarField) :: s_normalVel

contains

subroutine VariableSet_Init()

  !
  call GeometricField_Init(v_h, plMesh, "v_h", "fluid layer thickness", "m")
  call GeometricField_Init(v_hB, plMesh, "v_hB", "bottom topography", "m")
  call GeometricField_Init(v_div, plMesh, "v_div", "divergence", "s-1")
  call GeometricField_Init(s_normalVel, plMesh, "s_normalVel", "velocity normal to surface", "ms-1")
  
 v_hb = 0d0

end subroutine VariableSet_Init

subroutine VariableSet_Final()

  call GeometricField_Final(v_h)
  call GeometricField_Final(v_hb)
  call GeometricField_Final(v_div)
  call GeometricField_Final(s_normalVel)

end subroutine VariableSet_Final

end module VariableSet_mod
