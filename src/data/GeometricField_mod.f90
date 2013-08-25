module GeometricField_mod

  use pointScalarField_mod
  use pointVectorField_mod
  use surfaceScalarField_mod
  use surfaceVectorField_mod
  use volScalarField_mod
  use volVectorField_mod

  implicit none
  private

  ! Cascade
  public :: pointScalarField, pointVectorField
  public :: surfaceScalarField, surfaceVectorField
  public :: volScalarField, volVectorField

  public :: GeometricField_Init, GeometricField_Final
  public :: SetFieldAtitude
  public :: operator(+), assignment(=)

contains
end module GeometricField_mod
