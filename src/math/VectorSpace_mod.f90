module VectorSpace_mod
  use vector2d_mod
  use vector3d_mod

  use dotProd_vector2d_mod
  use dotProd_vector3d_mod

  use crossProd_mod

  implicit none 
  private

  ! Cascade
  public :: vector2d, vector3d
  public :: operator(+), operator(-), operator(*), operator(/), assignment(=)
  public :: toArray
  public :: operator(.dot.), normalizedVec, l2norm
  public :: operator(.cross.)

  public :: print

contains

end module VectorSpace_mod
