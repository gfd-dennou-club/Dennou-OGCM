module dotProd_vector3d_mod
  use dc_types
  use vector3d_mod
  implicit none
  private
  interface operator(.dot.)
     module procedure dotProduct
  end interface operator(.dot.)
  interface normalizedVec
     module procedure normalize_
  end interface normalizedVec
  public :: operator(.dot.), normalizedVec
contains
function normalize_(v) result(normVec)
  type(vector3d), intent(in) :: v
  type(vector3d) :: normVec
  normVec = v / sqrt(v .dot. v)
end function normalize_
function dotProduct(v1, v2) result(dot)
  type(vector3d), intent(in) :: v1, v2
  real(DP) :: dot
  dot = dot_product(v1%v_, v2%v_)
end function dotProduct
end module dotProd_vector3d_mod
