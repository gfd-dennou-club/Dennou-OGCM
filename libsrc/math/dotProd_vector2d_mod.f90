module dotProd_vector2d_mod
  use dc_types
  use vector2d_mod
  implicit none
  private
  interface operator(.dot.)
     module procedure dotProduct
  end interface operator(.dot.)
  interface normalizedVec
     module procedure normalize_
  end interface normalizedVec
  interface l2norm
    module procedure l2norm_
  end interface l2norm
  public :: operator(.dot.), normalizedVec, l2norm
contains
pure function normalize_(v) result(normVec)
  type(vector2d), intent(in) :: v
  type(vector2d) :: normVec
  normVec = v / sqrt(v .dot. v)
end function normalize_
pure function dotProduct(v1, v2) result(dot)
  type(vector2d), intent(in) :: v1, v2
  real(DP) :: dot
  dot = dot_product(v1%v_, v2%v_)
end function dotProduct
pure function l2norm_(v) result(nrm)
  type(vector2d), intent(in) :: v
  real(DP) :: nrm
  nrm = sqrt( dot_product(v%v_, v%v_ ) )
end function l2norm_
end module dotProd_vector2d_mod
