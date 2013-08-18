#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _moduleNameSuffix1(s) s ## _mod
#define _moduleNameSuffix2(s1, s2) s1 ## _ ## s2 ## _ ## mod
#define moduleName _eval2(_moduleNameSuffix2, dotProd, vectorspaceTypeName )
#define useModuleName _eval1(_moduleNameSuffix1, vectorspaceTypeName )

module moduleName
  use dc_types
  use useModuleName
  
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
  type(vectorspaceTypeName), intent(in) :: v
  type(vectorspaceTypeName) :: normVec

  normVec = v / sqrt(v .dot. v)

end function normalize_

function dotProduct(v1, v2) result(dot)
  type(vectorspaceTypeName), intent(in) :: v1, v2
  real(DP) :: dot

  dot = dot_product(v1%v_, v2%v_)

end function dotProduct


end module moduleName
