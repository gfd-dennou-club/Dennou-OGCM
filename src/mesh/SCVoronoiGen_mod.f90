module SCVoronoiGen_mod
  use dc_types
  use VectorSpace_mod
  use PolyMesh_mod
  use SVoronoiGen2_mod
  use SphericalCoord_mod

  implicit none
  private

  public :: SCVoronoiGen_Init, SCVoronoiGen_Final
  public :: SCVoroniDiagram_Generate, SCVoronoi_SetTopology


contains

subroutine SCVoronoiGen_Init

end subroutine SCVoronoiGen_Init

subroutine SCVoronoiGen_Final
  call SVoronoi2Gen_Final()
end subroutine SCVoronoiGen_Final

subroutine SCVoroniDiagram_Generate(pts, iniPtsId4, itrMax)
  type(Vector3d), intent(inout) :: pts(:)
  integer, intent(in) :: iniPtsId4(4)
  integer, intent(in) :: itrMax

  integer :: itr

  call SVoronoi2Gen_Init(size(pts))

  do itr=1, itrMax
     call SVoronoi2DiagramGen(pts, iniPtsId4)
     call moveSiteToMassCenter(pts)
  end do

end subroutine SCVoroniDiagram_Generate

subroutine SCVoronoi_SetTopology(mesh)
  type(PolyMesh), intent(inout) :: mesh

  call SVoronoi2_SetTopology(mesh)

end subroutine SCVoronoi_SetTopology

subroutine moveSiteToMassCenter(pts)
  type(Vector3d), intent(inout) :: pts(:)

  integer :: siteId, siteNum
  type(vector3d), allocatable :: lvrVxs(:)
  type(Vector3d) :: tmp, massPt, tri7Pt(7)
  integer :: lvrVxId, lvrVxNum, i
  real(DP) :: posVecs(7,3)
  real(DP) :: clusterEnergy, area

  siteNum = size(pts)
  clusterEnergy = 0d0

  do siteId=1, siteNum
     
     call GetVrRegionPoints(siteId, lvrVxs)
     lvrVxNum = size(lvrVxs)
     massPt = 0d0
     area = 0d0

!!$     tri7Pt(1) = pts(siteId)
!!$
     do lvrVxId=1, lvrVxNum
!!$        tri7Pt(2) = lvrVxs(lvrVxId)
!!$        tri7Pt(3) = lvrVxs(mod(lvrVxId,3)+1)
!!$        tri7Pt(4) = normalizedVec( 0.5d0*(tri7Pt(1)+tri7Pt(2)) )
!!$        tri7Pt(5) = normalizedVec( 0.5d0*(tri7Pt(2)+tri7Pt(3)) )
!!$        tri7Pt(6) = normalizedVec( 0.5d0*(tri7Pt(3)+tri7Pt(1)) )
!!$        tri7Pt(7) = normalizedVec( (tri7Pt(1)+tri7Pt(2)+tri7Pt(3))/3d0 )
!!$
!!$        do i=1, 7
!!$           posVecs(i,1:3) = toArray(tri7Pt(i))
!!$        end do
!!$
!!$        tmp = (/ volIntegrate(posVecs(:,1)), volIntegrate(posVecs(:,2)), volIntegrate(posVecs(:,3)) /)
        massPt = massPt + lvrVxs(lvrVxId)
        area = area + sphericalTriArea(pts(siteId), lvrVxs(lvrVxid), lvrVxs(mod(lvrVxId,lvrVxNum)+1))
     end do

     massPt = normalizedVec( massPt/dble(lvrVxNum) ) 
     clusterEnergy = clusterEnergy + l2norm( massPt - pts(siteId)) * area
     pts(siteId) = massPt
  end do

  write(*,*) "Clustering Energy=", clusterEnergy

end subroutine moveSiteToMassCenter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!

function volIntegrate(f) result(integ)
  real(DP), intent(in) :: f(7) 

  real(DP) :: integ

  integ = (  0.05d0*sum(f(1:3)) + 2d0/15d0*sum(f(4:6)) + 9d0/20d0*f(7) ) 
  
end function volIntegrate

function lIntrp(xi, eta, f7) result(val)
  real(DP), intent(in) :: xi ,eta
  real(DP), intent(in) :: f7(7)
  real(DP) :: val

  real(DP) :: c(7)
  
  c(1) = -3d0*xi*xi*eta -3d0*xi*eta*eta + 2d0*xi*xi + 2d0*eta*eta +7d0*xi*eta -3d0*xi -3d0*eta + 1d0
  c(2) = -3d0*xi*xi*eta -3d0*xi*eta*eta + 2d0*xi*xi +3d0*xi*eta -xi
  c(3) = -3d0*xi*xi*eta -3d0*xi*eta*eta + 2d0*eta*eta +3d0*xi*eta -eta
  c(4) = 12d0*xi*xi*eta + 12d0*xi*eta*eta -4d0*xi*xi -16d0*xi*eta + 4d0*xi
  c(5) = 12d0*xi*xi*eta + 12d0*xi*eta*eta -8d0*xi*eta
  c(6) = 12d0*xi*xi*eta + 12d0*xi*eta*eta -4d0*eta*eta -16d0*eta*xi +4d0*eta
  c(7) = -27d0*xi*xi*eta -27d0*xi*eta*eta + 27d0*xi*eta

  val = dot_product(c(:), f7(:))
  
end function lIntrp

function lIntrp_dy(xi, eta, f7) result(val)
  real(DP), intent(in) :: xi ,eta
  real(DP), intent(in) :: f7(7)
  real(DP) :: val

  real(DP) :: c(7)

  
  c(1) = -3d0*xi*xi -6d0*xi*eta + 4d0*eta +7d0*xi -3d0
  c(2) = -3d0*xi*xi -6d0*xi*eta + 3d0*xi
  c(3) = -3d0*xi*xi -6d0*xi*eta + 4d0*eta +3d0*xi -1d0
  c(4) = 12d0*xi*xi + 24d0*xi*eta -16d0*xi
  c(5) = 12d0*xi*xi + 24d0*xi*eta -8d0*xi
  c(6) = 12d0*xi*xi + 24d0*xi*eta -8d0*eta -16d0*xi +4d0
  c(7) = -27d0*xi*xi -54d0*xi*eta + 27d0*xi

  val = dot_product(c(:), f7(:))
  
end function lIntrp_dy

function lIntrp_dx(xi, eta, f7) result(val)
  real(DP), intent(in) :: xi ,eta
  real(DP), intent(in) :: f7(7)
  real(DP) :: val

  real(DP) :: c(7)
  
  c(1) = -6d0*xi*eta - 3d0*eta*eta + 4d0*xi +7d0*eta -3d0
  c(2) = -6d0*xi*eta -3d0*eta*eta + 4d0*xi +3d0*eta -1d0
  c(3) = -6d0*xi*eta -3d0*eta*eta + 3d0*eta 
  c(4) = 24d0*xi*eta + 12d0*eta*eta -8d0*xi -16d0*eta + 4d0
  c(5) = 24d0*xi*eta + 12d0*eta*eta -8d0*eta
  c(6) = 24d0*xi*eta + 12d0*eta*eta -16d0*eta
  c(7) = -54d0*xi*eta -27d0*eta*eta + 27d0*eta

  val = dot_product(c(:), f7(:))
  
end function lIntrp_dx

end module SCVoronoiGen_mod
