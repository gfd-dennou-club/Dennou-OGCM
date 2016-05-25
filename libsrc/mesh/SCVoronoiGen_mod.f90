module SCVoronoiGen_mod

  use dc_types
  use dc_message

  use VectorSpace_mod
  use PolyMesh_mod
  use SVoronoiGen2_mod
  use SphericalCoord_mod

  implicit none
  private

  public :: ClosedCurve_Init, ClosedCurve_Final
  public :: SCVoronoiGen_Init, SCVoronoiGen_Final
  public :: SCVoroniDiagram_Generate, SCVoronoi_SetTopology


  integer, parameter :: nShoreLinePts = 10000
  type(vector3d) :: shoreline_pts(nShoreLinePts)
  integer :: shoreLineCellIds(nShoreLinePts)
  integer :: nShoreLineCell

  type, public :: ClosedCurve
     type(vector3d), pointer :: pts(:) => null()
     integer, pointer :: boundaryCellIds(:) => null()
     integer :: nBoundaryCell
  end type ClosedCurve

  character(*), parameter :: module_name = 'SCVoronoiGen_mod'
contains

subroutine ClosedCurve_Init(self, nPts)
  type(ClosedCurve), intent(inout) :: self
  integer, intent(in) :: nPts

  call ClosedCurve_Final(self)

  allocate(self%pts(nPts), self%boundaryCellIds(nPts))

end subroutine ClosedCurve_Init

subroutine ClosedCurve_Final(self)
  type(ClosedCurve), intent(inout) :: self

  if(associated(self%pts)) then
     deallocate(self%pts, self%boundaryCellIds)
  end if
  
end subroutine ClosedCurve_Final

subroutine SCVoronoiGen_Init

end subroutine SCVoronoiGen_Init

subroutine SCVoronoiGen_Final
  call SVoronoi2Gen_Final()
end subroutine SCVoronoiGen_Final

subroutine SCVoroniDiagram_Generate(pts, iniPtsId4, itrMax, density_func, boundarys)
  type(Vector3d), intent(inout) :: pts(:)
  integer, intent(in) :: iniPtsId4(4)
  integer, intent(in) :: itrMax
  type(ClosedCurve), intent(inout), optional :: boundarys(:)
  interface 
     function density_func(x) result(density)
       use dc_types, only:DP
       use VectorSpace_mod, only: Vector3d
       type(vector3d), intent(in) :: x
       real(DP) :: density
     end function density_func
  end interface 

  integer :: itr, i
  real(DP), parameter :: convEPS = 1e-18
  real(DP) :: clusterEnergy, beforeClEnergy

  call SVoronoi2Gen_Init(size(pts))

  do itr=1, itrMax
     write(*,*) "* itr=", itr
     call SVoronoi2DiagramGen(pts, iniPtsId4)

     clusterEnergy = moveSiteToMassCenter(pts, density_func)
     if(present(boundarys)) then
        do i=1, size(boundarys)
           call moveSiteToCoastline(pts, boundarys(i))
        end do
     end if

     if( itr > 1 .and. abs(clusterEnergy - beforeClEnergy)  < convEPS ) exit;

     beforeClEnergy = clusterEnergy
  end do

  call MessageNotify('M', module_name, 'Finish generating a grid..')

end subroutine SCVoroniDiagram_Generate

subroutine SCVoronoi_SetTopology(mesh, boundarys)
  type(PolyMesh), intent(inout) :: mesh
  type(ClosedCurve), intent(in), optional :: boundarys(:)

  type(DomainBoundary) :: domain_boundary
  integer :: i

  call SVoronoi2_SetTopology(mesh)

  if(present(boundarys)) then
     do i=1, size(boundarys)
        call DomainBoundary_Init( domain_boundary, &
             & boundarys(i)%boundaryCellIds(1:boundarys(i)%nBoundaryCell) &
             & )
        call PolyMesh_setDomainBoundary(mesh, i, domain_boundary)
        call DomainBoundary_Final( domain_boundary )
     end do
  end if


end subroutine SCVoronoi_SetTopology

function moveSiteToMassCenter(pts, density_func) result(clusterEnergy)
  type(Vector3d), intent(inout) :: pts(:)
  real(DP) :: clusterEnergy
  interface 
     function density_func(x) result(density)
       use dc_types, only:DP
       use VectorSpace_mod, only: Vector3d
       type(vector3d), intent(in) :: x
       real(DP) :: density
     end function density_func
  end interface 

  integer :: siteId, siteNum
  type(vector3d), allocatable :: lvrVxs(:)
  type(Vector3d) :: tmp, massPt, tri7Pt(7)
  integer :: lvrVxId, lvrVxNum, i
  real(DP) :: posVecs(7,3)
  real(DP) :: triArea, polyArea

type(vector3d) :: newMassPt, tmpMassPt

  siteNum = size(pts)
  clusterEnergy = 0d0

  !$omp parallel do private(lvrVxs, polyArea, tmpMassPt, newMassPt, triArea, lvrVxId, lvrVxNum) reduction(+:clusterEnergy)
  do siteId=1, siteNum
     
     call GetVrRegionPoints(siteId, lvrVxs)
     lvrVxNum = size(lvrVxs)
!!$     massPt = 0d0
     polyArea = 0d0

     newMassPt = 0d0

     do lvrVxId=1, lvrVxNum
        triArea = sphericalTriArea(pts(siteId), lvrVxs(lvrVxid), lvrVxs(mod(lvrVxId,lvrVxNum)+1))
!        massPt = massPt + triArea * ( pts(siteId) + lvrVxs(lvrVxid) + lvrVxs(mod(lvrVxId,lvrVxNum)+1) ) / 3d0
        polyArea = polyArea + triArea


        call getTriMassCenter(tmpMassPt, &
             pts(siteId), lvrVxs(lvrVxid), lvrVxs(mod(lvrVxId,lvrVxNum)+1))
!!$if(siteId==243)then
!!$write(*,*) lvrVxId,triArea
!!$call print(RadToDegUnit(CartToSphPos(pts(siteId))))
!!$call print(RadToDegUnit(CartToSphPos( lvrVxs(lvrVxid))))
!!$call print(RadToDegUnit(CartToSphPos( lvrVxs(mod(lvrVxId,lvrVxNum)+1))))
!!$   call print(tmpMassPt)
!!$end if
        newMassPt = newMassPt + tmpMassPt
     end do

!!$     massPt = normalizedVec( massPt/polyArea ) 
!!$     clusterEnergy = clusterEnergy + l2norm( MassPt - pts(siteId) ) * PolyArea
!!$     pts(siteId) = massPt
!!$


     newMassPt = normalizedVec(newMassPt)
     clusterEnergy = clusterEnergy + l2norm( newMassPt - pts(siteId) ) * PolyArea

     pts(siteId) = newMassPt

if(isNan(clusterEnergy))then
write(*,*) siteId, ":", PolyArea
call print(newMassPt)
end if

!!$
!!$     call print(massPt)
!!$     call print(newMassPt)
  end do

  write(*,*) "Clustering Energy=", clusterEnergy

contains
subroutine getTriMassCenter(massCenter, x0, x1, x2) 
  type(vector3d), intent(in) :: x0, x1, x2
  type(vector3d), intent(out) :: massCenter
  
  type(vector3d) :: x, xp, b_1, b_2

  real(DP), parameter :: y1(6) = &
         & (/ 0.44594849091597, 0.44594849091597, 0.10810301816807, &
         &    0.09157621350977, 0.09157621350977, 0.81684757298046 /)
  real(DP), parameter :: y2(6) = &
         & (/ 0.44594849091597, 0.10810301816807, 0.44594849091597, &
         &    0.09157621350977, 0.81684757298046, 0.09157621350977 /)
  real(DP),  parameter :: sintWt(6) = &
         & (/ 0.22338158967801, 0.22338158967801, 0.22338158967801, &
         &    0.10995174365532, 0.10995174365532, 0.10995174365532 /)

  integer :: i
  real(DP) :: jacob, Wt(6)

  massCenter = 0d0
  do i=1, 6
     xp = get_xp(y1(i), y2(i), x0, x1, x2)
     b_1 = get_b_i(xp, x0, x1); b_2 = get_b_i(xp, x0, x2)
     jacob = sqrt((b_1.dot.b_1)*(b_2.dot.b_2) - (b_1.dot.b_2)**2)

     x = normalizedVec(xp)
     Wt(i) = 0.5d0*sintWt(i)*jacob*density_func(x)
     massCenter = massCenter + Wt(i)*x
  end do

!  massCenter = massCenter/sum(Wt)

end subroutine getTriMassCenter

function get_xp(y1, y2, x0, x1, x2) result(xp)
  real(DP), intent(in) :: y1, y2
  type(vector3d), intent(in) :: x0, x1, x2
  type(vector3d) :: xp
  
  xp = x0 + y1*(x1 - x0) + y2*(x2 - x0)
end function get_xp

function get_b_i(xp, x0, xi) result(b_i)
  type(vector3d), intent(in) :: xp, x0, xi
  type(vector3d) :: b_i
  
  real(DP) :: xp_abs

  xp_abs = l2norm(xp)
  b_i = (xi - x0) - (((xi - x0).dot.xp)/xp_abs**2)*xp
  b_i = b_i/xp_abs
end function get_b_i

end function moveSiteToMassCenter

subroutine moveSiteToCoastline(pts, boundary)
  use SVoronoiGen2_mod, only: &
       & search_VrRegionContainedPt

  type(Vector3d), intent(inout) :: pts(:)
  type(ClosedCurve), intent(inout) :: boundary

  integer :: i, siteId
  integer :: siteNum
  integer :: beforeVrRCId
  integer :: ShoreLinePt2VrRCId(size(boundary%pts))
  logical :: SiteNearBoundaryFlag(size(pts))
  integer :: shoreLinePtIdNearSite(size(pts))
  real(DP) :: dist_ShoreLinePt_Site(size(pts)), dist
  type(vector3d), allocatable :: lvrVxs(:)
  type(vector3d) :: movedNewSite
  integer :: j
  logical :: isSiteMoved
  integer :: boundaryCellId, nBoundaryCell

  siteNum = size(pts)  
  SiteNearBoundaryFlag(:) = .false.
  shoreLinePtIdNearSite(:) = -1
  dist_ShoreLinePt_Site(:) = 1d0
  beforeVrRCId = 1
  nShoreLineCell = 0
  nBoundaryCell = 0

  do i=1, size(boundary%pts)
     siteId = &
          & search_VrRegionContainedPt(boundary%pts(i), pts, beforeVrRCId)

     dist = geodesicArcLength(boundary%pts(i), pts(siteId))
     if( dist_ShoreLinePt_Site(siteId) > dist ) then
        dist_ShoreLinePt_Site(siteId) = dist
        shoreLinePtIdNearSite(siteId) = i
     end if

     if(.not. SiteNearBoundaryFlag(siteId)) then
        SiteNearBoundaryFlag(siteId) = .true.
        nboundaryCell = nboundaryCell + 1
        boundary%boundaryCellIds(nBoundaryCell) = siteId
     end if

     beforeVrRCId = siteId
     ShoreLinePt2VrRCId(i) = siteId
  end do
  
  boundary%nBoundaryCell = nBoundaryCell

!write(*,*) 'nBoundaryCell=', nBoundaryCell, "=>"
  do boundaryCellId=1, nBoundaryCell
     siteId = boundary%boundaryCellIds(boundaryCellId)
     movedNewSite = boundary%pts(shoreLinePtIdNearSite(siteId)) 
     call GetVrRegionPoints(siteId, lvrVxs)
     
     isSiteMoved = .true.
     do j=1, size(lvrVxs)
        if(SphericalTriArea(movedNewSite, lvrVxs(j), lvrVxs(1+mod(j,size(lvrVxs)))) < 1d-12 ) then
           isSiteMoved = .false.
           exit
        end if
     end do
        
     if(isSiteMoved) then
        pts(siteId) = movedNewSite
     else
        boundary%boundaryCellIds(boundaryCellId) = -1
     end if
  end do

  do boundaryCellId=1, nBoundaryCell
     if(boundary%boundaryCellIds(boundaryCellId) == -1) then
        boundary%boundaryCellIds(boundaryCellId:boundary%nBoundaryCell-1) &
             & = boundary%boundaryCellIds(boundaryCellId+1:boundary%nBoundaryCell)
        boundary%nBoundaryCell = boundary%nBoundaryCell - 1
     end if
  end do
!write(*,*) 'nBoundaryCell=', boundary%nBoundaryCell, ":"
  
end subroutine moveSiteToCoastline

end module SCVoronoiGen_mod
