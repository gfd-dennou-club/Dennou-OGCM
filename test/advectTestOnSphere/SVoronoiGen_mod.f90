module SVoronoiGen_mod

  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod

  implicit none
  private

  public :: SVoronoiGen_Init, SVoronoiGen_Final
  public :: SVoronoiDiagramGen

  integer, allocatable :: vorVx2VorVxId(:,:)
  integer, allocatable :: vorVx2VorRcId(:,:)
  type(Vector3d), allocatable :: vorVx(:)

  integer :: siteNum     ! The number of generators
  integer :: vorVxNum    ! The number of voronoi vetecies

  
contains
subroutine SVoronoiGen_Init(generatorsNum)
  integer, intent(in) :: generatorsNum

  siteNum = generatorsNum
  vorVxNum = 2*siteNum - 4

  if(allocated(vorVx)) deallocate(vorVx)
  if(allocated(vorVx2VorVxId)) deallocate(vorVx2VorVxId)
  if(allocated(vorVx2VorRcId)) deallocate(vorVx2VorRcId)

  allocate(vorVx(vorVxNum))
  allocate(vorVx2VorVxId(3, vorVxNum), vorVx2VorRcId(3, vorVxNum))

end subroutine SVoronoiGen_Init

subroutine SVoronoiDiagramGen(pts, ini4ptsIds_)
  type(Vector3d), intent(in) :: pts(:)
  integer, optional :: ini4ptsIds_(4)

  integer :: nowSiteId  
  integer :: ini4ptsIds(4)
  logical :: isAddedSite(size(pts))

  write(*,*) "= generate voronoi diagram.."
  write(*,*) "  site=", siteNum

  write(*,*) "* prepair.."

  if(present(ini4ptsIds_)) then
     ini4ptsIds(:) = ini4ptsIds_(:)
  else  
     ini4ptsIds(:) = (/ 1, 2, 3, 4 /) 
  end if
  
  isAddedSite = .false.
  isAddedSite(ini4ptsIds) = .true.
  call SVoronoiGen_prepair( pts(ini4ptsIds) )

  do nowSiteId=1, siteNum
     if( .not. isAddedSite(nowSiteId) ) then
        write(*,*) nowSiteId
     end if
  end do

end subroutine SVoronoiDiagramGen

subroutine SVoronoiGen_Final()

  if(allocated(vorVx)) deallocate(vorVx)
  if(allocated(vorVx2VorVxId)) deallocate(vorVx2VorVxId)
  if(allocated(vorVx2VorRcId)) deallocate(vorVx2VorRcId)

end subroutine SVoronoiGen_Final

subroutine SVoronoiGen_prepair(ini4pts)
  
  type(Vector3d), intent(in) :: ini4pts(1:4)

  ! Work variables
  !
  integer :: vxId, i
  integer :: neighVRIds(3) ! The index of three voronoi regions which are a neighbor of a voronoi vertex.  
  real(DP) :: dist_VorVx2unUsedSite

  ! Executable statements

  ! Define the mapping array to convert the index of a voronoi vertex
  ! into idecies of three voronoi regions. 
  vorVx2VorRcId(:,1) = (/ 2,3,4 /)
  vorVx2VorRcId(:,2) = (/ 1,3,4 /)
  vorVx2VorRcId(:,3) = (/ 1,2,4 /)
  vorVx2VorRcId(:,4) = (/ 1,2,3 /)

  ! Calculate the coordinates of four voronoi vetecies.
  do vxId=1, 4
     neighVRIds(:) = vorVx2VorRcId(:,vxId)

     VorVx(vxId) = calcUniSTriCenterPt( &
             &      ini4Pts(neighVRIds(1)), ini4Pts(neighVRIds(2)), ini4Pts(neighVRIds(3)) &
             &   ) 

     ! Check if vorVx2VorRcId is defined correctly.
     ! If the definition of id mapping is not correct, 
     ! the voronoi vertecies will be recomputed after the correction of index. 
     dist_VorVx2unUsedSite = abs( acos( ini4Pts(vxId) .dot. VorVx(vxId) ) )
     do i=1,3
        if( abs( acos( ini4Pts(neighVrIds(i)) .dot. VorVx(vxId) ) ) > dist_VorVx2unUsedSite ) then
           call swap(neighVrIds(2), neighVrIds(3))
           VorVx(vxId) = calcUniSTriCenterPt( &
             &      ini4Pts(neighVRIds(1)), ini4Pts(neighVRIds(2)), ini4Pts(neighVRIds(3)) &
             &   ) 
           
           vorVx2VorRcId(:,vxId) = neighVRIds(:)
           exit     
        end if
     end do

  end do

  ! Initialize the mapping array to convert the index of a voronoi vertex
  ! into the index of three vornoi vertecies which is a neighbor of it. 
  vorVx2VorVxId(:,1:4) = vorVx2VorRcId(:,1:4)

end subroutine SVoronoiGen_prepair

function calcUniSTriCenterPt(pt1, pt2, pt3) result(centerPt)
  type(Vector3d), intent(in) :: pt1, pt2, pt3
  type(Vector3d) :: centerPt

  centerPt = normalizedVec( (pt2 - pt1).cross.(pt3 - pt1) ) 

end function calcUniSTriCenterPt

subroutine swap(i1, i2)
  integer, intent(inout) :: i1, i2
  integer :: tmp

  tmp = i1
  i1 = i2
  i2 = tmp
end subroutine swap

end module SVoronoiGen_mod
