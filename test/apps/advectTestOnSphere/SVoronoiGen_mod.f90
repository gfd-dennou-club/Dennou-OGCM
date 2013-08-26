module SVoronoiGen_mod

  use dc_types
  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod

  implicit none
  private

  public :: SVoronoiGen_Init, SVoronoiGen_Final
  public :: SVoronoiDiagramGen

  public :: SVoronoi_setTopology

  integer, allocatable :: vorVx2VorVxId(:,:)
  integer, allocatable :: vorVx2VorRcId(:,:)
  real(DP), allocatable :: vorVxRadius(:)
  type(Vector3d), allocatable :: vorVx(:)
  type(Face), allocatable :: vorRCList(:)

  integer :: siteNum     ! The number of generators
  integer :: vorVxNum    ! The number of voronoi vetecies
  integer :: vorIdMapHead
  
contains
subroutine SVoronoiGen_Init(generatorsNum)
  integer, intent(in) :: generatorsNum

  type(Vector3d) :: zeroVec

  siteNum = generatorsNum
  vorVxNum = 2*siteNum - 4 

  if(allocated(vorVx)) deallocate(vorVx)
  if(allocated(vorVx2VorVxId)) deallocate(vorVx2VorVxId)
  if(allocated(vorVx2VorRcId)) deallocate(vorVx2VorRcId)
  if(allocated(vorVxRadius)) deallocate(vorVxRadius)
  if(allocated(vorRCList)) deallocate(vorRCList)

  allocate(vorVx(vorVxNum), vorVxRadius(vorVxNum))
  allocate(vorVx2VorVxId(3, vorVxNum), vorVx2VorRcId(3, vorVxNum))
  allocate(vorRCList(siteNum))

  vorVxRadius(:) = 0d0
  zeroVec = 0d0
  vorVx(:) = zeroVec
  vorIdMapHead = 0

end subroutine SVoronoiGen_Init

subroutine SVoronoiDiagramGen(pts, ini4ptsIds_)
  type(Vector3d), intent(in) :: pts(:)
  integer, optional :: ini4ptsIds_(4)

  integer :: nowSiteId, i, j, bVxId  
  integer :: ini4ptsIds(4)
  logical :: isAddedSite(size(pts))
  integer, allocatable :: brokenVxIdList(:)
  logical :: bVxFlag(vorVxNum)

  integer :: vorRgVx(6, siteNum)
  integer :: vxId, nextVxId, endVxId, cellVxNum
  

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
  call SVoronoiGen_prepair( pts, ini4ptsIds )

  !
  write(*,*) "* create voronoi mesh using increment method.."

  do nowSiteId=1, siteNum
    if( mod(nowSiteId, siteNum/10) == 0) then
      write(*,*) "..", nowSiteId*100/siteNum, "%"
    end if

    if( .not. isAddedSite(nowSiteId) ) then
 !       write(*,'(a,i4)') "----", nowSiteId

        ! Step1, 2: collect an information with the broken voronoi vertecies by adding new site.
        call getBrokenVorVxIdList(pts(nowSiteId), brokenVxIdList, bVxFlag)
 !       write(*,*) "broken Vx:", brokenVxIdList(:)

        ! Step3, 4: generate new voronoi region corresponding to new site
        !           using an information of broken voronoi vertecies.
        call construct_newVoronoiRegion(nowSiteId, brokenVxIdList, bVxFlag, pts)

     end if
  end do


!do i=1,vorVxNum
!   write(*,'(i2,a,3i7,a,3i7, a, 1f12.5, 5f12.5)') i, "*", vorVx2vorRcId(:,i), ":", vorVx2vorVxId(:,i), &
!        & ":",  VorVxRadius(i), & 
!        & radToDegUnit(cartToSphPos(vorVx(i))), vorVx(i)%v_(:)
!end do


  !
  !
  write(*,*) "* create an information of topology.." 

  do nowSiteId=1, siteNum

    do vxId=1, vorVxNum
      if(vorVx2vorRCId(1,vxId) == nowSiteId .or. vorVx2vorRCId(2,vxId) ==nowSiteId &
        & .or. vorVx2vorRCId(3,vxId) ==nowSiteId) exit;
    end do

    nextVxId = vxId; endVxId = vxId
    cellVxNum = 0
    do while(.true.)
      cellVxNum = cellVxNum + 1
      vorRCList(nowSiteId)%vertIdList(cellVxNum) = nextVxId

      do i=1,3
        if(nowSiteId==vorVx2vorRCId(i,nextVxId)) then
          nextVxId = vorVx2vorVxId(mod(i,3)+1,nextVxId)
          exit
        end if
      end do

      if(nextVxId == endVxId) exit;
    end do

    vorRCList(nowSiteId)%vertNum = cellVxNum
!    write(*,*) "site=",nowSiteId, ";", vorRCList(nowSiteId)%vertIdList(1:cellVxNum)
  end do


end subroutine SVoronoiDiagramGen

subroutine SVoronoiGen_Final()

  if(allocated(vorVx)) deallocate(vorVx)
  if(allocated(vorVx2VorVxId)) deallocate(vorVx2VorVxId)
  if(allocated(vorVx2VorRcId)) deallocate(vorVx2VorRcId)
  if(allocated(vorVxRadius)) deallocate(vorVxRadius)
  if(allocated(vorRCList)) deallocate(vorRCList)

end subroutine SVoronoiGen_Final


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SVoronoiGen_prepair(pts, siteIds)
  
  type(Vector3d), intent(in) :: pts(:)
  integer, intent(in) :: siteIds(4)

  ! Work variables
  !
  integer :: vxId, i
  integer :: neighVRIds(3) ! The index of three voronoi regions which are a neighbor of a voronoi vertex.  
  real(DP) :: dist_VorVx2unUsedSite

  ! Executable statements

  ! Define the mapping array to convert the index of a voronoi vertex
  ! into idecies of three voronoi regions. 

  ! Initialize the mapping array to convert the index of a voronoi vertex
  ! into the index of three vornoi vertecies which is a neighbor of it. 

  vorVx2VorVxId(:,1) = (/ 2,3,4 /)
  vorVx2VorVxId(:,2) = (/ 1,3,4 /)
  vorVx2VorVxId(:,3) = (/ 1,2,4 /)
  vorVx2VorVxId(:,4) = (/ 1,2,3 /)

  vorVx2VorRcId(:,1) = siteIds( vorVx2VorVxId(:,1) )
  vorVx2VorRcId(:,2) = siteIds( vorVx2VorVxId(:,2) )
  vorVx2VorRcId(:,3) = siteIds( vorVx2VorVxId(:,3) )
  vorVx2VorRcId(:,4) = siteIds( vorVx2VorVxId(:,4) )

  vorIdMapHead = 4

  ! Calculate the coordinates of four voronoi vetecies.
  do vxId=1, 4
     neighVRIds(:) = vorVx2VorRcId(:, vxId)

     VorVx(vxId) = calcUniSTriCenterPt( pts(neighVRIds(:)) ) 

     ! Check if vorVx2VorRcId is defined correctly.
     ! If the definition of id mapping is not correct, 
     ! the voronoi vertecies will be recomputed after the correction of index. 
     dist_VorVx2unUsedSite = geodesicArcLength( pts(siteIds(vxId)), VorVx(vxId) )
     do i=1,3
        if( geodesicArcLength(pts(neighVrIds(i)), VorVx(vxId)) > dist_VorVx2unUsedSite ) then
           call swap(neighVrIds(2), neighVrIds(3))
           call swap(vorVx2VorVxId(2,vxId), vorVx2VorVxId(3,vxId))
           VorVx(vxId) = calcUniSTriCenterPt( pts(neighVRIds(:)) ) 
           
           vorVx2VorRcId(:,vxId) = neighVRIds(:)
           exit     
        end if
     end do

     vorVxRadius(vxId) = geodesicArcLength( pts(neighVRIds(1)), VorVx(vxId) )

  end do
  
end subroutine SVoronoiGen_prepair

subroutine getBrokenVorVxIdList(newSitePos, brokenVxIdList, bVxFlag)
  type(Vector3d), intent(in) :: newSitePos
  integer, intent(inout), allocatable :: brokenVxIdList(:)
  logical, intent(inout) :: bVxFlag(:)

  integer :: startBrokenVxId
  integer :: parentVxId, childVxId
  integer :: tmpList(siteNum/3)
  integer :: lEndPtr
  integer :: searchPtr
  integer :: i

  
  bVxFlag(:) = .false.
  startBrokenVxId = -1
  call searchNeighVorVxId(newSitePos, vorIdMapHead, bVxFlag, startBrokenVxId)

  !
  bVxFlag(:) = .false.
  tmpList(:) = -1
  lEndPtr = 0; searchPtr = 0

  lEndPtr = lEndPtr + 1
  tmpList(lEndPtr) = startBrokenVxId; 
  bVxFlag(startBrokenVxId) = .true.
  
  do while(lEndPtr /= searchPtr)
    searchPtr = searchPtr + 1
    parentVxId = tmpList(searchPtr); 

     do i=1,3
        childVxId = vorVx2vorVxId(i,parentVxId)

        if ( (.not. bVxFlag(childVxId)) .and. &
           & geodesicArcLength(newSitePos, vorVx(childVxId)) < vorVxRadius(childVxId) ) then

           lEndPtr = lEndPtr + 1
           tmpList(lEndPtr) = childVxId; 
           bVxFlag(childVxId) = .true.

        end if
     end do
  end do

  !
  if(allocated(brokenVxIdList)) deallocate(brokenVxIdList)

  allocate(brokenVxIdList(lEndPtr))
!  write(*,*) "newSite:", radtodegunit(carttosphpos(newsitepos))
!  write(*,*) "tmpbroken:", startBrokenVxId, tmpList(1:lEndPtr)
  brokenVxIdList(1:lEndPtr) = tmpList(1:lEndPtr)

end subroutine getBrokenVorVxIdList

recursive subroutine searchNeighVorVxId(newSitePos, startVxId, judgedFlag, foundVxId)
  type(Vector3d), intent(in) :: newSitePos
  integer, intent(in) :: startVxId
  logical, intent(inout) :: judgedFlag(:)
  integer, intent(inout) :: foundVxId

  integer :: i

  if(.not. judgedFlag(startVxId)) then
     judgedFlag(startVxId) = .true.
     
     if ( geodesicArcLength(newSitePos, vorVx(startVxId)) < vorVxRadius(startVxId) ) then
        foundVxId = startVxId
        return
     else
        do i=1, 3
           call searchNeighVorVxId(newSitePos, vorVx2vorVxId(i, startVxId), judgedFlag, foundVxId)
        end do
        
     end if
  end if

end subroutine searchNeighVorVxId

subroutine construct_newVoronoiRegion( newSiteId, bVxIdList, bVxFlag, pts )
  integer, intent(in) :: newSiteId
  integer, intent(in) :: bvxIdList(:)
  logical, intent(in) :: bVxFlag(:)
  type(Vector3d), intent(in) :: pts(:)

  integer :: i, j, k, bVxId
  integer :: neighVxId, newVxId, newVorVxIds(size(bVxIdList)*3)
  integer :: newVxNum
  integer :: newVorVxIdLink(siteNum)
  integer :: endNewVxId, nextNewVxId, prevNewVxId
  integer :: rcids(3)
  integer :: tmpVx2RCId(3, size(bVxIdList)*3)
  integer :: tmpVx2VxId(3, size(bVxIdList)*3)
  integer :: neighVxCorrectedIdList(size(bVxIdList)*3)
  integer :: cids(3)

  newVxNum = 0 
  newVorVxIdLink(:) = -1

!do i=1,vorVxNum
!   write(*,'(i2,a,3i3,a,3i3, a, f12.5, 2f12.5)') i, "*", vorVx2vorRcId(:,i), ":", vorVx2vorVxId(:,i), &
!        & ":",  VorVxRadius(i), radToDegUnit(cartToSphPos(vorVx(i)))
!end do

  do i=1, size(bVxIdList)
     bVxId = bVxIdList(i)

     cids(:) = (/ 1, 2, 3 /)

     do j=1, 3
        neighVxId =  vorVx2vorVxId(j,bVxId)

        if( .not. bVxFlag(neighVxId) ) then

           newVxNum = newVxNum + 1

           if( newVxNum > size(bVxIdList)) then
              vorIdMapHead = vorIdMapHead + 1
              newVxId = vorIdMapHead
           else
              newVxId = bVxIdList(newVxNum)
           end if
           
           newVorVxIds(newVxNum) = newVxId

           tmpVx2RCId(:,newVxNum) = (/ &
                & newSiteId, vorVx2vorRCId(cids(2),bVxId), vorVx2vorRCId(cids(3),bVxId) /)

           VorVx(newVxId) = calcUniSTriCenterPt( pts(tmpVx2RCId(:, newVxNum)) ) 
           vorVxRadius(newVxId) = geodesicArcLength( pts(newSiteId), VorVx(newVxId) )
           
           tmpVx2VxId(1, newVxNum) = neighVxId
           do k=1, 3
             if(bVxId==vorVx2vorVxId(k,neighVxId)) neighVxCorrectedIdList(newVxNum) = k
           end do

           newVorVxIdLink(vorVx2vorRCId(cids(2),bVxId)) = newVxId

        end if

        cids = cshift(cids, 1)
     end do
  end do

  do i=1, newVxNum
    neighVxId = tmpVx2VxId(1, i)
    vorVx2vorVxId(neighVxCorrectedIdList(i), neighVxId) = newVorVxIds(i)
  end do

  vorVx2VorRCId(:,newVorVxIds) = tmpVx2RCId(:, 1:newVxNum)
  vorVx2VorVxId(:,newVorVxIds) = tmpVx2VxId(:, 1:newVxNum)

  !
  !write(*,*) "link:", newVorVxIdLink(:)
  endNewVxId = newVorVxIds(newVxNum)
  nextNewVxId = newVorVxIds(newVxNum)

  do while(.true.)

     prevNewVxId = nextNewVxId
     nextNewVxId = newVorVxIdLink(vorVx2vorRCId(3,prevNewVxId))

!write(*,*) prevNewVxId, "-.>", nextNewVxId
     if(nextNewVxId == -1) then
        write(*,*) "error"; stop
     end if
     vorVx2VorVxId(3, nextNewVxId) = prevNewVxId
     vorVx2VorVxId(2, prevNewVxId) = nextNewVxId
     neighVxId = vorVx2VorVxId(1, prevNewVxId)

     if ( nextNewVxId == endNewVxId ) exit
  end do

end subroutine construct_newVoronoiRegion

subroutine SVoronoi_setTopology(mesh)
  type(PolyMesh), intent(inout) :: mesh
  
  integer :: siteId
  type(Face) :: vorRC
  type(Face) :: tmpFace
  type(Cell) :: tmpCell
  integer :: edgeId, vId, vNum
  integer :: vId1, vId2, pairVorRCId
  integer :: m, n
  integer :: tmpVx2RCId(3)
  integer :: vIds(6), faceNum

  do vId=0, vorVxNum
    call PolyMesh_SetPoint(mesh, vId, vorVx(vId))
  end do
  
  faceNum = 0
  do siteId=1, siteNum
    vorRC = vorRCList(siteId)
    vNum = vorRC%vertNum
    tmpFace%vertNum = 2
    mesh%cellList(siteId)%faceNum = vNum

    do vId=1, vNum
      tmpVx2RCId(:) = vorVX2vorRCId(:,vorRC%vertIdList(vId))
      do m=1,3
        if(siteId == tmpVx2RCId(m)) then
          n = mod(m,3) + 1; 
          exit
        end if
      end do

      if(tmpVx2RCId(n) < siteId) then
!        write(*,*) "#reg v=", vId, "*tmp:", tmpVx2RCId

        tmpFace%vertIdList(1:2) = (/ vorVx2vorVxId(mod(n,3)+1, vorRC%vertIdList(vId)), vorRC%vertIdList(vId) /)
        tmpFace%ownCellId = siteId
        tmpFace%neighCellId = tmpVx2RCId(n)

!write(*,*) "Face; own=", siteId, "neigh=", tmpVx2RCId(n)

        faceNum = faceNum + 1
        call PolyMesh_SetFace(mesh, faceNum, tmpFace)

        mesh%cellList(siteId)%faceIdList(vId) = faceNum

        do m=1, vorRCList(tmpFace%neighCellId)%vertNum
          if( tmpFace%vertIdList(1) == vorRCList(tmpFace%neighCellId)%vertIdList(m) ) then
            mesh%cellList(tmpFace%neighCellId)%faceIdList(m) = faceNum
 !         write(*,*) "found: m=", m
            exit
          end if
          if(m == vorRCList(tmpFace%neighCellId)%vertNum) then
            write(*,*) "error:", tmpFace%vertIdList(1), vorRCList(tmpFace%neighCellId)%vertIdList(:), ":",  &
              & vorRCList(tmpFace%neighCellId)%vertNum
          end if
        end do

      end if

    end do
  end do

  write(*,*) "cell info--"
  do m=1, siteNum
    write(*,*) m, ":", mesh%cellList(m)%faceIdList(:), ":", mesh%cellList(m)%faceNum
  end do

end subroutine SVoronoi_setTopology

function calcUniSTriCenterPt(pts) result(centerPt)
  type(Vector3d), intent(in) :: pts(3)
  type(Vector3d) :: centerPt

  centerPt = normalizedVec( (pts(2) - pts(1)).cross.(pts(3) - pts(1)) ) 

end function calcUniSTriCenterPt

subroutine swap(i1, i2)
  integer, intent(inout) :: i1, i2
  integer :: tmp

  tmp = i1
  i1 = i2
  i2 = tmp
end subroutine swap

end module SVoronoiGen_mod
