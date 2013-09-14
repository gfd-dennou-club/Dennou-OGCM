program gridGen

  use dc_types
  use dc_message

  use VectorSpace_mod
  use SphericalCoord_mod

  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use netcdfDataWriter_mod
  use vtkDataWriter_mod


  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  type(netcdfDataWriter) :: ncwriter
  type(vtkDataWriter) :: vtkWriter

  integer, parameter :: glevel = 5
  integer, parameter :: maxItrNum = 2800
  logical, parameter :: outputVtkData = .true.
  character(STRING), parameter :: fileName = 'grid-glevel5'

  ! Setup
  write(*,*) 'hexgonal grid generation.. :glevel:', glevel
  call HexTriIcMesh_Init(htiMesh)
  call HexTriIcMesh_generate(htiMesh, glevel, maxItrNum)

  call fvMeshInfo_Init(fvInfo, htiMesh%mesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)

  !
  write(*,*) "* output netcdf .."
  call netcdfDataWriter_Init(ncwriter, trim(fileName)//".nc", htiMesh%mesh)
  call netcdfDataWriter_Regist(ncwriter, (/ fvInfo%v_CellVol /))
  call netcdfDataWriter_write(ncwriter, fvInfo%v_CellVol)
  call netcdfDataWriter_Final(ncwriter)

  if( outputVtkData ) then
     call vtkDataWriter_Init(vtkwriter, trim(fileName)//".vtk", htiMesh%mesh)
     call vtkDataWriter_Regist(vtkwriter, (/ fvInfo%v_CellVol /))
     call vtkDataWriter_write(vtkwriter)
     call vtkDataWriter_Final(vtkwriter)
  end if

  !
  call checkGridQuality()

  ! Finalization
  call HexTriIcMesh_Final(htiMesh)
  
contains


subroutine checkGridQuality()

  integer :: cellNum, faceNum, ptNum, i, j
  integer :: lfaceNum, neighCellId, faceId, ptId
  
  real(DP) :: dists(6)
  real(DP) :: sigMin, sigAvg, qMin, qAvg
  type(volScalarField) :: sigma
  type(pointScalarField) :: q
  type(PolyMesh), pointer :: plMesh

  plMesh => fvInfo%mesh

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)
  ptNUm = getPointListSize(plMesh)

  call GeometricField_Init(sigma, plMesh, "sigma")
  call GeometricField_Init(q, plMesh, "q")

  do i=1, cellNum
     lfaceNum = plMesh%cellList(i)%faceNum

     do j=1, lfaceNum
        faceId = fvInfo%Cell_FaceId(j,i)
        if( fvInfo%Face_CellId(1,faceId) == i) then
           neighCellId = fvInfo%Face_CellId(2,faceId)
        else
           neighCellId = fvInfo%Face_CellId(1,faceId)
        end if
!        dists(j) = geodesicArcLength( plMesh%CellPosList(i), plMesh%CellPosList(neighCellId) )
        dists(j) = l2norm(plMesh%CellPosList(i) - plMesh%CellPosList(neighCellId) )
     end do

     sigma%data%v_(1,i) = minval(dists(1:lfaceNum)) / maxval(dists(1:lfaceNum))
!write(*,*) i, sum(dists(1:lfaceNum))/lfaceNum
  end do

  do ptId=1, ptNum
     q%data%v_(1,ptId) = evalTriQuality(plMesh%CellPosList(fvInfo%Point_CellId(1:3, ptId)))
  end do


  sigMin = minval(sigma%data%v_ )
  sigAvg = sum(sigma%data%v_)/dble(cellNum)
  call MessageNotify("M", "SCVGridGen::CheckGridQuality", &
       & "minimum of local Uniformity for VoronoiMesh=%f", d=(/ sigMin /) )
  call MessageNotify("M", "SCVGridGen::CheckGridQuality", &
       & "average of local Uniformity for VoronoiMesh=%f", d=(/ sigAvg /) )

  qMin = minval(q%data%v_ )
  qAvg = sum(q%data%v_)/dble(ptNum)
  call MessageNotify("M", "SCVGridGen::CheckGridQuality", &
       & "minimum of local Uniformity for dual DelaunayTriMesh=%f", d=(/ qMin /) )
  call MessageNotify("M", "SCVGridGen::CheckGridQuality", &
       & "average of local Uniformity for dual DelaunayTriMesh=%f", d=(/ qAvg /) )

  call GeometricField_Final(sigma)
  call GeometricField_Final(q)

end subroutine checkGridQuality

function evalTriQuality(pts) result(ret)
  type(Vector3d), intent(in) :: pts(3)
  real(DP) :: ret

  real(DP) :: a, b, c, s

  a = geodesicArcLength(pts(1), pts(2))
  b = geodesicArcLength(pts(2), pts(3))
  c = geodesicArcLength(pts(3), pts(1))
  s = a + b + c
  ret = (s-2d0*a)*(s-2d0*b)*(s-2d0*c)/(a*b*c)

end function evalTriQuality



end program gridGen
