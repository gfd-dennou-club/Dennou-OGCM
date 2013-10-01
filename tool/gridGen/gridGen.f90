program gridGen

  use dc_types, only: &
       & DP, STRING
  use dc_message, only: &
       & MessageNotify
  use dc_string, only: &
       & StoI, StoA, CPrintf, Replace

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

  logical, parameter :: outputVtkData = .true.
  character(STRING), parameter :: fileName = 'grid-glevel5'

  integer :: glevel
  integer :: maxItrNum
  logical :: reExecuteFlag
  character(STRING) :: loadGridFileName
  character(STRING) :: genGridFileName

  type(PolyMesh) :: plMesh   ! This variable will be used only if re-executing is specified. 

  !****************************
  !* Executable statement
  !*****************************

  ! Setup
  !

  call argAnalysis()
  
  !

  ! Generate spherical constrained voronoi(SCV) grid
  !
  call MessageNotify("M", "SCVGridGen",  'hexgonal grid generation.. glevel=%d', i=(/ glevel /))

  if( reExecuteFlag ) then
     call loadGridData()
     call HexTriIcMesh_Init(htiMesh, globalMesh=plMesh)
  else
     call HexTriIcMesh_Init(htiMesh)
  end if
 
  call HexTriIcMesh_generate(htiMesh, glevel, maxItrNum)

  call fvMeshInfo_Init(fvInfo, htiMesh%globalMesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)


  ! Valify the quality of the generated grid
  !

  call checkGridQuality()

  ! Generate a grid data file
  !

  call MessageNotify("M",  "SCVGridGen", "* output grid data as a netcdf file ..")
  call netcdfDataWriter_Init(ncwriter, genGridFileName, htiMesh%globalMesh)
  call netcdfDataWriter_writeGlobalAttr(ncWriter, "planetRadius", 1d0)
  call netcdfDataWriter_Regist(ncwriter, (/ fvInfo%v_CellVol /))
  call netcdfDataWriter_write(ncwriter, fvInfo%v_CellVol)
  call netcdfDataWriter_Final(ncwriter)

  if( outputVtkData ) then
  call MessageNotify("M",  "SCVGridGen", "* output grid data as a vtk data file  ..")
     call vtkDataWriter_Init(vtkwriter, Replace(genGridFileName,".nc",".vtk"), htiMesh%globalMesh)
     call vtkDataWriter_Regist(vtkwriter, (/ fvInfo%v_CellVol /))
     call vtkDataWriter_write(vtkwriter)
     call vtkDataWriter_Final(vtkwriter)
  end if


  ! Finalization
  !
  call HexTriIcMesh_Final(htiMesh)
  
  if( reExecuteFlag ) then
     call PolyMesh_Final(plMesh)
  end if

contains
subroutine argAnalysis()

  use dc_args

  type(ARGS) :: arg

  logical :: optLoadFileName, optGenFileName
  character(STRING) :: valLoadFileName, valGenFileName
  character(STRING), pointer :: argStr(:) => null()

  call DCArgsOpen(arg)

  call DCArgsOption(arg, StoA('-r', '--reusedGridFile'), optLoadFileName, valLoadFileName, &
    &         help="Specify the grid file used in re-executing SCVGridGen." )

  call DCArgsOption(arg, StoA('-o', '--outputFile'), optGenFileName, valGenFileName, &
    &         help="Specify the output grid file." )


  call DCArgsDebug(arg)
  call DCArgsHelp(arg)
  call DCArgsStrict(arg)

  if ( DCArgsNumber(arg) /= 2 ) then
     call MessageNotify('E', 'SCVGridGen', 'Two aruguments must be specified. &
          & First argument is the "glevel of icosahedral grid". &
          & Second argument is the maximum number of iteration.')
  else
     call DCArgsGet(arg, argStr)
     glevel = StoI( argStr(1) )
     maxItrNum = StoI( argStr(2) )
     if(associated(argStr)) deallocate(argStr)
  end if

  if( optGenFileName ) then
     genGridFileName = valGenFileName
  else
     genGridFileName = cprintf( "grid-glevel%d.nc", i=(/ glevel /) )
  end if

  if( optLoadFileName ) then
     reExecuteFlag = .true. 
     loadGridFileName = valLoadFileName
     call MessageNotify("M", "SCVGridGen", "Re-executing mode.. '%a' is used as Original grid data", ca=(/ loadGridFileName /))
  else
     reExecuteFlag = .false.
  end if

end subroutine argAnalysis

subroutine loadGridData()
  use netcdfDataReader_mod
  
  type(netcdfDataReader) :: reader

  call netcdfDataReader_Init(reader, loadGridFileName, plMesh)
  call netcdfDataReader_Final(reader)

write(*,*) "listsize-",size(plMesh%CellPosList)

end subroutine loadGridData

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
