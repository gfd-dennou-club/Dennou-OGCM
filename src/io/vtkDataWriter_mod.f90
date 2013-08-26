module vtkDataWriter_mod
  use dc_types
  use VectorSpace_mod
  use polyMesh_mod
  use GeoMetricField_mod

  implicit none
  private

  type, public :: vtkDataWriter
    type(PolyMesh), pointer :: mesh => null()
    character(STRING) :: filePath
    integer :: outputNo
  end type vtkDataWriter

  interface write_GeometricField
    module procedure write_volScalarField
  end interface write_GeometricField
  
  public :: vtkDataWriter_Init, vtkDataWriter_Final
  public :: vtkDataWriter_writeHeader, write_GeometricField

contains
subroutine vtkDataWriter_Init(writer, fileName, mesh)
  type(vtkDataWriter), intent(inout) :: writer
  character(*), intent(in) :: fileName
  type(PolyMesh), intent(in), target :: mesh

  writer%mesh => mesh
  writer%filePath = fileName
  writer%outputNo = 17

end subroutine vtkDataWriter_Init

subroutine vtkDataWriter_writeHeader(writer)
  type(vtkDataWriter), intent(inout) :: writer

  open(writer%outputNo, file=writer%filePath, status="replace")

  write(writer%outputNo,'(a)') "# vtk DataFile Version 2.0"
  write(writer%outputNo,'(a)') "Spherical Voronoi Diagram"
  write(writer%outputNo,'(a)') "ASCII"
  write(writer%outputNo,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(writer%outputNo,*)

  call write_GridData(writer)

end subroutine vtkDataWriter_writeHeader

subroutine write_volScalarField(writer, v_scalar)
  type(vtkDataWriter), intent(inout) :: writer
  type(volScalarField), intent(in) :: v_scalar

 
end subroutine write_volScalarField

subroutine vtkDataWriter_Final(writer)
  type(vtkDataWriter), intent(inout) :: writer

  close(writer%outputNo)
  writer%mesh => null()

end subroutine vtkDataWriter_Final

!
subroutine write_GridData(writer)
  type(vtkDataWriter), intent(inout) :: writer

  integer :: ptListSize, cellListSize, cellDataInfoSize
  type(Vector3d) :: pt
  integer :: vxIds(6), faceLNum
  integer :: i, j
  type(PolyMesh), pointer :: mesh
  type(Face) :: face_

  mesh => writer%mesh

  ptListSize = getPointListSize(mesh)
  cellListSize = getCellListSize(mesh)

  write(writer%outputNo,'(a,i10,a)') "POINTS ", ptListSize, " float"
  do i=1, ptListSize
    pt = mesh%pointList(i)
    write(writer%outputNo,'(3f12.5)') pt%v_(1), pt%v_(2), pt%v_(3) 
  end do
  write(writer%outputNo,*)
  
  
  cellDataInfoSize = 0
  do i=1, cellListSize
    cellDataInfoSize = cellDataInfoSize + mesh%CellList(i)%faceNum + 1
  end do

  write(writer%outputNo,'(a,i10,i10)') "CELLS", cellListSize, cellDataInfoSize
  do i=1, cellListSize
    faceLNum = mesh%CellList(i)%faceNum
    do j=1, faceLNum
      face_ = mesh%faceList(mesh%CellList(i)%faceIdList(j))
      if(i == face_%ownCellId ) then
        vxIds(j) = face_%vertIdList(1)
      else
        vxIds(j) = face_%vertIdList(2)
      end if
    end do

    write(writer%outputNo,'(i0,6i7)') faceLNum, vxIds(1:faceLNum) - 1
  end do
  write(writer%outputNo,*)

  write(writer%outputNo,'(a,i10)') "CELL_TYPES", cellListSize
  do i=1, cellListSize
    write(writer%outputNo,'(i0)') 7
  end do
  write(writer%outputNo,*)

end subroutine write_GridData

end module vtkDataWriter_mod
