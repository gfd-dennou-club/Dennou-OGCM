program nc2vtk_main
  use dc_types
  use dc_message
  use dc_string

  use PolyMesh_mod
  use GeometricField_mod
  use vtkDataWriter_mod
  use netcdfDataReader_mod

  implicit none

  character(STRING) :: ncfileName, vtkfileName
  character(TOKEN), pointer :: p_FieldNames(:) => null()
  character(TOKEN), pointer :: v_FieldNames(:) => null()
  character(TOKEN), pointer :: s_FieldNames(:) => null()
  type(pointScalarField), allocatable :: p_Scalars(:)
  type(volScalarField), allocatable :: v_Scalars(:)
  type(surfaceScalarField), allocatable :: s_Scalars(:)
  real(DP) :: dataTime

  type(PolyMesh) :: plMesh

  call MessageNotify('M', "nc2vtk", "..")

  dataTime = 0d0

  call AnalizeArg()

  !
  call readDataFromNetCDF()

  !
  call outputVtkFile()

  !
  call releaseResource()


contains

subroutine AnalizeArg()

  use dc_args

  type(ARGS) :: arg
  
  logical :: optPFieldNames, optSFieldNames, optVFieldNames, optRange
  character(STRING) :: valPFieldNames, valSFieldNames, valVFieldNames, valRange
  character(STRING), pointer :: fileNames(:)
 
 !
  !

  call DCArgsOpen(arg)

  call DCArgsOption(arg, StoA('-p', '--p_FieldNames'), optPFieldNames, valPFieldNames, &
    &         help="Specify the name of some point scalar field." )

  call DCArgsOption(arg, StoA('-s', '--s_FieldNames'), optSFieldNames, valSFieldNames, &
    &         help="Specify the name of some surface scalar field." )

  call DCArgsOption(arg, StoA('-v', '--v_FieldNames'), optVFieldNames, valVFieldNames, &
    &         help="Specify the name of some volume scalar field." )

  call DCArgsOption(arg, StoA('-t', '--time'), optRange, valRange, &
    &         help="Specify the time of data." )


  call DCArgsDebug(arg)
  call DCArgsHelp(arg)
  call DCArgsStrict(arg)


  !
  !
  if ( DCArgsNumber(arg) /= 2 ) then
     call MessageNotify('E', 'nc2vtk', 'Two aruguments must be specified. &
          & First argument is a name of read netcdf file. &
          & Second argument is a name of outut vtk file.')
  else
     call DCArgsGet(arg, fileNames)
     ncfileName = fileNames(1)
     vtkfileName = fileNames(2)
     if(associated(fileNames)) deallocate(fileNames)
  end if

  if( optPFieldNames ) then
     call split(valPFieldNames, p_FieldNames, ",")
  else
     allocate(p_FieldNames(0))
  end if

  if( optSFieldNames ) then
     call split(valPFieldNames, s_FieldNames, ",")
  else
     allocate(s_FieldNames(0))
  end if

  if( optVFieldNames ) then
     call split(valVFieldNames, v_FieldNames, ",")
  else
     allocate(v_FieldNames(0))
  end if

  if( optRange ) then
     dataTime = StoD(valRange)
  end if

  !
  !

  call DCArgsClose(arg)

end subroutine AnalizeArg

subroutine readDataFromNetCDF()

  type(netcdfDataReader) :: reader
  integer :: i, fieldNum
  character(STRING) :: dataRange

  call netcdfDataReader_Init(reader, ncfileName, plMesh)

  dataRange = cprintf("time=%10f", d=(/ dataTime /))
  
  if(associated(p_FieldNames)) then
     fieldNum = size(p_FieldNames)
     allocate( p_Scalars(fieldNum) )
     do i=1, fieldNum
        call netcdfDataReader_get(reader, p_FieldNames(i), p_Scalars(i), dataRange )
     end do
  end if

  if(associated(v_FieldNames)) then
     fieldNum = size(v_FieldNames)
     allocate( v_Scalars(fieldNum) )
     do i=1, fieldNum
        call netcdfDataReader_get(reader, v_FieldNames(i), v_Scalars(i), dataRange )
     end do
  end if

  if(associated(s_FieldNames)) then
     fieldNum = size(s_FieldNames)
     allocate( s_Scalars(fieldNum) )
     do i=1, fieldNum
!        call netcdfDataReader_get(reader, s_FieldNames(i), s_Scalars(i), dataRange )
     end do
  end if

  call netcdfDataReader_Final(reader)

end subroutine readDataFromNetCDF


subroutine outputVtkFile()
  type(vtkDataWriter) :: writer

  call vtkDataWriter_Init(writer, vtkfileName, plMesh)

  call vtkDataWriter_Regist(writer, v_scalars, pointScalarFields=p_scalars)
  call vtkDataWriter_Write(writer)
  call vtkDataWriter_Final(writer)

end subroutine outputVtkFile


subroutine releaseResource()

  integer :: i

  if( associated(p_FieldNames) ) then
     do i=1, size(p_Scalars)
        call GeometricField_Final(p_Scalars(i))
     end do
     deallocate(p_Scalars)
     deallocate(p_FieldNames)
  end if

  if( associated(s_FieldNames) ) then
     do i=1, size(s_Scalars)
        call GeometricField_Final(s_Scalars(i))
     end do
     deallocate(s_Scalars)
     deallocate(s_FieldNames)
  end if

  if( associated(v_FieldNames) ) then
     do i=1, size(v_Scalars)
        call GeometricField_Final(v_Scalars(i))
     end do
     deallocate(v_Scalars)
     deallocate(v_FieldNames)
  end if

end subroutine releaseResource

end program nc2vtk_main
