program globalSWM

  use dc_types
  use dc_message
  use VectorSpace_mod
  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod
  use GovernEquationSolver_mod, only: &
       & GovernEquationSolver_Init, GovernEquationSolver_Final, &
       & Solve_GovernEquation

  use Exp_Williamson94_Case1, only: &
       & set_InitialCondition => setIniCond_ExpWS94Case1

  implicit none

  type(volScalarField) :: v_CellVol
  integer :: tStep

  

  ! Initialize
  !

  call GridSet_Init()
  
  call VariableSet_Init()

  call fvCalculus_Init(fvmInfo)

  call GovernEquationSolver_Init()

  ! Set initial condition
  !
  call MessageNotify( 'M', "globalSWM", "Set initial condition..")
  call Set_InitialCondition()
  
  call OutputData(0)

  !
  !

  do tStep=1, endTStep

     !
     if( mod(tStep, endTStep/20) == 0 ) then
        call MessageNotify( 'M', "globalSWM", "t=%d [s] (endTime=%d [s]) ..", i=(/ tStep*delTime, endTime /) )
     end if

     !
     if( mod(tStep*delTime, outputIntrVal) == 0 ) call OutputData(tStep)

     !
     call Solve_GovernEquation()

  end do


  ! Finalize
  !
  call GovernEquationSolver_Final()
  call fvCalculus_Final()
  call VariableSet_Final()
  call GridSet_Final()

contains

subroutine OutputData( tstep )

  use vtkDataWriter_mod
  use netcdfDataWriter_mod
  use dc_string
  use SphericalCoord_mod

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName
  type(vtkDataWriter) :: writer

integer :: maxId(1)
maxId = maxloc(v_height%data%v_)
write(*,*) "max height=", v_height.At.maxId(1), ":", RadToDegUnit(CartToSphPos(plmesh%cellPosList(maxId(1))))
  dataFileName = CPrintf("data-%08d.vtk", i=(/ tstep /))
  call vtkDataWriter_Init(writer, dataFileName, plMesh)
  call vtkDataWriter_Regist(writer, (/ v_height, fvmInfo%v_CellVol, v_div /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)

end subroutine OutputData

end program globalSWM
