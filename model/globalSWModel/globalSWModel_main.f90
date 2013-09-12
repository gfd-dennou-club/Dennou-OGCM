program globalSWM

  use dc_types
  use dc_message
  use VectorSpace_mod
  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use netcdfDataWriter_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  use GovernEquationSolver_mod, only: &
       & GovernEquationSolver_Init, GovernEquationSolver_Final, &
       & Solve_GovernEquation

  use OutputData_mod, only: &
       & OutputData_Init, OutputData_Final, &
       & OutputData


  !
  !
!!$  use Exp_Williamson94_Case1, only: &
!!$       & set_InitialCondition => setIniCond_ExpWS94Case1

  use Exp_Williamson94_Case2, only: &
       & set_InitialCondition => setIniCond_ExpWS94Case2, &
       & callBack_EndCurrentTimeStep

!!$  use Exp_Williamson94_Case5, only: &
!!$       & set_InitialCondition => setIniCond_ExpWS94Case5, &
!!$       & callBack_EndCurrentTimeStep



  implicit none

  integer :: tStep
  

  ! Initialize
  !

  call GridSet_Init()
  
  call VariableSet_Init()

  call fvCalculus_Init(fvmInfo)

  call GovernEquationSolver_Init()

  call OutputData_Init()

  ! Set initial condition
  !
  call MessageNotify( 'M', "globalSWM", "Set initial condition..")
  call Set_InitialCondition()
  
  call OutputData(0)
  call callBack_EndCurrentTimeStep(0, v_h, s_normalVel)
  
  !
  !

  do tStep=1, endTStep

     !
     if( mod(tStep, int(20) ) == 0 ) then
        call MessageNotify( 'M', "globalSWM", "t=%d [s] (endTime=%d [s]) ..", i=(/ tStep*delTime, endTime /) )
     end if

     !

     !
     call Solve_GovernEquation()

     call callBack_EndCurrentTimeStep(tstep, v_h, s_normalVel)

     if( mod(tStep*delTime, outputIntrVal) == 0 ) call OutputData(tStep)

  end do


  ! Finalize
  !
  call OutputData_Final()

  call GovernEquationSolver_Final()
  call fvCalculus_Final()
  call VariableSet_Final()
  call GridSet_Final()

end program globalSWM
