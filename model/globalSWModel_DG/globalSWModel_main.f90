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

  use LagrangePolyn_mod
  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod
  use DGHelper_mod, only: &
       & DGHelper_Init, DGHelper_Final, DGHelper_prepair

  use GovernEquationSolver_mod, only: &
       & GovernEquationSolver_Init, GovernEquationSolver_Final, &
       & Solve_GovernEquation

  use OutputData_mod, only: &
       & OutputData_Init, OutputData_Final, &
       & OutputData


  !
  !
!!$  use Exp_Williamson94_Case1, only: &
!!$       & set_InitialCondition => setIniCond_ExpWS94Case1, &
!!$       & callBack_EndCurrentTimeStep

!!$  use Exp_Williamson94_Case2, only: &
!!$       & set_InitialCondition => setIniCond_ExpWS94Case2, &
!!$       & callBack_EndCurrentTimeStep

!!$  use Exp_Williamson94_Case5, only: &
!!$       & set_InitialCondition => setIniCond_ExpWS94Case5, &
!!$       & callBack_EndCurrentTimeStep

  use Exp_WindDrivenCirc_mod, only: &
       & set_InitialCondition => setIniCond_ExpWindDrivenCirc, &
       & callBack_EndCurrentTimeStep


  implicit none

  integer :: tStep
  

  ! Initialize
  !
  call LagrangePolyn_Init()

  call GridSet_Init()
  call GridSet_prepair()

  call DGHelper_Init()
  call DGHelper_prepair()

  call VariableSet_Init()

!  call fvCalculus_Init(fvmInfo)

  call GovernEquationSolver_Init()

  call OutputData_Init()

  ! Set initial condition
  !
  call MessageNotify( 'M', "globalSWM", "Set initial condition..")
  call Set_InitialCondition()
  wc_hU1 = (meanDepth + wc_h)*wc_U1
  wc_hU2 = (meanDepth + wc_h)*wc_U2
  
   call OutputData(0)
   call callBack_EndCurrentTimeStep(0, wc_h, wc_hU1, wc_hU2)

  
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

     call callBack_EndCurrentTimeStep(tstep, wc_h, wc_hU1, wc_hU2)

     if( mod(tStep*delTime, outputIntrVal) == 0 ) call OutputData(tStep)

  end do


  ! Finalize
  !
  call OutputData_Final()

  call GovernEquationSolver_Final()
!  call fvCalculus_Final()
  call VariableSet_Final()
  call DGHelper_Final()
  call GridSet_Final()

end program globalSWM
