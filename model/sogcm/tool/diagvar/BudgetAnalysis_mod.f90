!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module BudgetAnalysis_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  use Constants_mod, only: &
       & RPlanet, Grav, RefDens, &
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use HydroBouEqSolverRHS_v2_mod, only: &
!!$  use HydroBouEqSolverRHS_old_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final
  
  use HydroBouEqSolverVImplProc_mod

  use DiagnoseUtil_mod
  use DiagVarEval_mod

  use DiagVarFileSet_mod, only: &
       & gtool_historyauto_info
    
!!$  use EnergyBudgetAnalysis_mod, only: &
!!$       & ENBUDGANAKEY_GLOBALMEANENERGY, ENBUDGANAKEY_ENERGYBUDGET, &
!!$       & EnergyBudgetAnalysis_Init, EnergyBudgetAnalysis_Final, &
!!$       & EnergyBudgetAnalysis_Perform

  use EnergyBudgetAnalysisT2013_mod, only: &
       & ENBUDGANAKEY_GLOBALMEANENERGY, ENBUDGANAKEY_ENERGYBUDGET,      &
       & EnergyBudgetAnalysis_Init => EnergyBudgetAnalysisT2013_Init,   &
       & EnergyBudgetAnalysis_Final => EnergyBudgetAnalysisT2013_Final, &
       & EnergyBudgetAnalysis_Perform => EnergyBudgetAnalysisT2013_Perform
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BudgetAnalysis_Init, BudgetAnalysis_Final
  public :: BudgetAnalysis_perform

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BudgetAnalysis_mod' !< Module Name



  character(*), parameter, public :: BUDGETANAKEY_ANGMOMBUDGET = 'AngMomBudget'
  character(*), parameter, public :: BUDGETANAKEY_ANGMOMAVG = 'AngMomAvg'


  type(gt_history) :: hst_angMomBudget
  logical :: angMomBudgAnaFlag

contains

  !>
  !!
  !!
  subroutine BudgetAnalysis_Init(diagVar_gthsInfo, budgetAnaName)

    ! モジュール引用; Use statements
    !
    
    ! 宣言文; Declare statements
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    character(*), intent(in) :: budgetAnaName(:)
    
    ! 作業変数
    ! Work variables
    !
    integer :: n
    logical :: globalMeanEnergyFlag
    logical :: energyBudgAnaFlag
    
    ! 実行文; Executable statements
    !

    globalMeanEnergyFlag = .false.
    energyBudgAnaFlag = .false.
    angMomBudgAnaFlag = .false.

    do n=1, size(budgetAnaName)
       select case(budgetAnaName(n))
       case (ENBUDGANAKEY_GLOBALMEANENERGY)
          globalMeanEnergyFlag = .true.
       case (ENBUDGANAKEY_ENERGYBUDGET) 
          energyBudgAnaFlag = .true.
       case (BUDGETANAKEY_ANGMOMBUDGET) 
          angMomBudgAnaFlag = .true.
       case ('')
       case Default
          call MessageNotify('E', module_name, &
               & "The specified type of budget analysis '%c' is invalid.", c1=trim(budgetAnaName(n)) )
       end select
    end do

    !
    call EnergyBudgetAnalysis_Init(diagVar_gthsInfo, globalMeanEnergyFlag, energyBudgAnaFlag)
    
    !
    call prepair_Output(diagVar_gthsInfo)

    !
    call HydroBouEqSolverRHS_Init()
    call HydroBouEqSolverVImplProc_Init()

  end subroutine BudgetAnalysis_Init

  !>
  !!
  !!
  subroutine BudgetAnalysis_Final()

    ! 実行文; Executable statements
    !
    if( angMomBudgAnaFlag ) then
       call HistoryClose(hst_angMomBudget)
    end if

    call HydroBouEqSolverRHS_Final()
    call HydroBouEqSolverVImplProc_Final()

  end subroutine BudgetAnalysis_Final

  subroutine BudgetAnalysis_perform( )

    ! 宣言文; Declare statements
    !
 
    ! 実行文; Executable statements
    !

    call EnergyBudgetAnalysis_perform()
    
    if( angMomBudgAnaFlag ) call analyze_angMomBudget()

  contains

    subroutine analyze_angMomBudget()

      ! モジュール引用; Use statements
      !
      use VariableSet_mod, only: xyz_UN
      
      use DiagVarEval_mod, only: &
           & eval_angularMomAvg

      ! 作業変数
      ! Work variables
      !
      real(DP) :: AngMomAvg

      ! 実行文; Executable statements
      !
      
      call MessageNotify("M", module_name, "Analyze angular momentum budget..")

      AngMomAvg = eval_angularMomAvg(xyz_UN)
      call HistoryPut( BUDGETANAKEY_ANGMOMAVG, AngMomAvg, hst_angMomBudget )

    end subroutine analyze_angMomBudget


  end subroutine BudgetAnalysis_perform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine prepair_Output(diagVar_gthsInfo)

    ! 宣言文; Declaration statement
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if ( angMomBudgAnaFlag ) then

       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'AngMomBudget.nc', title='energy budget analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(0), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_angMomBudget )  

       call HistoryAddVariable( &
            & varname=BUDGETANAKEY_ANGMOMAVG, dims=(/'t'/), & 
            & longname='global mean of relative angular momentum', units='kg*m2*s-1', xtype='double',&
             & history=hst_angMomBudget)
    end if

  end subroutine prepair_Output

end module BudgetAnalysis_mod

