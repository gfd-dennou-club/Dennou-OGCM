program ogcm_main

  ! モジュール引用; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use fvCalculus_mod, only: &
       & fvCalculus_Init, fvCalculus_Final

  use TemporalIntegSet_mod, only: &
       & CurrentTime, TotalIntegTime, &
       & TemporalIntegSet_AdvanceLongTStep

  use GridSet_mod, only: &
         & globalMesh, &
         & GridSet_getLocalMesh, GridSet_getLocalMeshInfo, nLocalPMesh, &
         & GridSet_getLocalFVMInfo

  use VariableSet_mod, only: &
       & VariableSet

  use DataFileSet_mod, only: &
       & DataFileSet, &
       & DataFileSet_OutputData

  use GovernEqSolverDriver_mod, only: &
       & GovernEqSolverDriver_AdvanceTStep

  use CGridFieldDataUtil_mod, only: &
       & CGridFieldDataUtil, &
       & CGridFieldDataUtil_Set

  use VGridFieldDataUtil_mod, only: &
       & VGridFieldDataUtil, &
       & VGridFieldDataUtil_Set

  ! 宣言文; Declaration statement
  !
  implicit none


  ! 局所変数
  ! Local variables
  !
  
  character(*), parameter :: PROGRAM_NAME = "ogcm_main"
  
  type(VariableSet), allocatable :: variables(:)
  type(CGridFieldDataUtil), allocatable :: cGridUtils(:)
  type(VGridFieldDataUtil), allocatable :: vGridUtils(:)
  type(DataFileSet), allocatable :: datFiles(:)

  integer :: lcMeshId

  ! 実行文; Executable statement
  ! 

  call MessageNotify("M", PROGRAM_NAME, "Start..")

  !***********************************
  ! Set up
  !***********************************

  call ogcm_setup()

  !***********************************
  ! Set initial condition
  !***********************************
  
  do lcMeshId=1, nLocalPMesh
     call DataFileSet_OutputData( datFiles(lcMeshId), variables(lcMeshId) )
  end do

  !***********************************
  ! The loop for temporal integration
  !************************************
  

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  do while(CurrentTime < TotalIntegTime)
     
     !
     do lcMeshId=1, nLocalPMesh
        call fvCalculus_Init( GridSet_getLocalFVMInfo(lcMeshId) )

        call CGridFieldDataUtil_Set( cGridUtils(lcMeshId) )
        call VGridFieldDataUtil_Set( vGridUtils(lcMeshId) )
        call GovernEqSolverDriver_AdvanceTStep( variables(lcMeshId), GridSet_getLocalMesh(lcMeshId) )

        call fvCalculus_Final()
     end do

     !
     call TemporalIntegSet_AdvanceLongTStep()

     !
     !
     do lcMeshId=1, nLocalPMesh
        call DataFileSet_OutputData( datFiles(lcMeshId), variables(lcMeshId) )
     end do

  end do

  call MessageNotify("M", PROGRAM_NAME, "[==== Finish temporal integration ====]")

  !*************************************
  ! Finalize 
  !*************************************

  call ogcm_finalize()

contains

  !> @brief 
  !!
  !!
  subroutine ogcm_setup()
    
    ! 引用文; Use statement
    !
    use PolyMesh_mod, only: &
         & PolyMesh

    use OptionParser_mod, only: &
         & OptionParser_Init, OptionParser_Final, &
         & OptionParser_getInfo

    use Constants_mod, only: &
         & Constants_Init
    
    use GridSet_mod, only: &
         & GridSet_Init

    use TemporalIntegSet_mod, only: &
         & TemporalIntegSet_Init

    use VariableSet_mod, only: &
         & VariableSet_Init

    use DataFileSet_mod, only: &
         & DataFileSet_Init

    use CGridFieldDataUtil_mod, only: &
         & CGridFieldDataUtil_Init

    use VGridFieldDataUtil_mod, only: &
         & VGridFieldDataUtil_Init

    use EOSDriver_mod, only: &
         & EOSDriver_Init

    use GovernEqSolverDriver_mod, only: &
         & GovernEqSolverDriver_Init

    ! 宣言文; Declaration statement
    !

    
    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile

    ! 実行文; Executable statement
    !

    call MessageNotify("M", PROGRAM_NAME, "[==== Set up some modules ====]")

    !
    call OptionParser_Init()
    call OptionParser_getInfo(configNmlFile)
    call OptionParser_Final()

    !
    call Constants_Init(configNmlFile)

    !
    call GridSet_Init(configNmlFile)

    !
    call TemporalIntegSet_Init(configNmlFile)

    !
    allocate( variables(nLocalPMesh), cGridUtils(nLocalPMesh), vGridUtils(nLocalPMesh), datFiles(nLocalPMesh) )
    do lcMeshId=1, nLocalPMesh
       call VariableSet_Init( variables(lcMeshId), GridSet_getLocalMesh(lcMeshId) )
       call CGridFieldDataUtil_Init( cGridUtils(lcMeshId), GridSet_getLocalFVMInfo(lcMeshId) )
       call VGridFieldDataUtil_Init( vGridUtils(lcMeshId), GridSet_getLocalFVMInfo(lcMeshId) )
       call DataFileSet_Init(datFiles(lcMeshId), configNmlFile, GridSet_getLocalMesh(lcMeshId), variables(lcMeshId) )
    end do


    !
    call EOSDriver_Init()

    !
    call GovernEqSolverDriver_Init()

  end subroutine ogcm_setup

  !> @brief 
  !!
  !!
  subroutine ogcm_finalize()

    ! 引用文; Use statement
    !
    use Constants_mod, only: &
         & Constants_Final

    use GridSet_mod, only: &
         & GridSet_Final, nLocalPMesh

    use TemporalIntegSet_mod, only: &
         & TemporalIntegSet_Final

    use VariableSet_mod, only: &
         & VariableSet_Final

    use DataFileSet_mod, only: &
         & DataFileSet_Final

    use CGridFieldDataUtil_mod, only: &
         & CGridFieldDataUtil_Final

    use VGridFieldDataUtil_mod, only: &
         & VGridFieldDataUtil_Final

    use EOSDriver_mod, only: &
         & EOSDriver_Final

    use GovernEqSolverDriver_mod, only: &
         & GovernEqSolverDriver_Final

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !

    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", PROGRAM_NAME, "[==== Finalize some modules ====]")

    !
    call GovernEqSolverDriver_Final()
    
    !
    call EOSDriver_Final()

    !
    call TemporalIntegSet_Final()

    !
    do lcMeshId=1, nLocalPMesh
       call VariableSet_Final( variables(lcMeshId) )
       call CGridFieldDataUtil_Final( cGridUtils(lcMeshId) )
       call VGridFieldDataUtil_Final( vGridUtils(lcMeshId) )
       call DataFileSet_Final( datFiles(lcMeshId) )
    end do
    deallocate(variables, cGridUtils, vGridUtils, datFiles)


    !
    call GridSet_Final()

    !
    call Constants_Final()


  end subroutine ogcm_finalize


end program ogcm_main
