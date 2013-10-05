program ogcm_main

  ! モジュール引用; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use Constants_mod
  use TemporalIntegSet_mod
  use GridSet_mod
  use VariableSet_mod

  ! 宣言文; Declaration statement
  !
  implicit none


  ! 局所変数
  ! Local variables
  !
  
  character(*), parameter :: PROGRAM_NAME = "ogcm_main"

  ! 実行文; Executable statement
  ! 

  call MessageNotify("M", PROGRAM_NAME, "Start..")


  call ogcm_setup()

  call ogcm_finalize()

  call MessageNotify("M", PROGRAM_NAME, "..End")


contains
  !> @brief 
  !!
  !!
  subroutine ogcm_setup()

    ! モジュール引用; Use statements
    !
    use OptionParser_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile

    ! 実行文; Executable statement
    !

    call OptionParser_Init()
    call OptionParser_GetInfo(configNmlFile)
    call OptionParser_Final()

    call Constants_Init(configNmlFile)

    call TemporalIntegSet_Init(configNmlFile)
    call GridSet_Init(configNmlFile)



    call VariableSet_Init()
    

  end subroutine ogcm_setup

  !> @brief 
  !!
  !!
  subroutine ogcm_finalize()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    call TemporalIntegSet_Final()
    call GridSet_Final()
    call VariableSet_Final()
    call Constants_Final()

  end subroutine ogcm_finalize

end program ogcm_main
