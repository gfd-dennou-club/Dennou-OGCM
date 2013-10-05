program ogcm_main

  ! �⥸�塼�����; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use Constants_mod
  use TemporalIntegSet_mod
  use GridSet_mod
  use VariableSet_mod

  ! ���ʸ; Declaration statement
  !
  implicit none


  ! �ɽ��ѿ�
  ! Local variables
  !
  
  character(*), parameter :: PROGRAM_NAME = "ogcm_main"

  ! �¹�ʸ; Executable statement
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

    ! �⥸�塼�����; Use statements
    !
    use OptionParser_mod

    ! ���ʸ; Declaration statement
    !
    
    
    ! �ɽ��ѿ�
    ! Local variables
    !
    character(STRING) :: configNmlFile

    ! �¹�ʸ; Executable statement
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
    
    ! ���ʸ; Declaration statement
    !
    
    
    ! �ɽ��ѿ�
    ! Local variables
    !
    
    
    ! �¹�ʸ; Executable statement
    !
    
    call TemporalIntegSet_Final()
    call GridSet_Final()
    call VariableSet_Final()
    call Constants_Final()

  end subroutine ogcm_finalize

end program ogcm_main
