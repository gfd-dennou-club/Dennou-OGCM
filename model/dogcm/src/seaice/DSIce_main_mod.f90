!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief Main part of sea ice model
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_main_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5

  use dc_types, only:               &
       & DP, STRING

  use dc_message, only:             &
       & MessageNotify

  !* Dennou-OGCM
  
  use DSIce_Admin_Constants_mod, only:         &
       & DSIce_Admin_Constants_Init,           &
       & DSIce_Admin_Constants_Final
  
  use DSIce_Admin_TInteg_mod, only:            &
       & DSIce_Admin_TInteg_Init,              &
       & DSIce_Admin_TInteg_Final,             &
       & DSIce_Admin_TInteg_AdvanceLongTStep,  &       
       & EndTemporalInteg,                     &
       & CurrentTimeStep,                      &
       & TLN => TIMELV_ID_N

!!$  use DSIce_Admin_BC_mod, only:                &
!!$       & DSIce_Admin_BC_Init,                  &
!!$       & DSIce_Admin_BC_Final

  use DSIce_Admin_Grid_mod, only:            &
       & DSIce_Admin_Grid_Init,              &
       & DSIce_Admin_Grid_Final
  
  use DSIce_Admin_GovernEq_mod, only:        &
       & DSIce_Admin_GovernEq_Init,          &
       & DSIce_Admin_GovernEq_Final
  
  use DSIce_Admin_Variable_mod, only:           &
       & DSIce_Admin_Variable_Init,             & 
       & DSIce_Admin_Variable_Final,            &
       & DSIce_Admin_Variable_AdvanceTStep,     &
       & DSIce_Admin_Variable_regist_OuputVars, &
       & DSIce_Admin_Variable_HistPut,          &
       & DSIce_Admin_Variable_HistGet,          &
       & DSIce_Admin_Variable_RestartPut,          &
       & xya_SIceCon, xya_IceThick, xya_SnowThick, &
       & xya_SIceSfcTemp, xyza_SIceTemp

  
  use DSIce_IO_History_mod, only:      &
       & DSIce_IO_History_Init,              &
       & DSIce_IO_History_Final,             &
       & DSIce_IO_History_Create,            &
       & DSIce_IO_History_Output,            &
       & DSIce_IO_History_HistPut
  
  use DSIce_IO_Restart_mod, only:            &
       & DSIce_IO_Restart_Init,              &
       & DSIce_IO_Restart_Final,             &
       & DSIce_IO_Restart_Create,            &
       & DSIce_IO_Restart_Input,             &
       & DSIce_IO_Restart_Output,            &
       & RestartFlag  

  use DSIce_Boundary_driver_mod, only: &
       & DSIce_Boundary_driver_Init,              &
       & DSIce_Boundary_driver_Final,             &
       & DSIce_Boundary_driver_UpdateBeforeTstep, &
       & DSIce_Boundary_driver_UpdateAfterTstep

  use DSIce_Boundary_vars_mod, only:  &
       & DSIce_Boundary_vars_HistPut
  
  use DSIce_TInt_driver_mod, only:   &
       & DSIce_TInt_driver_Init,           &
       & DSIce_TInt_driver_Final,          &
       & DSIce_TInt_driver_Do
  
  use ProfUtil_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_main_Init, DSIce_main_Final
  public :: DSIce_main_Setup, DSIce_main_Shutdown
  public :: DSIce_main_Advance_timestep

  ! For coupler
  public :: DSIce_main_update_OcnField
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !

  logical :: CoupledRunFlag
  character(*), parameter:: module_name = 'DSIce_main_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DSIce_main_Init( isCoupledRun  &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    logical, intent(in), optional :: isCoupledRun
    
    ! 実行文; Executable statements
    !

    CoupledRunFlag = .false.
    if(present(isCoupledRun)) CoupledRunFlag = isCoupledRun
    
  end subroutine DSIce_main_Init

  !>
  !!
  !!
  subroutine DSIce_main_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_main_Final

!-------------------------------------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DSIce_main_advance_timestep( &
       & tstep,                           & ! (in)
       & loop_end_flag,                   & ! (inout)
       & skip_flag                        & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: tstep
    logical, intent(inout) :: loop_end_flag
    logical, intent(in) :: skip_flag
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !

    
    if( EndTemporalInteg() ) then
       loop_end_flag = .true.; return
    else
       loop_end_flag = .false.
    end if

    call ProfUtil_RapStart('SIceComp', 1)
    
    if(tstep == 0) then
       call DSIce_IO_Restart_Input()
       if (RestartFlag) then
          call DSIce_Admin_Variable_HistGet()
       end if
    else
       !-  Advace one time step ----------------------------------------------------
       
       call DSIce_Boundary_driver_UpdateBeforeTstep( &
            & xya_SIceCon(:,:,TLN), xya_IceThick(:,:,TLN), xya_SnowThick(:,:,TLN), & ! (in)
            & xya_SIceSfcTemp(:,:,TLN), xyza_SIceTemp(:,:,:,TLN)                   & ! (in)
            & )    

       if ( .not. skip_flag ) then
          if(CurrentTimeStep == 1 .and. (.not. RestartFlag)) then
             call DSIce_TInt_driver_Do( isSelfStartSchemeUsed=.true. )
          else
             call DSIce_TInt_driver_Do( isSelfStartSchemeUsed=.false. )
          end if
          
          call DSIce_Admin_Variable_AdvanceTStep()
       end if

       call DSIce_Admin_TInteg_AdvanceLongTStep()
    end if

    call DSIce_Boundary_driver_UpdateAfterTstep( &
         & xya_SIceCon(:,:,TLN), xya_IceThick(:,:,TLN), xya_SnowThick(:,:,TLN), & ! (in)
         & xya_SIceSfcTemp(:,:,TLN), xyza_SIceTemp(:,:,:,TLN)                   & ! (in)
         & )    
    
    call ProfUtil_RapEnd('SIceComp', 1)       

    !- Output  --------------------------------------------------------------
    !

    ! history
    call DSIce_IO_History_Output()
    call DSIce_Admin_Variable_HistPut()
    call DSIce_Boundary_vars_HistPut()

    ! restart
    call DSIce_IO_Restart_Output()
    call DSIce_Admin_Variable_RestartPut()
    
  end subroutine DSIce_main_advance_timestep

!-------------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_main_setup( configNmlFileArg )


    ! モジュール引用; Use statements
    !
    use OptionParser_mod

    use DSIce_Admin_Grid_mod, only: &
         & DSIce_Admin_Grid_construct

#ifdef _OPENMP
    use omp_lib
#endif

    ! 宣言文; Declaration statement
    !
    character(STRING), intent(in), optional :: configNmlFileArg

    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile
    integer :: nThread

    ! 実行文; Executable statement
    !

    if (.not. present( configNmlFileArg )) then
       call OptionParser_Init()
       call OptionParser_GetInfo( configNmlFile )
       call OptionParser_Final()
    else
       configNmlFile = configNmlFileArg
    end if

    ! Initialize administrative modules 
    !

    ! -- Admin ------------

    call MessageNotify( 'M', module_name, 'Start to set up Dennou-SeaIce ..' ) 

    call MessageNotify( 'M', module_name, 'Set up Admin part ..' ) 
    
    call DSIce_Admin_Constants_Init( configNmlFile )
    call DSIce_Admin_TInteg_Init( configNmlFile )
    call DSIce_Admin_GovernEq_Init( configNmlFile )
    
    call DSIce_Admin_Grid_Init( configNmlFile )
    call DSIce_Admin_Grid_construct()

    ! -- IO ------------

    call MessageNotify( 'M', module_name, 'Set up IO part ..' ) 
    
    call DSIce_IO_History_Init( configNmlFile )
    call DSIce_IO_History_Create( configNmlFile )

    call DSIce_IO_Restart_Init( configNmlFile )
    call DSIce_IO_Restart_Create( configNmlFile )

    ! --------------

    call MessageNotify( 'M', module_name, 'Set up the remainder part ..' ) 
    
    call DSIce_Admin_Variable_Init()
    call DSIce_Admin_Variable_regist_OuputVars()
    
    call DSIce_Boundary_driver_Init( configNmlFile )        
    call DSIce_TInt_driver_Init( configNmlFile )

    call MessageNotify( 'M', module_name, "A sea ice component has been initialized.")
    
  end subroutine DSIce_main_setup

  
!-------------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_main_shutdown()

    ! モジュール引用; Use statements
    !
    use SpmlUtil_mod, only: &
         SpmlUtil_Final

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !

    call DSIce_Boundary_driver_Final()
    call DSIce_TInt_driver_Final()
    
    call DSIce_IO_History_Final()
    call DSIce_IO_Restart_Final()

    call DSIce_Admin_Variable_Final()
    call DSIce_Admin_Grid_Final()
!!$    call DSIce_Admin_BC_Final()
    call DSIce_Admin_TInteg_Final()
    call DSIce_Admin_GovernEq_Final()
    call DSIce_Admin_Constants_Final()

!!$    call SpmlUtil_Final()
    
  end subroutine DSIce_main_shutdown

  !-------------------------------------------------------------------------------------------------



  !-------------------------------------------------------------------------------------------------

  
  subroutine DSIce_main_update_OcnField( &
       & xy_SeaSfcTempRecv, xy_SeaSfcSaltRecv,     & ! (in)
       & xy_OcnMixLyrDepthRecv,                    & ! (in)
       & xy_UORecv, xy_VORecv                      & ! (in)
       & )

    ! モジュール引用; Use statements
    !    
    use DSIce_Admin_Constants_mod, only: &
         & IceMaskMin
    
    use DSIce_Admin_Variable_mod, only: &
         & xya_SIceCon

    use DSIce_Admin_Grid_mod, only: &
         & IS, IE, JS, JE, &
         & KS, IA, JA

    use DSIce_Boundary_vars_mod, only: &
         & xy_SeaSfcTemp, xy_SeaSfcSalt, &
         & xy_SeaSfcU, xy_SeaSfcV,       &
         & xy_OcnMixLyrDepth

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xy_SeaSfcTempRecv(IA,JA)
    real(DP), intent(in) :: xy_SeaSfcSaltRecv(IA,JA)
    real(DP), intent(in) :: xy_OcnMixLyrDepthRecv(IA,JA)
    real(DP), intent(in) :: xy_UORecv(IA,JA)
    real(DP), intent(in) :: xy_VORecv(IA,JA)

    ! 実行文; Executable statement
    !
    
    !$omp parallel
    !$omp workshare
    xy_SeaSfcTemp(:,:) = xy_SeaSfcTempRecv
    xy_SeaSfcSalt(:,:) = xy_SeaSfcSaltRecv
    xy_OcnMixLyrDepth(:,:) = xy_OcnMixLyrDepthRecv
    xy_SeaSfcU(:,:) = xy_UORecv
    xy_SeaSfcV(:,:) = xy_VORecv
    !$omp end workshare
    !$omp end parallel
    
  end subroutine DSIce_main_update_OcnField
  
end module DSIce_main_mod

