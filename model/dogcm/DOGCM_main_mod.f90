!-------------------------------------------------------------
! Copyright (c) 2015-2017 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module providing a few subroutines to run Dennou-OGCM. 
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_main_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5

  use dc_types, only:               &
       & DP, STRING

  use dc_message, only:             &
       & MessageNotify

  !* Dennou-OGCM
  
  use DOGCM_Admin_Constants_mod, only:         &
       & DOGCM_Admin_Constants_Init,           &
       & DOGCM_Admin_Constants_Final,          &
       & RPlanet
  
  use DOGCM_Admin_TInteg_mod, only:            &
       & DOGCM_Admin_TInteg_Init,              &
       & DOGCM_Admin_TInteg_Final,             &
       & DOGCM_Admin_TInteg_AdvanceLongTStep,  &       
       & EndTemporalInteg,                     &
       & CurrentTimeStep,                      &
       & TLN => TIMELV_ID_N,                   &
       & TLA => TIMELV_ID_A
  
  use DOGCM_Admin_BC_mod, only:                &
       & DOGCM_Admin_BC_Init,                  &
       & DOGCM_Admin_BC_Final

  use DOGCM_Admin_Grid_mod, only:            &
       & DOGCM_Admin_Grid_Init,              &
       & DOGCM_Admin_Grid_Final
  
  use DOGCM_Admin_GovernEq_mod, only:        &
       & DOGCM_Admin_GovernEq_Init,          &
       & DOGCM_Admin_GovernEq_Final,         &
       & SolverType,                         &
       & OCNGOVERNEQ_SOLVER_HSPM_VSPM,       &
       & OCNGOVERNEQ_SOLVER_HSPM_VFVM
  
  use DOGCM_Admin_Variable_mod, only:           &
       & DOGCM_Admin_Variable_Init,             & 
       & DOGCM_Admin_Variable_Final,            &
       & DOGCM_Admin_Variable_AdvanceTStep,     &
       & DOGCM_Admin_Variable_regist_OuputVars, &
       & DOGCM_Admin_Variable_HistPut,          &
       & DOGCM_Admin_Variable_HistGet,          &
       & DOGCM_Admin_Variable_RestartPut,       &
       & xyza_U, xyza_V, xyzaa_TRC, xya_SSH,    &
       & xyza_H

  use DOGCM_IO_History_mod, only:               &
       & DOGCM_IO_History_Init,                 &
       & DOGCM_IO_History_Final,                &
       & DOGCM_IO_History_Create
  
  use DOGCM_IO_Restart_mod, only:               &
       & DOGCM_IO_Restart_Init,                 &
       & DOGCM_IO_Restart_Final,                &
       & DOGCM_IO_Restart_Input,                &
       & DOGCM_IO_Restart_Create,               &
       & DOGCM_IO_Restart_Output,               &
       & RestartFlag  
  
  use DOGCM_Boundary_driver_mod, only: &
       & DOGCM_Boundary_driver_Init,              &
       & DOGCM_Boundary_driver_Final,             &
       & DOGCM_Boundary_driver_UpdateBeforeTstep, &
       & DOGCM_Boundary_driver_UpdateAfterTstep

  use DOGCM_Boundary_vars_mod, only: &
       & DOGCM_Boundary_vars_HistPut
  
  use DOGCM_TInt_driver_mod, only:     &
       & DOGCM_TInt_driver_Init,       &
       & DOGCM_TInt_driver_Final,      &
       & DOGCM_TInt_driver_Do
  
  use DOGCM_Exp_driver_mod, only:      &
       & DOGCM_Exp_driver_Init,        &
       & DOGCM_Exp_driver_Final

  use LPhys_RediGM_spm_mod, only: &
       & LPhys_RediGM_spm_Output

  use LPhys_RediGM_hspm_vfvm_mod, only: &
       & LPhys_RediGM_hspm_vfvm_Output

  use ProfUtil_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_main_Init, DOGCM_main_Final
  
  public :: DOGCM_main_setup, DOGCM_main_shutdown
  public :: DOGCM_main_advance_timestep
  public :: DOGCM_main_update_SIceField

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !

  logical :: CoupledRunFlag
  integer :: my_mpi_comm
  
  character(*), parameter:: module_name = 'DOGCM_main_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DOGCM_main_Init( mpi_comm, isCoupledRun  &  ! (in)
       &  )

    ! 宣言文; Declaration statement
    !
    integer, intent(in), optional :: mpi_comm
    logical, intent(in), optional :: isCoupledRun
    
    ! 実行文; Executable statements
    !

    my_mpi_comm = 1
    if(present(mpi_comm)) my_mpi_comm = mpi_comm
    
    CoupledRunFlag = .false.
    if(present(isCoupledRun)) CoupledRunFlag = isCoupledRun
    
  end subroutine DOGCM_main_Init

  !>
  !!
  !!
  subroutine DOGCM_main_Final()

    ! 実行文; Executable statements
    !
    
  end subroutine DOGCM_main_Final

!-------------------------------------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DOGCM_main_advance_timestep( &
       & tstep,                           & ! (in)
       & loop_end_flag,                   & ! (inout)
       & skip_flag                        & ! (in)
       & )
    
    ! モジュール引用; Use statements
    !    
!!$    use DOGCM_Admin_Constants_mod
!!$    use DOGCM_Admin_Grid_mod
    use DOGCM_Admin_Variable_mod, only: &
         & xyza_U, xyza_V, xyzaa_TRC, xyza_H, xya_SSH

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
    
    call ProfUtil_RapStart('OcnComp', 1)
    
    if(tstep == 0) then
       call DOGCM_IO_Restart_Input()
       if (RestartFlag) then
          call DOGCM_Admin_Variable_HistGet()
       end if
    else
       !-  Advace one time step ----------------------------------------------------
!!$       write(*,*) "Boundary UpdateBeforetstep"
       call DOGCM_Boundary_driver_UpdateBeforeTstep( &
            & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),  & ! (in)
            & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN)                            & ! (in)
            & )    

!!$       write(*,*) "TIntegMain"              
       if ( .not. skip_flag ) then
          if(CurrentTimeStep == 1 .and. (.not. RestartFlag)) then
             call DOGCM_TInt_driver_Do( isSelfStartSchemeUsed=.true. )
          else
             call DOGCM_TInt_driver_Do( isSelfStartSchemeUsed=.false. )
          end if
          call DOGCM_Admin_Variable_AdvanceTStep()
       end if

!!$       write(*,*) "TIntegAdvance"       
       call DOGCM_Admin_TInteg_AdvanceLongTStep()
!!$       write(*,*) "Boundary UpdateAftertstep"
    end if

    call DOGCM_Boundary_driver_UpdateAfterTstep( &
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),  & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN)                            & ! (in)
         & )
    
    call ProfUtil_RapEnd('OcnComp', 1)
    
    !- Output  --------------------------------------------------------------------

    call ProfUtil_RapStart('OcnIO', 1)
    
    ! history
    call DOGCM_Admin_Variable_HistPut()
    call DOGCM_Boundary_vars_HistPut()
!!$    call LPhys_RediGM_spm_Output()
    call LPhys_RediGM_hspm_vfvm_Output()
    
    ! restart
    call DOGCM_IO_Restart_Output()
    call DOGCM_Admin_Variable_RestartPut()

    call ProfUtil_RapEnd('OcnIO', 1)
    
  end subroutine DOGCM_main_advance_timestep

  !-------------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_main_setup( configNmlFileArg )

    ! モジュール引用; Use statements
    !
    use OptionParser_mod, only: &
         & OptionParser_Init, OptionParser_Final, &
         & OptionParser_GetInfo

    use GridIndex_mod, only: &
         & GridIndex_Init
    
    use DOGCM_Admin_Grid_mod, only: &
         & iMax, jMax, kMax,           &
         & nMax, tMax,                 &
         & IA, JA,                     &
         & KS, KE, KA, z_CDK,          &
         & DOGCM_Admin_Grid_construct

    use SpmlUtil_mod, only: &
         & SpmlUtil_Init

    use VFvmUtil_mod, only: &
         & VFvmUtil_Init

    use CalculusDriver_mod, only: &
         & CalculusDriver_Init
    
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

    call GridIndex_Init( configNmlFile )
    
    call DOGCM_Admin_Constants_Init( configNmlFile )
    call DOGCM_Admin_GovernEq_Init( configNmlFile )  
    call DOGCM_Admin_Grid_Init( configNmlFile )


    ! Initialize some helper modules to solve oceanic flow
    !
    
#ifdef _OPENMP
    !$omp parallel
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    call MessageNotify('M', module_name, "Execute as Thread Parallel Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np=nThread)
#else
    call MessageNotify('M', module_name, "Execute as Serial  Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
#endif
    
    !-- Grid -------------------------------------------
    
    call DOGCM_Admin_Grid_construct()
    

    !-- Solver -----------------------------------------

    select case ( SolverType )
    case( OCNGOVERNEQ_SOLVER_HSPM_VSPM )
       call CalculusDriver_Init('hspm_vspm')
    case( OCNGOVERNEQ_SOLVER_HSPM_VFVM )
       call VFvmUtil_Init(KS, KE, KA, z_CDK, IA, JA)
       call CalculusDriver_Init('hspm_vfvm')       
    end select
    
    call DOGCM_Admin_TInteg_Init( configNmlFile )
    call DOGCM_Admin_BC_Init( configNmlFile )
    
    !-- IO ---------------------------------------------
    
    call DOGCM_IO_History_Init( configNmlFile )
    call DOGCM_IO_History_Create( configNmlFile )

    call DOGCM_IO_Restart_Init( configNmlFile )
    call DOGCM_IO_Restart_Create()
    
    !-- Variable (admin)  -----------------------------
    
    call DOGCM_Admin_Variable_Init()
    call DOGCM_Admin_Variable_regist_OuputVars()

    
    !-- Boundary -------------------------------------
    
    call DOGCM_Boundary_driver_Init( configNmlFile )        

    !-- TInt     -------------------------------------
    
    call DOGCM_TInt_driver_Init( configNmlFile )
    
    !-- Exp      -------------------------------------
    
    call DOGCM_Exp_driver_Init( configNmlFile )
    
  end subroutine DOGCM_main_setup

  
!-------------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_main_shutdown()

    ! モジュール引用; Use statements
    !
    use GridIndex_mod, only: &
         & GridIndex_Final

    use SpmlUtil_mod, only: &
         SpmlUtil_Final

    use VFvmUtil_mod, only: &
         & VFvmUtil_Final

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !

    call DOGCM_Exp_driver_Final()
    call DOGCM_Boundary_driver_Final()
    call DOGCM_TInt_driver_Final()
    
    call DOGCM_IO_History_Final()
    call DOGCM_IO_Restart_Final()

    call DOGCM_Admin_Variable_Final()
    call DOGCM_Admin_Grid_Final()
    call DOGCM_Admin_BC_Final()
    call DOGCM_Admin_TInteg_Final()
    call DOGCM_Admin_Constants_Final()

    !    
    call SpmlUtil_Final()
    select case ( SolverType )
    case( OCNGOVERNEQ_SOLVER_HSPM_VFVM )
       call VFvmUtil_Final()
    end select

    call GridIndex_Final()
    
  end subroutine DOGCM_main_shutdown

  !-------------------------------------------------------------------------------------------------
  
  subroutine DOGCM_main_update_SIceField( &
       & xy_SIceMaskRecv,                                    &   ! (in)
       & xy_SIceConRecv, xy_SfcHFlxNsRecv, xy_SfcHFlxSrRecv, &   ! (in)
       & xy_FreshWtFlxSRecv                                  &   ! (in)
       & )

    ! モジュール引用; Use statements
    !    
    use DOGCM_Admin_Grid_mod, only: &
         & IS, IE, JS, JE, &
         & KS, IA, JA

    use DOGCM_Boundary_vars_mod, only: &
         & xy_OcnSfcCellMask,                       &
         & OCNCELLMASK_OCEAN, OCNCELLMASK_SICE,     &
         & xy_SfcHFlxIO_ns, xy_SfcHFlxIO_sr,        &
         & xy_FreshWtFlxSIO, xy_FreshWtFlxIO

    ! 宣言文; Declaration statement
    !    
    logical, intent(in) :: xy_SIceMaskRecv(IA,JA)
    real(DP), intent(in) :: xy_SIceConRecv(IA,JA)
    real(DP), intent(in) :: xy_SfcHFlxNsRecv(IA,JA)    
    real(DP), intent(in) :: xy_SfcHFlxSrRecv(IA,JA)
    real(DP), intent(in) :: xy_FreshWtFlxSRecv(IA,JA)

    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp workshare
    xy_SfcHFlxIO_ns(:,:)  = xy_SfcHFlxNsRecv
    xy_SfcHFlxIO_sr(:,:)  = xy_SfcHFlxSrRecv
    xy_FreshWtFlxSIO(:,:) = xy_FreshWtFlxSRecv
    
    where ( xy_SIceMaskRecv )
       xy_OcnSfcCellMask(:,:) = OCNCELLMASK_SICE
    elsewhere
       xy_OcnSfcCellMask(:,:) = OCNCELLMASK_OCEAN
    end where
    !$omp end workshare
    !$omp end parallel
    
    
  end subroutine DOGCM_main_update_SIceField

  !- private subroutines ---------------------------------------------------------------------------

  
end module DOGCM_main_mod

