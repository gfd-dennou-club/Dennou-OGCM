!-------------------------------------------------------------
! Copyright (c) 2017-2017 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
program ocndiag_main

  ! モジュール引用; Use statement
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use dc_calendar
  
  use gtool_history
  
  !* Dennou-OGCM / SeaIce

  use OptionParser_mod, only: &
       & OptionParser_Init,   &
       & OptionParser_Final,  &
       & OptionParser_GetInfo
  
  use DOGCM_main_mod, only: &
       & ogcm_main_Init => DOGCM_main_Init,    &
       & ogcm_main_Final => DOGCM_main_Final

!!$  use GridIndex_mod, only: &
!!$       & IS, IE, JS, JE, KS, KE

  use DOGCM_Admin_Grid_mod, only: &
       & DOGCM_Admin_Grid_Init,        &
       & DOGCM_Admin_Grid_Final,       &
       & IA, IS, IE, IM,               &
       & JA, JS, JE, JM, JBLOCK,       &
       & KA, KS, KE, KM,               &
       & xyz_Z, xy_Topo,               &
       & xyz_Lat,                      &
       & z_CDK,                        &
       & DOGCM_Admin_Grid_UpdateVCoord

  use DSIce_Admin_Grid_mod, only: &
       & DSIce_Admin_Grid_Init,              &
       & DSIce_Admin_Grid_Final,             &
       & ISS=>IS, IES=>IE, JSS=>JS, JES=>JE, &
       & KSS=>KS, KES=>KE
       
  use DOGCM_Admin_Constants_mod
  use DOGCM_Admin_BC_mod
  use DOGCM_Admin_TInteg_mod
  use DOGCM_Admin_Variable_mod
  use DOGCM_Admin_GovernEq_mod
  use DOGCM_Admin_Variable_mod
  use DOGCM_IO_History_mod
  use DOGCM_Boundary_driver_mod
  use DOGCM_Boundary_vars_mod
  use DOGCM_TInt_driver_mod, only: &
       & DOGCM_TInt_driver_Init, DOGCM_TInt_driver_Final, &
       & CurrentTime

  use DSIce_Admin_Constants_mod
  use DSIce_Admin_GovernEq_mod
  use DSIce_Admin_Variable_mod, only: &
       & DSIce_Admin_Variable_Init,   &
       & xya_IceThick, xya_SnowThick, xya_SIceU, xya_SIceV, xyza_SIceEn
  
  use DOGCM_Exp_driver_mod
  
  use OcnDiag_DiagUtil_mod, only: &
       & DiagUtil_GetBVFreq,      &
       & DiagUtil_GetDensPot
  
  use OcnDiag_HeatTransport_mod, only: &
       & OcnDiag_HeatTransport_prepair_Output,      &
       & OcnDiag_HTDiagnose, OcnDiag_SIceHTDiagnose

  ! 宣言文; Declareration statements
  !  
  implicit none

  character(*), parameter :: PROGRAM_NAME = 'ocndiag'  
  character(STRING) :: configNmlFile
  character(STRING) :: configNmlFileDOGCM
  real(DP) :: TimeStart
  real(DP) :: TimeEnd
  real(DP) :: TimeInt
  real(DP) :: TimeCurrent
  character(TOKEN) :: TimeUnits

  type(gt_history) :: hst_diagvar
  
  call OptionParser_Init()
  call OptionParser_GetInfo( configNmlFile )
  call OptionParser_Final()
  
  call read_nmlData( configNmlFile )

  !-----------------------------------------------------------
  call MessageNotify('M', PROGRAM_NAME, "Setup..")
  call ogcm_main_Init()
  call ocn_setup( configNmlFileDOGCM, configNmlFile )

!!$  call DOGCM_Exp_driver_SetInitCond()
  call prepair_Output()
  
  !-----------------------------------------------------------  

  TimeCurrent = TimeStart
  do while(TimeCurrent <= TimeEnd)
     CurrentTime = DCCalConvertByUnit(TimeCurrent, TimeUnits, 'sec')
     
     call get_ncdata( TimeCurrent )

     ! Calculate xyz_Z in DOGCM_Admin_Grid_mod
     call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_A),   & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                                   & ! (in)
         & )

    call diagnose_OcnStaticStability()
    call diagnose_OcnHeatTransport()
    
     TimeCurrent = TimeCurrent + TimeInt
  end do
  
  !-----------------------------------------------------------  

  call MessageNotify('M', PROGRAM_NAME, "Sutdown..")
  call ocn_shutdown()
  call ogcm_main_Final()
  
contains
  
  subroutine diagnose_OcnStaticStability()

    real(DP) :: xyz_BVFreq(IA,JA,KA)
    real(DP) :: xyz_DensPot(IA,JA,KA)
    
    call DiagUtil_GetDensPot( xyz_DensPot,           &
         & xyzaa_TRC(:,:,:,TRCID_PTEMP,TIMELV_ID_N), &
         & xyzaa_TRC(:,:,:,TRCID_SALT,TIMELV_ID_N),  &
         & xyz_Z )

    call DiagUtil_GetBVFreq( xyz_BVFreq,             &
         & xyzaa_TRC(:,:,:,TRCID_PTEMP,TIMELV_ID_N), &
         & xyzaa_TRC(:,:,:,TRCID_SALT,TIMELV_ID_N),  &
         & xyza_H(:,:,:,TIMELV_ID_N), xyz_Z )
    
    call DOGCM_IO_History_HistPut( 'DensPot', xyz_DensPot(IS:IE,JS:JE,KS:KE) )
    call DOGCM_IO_History_HistPut( 'BVFreq', xyz_BVFreq(IS:IE,JS:JE,KS:KE) )
    
  end subroutine diagnose_OcnStaticStability

  subroutine diagnose_OcnHeatTransport()

    real(DP) :: y_EulerHT(JA)
    real(DP) :: y_BolusHT(JA)
    real(DP) :: y_IsoDiffHT(JA)
    real(DP) :: y_OcnHT(JA)
    real(DP) :: y_SIceHT(JA)
    real(DP) :: y_TotHT(JA)
    real(DP) :: y_NumDiffTend(JA)
    
    real(DP), parameter :: Unit2PWFac = 1d-15

    call OcnDiag_HTDiagnose( &
         & y_EulerHT, y_BolusHT, y_IsoDiffHT, y_NumDiffTend,                                  &  ! (out)
         & xyza_U(:,:,:,TIMELV_ID_N), xyza_V(:,:,:,TIMELV_ID_N), xyza_OMG(:,:,:,TIMELV_ID_N), &  ! (in)
         & xyzaa_TRC(:,:,:,TRCID_PTEMP,TIMELV_ID_N), xyzaa_TRC(:,:,:,TRCID_SALT,TIMELV_ID_N), &  ! (in)
         & xyza_H(:,:,:,TIMELV_ID_N), xyz_Z, xy_Topo )                                           ! (in)

    call OcnDiag_SIceHTDiagnose( &
         & y_SIceHT,                                                      & ! (out)
         & xya_SIceU(:,:,TIMELV_ID_N), xya_SIceV(:,:,TIMELV_ID_N),        & ! (in)
         & xya_IceThick(:,:,TIMELV_ID_N), xya_SnowThick(:,:,TIMELV_ID_N), & ! (in)
         & xyza_SIceEn(:,:,:,TIMELV_ID_N) )                                 ! (in)

    y_OcnHT(:) = y_EulerHT + y_BolusHT + y_IsoDiffHT        
    y_TotHT(:) = y_OcnHT + y_SIceHT
    call DOGCM_IO_History_HistPut( 'TotHT', y_TotHT(JS:JE)*Unit2PWFac )

  end subroutine diagnose_OcnHeatTransport
  
  subroutine get_ncdata( time )

    use gtool_history, only: HistoryGet
    use dc_string, only: toChar

    real(DP), intent(in) :: time
    character(STRING) :: trange
    character(TOKEN), parameter :: SIceHistoryName = 'history_sice.nc'
    
    trange = trim( 'time='//toChar(time) )
    call MessageNotify('M', PROGRAM_NAME, "get_ncdata..(%a [%a])",  &
         & ca=(/ trange, TimeUnits /) )
    
    call HistoryGet('U.nc', 'U', range=trange, quiet=.true.,           &
         & array=xyza_U(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) )
    call HistoryGet('V.nc', 'V', range=trange, quiet=.true.,           &
         & array=xyza_V(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) )
    call HistoryGet('OMG.nc', 'OMG', range=trange, quiet=.true.,           &
         & array=xyza_OMG(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) )
    call HistoryGet('PTemp.nc', 'PTemp', range=trange, quiet=.true.,   &
         & array=xyzaa_TRC(IS:IE,JS:JE,KS:KE,TRCID_PTEMP,TIMELV_ID_N) )
    call HistoryGet('Salt.nc', 'Salt', range=trange, quiet=.true.,     &
         & array=xyzaa_TRC(IS:IE,JS:JE,KS:KE,TRCID_Salt,TIMELV_ID_N))
    call HistoryGet('H.nc', 'H', range=trange, quiet=.true.,           &
         & array=xyza_H(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) )

    xya_SSH(:,:,TIMELV_ID_N) = 0d0
    xy_Topo(:,:)             = 5.2d3
    
    call HistoryGet(SIceHistoryName, 'IceThick', range=trange, quiet=.true.,       &
         & array=xya_IceThick(ISS:IES,JSS:JES,TIMELV_ID_N) )
    call HistoryGet(SIceHistoryName, 'SnowThick', range=trange, quiet=.true.,      &
         & array=xya_SnowThick(ISS:IES,JSS:JES,TIMELV_ID_N) )
    call HistoryGet(SIceHistoryName, 'SIceEn', range=trange, quiet=.true.,         &
         & array=xyza_SIceEn(ISS:IES,JSS:JES,KSS:KES,TIMELV_ID_N) )
    call HistoryGet(SIceHistoryName, 'SIceV', range=trange, quiet=.true.,          &
         & array=xya_SIceV(ISS:IES,JSS-1:JES,TIMELV_ID_N) )
    
  end subroutine get_ncdata

  !--------------------------------------------------------------------------------

  subroutine prepair_Output()

    call DOGCM_IO_History_RegistVar( 'DensPot', 'IJKT', 'potential density', 'kg/m3' )
    call DOGCM_IO_History_RegistVar( 'BVFreq', 'IJKT', 'Brunt-Vaisala frequency', 's-1' )

    call OcnDiag_HeatTransport_prepair_Output()
    call DOGCM_IO_History_RegistVar( 'TotHT', 'JT', 'total heat transport(ocean and seaice)', 'PW' )
    
  end subroutine prepair_Output
  
  !--------------------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine ocn_setup( configNmlFileOGCM, configNmlFile )

    ! モジュール引用; Use statements
    !
    use GridIndex_mod, only: &
         & GridIndex_Init
    
    use DOGCM_Admin_Grid_mod, only: &
         & iMax, jMax, kMax,           &
         & nMax, tMax,                 &
         & IA, JA,                     &
         & KS, KE, KA, z_CDK,          &
         & DOGCM_Admin_Grid_construct

    use DSIce_Admin_Grid_mod, only: &
         & DSIce_Admin_Grid_construct
    
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
    character(*), intent(in) :: configNmlFileOGCM
    character(*), intent(in) :: configNmlFile

    ! 局所変数
    ! Local variables
    !
    integer :: nThread

    ! 実行文; Executable statement
    !
    
    ! Initialize administrative modules 
    !

    call GridIndex_Init( configNmlFileOGCM )
    
    call DOGCM_Admin_Constants_Init( configNmlFileOGCM )
    call DSIce_Admin_Constants_Init( configNmlFileOGCM )

    call DOGCM_Admin_GovernEq_Init( configNmlFileOGCM )  
    call DSIce_Admin_GovernEq_Init( configNmlFileOGCM )  

    call DOGCM_Admin_Grid_Init( configNmlFileOGCM )
    call DSIce_Admin_Grid_Init( configNmlFileOGCM )

    ! Initialize some helper modules to solve oceanic flow
    !
    
#ifdef _OPENMP
    !$omp parallel
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    call MessageNotify('M', PROGRAM_NAME, "Execute as Thread Parallel Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np=nThread)
#else
    call MessageNotify('M', PROGRAM_NAME, "Execute as Serial  Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
#endif
    
    !-- Grid -------------------------------------------
    
    call DOGCM_Admin_Grid_construct()
    call DSIce_Admin_Grid_construct()

    !-- Solver -----------------------------------------

    select case ( SolverType )
    case( OCNGOVERNEQ_SOLVER_HSPM_VSPM )
       call CalculusDriver_Init('hspm_vspm')
    case( OCNGOVERNEQ_SOLVER_HSPM_VFVM )
       call VFvmUtil_Init(KS, KE, KA, z_CDK, IA, JA)
       call CalculusDriver_Init('hspm_vfvm')       
    end select
    
    call DOGCM_Admin_TInteg_Init( configNmlFileOGCM )
    call DOGCM_Admin_BC_Init( configNmlFileOGCM )
    
    !-- IO ---------------------------------------------
    
    call DOGCM_IO_History_Init( configNmlFile )
    call DOGCM_IO_History_Create( configNmlFile )
!!$    call DOGCM_IO_Restart_Init( configNmlFileOGCM )
!!$    call DOGCM_IO_Restart_Create()
    
    !-- Variable (admin)  -----------------------------
    
    call DOGCM_Admin_Variable_Init()
!!$    call DOGCM_Admin_Variable_regist_OuputVars()
    call DSIce_Admin_Variable_Init()

    
    !-- Boundary -------------------------------------
    
    call DOGCM_Boundary_driver_Init( configNmlFileOGCM )        

    !-- TInt     -------------------------------------
    
    call DOGCM_TInt_driver_Init( configNmlFileOGCM )
    
    !-- Exp      -------------------------------------
    
    call DOGCM_Exp_driver_Init( configNmlFileOGCM )
    
  end subroutine ocn_setup
  
!-------------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine ocn_shutdown()

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
!!$    call DOGCM_IO_Restart_Final()

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
    
  end subroutine ocn_shutdown
  
  !--------------------------------------------------------------------------------
  
  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    !
    use dc_string, only: Split, Replace, StrInclude

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /ocndiag_nml/ &
         & configNmlFileDOGCM,                     &
         & TimeStart, TimeEnd, TimeInt, TimeUnits


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    configNmlFileDOGCM = ''
    TimeStart          = 0d0
    TimeEnd            = 0d0
    TimeInt            = 1d0
    TimeUnits          = 'day'
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', PROGRAM_NAME, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                                         ! (in)
            & nml = ocndiag_nml, iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    if (len(trim(configNmlFileDOGCM)) == 0) then
       call MessageNotify( 'E', PROGRAM_NAME, "Specify the value of 'configNmlFileDOGCM'.")
    end if
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, '  - configNmlFileDOGCM      = %a', ca=(/ trim(configNmlFileDOGCM) /)) 
    call MessageNotify( 'M', PROGRAM_NAME, '  - TimeStart=%f, TimeEnd=%f, TimeInt=%f [%c]', &
         &                   d=(/ TimeStart, TimeEnd, TimeInt /), ca=(/ trim(TimeUnits) /)  )
    
  end subroutine read_nmlData
  
end program ocndiag_main
