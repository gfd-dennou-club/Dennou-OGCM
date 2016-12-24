!-------------------------------------------------------------
! Copyright (c) 2013-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief Main program of Dennou-OGCM. 
!!
!! @author Yuta Kawai
!! @since 2013
!!
!!
program ogcm_main

  ! モジュール引用; Use statement
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  use OptionParser_mod, only: &
       & OptionParser_Init,   &
       & OptionParser_Final,  &
       & OptionParser_GetInfo
  
  use DOGCM_main_mod, only: &
       & ogcm_main_Init => DOGCM_main_Init,    &
       & ogcm_main_Final => DOGCM_main_Final,  &
       & ogcm_setup => DOGCM_main_setup,       &
       & ogcm_shutdown => DOGCM_main_shutdown, &
       & ogcm_advance_timestep => DOGCM_main_advance_timestep

  use DOGCM_Exp_driver_mod, only: &
       & DOGCM_Exp_driver_SetInitCond, &
       & DOGCM_Exp_driver_Do
  
  use DSIce_main_mod, only: &
       & sice_main_Init => DSIce_main_Init,    &
       & sice_main_Final => DSIce_main_Final,  &
       & sice_setup => DSIce_main_setup,       &
       & sice_shutdown => DSIce_main_shutdown, &
       & sice_advance_timestep => DSIce_main_advance_timestep

  ! 宣言文; Declaration statement
  !
  implicit none


  ! 局所変数
  ! Local variables
  !

  character(*), parameter :: PROGRAM_NAME = "ogcm_main"

  logical :: OCN_do
  logical :: SICE_do
  
  integer :: tstep_ocn
  integer :: tstep_sice
  
  logical :: loop_end_flag
  logical :: loop_end_flag_ocn
  logical :: loop_end_flag_sice

  character(STRING) :: configNmlFile

  !-------------------------------------------------------------
  
  ! 実行文; Executable statement
  !

  call MessageNotify("M", PROGRAM_NAME, "Start..")

  call OptionParser_Init()
  call OptionParser_GetInfo( configNmlFile )
  call OptionParser_Final()

  call read_nmlData( configNmlFile )
  
  call ogcm_main_Init()
  call sice_main_Init()
  
  !***********************************
  ! Set up
  !***********************************
  
  call ogcm_setup( configNmlFile )
  call sice_setup( configNmlFile )

  !- Set initial condition ----------------------------------------------------------------  

  call DOGCM_Exp_driver_SetInitCond()

  
  !***********************************
  ! The loop for temporal integration
  !************************************

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  loop_end_flag = .false.  

  tstep_ocn = 0
  tstep_sice = 0

  !- Time loop ---------------------------------------------------------------------------
  do while(.not. loop_end_flag)
          
     !* Sea ice component ************
!!$     call MessageNotify( 'M', PROGRAM_NAME, "SIce component tstep=%d", i=(/ tstep_sice /)) 
     if (SICE_do) call sice_advance_timestep(tstep_sice, loop_end_flag_sice, skip_flag=.false.)
     call pass_field_sice2ocn()
     tstep_sice = tstep_sice + 1
     
     !* Ocean component ************
!!$     call MessageNotify( 'M', PROGRAM_NAME, "OCN component tstep=%d", i=(/ tstep_ocn /))
     if (OCN_do) call ogcm_advance_timestep(tstep_ocn, loop_end_flag_ocn, skip_flag=.false.)
     call pass_field_ocn2sice()
     tstep_ocn = tstep_ocn + 1

     !* Call a subroutine defined by users.
     call DOGCM_Exp_driver_Do()
     
     !

     !
     loop_end_flag = ( loop_end_flag_ocn .and. loop_end_flag_sice )     
  end do

  !*************************************
  ! Finalize
  !*************************************

  call ogcm_shutdown()
  call sice_shutdown()
  
  call ogcm_main_Final()
  call sice_main_Final()
  
  call MessageNotify("M", PROGRAM_NAME, "..End")


contains

  subroutine pass_field_ocn2sice()

    ! モジュール引用; Use statements
    !
    
    use DSIce_main_mod, only: &
         & DSIce_main_update_OcnField

    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N

    use DOGCM_Admin_Grid_mod, only: &
         & IS, IE, JS, JE, &
         & KS, IA, JA,     &
         & z_KAXIS_Weight, xy_Topo
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyza_U, xyza_V,                     &
         & xyzaa_TRC, TRCID_PTEMP, TRCID_SALT, &
         & xyza_H

    ! 実行文; Executable statement
    !

    call DSIce_main_update_OcnField( &
         & xyzaa_TRC(:,:,KS,TRCID_PTEMP,TIMELV_ID_N),             & ! (in)
         & xyzaa_TRC(:,:,KS,TRCID_SALT,TIMELV_ID_N),              & ! (in)
         & z_KAXIS_Weight(KS)*xyza_H(:,:,KS,TIMELV_ID_N),         & ! (in)
         & xyza_U(:,:,KS,TIMELV_ID_N), xyza_V(:,:,KS,TIMELV_ID_N) & ! (in)
         & )
    
  end subroutine pass_field_ocn2sice

  !-----------------------------------------------------------
  
  subroutine pass_field_sice2ocn()

    ! モジュール引用; Use statements
    !
    
    use DOGCM_main_mod, only: &
         & DOGCM_main_update_SIceField

    use DSIce_Admin_Constants_mod, only: &
         & IceMaskMin
    
    use DSIce_Admin_TInteg_mod, only: &
         & TIMELV_ID_B, TIMELV_ID_N

    use DSIce_Admin_Grid_mod, only: &
         & IS, IE, JS, JE, &
         & KS, IA, JA
    
    use DSIce_Admin_Variable_mod, only: &
         & xya_SIceCon

    use DSIce_Boundary_vars_mod, only: &
         & xy_WindStressUAI, xy_WindStressVAI, &
         & xy_WindStressUIO, xy_WindStressVIO, &
         & xy_BtmHFlxIO, xy_FreshWtFlxS

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xy_BtmHFlxIO_sr(IA,JA)
    
    ! 実行文; Executable statement
    !

    xy_BtmHFlxIO_sr(:,:) = 0d0
    call DOGCM_main_update_SIceField( &
         & (xya_SIceCon(:,:,TIMELV_ID_B) >= IceMaskMin), & ! (in)
         & xya_SIceCon(:,:,TIMELV_ID_B),                 & ! (in)
         & xy_BtmHFlxIO, xy_BtmHFlxIO_sr,                & ! (in)
         & xy_FreshWtFlxS                                & ! (in)
         & )

  end subroutine pass_field_sice2ocn

  !-----------------------------------------------------------

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
    namelist /dogcm_nml/ &
         & OCN_do, SIce_do


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    OCN_do  = .true.
    SICE_do = .true.
    
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', PROGRAM_NAME, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                                         ! (in)
            & nml = dogcm_nml, iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    ! - Convert the type name into the corresponding ID ---------
    !

    ! Specify the governing equations used in thermodynamics model
    
    

    ! 印字 ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, '< DOGCM components             >')
    call MessageNotify( 'M', PROGRAM_NAME, '  - ocean         = %b', L = (/ OCN_do /)) 
    call MessageNotify( 'M', PROGRAM_NAME, '  - sea ice       = %b', L = (/ SICE_do /)) 

  end subroutine read_nmlData
  
end program ogcm_main
