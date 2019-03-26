!-------------------------------------------------------------
! Copyright (c) 2013-2019 Yuta Kawai. All rights reserved.
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

  use ProfUtil_mod, only: &
       & ProfUtil_Init, ProfUtil_Final,        &
       & ProfUtil_RapStart, ProfUtil_RapEnd,   &
       & ProfUtil_RapReport
  
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

  !---------------------------------------------------------------------------------------------
  
  !  実行文; Executable statement
  !

  call MessageNotify("M", PROGRAM_NAME, "Start..")

  call OptionParser_Init()
  call OptionParser_GetInfo( configNmlFile )
  call OptionParser_Final()
  
  call read_nmlData( configNmlFile )

  !*********************************************************************************************
  ! Set up
  !*********************************************************************************************  
  
  call ProfUtil_Init( configNmlFile )
  call ProfUtil_RapStart('Setup', 0) 
  
  call ogcm_main_Init()
  call sice_main_Init()
  
  call ogcm_setup( configNmlFile )
  call sice_setup( configNmlFile )

  !- Set initial condition ----------------------------------------------------------------  
  
  call DOGCM_Exp_driver_SetInitCond()

  call ProfUtil_RapEnd('Setup', 0) 
  
  !*********************************************************************************************
  ! The loop for temporal integration
  !*********************************************************************************************

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")
  
  call ProfUtil_RapStart('TimeLoop', 0) 
  
  loop_end_flag = .false.  
  tstep_ocn = 0
  tstep_sice = 0
  
  !- Time loop ---------------------------------------------------------------------------
  do while(.not. loop_end_flag)
          
     !* Sea ice component ************
!!$     call MessageNotify( 'M', PROGRAM_NAME, "SIce component tstep=%d", i=(/ tstep_sice /)) 
     if (SICE_do) then
       call sice_advance_timestep(tstep_sice, loop_end_flag_sice, skip_flag=.false.)
       if (OCN_do) call pass_field_sice2ocn()
       tstep_sice = tstep_sice + 1
     else
       loop_end_flag_sice = .true.
     end if

     !* Ocean component ************
!!$     call MessageNotify( 'M', PROGRAM_NAME, "OCN component tstep=%d", i=(/ tstep_ocn /))
     if (OCN_do) then
       call ogcm_advance_timestep(tstep_ocn, loop_end_flag_ocn, skip_flag=.false.)
       if (SICE_do) call pass_field_ocn2sice()
       tstep_ocn = tstep_ocn + 1
     else
          loop_end_flag_ocn = .true.   
     end if
     
     !* Call a subroutine defined by users.
     call DOGCM_Exp_driver_Do()

     !

     !
     loop_end_flag = ( loop_end_flag_ocn .and. loop_end_flag_sice )     
  end do

  call ProfUtil_RapEnd('TimeLoop', 0) 
  
  !*********************************************************************************************
  ! Finalize
  !*********************************************************************************************

  call ProfUtil_RapStart('Shutdown', 0)
  
  call ogcm_shutdown()
  call sice_shutdown()

  call ogcm_main_Final()
  call sice_main_Final()

  call ProfUtil_RapEnd('Shutdown', 0)

  call ProfUtil_RapReport()
  call ProfUtil_Final()
   
  !*************************************

  call MessageNotify("M", PROGRAM_NAME, "..End")
  
contains

  !-----------------------------------------------------------
  
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

    use DSIce_Admin_Variable_mod, only: &
         & xya_SIceCon, xya_IceThick, xya_SnowThick, &
         & xya_SIceSfcTemp, xyza_SIceTemp

    use DSIce_Boundary_driver_mod, only: &
         & DSIce_Boundary_driver_UpdateAfterTstep
    
!!$    use DOGCM_Phys_hspm_vfvm_mod, only: &
!!$         & reconstruct_sfctemp
    
!!$    use SpmlUtil_mod, only: w_xy, xy_w
!!$    use LPhys_DIFF_spm_mod, only: &
!!$         & w_Filter
!!$    
    real(DP) :: xy_SfcTemp(IA,JA)
    
    ! 実行文; Executable statement
    !

!!$    xy_SfcTemp(:,:) = xyzaa_TRC(:,:,KS,TRCID_PTEMP,TIMELV_ID_N)
!!$    call reconstruct_sfctemp( xy_SfcTemp )    
!!$    xy_SfcTemp(IS:IE,JS:JE) = xy_w(w_Filter*w_xy(xyzaa_TRC(IS:IE,JS:JE,KS,TRCID_PTEMP,TIMELV_ID_N)))
    
    call DSIce_main_update_OcnField( &
         & xyzaa_TRC(:,:,KS,TRCID_PTEMP,TIMELV_ID_N),             & ! (in)
!!$         & xy_SfcTemp, & ! (in)
         & xyzaa_TRC(:,:,KS,TRCID_SALT,TIMELV_ID_N),              & ! (in)
         & z_KAXIS_Weight(KS)*xyza_H(:,:,KS,TIMELV_ID_N),         & ! (in)
         & xyza_U(:,:,KS,TIMELV_ID_N), xyza_V(:,:,KS,TIMELV_ID_N) & ! (in)
         & )

    if (tstep_sice == 1) then
       !- Call a subroutine in sea-ice model again
       ! in order to set surface albedo of sea-ice (and ocean) grid
       call DSIce_Boundary_driver_UpdateAfterTstep( &
            & xya_SIceCon(:,:,TIMELV_ID_N), xya_IceThick(:,:,TIMELV_ID_N), xya_SnowThick(:,:,TIMELV_ID_N), & ! (in)
            & xya_SIceSfcTemp(:,:,TIMELV_ID_N), xyza_SIceTemp(:,:,:,TIMELV_ID_N)                           & ! (in)
            & )
    end if
    
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

  !-----------------------------------------------------------
  
  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !

    use dc_iounit, only: FileOpen
    use dc_types, only: STDOUT 
    use dc_string, only: Split, Replace, StrInclude

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml
    ! Unit number for NAMELIST file open

    integer:: iostat_nml 

    ! NAMELIST group name
    !
    namelist /dogcm_nml/ &
         & OCN_do, SIce_do


    ! 実行文; Executable statement
    !

    ! デフォルト値の設定
    ! Default values settings
    !

    OCN_do  = .true.
    SICE_do = .true.
    
    
    ! NAMELIST から入力
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
    
    

    ! 出力 ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, '< DOGCM components             >')
    call MessageNotify( 'M', PROGRAM_NAME, '  - ocean         = %b', L = (/ OCN_do /)) 
    call MessageNotify( 'M', PROGRAM_NAME, '  - sea ice       = %b', L = (/ SICE_do /)) 

  end subroutine read_nmlData  

end program ogcm_main
