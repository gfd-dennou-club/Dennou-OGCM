!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
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
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use DSOGCM_main_mod, only: &
       & ogcm_main_Init => DSOGCM_main_Init, &
       & ogcm_main_Final => DSOGCM_main_Final, &
       & ogcm_setup => DSOGCM_main_setup, &
       & ogcm_shutdown => DSOGCM_main_shutdown, &
       & ogcm_advance_timestep => DSOGCM_main_advance_timestep
  
  ! 宣言文; Declaration statement
  !
  implicit none


  ! 局所変数
  ! Local variables
  !

  character(*), parameter :: PROGRAM_NAME = "ogcm_main"
  integer :: tstep
  logical :: loop_end_flag
  
  ! 実行文; Executable statement
  !

  call MessageNotify("M", PROGRAM_NAME, "Start..")
  call ogcm_main_Init()
  
  !***********************************
  ! Set up
  !***********************************

  call ogcm_setup()

  !***********************************
  ! The loop for temporal integration
  !************************************

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  tstep = 0; loop_end_flag = .false.
  do while(.not. loop_end_flag)
     call ogcm_advance_timestep(tstep, loop_end_flag)
     tstep = tstep + 1
  end do

  !*************************************
  ! Finalize
  !*************************************

  call ogcm_shutdown()
  
  call ogcm_main_Final()
  call MessageNotify("M", PROGRAM_NAME, "..End")


end program ogcm_main
