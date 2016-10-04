!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_GridDef_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開変数
  ! Public variables
  !

  integer, public, parameter :: XDIR = 1
  integer, public, parameter :: YDIR = 2
  integer, public, parameter :: ZDIR = 3

  integer, public, parameter :: I_XY = 1
  integer, public, parameter :: I_UY = 2
  integer, public, parameter :: I_XV = 3
  integer, public, parameter :: I_UV = 4
  
end module DOGCM_Admin_GridDef_mod
