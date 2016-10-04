!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Dyn_driver_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/SIce

  use UnitConversion_mod, only: &
       & degC2K, K2degC
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, JA, KA, &
       & x_IAXIS_Weight, y_JAXIS_Weight

  use DSIce_Dyn_fvm_mod, only: &
       & DSIce_Dyn_fvm_Init, DSIce_Dyn_fvm_Final, &
       & DSIce_Dyn_fvm_SIceThickDiffRHS
  
    
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Dyn_driver_Init, DSIce_Dyn_driver_Final
  public :: DSIce_Dyn_driver_ADVRHS
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Dyn_driver_mod' !< Module Name


  logical :: initedFlag = .false.
  
contains

  !>
  !!
  !!
  subroutine DSIce_Dyn_driver_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName
    
    ! 実行文; Executable statements
    !

    call DSIce_Dyn_fvm_Init( configNmlName )
    
    initedFlag = .true.
    
  end subroutine DSIce_Dyn_driver_Init

  !>
  !!
  !!
  subroutine DSIce_Dyn_driver_Final()

    ! 実行文; Executable statements
    !

    call DSIce_Dyn_fvm_Final()
    
  end subroutine DSIce_Dyn_driver_Final

  !---------------------------------------------------------------
  
  subroutine DSIce_Dyn_driver_ADVRHS(                           & 
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS,     & ! (out)
       & xyz_SIceEn_RHS,                                        & ! (out)
       & xy_SIceU, xy_SIceV,                                    & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,              & ! (in)
       & xyz_SIceEn0                                            & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceCon_RHS(IA,JA)
    real(DP), intent(out) :: xy_IceThick_RHS(IA,JA)
    real(DP), intent(out) :: xy_SnowThick_RHS(IA,JA)
    real(DP), intent(out) :: xyz_SIceEn_RHS(IA,JA,KA)
    real(DP), intent(out) :: xy_SIceU(IA,JA)
    real(DP), intent(out) :: xy_SIceV(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xyz_SIceEn0(IA,JA,KA)

    !
    !
    call DSIce_Dyn_fvm_SIceThickDiffRHS(           & 
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS, & ! (out)
       & xyz_SIceEn_RHS,                                    & ! (out)
       & xy_SIceU, xy_SIceV,                                & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,          & ! (in)
       & xyz_SIceEn0                                        & ! (in)
       & )
    
  end subroutine DSIce_Dyn_driver_ADVRHS


end module DSIce_Dyn_driver_mod

