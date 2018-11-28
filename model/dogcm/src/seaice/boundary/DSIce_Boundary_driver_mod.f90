!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_Boundary_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DSIce_Admin_Constants_mod, only: &
       & Mu
  
  use DSIce_Boundary_vars_mod, only: &
       & DSIce_Boundary_vars_Init, DSIce_Boundary_vars_Final
  
  use DSIce_Boundary_common_mod, only: &
       & DSIce_Boundary_common_Init, DSIce_Boundary_common_Final, &
       & DSIce_Boundary_common_UpdateBeforeTstep,                 &
       & DSIce_Boundary_common_UpdateAfterTstep

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_Boundary_driver_Init, DSIce_Boundary_driver_Final

  interface DSIce_Boundary_driver_UpdateBeforeTstep
     module procedure DSIce_Boundary_common_UpdateBeforeTstep
  end interface DSIce_Boundary_driver_UpdateBeforeTstep
  public :: DSIce_Boundary_driver_UpdateBeforeTstep
  
  interface DSIce_Boundary_driver_UpdateAfterTstep
     module procedure DSIce_Boundary_common_UpdateAfterTstep
  end interface DSIce_Boundary_driver_UpdateAfterTstep
  public :: DSIce_Boundary_driver_UpdateAfterTstep
  
  public :: DSIce_Boundary_driver_ApplyBC
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DSIce_Boundary_driver_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DSIce_Boundary_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call DSIce_Boundary_vars_Init( configNmlName )
    call DSIce_Boundary_common_Init( configNmlName )

!    call read_nmlData(configNmlName)

    
  end subroutine DSIce_Boundary_driver_Init

  !>
  !!
  !!
  subroutine DSIce_Boundary_driver_Final()

    ! 実行文; Executable statements
    !

    call DSIce_Boundary_vars_Final()
    call DSIce_Boundary_common_Final()
    
  end subroutine DSIce_Boundary_driver_Final

  !-----------------------------------------
  
  subroutine DSIce_Boundary_driver_UpdateBC( &
       & xy_SIceCon, xy_SnowThick, xy_IceThick,                      & ! (in)
       & xy_SIceSfcTemp, xyz_SIceTemp                                & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in)  :: xy_SIceCon(IA,JA)
    real(DP), intent(in)  :: xy_SnowThick(IA,JA)    
    real(DP), intent(in)  :: xy_IceThick(IA,JA)    
    real(DP), intent(in)  :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in)  :: xyz_SIceTemp(IA,JA,KA)
    
    
    ! 実行文; Executable statements
    !

    
  end subroutine DSIce_Boundary_driver_UpdateBC

  !-----------------------------------------
  
  subroutine DSIce_Boundary_driver_ApplyBC( &
       & )

    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statements
    !

    
  end subroutine DSIce_Boundary_driver_ApplyBC
  
end module DSIce_Boundary_driver_mod
