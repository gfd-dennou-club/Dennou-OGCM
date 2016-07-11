!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Boundary_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use SpmlUtil_mod
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Boundary_vars_mod, only: &
       & DOGCM_Boundary_vars_Init, DOGCM_Boundary_vars_Final

  use DOGCM_Boundary_common_mod, only: &
       & DOGCM_Boundary_common_Init, DOGCM_Boundary_common_Final, &
       & DOGCM_Boundary_common_UpdateBeforeTstep,                 &
       & DOGCM_Boundary_common_UpdateAfterTstep
  
  use DOGCM_Boundary_spm_mod, only: &
       & DOGCM_Boundary_spm_Init, DOGCM_Boundary_spm_Final, &
       & DOGCM_Boundary_spm_ApplyBC

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_Boundary_driver_Init, DOGCM_Boundary_driver_Final

  interface DOGCM_Boundary_driver_UpdateBeforeTstep
     module procedure DOGCM_Boundary_common_UpdateBeforeTstep
  end interface DOGCM_Boundary_driver_UpdateBeforeTstep
  public :: DOGCM_Boundary_driver_UpdateBeforeTstep
  
  interface DOGCM_Boundary_driver_UpdateAfterTstep
     module procedure DOGCM_Boundary_common_UpdateAfterTstep
  end interface DOGCM_Boundary_driver_UpdateAfterTstep
  public :: DOGCM_Boundary_driver_UpdateAfterTstep
  
  public :: DOGCM_Boundary_driver_ApplyBC
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Boundary_driver_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Boundary_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_vars_Init( configNmlName )
    call DOGCM_Boundary_common_Init( configNmlName )
    call DOGCM_Boundary_spm_Init( configNmlName )
!    call read_nmlData(configNmlName)

    
  end subroutine DOGCM_Boundary_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_Boundary_driver_Final()

    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_spm_Final()
    call DOGCM_Boundary_common_Final()
    call DOGCM_Boundary_vars_Final()
    
  end subroutine DOGCM_Boundary_driver_Final

  
  subroutine DOGCM_Boundary_driver_ApplyBC( &
       & xyz_U, xyz_V, xyza_TRC,                    & ! (inout)
       & xyz_H, xyz_VViscCoef, xyz_VDiffCoef        & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_U(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V(IA,JA,KA)
    real(DP), intent(inout) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)
    
    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_spm_ApplyBC(    &
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE), & ! (inout)
         & xyza_TRC(IS:IE,JS:JE,KS:KE,:),                      & ! (inout)
         & xyz_H(IS:IE,JS:JE,KS:KE),                                           & ! (in)
         & xyz_VViscCoef(IS:IE,JS:JE,KS:KE), xyz_VDiffCoef(IS:IE,JS:JE,KS:KE)  & ! (in)
         & )
    
  end subroutine DOGCM_Boundary_driver_ApplyBC
  
end module DOGCM_Boundary_driver_mod
