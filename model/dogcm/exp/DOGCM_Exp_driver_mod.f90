!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Exp_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod, only: &
       & UNDEFVAL
  
  use DOGCM_Admin_TInteg_mod, only: &
       & TIMELV_ID_N

  use DOGCM_Admin_Constants_mod, only: &
       & vViscCoef, vDiffCoef

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM,       &
       & DOGCM_Admin_Grid_UpdateVCoord
  
  use DOGCM_Admin_Variable_mod, only: &
       & xya_SSH, xyza_H, xyza_U, xyza_V, xyza_OMG, xyzaa_TRC,    &
       & xyza_HydPres, xya_SfcPres,                               &
       & xyz_VViscCoef, xyz_VDiffCoef, xyz_ConvIndex,             &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Boundary_vars_mod, only: &
       & xy_SfcHFlx0_ns, xy_SfcHFlx0_sr,   &
       & xy_FreshWtFlx0, xy_FreshWtFlxS0,  &
       & xy_WindStressU, xy_WindStressV, &
       & xy_SeaSfcTemp, xy_SeaSfcSalt

  !* Dennou-SIce
  
  use DSIce_Admin_Variable_mod, only: &
       & xya_IceThick, xya_SnowThick,    &
       & xya_SIceSfcTemp, xyza_SIceTemp, &
       & xya_SIceCon

  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI,               &
       & xy_DSfcHFlxAIDTs
  
!------------------------------------------------  
!!$  use DOGCM_Exp_BarotRossby_mod, only:       &
!!$  use DOGCM_Exp_IGW_mod, only:               &
!!$  use DOGCM_Exp_WindDrivCirc_mod, only:      &
!!$  use DOGCM_Exp_EqJetAccel_mod, only:        &
!!$  use DOGCM_Exp_APEOGCirc_mod, only:            &
!!$  use DOGCM_Exp_APEOGCircSIce_mod, only:        &
  use DOGCM_Exp_APECoupleClimate_mod, only:        &
!!$  use DOGCM_Exp_APEOGCircI98BC_mod, only:        &
       & DOGCM_Exp_Init,                           &
       & DOGCM_Exp_Final,                          &
       & DOGCM_Exp_SetInitCond,                    &
       & DOGCM_Exp_Do
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_Exp_driver_Init, DOGCM_Exp_driver_Final

  public :: DOGCM_Exp_driver_SetInitCond
  public :: DOGCM_Exp_driver_Do
  
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Exp_driver_mod' !< Module Name

  real(DP), parameter :: DEFAULT_SSH         = 0d0    
  real(DP), parameter :: DEFAULT_VEL         = 0d0  
  real(DP), parameter :: DEFAULT_PTEMP       = 280d0
  real(DP), parameter :: DEFAULT_SALT        = 35d0
  real(DP), parameter :: DEFAULT_SfcPresInit = 0d0
  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Exp_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)
    call DOGCM_Exp_Init( configNmlName )

  end subroutine DOGCM_Exp_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_Exp_driver_Final()

    ! 実行文; Executable statements
    !
    
    call DOGCM_Exp_Final()


  end subroutine DOGCM_Exp_driver_Final

  !------------------------------------
  
  subroutine DOGCM_Exp_driver_SetInitCond()

    ! 実行文; Executable statements
    !

    !- Set default values of OGCM variables ----------

    xya_SSH(:,:,:)  = DEFAULT_SSH
    xyza_U(:,:,:,:) = DEFAULT_VEL
    xyza_V(:,:,:,:) = DEFAULT_VEL
    xyzaa_TRC(:,:,:,TRCID_PTEMP,:) = DEFAULT_PTEMP
    xyzaa_TRC(:,:,:,TRCID_SALT,:)  = DEFAULT_SALT
    xya_SfcPres(:,:,:) = DEFAULT_SfcPresInit
    
    xyz_VViscCoef(:,:,:) = vViscCoef
    xyz_VDiffCoef(:,:,:) = vDiffCoef
    xyz_ConvIndex(:,:,:) = 0d0
    
    xy_SfcHFlx0_sr(:,:)  = 0d0
    xy_SfcHFlx0_ns(:,:)  = 0d0
    xy_FreshWtFlx0(:,:)  = 0d0
    xy_FreshWtFlxS0(:,:) = 0d0
    xy_WindStressU(:,:) = 0d0
    xy_WindStressV(:,:) = 0d0
    xy_SeaSfcTemp(:,:) = DEFAULT_PTEMP
    xy_SeaSfcSalt(:,:) = DEFAULT_SALT

    !- Set default values of sea-ice variables ----------

    xya_IceThick(:,:,:) = 0d0
    xya_SnowThick(:,:,:) = 0d0
    xya_SIceSfcTemp(:,:,:) = UNDEFVAL
    xyza_SIceTemp(:,:,:,:) = UNDEFVAL
    xya_SIceCon(:,:,:) = 0d0

    xy_SfcHFlxAI(:,:) = 0d0
    xy_DSfcHFlxAIDTs(:,:) = 0d0
    
    !- Call subroutine to set intial values ---
    
    call DOGCM_Exp_SetInitCond()

    !------------------------------------------
    
    call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_N), & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                                & ! (in)
         & )
    
  end subroutine DOGCM_Exp_driver_SetInitCond
  
  !------------------------------------

  subroutine DOGCM_Exp_driver_Do()

    call DOGCM_Exp_Do()
    
  end subroutine DOGCM_Exp_driver_Do

  !------------------------------------

end module DOGCM_Exp_driver_mod
