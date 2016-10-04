!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_TInt_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SIce

  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM

  use TemporalIntegUtil_mod2, only: &
       & TimeIntMode_Euler, &
       & TimeIntMode_RK2,   &
       & TimeIntMode_LF,    &
       & TimeIntMode_LFAM3
  
  use DSIce_Admin_TInteg_mod, only: &
       & SIceTimeIntMode,                       &
       & TIMELV_ID_N, TIMELV_ID_B, TIMELV_ID_A, &
       & DelTime
  
  use DSIce_Admin_Variable_mod, only: &
       & xya_SIceCon, xya_IceThick, xya_SnowThick,   &
       & xya_SIceSfcTemp, xyza_SIceEn, xyza_SIceTemp, &
       & xy_Wice
  
  use DSIce_Dyn_driver_mod, only: &
       & DSIce_Dyn_driver_Init, DSIce_Dyn_driver_Final
       
!!$  use DSIce_Phys_driver_mod, only: &
!!$       & DSIce_Phys_driver_Init, DSIce_Phys_driver_Final

  use DSIce_TInt_common_mod, only: &
       & DSIce_TInt_common_Init, DSIce_TInt_common_Final, &
       & DSIce_TInt_common_advance_Dyn,                   &
       & DSIce_TInt_common_advance_ThermoDyn,             &
       & DSIce_TInt_common_advance_Phys
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_TInt_driver_mod' !< Module Name


  public :: DSIce_TInt_driver_Init, DSIce_TInt_driver_Final
  public :: DSIce_TInt_driver_Do


contains

  !>
  !!
  !!
  Subroutine DSIce_TInt_driver_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)

    call DSIce_Dyn_driver_Init( configNmlName )
    
!!$    call DSIce_Phys_driver_Init( configNmlName )

    call DSIce_TInt_common_Init( configNmlName )
    
  end subroutine DSIce_TInt_driver_Init

  !>
  !!
  !!
  subroutine DSIce_TInt_driver_Final()

    ! 実行文; Executable statements
    !

    call DSIce_TInt_common_Final()
    
    call DSIce_Dyn_driver_Final()
!!$    call DSIce_Phys_driver_Final()
    
  end subroutine DSIce_TInt_driver_Final

  !-------------------------------------

  
  subroutine DSIce_TInt_driver_Do( isSelfStartSchemeUsed )

    ! モジュール引用; Use statements
    !    

    ! 宣言文; Declaration statement
    !
    logical, intent(in) :: isSelfStartSchemeUsed

    ! 作業変数
    ! Work variables
    !
    
    ! 実行文; Executable statements
    !

    !-- Advance time step ---------------------------------------------------------
    
    if ( isSelfStartSchemeUsed ) then
          call DSIce_TInt_Euler_Do()
    else

       select case (SIceTimeIntMode)
       case (TimeIntMode_Euler)
          call DSIce_TInt_Euler_Do()
       case default
          call MessageNotify( 'E', module_name, &
               & "Unexcepted temporal scheme is specified. Check!" )
       end select
       
    end if
    
  end subroutine DSIce_TInt_driver_Do

  !-----------------------------------------

  subroutine DSIce_TInt_Euler_Do()

    ! 作業変数
    ! Work variables
    !
    integer :: TA
    integer :: TN

    real(DP) :: xy_SIceCon_thm(IA,JA)
    real(DP) :: xy_IceThick_thm(IA,JA)
    real(DP) :: xy_SnowThick_thm(IA,JA)
    real(DP) :: xyz_SIceEn_thm(IA,JA,KA)
    
    ! 実行文; Executable statement
    !

    TA = TIMELV_ID_A
    TN = TIMELV_ID_N
    
!    call MessageNotify('M', module_name, "TInt = Euler ..")

    xy_Wice(:,:) = 0d0

    call DSIce_TInt_common_advance_ThermoDyn(  &
         !----- Tendency due to thermodynamics process ------------------------
         & xy_SIceCon_thm, xy_IceThick_thm, xy_SnowThick_thm,                        & ! (out)
         & xyz_SIceEn_thm, xya_SIceSfcTemp(:,:,TA),                                  & ! (out) 
         & xy_Wice,                                                                  & ! (out)
         !----- Time level 0 --------------------------------------------------
         & xya_SIceCon(:,:,TN), xya_IceThick(:,:,TN), xya_SnowThick(:,:,TN),         & ! (in)
         & xya_SIceSfcTemp(:,:,TN), xyza_SIceTemp(:,:,:,TN), xyza_SIceEn(:,:,:,TN),  & ! (in)
         !---------------------------------------------------------------------
         & DelTime                                                                   & ! (in)
         & )

    call DSIce_TInt_common_advance_Dyn(  &
         !----- Time level A --------------------------------------------------
         & xya_SIceCon(:,:,TA), xya_IceThick(:,:,TA), xya_SnowThick(:,:,TA),         & ! (out)
         & xya_SIceSfcTemp(:,:,TA),                                                  & ! (inout)
         & xyza_SIceTemp(:,:,:,TA), xyza_SIceEn(:,:,:,TA),                           & ! (out)
         & xy_Wice,                                                                  & ! (inout)
         !----- Time level 0 --------------------------------------------------
         & xya_SIceCon(:,:,TN), xya_IceThick(:,:,TN), xya_SnowThick(:,:,TN),         & ! (out)
         & xya_SIceSfcTemp(:,:,TN), xyza_SIceTemp(:,:,:,TN), xyza_SIceEn(:,:,:,TN),  & ! (out)
         !---------------------------------------------------------------------
         & xy_SIceCon_thm, xy_IceThick_thm, xy_SnowThick_thm, xyz_SIceEn_thm,        & ! (in) 
         & DelTime                                                                   & ! (in)
         & )

  end subroutine DSIce_TInt_Euler_Do

  !-----------------------------------------------
  
end module DSIce_TInt_driver_mod
