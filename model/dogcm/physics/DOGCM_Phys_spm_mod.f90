!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Phys_spm_mod

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
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax, nMax

  use DOGCM_Admin_Constants_mod, only: &
       & PI, RPlanet,               &
       & hViscCoef, hHyperViscCoef, &
       & hDiffCoef, hHyperDiffCoef
       
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM,        &
       & TRCID_PTEMP, TRCID_SALT

  use SpmlUtil_mod

  use DOGCM_Admin_GovernEq_mod, only: &
       & isPhysicsCompActivated, &
       & OCNGOVERNEQ_LPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_LPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_LPHYS_REDIGM_NAME, &
       & OCNGOVERNEQ_VPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_VPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_VPHYS_CONVEC_NAME

  !-- Module for the parametrizations of lateral oceanic physics

  use LPhys_DIFF_spm_mod, only: &
       & LPhys_DIFF_spm_Init, LPhys_DIFF_spm_Final, &
       & LPhys_DIFF_spm_LMixMOMRHS,                       &
       & LPhys_DIFF_spm_LMixMOMRHSImpl,                   &
       & LPhys_DIFF_spm_LMixTRCRHS,                       &
       & LPhys_DIFF_spm_LMixTRCRHSImpl
       
  use LPhys_RediGM_spm_mod, only: &
       & LPhys_RediGM_spm_Init, LPhys_RediGM_spm_Final, &
       & LPhys_RediGM_spm_AddMixingTerm

  !-- Module for the parametrizations of vertical oceanic physics
  
  use DOGCM_VPhys_ConvAdjust_mod, only: &
       & DOGCM_VPhys_ConvAdjust_Init, DOGCM_VPhys_ConvAdjust_Final, &
       & DOGCM_VPhys_ConvAdjust_AddMixingTerm
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_Phys_spm_Init, DOGCM_Phys_spm_Final

  public :: DOGCM_Phys_spm_Do
  
  public :: DOGCM_Phys_spm_VMixMOMRHS
  public :: DOGCM_Phys_spm_VMixTRCRHS

  public :: DOGCM_Phys_spm_VImplUV
  public :: DOGCM_Phys_spm_VImplTRC
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Phys_spm_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Phys_spm_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !
!    call read_nmlData(configNmlName)

    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then
       call LPhys_DIFF_spm_Init( configNmlName = configNmlName )
    end if
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_spm_Init( confignmlFileName = configNmlName )
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Init()
    end if
    
  end subroutine DOGCM_Phys_spm_Init

  !>
  !!
  !!
  subroutine DOGCM_Phys_spm_Final()

    ! 実行文; Executable statements
    !
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_spm_Final()
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Final()
    end if
    
  end subroutine DOGCM_Phys_spm_Final

  !-----------------------------------------

  subroutine DOGCM_Phys_spm_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,        & ! (inout)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt                                                    & ! (in)
       & )

    use DOGCM_Admin_Grid_mod, only: &
         & z_KAXIS_Weight
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyz_ConvIndex

    use SpmlUtil_mod, only: &
         & AvrLonLat_xy, xy_IntSig_BtmToTop_xyz
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS_phy(0:iMax-1,jMax,0:kMax)
    real(DP), intent(inout) :: xyz_V_RHS_phy(0:iMax-1,jMax,0:kMax)
    real(DP), intent(inout) :: xyza_TRC_RHS_phy(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(inout) :: xyz_VViscCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(inout) :: xyz_VDiffCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(inout) :: xy_BtmFrictCoef(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_TOPO(0:iMax-1,jMax)
    real(DP), intent(in) :: dt


    real(DP) :: avr_ptemp_RHS_phys

    real(DP) :: xyz_UA(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_VA(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyza_TRCA(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    
    
    ! 実行文; Executable statements
    !

    !-- Horizontal momentum -----------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXMOM_NAME )  ) then
!!$       call LPhys_DIFF_spm_LMixMOMRHS( &
!!$            & xyz_U_RHS_phy, xyz_V_RHS_phy,                              & ! (inout)
!!$            & xyz_U, xyz_V, xyz_H, hViscCoef, hHyperViscCoef             & ! (in)
!!$            & )
       call LPhys_DIFF_spm_LMixMOMRHSImpl( &
            & xyz_U_RHS_phy, xyz_V_RHS_phy,                              & ! (inout)
            & xyz_U, xyz_V, xyz_H, hViscCoef, hHyperViscCoef,            & ! (in)
            & dt                                                         & ! (in)
            & )

    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXMOM_NAME )  ) then
       call DOGCM_Phys_spm_VMixMOMRHS(            & 
            & xyz_U_RHS_phy, xyz_V_RHS_phy,                              & ! (inout)
            & xyz_U, xyz_V, xyz_H, xyz_VViscCoef                         & ! (in)
            & )
    end if

    !-- Tracer ------------------------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then
!!$       call LPhys_DIFF_spm_LMixTRCRHS(            &
!!$            & xyza_TRC_RHS_phy,                                         & ! (inout)
!!$            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef                & ! (in)
!!$            & )
       
       call LPhys_DIFF_spm_LMixTRCRHSImpl(            &
            & xyza_TRC_RHS_phy,                                         & ! (inout)
            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef,               & ! (in)
            & dt )                                                        ! (in)

!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "->avr_ptemp_phys (+ LMixRHS): ", avr_ptemp_RHS_phys
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
       call DOGCM_Phys_spm_VMixTRCRHS(            &  
            & xyza_TRC_RHS_phy,                                         & ! (inout)
            & xyza_TRC, xyz_H, xyz_VDiffCoef                            & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ VMixRHS): ",  avr_ptemp_RHS_phys
!!$       write(*,*) "avr_ptemp_phys_local:", IntSig_BtmToTop(xyza_TRC_RHS_phy(0,jMax/2,:,TRCID_SALT))
    end if

    !-- Lateral mixing of tracers by eddy induced velocity -------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_spm_AddMixingTerm( &
            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP),                      & ! (inout)
            & xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),                       & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),  & ! (in)
            & xyz_H, xyz_Z, xy_Topo                                     & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ RediGMRHS): ",  avr_ptemp_RHS_phys
    end if

    !-- Vertical mixing of tracers by non-penetrative convection -------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_AddMixingTerm( &
            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP),                      & ! (inout)
            & xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),                       & ! (inout)
            & xyz_ConvIndex(IS:IE,JS:JE,KS:KE),                         & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),  & ! (in)
            & xyz_Z, z_KAXIS_Weight(KS:KE), dt                          & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ ConvAdjustRHS): ",  avr_ptemp_RHS_phys       
    end if

  end subroutine DOGCM_Phys_spm_Do

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VMixMOMRHS(     &
       & xyz_U_RHS, xyz_V_RHS,                            & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_VViscCoef               & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(inout) :: xyz_V_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VViscCoef(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !    
    integer :: n

    ! 実行文; Executable statements
    !

    xyz_U_RHS(:,:,:) = xyz_U_RHS + &
         & xyz_DSig_xyz( xyz_VViscCoef/xyz_H * xyz_DSig_xyz(xyz_U) )/xyz_H

    xyz_V_RHS(:,:,:) = xyz_V_RHS + &
         & xyz_DSig_xyz( xyz_VViscCoef/xyz_H * xyz_DSig_xyz(xyz_V) )/xyz_H
    
    
  end subroutine DOGCM_Phys_spm_VMixMOMRHS

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VMixTRCRHS(     &
       & xyza_TRC_RHS,                      & ! (out)
       & xyza_TRC, xyz_H, xyz_VDiffCoef     & ! (in)
       & )

    use SpmlUtil_mod, only: &
         DMat1 => tr_vDeriv1CoefMat
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyza_TRC_RHS(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VDiffCoef(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !    
    integer :: n

    ! 実行文; Executable statements
    !

    do n = 1, TRC_TOT_NUM
       xyza_TRC_RHS(:,:,:,n) = xyza_TRC_RHS(:,:,:,n) + &
            & xyz_DSig_xyz( xyz_VDiffCoef/xyz_H * xyz_DSig_xyz(xyza_TRC(:,:,:,n)) )/xyz_H
    end do
    
  end subroutine DOGCM_Phys_spm_VMixTRCRHS

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VImplTRC( xyza_TRCA,         & ! (out)
       & xyza_TRC0, xyza_HTRC_RHS,                       & ! (in)
       & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, alpha        & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & ThermBC_Surface, ThermBC_Bottom, &
         & SaltBC_Surface, SaltBC_Bottom
    
    use DOGCM_Boundary_spm_mod, only: &
         & DOGCM_Boundary_spm_InqVBCRHS_TRC

    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyza_TRCA(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC0(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)    
    real(DP), intent(in) :: xyza_HTRC_RHS(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VDiffCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xya_PTemp_VBCRHS(0:iMax-1,jMax,2)
    real(DP) :: xya_Salt_VBCRHS(0:iMax-1,jMax,2)
    real(DP) :: xyza_DHTRC(0:iMax-1,jMax,TRC_TOT_NUM)
    integer :: n
    integer :: k
    
    real(DP) :: avr_RHS
    real(DP) :: avr_TRC0_phys
    real(DP) :: avr_TRCA_phys
    real(DP) :: xyz_DzSalt(0:iMax-1,jMax,0:kMax)

    
    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_spm_InqVBCRHS_TRC( &
         & xya_PTemp_VBCRHS, xya_Salt_VBCRHS,                                & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP), xyza_TRC0(:,:,:,TRCID_SALT),        & ! (in)
         & xyz_H0, xyz_VDiffCoef                                             & ! (in)
         & )

    call calc_VDiffEq( xyza_TRCA(:,:,:,TRCID_PTEMP),                        & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_PTEMP)/xyz_H0,                         & ! (i1n)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                                & ! (in)
         & inquire_VBCSpecType(ThermBC_Surface), xya_PTemp_VBCRHS(:,:,1),   & ! (in)
         & inquire_VBCSpecType(ThermBC_Bottom),  xya_PTemp_VBCRHS(:,:,2),   & ! (in)
         & "PTemp"                                                          & ! (in)
         & )

    
    call calc_VDiffEq( xyza_TRCA(:,:,:,TRCID_SALT),                        & ! (out)
         & xyza_TRC0(:,:,:,TRCID_SALT),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_SALT)/xyz_H0,                         & ! (in)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                               & ! (in)
         & inquire_VBCSpecType(SaltBC_Surface), xya_Salt_VBCRHS(:,:,1),    & ! (in)
         & inquire_VBCSpecType(SaltBC_Bottom),  xya_Salt_VBCRHS(:,:,2),    & ! (in)
         & "Salt"                                                          & ! (in)
         & )

    
    do n = 1, TRC_TOT_NUM
       xyza_TRCA(:,:,:,n) = xyz_H0(:,:,:)*xyza_TRCA(:,:,:,n)/xyz_HA(:,:,:)
    end do
    
    
!!$    xyz_DzSalt = xyz_DSig_xyz(xyza_TRCA(:,:,:,TRCID_SALT))
!!$    avr_RHS = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_HTRC_RHS(:,:,:,TRCID_SALT)/xyz_H0 ))
!!$    avr_TRC0_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC0(:,:,:, TRCID_SALT) ))    
!!$    avr_TRCA_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRCA(:,:,:, TRCID_SALT) ))
!!$    write(*,*) "avr_dsalt (+ VIMPLTRC): ",  (avr_TRCA_phys - avr_TRC0_phys)/dt*dt, avr_RHS*dt 
    
  end subroutine DOGCM_Phys_spm_VImplTRC

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,          & ! (in)
       & xyz_H, xyz_VViscCoef, xy_BtmFrictCoef,         & ! (in)
       & dt, alpha                                      & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & DynBC_Surface, DynBC_Bottom
    
    use DOGCM_Boundary_spm_mod, only: &
         & DOGCM_Boundary_spm_InqVBCRHS_UV
    
    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyz_UA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_VA(0:iMax-1,jMax,0:kMax)        
    real(DP), intent(in) :: xyz_U0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V0(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_U_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VViscCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_BtmFrictCoef(0:iMax-1,jMax)    
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xya_U_VBCRHS(0:iMax-1,jMax,2)
    real(DP) :: xya_V_VBCRHS(0:iMax-1,jMax,2)

    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)    
    integer :: k

    ! 実行文; Executable statements
    !
    
    call DOGCM_Boundary_spm_InqVBCRHS_UV( &
         & xya_U_VBCRHS, xya_V_VBCRHS,                                     & ! (out)
         & xyz_U0, xyz_V0, xyz_H, xyz_VViscCoef                            & ! (in)
         & )
    

    call calc_VDiffEq( xyz_UA,                                             & ! (out)
         & xyz_U0, xyz_U_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
         & inquire_VBCSpecType(DynBC_Surface), xya_U_VBCRHS(:,:,1),        & ! (in)
         & inquire_VBCSpecType(DynBC_Bottom),  xya_U_VBCRHS(:,:,2),        & ! (in)
         & "U"                                                             & ! (in)
         & )

    call calc_VDiffEq( xyz_VA,                                             & ! (out)
         & xyz_V0, xyz_V_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
         & inquire_VBCSpecType(DynBC_Surface), xya_V_VBCRHS(:,:,1),        & ! (in)
         & inquire_VBCSpecType(DynBC_Bottom),  xya_V_VBCRHS(:,:,2),        & ! (in)
         & "V"                                                             & ! (in)
         & )

    
  end subroutine DOGCM_Phys_spm_VImplUV

  !-----------------------------------------------------------------------

  subroutine calc_VDiffEq( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_RHS,                     &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal,          &  ! (in)
       & qname                                  &  ! (in)
       & )

    use SpmlUtil_mod, only: &
         DMat1 => tr_vDeriv1CoefMat, &
         IntWt => g_Sig_WEIGHT
    
    ! 宣言文; Declaration statement
    !      
    real(DP), intent(out) :: xyz_qA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_q0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_q_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Kv(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    character, intent(in) :: SeaSfcBCType
    real(DP), intent(in) :: xy_SeaSfcBCVal(0:iMax-1,jMax)
    character, intent(in) :: SeaBtmBCType
    real(DP), intent(in) :: xy_SeaBtmBCVal(0:iMax-1,jMax)
    character(*), intent(in) :: qname

    ! 局所変数
    ! Local variables
    !
    real(DP) :: AMat(0:kMax,0:kMax)
    real(DP) :: z_Kv(0:KMax)
    real(DP) :: z_e3(0:kMax)

    real(DP) :: DSigMat(0:kMax,0:kMax)
    real(DP) :: IMat(0:kMax,0:kMax)
    real(DP) :: b(0:kMax)
    integer :: IPiv(0:kMax)
    integer :: info

    integer :: i
    integer :: j
    integer :: k
    integer :: l

    integer :: N
    integer :: Ns, Ne
   
   
    ! 実行文; Executable statements
    !

    N = size(IMat,1)
    IMat(:,:) = 0d0
    forAll(k=0:kMax) IMat(k,k) = 1d0
    do k = 0, kMax
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do

    DSigMat(:,:) = transpose(DMat1)

    !$omp parallel do collapse(2) private(z_Kv, z_e3, AMat, b, IPiv, info, i, k, l)
    do j = 1, jMax
       do i = 0, iMax-1

          ! Solve a linear equation system, 
          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 

          ! Set coefficient matrix A 
          !

          z_Kv(:) = xyz_Kv(i,j,:)
          z_e3(:) = xyz_H(i,j,:)
          
          do k=0, kMax
                AMat(k,:) =  IMat(k,:) &
                     &      -  alpha*dt/z_e3(k) * matmul( DSigMat(k,:) * z_Kv/z_e3, DSigMat(:,:) )
                
                b(k) = dt * xyz_q_RHS(i,j,k)
          end do

          ! Set boundary condition
          !

          select case(SeaSfcBCType)
          case('D')
             AMat(0,:) = IMat(0,:)
             b(0) = xy_SeaSfcBCVal(i,j) - xyz_q0(i,j,0)
          case('N')
                AMat(0,:) = DSigMat(0,:)
                b(0) = xy_SeaSfcBCVal(i,j) - sum(DSigMat(0,:)*xyz_q0(i,j,:))
          end select

          select case(SeaBtmBCType)
          case('D')
             AMat(kMax,:) = IMat(kMax,:)
             b(kMax) = xy_SeaBtmBCVal(i,j) - xyz_q0(i,j,kMax)
          case('N')
                AMat(kMax,:) = DSigMat(kMax,:)
                b(kMax) = xy_SeaBtmBCVal(i,j) - sum(DSigMat(kMax,:)*xyz_q0(i,j,:))
          end select
          
          call DGESV(N, 1, AMat, N, IPIV, b, N, info)

          xyz_qA(i,j,:) = xyz_q0(i,j,:) + b(:)

       end do
    end do

  end subroutine calc_VDiffEq
    
  !-----------------------------------------------------------------------

  subroutine calc_VDiffEq_v2( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_RHS,                     &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal,          &  ! (in)
       & qname                                  &  ! (in)
       & )

    use SpmlUtil_mod, only: &
         DMat1 => tr_vDeriv1CoefMat, &
         IntWt => g_Sig_WEIGHT
    
    ! 宣言文; Declaration statement
    !      
    real(DP), intent(out) :: xyz_qA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_q0(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_q_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Kv(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    character, intent(in) :: SeaSfcBCType
    real(DP), intent(in) :: xy_SeaSfcBCVal(0:iMax-1,jMax)
    character, intent(in) :: SeaBtmBCType
    real(DP), intent(in) :: xy_SeaBtmBCVal(0:iMax-1,jMax)
    character(*), intent(in) :: qname

    ! 局所変数
    ! Local variables
    !
    real(DP) :: AMat(0:kMax,0:kMax)
    real(DP) :: z_Kv(0:KMax)
    real(DP) :: z_e3(0:kMax)

    real(DP) :: DSigMat(0:kMax,0:kMax)
    real(DP) :: IMat(0:kMax,0:kMax)
    real(DP) :: b(0:kMax)
    integer :: IPiv(0:kMax)
    integer :: info

    integer :: i
    integer :: j
    integer :: k
    integer :: l

    integer :: N
    integer :: Ns, Ne

    real(DP) :: w(1:kMax+1)
    real(DP) :: rhow(1:kMax+1)
    real(DP) :: cb1stPt(1:kMax-1)
    real(DP) :: cglPt(1:kMax+1)
    real(DP) :: PMat(kMax-1,kMax+1)
    real(DP) :: PDSigMat(kMax-1,kMax+1)
    
    ! 実行文; Executable statements
    !

    N = size(IMat,1)

    rhow(:) = 2d0
    rhow(2:kMax) = 1d0
    w(1) = 1d0
    do k=2,kMax+1
       w(k) = - w(k-1)
    end do
    w(:) = w(:)/rhow(:)
    do k=1, N
       cglPt(k) = cos((k-1)*PI/dble(N-1))
    end do
    
    !---------
    
    do k=1, kMax-1
       cb1stPt(k) = cos((2*k-1)*PI/dble(2d0*(N-2)))
       PMat(k,:) = w(:) / (cb1stPt(k) - cglPt(:)) &
            &      / sum(w(:)/(cb1stPt(k) - cglPt(:)))
    end do

!!$    b(:) = cos(cglPt(:)*PI)
!!$    do k=1,N-2
!!$       write(*,*) "k=", k, "cb1stPT:", cb1stPt(k), "cglPT=", cglPt(k), ":", cos(cb1stPt(k)*PI), &
!!$            & sum(PMat(k,:)*b(:))
!!$    end do
!!$    write(*,*) sum( (cos(cb1stPt*PI) - matmul(PMat,b))**2 * PI/dble(kMax-1+1) )
!!$    stop
    
    !---------
    
    IMat(:,:) = 0d0
    forAll(k=0:kMax) IMat(k,k) = 1d0
    do k = 0, kMax
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do

    DSigMat(:,:) = transpose(DMat1)
    PDSigMat(:,:) = matmul(PMat,DSigMat)
    
    !$omp parallel do collapse(2) private(z_Kv, z_e3, AMat, b, IPiv, info, i, k, l)
    do j = 1, jMax
       do i = 0, iMax-1

          ! Solve a linear equation system, 
          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 

          ! Set coefficient matrix A 
          !

          z_Kv(:) = xyz_Kv(i,j,:)
          z_e3(:) = xyz_H(i,j,:)
          
          do k=0, kMax
                AMat(k,:) =  IMat(k,:) &
                     &      -  alpha*dt/z_e3(k) * matmul( DSigMat(k,:) * z_Kv/z_e3, DSigMat(:,:) )
                
                b(k) = dt * xyz_q_RHS(i,j,k)
          end do
          AMat(1:kMax-1,:) = matmul(PMat, AMat)
          b(1:kMax-1) = matmul(PMat,b)

          ! Set boundary condition
          !

          select case(SeaSfcBCType)
          case('D')
             AMat(0,:) = IMat(0,:)
             b(0) = xy_SeaSfcBCVal(i,j) - xyz_q0(i,j,0)
          case('N')
             AMat(0,:) = DSigMat(0,:)
             b(0) = xy_SeaSfcBCVal(i,j) - sum(DSigMat(0,:)*xyz_q0(i,j,:))
          end select

          select case(SeaBtmBCType)
          case('D')
             AMat(kMax,:) = IMat(kMax,:)
             b(kMax) = xy_SeaBtmBCVal(i,j) - xyz_q0(i,j,kMax)
          case('N')
             AMat(kMax,:) = DSigMat(kMax,:)
             b(kMax) = xy_SeaBtmBCVal(i,j) - sum(DSigMat(kMax,:)*xyz_q0(i,j,:))
          end select
          
          call DGESV(N, 1, AMat, N, IPIV, b, N, info)

          xyz_qA(i,j,:) = xyz_q0(i,j,:) + b(:)

       end do
    end do

  end subroutine calc_VDiffEq_v2
    
  
!!$  subroutine calc_VDiffEq_2( xyz_qA,              &  ! (inout)
!!$       & xyz_q0, xyz_q_RHS,                     &  ! (in)
!!$       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
!!$       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
!!$       & SeaBtmBCType, xy_SeaBtmBCVal,          &  ! (in)
!!$       & qname                                  &  ! (in)
!!$       & )
!!$
!!$    use SpmlUtil_mod, only: &
!!$         DMat1 => tr_vDeriv1CoefMat, &
!!$         IntWt => g_Sig_WEIGHT
!!$    
!!$    ! 宣言文; Declaration statement
!!$    !      
!!$    real(DP), intent(out) :: xyz_qA(0:iMax-1,jMax,0:kMax)
!!$    real(DP), intent(in) :: xyz_q0(0:iMax-1,jMax,0:kMax)
!!$    real(DP), intent(in) :: xyz_q_RHS(0:iMax-1,jMax,0:kMax)
!!$    real(DP), intent(in) :: xyz_Kv(0:iMax-1,jMax,0:kMax)
!!$    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
!!$    real(DP), intent(in) :: dt
!!$    real(DP), intent(in) :: alpha
!!$    character, intent(in) :: SeaSfcBCType
!!$    real(DP), intent(in) :: xy_SeaSfcBCVal(0:iMax-1,jMax)
!!$    character, intent(in) :: SeaBtmBCType
!!$    real(DP), intent(in) :: xy_SeaBtmBCVal(0:iMax-1,jMax)
!!$    character(*), intent(in) :: qname
!!$
!!$    ! 局所変数
!!$    ! Local variables
!!$    !
!!$    real(DP) :: AMat(0:kMax,0:kMax)
!!$    real(DP) :: z_Kv(0:KMax)
!!$    real(DP) :: z_e3(0:kMax)
!!$
!!$    real(DP) :: DSigMat(1:kMax-1,1:kMax-1)
!!$    real(DP) :: basisInv(2:kMax, 1:kMax-1)
!!$    real(DP) :: basis(1:kMax-1, 2:kMax)
!!$    real(DP) :: basis_x(1:kMax-1, 2:kMax)
!!$    real(DP) :: theta(0:kMax)
!!$    real(DP) :: Tcb(0:kMax,0:kMax)
!!$    real(DP) :: Tcb_x(0:kMax,0:kMax)
!!$    integer :: Nbasis
!!$    real(DP) :: workInv((kMax-1)*64)
!!$    real(DP) :: c(2:kMax)
!!$    
!!$    real(DP) :: IMat(0:kMax,0:kMax)
!!$    real(DP) :: b(0:kMax)
!!$    integer :: IPiv(0:kMax)
!!$    integer :: info
!!$
!!$    integer :: i
!!$    integer :: j
!!$    integer :: k
!!$    integer :: l
!!$
!!$    integer :: N
!!$    integer :: Ns, Ne
!!$
!!$    real(DP) :: sig(0:kMax)
!!$    real(DP) :: xy_B0(0:iMax-1,jMax)
!!$    real(DP) :: xy_B1(0:iMax-1,jMax)
!!$    real(DP) :: z_B(0:kMax)
!!$
!!$    real(DP) :: t_q(0:kMax)
!!$    real(DP) :: t_v(2:kMax)
!!$    
!!$    ! 実行文; Executable statements
!!$    !
!!$
!!$    N = kMax + 1
!!$    Nbasis = N - 2
!!$    
!!$    !---
!!$    IMat(:,:) = 0d0
!!$    forAll(k=0:kMax) IMat(k,k) = 1d0
!!$
!!$    do k=0, kMax
!!$       theta(k) = PI*k/dble(N-1)
!!$       sig(k) = cos(theta(k))
!!$       do l=0, kMax
!!$          Tcb(k,l) = cos(l*theta(k))
!!$          Tcb_x(k,l) = l*sin(l*theta(k))/sin(theta(k))
!!$       end do
!!$       Tcb_x(:,0) = 0d0
!!$       Tcb_x(:,1) = 1d0
!!$    end do
!!$
!!$    select case(SeaSfcBCType//SeaBtmBCType)
!!$    case('NN')
!!$       xy_B0 = 0d0
!!$       xy_B1 = 0d0
!!$       do l=0, kMax-2
!!$             c(l+2) = - (l/dble(l+2))**2
!!$             basis(1:kMax-1,l+2) = Tcb(1:kMax-1,l) + c(l+2)*Tcb(1:kMax-1,l+2)
!!$             basis_x(1:kMax-1,l+2) = Tcb_x(1:kMax-1,l) + c(l+2)*Tcb_x(1:kMax-1,l+2)
!!$       end do
!!$    case('ND')
!!$    case('DD')
!!$       do l=2, kMax
!!$          if (mod(l,2) == 0d0) then
!!$             basis(1:kMax-1,l) = Tcb(1:kMax-1,l) - Tcb(1:kMax-1,0)
!!$             basis_x(1:kMax-1,l) = Tcb_x(1:kMax-1,l)
!!$          else
!!$             basis(1:kMax-1,l) = Tcb(1:kMax-1,l) - Tcb(1:kMax-1,1)
!!$             basis_x(1:kMax-1,l) = Tcb_x(1:kMax-1,l) - 1d0          
!!$          end if
!!$       end do
!!$    end select
!!$
!!$    basisInv(:,:) = basis
!!$    call DGETRF(Nbasis, Nbasis, basisInv, Nbasis, IPiv(1:kMax-1), info)
!!$    call DGETRI(Nbasis, basisInv, Nbasis, IPiv(1:kMax-1), workInv, size(workInv), info)
!!$
!!$    !
!!$    DSigMat(:,:) = matmul(basis_x(1:kMax-1,2:kMax), basisInv(:,:))
!!$    
!!$    !$omp parallel do collapse(2) private(z_Kv, z_e3, AMat, b, IPiv, info, i, k, l, z_B, t_q, t_v)
!!$    do j = 1, jMax
!!$       do i = 0, iMax-1
!!$
!!$          ! Solve a linear equation system, 
!!$          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 
!!$
!!$          ! Set coefficient matrix A 
!!$          !
!!$
!!$          z_Kv(:) = xyz_Kv(i,j,:)
!!$          z_e3(:) = xyz_H(i,j,:)
!!$          z_B(:) = xy_B0(i,j) + sig(:)*xy_B1(i,j)
!!$          b(:) = xyz_q0(i,j,:)
!!$          
!!$          do k=1, kMax-1
!!$                AMat(k,1:kMax-1) =  IMat(k,1:kMax-1) &
!!$                     &      -  dt/z_e3(k) * matmul( DSigMat(k,:) * z_Kv(1:kMax-1)/z_e3(1:kMax-1), DSigMat(:,:) )
!!$                
!!$                b(k) = b(k) + dt * ( &
!!$                     &    xyz_q_RHS(i,j,k)                                                          &
!!$                     & -  1d0/z_e3(k)*sum(DSigMat(k,:)*z_Kv(1:kMax-1)/z_e3(1:kMax-1)*xy_B1(i,j))    &
!!$                     & -  z_B(k) )
!!$          end do
!!$          
!!$          call DGESV(Nbasis, 1, AMat(1:kMax-1,1:kMax-1), Nbasis, IPIV, b(1:kMax-1), Nbasis, info)
!!$          t_v(2:kMax) = matmul(basisInv, b(1:kMax-1))
!!$          
!!$          !
!!$          select case(SeaSfcBCType//SeaBtmBCType)
!!$          case('NN')
!!$             t_q(0:2) = t_v(2:4)
!!$             do l=3, kMax-2
!!$                t_q(l) = c(l)*t_v(l) + t_v(l+2)
!!$             end do
!!$             t_q(kMax-1:kMax) = c(kMax-1:kMax)*t_v(kMax-1:kMax)
!!$          case('DD')
!!$             t_q(2:kMax) = t_v(:)
!!$             t_q(0) = - sum(t_v(2:kMax:2))
!!$             t_q(1) = - sum(t_v(3:kMax:2))
!!$          end select
!!$
!!$          xyz_qA(i,j,0) = sum(Tcb(0,:)*t_q(:))
!!$          xyz_qA(i,j,kMax) = sum(Tcb(kMax,:)*t_q(:))
!!$          xyz_qA(i,j,1:kMax-1) = b(1:kMax-1) 
!!$          if (j==jMax/2) then
!!$             write(*,*) "j=",j, "qA:", xyz_qA(i,j,:)
!!$          end if
!!$          
!!$       end do
!!$    end do
!!$
!!$  end subroutine calc_VDiffEq_2
  
end module DOGCM_Phys_spm_mod
