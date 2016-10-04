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

  use DOGCM_LPhys_DIFF_spm_mod, only: &
       & DOGCM_LPhys_DIFF_spm_Init, DOGCM_LPhys_DIFF_spm_Final, &
       & DOGCM_LPhys_DIFF_spm_LMixMOMRHS,                       &
       & DOGCM_LPhys_DIFF_spm_LMixMOMRHSImpl,                   &
       & DOGCM_LPhys_DIFF_spm_LMixTRCRHS,                       &
       & DOGCM_LPhys_DIFF_spm_LMixTRCRHSImpl
       
  use DOGCM_LPhys_RediGM_spm_mod, only: &
       & DOGCM_LPhys_RediGM_spm_Init, DOGCM_LPhys_RediGM_spm_Final, &
       & DOGCM_LPhys_RediGM_spm_AddMixingTerm

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
       call DOGCM_LPhys_DIFF_spm_Init( configNmlName = configNmlName )
    end if
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call DOGCM_LPhys_RediGM_spm_Init( confignmlFileName = configNmlName )
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
       call DOGCM_LPhys_RediGM_spm_Final()
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Final()
    end if
    
  end subroutine DOGCM_Phys_spm_Final

  !-----------------------------------------

  subroutine DOGCM_Phys_spm_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef,                         & ! (in)
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
    real(DP), intent(in) :: xyz_VViscCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VDiffCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,0:kMax,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_TOPO(0:iMax-1,jMax)
    real(DP), intent(in) :: dt


    real(DP) :: avr_ptemp_RHS_phys
    
    ! 実行文; Executable statements
    !

    !-- Horizontal momentum -----------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXMOM_NAME )  ) then
!!$       call DOGCM_LPhys_DIFF_spm_LMixMOMRHS( &
!!$            & xyz_U_RHS_phy, xyz_V_RHS_phy,                              & ! (inout)
!!$            & xyz_U, xyz_V, xyz_H, hViscCoef, hHyperViscCoef             & ! (in)
!!$            & )
       call DOGCM_LPhys_DIFF_spm_LMixMOMRHSImpl( &
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
!!$       call DOGCM_LPhys_DIFF_spm_LMixTRCRHS(            &
!!$            & xyza_TRC_RHS_phy,                                         & ! (inout)
!!$            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef                & ! (in)
!!$            & )

       call DOGCM_LPhys_DIFF_spm_LMixTRCRHSImpl(            &
            & xyza_TRC_RHS_phy,                                         & ! (inout)
            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef,               & ! (in)
            & dt )                                                        ! (in)

!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_PTEMP) ))
!!$       write(*,*) "avr_ptemp_phys (+ LMixRHS): ",  avr_ptemp_RHS_phys
    end if


    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
       call DOGCM_Phys_spm_VMixTRCRHS(            &  
            & xyza_TRC_RHS_phy,                                         & ! (inout)
            & xyza_TRC, xyz_H, xyz_VDiffCoef                            & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_PTEMP) ))
!!$       write(*,*) "avr_ptemp_phys (+ VMixRHS): ",  avr_ptemp_RHS_phys
    end if

    !-- Lateral mixing of tracers by eddy induced velocity -------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call DOGCM_LPhys_RediGM_spm_AddMixingTerm( &
            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP),                      & ! (inout)
            & xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),                       & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),  & ! (in)
            & xyz_H, xyz_Z, xy_Topo                                     & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_PTEMP) ))
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
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_PTEMP) ))
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
         & xyz_DSig_xyz( xyz_VViscCoef*xyz_DSig_xyz(xyz_U)/xyz_H )/xyz_H

    xyz_V_RHS(:,:,:) = xyz_V_RHS + &
         & xyz_DSig_xyz( xyz_VViscCoef*xyz_DSig_xyz(xyz_V)/xyz_H )/xyz_H
    
    
  end subroutine DOGCM_Phys_spm_VMixMOMRHS

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VMixTRCRHS(     &
       & xyza_TRC_RHS,                      & ! (out)
       & xyza_TRC, xyz_H, xyz_VDiffCoef     & ! (in)
       & )
    
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
            & xyz_DSig_xyz( xyz_VDiffCoef*xyz_DSig_xyz(xyza_TRC(:,:,:,n))/xyz_H )/xyz_H
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

    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_spm_InqVBCRHS_TRC( &
         & xya_PTemp_VBCRHS, xya_Salt_VBCRHS,                                & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP), xyza_TRC0(:,:,:,TRCID_SALT),        & ! (in)
         & xyz_H0, xyz_VDiffCoef                                             & ! (in)
         & )

    call calc_VDiffEq( xyza_TRCA(:,:,:,TRCID_PTEMP),                        & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_PTEMP)/xyz_H0,                         & ! (in)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                                & ! (in)
         & inquire_VBCSpecType(ThermBC_Surface), xya_PTemp_VBCRHS(:,:,1),   & ! (in)
         & inquire_VBCSpecType(ThermBC_Bottom),  xya_PTemp_VBCRHS(:,:,2)    & ! (in)
         & )

    call calc_VDiffEq( xyza_TRCA(:,:,:,TRCID_SALT),                        & ! (out)
         & xyza_TRC0(:,:,:,TRCID_SALT),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_SALT)/xyz_H0,                         & ! (in)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                               & ! (in)
         & inquire_VBCSpecType(SaltBC_Surface), xya_Salt_VBCRHS(:,:,1),    & ! (in)
         & inquire_VBCSpecType(SaltBC_Bottom),  xya_Salt_VBCRHS(:,:,2)     & ! (in)
         & )


    !$omp parallel do
    do n = 1, TRC_TOT_NUM
       xyza_TRCA(:,:,:,n) = xyz_H0(:,:,:)*xyza_TRCA(:,:,:,n)/xyz_HA(:,:,:)
    end do

  end subroutine DOGCM_Phys_spm_VImplTRC

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_spm_VImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,          & ! (in)
       & xyz_H, xyz_VViscCoef, dt, alpha                & ! (in)
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
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xya_U_VBCRHS(0:iMax-1,jMax,2)
    real(DP) :: xya_V_VBCRHS(0:iMax-1,jMax,2)
    
    integer :: n

    ! 実行文; Executable statements
    !
    
    call DOGCM_Boundary_spm_InqVBCRHS_UV( &
         & xya_U_VBCRHS, xya_V_VBCRHS,                                     & ! (out)
         & xyz_U0, xyz_V0, xyz_H, xyz_VViscCoef                            & ! (in)
         & )
    

    call calc_VDiffEq( xyz_UA,                                             & ! (out)
         & xyz_U0, xyz_U_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
         & inquire_VBCSpecType(DynBC_Surface), xya_U_VBCRHS(:,:,1),        & ! (in)
         & inquire_VBCSpecType(DynBC_Bottom),  xya_U_VBCRHS(:,:,2)         & ! (in)
         & )

    call calc_VDiffEq( xyz_VA,                                             & ! (out)
         & xyz_V0, xyz_V_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
         & inquire_VBCSpecType(DynBC_Surface), xya_V_VBCRHS(:,:,1),        & ! (in)
         & inquire_VBCSpecType(DynBC_Bottom),  xya_V_VBCRHS(:,:,2)         & ! (in)
         & )
        
  end subroutine DOGCM_Phys_spm_VImplUV

  !-----------------------------------------------------------------------

  subroutine calc_VDiffEq( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_RHS,                     &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal           &  ! (in)
       & )
    
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
    integer :: N

    ! 実行文; Executable statements
    !

    N = size(IMat,1)
    IMat(:,:) = 0d0
    forAll(k=0:kMax) IMat(k,k) = 1d0
    do k = 0, kMax
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do

    !$omp parallel do collapse(2) private(z_Kv, z_e3, AMat, b, IPiv, info, i)
    do j = 1, jMax
       do i = 0, iMax-1

          ! Solve a linear equation system, 
          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 

          ! Set coefficient matrix A 
          !

          z_Kv(:) = xyz_Kv(i,j,:)
          z_e3(:) = xyz_H(i,j,:)

          do k=1, kMax-1
             AMat(k,:) =  z_e3(k)/dt * IMat(k,:) &
                  &      -  alpha * matmul( DSigMat(k,:) * z_Kv/z_e3, DSigMat(:,:) )

             b(k) = z_e3(k) * xyz_q_RHS(i,j,k)
          end do

          ! Set boundary condition
          !

          select case(SeaSfcBCType)
          case('D')
             AMat(0,:) = IMat(0,:)
          case('N')               
             AMat(0,:) = DSigMat(0,:)
          end select
          b(0) = 0d0!xy_SeaSfcBCVal(i,j)

          select case(SeaBtmBCType)
          case('D')
             AMat(kMax,:) = IMat(kMax,:)
          case('N')               
             AMat(kMax,:) = DSigMat(kMax,:)
          end select
          b(kMax) = 0d0!xy_SeaBtmBCVal(i,j)

          call DGESV(N, 1, AMat, N, IPIV, b, N, info)

          xyz_qA(i,j,:) = xyz_q0(i,j,:) + b(:)          
       end do
    end do

  end subroutine calc_VDiffEq
    
  !-----------------------------------------------------------------------
  
end module DOGCM_Phys_spm_mod
