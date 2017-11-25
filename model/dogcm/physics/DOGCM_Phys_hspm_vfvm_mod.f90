!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Phys_hspm_vfvm_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use ProfUtil_mod
  
  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax, nMax, &
       & z_CDK, z_FDK, z_RCDK, z_RFDK

  use DOGCM_Admin_Constants_mod, only: &
       & PI, RPlanet,               &
       & hViscCoef, hHyperViscCoef, &
       & hDiffCoef, hHyperDiffCoef
         
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM,        &
       & TRCID_PTEMP, TRCID_SALT

  use SpmlUtil_mod
  use VFvmUtil_mod
  
  use DOGCM_Admin_GovernEq_mod, only: &
       & isPhysicsCompActivated, &
       & OCNGOVERNEQ_LPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_LPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_LPHYS_REDIGM_NAME, &
       & OCNGOVERNEQ_VPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_VPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_VPHYS_CONVEC_NAME

  use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar, &
       & DOGCM_IO_History_HistPut
  
  !-- Module for the parametrizations of lateral oceanic physics

  use LPhys_DIFF_spm_mod, only: &
       & LPhys_DIFF_spm_Init, LPhys_DIFF_spm_Final, &
       & LPhys_DIFF_spm_LMixMOMRHS,                       &
       & LPhys_DIFF_spm_LMixMOMRHSImpl,                   &
       & LPhys_DIFF_spm_LMixTRCRHS,                       &
       & LPhys_DIFF_spm_LMixTRCRHSImpl,                   &
       & w_Filter, w_HDiffCoefH, w_HViscCoefH
       
  use LPhys_RediGM_hspm_vfvm_mod, only: &
       & LPhys_RediGM_hspm_vfvm_Init, LPhys_RediGM_hspm_vfvm_Final, &
       & LPhys_RediGM_hspm_vfvm_AddMixingTerm

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

  public :: DOGCM_Phys_hspm_vfvm_Init, DOGCM_Phys_hspm_vfvm_Final

  public :: DOGCM_Phys_hspm_vfvm_Do
  
  public :: DOGCM_Phys_hspm_vfvm_VMixMOMRHS
  public :: DOGCM_Phys_hspm_vfvm_VMixTRCRHS

  public :: DOGCM_Phys_hspm_vfvm_ImplUV
  public :: DOGCM_Phys_hspm_vfvm_ImplTRC

!!$  public :: reconstruct_sfctemp
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Phys_hspm_vfvm_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Phys_hspm_vfvm_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !
!    call read_nmlData(configNmlName)

    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then
       call LPhys_DIFF_spm_Init( configNmlName = configNmlName )
       call DOGCM_IO_History_RegistVar( 'PTemp_t_lphys_lmix', 'IJKT', 'PTemp tendency of lateral mixing', 'K/s' )
    end if
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_hspm_vfvm_Init( confignmlFileName = configNmlName )
       call DOGCM_IO_History_RegistVar( 'PTemp_t_lphys_RediGM', 'IJKT', 'PTemp tendency of RediGM', 'K/s' )
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Init()
       call DOGCM_IO_History_RegistVar( 'PTemp_t_vphys_CA', 'IJKT', 'PTemp tendency of convective adjustment', 'K/s' )
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
       call DOGCM_IO_History_RegistVar( 'PTemp_t_vphys_vmix', 'IJKT', 'PTemp tendency of vertical mixing', 'K/s' )
       call DOGCM_IO_History_RegistVar( 'PTemp_t_vphys_vmix_impl', 'IJKT', 'PTemp tendency of vertical mixing', 'K/s' )
    end if
    
  end subroutine DOGCM_Phys_hspm_vfvm_Init

  !>
  !!
  !!
  subroutine DOGCM_Phys_hspm_vfvm_Final()

    ! 実行文; Executable statements
    !
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_hspm_vfvm_Final()
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Final()
    end if
    
  end subroutine DOGCM_Phys_hspm_vfvm_Final

  !-----------------------------------------

  subroutine DOGCM_Phys_hspm_vfvm_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (inout)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,        & ! (inout)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt, lhst_tend                                         & ! (in)
       & )

    ! モジュール引用; Use statements
    !    
    use DOGCM_Admin_Grid_mod, only: &
         & z_KAXIS_Weight
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyz_ConvIndex

    use SpmlUtil_mod, only: &
         & AvrLonLat_xy, xy_IntSig_BtmToTop_xyz

    use DOGCM_Boundary_vars_mod, only: &
         & xy_SfcHFlx_sr, xy_SfcHFlx_ns
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(inout) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(inout) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(inout) :: xy_BtmFrictCoef(IA,JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_TOPO(IA,JA)
    real(DP), intent(in) :: dt
    logical,  intent(in) :: lhst_tend
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: avr_ptemp_RHS_phys
    integer :: k
    
!!$    real(DP) :: xyz_UA(0:iMax-1,jMax,KA)
!!$    real(DP) :: xyz_VA(0:iMax-1,jMax,KA)
!!$    real(DP) :: xyza_TRCA(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    
    ! 実行文; Executable statements
    !

    !-- Horizontal momentum -----------------------------------------------------

    call ProfUtil_RapStart('OcnPhys_Main', 3)
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXMOM_NAME )  ) then
       call ProfUtil_RapStart('OcnPhys_MOMLMix', 3)
!!$       call LPhys_DIFF_spm_LMixMOMRHS(            &
!!$            & xyz_U_RHS_phy, xyz_V_RHS_phy,               & ! (inout)
!!$            & xyz_U, xyz_V, xyz_H,                        & ! (in)
!!$            & hViscCoef, hHyperViscCoef                   & ! (in)
!!$            & )
!!$       call LPhys_DIFF_spm_LMixMOMRHSImpl( &
!!$            & xyz_U_RHS_phy, xyz_V_RHS_phy,               & ! (inout)
!!$            & xyz_U, xyz_V, xyz_H,                        & ! (in)
!!$            & hViscCoef, hHyperViscCoef, dt               & ! (in)
!!$            & )
      call ProfUtil_RapEnd('OcnPhys_MOMLMix', 3)       
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXMOM_NAME )  ) then
       call ProfUtil_RapStart('OcnPhys_MOMVMix', 3)
       call DOGCM_Phys_hspm_vfvm_VMixMOMRHS( xyz_U_RHS_phy, xyz_V_RHS_phy,   & ! (inout)
            & xyz_U, xyz_V, xyz_H, xyz_VViscCoef, xy_BtmFrictCoef            & ! (in)
            & )
       call ProfUtil_RapEnd('OcnPhys_MOMVMix', 3)

    end if

    !-- Tracer ------------------------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then

       call ProfUtil_RapStart('OcnPhys_TRCLMix', 3)
       xyza_TRC_RHS = 0d0
!!$       call LPhys_DIFF_spm_LMixTRCRHSImpl( xyza_TRC_RHS_phy,       & ! (inout)
!!$       call LPhys_DIFF_spm_LMixTRCRHSImpl( xyza_TRC_RHS,       & ! (inout)
!!$            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef, dt       &  ! (in)
!!$            & )
       xyza_TRC_RHS_phy = xyza_TRC_RHS_phy + xyza_TRC_RHS
!!$       if (lhst_tend) then
!!$          call DOGCM_IO_History_HistPut( 'PTemp_t_lphys_lmix', xyza_TRC_RHS(IS:IE,JS:JE,KS:KE,TRCID_PTEMP))
!!$          write(*,*) "MIXTRC:", xyza_TRC(IS,JS+19,KS,TRCID_PTEMP)*0d0 + xyza_TRC_RHS(IS,JS+43,KS,TRCID_PTEMP)*dt
!!$       end if
       call ProfUtil_RapEnd('OcnPhys_TRCLMix', 3)       
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(:,:,:, TRCID_SALT), xyz_H))
!!$       write(*,*) "->avr_ptemp_phys (+ LMixRHS): ", avr_ptemp_RHS_phys/35d0*dt*1d3
    end if

    !-- Lateral mixing of tracers by eddy induced velocity -------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then

       call ProfUtil_RapStart('OcnPhys_TRCRediGM', 3)              
       xyza_TRC_RHS = 0d0
       call LPhys_RediGM_hspm_vfvm_AddMixingTerm( &
            & xyza_TRC_RHS(:,:,:,TRCID_PTEMP), xyza_TRC_RHS(:,:,:,TRCID_SALT),  & ! (inout)
!!$            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP), xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),  & ! (inout)
            & xyz_VDiffCoef,                                                            & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),                  & ! (in)
            & xyz_H, xyz_Z, xy_Topo                                                     & ! (in)
            & )
       xyza_TRC_RHS_phy = xyza_TRC_RHS_phy + xyza_TRC_RHS
!!$       if( lhst_tend ) then
!!$          call DOGCM_IO_History_HistPut( 'PTemp_t_lphys_RediGM', xyza_TRC_RHS(IS:IE,JS:JE,KS:KE,TRCID_PTEMP))
!!$          write(*,*) "RediGM:", xyza_TRC(IS,JS+19,KS,TRCID_PTEMP)*0d0 + xyza_TRC_RHS(IS,JS+10:28,KS,TRCID_PTEMP)*43200d0
!!$       end if
       call ProfUtil_RapEnd('OcnPhys_TRCRediGM', 3)                     
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(:,:,:, TRCID_SALT), xyz_H))
!!$       write(*,*) "avr_ptemp_phys (+ RediGMRHS): ",  avr_ptemp_RHS_phys/35d0*dt*1d3
    end if

    !-- Vertical mixing of tracers by non-penetrative convection -------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then

       call ProfUtil_RapStart('OcnPhys_ConvAdjust', 3)                     
       xyza_TRC_RHS = 0d0
       call DOGCM_VPhys_ConvAdjust_AddMixingTerm( &
            & xyza_TRC_RHS(:,:,:,TRCID_PTEMP), xyza_TRC_RHS(:,:,:,TRCID_SALT),  & ! (inout)
!!$            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP), xyza_TRC_RHS_phy(:,:,:,TRCID_SALT), & ! (inout)
            & xyz_ConvIndex,                                                           & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),                 & ! (in)
            & xyz_H, xyz_Z, z_KAXIS_Weight, dt                                         & ! (in)
            & )
       xyza_TRC_RHS_phy = xyza_TRC_RHS_phy + xyza_TRC_RHS
!!$       if( lhst_tend ) then
!!$          call DOGCM_IO_History_HistPut( 'PTemp_t_vphys_CA', xyza_TRC_RHS(IS:IE,JS:JE,KS:KE,TRCID_PTEMP))
!!$          write(*,*) "Conv:", 0d0*xyza_TRC(IS,JS+19,KS,TRCID_PTEMP) + xyza_TRC_RHS(IS,JS+10:28,KS,TRCID_PTEMP)*43200d0
!!$       end if
       call ProfUtil_RapEnd('OcnPhys_ConvAdjust', 3)
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(:,:,:, TRCID_SALT), xyz_H))
!!$       write(*,*) "avr_ptemp_phys (+ ConvAdjustRHS): ",  avr_ptemp_RHS_phys/35d0*dt*1d3
    end if


    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
       call ProfUtil_RapStart('OcnPhys_VMixTRC', 3)
       xyza_TRC_RHS = 0d0
       call DOGCM_Phys_hspm_vfvm_VMixTRCRHS( xyza_TRC_RHS,   & ! (inout)
!!$       call DOGCM_Phys_hspm_vfvm_VMixTRCRHS( xyza_TRC_RHS_phy,   & ! (inout)
            & xyza_TRC, xyz_H, xyz_VDiffCoef                     & ! (in)
            & )
       xyza_TRC_RHS_phy = xyza_TRC_RHS_phy + xyza_TRC_RHS
       if( lhst_tend ) then
!          call DOGCM_IO_History_HistPut( 'PTemp_t_vphys_vmix', xyza_TRC_RHS(IS:IE,JS:JE,KS:KE,TRCID_PTEMP))
!          write(*,*) "VMIXTRC:", xyza_TRC_RHS(IS,JS+10:28,KS,TRCID_PTEMP)*43200d0
!!$          write(*,*) "PTempKS:", xyza_TRC(IS,JS+10:28,KS,TRCID_PTEMP)
!!$          write(*,*) "PTempKS+1:", xyza_TRC(IS,JS+10:28,KS+1,TRCID_PTEMP)
!!$          write(*,*) "SfcHFlx:", xy_SfcHFlx_ns(IS,JS+10:28) + xy_SfcHFlx_sr(IS,JS+10:28)

          !          write(*,*) "PTempP:", xyza_TRC(IS,JS+10:23,KS,TRCID_PTEMP) + xyza_TRC_RHS_phy(IS,JS+10:23,KS,TRCID_PTEMP)*dt
!          write(*,*) "Phys*:", xyza_TRC_RHS_phy(IS,JS+10:23,KS,TRCID_PTEMP)*dt
       end if
       call ProfUtil_RapEnd('OcnPhys_VMixTRC', 3)
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(:,:,:, TRCID_SALT), xyz_H))
!!$       write(*,*) "avr_ptemp_phys (+ VMixRHS): ",  avr_ptemp_RHS_phys/35d0*dt*1d3
!!$       write(*,*) "avr_ptemp_phys_local:", IntSig_BtmToTop(xyza_TRC_RHS_phy(0,jMax/2,:,TRCID_SALT))
    end if

    call ProfUtil_RapEnd('OcnPhys_Main', 3)
    
  end subroutine DOGCM_Phys_hspm_vfvm_Do

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_VMixMOMRHS(     &
       & xyz_U_RHS, xyz_V_RHS,                                  & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_VViscCoef, xy_BtmFrictCoef    & ! (in)
       & )

    ! モジュール引用; Use statements
    !        
    use DOGCM_Admin_Constants_mod, only: &
         & RefDens    

    use DOGCM_Admin_BC_mod, only: &
         & DynBC_Surface, DynBC_Bottom, &
         & DynBCTYPE_NoSlip
    
    use DOGCM_Boundary_vars_mod, only: &
         & xy_WindStressU, xy_WindStressV
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xy_BtmFrictCoef(IA,JA)
    
    ! 局所変数
    ! Local variables
    !

    real(DP) :: xy_Tmp(IA,JA)
    real(DP) :: xyr_DiffFlxU(IA,JA,KA)
    real(DP) :: xyr_DiffFlxV(IA,JA,KA)
    real(DP) :: m1
    real(DP) :: m2
    
    integer :: k
    
    ! 実行文; Executable statements
    !

    !$omp parallel do private(xy_Tmp, m1, m2)
    do k=KS, KE-1
       m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))       

       xy_Tmp(:,:) =   z_RFDK(k)*(m1*xyz_VViscCoef(:,:,k) + m2*xyz_VViscCoef(:,:,k+1)) &
            &                   /(m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1))
       
       xyr_DiffFlxU(:,:,k) = xy_Tmp(:,:)*(xyz_U(:,:,k) - xyz_U(:,:,k+1))
       xyr_DiffFlxV(:,:,k) = xy_Tmp(:,:)*(xyz_V(:,:,k) - xyz_V(:,:,k+1))
    end do
    
    select case(DynBC_Surface)
    case(DynBCTYPE_NoSlip)
       !$omp parallel
       !$omp workshare
       xyr_DiffFlxU(:,:,KS-1) = - z_RFDK(KS)*xyz_VViscCoef(:,:,KS)/xyz_H(:,:,KS)* 2d0*xyz_U(:,:,KS)
       xyr_DiffFlxV(:,:,KS-1) = - z_RFDK(KS)*xyz_VViscCoef(:,:,KS)/xyz_H(:,:,KS)* 2d0*xyz_V(:,:,KS)
       !$omp end workshare
       !$omp end parallel
    case default
       xyr_DiffFlxU(:,:,KS-1) = xy_WindStressU(:,:)/RefDens
       xyr_DiffFlxV(:,:,KS-1) = xy_WindStressV(:,:)/RefDens
    end select

    select case(DynBC_Bottom)
    case(DynBCTYPE_NoSlip)
       !$omp parallel
       !$omp workshare
       xyr_DiffFlxU(:,:,KE  ) = xy_BtmFrictCoef * xyz_U(:,:,KE)
       xyr_DiffFlxV(:,:,KE  ) = xy_BtmFrictCoef * xyz_V(:,:,KE)
       !$omp end workshare
       !$omp end parallel
    case default
       xyr_DiffFlxU(:,:,KE) = 0d0
       xyr_DiffFlxV(:,:,KE) = 0d0
    end select
    
    
    !$omp parallel do
    do k=KS, KE
       xyz_U_RHS(:,:,k) = xyz_U_RHS(:,:,k) + &
            & (xyr_DiffFlxU(:,:,k-1) - xyr_DiffFlxU(:,:,k))*z_RCDK(k)/xyz_H(:,:,k)

       xyz_V_RHS(:,:,k) = xyz_V_RHS(:,:,k) + &
            & (xyr_DiffFlxV(:,:,k-1) - xyr_DiffFlxV(:,:,k))*z_RCDK(k)/xyz_H(:,:,k)
    end do
    
  end subroutine DOGCM_Phys_hspm_vfvm_VMixMOMRHS

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_VMixTRCRHS(     &
       & xyza_TRC_RHS,                       & ! (out)
       & xyza_TRC, xyz_H, xyz_VDiffCoef      & ! (in)
       & )

    ! モジュール引用; Use statements
    !        
    use DOGCM_Admin_Constants_mod, only: &
         & RefDens, RefSalt, Cp0
    
    use DOGCM_Admin_BC_mod, only: &
         & ThermBC_Surface, ThermBC_Bottom, &
         & SaltBC_Surface, SaltBC_Bottom,   &
         & ThermBCTYPE_PrescTemp,           &
         & SaltBCTYPE_PrescSalt

    use DOGCM_Boundary_vars_mod, only: &
         & xy_SfcHFlx_sr, xy_SfcHFlx_ns, xy_SeaSfcTemp, &
         & xy_FreshWtFlxS, xy_SeaSfcSalt

    use SpmlUtil_mod, only: &
         & AvrLonLat_xy, xy_IntSig_BtmToTop_xyz
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)
    
    ! 局所変数
    ! Local variables
    !    

    real(DP) :: xy_Tmp(IA,JA)
    real(DP) :: xyra_DiffFlxTRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP) :: m1
    real(DP) :: m2
    
    integer :: k    
    integer :: n

    ! 実行文; Executable statements
    !

    !* Calculate the diffusive flux 
    !
    
    !$omp parallel do private(xy_Tmp, m1, m2, n)
    do k=KS, KE-1
       m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))       
       
       xy_Tmp(:,:) =   z_RFDK(k)*(m1*xyz_VDiffCoef(:,:,k) + m2*xyz_VDiffCoef(:,:,k+1)) &
            &                   /(m1*xyz_H(:,:,k) + m2*xyz_H(:,:,k+1))
       do n=1, TRC_TOT_NUM
          xyra_DiffFlxTRC(:,:,k,n) = xy_Tmp(:,:)*(xyza_TRC(:,:,k,n) - xyza_TRC(:,:,k+1,n))
       end do
    end do
    
    !* Set the diffusive fluxe at sea surface and bottom
    !

    ! Set zero as the default values for all tracers 
    xyra_DiffFlxTRC(:,:,KS-1,:) = 0d0
    xyra_DiffFlxTRC(:,:,KE  ,:) = 0d0

    ! for potential temperature
    k = KS; n = TRCID_PTEMP
    select case(ThermBC_Surface)
    case(ThermBCTYPE_PrescTemp)
       xyra_DiffFlxTRC(:,:,KS-1,n) = - z_RFDK(k)*xyz_VDiffCoef(:,:,KS)/xyz_H(:,:,KS)             &
            &                        * 2d0*(xyza_TRC(:,:,KS,n) - xy_SeaSfcTemp(:,:))
    case default
       xyra_DiffFlxTRC(:,:,KS-1,n) = - (xy_SfcHFlx_ns(:,:) + xy_SfcHFlx_sr(:,:)) &
            &                        / (RefDens * Cp0)
    end select

    ! for salinity    
    k = KS; n = TRCID_SALT
    select case(SaltBC_Surface)
    case(SaltBCTYPE_PrescSalt)
       xyra_DiffFlxTRC(:,:,KS-1,n) = - z_RFDK(k)*xyz_VDiffCoef(:,:,KS)/xyz_H(:,:,KS)             &
            &                        * 2d0*(xyza_TRC(:,:,KS,n) - xy_SeaSfcSalt(:,:))
    case default
       xyra_DiffFlxTRC(:,:,KS-1,n) = - xy_FreshWtFlxS(:,:) * RefSalt ! Virtual salinity flux (Set 35 psu as reference salinity)
    end select

    ! Calculate the divergence of diffusive flux
    !
    
    !$omp parallel do collapse(2)
    do n=1, TRC_TOT_NUM
       do k=KS, KE
          xyza_TRC_RHS(:,:,k,n) = xyza_TRC_RHS(:,:,k,n) + &
               & (xyra_DiffFlxTRC(:,:,k-1,n) - xyra_DiffFlxTRC(:,:,k,n))*z_RCDK(k)/xyz_H(:,:,k)       
       end do
    end do
    
  end subroutine DOGCM_Phys_hspm_vfvm_VMixTRCRHS

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_ImplTRC( xyza_TRCA,         & ! (out)
       & xyza_TRC0, xyza_HTRC_RHS,                       & ! (in)
       & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, alpha,       & ! (in)
       & lhst_tend                                       & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & ThermBC_Surface, ThermBC_Bottom, &
         & SaltBC_Surface, SaltBC_Bottom

    use DOGCM_Admin_Grid_mod, only: &
         & xyz_Z, z_KAXIS_Weight, xy_Lat
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyz_ConvIndex
    
    use DOGCM_Boundary_hspm_vfvm_mod, only: &
         & DOGCM_Boundary_hspm_vfvm_InqVBCRHS_TRC

    
    use LPhys_DIFF_spm_mod, only: &
         & LPhys_DIFF_AdaptiveFilter4SIce
    
    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyza_TRCA(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC0(IA,JA,KA,TRC_TOT_NUM)    
    real(DP), intent(in) :: xyza_HTRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_HA(IA,JA,KA)
    real(DP), intent(in) :: xyz_H0(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    logical, intent(in) :: lhst_tend

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xya_PTemp_VBCRHS(IA,JA,2)
    real(DP) :: xya_Salt_VBCRHS(IA,JA,2)
    integer :: n
    integer :: k
    integer :: j
    
    real(DP) :: avr_RHS
    real(DP) :: avr_TRC0_phys
    real(DP) :: avr_TRCA_phys

    real(DP) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    
    ! 実行文; Executable statements
    !

    
    call DOGCM_Boundary_hspm_vfvm_InqVBCRHS_TRC( &
         & xya_PTemp_VBCRHS, xya_Salt_VBCRHS,                                & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP), xyza_TRC0(:,:,:,TRCID_SALT),        & ! (in)
         & xyz_H0, xyz_VDiffCoef                                             & ! (in)
         & )

    call calc_VDiffTRCEq( xyza_TRCA(:,:,:,TRCID_PTEMP),                     & ! (out)
         & xyza_TRC0(:,:,:,TRCID_PTEMP),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_PTEMP),                                & ! (in)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                                & ! (in)
         & inquire_VBCSpecType(ThermBC_Surface), xya_PTemp_VBCRHS(:,:,1),   & ! (in)
         & inquire_VBCSpecType(ThermBC_Bottom),  xya_PTemp_VBCRHS(:,:,2),   & ! (in)
         & "PTemp"                                                          & ! (in)
         & )

    
    call calc_VDiffTRCEq( xyza_TRCA(:,:,:,TRCID_SALT),                     & ! (out)
         & xyza_TRC0(:,:,:,TRCID_SALT),                                    & ! (in)
         & xyza_HTRC_RHS(:,:,:,TRCID_SALT),                                & ! (in)
         & xyz_VDiffCoef, xyz_H0, dt, alpha,                               & ! (in)
         & inquire_VBCSpecType(SaltBC_Surface), xya_Salt_VBCRHS(:,:,1),    & ! (in)
         & inquire_VBCSpecType(SaltBC_Bottom),  xya_Salt_VBCRHS(:,:,2),    & ! (in)
         & "Salt"                                                          & ! (in)
         & )


!!$    if (lhst_tend ) then
!!$       write(*,*) "PTempA*:",  xyza_TRC0(IS,JS+10:28,KS,TRCID_PTEMP) &
!!$            &                + xyza_HTRC_RHS(IS,JS+10:28,KS,TRCID_PTEMP)/xyz_H0(IS,JS+10:28,KS)*dt
!!$       write(*,*) "PTempA :", xyza_TRCA(IS,JS+10:28,KS,TRCID_PTEMP)
!!$    end if

    !$omp parallel do private(n,k) collapse(2)
    do n = 1, TRC_TOT_NUM
    do k = KS, KE
       xyza_TRCA(IS:IE,JS:JE,k,n) = xy_w( &
            & w_Filter * w_xy(xyza_TRCA(IS:IE,JS:JE,k,n)) / (1d0 - dt*w_HDiffCoefH) &
            & )
       xyza_TRCA(:,:,k,n) = xyz_H0(:,:,k)*xyza_TRCA(:,:,k,n)/xyz_HA(:,:,k)
    end do
    end do    
    
!!$    if (lhst_tend ) then
!!$       write(*,*) "PTempAA :", xyza_TRCA(IS,JS+10:28,KS,TRCID_PTEMP)
!!$       
!!$!       call DOGCM_IO_History_HistPut( 'PTemp_t_vphys_vmix_impl', &
!!$!            & xyza_TRCA(IS:IE,JS:JE,KS:KE,TRCID_PTEMP)           &
!!$!            & - (xyza_TRC0(IS:IE,JS:JE,KS:KE,TRCID_PTEMP) + dt*xyza_HTRC_RHS(IS:IE,JS:JE,KS:KE,TRCID_PTEMP)/xyz_H0(IS:IE,JS:JE,KS:KE)) &
!!$!            & )
!!$    end if
!!$    avr_TRC0_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC0(IS:IE,JS:JE,:,TRCID_SALT), xyz_H0) )
!!$    avr_TRCA_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRCA(IS:IE,JS:JE,:,TRCID_SALT), xyz_H0) )
!!$    avr_RHS = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_HTRC_RHS(IS:IE,JS:JE,:,TRCID_SALT)/xyz_H0, xyz_H0) )
!!$    write(*,*) "Phys+Dyn, Phys+Dyn+Impl Salt:", &
!!$         & avr_RHS,                          &
!!$         & (avr_TRCA_phys - avr_TRC0_phys)/dt, & !/dt*dt
!!$         & AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRCA(IS:IE,JS:JE,:,TRCID_SALT) - xyza_TRC0(IS:IE,JS:JE,:,TRCID_SALT), xyz_H0) )/dt
    
  end subroutine DOGCM_Phys_hspm_vfvm_ImplTRC
  
  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_ImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,                & ! (in)
       & xyz_H, xyz_VViscCoef, xy_BtmFrictCoef, dt, alpha     & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & DynBC_Surface, DynBC_Bottom
    
    use DOGCM_Boundary_hspm_vfvm_mod, only: &
         & DOGCM_Boundary_hspm_vfvm_InqVBCRHS_UV

    use LPhys_DIFF_spm_mod, only: w_Filter
    
    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyz_UA(IA,JA,KA)
    real(DP), intent(out) :: xyz_VA(IA,JA,KA)        
    real(DP), intent(in) :: xyz_U0(IA,JA,KA)
    real(DP), intent(in) :: xyz_V0(IA,JA,KA)    
    real(DP), intent(in) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xy_BtmFrictCoef(IA,JA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xya_U_VBCRHS(IA,JA,2)
    real(DP) :: xya_V_VBCRHS(IA,JA,2)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    integer :: k

    ! 実行文; Executable statements
    !

    if (alpha > 0d0) then
       call DOGCM_Boundary_hspm_vfvm_InqVBCRHS_UV( &
            & xya_U_VBCRHS, xya_V_VBCRHS,                                     & ! (out)
            & xyz_U0, xyz_V0, xyz_H, xyz_VViscCoef                            & ! (in)
            & )
    
       call ProfUtil_RapStart('OcnPhys_ImplMom2', 3)
       call calc_VDiffMomEq( xyz_UA,                                          & ! (out)
            & xyz_U0, xyz_U_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
            & inquire_VBCSpecType(DynBC_Surface), xya_U_VBCRHS(:,:,1),        & ! (in)
            & inquire_VBCSpecType(DynBC_Bottom),  xya_U_VBCRHS(:,:,2),        & ! (in)
            & "U", xy_BtmFrictCoef                                            & ! (in)
            & )

       call calc_VDiffMomEq( xyz_VA,                                          & ! (out)
            & xyz_V0, xyz_V_RHS, xyz_VViscCoef, xyz_H, dt, alpha,             & ! (in)
            & inquire_VBCSpecType(DynBC_Surface), xya_V_VBCRHS(:,:,1),        & ! (in)
            & inquire_VBCSpecType(DynBC_Bottom),  xya_V_VBCRHS(:,:,2),        & ! (in)
            & "V",  xy_BtmFrictCoef                                           & ! (in)
            & )
       call ProfUtil_RapEnd('OcnPhys_ImplMom2', 3)
    end if

    !$omp parallel do private(w_Vor, w_Div)
    do k=KS, KE
       call calc_UVCosLat2VorDiv( &
            & xyz_UA(IS:IE,JS:JE,k)*xy_CosLat, xyz_VA(IS:IE,JS:JE,k)*xy_CosLat, &
            & w_Vor, w_Div                                                      &
            & )
       
       w_Vor(:) = w_Filter(:)*w_Vor/(1d0 - dt*w_HViscCoefH)
       w_Div(:) = w_Filter(:)*w_Div/(1d0 - dt*w_HViscCoefH)
       call calc_VorDiv2UV( w_Vor, w_Div,                     &
            & xyz_UA(IS:IE,JS:JE,k), xyz_VA(IS:IE,JS:JE,k) )
    end do

    
  end subroutine DOGCM_Phys_hspm_vfvm_ImplUV

  !-----------------------------------------------------------------------

  subroutine calc_VDiffTRCEq( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_HRHS,                       &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,                 &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,             &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal,             &  ! (in)
       & qname                                     &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(out) :: xyz_qA(IA,JA,KA)
    real(DP), intent(in) :: xyz_q0(IA,JA,KA)
    real(DP), intent(in) :: xyz_q_HRHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_Kv(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    character, intent(in) :: SeaSfcBCType
    real(DP), intent(in) :: xy_SeaSfcBCVal(IA,JA)
    character, intent(in) :: SeaBtmBCType
    real(DP), intent(in) :: xy_SeaBtmBCVal(IA,JA)
    character(*), intent(in) :: qname

    ! 局所変数
    ! Local variables
    !
    real(DP) :: z_Kv(KA)
    real(DP) :: z_e3(KA)
    real(DP) :: r_DFlxCoef(KA)
    real(DP) :: coef(KA)
    real(DP) :: m1
    real(DP) :: m2    

    real(DP) :: Al(KS:KE)
    real(DP) :: Ad(KS:KE)
    real(DP) :: Au(KS:KE)
    real(DP) :: b(KS:KE)
    integer :: info
    
    integer :: i
    integer :: j
    integer :: k
    integer :: N
    
   
    ! 実行文; Executable statements
    !

    N = KE - KS + 1

    !$omp parallel do collapse(2) private(z_Kv, z_e3, Al, Ad, Au, b, info, i, j, k, r_DFlxcoef, coef, m1, m2)
    do j = JS, JE
       do i = IS, IE

          ! Solve a linear equation system, 
          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 

          ! Set coefficient matrix A 
          !

          z_Kv(:) = xyz_Kv(i,j,:)
          z_e3(:) = xyz_H(i,j,:)
          coef(:) = dt/(z_e3(:)*z_CDK(:))

          do k=KS, KE-1
             m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
             m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))       
             r_DFlxCoef(k) =   (m1*z_Kv(k) + m2*z_Kv(k+1))            &
                  &          / (z_FDK(k)*(m1*z_e3(k) + m2*z_e3(k+1)))
          end do
          r_DFlxCoef(KS-1) = z_Kv(KS)/(z_FDK(KS-1)*z_e3(KS))
          r_DFlxCoef(KE  ) = z_Kv(KE)/(z_FDK(KE  )*z_e3(KE))
          
          do k=KS+1, KE-1
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k-1) + r_DFlxCoef(k))
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k)
          end do
          b(KS:KE) = dt * xyz_q_HRHS(i,j,KS:KE)/z_e3(KS:KE)

          ! Set boundary condition
          !
          
          k = KS
          select case(SeaSfcBCType)
          case('D')
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k) + 2d0*r_DFlxCoef(k-1))
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k) 
          case('N')
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k)                      )
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k) 
          end select

          k = KE
          select case(SeaBtmBCType)
          case('D')
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(2d0*r_DFlxCoef(k)    + r_DFlxCoef(k-1))
          case('N')
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(                  + r_DFlxCoef(k-1))
          end select
          
          call DGTSV(N, 1, Al(KS+1:KE), Ad(KS:KE), Au(KS:KE-1), b, N, info)

          xyz_qA(i,j,KS:KE) = xyz_q0(i,j,KS:KE) + b(KS:KE)          
       end do
    end do
    
  end subroutine calc_VDiffTRCEq

  !---------------------------------------------------------------------------
  
  subroutine calc_VDiffMOMEq( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_RHS,                     &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal,          &  ! (in)
       & qname, xy_BtmFrictCoef                 &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(out) :: xyz_qA(IA,JA,KA)
    real(DP), intent(in) :: xyz_q0(IA,JA,KA)
    real(DP), intent(in) :: xyz_q_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_Kv(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    character, intent(in) :: SeaSfcBCType
    real(DP), intent(in) :: xy_SeaSfcBCVal(IA,JA)
    character, intent(in) :: SeaBtmBCType
    real(DP), intent(in) :: xy_SeaBtmBCVal(IA,JA)
    character(*), intent(in) :: qname
    real(DP), intent(in) :: xy_BtmFrictCoef(IA,JA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: z_Kv(KA)
    real(DP) :: z_e3(KA)
    real(DP) :: r_DFlxCoef(KA)
    real(DP) :: coef(KA)
    real(DP) :: m1
    real(DP) :: m2    

    real(DP) :: Al(KS:KE)
    real(DP) :: Ad(KS:KE)
    real(DP) :: Au(KS:KE)
    real(DP) :: b(KS:KE)
    integer :: info
    
    integer :: i
    integer :: j
    integer :: k
    integer :: N
    
   
    ! 実行文; Executable statements
    !

    N = KE - KS + 1
    
    !$omp parallel do collapse(2) private(z_Kv, z_e3, Al, Ad, Au, b, info, i, j, k, r_DFlxcoef, coef, m1, m2)
    do j = JS, JE
       do i = IS, IE

          ! Solve a linear equation system, 
          ! [ e3*I/dt - DSig [ Kv/e3 DSig ] ] ( q^n+1 - q^0 ) =  e3 * RHS. 

          ! Set coefficient matrix A 
          !

          z_Kv(:) = xyz_Kv(i,j,:)
          z_e3(:) = xyz_H(i,j,:)
          coef(:) = dt/(z_e3(:)*z_CDK(:))

          do k=KS, KE-1
             m1 = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
             m2 = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))       
             r_DFlxCoef(k) =   (m1*z_Kv(k) + m2*z_Kv(k+1))            &
                  &          / (z_FDK(k)*(m1*z_e3(k) + m2*z_e3(k+1)))
          end do
          r_DFlxCoef(KS-1) = z_Kv(KS)/(z_FDK(KS-1)*z_e3(KS))
          r_DFlxCoef(KE  ) = z_Kv(KE)/(z_FDK(KE  )*z_e3(KE))
          
          do k=KS+1, KE-1
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k-1) + r_DFlxCoef(k))
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k)
          end do
          b(KS:KE) = dt * xyz_q_RHS(i,j,KS:KE)

          ! Set boundary condition
          !
          
          k = KS
          select case(SeaSfcBCType)
          case('D')
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k) + 2d0*r_DFlxCoef(k-1))
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k) 
          case('N')
             Ad(k) = 1d0 + alpha*coef(k)*(r_DFlxCoef(k)                      )
             Au(k) = - alpha*coef(k)*r_DFlxCoef(k) 
          end select

          k = KE
          select case(SeaBtmBCType)
          case('D')
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(xy_BtmFrictCoef(i,j) + r_DFlxCoef(k-1))
          case('N')
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(                  + r_DFlxCoef(k-1))
          end select
          
          call DGTSV(N, 1, Al(KS+1:KE), Ad(KS:KE), Au(KS:KE-1), b, N, info)

          xyz_qA(i,j,KS:KE) = xyz_q0(i,j,KS:KE) + b(KS:KE)          
       end do
    end do
    
  end subroutine calc_VDiffMOMEq

  !----

!!$  subroutine reconstruct_sfctemp( xy_SfcTemp )
!!$
!!$    use DOGCM_Admin_Grid_mod, only: &
!!$         & xy_Lat
!!$
!!$    use DOGCM_Boundary_Vars_mod, only: &
!!$         & xy_OcnSfcCellMask, OCNCELLMASK_SICE, OCNCELLMASK_OCEAN
!!$    
!!$    real(DP), intent(inout) :: xy_SfcTemp(IA,JA)
!!$    
!!$    real(DP) :: xy_SfcTempMod(IA,JA)
!!$    integer :: j
!!$    logical :: FlagUnderSIce
!!$
!!$    integer :: Mint
!!$    integer :: M
!!$    integer :: lambda
!!$    real(DP) :: SubDomLatStart
!!$    real(DP) :: SubDomLatEnd
!!$    logical :: do_reconstruct
!!$
!!$    FlagUnderSIce = .false.
!!$    if (xy_OcnSfcCellMask(IS,JS) == OCNCELLMASK_SICE) FlagUnderSIce = .true.
!!$    SubDomLatStart = - PI/2d0
!!$
!!$    do j=JS+1, JE
!!$       do_reconstruct = .false.
!!$       if (FlagUnderSIce) then
!!$          if (xy_OcnSfcCellMask(IS,j) == OCNCELLMASK_OCEAN) then
!!$             SubDomLatEnd = 0.5d0*(xy_Lat(IS,j-1) + xy_Lat(IS,j))
!!$             M = 10; lambda = 14; Mint = 60
!!$             do_reconstruct = .true.
!!$          else if(j==JE) then
!!$             SubDomLatEnd = PI/2d0
!!$             M = 10; lambda = 14; Mint = 60
!!$             do_reconstruct = .true.
!!$          end if
!!$       else
!!$          if (xy_OcnSfcCellMask(IS,j) == OCNCELLMASK_SICE) then
!!$             SubDomLatEnd = 0.5d0*(xy_Lat(IS,j-1) + xy_Lat(IS,j))
!!$             M = 18; lambda = 22; Mint = 60             
!!$             do_reconstruct = .true.
!!$          else if(j==JE) then
!!$             SubDomLatEnd = PI/2d0
!!$             M = 18; lambda = 22; Mint = 60             
!!$             do_reconstruct = .true.
!!$          end if
!!$       end if
!!$
!!$       if (do_reconstruct) then
!!$          write(*,*) "call freud_reconstruct:", &
!!$               & SubDomLatStart*180d0/PI, "-", SubDomLatEnd*180d0/PI, "M,lambda,Mint=", M, lambda, Mint
!!$          call freud_reconstruct( xy_SfcTempMod, &
!!$               & xy_SfcTemp(:,:), (/ SubDomLatStart, SubDomLatEnd /), M, lambda, Mint )
!!$          SubDomLatStart = SubDomLatEnd
!!$          FlagUnderSIce = .not.FlagUnderSIce
!!$       end if
!!$    end do
!!$
!    write(*,*) "en ori:", AvrLonLat_xy(xyz_PTemp(IS:IE,JS:JE,KS))
!!$    xy_SfcTemp(:,:) = xy_SfcTempMod(:,:) * &
!!$         & AvrLonLat_xy(xy_SfcTemp(IS:IE,JS:JE))/AvrLonLat_xy(xy_SfcTempMod(IS:IE,JS:JE))
!    do j=JS, JE
!       write(*,*) xy_Lat(IS,j), xyz_PTemp(IS,j,KS), xy_SfcTempMod(IS,j)
!    end do
!    write(*,*) "en check:", AvrLonLat_xy(xyz_PTemp(IS:IE,JS:JE,KS))
!!$    
!!$  end subroutine reconstruct_sfctemp
  
!!$  subroutine gegenbauer_reconstruct( &
!!$       xy_q, alat )
!!$
!!$    use w_zonal_module_sjpack, only: Interpolate_w, y_Lat
!!$    use DOGCM_Admin_Grid_mod, only: xy_Lat
!!$    real(DP), intent(inout) :: xy_q(IA,JA)
!!$    real(DP), intent(in) :: alat(:)
!!$
!!$    real(DP) :: a
!!$    real(DP) :: b
!!$    real(DP) :: eps
!!$    real(DP) :: del
!!$    real(DP) :: alpha
!!$    real(DP) :: beta
!!$
!!$    integer :: NLat
!!$    integer :: Mint
!!$    integer :: M
!!$    integer :: i
!!$    integer :: l
!!$    real :: lamb
!!$    integer :: j
!!$    
!!$    real(DP), allocatable :: xi(:)
!!$    real(DP), allocatable :: g(:)
!!$    real(DP), allocatable :: int_weight(:)
!!$    real(DP), allocatable :: lon_intp(:)
!!$    real(DP), allocatable :: lat_intp(:)
!!$    real(DP), allocatable :: q_intp(:)
!!$    real(DP), parameter :: PI = acos(-1d0)
!!$    real(DP), allocatable :: q_rc(:)
!!$    real(DP), allocatable :: xi_rc(:)
!!$    real(DP), allocatable :: Clamb_rc(:,:)
!!$
!!$    real(DP) :: w_q(lMax)
!!$    real(DP), allocatable  :: Clamb(:,:)
!!$    real(DP), allocatable :: hlamb(:)
!!$
!!$    NLat = size(alat)
!!$
!!$    M = 8 !6 !64/8
!!$    Mint = 40
!!$    lamb = 8.0 !6.0 !0.5 !8.0 !0.5 !1
!!$    allocate( xi(0:Mint), int_weight(0:Mint), q_intp(0:Mint), g(0:M) )
!!$    allocate( lon_intp(0:Mint), lat_intp(0:Mint), Clamb(0:M,0:Mint), hlamb(0:M) )
!!$    allocate( q_rc(JA), xi_rc(JA), Clamb_rc(0:M,JA) )
!!$
!!$    call calc_qaudrature(xi, int_weight, Mint+1)
!!$
!!$    a = alat(1)
!!$    b = alat(NLat)    
!!$    eps = 0.5d0*(b - a)
!!$    del = 0.5d0*(b + a)
!!$    w_q(:) = w_xy(xy_q(IS:IE,JS:JE))
!!$
!!$    do i = 0, Mint
!!$       lon_intp(i) = 0d0
!!$       lat_intp(i) = eps * xi(i)  + del
!!$       
!!$       q_intp(i) = Interpolate_w( w_q, lon_intp(i), lat_intp(i) )
!!$
!!$       Clamb(0,i) = 1d0; 
!!$       Clamb(1,i) = 2d0*lamb*xi(i)
!!$       do l=2, M
!!$          Clamb(l,i) = (2d0*xi(i)*(l + lamb - 1d0)*Clamb(l-1,i) - (l + 2d0*lamb - 2d0)*Clamb(l-2,i))/dble(l)
!!$       end do
!!$       write(*,*) "i=",i, ": Clamb", Clamb(:,i), "xi", xi(i), "intw=", int_weight(i)
!!$    end do
!!$    
!!$    do i=0,M       
!!$       hlamb(i) = sqrt(PI)* (gamma(i+2.0*lamb)/(gamma(i+1.0)*gamma(2.0*lamb)))*gamma(lamb + 0.5)/(dble(lamb + i)*gamma(lamb))
!!$!       write(*,* "Clamb_1:", gamma(real(2.0*lamb+
!       write(*,*) "lamb:", gamma(real(2.0*lamb)), gamma(real(lamb)), gamma(real(lamb+0.5)), "Clambn1", Clamb(i,M)
!!$    end do
!!$    write(*,*) "hlamb=", hlamb(:)
!!$    
!!$    g(:) = 0d0
!!$    do l=0, M
!!$       do i=0, Mint
!!$          g(l) = g(l) + &
!!$               &   int_weight(i) * (1d0 - xi(i)**2)**(lamb - 0.5d0) &
!!$               & * Clamb(l,i)                                       &
!!$               & * q_intp(i)
!!$          if(l==0) write(*,*) "i=", g(l)
!!$       end do
!!$       g(l) = g(l)/hlamb(l)
!!$    end do
!!$
!!$    q_rc(:) = 0d0
!!$    do j=JS, JE
!!$       xi_rc(j) = (xy_Lat(IS,j) - del)/eps
!!$       if ( abs(xi_rc(j)) <= 1d0 ) then
!!$          Clamb_rc(0,j) = 1d0
!!$          Clamb_rc(1,j) = 2d0*lamb*xi_rc(j)       
!!$          do l=2, M
!!$             Clamb_rc(l,j) = (2d0*xi_rc(j)*(l + lamb - 1d0)*Clamb_rc(l-1,j) - (l + 2d0*lamb - 2d0)*Clamb_rc(l-2,j))/dble(l)
!!$          end do
!          write(*,*) "j=", j , Clamb_rc(:,j), "xi_rc", xi_rc(j)
!!$          q_rc(j) = sum(g(:)*Clamb_rc(:,j))
!!$          write(*,*) j, xy_q(IS,j), q_rc(j)
!!$       end if
!!$    end do
!!$
!!$    write(*,*) "q_intp", q_intp
!!$    write(*,*) "g", g(:)
!    write(*,*) "q_rc", q_rc(JS:JE)
!    write(*,*) "q_ori", xy_q(IS,JS:JE)
!!$    stop
!!$       
!!$  end subroutine gegenbauer_reconstruct
!!$
!!$  subroutine freud_reconstruct( &
!!$       xy_q_rc, xy_q, alat, M, lamb, Mint )
!!$
!!$    use w_zonal_module_sjpack, only: Interpolate_w
!!$    use DOGCM_Admin_Grid_mod, only: xy_Lat
!!$    
!!$    real(DP), intent(out) :: xy_q_rc(IA,JA)
!!$    real(DP), intent(in) :: xy_q(IA,JA)
!!$    real(DP), intent(in) :: alat(2)
!!$    integer, intent(in) :: lamb
!!$    integer, intent(in) :: M
!!$    integer, intent(in) :: Mint
!!$
!!$    real(DP) :: a
!!$    real(DP) :: b
!!$    real(DP) :: eps
!!$    real(DP) :: del
!!$    real(DP) :: alpha
!!$    real(DP) :: beta
!!$
!!$    integer :: NLat
!!$    integer :: i
!!$    integer :: l
!!$    integer :: j
!!$    
!!$    real(DP), allocatable :: xi(:)
!!$    real(DP), allocatable :: g(:)
!!$    real(DP), allocatable :: int_weight(:)
!!$    real(DP), allocatable :: lon_intp(:)
!!$    real(DP), allocatable :: lat_intp(:)
!!$    real(DP), allocatable :: q_intp(:)
!!$    real(DP), parameter :: PI = acos(-1d0)
!!$    real(DP), allocatable :: xi_rc(:)
!!$    real(DP), allocatable :: Clamb_rc(:,:)
!!$
!!$    real(DP) :: w_q(lMax)
!!$    real(DP), allocatable  :: Clamb(:,:)
!!$    real(DP), allocatable :: hlamb(:)
!!$    real(DP), allocatable :: weightFunc(:)
!!$    real(DP) :: c
!!$    
!!$    NLat = size(alat)
!!$
!!$    c = 32d0
!!$    allocate( xi(0:Mint), int_weight(0:Mint), q_intp(0:Mint), weightFunc(0:Mint), g(0:M) )
!!$    allocate( lon_intp(0:Mint), lat_intp(0:Mint), Clamb(0:M,0:Mint), hlamb(0:M) )
!!$    allocate( xi_rc(JA), Clamb_rc(0:M,JA) )
!!$
!!$    call calc_qaudrature(xi, int_weight, Mint+1)
!!$
!!$    a = alat(1)
!!$    b = alat(NLat)    
!!$    eps = 0.5d0*(b - a)
!!$    del = 0.5d0*(b + a)
!!$    w_q(:) = w_xy(xy_q(IS:IE,JS:JE))
!!$    hlamb(:) = 0d0
!!$
!!$    do i=0, Mint
!!$       lon_intp(i) = 0d0
!!$       lat_intp(i) = eps * xi(i)  + del
!!$
!!$       q_intp(i) = Interpolate_w( w_q, lon_intp(i), lat_intp(i) )
!!$    end do
!!$
!!$    !-----------
!!$
!!$    weightFunc(:) = exp(-c*xi(:)**(2d0*lamb))
!!$    Clamb(0,:) = 1d0
!!$    Clamb(1,:) = xi(:)
!!$    do l=0, 1
!!$       hlamb(l) = sum( int_weight(:) * Clamb(l,:)**2 * weightFunc(:))
!!$    end do
!!$
!!$    do l=2, M
!!$       Clamb(l,:) = xi(:)*Clamb(l-1,:) - hlamb(l-1)/hlamb(l-2)*Clamb(l-2,:)
!!$       hlamb(l) = sum( int_weight(:) * Clamb(l,:)**2 * weightFunc )
!!$    end do
!!$
!    write(*,*) "xi:", xi(:)
!    do l=0, M
!       write(*,*) "l=", l, "Clamb:", Clamb(l,:), "hlamb", hlamb(l)
!    end do
!!$    !-----------
!!$        
!!$    do l=0, M
!!$       g(l) = sum( int_weight(:)*weightFunc(:)*Clamb(l,:)*q_intp(:) )/hlamb(l)
!!$    end do
!!$
!!$    !-------------
!!$    
!!$    do j=JS, JE
!!$       xi_rc(j) = (xy_Lat(IS,j) - del)/eps
!!$       if ( abs(xi_rc(j)) <= 1d0 ) then
!!$          Clamb_rc(0,j) = 1d0
!!$          Clamb_rc(1,j) = xi_rc(j)
!!$          do l=2, M
!!$             Clamb_rc(l,j) = xi_rc(j)*Clamb_rc(l-1,j) - hlamb(l-1)/hlamb(l-2)*Clamb_rc(l-2,j)
!!$          end do
!!$          xy_q_rc(IS,j) = sum(g(:)*Clamb_rc(:,j))
!          write(*,*) xy_Lat(IS,j), xy_q(IS,j), xy_q(IS,j)
!!$       end if
!!$    end do
!!$    
!    write(*,*) "q_intp", q_intp
!    write(*,*) "g", g(:)
!!$
!!$  end subroutine freud_reconstruct
!!$  
!!$  subroutine calc_qaudrature(xi, weight, M)
!!$    integer, intent(in) :: M
!!$    real(DP), intent(out) :: xi(M)
!!$    real(DP), intent(out) :: weight(M)
!!$
!!$    integer :: i
!!$    integer :: n
!!$    integer :: n_
!!$    real(DP) :: p, p1, p2, p3
!!$    real(DP) :: p_nm1, p_np1, p_nm2
!!$    real(DP), parameter :: PI = acos(-1d0)
!!$    real(DP) :: del_xi
!!$
!!$    xi(1) = -1d0
!!$    xi(M) =  1d0
!!$    weight(1) = 2d0/dble(M*(M-1))
!!$    weight(M) = weight(1)
!!$
!!$    do i=2, M-1
!!$       xi(i) = -(1d0 - 3d0*(M-2d0)/(8d0*(M-1d0)**3))*cos((4d0*i-3d0)/(4d0*(M-1d0)+1d0)*PI)
!!$       del_xi = 1d2
!!$       do while(abs(del_xi) > 1d-10)
!!$          p = xi(i); p_nm1 = 1d0
!!$          do n=1, M-1
!!$             p_np1 = ((2*n + 1)*xi(i)*p - n*p_nm1)/dble(n + 1d0)
!!$             p_nm2 = p_nm1; p_nm1 = p; p = p_np1; 
!!$          end do
!!$          n_ = M-1
!!$          p1 = n_*(p_nm2 - xi(i)*p_nm1)/(1d0 - xi(i)**2)
!!$          p2 = (2d0*xi(i)*p1 - n_*(n_+1)*p_nm1)/(1d0 - xi(i)**2)
!!$          p3 = (2d0*xi(i)*p2 - (n_*(n_+1)-2)*p1)/(1d0 - xi(i)**2)
!!$          del_xi = - 2d0*p1*p2/(2d0*p2**2 - p1*p3)
!!$          xi(i) = xi(i) + del_xi
!!$       end do
!!$       weight(i) = 2d0/dble(M*(M-1)*p_nm1**2)
!!$    end do
!!$
!!$    write(*,*) "xi", xi(:)
!!$    write(*,*) "w", weight(:)
!!$  end subroutine calc_qaudrature
  
end module DOGCM_Phys_hspm_vfvm_mod

