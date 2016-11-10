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

  !-- Module for the parametrizations of lateral oceanic physics

  use LPhys_DIFF_spm_mod, only: &
       & LPhys_DIFF_spm_Init, LPhys_DIFF_spm_Final, &
       & LPhys_DIFF_spm_LMixMOMRHS,                       &
       & LPhys_DIFF_spm_LMixMOMRHSImpl,                   &
       & LPhys_DIFF_spm_LMixTRCRHS,                       &
       & LPhys_DIFF_spm_LMixTRCRHSImpl
       
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

  public :: DOGCM_Phys_hspm_vfvm_VImplUV
  public :: DOGCM_Phys_hspm_vfvm_VImplTRC
  
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
    end if
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_hspm_vfvm_Init( confignmlFileName = configNmlName )
    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_Init()
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
    real(DP), intent(inout) :: xyz_U_RHS_phy(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyz_V_RHS_phy(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyza_TRC_RHS_phy(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(inout) :: xyz_VViscCoef(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyz_VDiffCoef(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xy_TOPO(0:iMax-1,jMax)
    real(DP), intent(in) :: dt


    real(DP) :: avr_ptemp_RHS_phys

!!$    real(DP) :: xyz_UA(0:iMax-1,jMax,KA)
!!$    real(DP) :: xyz_VA(0:iMax-1,jMax,KA)
!!$    real(DP) :: xyza_TRCA(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    
    
    ! 実行文; Executable statements
    !

    !-- Horizontal momentum -----------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXMOM_NAME )  ) then
       call LPhys_DIFF_spm_LMixMOMRHSImpl( &
            & xyz_U_RHS_phy(:,:,KS:KE), xyz_V_RHS_phy(:,:,KS:KE),        & ! (inout)
            & xyz_U(:,:,KS:KE), xyz_V(:,:,KS:KE), xyz_H(:,:,KS:KE),      & ! (in)
            & hViscCoef, hHyperViscCoef,                                 & ! (in)
            & dt                                                         & ! (in)
            & )

    end if

    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXMOM_NAME )  ) then
       call DOGCM_Phys_hspm_vfvm_VMixMOMRHS(            & 
            & xyz_U_RHS_phy, xyz_V_RHS_phy,                              & ! (inout)
            & xyz_U, xyz_V, xyz_H, xyz_VViscCoef                         & ! (in)
            & )
    end if

    !-- Tracer ------------------------------------------------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then
       
       call LPhys_DIFF_spm_LMixTRCRHSImpl(            &
            & xyza_TRC_RHS_phy(:,:,KS:KE,:),                                       & ! (inout)
            & xyza_TRC(:,:,KS:KE,:), xyz_H(:,:,KS:KE), hDiffCoef, hHyperDiffCoef,  & ! (in)
            & dt )                                                                   ! (in)

!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "->avr_ptemp_phys (+ LMixRHS): ", avr_ptemp_RHS_phys
    end if


    !-- Lateral mixing of tracers by eddy induced velocity -------------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
       call LPhys_RediGM_hspm_vfvm_AddMixingTerm( &
            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP),                      & ! (inout)
            & xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),                       & ! (inout)
            & xyz_VDiffCoef,                                            & ! (inout)
            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),  & ! (in)
            & xyz_H, xyz_Z, xy_Topo                                     & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ RediGMRHS): ",  avr_ptemp_RHS_phys
    end if

    !-- Vertical mixing of tracers by non-penetrative convection -------------------
    
    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
       call DOGCM_VPhys_ConvAdjust_AddMixingTerm( &
            & xyza_TRC_RHS_phy(:,:,KS:KE,TRCID_PTEMP),                         & ! (inout)
            & xyza_TRC_RHS_phy(:,:,KS:KE,TRCID_SALT),                          & ! (inout)
            & xyz_ConvIndex(IS:IE,JS:JE,KS:KE),                                & ! (inout)
            & xyza_TRC(:,:,KS:KE,TRCID_PTEMP), xyza_TRC(:,:,KS:KE,TRCID_SALT), & ! (in)
            & xyz_Z(:,:,KS:KE), z_KAXIS_Weight(KS:KE), dt                      & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ ConvAdjustRHS): ",  avr_ptemp_RHS_phys       
    end if


    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
       call DOGCM_Phys_hspm_vfvm_VMixTRCRHS(            &  
            & xyza_TRC_RHS_phy,                                         & ! (inout)
            & xyza_TRC, xyz_H, xyz_VDiffCoef                            & ! (in)
            & )
!!$       avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(:,:,:, TRCID_SALT) ))
!!$       write(*,*) "avr_ptemp_phys (+ VMixRHS): ",  avr_ptemp_RHS_phys
!!$       write(*,*) "avr_ptemp_phys_local:", IntSig_BtmToTop(xyza_TRC_RHS_phy(0,jMax/2,:,TRCID_SALT))
    end if
    
!!$    call DOGCM_Phys_hspm_vfvm_VImplTRC( xyza_TRCA,         & ! (out)
!!$       & xyza_TRC, xyza_TRC_RHS_phy*5.2d3,                 & ! (in)
!!$       & xyz_H, xyz_H, xyz_VDiffCoef, dt, 0.5d0            & ! (in)
!!$       & )
!!$
!!$    xyza_TRC_RHS_phy = (xyza_TRCA - xyza_TRC)/dt
!!$    
!!$    write(*,*) "Phys Salt:", AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(:,:,:,TRCID_SALT), xyz_H) )

  end subroutine DOGCM_Phys_hspm_vfvm_Do

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_VMixMOMRHS(     &
       & xyz_U_RHS, xyz_V_RHS,                            & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_VViscCoef               & ! (in)
       & )

    use DOGCM_Admin_Constants_mod, only: &
         & RefDens    

    use DOGCM_Admin_BC_mod, only: &
         & DynBC_Surface, DynBC_Bottom, &
         & DynBCTYPE_NoSlip
    
    use DOGCM_Boundary_vars_mod, only: &
         & xy_WindStressU, xy_WindStressV
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(inout) :: xyz_V_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VViscCoef(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    real(DP) :: xyr_DiffFlxU(0:iMax-1,jMax,KA)
    real(DP) :: xyr_DiffFlxV(0:iMax-1,jMax,KA)
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
       xyr_DiffFlxU(:,:,KS-1) = xy_WindStressU(IS:IE,JS:JE)/RefDens
       xyr_DiffFlxV(:,:,KS-1) = xy_WindStressV(IS:IE,JS:JE)/RefDens
    end select

    select case(DynBC_Bottom)
    case(DynBCTYPE_NoSlip)
       !$omp parallel
       !$omp workshare
       xyr_DiffFlxU(:,:,KE  ) = z_RFDK(KE)*xyz_VViscCoef(:,:,KE)/xyz_H(:,:,KE)* 2d0*xyz_U(:,:,KE)
       xyr_DiffFlxV(:,:,KE  ) = z_RFDK(KE)*xyz_VViscCoef(:,:,KE)/xyz_H(:,:,KE)* 2d0*xyz_V(:,:,KE)
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
       & xyza_TRC_RHS,                      & ! (out)
       & xyza_TRC, xyz_H, xyz_VDiffCoef     & ! (in)
       & )

    use DOGCM_Admin_Constants_mod, only: &
         & RefDens, Cp0

    
    use DOGCM_Admin_BC_mod, only: &
         & ThermBC_Surface, ThermBC_Bottom, &
         & SaltBC_Surface, SaltBC_Bottom,   &
         & ThermBCTYPE_PrescTemp,           &
         & SaltBCTYPE_PrescSalt

    use DOGCM_Boundary_vars_mod, only: &
         & xy_SfcHFlx_sr, xy_SfcHFlx_ns, xy_SeaSfcTemp, &
         & xy_FreshWtFlxS, xy_SeaSfcSalt
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyza_TRC_RHS(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !    

    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    real(DP) :: xyra_DiffFlxTRC(0:iMax-1,jMax,KA,TRC_TOT_NUM)
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
            &                        * 2d0*(xyza_TRC(:,:,KS,n) - xy_SeaSfcTemp(IS:IE,JS:JE))
    case default
       xyra_DiffFlxTRC(:,:,KS-1,n) = - (xy_SfcHFlx_ns(IS:IE,JS:JE) + xy_SfcHFlx_sr(IS:IE,JS:JE)) &
            &                        / (RefDens * Cp0)
    end select

    ! for salinity    
    k = KS; n = TRCID_SALT
    select case(ThermBC_Surface)
    case(ThermBCTYPE_PrescTemp)
       xyra_DiffFlxTRC(:,:,KS-1,n) = - z_RFDK(k)*xyz_VDiffCoef(:,:,KS)/xyz_H(:,:,KS)             &
            &                        * 2d0*(xyza_TRC(:,:,KS,n) - xy_SeaSfcSalt(IS:IE,JS:JE))
    case default
       xyra_DiffFlxTRC(:,:,KS-1,n) = - xy_FreshWtFlxS(IS:IE,JS:JE) * 35d0 ! Virtual salinity flux (Set 35 psu as reference salinity)
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
  
  subroutine DOGCM_Phys_hspm_vfvm_VImplTRC( xyza_TRCA,         & ! (out)
       & xyza_TRC0, xyza_HTRC_RHS,                       & ! (in)
       & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, alpha        & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & ThermBC_Surface, ThermBC_Bottom, &
         & SaltBC_Surface, SaltBC_Bottom

    use DOGCM_Admin_Grid_mod, only: &
         & xyz_Z, z_KAXIS_Weight
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyz_ConvIndex
    
    use DOGCM_Boundary_hspm_vfvm_mod, only: &
         & DOGCM_Boundary_hspm_vfvm_InqVBCRHS_TRC

    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyza_TRCA(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC0(0:iMax-1,jMax,KA,TRC_TOT_NUM)    
    real(DP), intent(in) :: xyza_HTRC_RHS(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H0(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(0:iMax-1,jMax,KA)
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

    real(DP) :: xyza_TRC_RHS_Conv(0:iMax-1,jMax,KA,2)
    
    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_hspm_vfvm_InqVBCRHS_TRC( &
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

    
!!$    avr_TRC0_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRC0(:,:,:,TRCID_SALT), xyz_H0) )
!!$    avr_TRCA_phys = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRCA(:,:,:,TRCID_SALT), xyz_H0) )
!!$    avr_RHS = AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_HTRC_RHS(:,:,:,TRCID_SALT)/xyz_H0, xyz_H0) )
!!$    write(*,*) "Phys+Dyn, Phys+Dyn+Impl Salt:", &
!!$         & avr_RHS,                          &
!!$         & (avr_TRCA_phys - avr_TRC0_phys)/dt, & !/dt*dt
!!$         & AvrLonLat_xy( VFvm_Int_BtmToTop(xyza_TRCA(:,:,:,TRCID_SALT) - xyza_TRC0(:,:,:,TRCID_SALT), xyz_H0) )/dt
    
  end subroutine DOGCM_Phys_hspm_vfvm_VImplTRC

  !-----------------------------------------------------------------------
  
  subroutine DOGCM_Phys_hspm_vfvm_VImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,                & ! (in)
       & xyz_H, xyz_VViscCoef, dt, alpha                      & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    use DOGCM_Admin_BC_mod, only: &
         & inquire_VBCSpecType, &
         & DynBC_Surface, DynBC_Bottom
    
    use DOGCM_Boundary_hspm_vfvm_mod, only: &
         & DOGCM_Boundary_hspm_vfvm_InqVBCRHS_UV
    
    ! 宣言文; Declaration statement
    !      

    real(DP), intent(out) :: xyz_UA(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_VA(0:iMax-1,jMax,KA)        
    real(DP), intent(in) :: xyz_U0(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V0(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_U_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VViscCoef(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
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
    
    call DOGCM_Boundary_hspm_vfvm_InqVBCRHS_UV( &
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
    
  end subroutine DOGCM_Phys_hspm_vfvm_VImplUV

  !-----------------------------------------------------------------------

  subroutine calc_VDiffEq( xyz_qA,              &  ! (inout)
       & xyz_q0, xyz_q_RHS,                     &  ! (in)
       & xyz_Kv, xyz_H, dt, alpha,              &  ! (in)
       & SeaSfcBCType, xy_SeaSfcBCVal,          &  ! (in)
       & SeaBtmBCType, xy_SeaBtmBCVal,          &  ! (in)
       & qname                                  &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(out) :: xyz_qA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_q0(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_q_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Kv(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
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
    real(DP) :: AMat_save(KS:KE,KS:KE)
    
   
    ! 実行文; Executable statements
    !


    N = KE - KS + 1
    
    !$omp parallel do collapse(2) private(z_Kv, z_e3, Al, Ad, Au, b, info, k, r_DFlxcoef, coef, m1, m2)
    do j = 1, jMax
       do i = 0, iMax-1

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
             Ad(k) = 1d0 + alpha*coef(k)*(2d0*r_DFlxCoef(k) + r_DFlxCoef(k-1))
          case('N')
             Al(k) = - alpha*coef(k)*r_DFlxCoef(k-1)
             Ad(k) = 1d0 + alpha*coef(k)*(                  + r_DFlxCoef(k-1))
          end select
          
          call DGTSV(N, 1, Al(KS+1:KE), Ad(KS:KE), Au(KS:KE-1), b, N, info)

          xyz_qA(i,j,KS:KE) = xyz_q0(i,j,KS:KE) + b(KS:KE)          
       end do
    end do

  end subroutine calc_VDiffEq
  
end module DOGCM_Phys_hspm_vfvm_mod
