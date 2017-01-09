!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to solve governing (partial differential) equation with line method 
!!
!! First we evaluate RHS of governing equation with some spatial discrete scheme, which obtains the
!! time tendency of prognostic variables. Then calculate values of prognostic variables at next time
!! step using specified temporal scheme.
!!
!! @author Yuta Kawai
!!
!!
module DOGCM_TInt_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_GovernEq_mod, only: &
       & DynEqType, EOSType, &
       & OCNGOVERNEQ_DYN_HYDROBOUSSINESQ, &
       & OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM,       &
       & xyz_Z, xy_Topo,       &
       & xyz_Lat

  use DOGCM_Admin_Constants_mod, only: &
       & Grav, RefDens, Omega
  
  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final, &
       & EOSDriver_Eval

  use TemporalIntegUtil_mod2, only: &
       & TimeIntMode_Euler, &
       & TimeIntMode_RK2,   &
       & TimeIntMode_LF,    &
       & TimeIntMode_LFAM3
       
  
  use DOGCM_Admin_TInteg_mod, only: &
       & TIMELV_ID_N, TIMELV_ID_B, TIMELV_ID_A, &
       & DelTime,                               &
       & BarocTimeIntMode,                      &
       & VDiffTermACoef, CoriolisTermACoef
  
  use DOGCM_Admin_Variable_mod, only: &
       & xya_SSH, xyza_H, xyza_U, xyza_V, xyza_OMG, xyzaa_TRC, &
       & xyza_HydPres, xya_SfcPres,                            &
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,        &
       & TRC_TOT_NUM,                                          &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Boundary_driver_mod, only: &
       & DOGCM_Boundary_driver_ApplyBC
  
  use DOGCM_Dyn_driver_mod, only: &
       & DOGCM_Dyn_driver_Init, DOGCM_Dyn_driver_Final
       
  use DOGCM_Phys_driver_mod, only: &
       & DOGCM_Phys_driver_Init, DOGCM_Phys_driver_Final

  use DOGCM_TInt_common_mod, only: &
       & DOGCM_TInt_common_Init, DOGCM_TInt_common_Final, &
       & DOGCM_TInt_common_advance_Dyn,                   &
       & DOGCM_TInt_common_advance_Phys

  use ProfUtil_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_TInt_driver_Init, DOGCM_TInt_driver_Final
  public :: DOGCM_TInt_driver_Do
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_TInt_driver_mod' !< Module Name

contains

  !>
  !!
  !!
  Subroutine DOGCM_TInt_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)

    call EOSDriver_Init( EOSType )
    call DOGCM_Dyn_driver_Init( configNmlName )
    call DOGCM_Phys_driver_Init( configNmlName )

    call DOGCM_TInt_common_Init( configNmlName )
    
  end subroutine DOGCM_TInt_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_TInt_driver_Final()

    ! 実行文; Executable statements
    !

    call DOGCM_TInt_common_Final()
    
    call EOSDriver_Final()
    call DOGCM_Dyn_driver_Final()
    call DOGCM_Phys_driver_Final()
    
  end subroutine DOGCM_TInt_driver_Final

  !-------------------------------------

  
  subroutine DOGCM_TInt_driver_Do( isSelfStartSchemeUsed )

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

!!$    call MessageNotify('M', module_name, "Do ..")

    !-- Advance time step ---------------------------------------------------------
    
    if ( isSelfStartSchemeUsed ) then
          call DOGCM_TInt_Euler_Do()
    else

       select case (BarocTimeIntMode)
       case (TimeIntMode_Euler)
          call DOGCM_TInt_Euler_Do()
       case (TimeIntMode_RK2)
          call DOGCM_TInt_RK2_Do()
       case (TimeIntMode_LF)
          call DOGCM_TInt_LF_Do()
       case (TimeIntMode_LFAM3)
          call DOGCM_TInt_LFAM3_Do()
       end select
       
    end if

    
  end subroutine DOGCM_TInt_driver_Do

  !-----------------------------------------

  subroutine DOGCM_TInt_Euler_Do()

    ! 作業変数
    ! Work variables
    !
    real(DP) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)    
    integer :: TLA
    integer :: TLN

    ! 実行文; Executable statement
    !

!!$    call MessageNotify('M', module_name, "TInt = Euler ..")
    
    TLA = TIMELV_ID_A
    TLN = TIMELV_ID_N    

!!$    call MessageNotify('M', module_name, "OCN_BC  ..")    
    call DOGCM_Boundary_driver_ApplyBC( &
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (inout)
         & xyza_H(:,:,:,TLN), xyz_VViscCoef, xyz_VDiffCoef                  & ! (in)
         & )

!!$    call MessageNotify('M', module_name, "OCN_Phys  ..")
    call DOGCM_TInt_common_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                  & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                   & ! (out)
       & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN),                            & ! (in)
       & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),     & ! (in)
       & xyz_Z, xy_Topo,                                                  & ! (in)
       & DelTime                                                          & ! (in)
       & )

!!$    call MessageNotify('M', module_name, "OCN_Dyn  ..")    
    call DOGCM_TInt_common_advance_Dyn(  &
         ! ------ Time lelve A ------------------------------------------------- 
         & xyza_U(:,:,:,TLA), xyza_V(:,:,:,TLA), xyza_OMG(:,:,:,TLA),      & ! (out)
         & xyza_H(:,:,:,TLA), xya_SSH(:,:,TLA), xyzaa_TRC(:,:,:,:,TLA),    & ! (out)
         & xya_SfcPres(:,:,TLA), xyza_HydPres(:,:,:,TLA),                  & ! (out)
         ! ------ Time lelve Dyn -------------------------------------------------
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyza_OMG(:,:,:,TLN),      & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (in)
         & xya_SfcPres(:,:,TLN), xyza_HydPres(:,:,:,TLN),                  & ! (in)
         ! ------ Time lelve 0 --------------------------------------------------
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyza_OMG(:,:,:,TLN),      & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (in)
         & xya_SfcPres(:,:,TLN), xyza_HydPres(:,:,:,TLN),                  & ! (in)
         ! ----------------------------------------------------------------------
         & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                 & ! (in)
         & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                  & ! (in)
         & DelTime,                                                        & ! (in)
         ! -- Note that when alpha=0.5 and gamma=1.0 in this temporal scheme,
         !    Crank-Nicolson method is applied for the Coriolis and surface pressure terms. 
         & alpha=0.5d0, gamma=1d0, lambda=VDiffTermACoef                   & ! (in)
         & )    


  end subroutine DOGCM_TInt_Euler_Do

  !-----------------------------------------------
  
  subroutine DOGCM_TInt_RK2_Do()

    ! 実行文; Executable statement
    !    
    call MessageNotify('M', module_name, "TInt = RK2 ..")
        
!!$    call DOGCM_TInt_common_advance( &
!!$         & TL_A    = TIMELV_ID_A,   &
!!$         & TL_0    = TIMELV_ID_N,   &
!!$         & TL_DYN  = TIMELV_ID_N,   &
!!$         & TL_PHYS = TIMELV_ID_N,   &
!!$         & dt      = DelTime        &
!!$         & )    
    
  end subroutine DOGCM_TInt_RK2_Do

  !-----------------------------------------------
  
  subroutine DOGCM_TInt_LF_Do()

    ! 作業変数
    ! Work variables
    !    
    real(DP), parameter :: filterCoef = 1d-3
    real(DP) :: coef1, coef2
    
    real(DP) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    integer :: TLA, TLN, TLB

    ! 実行文; Executable statement
    !
    
!!$    call MessageNotify('M', module_name, "TInt = LF..")
    
    TLA = TIMELV_ID_A
    TLN = TIMELV_ID_N
    TLB = TIMELV_ID_B

!!$    call DOGCM_Boundary_driver_ApplyBC( &
!!$         & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),   & ! (inout)
!!$         & xyza_H(:,:,:,TLB), xyz_VViscCoef, xyz_VDiffCoef                 & ! (in)
!!$         & )    
!!$    call DOGCM_Boundary_driver_ApplyBC( &
!!$         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),   & ! (inout)
!!$         & xyza_H(:,:,:,TLN), xyz_VViscCoef, xyz_VDiffCoef                 & ! (in)
!!$         & )
    
    call DOGCM_TInt_common_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                  & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                   & ! (out)
       & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB),                            & ! (in)
       & xyza_H(:,:,:,TLB), xya_SSH(:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),     & ! (in)
       & xyz_Z, xy_Topo,                                                  & ! (in)
       & 2d0*DelTime                                                      & ! (in)
       & )
    
    call DOGCM_TInt_common_advance_Dyn(  &
         ! ------ Time lelve A ------------------------------------------------- 
         & xyza_U(:,:,:,TLA), xyza_V(:,:,:,TLA), xyza_OMG(:,:,:,TLA),      & ! (out)
         & xyza_H(:,:,:,TLA), xya_SSH(:,:,TLA), xyzaa_TRC(:,:,:,:,TLA),    & ! (out)
         & xya_SfcPres(:,:,TLA), xyza_HydPres(:,:,:,TLA),                  & ! (out)
         ! ------ Time lelve Dyn -------------------------------------------------
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyza_OMG(:,:,:,TLN),      & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (in)
         & xya_SfcPres(:,:,TLN), xyza_HydPres(:,:,:,TLN),                  & ! (in)
         ! ------ Time lelve 0 --------------------------------------------------
         & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB), xyza_OMG(:,:,:,TLB),      & ! (in)
         & xyza_H(:,:,:,TLB), xya_SSH(:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),    & ! (in)
         & xya_SfcPres(:,:,TLB), xyza_HydPres(:,:,:,TLB),                  & ! (in)
         ! ----------------------------------------------------------------------
         & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                 & ! (in)
         & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                  & ! (in)
         & 2d0*DelTime,                                                    & ! (in)
         & alpha=1d0/3d0, gamma=1d0/3d0, lambda=VDiffTermACoef             & ! (in)
!         & alpha=CoriolisTermACoef, gamma=1d0-2d0*CoriolisTermACoef, lambda=VDiffTermACoef                   & ! (in)
         & )    

    
    ! Time filter by Asselin (1972)
    !
        
    coef1 = 1d0 - 2d0*filterCoef
    coef2 = filterCoef

    !$omp parallel
    !$omp workshare
    xyza_U(:,:,:,TLN) = coef1*xyza_U(:,:,:,TLN) + coef2*(xyza_U(:,:,:,TLA) + xyza_U(:,:,:,TLB))
    
    xyza_V(:,:,:,TLN) = coef1*xyza_V(:,:,:,TLN) + coef2*(xyza_V(:,:,:,TLA) + xyza_V(:,:,:,TLB))
    
    xyza_H(:,:,:,TLN) = coef1*xyza_H(:,:,:,TLN) + coef2*(xyza_H(:,:,:,TLA) + xyza_H(:,:,:,TLB))
    
    xyzaa_TRC(:,:,:,:,TLN) = coef1*xyzaa_TRC(:,:,:,:,TLN) + coef2*(xyzaa_TRC(:,:,:,:,TLA) + xyzaa_TRC(:,:,:,:,TLB))
    
    xya_SSH(:,:,TLN) = coef1*xya_SSH(:,:,TLN) + coef2*(xya_SSH(:,:,TLA) + xya_SSH(:,:,TLB))
    !$omp end workshare
    !$omp end parallel
    
  end subroutine DOGCM_TInt_LF_Do

  subroutine DOGCM_TInt_LFAM3_Do()

    use SpmlUtil_mod
    use VFvmUtil_mod
    use DOGCM_Boundary_Vars_mod
    
    ! 作業変数
    ! Work variables
    !
    real(DP) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)

    integer :: TLA
    integer :: TLN
    integer :: TLB

    real(DP) :: xyz_U(IA,JA,KA)
    real(DP) :: xyz_V(IA,JA,KA)
    real(DP) :: xyz_OMG(IA,JA,KA)
    real(DP) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP) :: xyz_H(IA,JA,KA)
    real(DP) :: xy_SSH(IA,JA)
    real(DP) :: xy_SfcPres(IA,JA)

    real(DP) :: avr_ptempN
    real(DP) :: avr_ptempA
    real(DP) :: avr_SfcHFlx
    real(DP) :: avr_ptemp_RHS_phys

    integer :: i
    integer :: j
    integer :: k
    integer :: n
    
    ! 実行文; Executable statement
    !
    
!!$    call MessageNotify('M', module_name, "TInt = LFAM3 ..")
    
    TLA = TIMELV_ID_A
    TLN = TIMELV_ID_N
    TLB = TIMELV_ID_B

    !
    !
    
!!$    call DOGCM_Boundary_driver_ApplyBC( &
!!$         & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),   & ! (inout)
!!$         & xyza_H(:,:,:,TLB), xyz_VViscCoef, xyz_VDiffCoef                 & ! (in)
!!$         & )
!!$

!!$    call DOGCM_Boundary_driver_ApplyBC( &
!!$         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),   & ! (inout)
!!$         & xyza_H(:,:,:,TLN), xyz_VViscCoef, xyz_VDiffCoef                 & ! (in)
!!$         & )
    
!!$    call MessageNotify('M', module_name, "OCN_Phys [1/2] ..")
    call DOGCM_TInt_common_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                  & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                   & ! (out)
       & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB),                            & ! (in)
       & xyza_H(:,:,:,TLB), xya_SSH(:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),     & ! (in)
       & xyz_Z, xy_Topo,                                                  & ! (in)
       & 2d0*DelTime                                                      & ! (in)
       & )

!!$    call MessageNotify('M', module_name, "OCN_Dyn [1/2] ..")    
    call DOGCM_TInt_common_advance_Dyn(  &
         ! ------ Time lelve A ------------------------------------------------- 
         & xyza_U(:,:,:,TLA), xyza_V(:,:,:,TLA), xyza_OMG(:,:,:,TLA),      & ! (out)
         & xyza_H(:,:,:,TLA), xya_SSH(:,:,TLA), xyzaa_TRC(:,:,:,:,TLA),    & ! (out)
         & xya_SfcPres(:,:,TLA), xyza_HydPres(:,:,:,TLA),                  & ! (out)
         ! ------ Time lelve Dyn -------------------------------------------------
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyza_OMG(:,:,:,TLN),      & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (in)
         & xya_SfcPres(:,:,TLN), xyza_HydPres(:,:,:,TLN),                  & ! (in)
         ! ------ Time lelve 0 --------------------------------------------------
         & xyza_U(:,:,:,TLB), xyza_V(:,:,:,TLB), xyza_OMG(:,:,:,TLB),      & ! (in)
         & xyza_H(:,:,:,TLB), xya_SSH(:,:,TLB), xyzaa_TRC(:,:,:,:,TLB),    & ! (in)
         & xya_SfcPres(:,:,TLB), xyza_HydPres(:,:,:,TLB),                  & ! (in)
         ! ----------------------------------------------------------------------
         & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                 & ! (in)
         & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                  & ! (in)
         & 2d0*DelTime,                                                    & ! (in)
!!$         & alpha=0.5d0, gamma=0d0, lambda=VDiffTermACoef                   & ! (in)    
         & alpha=CoriolisTermACoef, gamma=1d0-2d0*CoriolisTermACoef, lambda=VDiffTermACoef & ! (in)    
         & )    

    ! Calculate the values interpolated with 3rd-order Adams-Moulton formula.
    !
    
    call ProfUtil_RapStart('OcnAM3', 1)

    !$omp parallel private(i,j,k,n)
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IE
       xyz_U(i,j,k) = (5d0*xyza_U(i,j,k,TLA) + 8d0*xyza_U(i,j,k,TLN) - xyza_U(i,j,k,TLB))/12d0
       xyz_V(i,j,k) = (5d0*xyza_V(i,j,k,TLA) + 8d0*xyza_V(i,j,k,TLN) - xyza_V(i,j,k,TLB))/12d0
       xyz_H(i,j,k) = (5d0*xyza_H(i,j,k,TLA) + 8d0*xyza_H(i,j,k,TLN) - xyza_H(i,j,k,TLB))/12d0
    end do
    end do
    end do
    !$omp do collapse(3)
    do n=1, TRC_TOT_NUM
    do k=KS, KE
    do j=JS, JE
    do i=IS, IE
       xyza_TRC(i,j,k,n) = (5d0*xyzaa_TRC(i,j,k,n,TLA) + 8d0*xyzaa_TRC(i,j,k,n,TLN) - xyzaa_TRC(i,j,k,n,TLB))/12d0
    end do
    end do
    end do
    end do
    !$omp do collapse(1)
    do j=JS, JE
    do i=IS, IE
       xy_SSH(i,j)     = (5d0*xya_SSH(i,j,TLA) + 8d0*xya_SSH(i,j,TLN) - xya_SSH(i,j,TLB))/12d0
       xy_SfcPres(i,j) = (5d0*xya_SfcPres(i,j,TLA) + 8d0*xya_SfcPres(i,j,TLN) - xya_SfcPres(i,j,TLB))/12d0    
    end do
    end do
    !$omp end parallel
    
    call ProfUtil_RapEnd('OcnAM3', 1)

    !
    !
!!$    call MessageNotify('M', module_name, "OCN_Phys [2/2] ..")

    call DOGCM_TInt_common_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                  & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                   & ! (out)
       & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN),                            & ! (in)
       & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),     & ! (in)
       & xyz_Z, xy_Topo,                                                  & ! (in)
       & DelTime, .true.                                                  & ! (in)
       & )

!!$    call MessageNotify('M', module_name, "OCN_Dyn [2/2] ..")        
    call DOGCM_TInt_common_advance_Dyn(  &
         ! ------ Time lelve A ------------------------------------------------- 
         & xyza_U(:,:,:,TLA), xyza_V(:,:,:,TLA), xyza_OMG(:,:,:,TLA),      & ! (out)
         & xyza_H(:,:,:,TLA), xya_SSH(:,:,TLA), xyzaa_TRC(:,:,:,:,TLA),    & ! (out)
         & xya_SfcPres(:,:,TLA), xyza_HydPres(:,:,:,TLA),                  & ! (out)
         ! ------ Time lelve Dyn -------------------------------------------------
         & xyz_U, xyz_V, xyz_OMG,                                          & ! (in)
         & xyz_H, xy_SSH, xyza_TRC,                                        & ! (in)
         & xy_SfcPres, xyza_HydPres(:,:,:,TLN),                            & ! (in)
         ! ------ Time lelve 0 --------------------------------------------------
         & xyza_U(:,:,:,TLN), xyza_V(:,:,:,TLN), xyza_OMG(:,:,:,TLN),      & ! (in)
         & xyza_H(:,:,:,TLN), xya_SSH(:,:,TLN), xyzaa_TRC(:,:,:,:,TLN),    & ! (in)
         & xya_SfcPres(:,:,TLN), xyza_HydPres(:,:,:,TLN),                  & ! (in)
         ! ----------------------------------------------------------------------
         & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,                 & ! (in)
         & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,                  & ! (in)
         & DelTime,                                                        & ! (in)
!!$         & alpha=0.5d0, gamma=0d0, lambda=VDiffTermACoef,                  & ! (in)    
         & alpha=CoriolisTermACoef, gamma=1d0-2d0*CoriolisTermACoef, lambda=VDiffTermACoef, & ! (in)    
         & lhst_tend=.true. )    
    
!!$    avr_ptempN = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyzaa_TRC(IS:IE,JS:JE,KS:KE,TRCID_SALT,TLN)) )    
!!$    avr_ptempA = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyzaa_TRC(IS:IE,JS:JE,KS:KE,TRCID_SALT,TLA)) )
!!$    avr_ptempN = AvrLonLat_xy( VFvm_Int_BtmToTop(  (xyzaa_TRC(IS:IE,JS:JE,:,TRCID_SALT,TLN)), xyz_H(IS:IE,JS:JE,:)))
!!$    avr_ptempA = AvrLonLat_xy( VFvm_Int_BtmToTop(  (xyzaa_TRC(IS:IE,JS:JE,:,TRCID_SALT,TLA)), xyz_H(IS:IE,JS:JE,:)))
!!$    write(*,*) "FreshWt:", AvrLonLat_xy(xy_FreshWtFlxS(IS:IE,JS:JE))*1d3
!!$    write(*,*) "Salt:", (avr_ptempA - avr_ptempN)/deltime, -35d0*AvrLonLat_xy(xy_FreshWtFlxS(IS:IE,JS:JE))
!!$    avr_ptemp_RHS_phys = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( xyza_TRC_RHS_phy(IS:IE,JS:JE,KS:KE, TRCID_PTEMP) ))
!!$    
!!$    avr_SfcHFlx = AvrLonLat_xy( (xy_SfcHFlx_ns(IS:IE,JS:JE) + xy_SfcHFlx_sr(IS:IE,JS:JE) ) )*DelTime/(RefDens*Cp0)
!!$    write(*,*) "AvrPTempN=", avr_ptempN, ", AvrPTempA=", avr_ptempA, &
!!$         & ", avr_SfcHFlx*dt/(Rho0*Cp0)=", avr_SfcHFlx, ", DelAvrPTemp=", (avr_ptempA - avr_ptempN)*5.2d3, &
!!$         & ", RHS_phys=", avr_ptemp_RHS_phys*5.2d3*DelTime

  end subroutine DOGCM_TInt_LFAM3_Do
  
end module DOGCM_TInt_driver_mod
