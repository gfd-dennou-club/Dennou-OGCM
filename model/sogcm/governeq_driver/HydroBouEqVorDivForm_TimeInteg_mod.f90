!-------------------------------------------------------------
! Copyright (c) 2013-2015 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBoudEqVorDivForm_TimeInteg_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING 

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use Constants_mod, only: &
       & PI, Omega, Grav, RPlanet, Cp0, &
       & hViscCoef, vViscCoef, &
       & hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef, &
       & hHyperDiffCoef, vHyperDiffCoef, &
       & RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use GovernEqSet_mod, only: &
       & GOVERNEQSET_PHYSICS_EDDYMIX_NAME, &
       & GOVERNEQSET_PHYSICS_CONVADJUST_NAME, &       
       & isPhysicsCompActivated

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType, &
       & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress, &
       & ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed, & 
       & SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed, & 
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom

  
!!$  use TemporalIntegUtil_mod, only: &
!!$       & TemporalIntegUtil_Init, TemporalIntegUtil_Final, &
!!$       & timeIntMode_Euler, xy_timeIntEuler, wt_timeIntEuler, &
!!$       & timeIntMode_LFTR, xy_timeIntLFTR, wt_timeIntLFTR, &
!!$       & timeIntMode_LFAM3, xy_timeIntLFAM3, wt_timeIntLFAM3, &
!!$       & timeIntMode_RK2, xy_timeIntRK2, wt_timeIntRK2, &
!!$       & timeIntMode_RK4, xy_timeIntRK4, wt_timeIntRK4, &
!!$       & TemporalIntegUtil_GetDDtCoef, &
!!$       & TemporalIntegUtil_SetDelTime

  use TemporalIntegUtil_mod2
  
  use TemporalIntegSet_mod, only: &
       & CurrentTimeStep, SubCycleNum, &
       & nShortTimeLevel, &
       & CoriolisTermACoef, VDiffTermACoef

  use HydroBouEqRHS_DynProc_mod, only: &
       & HydroBouEqRHS_DynProc_Init, HydroBouEqRHS_DynProc_Final

  use HydroBouEqSolverVImplProc_mod, only: &
       & HydroBouEqSolverVImplProc_Init, HydroBouEqSolverVImplProc_Final, &
       & HydroBouEqSolverVImplProc_Prepare, &
       & Advance_VImplicitProc_DeltaForm_xyz

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEq_TimeInteg_Init, HydroBouEq_TimeInteg_Final
  public :: HydroBouEqSolver_AdvanceTStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBoudEq_TimeInteg_mod' !< Module Name
  real(DP), dimension(:,:), allocatable :: xy_CosLat

contains

  !>
  !!
  !!
  subroutine HydroBouEq_TimeInteg_Init()

    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    call HydroBouEqRHS_DynProc_Init()
    call TemporalIntegUtil_Init(0d0)! iMax, jMax, kMax, lMax, tMax, 0d0 )
    call HydroBouEqSolverVImplProc_Final()
    
    allocate(xy_CosLat(0:iMax-1,jMax))
    xy_CosLat = cos(xyz_Lat(:,:,0))

  end subroutine HydroBouEq_TimeInteg_Init

  !>
  !!
  !!
  subroutine HydroBouEq_TimeInteg_Final()

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    call HydroBouEqRHS_DynProc_Final()
    call TemporalIntegUtil_Final()
    call HydroBouEqSolverVImplProc_Final()
    deallocate(xy_CosLat)

  end subroutine HydroBouEq_TimeInteg_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTStep( &
       & DelTime, timeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt &
       & )
    
    ! モジュール引用; Use statements
    !

    use VariableSet_mod, only: &
         & xyz_UB, xyz_UN, xyz_UA, &
         & xyz_VB, xyz_VN, xyz_VA, &
         & xyz_PTempEddB, xyz_PTempEddN, xyz_PTempEddA, &
         & xyz_SaltB, xyz_SaltN, xyz_SaltA, &
         & xy_SurfHeightB, xy_SurfHeightN, xy_SurfHeightA, &
         & xy_SurfPressB, xy_SurfPressN, xy_SurfPressA, &
         & xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic, &
         & xyz_VViscCoefA, xyz_VViscCoefB, xyz_VViscCoefN, &
         & xyz_VDiffCoefA, xyz_VDiffCoefB, xyz_vDiffCoefN

    use BoundaryCondO_mod, only: &
         & apply_VBoundaryCondO

    use DiagnoseUtil_mod, only: &
         & Diagnose_SigDot, Diagnose_HydroPressEdd, Diagnose_GeoPot

    use EOSDriver_mod, only: &
         & EOSDriver_Eval
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTime
    integer, intent(in) :: timeIntMode
    integer, intent(in) :: nStage_BarocTimeInt
    logical, intent(in) :: isVarBUsed_BarocTimeInt

    ! 局所変数
    ! Local variables
    !

    !
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) ::   &
         & xyz_U, xyz_V, xyz_Vor, xyz_Div, &
         & xyz_PTempEdd, xyz_Salt,   &
         & xyz_PressBaroc, xyz_DensEdd, xyz_GeoPot
    
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHeight, xy_totDepth

    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_ForceUBaroc, xy_ForceVBaroc, &
         & xy_UBarot, xy_VBarot, xy_SurfPress, &
         & xy_UBarotOld, xy_VBarotOld, xy_SurfPressOld

    real(DP), dimension(lMax,0:kMax) :: &
         & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt, &
         & wz_VorN, wz_DivN, wz_PTempEddN, wz_SaltN, &
         & wz_VorB, wz_DivB, wz_PTempEddB, wz_SaltB
         

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PTempBasic
    
    
    ! Work variables for some temporal schemes, such as Runge=Kutta scheme etc. 
    !
    real(DP), dimension(lMax,0:kMax) :: &
         & wz_Vor_RHSTmp, wz_Div_RHSTmp, wz_PTemp_RHSTmp, wz_Salt_RHSTmp
    real(DP), dimension(lMax) :: &
         & w_SurfHeight_RHSTmp

    !
    !
    real(DP), dimension(lMax,0:kMax) :: &
         & wz_Vor_RHSEx, wz_Div_RHSEx, wz_PTemp_RHSEx, wz_Salt_RHSEx
    real(DP), dimension(0:iMax-1,jMax) :: &
         & w_SurfHeight_RHSEx

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U_RHSPhys, xyz_V_RHSPhys, xyz_PTemp_RHSPhys, xyz_Salt_RHSPhys

    
    real(DP) :: DelTauBaroc, DelTauBarot
    
    !
    !
    integer :: Stage
    character(TOKEN) :: TIntType_SurfPressTerm
    
    real(DP) :: xyz_CosLat(0:iMax-1, jMax, 0:kMax)
    integer :: j, k
    real(DP) :: xyz_Tmp(0:iMax-1,jMax,0:kMax)
!!$    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_VDiffCoef, xyz_VViscCoef
         
    logical :: isCoriTermSemiImplicit
    logical :: isVImplicitProc
    
    ! 実行文; Executable statement
    !
    
    ! * Preparation 


    !
    xyz_CosLat(:,:,:) = spread(xy_CosLat, 3, kMax+1)
    forAll(k=0:kMax) xyz_PTempBasic(:,:,k) = z_PTempBasic(k)

    if (CoriolisTermACoef > 0d0) then
       isCoriTermSemiImplicit = .true.
    else
       isCoriTermSemiImplicit = .false.
    end if

    isVImplicitProc = .false.
    

    ! Set some variables at the time level N. 

    call calc_VViscDiffCoef( &
            & xyz_VViscCoefN, xyz_VDiffCoefN,                           & ! (out)
            & xyz_UN, xyz_VN, xyz_PTempBasic + xyz_PTempEddN,  xyz_SaltN  & ! (in)
!            & xyz_U, xyz_V, xyz_PTempBasic + xyz_PTempEdd,  xyz_Salt  & ! (in)
            & )
       
    

    !
    !
    call apply_VBoundaryCondO( &
         & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_SaltN,     & ! (inout)
         & xyz_VViscCoefN, xyz_VDiffCoefN            & ! (in)
         & )    

    !
    !
    !

    !$omp parallel workshare
    xyz_U(:,:,:) = xyz_UN; xyz_V(:,:,:) = xyz_VN;
    xyz_PTempEdd(:,:,:) = xyz_PTempEddN; xyz_Salt(:,:,:) = xyz_SaltN
    xy_SurfHeight(:,:) = xy_SurfHeightN
    xy_SurfPress(:,:) = xy_SurfPressN
    xy_SurfPressOld(:,:) = xy_SurfPressN
    !$omp end parallel workshare

    if(isVarBUsed_BarocTimeInt) then
       call UV2VorDiv(xyz_UB, xyz_VB,    & ! (in)
            & wz_VorB, wz_DivB           & ! (out)
            & )
       wz_PTempEddB(:,:) = wz_xyz(xyz_PTempEddB)
       wz_SaltB(:,:) = wz_xyz(xyz_SaltB)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call TemporalIntegUtil_SetDelTime(DelTime)
    
    !
    
    !
    do Stage=1, nStage_BarocTimeInt


       != Preprocess of each temporal loop =================================

       !
       call UV2VorDiv(xyz_U, xyz_V,                 & ! (in)
            & wz_Vor, wz_Div, xyz_Vor, xyz_Div      & ! (out)
            & )
       wz_PTempEdd(:,:) = wz_xyz(xyz_PTempEdd)
       wz_Salt(:,:) = wz_xyz(xyz_Salt)

       if(Stage==1) then
          wz_VorN(:,:) = wz_Vor; wz_DivN(:,:) = wz_Div;
          wz_PTempEddN(:,:) = wz_PTempEdd; wz_SaltN(:,:) = wz_Salt
       end if
       
       !
       DelTauBaroc = TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime
       call HydroBouEqSolverVImplProc_Prepare( &
            & DelTauBaroc, CoriolisTermACoef, VDiffTermACoef, VDiffTermACoef &  ! (in)
            & )
       xy_ForceUBaroc(:,:) = 0d0; xy_ForceVBaroc(:,:) = 0d0
       xy_UBarot(:,:) = xy_IntSig_BtmToTop_xyz(xyz_U)
       xy_VBarot(:,:) = xy_IntSig_BtmToTop_xyz(xyz_V)

       !
       xy_totDepth(:,:) = xy_totDepthBasic + xy_SurfHeight
       
       
       != Core part of each temporal loop ====================================


       ! Set diagnosed variables
       !
       xyz_SigDot(:,:,:) = Diagnose_SigDot( xy_totDepth, xyz_U*xyz_CosLat, xyz_V*xyz_CosLat, xyz_Div)
       xyz_GeoPot(:,:,:) = Diagnose_GeoPot( xy_totDepth )

       call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                                        & ! (out)
            & theta=xyz_PTempBasic+xyz_PTempEdd, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)

       xyz_PressBaroc(:,:,:) = Diagnose_HydroPressEdd( xy_totDepth, xyz_DensEdd )
       
       !----------------------------------------------------------------------       
       !**********************************************************************
       !* Fast Physics 
       !**********************************************************************

       call calc_RHSFastPhysProc( &
           & xyz_U_RHSPhys, xyz_V_RHSPhys, xyz_PTemp_RHSPhys, xyz_Salt_RHSPhys,     & ! (out)
           & 'D'                                                                 )    ! (in)

       !**********************************************************************
       !* Slow Physics 
       !**********************************************************************

       call calc_RHSSlowPhysProc( &
           & xyz_U_RHSPhys, xyz_V_RHSPhys, xyz_PTemp_RHSPhys, xyz_Salt_RHSPhys,         & ! (out)
           & 'D'                                                                     )    ! (in)       
       

       !----------------------------------------------------------------------       
       !**********************************************************************
       !* Dynamics (For tracers)
       !**********************************************************************

       call calc_TracersRHSDynProc(                       &
            & wz_PTemp_RHSEx, wz_Salt_RHSEx,              &  ! (out)
            & xyz_PTemp_RHSPhys, xyz_Salt_RHSPhys         &  ! (in)  
            & )
       
       call integrate_ODE_wz( wz_PTempEdd,                & ! (out)
            & wz_PTemp_RHSEx, wz_PTempEddN, wz_PTempEddB, & ! (in)
            & wz_PTemp_RHSTmp,                            & ! (inout)
            & timeIntMode, Stage )                          ! (in)

       call integrate_ODE_wz( wz_Salt,                    & ! (out)
            & wz_Salt_RHSEx, wz_SaltN, wz_SaltB,          & ! (in)
            & wz_Salt_RHSTmp,                             & ! (inout)
            & timeIntMode, Stage )                          ! (in)


       
       !**********************************************************************
       !* Adjustment-typed Physics (For traces)
       !**********************************************************************
       
       xyz_PTempEdd(:,:,:) = xyz_wz(wz_PTempEdd)
       xyz_Salt(:,:,:) = xyz_wz(wz_Salt)
       
       !----------------------------------------------------------------------
       
       !**********************************************************************
       !* Dynamics (For momentums)
       !**********************************************************************

       call calc_MomBarocRHSDynProc(               &
           & wz_Vor_RHSEx, wz_Div_RHSEx,           &  ! (out)
           & xyz_U_RHSPhys, xyz_V_RHSPhys          &  ! (in)
           & )

       call integrate_ODE_wz( wz_Vor,                   & ! (out)
            & wz_Vor_RHSEx, wz_VorN, wz_VorB,           & ! (in)
            & wz_Vor_RHSTmp,                            & ! (inout)
            & timeIntMode, Stage )                        ! (in)

       call integrate_ODE_wz( wz_Div,                   & ! (out)
            & wz_Div_RHSEx, wz_DivN, wz_DivB,           & ! (in)
            & wz_Div_RHSTmp,                            & ! (inout)
            & timeIntMode, Stage )                        ! (in)
       
       !**********************************************************************
       !* Adjustment-typed Physics (For momentums)
       !**********************************************************************

       
       != Post process of each temporal loop =================================

       
       call VorDiv2UV(wz_Vor, wz_Div, &  ! (in)
            & xyz_U, xyz_V            &  ! (out)
            & )
       
    end do  ! End of do loop for a multi-stage temporal scheme.

    
    call apply_VBoundaryCondO( &
         & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt,     & ! (inout)
         & xyz_VViscCoefN, xyz_VDiffCoefN            & ! (in)
         & )
    
    ! ** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    xyz_UA(:,:,:) = xyz_U/xyz_CosLat;
    xyz_VA(:,:,:) = xyz_V/xyz_CosLat
    xyz_PTempEddA = xyz_PTempEdd
    xyz_SaltA = xyz_Salt
    xy_SurfHeightA(:,:) = xy_SurfHeight

  contains

      subroutine calc_RHSFastPhysProc( &
           & xyz_URHSPhys, xyz_VRHSPhys, xyz_PTempRHSPhys, xyz_SaltRHSPhys,   & ! (out)
           & EvalTLevel                                    )    ! (in)

        use HydroBouEqRHS_DynProc_mod, only: &
             & calc_HydroBouEqVViscRHS_xyz

        use TemporalIntegSet_mod, only: &
             & VDiffTermACoef
        
        real(DP), intent(out), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_URHSPhys, xyz_VRHSPhys, xyz_PTempRHSPhys, xyz_SaltRHSPhys
        character, intent(in) :: EvalTLevel

        real(DP) :: theta
        logical :: isRHSReplace
        
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_URHSTmp, xyz_VRHSTmp, xyz_PTempRHSTmp, xyz_SaltRHSTmp
        
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt
        real(DP), dimension(0:iMax-1,jMax) :: xy_SurfPressDummy
        

        theta = 1d0; isRHSReplace = .true.
        select case(EvalTLevel)
        case('B')
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHSTmp, xyz_VRHSTmp, xyz_PTempRHSTmp, xyz_SaltRHSTmp,     & ! (inout)  
                & xyz_UB, xyz_VB, xyz_PTempBasic/theta + xyz_PTempEddB, xyz_SaltB,                          & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,   & ! (in)
                & isRHSReplace=isRHSReplace )     
        case('N')
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHSTmp, xyz_VRHSTmp, xyz_PTempRHSTmp, xyz_SaltRHSTmp,     & ! (inout)  
                & xyz_UN, xyz_VN, xyz_PTempBasic/theta + xyz_PTempEddN, xyz_SaltN,                          & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,   & ! (in)
                & isRHSReplace=isRHSReplace )     
        case default
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHSTmp, xyz_VRHSTmp, xyz_PTempRHSTmp, xyz_SaltRHSTmp,     & ! (inout)  
                & xyz_U, xyz_V, xyz_PTempBasic/theta + xyz_PTempEdd, xyz_Salt,                              & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,   & ! (in)
                & isRHSReplace=isRHSReplace )     
        end select

        !
        !
        call Advance_VImplicitProc_DeltaForm_xyz( &
             & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt, xy_SurfPressDummy,              & ! (out)
             & xyz_URHSTmp, xyz_VRHSTmp, xyz_PTempRHSTmp, xyz_SaltRHSTmp,                                  & ! (inout)
             & xyz_VViscCoefN, xyz_VDiffCoefN, xy_totDepthBasic + xy_SurfHeight,                           & ! (in)
             & DynBC_Surface, DynBC_Bottom, ThermBC_Surface, ThermBC_Bottom, SaltBC_Surface, SaltBC_Bottom &  !(in)
             & )


        !
        !
        !$omp parallel workshare
        xyz_URHSPhys(:,:,:) = xyz_DU/DelTauBaroc
        xyz_VRHSPhys(:,:,:) = xyz_DV/DelTauBaroc
        xyz_PTempRHSPhys(:,:,:) = xyz_DPTempEdd/DelTauBaroc
        xyz_SaltRHSPhys(:,:,:) = xyz_DSalt/DelTauBaroc
        !$omp end parallel workshare

      end subroutine calc_RHSFastPhysProc

      subroutine calc_RHSSlowPhysProc( &
           & xyz_URHSPhys, xyz_VRHSPhys, xyz_PTempRHSPhys, xyz_SaltRHSPhys,   & ! (out)
           & EvalTLevel                                                    )    ! (in)

        use SGSEddyMixing_mod, only: SGSEddyMixing_AddMixingTerm
        
        real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_URHSPhys, xyz_VRHSPhys, xyz_PTempRHSPhys, xyz_SaltRHSPhys
        
        character, intent(in) :: EvalTLevel

        if (isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) then
           call SGSEddyMixing_AddMixingTerm(xyz_PTempRHSPhys, xyz_SaltRHSPhys,              & ! (inout)
                & xyz_PTempBasic + xyz_PTempEdd, xyz_Salt, xy_totDepthBasic + xy_SurfHeight & ! (in)
                & )
        end if
        
      end subroutine calc_RHSSlowPhysProc

      subroutine perform_AdjustTypePhysProc( &
           & xyz_PTempEdd, xyz_Salt )
        real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_PTempEdd, xyz_Salt

        
      end subroutine perform_AdjustTypePhysProc
      
      subroutine calc_TracersRHSDynProc(           &
           & wz_PTempRHS, wz_SaltRHS,              &  ! (out)
           & xyz_PTempRHSPhys, xyz_SaltRHSPhys     &  ! (in)  
           & )

        use HydroBouEqRHS_DynProc_mod, only: &
             & calc_TracersEqRHS_DynWithSrc, calc_HDiffRHS
        
        real(DP), intent(out), dimension(lMax,0:kMax) :: &
             & wz_PTempRHS, wz_SaltRHS
        real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_PTempRHSPhys, xyz_SaltRHSPhys

        integer, parameter :: nTracer = 2
        real(DP), dimension(lMax,0:kMax,nTracer) :: wza_TracersRHS
        real(DP), dimension(0:iMax-1,jMax,0:kMax,nTracer) :: &
             & xyza_Tracers, xyza_TracersRHSPhys
        integer :: n

        !
        !
        xyza_Tracers(:,:,:,1) = xyz_PTempEdd; xyza_TracersRHSPhys(:,:,:,1) = xyz_PTempRHSPhys
        xyza_Tracers(:,:,:,2) = xyz_Salt; xyza_TracersRHSPhys(:,:,:,2) = xyz_SaltRHSPhys

        call calc_TracersEqRHS_DynWithSrc( wza_TracersRHS,                      &
             & xyza_Tracers, xyza_TracersRHSPhys, nTracer,                      &
             & xyz_U*xyz_CosLat, xyz_V*xyz_CosLat, xyz_Div, xyz_SigDot          &
             & )

        wz_PTempRHS(:,:) = wza_TracersRHS(:,:,1)
        wz_SaltRHS(:,:) = wza_TracersRHS(:,:,2)

        !
        !
        call calc_HDiffRHS( wz_PTempRHS,                       & ! (inout)
             & wz_PTempEdd, hDiffCoef, hHyperDiffCoef, 'S',    & ! (in)
             & isRHSReplace=.false. )                            ! (in)

        call calc_HDiffRHS( wz_SaltRHS,                        & ! (inout)
             & wz_Salt, hDiffCoef, hHyperDiffCoef, 'S',        & ! (in)
             & isRHSReplace=.false. )                            ! (in)

      end subroutine calc_TracersRHSDynProc

      subroutine calc_MomBarocRHSDynProc(  &
           & wz_VorRHS, wz_DivRHS,         &  ! (out)
           & xyz_URHSPhys, xyz_VRHSPhys    &  ! (in)
           & )

        use HydroBouEqRHS_DynProc_mod, only: &
             & calc_VorDivEqRHS_DynWithSrc, calc_HDiffRHS
        
        real(DP), intent(out), dimension(lMax,0:kMax) :: &
             & wz_VorRHS, wz_DivRHS
        real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_URHSPhys, xyz_VRHSPhys

        real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: xyz_Press

        
        xyz_Press(:,:,:) = spread(xy_SurfPressN, 3, kMax+1) + xyz_PressBaroc
        call calc_VorDivEqRHS_DynWithSrc( wz_VorRHS, wz_DivRHS,                                       &
             & xyz_Vor, xyz_U, xyz_V, xy_SurfHeight, xyz_DensEdd, xyz_GeoPot, xyz_Press, xyz_SigDot,  &
             & xyz_U, xyz_V,                                                              &
             & xyz_URHSPhys, xyz_VRHSPhys )


        ! Horizontal mixing (with constant viscosity)
        !
        call calc_HDiffRHS(wz_VorRHS,                               &  ! (inout)
             & wz_Vor, hViscCoef, hHyperViscCoef, 'V',              &  ! (in)
             & isRHSReplace=.false. )                                  ! (in)

        call calc_HDiffRHS(wz_DivRHS,                               & ! (inout)
             & wz_Div, hViscCoef, hHyperViscCoef, 'V', 2d0,         & ! (in)
             & isRHSReplace=.false. )                                 ! (in)

        
      end subroutine calc_MomBarocRHSDynProc
      
      subroutine integrate_ODE_wz(wz_Ret, wz_RHS, wz_N, wz_B, wz_RHSTmp, odeIntMode, stage)
        
        real(DP), intent(out), dimension(lMax,0:kMax) :: wz_Ret
        real(DP), intent(in), dimension(lMax,0:kMax) :: wz_RHS, wz_N, wz_B
        real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_RHSTmp
        integer, intent(in) :: odeIntMode
        integer, intent(in) :: stage

        select case(odeIntMode)
        case(timeIntMode_Euler)
           wz_Ret(:,:) = timeIntEuler(wz_N, wz_RHS)
        case(timeIntMode_RK2)
           wz_Ret(:,:) = timeIntRK(wz_N, wz_RHS, 2, stage, wz_RHSTmp)
        case(timeIntMode_RK4)
           wz_Ret(:,:) = timeIntRK(wz_N, wz_RHS, 4, stage, wz_RHSTmp)
        case(timeIntMode_LFAM3)
           wz_Ret(:,:) = timeIntLFAM3(wz_N, wz_B, wz_RHS, stage)
        end select
        
      end subroutine integrate_ODE_wz

      subroutine calc_VViscDiffCoef( &
           & xyz_VViscCoef, xyz_VDiffCoef,      & ! (out)
           & xyz_U, xyz_V, xyz_PTemp, xyz_Salt  & ! (in)
           & )

         use EOSDriver_mod, only: EOSDriver_Eval

         ! 宣言文; Declaration statement
         !    
         real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
              & xyz_VViscCoef, xyz_VDiffCoef
         real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) ::  &
              & xyz_U, xyz_V, xyz_PTemp, xyz_Salt

         ! 局所変数
         ! Local variables
         !         
         real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
              & xyz_Depth, xyz_DensPot, xyz_RefPress, &
              & xyz_VViscCoefMixLyr, xyz_VDiffCoefMixLyr
         
         integer :: k

         ! 実行文; Executable statement
         !
         
         do k=0, kMax
            xyz_Depth(:,:,k) = (xy_totDepthBasic + xy_SurfHeight)*g_Sig(k)
         end do
         
         call calc_vViscDiffCoef_MixLyrSimple( xyz_VViscCoefMixLyr, xyz_VDiffCoefMixLyr, &
              & xyz_Depth )

!!$         xyz_VViscCoef(:,:,:) = vViscCoef; xyz_VDiffCoef(:,:,:) = vDiffCoef
         xyz_VViscCoef(:,:,:) = xyz_VViscCoefMixLyr + 1d-3!vViscCoef
         xyz_VDiffCoef(:,:,:) = xyz_VDiffCoefMixLyr + vDiffCoef
         return

         xyz_RefPress = 0d0
         call EOSDriver_Eval( rhoEdd=xyz_DensPot,               & ! (out)
              & theta=xyz_PTemp, S=xyz_Salt, p=xyz_RefPress )     ! (in)

         call calc_vViscDiffCoef_PP81( &
              & xyz_VViscCoef, xyz_VDiffCoef,                              & !(out)
              & xyz_U, xyz_V, xyz_DensPot, xy_totDepthBasic+xy_SurfHeight, & !(in)
              & vViscCoef, vDiffCoef                                       & !(in)
              & )

       end subroutine calc_VViscDiffCoef

       subroutine UV2VorDiv( xyz_U, xyz_V,          & ! (in)
            & wz_Vor, wz_Div, xyz_Vor, xyz_Div      & ! (out)
            & )
         real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax)  :: xyz_U, xyz_V
         real(DP), intent(out), dimension(lMax,0:kMax)  :: wz_Vor, wz_Div
         real(DP), optional, intent(out), dimension(0:iMax-1,jMax,0:kMax)  :: xyz_Vor, xyz_Div

         call wz_VectorCosLat2VorDiv(xyz_U*xyz_CosLat, xyz_V*xyz_CosLat, &
              & wz_Vor, wz_Div)
         
         if(present(xyz_Vor) .and. present(xyz_Div)) then
            xyz_Vor(:,:,:) = xyz_wz(wz_Vor)
            xyz_Div(:,:,:) = xyz_wz(wz_Div)
         end if
         
       end subroutine UV2VorDiv

       subroutine VorDiv2UV( wz_Vor, wz_Div, &  ! (in)
            & xyz_U, xyz_V, MultiCosFlag     &  ! (out)
            & )
         real(DP), optional, intent(in), dimension(lMax,0:kMax)  :: wz_Vor, wz_Div
         real(DP), intent(out), dimension(0:iMax-1,jMax,0:kMax)  :: xyz_U, xyz_V
         logical, intent(in), optional :: MultiCosFlag

         call wz_VorDiv2VectorCosLat(wz_Vor, wz_Div, &
              & xyz_U, xyz_V)

         if(present(MultiCosFlag) .and. MultiCosFlag) then
            return
         else
            xyz_U(:,:,:) = xyz_U/xyz_CosLat
            xyz_V(:,:,:) = xyz_V/xyz_CosLat
         end if
       end subroutine VorDiv2UV
       
     end subroutine HydroBouEqSolver_AdvanceTStep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_vViscDiffCoef_MixLyrSimple( &
       & xyz_VViscCoef, xyz_VDiffCoef, &
       & xyz_Depth &
       & )

    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_VViscCoef, xyz_VDiffCoef
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_Depth

    real(DP), parameter :: MixLyrDepth = 50d0
    real(DP), parameter :: LInv = 1d0/(0.1d0*MixLyrDepth)
    real(DP), parameter :: ViscfCoefMax = 2d-3
    real(DP), parameter :: DiffCoefMax = 2d-3

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Func

    !$omp parallel workshare
    xyz_Func(:,:,:) = 0.5d0 - atan((abs(xyz_Depth) - MixLyrDepth)*LInv)/PI
    xyz_VViscCoef(:,:,:) = ViscfCoefMax*xyz_Func
    xyz_VDiffCoef(:,:,:) = DiffCoefMax*xyz_Func
    !$omp end parallel workshare

  end subroutine calc_vViscDiffCoef_MixLyrSimple
       
    
  subroutine calc_vViscDiffCoef_PP81( &
       & xyz_VViscCoef, xyz_VDiffCoef, &
       & xyz_U, xyz_V, xyz_DensPot, xy_totDepth, vViscCoefBG, vDiffCoefBG )

    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_VViscCoef, xyz_VDiffCoef
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_U, xyz_V, xyz_DensPot
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: vViscCoefBG, vDiffCoefBG

    real(DP), parameter :: AvRic = 1d-2
    real(DP), parameter :: a = 5d0
    integer, parameter :: n = 2
    

    real(DP), dimension(0:iMax-1,jMax, 0:kMax) :: xyz_Ri
    integer :: i, j, k

    xyz_Ri(:,:,:) = diagnose_RicardsonNumber(xyz_U, xyz_V, xyz_DensPot, xy_totDepth)
    where(xyz_Ri < 0d0)
       xyz_Ri = 0d0
    end where
    
    !$omp parallel workshare
    xyz_VViscCoef(:,:,:) = AvRic/(1d0 + a*xyz_Ri)**n + 1d-3!vViscCoefBG
    xyz_VDiffCoef(:,:,:) = xyz_VViscCoef(:,:,:)/(1d0 + a*xyz_Ri) + vDiffCoefBG
    !$omp end parallel workshare

!!$    write(*,*) "=-------------"
!!$    write(*,*) "Av:", xyz_VViscCoef(0,1:32,0)
!!$    write(*,*) "MaxAv;", maxval(xyz_VViscCoef), maxloc(xyz_VViscCoef) 

  end subroutine calc_vViscDiffCoef_PP81



  function diagnose_RicardsonNumber(xyz_U, xyz_V, xyz_DensPot, xy_totDepth) result(xyz_Ri)

    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_U, xyz_V, xyz_DensPot
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Ri, xyz_N2

    
!!$    xyz_Ri(:,:,:) = ( &
!!$         &      -Grav/RefDens*xyz_Dz_xyz(xyz_DensPot)           &
!!$         &     /(xyz_Dz_xyz(xyz_U)**2 + xyz_Dz_xyz(xyz_V)**2) &
!!$         &   )

    xyz_Ri(:,:,:) = ( &
         &      -Grav/RefDens*xyz_Dz_xyz(xyz_DensPot)           &
         &     /(xyz_Dz_xyz(xyz_U)**2 + xyz_Dz_xyz(xyz_V)**2 + 1d-14) &
         &   )
    
!!$    xyz_N2 = -Grav/RefDens*xyz_Dz_xyz(xyz_DensPot)
!!$    write(*,*) "=-------------"
!!$    write(*,*) "Ri:", xyz_Ri(0,20,:)
!!$    write(*,*) "N2:", xyz_N2(0,20,:)
!!$    write(*,*) "(dUdz)^2:", xyz_N2(0,20,:)/xyz_Ri(0,20,:)
!!$    write(*,*) "DensPot:", xyz_Denspot(0,20,:)
!!$    write(*,*) 
    
  end function diagnose_RicardsonNumber

   function xyz_Dz_xyz(xyz) 

    use VariableSet_mod

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Dz_xyz(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: s(0:kMax), t(0:kMax), xyt(0:iMax-1,jMax,0:tMax)
    integer :: k

    ! 実行文; Executable statement
    !
    
    t(1:kMax-1) = g_Sig(0:kMax-2) - g_Sig(1:kMax-1)
    s(1:kMax-1) = g_Sig(1:kMax-1) - g_Sig(2:kMax)

    !$omp parallel do
    do k=1,kMax-1
       xyz_Dz_xyz(:,:,k) = &
            & (s(k)**2*xyz(:,:,k-1) - (s(k)**2-t(k)**2)*xyz(:,:,k) - t(k)**2*xyz(:,:,k+1)) &
            & /(s(k)*t(k)*(s(k) + t(k)))/xy_totDepthBasic
    end do
    xyz_Dz_xyz(:,:,0) = &
         & (xyz(:,:,0) - xyz(:,:,1))/(g_Sig(0) - g_Sig(1))/xy_totDepthBasic
    xyz_Dz_xyz(:,:,kMax) = &
         & (xyz(:,:,kMax-1) - xyz(:,:,kMax))/(g_Sig(kMax-1) - g_Sig(kMax))/xy_totDepthBasic

  end function xyz_Dz_xyz

end module HydroBoudEqVorDivForm_TimeInteg_mod

