!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBoudEq_TimeInteg_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING 

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use Constants_mod, only: &
       & Omega, Grav, RPlanet, &
       & hViscCoef, vViscCoef, &
       & hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef, &
       & RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use GovernEqSet_mod, only: &
       & SGSEddyMixType, &
       & isPhysicsCompActived

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType, &
       & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress, &
       & ThermBCTYPE_FluxFixed, ThermBCTYPE_Adiabat, ThermBCTYPE_TempFixed, ThermBCTYPE_TempRelaxed, & 
       & SaltBCTYPE_FluxFixed, SaltBCTYPE_Adiabat, SaltBCTYPE_SaltFixed, SaltBCTYPE_SaltRelaxed, & 
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom 

  use TemporalIntegUtil_mod, only: &
       & TemporalIntegUtil_Init, TemporalIntegUtil_Final, &
       & timeIntMode_Euler, xy_timeIntEuler, wt_timeIntEuler, &
       & timeIntMode_LFTR, xy_timeIntLFTR, wt_timeIntLFTR, &
       & timeIntMode_LFAM3, xy_timeIntLFAM3, wt_timeIntLFAM3, &
       & timeIntMode_RK2, xy_timeIntRK2, wt_timeIntRK2, &
       & timeIntMode_RK4, xy_timeIntRK4, wt_timeIntRK4, &
       & TemporalIntegUtil_GetDDtCoef, &
       & TemporalIntegUtil_SetDelTime

  use TemporalIntegSet_mod, only: &
       & CurrentTimeStep, SubCycleNum, &
       & SemiImplicitFlag, &
       & nShortTimeLevel

  use HydroBouEqSolverRHS_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final, &
       & calc_HydroBouEqInvisRHS, calc_HydroBouEqHViscRHS, calc_HydroBouEqVViscRHS, &
       & correct_DivEqRHSUnderRigidLid

  use HydroBouEqSolverVImplProc_mod, only: &
       & HydroBouEqSolverVImplProc_Init, HydroBouEqSolverVImplProc_Final, &
       & Advance_VImplicitProc

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

    call HydroBouEqSolverRHS_Init()
    call TemporalIntegUtil_Init( iMax, jMax, kMax, lMax, tMax, 0d0 )
    call HydroBouEqSolverVImplProc_Final()
    
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

    call HydroBouEqSolverRHS_Final()
    call TemporalIntegUtil_Final()
    call HydroBouEqSolverVImplProc_Final()

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
         & xy_WindStressU, xy_WindStressV, xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTime
    integer, intent(in) :: timeIntMode
    integer, intent(in) :: nStage_BarocTimeInt
    logical, intent(in) :: isVarBUsed_BarocTimeInt

    ! 局所変数
    ! Local variables
    !

    real(DP), dimension(lMax, 0:tMax) :: wt_Vor, wt_VorN, wt_VorB
    real(DP), dimension(lMax, 0:tMax) :: wt_Div, wt_DivN, wt_DivB
    real(DP), dimension(lMax, 0:tMax) :: wt_PTempEdd, wt_PTempEddN, wt_PTempEddB
    real(DP), dimension(lMax, 0:tMax) :: wt_Salt, wt_SaltN, wt_SaltB
    real(DP) :: xy_SurfHeight(0:iMax-1, jMax)

    real(DP) :: xy_SurfPress(0:iMax-1,jMax)


    real(DP), dimension(0:iMax-1, jMax, 0:kMax) :: xyz_Urf, xyz_Vrf
    real(DP), dimension(lMax, 0:kMax) ::  wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS
    real(DP) :: w_SurfHeightExplRHS(lMax)
    real(DP), dimension(lMax, 0:kMax) ::  wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp

    ! Work variables for some temporal schemes, such as Runge=Kutta scheme etc. 
    real(DP), dimension(lMax,0:tMax) :: wt_VorRHSTIntTmp, wt_DivRHSTIntTmp, wt_PTempRHSTIntTmp, &
         & wt_SaltRHSTIntTmp
    real(DP), dimension(0:iMax-1,jMax) :: xy_SurfHeightRHSTIntTmp

    integer :: Stage
    character(TOKEN) :: TIntType_SurfPressTerm
    
    real(DP) :: xyz_CosLat(0:iMax-1, jMax, 0:kMax)
    real(DP), dimension(lMax,2) :: wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS

    ! 実行文; Executable statement
    !
    
    ! * Preparation 


    xyz_CosLat = cos(xyz_Lat)

    ! Set some variables at the time level N. 
    xyz_Urf = xyz_UN*xyz_CosLat; xyz_Vrf = xyz_VN*xyz_CosLat
    call wt_VectorCosLat2VorDiv( xyz_Urf, xyz_Vrf, & ! (in)
         & wt_VorN, wt_DivN )                        ! (out)
    wt_PTempEddN = wt_xyz(xyz_PTempEddN)
    wt_SaltN = wt_xyz(xyz_SaltN)

    ! Set some variables at the time level B if they are used. 
    if( isVarBUsed_BarocTimeInt ) then
       call wt_VectorCosLat2VorDiv( xyz_UB*xyz_CosLat, xyz_VB*xyz_CosLat, & ! (in)
            & wt_VorB, wt_DivB )                        ! (out)
       wt_PTempEddB = wt_xyz(xyz_PTempEddB)
       wt_SaltB = wt_xyz(xyz_SaltB)
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call TemporalIntegUtil_SetDelTime(DelTime)

    wt_Div = wt_DivN;  wt_Vor = wt_VorN;  wt_PTempEdd = wt_PTempEddN; wt_Salt = wt_SaltN
    xy_SurfHeight = xy_SurfHeightN
    xy_SurfPress = xy_SurfPressN

    !
    do Stage=1, nStage_BarocTimeInt

       select case(timeIntMode) !======================================================================      
       case(timeIntMode_Euler)  !***** Using Euler scheme *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltExplRHS, 'D', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 
          call correct_DivEqRHS_RigidLid('CRANKNIC', wt_DivN)
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call update_VBCRHS( wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, wt_Vor, wt_Div, wt_PTempEdd, wt_Salt )

          call timeInt_Euler()

       case(timeIntMode_RK2)    !***** Using RK2 scheme   *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltExplRHS, 'D', 1d0, .false.); 
          call add_ImplRHS_into_RHS(); 
          call correct_DivEqRHS_RigidLid('CRANKNIC', wt_DivN)
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call update_VBCRHS( wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, wt_Vor, wt_Div, wt_PTempEdd, wt_Salt )

          call timeInt_RK2(Stage)

       case(timeIntMode_RK4)    !***** Using RK4 scheme   *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltExplRHS, 'D', 1d0, .false.); 
          call add_ImplRHS_into_RHS(); 
          call correct_DivEqRHS_RigidLid('CRANKNIC', wt_DivN)
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call update_VBCRHS( wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, wt_Vor, wt_Div, wt_PTempEdd, wt_Salt )

          call timeInt_RK4(Stage)

       case(timeIntMode_LFTR)   !***** Using LFTR scheme   *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltExplRHS, 'D', 1d0, .false.); 
          call add_ImplRHS_into_RHS(); 
          call correct_DivEqRHS_RigidLid('CRANKNIC', wt_DivN)
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call update_VBCRHS( wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, wt_Vor, wt_Div, wt_PTempEdd, wt_Salt )

          call timeInt_LFTR(Stage)

       case(timeIntMode_LFAM3)  !***** Using LFAM3 scheme  *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS,  w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )

          select case(Stage)
             case(1)
                call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'B', 0.5d0, .false.); 
                call add_ImplRHS_into_RHS(); 
                call correct_DivEqRHS_RigidLid('CRANKNIC_WithLF', wt_DivB)
             case(2)
                call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'N',0.5d0,.false.)!'N', 0.5d0, .false.); 
                call add_ImplRHS_into_RHS(); 
                call correct_DivEqRHS_RigidLid('CRANKNIC', wt_DivN)
             end select

             call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )

             call timeInt_LFAM3(Stage)

       end select              !=====================================================================================
       
       call apply_boundaryConditions2(wt_Vor, wt_Div, wt_PTempEdd, wt_Salt, &
            & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS )

       call wt_VorDiv2VectorCosLat( wt_Vor, wt_Div,      &  !(in)
            & xyz_Urf, xyz_Vrf                           &  !(out)
            & )

    end do  ! End of do loop for a multi-stage temporal scheme.


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_UA = xyz_Urf / xyz_CosLat;  xyz_VA = xyz_Vrf / xyz_CosLat 
    xyz_PTempEddA = xyz_wt(wt_PTempEdd)
    xyz_SaltA = xyz_wt(wt_Salt)
!    xy_SurfPressA = xy_SurfPress
    xy_SurfHeightA = xy_SurfHeight


    contains


      subroutine calc_InvisRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, & 
           & tLevel)
        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        real(DP), dimension(lMax), intent(inout) :: w_SurfHeightRHS
        character, intent(in), optional :: tLevel

        character :: tLvl

        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        select case(tLvl)
           case('B')
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wt(wt_Vor), xyz_wt(wt_Div), xyz_wt(wt_PTempEdd), xyz_wt(wt_Salt),    &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic                                &  ! (in)
                   & )
           case('N')
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wt(wt_Vor), xyz_wt(wt_Div), xyz_wt(wt_PTempEdd), xyz_wt(wt_Salt),    &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic                                &  ! (in)
                   & )
           case default
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wt(wt_Vor), xyz_wt(wt_Div), xyz_wt(wt_PTempEdd), xyz_wt(wt_Salt),    &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic                                &  ! (in)
                   & )
        end select
      end subroutine calc_InvisRHS

      subroutine calc_HViscRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, & 
           & tLevel, isRHSAppend)

        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        logical, intent(in), optional :: isRHSAppend

        character :: tLvl
        logical :: isRHSReplace

        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        isRHSReplace = .true.
        if (present(isRHSAppend)) isRHSReplace = (.not. isRHSAppend)

        select case(tLvl)
        case('B')
           call calc_HydroBouEqHViscRHS( &
                & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,                          & ! (inout)  
                & wz_wt(wt_VorB), wz_wt(wt_DivB), wz_wt(wt_PTempEddB), wz_wt(wt_SaltB),   & ! (in)
                & hViscCoef, hHyperViscCoef, hDiffCoef,                                   & ! (in)
                & isRHSReplace=isRHSReplace )

        case default
           call calc_HydroBouEqHViscRHS( &
                & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,                        & ! (inout)  
                & wz_wt(wt_Vor), wz_wt(wt_Div), wz_wt(wt_PTempEdd), wz_wt(wt_Salt),     & ! (in)
                & hViscCoef, hHyperViscCoef, hDiffCoef,                                 & ! (in)
                & isRHSReplace=isRHSReplace )
        end select
      end subroutine calc_HViscRHS

      ! Evaluation of the vertical viscid term. 
      ! \[
      !   dq/dt = \theta F(q^m), 
      ! \]
      ! where $F(q^m)$ is vertical viscid term, q is an arbitary physical quantity, 
      ! and \theta is a coefficent of $F(q^m)$.
      ! 
      subroutine calc_VViscRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, & 
           & tLevel, vViscTermCoef, isRHSAppend)

        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        real(DP), intent(in), optional :: vViscTermCoef
        logical, intent(in), optional :: isRHSAppend

        character :: tLvl
        real(DP) :: theta
        logical :: isRHSReplace
        real(DP) :: xyz_PTempBasic(0:iMax-1,jMax,0:kMax)
        integer :: k

        !
        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        theta = 1d0
        if (present(vViscTermCoef)) theta = vViscTermCoef

        isRHSReplace = .true.
        if (present(isRHSAppend)) isRHSReplace = (.not. isRHSAppend)

        !
        forAll(k=0:kMax) xyz_PTempBasic(:,:,k) = z_PTempBasic(k)

        !
        select case(tLvl)
        case('B')
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,                           & ! (inout)  
                & wz_wt(wt_VorB), wz_wt(wt_DivB), wz_xyz(xyz_PTempBasic/theta + xyz_PTempEddB), wz_wt(wt_SaltB),  & ! (in)
                & theta*vViscCoef, vHyperViscCoef, theta*vDiffCoef,                                               & ! (in)
                & isRHSReplace=isRHSReplace )     
        case('N')
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,                            & ! (inout)  
                & wz_wt(wt_VorN), wz_wt(wt_DivN), wz_xyz(xyz_PTempBasic/theta + xyz_PTempEddN), wz_wt(wt_SaltN), & ! (in)
                & theta*vViscCoef, vHyperViscCoef, theta*vDiffCoef,                                                & ! (in)
                & isRHSReplace=isRHSReplace )     
        case default 
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,                            & ! (inout)  
                & wz_wt(wt_Vor), wz_wt(wt_Div), wz_xyz(xyz_PTempBasic/theta) + wz_wt(wt_PTempEdd), wz_wt(wt_Salt), & ! (in)
                & theta*vViscCoef, vHyperViscCoef, theta*vDiffCoef,                                                & ! (in)
                & isRHSReplace=isRHSReplace )
        end select

      end subroutine calc_VViscRHS

      subroutine calc_ExplTermWithPhysicsRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, & 
           & tLevel, isRHSAppend)

        use SGSEddyMixing_mod, only: SGSEddyMixing_AddMixingTerm

        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        logical, intent(in), optional :: isRHSAppend

        character :: tLvl
        logical :: isRHSReplace
        integer :: k
        real(DP) :: xyz_PTempBasic(0:iMax-1,jMax,0:kMax)

        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        forAll(k=0:kMax) xyz_PTempBasic(:,:,k) = z_PTempBasic(k)

        if(isPhysicsCompActived(SGSEddyMixType)) then
           call SGSEddyMixing_AddMixingTerm(wz_PTempRHS, wz_SaltRHS, &
                & wz_wt(wt_PTempEdd)+wz_xyz(xyz_PTempBasic), wz_wt(wt_Salt), xy_totDepthBasic+xy_SurfHeight)
        end if

      end subroutine calc_ExplTermWithPhysicsRHS

      subroutine add_ImplRHS_into_RHS()

        wz_VorExplRHS = wz_VorExplRHS + wz_VorImplRHSTmp
        wz_DivExplRHS = wz_DivExplRHS + wz_DivImplRHSTmp
        wz_PTempExplRHS = wz_PTempExplRHS + wz_PTempImplRHSTmp
        wz_SaltExplRHS  = wz_SaltExplRHS + wz_SaltImplRHSTmp

      end subroutine add_ImplRHS_into_RHS

      subroutine correct_DivEqRHS_RigidLid( &
           & TIntType_SurfPressTerm, wt_Div)

        character(*), intent(in) :: TIntType_SurfPressTerm
        real(DP), intent(in) :: wt_Div(lMax,0:tMax)

        call correct_DivEqRHSUnderRigidLid( wz_DivExplRHS, xy_SurfPressA, &
             & wz_wt(wt_Div), xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
             & 1d0 / (TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime), &
             & TIntType_SurfPressTerm )
        xy_SurfPress = xy_SurfPressN
      end subroutine correct_DivEqRHS_RigidLid


      !* Interfaces for calling each implicit scheme. 
      !

      function wt_implicitSolver_Vor(wt_val) result(wt_ret)
         real(DP), intent(in) :: wt_val(lMax,0:tMax)
         real(DP) :: wt_ret(lMax,0:tMax)

         wt_ret = wt_val
         
         call Update_VBCRHS(wa_VorBCRHS=wa_VorBCRHS, wt_Vor=wt_ret)
         call Advance_VImplicitProc(wt_Vor=wt_ret, &
              & wa_VorBCRHS=wa_VorBCRHS, xy_totDepth=xy_totDepthBasic+xy_SurfHeight, &
              & vViscDiffTermCoef=0.5d0, dt=TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime,  &
              & DynBCSurf=DynBC_Surface, DynBCBottom=DynBC_Bottom )
       end function wt_implicitSolver_Vor

       function wt_implicitSolver_Div(wt_val) result(wt_ret)
         real(DP), intent(in) :: wt_val(lMax,0:tMax)
         real(DP) :: wt_ret(lMax,0:tMax)

         wt_ret = wt_val

         call Update_VBCRHS(wa_DivBCRHS=wa_DivBCRHS, wt_Div=wt_ret)
         call Advance_VImplicitProc(wt_Div=wt_ret, &
              & wa_DivBCRHS=wa_DivBCRHS, xy_totDepth=xy_totDepthBasic+xy_SurfHeight, &
              & vViscDiffTermCoef=0.5d0, dt=TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime,  &
              & DynBCSurf=DynBC_Surface, DynBCBottom=DynBC_Bottom )
       end function wt_implicitSolver_Div

       function wt_implicitSolver_PTemp(wt_val) result(wt_ret)
         real(DP), intent(in) :: wt_val(lMax,0:tMax)
         real(DP) :: wt_ret(lMax,0:tMax)

         wt_ret = wt_val
         call Update_VBCRHS(wa_PTempEddBCRHS=wa_PTempEddBCRHS, wt_PTempEdd=wt_ret)
         call Advance_VImplicitProc(wt_PTempEdd=wt_ret, &
              & wa_PTempEddBCRHS=wa_PTempEddBCRHS, xy_totDepth=xy_totDepthBasic+xy_SurfHeight, &
              & vViscDiffTermCoef=0.5d0, dt=TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime,  &
              & ThermBCSurf=ThermBC_Surface, ThermBCBottom=ThermBC_Bottom )

       end function wt_implicitSolver_PTemp

       function wt_implicitSolver_Salt(wt_val) result(wt_ret)
         real(DP), intent(in) :: wt_val(lMax,0:tMax)
         real(DP) :: wt_ret(lMax,0:tMax)

         wt_ret = wt_val
         call Update_VBCRHS(wa_SaltBCRHS=wa_SaltBCRHS, wt_Salt=wt_ret)
         call Advance_VImplicitProc(wt_Salt=wt_ret, &
              & wa_SaltBCRHS=wa_SaltBCRHS, xy_totDepth=xy_totDepthBasic+xy_SurfHeight, &
              & vViscDiffTermCoef=0.5d0, dt=TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime,  &
              & SaltBCSurf=SaltBC_Surface, SaltBCBottom=SaltBC_Bottom )

       end function wt_implicitSolver_Salt

       function xy_implicitSolver_SurfHeight(xy_val) result(xy_ret)
         real(DP), intent(in) :: xy_val(0:iMax-1,jMax)
         real(DP) :: xy_ret(0:iMax-1,jMax)

       end function xy_implicitSolver_SurfHeight

       ! * Interfaces for calling a subroutine to perform a temporal 
       !   integration of ODE for prognostic variables. 
       !

       subroutine timeInt_Euler()
         wt_Vor = wt_timeIntEuler( wt_VorN, wt_wz(wz_VorExplRHS) )
         wt_Div = wt_timeIntEuler( wt_DivN, wt_wz(wz_DivExplRHS) )
         wt_PTempEdd = wt_timeIntEuler( wt_PTempEddN, wt_wz(wz_PTempExplRHS) )
         wt_Salt = wt_timeIntEuler( wt_SaltN, wt_wz(wz_SaltExplRHS) )
         xy_SurfHeight = xy_timeIntEuler( xy_SurfHeightN, xy_w(w_SurfHeightExplRHS) )
       end subroutine timeInt_Euler
       
       subroutine timeInt_RK2(RKStage)
         integer, intent(in) :: RKStage
         wt_Vor = wt_timeIntRK2( wt_VorN, wt_wz(wz_VorExplRHS), RKStage, wt_VorRHSTIntTmp )
         wt_Div = wt_timeIntRK2( wt_DivN, wt_wz(wz_DivExplRHS), RKStage, wt_DivRHSTIntTmp )
         wt_PTempEdd = wt_timeIntRK2( wt_PTempEddN, wt_wz(wz_PTempExplRHS), RKStage, wt_PTempRHSTIntTmp )
         wt_Salt = wt_timeIntRK2( wt_SaltN, wt_wz(wz_SaltExplRHS), RKStage, wt_SaltRHSTIntTmp )
         xy_SurfHeight = xy_timeIntRK2( xy_SurfHeightN, xy_w(w_SurfHeightExplRHS), RKStage, xy_SurfHeightRHSTIntTmp )
       end subroutine timeInt_RK2

       subroutine timeInt_RK4(RKStage)
         integer, intent(in) :: RKStage
         wt_Vor = wt_timeIntRK4( wt_VorN, wt_wz(wz_VorExplRHS), RKStage, wt_VorRHSTIntTmp )
         wt_Div = wt_timeIntRK4( wt_DivN, wt_wz(wz_DivExplRHS), RKStage, wt_DivRHSTIntTmp )
         wt_PTempEdd = wt_timeIntRK4( wt_PTempEddN, wt_wz(wz_PTempExplRHS), RKStage, wt_PTempRHSTIntTmp )
         wt_Salt = wt_timeIntRK4( wt_SaltN, wt_wz(wz_SaltExplRHS), RKStage, wt_SaltRHSTIntTmp )
         xy_SurfHeight = xy_timeIntRK4( xy_SurfHeightN, xy_w(w_SurfHeightExplRHS), RKStage, xy_SurfHeightRHSTIntTmp )
       end subroutine timeInt_RK4

       subroutine timeInt_LFTR(Stage)
         integer, intent(in) :: Stage
         wt_Vor = wt_timeIntLFTR( wt_VorN, wt_VorB, wt_wz(wz_VorExplRHS), Stage )
         wt_Div = wt_timeIntLFTR( wt_DivN, wt_DivB, wt_wz(wz_DivExplRHS), Stage )
         wt_PTempEdd = wt_timeIntLFTR( wt_PTempEddN, wt_PTempEddB, wt_wz(wz_PTempExplRHS), Stage )
         wt_Salt = wt_timeIntLFAM3( wt_SaltN, wt_SaltB, wt_wz(wz_SaltExplRHS), Stage )
         xy_SurfHeight = xy_timeIntLFTR( xy_SurfHeightN, xy_SurfHeightB, xy_w(w_SurfHeightExplRHS), Stage )
       end subroutine timeInt_LFTR

       subroutine timeInt_LFAM3(Stage)
         integer, intent(in) :: Stage
         wt_Vor = wt_timeIntLFAM3( wt_VorN, wt_VorB, wt_wz(wz_VorExplRHS), Stage, wt_implicitSolver_Vor )
         wt_Div = wt_timeIntLFAM3( wt_DivN, wt_DivB, wt_wz(wz_DivExplRHS), Stage, wt_implicitSolver_Div )
         wt_PTempEdd = wt_timeIntLFAM3( wt_PTempEddN, wt_PTempEddB, wt_wz(wz_PTempExplRHS), Stage, wt_implicitSolver_PTemp )
         wt_Salt = wt_timeIntLFAM3( wt_SaltN, wt_SaltB, wt_wz(wz_SaltExplRHS), Stage, wt_implicitSolver_Salt )
         xy_SurfHeight = xy_timeIntLFAM3( xy_SurfHeightN, xy_SurfHeightB, xy_w(w_SurfHeightExplRHS), Stage )
       end subroutine timeInt_LFAM3


       subroutine replace_RHS_with_VBCTIntRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS )

         ! モジュール引用; Use statements
         !
         use VariableSet_mod, only: &
              & xy_SeaSurfTemp, xy_SeaSurfSalt

         use BoundCondSet_mod, only: &
              & SurfTempRelaxedTime, SurfSaltRelaxedTime

         ! 宣言文; Declaration statement
         !
         real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
         

         ! 作業変数
         ! Work variables
         !
         real(DP), dimension(lMax,0:kMax) :: wz_PTemp, wz_Salt

         ! 実行文; Executable statement
         !

         wz_PTemp = wz_wt(wt_PTempEdd)  + wz_xyz( spread(spread(z_PTempBasic,1,jMax),1,iMax) )

         if(ThermBC_Surface == ThermBCTYPE_TempRelaxed) then
            wz_PTempRHS(:,0) = - ( wz_PTemp(:,0) - w_xy(xy_SeaSurfTemp) )/SurfTempRelaxedTime
         end if

         if(ThermBC_Bottom == ThermBCTYPE_TempRelaxed) then
            wz_PTempRHS(:,kMax) = 0d0
         end if

         !
         wz_Salt = wz_wt(wt_Salt)
         if(SaltBC_Surface == SaltBCTYPE_SaltRelaxed) then
            wz_SaltRHS(:,0) = - ( wz_Salt(:,0) - w_xy(xy_SeaSurfSalt) )/SurfSaltRelaxedTime
         end if

         if(SaltBC_Bottom == SaltBCTYPE_SaltRelaxed) then
            wz_SaltRHS(:,kMax) = 0d0
         end if

       end subroutine replace_RHS_with_VBCTIntRHS

  end subroutine HydroBouEqSolver_AdvanceTStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine update_VBCRHS( &
       & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, &
       & wt_Vor, wt_Div, wt_PTempEdd, wt_Salt )

    
    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & xy_WindStressU, xy_WindStressV, &
         & xy_SeaSurfTemp, xy_SurfTempFlux, &
         & xy_SeaSurfSalt, xy_SurfSaltFlux, &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax, 1:2), intent(inout), optional :: &
         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
    real(DP), dimension(lMax,0:tMax), intent(in), optional :: wt_Vor, wt_Div, wt_PTempEdd, wt_Salt

    ! 作業変数
    ! Work variables
    !
    real(DP) :: xy_Coef(0:iMax-1,jMax)
    real(DP) :: xyz_PTempTmp(0:iMax-1,jMax,0:kMax)
    real(DP) :: wz_PTempTmp(lMax,0:kMax), wz_SaltTmp(lMax,0:kMax)
    integer :: k
    logical :: isVorUpdate, isDivUpdate, isPTempEddUpdate, isSaltUpdate

    ! 実行文; Executable statement
    !

    xy_Coef = xy_totDepthBasic*cos(xyz_Lat(:,:,1))

    isVorUpdate = .false.; isDivUpdate = .false.; isPTempEddUpdate = .false.; isSaltUpdate = .false.
    if(present(wa_VorBCRHS)) isVorUpdate = .true.
    if(present(wa_DivBCRHS)) isDivUpdate = .true.
    if(present(wa_PTempEddBCRHS)) isPTempEddUpdate = .true.
    if(present(wa_SaltBCRHS)) isSaltUpdate = .true.

    if(isVorUpdate .or. isDivUpdate ) then
       select case(DynBC_Surface)
       case(DynBCTYPE_NoSlip)
          if(isVorUpdate) wa_VorBCRHS(:,1) = 0d0
          if(isDivUpdate) wa_DivBCRHS(:,1) = 0d0
       case(DynBCTYPE_SpecStress)
          if(isVorUpdate) &
               & wa_VorBCRHS(:,1) = w_AlphaOptr_xy(xy_WindStressV*xy_Coef, -xy_WindStressU*xy_Coef)/(RefDens*vViscCoef)
          if(isDivUpdate) &
               & wa_DivBCRHS(:,1) = w_AlphaOptr_xy(xy_WindStressU*xy_Coef,  xy_WindStressV*xy_Coef)/(RefDens*vViscCoef)
       case(DynBCTYPE_Slip)
          if(isVorUpdate) wa_VorBCRHS(:,1) = 0d0
          if(isDivUpdate) wa_DivBCRHS(:,1) = 0d0
       end select

       select case(DynBC_Bottom)
       case(DynBCTYPE_NoSlip)
          if(isVorUpdate) wa_VorBCRHS(:,2) = 0d0
          if(isDivUpdate) wa_DivBCRHS(:,2) = 0d0
       case(DynBCTYPE_SpecStress)
       case(DynBCTYPE_Slip)
          if(isVorUpdate) wa_VorBCRHS(:,2) = 0d0
          if(isDivUpdate) wa_DivBCRHS(:,2) = 0d0
       end select
    end if

    if (isPTempEddUpdate) then
       select case(ThermBC_Surface)
       case(ThermBCTYPE_Adiabat)
          forAll(k=0:kMax) xyz_PTempTmp(:,:,k) = z_PTempBasic(k)
          wz_PTempTmp = wz_wt(wt_DSig_wt(wt_xyz(xyz_PTempTmp)))
          wa_PTempEddBCRHS(:,1) = - wz_PTempTmp(:,0)
       case(ThermBCTYPE_FluxFixed)
          wa_PTempEddBCRHS(:,1) = w_xy(xy_SurfTempFlux)
       case(ThermBCTYPE_TempFixed)
          wa_PTempEddBCRHS(:,1) = w_xy(xy_SeaSurfTemp)
       case(ThermBCTYPE_TempRelaxed)
          wz_PTempTmp = wz_wt(wt_PTempEdd)
          wa_PTempEddBCRHS(:,1) = wz_PTempTmp(:,0)
       end select

       select case(ThermBC_Bottom)
       case(ThermBCTYPE_Adiabat)
          forAll(k=0:kMax) xyz_PTempTmp(:,:,k) = z_PTempBasic(k)
          wz_PTempTmp = wz_wt(wt_DSig_wt(wt_xyz(xyz_PTempTmp)))
          wa_PTempEddBCRHS(:,2) = - wz_PTempTmp(:,kMax)
       case(ThermBCTYPE_FluxFixed)
       case(ThermBCTYPE_TempFixed)
       case(ThermBCTYPE_TempRelaxed)
       end select
    end if

    if (isSaltUpdate) then
       select case(SaltBC_Surface)
       case(SaltBCTYPE_Adiabat)
          wa_SaltBCRHS(:,1) = 0d0
       case(SaltBCTYPE_FluxFixed)
          wa_SaltBCRHS(:,1) = w_xy(xy_SurfSaltFlux)
       case(SaltBCTYPE_SaltFixed)
          wa_SaltBCRHS(:,1) = w_xy(xy_SeaSurfSalt)
       case(SaltBCTYPE_SaltRelaxed)
          wz_SaltTmp = wz_wt(wt_Salt)
          wa_SaltBCRHS(:,1) = wz_SaltTmp(:,0)
       end select

       select case(SaltBC_Bottom)
       case(SaltBCTYPE_Adiabat)
          wa_SaltBCRHS(:,2) = 0d0
       case(SaltBCTYPE_FluxFixed)
       case(SaltBCTYPE_SaltFixed)
       case(SaltBCTYPE_SaltRelaxed)
       end select
    end if

  end subroutine Update_VBCRHS


  !> @brief 
  !!
  !!
  subroutine apply_boundaryConditions2( &
       & wt_Vor, wt_Div, wt_PTempEdd, wt_Salt, &
       & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS )
    
    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:tMax), intent(inout) :: wt_Vor, wt_Div, wt_PTempEdd, wt_Salt
    real(DP), dimension(lMax, 1:2), intent(in) :: &
         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
    
    ! 実行文; Executable statement
    !

    call apply_ZBoundaryCond( wt_Vor, &
         & inquire_VBCSpecType(DynBC_Surface), inquire_VBCSpecType(DynBC_Bottom), &
         & w_SurfBCWork=wa_VorBCRHS(:,1), w_BtmBCWork=wa_VorBCRHS(:,2) )

    call apply_ZBoundaryCond( wt_Div, &
         & inquire_VBCSpecType(DynBC_Surface), inquire_VBCSpecType(DynBC_Bottom), &
         & w_SurfBCWork=wa_DivBCRHS(:,1), w_BtmBCWork=wa_DivBCRHS(:,2) )

    call apply_ZBoundaryCond( wt_PTempEdd, &
         & inquire_VBCSpecType(ThermBC_Surface), inquire_VBCSpecType(ThermBC_Bottom), &
         & w_SurfBCWork=wa_PTempEddBCRHS(:,1), w_BtmBCWork=wa_PTempEddBCRHS(:,2) )

    call apply_ZBoundaryCond( wt_Salt, &
         & inquire_VBCSpecType(SaltBC_Surface), inquire_VBCSpecType(SaltBC_Bottom), &
         & w_SurfBCWork=wa_SaltBCRHS(:,1), w_BtmBCWork=wa_SaltBCRHS(:,2) )

  end subroutine apply_boundaryConditions2

end module HydroBoudEq_TimeInteg_mod

