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
       & ThermBCTYPE_FluxFixed, ThermBCTYPE_Adiabat, ThermBCTYPE_TempFixed, ThermBCTYPE_TempRelaxed, & 
       & SaltBCTYPE_FluxFixed, SaltBCTYPE_Adiabat, SaltBCTYPE_SaltFixed, SaltBCTYPE_SaltRelaxed, & 
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

  use HydroBouEqSolverRHS_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final, &
       & calc_HydroBouEqInvisRHS, calc_HydroBouEqHViscRHS, calc_HydroBouEqVViscRHS, &
       & correct_DivEqRHSUnderRigidLid, correct_DivEqRHSUnderRigidLid2, correct_DivVorEqRHSUnderRigidLid

  use HydroBouEqSolverVImplProc_mod, only: &
       & HydroBouEqSolverVImplProc_Init, HydroBouEqSolverVImplProc_Final, &
       & HydroBouEqSolverVImplProc_Prepare, &
       & Advance_VImplicitProc, Advance_VImplicitProc_DeltaForm

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

    call HydroBouEqSolverRHS_Init()
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

    call HydroBouEqSolverRHS_Final()
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

    real(DP), dimension(lMax, 0:kMax) :: wz_Vor, wz_VorN, wz_VorB
    real(DP), dimension(lMax, 0:kMax) :: wz_Div, wz_DivN, wz_DivB
    real(DP), dimension(lMax, 0:kMax) :: wz_PTempEdd, wz_PTempEddN, wz_PTempEddB
    real(DP), dimension(lMax, 0:kMax) :: wz_Salt, wz_SaltN, wz_SaltB
    real(DP) :: xy_SurfHeight(0:iMax-1, jMax)
    real(DP) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP) :: xyz_PTempBasic(0:iMax-1,jMax,0:kMax)

    real(DP), dimension(0:iMax-1, jMax, 0:kMax) :: xyz_Urf, xyz_Vrf
    real(DP), dimension(lMax, 0:kMax) ::  wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS
    real(DP) :: w_SurfHeightExplRHS(lMax)
    real(DP), dimension(lMax, 0:kMax) ::  wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp

    ! Work variables for some temporal schemes, such as Runge=Kutta scheme etc. 
    real(DP), dimension(lMax,0:kMax) :: wz_VorRHSTIntTmp, wz_DivRHSTIntTmp, wz_PTempRHSTIntTmp, &
         & wz_SaltRHSTIntTmp
    real(DP), dimension(0:iMax-1,jMax) :: xy_SurfHeightRHSTIntTmp

    integer :: Stage
    character(TOKEN) :: TIntType_SurfPressTerm
    
    real(DP) :: xyz_CosLat(0:iMax-1, jMax, 0:kMax)
    real(DP), dimension(lMax,2) :: wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
    integer :: j, k

    logical :: isCoriTermSemiImplicit
    
    ! 実行文; Executable statement
    !
    
    ! * Preparation 


    xyz_CosLat(:,:,:) = spread(xy_CosLat, 3, kMax+1)
    forAll(k=0:kMax) xyz_PTempBasic(:,:,k) = z_PTempBasic(k)

    if (CoriolisTermACoef > 0d0) then
       isCoriTermSemiImplicit = .true.
    else
       isCoriTermSemiImplicit = .false.
    end if
    
    !
    !

    ! Set some variables at the time level N. 

    xyz_Urf(:,:,:) = xyz_UN*xyz_CosLat; xyz_Vrf(:,:,:) = xyz_VN*xyz_CosLat
    call wz_VectorCosLat2VorDiv(xyz_Urf, xyz_Vrf, & !(in) 
         & wz_VorN, wz_DivN                     )   !(out)
    wz_PTempEddN(:,:) = wz_xyz(xyz_PTempEddN)
    wz_SaltN(:,:) = wz_xyz(xyz_SaltN)

    ! Set some variables at the time level B if they are used. 
    if( isVarBUsed_BarocTimeInt ) then
       call wz_VectorCosLat2VorDiv(xyz_UB*xyz_CosLat, xyz_VB*xyz_CosLat, & ! (in)
            & wz_VorB, wz_DivB                                         )   ! (out)
       wz_PTempEddB(:,:) = wz_xyz(xyz_PTempEddB)
       wz_SaltB(:,:) = wz_xyz(xyz_SaltB)
    end if

    !

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call TemporalIntegUtil_SetDelTime(DelTime)

    !

    !$omp parallel workshare
    wz_Div(:,:) = wz_DivN; wz_Vor(:,:) = wz_VorN; wz_PTempEdd(:,:) = wz_PTempEddN; wz_Salt(:,:) = wz_SaltN;
    xy_SurfHeight(:,:) = xy_SurfHeightN; xy_SurfPress(:,:) = xy_SurfPressN;
    !$omp end parallel workshare

    !
    do Stage=1, nStage_BarocTimeInt
       
       call HydroBouEqSolverVImplProc_Prepare( &
            & TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime, &
            & CoriolisTermACoef, VDiffTermACoef, VDiffTermACoef)

       select case(timeIntMode) !======================================================================      
       case(timeIntMode_Euler)  !***** Using Euler scheme *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'D', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )

          call timeInt_Euler()
          call correct_DivEqRHS_RigidLid2('CRANKNIC')
          
       case(timeIntMode_RK2)    !***** Using RK2 scheme   *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'D', 1d0, .false.); 
          call add_ImplRHS_into_RHS(); 
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )

          call timeInt_RK(Stage,2)
          call correct_DivEqRHS_RigidLid2('CRANKNIC')
          
       case(timeIntMode_RK4)    !***** Using RK4 scheme   *******************

          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, w_SurfHeightExplRHS, 'D' ); 
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'D', 1d0, .false.); 

          call add_ImplRHS_into_RHS(); 
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )

          call timeInt_RK(Stage,4)
          call correct_DivEqRHS_RigidLid2('BEuler1')

       case(timeIntMode_LF)  !***** Using Leap-Frog scheme  *******************

          xy_SurfPress = 1.5d0*xy_SurfPressN - 0.5d0*xy_SurfPressB !1d0/3d0*xy_SurfPressN + 1d0/3d0*xy_SurfPressB
          call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS,  w_SurfHeightExplRHS, 'D', 'B');
          call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
          call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
          call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'B', 1d0, .false.); 
          call add_ImplRHS_into_RHS();

          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call timeInt_LF(Stage)

       case(timeIntMode_LFAM3)  !***** Using LFAM3 scheme  *******************

          select case(Stage)
          case(1)
             xy_SurfPress = 1.5d0*xy_SurfPressN - 0.5d0*xy_SurfPressB !1d0/3d0*xy_SurfPressN + 1d0/3d0*xy_SurfPressB
             call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS,  w_SurfHeightExplRHS, 'D', 'B');
             call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
             call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
             call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'B', 1d0, .false.); 
          case(2)
             xy_SurfPress = 0.5d0*(xy_SurfPressA + xy_SurfPressN)!2d0*xy_SurfPressN - xy_SurfPressB
             call calc_InvisRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS,  w_SurfHeightExplRHS, 'D', 'N');
             call calc_HViscRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. ); 
             call calc_ExplTermWithPhysicsRHS(wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, 'D', .true. )
             call calc_VViscRHS(wz_VorImplRHSTmp, wz_DivImplRHSTmp, wz_PTempImplRHSTmp, wz_SaltImplRHSTmp, 'N', 1d0,.false.)!'N', 0.5d0, .false.); 
          end select
          
          call add_ImplRHS_into_RHS();
          call replace_RHS_with_VBCTIntRHS( wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS )
          call timeInt_LFAM3(Stage)
       end select              !=====================================================================================
       
       call update_VBCRHS( wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, wz_Vor, wz_Div, wz_PTempEdd, wz_Salt )
       call apply_boundaryConditions2(wz_Vor, wz_Div, wz_PTempEdd, wz_Salt, &
            & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS )

       call wz_VorDiv2VectorCosLat( wz_Vor, wz_Div,      &  !(in)
            & xyz_Urf, xyz_Vrf                           &  !(out)
            & )


    end do  ! End of do loop for a multi-stage temporal scheme.



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xyz_PTempEddA(:,:,:) = xyz_wz(wz_PTempEdd)
    xyz_SaltA(:,:,:) = xyz_wz(wz_Salt)

    !$omp parallel workshare
!    xy_SurfPressA = xy_SurfPress
    xyz_UA(:,:,:) = xyz_Urf/xyz_CosLat
    xyz_VA(:,:,:) = xyz_Vrf/xyz_CosLat
    !$omp end parallel workshare
    xy_SurfHeightA(:,:) = xy_SurfHeight

    !
    !

    contains
      
      subroutine calc_InvisRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, & 
           & tLevel, SemiImplicit_tLevel &
           & )

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        real(DP), dimension(lMax), intent(inout) :: w_SurfHeightRHS
        character, intent(in), optional :: tLevel, SemiImplicit_tLevel

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_UrfCori, xyz_VrfCori

        ! 実行文; Executable statement
        !
        
        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        if (isCoriTermSemiImplicit .and. present(SemiImplicit_tLevel)) then
           select case(SemiImplicit_tLevel)
           case('B')
              xyz_UrfCori(:,:,:) = xyz_UB*xyz_CosLat; xyz_VrfCori(:,:,:) = xyz_VB*xyz_CosLat; 
           case('N')
              xyz_UrfCori(:,:,:) = xyz_UN*xyz_CosLat; xyz_VrfCori(:,:,:) = xyz_VN*xyz_CosLat; 
           end select
        else
           xyz_UrfCori(:,:,:) = xyz_Urf; xyz_VrfCori(:,:,:) = xyz_Vrf
        end if

        !
        select case(tLvl)
           case('B')
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wz(wz_Vor), xyz_wz(wz_Div), xyz_wz(wz_PTempEdd), xyz_wz(wz_Salt),     &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic,                                &  ! (in)
                   & xyz_UrfCori, xyz_VrfCori &
                   & )
           case('N')
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wz(wz_Vor), xyz_wz(wz_Div), xyz_wz(wz_PTempEdd), xyz_wz(wz_Salt),    &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic,                                &  ! (in)
                   & xyz_UrfCori, xyz_VrfCori &
                   & )
           case default
              call calc_HydroBouEqInvisRHS( &
                   & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &                  ! (out)
                   & xyz_Urf, xyz_Vrf, xyz_wz(wz_Vor), xyz_wz(wz_Div), xyz_wz(wz_PTempEdd), xyz_wz(wz_Salt),    &  ! (in)
                   & xy_SurfHeight, xy_SurfPress, xy_totDepthBasic, z_PTempBasic,                                &  ! (in)
                   & xyz_UrfCori, xyz_VrfCori &
                   & )
        end select
      end subroutine calc_InvisRHS

      subroutine calc_HViscRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, & 
           & tLevel, isRHSAppend)

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        logical, intent(in), optional :: isRHSAppend

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        logical :: isRHSReplace

        ! 実行文; Executable statement
        !
        
        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        isRHSReplace = .true.
        if (present(isRHSAppend)) isRHSReplace = (.not. isRHSAppend)

        select case(tLvl)
        case('B')
           call calc_HydroBouEqHViscRHS( &
                & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,           &  ! (inout)  
                & wz_VorB, wz_DivB, wz_PTempEddB, wz_SaltB,                &  ! (in)
                & hViscCoef, hHyperViscCoef, hDiffCoef, hHyperDiffCoef,    &  ! (in)
                & isRHSReplace=isRHSReplace )

        case default
           call calc_HydroBouEqHViscRHS( &
                & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,          & ! (inout)  
                & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt,                   & ! (in)
                & hViscCoef, hHyperViscCoef, hDiffCoef, hHyperDiffCoef,   & ! (in)
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

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        real(DP), intent(in), optional :: vViscTermCoef
        logical, intent(in), optional :: isRHSAppend

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        real(DP) :: theta
        logical :: isRHSReplace

        ! 実行文; Executable statement
        !
        
        !
        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        theta = 1d0
        if (present(vViscTermCoef)) theta = vViscTermCoef

        isRHSReplace = .true.
        if (present(isRHSAppend)) isRHSReplace = (.not. isRHSAppend)

        !

        !
        select case(tLvl)
        case('B')
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,              & ! (inout)  
                & wz_VorB, wz_DivB, wz_xyz(xyz_PTempBasic/theta + xyz_PTempEddB), wz_SaltB,          & ! (in)
                & theta*vViscCoef, theta*vHyperViscCoef, theta*vDiffCoef, theta*vHyperDiffCoef,      & ! (in)
                & isRHSReplace=isRHSReplace )     
        case('N')
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,          & ! (inout)  
                & wz_VorN, wz_DivN, wz_xyz(xyz_PTempBasic/theta + xyz_PTempEddN), wz_SaltN,      & ! (in)
                & theta*vViscCoef, theta*vHyperViscCoef, theta*vDiffCoef, theta*vHyperDiffCoef,  & ! (in)
                & isRHSReplace=isRHSReplace )     
        case default 
           call calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS,          & ! (inout)  
                & wz_Vor, wz_Div, wz_xyz(xyz_PTempBasic/theta) + wz_PTempEdd, wz_Salt,           & ! (in)
                & theta*vViscCoef, theta*vHyperViscCoef, theta*vDiffCoef, theta*vHyperDiffCoef,  & ! (in)
                & isRHSReplace=isRHSReplace )
        end select

      end subroutine calc_VViscRHS

      subroutine calc_ExplTermWithPhysicsRHS( &
           & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, & 
           & tLevel, isRHSAppend)

        use SGSEddyMixing_mod, only: SGSEddyMixing_AddMixingTerm

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
        character, intent(in), optional :: tLevel
        logical, intent(in), optional :: isRHSAppend

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        logical :: isRHSReplace
        integer :: k

        ! 実行文; Executable statement
        !

        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) then
           call SGSEddyMixing_AddMixingTerm(wz_PTempRHS, wz_SaltRHS, &
                & wz_PTempEdd + wz_xyz(xyz_PTempBasic), wz_Salt, xy_totDepthBasic+xy_SurfHeight)
        end if

      end subroutine calc_ExplTermWithPhysicsRHS

      subroutine perform_adjustmentProcess(wz_PTempEdd, wz_Salt, xy_SurfHeight)

        use TemporalIntegSet_mod, only: DelTime, CurrentTime
        use VariableSet_mod, only: xyz_ConvIndex

        use SGSConvAdjust_mod, only: &
             & SGSConvAdjust_perform

        use SGSSlowConvAdjust_mod, only: &
             & SGSSlowConvAdjust_perform

        ! 宣言文; Declaration statement
        !        
        real(DP), dimension(lMax, 0:kMax), intent(inout) :: wz_PTempEdd, wz_Salt
        real(DP), dimension(0:iMax-1, jMax), intent(in) :: xy_SurfHeight

        ! 局所変数
        ! Local variables
        !
        integer :: k
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTemp, xyz_Salt
        logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_adjustedFlag
        integer :: nTStep

        ! 実行文; Executable statement
        !

        !
        if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_CONVADJUST_NAME)) then

           xyz_PTemp(:,:,:) = xyz_wz(wz_PTempEdd) + xyz_PTempBasic
           xyz_Salt(:,:,:) = xyz_wz(wz_Salt)

           call SGSConvAdjust_perform( xyz_PTemp, xyz_Salt, &
                & xy_totDepthBasic+xy_SurfHeight, xyz_adjustedFlag )


           !
           nTStep = CurrentTime/DelTime
           where(xyz_adjustedFlag)
              xyz_ConvIndex = (xyz_ConvIndex*nTStep + 1d0)/(nTStep + 1d0)
           elsewhere
              xyz_ConvIndex = xyz_ConvIndex*nTStep/(nTStep + 1d0)              
           end where

           wz_PTempEdd(:,:) = wz_xyz(xyz_PTemp - xyz_PTempBasic)
           wz_Salt(:,:) = wz_xyz(xyz_Salt)
        end if
!!$        call SGSSlowConvAdjust_perform( xyz_PTemp, xyz_Salt, &
!!$             & xy_totDepthBasic+xy_SurfHeight, xy_adjustedFlag )
!!$

      end subroutine perform_adjustmentProcess

      subroutine add_ImplRHS_into_RHS()

        ! 実行文; Executable statement
        !

        !$omp parallel workshare
        wz_VorExplRHS(:,:) = wz_VorExplRHS + wz_VorImplRHSTmp
        wz_DivExplRHS(:,:) = wz_DivExplRHS + wz_DivImplRHSTmp
        wz_PTempExplRHS(:,:) = wz_PTempExplRHS + wz_PTempImplRHSTmp
        wz_SaltExplRHS(:,:)  = wz_SaltExplRHS + wz_SaltImplRHSTmp
        !$omp end parallel workshare

      end subroutine add_ImplRHS_into_RHS

      subroutine correct_DivEqRHS_RigidLid2( &
           & TIntType_SurfPressTerm )

        ! 宣言文; Declaration statement
        !
        character(*), intent(in) :: TIntType_SurfPressTerm

        ! 実行文; Executable statement
        !

        call correct_DivEqRHSUnderRigidLid2( wz_Div, xy_SurfPressA, &
             & xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
             & 1d0 / (TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime), &
             & TIntType_SurfPressTerm )
        xy_SurfPress = xy_SurfPressN
      end subroutine correct_DivEqRHS_RigidLid2

      subroutine correct_DivVor_RigidLid( &
           & wz_Div, wz_Vor, &
           & CoriTermCoef, TIntType_SurfPressTerm )

        ! 宣言文; Declaration statement
        !        
        real(DP), intent(inout) :: wz_Div(lMax,0:kMax), wz_Vor(lMax,0:kMax)
        real(DP), intent(in) :: CoriTermCoef
        character(*), intent(in) :: TIntType_SurfPressTerm

        !
        real(DP) :: w_DDivCorrect(lMax), w_DVorCorrect(lMax)
        
        ! 実行文; Executable statement
        !
        
        call correct_DivVorEqRHSUnderRigidLid( w_DDivCorrect, w_DVorCorrect, xy_SurfPressA, &
             & wz_Div, xy_SurfPress, xy_SurfPressN, xy_SurfPressB, &
             & TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime, CoriTermCoef, &
             & TIntType_SurfPressTerm )
        xy_SurfPress = xy_SurfPressN

        wz_Div(:,:) = wz_Div + spread(w_DDivCorrect, 2, kMax+1)
        wz_Vor(:,:) = wz_Vor + spread(w_DVorCorrect, 2, kMax+1)
        
      end subroutine correct_DivVor_RigidLid

        
      !* Interfaces for calling implicit scheme. 
      !
      subroutine ImplicitProc(wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt)

        ! 宣言文; Declaration statement
        !
        real(DP), intent(out), dimension(lMax,0:kMax) :: &
             & wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt

        ! 実行文; Executable statement
        !

        call Advance_VImplicitProc_DeltaForm(wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt, xy_SurfPressA, &
             & wz_VorExplRHS, wz_DivExplRHS, wz_PTempExplRHS, wz_SaltExplRHS, &
             & xy_totDepthBasic+xy_SurfHeight,  &
             & DynBC_Surface, DynBC_Bottom, ThermBC_Surface, ThermBC_Bottom, SaltBC_Surface, SaltBC_Bottom &  !(in)
             & )
    
      end subroutine ImplicitProc

       ! * Interfaces for calling a subroutine to perform a temporal 
       !   integration of ODE for prognostic variables. 
       !

       subroutine timeInt_Euler()
         wz_Vor(:,:) = timeIntEuler( wz_VorN, wz_VorExplRHS )
         wz_Div(:,:) = timeIntEuler( wz_DivN, wz_DivExplRHS )
         wz_PTempEdd(:,:) = timeIntEuler( wz_PTempEddN, wz_PTempExplRHS )
         wz_Salt(:,:) = timeIntEuler( wz_SaltN, wz_SaltExplRHS )
         xy_SurfHeight(:,:) = timeIntEuler( xy_SurfHeightN, xy_w(w_SurfHeightExplRHS) )

         !
         call perform_adjustmentProcess(wz_PTempEdd, wz_Salt, &   ! (inout)
              & xy_SurfHeight)                                    ! (in)
         
       end subroutine timeInt_Euler
       
       subroutine timeInt_RK(RKStage, RKOrder)
         integer, intent(in) :: RKStage, RKOrder

         wz_Vor = timeIntRK( wz_VorN, wz_VorExplRHS, RKOrder, RKStage, wz_VorRHSTIntTmp )
         wz_Div = timeIntRK( wz_DivN, wz_DivExplRHS, RKOrder, RKStage, wz_DivRHSTIntTmp )
         wz_PTempEdd = timeIntRK( wz_PTempEddN, wz_PTempExplRHS, RKOrder, RKStage, wz_PTempRHSTIntTmp )
         wz_Salt = timeIntRK( wz_SaltN, wz_SaltExplRHS, RKOrder,  RKStage, wz_SaltRHSTIntTmp )
         xy_SurfHeight = timeIntRK( xy_SurfHeightN, xy_w(w_SurfHeightExplRHS), RKOrder, RKStage, xy_SurfHeightRHSTIntTmp )

         !
         call perform_adjustmentProcess(wz_PTempEdd, wz_Salt, &   ! (inout)
              & xy_SurfHeight)                                    ! (in)
         
       end subroutine timeInt_RK

       subroutine timeInt_LFAM3(Stage)

         use VariableSet_mod, only: z_PTempBasic
         ! 宣言文; Declaration statement
         !
         integer, intent(in) :: Stage

         ! 局所変数
         ! Local variables
         !
         real(DP), dimension(lMax,0:kMax) :: wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt
         
         !
         call ImplicitProc(wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt)

         !$omp parallel sections
         !$omp section
         wz_Vor(:,:)= timeIntLFAM3_IMEX(wz_VorN, wz_VorB, wz_DVor, Stage, AM3=.false.)
         !$omp section
         wz_Div(:,:) = timeIntLFAM3_IMEX(wz_DivN, wz_DivB, wz_DDiv, Stage, AM3=.false.)
         !$omp section
         wz_PTempEdd(:,:) = timeIntLFAM3_IMEX(wz_PTempEddN, wz_PTempEddB, wz_DPTempEdd, Stage, AM3=.false.)
         !$omp section
         wz_Salt(:,:) = timeIntLFAM3_IMEX(wz_SaltN, wz_SaltB, wz_DSalt, Stage, AM3=.false.)
         !$omp section
         xy_SurfHeight(:,:) = 0d0
         !$omp end parallel sections

         !
         call correct_DivVor_RigidLid( wz_Div, wz_Vor, &  ! (inout)
              & CoriolisTermACoef, 'BEuler1')             ! (in)
         
         !
         call perform_adjustmentProcess(wz_PTempEdd, wz_Salt, &   ! (inout)
              & xy_SurfHeight)                                    ! (in)

         if(Stage == 1) then
            !$omp parallel sections
            !$omp section
            wz_Vor(:,:)= timeIntLFAM3_IMEX(wz_VorN, wz_VorB, wz_Vor, Stage, AM3=.true.)
            !$omp section
            wz_Div(:,:) = timeIntLFAM3_IMEX(wz_DivN, wz_DivB, wz_Div, Stage, AM3=.true.)
            !$omp section
            wz_PTempEdd(:,:) = timeIntLFAM3_IMEX(wz_PTempEddN, wz_PTempEddB, wz_PTempEdd, Stage, AM3=.true.)
            !$omp section
            wz_Salt(:,:) = timeIntLFAM3_IMEX(wz_SaltN, wz_SaltB, wz_Salt, Stage, AM3=.true.)
            !$omp section
            xy_SurfHeight(:,:) = 0d0
            !$omp end parallel sections

!!$            write(*,*) "End Stage1"
         end if

       end subroutine timeInt_LFAM3
       
       subroutine timeInt_LF(Stage)

         ! 宣言文; Declaration statement
         !
         integer, intent(in) :: Stage

         ! 局所変数
         ! Local variables
         !

         real(DP), dimension(lMax,0:kMax) :: wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt

         call ImplicitProc(wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt)
         call correct_DivVor_RigidLid(wz_DDiv, wz_DVor, CoriolisTermACoef, 'BEuler1')

         !$omp parallel sections
         !$omp section
         wz_Vor(:,:)= timeIntLF_IMEX(wz_VorB, wz_DVor, Stage)
         !$omp section
         wz_Div(:,:) = timeIntLF_IMEX(wz_DivB, wz_DDiv, Stage)
         !$omp section
         wz_PTempEdd(:,:) = timeIntLF_IMEX(wz_PTempEddB, wz_DPTempEdd, Stage)
         !$omp section
         wz_Salt(:,:) = timeIntLF_IMEX(wz_SaltB, wz_DSalt, Stage)
         !$omp section
         xy_SurfHeight(:,:) = 0d0
         !$omp end parallel sections

         !
         call perform_adjustmentProcess(wz_PTempEdd, wz_Salt, &   ! (inout)
              & xy_SurfHeight)                                    ! (in)
         
       end subroutine timeInt_LF

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
         integer :: k
         real(DP), dimension(lMax,0:kMax) :: wz_PTemp

         ! 実行文; Executable statement
         !


         wz_PTemp(:,:) = wz_PTempEdd  + wz_xyz(xyz_PTempBasic)

         if(ThermBC_Surface == ThermBCTYPE_TempRelaxed) then
            wz_PTempRHS(:,1) =  wz_PTempRHS(:,1) &
                 & - ( wz_PTemp(:,1) - w_xy(xy_SeaSurfTemp) )/SurfTempRelaxedTime
         end if

         if(ThermBC_Bottom == ThermBCTYPE_TempRelaxed) then
            wz_PTempRHS(:,kMax) = 0d0
         end if

         !
         if(SaltBC_Surface == SaltBCTYPE_SaltRelaxed) then
            wz_SaltRHS(:,1) = wz_SaltRHS(:,1) &
                 & - ( wz_Salt(:,1) - w_xy(xy_SeaSurfSalt) )/SurfSaltRelaxedTime
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
       & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt )

    
    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & xy_WindStressU, xy_WindStressV, &
         & xy_SeaSurfTemp, xy_SurfTempFlux, &
         & xy_SeaSurfSalt, xy_SurfSaltFlux, &
         & xy_totDepthBasic, z_PTempBasic

    use BoundCondSet_mod, only: &
         & SurfTempRelaxedTime, SurfSaltRelaxedTime

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax, 1:2), intent(inout), optional :: &
         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
    real(DP), dimension(lMax,0:kMax), intent(in), optional :: wz_Vor, wz_Div, wz_PTempEdd, wz_Salt
    

    ! 作業変数
    ! Work variables
    !
    real(DP) :: xy_Coef(0:iMax-1,jMax)
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTempTmp, xyz_PTempBasic
    real(DP) :: wz_PTempTmp(lMax,0:kMax), wz_SaltTmp(lMax,0:kMax)
    integer :: k
    logical :: isVorUpdate, isDivUpdate, isPTempEddUpdate, isSaltUpdate

    ! 実行文; Executable statement
    !

    xy_Coef(:,:) = xy_totDepthBasic*xy_CosLat
    forAll(k=0:kMax) xyz_PTempBasic(:,:,k) = z_PTempBasic(k)

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
          wz_PTempTmp = wz_wt(wt_DSig_wt(wt_xyz(xyz_PTempBasic)))
          wa_PTempEddBCRHS(:,1) = - wz_PTempTmp(:,0)
       case(ThermBCTYPE_FluxFixed)
          wa_PTempEddBCRHS(:,1) = w_xy(xy_SurfTempFlux)
       case(ThermBCTYPE_TempFixed)
          wa_PTempEddBCRHS(:,1) = w_xy(xy_SeaSurfTemp  - z_PTempBasic(0))
       case(ThermBCTYPE_TempRelaxed)
          wa_PTempEddBCRHS(:,1) = w_xy(xy_SeaSurfTemp - z_PTempBasic(0))
       end select

       select case(ThermBC_Bottom)
       case(ThermBCTYPE_Adiabat)
          wz_PTempTmp = wz_wt(wt_DSig_wt(wt_xyz(xyz_PTempBasic)))
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
          wa_SaltBCRHS(:,1) = w_xy(xy_SeaSurfSalt)
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
       & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt, &
       & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS )
    
    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_Vor, wz_Div, wz_PTempEdd, wz_Salt
    real(DP), dimension(lMax, 1:2), intent(in) :: &
         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS

    !
    !
    real(DP), dimension(lMax,0:tMax) :: wt_Vor, wt_Div, wt_PTempEdd, wt_Salt

    
    ! 実行文; Executable statement
    !

    wt_Vor = wt_wz(wz_Vor)
    wt_Div = wt_wz(wz_Div)
    wt_PTempEdd = wt_wz(wz_PTempEdd)
    wt_Salt = wt_wz(wz_Salt)

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


    wz_Vor = wz_wt(wt_Vor)
    wz_Div = wz_wt(wt_Div)
    wz_PTempEdd = wz_wt(wt_PTempEdd)
    wz_Salt = wz_wt(wt_Salt)

  end subroutine apply_boundaryConditions2

end module HydroBoudEq_TimeInteg_mod

