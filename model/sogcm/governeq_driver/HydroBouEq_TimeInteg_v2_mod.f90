!-------------------------------------------------------------
! Copyright (c) 2013-2015 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBoudEq_TimeInteg_v2_mod 

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
       & xyz_Lat, xyz_Lon, z_LyrThickSig

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

  use HydroBouEqSolverRHS_v2_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final, &
       & calc_HydroBouEqRHSExpl_xyz, calc_HydroBouEqVViscRHS_xyz, &
       & Advance_BarotEqStep, Advance_BarotEqStep2

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
         & xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic, &
         & xyz_VViscCoefA, xyz_VViscCoefB, xyz_VViscCoefN, &
         & xyz_VDiffCoefA, xyz_VDiffCoefB, xyz_vDiffCoefN

    use BoundaryCondO_mod, only: &
         & apply_VBoundaryCondO

    use TemporalIntegSet_mod, only: &
         & CurrentTime
    
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
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHeight
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_ForceUBaroc, xy_ForceVBaroc, &
         & xy_UBarot, xy_VBarot, xy_SurfPress, &
         & xy_UBarotOld, xy_VBarotOld, xy_SurfPressOld
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTempBasic
    
    ! Work variables for some temporal schemes, such as Runge=Kutta scheme etc. 
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U_RHSTmp, xyz_V_RHSTmp, xyz_PTemp_RHSTmp, xyz_Salt_RHSTmp
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHeight_RHSTmp

    !
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHeight_RHSEx

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHeight_RHSIm

    real(DP), dimension(lMax,0:kMax) :: wz_Vor, wz_Div
    
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
           

    call apply_VBoundaryCondO( &
         & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_SaltN,     & ! (inout)
         & xyz_VViscCoefN, xyz_VDiffCoefN                & ! (in)
         & )
    
    !
    !
    !
    !$omp parallel
    !$OMP workshare
    xyz_U(:,:,:) = xyz_UN; xyz_V(:,:,:) = xyz_VN;
    xyz_PTempEdd(:,:,:) = xyz_PTempEddN; xyz_Salt(:,:,:) = xyz_SaltN
    xy_SurfHeight(:,:) = xy_SurfHeightN
    xy_SurfPress(:,:) = xy_SurfPressN
    xy_SurfPressOld(:,:) = xy_SurfPressN
!    xyz_VViscCoef(:,:,:) = xyz_VViscCoefN
!    xyz_VDiffCoef(:,:,:) = xyz_VDiffCoefN    
    !$omp end workshare
    !$omp end parallel

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call TemporalIntegUtil_SetDelTime(DelTime)
    
    
    !
    do Stage=1, nStage_BarocTimeInt

       DelTauBaroc = TemporalIntegUtil_GetDDtCoef(timeIntMode, Stage)*DelTime
       call HydroBouEqSolverVImplProc_Prepare( &
            & DelTauBaroc, CoriolisTermACoef, VDiffTermACoef, VDiffTermACoef &  ! (in)
            & )

       !$omp parallel
       !$omp workshare
       xyz_U_RHSEx = 0d0; xyz_V_RHSEx = 0d0; xyz_PTemp_RHSEx = 0d0; xyz_Salt_RHSEx = 0d0
       xyz_U_RHSIm = 0d0; xyz_V_RHSIm = 0d0; xyz_PTemp_RHSIm = 0d0; xyz_Salt_RHSIm = 0d0
       !$omp end workshare
       !$omp end parallel 
       xy_ForceUBaroc(:,:) = 0d0; xy_ForceVBaroc(:,:) = 0d0


       call calc_IntSig_BtmToTop(xy_UBarot, xyz_U)
       call calc_IntSig_BtmToTop(xy_VBarot, xyz_V)
       
       select case(timeIntMode) !======================================================================      
       case(timeIntMode_Euler)  !***** Using Euler scheme *******************

          call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D' ); 
          call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
          call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'D', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 
          call timeInt_Euler()
          call perform_adjustmentProcess(xyz_PTempEdd, xyz_Salt, xy_SurfHeight)
          call BarotEqStep(xyz_U, xyz_V)
          
       case(timeIntMode_RK2)    !***** Using RK2 scheme   *******************
          xy_SurfPressOld(:,:) = xy_SurfPressN
          call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D' ); 
          call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
          call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'D', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 
          call timeInt_RK(Stage,2)
          call perform_adjustmentProcess(xyz_PTempEdd, xyz_Salt, xy_SurfHeight)
          call BarotEqStep(xyz_U, xyz_V)
          
       case(timeIntMode_RK4)    !***** Using RK4 scheme   *******************
          xy_SurfPressOld(:,:) = xy_SurfPressN
          xy_UBarotOld(:,:) = xy_IntSig_BtmToTop_xyz(xyz_UN)
          xy_VBarotOld(:,:) = xy_IntSig_BtmToTop_xyz(xyz_VN)

          call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D', 'N' );
          call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
          call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'D', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 

          call timeInt_RK(Stage,4)
          call perform_adjustmentProcess(xyz_PTempEdd, xyz_Salt, xy_SurfHeight)
          call BarotEqStep(xyz_U, xyz_V)
          
       case(timeIntMode_LF)  !***** Using Leap-Frog scheme  *******************

          xy_SurfPress = 1.5d0*xy_SurfPressN - 0.5d0*xy_SurfPressB !1d0/3d0*xy_SurfPressN + 1d0/3d0*xy_SurfPressB
          
          call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D', 'B'); 
          call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
          call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'B', 1d0, .false.);
          call add_ImplRHS_into_RHS(); 

          call timeInt_LF(Stage)
          call perform_adjustmentProcess(xyz_PTempEdd, xyz_Salt, xy_SurfHeight)
          call BarotEqStep(xyz_U, xyz_V)
          
       case(timeIntMode_LFAM3)  !***** Using LFAM3 scheme  *******************

          isVImplicitProc = .true.
          
          select case(Stage)
          case(1)
             call calc_IntSig_BtmToTop(xy_UBarotOld, xyz_UB)
             call calc_IntSig_BtmToTop(xy_VBarotOld, xyz_VB)

             xy_SurfPressOld(:,:) = xy_SurfPressB
!             xy_SurfPress = 1.5d0*xy_SurfPressN - 0.5d0*xy_SurfPressB !1d0/3d0*xy_SurfPressN + 1d0/3d0*xy_SurfPressB
             
             call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D', 'B'); 
             call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
             call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'B', 1d0, .false.);
          case(2)
             call calc_IntSig_BtmToTop(xy_UBarotOld, xyz_UN)
             call calc_IntSig_BtmToTop(xy_VBarotOld, xyz_VN)
             xy_SurfPressOld(:,:) = xy_SurfPressN
!             xy_SurfPress = 0.5d0*(xy_SurfPressA + xy_SurfPressN)!2d0*xy_SurfPressN - xy_SurfPressB
             
             call calc_ExplRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, xy_SurfHeight_RHSEx, 'D', 'N'); 
!!$          if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$             write(*,*) "RHSEx=",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTemp_RHSEx)*xy_totDepthBasic), &
!!$                  AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSEx)*xy_totDepthBasic)
!!$          end if
             call calc_ExplTermWithPhysicsRHS(xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx, 'D', .true. )
!!$          if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$             write(*,*) "RHSEx=",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTemp_RHSEx)*xy_totDepthBasic), &
!!$                  AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSEx)*xy_totDepthBasic)
!!$          end if
             call calc_VViscRHS(xyz_U_RHSIm, xyz_V_RHSIm, xyz_PTemp_RHSIm, xyz_Salt_RHSIm, 'N', 1d0, .false.);
          end select          
!!$          if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$             write(*,*) "--> RHSEx=",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTemp_RHSEx)*xy_totDepthBasic), &
!!$                                    AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSEx)*xy_totDepthBasic)
!!$             write(*,*) "--> RHSIm=",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTemp_RHSIm)*xy_totDepthBasic), &
!!$                                    AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSIm)*xy_totDepthBasic)
!!$          end if

          call add_ImplRHS_into_RHS();
          call timeInt_LFAM3(Stage)
       end select              !=====================================================================================

       !

!!$    if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$       write(*,*) "Before VBC -----------------------"
!!$       write(*,*) "PTempEdd=", &
!!$            RefDens*Cp0*AvrLonLat_xy( (xyz_PTempEdd(:,:,0)-xyz_PTempEddN(:,:,0))*z_LyrThickSig(0)*xy_totDepthBasic(:,:))/DelTime, &
!!$            RefDens*Cp0*AvrLonLat_xy( (xyz_PTempEdd(:,:,kMax)-xyz_PTempEddN(:,:,kMax))*z_LyrThickSig(kMax)*xy_totDepthBasic(:,:))/DelTime
!!$       call check_CpT()
!!$    end if

       call DealiasingStep()

       call apply_VBoundaryCondO( &
            & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt,     & ! (inout)
            & xyz_VViscCoefN, xyz_VDiffCoefN            & ! (in)
            & )


!!$    if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$       write(*,*) "After VBC -----------------------"
!!$       write(*,*) "PTempEdd=", &
!!$            RefDens*Cp0*AvrLonLat_xy( (xyz_PTempEdd(:,:,0)-xyz_PTempEddN(:,:,0))*z_LyrThickSig(0)*xy_totDepthBasic(:,:))/DelTime, &
!!$            RefDens*Cp0*AvrLonLat_xy( (xyz_PTempEdd(:,:,kMax)-xyz_PTempEddN(:,:,kMax))*z_LyrThickSig(kMax)*xy_totDepthBasic(:,:))/DelTime
!!$       call check_surfflx(xyz_PTempEdd)
!!$       call check_CpT()
!!$    end if

       !
       
    end do  ! End of do loop for a multi-stage temporal scheme.

    
    ! ** !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !
    !
    !$omp parallel
    !$omp workshare
    xyz_UA(:,:,:) = xyz_U
    xyz_VA(:,:,:) = xyz_V
    xyz_PTempEddA(:,:,:) = xyz_PTempEdd
    xyz_SaltA(:,:,:) = xyz_Salt
    xy_SurfHeightA(:,:) = xy_SurfHeight
    !$omp end workshare
    !$omp end parallel
    
!    xyz_VViscCoefA(:,:,:) = xyz_VViscCoef
    !    xyz_VDiffCoefA(:,:,:) = xyz_VDiffCoef

    if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
       write(*,*) "After -----------------------"
!       call check_surfflx(xyz_PTempEdd)
       call check_CpT()
       call check_Salt()
!       write(*,*) "************************************"
    end if
    
    contains
      subroutine check_CpT()
        use BoundaryCondO_mod, only: xy_SurfHFlxO

        !        real(DP), intent(in) :: xyz_PTempEdd(0:iMax-1,jMax,0:kMax)
        real(DP) :: heatbudget

        heatbudget = RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTempEdd - xyz_PTempEddN)*xy_totDepthBasic)
        write(*,*) "check CpT=", heatbudget/DelTime, -AvrLonLat_xy(xy_SurfHFlxO)
      end subroutine check_CpT
      subroutine check_Salt()
        use BoundaryCondO_mod, only: xy_SurfFwFlxO

        !        real(DP), intent(in) :: xyz_PTempEdd(0:iMax-1,jMax,0:kMax)
        real(DP) :: saltbudget
        real(DP), parameter :: RefSalt = 35d0

        saltbudget = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt - xyz_SaltN)*xy_totDepthBasic)
        write(*,*) "check Salt=", saltbudget/DelTime, -AvrLonLat_xy(xy_SurfFwFlxO)*RefSalt
      end subroutine check_Salt

      subroutine check_surfflx(xyz_PTempEdd)


        use BoundaryCondO_mod, only: xy_SurfHFlxO
        
        real(DP), intent(in) :: xyz_PTempEdd(0:iMax-1,jMax,0:kMax)
         
        real(DP) :: xyz_DSigPTemp(0:iMax-1,jMax,0:kMax)
        
        xyz_DSigPTemp = xyz_DSig_xyz(xyz_PTempEdd)
         write(*,*) "--", &
              & AvrLonLat_xy(RefDens*Cp0*xyz_VDiffCoefN(:,:,0)*xyz_DSigPTemp(:,:,0)/xy_totDepthBasic), &
              & "glmean given surfHFlx:", AvrLonLat_xy(xy_SurfHFlxO)

      end subroutine check_surfflx
      
      subroutine calc_ExplRHS( &
           & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, xy_SurfHeightRHS, & 
           & tLevel, SemiImplicit_tLevel &
           & )

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
             & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
        real(DP), dimension(0:iMax-1,jMax), intent(inout) :: xy_SurfHeightRHS
        character, intent(in), optional :: tLevel, SemiImplicit_tLevel

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        integer :: k
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_UrfCori, xyz_VrfCori
        real(DP), dimension(0:iMax-1,jMax) :: xy_UBarotCori, xy_VBarotCori
        
        ! 実行文; Executable statement
        !
        
        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        if (isCoriTermSemiImplicit .and. present(SemiImplicit_tLevel)) then
           select case(SemiImplicit_tLevel)
           case('B')
             call calc_IntSig_BtmToTop( xy_UBarotCori, xyz_UB)
             call calc_IntSig_BtmToTop( xy_VBarotCori, xyz_VB)
              xyz_UrfCori(:,:,:) = xyz_UB; xyz_VrfCori(:,:,:) = xyz_VB
           case('N')
             call calc_IntSig_BtmToTop( xy_UBarotCori, xyz_UN)
             call calc_IntSig_BtmToTop( xy_VBarotCori, xyz_VN)
              xyz_UrfCori(:,:,:) = xyz_UN; xyz_VrfCori(:,:,:) = xyz_VN
           case default
              stop
           end select
        else
           call calc_IntSig_BtmToTop( xy_UBarotCori, xyz_U)
           call calc_IntSig_BtmToTop( xy_VBarotCori, xyz_V)
           xyz_UrfCori(:,:,:) = xyz_U; xyz_VrfCori(:,:,:) = xyz_V
        end if
        
        xyz_UrfCori = (xyz_UrfCori - spread(xy_UBarotCori,3,kMax+1))*xyz_CosLat
        xyz_VrfCori = (xyz_VrfCori - spread(xy_VBarotCori,3,kMax+1))*xyz_CosLat
        xy_ForceUBaroc = 0d0!- 2d0*Omega*sin(xyz_Lat(:,:,0))*xy_VBarotCori
        xy_ForceVBaroc = 0d0!  2d0*Omega*sin(xyz_Lat(:,:,0))*xy_UBarotCori
        
        !      
        !
        select case(tLvl)
           case('B')
              stop
           case('N')
              stop
           case default
              call calc_HydroBouEqRHSExpl_xyz( &
                   & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, xy_SurfHeightRHS, &  ! (out)
                   & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt, xy_SurfHeight,             &  ! (in)
                   & 0d0*xy_SurfPress, xy_totDepthBasic, z_PTempBasic,                    &  ! (in)
                   & hViscCoef, hHyperViscCoef, hDiffCoef, hHyperDiffCoef,            &  ! (in)
                   & xyz_UrfCori, xyz_VrfCori                                         &  ! (in)
                   & )
       end select

      end subroutine calc_ExplRHS


      ! Evaluation of the vertical viscid term. 
      ! \[
      !   dq/dt = \theta F(q^m), 
      ! \]
      ! where $F(q^m)$ is vertical viscid term, q is an arbitary physical quantity, 
      ! and \theta is a coefficent of $F(q^m)$.
      ! 
      subroutine calc_VViscRHS( &
           & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, & 
           & tLevel, vViscTermCoef, isRHSAppend)

        !

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
             & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
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
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS,                 & ! (inout)  
                & xyz_UB, xyz_VB, xyz_PTempBasic/theta + xyz_PTempEddB, xyz_SaltB,                          & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,   & ! (in)
                & isRHSReplace=isRHSReplace )     
        case('N')
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS,                 & ! (inout)  
                & xyz_UN, xyz_VN, xyz_PTempBasic/theta + xyz_PTempEddN, xyz_SaltN,                          & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,   & ! (in)
                & isRHSReplace=isRHSReplace )     
        case default 
           call calc_HydroBouEqVViscRHS_xyz( xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS,                 & ! (inout)  
                & xyz_U, xyz_V, xyz_PTempBasic/theta + xyz_PTempEdd, xyz_Salt,                              & ! (in)
                & theta*xyz_VViscCoefN, theta*vHyperViscCoef, theta*xyz_VDiffCoefN, theta*vHyperDiffCoef,     & ! (in)
                & isRHSReplace=isRHSReplace )     
        end select
        
                
      end subroutine calc_VViscRHS

      subroutine calc_ExplTermWithPhysicsRHS( &
           & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, &  ! (inout)
           & tLevel, isRHSAppend                            &  ! (in)
           & )

        use SGSEddyMixing_mod, only: SGSEddyMixing_AddMixingTerm
        use TemporalIntegSet_mod, only: DelTime, CurrentTime
        use VariableSet_mod, only: xyz_ConvIndex
        use SGSConvAdjust_mod, only: &
             & SGSConvAdjust_perform

        ! 宣言文; Declaration statement
        !
        real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
             & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
        character, intent(in), optional :: tLevel
        logical, intent(in), optional :: isRHSAppend

        ! 局所変数
        ! Local variables
        !
        character :: tLvl
        logical :: isRHSReplace
        integer :: k
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTempTmp, xyz_SaltTmp
        logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_adjustedFlag
        integer :: nTStep

        ! 実行文; Executable statement
        !

        tLvl = 'D'
        if (present(tLevel)) tLvl = tLevel

        if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) then
           call SGSEddyMixing_AddMixingTerm(xyz_PTempRHS, xyz_SaltRHS, &
                & xyz_PTempBasic + xyz_PTempEdd, xyz_Salt, xy_totDepthBasic + xy_SurfHeight)
        end if

!!$
        if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_CONVADJUST_NAME)) then

           xyz_PTempTmp = xyz_PTempEdd + xyz_PTempBasic
           xyz_SaltTmp = xyz_Salt
           call SGSConvAdjust_perform( xyz_PTempTmp, xyz_SaltTmp, &
                & xy_totDepthBasic+xy_SurfHeight, xyz_adjustedFlag )

           !$omp parallel
           !$omp workshare
           xyz_PTempRHS = xyz_PTempRHS + &
                (xyz_PTempTmp - xyz_PTempEdd - xyz_PTempBasic)/DelTauBaroc
           xyz_SaltRHS = xyz_SaltRHS + &
                (xyz_SaltTmp - xyz_Salt)/DelTauBaroc
           !$omp end workshare
           !$omp end parallel
        end if

      end subroutine calc_ExplTermWithPhysicsRHS

      subroutine perform_adjustmentProcess( &
           & xyz_PTempEdd, xyz_Salt,        & ! (inout)
           & xy_SurfHeight                  & ! (in)
           & )

        use TemporalIntegSet_mod, only: DelTime, CurrentTime
        use VariableSet_mod, only: xyz_ConvIndex

        use SGSConvAdjust_mod, only: &
             & SGSConvAdjust_perform

        use SGSSlowConvAdjust_mod, only: &
             & SGSSlowConvAdjust_perform

        ! 宣言文; Declaration statement
        !        
        real(DP), dimension(0:iMax-1,jMax, 0:kMax), intent(inout) :: &
             & xyz_PTempEdd, xyz_Salt
        real(DP), dimension(0:iMax-1, jMax), intent(in) :: xy_SurfHeight

        ! 局所変数
        ! Local variables
        !
        integer :: k
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTemp
        logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_adjustedFlag
        integer :: nTStep

        ! 実行文; Executable statement
        !

        !
        if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_CONVADJUST_NAME)) then

           xyz_PTemp(:,:,:) = xyz_PTempEdd + xyz_PTempBasic

           call SGSConvAdjust_perform( xyz_PTemp, xyz_Salt, &
                & xy_totDepthBasic+xy_SurfHeight, xyz_adjustedFlag )

!!$           call SGSSlowConvAdjust_perform( xyz_PTemp, xyz_Salt, &
!!$                & xy_totDepthBasic+xy_SurfHeight, xyz_adjustedFlag )
!!$
           !
           nTStep = CurrentTime/DelTime
           !$omp parallel
           !$omp workshare
           where(xyz_adjustedFlag)
              xyz_ConvIndex = (xyz_ConvIndex*nTStep + 1d0)/(nTStep + 1d0)
           elsewhere
              xyz_ConvIndex = xyz_ConvIndex*nTStep/(nTStep + 1d0)              
           end where

           xyz_PTempEdd = xyz_PTemp - xyz_PTempBasic
           !$omp end workshare
           !$omp end parallel
        end if

      end subroutine perform_adjustmentProcess

      subroutine add_ImplRHS_into_RHS()

        ! 実行文; Executable statement
        !

        !$omp parallel
        !$omp workshare
        xyz_U_RHSEx(:,:,:) = xyz_U_RHSEx + xyz_U_RHSIm
        xyz_V_RHSEx(:,:,:) = xyz_V_RHSEx + xyz_V_RHSIm
        xyz_PTemp_RHSEx(:,:,:) = xyz_PTemp_RHSEx + xyz_PTemp_RHSIm
        xyz_Salt_RHSEx(:,:,:)  = xyz_Salt_RHSEx + xyz_Salt_RHSIm
        !$omp end workshare
        !$omp end parallel

!        if( .not. isVImplicitProc ) then
           xy_ForceUBaroc(:,:) = xy_ForceUBaroc + xy_IntSig_BtmToTop_xyz(xyz_U_RHSEx)
           xy_ForceVBaroc(:,:) = xy_ForceVBaroc + xy_IntSig_BtmToTop_xyz(xyz_V_RHSEx)
!        end if
        
      end subroutine add_ImplRHS_into_RHS

      subroutine BarotEqStep( &
           & xyz_U, xyz_V                           & ! (inout)
           & )

        ! 宣言文; Declaration statement
        !        
        real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_U, xyz_V

        ! 局所変数
        ! Local variables
        !
        real(DP), dimension(0:iMax-1,jMax) :: &
             & xy_UBarotA, xy_VBarotA, xy_CoriImpFac, xy_UTmp
        integer :: k
        
        ! 実行文; Executable statement
        !
        
        call Advance_BarotEqStep( &
             & xy_UBarotA, xy_VBarotA, xy_SurfPressA,                     &  ! (out)
             & xy_UBarot, xy_VBarot, xy_SurfPress,                        &
             & xy_UBarotOld, xy_VBarotOld, xy_SurfPressOld,               &  ! (in)
             & xy_ForceUBaroc, xy_ForceVBaroc,                            &  ! (in)             
             & DelTauBaroc, CoriolisTermACoef                             &  ! (in)
             & )

        xy_UBarotA = xy_UBarotA - xy_IntSig_BtmToTop_xyz(xyz_U)
        xy_VBarotA = xy_VBarotA - xy_IntSig_BtmToTop_xyz(xyz_V)
        !$omp parallel
        !$omp do
        do k=0, kMax
           xyz_U(:,:,k) = xyz_U(:,:,k) + xy_UBarotA
        end do
        !$omp do
        do k=0, kMax
           xyz_V(:,:,k) = xyz_V(:,:,k) + xy_VBarotA
        end do
        !$omp end parallel

!!$        call Advance_BarotEqStep2( &
!!$             & xyz_U, xyz_V, xy_SurfPressA, xy_SurfPress, DelTauBaroc)
      end subroutine BarotEqStep

        
      !* Interfaces for calling implicit scheme. 
      !
      subroutine ImplicitProc( &
           & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt &  ! (out)
           & )

        ! 宣言文; Declaration statement
        !
        real(DP), intent(out), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt
        

        !
        !
        real(DP), dimension(0:iMax-1,jMax) :: xy_DUBarot, xy_DVBarot, xy_DUTmp, xy_CoriImpFac
        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_U_RHSEx_, xyz_V_RHSEx_
        integer :: k
        real(DP) :: xyz_Tmp(0:iMax-1,jMax,0:kMax) 

        ! 実行文; Executable statement
        !
        if( .not. isVImplicitProc ) then
           call MessageNotify('E', module_name, &
                & 'VImplicitProcFlag is set to .false., but ImplicitProc has been called.')
        end if

!!$          if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$             write(*,*) "BEFORE Impl."
!!$             write(*,*) "--> RHS=",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_PTemp_RHSEx)*xy_totDepthBasic),  &
!!$                  AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSEx)*xy_totDepthBasic)
!!$             
!!$          end if

        call Advance_VImplicitProc_DeltaForm_xyz( &
             & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt, xy_SurfPressA,             & ! (out)
             & xyz_U_RHSEx, xyz_V_RHSEx, xyz_PTemp_RHSEx, xyz_Salt_RHSEx,           & ! (inout)
             & xyz_VViscCoefN, xyz_VDiffCoefN, xy_totDepthBasic + xy_SurfHeight,                           & ! (in)
             & DynBC_Surface, DynBC_Bottom, ThermBC_Surface, ThermBC_Bottom, SaltBC_Surface, SaltBC_Bottom &  !(in)
             & )

        call calc_IntSig_BtmToTop( xy_DUBarot, xyz_DU )
        call calc_IntSig_BtmToTop( xy_DVBarot, xyz_DV )

        xy_ForceUBaroc(:,:) = xy_DUBarot / DelTauBaroc 
        xy_ForceVBaroc(:,:) = xy_DVBarot / DelTauBaroc 
        xy_CoriImpFac(:,:) = CoriolisTermACoef*DelTauBaroc*2d0*Omega*sin(xyz_Lat(:,:,0))        

        !$omp parallel do private(xy_DUTmp)
        do k=0, kMax
           xyz_DU(:,:,k) = xyz_DU(:,:,k) - xy_DUBarot
           xyz_DV(:,:,k) = xyz_DV(:,:,k) - xy_DVBarot
           
           xy_DUTmp = xyz_DU(:,:,k)
           xyz_DU(:,:,k) = (xyz_DU(:,:,k) + xy_CoriImpFac*xyz_DV(:,:,k))/(1d0 + xy_CoriImpFac**2)
           xyz_DV(:,:,k) = (xyz_DV(:,:,k) - xy_CoriImpFac*xy_DUTmp(:,:))/(1d0 + xy_CoriImpFac**2)
        end do
        
      end subroutine ImplicitProc

       ! * Interfaces for calling a subroutine to perform a temporal 
       !   integration of ODE for prognostic variables. 
       !

       subroutine timeInt_Euler()
         xyz_U(:,:,:)        = timeIntEuler( xyz_UN, xyz_U_RHSEx )
         xyz_V(:,:,:)        = timeIntEuler( xyz_VN, xyz_V_RHSEx )
         xyz_PTempEdd(:,:,:) = timeIntEuler( xyz_PTempEddN, xyz_PTemp_RHSEx )
         xyz_Salt(:,:,:)     = timeIntEuler( xyz_SaltN, xyz_Salt_RHSEx )
         xy_SurfHeight(:,:)  = timeIntEuler( xy_SurfHeightN, xy_SurfHeight_RHSEx )
       end subroutine timeInt_Euler
       
       subroutine timeInt_RK(RKStage, RKOrder)
         integer, intent(in) :: RKStage, RKOrder

         xyz_U(:,:,:)        = timeIntRK( xyz_UN, xyz_U_RHSEx, RKOrder, RKStage, xyz_U_RHSTmp )
         xyz_V(:,:,:)        = timeIntRK( xyz_VN, xyz_V_RHSEx, RKOrder, RKStage, xyz_V_RHSTmp )
         xyz_PTempEdd(:,:,:) = timeIntRK( xyz_PTempEddN, xyz_PTemp_RHSEx, RKOrder, RKStage, xyz_PTemp_RHSTmp )
         xyz_Salt(:,:,:)     = timeIntRK( xyz_SaltN, xyz_Salt_RHSEx, RKOrder,  RKStage, xyz_Salt_RHSTmp )
         xy_SurfHeight(:,:)  = timeIntRK( xy_SurfHeightN, xy_SurfHeight_RHSEx, RKOrder, RKStage, xy_SurfHeight_RHSTmp )

       end subroutine timeInt_RK

       subroutine timeInt_LFAM3(Stage)

         ! 宣言文; Declaration statement
         !
         integer, intent(in) :: Stage

         ! 局所変数
         ! Local variables
         !
         real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
              & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt
         integer :: k

         !
         call ImplicitProc(xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt)
         
         select case (Stage)
         case (1)
            !$omp parallel 
            !$omp workshare
            xyz_U = xyz_UB + xyz_DU
            xyz_V = xyz_VB + xyz_DV
            xyz_PTempEdd = xyz_PTempEddB + xyz_DPTempEdd
            xyz_Salt = xyz_SaltB + xyz_DSalt
            !$omp end workshare
            !$omp end parallel
            xy_SurfHeight = 0d0
         case(2)
            !$omp parallel 
            !$omp workshare
            xyz_U = xyz_UN + xyz_DU
            xyz_V = xyz_VN + xyz_DV
            xyz_PTempEdd = xyz_PTempEddN + xyz_DPTempEdd
            xyz_Salt = xyz_SaltN + xyz_DSalt
            !$omp end workshare
            !$omp end parallel
            xy_SurfHeight = 0d0
         end select
         
         call BarotEqStep(xyz_U, xyz_V)

!!$        if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$           write(*,*) "After Implicit -----------------------"
!!$           write(*,*) "CpT==",  RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_DPTempEdd/DelTauBaroc)*xy_totDepthBasic), &
!!$                AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_DSalt/DelTauBaroc)*xy_totDepthBasic)
!!$           write(*,*) "CpT check2==",  &
!!$                RefDens*Cp0*AvrLonLat_xy(xy_IntSig_BtmToTop_xyz( xyz_PTemp_RHSex )*xy_totDepthBasic ),  &
!!$                AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_Salt_RHSex/DelTauBaroc)*xy_totDepthBasic)
!!$        end if
!!$
!!$
!!$        if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$           write(*,*) "Before Conv -----------------------"
!!$           call check_CpT()
!!$        end if
!!$         write(*,*) "Skip conv adjustment.."
!!$         call perform_adjustmentProcess(xyz_PTempEdd, xyz_Salt, &   ! (inout)
!!$              & xy_SurfHeight)                                      ! (in)

!!$        if ( mod(CurrentTime, 3600d0*24d0*10d0) == 0) then
!!$           write(*,*) "After Conv -----------------------"
!!$           call check_CpT()
!!$        end if

         
         !
         !
         if(Stage == 1) then
            !$omp parallel
            !$omp workshare
            xyz_U = (5d0*xyz_U + 8d0*xyz_UN - xyz_UB)/12d0
            xyz_V = (5d0*xyz_V + 8d0*xyz_VN - xyz_VB)/12d0
            xyz_PTempEdd = (5d0*xyz_PTempEdd + 8d0*xyz_PTempEddN - xyz_PTempEddB)/12d0
            xyz_Salt = (5d0*xyz_Salt + 8d0*xyz_SaltN - xyz_SaltB)/12d0
            !$omp end workshare
            !$omp end parallel
            xy_SurfHeight(:,:)  = 0d0
         end if
       end subroutine timeInt_LFAM3
       
       subroutine timeInt_LF(Stage)

         ! 宣言文; Declaration statement
         !
         integer, intent(in) :: Stage

         ! 局所変数
         ! Local variables
         !
         real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
              & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt
         
         !
         call ImplicitProc(xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt)

         !$omp parallel
         !$omp sections
         !$omp section
         xyz_U(:,:,:)        = timeIntLF_IMEX(xyz_UB, xyz_DU, Stage)
         !$omp section
         xyz_V(:,:,:)        = timeIntLF_IMEX(xyz_VB, xyz_DV, Stage)
         !$omp section
         xyz_PTempEdd(:,:,:) = timeIntLF_IMEX(xyz_PTempEddB, xyz_DPTempEdd, Stage)
         !$omp section
         xyz_Salt(:,:,:)     = timeIntLF_IMEX(xyz_SaltB, xyz_DSalt, Stage)
         !$omp section
         xy_SurfHeight(:,:)  = 0d0
         !$omp end sections
         !$omp end parallel

       end subroutine timeInt_LF


       subroutine DealiasingStep
         
         use GridSet_mod, only: nMax

         use w_zonal_module_sjpack

         integer :: k
         real(DP) :: wz_PTempEdd(lMax,0:kMax)
         real(DP) :: wz_Salt(lMax,0:kMax)
         real(DP) :: w_LaplaEigVal(lMax)
         real(DP) :: w_HDifCoefH(lMax)
         real(DP) :: w_HDifCoefM(lMax)
         integer, parameter :: HDOrder = 8
         real(DP), parameter :: HDEFoldTimeSec = 5*86400d0
         real(DP) :: VisCoef
         

         w_LaplaEigVal = rn(:,1) / RPlanet**2
         VisCoef = ( nMax*(NMax + 1) / RPlanet**2 )**(-HDOrder/2) / HDEFoldTimeSec
         w_HDifCoefH = - VisCoef * ( ( - W_LaplaEigVal )**(HDOrder/2) )
         w_HDifCoefM = w_HDifCoefH - VisCoef * ( - (2d0/RPlanet**2)**(HDOrder/2) )

         call wz_VectorCosLat2VorDiv( xyz_U*xyz_CosLat, xyz_V*xyz_CosLat, &
              wz_Vor, wz_Div )
         !$omp parallel do
         do k=0, kMax
            wz_Vor(:,k) = (1d0 / (1d0 - DelTauBaroc * w_HDifCoefM))*wz_Vor(:,k)
            wz_Div(:,k) = (1d0 / (1d0 - DelTauBaroc * w_HDifCoefM))*wz_Div(:,k)
         end do
         call wz_VorDiv2VectorCosLat( wz_Vor, wz_Div, &
              xyz_U, xyz_V )

         !$omp parallel do
         do k=0 ,kMax
            xyz_U(:,:,k) = xyz_U(:,:,k)/xyz_CosLat(:,:,k)
            xyz_V(:,:,k) = xyz_V(:,:,k)/xyz_CosLat(:,:,k)
         end do

         wz_PTempEdd = wz_xyz(xyz_PTempEdd)
         wz_Salt = wz_xyz(xyz_Salt)
         !$omp parallel do
         do k=0, kMax
            wz_PTempEdd(:,k) = (1d0 / (1d0 - DelTauBaroc * w_HDifCoefH))*wz_PTempEdd(:,k)
            wz_Salt(:,k) = (1d0 / (1d0 - DelTauBaroc * w_HDifCoefH))*wz_Salt(:,k)
         end do
         xyz_PTempEdd = xyz_wz(wz_PTempEdd)
         xyz_Salt = xyz_wz(wz_Salt)

       end subroutine DealiasingStep


!!$
!!$       subroutine replace_RHS_with_VBCTIntRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS )
!!$
!!$         ! モジュール引用; Use statements
!!$         !
!!$         use VariableSet_mod, only: &
!!$              & xy_totDepthBasic
!!$
!!$         use BoundCondSet_mod, only: &
!!$              & SurfTempRelaxedTime, SurfSaltRelaxedTime
!!$
!!$         use BoundaryCondO_mod, only: &
!!$              & xy_SeaSurfTemp, xy_SeaSurfSalt, &
!!$              & xy_SurfHFlxO, xy_SurfFwFlxO
!!$         
!!$         use SpmlUtil_mod
!!$         
!!$         ! 宣言文; Declaration statement
!!$         !
!!$         real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
!!$         
!!$
!!$         ! 作業変数
!!$         ! Work variables
!!$         !
!!$         integer :: k
!!$         real(DP) :: restoreTimeScale, dz1
!!$         
!!$         ! 実行文; Executable statement
!!$         !
!!$
!!$         dz1 = xy_totDepthBasic(0,1)*(-g_Sig(1)*2d0)
!!$
!!$         if(ThermBC_Surface == ThermBCTYPE_TempRelaxed) then
!!$
!!$            wz_PTempRHS(:,1) = wz_PTempRHS(:,1) &
!!$                 & - ( wz_PTempEdd(:,1) - w_xy(xy_SeaSurfTemp-xyz_PTempBasic(:,:,1)) )/SurfTempRelaxedTime &
!!$                 & - 0d0*w_xy(xy_SurfHFlxO)/(RefDens*3986d0*dz1)
!!$         end if
!!$
!!$         if(ThermBC_Bottom == ThermBCTYPE_TempRelaxed) then
!!$            wz_PTempRHS(:,kMax) = 0d0
!!$         end if
!!$
!!$         !
!!$         if(SaltBC_Surface == SaltBCTYPE_SaltRelaxed) then
!!$            wz_SaltRHS(:,1) = wz_SaltRHS(:,1) &
!!$                 & - ( wz_Salt(:,1) - w_xy(xy_SeaSurfSalt) )/SurfSaltRelaxedTime &
!!$                 & + 0d0*w_xy(xy_SurfFwFlxO*35d0)/dz1
!!$         end if
!!$
!!$         if(SaltBC_Bottom == SaltBCTYPE_SaltRelaxed) then
!!$            wz_SaltRHS(:,kMax) = 0d0
!!$         end if
!!$
!!$       end subroutine replace_RHS_with_VBCTIntRHS

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

         xyz_VViscCoef(:,:,:) = xyz_VViscCoefMixLyr + 1d-3!vViscCoef
         xyz_VDiffCoef(:,:,:) = xyz_VDiffCoefMixLyr + vDiffCoef

!!$         xyz_VViscCoef(:,:,:) = vViscCoef;
!!$         xyz_VDiffCoef(:,:,:) = vDiffCoef
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
  end subroutine HydroBouEqSolver_AdvanceTStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_vViscDiffCoef_MixLyrSimple( &
       & xyz_VViscCoef, xyz_VDiffCoef,           &
       & xyz_Depth      &
       & )

    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_VViscCoef, xyz_VDiffCoef
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_Depth

    real(DP), parameter :: MixLyrDepth = 40d0
    real(DP), parameter :: LInv = 1d0/(0.1d0*MixLyrDepth)
    real(DP), parameter :: ViscfCoefMax = 2d-3
    real(DP), parameter :: DiffCoefMax = 1d-3

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Func

    !$omp parallel
    !$omp workshare
    xyz_Func(:,:,:) = 0.5d0 - atan((abs(xyz_Depth) - MixLyrDepth)*LInv)/PI
    xyz_VViscCoef(:,:,:) = ViscfCoefMax*xyz_Func
    xyz_VDiffCoef(:,:,:) = DiffCoefMax*xyz_Func
    !$omp end workshare
    !$omp end parallel

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
    
    !$omp parallel
    !$omp workshare
    xyz_VViscCoef(:,:,:) = AvRic/(1d0 + a*xyz_Ri)**n + 1d-3!vViscCoefBG
    xyz_VDiffCoef(:,:,:) = xyz_VViscCoef(:,:,:)/(1d0 + a*xyz_Ri) + vDiffCoefBG
    !$omp end workshare
    !$omp end parallel


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

    !$omp parallel
    !$omp do
    do k=1,kMax-1
       xyz_Dz_xyz(:,:,k) = &
            & (s(k)**2*xyz(:,:,k-1) - (s(k)**2-t(k)**2)*xyz(:,:,k) - t(k)**2*xyz(:,:,k+1)) &
            & /(s(k)*t(k)*(s(k) + t(k)))/xy_totDepthBasic
    end do
    
    !$omp workshare
    xyz_Dz_xyz(:,:,0) = &
         & (xyz(:,:,0) - xyz(:,:,1))/(g_Sig(0) - g_Sig(1))/xy_totDepthBasic
    xyz_Dz_xyz(:,:,kMax) = &
         & (xyz(:,:,kMax-1) - xyz(:,:,kMax))/(g_Sig(kMax-1) - g_Sig(kMax))/xy_totDepthBasic
    !$omp end workshare
    !$omp end parallel

  end function xyz_Dz_xyz

end module HydroBoudEq_TimeInteg_v2_mod

