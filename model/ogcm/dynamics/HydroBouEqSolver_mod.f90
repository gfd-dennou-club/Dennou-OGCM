!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolver_mod 

  ! モジュール引用; Use statements
  !

  use dc_types, only: &
       & DP, STRING 

  use VectorSpace_mod, only: &
       & Vector3d

  use SphericalCoord_mod, only: &
       & CartToSphPos

  use PolyMesh_mod, only: &
       & PolyMesh

  use GeometricField_mod, only: &
       & volScalarField, pointScalarField, surfaceScalarField, &
       & GeometricField_Init, GeometricField_Final, &
       & assignment(=), operator(+), operator(-), operator(*), operator(/), &
       & DeepCopy
  
  use fvCalculus_mod, only: &
       & curl, grad, div

  use Constants_mod, only: &
       & Omega, Grav

  use GridSet_mod, only: &
       & GridSet_getLocalMeshInfo, &
       & nVzLyr, nVrLyr, vHaloSize

  use TemporalIntegSet_mod, only: &
       & DelTime, SubCycleNum, &
       & nShortTimeLevel

  use EOSDriver_mod, only: &
       & EOSDriver_Eval

  use CGridFieldDataUtil_mod, only: &
       & e_c, ToTangentVel

  use VGridFieldDataUtil_mod, only: &
       & z_r, r_z, verticalInt

  use VariableSet_mod, only: &
       & VariableSet, &
       & SaltTracerID, PTempTracerID, TracerNum, &       
       & refDens

  use HydroBouEqSolverRHS_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final, &
       & calc_RHShNormVel, calc_RHSTracer, calc_RHSLyrThick

  use BarotModeTimeFilter_mod, only: &
       & BarotModeTimeFilter_Init, BarotModeTimeFilter_Final


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolver_Init, HydroBouEqSolver_Final
  public :: HydroBouEqSolver_AdvanceTStep
  public :: diagnosePress, diagnosevNormVel

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolver_mod' !< Module Name

  type(PolyMesh), pointer :: mesh

  type(volScalarField) :: zc_lyrThickMl, zc_lyrThickBMl
  type(surfaceScalarField) :: ze_hNormVelMl
  type(volScalarField), allocatable :: zc_TracersMl(:)
  type(surfaceScalarField) :: e_hMFluxTAvg1
  type(surfaceScalarField) :: e_hMFluxTAvg2
  type(volScalarField) :: c_surfHeightTAvg1
  type(volScalarField) :: c_surfHeight(nShortTimeLevel), c_surfHeightMs
  type(surfaceScalarField) :: e_hNormVel(nShortTimeLevel), e_hNormVelMs

  type(surfaceScalarField) :: e_RHShNormVelMl

  real(DP), parameter :: LFAM3_GAM = 1d0/12d0
  real(DP), parameter :: LFAM3_BETA = 0d0
  real(DP), parameter :: LFAM3_EPS = 0d0
  real(DP), parameter :: GFB_AM3BETA = 0.281105d0
  real(DP), parameter :: GFB_DELTA = 0.614d0
  real(DP), parameter :: GFB_EPS = 0.013d0
  real(DP), parameter :: GFB_GAM = 0.088d0
  

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Init(plMesh)

    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    type(PolyMesh), intent(in), target :: plMesh

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !
    
    call GeometricField_Init(zc_lyrThickMl, plMesh, "zc_lyrThickMl", vHaloSize=vHaloSize)
    call GeometricField_Init(zc_lyrThickBMl, plMesh, "zc_lyrThickBMl", vHaloSize=vHaloSize)
    call GeometricField_Init(ze_hNormVelMl, plMesh, "ze_hNormVelMl", vHaloSize=vHaloSize)
    
    allocate( zc_TracersMl(TracerNum) )
    do n=1, TracerNum
       call GeometricField_Init(zc_TracersMl(n), plMesh, "zc_TracersMl", vHaloSize=vHaloSize)
    end do

    call GeometricField_Init(e_hMFluxTAvg1, plMesh, "e_hMFluxTAvg1", vLayerNum=1)
    call GeometricField_Init(e_hMFluxTAvg2, plMesh, "e_hMFluxTAvg2", vLayerNum=1)
    call GeometricField_Init(c_surfHeightTAvg1, plMesh, "c_surfHeightTAvg1", vLayerNum=1)

    do tl=1, nShortTimeLevel
       call GeometricField_Init(e_hNormVel(tl), plMesh, "e_hNormVel", vLayerNum=1)
       call GeometricField_Init(c_surfHeight(tl), plMesh, "c_surfHeight", vLayerNum=1)
    end do

    call GeometricField_Init(e_hNormVelMs, plMesh, "e_hNormVelMs", vLayerNum=1)
    call GeometricField_Init(c_surfHeightMs, plMesh, "c_surfHeightMs", vLayerNum=1)

    call GeometricField_Init(e_RHShNormVelMl, plMesh, "e_RHShNormVelMl", vLayerNum=1)

    mesh => plMesh

    !
    call HydroBouEqSolverRHS_Init()
    call BarotModeTimeFilter_Init( DelTimeLong=DelTime, DelTimeShort=DelTime/dble(SubCycleNum) )

  end subroutine HydroBouEqSolver_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Final()

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !


    !
    call HydroBouEqSolverRHS_Final()
    call BarotModeTimeFilter_Final()

    call GeometricField_Final(zc_lyrThickMl)
    call GeometricField_Final(zc_lyrThickBMl)
    call GeometricField_Final(ze_hNormVelMl)
    
    do n=1, TracerNum
       call GeometricField_Final(zc_TracersMl(n))
    end do
    deallocate(zc_TracersMl)

    call GeometricField_Final(e_hMFluxTAvg1)
    call GeometricField_Final(e_hMFluxTAvg2)
    call GeometricField_Final(c_surfHeightTAvg1)

    do tl=1, nShortTimeLevel
       call GeometricField_Final( e_hNormVel(tl) )
       call GeometricField_Final( c_surfHeight(tl) )
    end do
    call GeometricField_Final(e_hNormVelMs)
    call GeometricField_Final(c_surfHeightMs)

    call GeometricField_Final(e_RHShNormVelMl)

    mesh => null()

  end subroutine HydroBouEqSolver_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTStep(var)
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & Al, Nl, Bl

    ! 宣言文; Declaration statement
    !
    type(VariableSet), intent(inout) :: var
    
    ! 局所変数
    ! Local variables
    !

    ! 実行文; Executable statement
    !
    

    call AdvanceBarocStep( e_RHShNormVelMl, var%ze_hNormVel(Al), &
         & var%ze_hNormVel(Nl), var%ze_hNormVel(Bl), var%zc_lyrThick(Nl), var%zc_lyrThick(Bl), &
         & var%zc_Tracers(:,Nl), var%zc_Tracers(:,Bl), &
         & var%rc_vNormVel, var%zc_Dens, var%zc_Press, var%c_surfPress, var%zc_ZMid )

    call AdvanceBarotStep( e_RHShNormVelMl, var%zc_Dens, var%c_totDepthBasic )
    
    call FinalizeBarocStep( var%ze_hNormVel(Al), var%zc_lyrThick(Al), var%zc_Tracers(:,Al), &
       & var%ze_hNormVel(Nl), var%ze_hNormVel(Bl), var%zc_lyrThick(Nl), var%zc_Tracers(:,Nl), &
       & var%rc_vNormVel, var%c_totDepthBasic, var%zc_lyrThickBasic )

    
  end subroutine HydroBouEqSolver_AdvanceTStep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine AdvanceBarocStep( e_RHShNormVelMl, ze_hNormVelAl, &
       & ze_hNormVelNl, ze_hNormVelBl, zc_lyrThickNl, zc_lyrThickBl, zc_TracersNl, zc_TracersBl, &
       & rc_vNormVel, zc_Press, zc_Dens, c_surfPress, zc_ZMid )

    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & DelTime
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: e_RHShNormVelMl
    type(surfaceScalarField), intent(inout) :: ze_hNormVelAl
    type(surfaceScalarField), intent(in) :: ze_hNormVelNl, ze_hNormVelBl
    type(volScalarField), intent(in) :: zc_lyrThickNl, zc_lyrThickBl
    type(volScalarField), intent(in) :: zc_TracersNl(:), zc_TracersBl(:)
    type(volScalarField), intent(inout) :: rc_vNormVel
    type(volScalarField), intent(inout) :: zc_Dens
    type(volScalarField), intent(inout) :: zc_Press
    type(volScalarField), intent(in) :: c_surfPress
    type(volScalarField), intent(in) :: zc_ZMid


    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: e
    integer :: n
    integer :: nCell, nEdge

    type(surfaceScalarField) :: ze_RHShNormVel
    type(volScalarField) :: zc_RHSlyrThick
    type(volScalarField) :: zc_RHSTracer
    type(surfaceScalarField) :: ze_massHFlux
    
    ! 実行文; Executable statement
    !

    call GridSet_getLocalMeshInfo(mesh, nCell=nCell, nEdge=nEdge)

    !******************************************************************
    ! Stage1: Predicotr update for the 3D momentum
    !******************************************************************

    !
    ! Advance the 3D momentum to n+1/2

    call calc_RHShNormVel( ze_RHShNormVel, &
         & ze_hNormVelNl, rc_vNormVel, zc_lyrThickNl, zc_ZMid, zc_Press, zc_Dens )

    !$omp parallel do
    do e=1, nEdge
       ze_hNormVelMl%data%v_(1:nVzLyr,e) = &
            &   (5d-1 - 2d0*LFAM3_GAM)*ze_hNormVelBl%data%v_(1:nVzLyr,e) &
            & + (5d-1 + 2d0*LFAM3_GAM)*ze_hNormVelNl%data%v_(1:nVzLyr,e) &
            & + ( 1d0 - 2d0*LFAM3_GAM)*DelTime*ze_RHShNormVel%data%v_(1:nVzLyr,e) 
    end do

    !
    ! Advance  the thickness of vertical layers to n+1/2 with pseudo-compressible algorithm

    ze_massHFlux = ze_hNormVelNl*e_c(zc_lyrThickNl)
    call calc_RHSLyrThick( zc_RHSlyrThick, &
       & ze_massHFlux, rc_vNormVel )

    !$omp parallel do
    do i=1, nCell
       zc_lyrThickBl%data%v_(1:nVzLyr,i) = zc_lyrThickNl%data%v_(1:nVzLyr,i) &
            & - (5d-1 - LFAM3_GAM)*DelTime*zc_RHSlyrThick%data%v_(1:nVzLyr,i)

       zc_lyrThickMl%data%v_(1:nVzLyr,i) = zc_lyrThickNl%data%v_(1:nVzLyr,i) &
            & + (5d-1 -LFAM3_GAM)*DelTime*zc_RHSlyrThick%data%v_(1:nVzLyr,i)
    end do

    !*
    call adjust_hVelocity(ze_hNormVelMl, e_hMFluxTAvg1, zc_lyrThickMl)


    !*******************************************************************
    ! Stage2: Predicotr update for the tracers
    !*******************************************************************

    ze_massHFlux = e_c(zc_lyrThickMl)*( &
         &   (1d0 - 2d0*LFAM3_GAM)*ze_hNormVelNl &
         & + LFAM3_BETA*(2d0*ze_hNormVelMl - 3d0*ze_hNormVelNl + ze_hNormVelAl) &
         & )

    do n=1, TracerNum

       call calc_RHSTracer( zc_RHSTracer, &
            & zc_TracersNl(n), ze_massHFlux, rc_vNormVel )

      !$omp parallel do
       do i=1, nCell
          zc_TracersMl(n)%data%v_(1:nVzLyr,i) = &
               &   (   (5d-1 - 2d0*LFAM3_GAM)*zc_TracersBl(n)%data%v_(1:nVzLyr,i) &
               &     + (5d-1 + 2d0*LFAM3_GAM)*zc_TracersNl(n)%data%v_(1:nVzLyr,i) & 
               &   )*zc_lyrThickBl%data%v_(1:nVzLyr,i)                               &
               & + (1d0 - 2d0*LFAM3_GAM)*DelTime*zc_RHSTracer%data%v_(1:nVzLyr,i)    &
               &   / zc_lyrThickMl%data%v_(1:nVzLyr,i)
       end do
    end do

    !******************************************************************
    ! Stage3: Corrector update for the 3D momentum
    !******************************************************************

    call EOSDriver_Eval( rho=zc_Dens, &
         & theta=zc_TracersMl(PTempTracerID), S=zc_TracersMl(SaltTracerID), p=zc_Press)

    call diagnosePress( zc_Press, &
         & zc_Dens, zc_lyrThickMl, c_surfPress)

    call calc_RHShNormVel( ze_RHShNormVel, &
         & ze_hNormVelMl, rc_vNormVel, zc_lyrThickMl, zc_ZMid, zc_Press, zc_Dens )

    !$omp parallel do
    do e=1, nEdge
       ze_hNormVelAl%data%v_(1:nVzLyr,e) = ze_hNormVelNl%data%v_(1:nVzLyr,e) &
            & + DelTime*ze_RHShNormVel%data%v_(1:nVzLyr,e) 
    end do

    e_RHShNormVelMl = verticalInt( ze_RHShNormVel, e_c(zc_lyrThickMl), avgFlag=.true. )

    !*******************************************************************
    ! Finalize local variables
    !*******************************************************************

    call GeometricField_Final(ze_RHShNormVel)
    call GeometricField_Final(ze_massHFlux)
    call GeometricField_Final(zc_RHSlyrThick)
    call GeometricField_Final(zc_RHSTracer)

  end subroutine AdvanceBarocStep


  !> @brief 
  !!
  !!
  subroutine AdvanceBarotStep( &
       & e_RHShNormVelMl, zc_Dens, c_totDepthBasic )
    
    ! モジュール引用; Use statements
    !
    use BarotModeTimeFilter_mod, only: &
         & nTotBarotStep, Am, Bm

    use TemporalIntegSet_mod, only: &
         & TemporalIntegSet_AdvanceShortTStep, As, Ns, Bs

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: e_RHShNormVelMl
    type(volScalarField), intent(in) :: zc_Dens
    type(volScalarField), intent(in) :: c_totDepthBasic

    ! 局所変数
    ! Local variables
    !
    type(volScalarField) :: c_DensVAvg1
    type(volScalarField) :: c_DensVAvg2
    type(volScalarField) :: rc_Dens
    type(surfaceScalarField) :: e_totDepthMs
    type(surfaceScalarField) :: e_totDepthAs
    type(surfaceScalarField) :: e_PressForce
    type(volScalarField) :: c_SurfHeightTmp
    type(volScalarField) :: c_f
    type(surfaceScalarField) :: e_f
    integer :: i
    integer :: e
    integer :: kz
    real(DP) :: tmp(nVzLyr)
    integer :: nCell

    integer :: m
    real(DP) :: DelTimeShort
    type(Vector3d) :: geoPos


    ! 実行文; Executable statement
    !

    call GridSet_getLocalMeshInfo(mesh, nCell=nCell)

    call GeometricField_Init(c_DensVAvg2, mesh, "c_DensVAvg2", vLayerNum=1)
    call GeometricField_Init(e_PressForce, mesh, "e_PressForce", vLayerNum=1)
    call GeometricField_Init(c_SurfHeightTmp, mesh, "c_SurfHeightTmp", vLayerNum=1)
    call GeometricField_Init(c_f, mesh, "c_f", vLayerNum=1)
    call GeometricField_Init(e_f, mesh, "e_f", vLayerNum=1)
 
    c_DensVAvg1 = verticalInt(zc_Dens, zc_lyrThickMl, avgFlag=.true.)
    
    rc_Dens = r_z(zc_Dens)

    !$omp parallel do private(kz, tmp, geoPos)
    do i=1, nCell
       tmp(1) = 0d0
       do kz=2, nVzLyr
          tmp(kz) = sum(zc_Dens%data%v_(1:kz-1,i)*zc_lyrThickMl%data%v_(1:kz-1,i))*zc_lyrThickMl%data%v_(kz,i)
       end do

       c_DensVAvg2%data%v_(1,i) = sum( tmp + &
            &      0.5d0*zc_lyrThickMl%data%v_(1:nVzLyr,i)**2 &
            &         *( zc_Dens%data%v_(1:nVzLyr,i) + (rc_Dens%data%v_(1:nVzLyr,i) - rc_Dens%data%v_(2:nVzLyr+1,i))/6d0 ) &
            &    )

       geoPos = CartToSphPos(mesh%CellPosList(i))
       c_f%data%v_(1,i) = 2d0*Omega*sin(geoPos%v_(2))
    end do
    
    e_f = e_c(c_f)
    call GeometricField_Final(c_f)

    !
    DelTimeShort = DelTime/dble(SubCycleNum)

    call DeepCopy(e_hNormVel(Ns), e_hMFluxTAvg1)
    call DeepCOpy(c_surfHeight(Ns), c_surfHeightTAvg1)

    e_hMFluxTAvg1 = 0d0; e_hMFluxTAvg2 = 0d0
    c_surfHeightTAvg1 = 0d0

    do m=1, nTotBarotStep

       e_totDepthMs = e_c(c_totDepthBasic + c_surfHeightMs)
       c_surfHeight(As) = c_surfHeight(Ns) - DelTimeShort*div(e_totDepthMs*e_hNormVelMs)
       
       c_SurfHeightTmp = GFB_DELTA*c_surfHeight(As) &
            & + (1d0 - GFB_DELTA - GFB_GAM - GFB_EPS)*c_surfHeight(Ns) &
            & + GFB_GAM*c_surfHeight(Bs) + GFB_EPS*c_SurfHeightTmp

       call calc_pressForce(c_SurfHeightTmp)
       e_totDepthAs = e_c( c_totDepthBasic + c_surfHeight(As) )
       e_hNormVel(As) = ( &
            &      e_c( c_surfHeight(Ns) )*e_hNormVel(Ns) &
            &    + DelTimeShort*( e_PressForce  + e_RHShNormVelMl - e_totDepthMs*e_f*ToTangentVel(e_hNormVelMs) ) &
            &  )/e_totDepthAs

       e_hMFluxTAvg1 =  e_hMFluxTAvg1 + Am(m)*e_hNormVel(As)*e_totDepthAs
       e_hMFluxTAvg2 =  e_hMFluxTAvg2 + Bm(m)*e_hNormVelMs*e_totDepthMs
       c_surfHeightTAvg1 = c_surfHeightTAvg1 + Am(m)*c_surfHeight(As)


       c_surfHeightMs = (1.5d0 + GFB_AM3BETA)*c_surfHeight(As) - (0.5d0 + 2d0*GFB_AM3BETA)*c_surfHeight(Ns) &
            &         + GFB_AM3BETA*c_surfHeight(Bs)
       e_hNormVelMs = (1.5d0 + GFB_AM3BETA)*e_hNormVel(As) - (0.5d0 + 2d0*GFB_AM3BETA)*e_hNormVel(Ns) &
            &         + GFB_AM3BETA*e_hNormVel(Bs)

       call DeepCopy( c_SurfHeightTmp, c_surfHeight(Bs) )
       call TemporalIntegSet_AdvanceShortTStep()
    end do

    !*******************************************************************
    ! Finalize local variables
    !*******************************************************************

    call GeometricField_Final(c_DensVAvg1)
    call GeometricField_Final(c_DensVAvg2)
    call GeometricField_Final(rc_Dens)
    call GeometricField_Final(e_PressForce)
    call GeometricField_Final(c_SurfHeightTmp)
    call GeometricField_Final(e_f)
    call GeometricField_Final(e_totDepthMs)
    call GeometricField_Final(e_totDepthAs)

  contains
    subroutine calc_pressForce(c_surfHeight)
      type(volScalarField), intent(in) :: c_surfHeight
      e_PressForce = (-Grav/refDens)*( &
           & e_c(c_totDepthBasic)*grad(c_DensVAvg2*c_surfHeight) &
           & + grad(0.5d0*c_DensVAvg2*c_surfHeight*c_surfHeight) &
           & + e_c( (c_DensVAvg2 - c_DensVAvg1)*c_surfHeight )*grad(c_totDepthBasic) &
           & )
    end subroutine calc_pressForce
  end subroutine AdvanceBarotStep

  !> @brief 
  !!
  !!
  subroutine FinalizeBarocStep( ze_hNormVelAl, zc_lyrThickAl, zc_TracersAl, &
       & ze_hNormVelNl, ze_hNormVelBl, zc_lyrThickNl, zc_TracersNl, &
       & rc_vNormVel, c_totDepthBasic, zc_lyrThickBasic )
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & DelTime

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: ze_hNormVelAl
    type(volScalarField), intent(inout) :: zc_lyrThickAl
    type(volScalarField), intent(inout) :: zc_TracersAl(:)
    type(surfaceScalarField), intent(in) :: ze_hNormVelNl
    type(surfaceScalarField), intent(in) :: ze_hNormVelBl
    type(volScalarField), intent(in) :: zc_lyrThickNl
    type(volScalarField), intent(in) :: zc_TracersNl(:)
    type(volScalarField), intent(inout) :: rc_vNormVel
    type(volScalarField), intent(in) :: c_totDepthBasic
    type(volScalarField), intent(in) :: zc_lyrThickBasic
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: n
    integer :: nCell
    type(volScalarField) :: zc_RHSTracer
    type(surfaceScalarField) :: ze_massHFluxTmp

    ! 実行文; Executable statement
    !
    
    call GridSet_getLocalMeshInfo(mesh, nCell=nCell)

    ! Stage 4
    !

    !* Update vertical coordinate

    !$omp parallel do
    do i=1, nCell
       zc_lyrThickAl%data%v_(1:nVzLyr,i) = zc_lyrThickBasic%data%v_(1:nVzLyr,i) &
            & * ( 1d0 + c_surfHeightTAvg1%data%v_(1,i)/c_totDepthBasic%data%v_(1,i) )
    end do

    !*
    call adjust_hVelocity(ze_hNormVelAl, e_hMFluxTAvg2, zc_lyrThickAl)

    ze_massHFluxTmp = e_c(zc_lyrThickMl)*( &
         &   (1d0 - LFAM3_EPS)*ze_hNormVelMl &
         & + LFAM3_EPS*((0.5d0-LFAM3_GAM)*ze_hNormVelAl + (0.5d0+2d0*LFAM3_GAM)*ze_hNormVelNl - LFAM3_GAM*ze_hNormVelBl) &
         & )

    do n=1, TracerNum
       call calc_RHSTracer( zc_RHSTracer, &
            & zc_TracersMl(n), ze_massHFluxTmp, rc_vNormVel )

       zc_TracersAl(n) = ( zc_TracersNl(n)*zc_lyrThickNl + DelTime*zc_RHSTracer )/zc_lyrThickAl
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call GeometricField_Final(ze_massHFluxTmp)
    call GeometricField_Final(zc_RHSTracer)

  end subroutine FinalizeBarocStep

  !> @brief 
  !!
  !!
  subroutine adjust_hVelocity(ze_hNormVel, e_hNormVel, zc_lyrThick)
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: ze_hNormVel
    type(surfaceScalarField), intent(in) :: e_hNormVel
    type(volScalarField), intent(in) :: zc_lyrThick
    
    ! 局所変数
    ! Local variables
    !
    type(surfaceScalarField) :: e_correction
    integer :: e
    integer :: nEdge
    real(DP) :: correction

    ! 実行文; Executable statement
    !
    call GridSet_getLocalMeshInfo(mesh, nEdge=nEdge)

    call GeometricField_Init(e_correction, mesh, "e_correction", vLayerNum=1)
    
    e_correction = e_hNormVel - verticalInt(ze_hNormVel, e_c(zc_lyrThick), avgFlag=.true.)

    !$omp parallel do private(correction)
    do e=1, nEdge
       ze_hNormVel%data%v_(1:nVzLyr,e) = ze_hNormVel%data%v_(1:nVzLyr,e) + e_correction%data%v_(1,e)  
    end do

    call GeometricField_Final(e_correction)

  end subroutine adjust_hVelocity


  !> @brief 
  !!
  !!
  subroutine diagnosevNormVel( rc_vNormVel, &
       & ze_hNormVel, zc_lyrThick, dt )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: rc_vNormVel
    type(surfaceScalarField), intent(in) :: ze_hNormVel
    type(volScalarField), intent(in) :: zc_lyrThick
    real(DP), intent(in) :: dt
    
    ! 局所変数
    ! Local variables
    !
    type(volScalarField) :: zc_hdivMassFlux
    integer :: k, i
    integer :: nCell

    ! 実行文; Executable statement
    !
    call GridSet_getLocalMeshInfo(mesh, nCell=nCell)
    call GeometricField_Init(zc_hdivMassFlux, mesh, "zc_hdivMassFlux")
    
    zc_hdivMassFlux = div(zc_lyrThick, ze_hNormVel)
    
    !$omp parallel do 
    do i=1, nCell
       rc_vNormVel%data%v_(nVrLyr,i) = 0d0
       do k=nVrLyr-1, 1, -1
          rc_vNormVel%data%v_(k,i) = rc_vNormVel%data%v_(k+1,i) + zc_hdivMassFlux%data%v_(k,i)
       end do
    end do
    
    call GeometricField_Final(zc_hdivMassFlux)
    
  end subroutine diagnosevNormVel


  !> @brief 
  !!
  !!
  subroutine diagnosePress( zc_Press, &
       & zc_Dens, zc_lyrThick, c_surfPress )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: zc_Press
    type(volScalarField), intent(in) :: zc_Dens
    type(volScalarField), intent(in) :: zc_lyrThick
    type(volScalarField), intent(in) :: c_surfPress

    ! 局所変数
    ! Local variables
    !
    type(volScalarField) :: rc_Dens
    real(DP) :: r_Press(nVrLyr)
    integer :: kr
    integer :: i
    integer :: nCell

    ! 実行文; Executable statement
    !

    call GridSet_getLocalMeshInfo(zc_Press%mesh, nCell=nCell)

    rc_Dens = r_z(zc_Dens)

    !$omp parallel do private(kr, r_Press)
    do i=1, nCell
       r_Press(1) = c_SurfPress%data%v_(1,i)
       do kr=2, nVrLyr
          r_Press(kr) = r_Press(kr-1)  + Grav*zc_Dens%data%v_(kr-1,i)*zc_lyrThick%data%v_(kr-1,i)
       end do

       zc_Press%data%v_(1:nVzLyr,i) = 0.5d0*(r_Press(1:nVrLyr-1) + r_Press(2:nVrLyr)) &
            & + Grav/12d0 *zc_lyrThick%data%v_(1:nVzLyr,i) &
            &   *(rc_Dens%data%v_(1:nVzLyr,i) - rc_Dens%data%v_(2:nVrLyr+1,i)) 
    end do
    
    call GeometricField_Final(rc_Dens)

  end subroutine diagnosePress


end module HydroBouEqSolver_mod

