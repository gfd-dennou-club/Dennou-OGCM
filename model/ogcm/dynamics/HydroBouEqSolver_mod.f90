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
       & plMesh, fvmInfo, &
       & nEdge, nCell, nVertex, &
       & nVzLyr, nVrLyr, vHaloSize

  use EOSDriver_mod, only: &
       & EOSDriver_Eval

  use CGridFieldDataUtil_mod, only: &
       & e_c, ToTangentVel

  use VGridFieldDataUtil_mod, only: &
       & z_r, r_z, verticalInt

  use VariableSet_mod, only: &
       & zc_lyrThickB, zc_lyrThickN, zc_lyrThickA, &
       & ze_hNormVelB, ze_hNormVelN, ze_hNormVelA, &
       & zc_TracersB, zc_TracersN, zc_TracersA, &
       & SaltTracerID, PTempTracerID, TracerNum, &
       & rc_vNormVel, zc_Dens, zc_Press, zc_ZMid, &
       & c_SurfPress, c_totDepthBasic, zc_lyrThickBasic, &
       & refDens

  use HydroBouEqSolverRHS_mod, only: &
       & HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final, &
       & calc_RHShNormVel, calc_RHSTracer, calc_RHSLyrThick

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolver_Init, HydroBouEqSolver_Final
  public :: HydroBouEqSolver_AdvanceTime

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolver_mod' !< Module Name

  type(volScalarField) :: zc_lyrThickMl, zc_lyrThickBMl
  type(surfaceScalarField) :: ze_hNormVelMl
  type(volScalarField), allocatable :: zc_TracersMl(:)
  type(surfaceScalarField) :: e_hMFluxTAvg1
  type(surfaceScalarField) :: e_hMFluxTAvg2
  type(volScalarField) :: c_surfHeightTAvg1
  type(volScalarField) :: c_surfHeightAs, c_surfHeightNs, c_surfHeightBs, c_surfHeightMs
  type(surfaceScalarField) :: e_hNormVelAs, e_hNormVelNs, e_hNormVelBs, e_hNormVelMs

  type(surfaceScalarField) :: e_RHShNormVelMl

  real(DP), parameter :: LFAM3_GAM = 1d0/12d0
  real(DP), parameter :: GFB_AM3BETA = 0.281105d0
  real(DP), parameter :: GFB_Delta = 0.614d0
  real(DP), parameter :: GFB_EPS = 0.013d0
  real(DP), parameter :: GFB_GAM = 0.088d0
  

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Init()

    ! モジュール引用; Use statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: n

    ! 実行文; Executable statements
    !
    
    call GeometricField_Init(zc_lyrThickMl, plMesh, "zc_lyrThickMl", vHaloSize=vHaloSize)
    call GeometricField_Init(zc_lyrThickBMl, plMesh, "zc_lyrThickBMl", vHaloSize=vHaloSize)
    call GeometricField_Init(ze_hNormVelMl, plMesh, "ze_hNormVelMl", vHaloSize=vHaloSize)
    
    allocate( zc_TracersMl(TracerNum) )
    do n=1, TracerNum
       call GeometricField_Init(zc_TracersMl(n), plMesh, trim(zc_TracersN(n)%name)//"Ml")
    end do

    call GeometricField_Init(e_hMFluxTAvg1, plMesh, "e_hMFluxTAvg1", vLayerNum=1)
    call GeometricField_Init(e_hMFluxTAvg2, plMesh, "e_hMFluxTAvg2", vLayerNum=1)
    call GeometricField_Init(c_surfHeightTAvg1, plMesh, "c_surfHeightTAvg1", vLayerNum=1)

    call GeometricField_Init(e_hNormVelAs, plMesh, "e_hNormVelAs", vLayerNum=1)
    call GeometricField_Init(e_hNormVelNs, plMesh, "e_hNormVelNs", vLayerNum=1)
    call GeometricField_Init(e_hNormVelBs, plMesh, "e_hNormVelBs", vLayerNum=1)
    call GeometricField_Init(e_hNormVelMs, plMesh, "e_hNormVelMs", vLayerNum=1)

    call GeometricField_Init(c_surfHeightAs, plMesh, "c_surfHeightAs", vLayerNum=1)
    call GeometricField_Init(c_surfHeightNs, plMesh, "c_surfHeightNs", vLayerNum=1)
    call GeometricField_Init(c_surfHeightBs, plMesh, "c_surfHeightBs", vLayerNum=1)
    call GeometricField_Init(c_surfHeightMs, plMesh, "c_surfHeightMs", vLayerNum=1)

    call GeometricField_Init(e_RHShNormVelMl, plMesh, "e_RHShNormVelMl", vLayerNum=1)

    !
    call HydroBouEqSolverRHS_Init()

  end subroutine HydroBouEqSolver_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Final()

    ! 局所変数
    ! Local variable
    !
    integer :: n

    ! 実行文; Executable statements
    !


    !
    call HydroBouEqSolverRHS_Final()

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

    call GeometricField_Final(e_hNormVelAs)
    call GeometricField_Final(e_hNormVelNs)
    call GeometricField_Final(e_hNormVelBs)
    call GeometricField_Final(e_hNormVelMs)
    call GeometricField_Final(e_RHShNormVelMl)

    call GeometricField_Final(c_surfHeightAs)
    call GeometricField_Final(c_surfHeightNs)
    call GeometricField_Final(c_surfHeightBs)
    call GeometricField_Final(c_surfHeightMs)

  end subroutine HydroBouEqSolver_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTime()
    
    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !

    ! 実行文; Executable statement
    !
    

    call AdvanceBarocStep( e_RHShNormVelMl )
    call AdvanceBarotStep( e_RHShNormVelMl )
    call FinalizeBarocStep()


  end subroutine HydroBouEqSolver_AdvanceTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine AdvanceBarocStep(e_RHShNormVelMl)

    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & DelTime
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: e_RHShNormVelMl


    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: e
    integer :: n

    type(surfaceScalarField) :: ze_RHShNormVel
    type(volScalarField) :: zc_RHSlyrThick
    type(volScalarField) :: zc_RHSTracer
    type(surfaceScalarField) :: ze_massHFlux

    ! 実行文; Executable statement
    !

    !******************************************************************
    ! Stage1: Predicotr update for the 3D momentum
    !******************************************************************

    !
    ! Advance the 3D momentum to n+1/2

    call calc_RHShNormVel( ze_RHShNormVel, &
         & ze_hNormVelN, rc_vNormVel, zc_lyrThickN, zc_ZMid, zc_Press, zc_Dens )

    !$omp parallel do
    do e=1, nEdge
       ze_hNormVelMl%data%v_(1:nVzLyr,e) = &
            &   (5d-1 - 2d0*LFAM3_GAM)*ze_hNormVelA%data%v_(1:nVzLyr,e) &
            & + (5d-1 + 2d0*LFAM3_GAM)*ze_hNormVelN%data%v_(1:nVzLyr,e) &
            & + ( 1d0 - 2d0*LFAM3_GAM)*DelTime*ze_RHShNormVel%data%v_(1:nVzLyr,e) 
    end do
    
    !
    ! Advance  the thickness of vertical layers to n+1/2 with pseudo-compressible algorithm

    ze_massHFlux = ze_hNormVelN*e_c(zc_lyrThickN)
    call calc_RHSLyrThick( zc_RHSlyrThick, &
       & ze_massHFlux, rc_vNormVel )

    !$omp parallel do
    do i=1, nCell
       zc_lyrThickB%data%v_(1:nVzLyr,i) = zc_lyrThickN%data%v_(1:nVzLyr,i) &
            & - (5d-1 - LFAM3_GAM)*DelTime*zc_RHSlyrThick%data%v_(1:nVzLyr,i)

       zc_lyrThickMl%data%v_(1:nVzLyr,i) = zc_lyrThickN%data%v_(1:nVzLyr,i) &
            & + (5d-1 -LFAM3_GAM)*DelTime*zc_RHSlyrThick%data%v_(1:nVzLyr,i)
    end do

    !*
    call adjust_hVelocity(ze_hNormVelMl, e_hMFluxTAvg1, zc_lyrThickMl)


    !*******************************************************************
    ! Stage2: Predicotr update for the tracers
    !*******************************************************************

    do n=1, TracerNum

       call calc_RHSTracer( zc_RHSTracer, &
            & zc_TracersN(n), ze_massHFlux, rc_vNormVel )

      !$omp parallel do
       do i=1, nCell
          zc_TracersMl(n)%data%v_(1:nVzLyr,i) = &
               &   (   (5d-1 - 2d0*LFAM3_GAM)*zc_TracersA(n)%data%v_(1:nVzLyr,i) &
               &     + (5d-1 + 2d0*LFAM3_GAM)*zc_TracersN(n)%data%v_(1:nVzLyr,i) & 
               &   )*zc_lyrThickB%data%v_(1:nVzLyr,i)                               &
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
         & zc_Dens, zc_lyrThickMl)

    call calc_RHShNormVel( ze_RHShNormVel, &
         & ze_hNormVelMl, rc_vNormVel, zc_lyrThickMl, zc_ZMid, zc_Press, zc_Dens )

    !$omp parallel do
    do e=1, nEdge
       ze_hNormVelA%data%v_(1:nVzLyr,e) = ze_hNormVelN%data%v_(1:nVzLyr,e) &
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
  subroutine AdvanceBarotStep( e_RHShNormVelMl )
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & DelTime, SubCycleNum

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: e_RHShNormVelMl

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

    integer :: m
    integer :: mMax
    real(DP) :: DelTimeShort
    real(DP), allocatable :: avgCoef1(:)
    real(DP), allocatable :: avgCoef2(:)
    type(Vector3d) :: geoPos

    ! 実行文; Executable statement
    !

    call GeometricField_Init(c_DensVAvg2, plMesh, "c_DensVAvg2", vLayerNum=1)
    call GeometricField_Init(e_PressForce, plMesh, "e_PressForce", vLayerNum=1)
    call GeometricField_Init(c_SurfHeightTmp, plMesh, "c_SurfHeightTmp", vLayerNum=1)
    call GeometricField_Init(c_f, plMesh, "c_f", vLayerNum=1)
    call GeometricField_Init(e_f, plMesh, "e_f", vLayerNum=1)
 
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

       geoPos = CartToSphPos(plMesh%CellPosList(i))
       c_f%data%v_(1,i) = 2d0*Omega*sin(geoPos%v_(2))
    end do
    
    e_f = e_c(c_f)
    call GeometricField_Final(c_f)

    !
    DelTimeShort = DelTime/dble(SubCycleNum)
    mMax = int(1.25d0*SubCycleNum)
    allocate( avgCoef1(mMax), avgCoef2(mMax) )

    call DeepCopy(e_hNormVelNs, e_hMFluxTAvg1)
    call DeepCOpy(c_surfHeightNs, c_surfHeightTAvg1)

    e_hMFluxTAvg1 = 0d0
    e_hMFluxTAvg2 = 0d0
    c_surfHeightTAvg1 = 0d0

    do m=1, mMax

       e_totDepthMs = e_c(c_totDepthBasic + c_surfHeightMs)
       c_surfHeightAs = c_surfHeightNs - DelTimeShort*div(e_totDepthMs*e_hNormVelMs)
       
       c_SurfHeightTmp = GFB_EPS*c_surfHeightAs &
            & + (1d0 - GFB_Delta - GFB_GAM - GFB_EPS)*c_surfHeightNs &
            & + GFB_GAM*c_surfHeightBs + GFB_EPS*c_SurfHeightTmp

       call calc_pressForce(c_SurfHeightTmp)
       e_totDepthAs = e_c(c_totDepthBasic + c_surfHeightAs)
       e_hNormVelAs = ( &
            &      e_c(c_surfHeightNs)*e_hNormVelNs &
            &    + DelTimeShort*( e_PressForce  + e_RHShNormVelMl - e_totDepthMs*e_f*ToTangentVel(e_hNormVelMs) ) &
            &  )/e_totDepthAs

       e_hMFluxTAvg1 =  e_hMFluxTAvg1 + avgCoef1(m)*e_hNormVelAs*e_totDepthAs
       e_hMFluxTAvg2 =  e_hMFluxTAvg2 + avgCoef2(m)*e_hNormVelMs*e_totDepthMs
       c_surfHeightTAvg1 = c_surfHeightTAvg1 + avgCoef1(m)*c_surfHeightAs


       c_surfHeightMs = (1.5d0 + GFB_AM3BETA)*c_surfHeightAs - (0.5d0 + 2d0*GFB_AM3BETA)*c_surfHeightNs &
            &         + GFB_AM3BETA*c_surfHeightBs
       e_hNormVelMs = (1.5d0 + GFB_AM3BETA)*e_hNormVelAs - (0.5d0 + 2d0*GFB_AM3BETA)*e_hNormVelNs &
            &         + GFB_AM3BETA*e_hNormVelBs

       call DeepCopy(c_SurfHeightTmp, c_surfHeightBs)
       call DeepCopy(c_surfHeightBs, c_surfHeightNs)
       call DeepCopy(c_surfHeightNs, c_surfHeightAs)
       call DeepCopy(e_hNormVelBs, e_hNormVelNs)
       call DeepCopy(e_hNormVelNs, e_hNormVelAs)

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
  subroutine FinalizeBarocStep()
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
         & DelTime

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: n
    type(volScalarField) :: zc_RHSTracer
    type(surfaceScalarField) :: ze_massHFluxMl

    ! 実行文; Executable statement
    !
    
    ! Stage 4
    !

    !* Update vertical coordinate

    !$omp parallel do
    do i=1, nCell
       zc_lyrThickA%data%v_(1:nVzLyr,i) = zc_lyrThickBasic%data%v_(1:nVzLyr,i) &
            & * ( 1d0 + c_surfHeightTAvg1%data%v_(1,i)/c_totDepthBasic%data%v_(1,i) )
    end do

    !*
    call adjust_hVelocity(ze_hNormVelA, e_hMFluxTAvg2, zc_lyrThickA)

    ze_massHFluxMl = ze_hNormVelMl*e_c(zc_lyrThickMl)

    do n=1, TracerNum
       call calc_RHSTracer( zc_RHSTracer, &
            & zc_TracersMl(n), ze_massHFluxMl, rc_vNormVel )

       zc_TracersA(n) = ( zc_TracersN(n)*zc_lyrThickN + DelTime*zc_RHSTracer )/zc_lyrThickA
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call GeometricField_Final(ze_massHFluxMl)
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
    real(DP) :: correction

    ! 実行文; Executable statement
    !

    call GeometricField_Init(e_correction, plMesh, "e_correction", vLayerNum=1)
    
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
    
    ! 実行文; Executable statement
    !
    call GeometricField_Init(zc_hdivMassFlux, plMesh, "zc_hdivMassFlux")
    
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
       & zc_Dens, zc_lyrThick )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: zc_Press
    type(volScalarField), intent(in) :: zc_Dens
    type(volScalarField), intent(in) :: zc_lyrThick
    
    ! 局所変数
    ! Local variables
    !
    type(volScalarField) :: rc_Dens
    real(DP) :: r_Press(nVrLyr)
    integer :: kr
    integer :: i

    ! 実行文; Executable statement
    !


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

