!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolverSelfStart_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING 

  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod, only: &
       & Vector3d

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

  use HydroBouEqSolver_mod, only: &
       & diagnosePress, diagnosevNormVel

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
  public :: HydroBouEqSolverSelfStart_Init, HydroBouEqSolverSelfStart_Final
  public :: HydroBouEqSolverSelfStart_AdvanceTStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverSelfStart_mod' !< Module Name
  type(PolyMesh), pointer :: mesh

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverSelfStart_Init(plMesh)

    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !
    type(PolyMesh), intent(in), target :: plMesh

    ! 実行文; Executable statements
    !

    mesh => plMesh

  end subroutine HydroBouEqSolverSelfStart_Init


  !>
  !!
  !!
  subroutine HydroBouEqSolverSelfStart_AdvanceTStep(var)
    
    ! モジュール引用; Use statements
    !
    use TemporalIntegSet_mod, only: &
       & DelTime, SubCycleNum, &
       & Al, Nl, Bl, &
       & As, Ns, Bs

    ! 宣言文; Declaration statement
    !
    type(VariableSet), intent(inout) :: var
    
    ! 局所変数
    ! Local variables
    !
    
    real(DP) :: dt
    integer :: n

    type(surfaceScalarField) :: ze_hNormVel, ze_hNormVelTmp, ze_hNormVel0, ze_RHShNormVel
    type(volScalarField) :: zc_lyrThick, zc_lyrThickTmp, zc_lyrThick0, zc_RHSlyrThick
    type(surfaceScalarField) :: ze_massHFluxTmp
    type(volScalarField) :: zc_TracersTmp(TracerNum), zc_Tracers0(TracerNum), zc_Tracers(TracerNum), &
         & zc_RHSTracers(TracerNum)

    integer :: nIntrnalLoop
    integer :: internalStep

    ! 実行文; Executable statements
    !
    
    nIntrnalLoop = SubCycleNum*2
    dt = DelTime/dble(nIntrnalLoop)

    call GeometricField_Init(ze_hNormVel0, mesh, "ze_hNormVel0", vHaloSize=vHaloSize)
    call GeometricField_Init(ze_RHShNormVel, mesh, "ze_RHShNormVel", vHaloSize=vHaloSize)
    ze_hNormVel = var%ze_hNormVel(Nl)

    call GeometricField_Init(zc_lyrThick0, mesh, "zc_lyrThick0", vHaloSize=vHaloSize)
    call GeometricField_Init(zc_RHSlyrThick, mesh, "zc_RHSlyrThick", vHaloSize=vHaloSize)
    zc_lyrThick = var%zc_lyrThick(Nl)

    do n=1, TracerNum
       call GeometricField_Init(zc_Tracers0(n), mesh, "zc_Tracers0", vHaloSize=vHaloSize)
       call GeometricField_Init(zc_RHSTracers(n), mesh, "zc_RHSTracers", vHaloSize=vHaloSize)
       zc_Tracers(n) = var%zc_Tracers(n,Nl)
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do internalStep=1, nIntrnalLoop
       call MessageNotify("M", module_name, "Self starting.. %d/%d", i=(/ internalStep, nIntrnalLoop /))
       call advanceIntrnalStep()
    end do

    var%ze_hNormVel(Al) = ze_hNormVel
    var%zc_lyrThick(Al) = zc_lyrThick
    do n=1, TracerNum
       var%zc_Tracers(n,Al) = zc_Tracers(n)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call GeometricField_Final(ze_RHShNormVel)
    call GeometricField_Final(ze_hNormVel)
    call GeometricField_Final(ze_hNormVel0)
    call GeometricField_Final(ze_hNormVelTmp)

    call GeometricField_Final(zc_RHSlyrThick)
    call GeometricField_Final(zc_lyrThick)
    call GeometricField_Final(zc_lyrThick0)
    call GeometricField_Final(zc_lyrThickTmp)

    call GeometricField_Final(ze_massHFluxTmp)

    do n=1, TracerNum
       call GeometricField_Final(zc_Tracers0(n))
       call GeometricField_Final(zc_TracersTmp(n))
       call GeometricField_Final(zc_Tracers(n))
       call GeometricField_Final(zc_RHSTracers(n))
    end do

  contains 
    subroutine advanceIntrnalStep()
      integer :: RKStep
      real(DP), parameter :: RKCoef(4) = (/ 1d0, 2d0, 2d0, 1d0 /)
      integer :: tmpBl


      ze_RHShNormVel = 0d0
      call DeepCopy(ze_hNormVel0, ze_hNormVel)

      zc_RHSlyrThick = 0d0
      do n=1, TracerNum
         zc_RHSTracers(n) = 0d0
      end do

      do RKStep=1, size(RKCoef)
         
           ze_hNormVelTmp = ze_hNormVel0 + (dt/RKCoef(RKStep))*ze_RHShNormVel
           zc_lyrThickTmp = zc_lyrThick0 + (dt/RKCoef(RKStep))*zc_RHSlyrThick
           
           do n=1, TracerNum
              zc_TracersTmp(n) = zc_Tracers0(n) + (-dt/RKCoef(RKStep))*zc_RHSTracers(n)
           end do

           call EOSDriver_Eval( rho=var%zc_Dens, &
                & theta=zc_TracersTmp(PTempTracerID), S=zc_TracersTmp(SaltTracerID), p=var%zc_Press)

           call diagnosePress( var%zc_Press, &
                & var%zc_Dens, zc_lyrThickTmp, var%c_surfPress)

           call calc_RHShNormVel( ze_RHShNormVel, &
                & ze_hNormVelTmp, var%rc_vNormVel, zc_lyrThickTmp, var%zc_ZMid, var%zc_Press, var%zc_Dens )

           ze_massHFluxTmp = ze_hNormVelTmp*e_c(zc_lyrThickTmp)
           call calc_RHSLyrThick( zc_RHSlyrThick, &
                & ze_massHFluxTmp, var%rc_vNormVel )

           do n=1, TracerNum
              call calc_RHSTracer( zc_RHSTracers(n), &
                   & zc_TracersTmp(n), ze_massHFluxTmp, var%rc_vNormVel )
           end do

           var%ze_hNormVel(Al) = var%ze_hNormVel(Al) + (dt/RKCoef(RKStep))*ze_RHShNormVel
           var%zc_lyrThick(Al) = var%zc_lyrThick(Al) + (dt/RKCoef(RKStep))*zc_RHSlyrThick
           do n=1, TracerNum
              var%zc_Tracers(n,Al) = var%zc_Tracers(n,Al) + (dt/RKCoef(RKStep))*zc_RHSTracers(n)
           end do

        end do

        tmpBl = Bl
        Bl = Nl
        Nl = Al
        Al = tmpBl

      end subroutine advanceIntrnalStep

  end subroutine HydroBouEqSolverSelfStart_AdvanceTStep

  !>
  !!
  !!
  subroutine HydroBouEqSolverSelfStart_Final()

    ! 実行文; Executable statements
    !

    mesh => null()

  end subroutine HydroBouEqSolverSelfStart_Final

end module HydroBouEqSolverSelfStart_mod
