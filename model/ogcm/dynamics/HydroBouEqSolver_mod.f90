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
       & EdgeNum, CellNum, VertexNum, &
       & vzLayerNum, vrLayerNum

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

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Init()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolver_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolver_Final()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolver_Final

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolver_AdvanceTime()
    
    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & ze_hNormVel, rc_vNormVel, zc_lyrThick, c_totDepthBasic, &
         & zc_Press, zc_Zmid, zc_Dens

    use TemporalIntegSet_mod, only: &
         & DelTime, SubCycleNum

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    type(surfaceScalarField) :: ze_hNormVelBaroc
    type(surfaceScalarField) :: e_hNormVelBarot
    type(surfaceScalarField) :: e_hNormVelBarotAvg
    type(volScalarField) :: c_SurfHeight
    type(surfaceScalarField) :: e_G               !< Forcing that eliminates barotoropic velocity components. 
    
    ! 実行文; Executable statement
    !

    !
    call GeometricField_Init(ze_hNormVelBaroc, plMesh, "ze_hNormVelBaroc")
    call GeometricField_Init(e_hNormVelBarot, plMesh, "e_hNormVelBarot", vLayerNum=1)
    call GeometricField_Init(e_hNormVelBarotAvg, plMesh, "e_hNormVelBarotAvg", vLayerNum=1)
    call GeometricField_Init(c_SurfHeight, plMesh, "c_SurfHeight", vLayerNum=1)
    call GeometricField_Init(e_G, plMesh, "e_G", vLayerNum=1)



    call splitToBarotroBaroCliComp( &
         & e_hNormVelBarot, ze_hNormVelBaroc, c_SurfHeight, &     ! (inout)
         & ze_hNormVel, zc_lyrThick, c_totDepthBasic         &    ! (in)            
         & )

    !********************************************************************
    ! Calculate baroclinic velocity (3D) 
    !*********************************************************************

    call advanceBaroclinicCompStep( &
       & ze_hNormVelBaroc, e_G,        &  !(inout) 
       & ze_hNormVelBaroc, e_hNormVelBarot, ze_hNormVel,        &  !(in)
       & rc_vNormVel, zc_lyrThick, c_surfHeight, zc_Zmid, zc_Press, zc_Dens, &  !(in)
       & DelTime )

    !*********************************************************************
    ! Calculate bartropic velocity and surface height perturbation (2D) 
    !**********************************************************************
    !
    call advanceBarotropicCompStep( &
         & e_hNormVelBarot, e_hNormVelBarotAvg, c_SurfHeight, &
         & e_hNormVelBarot, c_SurfHeight, e_G, &
         & DelTime, SubCycleNum )

    !
    call GeometricField_Final(ze_hNormVelBaroc)
    call GeometricField_Final(e_hNormVelBarot)
    call GeometricField_Final(c_SurfHeight)
    call GeometricField_Final(e_G)

  end subroutine HydroBouEqSolver_AdvanceTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


  !> @brief 
  !!
  !!
  subroutine splitToBarotroBaroCliComp( &
       & e_hNormVelBarot, ze_hNormVelBaroc, c_SurfHeight,  &  ! (inout)
       & ze_NormVel, zc_lyrThick, c_totDepthBasic          &  ! (in)            
       & )
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: e_hNormVelBarot
    type(surfaceScalarField), intent(inout) :: ze_hNormVelBaroc
    type(volScalarField), intent(inout) :: c_SurfHeight
    type(surfaceScalarField), intent(in) :: ze_NormVel
    type(volScalarField), intent(in) :: zc_lyrThick
    type(volScalarField), intent(in) :: c_totDepthBasic
    
    ! 局所変数
    ! Local variables
    !
    integer :: e, i, k
    real(DP) :: z_lyrThickTmp(vzLayerNum)

    ! 実行文; Executable statement
    !

    !$omp parallel do private(z_lyrThickTmp)
    do e=1, EdgeNum
       z_lyrThickTmp(1:vzLayerNum) = 0.5d0*sum(zc_lyrThick%data%v_(1:vzLayerNum, fvmInfo%Face_CellId(1:2,e)), 2)
       e_hNormVelBarot%data%v_(1,e) = &
            & sum(z_lyrThickTmp(1:vzLayerNum)*ze_NormVel%data%v_(1:vzLayerNum,e)) &
            & / sum(z_lyrThickTmp)

       ze_hNormVelBaroc%data%v_(1:vzLayerNum,e) = ze_NormVel%data%v_(1:vzLayerNum,e) - e_hNormVelBarot%data%v_(1,e)
    end do

    !$omp parallel do
    do i=1, cellNum
       c_SurfHeight%data%v_(1,i) = sum(zc_lyrThick%data%v_(1:vzLayerNum,i)) - c_totDepthBasic%data%v_(1,i)
    end do
 
 end subroutine splitToBarotroBaroCliComp

  !> @brief 
  !!
  !!
  subroutine advanceBaroclinicCompStep( &
       & ze_hNormVelBarocA, e_G,        &  !(inout) 
       & ze_hNormVelBarocN, e_hNormVelBarotN, ze_hNormVelN, &
       & rc_vNormVelN, zc_lyrThickN, c_surfHeightN, zc_ZmidN, zc_PressN, zc_DensN, &  !(in)
       & dt &
       & )

    ! モジュール引用; Use statement
    !


    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: ze_hNormVelBarocA
    type(surfaceScalarField), intent(inout) :: e_G
    type(surfaceScalarField), intent(in) :: ze_hNormVelBarocN
    type(surfaceScalarField), intent(in) :: e_hNormVelBarotN
    type(surfaceScalarField), intent(in) :: ze_hNormVelN
    type(volScalarField), intent(in) :: rc_vNormVelN
    type(volScalarField), intent(in) :: zc_lyrThickN
    type(volScalarField), intent(in) :: c_surfHeightN
    type(volScalarField), intent(in) :: zc_ZmidN
    type(volScalarField), intent(in) :: zc_PressN
    type(volScalarField), intent(in) :: zc_DensN
    real(DP), intent(in) :: dt
    
    ! 局所変数
    ! Local variables
    !
    integer :: e
    type(Vector3d) :: geoPos
    real(DP) :: f
    real(DP) :: z_lyrThickTmp(vzLayerNum)
    type(surfaceScalarField) :: e_grad_surfHeightN
    type(surfaceScalarField) :: ze_hNormVelTendN

    ! 実行文; Executable statement
    !
    
    call GeometricField_Init(e_grad_surfHeightN, plMesh, "e_grad_surfHeight", vLayerNum=1)
    call GeometricField_Init(ze_hNormVelTendN, plMesh, "ze_hNormVelTendN")

    e_grad_surfHeightN = grad(c_surfHeightN)

    call calchNormVelTendency(ze_hNormVelTendN, &
         & ze_hNormVelN, rc_vNormVelN, zc_lyrThickN, zc_ZmidN, &
         & zc_PressN, zc_DensN &
         & )


    !$omp parallel do private(geoPos, f, z_lyrThickTmp)
    do e=1, EdgeNum

       !
       geoPos = CartToSphPos(fvmInfo%s_faceCenter%data%v_(1,e))
       f = 2d0*Omega*sin(geoPos%v_(2))

       ze_hNormVelBarocA%data%v_(1:vzLayerNum,e) = ze_hNormVelBarocN%data%v_(1:vzLayerNum,e) &
            & + dt*( &
            &   - f*( ze_hNormVelBarocN%data%v_(1:vzLayerNum,e) - e_hNormVelBarotN%data%v_(1,e) ) &
            &   + ze_hNormVelTendN%data%v_(1:vzLayerNum,e)         &
            &   + Grav*e_grad_surfHeightN%data%v_(1:vzLayerNum,e)  & 
            & )

       !
       z_lyrThickTmp(1:vzLayerNum) = 0.5d0*sum(zc_lyrThickN%data%v_(1:vzLayerNum, fvmInfo%Face_CellId(1:2,e)), 2)
       e_G%data%v_(1,e) = sum(z_lyrThickTmp(:)*ze_hNormVelBarocA%data%v_(1:vzLayerNum,e)) / sum(z_lyrThickTmp) / dt

       ze_hNormVelBarocA%data%v_(1:vzLayerNum,e) = ze_hNormVelBarocA%data%v_(1:vzLayerNum,e) -  dt*e_G%data%v_(1:vzLayerNum,e)

    end do
    
    call GeometricField_Final(e_grad_surfHeightN)
    call GeometricField_Final(ze_hNormVelTendN)

  end subroutine advanceBaroclinicCompStep

  !> @brief 
  !!
  !!
  subroutine advanceBarotropicCompStep( &
       & e_hNormVelBarotA, e_hNormVelBarotAvg, c_SurfHeightA, &
       & e_hNormVelBarotN, c_SurfHeightN, e_G, &       
       & dt, J )
    
    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & c_totDepthBasic

    use CGridFieldDataUtil_mod, only: &
         & e_c

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: e_hNormVelBarotA
    type(surfaceScalarField), intent(inout) :: e_hNormVelBarotAvg
    type(volScalarField), intent(inout) :: c_SurfHeightA
    type(surfaceScalarField), intent(in) :: e_hNormVelBarotN
    type(volScalarField), intent(in) :: c_SurfHeightN
    type(surfaceScalarField), intent(in) :: e_G    
    real(DP), intent(in) :: dt
    integer, intent(in) :: J
    
    ! 局所変数
    ! Local variables
    !
    integer :: subCycleCount, e, i
    type(Vector3d) :: geoPos
    real(DP) :: f
    type(surfaceScalarField) :: e_hTangenVel
    type(surfaceScalarField) :: e_grad_surfHeight
    type(surfaceScalarField) :: e_hNormVelBarotTmp
    type(volScalarField) :: c_SurfHeightTmp

    ! 実行文; Executable statement
    !
    
    call GeometricField_Init(e_hNormVelBarotTmp, plMesh, "e_hNormVelBarotTmp", vLayerNum=1)
    call GeometricField_Init(c_SurfHeightTmp, plMesh, "c_SurfHeightTmp", vLayerNum=1)
    call GeometricField_Init(e_hTangenVel, plMesh, "e_hTangenVel", vLayerNum=1)
    call GeometricField_Init(e_grad_surfHeight, plMesh, "e_grad_surfHeight", vLayerNum=1)
    
    call DeepCopy(e_hNormVelBarotTmp, e_hNormVelBarotN)
    call DeepCopy(e_hNormVelBarotAvg, e_hNormVelBarotN)
    call DeepCopy(c_SurfHeightTmp, c_SurfHeightN)
    
    do subCycleCount=1, J
       e_grad_surfHeight = grad(c_surfHeightTmp)

       c_SurfHeightTmp  = c_SurfHeightTmp &
            & - (dt/dble(J)) * div(  c_SurfHeightTmp + c_totDepthBasic, e_hNormVelBarotTmp )

       !$omp parallel do private(geoPos, f)
       do e=1, EdgeNum
          !
          geoPos = CartToSphPos(fvmInfo%s_faceCenter%data%v_(1,e))
          f = 2d0*Omega*sin(geoPos%v_(2))

          e_hNormVelBarotTmp%data%v_(1,e) = e_hNormVelBarotTmp%data%v_(1,e) + dt/dble(J)*( &
               & - f*e_hTangenVel%data%v_(1,e) - Grav*e_grad_surfHeight%data%v_(1,e) &
               & + e_G%data%v_(1,e) &
               & )
       end do

       e_hNormVelBarotAvg = e_hNormVelBarotAvg + e_hNormVelBarotTmp
    end do

    e_hNormVelBarotAvg = e_hNormVelBarotAvg/dble(J)

    call DeepCopy(e_hNormVelBarotA, e_hNormVelBarotTmp)
    call DeepCopy(c_SurfHeightA, c_SurfHeightTmp)

    call GeometricField_Final(e_hTangenVel)
    call GeometricField_Final(e_grad_surfHeight)
    call GeometricField_Final(e_hNormVelBarotTmp)
    call GeometricField_Final(c_SurfHeightTmp)

  end subroutine advanceBarotropicCompStep

  !> @brief 
  !!
  !!
  subroutine calchNormVelTendency( ze_hNormVelTend, &
       & ze_hNormVel, rc_vNormVel, zc_lyrThick, zc_ZMid, zc_Press, zc_Dens )
    
    ! モジュール引用; Use statement
    !
    use VariableSet_mod, only: &
         & refDens

    use CGridFieldDataUtil_mod, only: &
         & e_c, CalcKineticEnergy, CalcPVFlux
    use VGridFieldDataUtil_mod, only: &
         & z_r, delz

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: ze_hNormVelTend
    type(surfaceScalarField), intent(in) :: ze_hNormVel
    type(volScalarField), intent(in) :: rc_vNormVel
    type(volScalarField), intent(in) :: zc_lyrThick
    type(volScalarField), intent(in) :: zc_Zmid
    type(volScalarField), intent(in) :: zc_Press
    type(volScalarField), intent(in) :: zc_Dens

    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    ! ddt(u) = - (1/rho0)*( grad(P) + g*rho*grad(zmid) )
    !         - grad(K)
    !         - PVFlux
    !         - w*ddz(u)
    !
    ze_hNormVelTend = &
         & (-1d0/refDens)*( grad(zc_Press)  + Grav*e_c(zc_Dens)*grad(zc_Zmid) ) &
         & - grad( CalcKineticEnergy(ze_hNormVel) ) &
         & - CalcPVFlux(zc_lyrThick, ze_hNormVel)   &
         & - z_r( e_c(rc_vNormVel) * delz(ze_hNormVel, zc_lyrThick) )

  end subroutine calchNormVelTendency

  !> @brief 
  !!
  !!
  subroutine calchLyrThickTendency( zc_lyrThickTend &
       & )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: zc_lyrThickTend 
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
  end subroutine calchLyrThickTendency

end module HydroBouEqSolver_mod

