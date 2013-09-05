module GovernEquationSolver_mod

  !
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify

  !
  !

  use VectorSpace_mod

  use SphericalCoord_mod, only: &
       & SphericalTriArea

  use PolyMesh_mod

  use HexTriIcMesh_mod

  use GeometricField_mod

  use fvMeshInfo_mod

  use fvCalculus_mod, only: &
       div, curl, grad

  !
  !
  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod, only: &
       & v_height, s_normalVel

  implicit none
  private

  real(DP), allocatable :: R_iv(:,:)

  public :: GovernEquationSolver_Init, GovernEquationSolver_Final
  public :: Solve_GovernEquation
  
contains
subroutine GovernEquationSolver_Init()
  
  call prepairIntrpWeight()
  
end subroutine GovernEquationSolver_Init

subroutine GovernEquationSolver_Final()

  deallocate(R_iv)

end subroutine GovernEquationSolver_Final

subroutine Solve_GovernEquation()

  type(volScalarField) :: v_mdhdt, v_heightN, v_heightTmp
  type(surfaceScalarField) :: s_mdvdt, s_normalVelN, s_normalVelTmp

  integer :: i , RKStep
  real(DP) :: dt
  real(DP), parameter :: RKCoef(4) = (/ 1d0, 2d0, 2d0, 1d0 /)

  dt = dble(delTime)
  call GeometricField_Init(v_mdhdt, plmesh, "v_Tend_height", "tendency of height", "ms-1")
  call GeometricField_Init(v_heightN, plMesh, "v_heightN", "height of current tstep", "m")
  call GeometricField_Init(v_heightTmp, plMesh, "v_heightTmp", "height of current tstep", "m")
  call GeometricField_Init(s_mdvdt, plmesh, "s_Tend_normalVel", "tendency of velocity", "ms-2")
  call GeometricField_Init(s_normalVelN, plMesh, "s_normalVelN", "normal velocity of current tstep", "ms-1")
  call GeometricField_Init(s_normalVelTmp, plMesh, "s_normalVelTmp", "normal velocity of current tstep", "ms-1")


  v_mdhdt = 0d0; s_mdvdt = 0d0
  call DeepCopy(v_heightN,  v_height)
  call DeepCopy(s_normalVelN, s_normalVel)
  call DeepCopy(v_heightTmp,  v_height)
  call DeepCopy(s_normalVelTmp, s_normalVel)

  do RKStep=1, 4
     v_heightTmp = v_heightN + (-dt/RKCoef(RKStep))*v_mdhdt
     s_normalVelTmp = s_normalVelN + (-dt/RKCoef(RKStep))*s_mdvdt
     call calcTendency( v_mdhdt, s_mdvdt, v_heightTmp, s_normalVelTmp)

     v_height = v_height + (-dt*RKCoef(RKStep)/6d0) * v_mdhdt
     s_normalVel = s_normalVel + (-dt*RKCoef(RKStep)/6d0) * s_mdvdt
  end do

  
  call GeometricField_Final(v_mdhdt)
  call GeometricField_Final(s_mdvdt)
  call GeometricField_Final(v_heightN)
  call GeometricField_Final(v_heightTmp)
  call GeometricField_Final(s_normalVelN)
  call GeometricField_Final(s_normalVelTmp)

end subroutine Solve_GovernEquation

subroutine calcTendency(dmhdt, dmvdt, h, v)
  type(volScalarField), intent(inout) :: dmhdt
  type(surfaceScalarField), intent(inout) :: dmvdt
  type(volScalarField), intent(in) :: h
  type(surfaceScalarField), intent(in) :: v

  dmhdt = div(h, v)
  dmvdt = 0d0!PVFlux() + grad( Grav*h + KEnergy(v) )

end subroutine calcTendency

function KEnergy(s_normalVel) result(KE)
  type(surfaceScalarField), intent(in) :: s_normalVel
  type(volScalarField) :: KE
  
  integer :: cellId, faceId
  integer :: lfaceNum, faceNum
  real(DP), allocatable :: s_localKE(:)

  call GeometricField_Init(KE, plMesh, "KE")
  KE%TempDataFlag = .true.

  faceNum = getFaceListSize(plMesh)
  allocate(s_localKE(faceNum))

  do faceId=1, faceNum
     s_localKE(faceId) = 0.25d0 * l2norm( fvmInfo%s_faceAreaVec.At.faceId ) * (fvmInfo%s_dualMeshFaceArea.At.FaceId) &
          &              * (s_normalVel.At.faceId)**2
  end do

  do cellId=1, getCellListSize(plMesh)
     lfaceNum = plMesh%CellList(cellId)%faceNum
     KE%data%v_(cellId) = sum( s_localKE(fvmInfo%Cell_PointId(1:lfaceNum, cellId)) )/(fvmInfo%v_CellVol.At.cellId)
  end do

end function KEnergy

function PVFlux(s_normalVel) result(flux)
  type(surfaceScalarField), intent(in) :: s_normalVel
  type(surfaceScalarField) :: flux

  type(PointScalarField) :: p_eta ! absolute vorticity
  integer :: faceId, pairCellfaceId, faceNum, k
  real(DP) :: weight(MAX_CELL_FACE_NUM*2)

  call GeometricField_Init(flux, plMesh, "PVFlux")
  flux%TempDataFlag = .true.

  call GeometricField_Init(p_eta, plMesh, "p_eta")


  p_eta = (2d0*Omega) + curl(s_normalVel)
  flux = 0d0

  faceNum = getFaceListSize(plMesh)

  do faceId=1, faceNum
     
!!$     weight = 0d0
!!$     do while( .true. )
!!$        pairCellFaceId = fvmInfo%Face_PairCellFaceId(k, faceId)
!!$     end do
  end do

  call GeometricField_Final(p_eta)

end function PVFlux

subroutine prepairIntrpWeight()

  integer :: cellId
  integer :: lptNum, lptId, m
  real(DP) :: areas(2)
  type(Vector3d) :: ptPos
  integer :: cellIds(2)

  allocate( R_iv(MAX_CELL_FACE_NUM, getCellListSize(plMesh)) )

  do cellId=1, getCellListSize(plMesh)
     lptNum = plMesh%CellList(cellId)%faceNum
     do lptId=1, lptNum
        ptPos = plMesh%PointList(fvmInfo%Cell_PointId(lptId,cellId))
        do m=1, 2
           cellIds(:) = fvmInfo%Face_CellId(1:2, fvmInfo%CellPoint_PairFaceId(m, lptId, cellId))
           areas(m) = 0.5d0*sphericalTriArea( ptPos, plMesh%CellPosList(cellIds(1)), plMesh%CellPosList(cellIds(2)) )
        end do
        R_iv(lptId, cellId) = sum(areas)
     end do
  end do

end subroutine prepairIntrpWeight

end module GovernEquationSolver_mod
