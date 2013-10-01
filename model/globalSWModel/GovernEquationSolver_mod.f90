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
       & SphericalTriArea, &
       & CartToSphPos, radToDegUnit

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
       & v_h, s_normalVel, v_hb

  implicit none
  private

  real(DP), allocatable :: R_iv(:,:)
  real(DP), allocatable :: Wt_ff(:,:)
  type(volScalarField) :: v_mdhdt, v_hN, v_hTmp
  type(surfaceScalarField) :: s_mdvdt, s_normalVelN, s_normalVelTmp

  public :: GovernEquationSolver_Init, GovernEquationSolver_Final
  public :: Solve_GovernEquation
  
contains
subroutine GovernEquationSolver_Init()
  
  call prepairIntrpWeight()

  call GeometricField_Init(v_mdhdt, plmesh, "v_Tend_height", "tendency of height", "ms-1")
  call GeometricField_Init(v_hN, plMesh, "v_hN", "height of current tstep", "m")
  call GeometricField_Init(v_hTmp, plMesh, "v_hTmp", "height of current tstep", "m")
  call GeometricField_Init(s_mdvdt, plmesh, "s_Tend_normalVel", "tendency of velocity", "ms-2")
  call GeometricField_Init(s_normalVelN, plMesh, "s_normalVelN", "normal velocity of current tstep", "ms-1")
  call GeometricField_Init(s_normalVelTmp, plMesh, "s_normalVelTmp", "normal velocity of current tstep", "ms-1")
  
end subroutine GovernEquationSolver_Init

subroutine GovernEquationSolver_Final()

  deallocate(R_iv)
  deallocate(Wt_ff)

  
  call GeometricField_Final(v_mdhdt)
  call GeometricField_Final(s_mdvdt)
  call GeometricField_Final(v_hN)
  call GeometricField_Final(v_hTmp)
  call GeometricField_Final(s_normalVelN)
  call GeometricField_Final(s_normalVelTmp)

end subroutine GovernEquationSolver_Final

subroutine Solve_GovernEquation()


  integer :: i , RKStep
  real(DP) :: dt
  real(DP), parameter :: RKCoef(4) = (/ 1d0, 2d0, 2d0, 1d0 /)

  dt = dble(delTime)


  v_mdhdt = 0d0; s_mdvdt = 0d0
  call DeepCopy(v_hN,  v_h)
  call DeepCopy(s_normalVelN, s_normalVel)
  call DeepCopy(v_hTmp,  v_h)
  call DeepCopy(s_normalVelTmp, s_normalVel)

  do RKStep=1, 4
     v_hTmp = v_hN + (-dt/RKCoef(RKStep))*v_mdhdt
     s_normalVelTmp = s_normalVelN + (-dt/RKCoef(RKStep))*s_mdvdt
     call calcTendency( v_mdhdt, s_mdvdt, v_hTmp, s_normalVelTmp)

     v_h = v_h + (-dt*RKCoef(RKStep)/6d0) * v_mdhdt
     s_normalVel = s_normalVel + (-dt*RKCoef(RKStep)/6d0) * s_mdvdt
  end do


end subroutine Solve_GovernEquation

subroutine calcTendency(dmhdt, dmvdt, h, v)
  type(volScalarField), intent(inout) :: dmhdt
  type(surfaceScalarField), intent(inout) :: dmvdt
  type(volScalarField), intent(in) :: h
  type(surfaceScalarField), intent(in) :: v


  dmhdt = div(h, v)
  dmvdt = grad( Grav*(h + v_hb) + KEnergy(v) ) + PVFlux(h, v) 

end subroutine calcTendency



!
!
!
function KEnergy(s_normalVel) result(KE)
  type(surfaceScalarField), intent(in) :: s_normalVel
  type(volScalarField) :: KE
  
  integer :: cellId, faceId
  integer :: lfaceNum, faceNum, cellNum
  real(DP) :: s_localKE(getVLayerSize(plMesh), getFaceListSize(plMesh))

  call GeometricField_Init(KE, plMesh, "KE")
  KE%TempDataFlag = .true.

  faceNum = getFaceListSize(plMesh)
  cellNum = getCellListSize(plMesh)

  !$omp parallel do
  do faceId=1, faceNum
     s_localKE(:,faceId) = 0.25d0 * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId) ) &
          & * fvmInfo%s_dualMeshFaceArea%data%v_(1,FaceId) &
          & * s_normalVel%data%v_(:,FaceId)**2
  end do

  !$omp parallel do private(lfaceNum)
  do cellId=1, cellNum
     lfaceNum = plMesh%CellList(cellId)%faceNum
     KE%data%v_(:,cellId) = sum( s_localKE(:,fvmInfo%Cell_FaceId(1:lfaceNum, cellId)), 2)/ fvmInfo%v_CellVol%data%v_(:,cellId)
  end do

end function KEnergy

function PVFlux(v_h, s_normalVel) result(s_PVflux)
  type(volScalarField), intent(in) :: v_h
  type(surfaceScalarField), intent(in) :: s_normalVel
  type(surfaceScalarField) :: s_PVFlux

  type(PointScalarField) :: p_planetVor
  type(PointScalarField) :: p_height, p_PV
 
  integer :: faceId, faceId2, faceNum, k, i
  real(DP) :: weight, tmp(1), s2_MomFlux(1), s_PV(1), s2_PV(1)
  integer :: ptNum, ptId, cellIds(3), cellId
  type(Vector3d) :: geoPos


  call GeometricField_Init(p_planetVor, plMesh, "p_planetVor")
  call GeometricField_Init(p_PV, plMesh, "p_eta")
  call GeometricField_Init(p_height, plMesh, "p_height")

  ptNum = getPointListSize(plMesh)

  !$omp parallel do private(geoPos, cellIds)
  do ptId=1, ptNum
     geoPos = CartToSphPos( plMesh%PointPosList(ptId) )
     p_planetVor%data%v_(1,ptId) = (2d0*Omega)*sin(geoPos%v_(2))

     cellIds(:) = fvmInfo%Point_CellId(1:3, ptId)
    
     p_height%data%v_(:,ptId) = sum( &
          & (spread(R_iv(1:3,ptId),1,1)) * fvmInfo%v_cellVol%data%v_(:,cellIds(:)) * v_h%data%v_(:,cellIds(:)) &
          & , 2) / fvmInfo%p_dualMeshCellVol%data%v_(:, ptId)

  end do

  p_PV =  ( p_planetVor + curl(s_normalVel) )/ p_height

  !
  !

  call GeometricField_Init(s_PVFlux, plMesh, "s_PVFlux")
  s_PVFlux%TempDataFlag = .true.

  faceNum = getFaceListSize(plMesh)
  
  !$omp parallel do private(s_PV, s2_PV, s2_MomFlux, tmp, k, faceId2) 
  do faceId=1, faceNum
     s_PV(:) = 0.5d0*sum( p_PV%data%v_(:,fvmInfo%Face_PointId(1:2, faceId)), 2)
     tmp = 0d0

     do k=1, size(fvmInfo%Face_PairCellFaceId,1)
        faceId2 = fvmInfo%Face_PairCellFaceId(k, faceId)
        if( faceId2 == -1) exit;

        s2_MomFlux(:) = 0.5d0*sum( v_h%data%v_(:, fvmInfo%Face_CellId(1:2, faceId2)), 2) * s_normalVel%data%v_(:,faceId2)
        s2_PV(:) = 0.5d0*sum( p_PV%data%v_(:, fvmInfo%Face_PointId(1:2, faceId2)), 2)

        tmp(:) =    tmp(:) &
             & + wt_ff(k,faceId) * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId2) ) * s2_MomFlux * 0.5d0*(s_PV + s2_PV) 


     end do ! do while

     s_PVFlux%data%v_(:,faceId) = tmp(:) / fvmInfo%s_dualMeshFaceArea%data%v_(:,faceId)
    
  end do

  !
  !
  call GeometricField_Final(p_planetVor)
  call GeometricField_Final(p_PV)
  call GeometricField_Final(p_height)

end function PVFlux

subroutine prepairIntrpWeight()

  integer :: cellNum, ptNum, faceNum, lfaceNum
  integer :: cellId, cellIds(3), ptId, faceId, lfaceId
  integer :: m, n, k
  real(DP) :: areas(3)
  type(Vector3d) :: ptPos
  real(DP), allocatable :: R_vi(:,:)
  integer :: startPtId, startlPtId, startlfaceId, t_fv2
  real(DP) :: lR_vi(MAX_CELL_FACE_NUM), ln_fv(MAX_CELL_FACE_NUM)

  cellNum = getCellListSize(plMesh)
  ptNum = getPointListSize(plMesh)
  faceNum = getFaceListSize(plMesh)
  allocate( R_iv(3, ptNum) )
  allocate( Wt_ff(MAX_CELL_FACE_NUM*2, faceNum) )
  allocate(R_vi(MAX_CELL_FACE_NUM, cellNum))

  R_vi = 0d0
  do ptId=1, ptNum
     ptPos = plMesh%PointPosList(ptId)
     cellIds = cshift(fvmInfo%Point_CellId(1:3,ptId), -1)

     do m=1, 3
        areas(m) = sphericalTriArea(ptPos, plMesh%CellPosList(cellIds(m)), plMesh%CellPosList(cellIds(mod(m,3)+1)) )
     end do

     do m=1, 3
        cellId = fvmInfo%Point_CellId(m, ptId) 
        R_iv(m, ptId) = 0.5d0*( areas(m) + areas(mod(m,3)+1) )/ At(fvmInfo%v_CellVol, cellId) 

        do n=1, MAX_CELL_FACE_NUM
           if( fvmInfo%Cell_PointId(n,cellId) == ptId) then
              R_vi(n, cellId) = R_iv(m, ptId); exit
           end if
        end do

     end do
  end do

  Wt_ff = 0d0
  do faceId=1, faceNum   

     k=2
     do m=1, 2
        cellId = fvmInfo%Face_CellId(m, faceId)
        startPtId = fvmInfo%Face_PointId(mod(m,2)+1, faceId)
        lfaceNum = plMesh%CellList(cellId)%faceNum
        do n=1, lfaceNum
           if( fvmInfo%Cell_PointId(n,cellId) == startPtId ) then
              startlPtId = n; exit
           end if
        end do

        do n=1, lfaceNum
           if( fvmInfo%Cell_FaceId(n,cellId) == faceId ) then
              startlFaceId = n; exit
           end if

        end do


        !
        lR_vi(1:lfaceNum) = cshift(R_vi(1:lfaceNum,cellId), startlPtId-1)
        ln_fv(1:lfaceNum) = cshift(fvmInfo%n_fv(1:lfaceNum,cellId), startlFaceId-1 )
        do lfaceId=1, 3
           if( fvmInfo%Point_FaceId(lfaceId,startPtId) == faceId) then
              t_fv2 = fvmInfo%t_fv(lfaceId, startPtId); exit
           end if
        end do

        do n=2, lfaceNum
           Wt_ff(k, faceId) =  dble(ln_fv(n)*t_fv2) * (sum(lR_vi(1:n-1)) - 0.5d0)
           k = k + 1
        end do

     end do

! write(*,'(a, i6, a, 13f12.5)') "wt_ff", faceId, ":", wt_ff(:,faceId), sum(wt_ff(:,faceId))
  end do

end subroutine prepairIntrpWeight


end module GovernEquationSolver_mod
