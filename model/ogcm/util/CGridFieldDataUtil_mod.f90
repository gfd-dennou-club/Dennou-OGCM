!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module CGridFieldDataUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod, only: &
       & Vector3d, l2Norm

  use GridSet_mod, only: &
       & plMesh, fvmInfo, &
       & nCell, nEdge, nVertex, &
       & nVzLyr, nVrLyr, vHaloSize

  use fvCalculus_mod, only: &
       div, curl, grad, t_fv, n_fv

  use GeometricField_mod, only: &
       & volScalarField, surfaceScalarField, pointScalarField, &
       & GeometricField_Init, GeometricField_Final, Release, &
       & At, assignment(=), operator(+), operator(/)

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  interface e_c
     module procedure ze_zc
  end interface e_c

  public :: CGridFieldDataUtil_Init, CGridFieldDataUtil_Final
  public :: e_c
  public :: CalcKineticEnergy, CalcPVFlux, ToTangentVel

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'CGridFieldDataUtil_mod' !< Module Name
  integer, parameter :: MAX_CELL_FACE_NUM = 6
  real(DP), allocatable :: R_iv(:,:)
  real(DP), allocatable :: Wt_ff(:,:)

contains

  !>
  !!
  !!
  subroutine CGridFieldDataUtil_Init()

    ! 実行文; Executable statements
    !

    
    call prepairIntrpWeight()

  end subroutine CGridFieldDataUtil_Init

  !>
  !!
  !!
  subroutine CGridFieldDataUtil_Final()

    ! 実行文; Executable statements
    !

    deallocate(R_iv)
    deallocate(Wt_ff)

  end subroutine CGridFieldDataUtil_Final


  !> @brief 
  !!
  !! @return 
  !!
  function ze_zc(zc_field) result(ze_field)
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(in) :: zc_field
    type(surfaceScalarField) :: ze_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: e, lyrNum
    
    ! 実行文; Executable statement
    !
    
    lyrNum = zc_field%vLayerNum
    call GeometricField_Init(ze_field, plMesh, "ze_field", vLayerNum=lyrNum)
    ze_field%TempDataFlag = .true.

    !$omp parallel do
    do e=1, nEdge
       ze_field%data%v_(1:lyrNum,e) = 0.5d0*sum(zc_field%data%v_(1:lyrNum,fvmInfo%Face_CellId(1:2,e)), 2)
    end do

    if(zc_field%TempDataFlag) call Release(zc_field)
    
  end function ze_zc



  !
  !
  !
  function calcKineticEnergy(ze_normalVel) result(zc_KE)
    type(surfaceScalarField), intent(in) :: ze_normalVel
    type(volScalarField) :: zc_KE

    integer :: cellId, faceId
    integer :: lfaceNum
    real(DP) :: ze_localKE(nVzLyr, nEdge)

    call GeometricField_Init(zc_KE, plMesh, "KE")
    zc_KE%TempDataFlag = .true.

    !$omp parallel do
    do faceId=1, nEdge
       ze_localKE(1:nVzLyr,faceId) = 0.25d0 * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId) ) &
            & * fvmInfo%s_dualMeshFaceArea%data%v_(1,FaceId) &
            & * ze_normalVel%data%v_(1:nVzLyr,FaceId)**2
    end do

    !$omp parallel do private(lfaceNum)
    do cellId=1, nCell
       lfaceNum = plMesh%CellList(cellId)%faceNum
       zc_KE%data%v_(1:nVzLyr,cellId) = &
            & sum( ze_localKE(1:nVzLyr,fvmInfo%Cell_FaceId(1:lfaceNum, cellId)), 2) &
            & / fvmInfo%v_CellVol%data%v_(1:nVzLyr,cellId)
    end do

  end function CalcKineticEnergy

  function CalcPVFlux(zc_h, ze_normalVel) result(ze_PVflux)

    use SphericalCoord_mod, only: &
         & CartToSphPos

    use Constants_mod, only: &
         & Omega


    type(volScalarField), intent(in) :: zc_h
    type(surfaceScalarField), intent(in) :: ze_normalVel
    type(surfaceScalarField) :: ze_PVFlux

    type(PointScalarField) :: zv_planetVor
    type(PointScalarField) :: zv_height
    type(PointScalarField) :: zv_PV

    integer :: faceId, faceId2, k, i
    real(DP) :: weight, tmp(nVzLyr), s2_MomFlux(nVzLyr), s_PV(nVzLyr), s2_PV(nVzLyr)
    integer :: ptId, cellIds(3), cellId
    type(Vector3d) :: geoPos


    call GeometricField_Init(zv_planetVor, plMesh, "p_planetVor")
    call GeometricField_Init(zv_PV, plMesh, "p_eta")
    call GeometricField_Init(zv_height, plMesh, "p_height")

    !$omp parallel do private(geoPos, cellIds)
    do ptId=1, nVertex
       geoPos = CartToSphPos( plMesh%PointPosList(ptId) )
       zv_planetVor%data%v_(1:nVzLyr,ptId) = (2d0*Omega)*sin(geoPos%v_(2))

       cellIds(:) = fvmInfo%Point_CellId(1:3, ptId)

       forAll(k=1:nVzLyr) &
            & zv_height%data%v_(k,ptId) = &
            & sum( R_iv(1:3,ptId) * fvmInfo%v_cellVol%data%v_(k,cellIds(:)) * zc_h%data%v_(k,cellIds(:)) ) &
            & / fvmInfo%p_dualMeshCellVol%data%v_(k, ptId)

    end do

    zv_PV = ( zv_planetVor + curl(ze_normalVel) )/ zv_height

    !
    !

    call GeometricField_Init(ze_PVFlux, plMesh, "s_PVFlux")
    ze_PVFlux%TempDataFlag = .true.

    !$omp parallel do private(s_PV, s2_PV, s2_MomFlux, tmp, k, faceId2) 
    do faceId=1, nEdge
       s_PV(:) = 0.5d0*sum( zv_PV%data%v_(1:nVzLyr,fvmInfo%Face_PointId(1:2, faceId)), 2)
       tmp = 0d0

       do k=1, size(fvmInfo%Face_PairCellFaceId,1)
          faceId2 = fvmInfo%Face_PairCellFaceId(k, faceId)
          if( faceId2 == -1) exit;

          s2_MomFlux(:) = 0.5d0*sum( zc_h%data%v_(1:nVzLyr, fvmInfo%Face_CellId(1:2, faceId2)), 2) &
               &           * ze_normalVel%data%v_(1:nVzLyr, faceId2)
          s2_PV(:) = 0.5d0*sum( zv_PV%data%v_(1:nVzLyr, fvmInfo%Face_PointId(1:2, faceId2)), 2)

          tmp(:) =    tmp(:) &
               & + wt_ff(k,faceId) * l2norm( fvmInfo%s_faceAreaVec%data%v_(1,faceId2) ) * s2_MomFlux * 0.5d0*(s_PV + s2_PV) 


       end do ! do while

       ze_PVFlux%data%v_(1:nVzLyr,faceId) = tmp(:) / fvmInfo%s_dualMeshFaceArea%data%v_(1,faceId)

    end do

    !
    !
    call GeometricField_Final(zv_planetVor)
    call GeometricField_Final(zv_height)
    call GeometricField_Final(zv_PV)

  end function CalcPVFlux

  function ToTangentVel(ze_normalVel) result(ze_tangentVel)

    use SphericalCoord_mod, only: &
         & CartToSphPos

    use Constants_mod, only: &
         & Omega

    type(surfaceScalarField), intent(in) :: ze_normalVel
    type(surfaceScalarField) :: ze_tangentVel

    integer :: faceId, faceId2, k, i
    real(DP) :: weight, tmp(nVzLyr)
    integer :: ptId, cellIds(3), cellId


    !
    !

    call GeometricField_Init(ze_tangentVel, plMesh, "ze_tangentVel")
    ze_tangentVel%TempDataFlag = .true.

    !$omp parallel do private(tmp, k, faceId2) 
    do faceId=1, nEdge
       tmp = 0d0

       do k=1, size(fvmInfo%Face_PairCellFaceId,1)
          faceId2 = fvmInfo%Face_PairCellFaceId(k, faceId)
          if( faceId2 == -1) exit;

          tmp(:) =    tmp(:)  + &
               & wt_ff(k,faceId) * l2norm(fvmInfo%s_faceAreaVec%data%v_(1,faceId2)) * ze_normalVel%data%v_(1:nVzLyr, faceId2)

       end do ! do while

       ze_tangentVel%data%v_(1:nVzLyr,faceId) = tmp(:) / fvmInfo%s_dualMeshFaceArea%data%v_(1,faceId)
    end do

    !
    !

  end function ToTangentVel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepairIntrpWeight()
    
    ! モジュール引用; Use statements
    !
    use SphericalCoord_mod, only: &
         & SphericalTriArea

    ! 局所変数
    ! Local variables
    !

    integer :: lfaceNum
    integer :: cellId, cellIds(3), ptId, faceId, lfaceId
    integer :: m, n, k
    real(DP) :: areas(3)
    type(Vector3d) :: ptPos
    real(DP) :: R_vi(MAX_CELL_FACE_NUM, nCell)
    integer :: startPtId, startlPtId, startlfaceId, t_fv2
    real(DP) :: lR_vi(MAX_CELL_FACE_NUM), ln_fv(MAX_CELL_FACE_NUM)

    ! 実行文; Executable statements
    !

    call MessageNotify('M', module_name, "prepare some data about weight required interpolation..")

    allocate( R_iv(3, nVertex) )
    allocate( Wt_ff(MAX_CELL_FACE_NUM*2, nEdge) )

    R_vi = 0d0
    do ptId=1, nVertex
       ptPos = plMesh%PointPosList(ptId)
       cellIds = cshift(fvmInfo%Point_CellId(1:3,ptId), -1)

       do m=1, 3
          areas(m) = sphericalTriArea(ptPos, plMesh%CellPosList(cellIds(m)), plMesh%CellPosList(cellIds(mod(m,3)+1)) )
       end do

       do m=1, 3
          cellId = fvmInfo%Point_CellId(m, ptId) 
          R_iv(m, ptId) = 0.5d0*( areas(m) + areas(mod(m,3)+1) )/ At(fvmInfo%v_CellVol, 1, cellId) 

          do n=1, MAX_CELL_FACE_NUM
             if( fvmInfo%Cell_PointId(n,cellId) == ptId) then
                R_vi(n, cellId) = R_iv(m, ptId); exit
             end if
          end do

       end do
    end do

    Wt_ff = 0d0
    do faceId=1, nEdge   

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
          ln_fv(1:lfaceNum) = cshift(n_fv(1:lfaceNum,cellId), startlFaceId-1 )
          do lfaceId=1, 3
             if( fvmInfo%Point_FaceId(lfaceId,startPtId) == faceId) then
                t_fv2 = t_fv(lfaceId, startPtId); exit
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


end module CGridFieldDataUtil_mod

