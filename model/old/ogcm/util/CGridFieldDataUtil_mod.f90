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

  use PolyMesh_mod, only: &
       & PolyMesh

  use fvMeshInfo_mod, only: &
       & fvMeshInfo

  use GridSet_mod, only: &
       & GridSet_getLocalMeshInfo, &
       & nVzLyr, nVrLyr, vHaloSize

  use fvCalculus_mod, only: &
       div, curl, grad

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

  type, public :: CGridFieldDataUtil
     real(DP), allocatable :: R_iv(:,:)
     real(DP), allocatable :: Wt_ff(:,:)
     type(PolyMesh), pointer :: mesh => null()
     type(fvMeshInfo), pointer :: fvmInfo => null() 
  end type CGridFieldDataUtil

  public :: CGridFieldDataUtil_Init, CGridFieldDataUtil_Final
  public :: CGridFieldDataUtil_Set
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

  type(CGridFieldDataUtil), pointer :: defaultUtilObj => null()

contains

  !>
  !!
  !!
  subroutine CGridFieldDataUtil_Init(this, fvmInfo)

    ! 宣言文; Declaration statement
    !
    type(CGridFieldDataUtil), intent(inout) :: this
    type(fvMeshInfo), intent(in), target :: fvmInfo

    ! 実行文; Executable statements
    !

    this%fvmInfo => fvmInfo
    this%mesh => fvmInfo%mesh
    call prepairIntrpWeight(this)

  end subroutine CGridFieldDataUtil_Init

  !>
  !!
  !!
  subroutine CGridFieldDataUtil_Final(this)

    type(CGridFieldDataUtil), intent(inout) :: this

    ! 実行文; Executable statements
    !

    deallocate(this%R_iv)
    deallocate(this%Wt_ff)

  end subroutine CGridFieldDataUtil_Final


  !> @brief 
  !!
  !!
  subroutine CGridFieldDataUtil_Set(setUtilObj)
    
    ! 宣言文; Declaration statement
    !
    type(CGridFieldDataUtil), intent(in), target :: setUtilObj
    
    ! 実行文; Executable statement
    !
    
    defaultUtilObj => setUtilObj

#ifdef DEBUG
    call MessageNotify("M", module_name, "Set default object used in the subroutines in this module.")
#endif
    
  end subroutine CGridFieldDataUtil_Set

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
    integer :: e, lyrNum, nEdge
    type(PolyMesh), pointer :: mesh
    type(fvMeshInfo), pointer :: fvm

    ! 実行文; Executable statement
    !

#ifdef DEBUG
    if( .not. associated(defaultUtilObj) ) &
         & call MessageNotify("E", module_name, "Set the pointer of default object to actual object.")
#endif

    lyrNum = zc_field%vLayerNum
    mesh => defaultUtilObj%mesh
    fvm => defaultUtilObj%fvmInfo

    call GridSet_getLocalMeshInfo(mesh, nEdge=nEdge)
    call GeometricField_Init(ze_field, mesh, "ze_field", vLayerNum=lyrNum)
    ze_field%TempDataFlag = .true.

    !$omp parallel do
    do e=1, nEdge
       ze_field%data%v_(1:lyrNum,e) = 0.5d0*sum(zc_field%data%v_(1:lyrNum, fvm%Face_CellId(1:2,e)), 2)
    end do

    if(zc_field%TempDataFlag) call Release(zc_field)
    
  end function ze_zc



  !
  !
  !
  function calcKineticEnergy(ze_normalVel) result(zc_KE)

    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: ze_normalVel
    type(volScalarField) :: zc_KE

    ! 局所変数
    ! Local variables
    !
    integer :: i, e
    integer :: nEdge, nCell
    integer :: lfaceNum
    real(DP), allocatable :: ze_localKE(:,:)
    type(PolyMesh), pointer :: mesh
    type(fvMeshInfo), pointer :: fvm

    ! 実行文; Executable statement
    !

#ifdef DEBUG
    if( .not. associated(defaultUtilObj) ) &
         & call MessageNotify("E", module_name, "Set the pointer of default object to actual object.")
#endif

    mesh => defaultUtilObj%mesh
    fvm => defaultUtilObj%fvmInfo

    call GridSet_getLocalMeshInfo(mesh, nEdge=nEdge, nCell=nCell)
    call GeometricField_Init(zc_KE, mesh, "KE")
    zc_KE%TempDataFlag = .true.

    allocate(ze_localKE(nVzLyr, nEdge))
    !$omp parallel do
    do e=1, nEdge
       ze_localKE(1:nVzLyr,e) = 0.25d0 * l2norm( fvm%s_faceAreaVec%data%v_(1,e) ) &
            & * fvm%s_dualMeshFaceArea%data%v_(1,e) * ze_normalVel%data%v_(1:nVzLyr,e)**2
    end do

    !$omp parallel do private(lfaceNum)
    do i=1, nCell
       lfaceNum = mesh%CellList(i)%faceNum
       zc_KE%data%v_(1:nVzLyr,i) = &
            & sum( ze_localKE(1:nVzLyr, fvm%Cell_FaceId(1:lfaceNum,i)), 2) &
            & / fvm%v_CellVol%data%v_(1:nVzLyr,i)
    end do

  end function CalcKineticEnergy

  function CalcPVFlux(zc_h, ze_normalVel) result(ze_PVflux)

    use SphericalCoord_mod, only: &
         & CartToSphPos

    use Constants_mod, only: &
         & Omega
    
    ! 宣言文; Declaration statement
    !

    type(volScalarField), intent(in) :: zc_h
    type(surfaceScalarField), intent(in) :: ze_normalVel
    type(surfaceScalarField) :: ze_PVFlux

    ! 局所変数
    ! Local variables
    !
    type(PointScalarField) :: zv_planetVor
    type(PointScalarField) :: zv_height
    type(PointScalarField) :: zv_PV

    integer :: e, e2, k, i
    real(DP) :: weight, tmp(nVzLyr), s2_MomFlux(nVzLyr), s_PV(nVzLyr), s2_PV(nVzLyr)
    integer :: v, cellIds(3), cellId
    type(Vector3d) :: geoPos
    integer :: nVertex, nEdge

    type(CGridFieldDataUtil), pointer :: this
    type(PolyMesh), pointer :: mesh
    type(fvMeshInfo), pointer :: fvm

    ! 実行文; Executable statement
    !

#ifdef DEBUG
    if( .not. associated(defaultUtilObj) ) &
         & call MessageNotify("E", module_name, "Set the pointer of default object to actual object.")
#endif

    this => defaultUtilObj
    mesh => this%mesh
    fvm => this%fvmInfo

    call GridSet_getLocalMeshInfo(mesh, nEdge=nEdge, nVertex=nVertex)
    call GeometricField_Init(zv_planetVor, mesh, "p_planetVor")
    call GeometricField_Init(zv_PV, mesh, "p_eta")
    call GeometricField_Init(zv_height, mesh, "p_height")

    !$omp parallel do private(geoPos, cellIds)
    do v=1, nVertex
       geoPos = CartToSphPos( mesh%PointPosList(v) )
       zv_planetVor%data%v_(1:nVzLyr,v) = (2d0*Omega)*sin(geoPos%v_(2))

       cellIds(:) = this%fvmInfo%Point_CellId(1:3,v)

       forAll(k=1:nVzLyr) &
            & zv_height%data%v_(k,v) = &
            & sum( this%R_iv(1:3,v) * fvm%v_cellVol%data%v_(k,cellIds(:)) * zc_h%data%v_(k,cellIds(:)) ) &
            & / fvm%p_dualMeshCellVol%data%v_(k,v)

    end do

    zv_PV = ( zv_planetVor + curl(ze_normalVel) )/ zv_height

    !
    !

    call GeometricField_Init(ze_PVFlux, mesh, "s_PVFlux")
    ze_PVFlux%TempDataFlag = .true.

    !$omp parallel do private(s_PV, s2_PV, s2_MomFlux, tmp, k, e2) 
    do e=1, nEdge
       s_PV(:) = 0.5d0*sum( zv_PV%data%v_(1:nVzLyr,fvm%Face_PointId(1:2,e)), 2)
       tmp = 0d0

       do k=1, size(fvm%Face_PairCellFaceId,1)
          e2 = fvm%Face_PairCellFaceId(k, e)
          if( e2 == -1) exit;

          s2_MomFlux(:) = 0.5d0*sum( zc_h%data%v_(1:nVzLyr, fvm%Face_CellId(1:2,e2)), 2) &
               &           * ze_normalVel%data%v_(1:nVzLyr,e2)
          s2_PV(:) = 0.5d0*sum( zv_PV%data%v_(1:nVzLyr, fvm%Face_PointId(1:2,e2)), 2)

          tmp(:) =    tmp(:) &
               & + this%wt_ff(k,e) * l2norm(fvm%s_faceAreaVec%data%v_(1,e2)) * s2_MomFlux * 0.5d0*(s_PV + s2_PV) 


       end do ! do while

       ze_PVFlux%data%v_(1:nVzLyr,e) = tmp(:) / fvm%s_dualMeshFaceArea%data%v_(1,e)

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

    integer :: e, e2, k
    real(DP) :: weight, tmp(nVzLyr)
    integer :: nEdge
    type(CGridFieldDataUtil), pointer :: this
    type(PolyMesh), pointer :: mesh
    type(fvMeshInfo), pointer :: fvm

    ! 実行文; Executable statement
    !

#ifdef DEBUG
    if( .not. associated(defaultUtilObj) ) &
         & call MessageNotify("E", module_name, "Set the pointer of default object to actual object.")
#endif

    this => defaultUtilObj
    mesh => this%mesh
    fvm => this%fvmInfo


    call GridSet_getLocalMeshInfo(mesh, nEdge=nEdge)
    call GeometricField_Init(ze_tangentVel, mesh, "ze_tangentVel")
    ze_tangentVel%TempDataFlag = .true.

    !$omp parallel do private(tmp, k, e2) 
    do e=1, nEdge
       tmp = 0d0

       do k=1, size(fvm%Face_PairCellFaceId,1)
          e2 = fvm%Face_PairCellFaceId(k, e)
          if( e2 == -1) exit;

          tmp(:) =    tmp(:)  + &
               & this%wt_ff(k,e) * l2norm(fvm%s_faceAreaVec%data%v_(1,e2)) * ze_normalVel%data%v_(1:nVzLyr, e2)

       end do ! do while

       ze_tangentVel%data%v_(1:nVzLyr,e) = tmp(:) / fvm%s_dualMeshFaceArea%data%v_(1,e)
    end do

    !
    !

  end function ToTangentVel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepairIntrpWeight(this)
    
    ! モジュール引用; Use statements
    !
    use SphericalCoord_mod, only: &
         & SphericalTriArea

    ! 宣言文; Declaration statement
    !
    type(CGridFieldDataUtil), intent(inout) :: this
        
    ! 局所変数
    ! Local variables
    !

    integer :: lfaceNum
    integer :: cellId, cellIds(3), ptId, faceId, lfaceId
    integer :: m, n, k
    real(DP) :: areas(3)
    type(Vector3d) :: ptPos
    real(DP), allocatable :: R_vi(:,:)
    integer :: startPtId, startlPtId, startlfaceId, t_fv2
    real(DP) :: lR_vi(MAX_CELL_FACE_NUM), ln_fv(MAX_CELL_FACE_NUM)
    type(fvMeshInfo), pointer :: fvmInfo
    integer :: nVertex, nEdge, nCell

    ! 実行文; Executable statements
    !

    call MessageNotify('M', module_name, "prepare some data about weight required interpolation..")

    fvmInfo => this%fvmInfo
    call GridSet_getLocalMeshInfo(this%mesh, nCell=nCell, nEdge=nEdge, nVertex=nVertex)
    allocate( this%R_iv(3, nVertex) )
    allocate( this%Wt_ff(MAX_CELL_FACE_NUM*2, nEdge) )
    allocate( R_vi(MAX_CELL_FACE_NUM, nCell) )

    R_vi = 0d0
    do ptId=1, nVertex
       ptPos = this%mesh%PointPosList(ptId)
       cellIds = cshift(fvmInfo%Point_CellId(1:3,ptId), -1)

       do m=1, 3
          areas(m) = sphericalTriArea(ptPos, this%mesh%CellPosList(cellIds(m)), this%mesh%CellPosList(cellIds(mod(m,3)+1)) )
       end do

       do m=1, 3
          cellId = fvmInfo%Point_CellId(m, ptId) 
          this%R_iv(m, ptId) = 0.5d0*( areas(m) + areas(mod(m,3)+1) )/ At(fvmInfo%v_CellVol, 1, cellId) 

          do n=1, MAX_CELL_FACE_NUM
             if( fvmInfo%Cell_PointId(n,cellId) == ptId) then
                R_vi(n, cellId) = this%R_iv(m, ptId); exit
             end if
          end do

       end do
    end do

    this%Wt_ff = 0d0
    do faceId=1, nEdge   

       k=2
       do m=1, 2
          cellId = fvmInfo%Face_CellId(m, faceId)
          startPtId = fvmInfo%Face_PointId(mod(m,2)+1, faceId)
          lfaceNum = this%mesh%CellList(cellId)%faceNum
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
             this%Wt_ff(k, faceId) =  dble(ln_fv(n)*t_fv2) * (sum(lR_vi(1:n-1)) - 0.5d0)
             k = k + 1
          end do

       end do

       ! write(*,'(a, i6, a, 13f12.5)') "wt_ff", faceId, ":", wt_ff(:,faceId), sum(wt_ff(:,faceId))
    end do

  end subroutine prepairIntrpWeight


end module CGridFieldDataUtil_mod

