!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
module HydroBouEqSolverVDiffProc_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP

  use dc_message, only: &
       & MessageNotify

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use BoundCondSet_mod, only: &
       & DynBCTYPE_Slip, DynBCTYPE_NoSlip, &
       & KinBCTYPE_FreeSurf, KinBCTYPE_RigidLid

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolverVDiffProc_Init, HydroBouEqSolverVDiffProc_Final
  public :: HydroBouEqSolverVDiffProc_constructMat
  public :: HydroBouEqSolverVDiffProc_Advance

  public :: Advance_VDiffProc

  public :: construct_vDiffProcMat, Solve

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverVDiffProc_mod' !< Module Name

  real(DP), allocatable :: vDiffProcMatVor(:,:,:), vDiffProcMatDiv(:,:,:), vDiffProcMatHeat(:,:,:)
  integer, allocatable :: vDiffProcMatVorKp(:,:), vDiffProcMatDivKp(:,:), vDiffProcMatHeatKp(:,:)

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverVDiffProc_Init()

    !
    !

    ! 実行文; Executable statements
    !

    
  end subroutine HydroBouEqSolverVDiffProc_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolverVDiffProc_Final()

    ! 実行文; Executable statements
    !

    if(allocated(vDiffProcMatVorKp)) deallocate(vDiffProcMatVorKp, vDiffProcMatVor)
    if(allocated(vDiffProcMatDivKp)) deallocate(vDiffProcMatDivKp, vDiffProcMatDiv)
    if(allocated(vDiffProcMatHeatKp)) deallocate(vDiffProcMatHeatKp, vDiffProcMatHeat)

  end subroutine HydroBouEqSolverVDiffProc_Final


  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolverVDiffProc_constructMat( & 
       & vViscTermCoef, vDiffTermCoef, dt, &
       & xy_totDepth, DynBCSurf, DynBCBottom )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: vViscTermCoef, vDiffTermCoef
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    integer, intent(in) :: DynBCSurf
    integer, intent(in) :: DynBCBottom    
    
    ! 局所変数
    ! Local variables
    !
    character :: BCKindUpper, BCKindBottom
    
    ! 実行文; Executable statement
    !

    if (DynBCSurf == DynBCTYPE_NoSlip) then
       BCKindUpper = 'N'
    else
       BCKindUpper = 'N'
    end if
    if (DynBCBottom == DynBCTYPE_NoSlip) then
       BCKindBottom = 'D'
    else
       BCKindBottom = 'N'
    end if

    call construct_vDiffProcMat(vDiffProcMatVor, vDiffProcMatVorKp, &   !(inout)
         & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )             !(in)

    call construct_vDiffProcMat2(vDiffProcMatDiv, vDiffProcMatDivKp, &   !(inout)
         & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )              !(in)
!!$    call construct_vDiffProcMat(vDiffProcMatDiv, vDiffProcMatDivKp, &   !(inout)
!!$         & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )             !(in)
!!$    wt_Div = solve(vDiffProcMatDiv, vDiffProcMatDivKp, xyz_Work(:,:,0:kMax))

    !
    call construct_vDiffProcMat(vDiffProcMatHeat, vDiffProcMatHeatKp,&
         & vDiffTermCoef, dt, xy_totDepth, 'N', 'N')

    
  end subroutine HydroBouEqSolverVDiffProc_constructMat

  !> @brief 
  !!
  !!
  subroutine HydroBouEqSolverVDiffProc_Advance(wt_Vor, wt_Div, wt_PTempEdd, &
    & xy_WindStressU, xy_WindStressV, xy_totDepth, &
    & vViscCoef, &
    & DynBCSurf, DynBCBottom )
    
    !
    !
    use Constants_mod, only: refDens
use at_module_omp
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), optional :: wt_Vor(lMax,0:tMax)
    real(DP), intent(inout), optional :: wt_Div(lMax,0:tMax)
    real(DP), intent(inout), optional :: wt_PTempEdd(lMax,0:tMax)
    real(DP), intent(in) :: xy_WindStressU(0:iMax-1,jMax), xy_WindStressV(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: vViscCoef
    integer, intent(in) :: DynBCSurf
    integer, intent(in) :: DynBCBottom
!    logical, intent(in) :: isVMatReconst

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_Work(0:iMax-1,jMax,0:kMax+1)
    real(DP),dimension(iMax*jMax,0:kMax) :: workRHS
    integer :: l, k
    real(DP) :: xy_CosLat(0:iMax-1,jMax)
    real(dp) :: xy_DivSig(0:iMax-1,jMax)
    real(dp) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    logical :: isVMatReconstFlag

    ! 実行文; Executable statement
    !

    xy_CosLat = cos(xyz_Lat(:,:,0))
    
   
    !
    !


    !
    !
    if( present(wt_Vor) ) then
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_Vor)
       if (DynBCSurf == DynBCTYPE_NoSlip) then
          xyz_Work(:,:,0) = xy_w( &
               & w_AlphaOptr_xy(xy_WindStressV*xy_CosLat, -xy_WindStressU*xy_CosLat) &
               & )/(RefDens*vViscCoef)
       else
          xyz_Work(:,:,0) = 0d0
       end if
       if (DynBCBottom == DynBCTYPE_NoSlip) then
          xyz_Work(:,:,kMax) = 0d0 
       else
          xyz_Work(:,:,kMax) = 0d0
       end if

       wt_Vor(:,:) = solve(vDiffProcMatVor, vDiffProcMatVorKp, xyz_Work(:,:,0:kMax))
    end if

    !
    if( present(wt_Div) ) then
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_Div)
       xyz_Work(:,:,kMax+1) = 0d0!xy_IntSig_BtmToTop_xyz(xyz_Work(:,:,0:kMax)) !0d0

       if (DynBCSurf == DynBCTYPE_NoSlip) then 
          xyz_Work(:,:,0) = xy_w( &
               & w_AlphaOptr_xy(xy_WindStressU*xy_CosLat, xy_WindStressV*xy_CosLat) &
               & )/(RefDens*vViscCoef)
       else
          xyz_Work(:,:,0) = 0d0
       end if
       if (DynBCBottom == DynBCTYPE_NoSlip) then
          xyz_Work(:,:,kMax) = 0d0 
       else
          xyz_Work(:,:,kMax) = 0d0
       end if

       wt_Div(:,:) = solve(vDiffProcMatDiv, vDiffProcMatDivKp, xyz_Work(:,:,0:kMax+1))
    end if

    if( present(wt_PTempEdd) ) then
       !
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_PTempEdd)
       xyz_Work(:,:,0) = 0d0 
       xyz_Work(:,:,kMax) = 0d0 
       wt_PTempEdd(:,:) = solve(vDiffProcMatHeat, vDiffProcMatHeatKp, xyz_Work(:,:,0:kMax))
    end if

  end subroutine HydroBouEqSolverVDiffProc_Advance

  !> @brief 
  !!
  !!
  subroutine Advance_VDiffProc(wt_Vor, wt_Div, wt_PTempEdd, &
    & xy_WindStressU, xy_WindStressV, xy_totDepth, &
    & vViscTermCoef, vDiffTermCoef, vViscCoef, dt,  &
    & DynBCSurf, DynBCBottom )
    
    !
    !
    use Constants_mod, only: refDens
use at_module_omp
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), optional :: wt_Vor(lMax,0:tMax)
    real(DP), intent(inout), optional :: wt_Div(lMax,0:tMax)
    real(DP), intent(inout), optional :: wt_PTempEdd(lMax,0:tMax)
    real(DP), intent(in) :: xy_WindStressU(0:iMax-1,jMax), xy_WindStressV(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: vViscTermCoef, vDiffTermCoef, vViscCoef
    real(DP), intent(in) :: dt
    integer, intent(in) :: DynBCSurf
    integer, intent(in) :: DynBCBottom
!    logical, intent(in) :: isVMatReconst

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_Work(0:iMax-1,jMax,0:kMax+1)
    real(DP),dimension(iMax*jMax,0:kMax) :: workRHS
    integer :: l, k
    real(DP) :: xy_CosLat(0:iMax-1,jMax)
    character :: BCKindUpper, BCKindBottom
    real(dp) :: xy_DivSig(0:iMax-1,jMax)
    real(dp) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    logical :: isVMatReconstFlag

    ! 実行文; Executable statement
    !

    xy_CosLat = cos(xyz_Lat(:,:,0))
    
   
    !
    !


    !
    !
    if( present(wt_Vor) ) then
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_Vor)
       if (DynBCSurf == DynBCTYPE_NoSlip) then
          BCKindUpper = 'N'
          xyz_Work(:,:,0) = xy_w( &
               & w_AlphaOptr_xy(xy_WindStressV*xy_CosLat, -xy_WindStressU*xy_CosLat) &
               & )/(RefDens*vViscCoef)
       else
          BCKindUpper = 'N'
          xyz_Work(:,:,0) = 0d0
       end if
       if (DynBCBottom == DynBCTYPE_NoSlip) then
          BCKindBottom = 'D'
          xyz_Work(:,:,kMax) = 0d0 
       else
          BCKindBottom = 'N'
          xyz_Work(:,:,kMax) = 0d0
       end if

       call construct_vDiffProcMat(vDiffProcMatVor, vDiffProcMatVorKp, &   !(inout)
            & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )             !(in)
       wt_Vor(:,:) = solve(vDiffProcMatVor, vDiffProcMatVorKp, xyz_Work(:,:,0:kMax))
    end if

    !
    if( present(wt_Div) ) then
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_Div)
       xyz_Work(:,:,kMax+1) = xy_IntSig_BtmToTop_xyz(xyz_Work(:,:,0:kMax)) !0d0

       if (DynBCSurf == DynBCTYPE_NoSlip) then 
          BCKindUpper = 'N'
          xyz_Work(:,:,0) = xy_w( &
               & w_AlphaOptr_xy(xy_WindStressU*xy_CosLat, xy_WindStressV*xy_CosLat) &
               & )/(RefDens*vViscCoef)
       else
          BCKindUpper = 'N'
          xyz_Work(:,:,0) = 0d0
       end if
       if (DynBCBottom == DynBCTYPE_NoSlip) then
          BCKindUpper = 'D'
          xyz_Work(:,:,kMax) = 0d0 
       else
          BCKindUpper = 'N'
          xyz_Work(:,:,kMax) = 0d0
       end if

       call construct_vDiffProcMat2(vDiffProcMatDiv, vDiffProcMatDivKp, &   !(inout)
            & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )              !(in)
       wt_Div(:,:) = solve(vDiffProcMatDiv, vDiffProcMatDivKp, xyz_Work(:,:,0:kMax+1))
!!$    call construct_vDiffProcMat(vDiffProcMatDiv, vDiffProcMatDivKp, &   !(inout)
!!$         & vViscTermCoef, dt, xy_totDepth, BCKindUpper, BCKindBottom )             !(in)
!!$    wt_Div = solve(vDiffProcMatDiv, vDiffProcMatDivKp, xyz_Work(:,:,0:kMax))

    end if

    if( present(wt_PTempEdd) ) then
       !
       call construct_vDiffProcMat(vDiffProcMatHeat, vDiffProcMatHeatKp,&
            & vDiffTermCoef, dt, xy_totDepth, 'N', 'N')
       xyz_Work(:,:,0:kMax) = xyz_wt(wt_PTempEdd)
       xyz_Work(:,:,0) = 0d0 
       xyz_Work(:,:,kMax) = 0d0 
       wt_PTempEdd(:,:) = solve(vDiffProcMatHeat, vDiffProcMatHeatKp, xyz_Work(:,:,0:kMax))
    end if

  end subroutine Advance_VDiffProc


  function Solve(vDiffProcMat, vDiffProcMatKp, xyz_RHS) result(wt_ret)

    use lumatrix, only: lusolve
    use omp_lib
use TemporalIntegSet_mod, only: CurrentTime
use Constants_mod, only: RefDens

    real(DP), intent(in) :: vDiffProcMat(:,:,:)
    integer, intent(in) :: vDiffProcMatKp(:,:)
    real(DP), intent(in) :: xyz_RHS(:,:,:)
    real(DP) :: xyt_retTmp(0:iMax-1,jMax,0:size(vDiffProcMat,3)-1), wt_ret(lMax,0:tMax)

    real(DP) :: az_RHS(size(xyz_RHS,1)*size(xyz_RHS,2),size(xyz_RHS,3))
!!$    real(DP) :: xy_SurfPress(0:iMax-1,jMax)

    !$omp parallel workshare
    az_RHS = reshape(xyz_RHS, shape(az_RHS))
    !$omp end parallel workshare

    xyt_retTmp(:,:,:) = reshape( &
         & lusolve( vDiffProcMat, vDiffProcMatKp,  az_RHS ), &
         & shape(xyz_RHS) )

    wt_ret = wa_xya( xyt_retTmp(:,:,0:tMax) )


!!$if(size(xyz_RHS,3)==kMax+2.and.mod(int(CurrentTime),43200)==0)then
!!$   xy_SurfPress = xy_w(w_InvLapla2D_w ( w_xy(xyt_retTmp(:,:,tMax+1)) ))*RefDens
!!$   write(*,*) "surfPress:", xy_SurfPress(0,:)
!!$   write(*,*) "intSIg(D):", xy_IntSig_BtmToTop_xyz(xyz_wt(wt_ret))
!!$endif

  end function Solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine construct_vDiffProcMat2(vDiffProcMat, vDiffProcMatKp, &
    & vDiffTermCoef, dt, xy_totDepth, &
    & BCKindUpper, BCKindBottom )
    
    !
    !
    use lumatrix
    use at_module_omp
    use omp_lib

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), allocatable :: vDiffProcMat(:,:,:)
    integer, intent(inout), allocatable :: vDiffProcMatKp(:,:)
    real(DP), intent(in) :: vDiffTermCoef
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    character, intent(in) :: BCKindUpper, BCKindBottom ! N or D
    
    ! 局所変数
    ! Local variables
    !
    integer :: t, k
    real(DP) :: tt_data(0:tMax,0:tMax)
    real(DP) :: tg_data(0:tMax,0:kMax)
    real(DP) :: dwork_tg_data(0:tMax,0:kMax), iwork_tg_data(0:tMax)
    real(DP) :: tt_I(0:tMax,0:tMax)
    real(DP) :: Depth(iMax*jMax)


    ! 実行文; Executable statement
    !
    
    if(.not. allocated(vDiffProcMat)) then
       allocate( vDiffProcMat(iMax*jMax,0:kMax+1,0:tMax+1), vDiffProcMatKp(iMax*jMax,0:kMax+1) )
!   else
!       return
    end if

    vDiffProcMat = 0d0
    Depth(:) = reshape( xy_totDepth, shape(Depth) )
    
    
    tt_I = 0d0
    do t=0,tMax
       tt_I(t,t) = 1d0
    end do
    tg_data = ag_at(tt_I)

    dwork_tg_data = ag_at( at_Dx_at(at_Dx_at(tt_I)) )
    do k=1,kMax-1
       !$omp parallel workshare
       forAll(t=0:tMax) &
            & vDiffProcMat(:,k,t) = tg_data(t,k) - dt*vDiffTermCoef/Depth**2 * dwork_tg_data(t,k)
       !$omp end parallel workshare
       vDiffProcMat(:,k,tMax+1) = 0.5d0*dt
    end do

    dwork_tg_data = ag_at( at_Dx_at(tt_I) )
    if(BCKindUpper == 'N') then
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,0,t) = dwork_tg_data(t,0)/Depth
       !$omp end parallel workshare
    else
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,0,t) = tg_data(t,0)
       !$omp end parallel workshare
    end if

    if(BCKindBottom == 'N') then
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = dwork_tg_data(t,kMax)/Depth
       !$omp end parallel workshare
    else
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = tg_data(t,kMax)
       !$omp end parallel workshare
    end if

    iwork_tg_data = matmul(g_X_WEIGHT, tg_data)
    forall(t=0:tMax) vDiffProcMat(:,kMax+1,t) = dot_product(g_X_WEIGHT, tg_data(t,:))

!!$    forall(t=0:tMax) &
!!$         & vDiffProcMat(:,kMax+1,t) = - dt*vDiffTermCoef/Depth**2 * (dwork_tg_data(t,0) - dwork_tg_data(t,kMax))
!!$    vDiffProcMat(:,kMax+1,tMax+1) = 0.5d0*dt

    call ludecomp(vDiffProcMat, vDiffProcMatKp)

  end subroutine construct_vDiffProcMat2


  !> @brief 
  !!
  !!
  subroutine construct_vDiffProcMat(vDiffProcMat, vDiffProcMatKp, &
    & vDiffTermCoef, dt, xy_totDepth, &
    & BCKindUpper, BCKindBottom )

    !
    !
    use lumatrix
    use at_module_omp!_omp

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), allocatable :: vDiffProcMat(:,:,:)
    integer, intent(inout), allocatable :: vDiffProcMatKp(:,:)
    real(DP), intent(in) :: vDiffTermCoef
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    character, intent(in) :: BCKindUpper, BCKindBottom ! N or D
    
    ! 局所変数
    ! Local variables
    !
    integer :: t, k
    real(DP) :: tt_data(0:tMax,0:tMax)
    real(DP) :: tg_data(0:tMax,0:kMax), dwork_tg_data(0:tMax,0:kMax)
    real(DP) :: tt_I(0:tMax,0:tMax)
    real(DP) :: Depth(iMax*jMax)


    ! 実行文; Executable statement
    !
    
    if(.not. allocated(vDiffProcMat)) then
       allocate( vDiffProcMat(iMax*jMax,0:kMax,0:tMax), vDiffProcMatKp(iMax*jMax,0:kMax) )
!    else
!       return
    end if

    vDiffProcMat = 0d0
    Depth(:) = reshape( xy_totDepth, shape(Depth) )

    tt_I = 0d0
    do t=0,tMax
       tt_I(t,t) = 1d0
    end do
    tg_data = ag_at(tt_I)

    dwork_tg_data = ag_at( at_Dx_at(at_Dx_at(tt_I)) )
    do k=1,kMax-1
       !$omp parallel workshare
       forAll(t=0:tMax) &
            & vDiffProcMat(:,k,t) = tg_data(t,k) - dt*vDiffTermCoef/Depth**2 * dwork_tg_data(t,k)
       !$omp end parallel workshare
    end do

    dwork_tg_data = ag_at( at_Dx_at(tt_I) )
    if(BCKindUpper == 'N') then
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,0,t) = dwork_tg_data(t,0)/Depth
       !$omp end parallel workshare
    else
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,0,t) = tg_data(t,0)
       !$omp end parallel workshare
    end if

    if(BCKindBottom == 'N') then
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = dwork_tg_data(t,kMax)/Depth
       !$omp end parallel workshare
    else
       !$omp parallel workshare
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = tg_data(t,kMax)
       !$omp end parallel workshare
    end if

    call ludecomp(vDiffProcMat, vDiffProcMatKp)
  end subroutine construct_vDiffProcMat

end module HydroBouEqSolverVDiffProc_mod


