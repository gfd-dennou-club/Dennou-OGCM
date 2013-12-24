!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
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
  public :: Advance_VDiffProc
  public :: construct_vDiffProcMat, Solve

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverVDiffProc_mod' !< Module Name
  real(DP), allocatable :: vDiffProcMat(:,:,:)
  integer, allocatable :: vDiffProcMatKp(:,:)


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

    deallocate(vDiffProcMatKp, vDiffProcMat)

  end subroutine HydroBouEqSolverVDiffProc_Final


  !> @brief 
  !!
  !!
  subroutine Advance_VDiffProc(wt_Vor, wt_Div, wt_PTempEdd, &
    & xy_WindStressU, xy_WindStressV, xy_totDepth, Av, dt,  &
    & DynBCSurf, DynBCBottom )
    
    !
    !
    use Constants_mod, only: refDens

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wt_Vor(lMax,0:tMax)
    real(DP), intent(inout) :: wt_Div(lMax,0:tMax)
    real(DP), intent(inout) :: wt_PTempEdd(lMax,0:tMax)
    real(DP), intent(in) :: xy_WindStressU(0:iMax-1,jMax), xy_WindStressV(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: Av
    real(DP), intent(in) :: dt
    integer, intent(in) :: DynBCSurf
    integer, intent(in) :: DynBCBottom
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_Work(0:iMax-1,jMax,0:kMax)
    real(DP),dimension(iMax*jMax,0:kMax) :: workRHS
    integer :: l, k
    real(DP) :: xy_CosLat(0:iMax-1,jMax)
    character :: BCKindUpper, BCKindBottom
    real(dp) :: xy_DivSig(0:iMax-1,jMax)
    real(dp) :: xyz_Div(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !

    xy_CosLat = cos(xyz_Lat(:,:,1))
    

    !
    !


    !
    !

    xyz_Work = xyz_wt(wt_Vor)
    if (DynBCSurf == DynBCTYPE_NoSlip) then
       BCKindUpper = 'N'
       xyz_Work(:,:,0) = xy_w( &
            & w_AlphaOptr_xy(xy_WindStressV*xy_CosLat, -xy_WindStressU*xy_CosLat) &
            & )/(RefDens*Av)
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

    call construct_vDiffProcMat(Av, dt, xy_totDepth, BCKindUpper, BCKindBottom)
    wt_Vor = solve(xyz_Work)

    xyz_Work = xyz_wt(wt_Div)
    if (DynBCSurf == DynBCTYPE_NoSlip) then 
       xyz_Work(:,:,0) = xy_w( &
            & w_AlphaOptr_xy(xy_WindStressU*xy_CosLat, xy_WindStressV*xy_CosLat) &
            & )/(RefDens*Av)
    else
       xyz_Work(:,:,0) = 0d0
    end if
    if (DynBCBottom == DynBCTYPE_NoSlip) then
       xyz_Work(:,:,kMax) = 0d0 
    else
       xyz_Work(:,:,kMax) = 0d0
    end if
    wt_Div = solve(xyz_Work)
!!$    xyz_Div = xyz_wt(wt_Div)
!!$    xy_DivSig = xy_IntSig_BtmToTop_xyz(xyz_Div)
!!$    forAll(k=0:kMax) xyz_Div(:,:,k) = xyz_Div(:,:,k) - xy_DivSig
!!$    wt_Div = wt_xyz(xyz_Div)

    call construct_vDiffProcMat(Av, dt, xy_totDepth, 'N', 'N')
    xyz_Work = xyz_wt(wt_PTempEdd)
    xyz_Work(:,:,0) = 0d0 
    xyz_Work(:,:,kMax) = 0d0 
    wt_PTempEdd = solve(xyz_Work)

  end subroutine Advance_VDiffProc

  function Solve(xyz_RHS) result(wt_ret)

    use lumatrix, only: lusolve

    real(DP), intent(in) :: xyz_RHS(0:iMax-1,jMax,0:kMax)
    real(DP) :: wt_ret(lMax,0:tMax)

    wt_ret = wa_xya( reshape( &
         & lusolve( vDiffProcMat, vDiffProcMatKp, reshape(xyz_RHS, (/ iMax*jMax,kMax+1 /)) ), &
         & shape(xyz_RHS) ))
  end function Solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine construct_vDiffProcMat(Av, dt, xy_totDepth, &
       & BCKindUpper, BCKindBottom )
    
    !
    !
    use at_module
    use lumatrix

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: Av
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
    end if

    Depth(:) = reshape( xy_totDepth, shape(Depth) )

    tt_I = 0d0
    do t=0,tMax
       tt_I(t,t) = 1d0
    end do
    tg_data = ag_at(tt_I)

    dwork_tg_data = ag_at( at_Dx_at(at_Dx_at(tt_I)) )
    do k=1,kMax-1
       forAll(t=0:tMax) &
            & vDiffProcMat(:,k,t) = tg_data(t,k) - dt*Av/Depth**2 * dwork_tg_data(t,k)
    end do

    dwork_tg_data = ag_at( at_Dx_at(tt_I) )
    if(BCKindUpper == 'N') then
       forall(t=0:tMax) vDiffProcMat(:,0,t) = dwork_tg_data(t,0)/Depth
    else
       forall(t=0:tMax) vDiffProcMat(:,0,t) = tg_data(t,0)
    end if

    if(BCKindBottom == 'N') then
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = dwork_tg_data(t,kMax)/Depth
    else
       forall(t=0:tMax) vDiffProcMat(:,kMax,t) = tg_data(t,kMax)
    end if

    call ludecomp(vDiffProcMat, vDiffProcMatKp)
 
  end subroutine construct_vDiffProcMat

end module HydroBouEqSolverVDiffProc_mod


