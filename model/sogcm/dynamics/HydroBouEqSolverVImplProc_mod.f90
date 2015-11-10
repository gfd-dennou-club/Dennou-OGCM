!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolverVImplProc_mod 

  ! モジュール引用; Use statements
  !

  use dc_types, only: DP

  use dc_message, only: MessageNotify
  
  use GridSet_mod, only: &
       & iMax, jMax, kMax, tMax, lMax, nMax, &
       & xyz_Lat

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType

  use Constants_mod, only: &
       & refDens, &
       & Omega,   &
       & vViscCoef, vHyperViscCoef, vDiffCoef, vHyperDiffCoef

  use SpmlUtil_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolverVImplProc_Init, HydroBouEqSolverVImplProc_Final
  public :: HydroBouEqSolverVImplProc_Prepare
!!$  public :: Advance_VImplicitProc
  public :: Advance_VImplicitProc_DeltaForm
  public :: Advance_VImplicitProc_DeltaForm_xyz

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverVImplProc_mod' !< Module Name
  
  real(DP), dimension(:,:), allocatable :: gt_Mat, gt_DSig1Mat, gt_DSig2Mat, gt_DSig4Mat
  real(DP), dimension(:), allocatable :: t_IntSigMat

  logical :: isDifferentialMatInit

  real(DP) :: vViscDiffTermCoef
  real(DP) :: vHyperViscDiffTermCoef
  real(DP) :: CoriTermCoef
  real(DP) :: dt
  
contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverVImplProc_Init()

    ! 実行文; Executable statements
    !

    isDifferentialMatInit = .false.
    vViscDiffTermCoef = 0d0
    vHyperViscDiffTermCoef = 0d0
    CoriTermCoef = 0d0
    
  end subroutine HydroBouEqSolverVImplProc_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolverVImplProc_Final()

    ! 実行文; Executable statements
    !
    
    if(isDifferentialMatInit) then
       deallocate(gt_Mat, gt_DSig1Mat, gt_DSig2Mat, gt_DSig4Mat, t_IntSigMat)
    end if

    
  end subroutine HydroBouEqSolverVImplProc_Final

  subroutine HydroBouEqSolverVImplProc_Prepare( &
       & DelTimeImpl, CoriTermTAvgCoefA, vViscDiffTermTAvgCoefA, vHyperViscDiffTermTAvgCoefA )
    
    use at_module_omp
    real(DP), intent(in) :: DelTimeImpl
    real(DP), intent(in) :: CoriTermTAvgCoefA, vViscDiffTermTAvgCoefA, vHyperViscDiffTermTAvgCoefA
    
    real(DP), dimension(0:tMax,0:tMax) :: tt_I
    integer :: t

    dt = DelTimeImpl
    CoriTermCoef = CoriTermTAvgCoefA
    vViscDiffTermCoef = vViscDiffTermTAvgCoefA
    vHyperViscDiffTermCoef = vHyperViscDiffTermTAvgCoefA
    
    if(isDifferentialMatInit) then
       deallocate(gt_Mat, gt_DSig1Mat, gt_DSig2Mat, gt_DSig4Mat, t_IntSigMat)
    end if

    allocate(gt_Mat(0:kMax,0:tMax), gt_DSig1Mat(0:kMax,0:tMax), gt_DSig2Mat(0:kMax,0:tMax), gt_DSig4Mat(0:kMax,0:tMax), &
         &   t_IntSigMat(0:tMax))

    tt_I = 0d0
    do t=0, tMax
       tt_I(t,t) = 1d0
    end do
    gt_Mat = tt_I
    
    tt_I = at_ag(tt_I)
    gt_DSig1Mat = ag_at( at_Dx_at(tt_I) )
    gt_DSig2Mat = ag_at( at_Dx_at(at_Dx_at(tt_I)) )
    gt_DSig4Mat = ag_at( at_Dx_at(at_Dx_at(at_Dx_at(at_Dx_at(tt_I)))) )
    forAll(t=0:tMax) t_IntSigMat(t) = dot_product(g_X_WEIGHT, gt_Mat(t,:))

    isDifferentialMatInit = .true.

  end subroutine HydroBouEqSolverVImplProc_Prepare

 
!!$  !> @brief 
!!$  !!
!!$  !!
!!$  subroutine Advance_VImplicitProc(wz_Vor, wz_Div, wz_PTempEdd, wz_Salt, xy_SurfPress, &
!!$    & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, xy_totDepth, &
!!$    & DynBCSurf, DynBCBottom, ThermBCSurf, ThermBCBottom, SaltBCSurf, SaltBCBottom )
!!$    
!!$    !
!!$    !
!!$
!!$use at_module_omp
!!$
!!$    ! 宣言文; Declaration statement
!!$    !
!!$    real(DP), dimension(lMax,0:kMax), intent(inout)  :: &
!!$         & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt
!!$    real(DP), intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
!!$    real(DP), dimension(lMax,2), intent(in) :: &
!!$         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
!!$    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
!!$    integer, intent(in) :: DynBCSurf, DynBCBottom
!!$    integer, intent(in) :: ThermBCSurf, ThermBCBottom
!!$    integer, intent(in) :: SaltBCSurf, SaltBCBottom
!!$
!!$
!!$    ! 局所変数
!!$    ! Local variables
!!$    !
!!$    real(DP) :: xyz_RHSWork(0:iMax-1,jMax,0:kMax+1)
!!$
!!$    ! 実行文; Executable statement
!!$    !
!!$    
!!$    if(.not. isDifferentialMatInit) then
!!$       call MessageNotify('E', module_name, "HydroBouEqSolverVImplProc_Prepare has not been called.")
!!$    end if
!!$
!!$    !
!!$
!!$    !
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_Vor)
!!$    xyz_RHSWork(:,:,0) = xy_w(wa_VorBCRHS(:,1))
!!$    xyz_RHSWork(:,:,kMax) = xy_w(wa_VorBCRHS(:,2))
!!$
!!$    wz_Vor(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
!!$         & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
!!$         & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, .false. )
!!$
!!$    !
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_Div)
!!$    xyz_RHSWork(:,:,0) = xy_w(wa_DivBCRHS(:,1))
!!$    xyz_RHSWork(:,:,kMax) = xy_w(wa_DivBCRHS(:,2))
!!$    xyz_RHSWork(:,:,kMax+1) = 0d0
!!$
!!$    wz_Div(:,:) = solve( xyz_RHSWork(:,:,0:kMax+1), &
!!$         & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
!!$         & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, .true., xy_SurfPress )
!!$ !      wt_Div(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
!!$ !           & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
!!$ !           & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, .false. )
!!$
!!$
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_PTempEdd)
!!$    xyz_RHSWork(:,:,0) = xy_w(wa_PTempEddBCRHS(:,1))
!!$    xyz_RHSWork(:,:,kMax) = xy_w(wa_PTempEddBCRHS(:,2))
!!$
!!$    wz_PTempEdd(:,:) = solve( xyz_RHSWork(:,:,0:kMax), & 
!!$         & inquire_VBCSpecType(ThermBCSurf), inquire_VBCSpecType(ThermBCBottom), &
!!$         & vViscDiffTermCoef*vDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, .false. )
!!$
!!$
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_Salt)
!!$    xyz_RHSWork(:,:,0) = xy_w(wa_SaltBCRHS(:,1))
!!$    xyz_RHSWork(:,:,kMax) = xy_w(wa_SaltBCRHS(:,2))
!!$
!!$    wz_Salt(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
!!$         & inquire_VBCSpecType(SaltBCSurf), inquire_VBCSpecType(SaltBCBottom), &
!!$         & vViscDiffTermCoef*vDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, .false. )
!!$
!!$  end subroutine Advance_VImplicitProc

     
  !> @brief 
  !!
  !!
  subroutine Advance_VImplicitProc_DeltaForm(wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt, xy_SurfPress, &
       & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, &
       & xyz_VViscCoef, xyz_VDiffCoef, xy_totDepth, &
       & DynBCSurf, DynBCBottom, ThermBCSurf, ThermBCBottom, SaltBCSurf, SaltBCBottom )
    
    !
    !
    use Constants_mod, only: Cp0, RefDens

    use BoundaryCondO_mod, only: &
         & xy_SurfFwFlxO
    
use at_module_omp

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(out)  :: &
         & wz_DVor, wz_DDiv, wz_DPTempEdd, wz_DSalt
    real(dp), intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP), dimension(lMax,0:kMax), intent(inout) :: &
         & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_VViscCoef, xyz_VDiffCoef
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    integer, intent(in) :: DynBCSurf, DynBCBottom
    integer, intent(in) :: ThermBCSurf, ThermBCBottom
    integer, intent(in) :: SaltBCSurf, SaltBCBottom

    ! 局所変数
    ! Local variables
    !

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_UrfRHS, xyz_VrfRHS, xyz_fCoef
    real(DP) :: xyz_RHSWork(0:iMax-1,jMax,0:kMax+1)
    real(DP), dimension(2,0:iMax-1,jMax) :: axy_RobinBCUpCoef, axy_RobinBCBtmCoef
    
    ! 実行文; Executable statement
    !

    if(.not. isDifferentialMatInit) then
       call MessageNotify('E', module_name, "HydroBouEqSolverVImplProc_Prepare has not been called.")
    end if

!!$    !
    xyz_fCoef = CoriTermCoef*dt*(2d0*Omega*spread(sin(xyz_Lat(:,:,0)),3,kMax+1))
    
    call wz_VorDiv2VectorCosLat( wz_VorRHS, wz_DivRHS, &
         & xyz_UrfRHS, xyz_VrfRHS )

    call wz_VectorCosLat2VorDiv( &
         & (xyz_UrfRHS + xyz_fCoef*xyz_VrfRHS)/(1d0 + xyz_fCoef**2), &
         & (xyz_VrfRHS - xyz_fCoef*xyz_UrfRHS)/(1d0 + xyz_fCoef**2), &
         & wz_VorRHS, wz_DivRHS )

!!$    call invert_CoriTermOptr(wz_DivRHS, wz_VorRHS)    

    !
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(DynBCSurf))
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(DynBCBottom))

    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_VorRHS*dt)
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
    wz_DVor(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

    !
    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_DivRHS*dt)
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax:kMax+1) = 0d0;
    wz_DDiv(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

    
    !
    !
!!$    axy_RobinBCUpCoef(1,:,:) = 0d0!vDiffCoef/xy_totDepth
!!$    axy_RobinBCUpCoef(2,:,:) = 1d0!40d0/(RefDens*Cp0)
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(ThermBCSurf))        
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(ThermBCBottom))    
    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_PTempRHS*dt)
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;

    wz_DPTempEdd(:,:) = solve( xyz_RHSWork(:,:,0:kMax), & 
         & vViscDiffTermCoef*xyz_VDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef  )

    !
!!$    axy_RobinBCUpCoef(1,:,:) = vDiffCoef/xy_totDepth
!!$    axy_RobinBCUpCoef(2,:,:) = 0d0!(27.7d0*1d-3/86400d0)*35d0
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(SaltBCSurf))            
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(SaltBCBottom))        
    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_SaltRHS*dt)
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
    wz_DSalt(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

  end subroutine Advance_VImplicitProc_DeltaForm

  !> @brief 
  !!
  !!
  subroutine Advance_VImplicitProc_DeltaForm_xyz( &
       & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt, xy_SurfPress,                      &  ! (out)
       & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS,                               &  ! (inout)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_totDepth,                                   &  ! (in)
       & DynBCSurf, DynBCBottom, ThermBCSurf, ThermBCBottom, SaltBCSurf, SaltBCBottom &  ! (in)
       & )
    
    !
    !
    use Constants_mod, only: Cp0, RefDens

    use BoundaryCondO_mod, only: &
         & xy_SurfFwFlxO
    
use at_module_omp

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out)  :: &
         & xyz_DU, xyz_DV, xyz_DPTempEdd, xyz_DSalt
    real(dp), intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
         & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_VViscCoef, xyz_VDiffCoef
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    integer, intent(in) :: DynBCSurf, DynBCBottom
    integer, intent(in) :: ThermBCSurf, ThermBCBottom
    integer, intent(in) :: SaltBCSurf, SaltBCBottom

    ! 局所変数
    ! Local variables
    !

    real(DP), dimension(0:iMax-1,jMax) :: xy_fCoef
    real(DP) :: xyz_RHSWork(0:iMax-1,jMax,0:kMax+1), xy_URHSTmp(0:iMax-1,jMax)
    real(DP), dimension(2,0:iMax-1,jMax) :: axy_RobinBCUpCoef, axy_RobinBCBtmCoef
    integer :: k

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_CosLat
    real(DP), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_DVor, wz_DDiv
    ! 実行文; Executable statement
    !

    if(.not. isDifferentialMatInit) then
       call MessageNotify('E', module_name, "HydroBouEqSolverVImplProc_Prepare has not been called.")
    end if

    !
    !
    
    !
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(DynBCSurf))
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(DynBCBottom))

!!$    xyz_CosLat = cos(xyz_Lat)
!!$    call wz_VectorCosLat2VorDiv(xyz_URHS*xyz_CosLat, xyz_VRHS*xyz_CosLat, &
!!$         & wz_VorRHS ,wz_DivRHS)
!!$    
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_VorRHS*dt)
!!$    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
!!$    wz_DVor(:,:) = wz_xyz(solve_xyz( xyz_RHSWork(:,:,0:kMax), &
!!$         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
!!$         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef ))
!!$
!!$    !
!!$    xyz_RHSWork(:,:,0:kMax) = xyz_wz(wz_DivRHS*dt)
!!$    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax:kMax+1) = 0d0;
!!$    wz_DDiv(:,:) = wz_xyz(solve_xyz( xyz_RHSWork(:,:,0:kMax), &
!!$         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
!!$         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef ))
!!$
!!$    call wz_VorDiv2VectorCosLat(wz_DVor, wz_DDiv, &
!!$         & xyz_DU, xyz_DV)
!!$    xyz_DU = xyz_DU/xyz_CosLat; xyz_DV = xyz_DV/xyz_CosLat
    
!!$    
    xyz_RHSWork(:,:,0:kMax) = xyz_URHS*dt
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
    xyz_DU(:,:,:) = solve_xyz( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

    !
    xyz_RHSWork(:,:,0:kMax) = xyz_VRHS*dt
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax:kMax+1) = 0d0;
    xyz_DV(:,:,:) = solve_xyz( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VViscCoef, vViscDiffTermCoef*vHyperViscCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

    
!!$    if( CoriTermCoef > 0d0) then
!!$       
!!$       xy_fCoef(:,:) = CoriTermCoef*dt*( 2d0*Omega*sin(xyz_Lat(:,:,0)) )
!!$       !$omp parallel do private(xy_URHSTmp)
!!$       do k=0, kMax
!!$          xy_URHSTmp(:,:) = xyz_URHS(:,:,k)
!!$          xyz_URHS(:,:,k) = ( xyz_URHS(:,:,k) + xy_fCoef(:,:)*xyz_VRHS(:,:,k)  )/( 1d0 + xy_fCoef(:,:)**2 )
!!$          xyz_VRHS(:,:,k) = ( xyz_VRHS(:,:,k) - xy_fCoef(:,:)*xy_URHSTmp(:,:)  )/( 1d0 + xy_fCoef(:,:)**2 )       
!!$       end do
!!$    end if
    
    
    !
    !
!!$    axy_RobinBCUpCoef(1,:,:) = 0d0!vDiffCoef/xy_totDepth
!!$    axy_RobinBCUpCoef(2,:,:) = 1d0!40d0/(RefDens*Cp0)
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(ThermBCSurf))        
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(ThermBCBottom))    
    xyz_RHSWork(:,:,0:kMax) = xyz_PTempRHS*dt
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
    xyz_DPTempEdd(:,:,:) = solve_xyz( xyz_RHSWork(:,:,0:kMax), & 
         & vViscDiffTermCoef*xyz_VDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef  )

    !
!!$    axy_RobinBCUpCoef(1,:,:) = vDiffCoef/xy_totDepth
!!$    axy_RobinBCUpCoef(2,:,:) = 0d0!(27.7d0*1d-3/86400d0)*35d0
    call set_RobinBCCoef(axy_RobinBCUpCoef, inquire_VBCSpecType(SaltBCSurf))            
    call set_RobinBCCoef(axy_RobinBCBtmCoef, inquire_VBCSpecType(SaltBCBottom))        
    xyz_RHSWork(:,:,0:kMax) = xyz_SaltRHS*dt
    xyz_RHSWork(:,:,0) = 0d0; xyz_RHSWork(:,:,kMax) = 0d0;
    xyz_DSalt(:,:,:) = solve_xyz( xyz_RHSWork(:,:,0:kMax), &
         & vViscDiffTermCoef*xyz_VDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, xy_totDepth, &
         & .false., axy_RobinBCUpCoef, axy_RobinBCBtmCoef )

  end subroutine Advance_VImplicitProc_DeltaForm_xyz
  
  subroutine set_RobinBCCoef(axy_RobinBCCoef, BCType)
    real(DP), intent(inout) :: axy_RobinBCCoef(2,0:iMax-1,jMax)
    character :: BCType

    select case(BCType)
    case ('D')
       axy_RobinBCCoef(1,:,:) = 0d0;  axy_RobinBCCoef(2,:,:) = 1d0;
    case ('N')
       axy_RobinBCCoef(1,:,:) = 1d0;  axy_RobinBCCoef(2,:,:) = 0d0;
    end select
  end subroutine set_RobinBCCoef
  
  function solve( &
       & xyz_RHS, xyz_vDiffTermCoef, vHyperDiffTermCoef, xy_totDepth, isDivEq, &
       & axy_RobinBCUpCoef, axy_RobinBCBtmCoef, xy_SurfPress) result(wz_ret)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_RHS(0:,:,0:)
    real(DP), intent(in) :: xyz_vDiffTermCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: vHyperDiffTermCoef
    logical, intent(in) :: isDivEq
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in), dimension(2,0:iMax-1,jMax) :: axy_RobinBCUpCoef, axy_RobinBCBtmCoef
    real(DP), optional, intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP) :: wz_ret(lMax,0:kMax)


    ! 局所変数
    ! Local variables
    !
    real(DP) :: gz_VImplMat(0:kMax+1,0:kMax+1)
    real(DP) :: xyz_Tmp(0:iMax-1,jMax,0:size(xyz_RHS,3)-1)
    integer :: i, j

    ! 実行文; Executable statement
    !

    !$omp parallel do private(i, gz_VImplMat)
    do j=1,jMax
       do i=0,iMax-1
          call constrct_VImplMat(gz_VImplMat,                                &  ! (out)
               & dt, xy_totDepth(i,j), xyz_vDiffTermCoef(i,j,:), vHyperDiffCoef, & ! (in)
               & isDivEq, axy_RobinBCUpCoef(:,i,j), axy_RobinBCBtmCoef(:,i,j)    & ! (in)
               & )

          xyz_Tmp(i,j,:) = solve_LinearEq(gz_VImplMat(0:size(xyz_RHS,3)-1,0:size(xyz_RHS,3)-1), xyz_RHS(i,j,:))
       end do
    end do

    wz_ret(:,:) = wz_xyz(xyz_Tmp(:,:,0:tMax))

  end function solve

  function solve_xyz( &
       & xyz_RHS, xyz_vDiffTermCoef, vHyperDiffTermCoef, xy_totDepth, isDivEq, &
       & axy_RobinBCUpCoef, axy_RobinBCBtmCoef, xy_SurfPress) result(xyz_ret)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_RHS(0:,:,0:)
    real(DP), intent(in) :: xyz_vDiffTermCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: vHyperDiffTermCoef
    logical, intent(in) :: isDivEq
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in), dimension(2,0:iMax-1,jMax) :: axy_RobinBCUpCoef, axy_RobinBCBtmCoef
    real(DP), optional, intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP) :: xyz_ret(0:iMax-1,jMax,0:kMax)


    ! 局所変数
    ! Local variables
    !
    real(DP) :: gz_VImplMat(0:kMax+1,0:kMax+1)
    integer :: i, j

    ! 実行文; Executable statement
    !

    !$omp parallel do private(i, gz_VImplMat)
    do j=1,jMax
       do i=0,iMax-1
          call constrct_VImplMat(gz_VImplMat,                                    &  ! (out)
               & dt, xy_totDepth(i,j), xyz_vDiffTermCoef(i,j,:), vHyperDiffCoef, &  ! (in)
               & isDivEq, axy_RobinBCUpCoef(:,i,j), axy_RobinBCBtmCoef(:,i,j)    &  ! (in)
               & )

          xyz_ret(i,j,:) = solve_LinearEq(gz_VImplMat(0:size(xyz_RHS,3)-1,0:size(xyz_RHS,3)-1), xyz_RHS(i,j,:))
       end do
    end do

  end function solve_xyz
  
  function solve_LinearEq(A, b) result(x)

    use lumatrix

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: A(:,:)
    real(DP), intent(in) :: b(size(A,1))
    real(DP) :: x(size(A,1))

    ! 局所変数
    ! Local variables
    !
    integer :: N
    integer :: IPIV(size(A,1))
    integer :: INFO

!!$    integer :: kp(size(A,1))

    ! 実行文; Executable statement
    !

    N = size(A,1)
    x(:) = b
    call DGESV(N, 1, A, N, IPIV, x, N, INFO)

    !
!!$    call ludecomp(A, kp)
!!$    x = lusolve(A, kp, b)

  end function solve_LinearEq

  subroutine constrct_VImplMat(gt_VImplMat, &
       & dt, totDepth, z_vDiffTermCoef, vHyperDiffTermCoef, &
       & isDivEq, RobinBCUpCoef, RobinBCBtmCoef )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: gt_VImplMat(0:kMax+1,0:tMax+1)
    real(DP), intent(in) :: dt, totDepth, z_vDiffTermCoef(0:kMax), vHyperDiffTermCoef
    logical, intent(in) :: isDivEq
    real(DP), optional, intent(in) :: RobinBCUpCoef(2), RobinBCBtmCoef(2)

    ! 局所変数
    ! Local variables
    !
    integer :: k, t

    ! 実行文; Executable statement
    !

    gt_VImplMat(:,:) = 0d0
    do k=1,kMax-1
       do t=0,tMax
          gt_VImplMat(k,t) = &
               &   gt_Mat(t,k) &
               & - dt/totDepth**2 * dot_product(gt_DSig1Mat(t,:),z_vDiffTermCoef(:)*gt_DSig1Mat(:,k)) &
               & + dt*vHyperDiffTermCoef/totDepth**4 * gt_DSig4Mat(t,k)
       end do
    end do

    gt_VImplMat(0,0:kMax) = &
         & RobinBCUpCoef(1)*gt_DSig1Mat(:,0) + RobinBCUpCoef(2)*gt_Mat(:,0)

    gt_VImplMat(kMax,0:kMax) = &
         & RobinBCBtmCoef(1)*gt_DSig1Mat(:,kMax) + RobinBCBtmCoef(2)*gt_Mat(:,kMax)
    
    if(isDivEq) then
       gt_VImplMat(1:kMax-1,tMax+1) = dt
       gt_VImplMat(kMax+1,0:tMax) = t_IntSigMat
    end if

  end subroutine constrct_VImplMat

  subroutine invert_CoriTermOptr(wz_DivRHS, wz_VorRHS)
    real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_DivRHS, wz_VorRHS
#ifdef DSOGCM_MODE_AXISYM
    complex(DP), dimension(0:0,-2:nMax+2) :: A, B, C
#else
    complex(DP), dimension(0:nMax,-2:nMax+2) :: A, B, C
#endif
    integer :: k, m, MM


    !
    !
    
#ifdef DSOGCM_MODE_AXISYM
     MM = 0 
#else
     MM = nMax
#endif

    call construct_coef()
    do k=0,kMax
       do m=0, MM
          call solve_for_m(m, wz_VorRHS(:,k), wz_DivRHS(:,k))
       end do
    end do
    
  contains
    subroutine construct_coef()
      integer :: m, n
      real(DP) :: tau
      
      tau = CoriTermCoef*2d0*Omega
      A(:,:) = 1d0; B(:,:) = 0d0; C(:,:) = 0d0
      do m=0,MM
         do n=m,nMax
            A(m,n) = cmplx(1d0/dt, -m*tau/dble(n*(n + 1)))
            B(m,n) = cmplx(tau*(n + 1)/dble(n)*e_(m,n), 0d0)
            C(m,n) = cmplx(tau*n/dble(n+1)*e_(m,n+1), 0d0)
         end do
      end do
      A(0,0) = 1d0/dt; B(0,0) = 0d0
      
    end subroutine construct_coef

    subroutine solve_for_m(m, w_VorRHS, w_DivRHS)
      integer, intent(in) :: m
      real(DP), intent(inout) :: w_VorRHS(lMax), w_DivRHS(lMax)
      
      complex(DP) :: Div(m-1:nMax+1), RHS(m:nMax), Zeta(m:nMax)
      complex(DP) :: CoefMat(m:nMax,m-2:nMax+2)
      integer :: n
      integer :: info, ipiv(nMax-m+1)
      
      CoefMat(:,:) = 0d0
      Div(:) = 0d0
      
      do n=m, nMax
         if (n==m) then
            CoefMat(n,n-2) = 0d0
            CoefMat(n,n) = A(m,n) + B(m,n+1)*C(m,n)/A(m,n+1)
            CoefMat(n,n+2) = C(m,n)*C(m,n+1)/A(m,n+1)
            Div(n) =  w2fc(w_DivRHS,m,n) + C(m,n)/A(m,n+1)*w2fc(w_VorRHS,m,n+1)
         else if(n==nMax) then
            CoefMat(n,n-2) = B(m,n-1)*B(m,n)/A(m,n-1)
            CoefMat(n,n) = A(m,n) + B(m,n)*C(m,n-1)/A(m,n-1)
            CoefMat(n,n+2) = 0d0!C(m,n)*C(m,n+1)/A(m,n+1)
            Div(n) =  w2fc(w_DivRHS,m,n) + B(m,n)/A(m,n-1)*w2fc(w_VorRHS,m,n-1)
         else
            CoefMat(n,n-2) = B(m,n-1)*B(m,n)/A(m,n-1)
            CoefMat(n,n) = A(m,n) + B(m,n)*C(m,n-1)/A(m,n-1) + B(m,n+1)*C(m,n)/A(m,n+1)
            CoefMat(n,n+2) = C(m,n)*C(m,n+1)/A(m,n+1)
            Div(n) =  w2fc(w_DivRHS,m,n) &
              &  + B(m,n)/A(m,n-1)*w2fc(w_VorRHS,m,n-1) + C(m,n)/A(m,n+1)*w2fc(w_VorRHS,m,n+1)
         end if
      end do

!!$      write(*,*) Div(m:nMax)
      call zgesv(nMax-m+1, 1, CoefMat(m:nMax,m:nMax), nMax-m+1, ipiv, Div(m:nMax), nMax-m+1, info)
!!$      write(*,*) "->", Div(m:nMax)
!!$      write(*,*) "--------"

      do n=m, nMax
         Zeta(n) = (w2fc(w_VorRHS,m,n) - B(m,n)*Div(n-1) - C(m,n)*Div(n+1))/A(m,n)
      end do


      !
      !      write(*,*) w_VorRHS
      
      call fc2w(Div(m:nMax)/dt, w_DivRHS, m)
      call fc2w(Zeta(m:nMax)/dt, w_VorRHS, m)
!      write(*,*) "->>", w_VorRHS
!!$      stop
    end subroutine solve_for_m

    complex function w2fc(w,m,n)
      real(DP), intent(in) :: w(lMax)
      integer, intent(in) :: m, n

      w2fc = cmplx( w(l_nm(n,m)), w(l_nm(n,-m)) )
    end function w2fc

    subroutine fc2w(fc, w, m)
      complex(DP), intent(in) :: fc(m:nMax)
      real(DP), intent(inout) :: w(lMax)
      integer, intent(in) :: m
      integer :: n

      do n=m, nMax
         w(l_nm(n,m)) = real(fc(n))
         w(l_nm(n,-m)) = aimag(fc(n))
      end do
    end subroutine fc2w

    real(DP) function  e_(m,n)
      integer, intent(in) :: m, n
      e_ = sqrt((n**2 - m**2)/(4d0*n**2 - 1d0))
    end function e_

  end subroutine invert_CoriTermOptr


end module HydroBouEqSolverVImplProc_mod

