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

  use GridSet_mod, only: &
       & iMax, jMax, kMax, tMax, lMax, &
       & xyz_Lat

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType

  use Constants_mod, only: &
       & refDens, &
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
  public :: Advance_VImplicitProc

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

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverVImplProc_Init()

    ! 実行文; Executable statements
    !

    isDifferentialMatInit = .false.

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

  subroutine HydroBouEqSolverVImplProc_Prepare()
    
    use at_module_omp

    real(DP), dimension(0:tMax,0:tMax) :: tt_I
    integer :: t

    if(isDifferentialMatInit) then
       deallocate(gt_Mat, gt_DSig1Mat, gt_DSig2Mat, gt_DSig4Mat, t_IntSigMat)
    end if

    allocate(gt_Mat(0:kMax,0:tMax), gt_DSig1Mat(0:kMax,0:tMax), gt_DSig2Mat(0:kMax,0:tMax), gt_DSig4Mat(0:kMax,0:tMax), &
         &   t_IntSigMat(0:tMax))

    tt_I = 0d0
    do t=0, tMax
       tt_I(t,t) = 1d0
    end do
    gt_Mat = ag_at(tt_I)
    gt_DSig1Mat = ag_at( at_Dx_at(tt_I) )
    gt_DSig2Mat = ag_at( at_Dx_at(at_Dx_at(tt_I)) )
    gt_DSig4Mat = ag_at( at_Dx_at(at_Dx_at(at_Dx_at(at_Dx_at(tt_I)))) )
    forAll(t=0:tMax) t_IntSigMat(t) = dot_product(g_X_WEIGHT, gt_Mat(t,:))

    isDifferentialMatInit = .true.

  end subroutine HydroBouEqSolverVImplProc_Prepare

  !> @brief 
  !!
  !!
  subroutine Advance_VImplicitProc(wt_Vor, wt_Div, wt_PTempEdd, wt_Salt, &
    & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS, xy_totDepth, &
    & vViscDiffTermCoef, dt,  &
    & DynBCSurf, DynBCBottom, ThermBCSurf, ThermBCBottom, SaltBCSurf, SaltBCBottom )
    
    !
    !

use at_module_omp

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:tMax), intent(inout), optional :: &
         & wt_Vor, wt_Div, wt_PTempEdd, wt_Salt
    real(DP), dimension(lMax,2), intent(in), optional :: &
         & wa_VorBCRHS, wa_DivBCRHS, wa_PTempEddBCRHS, wa_SaltBCRHS
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: vViscDiffTermCoef
    real(DP), intent(in) :: dt
    integer, intent(in), optional :: DynBCSurf, DynBCBottom
    integer, intent(in), optional :: ThermBCSurf, ThermBCBottom
    integer, intent(in), optional :: SaltBCSurf, SaltBCBottom


    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_RHSWork(0:iMax-1,jMax,0:kMax+1)


    ! 実行文; Executable statement
    !

    !
    !
    if( present(wt_Vor) ) then

       xyz_RHSWork(:,:,0:kMax) = xyz_wt(wt_Vor)
       xyz_RHSWork(:,:,0) = xy_w(wa_VorBCRHS(:,1))
       xyz_RHSWork(:,:,kMax) = xy_w(wa_VorBCRHS(:,2))

       wt_Vor(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
            & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
            & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, .false. )
    end if

    !
    if( present(wt_Div) ) then

       xyz_RHSWork(:,:,0:kMax) = xyz_wt(wt_Div)
       xyz_RHSWork(:,:,0) = xy_w(wa_DivBCRHS(:,1))
       xyz_RHSWork(:,:,kMax) = xy_w(wa_DivBCRHS(:,2))
       xyz_RHSWork(:,:,kMax+1) = 0d0

       wt_Div(:,:) = solve( xyz_RHSWork(:,:,0:kMax+1), &
            & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
            & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, .true. )
!!$       wt_Div(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
!!$            & inquire_VBCSpecType(DynBCSurf), inquire_VBCSpecType(DynBCBottom), &
!!$            & vViscDiffTermCoef*vViscCoef, vViscDiffTermCoef*vHyperViscCoef, .false. )

    end if

    if( present(wt_PTempEdd) ) then

       xyz_RHSWork(:,:,0:kMax) = xyz_wt(wt_PTempEdd)
       xyz_RHSWork(:,:,0) = xy_w(wa_PTempEddBCRHS(:,1))
       xyz_RHSWork(:,:,kMax) = xy_w(wa_PTempEddBCRHS(:,2))

       wt_PTempEdd(:,:) = solve( xyz_RHSWork(:,:,0:kMax), & 
            & inquire_VBCSpecType(ThermBCSurf), inquire_VBCSpecType(ThermBCBottom), &
            & vViscDiffTermCoef*vDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, .false. )
    end if

    if( present(wt_Salt) ) then

       xyz_RHSWork(:,:,0:kMax) = xyz_wt(wt_Salt)
       xyz_RHSWork(:,:,0) = xy_w(wa_SaltBCRHS(:,1))
       xyz_RHSWork(:,:,kMax) = xy_w(wa_SaltBCRHS(:,2))

       wt_Salt(:,:) = solve( xyz_RHSWork(:,:,0:kMax), &
            & inquire_VBCSpecType(SaltBCSurf), inquire_VBCSpecType(SaltBCBottom), &
            & vViscDiffTermCoef*vDiffCoef, vViscDiffTermCoef*vHyperDiffCoef, .false. )
    end if

    contains

      function solve(xyz_RHS, BCUpper, BCBottom, vDiffTermCoef, vHyperDiffTermCoef, isDivEq) result(wt_ret)

        ! 宣言文; Declaration statement
        !
        real(DP), intent(in) :: xyz_RHS(0:,:,0:)
        character, intent(in) :: BCUpper, BCBottom
        real(DP), intent(in) :: vDiffTermCoef, vHyperDiffTermCoef
        logical, intent(in) :: isDivEq
        real(DP) :: wt_ret(lMax,0:tMax)

        ! 局所変数
        ! Local variables
        !
        real(DP) :: gt_VImplMat(0:kMax+1,0:tMax+1)
        real(DP) :: xyt_Tmp(0:iMax-1,jMax,0:size(xyz_RHS,3)-1)
        integer :: i, j

        ! 実行文; Executable statement
        !

        !$omp parallel do private(i, gt_VImplMat)
        do j=1,jMax
           do i=0,iMax-1
              call constrct_VImplMat(gt_VImplMat,                                &  ! (out)
                   & dt, xy_totDepth(i,j), vDiffTermCoef, vHyperDiffCoef, BCUpper, BCBottom, isDivEq  & !(in)
                   & )

              xyt_Tmp(i,j,:) = solve_LinearEq(gt_VImplMat(0:size(xyz_RHS,3)-1,0:size(xyz_RHS,3)-1), xyz_RHS(i,j,:))
           end do
        end do

        wt_ret(:,:) = wt_xyt(xyt_Tmp(:,:,0:tMax))

      end function solve

    end subroutine Advance_VImplicitProc

  function solve_LinearEq(A, b) result(x)

    use lumatrix

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: A(:,:)
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
       & dt, totDepth, vDiffTermCoef, vHyperDiffTermCoef, BCKindUpper, BCKindBottom, isDivEq )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: gt_VImplMat(0:kMax+1,0:tMax+1)
    real(DP), intent(in) :: dt, totDepth, vDiffTermCoef, vHyperDiffTermCoef
    character, intent(in) :: BCKindUpper, BCKindBottom
    logical, intent(in) :: isDivEq

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
               & - dt*vDiffTermCoef/totDepth**2 * gt_DSig2Mat(t,k) &
               & + dt*vHyperDiffTermCoef/totDepth**4 * gt_DSig4Mat(t,k)
       end do
    end do

    select case(BCKindUpper)
       case ('N') 
          gt_VImplMat(0,0:tMax) = gt_DSig1Mat(:,0)
       case ('D')
          gt_VImplMat(0,0:tMax) = gt_Mat(:,0)
    end select

    select case(BCKindBottom)
       case ('N') 
          gt_VImplMat(kMax,0:tMax) = gt_DSig1Mat(:,kMax)
       case ('D')
          gt_VImplMat(kMax,0:tMax) = gt_Mat(:,kMax)
    end select
    
    if(isDivEq) then
       gt_VImplMat(1:kMax-1,tMax+1) = 0.5d0*dt
       gt_VImplMat(kMax+1,:) = t_IntSigMat
    end if

  end subroutine constrct_VImplMat

end module HydroBouEqSolverVImplProc_mod

