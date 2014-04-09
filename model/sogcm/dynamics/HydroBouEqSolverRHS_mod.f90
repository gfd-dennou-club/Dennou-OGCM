!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolverRHS_mod 

  ! モジュール引用; Use statements
  !

  use dc_types, only: &
       & DP, STRING 

  use Constants_mod, only: &
       & Omega, Grav, RPlanet, RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final
  public :: calc_VorEqDivEqInvisRHS, calc_VorEqDivEqDiffRHS
  public :: calc_TracerEqInvisRHS, calc_TracerEqDiffRHS
  public :: calc_SurfHeightRHS
  public :: correct_DivEqRHSUnderRigidLid

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverRHS_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverRHS_Init()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolverRHS_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolverRHS_Final()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolverRHS_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine calc_VorEqDivEqInvisRHS(wz_RHSVor, wz_RHSDiv, &
       & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_PEdd, xyz_GeoPot, xyz_SigDot &
       & )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: wz_RHSVor(lMax, 0:kMax)
    real(DP), intent(out) :: wz_RHSDiv(lMax, 0:kMax)
    real(DP), intent(in) :: xyz_Vor(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_SurfHeight(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_PEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_A(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_B(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_KinEngy(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_AbsVor(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_GeoPotGradCoef(0:iMax-1,jMax,0:kMax)
    real(DP) :: wz_GeoPot(lMax, 0:kMax)

    ! 実行文; Executable statement
    !

    !$omp parallel

    !$omp workshare
    xyz_KinEngy = (xyz_Urf**2 + xyz_Vrf**2)/(2d0*cos(xyz_Lat)**2)
    !$omp end workshare

    !$omp workshare
    xyz_AbsVor = xyz_Vor + 2d0*Omega*sin(xyz_Lat)
    !$omp end workshare

    !$omp workshare
    xyz_GeoPotGradCoef = xyz_DensEdd/(RefDens*RPlanet)
    !$omp end workshare

    !$omp end parallel

    wz_GeoPot = wz_xyz(xyz_GeoPot)

    xyz_A = &
         &   xyz_AbsVor*xyz_Urf + xyz_SigDot*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Vrf))) &
         & + xyz_GeoPotGradCoef*xya_GradMu_wa(wz_GeoPot)
    xyz_B = &
         &   xyz_AbsVor*xyz_Vrf - xyz_SigDot*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Urf))) &
         & - xyz_GeoPotGradCoef*xya_GradLambda_wa(wz_GeoPot)

    wz_RHSVor = - wz_AlphaOptr_xyz( xyz_A, xyz_B )

    wz_RHSDiv = wz_AlphaOptr_xyz( xyz_B, -xyz_A ) &
         &      - wz_Lapla2D_wz(wz_xyz( &
         &             xyz_KinEngy + Grav*spread(xy_SurfHeight,3,kMax+1) &
         &           + xyz_PEdd/RefDens &
         &        ))

  end subroutine calc_VorEqDivEqInvisRHS

  !> @brief 
  !!
  !!
  subroutine calc_VorEqDivEqDiffRHS(wz_RHSVor, wz_RHSDiv, &
       & wz_Vor, wz_Div,&
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & xyz_totDepth, &
       & overwrite )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSVor(lMax, 0:kMax)
    real(DP), intent(inout) :: wz_RHSDiv(lMax, 0:kMax)
    real(DP), intent(in) :: wz_Vor(lMax,0:kMax)
    real(DP), intent(in) :: wz_Div(lMax,0:kMax)
    real(DP), intent(in) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: hViscCoef, vViscCoef
    real(DP), intent(in) :: hHyperViscCoef, vHyperViscCoef
    logical, optional, intent(in) :: overWrite
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyt_VorDIVDep2(0:iMax-1, jMax, 0:tMax)
    real(DP) :: xyt_DivDIVDep2(0:iMax-1, jMax, 0:tMax)

    ! 実行文; Executable statement
    !
    
    if( present(overwrite) .and. overwrite ) then
       wz_RHSVor = 0d0; wz_RHSDiv = 0d0
    end if

    xyt_VorDIVDep2 = xyt_xyz(xyz_wz(wz_Vor))/xyz_totDepth**2
    wz_RHSVor = wz_RHSVor &
!!$         &      + hViscCoef*wz_Vor/RPlanet**2 &
         &      + wz_Lapla2D_wz(hViscCoef*wz_Vor - hHyperViscCoef*wz_Lapla2D_wz(wz_Vor))  &
         &      + wz_xyz(xyz_xyt( &
         &          xyt_DSigDSig_xyt( vViscCoef*xyt_VorDIVDep2 - vHyperViscCoef*xyt_DSigDSig_xyt(xyt_VorDIVDep2/xyz_totDepth**2) ) &
         &      ))

    xyt_DivDIVDep2 = xyt_xyz(xyz_wz(wz_Div))/xyz_totDepth**2
    wz_RHSDiv = wz_RHSDiv &
!!$         &      + hViscCoef*wz_Div/RPlanet**2 &
         &      + wz_Lapla2D_wz(hViscCoef*wz_Div - hHyperViscCoef*wz_Lapla2D_wz(wz_Div))  &
         &      + wz_xyz(xyz_xyt( &
         &          xyt_DSigDSig_xyt( vViscCoef*xyt_DivDIVDep2 - vHyperViscCoef*xyt_DSigDSig_xyt(xyt_DivDIVDep2/xyz_totDepth**2) ) &
         &      ))


  end subroutine calc_VorEqDivEqDiffRHS

  subroutine calc_TracerEqInvisRHS( wz_RHSTracer, xyz_Tracer, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    real(DP), intent(out) :: wz_RHSTracer(lMax, 0:kMax)
    real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    wz_RHSTracer =  &
         & - wz_AlphaOptr_xyz(xyz_Tracer*xyz_Urf, xyz_Tracer*xyz_Vrf) &
         & + wz_xyz(   xyz_Tracer*xyz_Div & 
         &           - xyz_SigDot*xyz_wt(wt_DSig_wt(wt_xyz(xyz_Tracer))) )

!!$    wz_RHSTracer =  &
!!$         & - wz_AlphaOptr_xyz(xyz_Tracer*xyz_Urf, xyz_Tracer*xyz_Vrf) &
!!$         & + wz_xyz(   xyz_Tracer*xyz_Div &
!!$         &           - xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot*xyz_Tracer))) &
!!$         &           + xyz_Tracer*xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot))) )

  end subroutine calc_TracerEqInvisRHS

  !> @brief 
  !!
  !!
  subroutine calc_TracerEqDiffRHS(wz_RHSTracer, &
       & wz_Tracer, Kh, Kv, xyz_totDepth, &
       & overwrite )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSTracer(lMax, 0:kMax)
    real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
    real(DP), intent(in) :: Kh, Kv
    real(DP), intent(in) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    logical, optional, intent(in) :: overWrite
    
    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !
    
    if( present(overwrite) .and. overwrite ) then
       wz_RHSTracer = 0d0
    end if

    wz_RHSTracer = wz_RHSTracer & 
         & + Kh * wz_Lapla2D_wz(wz_Tracer) &
         & + Kv * wz_wt( &
         &          wt_DSig_wt(wt_DSig_wt( wt_xyz(xyz_wz(wz_Tracer)/xyz_totDepth**2) )) &
         &        )

  end subroutine calc_TracerEqDiffRHS


  subroutine calc_SurfHeightRHS( w_RHSSurfHeight, xyz_Urf, xyz_Vrf, xy_totDepth )
    real(DP), intent(out) :: w_RHSSurfHeight(lMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)

    w_RHSSurfHeight =  0d0!- w_AlphaOptr_xy( &
!         & xy_totDepth*xy_IntSig_BtmToTop_xyz(xyz_Urf), xy_totDepth*xy_IntSig_BtmToTop_xyz(xyz_Vrf) )

  end subroutine calc_SurfHeightRHS


  !> @brief 
  !!
  !!
  subroutine correct_DivEqRHSUnderRigidLid(wz_RHSDivEqN, &
       & xy_SurfPress, wz_DivN, xy_totDepth, Av, dt )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSDivEqN(lMax, 0:kMax)
    real(DP), intent(inout) :: xy_SurfPress(0:iMax-1,jMax)
    real(DP), intent(in) :: wz_DivN(lMax,0:kMax), xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: Av, dt

    ! 局所変数
    ! Local variables
    !
    real(DP) :: w_CorrectTerm(lMax)
    real(DP) :: xyz_tmp(0:iMax-1,jMax,0:kMax)
    integer :: k
    real(DP) :: wz_Av_DSig_Div(lMax,0:kMax)
    
    real(DP) :: xyz_Av_DSig_Div(0:iMax-1,jMax,0:kMax)
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Div, xyz_RHSDiv
    real(DP) :: xy_Correct(0:iMax-1,jMax)
    real(DP) :: xy_DivSig(0:iMax-1,jMax), xy_RHSDivSig(0:iMax-1,jMax)

    ! 実行文; Executable statement
    !

!!$    wz_Av_Dsig_Div = Av*wz_wt( &
!!$         & wt_DSig_wt(wt_xyz(xyz_wz(wz_DivN)/spread(xy_totDepth**2,3,kMax+1))) )
    w_CorrectTerm = - w_IntSig_BtmToTop_wz( wz_RHSDivEqN + wz_DivN/dt ) !&
!         &          + wz_Av_DSig_Div(:,kMax)
!write(*,*) "Av*DSig_Div:", wz_Av_DSig_Div(:,kMax)


    xy_SurfPress = - xy_w( w_InvLapla2D_w( w_CorrectTerm*RefDens ) )

!!$    xyz_Av_DSig_Div = xyz_wz(wz_Av_DSig_Div)
!!$    xy_Correct = xy_w(w_CorrectTerm)
!!$    xyz_Div = xyz_wz(wz_DivN)
!!$    xyz_RHSDiv = xyz_wz(wz_RHSDivEqN)
!!$    xy_DivSig = xy_IntSig_BtmToTop_xyz(xyz_Div)
!!$    xy_RHSDivSig = xy_IntSig_BtmToTop_xyz(xyz_wz(wz_RHSDivEqN))
!!$write(*,'(13(1x,e12.5),a,13(1x,e12.5),a,e12.5)') xyz_RHSDiv(1,16,0:12), "*", xyz_Div(1,16,0:12)/dt, &
!!$     & "*", -xyz_Av_DSig_Div(1,16,kMax) !xy_DivSig(1,:)
!!$write(*,'(3(1x,e12.5))') xy_Correct(1,16), xy_RHSDivSig(1,16), xy_DivSig(1,16)/dt
!!$write(*,*) '-----'

    !$omp parallel do
    do k=0,kMax
       wz_RHSDivEqN(:,k) = wz_RHSDivEqN(:,k) + w_CorrectTerm
    end do

  end subroutine correct_DivEqRHSUnderRigidLid

end module HydroBouEqSolverRHS_mod

