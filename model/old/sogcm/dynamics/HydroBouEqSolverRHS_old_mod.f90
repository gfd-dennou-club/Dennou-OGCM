!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief This module provides the subroutines to calculate the tendency of each time evolution equation 
!!that hydrostaic boussinesq equation consist of.  
!!
!! @attention This module require that the initialization of EOSDriver_mod, DiagnoseUtil_mod has been finished.  
!! 
!! @author Yuta Kawai
!! @since 2013
!!
!!
module HydroBouEqSolverRHS_old_mod 

  ! モジュール引用; Use statements
  !

  !* gtool 
  
  use dc_types, only: &
       & DP, STRING 

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use Constants_mod, only: &
       & Omega, Grav, RPlanet, RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use EOSDriver_mod, only: &
       & EOSDriver_Eval
  
  use DiagnoseUtil_mod, only: &
       & Diagnose_SigDot, Diagnose_HydroPressEdd, Diagnose_GeoPot

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  ! Driver routines
  public :: HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final
  public :: calc_HydroBouEqInvisRHS
  public :: calc_HydroBouEqHViscRHS
  public :: calc_HydroBouEqVViscRHS, calc_HydroBouEqVViscRHS_xyz

  ! Calculation routines
  public :: calc_VorEqDivEqInvisRHS, calc_TracerEqInvisRHS
  public :: calc_HDiffRHS
  public :: calc_VDiffRHS_xyz
  public :: calc_SurfHeightRHS

  public :: correct_DivEqRHSUnderRigidLid, correct_DivEqRHSUnderRigidLid2
  public :: correct_DivVorEqRHSUnderRigidLid
  
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
!
! * Driver routines in which some RHS calculation routines are called. 
!   The called routines calculate the tendency of each process such as advection, diffusion and etc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> @brief 
  !!
  !!
  subroutine calc_HydroBouEqInvisRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &
       & xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt, xy_SurfHeight, &
       & xy_SurfPress, xy_totDepthBasic, z_PTempBasic, &
       & xyz_UrfCori, xyz_VrfCori )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(out), dimension(lMax) :: w_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight, xy_SurfPress, xy_totDepthBasic
    real(DP), intent(in), dimension(0:kMax) :: z_PTempBasic
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_UrfCori, xyz_VrfCori

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_SigDot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_GeoPot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_Press(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)    
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)
    integer :: i, j, k

    ! 実行文; Executable statement
    !

    xy_totDepth(:,:)   = xy_totDepthBasic! + xy_SurfHeightN
    xyz_SigDot(:,:,:)  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )
    xyz_GeoPot(:,:,:)  = Diagnose_GeoPot( xy_totDepth )

    forAll(k=0:kMax) &
         & xyz_PTemp(:,:,k) = xyz_PTempEdd(:,:,k) + z_PTempBasic(k)

    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)

    ! Calculate the pressure which is deviation from -RefDens*Grav*z).
    !
    xyz_Press(:,:,:) =   spread(xy_SurfPress, 3, kMax+1)                 &    ! barotropic component
                &      + Diagnose_HydroPressEdd( xy_totDepth, xyz_DensEdd )      ! baroclinic component

    !
    !
    call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot, &
         & xyz_UrfCori, xyz_VrfCori )

    Call calc_TracerEqInvisRHS(wz_PTempRHS, &
         & xyz_PTemp, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )

    call calc_TracerEqInvisRHS(wz_SaltRHS, &
         & xyz_Salt, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    
    call calc_SurfHeightRHS(w_SurfHeightRHS, &
         & xyz_Urf, xyz_Vrf, xy_totDepth ) 
    
  end subroutine calc_HydroBouEqInvisRHS

  !> @brief 
  !!
  !!
  subroutine calc_HydroBouEqInvisRHS_xyz( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, w_SurfHeightRHS, &
       & xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt, xy_SurfHeight, &
       & xy_SurfPress, xy_totDepthBasic, z_PTempBasic, &
       & xyz_UrfCori, xyz_VrfCori )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(out), dimension(lMax) :: w_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight, xy_SurfPress, xy_totDepthBasic
    real(DP), intent(in), dimension(0:kMax) :: z_PTempBasic
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_UrfCori, xyz_VrfCori

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_SigDot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_GeoPot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_Press(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)    
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)
    integer :: i, j, k

    ! 実行文; Executable statement
    !

    xy_totDepth(:,:)   = xy_totDepthBasic! + xy_SurfHeightN
    xyz_SigDot(:,:,:)  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )
    xyz_GeoPot(:,:,:)  = Diagnose_GeoPot( xy_totDepth )

    forAll(k=0:kMax) &
         & xyz_PTemp(:,:,k) = xyz_PTempEdd(:,:,k) + z_PTempBasic(k)

    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)

    ! Calculate the pressure which is deviation from -RefDens*Grav*z).
    !
    xyz_Press(:,:,:) =   spread(xy_SurfPress, 3, kMax+1)                 &    ! barotropic component
                &      + Diagnose_HydroPressEdd( xy_totDepth, xyz_DensEdd )      ! baroclinic component

    !
    !
    call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot, &
         & xyz_UrfCori, xyz_VrfCori )

    Call calc_TracerEqInvisRHS(wz_PTempRHS, &
         & xyz_PTemp, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )

    call calc_TracerEqInvisRHS(wz_SaltRHS, &
         & xyz_Salt, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    
    call calc_SurfHeightRHS(w_SurfHeightRHS, &
         & xyz_Urf, xyz_Vrf, xy_totDepth ) 
    
  end subroutine calc_HydroBouEqInvisRHS_xyz
  
  !> @brief 
  !!
  !!
  subroutine calc_HydroBouEqHViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, &
       & wz_Vor, wz_Div, wz_PTempEdd, wz_Salt, &
       & hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef, hHyperDiffTermCoef, &
       & isRHSReplace )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTempEdd, wz_Salt
    real(DP), intent(in) :: hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef, hHyperDiffTermCoef
    logical, intent(in) :: isRHSReplace

    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !

    call calc_HDiffRHS(wz_VorRHS,                          &  !(inout)
         & wz_Vor, hViscTermCoef, hHyperViscTermCoef, 'V', &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_HDiffRHS(wz_DivRHS,                                &  !(inout)
         & wz_Div, hViscTermCoef, hHyperViscTermCoef, 'V', 2d0,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_HDiffRHS(wz_PTempRHS,                               &  !(inout)
         & wz_PTempEdd, hDiffTermCoef, hHyperDiffTermCoef,  'S',  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_HDiffRHS(wz_SaltRHS,                             &  !(inout)
         & wz_Salt, hDiffTermCoef, hHyperDiffTermCoef,  'S',   &  !(in)
         & isRHSReplace=isRHSReplace )

  end subroutine calc_HydroBouEqHViscRHS


  !> @brief
  !!
  !!
  subroutine calc_HydroBouEqVViscRHS( wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS, &
       & wz_Vor, wz_Div, wz_PTemp, wz_Salt, &
       & xyz_VViscTermCoef, vHyperViscTermCoef, xyz_VDiffTermCoef, vHyperDiffTermCoef, &
       & isRHSReplace )

    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTemp, wz_Salt
    real(DP), intent(in) :: xyz_VViscTermCoef(0:iMax-1,jMax,0:kMax), vHyperViscTermCoef
    real(DP), intent(in) :: xyz_VDiffTermCoef(0:iMax-1,jMax,0:kMax), vHyperDiffTermCoef
    logical, intent(in) :: isRHSReplace

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !


    xyz_totDepth(:,:,:) = spread(xy_totDepthBasic, 3, kMax+1)

    call calc_VDiffRHS_wz(wz_VorRHS,                                        &  !(inout)
         & wz_Vor, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_wz(wz_DivRHS,                                        &  !(inout)
         & wz_Div, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_wz(wz_PTempRHS,                               &  !(inout)
         & wz_PTemp, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,   &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_wz(wz_SaltRHS,                               &  !(inout)
         & wz_Salt, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,     &  !(in)
         & isRHSReplace=isRHSReplace )

  end subroutine calc_HydroBouEqVViscRHS
  
  !> @brief 
  !!
  !!
  subroutine calc_HydroBouEqVViscRHS_xyz( xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, &
       & xyz_U, xyz_V, xyz_PTemp, xyz_Salt, &
       & xyz_VViscTermCoef, vHyperViscTermCoef, xyz_VDiffTermCoef, vHyperDiffTermCoef, &
       & isRHSReplace )
    
    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_U, xyz_V, xyz_PTemp, xyz_Salt
    real(DP), intent(in) :: xyz_VViscTermCoef(0:iMax-1,jMax,0:kMax), vHyperViscTermCoef
    real(DP), intent(in) :: xyz_VDiffTermCoef(0:iMax-1,jMax,0:kMax), vHyperDiffTermCoef
    logical, intent(in) :: isRHSReplace

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !


    xyz_totDepth(:,:,:) = spread(xy_totDepthBasic, 3, kMax+1)

    call calc_VDiffRHS_xyz(xyz_URHS,                                        &  !(inout)
         & xyz_U, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_VRHS,                                        &  !(inout)
         & xyz_V, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_PTempRHS,                               &  !(inout)
         & xyz_PTemp, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,   &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_SaltRHS,                               &  !(inout)
         & xyz_Salt, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,     &  !(in)
         & isRHSReplace=isRHSReplace )

  end subroutine calc_HydroBouEqVViscRHS_xyz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! * Subroutines to calculate the right-side hand of each time evolution equation.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief 
  !!
  !!
  subroutine calc_VorEqDivEqInvisRHS(wz_RHSVor, wz_RHSDiv, &
       & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot, &
       & xyz_UrfCori, xyz_VrfCori  &
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(out) :: wz_RHSVor, wz_RHSDiv
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot, &
         & xyz_UrfCori, xyz_VrfCori
    real(DP), intent(in) :: xy_SurfHeight(0:iMax-1,jMax)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_A(0:iMax-1,jMax)
    real(DP) :: xy_B(0:iMax-1,jMax)
    real(DP) :: w_GeoPot(lMax)
    real(DP) :: xy_SinLat(0:iMax-1,jMax)
    real(DP) :: xyz_DSigUrf(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DSigVrf(0:iMax-1,jMax,0:kMax)
    integer :: k

    ! 実行文; Executable statement
    !


    xy_SinLat(:,:) = sin(xyz_Lat(:,:,0))
!!$    xyz_DSigUrf(:,:,:) = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Urf)))
!!$    xyz_DSigVrf(:,:,:) = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Vrf)))
    xyz_DSigUrf(:,:,:) = xyz_DSig_xyz(xyz_Urf)
    xyz_DSigVrf(:,:,:) = xyz_DSig_xyz(xyz_Vrf)
   
    !$omp parallel do private(w_GeoPot, xy_A, xy_B)
    do k=0, kMax

       w_GeoPot(:) = w_xy(xyz_GeoPot(:,:,k))

       !
       xy_A(:,:) = &
            &   xyz_Vor(:,:,k)*xyz_Urf(:,:,k) + 2d0*Omega*xy_SinLat*xyz_UrfCori(:,:,k) &
            & + xyz_SigDot(:,:,k)*xyz_DSigVrf(:,:,k) &
            & + xyz_DensEdd(:,:,k)*xy_GradMu_w(w_GeoPot)/(RefDens*RPlanet)

       xy_B(:,:) = &
            &   xyz_Vor(:,:,k)*xyz_Vrf(:,:,k) + 2d0*Omega*xy_SinLat*xyz_VrfCori(:,:,k) &
            & - xyz_SigDot(:,:,k)*xyz_DSigUrf(:,:,k) &
            & + xyz_DensEdd(:,:,k)*xy_GradLambda_w(w_GeoPot)/(RefDens*RPlanet)

       !
       wz_RHSVor(:,k) = - w_AlphaOptr_xy(xy_A, xy_B)
       
       wz_RHSDiv(:,k) = &
            &    w_AlphaOptr_xy(xy_B, -xy_A) &
            &  - w_Lapla2D_w(w_xy( &
            &        0.5d0*(xyz_Urf(:,:,k)**2 + xyz_Vrf(:,:,k)**2)/(1d0 - xy_SinLat**2) &
            &      + Grav*xy_SurfHeight + xyz_Press(:,:,k)/RefDens                      &
            &    ))
    end do

  end subroutine calc_VorEqDivEqInvisRHS


  !> @brief 
  !!
  !!
  subroutine calc_HDiffRHS(wz_RHSQuant, &
       & wz_Quant, hDiffCoef, hHyperDiffCoef, LhOptrType, LhOptrCoef, &
       & isRHSReplace )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSQuant(lMax, 0:kMax)
    real(DP), intent(in) :: wz_Quant(lMax,0:kMax)
    real(DP), intent(in) :: hDiffCoef, hHyperDiffCoef
    character(1), intent(in) :: LhOptrType              ! 'V': Vector, 'S': Scalar, Default setteing is 'S'
    real(DP), optional, intent(in) :: LhOptrCoef 
    logical, optional, intent(in) :: isRHSReplace
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: ScaLaplaOptrCoef
    integer :: k

    ! 実行文; Executable statement
    !
    
    if( present(isRHSReplace) .and. isRHSReplace ) then
       wz_RHSQuant = 0d0
    end if

    ScaLaplaOptrCoef = 1d0
    if( present(LhOptrCoef) ) ScaLaplaOptrCoef = LhOptrCoef

    if( LhOptrType == 'V' ) then

       !$omp parallel do
       do k=0, kMax
          wz_RHSQuant(:,k) = wz_RHSQuant(:,k)  &
               &      +  hDiffCoef* 2d0*wz_Quant(:,k)/RPlanet**2 &
               &      +  w_Lapla2D_w(  &
               &              ScaLaplaOptrCoef*(hDiffCoef - 2d0*hHyperDiffCoef/RPlanet**2)*wz_Quant(:,k) & 
               &            - ScaLaplaOptrCoef*hHyperDiffCoef*w_Lapla2D_w(wz_Quant(:,k))               &
               &         )
       end do
    else
       !$omp parallel do
       do k=0, kMax
          wz_RHSQuant(:,k) = wz_RHSQuant(:,k) &
            &      + w_Lapla2D_w( hDiffCoef*wz_Quant(:,k) - hHyperDiffCoef*w_Lapla2D_w(wz_Quant(:,k)) )
       end do
    end if

  end subroutine calc_HDiffRHS

  !> @brief 
  !!
  !!
  subroutine calc_VDiffRHS_xyz(xyz_RHSQuant,           & !(inout)
       & xyz_Quant, xyz_vDiffCoef, vHyperDiffCoef, & !(in)
       & xyz_totDepth, isRHSReplace                & !(in)
       & )
    

    use at_module_omp

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_RHSQuant(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Quant(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_vDiffCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: vHyperDiffCoef
    logical, optional, intent(in) :: isRHSReplace
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_VDiffTmp, xyz_QuantDSig4
    integer :: i, j, k

    ! 実行文; Executable statement
    !
    
    if( present(isRHSReplace) .and. isRHSReplace ) then
       xyz_RHSQuant(:,:,:) = 0d0
    end if

    xyz_VDiffTmp(:,:,:) = xyz_DSig_xyz( xyz_vDiffCoef*xyz_DSig_xyz(xyz_Quant) )

    if(vHyperDiffCoef > 0d0) then
       xyz_QuantDSig4(:,:,:) = xyz_DSigDSig_xyz(xyz_DSigDSig_xyz(xyz_Quant))
    else
       xyz_QuantDSig4(:,:,:) = 0d0
    end if
    
    !$omp parallel do
    do k=0, kMax
       xyz_RHSQuant(:,:,k) = &
            &   (   xyz_VDiffTmp(:,:,k)                                         &
            &     - vHyperDiffCoef*xyz_QuantDSig4(:,:,k)/xyz_totDepth(:,:,k)**2 &
            &   )/xyz_totDepth(:,:,k)**2
    end do

  end subroutine calc_VDiffRHS_xyz

  !> @brief 
  !!
  !!
  subroutine calc_VDiffRHS_wz(wz_RHSQuant,           & !(inout)
       & wz_Quant, xyz_vDiffCoef, vHyperDiffCoef, & !(in)
       & xyz_totDepth, isRHSReplace                & !(in)
       & )
    

    use at_module_omp

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSQuant(lMax,0:kMax)
    real(DP), intent(in) :: wz_Quant(lMax,0:kMax)
    real(DP), intent(in) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_vDiffCoef(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: vHyperDiffCoef
    logical, optional, intent(in) :: isRHSReplace
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_RHSQuant(0:iMax-1,jMax,0:kMax)
    
    ! 実行文; Executable statement
    !

    call calc_VDiffRHS_xyz(xyz_RHSQuant, &
         & xyz_wz(wz_Quant), xyz_VDiffCoef, vHyperDiffCoef, &
         & xyz_totDepth, isRHSReplace                       &
         & )
    wz_RHSQuant(:,:) = wz_xyz(xyz_RHSQuant)
    
  end subroutine calc_VDiffRHS_wz
  

  subroutine calc_TracerEqInvisRHS( wz_RHSTracer, xyz_Tracer, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    real(DP), intent(out) :: wz_RHSTracer(lMax, 0:kMax)
    real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    integer :: k
    real(DP) :: xyz_DSigTracer(0:iMax-1,jMax,0:kMax)

    xyz_DSigTracer(:,:,:) = xyz_DSig_xyz(xyz_Tracer)
    
    !$omp parallel do
    do k=0, kMax
       wz_RHSTracer(:,k) =  &
            & - w_AlphaOptr_xy(xyz_Tracer(:,:,k)*xyz_Urf(:,:,k), xyz_Tracer(:,:,k)*xyz_Vrf(:,:,k)) &
            & + w_xy(xyz_Tracer(:,:,k)*xyz_Div(:,:,k) - xyz_SigDot(:,:,k)*xyz_DSigTracer(:,:,k))
    end do

    
!!$    wz_RHSTracer =  &
!!$         & - wz_AlphaOptr_xyz(xyz_Tracer*xyz_Urf, xyz_Tracer*xyz_Vrf) &
!!$         & + wz_xyz(   xyz_Tracer*xyz_Div &
!!$         &           - xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot*xyz_Tracer))) &
!!$         &           + xyz_Tracer*xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot))) )

  end subroutine calc_TracerEqInvisRHS

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
  subroutine correct_DivEqRHSUnderRigidLid( wz_DivRHS, xy_SurfPressA, &
       & wz_DivOld, xy_SurfPressRef, xy_SurfPressN, xy_SurfPressB, dDivdt_coef, &
       & TIntType_SurfPressTerm )
    


    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_DivRHS(lMax, 0:kMax)
    real(DP), intent(inout) :: xy_SurfPressA(0:iMax-1,jMax)
    real(DP), intent(in) :: wz_DivOld(lMax,0:kMax)
    real(DP), intent(in) :: xy_SurfPressRef(0:iMax-1,jMax), xy_SurfPressN(0:iMax-1,jMax), xy_SurfPressB(0:iMax-1,jMax)
    real(DP), intent(in) :: dDivdt_coef
    character(*), intent(in) :: TIntType_SurfPressTerm

    ! 局所変数
    ! Local variables
    !

    !> Horizontal divergence of velcoity potential used in Pressure Correction method
    !! to keep velocity divergence-free. 
    real(DP) :: w_HDivVelPot(lMax)              

    !> Coefficients used to exptrapolate new value of surface pressure.  
    real(DP) :: ecoefPA, ecoefPB, ecoefPN

    integer :: k

    ! 実行文; Executable statement
    !


    ! 
    ! 

    select case(TIntType_SurfPressTerm)
       case("BEuler1")
          ecoefPA = 1d0; ecoefPN = 0d0; ecoefPB = 0d0;
       case("CRANKNIC")
          ecoefPA = 0.5d0; ecoefPN = 0.5d0; ecoefPB = 0d0;
       case("CRANKNIC_WithLF")
          ecoefPA = 0.5d0; ecoefPN = 0d0; ecoefPB = 0.5d0;
       case("AM3")
          ecoefPA = 5d0/12d0; ecoefPN = 8d0/12d0; ecoefPB = - 1d0/12d0;
       case default
          call MessageNotify('E', module_name, &
               & "'%c' is not allowable as temporal integration method for surface pressure term", &
               & c1 = TIntType_SurfPressTerm )
    end select

    w_HDivVelPot = RefDens * w_IntSig_BtmToTop_wz( wz_DivOld*dDivdt_coef + wz_DivRHS )


    xy_SurfPressA = ( &
         & xy_w( w_InvLapla2D_w( w_HDivVelPot ) ) + xy_SurfPressRef &
         & - ecoefPN*xy_SurfPressN - ecoefPB*xy_SurfPressB &
         & )/ecoefPA


    wz_DivRHS(:,:) =   wz_DivRHS(:,:) &
         &      - spread(w_HDivVelPot, 2, kMax+1)/RefDens

  end subroutine correct_DivEqRHSUnderRigidLid

  !> @brief 
  !! 
  !!
  subroutine correct_DivEqRHSUnderRigidLid2( wz_DivA, xy_SurfPressA, &
       & xy_SurfPressRef, xy_SurfPressN, xy_SurfPressB, dDivdt_coef, &
       & TIntType_SurfPressTerm )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_DivA(lMax, 0:tMax)
    real(DP), intent(out) :: xy_SurfPressA(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SurfPressRef(0:iMax-1,jMax), xy_SurfPressN(0:iMax-1,jMax), xy_SurfPressB(0:iMax-1,jMax)
    real(DP), intent(in) :: dDivdt_coef
    character(*), intent(in) :: TIntType_SurfPressTerm

    ! 局所変数
    ! Local variables
    !

    !> Horizontal divergence of velcoity potential used in Pressure Correction method
    !! to keep velocity divergence-free. 
    real(DP) :: w_HDivPressInvRho(lMax)              

    !> Coefficients used to exptrapolate new value of surface pressure.  
    real(DP) :: ecoefPA, ecoefPB, ecoefPN

    integer :: k

    ! 実行文; Executable statement
    !


    ! 
    ! 

    select case(TIntType_SurfPressTerm)
       case("BEuler1")
          ecoefPA = 1d0; ecoefPN = 0d0; ecoefPB = 0d0;
       case("CRANKNIC")
          ecoefPA = 0.5d0; ecoefPN = 0.5d0; ecoefPB = 0d0;
       case("CRANKNIC_WithLF")
          ecoefPA = 0.5d0; ecoefPN = 0d0; ecoefPB = 0.5d0;
       case("AM3")
          ecoefPA = 5d0/12d0; ecoefPN = 8d0/12d0; ecoefPB = - 1d0/12d0;
       case default
          call MessageNotify('E', module_name, &
               & "'%c' is not allowable as temporal integration method for surface pressure term", &
               & c1 = TIntType_SurfPressTerm )
    end select


    w_HDivPressInvRho = w_IntSig_BtmToTop_wz( wz_DivA*dDivdt_coef )


    xy_SurfPressA = ( &
         & xy_w( w_InvLapla2D_w( w_HDivPressInvRho*RefDens ) ) + xy_SurfPressRef &
         & - ecoefPN*xy_SurfPressN - ecoefPB*xy_SurfPressB &
         & )/ecoefPA



    wz_DivA(:,:) = &
         &        wz_DivA(:,:) &
         &      - spread(w_HDivPressInvRho, 2, kMax+1)/dDivdt_coef

  end subroutine correct_DivEqRHSUnderRigidLid2

  !> @brief 
  !! 
  !!
  subroutine correct_DivVorEqRHSUnderRigidLid( w_DDivCorrect, w_DVorCorrect, xy_SurfPressA,    &
       & wz_DivOld, xy_SurfPressRef, xy_SurfPressN, xy_SurfPressB, dt, CoriTermCoef,  &
       & TIntType_SurfPressTerm )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax) :: w_DDivCorrect, w_DVorCorrect
    real(DP), intent(out) :: xy_SurfPressA(0:iMax-1,jMax)
    real(DP), intent(in) :: wz_DivOld(lMax,0:kMax)
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfPressRef, xy_SurfPressN, xy_SurfPressB
    real(DP), intent(in) :: dt, CoriTermCoef
    character(*), intent(in) :: TIntType_SurfPressTerm

    
    ! 局所変数
    ! Local variables
    !

    !> Horizontal divergence of velcoity potential used in Pressure Correction method
    !! to keep velocity divergence-free. 
    real(DP) :: w_HDivPressInvRho(lMax)              

    real(DP) :: w_PhiInvRho(lMax)
    real(DP), dimension(0:iMax-1,jMax) :: xy_UrfCorrect, xy_VrfCorrect
    
    !> Coefficients used to exptrapolate new value of surface pressure.  
    real(DP) :: ecoefPA, ecoefPB, ecoefPN

    integer :: k, itr
    real(DP) :: theta
    
    ! 実行文; Executable statement
    !


    ! 
    !
    theta = 1d0    
    select case(TIntType_SurfPressTerm)
       case("BEuler1")
          ecoefPA = 1d0; ecoefPN = 0d0; ecoefPB = 0d0;
       case("CRANKNIC")
          ecoefPA = theta; ecoefPN = 1d0 - theta; ecoefPB = 0d0;
       case("CRANKNIC_WithLF")
          ecoefPA = theta; ecoefPN = 0d0; ecoefPB = 1d0 - theta;
       case("AM3")
          ecoefPA = 5d0/12d0; ecoefPN = 8d0/12d0; ecoefPB = - 1d0/12d0;
       case default
          call MessageNotify('E', module_name, &
               & "'%c' is not allowable as temporal integration method for surface pressure term", &
               & c1 = TIntType_SurfPressTerm )
    end select


    w_HDivPressInvRho = w_IntSig_BtmToTop_wz( wz_DivOld/dt )
    w_PhiInvRho = w_InvLapla2D_w(w_HDivPressInvRho)

    w_DDivCorrect(:) = - dt*w_HDivPressInvRho
    xy_VrfCorrect(:,:) = - dt*xy_GradMu_w(w_PhiInvRho)/RPlanet
    do itr=1, 2
       w_DVorCorrect(:) = - w_xy( &
              & CoriTermCoef*dt*2d0*Omega*( &
              &     sin(xyz_Lat(:,:,0))*xy_w(w_DDivCorrect) &
              &   + xy_VrfCorrect/RPlanet                  &
              & ) &
            & )
       call w_VorDiv2VectorCosLat(w_DVorCorrect, w_DDivCorrect, &
            & xy_UrfCorrect, xy_VrfCorrect)
    end do


    xy_SurfPressA = ( &
         &   RefDens*xy_w( &
         &       w_PhiInvRho &
         &     - w_InvLapla2D_w(CoriTermCoef*2d0*Omega*w_xy(    &
         &       - sin(xyz_Lat(:,:,0))*xy_w(w_DVorCorrect)         &
         &       + xy_UrfCorrect/RPlanet                           &
         &     )) &
         &   ) &
         & + xy_SurfPressRef &
         & - ecoefPN*xy_SurfPressN - ecoefPB*xy_SurfPressB &
         & )/ecoefPA


  end subroutine correct_DivVorEqRHSUnderRigidLid

end module HydroBouEqSolverRHS_old_mod

