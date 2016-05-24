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
module HydroBouEqSolverRHS_v2_mod 

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
  public :: calc_HydroBouEqRHSExpl_xyz
!  public :: calc_HydroBouEqHViscRHS
  public :: calc_HydroBouEqVViscRHS_xyz

  ! Calculation routines
  public :: calc_VorEqDivEqInvisRHS
  public :: calc_TracerEqInvisRHS_wz, calc_TracerEqInvisRHS_xyz
  public :: calc_HDiffRHS, calc_VDiffRHS_xyz
  public :: calc_SurfHeightRHS_xy

  public :: Advance_BarotEqStep, Advance_BarotEqStep2

  
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
  subroutine calc_HydroBouEqRHSExpl_xyz( xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS, xy_SurfHeightRHS, &
       & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt, xy_SurfHeight, &
       & xy_SurfPress, xy_totDepthBasic, z_PTempBasic, &
       & hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef, hHyperDiffTermCoef, &       
       & xyz_UrfCori, xyz_VrfCori )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax,0:kMax) :: xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: xy_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight, xy_SurfPress, xy_totDepthBasic
    real(DP), intent(in), dimension(0:kMax) :: z_PTempBasic
    real(DP), intent(in) :: hViscTermCoef, hHyperViscTermCoef, hDiffTermCoef, hHyperDiffTermCoef    
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_UrfCori, xyz_VrfCori

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, &
         & xyz_SigDot, xyz_GeoPot, xyz_DensEdd, xyz_Press, xyz_PTemp

    real(DP), dimension(lMax,0:kMax) :: &
         & wz_Vor, wz_Div, &
         & wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    
    real(DP), dimension(0:iMax-1,jMax) :: xy_totDepth, xy_CosLat
    integer :: i, j, k

    ! 実行文; Executable statement
    !

    xy_CosLat(:,:) = cos(xyz_Lat(:,:,0))
    !$omp parallel do
    do k=0, kMax
       xyz_PTemp(:,:,k) = xyz_PTempEdd(:,:,k) + z_PTempBasic(k)
       xyz_Urf(:,:,k) = xyz_U(:,:,k)*xy_CosLat(:,:)
       xyz_Vrf(:,:,k) = xyz_V(:,:,k)*xy_CosLat(:,:)       
    end do

    call wz_VectorCosLat2VorDiv( xyz_Urf, xyz_Vrf, & ! (in)
         & wz_Vor, wz_Div                          & ! (out)
         & )
    xyz_Vor(:,:,:) = xyz_wz(wz_Vor); xyz_Div(:,:,:) = xyz_wz(wz_Div)

    !
    !
    
    xy_totDepth(:,:)   = xy_totDepthBasic! + xy_SurfHeightN
    xyz_SigDot(:,:,:)  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )
    xyz_GeoPot(:,:,:)  = Diagnose_GeoPot( xy_totDepth )
    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*xyz_GeoPot )     ! (in)


    ! Calculate the pressure which is deviation from -RefDens*Grav*z).
    !
    xyz_Press(:,:,:) =   spread(xy_SurfPress, 3, kMax+1)                    &    ! barotropic component
                &      + Diagnose_HydroPressEdd( xy_totDepth, xyz_DensEdd )      ! baroclinic component

    !
    !
    call calc_VorEqDivEqInvisRHS(wz_VorRHS, wz_DivRHS, &
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot, &
         & xyz_UrfCori, xyz_VrfCori )

    Call calc_TracerEqInvisRHS_wz(wz_PTempRHS, &
         & xyz_PTemp, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )

    call calc_TracerEqInvisRHS_wz(wz_SaltRHS, &
         & xyz_Salt, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    
    call calc_SurfHeightRHS_xy(xy_SurfHeightRHS, &
         & xyz_Urf, xyz_Vrf, xy_totDepth ) 

    ! Horizontal mixing (with constant viscosity or diffusivity)
    !
    call calc_HDiffRHS(wz_VorRHS,                          &  !(inout)
         & wz_Vor, hViscTermCoef, hHyperViscTermCoef, 'V', &  !(in)
         & isRHSReplace=.false. )

    call calc_HDiffRHS(wz_DivRHS,                                &  !(inout)
         & wz_Div, hViscTermCoef, hHyperViscTermCoef, 'V', 2d0,  &  !(in)
         & isRHSReplace=.false. )

    call calc_HDiffRHS(wz_PTempRHS,                                        &  !(inout)
         & wz_xyz(xyz_PTempEdd), hDiffTermCoef, hHyperDiffTermCoef,  'S',  &  !(in)
         & isRHSReplace=.false. )

    call calc_HDiffRHS(wz_SaltRHS,                                      &  !(inout)
         & wz_xyz(xyz_Salt), hDiffTermCoef, hHyperDiffTermCoef,  'S',   &  !(in)
         & isRHSReplace=.false. )

    
    ! Spectral space -> Grid space
    !
    
    call wz_VorDiv2VectorCosLat( wz_VorRHS, wz_DivRHS, & ! (in)
         & xyz_URHS, xyz_VRHS                        )   ! (out)
    do k=0, kMax
       xyz_URHS(:,:,k) = xyz_URHS(:,:,k)/xy_CosLat(:,:)
       xyz_VRHS(:,:,k) = xyz_VRHS(:,:,k)/xy_CosLat(:,:)
    end do
    xyz_PTempRHS(:,:,:) = xyz_wz(wz_PTempRHS)    
    xyz_SaltRHS(:,:,:) = xyz_wz(wz_SaltRHS)
    
  end subroutine calc_HydroBouEqRHSExpl_xyz
  


  !> @brief 
  !!
  !!
  subroutine calc_HydroBouEqVViscRHS_xyz( &
       & xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS,                                 &  ! (inout) 
       & xyz_U, xyz_V, xyz_PTemp, xyz_Salt,                                             &  ! (in)
       & xyz_VViscTermCoef, vHyperViscTermCoef, xyz_VDiffTermCoef, vHyperDiffTermCoef,  &  ! (in)
       & isRHSReplace )                                                                    ! (in)
    
    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: xyz_URHS, xyz_VRHS, xyz_PTempRHS, xyz_SaltRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_U, xyz_V, xyz_PTemp, xyz_Salt
    real(DP), intent(in) :: xyz_VViscTermCoef(0:iMax-1,jMax,0:kMax), vHyperViscTermCoef
    real(DP), intent(in) :: xyz_VDiffTermCoef(0:iMax-1,jMax,0:kMax), vHyperDiffTermCoef
    logical, intent(in) :: isRHSReplace

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: xyz_totDepth
    
    ! 実行文; Executable statement
    !


    xyz_totDepth(:,:,:) = spread(xy_totDepthBasic, 3, kMax+1)
    
    call calc_VDiffRHS_xyz(xyz_URHS,                                        &  !(inout)
         & xyz_U, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,      &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_VRHS,                                        &  !(inout)
         & xyz_V, xyz_vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,      &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_PTempRHS,                                    &  !(inout)
         & xyz_PTemp, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS_xyz(xyz_SaltRHS,                                     &  !(inout)
         & xyz_Salt, xyz_vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,   &  !(in)
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
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_VDiffTmp, xyz_QuantDSig4
    integer :: k

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


  subroutine calc_TracerEqInvisRHS_xyz( &
       & xyz_RHSTracer,                                     & ! (out)
       & xyz_Tracer, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot  & ! (in)
       & )

    real(DP), intent(out) :: xyz_RHSTracer(0:iMax-1,jMax, 0:kMax)
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
       xyz_RHSTracer(:,:,k) =  &
            & - xy_w(w_AlphaOptr_xy(xyz_Tracer(:,:,k)*xyz_Urf(:,:,k), xyz_Tracer(:,:,k)*xyz_Vrf(:,:,k))) &
            & + xyz_Tracer(:,:,k)*xyz_Div(:,:,k) - xyz_SigDot(:,:,k)*xyz_DSigTracer(:,:,k)
    end do

    
!!$    wz_RHSTracer =  &
!!$         & - wz_AlphaOptr_xyz(xyz_Tracer*xyz_Urf, xyz_Tracer*xyz_Vrf) &
!!$         & + wz_xyz(   xyz_Tracer*xyz_Div &
!!$         &           - xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot*xyz_Tracer))) &
!!$         &           + xyz_Tracer*xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot))) )

  end subroutine calc_TracerEqInvisRHS_xyz

  subroutine calc_TracerEqInvisRHS_wz( &
       & wz_RHSTracer,                                     & ! (out)
       & xyz_Tracer, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot & ! (in)
       & )
    
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

  end subroutine calc_TracerEqInvisRHS_wz

  
  subroutine calc_SurfHeightRHS_xy( xy_RHSSurfHeight, xyz_Urf, xyz_Vrf, xy_totDepth )
    real(DP), intent(out) :: xy_RHSSurfHeight(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)

    xy_RHSSurfHeight =  0d0!- w_AlphaOptr_xy( &
!         & xy_totDepth*xy_IntSig_BtmToTop_xyz(xyz_Urf), xy_totDepth*xy_IntSig_BtmToTop_xyz(xyz_Vrf) )

  end subroutine calc_SurfHeightRHS_xy

  subroutine Advance_BarotEqStep2( &
       & xyz_U, xyz_V, xy_SurfPressA, xy_SurfPressRef, dt)
    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: xy_SurfPressA
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: xy_SurfPressRef
    real(DP), intent(in) :: dt

    real(DP), dimension(lMax) :: w_HDivPressInvRho, w_PhiInvRho, w_DDiv, w_DVor
    real(DP), dimension(lMax,0:kMax) :: wz_Div, wz_Vor
    real(DP), dimension(0:iMax-1,jMax) :: xy_DUrf, xy_DVrf, xy_CosLat

    xy_CosLat(:,:) = cos(xyz_Lat(:,:,0))
    
    call wz_VectorCosLat2VorDiv(xyz_U*spread(xy_CosLat,3,kMax+1), xyz_V*spread(xy_CosLat,3,kMax+1), &
         & wz_Vor, wz_Div)

    w_HDivPressInvRho = w_IntSig_BtmToTop_wz(wz_Div/dt)
    w_PhiInvRho = w_InvLapla2D_w(w_HDivPressInvRho)

    w_DDiv(:) = -dt*w_HDivPressInvRho; w_DVor(:) = 0d0
    call w_VorDiv2VectorCosLat(w_DVor, w_DDiv, &
         & xy_DUrf, xy_DVrf)

    xy_SurfPressA(:,:) = xy_SurfPressRef + RefDens*xy_w(w_PhiInvRho)
    xyz_U(:,:,:) = xyz_U + spread(xy_DUrf/xy_CosLat, 3, kMax+1)
    xyz_V(:,:,:) = xyz_V + spread(xy_DVrf/xy_CosLat, 3, kMax+1)
    
  end subroutine Advance_BarotEqStep2
  !> @brief 
  !! 
  !!
  subroutine Advance_BarotEqStep( &
       & xy_UBarotA, xy_VBarotA, xy_SurfPressA,                     &  ! (out)
       & xy_UBarot, xy_VBarot, xy_SurfPress,                        &  ! (in)
       & xy_UBarotOld, xy_VBarotOld, xy_SurfPressOld,               &  ! (in)
       & xy_ForceUBaroc, xy_ForceVBaroc,                            &  ! (in)
       & dt, CoriTermCoefA                                          &  ! (in)
       & )
        
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: &
         & xy_UBarotA, xy_VBarotA, xy_SurfPressA
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_UBarot, xy_VBarot, xy_SurfPress, &
         & xy_UBarotOld, xy_VBarotOld, xy_SurfPressOld, &
         & xy_ForceUBaroc, xy_ForceVBaroc
         
    real(DP), intent(in) :: dt, CoriTermCoefA

    
    ! 局所変数
    ! Local variables
    !

    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_UBarot_RHSEx, xy_VBarot_RHSEx

    real(DP), dimension(lMax) :: &
         & w_DivA, w_VorA, w_SurfPress, &
         & w_DDiv, w_DVor, w_DSurfPress
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_DUrf, xy_DVrf
    
    
    real(DP) :: dt_
    real(DP), dimension(0:iMax-1,jMax) :: xy_CosLat, xy_SinLat
    real(DP), dimension(0:iMax-1,jMax) :: xy_CoriParam, xy_CoriImpFac, xy_DSurfPress
    real(DP) :: PressTermCoefA, UVPAvgCoefA, UVPAvgCoefN, UVPAvgCoefB
    real(DP) :: theta

    integer :: itr
    
    ! 実行文; Executable statement
    !


    ! 
    !
    dt_ = dt
    theta = 1d0 - sqrt(2d0)/2d0

    PressTermCoefA = 0.5d0
!    if(CoriTermCoefA > 0d0) PressTermCoefA = CoriTermCoefA
    UVPAvgCoefA = 0.53d0!1d0/3d0
    UVPAvgCoefN = 0d0!1d0/3d0
    UVPAvgCoefB = 1d0 - UVPAvgCoefA - UVPAvgCoefN
    
    !
    !
    

    xy_SinLat(:,:) = sin(xyz_Lat(:,:,0))
    xy_CosLat(:,:) = cos(xyz_Lat(:,:,0))
    xy_CoriParam(:,:) = 2d0*Omega*xy_SinLat

    w_SurfPress(:) = w_xy( UVPAvgCoefN*xy_SurfPress + (1d0 - UVPAvgCoefN)*xy_SurfPressOld )
    
    !
    !    
    xy_UBarot_RHSEx(:,:) = &
!         &   xy_CoriParam*xy_VBarot                             &
         &   xy_CoriParam*( UVPAvgCoefN*xy_VBarot  + (1d0 - UVPAvgCoefN)*xy_VBarotOld ) &
         & - xy_GradLon_w(w_SurfPress)/(RefDens*RPlanet)                                &
         & + xy_ForceUBaroc
    xy_VBarot_RHSEx(:,:) = &
!         & - xy_CoriParam*xy_UBarot &
         & - xy_CoriParam*( UVPAvgCoefN*xy_UBarot  + (1d0 - UVPAvgCoefN)*xy_UBarotOld ) &
         & - xy_GradLat_w(w_SurfPress)/(RefDens*RPlanet)                                &
         & + xy_ForceVBaroc

    xy_CoriImpFac(:,:) = UVPAvgCoefA*dt_*xy_CoriParam(:,:)
    xy_UBarotA(:,:) = xy_UBarotOld + dt_*( &
         & (xy_UBarot_RHSEx + xy_CoriImpFac*xy_VBarot_RHSEx)/(1d0 + xy_CoriImpFac**2) &
         & )
    xy_VBarotA(:,:) = xy_VBarotOld + dt_*( &
         & (xy_VBarot_RHSEx - xy_CoriImpFac*xy_UBarot_RHSEx)/(1d0 + xy_CoriImpFac**2) &
         & )

    !* = Projection step ==
    !  In addition to Prresure gradient term, coriolis term is also considered in this
    ! projection step.
    !
    call w_VectorCosLat2VorDiv_2( xy_UBarotA*xy_CosLat, xy_VBarotA*xy_CosLat, & ! (in)
         & w_VorA, w_DivA                                                     & ! (out)
         & )

    w_DDiv(:) = - w_DivA(:); w_DVor(:) = 0d0
    xy_DUrf(:,:) = 0d0
    xy_DVrf(:,:) = - xy_GradMu_w(w_InvLapla2D_w(w_DDiv))/RPlanet
    do itr=1, 2
       w_DVor(:) = - w_xy( &
            &   UVPAvgCoefA*dt_*2d0*Omega*( &
            &     xy_SinLat*xy_w(w_DDiv) + xy_DVrf/RPlanet &
            &   ) &
            & )
       call w_VorDiv2VectorCosLat(w_DVor, w_DDiv, &
            & xy_DUrf, xy_DVrf)
    end do

    w_DSurfPress(:) = RefDens*w_InvLapla2D_w( &
         &      w_DivA/(UVPAvgCoefA*dt_)                                       &
         &    + w_xy( 2d0*Omega*(xy_SinLat*xy_w(w_DVor) - xy_DUrf/RPlanet) )   &
         & )
    
    xy_SurfPressA(:,:) = xy_SurfPressOld + xy_w(w_DSurfPress)    
    xy_UBarotA(:,:) = xy_UBarotA + xy_DUrf/xy_CosLat
    xy_VBarotA(:,:) = xy_VBarotA + xy_DVrf/xy_CosLat

!***** If coriolis term is explicitly treated, above code for presure correction become very
!**** simple as follows.
!!$
!!$    UVPAvgCoefA = 1d0
!!$    w_DSurfPress(:) = RefDens*w_InvLapla2D_w( w_DivA/(UVPAvgCoefA*dt_) )
!!$    xy_SurfPressA(:,:) = xy_SurfPress + xy_w(w_DSurfPress)        
!!$    xy_UBarotA(:,:) = xy_UBarotA  &
!!$         & - UVPAvgCoefA*dt_ * xy_GradLon_w(w_DSurfPress)/(RefDens*RPlanet)
!!$    
!!$    xy_VBarotA(:,:) = xy_VBarotA  &
!!$         & - UVPAvgCoefA*dt_ * xy_GradLat_w(w_DSurfPress)/(RefDens*RPlanet)
!*****************************************************************************************
    
  end subroutine Advance_BarotEqStep

end module HydroBouEqSolverRHS_v2_mod

