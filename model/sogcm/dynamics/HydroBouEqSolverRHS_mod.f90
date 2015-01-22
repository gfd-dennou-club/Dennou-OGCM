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
module HydroBouEqSolverRHS_mod 

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
  public :: calc_HydroBouEqHViscRHS, calc_HydroBouEqVViscRHS

  ! Calculation routines
  public :: calc_VorEqDivEqInvisRHS, calc_TracerEqInvisRHS
  public :: calc_HDiffRHS, calc_VDiffRHS
  public :: calc_SurfHeightRHS

  public :: correct_DivEqRHSUnderRigidLid, correct_DivEqRHSUnderRigidLid2

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
       & xy_SurfPress, xy_totDepthBasic, z_PTempBasic )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(out), dimension(lMax) :: w_SurfHeightRHS
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Urf, xyz_Vrf, xyz_Vor, xyz_Div, xyz_PTempEdd, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_SurfHeight, xy_SurfPress, xy_totDepthBasic
    real(DP), intent(in), dimension(0:kMax) :: z_PTempBasic

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
         & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot)

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
       & vViscTermCoef, vHyperViscTermCoef, vDiffTermCoef, vHyperDiffTermCoef, &
       & isRHSReplace )
    
    !
    !
    use VariableSet_mod, only: &
         & xy_totDepthBasic, z_PTempBasic

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(lMax,0:kMax) :: wz_VorRHS, wz_DivRHS, wz_PTempRHS, wz_SaltRHS
    real(DP), intent(in), dimension(lMax,0:kMax) :: wz_Vor, wz_Div, wz_PTemp, wz_Salt
    real(DP), intent(in) :: vViscTermCoef, vHyperViscTermCoef, vDiffTermCoef, vHyperDiffTermCoef
    logical, intent(in) :: isRHSReplace

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !


    xyz_totDepth = spread(xy_totDepthBasic, 3, kMax+1)

    call calc_VDiffRHS(wz_VorRHS,                                    &  !(inout)
         & wz_Vor, vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS(wz_DivRHS,                                    &  !(inout)
         & wz_Div, vViscTermCoef, vHyperViscTermCoef, xyz_totDepth,  &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS(wz_PTempRHS,                               &  !(inout)
         & wz_PTemp, vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,            &  !(in)
         & isRHSReplace=isRHSReplace )

    call calc_VDiffRHS(wz_SaltRHS,                               &  !(inout)
         & wz_Salt, vDiffTermCoef, vHyperDiffTermCoef, xyz_totDepth,            &  !(in)
         & isRHSReplace=isRHSReplace )

  end subroutine calc_HydroBouEqVViscRHS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! * Subroutines to calculate the right-side hand of each time evolution equation.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief 
  !!
  !!
  subroutine calc_VorEqDivEqInvisRHS(wz_RHSVor, wz_RHSDiv, &
       & xyz_Vor, xyz_Urf, xyz_Vrf, xy_SurfHeight, xyz_DensEdd, xyz_Press, xyz_GeoPot, xyz_SigDot &
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
    real(DP), intent(in) :: xyz_Press(0:iMax-1,jMax,0:kMax)
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

!real(DP) :: w(lMax)

    ! 実行文; Executable statement
    !

    !$omp parallel workshare
    xyz_KinEngy = (xyz_Urf**2 + xyz_Vrf**2)/(2d0*cos(xyz_Lat)**2)
    xyz_AbsVor = xyz_Vor + 2d0*Omega*sin(xyz_Lat)
    xyz_GeoPotGradCoef = xyz_DensEdd/(RefDens)
    !$omp end parallel workshare

    wz_GeoPot(:,:) = wz_xyz(xyz_GeoPot)

    xyz_A(:,:,:) = &
         &   xyz_AbsVor*xyz_Urf + xyz_SigDot*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Vrf))) &
         & + xyz_GeoPotGradCoef*xyz_GradMu_wz(wz_GeoPot)
    xyz_B(:,:,:) = &
         &   xyz_AbsVor*xyz_Vrf - xyz_SigDot*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Urf))) &
         & - xyz_GeoPotGradCoef*xyz_GradLambda_wz(wz_GeoPot)


    wz_RHSVor(:,:) = - wz_AlphaOptr_xyz( xyz_A, xyz_B )

    wz_RHSDiv(:,:) = wz_AlphaOptr_xyz( xyz_B, -xyz_A ) &
         &      - wz_Lapla2D_wz(wz_xyz( &
         &             xyz_KinEngy + Grav*spread(xy_SurfHeight,3,kMax+1) &
         &           + xyz_Press/RefDens &
         &        ))

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

    ! 実行文; Executable statement
    !
    
    if( present(isRHSReplace) .and. isRHSReplace ) then
       wz_RHSQuant = 0d0
    end if

    ScaLaplaOptrCoef = 1d0
    if( present(LhOptrCoef) ) ScaLaplaOptrCoef = LhOptrCoef

    if( LhOptrType == 'V' ) then
       wz_RHSQuant = wz_RHSQuant  &
            &      +  hDiffCoef* 2d0*wz_Quant/RPlanet**2 &
            &      +  wz_Lapla2D_wz(  &
            &              ScaLaplaOptrCoef*(hDiffCoef - 2d0*hHyperDiffCoef/RPlanet**2)*wz_Quant & 
            &            - ScaLaplaOptrCoef*hHyperDiffCoef*wz_Lapla2D_wz(wz_Quant)               &
            &         )
    else
       wz_RHSQuant = wz_RHSQuant &
            &      + wz_Lapla2D_wz( hDiffCoef*wz_Quant - hHyperDiffCoef*wz_Lapla2D_wz(wz_Quant) )
    end if

  end subroutine calc_HDiffRHS

  !> @brief 
  !!
  !!
  subroutine calc_VDiffRHS(wz_RHSQuant,       & !(inout)
       & wz_Quant, vDiffCoef, vHyperDiffCoef, & !(in)
       & xyz_totDepth, isRHSReplace           & !(in)
       & )
    

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: wz_RHSQuant(lMax, 0:kMax)
    real(DP), intent(in) :: wz_Quant(lMax,0:kMax)
    real(DP), intent(in) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: vDiffCoef, vHyperDiffCoef
    logical, optional, intent(in) :: isRHSReplace
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyt_QuantDIVDep2(0:iMax-1, jMax, 0:tMax)
    real(DP) :: xyt_QuantDIVDep4(0:iMax-1, jMax, 0:tMax)

    ! 実行文; Executable statement
    !
    
    if( present(isRHSReplace) .and. isRHSReplace ) then
       wz_RHSQuant = 0d0
    end if

    xyt_QuantDIVDep2(:,:,:) = xyt_xyz(xyz_wz(wz_Quant)/xyz_totDepth**2)   
    xyt_QuantDIVDep4(:,:,:) = xyt_xyz(xyz_wz(wz_Quant)/xyz_totDepth**4)   

    wz_RHSQuant(:,:) = wz_RHSQuant  + wz_xyz(xyz_xyt( &
         &          xyt_DSigDSig_xyt( &
         &              vDiffCoef*xyt_QuantDIVDep2                                     &
         &            - vHyperDiffCoef*xyt_DSigDSig_xyt(xyt_QuantDIVDep4)              &
         &          )                                                                  &
         &      ))

  end subroutine calc_VDiffRHS


  subroutine calc_TracerEqInvisRHS( wz_RHSTracer, xyz_Tracer, xyz_Urf, xyz_Vrf, xyz_Div, xyz_SigDot )
    real(DP), intent(out) :: wz_RHSTracer(lMax, 0:kMax)
    real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    wz_RHSTracer(:,:) =  &
         & - wz_AlphaOptr_xyz(xyz_Tracer*xyz_Urf, xyz_Tracer*xyz_Vrf) &
         & + wz_xyz(   xyz_Tracer*xyz_Div & 
         &           - xyz_SigDot*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer))) )

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
  subroutine correct_DivEqRHSUnderRigidLid2( wt_DivA, xy_SurfPressA, &
       & wz_DivOld, xy_SurfPressRef, xy_SurfPressN, xy_SurfPressB, dDivdt_coef, &
       & TIntType_SurfPressTerm )
    


    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: wt_DivA(lMax, 0:tMax)
    real(DP), intent(out) :: xy_SurfPressA(0:iMax-1,jMax)
    real(DP), intent(in) :: wz_DivOld(lMax,0:kMax)
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


    w_HDivPressInvRho = w_IntSig_BtmToTop_wz( wz_DivOld*dDivdt_coef )


    xy_SurfPressA = ( &
         & xy_w( w_InvLapla2D_w( w_HDivPressInvRho*RefDens ) ) + xy_SurfPressRef &
         & - ecoefPN*xy_SurfPressN - ecoefPB*xy_SurfPressB &
         & )/ecoefPA



    wt_DivA(:,:) = wt_wz( &
         &        wz_DivOld(:,:) &
         &      - spread(w_HDivPressInvRho, 2, kMax+1)/dDivdt_coef ) 

  end subroutine correct_DivEqRHSUnderRigidLid2

end module HydroBouEqSolverRHS_mod

