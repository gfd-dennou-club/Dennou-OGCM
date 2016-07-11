!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module NonDynMixedLyr_TimeInteg_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING 

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use Constants_mod, only: &
       & PI, Omega, Grav, RPlanet, Cp0, &
       & hViscCoef, vViscCoef, &
       & hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef, &
       & hHyperDiffCoef, vHyperDiffCoef, &
       & RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon, z_LyrThickSig

  use SpmlUtil_mod

  use GovernEqSet_mod, only: &
       & GOVERNEQSET_PHYSICS_EDDYMIX_NAME, &
       & GOVERNEQSET_PHYSICS_CONVADJUST_NAME, &       
       & isPhysicsCompActivated

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType, &
       & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress, &
       & ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed, & 
       & SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed, & 
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom
  
!!$  use TemporalIntegUtil_mod, only: &
!!$       & TemporalIntegUtil_Init, TemporalIntegUtil_Final, &
!!$       & timeIntMode_Euler, xy_timeIntEuler, wt_timeIntEuler, &
!!$       & timeIntMode_LFTR, xy_timeIntLFTR, wt_timeIntLFTR, &
!!$       & timeIntMode_LFAM3, xy_timeIntLFAM3, wt_timeIntLFAM3, &
!!$       & timeIntMode_RK2, xy_timeIntRK2, wt_timeIntRK2, &
!!$       & timeIntMode_RK4, xy_timeIntRK4, wt_timeIntRK4, &
!!$       & TemporalIntegUtil_GetDDtCoef, &
!!$       & TemporalIntegUtil_SetDelTime

  use TemporalIntegUtil_mod2
  
  use TemporalIntegSet_mod, only: &
       & CurrentTimeStep, SubCycleNum, &
       & nShortTimeLevel, &
       & CoriolisTermACoef, VDiffTermACoef


  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: NonDynMixedLyr_TimeInteg_Init, NonDynMixedLyr_TimeInteg_Final
  public :: NonDynMixedLyr_AdvanceTStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'NonDynMixedLyr_mod' !< Module Name
  real(DP), dimension(:,:), allocatable :: xy_CosLat

contains

  !>
  !!
  !!
  subroutine NonDynMixedLyr_TimeInteg_Init()

    ! モジュール引用; Use statement
    !

    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    call TemporalIntegUtil_Init(0d0)! iMax, jMax, kMax, lMax, tMax, 0d0 )

  end subroutine NonDynMixedLyr_TimeInteg_Init

  !>
  !!
  !!
  subroutine NonDynMixedLyr_TimeInteg_Final()

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !


    call TemporalIntegUtil_Final()

  end subroutine NonDynMixedLyr_TimeInteg_Final

  !>
  !!
  !!
  subroutine NonDynMixedLyr_AdvanceTStep( &
       & DelTime, timeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt &
       & )

    ! モジュール引用; Use statements
    !

    use VariableSet_mod, only: &
         & xyz_UB, xyz_UN, xyz_UA, &
         & xyz_VB, xyz_VN, xyz_VA, &
         & xyz_PTempEddB, xyz_PTempEddN, xyz_PTempEddA, &
         & xyz_SaltB, xyz_SaltN, xyz_SaltA, &
         & xy_SurfHeightB, xy_SurfHeightN, xy_SurfHeightA, &
         & xy_SurfPressB, xy_SurfPressN, xy_SurfPressA, &
         & xy_totDepthBasic, &
         & xyz_SigDot, z_PTempBasic, &
         & xyz_VViscCoefA, xyz_VViscCoefB, xyz_VViscCoefN, &
         & xyz_VDiffCoefA, xyz_VDiffCoefB, xyz_vDiffCoefN

    use BoundaryCondO_mod, only: &
         & xy_SurfHFlxO
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTime
    integer, intent(in) :: timeIntMode
    integer, intent(in) :: nStage_BarocTimeInt
    logical, intent(in) :: isVarBUsed_BarocTimeInt

    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !

    xyz_PTempEddA(:,:,0) = xyz_PTempEddN(:,:,0) &
         - DelTime / (Cp0 * RefDens * xy_totDepthBasic )*xy_SurfHFlxO
    
    xyz_UA = 0d0
    xyz_VA = 0d0
    xyz_SaltA = 0d0
    xy_SurfHeightA = 0d0
    
  end subroutine NonDynMixedLyr_AdvanceTStep
  
end module NonDynMixedLyr_TimeInteg_mod
