!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module VariableSet_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, STRING, TOKEN 

  !* Dennou-OGCM

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax

  use TemporalIntegSet_mod, only: &
       & nLongTimeLevel

  use BoundCondSet_mod, only: &
       & ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, &
       & ThermBCTYPE_PrescTemp, ThermBCTYPELBL_TempRelaxed, &
       & ThermBC_Surface, &
       & SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, &
       & SaltBCTYPE_PrescSalt, SaltBCTYPELBL_SaltRelaxed, &
       & SaltBC_Surface


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: VariableSet_Init, VariableSet_Final
  public :: VariableSet_AdvanceTStep

  ! 公開変数
  ! Public variable
  !

  integer, public, save :: TracerNum
  integer, public, save :: SaltTracerID
  integer, public, save :: PTempTracerID

  real(DP), public, save, allocatable :: xyz_uA(:,:,:), xyz_uN(:,:,:), xyz_uB(:,:,:)
  real(DP), public, save, allocatable :: xyz_vA(:,:,:), xyz_vN(:,:,:), xyz_vB(:,:,:)
  real(DP), public, save, allocatable :: xyz_SigDot(:,:,:)
  real(DP), public, save, allocatable :: xyz_PTempEddA(:,:,:), xyz_PTempEddN(:,:,:), xyz_PTempEddB(:,:,:)
  real(DP), public, save, allocatable :: xyz_SaltA(:,:,:), xyz_SaltN(:,:,:), xyz_SaltB(:,:,:)
  real(DP), public, save, allocatable :: xy_SurfHeightA(:,:), xy_SurfHeightN(:,:), xy_SurfHeightB(:,:)
  real(DP), public, save, allocatable :: xy_totDepthBasic(:,:)
  real(DP), public, save, allocatable :: z_PTempBasic(:)
  real(DP), public, save, allocatable :: xy_SurfPressA(:,:), xy_SurfPressN(:,:), xy_SurfPressB(:,:)

  real(DP), public, save, allocatable :: xyz_ConvIndex(:,:,:)
  real(DP), public, save, allocatable :: xyz_VViscCoefA(:,:,:), xyz_VViscCoefN(:,:,:), xyz_VViscCoefB(:,:,:)
  real(DP), public, save, allocatable :: xyz_VDiffCoefA(:,:,:), xyz_VDiffCoefN(:,:,:), xyz_VDiffCoefB(:,:,:)

!!$  real(DP), public, save, allocatable :: wz_VorA(:,:), wz_VorN(:,:), wz_VorB(:,:)
!!$  real(DP), public, save, allocatable :: wz_DivA(:,:), wz_DivN(:,:), wz_DivB(:,:)
!!$  real(DP), public, save, allocatable :: wz_SaltA(:,:), wz_SaltN(:,:), wz_SaltB(:,:)
!!$  real(DP), public, save, allocatable :: wz_PTempEddA(:,:), wz_PTempEddN(:,:), wz_PTempEddB(:,:)
  
  character(TOKEN), public, parameter :: VARSET_KEY_U  = 'U'
  character(TOKEN), public, parameter :: VARSET_KEY_UB = 'UB'
  character(TOKEN), public, parameter :: VARSET_KEY_V  = 'V'
  character(TOKEN), public, parameter :: VARSET_KEY_VB = 'VB'
  character(TOKEN), public, parameter :: VARSET_KEY_PTEMPEDD  = 'PTempEdd'
  character(TOKEN), public, parameter :: VARSET_KEY_PTEMPEDDB = 'PTempEddB'
  character(TOKEN), public, parameter :: VARSET_KEY_SALT  = 'Salt'
  character(TOKEN), public, parameter :: VARSET_KEY_SALTB = 'SaltB'
  character(TOKEN), public, parameter :: VARSET_KEY_SURFHEIGHT = 'SurfHeight'
  character(TOKEN), public, parameter :: VARSET_KEY_SURFHEIGHTB = 'SurfHeightB'

  character(TOKEN), public, parameter :: VARSET_KEY_SIGDOT = 'SigDot'
  character(TOKEN), public, parameter :: VARSET_KEY_SURFPRESS = 'SurfPress'
  character(TOKEN), public, parameter :: VARSET_KEY_HYDROPRESSEDD = 'HydroPressEdd'

  character(TOKEN), public, parameter :: VARSET_KEY_TOTDEPTHBASIC = 'TotDepthBasic'
  character(TOKEN), public, parameter :: VARSET_KEY_PTEMPBASIC = 'PTempBasic'


  character(TOKEN), public, parameter :: VARSET_KEY_CONVINDEX = 'ConvIndex'

  character(TOKEN), public, parameter :: VARSET_KEY_VVISCCOEF = 'VViscCoef'
  character(TOKEN), public, parameter :: VARSET_KEY_VDIFFCOEF = 'VDiffCoef'
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VariableSet_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine VariableSet_Init()

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !

    ! 実行文; Executable statements
    !

    call malloc3DVar(xyz_UA); call malloc3DVar(xyz_UN); call malloc3DVar(xyz_UB); 
    call malloc3DVar(xyz_VA); call malloc3DVar(xyz_VN); call malloc3DVar(xyz_VB);
    call malloc3DVar(xyz_SigDot)
    call malloc3DVar(xyz_PTempEddA); call malloc3DVar(xyz_PTempEddN); call malloc3DVar(xyz_PTempEddB); 
    call malloc3DVar(xyz_SaltA); call malloc3DVar(xyz_SaltN); call malloc3DVar(xyz_SaltB); 
    call malloc2DVar(xy_SurfHeightA); call malloc2DVar(xy_SurfHeightN); call malloc2DVar(xy_SurfHeightB);
    call malloc2DVar(xy_totDepthBasic)
    call malloc1DVar(z_PTempBasic)
    call malloc2DVar(xy_SurfPressA);  call malloc2DVar(xy_SurfPressN); call malloc2DVar(xy_SurfPressB);

    call malloc3DVar(xyz_VViscCoefA); call malloc3DVar(xyz_VViscCoefN); call malloc3DVar(xyz_VViscCoefB);
    call malloc3DVar(xyz_VDiffCoefA); call malloc3DVar(xyz_VDiffCoefN); call malloc3DVar(xyz_VDiffCoefB);
    
    call malloc3DVar(xyz_ConvIndex)

!!$    call malloc2DVar_wz(wz_VorA); call malloc2DVar_wz(wz_VorN); call malloc2DVar_wz(wz_VorB)
!!$    call malloc2DVar_wz(wz_DivA); call malloc2DVar_wz(wz_DivN); call malloc2DVar_wz(wz_DivB)
!!$    call malloc2DVar_wz(wz_PTempEddA); call malloc2DVar_wz(wz_PTempEddN); call malloc2DVar_wz(wz_PTempEddB)
!!$    call malloc2DVar_wz(wz_SaltA); call malloc2DVar_wz(wz_SaltN); call malloc2DVar_wz(wz_SaltB)
    
    !
    TracerNum = 2
    SaltTracerID = 1
    PTempTracerID = 2
    
  contains
    subroutine malloc3DVar( var )
      real(DP), intent(inout), allocatable :: var(:,:,:)
      allocate( var(0:iMax-1, 1:jMax, 0:kMax) )
    end subroutine malloc3DVar
    subroutine malloc2DVar( var )
      real(DP), intent(inout), allocatable :: var(:,:)
      allocate( var(0:iMax-1, 1:jMax) )
    end subroutine malloc2DVar
    subroutine malloc2DVar_wz( var )
      real(DP), intent(inout), allocatable :: var(:,:)
      allocate( var(lMax, 0:kMax) )
    end subroutine malloc2DVar_wz    
    subroutine malloc1DVar( var )
      real(DP), intent(inout), allocatable :: var(:)
      allocate( var(0:kMax) )
    end subroutine malloc1DVar
  end subroutine VariableSet_Init

  !>
  !!
  !!
  subroutine VariableSet_Final()

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variable
    !
    
    ! 実行文; Executable statements
    !

    if( allocated(xyz_UA) ) then
       deallocate( xyz_UA, xyz_UN, xyz_UB )
       deallocate( xyz_VA, xyz_VN, xyz_VB )
       deallocate( xyz_SigDot )
       deallocate( xyz_PTempEddA, xyz_PTempEddN, xyz_PTempEddB )
       deallocate( xyz_SaltA, xyz_SaltN, xyz_SaltB )
       deallocate( xy_SurfHeightA, xy_SurfHeightN, xy_SurfHeightB )
       deallocate( xy_totDepthBasic )
       deallocate( xy_SurfPressA, xy_SurfPressN, xy_SurfPressB )
       deallocate( xyz_ConvIndex )
       deallocate( xyz_VViscCoefA, xyz_VViscCoefN, xyz_VViscCoefB )
       deallocate( xyz_VDiffCoefA, xyz_VDiffCoefN, xyz_VDiffCoefB )

!!$       deallocate( wz_VorA, wz_VorN, wz_VorB )
!!$       deallocate( wz_DivA, wz_DivN, wz_DivB )
!!$       deallocate( wz_PTempEddA, wz_PTempEddN, wz_PTempEddB )
!!$       deallocate( wz_SaltA, wz_SaltN, wz_SaltB )
    end if

  end subroutine VariableSet_Final


  !> @brief 
  !!
  !!
  subroutine VariableSet_AdvanceTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    integer :: j, k
    
    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp workshare
    xyz_UB(:,:,:) = xyz_UN; xyz_UN(:,:,:) = xyz_UA; xyz_UA(:,:,:) = 0d0
    xyz_VB(:,:,:) = xyz_VN; xyz_VN(:,:,:) = xyz_VA; xyz_VA(:,:,:) = 0d0
    xyz_PTempEddB(:,:,:) = xyz_PTempEddN; xyz_PTempEddN(:,:,:) = xyz_PTempEddA; xyz_PTempEddA(:,:,:) = 0d0
    xyz_SaltB(:,:,:) = xyz_SaltN; xyz_SaltN(:,:,:) = xyz_SaltA; xyz_SaltA(:,:,:) = 0d0

    xyz_VViscCoefB(:,:,:) = xyz_VViscCoefN; !xyz_VViscCoefN(:,:,:) = xyz_VViscCoefA; xyz_VViscCoefA(:,:,:) = 0d0 
    xyz_VDiffCoefB(:,:,:) = xyz_VDiffCoefN; !xyz_VDiffCoefN(:,:,:) = xyz_VDiffCoefA; xyz_VDiffCoefA(:,:,:) = 0d0 
    !$omp end workshare

    !$omp workshare
    xy_SurfHeightB(:,:) = xy_SurfHeightN; xy_SurfHeightN(:,:) = xy_SurfHeightA; xy_SurfHeightA(:,:) = 0d0
    xy_SurfPressB(:,:) = xy_SurfPressN; xy_SurfPressN(:,:) = xy_SurfPressA; xy_SurfPressA(:,:) = 0d0
    !$omp end workshare
    !$omp end parallel

!!$    !$omp parallel workshare
!!$    wz_VorB(:,:) = wz_VorN; wz_VorN(:,:) = wz_VorA; wz_VorA(:,:) = 0d0
!!$    wz_DivB(:,:) = wz_VorN; wz_DivN(:,:) = wz_DivA; wz_DivA(:,:) = 0d0
!!$    wz_PTempEddB(:,:) = wz_PTempEddN; wz_PTempEddN(:,:) = wz_PTempEddA; wz_PTempEddA(:,:) = 0d0
!!$    wz_SaltB(:,:) = wz_SaltN; wz_SaltN(:,:) = wz_SaltA; wz_SaltA(:,:) = 0d0
!!$    !$omp end parallel workshare

  end subroutine VariableSet_AdvanceTStep

end module VariableSet_mod

