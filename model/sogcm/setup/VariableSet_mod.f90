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

  use dc_types, only: &
       & DP, STRING 

  use TemporalIntegSet_mod, only: &
       & nLongTimeLevel

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
  real(DP), public, save, allocatable :: xy_SurfPress(:,:)
  real(DP), public, save, allocatable :: xy_WindStressU(:,:)
  real(DP), public, save, allocatable :: xy_WindStressV(:,:)

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
    use PolyMesh_mod, only: &
         & PolyMesh

    use GridSet_mod, only: &
         & iMax, jMax, kMax

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
    call malloc2DVar(xy_SurfPress)
    call malloc2DVar(xy_WindStressU); call malloc2DVar(xy_WindStressV)
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
       deallocate( xy_SurfPress )
       deallocate( xy_WindStressU, xy_WindStressV )
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
    
    
    ! 実行文; Executable statement
    !

    xyz_UB = xyz_UN; xyz_UN = xyz_UA; xyz_UA = 0d0
    xyz_VB = xyz_VN; xyz_VN = xyz_VA; xyz_VA = 0d0
    xyz_PTempEddB = xyz_PTempEddN; xyz_PTempEddN = xyz_PTempEddA; xyz_PTempEddA = 0d0
    xyz_SaltB = xyz_SaltN; xyz_SaltN = xyz_SaltA; xyz_SaltA = 0d0
    xy_SurfHeightB = xy_SurfHeightN; xy_SurfHeightN = xy_SurfHeightA; xy_SurfHeightA = 0d0

  end subroutine VariableSet_AdvanceTStep

end module VariableSet_mod

