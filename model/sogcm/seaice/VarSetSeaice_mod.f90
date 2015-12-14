!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module VarSetSeaice_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  
  use GridSet_mod, only: &
       & iMax, jMax

  use SeaIceConstants_mod, only: &
       & Mu, SaltSeaIce
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: VarSetSeaice_Init, VarSetSeaice_Final
  public :: VarSetSeaice_SetDefualtValue
  public :: VarSetSeaice_AdvanceTStep
  
  ! 公開変数
  ! Public variable
  !
  real(DP), allocatable, public :: xy_SIceConA(:,:)
  real(DP), allocatable, public :: xy_SIceConN(:,:)
  real(DP), allocatable, public :: xy_SIceConB(:,:)
  
  real(DP), allocatable, public :: xy_SIceSurfTempA(:,:)
  real(DP), allocatable, public :: xy_SIceSurfTempN(:,:)
  real(DP), allocatable, public :: xy_SIceSurfTempB(:,:)

  real(DP), allocatable, public :: xyz_SIceTempA(:,:,:)
  real(DP), allocatable, public :: xyz_SIceTempN(:,:,:)
  real(DP), allocatable, public :: xyz_SIceTempB(:,:,:)

  real(DP), allocatable, public :: xy_IceThickA(:,:)
  real(DP), allocatable, public :: xy_IceThickN(:,:)
  real(DP), allocatable, public :: xy_IceThickB(:,:)
  
  real(DP), allocatable, public :: xy_SnowThickA(:,:)
  real(DP), allocatable, public :: xy_SnowThickN(:,:)
  real(DP), allocatable, public :: xy_SnowThickB(:,:)

  character(*), parameter, public :: VARSET_KEY_SICECON = 'SIceCon'
  character(*), parameter, public :: VARSET_KEY_SNOWTHICK = 'SnowThick'  
  character(*), parameter, public :: VARSET_KEY_ICETHICK = 'IceThick'
  character(*), parameter, public :: VARSET_KEY_SICETEMP = 'SIceTemp'
  character(*), parameter, public :: VARSET_KEY_SICESURFTEMP = 'SIceSurfTemp'  

  character(*), parameter, public :: VARSET_KEY_SICECONB = 'SIceConB'
  character(*), parameter, public :: VARSET_KEY_SNOWTHICKB = 'SnowThickB'  
  character(*), parameter, public :: VARSET_KEY_ICETHICKB = 'IceThickB'
  character(*), parameter, public :: VARSET_KEY_SICETEMPB = 'SIceTempB'
  character(*), parameter, public :: VARSET_KEY_SICESURFTEMPB = 'SIceSurfTempB'  
  

  real(DP), allocatable, public :: xy_Wice(:,:)
  
  ! 非公開手続き
  ! Private procedure
  !
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VarSetSeaice_mod' !< Module Name
  
contains

  !>
  !!
  !!
  subroutine VarSetSeaice_Init(nIceLyr)

    !
    !
    integer, intent(in) :: nIceLyr
    
    ! 実行文; Executable statements
    !

    call allocXY(xy_SIceConA, xy_SIceConN, xy_SIceConB)
    call allocXY(xy_SIceSurfTempA, xy_SIceSurfTempN, xy_SIceSurfTempB)
    call allocXY(xy_IceThickA, xy_IceThickN, xy_IceThickB)
    call allocXY(xy_SnowThickA, xy_SnowThickN, xy_SnowThickB)
    call allocXYZ(xyz_SIceTempA, xyz_SIceTempN, xyz_SIceTempB, nIceLyr)

    allocate( xy_Wice(0:iMax-1,jMax) )
    
  contains
    subroutine allocXY(xyA, xyN, xyB)
      real(DP), dimension(:,:), allocatable :: xyA, xyN, xyB

      allocate(xyA(0:iMax-1,jMax), xyN(0:iMax-1,jMax), xyB(0:iMax-1,jMax))
    end subroutine allocXY

    subroutine allocXYZ(xyzA, xyzN, xyzB, nZ)

      real(DP), dimension(:,:,:), allocatable :: xyzA, xyzN, xyzB
      integer, intent(in) :: nZ
      
      allocate( &
           & xyzA(0:iMax-1,jMax,nZ), xyzN(0:iMax-1,jMax,nz), xyzB(0:iMax-1,jMax,nz) &
           & )
    end subroutine allocXYZ
    
  end subroutine VarSetSeaice_Init

  !> @brief 
  !!
  !!
  subroutine VarSetSeaice_SetDefualtValue()
    
    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statement
    !

    xy_SIceConN = 0d0
    xy_SIceSurfTempN = 0d0
    xyz_SIceTempN = 0d0!- Mu*SaltSeaIce
    xyz_SIceTempA = 0d0!- Mu*SaltSeaIce    
    xy_SnowThickN = 0d0
    xy_IceThickN = 0d0

    xy_Wice = 0d0
    
  end subroutine VarSetSeaice_SetDefualtValue

  !>
  !!
  !!
  subroutine VarSetSeaice_Final()

    ! 実行文; Executable statements
    !

    deallocate( xy_SIceConA, xy_SIceConN, xy_SIceConB )
    deallocate( xy_SIceSurfTempA, xy_SIceSurfTempN, xy_SIceSurfTempB )
    deallocate( xyz_SIceTempA, xyz_SIceTempN, xyz_SIceTempB )
    deallocate( xy_SnowThickA, xy_SnowThickN, xy_SnowThickB )
    deallocate( xy_IceThickA, xy_IceThickN, xy_IceThickB )

    deallocate( xy_Wice )
    
  end subroutine VarSetSeaice_Final


  !> @brief 
  !!
  !!
  subroutine VarSetSeaice_AdvanceTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp workshare
    xy_SIceConB = xy_SIceConN; xy_SIceConN = xy_SIceConA; xy_SIceConA = 0d0

    xy_IceThickB = xy_IceThickN; xy_IceThickN = xy_IceThickA; xy_IceThickA = 0d0
    xy_SnowThickB = xy_SnowThickN; xy_SnowThickN = xy_SnowThickA; xy_SnowThickA = 0d0

    xy_SIceSurfTempB = xy_SIceSurfTempN; xy_SIceSurfTempN = xy_SIceSurfTempA; xy_SIceSurfTempA = 0d0
    xyz_SIceTempB = xyz_SIceTempN; xyz_SIceTempN = xyz_SIceTempA; xyz_SIceTempA = 0d0
    !$omp end workshare
    !$omp end parallel
    
  end subroutine VarSetSeaice_AdvanceTStep

end module VarSetSeaice_mod

