!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DiagVarSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  use GridSet_mod, only: &
       & iMax, jMax, kMax


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DiagVarSet_Init, DiagVarSet_Final

  ! 公開変数
  ! Public variable
  !
  real(DP), dimension(:,:,:), allocatable, public :: &
       & xyz_Div, xyz_Vor, xyz_BarocPress, xyz_PressEdd, xyz_DensEdd

  real(DP), dimension(:,:), allocatable, public :: &
       & xy_totDepth, yz_MassStreamFunc

  character(*), parameter, public :: DVARKEY_VOR = 'Vor'
  character(*), parameter, public :: DVARKEY_DIV = 'Div'
  character(*), parameter, public :: DVARKEY_PRESSEDD = 'PressEdd'
  character(*), parameter, public :: DVARKEY_DENSEDD = 'DensEdd'
  character(*), parameter, public :: DVARKEY_MASSSTREAMFUNC = 'MassStreamFunc'
  character(*), parameter, public :: DVARKEY_PTEMP = 'PTemp'
  character(*), parameter, public :: DVARKEY_STATICSTABILITY = 'StaticStability'
  

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DiagVarSet_mod' !< Module Name


contains

  !>
  !!
  !!
  subroutine DiagVarSet_Init( diagVarsName )

    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: diagVarsName(:)

    ! 作業変数
    ! Work variables
    !
    integer :: varID

    ! 実行文; Executable statements
    !
    allocate(xy_totDepth(0:iMax-1,jMax))
    allocate(xyz_BarocPress(0:iMax-1,jMax,0:kMax))
    allocate(xyz_PressEdd(0:iMax-1,jMax,0:kMax))
    allocate(xyz_DensEdd(0:iMax-1,jMax,0:kMax))

    do varID=1, size(diagVarsName)
       write(*,*) diagVarsName(varID)

       select case( diagVarsName(varID) )
          case( DVARKEY_DIV )
             allocate( xyz_Div(0:iMax-1,jMax,0:kMax) )
          case( DVARKEY_VOR )
             allocate( xyz_Vor(0:iMax-1,jMax,0:kMax) )
          case ( DVARKEY_PRESSEDD )
          case ( DVARKEY_DENSEDD )
          case (DVARKEY_MASSSTREAMFUNC )
             allocate( yz_MassStreamFunc(jMax, 0:kMax) )
          case (DVARKEY_PTEMP)
          case (DVARKEY_STATICSTABILITY)
          case Default
             call MessageNotify('E', module_name, &
                  & "The name of specified diagnostic variable '%c'is invalid.", c1=trim(diagVarsName(varID)) )
       end select
    end do
  end subroutine DiagVarSet_Init

  !>
  !!
  !!
  subroutine DiagVarSet_Final()

    ! 実行文; Executable statements
    !

    deallocate( xy_totDepth, xyz_BarocPress, xyz_PressEdd, xyz_DensEdd ) 
    if( allocated(xyz_Div) ) deallocate(xyz_Div)
    if( allocated(xyz_Vor) ) deallocate(xyz_Vor)
    if( allocated(xyz_BarocPress) ) deallocate(xyz_BarocPress)
    if( allocated(xyz_PressEdd) ) deallocate(xyz_PressEdd)
    if( allocated(yz_MassStreamFunc) ) deallocate(yz_MassStreamFunc)

  end subroutine DiagVarSet_Final

end module DiagVarSet_mod
