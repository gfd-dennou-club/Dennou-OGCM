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
       & xyz_Div, xyz_Vor, xyz_PressBaroc, xyz_Press

  character(*), parameter, public :: DVARKEY_VOR = 'Vor'
  character(*), parameter, public :: DVARKEY_DIV = 'Div'
  character(*), parameter, public :: DVARKEY_PRESSBAROC = 'PressBaroc'
  character(*), parameter, public :: DVARKEY_PRESS = 'Press'
  character(*), parameter, public :: DVARKEY_DENS = 'DENS'

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
    
    do varID=1, size(diagVarsName)
       write(*,*) diagVarsName(varID)

       select case( diagVarsName(varID) )
          case( DVARKEY_DIV )
             allocate( xyz_Div(0:iMax-1,jMax,0:kMax) )
          case( DVARKEY_VOR )
             allocate( xyz_Vor(0:iMax-1,jMax,0:kMax) )
          case( DVARKEY_PRESSBAROC )
             allocate( xyz_PressBaroc(0:iMax-1,jMax,0:kMax) )
          case( DVARKEY_PRESS )
             allocate( xyz_Press(0:iMax-1,jMax,0:kMax) )
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

    if( allocated(xyz_Div) ) deallocate(xyz_Div)
    if( allocated(xyz_Vor) ) deallocate(xyz_Vor)
    if( allocated(xyz_PressBaroc) ) deallocate(xyz_PressBaroc)
    if( allocated(xyz_Press) ) deallocate(xyz_Press)


  end subroutine DiagVarSet_Final

end module DiagVarSet_mod
