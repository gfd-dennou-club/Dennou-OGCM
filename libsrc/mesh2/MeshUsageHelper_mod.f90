!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module MeshUsageHelper_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING
  use dc_message, only: MessageNotify

  use PolyMesh_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  type, public :: MeshUsage
     integer :: usageTypeId
     character(STRING) :: usageTypeName
  end type MeshUsage

  type, public :: MeshUsageHelper
     type(MeshUsage), pointer :: usageList(:) => null()
     type(PolyMesh), pointer :: mesh => null()
  end type MeshUsageHelper



  public :: MeshUsageHelper_Init, MeshUsageHelper_Final
  public :: MeshUsageHelper_GetUsage, MeshUsageHelper_SetUsage

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'MeshUsageHelper_mod' !< Module Name

contains
  

  subroutine MeshUsage_Init(self, &
       & id, usageTypeName )
    type(MeshUsage), intent(inout) :: self
    integer, intent(in) :: id
    character(*), intent(in) :: usageTypeName
    
    self%usageTypeId = id
    self%usageTypeName = usageTypeName
  end subroutine MeshUsage_Init

  !>
  !!
  !!
  subroutine MeshUsageHelper_Init(self, mesh, nUsageType)

    type(MeshUsageHelper), intent(inout) :: self
    type(PolyMesh), intent(in), target :: mesh
    integer, intent(in) :: nUsageType

    ! 実行文; Executable statements
    !

    call MeshUsageHelper_Final(self)

    self%mesh => mesh
    allocate(self%usageList(nUsageType))
    
  end subroutine MeshUsageHelper_Init

  !>
  !!
  !!
  subroutine MeshUsageHelper_Final(self)

    type(MeshUsageHelper), intent(inout) :: self

    ! 実行文; Executable statements
    !

    if(associated(self%usageList)) deallocate(self%usageList)
    self%mesh => null()

  end subroutine MeshUsageHelper_Final

  !> @brief 
  !!
  !!
  subroutine MeshUsageHelper_SetUsage(self, id, usage)
    
    ! 宣言文; Declaration statement
    !
    type(MeshUsageHelper), intent(inout) :: self
    integer, intent(in) :: id
    type(MeshUsage), intent(in) :: usage
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    if(id > size(self%usageList)) then
       call MessageNotify('E', module_name, &
            & "Specified id of MeshUsage is invalid")
    end if
    
    self%usageList(id) = usage
    
  end subroutine MeshUsageHelper_SetUsage

  !> @brief 
  !!
  !! @return 
  !!
  function MeshUsageHelper_GetUsage(self, id) result(usage)
    
    ! 宣言文; Declaration statement
    !
    type(MeshUsageHelper), intent(in) :: self
    integer, intent(in) :: id
    type(MeshUsage) :: usage

    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    if(id > size(self%usageList)) then
       call MessageNotify('E', module_name, &
            & "Specified id of MeshUsage is invalid")
    end if

    usage = self%usageList(id)

  end function MeshUsageHelper_GetUsage

  
  
end module MeshUsageHelper_mod

