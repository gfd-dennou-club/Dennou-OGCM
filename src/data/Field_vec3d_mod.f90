module Field_vec3d_mod
  !
  use dc_types
  use VectorSpace_mod
  !
  implicit none
  private
  !
  type, public :: Field_vec3d
     type(vector3d), pointer :: v_(:)
     character(TOKEN) :: name
     character(TOKEN) :: units
  end type Field_vec3d
  interface setAtitude
     module procedure Field_vec3d_setAtitude
  end interface setAtitude
  ! Public procedures
  public :: Field_vec3d_Init, Field_vec3d_Final, setAtitude
contains
!
!
subroutine Field_vec3d_Init (list, list_size)
  type(Field_vec3d), intent(inout) :: list
  integer, intent(in) :: list_size
  allocate(list%v_(list_size))
end subroutine Field_vec3d_Init
!
!
!
subroutine Field_vec3d_setAtitude (list, name, units)
  type(Field_vec3d), intent(inout) :: list
  character(*), intent(in) :: name
  character(*), intent(in) :: units
  list%name = name
  list%units = units
end subroutine Field_vec3d_setAtitude
!
!
subroutine Field_vec3d_Final (list)
  type(Field_vec3d), intent(inout) :: list
  deallocate(list%v_)
end subroutine Field_vec3d_Final
end module Field_vec3d_mod
