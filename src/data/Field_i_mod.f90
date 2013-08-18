module Field_i_mod
  !
  use dc_types
  use VectorSpace_mod
  !
  implicit none
  private
  !
  type, public :: Field_i
     integer, pointer :: v_(:)
     character(TOKEN) :: name
     character(TOKEN) :: units
  end type Field_i
  interface setAtitude
     module procedure Field_i_setAtitude
  end interface setAtitude
  ! Public procedures
  public :: Field_i_Init, Field_i_Final, setAtitude
contains
!
!
subroutine Field_i_Init (list, list_size)
  type(Field_i), intent(inout) :: list
  integer, intent(in) :: list_size
  allocate(list%v_(list_size))
end subroutine Field_i_Init
!
!
!
subroutine Field_i_setAtitude (list, name, units)
  type(Field_i), intent(inout) :: list
  character(*), intent(in) :: name
  character(*), intent(in) :: units
  list%name = name
  list%units = units
end subroutine Field_i_setAtitude
!
!
subroutine Field_i_Final (list)
  type(Field_i), intent(inout) :: list
  deallocate(list%v_)
end subroutine Field_i_Final
end module Field_i_mod
