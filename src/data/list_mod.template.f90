#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _moduleName(s) s ## _mod
#define moduleName _eval1(_moduleName, ListTypeName)
#define _funcNameSuffix(pref, suff) pref ## _ ## suff
#define funcNameSuffix(suff) _eval2(_funcNameSuffix, ListTypeName, suff)

module moduleName
  
  !
  use dc_types
  use VectorSpace_mod

  !
  implicit none
  private

  !
  type, public :: ListTypeName
     ListElemType, pointer :: v_(:)
     character(TOKEN) :: name
     character(TOKEN) :: units
  end type ListTypeName

  interface setAtitude
     module procedure funcNameSuffix(setAtitude)
  end interface setAtitude
  
  ! Public procedures
  public :: funcNameSuffix(Init), funcNameSuffix(Final), setAtitude

  
contains

!
!
subroutine funcNameSuffix(Init) (list, list_size)
  type(ListTypeName), intent(inout) :: list
  integer, intent(in) :: list_size
  
  allocate(list%v_(list_size))

end subroutine funcNameSuffix(Init)

!
!
!
subroutine funcNameSuffix(setAtitude) (list, name, units)  
  type(ListTypeName), intent(inout) :: list
  character(*), intent(in) :: name
  character(*), intent(in) :: units

  list%name = name
  list%units = units

end subroutine funcNameSuffix(setAtitude)

!
!
subroutine funcNameSuffix(Final) (list)
  type(ListTypeName), intent(inout) :: list

  deallocate(list%v_)

end subroutine funcNameSuffix(Final)

end module moduleName

