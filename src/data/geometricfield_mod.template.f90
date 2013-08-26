#define _eval(s) s
#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _moduleName(s) s ## _mod
#define moduleName _eval1(_moduleName, FieldTypeName)
#define _funcNameSuffix(pref, suff) pref ## _ ## suff
#define funcNameSuffix(suff) _eval2(_funcNameSuffix, FieldTypeName, suff)

module moduleName
  
  !
  use dc_types
  use VectorSpace_mod
  use List_mod
  use PolyMesh_mod

  !
  implicit none
  private

  !
  type, public :: FieldTypeName
    type(FieldDataType), pointer :: data => null()
    type(PolyMesh), pointer :: mesh => null()
    character(TOKEN) :: name
    character(STRING) :: long_name
    character(TOKEN) :: units
    logical :: tempDataFlag = .false.
  end type FieldTypeName

  !
  !
  interface GeometricField_Init
    module procedure funcNameSuffix(Init)
  end interface GeometricField_Init

  !
  interface GeometricField_Final
    module procedure funcNameSuffix(Final)
  end interface GeometricField_Final

  interface SetFieldAtitude
    module procedure funcNameSuffix(SetFieldAtitude)
  end interface SetFieldAtitude
  
  !
  interface assignment(=)
    module procedure funcNameSuffix(assignFieldObj)
    module procedure funcNameSuffix(assignFDElemType)
    module procedure funcNameSuffix(assignFDElemTypeArray)
  end interface assignment(=)

!!$  interface assignment(:=)
!!$    module procedure funcNameSuffix(assignForExpr)
!!$  end interface assignment(:=)
  
  interface operator(+)
    module procedure funcNameSuffix(addOptr1)
!    module procedure funcNameSuffix(addOptr2)
!    module procedure funcNameSuffix(addOptr3)
  end interface operator(+)

  interface operator(-)
    module procedure funcNameSuffix(subOptr1)
!    module procedure funcNameSuffix(subOptr2)
!    module procedure funcNameSuffix(subOptr3)
  end interface operator(-)
  
  interface operator(.At.)
    module procedure funcNameSuffix(ReturnElemRef)
  end interface operator(.At.)
  
  !
  public :: GeometricField_Init, GeometricField_Final
  public :: SetFieldAtitude
  public :: operator(+), operator(-), assignment(=), operator(.At.)

  
  character(TOKEN), parameter :: DATADEFLOC_SURFACE = "surface"
  character(TOKEN), parameter :: DATADEFLOC_POINT = "point"
  character(TOKEN), parameter :: DATADEFLOC_VOL = "vol"

contains

!
subroutine funcNameSuffix(Init) (field, mesh, name, long_name, units)
  type(FieldTypeName), intent(inout) :: field
  type(PolyMesh), target, intent(in) :: mesh
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units

  ! Work variables
  !
  integer :: dataSize
  
  ! Executable statements
  !

  !
  !

  !
  call GeometricField_Final(field)
  
  !
  field%mesh => mesh
  call setFieldAtitude(field, "field", " ", "1")

  if(present(name)) call setFieldAtitude(field, name=name)
  if(present(long_name)) call setFieldAtitude(field, long_name=long_name)
  if(present(units)) call setFieldAtitude(field, units=units)

  !
  select case( DataDefLoc )
    case (DATADEFLOC_SURFACE)
      dataSize = getFaceListSize(field%mesh)
    case (DATADEFLOC_VOL)
      dataSize = getCellListSize(field%mesh)
    case (DATADEFLOC_POINT)
      dataSize = getPointListSize(field%mesh)
  end select

  !
  allocate(field%data)
  call List_Init(field%data, dataSize)
  write(*,*) "allocated..", ", dataSize=", dataSize

end subroutine funcNameSuffix(Init)

subroutine funcNameSuffix(SetFieldAtitude) (field, name, long_name, units)
  type(FieldTypeName), intent(inout) :: field
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units

  if(present(name)) field%name = name
  if(present(long_name)) field%long_name = long_name
  if(present(units)) field%units = units

end subroutine funcNameSuffix(SetFieldAtitude)

!
subroutine funcNameSuffix(Final) (field)
  type(FieldTypeName), intent(inout) :: field

  write(*,*) "Call releasedataRef"
  call funcNameSuffix(releaseDataRef) (field)
  
end subroutine funcNameSuffix(Final)


!
subroutine funcNameSuffix(releaseDataRef) (field)
  type(FieldTypeName) :: field
  
  logical :: refIs0

  if( .not. associated(field%data) ) return

  call decRef(field%data, refIs0)
  
  if( refIs0 ) then
    call List_Final(field%data)
    write(*,*) "FieldName=", field%name
    if( field%tempDataFlag ) write(*,*) "This object is a temporary data. ^"
    deallocate(field%data)
  end if

  field%data => null()
  
end subroutine funcNameSuffix(releaseDataRef)


! operation
!
subroutine funcNameSuffix(assignFieldObj) (this, other)
  type(FieldTypeName), intent(inout) :: this
  type(FieldTypeName), intent(in) :: other

!  call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  
  if( .not. associated(other%data) ) then
    write(*,*) "Error: try to assignment to uninitialized list data!" 
    stop
  end if

write(*,*) "In assignment.."
  !
  call funcNameSuffix(releaseDataRef) (this)
write(*,*) ".."

  this%data => other%data
  call incRef(this%data)

  if( other%tempDataFlag ) then
    call funcNameSuffix(releaseDataRef) (other)
  end if

end subroutine funcNameSuffix(assignFieldObj)

subroutine funcNameSuffix(assignFDElemType) (this, elemData)
  type(FieldTypeName), intent(inout) :: this
  FieldDataElemType, intent(in) :: elemData

!  call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  this%data%v_(:) = elemData

end subroutine funcNameSuffix(assignFDElemType)

subroutine funcNameSuffix(assignFDElemTypeArray) (this, elemArrayData)
  type(FieldTypeName), intent(inout) :: this
  FieldDataElemType, intent(in) :: elemArrayData(:)


  if( getListSize(this%data) /= size(elemArrayData) ) then
    write(*,*) "Error in assignment:"
  end if

  this%data%v_(:) = elemArrayData(:)

end subroutine funcNameSuffix(assignFDElemTypeArray)

function funcNameSuffix(ReturnElemRef) (this, id) result(elemRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: id
  FieldDataElemType :: elemRef

  elemRef = this%data%v_(id)

end function funcNameSuffix(ReturnElemRef)

!*
!* 

! Define add operators
#define op +
#define opname addOptr
#include "geometricfield_binaryOptr.template.f90"
#undef opname
#undef op


! Define sub operators
#define op -
#define opname subOptr
#include "geometricfield_binaryOptr.template.f90"
#undef op
#undef opname


!
subroutine checkDataSizeEquality(data1, data2, info)
  type(FieldDataType), intent(in) :: data1, data2
  character(*), intent(in) :: info

  if( getListSize(data1) /= getListSize(data2) ) then
    write(*,*) info, ": The size of list data in binary operation is not equal to other."
    stop
  end if
  
end subroutine checkDataSizeEquality


end module moduleName
