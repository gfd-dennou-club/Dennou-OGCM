#define _eval(s) s
#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _moduleName(s) s ## _mod
#define moduleName _eval1(_moduleName, FieldTypeName)
#define _funcNameSuffix(pref, suff) pref ## _ ## suff
#define funcNameSuffix(suff) _eval2(_funcNameSuffix, FieldTypeName, suff)

#if DataTypeId == 1
#define DATATYPE_SCALAR
#undef DATATYPE_VECTOR
#elif DataTypeId == 2
#define DATATYPE_VECTOR
#undef DATATYPE_SCALAR 
#endif

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
    integer :: vlayerNum
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
  
  interface Release
     module procedure funcNameSuffix(releaseDataRef)
  end interface Release

  interface DeepCopy
     module procedure funcNameSuffix(DeepCopy)
  end interface DeepCopy

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
    module procedure funcNameSuffix(addOptr2)
    module procedure funcNameSuffix(addOptr3)
  end interface operator(+)

  interface operator(-)
    module procedure funcNameSuffix(subOptr1)
    module procedure funcNameSuffix(subOptr2)
    module procedure funcNameSuffix(subOptr3)
  end interface operator(-)
  
  interface operator(*)
#ifdef DATATYPE_SCALAR
     module procedure funcNameSuffix(mulOptr1)
#endif
     module procedure funcNameSuffix(mulOptr2)
  end interface operator(*)

  interface operator(/)
#ifdef DATATYPE_SCALAR
     module procedure funcNameSuffix(divOptr1)
     module procedure funcNameSuffix(divOptr2)
#endif
     module procedure funcNameSuffix(divOptr3)
  end interface operator(/)

  interface At
     module procedure funcNameSuffix(At1)
     module procedure funcNameSuffix(At2)
  end interface At
  
  interface hSlice
     module procedure funcNameSuffix(hSlice1)
     module procedure funcNameSuffix(hSlice2)
  end interface hSlice

  interface vSlice
     module procedure funcNameSuffix(vSlice1)
     module procedure funcNameSuffix(vSlice2)
  end interface vSlice

  !
  public :: GeometricField_Init, GeometricField_Final, Release
  public :: SetFieldAtitude, DeepCopy
  public :: operator(+), operator(/), operator(-), operator(*), assignment(=)
  public :: At, hSlice, vSlice

  
  character(TOKEN), parameter :: DATADEFLOC_SURFACE = "surface"
  character(TOKEN), parameter :: DATADEFLOC_POINT = "point"
  character(TOKEN), parameter :: DATADEFLOC_VOL = "vol"

contains

!
subroutine funcNameSuffix(Init) (field, mesh, name, long_name, units, vlayerNum)
  type(FieldTypeName), intent(inout) :: field
  type(PolyMesh), target, intent(in) :: mesh
  character(*), optional, intent(in) :: name
  character(*), optional, intent(in) :: long_name
  character(*), optional, intent(in) :: units
  integer, optional, intent(in) :: vlayerNum

  ! Work variables
  !
  integer :: hdataSize
  
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
  
  if(present(vlayerNum)) then
     field%vlayerNum = vlayerNum
  else
     field%vlayerNum = mesh%vlayerNum
  end if

  !
  select case( DataDefLoc )
    case (DATADEFLOC_SURFACE)
       hdataSize = getFaceListSize(mesh)
    case (DATADEFLOC_VOL)
      hdataSize = getCellListSize(mesh)
    case (DATADEFLOC_POINT)
      hdataSize = getPointListSize(mesh)
  end select

  !
  allocate(field%data)
  call List_Init(field%data, field%vLayerNum, hdataSize)
!  write(*,*) field%name, "allocated..", ", hdataSize=", hdataSize, "vlayerNum=", field%vlayerNum

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

  !write(*,*) "Call releasedataRef"
  call funcNameSuffix(releaseDataRef) (field)
  
end subroutine funcNameSuffix(Final)


!
subroutine funcNameSuffix(DeepCopy) (self, field)
  type(FieldTypeName), intent(inout) :: self
  type(FieldTypeName), intent(in) :: field
  

  self%data%v_(:,:) = field%data%v_(:,:) 
  
end subroutine funcNameSuffix(DeepCopy)


!
subroutine funcNameSuffix(releaseDataRef) (field)
  type(FieldTypeName) :: field
  
  logical :: refIs0

  if( .not. associated(field%data) ) return

  call decRef(field%data, refIs0)
  
  if( refIs0 ) then
    call List_Final(field%data)
    !write(*,*) "FieldName=", field%name
    !if( field%tempDataFlag ) write(*,*) "This object is a temporary data. ^"
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

  !
  call funcNameSuffix(releaseDataRef) (this)

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
  this%data%v_ = elemData

end subroutine funcNameSuffix(assignFDElemType)

subroutine funcNameSuffix(assignFDElemTypeArray) (this, elemArrayData)
  type(FieldTypeName), intent(inout) :: this
  FieldDataElemType, intent(in) :: elemArrayData(:,:)


!!$  if( getListSize(this%data) /= size(elemArrayData) ) then
!!$    write(*,*) "Error in assignment:"
!!$  end if

  this%data%v_(:,:) = elemArrayData(:,:)

end subroutine funcNameSuffix(assignFDElemTypeArray)

pure function funcNameSuffix(At1) (this, k, hid) result(elemRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: hid, k
  FieldDataElemType :: elemRef

  elemRef = this%data%v_(k,hid)

end function funcNameSuffix(At1)

pure function funcNameSuffix(At2) (this, hid) result(elemRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: hid
  FieldDataElemType :: elemRef

  elemRef = this%data%v_(1,hid)

end function funcNameSuffix(At2)

pure function funcNameSuffix(hSlice1) (this, vid) result(elemRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: vid
  FieldDataElemType :: elemRef(getHListSize(this%data))

  elemRef(:) = this%data%v_(vId,:)

end function funcNameSuffix(hSlice1)

pure function funcNameSuffix(hSlice2) (this, vids) result(elemsRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: vids(:)
  FieldDataElemType :: elemsRef(size(vids), getHListSize(this%data))

  elemsRef(:,:) = this%data%v_(vids(:),:)

end function funcNameSuffix(hSlice2)

pure function funcNameSuffix(vSlice1) (this, hid) result(elemRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: hid
  FieldDataElemType :: elemRef(this%vlayerNum)

  elemRef(:) = this%data%v_(:, hid)

end function funcNameSuffix(vSlice1)

pure function funcNameSuffix(vSlice2) (this, hids) result(elemsRef)
  type(FieldTypeName), intent(in) :: this
  integer, intent(in) :: hids(:)
  FieldDataElemType :: elemsRef(this%vlayerNum, size(hids))

  elemsRef(:,:) = this%data%v_(:, hids(:))

end function funcNameSuffix(vSlice2)

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

! Define mul operators
#define op *
#define opname mulOptr

#ifdef DATATYPE_VECTOR
#define DISABLE_OPTRGENTYPE1
#endif
#define DISABLE_OPTRGENTYPE3

#include "geometricfield_binaryOptr.template.f90"

#ifdef DATATYPE_VECTOR
#undef DISABLE_OPTRGENTYPE1
#endif
#undef DISABLE_OPTRGENTYPE3
#undef op
#undef opname

!*** Define div operators
!*
#define op /
#define opname divOptr

#ifdef DATATYPE_VECTOR
#define DISABLE_OPTRGENTYPE1
#define DISABLE_OPTRGENTYPE2
#endif


#include "geometricfield_binaryOptr.template.f90"

#ifdef DATATYPE_VECTOR
#undef DISABLE_OPTRGENTYPE1
#undef DISABLE_OPTRGENTYPE2
#endif

#undef op
#undef opname


!
subroutine checkDataSizeEquality(data1, data2, info)
  type(FieldDataType), intent(in) :: data1, data2
  character(*), intent(in) :: info

  if( getHListSize(data1) /= getHListSize(data2) &
       & .or. getVListSize(data1) /= getVListSize(data2) ) then
    write(*,*) info, ": The size of list data in binary operation is not equal to other."
    stop
  end if
  
end subroutine checkDataSizeEquality


end module moduleName
