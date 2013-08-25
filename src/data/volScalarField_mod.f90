module volScalarField_mod
  !
  use dc_types
  use VectorSpace_mod
  use List_mod
  use PolyMesh_mod
  !
  implicit none
  private
  !
  type, public :: volScalarField
    type(List_d), pointer :: data => null()
    type(PolyMesh), pointer :: mesh => null()
    character(TOKEN) :: name
    character(STRING) :: long_name
    character(TOKEN) :: units
    logical :: tempDataFlag = .false.
  end type volScalarField
  !
  !
  interface GeometricField_Init
    module procedure volScalarField_Init
  end interface GeometricField_Init
  !
  interface GeometricField_Final
    module procedure volScalarField_Final
  end interface GeometricField_Final
  interface SetFieldAtitude
    module procedure volScalarField_SetFieldAtitude
  end interface SetFieldAtitude
  !
  interface assignment(=)
    module procedure volScalarField_assignOptr
  end interface assignment(=)
  interface operator(+)
    module procedure volScalarField_addOptr1
! module procedure volScalarField_addOptr2
! module procedure volScalarField_addOptr3
  end interface operator(+)
  interface operator(-)
    module procedure volScalarField_subOptr1
! module procedure volScalarField_subOptr2
! module procedure volScalarField_subOptr3
  end interface operator(-)
  !
  public :: GeometricField_Init, GeometricField_Final
  public :: SetFieldAtitude
  public :: operator(+), assignment(=)
  character(TOKEN), parameter :: DATADEFLOC_SURFACE = "surface"
  character(TOKEN), parameter :: DATADEFLOC_POINT = "point"
  character(TOKEN), parameter :: DATADEFLOC_VOL = "vol"
contains
!
subroutine volScalarField_Init (field, mesh, name, long_name, units)
  type(volScalarField), intent(inout) :: field
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
  call GeometrcField_Final(field)
  !
  allocate(field%data)
  call List_Init(field%data, dataSize)
  !
  field%mesh => mesh
  call setFieldAtitude(field, "field", " ", "1")
  if(present(name)) call setFieldAtitude(field, name=name)
  if(present(long_name)) call setFieldAtitude(field, long_name=long_name)
  if(present(units)) call setFieldAtitude(field, units=units)
  !
  select case( DATADEFLOC_vol )
    case (DATADEFLOC_SURFACE)
      dataSize = getFaceListSize(field%mesh)
    case (DATADEFLOC_VOL)
      dataSize = getCellListSize(field%mesh)
    case (DATADEFLOC_POINT)
      dataSize = getPointListSize(field%mesh)
  end select
end subroutine volScalarField_Init
subroutine volScalarField_SetFieldAtitude (field, name, long_name, units)
  type(volScalarField), intent(inout) :: field
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units
  if(present(name)) field%name = name
  if(present(long_name)) field%long_name = long_name
  if(present(units)) field%units = units
end subroutine volScalarField_SetFieldAtitude
!
subroutine volScalarField_Final (field)
  type(volScalarField), intent(inout) :: field
  call volScalarField_releaseDataRef (field)
end subroutine volScalarField_Final
!
subroutine volScalarField_releaseDataRef (field)
  type(volScalarField), intent(inout) :: field
  logical :: refIs0
  if( .not. associated(field%data) ) return
  call decRef(field%data, refIs0)
  if( refIs0 ) then
    write(*,*) "Safe to release a resource of list data."
    call List_Final(field%data)
    deallocate(field%data)
  end if
  field%data => null()
end subroutine volScalarField_releaseDataRef
! operation
!
subroutine volScalarField_assignOptr (this, other)
  type(volScalarField), intent(inout) :: this
  type(volScalarField), intent(in) :: other
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  if( .not. associated(other%data) ) then
    write(*,*) "Error: try to assignment to uninitialized list data!"
    stop
  end if
  !
  call volScalarField_releaseDataRef (this)
  this%data => other%data
  call incRef(other%data)
end subroutine volScalarField_assignOptr
!*
!*
! Define add operators
function volScalarField_addOptr1 (field1, field2) result(tmpfield)
  type(volScalarField), intent(in) :: field1
  type(volScalarField), intent(in) :: field2
  type(volScalarField) :: tmpField
  integer :: dataSize, i
  dataSize = getListSize(field1%data)
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data")
  tmpField%tempDataFlag = .true.
  do i=1,dataSize
    tmpField%data%v_(i) = field1%data%v_(i) + field2%data%v_(i)
  end do
  if( field1%tempDataFlag ) then
    call List_Final(field1%data)
  end if
! if( field2%tempDataFlag ) call volScalarField_releaseDataRef(field2)
end function volScalarField_addOptr1
! Define sub operators
function volScalarField_subOptr1 (field1, field2) result(tmpfield)
  type(volScalarField), intent(in) :: field1
  type(volScalarField), intent(in) :: field2
  type(volScalarField) :: tmpField
  integer :: dataSize, i
  dataSize = getListSize(field1%data)
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data")
  tmpField%tempDataFlag = .true.
  do i=1,dataSize
    tmpField%data%v_(i) = field1%data%v_(i) - field2%data%v_(i)
  end do
  if( field1%tempDataFlag ) then
    call List_Final(field1%data)
  end if
! if( field2%tempDataFlag ) call volScalarField_releaseDataRef(field2)
end function volScalarField_subOptr1
!
subroutine checkDataSizeEquality(data1, data2, info)
  type(List_d), intent(in) :: data1, data2
  character(*), intent(in) :: info
  if( getListSize(data1) /= getListSize(data2) ) then
    write(*,*) info, ": The size of list data in binary operation is not equal to other."
    stop
  end if
end subroutine checkDataSizeEquality
end module volScalarField_mod
