module pointScalarField_mod
  !
  use dc_types
  use VectorSpace_mod
  use List_mod
  use PolyMesh_mod
  !
  implicit none
  private
  !
  type, public :: pointScalarField
    type(List_d), pointer :: data => null()
    type(PolyMesh), pointer :: mesh => null()
    character(TOKEN) :: name
    character(STRING) :: long_name
    character(TOKEN) :: units
    logical :: tempDataFlag = .false.
    integer :: vlayerNum
    integer :: vHaloSize = 0
  end type pointScalarField
  !
  !
  interface GeometricField_Init
    module procedure pointScalarField_Init
  end interface GeometricField_Init
  !
  interface GeometricField_Final
    module procedure pointScalarField_Final
  end interface GeometricField_Final
  interface SetFieldAtitude
    module procedure pointScalarField_SetFieldAtitude
  end interface SetFieldAtitude
  interface Release
     module procedure pointScalarField_releaseDataRef
  end interface Release
  interface DeepCopy
     module procedure pointScalarField_DeepCopy
  end interface DeepCopy
  !
  interface assignment(=)
    module procedure pointScalarField_assignFieldObj
    module procedure pointScalarField_assignFDElemType
    module procedure pointScalarField_assignFDElemTypeArray
  end interface assignment(=)
  interface operator(+)
    module procedure pointScalarField_addOptr1
    module procedure pointScalarField_addOptr2
    module procedure pointScalarField_addOptr3
  end interface operator(+)
  interface operator(-)
    module procedure pointScalarField_subOptr1
    module procedure pointScalarField_subOptr2
    module procedure pointScalarField_subOptr3
  end interface operator(-)
  interface operator(*)
     module procedure pointScalarField_mulOptr1
     module procedure pointScalarField_mulOptr2
  end interface operator(*)
  interface operator(/)
     module procedure pointScalarField_divOptr1
     module procedure pointScalarField_divOptr2
     module procedure pointScalarField_divOptr3
  end interface operator(/)
  interface At
     module procedure pointScalarField_At1
     module procedure pointScalarField_At2
  end interface At
  interface hSlice
     module procedure pointScalarField_hSlice1
     module procedure pointScalarField_hSlice2
  end interface hSlice
  interface vSlice
     module procedure pointScalarField_vSlice1
     module procedure pointScalarField_vSlice2
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
subroutine pointScalarField_Init (field, mesh, name, long_name, units, vlayerNum, vHaloSize)
  type(pointScalarField), intent(inout) :: field
  type(PolyMesh), target, intent(in) :: mesh
  character(*), optional, intent(in) :: name
  character(*), optional, intent(in) :: long_name
  character(*), optional, intent(in) :: units
  integer, optional, intent(in) :: vlayerNum
  integer, optional, intent(in) :: vHaloSize
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
  if(present(vHaloSize)) field%vHaloSize = vHaloSize
  !
  select case( DATADEFLOC_point )
    case (DATADEFLOC_SURFACE)
      hdataSize = getFaceListSize(mesh)
    case (DATADEFLOC_VOL)
      hdataSize = getCellListSize(mesh)
    case (DATADEFLOC_POINT)
      hdataSize = getPointListSize(mesh)
  end select
  !
  allocate(field%data)
  call List_Init(field%data, field%vLayerNum, hdataSize, field%vHaloSize)
! write(*,*) field%name, "allocated..", ", hdataSize=", hdataSize, "vlayerNum=", field%vlayerNum
end subroutine pointScalarField_Init
subroutine pointScalarField_SetFieldAtitude (field, name, long_name, units)
  type(pointScalarField), intent(inout) :: field
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units
  if(present(name)) field%name = name
  if(present(long_name)) field%long_name = long_name
  if(present(units)) field%units = units
end subroutine pointScalarField_SetFieldAtitude
!
subroutine pointScalarField_Final (field)
  type(pointScalarField), intent(inout) :: field
  !write(*,*) "Call releasedataRef"
  call pointScalarField_releaseDataRef (field)
end subroutine pointScalarField_Final
!
subroutine pointScalarField_DeepCopy (self, field)
  type(pointScalarField), intent(inout) :: self
  type(pointScalarField), intent(in) :: field
  self%data%v_(:,:) = field%data%v_(:,:)
end subroutine pointScalarField_DeepCopy
!
subroutine pointScalarField_releaseDataRef (field)
  type(pointScalarField) :: field
  logical :: refIs0
  if( .not. associated(field%data) ) return
  call decRef(field%data, refIs0)
  if( refIs0 ) then
    call List_Final(field%data)
    !write(*,*) "Delete data FieldName=", field%name
    !if( field%tempDataFlag ) write(*,*) "This object is a temporary data. ^"
    deallocate(field%data)
  end if
  field%data => null()
end subroutine pointScalarField_releaseDataRef
! operation
!
subroutine pointScalarField_assignFieldObj (this, other)
  type(pointScalarField), intent(inout) :: this
  type(pointScalarField), intent(in) :: other
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  if( .not. associated(other%data) ) then
    write(*,*) "Error: try to assignment to uninitialized list data!"
    stop
  end if
  !
  !write(*,*) "assign:", this%name, "=", other%name
  call pointScalarField_releaseDataRef (this)
  this%data => other%data
  call incRef(this%data)
  this%mesh => other%mesh
  this%vlayerNum = other%vlayerNum
  this%vHaloSize = other%vHaloSize
  if( other%tempDataFlag ) then
    call pointScalarField_releaseDataRef (other)
  end if
end subroutine pointScalarField_assignFieldObj
subroutine pointScalarField_assignFDElemType (this, elemData)
  type(pointScalarField), intent(inout) :: this
  real(DP), intent(in) :: elemData
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  this%data%v_ = elemData
end subroutine pointScalarField_assignFDElemType
subroutine pointScalarField_assignFDElemTypeArray (this, elemArrayData)
  type(pointScalarField), intent(inout) :: this
  real(DP), intent(in) :: elemArrayData(:,:)
!!$ if( getListSize(this%data) /= size(elemArrayData) ) then
!!$ write(*,*) "Error in assignment:"
!!$ end if
  this%data%v_(:,:) = elemArrayData(:,:)
end subroutine pointScalarField_assignFDElemTypeArray
pure function pointScalarField_At1 (this, k, hid) result(elemRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: hid, k
  real(DP) :: elemRef
  elemRef = this%data%v_(k,hid)
end function pointScalarField_At1
pure function pointScalarField_At2 (this, hid) result(elemRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: hid
  real(DP) :: elemRef
  elemRef = this%data%v_(1,hid)
end function pointScalarField_At2
pure function pointScalarField_hSlice1 (this, vid) result(elemRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: vid
  real(DP) :: elemRef(getHListSize(this%data))
  elemRef(:) = this%data%v_(vId,:)
end function pointScalarField_hSlice1
pure function pointScalarField_hSlice2 (this, vids) result(elemsRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: vids(:)
  real(DP) :: elemsRef(size(vids), getHListSize(this%data))
  elemsRef(:,:) = this%data%v_(vids(:),:)
end function pointScalarField_hSlice2
pure function pointScalarField_vSlice1 (this, hid) result(elemRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: hid
  real(DP) :: elemRef(this%vlayerNum)
  elemRef(:) = this%data%v_(:, hid)
end function pointScalarField_vSlice1
pure function pointScalarField_vSlice2 (this, hids) result(elemsRef)
  type(pointScalarField), intent(in) :: this
  integer, intent(in) :: hids(:)
  real(DP) :: elemsRef(this%vlayerNum, size(hids))
  elemsRef(:,:) = this%data%v_(:, hids(:))
end function pointScalarField_vSlice2
!*
!*
! Define add operators
function pointScalarField_addOptr1 (field1, field2) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_addOptr1
function pointScalarField_addOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val + field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_addOptr2
function pointScalarField_addOptr3 (field1, val) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + val
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
end function pointScalarField_addOptr3
! Define sub operators
function pointScalarField_subOptr1 (field1, field2) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_subOptr1
function pointScalarField_subOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val - field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_subOptr2
function pointScalarField_subOptr3 (field1, val) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - val
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
end function pointScalarField_subOptr3
! Define mul operators
function pointScalarField_mulOptr1 (field1, field2) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) * field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_mulOptr1
function pointScalarField_mulOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val * field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_mulOptr2
!*** Define div operators
!*
function pointScalarField_divOptr1 (field1, field2) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) / field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_divOptr1
function pointScalarField_divOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointScalarField), intent(in) :: field2
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val / field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointScalarField_releaseDataRef (field2)
end function pointScalarField_divOptr2
function pointScalarField_divOptr3 (field1, val) result(tmpfield)
  type(pointScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) / val
     end do
  end do
  if( field1%tempDataFlag ) call pointScalarField_releaseDataRef (field1)
end function pointScalarField_divOptr3
!
subroutine checkDataSizeEquality(data1, data2, info)
  type(List_d), intent(in) :: data1, data2
  character(*), intent(in) :: info
  if( getHListSize(data1) /= getHListSize(data2) &
       & .or. getVListSize(data1) /= getVListSize(data2) ) then
    write(*,*) info, ": The size of list data in binary operation is not equal to other."
    stop
  end if
end subroutine checkDataSizeEquality
end module pointScalarField_mod
