module pointVectorField_mod
  !
  use dc_types
  use VectorSpace_mod
  use List_mod
  use PolyMesh_mod
  !
  implicit none
  private
  !
  type, public :: pointVectorField
    type(List_vec3d), pointer :: data => null()
    type(PolyMesh), pointer :: mesh => null()
    character(TOKEN) :: name
    character(STRING) :: long_name
    character(TOKEN) :: units
    logical :: tempDataFlag = .false.
    integer :: vlayerNum
    integer :: vHaloSize = 0
  end type pointVectorField
  !
  !
  interface GeometricField_Init
    module procedure pointVectorField_Init
  end interface GeometricField_Init
  !
  interface GeometricField_Final
    module procedure pointVectorField_Final
  end interface GeometricField_Final
  interface SetFieldAtitude
    module procedure pointVectorField_SetFieldAtitude
  end interface SetFieldAtitude
  interface Release
     module procedure pointVectorField_releaseDataRef
  end interface Release
  interface DeepCopy
     module procedure pointVectorField_DeepCopy
  end interface DeepCopy
  !
  interface assignment(=)
    module procedure pointVectorField_assignFieldObj
    module procedure pointVectorField_assignFDElemType
    module procedure pointVectorField_assignFDElemTypeArray
  end interface assignment(=)
  interface operator(+)
    module procedure pointVectorField_addOptr1
    module procedure pointVectorField_addOptr2
    module procedure pointVectorField_addOptr3
  end interface operator(+)
  interface operator(-)
    module procedure pointVectorField_subOptr1
    module procedure pointVectorField_subOptr2
    module procedure pointVectorField_subOptr3
  end interface operator(-)
  interface operator(*)
     module procedure pointVectorField_mulOptr2
  end interface operator(*)
  interface operator(/)
     module procedure pointVectorField_divOptr3
  end interface operator(/)
  interface At
     module procedure pointVectorField_At1
     module procedure pointVectorField_At2
  end interface At
  interface hSlice
     module procedure pointVectorField_hSlice1
     module procedure pointVectorField_hSlice2
  end interface hSlice
  interface vSlice
     module procedure pointVectorField_vSlice1
     module procedure pointVectorField_vSlice2
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
subroutine pointVectorField_Init (field, mesh, name, long_name, units, vlayerNum, vHaloSize)
  type(pointVectorField), intent(inout) :: field
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
end subroutine pointVectorField_Init
subroutine pointVectorField_SetFieldAtitude (field, name, long_name, units)
  type(pointVectorField), intent(inout) :: field
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units
  if(present(name)) field%name = name
  if(present(long_name)) field%long_name = long_name
  if(present(units)) field%units = units
end subroutine pointVectorField_SetFieldAtitude
!
subroutine pointVectorField_Final (field)
  type(pointVectorField), intent(inout) :: field
  !write(*,*) "Call releasedataRef"
  call pointVectorField_releaseDataRef (field)
end subroutine pointVectorField_Final
!
subroutine pointVectorField_DeepCopy (self, field)
  type(pointVectorField), intent(inout) :: self
  type(pointVectorField), intent(in) :: field
  self%data%v_(:,:) = field%data%v_(:,:)
end subroutine pointVectorField_DeepCopy
!
subroutine pointVectorField_releaseDataRef (field)
  type(pointVectorField) :: field
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
end subroutine pointVectorField_releaseDataRef
! operation
!
subroutine pointVectorField_assignFieldObj (this, other)
  type(pointVectorField), intent(inout) :: this
  type(pointVectorField), intent(in) :: other
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  if( .not. associated(other%data) ) then
    write(*,*) "Error: try to assignment to uninitialized list data!"
    stop
  end if
  !
  !write(*,*) "assign:", this%name, "=", other%name
  call pointVectorField_releaseDataRef (this)
  this%data => other%data
  call incRef(this%data)
  this%mesh => other%mesh
  this%vlayerNum = other%vlayerNum
  this%vHaloSize = other%vHaloSize
  if( other%tempDataFlag ) then
    call pointVectorField_releaseDataRef (other)
  end if
end subroutine pointVectorField_assignFieldObj
subroutine pointVectorField_assignFDElemType (this, elemData)
  type(pointVectorField), intent(inout) :: this
  type(vector3d), intent(in) :: elemData
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  this%data%v_ = elemData
end subroutine pointVectorField_assignFDElemType
subroutine pointVectorField_assignFDElemTypeArray (this, elemArrayData)
  type(pointVectorField), intent(inout) :: this
  type(vector3d), intent(in) :: elemArrayData(:,:)
!!$ if( getListSize(this%data) /= size(elemArrayData) ) then
!!$ write(*,*) "Error in assignment:"
!!$ end if
  this%data%v_(:,:) = elemArrayData(:,:)
end subroutine pointVectorField_assignFDElemTypeArray
pure function pointVectorField_At1 (this, k, hid) result(elemRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: hid, k
  type(vector3d) :: elemRef
  elemRef = this%data%v_(k,hid)
end function pointVectorField_At1
pure function pointVectorField_At2 (this, hid) result(elemRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: hid
  type(vector3d) :: elemRef
  elemRef = this%data%v_(1,hid)
end function pointVectorField_At2
pure function pointVectorField_hSlice1 (this, vid) result(elemRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: vid
  type(vector3d) :: elemRef(getHListSize(this%data))
  elemRef(:) = this%data%v_(vId,:)
end function pointVectorField_hSlice1
pure function pointVectorField_hSlice2 (this, vids) result(elemsRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: vids(:)
  type(vector3d) :: elemsRef(size(vids), getHListSize(this%data))
  elemsRef(:,:) = this%data%v_(vids(:),:)
end function pointVectorField_hSlice2
pure function pointVectorField_vSlice1 (this, hid) result(elemRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: hid
  type(vector3d) :: elemRef(this%vlayerNum)
  elemRef(:) = this%data%v_(:, hid)
end function pointVectorField_vSlice1
pure function pointVectorField_vSlice2 (this, hids) result(elemsRef)
  type(pointVectorField), intent(in) :: this
  integer, intent(in) :: hids(:)
  type(vector3d) :: elemsRef(this%vlayerNum, size(hids))
  elemsRef(:,:) = this%data%v_(:, hids(:))
end function pointVectorField_vSlice2
!*
!*
! Define add operators
function pointVectorField_addOptr1 (field1, field2) result(tmpfield)
  type(pointVectorField), intent(in) :: field1
  type(pointVectorField), intent(in) :: field2
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointVectorField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointVectorField_releaseDataRef (field2)
end function pointVectorField_addOptr1
function pointVectorField_addOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointVectorField), intent(in) :: field2
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val + field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointVectorField_releaseDataRef (field2)
end function pointVectorField_addOptr2
function pointVectorField_addOptr3 (field1, val) result(tmpfield)
  type(pointVectorField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + val
     end do
  end do
  if( field1%tempDataFlag ) call pointVectorField_releaseDataRef (field1)
end function pointVectorField_addOptr3
! Define sub operators
function pointVectorField_subOptr1 (field1, field2) result(tmpfield)
  type(pointVectorField), intent(in) :: field1
  type(pointVectorField), intent(in) :: field2
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call pointVectorField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call pointVectorField_releaseDataRef (field2)
end function pointVectorField_subOptr1
function pointVectorField_subOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointVectorField), intent(in) :: field2
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val - field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointVectorField_releaseDataRef (field2)
end function pointVectorField_subOptr2
function pointVectorField_subOptr3 (field1, val) result(tmpfield)
  type(pointVectorField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - val
     end do
  end do
  if( field1%tempDataFlag ) call pointVectorField_releaseDataRef (field1)
end function pointVectorField_subOptr3
! Define mul operators
function pointVectorField_mulOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(pointVectorField), intent(in) :: field2
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val * field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call pointVectorField_releaseDataRef (field2)
end function pointVectorField_mulOptr2
!*** Define div operators
!*
function pointVectorField_divOptr3 (field1, val) result(tmpfield)
  type(pointVectorField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(pointVectorField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) / val
     end do
  end do
  if( field1%tempDataFlag ) call pointVectorField_releaseDataRef (field1)
end function pointVectorField_divOptr3
!
subroutine checkDataSizeEquality(data1, data2, info)
  type(List_vec3d), intent(in) :: data1, data2
  character(*), intent(in) :: info
  if( getHListSize(data1) /= getHListSize(data2) &
       & .or. getVListSize(data1) /= getVListSize(data2) ) then
    write(*,*) info, ": The size of list data in binary operation is not equal to other."
    stop
  end if
end subroutine checkDataSizeEquality
end module pointVectorField_mod
