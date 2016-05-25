module surfaceScalarField_mod
  !
  use dc_types
  use VectorSpace_mod
  use List_mod
  use PolyMesh_mod
  !
  implicit none
  private
  !
  type, public :: surfaceScalarField
    type(List_d), pointer :: data => null()
    type(PolyMesh), pointer :: mesh => null()
    character(TOKEN) :: name
    character(STRING) :: long_name
    character(TOKEN) :: units
    logical :: tempDataFlag = .false.
    integer :: vlayerNum
    integer :: vHaloSize = 0
  end type surfaceScalarField
  !
  !
  interface GeometricField_Init
    module procedure surfaceScalarField_Init
  end interface GeometricField_Init
  !
  interface GeometricField_Final
    module procedure surfaceScalarField_Final
  end interface GeometricField_Final
  interface SetFieldAtitude
    module procedure surfaceScalarField_SetFieldAtitude
  end interface SetFieldAtitude
  interface Release
     module procedure surfaceScalarField_releaseDataRef
  end interface Release
  interface DeepCopy
     module procedure surfaceScalarField_DeepCopy
  end interface DeepCopy
  !
  interface assignment(=)
    module procedure surfaceScalarField_assignFieldObj
    module procedure surfaceScalarField_assignFDElemType
    module procedure surfaceScalarField_assignFDElemTypeArray
  end interface assignment(=)
  interface operator(+)
    module procedure surfaceScalarField_addOptr1
    module procedure surfaceScalarField_addOptr2
    module procedure surfaceScalarField_addOptr3
  end interface operator(+)
  interface operator(-)
    module procedure surfaceScalarField_subOptr1
    module procedure surfaceScalarField_subOptr2
    module procedure surfaceScalarField_subOptr3
  end interface operator(-)
  interface operator(*)
     module procedure surfaceScalarField_mulOptr1
     module procedure surfaceScalarField_mulOptr2
  end interface operator(*)
  interface operator(/)
     module procedure surfaceScalarField_divOptr1
     module procedure surfaceScalarField_divOptr2
     module procedure surfaceScalarField_divOptr3
  end interface operator(/)
  interface At
     module procedure surfaceScalarField_At1
     module procedure surfaceScalarField_At2
  end interface At
  interface hSlice
     module procedure surfaceScalarField_hSlice1
     module procedure surfaceScalarField_hSlice2
  end interface hSlice
  interface vSlice
     module procedure surfaceScalarField_vSlice1
     module procedure surfaceScalarField_vSlice2
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
subroutine surfaceScalarField_Init (field, mesh, name, long_name, units, vlayerNum, vHaloSize)
  type(surfaceScalarField), intent(inout) :: field
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
  select case( DATADEFLOC_surface )
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
end subroutine surfaceScalarField_Init
subroutine surfaceScalarField_SetFieldAtitude (field, name, long_name, units)
  type(surfaceScalarField), intent(inout) :: field
  character(*), optional :: name
  character(*), optional :: long_name
  character(*), optional :: units
  if(present(name)) field%name = name
  if(present(long_name)) field%long_name = long_name
  if(present(units)) field%units = units
end subroutine surfaceScalarField_SetFieldAtitude
!
subroutine surfaceScalarField_Final (field)
  type(surfaceScalarField), intent(inout) :: field
  !write(*,*) "Call releasedataRef"
  call surfaceScalarField_releaseDataRef (field)
end subroutine surfaceScalarField_Final
!
subroutine surfaceScalarField_DeepCopy (self, field)
  type(surfaceScalarField), intent(inout) :: self
  type(surfaceScalarField), intent(in) :: field
  self%data%v_(:,:) = field%data%v_(:,:)
end subroutine surfaceScalarField_DeepCopy
!
subroutine surfaceScalarField_releaseDataRef (field)
  type(surfaceScalarField) :: field
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
end subroutine surfaceScalarField_releaseDataRef
! operation
!
subroutine surfaceScalarField_assignFieldObj (this, other)
  type(surfaceScalarField), intent(inout) :: this
  type(surfaceScalarField), intent(in) :: other
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  if( .not. associated(other%data) ) then
    write(*,*) "Error: try to assignment to uninitialized list data!"
    stop
  end if
  !
  !write(*,*) "assign:", this%name, "=", other%name
  call surfaceScalarField_releaseDataRef (this)
  this%data => other%data
  call incRef(this%data)
  this%mesh => other%mesh
  this%vlayerNum = other%vlayerNum
  this%vHaloSize = other%vHaloSize
  if( other%tempDataFlag ) then
    call surfaceScalarField_releaseDataRef (other)
  end if
end subroutine surfaceScalarField_assignFieldObj
subroutine surfaceScalarField_assignFDElemType (this, elemData)
  type(surfaceScalarField), intent(inout) :: this
  real(DP), intent(in) :: elemData
! call checkDataSizeEquality(field1%data, field2%data, "assign operator")
  this%data%v_ = elemData
end subroutine surfaceScalarField_assignFDElemType
subroutine surfaceScalarField_assignFDElemTypeArray (this, elemArrayData)
  type(surfaceScalarField), intent(inout) :: this
  real(DP), intent(in) :: elemArrayData(:,:)
!!$ if( getListSize(this%data) /= size(elemArrayData) ) then
!!$ write(*,*) "Error in assignment:"
!!$ end if
  this%data%v_(:,:) = elemArrayData(:,:)
end subroutine surfaceScalarField_assignFDElemTypeArray
pure function surfaceScalarField_At1 (this, k, hid) result(elemRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: hid, k
  real(DP) :: elemRef
  elemRef = this%data%v_(k,hid)
end function surfaceScalarField_At1
pure function surfaceScalarField_At2 (this, hid) result(elemRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: hid
  real(DP) :: elemRef
  elemRef = this%data%v_(1,hid)
end function surfaceScalarField_At2
pure function surfaceScalarField_hSlice1 (this, vid) result(elemRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: vid
  real(DP) :: elemRef(getHListSize(this%data))
  elemRef(:) = this%data%v_(vId,:)
end function surfaceScalarField_hSlice1
pure function surfaceScalarField_hSlice2 (this, vids) result(elemsRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: vids(:)
  real(DP) :: elemsRef(size(vids), getHListSize(this%data))
  elemsRef(:,:) = this%data%v_(vids(:),:)
end function surfaceScalarField_hSlice2
pure function surfaceScalarField_vSlice1 (this, hid) result(elemRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: hid
  real(DP) :: elemRef(this%vlayerNum)
  elemRef(:) = this%data%v_(:, hid)
end function surfaceScalarField_vSlice1
pure function surfaceScalarField_vSlice2 (this, hids) result(elemsRef)
  type(surfaceScalarField), intent(in) :: this
  integer, intent(in) :: hids(:)
  real(DP) :: elemsRef(this%vlayerNum, size(hids))
  elemsRef(:,:) = this%data%v_(:, hids(:))
end function surfaceScalarField_vSlice2
!*
!*
! Define add operators
function surfaceScalarField_addOptr1 (field1, field2) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_addOptr1
function surfaceScalarField_addOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val + field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_addOptr2
function surfaceScalarField_addOptr3 (field1, val) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) + val
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
end function surfaceScalarField_addOptr3
! Define sub operators
function surfaceScalarField_subOptr1 (field1, field2) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_subOptr1
function surfaceScalarField_subOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val - field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_subOptr2
function surfaceScalarField_subOptr3 (field1, val) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) - val
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
end function surfaceScalarField_subOptr3
! Define mul operators
function surfaceScalarField_mulOptr1 (field1, field2) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) * field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_mulOptr1
function surfaceScalarField_mulOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val * field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_mulOptr2
!*** Define div operators
!*
function surfaceScalarField_divOptr1 (field1, field2) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) / field2%data%v_(i,j)
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_divOptr1
function surfaceScalarField_divOptr2 (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(surfaceScalarField), intent(in) :: field2
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val / field2%data%v_(i,j)
     end do
  end do
  if( field2%tempDataFlag ) call surfaceScalarField_releaseDataRef (field2)
end function surfaceScalarField_divOptr2
function surfaceScalarField_divOptr3 (field1, val) result(tmpfield)
  type(surfaceScalarField), intent(in) :: field1
  real(DP), intent(in) :: val
  type(surfaceScalarField) :: tmpField
  integer :: i, j
  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.
  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) / val
     end do
  end do
  if( field1%tempDataFlag ) call surfaceScalarField_releaseDataRef (field1)
end function surfaceScalarField_divOptr3
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
end module surfaceScalarField_mod
