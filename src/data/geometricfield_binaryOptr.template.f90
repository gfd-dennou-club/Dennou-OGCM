#define addId(name, id) name ## id
#define optrFuncNameSuffix(name,id) _eval2(addId, _eval2(_funcNameSuffix, FieldTypeName, name), id)

#ifndef DISABLE_OPTRGENTYPE1
function optrFuncNameSuffix(opname,1) (field1, field2) result(tmpfield)
  type(FieldTypeName), intent(in) :: field1
  type(FieldTypeName), intent(in) :: field2
  type(FieldTypeName) :: tmpField
                                
  integer :: i, j 

  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.

  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j)  op  field2%data%v_(i,j)
     end do
  end do

  if( field1%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field1)
  if( field2%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field2)
  
end function optrFuncNameSuffix(opname,1)
#endif

#ifndef DISABLE_OPTRGENTYPE2
function optrFuncNameSuffix(opname,2) (val, field2) result(tmpfield)
  real(DP), intent(in) :: val
  type(FieldTypeName), intent(in) :: field2
  type(FieldTypeName) :: tmpField
                                
  integer :: i, j

  call GeometricField_Init(tmpField, field2%mesh, name="temporay data", vlayerNum=field2%vlayerNum)
  tmpField%tempDataFlag = .true.

  !$omp parallel do private(i)
  do j=1,getHListSize(field2%data)
     do i=1, getVListSize(field2%data)
        tmpField%data%v_(i,j) = val  op  field2%data%v_(i,j)
     end do
  end do

  if( field2%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field2)
  
end function optrFuncNameSuffix(opname,2)
#endif

#ifndef DISABLE_OPTRGENTYPE3
function optrFuncNameSuffix(opname,3) (field1, val) result(tmpfield)
  type(FieldTypeName), intent(in) :: field1
  real(DP), intent(in) :: val
  type(FieldTypeName) :: tmpField
                                
  integer :: i, j


  call GeometricField_Init(tmpField, field1%mesh, name="temporay data", vlayerNum=field1%vlayerNum)
  tmpField%tempDataFlag = .true.

  !$omp parallel do private(i)
  do j=1,getHListSize(field1%data)
     do i=1, getVListSize(field1%data)
        tmpField%data%v_(i,j) = field1%data%v_(i,j) op val
     end do
  end do

  if( field1%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field1)
  
end function optrFuncNameSuffix(opname,3)
#endif
