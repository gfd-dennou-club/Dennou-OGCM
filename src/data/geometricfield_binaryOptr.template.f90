#define addId(name, id) name ## id
#define optrFuncNameSuffix(name,id) _eval2(addId, _eval2(_funcNameSuffix, FieldTypeName, name), id)

#ifndef DISABLE_OPTRGENTYPE1
function optrFuncNameSuffix(opname,1) (field1, field2) result(tmpfield)
  type(FieldTypeName), intent(in) :: field1
  type(FieldTypeName), intent(in) :: field2
  type(FieldTypeName) :: tmpField
                                
  integer :: dataSize, i 

  dataSize = getListSize(field1%data)

  call GeometricField_Init(tmpField, field1%mesh, name="temporay data")
  tmpField%tempDataFlag = .true.

  do i=1,dataSize 
    tmpField%data%v_(i) = field1%data%v_(i)  op  field2%data%v_(i) 
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
                                
  integer :: dataSize, i 

  dataSize = getListSize(field2%data)

  call GeometricField_Init(tmpField, field2%mesh, name="temporay data")
  tmpField%tempDataFlag = .true.

  do i=1,dataSize 
    tmpField%data%v_(i) = val  op  field2%data%v_(i) 
  end do

  if( field2%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field2)
  
end function optrFuncNameSuffix(opname,2)
#endif

#ifndef DISABLE_OPTRGENTYPE3
function optrFuncNameSuffix(opname,3) (field1, val) result(tmpfield)
  type(FieldTypeName), intent(in) :: field1
  real(DP), intent(in) :: val
  type(FieldTypeName) :: tmpField
                                
  integer :: dataSize, i 

  dataSize = getListSize(field1%data)

  call GeometricField_Init(tmpField, field1%mesh, name="temporay data")
  tmpField%tempDataFlag = .true.

  do i=1,dataSize 
    tmpField%data%v_(i) = field1%data%v_(i) op val 
  end do

  if( field1%tempDataFlag ) call funcNameSuffix(releaseDataRef) (field1)
  
end function optrFuncNameSuffix(opname,3)
#endif
