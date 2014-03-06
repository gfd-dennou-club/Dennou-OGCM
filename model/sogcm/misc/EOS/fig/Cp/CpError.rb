#!/usr/bin/ruby

require 'numru/ggraph'
include NumRu

NCFILE = "Cp.nc"

cp_JM95 = GPhys::IO.open(NCFILE, 'EOS_JM95')
cp_Linear = GPhys::IO.open(NCFILE, 'EOS_LINEAR')
cp_NonLinear = GPhys::IO.open(NCFILE, 'EOS_SIMPLENONLINEAR')

cpError_Linear = (cp_Linear - cp_JM95)/(cp_JM95)
cpError_NonLinear = (cp_NonLinear - cp_JM95)/(cp_JM95)

cpError_Linear.name = "cpError_Linear"
cpError_Linear.long_name = "relative error of EOS_Linear to EOS_JM95"
cpError_NonLinear.name = "cpError_SimpleNonLinear"
cpError_NonLinear.long_name = "relative error of EOS_SimpleNonLinear to EOS_JM95"

p cpError_Linear
p cpError_NonLinear

newncfile = NetCDF.create(File.dirname(NCFILE) + "/cpError.nc")
GPhys::IO.write(newncfile, cpError_Linear)
GPhys::IO.write(newncfile, cpError_NonLinear)
newncfile.close

rmsError = ( (cpError_Linear**2).mean(2).mean(1).mean(0) ).sqrt
p "Linear= #{rmsError.to_f}"
rmsError = ( (cpError_NonLinear**2).mean(2).mean(1).mean(0) ).sqrt
p "SimpleNonLinear= #{rmsError.to_f}"
