#!/usr/bin/ruby

require 'numru/ggraph'
include NumRu

NCFILE = "dens.nc"

dens_JM95 = GPhys::IO.open(NCFILE, 'EOS_JM95')
dens_Linear = GPhys::IO.open(NCFILE, 'EOS_LINEAR')
dens_NonLinear = GPhys::IO.open(NCFILE, 'EOS_SIMPLENONLINEAR')

RefDens = UNumeric[ 1.027e03, 'kg.m-3' ]
densError_Linear = (dens_Linear - dens_JM95)/(dens_JM95 + RefDens)
densError_NonLinear = (dens_NonLinear - dens_JM95)/(dens_JM95 + RefDens)

densError_Linear.name = "densError_Linear"
densError_Linear.long_name = "relative error of EOS_Linear to EOS_JM95"
densError_NonLinear.name = "densError_SimpleNonLinear"
densError_NonLinear.long_name = "relative error of EOS_SimpleNonLinear to EOS_JM95"

p densError_Linear
p densError_NonLinear

newncfile = NetCDF.create(File.dirname(NCFILE) + "/densError.nc")
GPhys::IO.write(newncfile, densError_Linear)
GPhys::IO.write(newncfile, densError_NonLinear)
newncfile.close

rmsError = ( (densError_Linear**2).mean(2).mean(1).mean(0) ).sqrt
p "Linear= #{rmsError.to_f}"
rmsError = ( (densError_NonLinear**2).mean(2).mean(1).mean(0) ).sqrt
p "SimpleNonLinear= #{rmsError.to_f}"
