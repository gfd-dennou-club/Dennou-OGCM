#!/usr/bin/ruby

require 'numru/ggraph'
include NumRu

NCFILE = "PTempTempDiff.nc"

ptemptempdiff_JM95 = GPhys::IO.open(NCFILE, 'EOS_JM95')
ptemptempdiff_Linear = GPhys::IO.open(NCFILE, 'EOS_LINEAR')
ptemptempdiff_NonLinear = GPhys::IO.open(NCFILE, 'EOS_SIMPLENONLINEAR')

temp = ptemptempdiff_JM95.axis("temp").to_gphys + 273.15
ptemptempdiffError_Linear = (ptemptempdiff_Linear - ptemptempdiff_JM95)
ptemptempdiffError_NonLinear = (ptemptempdiff_NonLinear - ptemptempdiff_JM95)

ptemptempdiffError_Linear.name = "ptemptempdiffError_Linear"
ptemptempdiffError_Linear.long_name = "error of EOS_Linear to EOS_JM95"
ptemptempdiffError_NonLinear.name = "ptemptempdiffError_SimpleNonLinear"
ptemptempdiffError_NonLinear.long_name = "error of EOS_SimpleNonLinear to EOS_JM95"

p ptemptempdiffError_Linear
p ptemptempdiffError_NonLinear

newncfile = NetCDF.create(File.dirname(NCFILE) + "/PTempTempDiffError.nc")
GPhys::IO.write(newncfile, ptemptempdiffError_Linear)
GPhys::IO.write(newncfile, ptemptempdiffError_NonLinear)
newncfile.close

rmsError = ( (ptemptempdiffError_Linear**2).mean(2).mean(1).mean(0) ).sqrt
p "Linear= #{rmsError.to_f}"
rmsError = ( (ptemptempdiffError_NonLinear**2).mean(2).mean(1).mean(0) ).sqrt
p "SimpleNonLinear= #{rmsError.to_f}"
