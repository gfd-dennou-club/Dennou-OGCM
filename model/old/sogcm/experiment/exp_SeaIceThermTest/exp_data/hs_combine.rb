#!/usr/bin/ruby

require 'numru/ggraph'
include NumRu

ncfile = './ctrl.nc'

gp_hs = GPhys::IO.open(ncfile, 'hs')
gp_hi  = GPhys::IO.open(ncfile, 'hi')              

gp_hs = gp_hs + gp_hi

newncfile = NetCDF.create("ctrl_hs_comb.nc")
GPhys::IO.write(newncfile, gp_hs)
newncfile.close
