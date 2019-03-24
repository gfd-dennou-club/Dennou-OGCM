#!/usr/bin ruby

require 'numru/ggraph'
include NumRu

DIRS=["./exp1_result/", "./exp2_result/", "./exp3_result/"] 
PERIODS= [24.0, 72.0, 36.0]
FILENAME="Run1_k20_Psi.nc"
NEWFILENAME="Run1_k20_Psi_2.nc"
VARNAME="Psi"


def normalize_taxis(fname, newfname, varname, period)
  gp = GPhys::IO.open(fname, varname)
  t = GPhys::IO.open(fname, 't')
  t2 = t.copy
  t2 = t / period
  t2.long_name = "time normalized by wave period"
  t2.units = "(1)"

  newgrid = gp.grid_copy.change_axis(3,Axis.new.set_pos(t2.data))

  ofile = NetCDF.create(newfname)

  gp2 = GPhys.new(newgrid, gp.data.rename(varname + "_2"))
  GPhys::IO.write(ofile, gp2)
  ofile.close
end


DIRS.each_with_index{|dir,i|

  fname = ""; newfname = "";
  fname << dir << FILENAME; 
  newfname << dir << NEWFILENAME;
  p "#{fname},  #{newfname}"
  normalize_taxis(fname, newfname, VARNAME, PERIODS[i])

  p "ok!"
}
