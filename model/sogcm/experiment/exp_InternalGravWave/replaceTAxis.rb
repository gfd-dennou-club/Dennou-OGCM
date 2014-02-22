#!/usr/bin ruby

require 'numru/ggraph'
include NumRu

DIRS=["./exp1_result/", "./exp2_result/", "./exp3_result/", "./exp4_result/", "./exp5_result/"] 
PERIODS= [20.6, 11.9, 11.9, 41.2, 205.8]
VARNAME="Div"
FILENAME="Run1_k20_#{VARNAME}.nc"
NEWFILENAME="Run1_k20_#{VARNAME}_2.nc"



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
