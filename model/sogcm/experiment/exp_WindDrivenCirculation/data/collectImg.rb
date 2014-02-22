#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

# * exp0: wind stress provided with analystic function.
# * exp1: wind stress provided with the result of aquaplanet numerical simulation in Marshall et al(2007). 
#    - exp1-1: The number of vertical grid points is 20(k=20).
#    - exp1-2: The number of vertical grid points is 40(k=40).
#    - exp1-3: The cofficient of horizontal eddy diffusion is twice larger than exp1-1. 

require "./../../collectImgLib"
include CollectImgLib

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"


VarNoAnimFig = Struct.new("VarNoAnimFig", :ncFileSuffix, :name, :cutpos, :intrv, :range, :prefix, :suffix, :gpOpts)
VarAnimFig = Struct.new("VarAnimFig", :name, :cutpos, :intrv, :range, :animTimeInfo, :prefix, :suffix, :gpOpts)


LAST_DAY=11300
xy_Div_SeaSurf_init = VarNoAnimFig.new("", "Div", "t=0,sig=0",  "1e-13", "-8e-13:8e-13", "xy", "SeaSurf_init","") 
xz_Div_lat0_init = VarNoAnimFig.new("", "Div", "t=0,lat=0",  "1e-13", "-8e-13:8e-13", "xz", "lat0_init","") 
xz_Div_lat45_init = VarNoAnimFig.new("", "Div", "t=0,lat=45",  "1e-13", "-8e-13:8e-13", "xz", "lat45_init","") 

#tz_SigDot_lat0lon180 = VarNoAnimFig.new("SigDot", "lat=0,lon=180",  "2.5e-14", "-3e-13:3e-13", "tz", "lat0lon180", "--exch")
xy_U_SeaSurf = VarNoAnimFig.new("", "U", "lon=0,sig=0", "0.1", "-1.8:1.8", "xy", "SeaSurf", "--clrmap 10 --title 'U' ")
yz_U_mplane = VarNoAnimFig.new("", "U", "lon=0,t=#{LAST_DAY}", "0.1", "-1.8:1.8", "yz", "mplane", "--clrmap 10 --title 'U' ")
yz_MassStrFunc_mplane = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY}", "0.05", "-0.35:0.35", "yz", "mplane", "--title 'MassStreamFunc'")
yz_PTempEdd_mplane = VarNoAnimFig.new("", "PTempEdd", "lon=0,t=#{LAST_DAY}", "0.5", "-5:5", "yz", "mplane", "--title 'PTempEdd'")
t_KEAvg = VarNoAnimFig.new("EnergyBudget", "KEAvg", "", "", "0:0.06", "t", "", "--title 'K.E.(global mean)'" )


exp_noAnimFigs =  [ t_KEAvg ] #xy_U_SeaSurf, yz_U_mplane, yz_MassStrFunc_mplane, yz_PTempEdd_mplane ]
noAnimFigVars = [
# [xy_Div_SeaSurf_init, xz_Div_lat0_init, tz_SigDot_lat0lon180, tz_Div_lat0lon180], 
    exp_noAnimFigs, exp_noAnimFigs, exp_noAnimFigs, exp_noAnimFigs
]

DATASTORAGE_DIR="/home/ykawai/exp_backup/current/"
exps = [
        Exp.new("exp_Ah800L20", "#{DATASTORAGE_DIR}/exp_Ah800L20/", "Run1_k20_"), 
        Exp.new("exp_Ah800L60", "#{DATASTORAGE_DIR}/exp_Ah800L60/", "Run1_k60_"), 
        Exp.new("exp_Ah400L20", "#{DATASTORAGE_DIR}/exp_Ah400L20/", "Run1_k20_"),
        Exp.new("exp_Ah1600L20", "#{DATASTORAGE_DIR}/exp_Ah1600L20/", "Run1_k20_")
]

animTimeInfo = AnimTimeInfo.new(0, 10000, 100.0, "day")
animFigFlag = false
animFigVars = [
# [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=45", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ], 
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ]
]

#
exp_common = 

extra_exps = [ Exp.new("exp_comm", "./common/"), 
               Exp.new("exp_LComp", "./LCompare/"),
               Exp.new("exp_AhComp", "./AhCompare/")]
extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png")
}

####################################################

exps.each_with_index{|exp, i|
  puts "= experiment: #{exp.name}"

  noAnimFigVars[i].each{|figvar|
      exp.add_NoAnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.prefix, figvar.suffix, figvar.gpOpts, "", figvar.ncFileSuffix) 
  }

if animFigFlag then
  animFigVars[i].each{|figvar|
    exp.add_AnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.animTimeInfo, figvar.prefix, figvar.suffix, figvar.gpOpts)
p figvar.animTimeInfo 
  }
end if

  exp.create_allFig
#  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0 and animFigFlag
}

