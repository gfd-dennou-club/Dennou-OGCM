#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

require "./../collectImgLib"
include CollectImgLib

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"


VarNoAnimFig = Struct.new("VarNoAnimFig", :name, :cutpos, :intrv, :range, :prefix, :suffix, :gpOpts)
VarAnimFig = Struct.new("VarAnimFig", :name, :cutpos, :intrv, :range, :animTimeInfo, :prefix, :suffix, :gpOpts)


animTimeInfos = [
  AnimTimeInfo.new(0, 42, 1.0, "day"), 
  AnimTimeInfo.new(0, 24, 0.5, "day"), 
  AnimTimeInfo.new(0, 24, 0.5, "day"),  
  AnimTimeInfo.new(0, 90, 2.0, "day"),  
  AnimTimeInfo.new(0, 420, 10.0, "day") 
]

xy_Div_SeaSurf_init = VarNoAnimFig.new("Div", "t=0,sig=0",  "1e-13", "-8e-13:8e-13", "xy", "SeaSurf_init","") 
xz_Div_lat0_init = VarNoAnimFig.new("Div", "t=0,lat=0",  "1e-13", "-8e-13:8e-13", "xz", "lat0_init","") 
xz_Div_lat45_init = VarNoAnimFig.new("Div", "t=0,lat=45",  "1e-13", "-8e-13:8e-13", "xz", "lat45_init","") 

tz_SigDot_lat0lon180 = VarNoAnimFig.new("SigDot", "lat=0,lon=180",  "2.5e-14", "-3e-13:3e-13", "tz", "lat0lon180", "--exch")
tz_SigDot_lat45lon180 = VarNoAnimFig.new("SigDot", "lat=45,lon=180",  "2.5e-14", "-3e-13:3e-13", "tz", "lat45lon180", "--exch") 
tz_Div_lat0lon180 = VarNoAnimFig.new("Div", "lat=0,lon=180",  "1e-13", "-8e-13:8e-13", "tz", "lat0lon180", "--exch") 
tz_Div_lat45lon180 = VarNoAnimFig.new("Div", "lat=45,lon=180",  "1e-13", "-8e-13:8e-13", "tz", "lat0lon180", "--exch") 

noAnimFigVars = [
 [xy_Div_SeaSurf_init, xz_Div_lat0_init, tz_SigDot_lat0lon180, tz_Div_lat0lon180], 
 [xy_Div_SeaSurf_init, xz_Div_lat45_init, tz_SigDot_lat45lon180, tz_Div_lat45lon180], 
 [xy_Div_SeaSurf_init, xz_Div_lat0_init, tz_SigDot_lat0lon180, tz_Div_lat0lon180], 
 [xy_Div_SeaSurf_init, xz_Div_lat0_init, tz_SigDot_lat0lon180, tz_Div_lat0lon180], 
 [xy_Div_SeaSurf_init, xz_Div_lat0_init, tz_SigDot_lat0lon180, tz_Div_lat0lon180], 
]

exps = [
 Exp.new("exp1", "./exp1_result/", "Run1_k20_"), 
 Exp.new("exp2", "./exp2_result/", "Run1_k20_"), 
 Exp.new("exp3", "./exp3_result/", "Run1_k20_"), 
 Exp.new("exp4", "./exp4_result/", "Run1_k20_"),
 Exp.new("exp5", "./exp5_result/", "Run1_k20_")
]

animFigVars = [
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfos[0], "xz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=45", "2.5e-14", "-3e-13:3e-13", animTimeInfos[1], "xz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfos[2], "xz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfos[3], "xz", "anim", "") ],
 [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfos[4], "xz", "anim", "") ]
]

####################################################

exps.each_with_index{|exp, i|
  puts "= experiment: #{exp.name}"

  noAnimFigVars[i].each{|figvar|
    exp.add_NoAnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.prefix, figvar.suffix, figvar.gpOpts) 
  }

  animFigVars[i].each{|figvar|
    exp.add_AnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.animTimeInfo, figvar.prefix, figvar.suffix, figvar.gpOpts)
p figvar.animTimeInfo 
  }

  exp.create_allFig
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0
}

