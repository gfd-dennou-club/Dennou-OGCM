#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

require "./../collectImgLib"
include CollectImgLib

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"


VarNoAnimFig = Struct.new("VarNoAnimFig", :name, :cutpos, :intrv, :range, :prefix, :suffix)
VarAnimFig = Struct.new("VarAnimFig", :name, :cutpos, :intrv, :range, :animTimeInfo, :prefix, :suffix)


animTimeInfos = [
  AnimTimeInfo.new(0, 72, 2, "hour"), 
  AnimTimeInfo.new(0, 192, 3, "hour"), 
  AnimTimeInfo.new(0, 96, 1.5, "hour") 
]

xy_Psi_SeaSurf_init = VarNoAnimFig.new("Psi", "t=0,sig=0",  "1", "-8:8", "xy", "SeaSurf_init") 
xy_Psi_SeaBtm_init = VarNoAnimFig.new("Psi", "t=0,sig=-1",  "1", "-8:8", "xy", "SeaBtm_init") 
t_Psi_SeaSurf_x0 = VarNoAnimFig.new("Psi", "lat=0,lon=0,sig=0",  "1", "-8:8", "t", "SeaSurf") 
t_Psi2_SeaSurf_x0 = VarNoAnimFig.new("Psi_2", "lat=0,lon=0,sig=0,t=0:3",  "1", "-8:8", "t", "SeaSurf") 
t_Psi_SeaBtm_x0 = VarNoAnimFig.new("Psi", "lat=0,lon=0,sig=-1", "1", "-8:8", "t", "SeaBtm") 
t_Psi_SeaSurf_x45 = VarNoAnimFig.new("Psi", "lat=45,lon=0,sig=0",  "1", "-8:8", "t", "SeaSurf") 
t_Psi2_SeaSurf_x45 = VarNoAnimFig.new("Psi_2", "lat=45,lon=0,sig=0,t=0:3",  "1", "-8:8", "t", "SeaSurf") 
t_Psi_SeaBtm_x45 = VarNoAnimFig.new("Psi", "lat=45,lon=0,sig=-1", "1", "-8:8", "t", "SeaBtm") 

noAnimFigVars = [
 [ xy_Psi_SeaSurf_init, xy_Psi_SeaBtm_init, t_Psi_SeaSurf_x0, t_Psi_SeaBtm_x0, t_Psi2_SeaSurf_x0 ], 
 [ xy_Psi_SeaSurf_init, xy_Psi_SeaBtm_init, t_Psi_SeaSurf_x45, t_Psi_SeaBtm_x45, t_Psi2_SeaSurf_x45], 
 [ xy_Psi_SeaSurf_init, xy_Psi_SeaBtm_init, t_Psi_SeaSurf_x0, t_Psi_SeaBtm_x0, t_Psi2_SeaSurf_x0] 
]

exps = [
 Exp.new("exp1", "./exp1_result/", "Run1_k20_"), 
 Exp.new("exp2", "./exp2_result/", "Run1_k20_"), 
 Exp.new("exp3", "./exp3_result/", "Run1_k20_")
]

animFigFlag = false
animFigVars = [
 [VarAnimFig.new("Psi", "sig=-0.5", "1", "-8:8", animTimeInfos[0], "xy", "anim")], 
 [VarAnimFig.new("Psi", "sig=-0.5", "1", "-8:8", animTimeInfos[1], "xy", "anim")], 
 [VarAnimFig.new("Psi", "sig=-0.5", "1", "-8:8", animTimeInfos[2], "xy", "anim")]
]

####################################################

exps.each_with_index{|exp, i|
  puts "= experiment: #{exp.name}"

  noAnimFigVars[i].each{|figvar|
    exp.add_NoAnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.prefix, figvar.suffix) 
  }

if animFigFlag then
  animFigVars[i].each{|figvar|
    exp.add_AnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.animTimeInfo, figvar.prefix, figvar.suffix)
p figvar.animTimeInfo 
  }
end

  exp.create_allFig
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif")
}

