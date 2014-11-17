#!/usr/bin/env ruby
# -*- coding: euc-jp -*-
#
#     

##########################
## Configuration
#########################

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"

# End day of time integration
#LAST_DAY=109600  #[day]
LAST_DAY=102000  #[day]

# Path of parent directory where the simulation results are saved.
DATASTORAGE_DIR="/home/ykawai/exp_backup/APEOGCirc/"

# Path of parent directory where the figures of simulation results are saved.
FIGSTORAGE_DIR="./"

#################################################################

require "./../../collectImgLib"
include CollectImgLib


VarNoAnimFig = Struct.new("VarNoAnimFig", :ncFileSuffix, :name, :cutpos, :intrv, :range, :prefix, :suffix, :gpOpts)
VarAnimFig = Struct.new("VarAnimFig", :name, :cutpos, :intrv, :range, :animTimeInfo, :prefix, :suffix, :gpOpts)


xy_Div_SeaSurf_init = VarNoAnimFig.new("", "Div", "t=0,sig=0",  "1e-13", "-8e-13:8e-13", "xy", "SeaSurf_init","") 
xz_Div_lat0_init = VarNoAnimFig.new("", "Div", "t=0,lat=0",  "1e-13", "-8e-13:8e-13", "xz", "lat0_init","") 
xz_Div_lat45_init = VarNoAnimFig.new("", "Div", "t=0,lat=45",  "1e-13", "-8e-13:8e-13", "xz", "lat45_init","") 

yt_U_SeaSurf = VarNoAnimFig.new("", "U", "lon=0,sig=0", "0.05", "-1.2:1.2", "xy", "SeaSurf", "--clrmap 10 --title 'U' ")
yz_U_mplane = VarNoAnimFig.new("", "U", "lon=0,t=#{LAST_DAY}", "0.05", "-1.2:1.2", "yz", "mplane_300yr", "--clrmap 10 --title 'U' ")
yz_U_mplane_eq = VarNoAnimFig.new("", "U", "lat=-10:10,lon=0,t=#{LAST_DAY}", "0.05", "-1.2:1.2", "yz", "mplane_eq_300yr", "--clrmap 10 --title 'U' ")

yz_MassStrFunc_mplane = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-90:90", "10", "-100:100", "yz", "mplane_300yr", "--title 'MassStreamFunc'")
yz_MassStrFunc_mplane_eq = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-10:10", "10", "-100:100", "yz", "mplane_eq_300yr", "--title 'MassStreamFunc'")
yt_PTemp_SeaSurf = VarNoAnimFig.new("", "PTemp", "lon=0,sig=0,t=0:#{LAST_DAY}", "2", "270:310", "tz", "SeaSurf", "--title 'Potential Temperature'")
yz_PTemp_mplane = VarNoAnimFig.new("", "PTemp", "lon=0,sig=-1:0,t=#{LAST_DAY}", "2", "270:310", "yz", "mplane_300yr", "--title 'Potential Temperature'")
yt_Salt_SeaSurf = VarNoAnimFig.new("", "Salt", "lon=0,sig=0,t=0:#{LAST_DAY}", "0.2", "33:37", "tz", "SeaSurf", "--title 'Salinity'")
yz_Salt_mplane = VarNoAnimFig.new("", "Salt", "lon=0,sig=-1:0,t=#{LAST_DAY}", "0.2", "33:37", "yz", "mplane_300yr", "--title 'Salinity'")

yz_PressEdd_mplane = VarNoAnimFig.new("", "PressEdd", "lon=0,sig=-1:0,t=#{LAST_DAY}", "2.5e4", "-2e5:6e5", "yz", "mplane_300yr", "--title 'Pressure deviation from static pressure'")
yz_DensEdd_mplane = VarNoAnimFig.new("", "DensEdd", "lon=0,sig=-1:0,t=#{LAST_DAY}", "1.5", "-6:30", "yz", "mplane_300yr", "--title 'Density deviation from refrence density'")
yz_DensPot_mplane = VarNoAnimFig.new("", "DensPot", "lon=0,sig=-1:0,t=#{LAST_DAY}", "0.25", "-3:3", "yz", "mplane_300yr", "--title 'Potential density deviation from refrence density'")

t_KEAvg = VarNoAnimFig.new("", "KEAvg", "", "", "0:0.02", "t", "", "--title 'K.E.(global mean)'" )


exps = [
#        Exp.new("exp_constSalt",     "#{DATASTORAGE_DIR}/exp_APEOGCirc_constSalt/", ""),      
#        Exp.new("exp_Ah1e3Prh1Prv1", 
#                "#{DATASTORAGE_DIR}/exp_APEOGCirc_noConstSalt/exp_Ah1e3Prh1Prv1_Pl171L60/", ""),      
#        Exp.new("exp_Ah1e4Prh1Prv1", 
#                "#{DATASTORAGE_DIR}/exp_APEOGCirc_noConstSalt/exp_Ah1e4Prh1Prv1_Pl171L60/", ""),      
#        Exp.new("exp_Ah1e5Prh1Prv1", 
#                "#{DATASTORAGE_DIR}/exp_APEOGCirc_noConstSalt/exp_Ah1e5Prh1Prv1_Pl171L60/", ""),      
#        Exp.new("exp_Ah1e4Prh10Prv10", 
#                "#{DATASTORAGE_DIR}/exp_APEOGCirc_noConstSalt/exp_Ah1e4Prh10Prv10Pl171L60/", ""),      
        Exp.new("exp_Ah1e4Prh10Prv10GM", 
                "#{DATASTORAGE_DIR}/exp_APEOGCirc_noConstSalt/exp_Ah1e4Prh10Prv10GMPl171L60/", ""),      ]

noAnimFigVars = []
exp_noAnimFigs =  [ \
                    t_KEAvg, \
                    yt_U_SeaSurf, yz_U_mplane, yz_U_mplane_eq, \
                    yz_MassStrFunc_mplane, yz_MassStrFunc_mplane_eq, \
                    yt_PTemp_SeaSurf,  yz_PTemp_mplane, \
                    yt_Salt_SeaSurf, yz_Salt_mplane, \
#                    yz_PressEdd_mplane, 
                    yz_DensEdd_mplane, yz_DensPot_mplane, 
                  ]

# 
animTimeInfo = AnimTimeInfo.new(0, LAST_DAY, 1500.0, "day")
animFigFlag = false

animFigVars = []
exp_AnimFigs =  [ 
                 VarAnimFig.new("U", "lon=0", "0.05", "-1.2:1.2", animTimeInfo, "yz", "anim", ""),  
                 VarAnimFig.new("PTemp", "lon=0,sig=-1:0", "2", "270:310", animTimeInfo, "yz", "anim", ""),
                 VarAnimFig.new("Salt", "lon=0,sig=-1:0", "0.2", "33:37", animTimeInfo, "yz", "anim", ""),
                 VarAnimFig.new("MassStreamFunc", "lat=-90:90,sig=-1:0", "10", "-100:100", animTimeInfo, "yz", "anim", "") 
                ]




#################################

extra_exps = [ 
#              Exp.new("exp_comm", "./common/"), 
#              Exp.new("exp_HViscDiffComp", "./HViscDiffComp/"), 
             ]


####################################

# Register figure objects for each experiment
@expsHash = {}

exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
  noAnimFigVars[i] = exp_noAnimFigs
  animFigVars[i] = exp_AnimFigs
}

extra_exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
}


####################################

def getExpDirPath(expHashKey)
  return @expsHash[expHashKey].dirPath
end

=begin
t_KECompari_LComp = NoAnimOverplotFig.new("KEAvg_LCompari", 
                     "#{getExpDirPath("Ah1e4Pl341L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e4Pl341L30")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e4Pl341L120")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "t=0:4000", "", "0:0.006", "--title 'K.E.(global mean)'") 
t_KECompari_LComp.createFigure(getExpDirPath("LComp"))

t_KECompari_HComp = NoAnimOverplotFig.new("KEAvg_HCompari", 
                     "#{getExpDirPath("Ah1e4Pl341L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e4Pl170L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e4Pl682L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "t=0:4000", "", "0:0.006", "--title 'K.E.(global mean)'") 
t_KECompari_HComp.createFigure(getExpDirPath("HComp"))


t_KECompari_HViscComp = NoAnimOverplotFig.new("KEAvg_HViscDiffCompari", 
                     "#{getExpDirPath("Ah1e4Prh1Prv1")}/KEAvg.nc, \
                      #{getExpDirPath("Ah1e3Prh1Prv1")}/KEAvg.nc, \
                      #{getExpDirPath("Ah1e5Prh1Prv1")}/KEAvg.nc", 
                     "KEAvg,KEAvg,KEAvg", "t=0:#{LAST_DAY}", "", "0:0.025", "--title 'K.E.(global mean)'") 
t_KECompari_HViscComp.createFigure(getExpDirPath("HViscDiffComp"))
=end

####################################################

extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
#  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png")
}
#exit

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
  end

  exp.create_allFig
  # exp.create_thumb(DCMODEL_THUM_ORIGIN, "png") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0 and animFigFlag

  p "Move created figures from '#{exp.dirPath}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.mv(Dir.glob("#{exp.dirPath}/*.{png,jpg,gif}"), "#{FIGSTORAGE_DIR}/#{exp.name}/")


  config_nml = Dir.glob("#{exp.dirPath}/config*.nml")[0]
  p "Move c '#{config_nml}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.cp(config_nml, "#{FIGSTORAGE_DIR}/#{exp.name}/config.nml")

}

