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
LAST_DAY=10950

# Path of parent directory where the simulation results are saved.
DATASTORAGE_DIR="/home/ykawai/exp_backup/current/"

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

#tz_SigDot_lat0lon180 = VarNoAnimFig.new("SigDot", "lat=0,lon=180",  "2.5e-14", "-3e-13:3e-13", "tz", "lat0lon180", "--exch")
xy_U_SeaSurf = VarNoAnimFig.new("", "U", "lon=0,sig=0", "0.025", "-0.5:0.5", "xy", "SeaSurf", "--clrmap 10 --title 'U' ")
yz_U_mplane = VarNoAnimFig.new("", "U", "lon=0,t=#{LAST_DAY}", "0.025", "-0.5:0.5", "yz", "mplane", "--clrmap 10 --title 'U' ")
yz_U_mplane_eq = VarNoAnimFig.new("", "U", "lat=-10:10,lon=0,t=#{LAST_DAY}", "0.01", "-0.2:0.2", "yz", "mplane_eq", "--clrmap 10 --title 'U' ")
yz_PressEdd_mplane = VarNoAnimFig.new("", "PressEdd", "lon=0,t=#{LAST_DAY}", "5e04", "-1e06:1e05", "yz", "mplane", "--clrmap 10 --title 'PressEdd' ")
yz_MassStrFunc_mplane = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-90:90", "10", "-120:120", "yz", "mplane", "--title 'MassStreamFunc'")
yz_MassStrFunc_mplane_eq = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-10:10", "20", "-220:220", "yz", "mplane_eq", "--title 'MassStreamFunc'")
tz_PTempEdd_SeaSurf = VarNoAnimFig.new("", "PTemp", "lon=0,sig=0", "0.1", "275.5:277.5", "tz", "SeaSurf", "--title 'PTemp'")
yz_PTempEdd_mplane_30yr = VarNoAnimFig.new("", "PTemp", "lon=0,t=10950", "0.1", "275.5:277.5", "yz", "mplane_30yr", "--title 'PTemp'")
yz_PTempEdd_mplane_100yr = VarNoAnimFig.new("", "PTemp", "lon=0,t=36500", "0.1", "275.5:277.5", "yz", "mplane_100yr", "--title 'PTemp'")
yz_PTempEdd_mplane = VarNoAnimFig.new("", "PTemp", "lon=0,t=#{LAST_DAY}", "0.1", "275.5:277.5", "yz", "mplane_300yr", "--title 'PTemp'")
t_KEAvg = VarNoAnimFig.new("EnergyBudget", "KEAvg", "", "", "0:0.01", "t", "", "--title 'K.E.(global mean)'" )


exps = [
        Exp.new("exp_Ah1e4Pl170L60",     "#{DATASTORAGE_DIR}/exp_Ah1e4Pl170L60/", "Run1_"),      
        Exp.new("exp_Ah1e4Pl341L60",     "#{DATASTORAGE_DIR}/exp_Ah1e4Pl341L60/", "Run1_"),      
        Exp.new("exp_Ah1e4Pl682L60",     "#{DATASTORAGE_DIR}/exp_Ah1e4Pl682L60/", "Run1_"),      
        Exp.new("exp_Ah1e4Pl341L30",     "#{DATASTORAGE_DIR}/exp_Ah1e4Pl341L30/", "Run1_"),      
        Exp.new("exp_Ah1e4Pl341L120",     "#{DATASTORAGE_DIR}/exp_Ah1e4Pl341L120/", "Run1_"),      
        Exp.new("exp_Ah1e3Pl341L60",     "#{DATASTORAGE_DIR}/exp_Ah1e3Pl341L60/", "Run1_"),      
        Exp.new("exp_Ah1e5Pl341L60",     "#{DATASTORAGE_DIR}/exp_Ah1e5Pl341L60/", "Run1_"),      
]

noAnimFigVars = []
exp_noAnimFigs =  [ \
                    t_KEAvg, xy_U_SeaSurf, yz_U_mplane, yz_MassStrFunc_mplane, yz_PressEdd_mplane, \
                     yz_U_mplane, yz_U_mplane_eq, \
                      yz_MassStrFunc_mplane, yz_MassStrFunc_mplane_eq
                  ]

# 
animTimeInfo = AnimTimeInfo.new(0, LAST_DAY, 500.0, "day")
animFigFlag = false

animFigVars = []
exp_AnimFigs =  [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ]




#################################

extra_exps = [ Exp.new("exp_comm", "./common/"), 
               Exp.new("exp_LComp", "./LCompare/"),
               Exp.new("exp_HComp", "./HCompare/"), 
               Exp.new("exp_HViscComp", "./HViscCompare/"), 
               Exp.new("exp_RefSolComp", "./exp_RefSolComp")
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

t_KECompari_HViscComp = NoAnimOverplotFig.new("KEAvg_HViscCompari", 
                     "#{getExpDirPath("Ah1e4Pl341L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e3Pl341L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Ah1e5Pl341L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "t=0:4000", "", "0:0.006", "--title 'K.E.(global mean)'") 
t_KECompari_HViscComp.createFigure(getExpDirPath("HViscComp"))

extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
}

exit

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
  end

  exp.create_allFig
  # exp.create_thumb(DCMODEL_THUM_ORIGIN, "png") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0 and animFigFlag

  p "Move created figures from '#{exp.dirPath}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.mv(Dir.glob("#{exp.dirPath}/*.{png,jpg}"), "#{FIGSTORAGE_DIR}/#{exp.name}/")


  config_nml = Dir.glob("#{exp.dirPath}/config*.nml")[0]
  p "Move c '#{config_nml}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.cp(config_nml, "#{FIGSTORAGE_DIR}/#{exp.name}/config.nml")

}

