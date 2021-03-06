#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

# (* exp0: wind stress provided with analystic function.)
# * exp1: wind stress provided with the result of aquaplanet numerical simulation in Marshall et al(2007). 
#    - exp_Ah800L20 (Control Run): The number of vertical grid points is 20(k=20). 
#    - exp_Ah800L60: The number of vertical grid points is 60(k=60).
#    - exp_Ah400L20: The cofficient of horizontal eddy diffusion is a half of one in a control experiement. 
#    - exp_Ah1600L20: The cofficient of horizontal eddy diffusion is twice one in a control experiement. 

##########################
## Configuration
#########################

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"

# End day of time integration
LAST_DAY=109500

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
xy_U_SeaSurf = VarNoAnimFig.new("", "U", "lon=0,sig=0", "0.05", "-1:1", "xy", "SeaSurf", "--clrmap 10 --title 'U' ")
yz_U_mplane = VarNoAnimFig.new("", "U", "lon=0,t=#{LAST_DAY}", "0.05", "-1:1", "yz", "mplane", "--clrmap 10 --title 'U' ")
yz_PressEdd_mplane = VarNoAnimFig.new("", "PressEdd", "lon=0,t=#{LAST_DAY}", "5e04", "-1e06:1e05", "yz", "mplane", "--clrmap 10 --title 'PressEdd' ")
yz_MassStrFunc_mplane = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-90:90", "0.05", "-0.5:0.5", "yz", "mplane", "--title 'MassStreamFunc'")
tz_PTempEdd_SeaSurf = VarNoAnimFig.new("", "PTemp", "lon=0,sig=0", "0.1", "275.5:277.5", "tz", "SeaSurf", "--title 'PTemp'")
yz_PTempEdd_mplane_30yr = VarNoAnimFig.new("", "PTemp", "lon=0,t=10950", "0.1", "275.5:277.5", "yz", "mplane_30yr", "--title 'PTemp'")
yz_PTempEdd_mplane_100yr = VarNoAnimFig.new("", "PTemp", "lon=0,t=36500", "0.1", "275.5:277.5", "yz", "mplane_100yr", "--title 'PTemp'")
yz_PTempEdd_mplane = VarNoAnimFig.new("", "PTemp", "lon=0,t=#{LAST_DAY}", "0.1", "275.5:277.5", "yz", "mplane_300yr", "--title 'PTemp'")
t_KEAvg = VarNoAnimFig.new("EnergyBudget", "KEAvg", "", "", "0:0.05", "t", "", "--title 'K.E.(global mean)'" )


exps = [
#        Exp.new("exp_Kh800T21L20", "#{DATASTORAGE_DIR}/exp_Kh800T21L20/", "Run1_"),      # T21
#        Exp.new("exp_Kh800T21L60", "#{DATASTORAGE_DIR}/exp_Kh800T21L60/", "Run1_"), 
#        Exp.new("exp_Kh400T21L20", "#{DATASTORAGE_DIR}/exp_Kh400T21L20/", "Run1_"),
#        Exp.new("exp_Kh1600T21L20", "#{DATASTORAGE_DIR}/exp_Kh1600T21L20/", "Run1_"),
\
\
##        Exp.new("exp_Kh800T42L20", "#{DATASTORAGE_DIR}/exp_Kh800T42L20/", "Run1_"),      # T42 
        Exp.new("exp_Kh800T42L60", "#{DATASTORAGE_DIR}/exp_Kh800T42L60/", "Run1_"), 
        Exp.new("exp_Kh800Pr10T42L60", "#{DATASTORAGE_DIR}/exp_Kh800Pr10T42L60/", "Run1_"), 
        Exp.new("exp_Kh800HAhT42L60", "#{DATASTORAGE_DIR}/exp_Kh800HAhT42L60/", "Run1_"), 
##        Exp.new("exp_Kh400T42L20", "#{DATASTORAGE_DIR}/exp_Kh400T42L20/", "Run1_"),      
##        Exp.new("exp_Kh1600T42L20", "#{DATASTORAGE_DIR}/exp_Kh1600T42L20/", "Run1_")
\
\
        Exp.new("exp_Kh800T85L60",     "#{DATASTORAGE_DIR}/exp_Kh800T85L60/", "Run1_"),        # T85
        Exp.new("exp_Kh800T85L40",     "#{DATASTORAGE_DIR}/exp_Kh800T85L40/", "Run1_"),        
        Exp.new("exp_Kh800T85L80",     "#{DATASTORAGE_DIR}/exp_Kh800T85L80/", "Run1_"),        
        Exp.new("exp_Kh800Pr10T85L60", "#{DATASTORAGE_DIR}/exp_Kh800Pr10T85L60/", "Run1_"), 
        Exp.new("exp_Kh800HAhT85L60",  "#{DATASTORAGE_DIR}/exp_Kh800HAhT85L60/", "Run1_"), 
\
\
        Exp.new("exp_Kh800T170L60",     "#{DATASTORAGE_DIR}/exp_Kh800T170L60/", "Run1_"),       # T170
        Exp.new("exp_Kh800Pr10T170L60", "#{DATASTORAGE_DIR}/exp_Kh800Pr10T170L60/", "Run1_"), 
        Exp.new("exp_Kh800HAhT170L60",  "#{DATASTORAGE_DIR}/exp_Kh800HAhT170L60/", "Run1_") 
]

@expsHash = {}

noAnimFigVars = []
exp_noAnimFigs =  [ t_KEAvg, xy_U_SeaSurf, yz_U_mplane, yz_MassStrFunc_mplane, 
                    yz_PTempEdd_mplane_30yr, yz_PTempEdd_mplane_100yr, yz_PTempEdd_mplane, tz_PTempEdd_SeaSurf, 
                    yz_PressEdd_mplane ]

########

animTimeInfo = AnimTimeInfo.new(0, LAST_DAY, 500.0, "day")
animFigFlag = false

animFigVars = []
exp_AnimFigs =  [ VarAnimFig.new("SigDot", "lat=0", "2.5e-14", "-3e-13:3e-13", animTimeInfo, "yz", "anim", "") ]


##########

# Register figure objects for each experiment

exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
  noAnimFigVars[i] = exp_noAnimFigs
  animFigVars[i] = exp_AnimFigs
}


#################################

#=begin
extra_exps = [ Exp.new("exp_comm", "./common/"), 
               Exp.new("exp_LComp", "./LCompare/"),
               Exp.new("exp_KhComp", "./KhCompare/"),
               Exp.new("exp_HComp", "./HCompare/"), 
               Exp.new("exp_HViscComp", "./HViscCompare/")
             ]


def getExpDirPath(expHashKey)
  return @expsHash[expHashKey].dirPath
end

t_KECompari_LComp = NoAnimOverplotFig.new("KEAvg_LCompari", 
                     "#{getExpDirPath("Kh800T85L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T85L40")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T85L80")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.01", "--title 'K.E.(global mean)'") 
t_KECompari_KhComp = NoAnimOverplotFig.new("KEAvg_KhCompari", 
                     "#{exps[4].dirPath}/Run1_EnergyBudget.nc,#{exps[6].dirPath}/Run1_EnergyBudget.nc,#{exps[7].dirPath}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.04", "--title 'K.E.(global mean)'") 
t_KECompari_HComp = NoAnimOverplotFig.new("KEAvg_HCompari", 
                     "#{getExpDirPath("Kh800T85L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T42L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T170L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.01", "--title 'K.E.(global mean)'") 

t_KECompari_HViscCompT42 = NoAnimOverplotFig.new("KEAvg_HCompari_T42", 
                     "#{getExpDirPath("Kh800Pr10T42L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800HAhT42L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T42L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.01", "--title 'K.E.(global mean)'") 
t_KECompari_HViscCompT85 = NoAnimOverplotFig.new("KEAvg_HCompari_T85", 
                     "#{getExpDirPath("Kh800Pr10T85L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800HAhT85L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T85L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.01", "--title 'K.E.(global mean)'") 
t_KECompari_HViscCompT170 = NoAnimOverplotFig.new("KEAvg_HCompari_T170", 
                     "#{getExpDirPath("Kh800Pr10T170L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800HAhT170L60")}/Run1_EnergyBudget.nc, \
                      #{getExpDirPath("Kh800T170L60")}/Run1_EnergyBudget.nc", 
                     "KEAvg,KEAvg,KEAvg", "", "", "0:0.01", "--title 'K.E.(global mean)'") 


t_KECompari_LComp.createFigure(extra_exps[1].dirPath)
t_KECompari_KhComp.createFigure(extra_exps[2].dirPath)
t_KECompari_HComp.createFigure(extra_exps[3].dirPath)
t_KECompari_HViscCompT42.createFigure(extra_exps[4].dirPath)
t_KECompari_HViscCompT85.createFigure(extra_exps[4].dirPath)
t_KECompari_HViscCompT170.createFigure(extra_exps[4].dirPath)

[ "Kh800T42L60", 
  "Kh800T85L40", "Kh800T85L60", "Kh800T85L80", 
  "Kh800T170L60", 
  "Kh800Pr10T42L60", "Kh800Pr10T85L60", "Kh800Pr10T170L60", 
  "Kh800HAhT42L60", "Kh800HAhT85L60", "Kh800HAhT170L60" 
].each{|expName|
  exp = @expsHash[expName]

  ncFilePaths = ""; varNames = ""
  ["KEGenNet","KEInputSurfAvg","VDiffDispAvg","PE2KEAvg","HDiffDispAvg","AdvectWorkAvg"].each{|var|
    ncFilePaths << "#{exp.dirPath}/Run1_EnergyBudget.nc,"
    varNames << "#{var},"
  }
  ncFilePaths.chop!; varNames.chop!;

  t_KEBudget = NoAnimOverplotFig.new("KEBudget", ncFilePaths, varNames, "", "", "-5e-09:5e-09", "--title 'K.E. Budget'")
  t_KEBudget.createFigure(exp.dirPath)

  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
  p "Move created figures from '#{exp.dirPath}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.mv(Dir.glob("#{exp.dirPath}/*.{png,jpg}"), "#{FIGSTORAGE_DIR}/#{exp.name}/")
}

extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
}
exit
#=end

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
  # exp.create_thumb(DCMODEL_THUM_ORIGIN, "png") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0 and animFigFlag

  p "Move created figures from '#{exp.dirPath}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.mv(Dir.glob("#{exp.dirPath}/*.{png,jpg}"), "#{FIGSTORAGE_DIR}/#{exp.name}/")


  config_nml = Dir.glob("#{exp.dirPath}/config*.nml")[0]
  p "Move c '#{config_nml}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
  FileUtils.cp(config_nml, "#{FIGSTORAGE_DIR}/#{exp.name}/config.nml")

}

