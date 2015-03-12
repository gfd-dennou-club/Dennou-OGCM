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
LAST_DAY= 3300000  #[day]
#LAST_DAY=730000  #[day]
LASTDATE_SUFFIX="_10000yr"
YT_SUFFIX="_0-10000yr"

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

yt_U_SeaSurf = VarNoAnimFig.new("", "U", "lon=0,sig=0", "0.05", "-1.2:1.2", "xy", "SeaSurf#{YT_SUFFIX}", "--clrmap 10 --title 'U' ")
yt_U_SeaBtm = VarNoAnimFig.new("", "U", "lon=0,sig=-0.995", "0.05", "-1.2:1.2", "xy", "SeaBtm#{YT_SUFFIX}", "--clrmap 10 --title 'U' ")
yz_U_mplane = VarNoAnimFig.new("", "U", "lon=0,t=#{LAST_DAY}", "0.05", "-1.2:1.2", "yz", "mplane#{LASTDATE_SUFFIX}", "--clrmap 10 --title 'U' ")
yz_U_mplane_eq = VarNoAnimFig.new("", "U", "lat=-10:10,lon=0,t=#{LAST_DAY}", "0.05", "-1.2:1.2", "yz", "mplane_eq#{LASTDATE_SUFFIX}", "--clrmap 10 --title 'U' ")

yz_MassStrFunc_mplane = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-90:90", "10", "-100:100", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'MassStreamFunc'")
yz_MassStrFunc_mplane_eq = VarNoAnimFig.new("", "MassStreamFunc", "t=#{LAST_DAY},lat=-10:10", "10", "-100:100", "yz", "mplane_eq#{LASTDATE_SUFFIX}", "--title 'MassStreamFunc'")
yt_PTemp_SeaSurf = VarNoAnimFig.new("", "PTemp", "lon=0,sig=0,t=0:#{LAST_DAY}", "2", "270:310", "tz", "SeaSurf#{YT_SUFFIX}", "--title 'Potential Temperature'")
yt_PTemp_SeaBtm = VarNoAnimFig.new("", "PTemp", "lon=0,sig=-0.995,t=0:#{LAST_DAY}", "2", "270:310", "tz", "SeaBtm#{YT_SUFFIX}", "--title 'Potential Temperature'")
yz_PTemp_mplane = VarNoAnimFig.new("", "PTemp", "lon=0,sig=-1:0,t=#{LAST_DAY}", "2", "270:310", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'Potential Temperature'")
yt_Salt_SeaSurf = VarNoAnimFig.new("", "Salt", "lon=0,sig=0,t=0:#{LAST_DAY}", "0.2", "33:37", "tz", "SeaSurf#{YT_SUFFIX}", "--title 'Salinity'")
yt_Salt_SeaBtm = VarNoAnimFig.new("", "Salt", "lon=0,sig=-0.995,t=0:#{LAST_DAY}", "0.2", "33:37", "tz", "SeaBtm#{YT_SUFFIX}", "--title 'Salinity'")
yz_Salt_mplane = VarNoAnimFig.new("", "Salt", "lon=0,sig=-1:0,t=#{LAST_DAY}", "0.2", "33:37", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'Salinity'")

yz_PressEdd_mplane = VarNoAnimFig.new("", "PressEdd", "lon=0,sig=-1:0,t=#{LAST_DAY}", "2.5e4", "-2e5:6e5", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'Pressure deviation from static pressure'")
yz_DensEdd_mplane = VarNoAnimFig.new("", "DensEdd", "lon=0,sig=-1:0,t=#{LAST_DAY}", "1.5", "-6:30", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'Density deviation from refrence density'")
yz_DensPot_mplane = VarNoAnimFig.new("", "DensPot", "lon=0,sig=-1:0,t=#{LAST_DAY}", "0.25", "-4:4", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'Potential density deviation from refrence density'")

t_KEAvg = VarNoAnimFig.new("", "KEAvg", "", "", "0:0.050", "t", "", "--title 'K.E.(global mean)'" )

yz_BoulusMassStrFunc_mplane = VarNoAnimFig.new("", "BolusMStreamFunc", "t=#{LAST_DAY},lat=-90:90", "10", "-100:100", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'BolusMassStreamFunc'")
yz_ResMassStrFunc_mplane = VarNoAnimFig.new("", "ResMStreamFunc", "t=#{LAST_DAY},lat=-90:90", "10", "-100:100", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'ResidualMassStreamFunc'")

yz_ConvIndex_mplane = VarNoAnimFig.new("", "ConvIndex", "t=#{LAST_DAY},lat=-90:90,lon=0", "0.1", "0:1", "yz", "mplane#{LASTDATE_SUFFIX}", "--title 'ConvectiveIndex'")

exps = [
        Exp.new("exp_EOSL_IDIFF_CA_GM", 
                "#{DATASTORAGE_DIR}/exp_EOSL_IDIFF_CA_GM2/", ""), 
        Exp.new("exp_EOSQ_IDIFF_CA_GM", 
                "#{DATASTORAGE_DIR}/exp_EOSQ_IDIFF_CA_GM2/", ""), 
        Exp.new("exp_EOSJM95_IDIFF_CA_GM", 
                "#{DATASTORAGE_DIR}/exp_EOSJM95_IDIFF_CA_GM2/", ""), 
]

noAnimFigVars = []
exp_noAnimFigs =  [ \
                    t_KEAvg, yt_U_SeaSurf, 
                    yz_U_mplane, yz_U_mplane_eq, \
                    yz_MassStrFunc_mplane, yz_MassStrFunc_mplane_eq, \
                    yt_PTemp_SeaSurf,  yz_PTemp_mplane, \
                    yt_Salt_SeaSurf, yz_Salt_mplane, \
#                    yz_PressEdd_mplane, 
                    yz_DensPot_mplane, 
                    yt_U_SeaBtm, yt_PTemp_SeaBtm, yt_Salt_SeaBtm, 
                    yz_BoulusMassStrFunc_mplane, yz_ResMassStrFunc_mplane,
                    yz_ConvIndex_mplane
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
              Exp.new("exp_comm", "./common/"), 
              Exp.new("exp_EOSComp_CARediGM", "./EOSComp/CARediGM/"), 
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
p expHashKey
  return @expsHash[expHashKey].dirPath
end


t_KECompari_EOSComp = NoAnimOverplotFig.new("KEAvg_EOSCompari", 
                     "#{getExpDirPath("EOSQ_IDIFF_CA_GM")}/KEAvg.nc, \
                      #{getExpDirPath("EOSL_IDIFF_CA_GM")}/KEAvg.nc, \
                      #{getExpDirPath("EOSJM95_IDIFF_CA_GM")}/KEAvg.nc", 
                     "KEAvg,KEAvg,KEAvg", "t=0:1e10", "", "0:0.05", "--title 'K.E.(global mean)'") 
t_KECompari_EOSComp.createFigure(getExpDirPath("EOSComp_CARediGM"))

#t_KECompari_VDiffComp = NoAnimOverplotFig.new("KEAvg_VDiffCompari", 
#                     "#{getExpDirPath("EOSL_HDIFF")}/KEAvg.nc, \
#                      #{getExpDirPath("EOSL_HDIFF_VDIFF100")}/KEAvg.nc", 
#                     "KEAvg,KEAvg", "t=0:1e10", "", "0:0.05", "--title 'K.E.(global mean)'") 
#t_KECompari_VDiffComp.createFigure(getExpDirPath("VDiffComp_noCARediGM"))


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

