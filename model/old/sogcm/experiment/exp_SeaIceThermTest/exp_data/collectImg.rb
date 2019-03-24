#!/usr/bin/env ruby
# -*- coding: euc-jp -*-
#
#     

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"

#################################################################

require "./../../collectImgLib"
include CollectImgLib
require "numru/ggraph"
require "../../GpviewLib.rb"
require "./MyFigGen.rb"

TRANGE="#{99*360}:36000"

#`gpview --range 0:7 --overplot 5 --type 1 --index 12 --aspect 0.9 --var hi,lon=0,lat=0,t=#{TRANGE}  Fbx0.nc Fbx0.5.nc ctrl.nc Fbx2.nc Fbx3.nc`


#################################

exps = [ 
  Exp.new("exp_comm", "./common/"),
              Exp.new("exp_ctrl", "./ctrl/"), 
              Exp.new("exp_FbComp", "./FbComp/"),
              Exp.new("exp_i0Comp", "./i0Comp/"),
              Exp.new("exp_SnowComp", "./SnowComp/"),                             
             ]


####################################

# Register figure objects for each experiment
@expsHash = {}

exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
}


####################################

def getExpDirPath(expHashKey)
p expHashKey
  return @expsHash[expHashKey].dirPath
end


=begin

###########################################

t_SSFlux = NoAnimOverplotFig.new("SSFlux", 
                                 "ctrl.nc, ctrl.nc, ctrl.nc, ctrl.nc",
                                 "LW,SW,SI,LI", "t=#{TRANGE},lon=0,lat=0", "", "50:-400",
                                 "--title 'Forcing at sea ice surface'")
t_SSFlux.createFigure( getExpDirPath("comm"),
                       MyFigGen.new(
                         'Forcing at sea ice surface', 
                         t_SSFlux.varNames.split(","), [1,1,1,1], [12,22,32,42]
                       )
                     )

############ Control run #################

t_hihs = NoAnimOverplotFig.new("hshi",
                               "ctrl_hs_comb.nc, ctrl.nc",
                               "hs,hi", "t=#{TRANGE},lon=0,lat=0", "", "0:4",
                               "--title 'Snow and ice thickness'")

t_hihs.createFigure( getExpDirPath("ctrl"),
                     MyFigGen.new('Snow and ice thickness', ["hs+hi","hi"], [1,2], [12, 12])
                   )


t_IceTemp = NoAnimOverplotFig.new("T1T2",
                                  "ctrl.nc, ctrl.nc",
                                  "T1,T2", "t=#{TRANGE},lon=0,lat=0", "", "-15:0",
                                  "--title 'Ice temperature'")
t_IceTemp.createFigure( getExpDirPath("ctrl"),
                        MyFigGen.new('Ice temperature', t_IceTemp.varNames.split(","), [1,2], [12,12])
                      )


############# Fb dependence #################

=end

t_hi_FbComp = NoAnimOverplotFig.new("hi_FbComp",
                                    "Fbx0.nc, Fbx0.5.nc, ctrl.nc, Fbx2.nc, Fbx3.nc, Fbx4.nc",
                                    "hi,hi,hi,hi,hi,hi", "t=#{TRANGE},lon=0,lat=0", "", "0:7",
                                    "--title 'Ice thickness' --type 1 --index 12 --aspect 0.9"
                                   )
t_hi_FbComp.createFigure( getExpDirPath("FbComp"),
                          MyFigGen.new(
                            'Ice thickness', 
                            ["Fb=0.0", "Fb=1.0", "Fb=2.0", "Fb=4.0", "Fb=6.0", "Fb=8.0"],
                            [1,1,1,1,1,1], [22,32,12,42,52,62], 
                            0.85, {"ylabelint"=>100.0, "ytickint"=>100.0}
                          )
                        )

axis_Fb = Axis.new.set_pos(VArray.new(NArray.to_na([0.0,1.0,2.0,4.0,6.0,8.0]),
             {"long_name"=>"Ocean heat flux", "units"=>"W/m2"}, "Fb")
)
t_hi_FbComp.name = "hi_FbS76Comp"
t_hi_FbComp.createFigure(
  getExpDirPath("FbComp"),
  S76CompariFig.new('Ice thickness',
                    [619.0, 415.0, 287.0, 146.0, 61,0],
                    axis_Fb)
)

############# i0 dependence #################

t_hi_i0Comp = NoAnimOverplotFig.new("hi_i0Comp",
                                    "i0x0.5.nc, ctrl.nc, i0x1.5.nc, i0x2.nc",
                                    "hi,hi,hi,hi", "t=#{TRANGE},lon=0,lat=0", "", "0:5",
                                    "--title 'Ice thickness' --type 1 --index 12 --aspect 1.2"
                                   )
t_hi_i0Comp.createFigure(getExpDirPath("i0Comp"),
                         MyFigGen.new(
                           'Ice thickness', 
                           ["i0=0.085", "i0=0.17", "i0=0.235", "i0=0.34"],
                           [1,1,1,1], [22,12,32,42],                            
                           1.2, {"ylabelint"=>100.0,"ytickint"=>100.0}
                         )
                        )


t_hi_i0Comp.name = "hi_i0S76Comp"
t_hi_i0Comp.createFigure(
  getExpDirPath("i0Comp"), 
  S76CompariFig.new(
    'Ice thickness', [262.0, 287.0, 320.0, 351.0],
    Axis.new.set_pos( VArray.new(NArray.to_na([0.085, 0.17, 0.235, 0.34]),
                      {"long_name"=>"Fraction of penetrating shortwave radiation", "units"=>"(1)"}, "i0") )
  )
)


####### Snowfall dependency

t_hi_SnowComp = NoAnimOverplotFig.new("hi_SnowComp",
                                    "Snowx0.5.nc, ctrl.nc, Snowx1.5.nc, Snowx2.0.nc, Snowx2.5.nc, Snowx3.0.nc",
                                    "hi,hi,hi,hi,hi,hi", "t=#{TRANGE},lon=0,lat=0", "", "0:10",
                                    "--title 'Ice thickness' --type 1 --index 12 --aspect 1.2"
                                   )
t_hi_SnowComp.createFigure(getExpDirPath("SnowComp"),
                         MyFigGen.new(
                           'Ice thickness', 
                           ["sf=20", "sf=40", "sf=60", "sf=80", "sf=100", "sf=120"],
                           [1,1,1,1,1,1], [22,12,32,42,52,62],                            
                           1.2, {"ylabelint"=>100.0,"ytickint"=>100.0}
                         )
                        )

t_hi_SnowComp.name = 'hi_SnowS76Comp'
t_hi_SnowComp.createFigure(
  getExpDirPath("SnowComp"),
  S76CompariFig.new(
    'Ice thickness', [338.0, 287.0, 266.0, 277.0, 395.0, 660.0],
    Axis.new.set_pos( VArray.new(NArray.to_na([20.0, 40.0, 60.0, 80.0, 100.0, 120.0]),
                                 {"long_name"=>"Total yearly snowfall", "units"=>"cm"}, "sf") )
  )
)

##############################################

# Generate thumnails

exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
}
