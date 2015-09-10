#!/usr/bin/env ruby
#
require "numru/ggraph"
include NumRu

DATADIR = "./data_T42/"
BCDATADIR="./data_T42/ogcm_BC/"

TauXVarName   = "TauX"
TauYVarName   = "TauY"
RadSDWFlxName = "RadSDWFLXA"
RadLDWFlxName = "RadLDWFLXA"
SensHeatFlxName = "Sens"
LatentHeatFlxName = "Evap"
PrcpName = "PRCP"
SurfTempName = "SurfTemp"

DensFreshWater = UNumeric.new(1e3, "kg.m-3")
LatentHeat = UNumeric.new(2.501e6, "J.kg-1")

###########################

#
varNames = [ \
             TauXVarName, TauYVarName, \
             RadSDWFlxName, RadLDWFlxName, SensHeatFlxName, LatentHeatFlxName, \
             PrcpName, SurfTempName
           ]

#
def output_ncfile(dir, varName, gphys)
  ofile = NetCDF.create("#{dir}#{varName}.nc")
  GPhys::IO.write(ofile, gphys)
  ofile.close
end

#
`mkdir -p #{BCDATADIR}` unless FileTest.exist?("#{BCDATADIR}")

# Time mean

varNames.each{|varName|
  puts "Take average of #{varName}.." 
  var = GPhys::IO.open_gturl("#{DATADIR}#{varName}.nc@#{varName}")
  var = var.cut("sigm"=> 1.0) if var.axnames.include?("sigm")
  var = var.mean("lon").mean("time")

  output_ncfile(BCDATADIR, varName, var)
}

# Heat flux and freshwater flux

puts "Output non-penetrate and penetarate heat flux.."
radSDWFlx = GPhys::IO.open_gturl("#{BCDATADIR}#{RadSDWFlxName}.nc@#{RadSDWFlxName}")
radLDWFlx = GPhys::IO.open_gturl("#{BCDATADIR}#{RadLDWFlxName}.nc@#{RadLDWFlxName}")
sensHFlx = GPhys::IO.open_gturl("#{BCDATADIR}#{SensHeatFlxName}.nc@#{SensHeatFlxName}")
latentHFlx = GPhys::IO.open_gturl("#{BCDATADIR}#{LatentHeatFlxName}.nc@#{LatentHeatFlxName}")

q_ns = - radLDWFlx + sensHFlx + latentHFlx
q_ns.rename("Qns")
q_ns.long_name = "Non penetrtive heat flux(the sum of sensible, latent and long wave heat fluxes"
output_ncfile(BCDATADIR, "Qns", q_ns)

q_sr = - radSDWFlx.copy
q_sr.rename("Qsr")
q_sr.long_name = "Penetrative heat flux"
output_ncfile(BCDATADIR, "Qsr", q_sr)

puts "Output freshwater"
prcp = GPhys::IO.open_gturl("#{BCDATADIR}#{PrcpName}.nc@#{PrcpName}")
freshWaterFlx =  (- prcp - (-latentHFlx/LatentHeat))/DensFreshWater
freshWaterFlx.rename("Fw")
freshWaterFlx.long_name = "Freshwater flux"
output_ncfile(BCDATADIR, "Fw", freshWaterFlx)
