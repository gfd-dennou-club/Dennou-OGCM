##!/usr/bin/env ruby
#

###############################################################

CONV_HRES_PROG = "../util/bin/conv_hres"
CONF_FILENAME  = "conv_hres.conf"

CONV_HRES_VARS = [ 'PRCP', 'TauX', 'TauY', 'RadSDWFLXA', 'RadLDWFLXA', 'Sens', 'Evap' ]
POSITIVE_VARS  = [ 'PRCP' ]

HRES_IMAX  = 128
TIME_START = 730
TIME_END   = 1096

INNCFILE_DIR  = "./data_T21/"
OUTNCFILE_DIR = "./data_T42/"


#############################################################


def create_namelist(inncfn, outncfn, positVarNames, nmlName=CONF_FILENAME)

  positVarNamesStr = ''
  positVarNames.each{|var|
    positVarNamesStr << "'#{var.to_s}', "
  }

  nmlFile = File.open(nmlName, "w")
  nmlFile.print <<"EOS"
&file inncfn = '#{inncfn}', outncfn = '#{outncfn}' /
&time ts = #{TIME_START}, te=#{TIME_END}  /
&hres imax=#{HRES_IMAX} /
&varsign PositiveVar=#{positVarNamesStr} /
EOS
  nmlFile.close
  
end


CONV_HRES_VARS.each{|var|
  inncfile="#{INNCFILE_DIR}#{var}.nc"
  outncfile="#{OUTNCFILE_DIR}#{var}.nc"

  puts "convert #{inncfile} ----> #{outncfile}"

  positVarNames = []
  positVarNames << ["#{var}"] if POSITIVE_VARS.include?(var)
  create_namelist(inncfile, outncfile, positVarNames, CONF_FILENAME)
  
  puts `#{CONV_HRES_PROG}`
  `rm #{CONF_FILENAME}`
}


