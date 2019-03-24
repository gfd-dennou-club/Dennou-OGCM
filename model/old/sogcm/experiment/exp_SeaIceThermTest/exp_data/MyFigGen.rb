#!/usr/bin/env ruby
# -*- coding: euc-jp -*-
#
#     
require "numru/ggraph"
include NumRu

def DCLInitialize(iws=4, t0, t1, minv, maxv)
  
  # Setting same as gpview
  DCL.swlstx('lwnd', false)    
  DCL.gropn(iws)
  DCL::gllset('LMISS', true)
  GGraph.set_fig('viewport'=> GpviewLib.set_vpsize([0.15,0.85,0.2,0.55], @aspect),
                   'window'=>[t0, t1, minv, maxv])
  DCL.sgpset('lcntl', false) ;
  DCL.sgpset('lfull', true)               # use full area in the window
  DCL.sgpset('lfprop',true)               # use proportional font
  DCL.uscset('cyspos', 'B' )              # move unit y axis
  #    DCL.sgiset('ifont', 2)
end

class MyFigGen < FigureGenCmd

  include NumRu

  attr_accessor :aspect, :axesOpts
  def initialize(title, labels, lineTypes, lineIndexs, aspectRatio=2.0, axesOpts=Hash.new)
    @title = title
    @labels = labels
    @lineTypes = lineTypes
    @lineIndexs = lineIndexs
    @aspect = aspectRatio
    @axesOpts = axesOpts
    @semtner76Sols = []
  end

  def execute(fig)

    #
    ncFilePaths =  fig.ncFilePaths.gsub(" ", "").split(",")
    varNames = fig.varNames.gsub(" ", "").split(",")
    minv, maxv = fig.range.split(":")
    minv = minv.to_f; maxv = maxv.to_f
    t0str, t1str = TRANGE.split(":")
    t0 = t0str.to_f; t1 = t1str.to_f
    
    
    create_SeasonalCycleFig(fig, ncFilePaths, varNames, t0+15.0, t1-15.0, minv, maxv)

    #
    fig.convertOpt = "#{fig.convertOpt} -quality 50"
  end

  def create_SeasonalCycleFig(fig, ncFilePaths, varNames, t0, t1, minv, maxv)
    
    # setting for lateral axes
    nx1 = 12 # the number of tick
    nx2 = 12 # the number of label
    rx1 = NArray.sfloat(nx1).indgen(0.0, 30.0) + t0
    rx2 = NArray.sfloat(nx2).indgen(0.0, 30.0) + t0
    cx2 = "J F M A M J J A S O N D".split(" ")

    #
    DCLInitialize(4, t0, t1, minv, maxv)    
    DCL.uzfact(0.7)
    GGraph.set_axes({'xside'=>'t', 'yside'=>'lr'}.merge(@axesOpts)) # axes at bottom is not drawn

    newframe = true
    for i in 0..ncFilePaths.length-1
      gp = GPhys::IO.open_gturl("#{ncFilePaths[i]}@#{varNames[i]},t=#{t0}:#{t1},lon=0,lat=0").copy
      gp.name = @labels[i]
      if(gp.units.to_s == "m") then
        gp = gp.convert_units("cm")
        minv *= 100.0; maxv *=100.0
      end
      lineOpts = {"title"=>@title, "min"=>minv, "max"=>maxv, "index"=>@lineIndexs[i], "type"=>@lineTypes[i], "legend"=>true, "annotate"=>false}
      GGraph.line(gp, newframe, lineOpts)
      newframe = false
    end

    # draw axes at bottom 
    DCL::uxaxlb('B', rx1, rx2, cx2, 1)
    DCL::uxsttl('B', 'time', 0.0)

    DCL.uzfact(10.0/7.0)
    DCL.grcls
  end
end

class S76CompariFig < FigureGenCmd
  def initialize(title, s76sols, axis_x, aspectRatio=2.0)
    @title = title
    @s76sols = s76sols
    @axis_x = axis_x
    @aspect = aspectRatio
  end

  def execute(fig)
    #
    ncFilePaths =  fig.ncFilePaths.gsub(" ", "").split(",")
    varNames = fig.varNames.gsub(" ", "").split(",")
    minv, maxv = fig.range.split(":")
    minv = minv.to_f; maxv = maxv.to_f
    t0str, t1str = TRANGE.split(":")
    t0 = t0str.to_f; t1 = t1str.to_f

    #
    mysols = []; units = nil
    ncFilePaths.each_with_index{|ncFilePath, i|
      gp = GPhys::IO.open_gturl("#{ncFilePaths[i]}@#{varNames[i]},t=#{t0}:#{t1},lon=0,lat=0").copy
      units = gp.units.to_s
      if units == "m" then
        gp = gp.convert_units("cm"); units = "cm"
      end
      mysols << gp.mean('t').val
    }
    
    mysols_gp = ary2GPhys(mysols, varNames[0], units)
    s76sols_gp = ary2GPhys(@s76sols, "#{varNames[0]}S76", units)

    minv *= 100.0; maxv *= 100.0 if units == "cm"
    DCLInitialize(4, @axis_x.pos.min.to_f, @axis_x.pos.max.to_f, minv, maxv)
    DCL.uzfact(0.7)

    GGraph.set_axes({'xside'=>'tb', 'yside'=>'lr'})    

    GGraph.line(mysols_gp, true, {"type"=>1, "index"=>12, "title"=>@title})
    GGraph.mark(mysols_gp, false, {"type"=>4, "index"=>12})
    GGraph.line(s76sols_gp, false, {"type"=>2, "index"=>12})
    GGraph.mark(s76sols_gp, false, {"type"=>5, "index"=>12})
    DCL.uzfact(10.0/7.0)
    DCL.grcls    

  end

  def ary2GPhys(ary, varName, units)
    na = NArray.to_na(ary)
    vary = VArray.new(na, {"long_name"=>varName, "units"=>units}, varName)
    return GPhys.new(Grid.new(@axis_x), vary)
  end

end
  
