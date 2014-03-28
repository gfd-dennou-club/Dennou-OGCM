#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

include Math

RPalnet = 6.371e06
PI = acos(-1.0)
Depth = 5e03
CourantNum = 0.2
DelX_min = 300e03
DelY_min = 500e03

#############################

class Exp
 def initialize(exp_name, bvFreq, m, n, l)
   @name = exp_name
   @bvFreq = bvFreq
   @m = m
   @n = n 
   @l = l
 end

 def print_IGravWaveProp()
   omega = sprintf("%7.2e", get_IGravWaveFreq())
   period = sprintf("%7.2e", get_IGravWavePeriod())
   delTime = sprintf("%7.2e", get_validDelTime(CourantNum))
   puts "[#{@name}]: N=#{@bvFreq}, m=#{@m}, n=#{@n}, l=#{@l} ::=> frequency=#{omega}, Period=#{period}, DelTime=#{delTime}"
   period_day = get_IGravWavePeriod()/86400.0
   delTime_hr = get_validDelTime(CourantNum)/3600.0
   puts "------------> Period:#{period_day}[day], DelTime:#{delTime_hr}[hr]"
 end

 def get_validDelTime(courantNum)
   cx_wave = RPalnet/(@m * get_IGravWavePeriod())
   cy_wave = RPalnet/(@n * get_IGravWavePeriod())

   return CourantNum* [DelX_min/cx_wave, DelY_min/cy_wave].min
 end

 def get_IGravWaveFreq
   return sqrt(@n*(@n+1))/(RPalnet*PI*@l/Depth) * @bvFreq
 end

 def get_IGravWavePeriod
   return 2*PI/get_IGravWaveFreq()
 end
 

end

exps = []
exps[0] = Exp.new("exp1", 1e-02, 1, 1, 1)
exps[1] = Exp.new("exp2", 1e-02, 1, 2, 1)
exps[2] = Exp.new("exp3", 1e-02, 2, 2, 1)
exps[3] = Exp.new("exp4", 1e-02, 1, 1, 2)
exps[4] = Exp.new("exp5", 1e-03, 1, 1, 1)
exps[5] = Exp.new("exp6", 1e-02, 0, 1, 2)

exps.each{|exp|
 exp.print_IGravWaveProp
}
