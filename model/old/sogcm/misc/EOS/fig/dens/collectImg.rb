#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

require "fileutils"

DENS_NCFILE="dens.nc"
DENSERROR_NCFILE="densError.nc"

FIXSAL_CUT = [ 0, 10, 20, 30, 35, 40 ]
FIXPTEMP_CUT = [ 0, 10, 20, 30, 40 ]
FIXPRESS_CUT = [ 0, 100, 200, 300, 400, 500 ]

DENS_RANGE = "-10:80"
DENS_INT = "2.5"

DENSERROR_RANGE = "-1.1e-02:1.1e-02"
DENSERROR_INT = "5e-04"

CutFigList = { "sal"=>FIXSAL_CUT, "ptemp"=>FIXPTEMP_CUT, "press"=>FIXPRESS_CUT }

#ENV.store("DCLENVCHAR","_") # ":"をbashでも使えるようにする
#ENV.store("SW_LDUMP","true") # trueならXのダンプファイルを裏画面で作成
#ENV.store("SW_LWAIT","false") # falseなら改ページのタイミングで一時停止しない
#ENV.store("SW_LWAIT1","false") # falseならデバイスをクローズするときに一時停止しない

def createFig(ncFile, varName, varRange, varInt, title, gpOpt="")
  CutFigList.each_pair{ |fixAxis, cut|
    cut.each{|fixVal|
      p "ncFile=#{ncFile}, var=#{varName} : fixAxis=#{fixAxis}, fixVal=#{fixVal}"
      options = "--wsn 2 --range #{varRange} --int #{varInt} --title '#{title}'"
      options << " #{gpOpt}" if gpOpt.length > 0
      `gpview #{ncFile}@#{varName},#{fixAxis}=#{fixVal} #{options}`
      `convert -trim -rotate 90 -density 400x400 -scale 600x600 -units PixelsPerInch dcl.ps #{varName}_#{fixAxis}_#{fixVal}.jpg`
    }
  }
end

["JM95", "LINEAR", "SIMPLENONLINEAR"].each{|eos|
#  createFig(DENS_NCFILE, "EOS_#{eos}", DENS_RANGE, DENS_INT, "density anomaly(sigma)")
}
["Linear", "SimpleNonLinear"].each{|eos|
  createFig(DENSERROR_NCFILE, "densError_#{eos}", DENSERROR_RANGE, DENSERROR_INT, "error relative to density by JM95 EOS", "--clrmap 14")
}


p "create thumbnail.."
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"
FileUtils.mkdir(["figdir", "thum-src"])
FileUtils.cp(DCMODEL_THUM_ORIGIN, "thum-src/")
FileUtils.cp(Dir.glob("{EOS,densError}_*.jpg"), "figdir")
Dir.chdir("./thum-src/"){
`ruby1.8 ./dcmodel-thum.rb`
`ruby1.8 ./dcmodel-thum-make.rb`
}
FileUtils.mv(Dir.glob("./thumbdir/*.png"), ".")
FileUtils.rm_r(["thumbdir", "figdir", "thum-src"])
FileUtils.rm(Dir.glob("sample_thum.htm*"))
FileUtils.rm("thumbdir.SIGEN")

 
