#!/usr/bin/env ruby
# -*- coding: euc-jp -*-

require "fileutils"


Dir.glob("{expFunc,waveFunc}*.eps"){|f|
  jpgFile = File.basename(f, "eps") + "jpg"
  p "Convert #{f} -> #{jpgFile}"
  `convert -density 300x300 #{f} #{jpgFile}`
}

p "create thumbnail.."
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"
FileUtils.mkdir(["figdir", "thum-src"])
FileUtils.cp(DCMODEL_THUM_ORIGIN, "thum-src/")
FileUtils.cp(Dir.glob("*.jpg"), "figdir")
Dir.chdir("./thum-src/"){
`ruby1.8 ./dcmodel-thum.rb`
`ruby1.8 ./dcmodel-thum-make.rb`
}
FileUtils.mv(Dir.glob("./thumbdir/*.png"), ".")
FileUtils.rm_r(["thumbdir", "figdir", "thum-src"])
FileUtils.rm(Dir.glob("sample_thum.htm*"))
FileUtils.rm("thumbdir.SIGEN")

 
