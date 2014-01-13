module CollectImgLib

require "fileutils"

class Figure
  attr_accessor :name, :ncFilePath, :varName, :cutpos, :gpviewOpt, :convertOpt, :figPicExt

  def initialize(figname, ncFilePath, varName, cutpos, intrv="", range="", gpOpts="")
    @name = figname
    @ncFilePath = ncFilePath
    @varName = varName
    @cutpos = cutpos

    @gpviewOpt = "--wsn 2 "
    @gpviewOpt << " --range #{range}" if range.length != 0
    @gpviewOpt << " --int #{intrv}" if intrv.length != 0
    @gpviewOpt << " #{gpOpts}" if gpOpts.length != 0

    @convertOpt = @@DEFAULT_CONVERT_OPT
  end

  def createFigure(dirPath)
  end

  @@DEFAULT_CONVERT_OPT = "-rotate 90 -trim -density 400x400 -scale 600x600 -units PixelsPerInch"
end


AnimTimeInfo = Struct.new("AnimTimeInfo", :init, :end, :intrv, :unit)

class AnimFig < Figure
  attr_accessor :animTimeInfo

  def initialize(figname, ncFilePath, varName, cutpos, animTimeInfo, intrv="", range="", gpOpts="")
    super(figname, ncFilePath, varName, cutpos, intrv, range, gpOpts)
    @figPicExt = "gif"
    @animTimeInfo = animTimeInfo
  end

  def createFigure(dirPath)
    for t in 0..(@animTimeInfo.end - @animTimeInfo.init)/@animTimeInfo.intrv
      time = @animTimeInfo.init + @animTimeInfo.intrv*t
      p "command gpview #{@ncFilePath}@#{@varName},#{@cutpos},t=#{time} #{@gpviewOpt}`"
      `gpview #{@ncFilePath}@#{@varName},#{@cutpos},t=#{time} #{@gpviewOpt}`
      
      timePadd = format("%04d", "#{time}".to_i)
      `mv dcl.ps #{@name}_#{timePadd}.ps`
    end

    animFigPath = "#{dirPath}#{@name}_#{animTimeInfo.init}-#{@animTimeInfo.end}.#{@figPicExt}"
    p "convert.. ps => gif animation (options #{@convertOpt})"
    `convert #{@convertOpt} #{@name}_*.ps #{animFigPath}`
    `rm #{@name}_*.ps` 
  end 

end

class NoAnimFig < Figure
  def initialize(figname, ncFilePath, varName, cutpos, intrv="", range="", gpOpts="")
    super
    @figPicExt = "jpg"
  end

  def createFigure(dirPath)

    p "command: gpview #{@ncFilePath}@#{@varName},#{@cutpos} #{@gpviewOpt}"
    `gpview #{@ncFilePath}@#{@varName},#{@cutpos} #{@gpviewOpt}`
    p "convert.. ps => #{figPicExt} (options #{@convertOpt}"    
    `convert #{@convertOpt} dcl.ps #{dirPath}#{@name}.#{figPicExt}`
    `rm dcl.ps`
  end 
end

class Exp
  attr_accessor :name, :dirPath, :ncFilePrefix

  def initialize(name="", dirPath="./", ncFilePrefix="")
    @name = name
    @dirPath = dirPath
    @ncFilePrefix = ncFilePrefix
    @figList = []
  end

  def add_NoAnimFig(varName, cutPos, figIntrv, figRange, prefix="", suffix="", gpOpts="", figPicExt="" )

    figName, ncFilePath = create_FigMetaInfo(varName, prefix, suffix)

    fig = NoAnimFig.new(figName, ncFilePath, varName, cutPos, figIntrv, figRange, gpOpts)
    fig.figPicExt = figPicExt if figPicExt.length != 0
    @figList.push(fig)
  end

  def add_AnimFig(varName, cutPos, figIntrv, figRange, animTimeInfo, prefix="", suffix="", gpOpts="")

    figName, ncFilePath = create_FigMetaInfo(varName, prefix, suffix)

    fig = AnimFig.new(figName, ncFilePath, varName, cutPos, animTimeInfo, figIntrv, figRange, gpOpts)
    @figList.push(fig)
    return figName
  end

  def create_allFig(dirPath=@dirPath)
    @figList.each{|fig|
      puts "Create figure `#{fig.name}`..." 
      fig.createFigure(dirPath)
    }
  end

  def create_FigMetaInfo(varName, prefix, suffix)
    figName = ""
    figName << "#{prefix}_" if prefix.length != 0
    figName << varName
    figName << "_#{suffix}" if suffix.length != 0
    
    ncFilePath = "#{@dirPath}#{ncFilePrefix}#{varName}.nc"

    return figName, ncFilePath
  end

  def create_thumb(dcmodel_thumb_rb, targetFileExt, rb18cmd="ruby1.8")

    p "chdir experiment dir `#{@dirPath}`.."

    FileUtils.chdir(@dirPath){ 
      p "create thumbnail.."

      FileUtils.mkdir(["figdir", "thum-src"])
      FileUtils.cp(dcmodel_thumb_rb, "thum-src/dcmodel-thum.rb")
      FileUtils.cp(Dir.glob("*.#{targetFileExt}"), "figdir")
      Dir.chdir("./thum-src/"){
        `#{rb18cmd} ./dcmodel-thum.rb`
        `#{rb18cmd} ./dcmodel-thum-make.rb`
      }
      FileUtils.mv(Dir.glob("./thumbdir/*.png"), ".")
      FileUtils.rm_r(["thumbdir", "figdir", "thum-src"])
      FileUtils.rm(Dir.glob("sample_thum.htm*"))
      FileUtils.rm("thumbdir.SIGEN")
    }
  end


end

end 
