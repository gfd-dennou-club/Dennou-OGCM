module CollectImgLib

require "fileutils"

class Figure
  attr_accessor :name, :ncFilePath, :varName, :cutpos, :gpviewOpt, :convertOpt, :figPicExt

  def initialize(figname, ncFilePath, varName, cutpos, intrv="", range="", gpOpts="")
    @name = figname
    @ncFilePath = ncFilePath
    @varName = varName
    @cutpos = cutpos

    @gpviewOpt = "--wsn 4 -sw:lwnd=f"
    @gpviewOpt << " --range #{range}" if range.length != 0
    @gpviewOpt << " --int #{intrv}" if intrv.length != 0
    @gpviewOpt << " #{gpOpts}" if gpOpts.length != 0

    @convertOpt = @@DEFAULT_CONVERT_OPT
  end

  def createFigure(dirPath)
  end

  @@DEFAULT_CONVERT_OPT = "-trim -delay 20"
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
      
      endTimeNum = Math.log10(@animTimeInfo.end).to_i
      timePadd = format("%0#{endTimeNum+1}d", "#{time}".to_i)
      `mv dcl_001.png #{@name}_#{timePadd}.png`
#p      "mv dcl.ps #{@name}_#{timePadd}.ps"

    end

    animFigPath = "#{dirPath}#{@name}_#{animTimeInfo.init}-#{@animTimeInfo.end}.#{@figPicExt}"
    p "convert.. png => gif animation (options #{@convertOpt})"
    `convert #{@convertOpt} #{@name}_*.png #{animFigPath}`
    `rm #{@name}_*.png` 
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
    `convert #{@convertOpt} dcl_001.png #{dirPath}#{@name}.#{figPicExt}`
    `rm dcl_001.png`
  end 
end

class NoAnimOverplotFig < Figure
  attr_accessor :ncFilePaths, :varNames
  def initialize(figname, ncFilePaths, varNames, cutpos, intrv="", range="", gpOpts="")
    super(figname, "", "", cutpos, intrv, range, gpOpts)
    @ncFilePaths = ncFilePaths
    @varNames = varNames
    @figPicExt = "jpg"
  end

  def createFigure(dirPath)

    ncFilePaths_ary = @ncFilePaths.split(",")
    
    targets = "--overplot #{ncFilePaths_ary.length} "
    @varNames.split(",").each_with_index{|varName,i|
      targets << " #{ncFilePaths_ary[i]}@#{varName}"
      targets << ",#{@cutpos}" if @cutpos.length > 0
    }

    p "command: gpview #{targets} #{@gpviewOpt}"
    `gpview #{targets} #{@gpviewOpt}`
    p "convert.. ps => #{dirPath}#{@name}.#{figPicExt} (options #{@convertOpt})"    
    `convert #{@convertOpt} dcl_001.png #{dirPath}#{@name}.#{figPicExt}`
    `rm dcl_001.png`
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


  def add_NoAnimFig(varName, cutPos, figIntrv, figRange, prefix="", suffix="", gpOpts="", figPicExt="", ncFileSuffix="" )

    if ncFileSuffix.length > 0 then
      figName, ncFilePath = create_FigMetaInfo(ncFileSuffix, prefix, suffix)
    else
      figName, ncFilePath = create_FigMetaInfo(varName, prefix, suffix)
    end

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
      p "create thumbnail.. target files: *.#{targetFileExt}"

      FileUtils.mkdir(["figdir", "thum-src"])
      FileUtils.cp(dcmodel_thumb_rb, "thum-src/dcmodel-thum.rb")
      FileUtils.cp(Dir.glob("*.#{targetFileExt}"), "figdir")
      Dir.chdir("./figdir/"){
        FileUtils.rm_r(Dir.glob("*_thumb.png"))
      }
      if ( Dir.glob("./figdir/*.#{targetFileExt}").count > 0 ) then
        Dir.chdir("./thum-src/"){
          `#{rb18cmd} ./dcmodel-thum.rb`
          `#{rb18cmd} ./dcmodel-thum-make.rb`
        }
        FileUtils.mv(Dir.glob("./thumbdir/*_thumb.png"), ".")
      end
      FileUtils.rm_r(["thumbdir", "figdir", "thum-src"])
      FileUtils.rm(Dir.glob("sample_thum.htm*"))
      FileUtils.rm("thumbdir.SIGEN")
    }
  end


end

end 
