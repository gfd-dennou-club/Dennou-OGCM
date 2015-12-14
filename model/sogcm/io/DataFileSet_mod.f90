!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DataFileSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_historyauto

  use Constants_mod

  use VariableSet_mod

  use VarSetSeaice_mod

  use BoundaryCondO_mod, only: &
       & VARSET_KEY_WINDSTRESSLON, VARSET_KEY_WINDSTRESSLAT, &
       & xy_WindStressU, xy_WindStressV
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  type, public :: DataFileSet

     character(String) :: FilePrefix
     real(DP) :: outputIntrvalSec

  end type DataFileSet

  public :: DataFileSet_Init, DataFileSet_Final
  public :: DataFileSet_OutputData, DataFileSet_OutputBasicData
  public :: DataFileSet_isOutputTiming

  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DataFileSet_mod' !< Module Name
  real(DP) :: outputIntTimeSec
  character(TOKEN) :: outputIntUnit
  
contains

  !>
  !!
  !!
  subroutine DataFileSet_Init(this, configNmlFileName)

    ! モジュール引用; Use statement
    !
    use TemporalIntegSet_mod, only: &
         & Nl, RestartTime, IntegTime, EndTime

    use GridSet_mod, only: &
         & iMax, jMax, kMax

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variable
    !
    character(STRING) :: outputFileName

    ! 実行文; Executable statements
    !

    !
    call read_nmlData(this, configNmlFileName, outputFileName)
    
    !
    !

    !
    call HistoryAutoCreate( &                            ! ヒストリー作成
         & title='OGCM Output',             &
         & source='OGCM Output',                                        &
         & institution='GFD_Dennou Club OGCM project',                  &
         & dims=(/'lon ', 'lat ', 'sig ', 'sig2', 'time'/),             &
         & dimsizes=(/iMax, jMax, kMax+1, 2, 0 /),                      &
         & longnames=(/'longitude   ', 'latitude    ', 'sigma       ',  &
         &             'sigma-seaice', 'time        '   /),             &
         & units=(/'degree_east ', 'degree_north', '1              ',   &
         &         '1           ', 'sec.        ' /),                   &
         & origin=RestartTime, interval=this%outputIntrvalSec, terminus=EndTime,          &
         & namelist_filename=configNmlFileName )    

    ! Regist the axises and variables which will be output. 
    call regist_OutputAxisAndVar()
 
  end subroutine DataFileSet_Init

  !>
  !!
  !!
  subroutine DataFileSet_Final(this)


    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this

    ! 実行文; Executable statements
    !

    call HistoryAutoClose()

  end subroutine DataFileSet_Final

  logical function DataFileSet_isOutputTiming(this, CurrentTimeSec)

    type(DataFileSet), intent(in) :: this
    real(DP), intent(in) :: CurrentTimeSec

    DataFileSet_isOutputTiming &
         & = ( mod(CurrentTimeSec, this%outputIntrvalSec) == 0 )

  end function DataFileSet_isOutputTiming

  !> @brief 
  !!
  !!
  subroutine DataFileSet_OutputData(this)

    ! モジュール引用; Use statement
    !
    use dc_calendar, only: &
         & DCCalConvertByUnit

    use TemporalIntegSet_mod, only: &
         & CurrentTime, Nl

    use GridSet_mod, only: &
         & iMax, jMax, kMax, lMax, &
         & xyz_Lat

    use SpmlUtil_mod

    use DiagnoseUtil_mod

    use EOSDriver_mod, only: &
         & EOSDriver_Eval

    
    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    
    ! 局所変数
    ! Local variables
    !
    
    real(DP) :: wz_Vor(lMax, 0:kMax)
    real(DP) :: wz_Div(lMax, 0:kMax)
    real(DP) :: xyz_Psi(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_Chi(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_CosLat(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xy_totDepth(0:iMax-1, jMax)    
    real(DP) :: xyz_GeoPot(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_HydroPressEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)

    ! 実行文; Executable statement
    !

    if( .not. DataFileSet_isOutputTiming(this, CurrentTime) ) return 

    call MessageNotify("M", module_name, "Output data of some field at %f [%c] ..", &
         & d=(/ DCCalConvertByUnit(CurrentTime, 'sec', outputIntUnit) /), c1=trim(outputIntUnit) )


    xy_totDepth = xy_totDepthBasic + xy_SurfHeightN

    ! Output variables in OGCM
    !
    
    call HistoryAutoPut(CurrentTime, VARSET_KEY_U, xyz_UN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_V, xyz_VN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFHEIGHT, xy_SurfHeightN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_PTEMPEDD, xyz_PTempEddN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SALT, xyz_SaltN)

    call HistoryAutoPut(CurrentTime, VARSET_KEY_UB, xyz_UB)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_VB, xyz_VB)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFHEIGHTB, xy_SurfHeightB)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_PTEMPEDDB, xyz_PTempEddB)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SALTB, xyz_SaltB)
    
    xyz_CosLat = cos(xyz_Lat)
    wz_Vor = wz_AlphaOptr_xyz(xyz_VN*xyz_CosLat, -xyz_UN*xyz_CosLat) 
    wz_Div = wz_AlphaOptr_xyz(xyz_UN*xyz_CosLat,  xyz_VN*xyz_CosLat) 
    xyz_Psi = xyz_wz( wz_InvLapla2D_wz( wz_Vor ) )
    xyz_Chi = xyz_wz( wz_InvLapla2D_wz( wz_Div ) )
    xyz_SigDot = Diagnose_SigDot( xy_totDepth, xyz_UN*xyz_CosLat, xyz_VN*xyz_CosLat, xyz_wz(wz_Div) )

    xyz_PTemp = xyz_PTempEddN + spread(spread(z_PTempBasic,1,jMax), 1, iMax)
    xyz_GeoPot = Diagnose_GeoPot( xy_totDepth ) 
    call EOSDriver_Eval( rhoEdd=xyz_DensEdd,                      & ! (out)
         & theta=xyz_PTemp, S=xyz_SaltN, p=-RefDens*xyz_GeoPot )     ! (in)

    xyz_HydroPressEdd = Diagnose_HydroPressEdd(xy_totDepth, xyz_DensEdd)

    !
    !
    call HistoryAutoPut(CurrentTime, "Psi", xyz_Psi)
    call HistoryAutoPut(CurrentTime, "Chi", xyz_Chi)
    call HistoryAutoPut(CurrentTime, "Div", xyz_wz(wz_Div))
    call HistoryAutoPut(CurrentTime, "Vor", xyz_wz(wz_Vor))
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFPRESS, xy_SurfPressN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SIGDOT, xyz_SigDot)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_HYDROPRESSEDD, xyz_HydroPressEdd)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_WINDSTRESSLON, xy_WindStressU)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_WINDSTRESSLAT, xy_WindStressV)

    call HistoryAutoPut(CurrentTime, VARSET_KEY_CONVINDEX, xyz_ConvIndex)

    call HistoryAutoPut(CurrentTime, VARSET_KEY_VVISCCOEF, xyz_VViscCoefN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_VDIFFCOEF, xyz_VDiffCoefN)
    
   
    ! Output variables in sea-ice model
    !
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SICECON, xy_SIceConN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_ICETHICK, xy_IceThickN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SNOWTHICK, xy_SnowThickN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SICETEMP, xyz_SIceTempN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SICESURFTEMP, xy_SIceSurfTempN)
    
  end subroutine DataFileSet_OutputData

  !> @brief 
  !!
  !!
  subroutine DataFileSet_OutputBasicData()

    ! モジュール引用; Use statement
    !
    use GridSet_mod, only: &
         & GRIDSET_KEY_LYRTHICKSIG, &
         z_LyrThickSig
    
    use TemporalIntegSet_mod, only: &
         & CurrentTime
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    call HistoryAutoPut(CurrentTime, GRIDSET_KEY_LYRTHICKSIG, z_LyrThickSig )
    call HistoryAutoPut(CurrentTime, VARSET_KEY_TOTDEPTHBASIC, xy_totDepthBasic )
    call HistoryAutoPut(CurrentTime, VARSET_KEY_PTEMPBASIC, z_PTempBasic )

  end subroutine DataFileSet_OutputBasicData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  !> @brief 
  !!
  !!
  subroutine regist_OutputAxisAndVar()

    ! モジュール引用;
    ! Use statements
    use Constants_mod, only: PI

    use GridSet_mod, only: &
         & iMax, jMax, kMax, &
         & xyz_Lon, xyz_Lat, z_Sig,                                      &
         & x_Lon_Weight, y_Lat_Weight, z_Sig_Weight,                     &
         & lonName => GRIDSET_KEY_XAXIS, latName  => GRIDSET_KEY_YAXIS,  &
         & sigName => GRIDSET_KEY_ZAXIS, sig2Name => GRIDSET_KEY_ZAXIS2, &
         & timeName => GRIDSET_KEY_TAXIS, &
         & GRIDSET_KEY_LYRTHICKSIG

    
    ! 宣言文; Declaration statement
    !

    
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: dims_Z(1), dims_XY(2), dims_ZT(2), dims_XYT(3), dims_XYZT(4), dims_XYZ2T(4)
    
    ! 実行文; Executable statement
    !

    ! Set arrays storing the name of axises
    !
    dims_Z = (/ sigName /)
    dims_XY = (/ lonName, latName /)
    dims_ZT = (/ sigName, timeName /)
    dims_XYT = (/ lonName, latName, timeName /)
    dims_XYZT = (/ lonName, latName, sigName, timeName /)
    dims_XYZ2T = (/ lonName, latName, sig2Name, timeName /)
    
    ! Regist coordinates
    !
    call HistoryAutoAddAttr(lonName, 'standard_name', 'longitude')
    call HistoryAutoAddAttr(latName, 'standard_name', 'latitude')
    call HistoryAutoAddAttr(sigName, 'standard_name', 'ocean_sigma_coordinate')
    call HistoryAutoAddAttr(sig2Name, 'standard_name', 'seaice_sigma_coordinate')

    call HistoryAutoPutAxis(lonName, xyz_Lon(:,1,0)*180/PI)
    call HistoryAutoAddAttr(lonName, 'topology', 'circular')
    call HistoryAutoAddAttr(lonName, 'modulo', 360.0)
    call HistoryAutoPutAxis(latName, xyz_Lat(0,:,0)*180/PI)
    call HistoryAutoPutAxis(sigName, z_Sig)
    call HistoryAutoPutAxis(sig2Name, (/ -0.25d0, -0.75d0 /))

    call HistoryAutoAddWeight(lonName, x_Lon_Weight, 'radian', xtype='double')
    call HistoryAutoAddWeight(latName, y_Lat_Weight, 'radian', xtype='double')
    call HistoryAutoAddWeight(sigName, z_Sig_Weight, '1',      xtype='double')
    
    ! Regist prognostic variables
    !

    ! Regist variables in OGCM
    
    call HistoryAutoAddVariable( varname=VARSET_KEY_U, &
         & dims=dims_XYZT, longname='velocity(longitude) ', units='m/s')

    call HistoryAutoAddVariable( varname=VARSET_KEY_UB, &
         & dims=dims_XYZT, longname='velocity(longitude) ', units='m/s')

    call HistoryAutoAddVariable( varname=VARSET_KEY_V, &
         & dims=dims_XYZT, longname='velocity(latitude) ', units='m/s')
    call HistoryAutoAddVariable( varname=VARSET_KEY_VB, &
         & dims=dims_XYZT, longname='velocity(latitude) ', units='m/s')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SURFHEIGHT, &
         & dims=dims_XYT, longname='surface height ', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SURFHEIGHTB, &
         & dims=dims_XYT, longname='surface height ', units='m')
    
    call HistoryAutoAddVariable( varname=VARSET_KEY_PTEMPEDD, &
         & dims=dims_XYZT, longname='eddy component of potential temperature ', units='K')

    call HistoryAutoAddVariable( varname=VARSET_KEY_PTEMPEDDB, &
         & dims=dims_XYZT, longname='eddy component of potential temperature ', units='K')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SALT, &
         & dims=dims_XYZT, longname='salinity', units='psu')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SALTB, &
         & dims=dims_XYZT, longname='salinity', units='psu')

    call HistoryAutoAddVariable( varname=VARSET_KEY_VVISCCOEF, &
         & dims=dims_XYZT, longname='Vertical eddy viscosity', units='m2/s')

    call HistoryAutoAddVariable( varname=VARSET_KEY_VDIFFCOEF, &
         & dims=dims_XYZT, longname='Vertical eddy diffusivity', units='m2/s')

    
    !
    ! Regist variables in seaice model
    !
    
    call HistoryAutoAddVariable( varname=VARSET_KEY_SICECON, &
         & dims=dims_XYT, longname='seaice concentration ', units='1')

    call HistoryAutoAddVariable( varname=VARSET_KEY_ICETHICK, &
         & dims=dims_XYT, longname='effective sea-ice thickness ', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_ICETHICKB, &
         & dims=dims_XYT, longname='effective sea-ice thickness ', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SNOWTHICK, &
         & dims=dims_XYT, longname='effective snow depth ', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SNOWTHICKB, &
         & dims=dims_XYT, longname='effective snow depth ', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SICETEMP, &
         & dims=dims_XYZ2T, longname='sea-ice temperature ', units='K')  

    call HistoryAutoAddVariable( varname=VARSET_KEY_SICETEMPB, &
         & dims=dims_XYZ2T, longname='sea-ice temperature ', units='K')  
    
    ! 
    ! For variables in ocean model
    ! Regist diagnostic variables
    !

    call HistoryAutoAddVariable( varname=VARSET_KEY_SIGDOT, &
         & dims=dims_XYZT, longname='vertical velocity in Sigma coordinate ', units='s-1')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SURFPRESS, &
         & dims=dims_XYT, longname='surface(barotropic) pressure ', units='Pa')

    call HistoryAutoAddVariable( varname=VARSET_KEY_HYDROPRESSEDD, &
         & dims=dims_XYZT, longname='deviation of hydrostatic pressure ', units='Pa')

    ! 
    ! For variables in sea-ice model
    ! Regist diagnostic variables
    !
    call HistoryAutoAddVariable( varName=VARSET_KEY_SICESURFTEMP, &
         & dims=dims_XYT, longname='Surface temperature of snow or ice layer', units='degC')
    
    ! Regist accessory variables
    !
    call HistoryAutoAddVariable( varname=GRIDSET_KEY_LYRTHICKSIG, &
         & dims=dims_Z, longname='nomdimensional layer thickness coressponding to each vertiocal grid point.', units='1')

    call HistoryAutoAddVariable( varname=VARSET_KEY_TOTDEPTHBASIC, &
         & dims=dims_XY, longname='basic state of total depth', units='m')

    call HistoryAutoAddVariable( varname=VARSET_KEY_PTEMPBASIC, &
         & dims=dims_Z, longname='basic state of potential temperature', units='K')

    call HistoryAutoAddVariable( varname='Chi', &
         & dims=dims_XYZT, longname='velocity potential ', units='m2/s')

    call HistoryAutoAddVariable( varname='Div', &
         & dims=dims_XYZT, longname='divergence ', units='s-1')

    call HistoryAutoAddVariable( varname='Psi', &
         & dims=dims_XYZT, longname='stream function', units='m2/s')

    call HistoryAutoAddVariable( varname='Vor', &
         & dims=dims_XYZT, longname='vorcity', units='s-1')


    call HistoryAutoAddVariable( varname=VARSET_KEY_WINDSTRESSLAT, &
         & dims=dims_XYT, longname='wind stress(latitude)', units='kg.m-1.s-2')

    call HistoryAutoAddVariable( varname=VARSET_KEY_WINDSTRESSLON, &
         & dims=dims_XYT, longname='wind stress(longitude)', units='kg.m-1.s-2')

    !
    call HistoryAutoAddVariable( varname=VARSET_KEY_CONVINDEX, &
         & dims=dims_XYZT, &
         & longname='convective index(The number of calling a routine for convective adjustment per time step)', &
         & units='times per time step')

  end subroutine regist_OutputAxisAndVar


  subroutine read_nmlData( this, configNmlFileName, &
       & outputFileName )

    ! モジュール引用; Use statement
    !
    use dc_calendar, only: &
         & DCCalConvertByUnit

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    character(*), intent(in) :: configNmlFileName
    character(STRING), intent(out) :: outputFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
                              ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
                              ! IOSTAT of NAMELIST read

    character(TOKEN) :: pos_nml

    real(DP) :: IntValue
    character(TOKEN) :: IntUnit
    character(STRING) :: Name, FilePrefix

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    IntValue = 0d0
    IntUnit = "sec"

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       pos_nml = ''; iostat_nml = 0
       do while ( trim(pos_nml) /= 'APPEND' .and. iostat_nml == 0 ) 
          read( unit_nml, &           ! (in)
               & nml = gtool_historyauto_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, &      !(in)
               & position=pos_nml ) !(out)
       end do
       close( unit_nml )
    end if

    outputIntTimeSec = DCCalConvertByUnit(IntValue, IntUnit, "sec")
    outputIntUnit = IntUnit
    this%FilePrefix = FilePrefix
    this%outputIntrvalSec = outputIntTimeSec

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, 'Output Interval %f [%c]', d=(/IntValue/), c1=trim(IntUnit) )

  end subroutine read_nmlData

end module DataFileSet_mod

