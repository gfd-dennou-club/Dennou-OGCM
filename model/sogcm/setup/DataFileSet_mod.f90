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

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  type, public :: DataFileSet
     
     real(DP) :: outputIntrval

  end type DataFileSet

  public :: DataFileSet_Init, DataFileSet_Final
  public :: DataFileSet_OutputData, DataFileSet_OutputBasicData


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DataFileSet_mod' !< Module Name
  

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
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/'lon','lat','sig','t  '/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', 'sec.        '/), &
         & origin=RestartTime, interval=this%outputIntrval, terminus=EndTime,          &
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

  !> @brief 
  !!
  !!
  subroutine DataFileSet_OutputData(this)

    ! モジュール引用; Use statement
    !
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
    real(DP) :: xyz_PressBaroc(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1, jMax, 0:kMax)
    real(DP) :: xyz_PTemp(0:iMax-1, jMax, 0:kMax)

    ! 実行文; Executable statement
    !

    if( mod(CurrentTime, this%outputIntrval) /= 0 ) return 


    xy_totDepth = xy_totDepthBasic + xy_SurfHeightN

    call MessageNotify("M", module_name, "Output data of some field at %d [sec] ..", i=(/ int(CurrentTime) /))
    call HistoryAutoPut(CurrentTime, VARSET_KEY_U, xyz_UN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_V, xyz_VN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFHEIGHT, xy_SurfHeightN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_PTEMPEDD, xyz_PTempEddN)

    call HistoryAutoPut(CurrentTime, VARSET_KEY_UB, xyz_UB)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_VB, xyz_VB)
!    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFHEIGHT, xy_SurfHeightN)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_PTEMPEDDB, xyz_PTempEddB)
    
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

    xyz_PressBaroc = Diagnose_PressBaroc(xy_totDepth, xyz_DensEdd)

    call HistoryAutoPut(CurrentTime, "Psi", xyz_Psi)
    call HistoryAutoPut(CurrentTime, "Chi", xyz_Chi)
    call HistoryAutoPut(CurrentTime, "Div", xyz_wz(wz_Div))
    call HistoryAutoPut(CurrentTime, "Vor", xyz_wz(wz_Vor))
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFPRESS, xy_SurfPress)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SIGDOT, xyz_SigDot)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_BAROCPRESS, xyz_PressBaroc)

  end subroutine DataFileSet_OutputData

  !> @brief 
  !!
  !!
  subroutine DataFileSet_OutputBasicData()

    ! モジュール引用; Use statement
    !
    use TemporalIntegSet_mod, only: &
         & CurrentTime
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
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
         & x_Lon, y_Lat

    use SpmlUtil_mod, only: g_Sig
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: lonName
    character(TOKEN) :: latName
    character(TOKEN) :: sigName
    character(TOKEN) :: timeName

    character(TOKEN) :: dims_Z(1), dims_XY(2), dims_ZT(2), dims_XYT(3), dims_XYZT(4)
    
    ! 実行文; Executable statement
    !
 
    ! Regist coordinates
    !
    lonName = 'lon'; latName='lat'; sigName='sig'; timeName='t'

    call HistoryAutoPutAxis(lonName, x_Lon*180/PI)
    call HistoryAutoAddAttr(lonName, 'topology', 'circular')
    call HistoryAutoAddAttr(lonName, 'modulo', 360.0)
    call HistoryAutoPutAxis(latName, y_Lat*180/PI)
    call HistoryAutoPutAxis(sigName, g_Sig)

    !
    !
    dims_Z = (/ sigName /)
    dims_XY = (/ lonName, latName /)
    dims_ZT = (/ sigName, timeName /)
    dims_XYT = (/ lonName, latName, timeName /)
    dims_XYZT = (/ lonName, latName, sigName, timeName /)

    ! Regist prognostic variables
    !

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

    call HistoryAutoAddVariable( varname=VARSET_KEY_PTEMPEDD, &
         & dims=dims_XYZT, longname='eddy component of potential temperature ', units='K')

    call HistoryAutoAddVariable( varname=VARSET_KEY_PTEMPEDDB, &
         & dims=dims_XYZT, longname='eddy component of potential temperature ', units='K')

    ! Regist diagnostic variables
    !

    call HistoryAutoAddVariable( varname=VARSET_KEY_SIGDOT, &
         & dims=dims_XYZT, longname='vertical velocity in Sigma coordinate ', units='s-1')

    call HistoryAutoAddVariable( varname=VARSET_KEY_SURFPRESS, &
         & dims=dims_XYT, longname='surface(barotropic) pressure ', units='Pa')

    call HistoryAutoAddVariable( varname=VARSET_KEY_BAROCPRESS, &
         & dims=dims_XYZT, longname='baroclinic pressure ', units='Pa')

    ! Regist accessory variables
    !

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
   
  end subroutine regist_OutputAxisAndVar


  subroutine read_nmlData( this, configNmlFileName, &
       & outputFileName )

    ! モジュール引用; Use statement
    !

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

    real(DP) :: outputIntrval

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /datafile_nml/ &
         & outputFileName, outputIntrval

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    outputFileName = "data.nc"
    outputIntrval  = -1d0

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = datafile_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    this%outputIntrval = outputIntrval

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  outputFileName        = %a', ca=(/ outputFileName /))
    call MessageNotify( 'M', module_name, '  outputIntrval        = %d [sec]', i=(/ int(this%outputIntrval) /))

  end subroutine read_nmlData

end module DataFileSet_mod

