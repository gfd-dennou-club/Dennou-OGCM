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
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_historyauto
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
  public :: DataFileSet_OutputData

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

    use Constants_mod, only: RPlanet, PI

    use TemporalIntegSet_mod, only: &
         & Nl, StartTime, TotalIntegTime

    use GridSet_mod, only: &
         & iMax, jMax, kMax, &
         & x_Lon, y_Lat

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
    call HistoryAutoCreate( &                            ! ヒストリー作成
         & title='OGCM Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/'lon','lat','sig','t  '/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', 'sec.        '/), &
         & origin=StartTime, interval=this%outputIntrval, terminus=TotalIntegTime,          &
         & namelist_filename=configNmlFileName )    

    call HistoryAutoPutAxis('lon', x_Lon*180/PI)
    call HistoryAutoAddAttr('lon', 'topology', 'circular')
    call HistoryAutoAddAttr('lon', 'modulo', 360.0)
    call HistoryAutoPutAxis('lat', y_Lat*180/PI)

    call HistoryAutoAddVariable( &
         varname='u', dims=(/'lon','lat','sig','t  '/), & 
         longname='velocity(longitude) ', units='m/s')

    call HistoryAutoAddVariable( &                  
         varname='v', dims=(/'lon','lat','sig','t  '/), & 
         longname='velocity(latitude) ', units='m/s')

    call HistoryAutoAddVariable( &                  
         varname='eta', dims=(/'lon','lat', 't  '/), & 
         longname='surface height ', units='m')

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

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if( mod(CurrentTime, this%outputIntrval) /= 0 ) return 


    call MessageNotify("M", module_name, "Output data of some field at %d [sec] ..", i=(/ int(CurrentTime) /))
    call HistoryAutoPut(CurrentTime, "u", xyz_UN)
    call HistoryAutoPut(CurrentTime, "v", xyz_VN)
    call HistoryAutoPut(CurrentTime, "eta", xy_SurfHeightN)

  end subroutine DataFileSet_OutputData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

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
       call FileOpen( unit_nml, &          ! (out)
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

