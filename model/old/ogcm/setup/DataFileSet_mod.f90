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

  use PolyMesh_mod, only: &
       & PolyMesh

  use netcdfDataWriter_mod, only: &
       & netcdfDataWriter, &
       & netcdfDataWriter_Init, netcdfDataWriter_Final, &
       & netcdfDataWriter_writeGlobalAttr, netcdfDataWriter_Regist, &
       & netcdfDataWriter_write, &
       & netcdfDataWriter_AdvanceTimeStep


  use VariableSet_mod, only: &
       & VariableSet

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  type, public :: DataFileSet
     
     type(netcdfDataWriter) :: dataWriter  
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
  subroutine DataFileSet_Init(this, configNmlFileName, plMesh, var)

    ! モジュール引用; Use statement
    !

    use Constants_mod, only: RPlanet

    use TemporalIntegSet_mod, only: &
         & Nl

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    character(*), intent(in) :: configNmlFileName
    type(PolyMesh), intent(in) :: plMesh
    type(VariableSet), intent(in) :: var

    ! 局所変数
    ! Local variable
    !
    character(STRING) :: outputFileName


    ! 実行文; Executable statements
    !

    !
    call read_nmlData(this, configNmlFileName, outputFileName)
    
    !
    call netcdfDataWriter_Init(this%dataWriter, outputFileName, plMesh)
    
    call netcdfDataWriter_writeGlobalAttr(this%dataWriter, "RPlanet", RPlanet)
    call netcdfDataWriter_Regist(this%dataWriter, (/ var%zc_lyrThick(Nl) /))
    
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

    call netcdfDataWriter_Final(this%dataWriter)

  end subroutine DataFileSet_Final

  !> @brief 
  !!
  !!
  subroutine DataFileSet_OutputData(this, var)

    ! モジュール引用; Use statement
    !
    use TemporalIntegSet_mod, only: &
         & CurrentTime, Nl

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(inout) :: this
    type(VariableSet), intent(in) :: var
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if( mod(CurrentTime, this%outputIntrval) /= 0 ) return 


    call MessageNotify("M", module_name, "Output data of some field at %d [sec] ..", i=(/ int(CurrentTime) /))
    call netcdfDataWriter_write(this%dataWriter, var%zc_lyrThick(Nl))
    
    !
    call netcdfDataWriter_AdvanceTimeStep(this%dataWriter, CurrentTime)

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

