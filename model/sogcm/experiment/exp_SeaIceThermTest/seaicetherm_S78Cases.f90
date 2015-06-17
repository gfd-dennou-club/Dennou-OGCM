!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
program SeaIceThermDyn_main

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  use dc_calendar
  use gtool_history

  use Constants_mod
  use SeaIceConstants_mod
  use SeaIceThermDyn_Winton2000_mod

  ! 宣言文; Declareration statements
  !  
  implicit none

  !**        *************

  character(*), parameter :: ConfigNmlName = 'config.nml'
  real(DP), parameter :: dt=86400d0
  real(DP), parameter :: EndTimeSec=86400d0*360d0*100d0
  real(DP), parameter :: MonitorIntrv = 86400d0*360d0
  
  real(DP), parameter :: OutputIntrv = 86400d0*3d0!* 60d0*10d0
  integer, parameter :: OutputChunkSize = 12001
  
  !　Variables read from namelist
  !

  real(DP) :: Fb
  real(DP) :: SnowFallFactor
  real(DP) :: init_hi
  character(STRING) :: ncFileName

  ! 
  !
  
  type(SeaIceThermDynVarSet), save :: varset
  real(DP) :: excessMeltEnergy
  real(DP) :: Q_SI, Q_LI, SW, LW
  
  integer, parameter :: iMax = 1, jMax=1
  
  type(DC_CAL_DATE) :: start_date, end_date
  real(DP) :: CurrentTimeSec, dhs

  type(gt_history) :: hst_seaice
  integer :: OutputChunkCounter
  real(DP), dimension(OutputChunkSize) :: &
       & chnk_hi, chnk_hs, chnk_Ts, chnk_T1, chnk_T2, &
       & chnk_SI, chnk_LI, chnk_SW, chnk_LW
  character(TOKEN) :: chnk_datetime(OutputChunkSize)

  integer :: n
  
  !
  character(*), parameter :: PROGRAM_NAME = 'seaicetherm_S78Cases'
  

  ! 実行文; Executable statements
  !
  
  write(*,*) "***** Initialization *******"  
  call initialize()
  
  !
  varset%hsN = 0.29d0
  varset%hiN = init_hi
  varset%TsN = -30d0!-29.2d0
  varset%T1N = -12d0!-13.3d0
  varset%T2N = -5d0!5.95d0

  !
  call get_Forcing(Q_SI, Q_LI, SW, LW, &
       & dt*0d0 )
  call Output_vars(dt*0d0, Q_SI, Q_LI, SW, LW)
 

  do while(CurrentTimeSec < EndTimeSec)
     call print(CurrentTimeSec)
     
     call get_SnowFallAccum(dhs, &
          & CurrentTimeSec, varset%TsN, varset%hsN > 0d0)
     varset%hsN = varset%hsN + dhs

     call get_Forcing(Q_SI, Q_LI, SW, LW, &
          & CurrentTimeSec )     

     call advance_SeaIceThermDynProc(varset, & ! (inout)
          & DelTime=dt,                                & ! (in)
          & Q_SI=Q_SI, Q_LI=Q_LI, SW=SW, LW=LW, Fb=Fb, & ! (in)
          & excessMeltEnergy=excessMeltEnergy          & ! (out)
          & )

     CurrentTimeSec = CurrentTimeSec + dt
     
     call update_SeaIceThermDynProcVar(varset)
     call Output_vars(CurrentTimeSec, Q_SI, Q_LI, SW, LW, &
          & isForcedChunkClear=(CurrentTimeSec==EndTimeSec) &
          & )
  end do
  
  !
  write(*,*) "***** Finalization *******"
  call SeaIceThermDyn_Winton2000_Final()
  call SeaIceConstants_Final()

contains

  subroutine initialize()


    
    ! 宣言文; Declareration statements
    !

    ! 作業変数
    ! Work variables
    !
    integer :: day_in_month(12)    

    ! 実行文; Executable statements
    !

    
    ! Initialize a thermodynamics component of sea ice model. 
    !
    call SeaIceConstants_Init(ConfigNmlName)
    call SeaIceThermDyn_Winton2000_Init()

    ! Read namelist to get the values of some parameters. 
    call read_nml()
    
    ! Initialize some variables to manage the passage of time.
    !
    CurrentTimeSec = 0d0
    
    day_in_month(:) = 30
    call DCCalCreate(12, day_in_month, 24, 60, 60d0)
    call DCCalDateCreate(year=2000,month=1,day=1,hour=0,min=0,sec=0d0, &
         & date=start_date)
   
   ! Preparation to output data
   !
   call prepair_Output(OriginTime=CurrentTimeSec, EndTimeSec=EndTimeSec, &
         & IntrvSec=OutputIntrv)
    OutputChunkCounter = 0
   
  end subroutine initialize

  subroutine print(CurrentSec)

    ! 宣言文; Declareration statements
    !    
    
    real(DP), intent(in) :: CurrentSec
    integer :: year, month, day
    real(DP) ::Sec

    character(STRING) :: date_str, enddate_str

    ! 実行文; Executable statements
    !
    
    if(mod(CurrentTimeSec,MonitorIntrv)/=0) return
    
    call DCCalDateInquire(date_str, &
         & elapse_sec=CurrentSec, date=start_date)

    call DCCalDateInquire(enddate_str, &
         & elapse_sec=EndTimeSec, date=start_date)

    write(*,*) CurrentSec / EndTimeSec
    write(*,*) "[", trim(date_str), "] (--[", trim(enddate_str), "] *********************"
    write(*,*) "hs=",varset%hsN,"hi=",varset%hiN, &
         & "Ts=", varset%TsN, "T1=", varset%T1N, "T2=", varset%T2N
  end subroutine print

  subroutine get_Forcing(Q_SI, Q_LI, SW, LW, &
       & CurrentSec )

    ! 宣言文; Declareration statements
    !    
    real(DP), intent(out) :: Q_SI, Q_LI, SW, LW
    real(DP), intent(in) :: CurrentSec

    ! 作業変数
    ! Work variables
    !    
    integer :: year, month, day
    real(DP) ::Sec

    real(DP), parameter :: SI_DAT(12) = &
         & (/ -1.18d0, -0.76d0, -0.72d0, -0.29d0, 0.45d0, 0.39d0, 0.30d0, 0.40d0, 0.17d0, -0.10d0, -0.56d0, -0.79d0 /) &
         & *4.184d7/(30d0*86400d0)

    real(DP), parameter :: LI_DAT(12) = &
         & (/ 0d0, 0.02d0, 0.03d0, 0.09d0, 0.46d0, 0.70d0, 0.64d0, 0.66d0, 0.39d0, 0.19d0, 0.01d0, 0.01d0 /) &
         & *4.184d7/(30d0*86400d0)

    real(DP), parameter :: SW_DAT(12) = &
         & (/ 0d0, 0d0, -1.9d0, -9.9d0, -17.7d0, -19.2d0, -13.6d0, -9.0d0, -3.7d0, -0.4d0, 0d0, 0d0 /) &
         & *4.184d7/(30d0*86400d0)

    real(DP), parameter :: LW_DAT(12) = &
         & (/ -10.4d0, -10.3d0, -10.3d0, -11.6d0, -15.1d0, -18.0d0, &
         &    -19.1d0, -18.7d0, -16.5d0, -13.9d0, -11.2d0, -10.9d0 /) &
         & *4.184d7/(30d0*86400d0)


    ! 実行文; Executable statements
    !
    
    call DCCalDateInquire(year=year,month=month,day=day, sec=sec, &
         & elapse_sec=CurrentSec, date=start_date)

    Q_SI = interpolate_DAT12(SI_DAT, month, day, sec)
    Q_LI = interpolate_DAT12(LI_DAT, month, day, sec)
    SW = interpolate_DAT12(SW_DAT, month, day, sec)
    Lw = interpolate_DAT12(LW_DAT, month, day, sec)
    
!!$    write(*,*) year, month, day, sec, ":", Q_SI, Q_LI, SW, LW

  end subroutine get_Forcing

  subroutine get_SnowFallAccum(dhs, &
       & CurrentSec, Ts, isSnowLyrExist )
    ! 宣言文; Declareration statements
    !        
    real(DP), intent(out) :: dhs
    real(DP), intent(in) :: CurrentSec, Ts
    logical, intent(in) :: isSnowLyrExist

    !
    ! 局所変数
    ! Local variables
    !    
    integer :: year, month, day
    real(DP) ::Sec
    integer :: dayOfYear

    ! 実行文; Executable statements
    !
    
    call DCCalDateInquire(year=year,month=month,day=day, sec=sec, &
         & elapse_sec=CurrentSec, date=start_date)

    dayOfYear = (month-1)*30  + day
    if((8-1)*30 + 20 <= dayOfYear .and. dayOfYear <= 10*30) then
       if((isSnowLyrExist.and.Ts<0d0).or.(.not.isSnowLyrExist .and. Ts<-Mu*SaltSeaIce))then
          dhs  = 0.3d0/(70d0*86400d0)*dt
       else
          dhs = 0d0
       end if
    else if(10*30 < dayOfYear .or. dayOfYear <= 4*30) then
       dhs  = 0.05d0/(180d0*86400d0)*dt
    else if(4*30<dayOfYear .and. dayOfYear<=5*30 .and. (isSnowLyrExist.and.Ts < 0d0))then
       dhs = 0.05/(30*86400d0)*dt
    else
       dhs  = 0.0d0
     end if

     dhs = dhs*SnowFallFactor
     
  end subroutine get_SnowFallAccum

  function interpolate_DAT12(DAT12, month, day, sec) result(val)

    ! 宣言文; Declareration statements
    !        
    real(DP), intent(in) :: DAT12(12)
    integer, intent(in) :: month, day
    real(DP), intent(in) ::sec
    real(DP) :: val

    !
    ! 局所変数
    ! Local variables
    !    
    real(DP), parameter :: SecOfDay = 86400d0
    integer, parameter :: ids(-1:14) = (/ 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2 /)
    real(DP) :: t0, t1, t2, t3
    integer :: m0
    real(DP) :: t, Work

    ! 実行文; Executable statements
    !
    
    t = (day - 1)*SecOfDay + sec !- (15d0*86400d0)

    if(t <= 15d0*SecOfDay) then
       m0 = month - 2
       t2 = 15d0*SecOfDay;
       t0 = t2 - 60d0*SecOfDay; t1 = t2 - 30d0*SecOfDay; t3 = t2 + 30d0*SecOfDay
    else
       m0 = month - 1        
       t1 = 15d0*SecOfDay;
       t0 = t1 - 30d0*SecOfDay; t2 = t1 + 30d0*SecOfDay; t3 = t1 + 60d0*SecOfDay       
    end if
    
    val = &
      &   DAT12(ids(m0))   * (t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3))     &
      & + DAT12(ids(m0+1)) * (t-t0)*(t-t2)*(t-t3)/((t1-t0)*(t1-t2)*(t1-t3))     &
      & + DAT12(ids(m0+2)) * (t-t0)*(t-t1)*(t-t3)/((t2-t0)*(t2-t1)*(t2-t3))     &
      & + DAT12(ids(m0+3)) * (t-t0)*(t-t1)*(t-t2)/((t3-t0)*(t3-t1)*(t3-t2))
    
  end function interpolate_DAT12

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine Output_vars( &
       & CurrentTimeSec, Q_SI, Q_LI, SW, LW, &
       & isForcedChunkClear )

    use dc_string, only: CPrintf

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: CurrentTimeSec
    real(DP), intent(in) :: Q_SI, Q_LI, SW, LW
    logical, optional  :: isForcedChunkClear

    !
    ! 局所変数
    ! Local variables
    !    
    character(TOKEN) :: range_str, date_str
    integer :: n, i
    real :: CurrentDay, IntrvDay
    character(TOKEN*OutputChunkSize) :: datestr1D
    
    ! 実行文; Executable statements
    !
    
    if(mod(CurrentTimeSec,OutputIntrv)/=0) return

    OutputChunkCounter = OutputChunkCounter + 1
    chnk_hs(OutputChunkCounter) = varset%hsN
    chnk_hi(OutputChunkCounter) = varset%hiN
    chnk_Ts(OutputChunkCounter) = varset%TsN
    chnk_T1(OutputChunkCounter) = varset%T1N
    chnk_T2(OutputChunkCounter) = varset%T2N
    chnk_SI(OutputChunkCounter) = Q_SI
    chnk_LI(OutputChunkCounter) = Q_LI
    chnk_SW(OutputChunkCounter) = SW
    chnk_LW(OutputChunkCounter) = LW

    call DCCalDateInquire(date_str, &
         & elapse_sec=CurrentTimeSec, date=start_date)
    chnk_datetime(OutputChunkCounter) = date_str
    
    if( present(isForcedChunkClear) )  then
       if(.not. isForcedChunkClear .and. OutputChunkCounter < OutputChunkSize) return
    elseif(OutputChunkCounter < OutputChunkSize) then
       return
    end if
    
    CurrentDay = CurrentTimeSec/86400d0
    IntrvDay = OutputIntrv/86400d0
    range_str = CPrintf("t=%r:%r", r=(/ CurrentDay-IntrvDay*(OutputChunkCounter-1), CurrentDay /))
    n = OutputChunkCounter
    call HistoryPut('hs', chnk_hs(1:n), hst_seaice, range_str)
    call HistoryPut('hi', chnk_hi(1:n), hst_seaice, range_str)
    call HistoryPut('Ts', chnk_Ts(1:n), hst_seaice, range_str)
    call HistoryPut('T1', chnk_T1(1:n), hst_seaice, range_str)
    call HistoryPut('T2', chnk_T2(1:n), hst_seaice, range_str)

    call HistoryPut('SI', chnk_SI(1:n), hst_seaice, range_str)    
    call HistoryPut('LI', chnk_LI(1:n), hst_seaice, range_str)
    call HistoryPut('SW', chnk_SW(1:n), hst_seaice, range_str)
    call HistoryPut('LW', chnk_LW(1:n), hst_seaice, range_str)

    do i=1, n
       dateStr1d(TOKEN*(i-1)+1:TOKEN*i) = trim(chnk_datetime(i))
    end do
    call HistoryPut('datetime', dateStr1d(1:TOKEN*n), hst_seaice, range_str)

    !
    OutputChunkCounter = 0
    
  end subroutine Output_vars
  
  subroutine prepair_Output(OriginTime, EndTimeSec, IntrvSec)

    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: OriginTime, EndTimeSec, IntrvSec

    !
    ! 局所変数
    ! Local variables
    !
    Character(TOKEN) :: lonName, latName, sigName, timeName, timeStrName
    real(DP) :: xy_Lon(0:0, 1), xy_Lat(0:0,1)
    real(DP), allocatable :: a_TimeDay(:)
    real(DP) :: IntrvDay
    integer :: tMax, n
    
    ! 実行文; Executable statement
    !

    !
    lonName = 'lon'; latName='lat'; timeStrName = 'timestr'; timeName='t'
    xy_Lon = 0d0; xy_Lat = 0d0

    tMax = int(EndTimeSec/IntrvSec) + 1
    allocate(a_TimeDay(tMax))
    IntrvDay = IntrvSec/86400d0    
    do n=1, tMax
       a_TimeDay(n) = (n-1)*IntrvDay
    end do
    
    call HistoryCreate( &                            ! ヒストリー作成
         & file=trim(ncFileName), title='OGCM Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/lonName,latName,timeStrName, timeName/),&
         & dimsizes=(/iMax, jMax, TOKEN, tMax/),       &
         & longnames=(/ 'longitude                        ', &
         &              'latitude                         ', &
         &              'number of characters for datetime', &
         &              'time                             '/),      &
         & units    =(/'degree_east           ', &
         &             'degree_north          ',    &
         &             '1                     ', &
         &             'days since 2000-01-01 '/),  &
         & xtypes   =(/'float ', 'float ', 'int   ', 'double'/), &
         & origin=real(OriginTime), interval=real(IntrvDay),  &        
         & history=hst_seaice  )

    
    call HistoryPut(lonName, xy_Lon(:,1)*180d0/PI, hst_seaice)
    call HistoryAddAttr(lonName, 'topology', 'circular', hst_seaice)
    call HistoryAddAttr(lonName, 'modulo', 360.0, hst_seaice)
    call HistoryPut(latName, xy_Lat(0,:)*180d0/PI, hst_seaice)
    call HistoryPut(timeName, a_TimeDay, hst_seaice)
    
    call regist_XYTVariable('Ts', 'temperature at surface', 'deg C')
    call regist_XYTVariable('hs', 'thickness of snow layer', 'm')
    call regist_XYTVariable('hi', 'thickness of ice layer', 'm')
    call regist_XYTVariable('T1', 'temperature at upper ice layer', 'deg C')
    call regist_XYTVariable('T2', 'temperature at lower ice layer', 'deg C')

    call regist_XYTVariable('SI', 'flux of sensible heat', 'W/m2')
    call regist_XYTVariable('LI', 'flux of latent heat', 'W/m2')
    call regist_XYTVariable('SW', 'shortwave radiation', 'W/m2')
    call regist_XYTVariable('LW', 'longwave radiation', 'W/m2')

    call HistoryAddVariable( &
         & 'datetime', (/ 'timestr', 't      '/), &
         & 'time represented as strings', '1', xtype='char', history=hst_seaice &
         & )

  end subroutine prepair_Output

  subroutine regist_XYTVariable(varName, long_name, units)
    character(*), intent(in) :: varName, long_name, units

    character(TOKEN) :: dims_XYT(3)

    dims_XYT = (/ 'lon', 'lat', 't  ' /)      
    call HistoryAddVariable(varName, dims_XYT, &
         & long_name, units, xtype='float', history=hst_seaice)
    
  end subroutine regist_XYTVariable

  subroutine read_nml()
    ! 引用文; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! 作業変数
    ! Work variables
    !

    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
                              ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
                              ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /SeaIceThermTest_nml/ &
         & ncFileName, Fb, SnowFallFactor, init_hi


    ! 実行文; Executable statements
    !

    
    ! Set default values
    Fb = 2d0
    ncFileName = 'Run_SeaIceThermTest.nc'
    init_hi = 2.8d0
    SnowFallFactor = 1d0
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = SeaIceThermTest_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
   end if

   call MessageNotify('M', PROGRAM_NAME, 'ncFileName = %c', c1=trim(ncFileName))
   call MessageNotify('M', PROGRAM_NAME, 'Fb              = %f  [W/m2]', d=(/ Fb /)) 
   call MessageNotify('M', PROGRAM_NAME, 'init_hi         = %f  [W/m2]', d=(/ init_hi /)) 
   call MessageNotify('M', PROGRAM_NAME, 'SnowFallFactor  = %f  [(1)]', d=(/ SnowFallFactor /))
   
  end subroutine read_nml
  
end program SeaIceThermDyn_main
