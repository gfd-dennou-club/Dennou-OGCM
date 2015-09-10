program SeaIceThermDyn_main

  ! モジュール引用
  !
  
  use dc_types, only: &
       & DP, TOKEN, STRING
  
  use dc_calendar
  use gtool_history

  use Constants_mod
  use SeaIceConstants_mod
  use SeaIceThermDyn_Winton2000_mod
  use ice_thm_mod
  
  implicit none

  !**        *************
  real(DP), parameter :: dt=86400d0
  real(DP), parameter :: EndTimeSec=86400d0*360d0*60d0  
  real(DP), parameter :: OutputIntrv = 86400d0*6d0
  character(*), parameter :: FilePrefix = 'Run2_'
  
  !
  
  type(SeaIceThermDynVarSet), save :: varset
  integer :: n
  real(DP) :: Q_SI, Q_LI, SW, LW

  integer, parameter :: iMax = 1, jMax=1
  
  type(DC_CAL_DATE) :: start_date, end_date
  real(DP) :: CurrentTimeSec, dhs

  type(gt_history) :: hst_seaice
 
  real(DP) :: LW_yr, SW_yr, SI_yr, LI_yr
  real :: hs, hi, ts, t1, t2, A, B, I, tfw, fb, tmelt, bmelt
  real :: alb, pen, trn
  real :: heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice
  
  ! 実行文
  !
  call initialize()
  
  call SeaIceConstants_Init()

  !call SeaIceThermDyn_Winton2000_Init()

  !
  hs = 0.316d0
  hi = 2.8d0
  ts = -29.2d0
  t1 = -13.3d0
  t2 = -5.95d0
  tfw = FreezeTempSeaWater
  fb = 0d0

  !
  varset%hsN = 0.316d0
  varset%hiN = 2.8d0
  varset%TsN = -29.2d0
  varset%T1N = -13.3d0
  varset%T2N = -5.95d0
  
  
  !
  LW_yr = 0d0; SW_yr = 0d0; SI_yr=0d0; LI_yr = 0d0
  call get_Forcing(Q_SI, Q_LI, SW, LW, &
       & dt*0d0 )
  call Output_vars(dt*0d0, Q_SI, Q_LI, SW, LW)

  
  do n=1, int(EndTimeSec/dt)
     call print(CurrentTimeSec)
     
     call get_SnowFallAccum(dhs, &
          & CurrentTimeSec, varset%TsN, varset%hsN > 0d0)
     varset%hsN = varset%hsN + dhs
     
     call get_Forcing(Q_SI, Q_LI, SW, LW, &
          & CurrentTimeSec )     

     call ice_optics(alb, pen, trn, &
          & hs, hi, ts, tfw )
     if(hs>0.0) then
        alb = 0.8; pen = 0.0
        if(ts==0.0) alb =0.75
     else
        pen = 0.3
        alb = 0.65
     end if
     write(*,*) "SIS alb,pen,trn:", alb,pen,trn


     B = 4.0*SBConst*(ts+273.15)**3     
     A =    Q_LI + Q_SI + LW + (1.0 - alb)*(1.0 - pen)*SW + SBConst*(ts+273.15)**4 &
          & - ts*B
     
     I = -(1.0 - alb)*pen*SW
     write(*,*) "*SIS I", I
     tmelt = 0d0; bmelt = 0d0
     call ice3lay_temp( hs,  hi, &
          & t1, t2, ts, &
          & A, B, I, tfw, fb, real(dt), &
          & tmelt, bmelt )
     write(*,*) "*SIS ts, t1, t2:", ts, t1, t2

     call ice3lay_resize(hs, hi, t1, t2, &
          & real(dhs*DensSnow), 0.0, 0.0, tmelt, bmelt, &
          & tfw, heat_to_ocn, h2o_to_ocn, h2o_from_ocn, snow_to_ice )
     write(*,*) "*SIS hs hi:",  hs, hi
     
     call advance_SeaIceThermDynProc(varset, & !(inout)
          & DelTime=dt, &
          & Q_SI=Q_SI, Q_LI=Q_LI, SW=SW, LW=LW, Fb=2d0*0d0 )


     write(*,*) "--------------------------------"
     CurrentTimeSec = CurrentTimeSec + dt
     call update_SeaIceThermDynProcVar(varset)
!!$     call Output_vars(dt*n, Q_SI, Q_LI, SW, LW)
  end do
  
  !
  !call SeaIceThermDyn_Winton2000_Final()
  call SeaIceConstants_Final()
  
contains

  subroutine initialize()

    integer :: day_in_month(12)
    
    CurrentTimeSec = 0d0

    !
    day_in_month(:) = 30
    call DCCalCreate(12, day_in_month, 24, 60, 60d0)
    call DCCalDateCreate(year=2000,month=1,day=1,hour=0,min=0,sec=0d0, &
         & date=start_date)

    !
    call prepair_Output(OriginTime=CurrentTimeSec, EndTime=EndTimeSec, &
         & Intrv=OutputIntrv, FilePrefix=FilePrefix)

  end subroutine initialize

  subroutine print(CurrentSec)

    real(DP), intent(in) :: CurrentSec
    integer :: year, month, day
    real(DP) ::Sec

    character(STRING) :: date_str
    
    call DCCalDateInquire(year=year,month=month,day=day, sec=sec, &
         & elapse_sec=CurrentSec, date=start_date)

    call DCCalDateInquire(date_str, &
         & elapse_sec=CurrentSec, date=start_date)
    write(*,*) "[", trim(date_str), "]"
    
  end subroutine print

  subroutine get_Forcing(Q_SI, Q_LI, SW, LW, &
       & CurrentSec )
    real(DP), intent(out) :: Q_SI, Q_LI, SW, LW
    real(DP), intent(in) :: CurrentSec

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
    
    
    call DCCalDateInquire(year=year,month=month,day=day, sec=sec, &
         & elapse_sec=CurrentSec, date=start_date)

    Q_SI = interpolate_DAT12(SI_DAT, month, day, sec)
    Q_LI = interpolate_DAT12(LI_DAT, month, day, sec)
    SW = interpolate_DAT12(SW_DAT, month, day, sec)
    Lw = interpolate_DAT12(LW_DAT, month, day, sec)

    LW_yr = LW_yr + LW*dt
    SW_yr = SW_yr + SW*dt
    SI_yr = SI_yr + Q_SI*dt
    LI_yr = LI_yr + Q_LI*dt
    
    write(*,*) year, month, day, sec, ":", &
         & Q_SI, Q_LI, SW, LW
    write(*,*) LW_yr/4.184d7, SW_yr/4.184d7, LI_yr/4.184d7, SI_yr/4.184d7
  end subroutine get_Forcing

  subroutine get_SnowFallAccum(dhs, &
       & CurrentSec, Ts, isSnowLyrExist )
    real(DP), intent(out) :: dhs
    real(DP), intent(in) :: CurrentSec, Ts
    logical, intent(in) :: isSnowLyrExist
    
    integer :: year, month, day
    real(DP) ::Sec
    integer :: dayOfYear


    call DCCalDateInquire(year=year,month=month,day=day, sec=sec, &
         & elapse_sec=CurrentSec, date=start_date)

    dayOfYear = (month-1)*30  + day
    if((8-1)*30 + 20 <= dayOfYear .and. dayOfYear <= 10*30) then
!       if((isSnowLyrExist.and.Ts<0d0).or.(.not.isSnowLyrExist .and. Ts<-Mu*SaltSeaIce))then
          dhs  = 0.3d0/(70d0*86400d0)*dt
!       else
!          dhs = 0d0
!       end if
    else if(10*30 < dayOfYear .or. dayOfYear <= 4*30) then
       dhs  = 0.05d0/(180d0*86400d0)*dt
    else if(4*30<dayOfYear .and. dayOfYear<=5*30 .and. (isSnowLyrExist.and.Ts < 0d0))then
       dhs = 0.05/(30*86400d0)*dt
    else
       dhs  = 0.0d0
     end if
    
  end subroutine get_SnowFallAccum

  function interpolate_DAT12(DAT12, month, day, sec) result(val)
    real(DP), intent(in) :: DAT12(12)
    integer, intent(in) :: month, day
    real(DP), intent(in) ::sec
    real(DP) :: val

    real(DP), parameter :: SecOfDay = 86400d0
    integer, parameter :: ids(-1:14) = (/ 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2 /)
    real(DP) :: t0, t1, t2, t3
    integer :: m0
    real(DP) :: t, Work

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
       &   CurrentTimeSec, Q_SI, Q_LI, SW, LW )

    real(DP), intent(in) :: CurrentTimeSec
    real(DP), intent(in) :: Q_SI, Q_LI, SW, LW
    
    if(mod(int(CurrentTimeSec),int(OutputIntrv))/=0) return

    call HistoryPut('hs', varset%hsN, hst_seaice)
    call HistoryPut('hi', varset%hiN, hst_seaice)
    call HistoryPut('Ts', varset%TsN, hst_seaice)
    call HistoryPut('T1', varset%T1N, hst_seaice)
    call HistoryPut('T2', varset%T2N, hst_seaice)

    call HistoryPut('SI', Q_SI, hst_seaice)    
    call HistoryPut('LI', Q_LI, hst_seaice)
    call HistoryPut('SW', SW, hst_seaice)
    call HistoryPut('LW', LW, hst_seaice)
    
    
  end subroutine Output_vars
  
  subroutine prepair_Output(OriginTime, EndTime, Intrv, FilePrefix)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: OriginTime, EndTime, Intrv
    character(*), intent(in) :: FilePrefix

    !
    ! 局所変数
    ! Local variables
    !
    Character(TOKEN) :: lonName, latName, sigName, timeName
    real(DP) :: xy_Lon(0:0, 1), xy_Lat(0:0,1)

    ! 実行文; Executable statement
    !

    !
    lonName = 'lon'; latName='lat'; timeName='t'
    xy_Lon = 0d0; xy_Lat = 0d0
    
    call HistoryCreate( &                            ! ヒストリー作成
         & file=trim(FilePrefix) // 'SeaiceOutput.nc', title='OGCM Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/lonName,latName,timeName/), dimsizes=(/iMax,jMax,0/),       &
         & longnames=(/'longitude','latitude ', 'time     '/),      &
         & units=(/'degree_east ','degree_north', 'sec         '/),  &
         & origin=real(OriginTime), interval=real(Intrv),  &        
         & history=hst_seaice  )

    
    call HistoryPut(lonName, xy_Lon(:,1)*180d0/PI, hst_seaice)
    call HistoryAddAttr(lonName, 'topology', 'circular', hst_seaice)
    call HistoryAddAttr(lonName, 'modulo', 360.0, hst_seaice)
    call HistoryPut(latName, xy_Lat(0,:)*180d0/PI, hst_seaice)

    call regist_XYTVariable('Ts', 'temperature at surface', 'deg C')
    call regist_XYTVariable('hs', 'thickness of snow layer', 'm')
    call regist_XYTVariable('hi', 'thickness of ice layer', 'm')
    call regist_XYTVariable('T1', 'temperature at upper ice layer', 'deg C')
    call regist_XYTVariable('T2', 'temperature at lower ice layer', 'deg C')

    call regist_XYTVariable('SI', 'flux of sensible heat', 'W/m2')
    call regist_XYTVariable('LI', 'flux of latent heat', 'W/m2')
    call regist_XYTVariable('SW', 'shortwave radiation', 'W/m2')
    call regist_XYTVariable('LW', 'longwave radiation', 'W/m2')
    
  end subroutine prepair_Output

  subroutine regist_XYTVariable(varName, long_name, units)
    character(*), intent(in) :: varName, long_name, units

    character(TOKEN) :: dims_XYT(3)

    dims_XYT = (/ 'lon', 'lat', 't  ' /)      
    call HistoryAddVariable(varName, dims_XYT, &
         & long_name, units, xtype='float', history=hst_seaice)
    
  end subroutine regist_XYTVariable
  
end program SeaIceThermDyn_main
