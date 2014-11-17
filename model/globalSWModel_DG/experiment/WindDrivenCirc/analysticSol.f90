program analysticSol
  use dc_types
  use l_module

  implicit none

  integer, parameter :: im=1024
  integer, parameter :: jm=256
  integer, parameter :: nm = 171
  integer, parameter :: lm = (nm+1)!**2

  real(DP), dimension(0:im-1,jm), save  :: &
       & xy_Psi, xy_TauLon, xy_TauLat, &
       & xy_CosLat, xy_Lat, xy_Lon, xy_RotWindStress, &
       & xy_Ucos, xy_Vcos

  real(DP), save :: l_Psi(lm), l_RotSrc(lm)
  real(DP) :: PI = acos(-1d0)

  real(DP), parameter :: radius  = 6.37122d06
  real(DP), parameter :: Omega = 7.292d-05
  real(DP), parameter :: Grav = 9.80616
  real(DP), parameter :: meanDepth = 3000d0
  real(DP), parameter :: refDens = 1000d0
  real(DP), parameter :: LinearDragCoef = 3d-6
  real(DP), parameter :: WindStessMag = 1d-1
  
  integer :: l, m, n, i, j, k
  integer :: nm_(2)
  real(DP) :: r

  !

  call l_Initial(nm, jm)

  !
  do i=0, im-1
     xy_Lat(i,:) = y_Lat
     xy_Lon(i,:) = -PI/2d0 + PI/(im-1)*i
  end do

  xy_CosLat = cos(xy_Lat)

  r = LinearDragCoef/Omega
  xy_TauLon = -(WindStessMag/(refDens*Radius*Omega**2*meanDepth))&
       &*cos(3d0*xy_Lat)
  xy_TauLat = 0d0
  
  l_RotSrc = w_AlphaOptr_xy(xy_TauLat(0,:)*xy_CosLat(0,:), -xy_TauLon(0,:)*xy_CosLat(0,:))
  xy_RotWindStress = spread(y_l(l_RotSrc), 1, im)

  !
  call calc_streamFuncField()

  !
  call Output_data()

contains
  subroutine calc_streamFuncField()
    integer :: l, n, i, j
    real(DP) :: delM, gam, C1, C2, C0, r1, r2
    real(DP) :: lonW, lonE, lon, dlon

    real(DP) :: y_Pn(jm,lm), l_Work(lm)

    lonW = -PI/2d0; lonE = PI/2d0
    dlon = lonE - lonW

    do l=1,lm
       l_Work = 0d0; l_Work(l) = 1d0
       y_Pn(:,l) = y_l(l_Work)
    end do
    
    do j=1, jm
       delM = (LinearDragCoef/Omega)/(2d0*cos(y_Lat(j))**2 )

       xy_Psi(:,j) = 0d0
       xy_Ucos(:,j) = 0d0
       xy_Vcos(:,j) = 0d0

       do l=1, lm
          n = l-1
          gam = (LinearDragCoef/Omega)*(n*(n+1)/2d0)
          
          r1 = (- 1d0 + sqrt(1d0 + 4d0*gam*delM))/(2d0*delM)
          r2 = (- 1d0 - sqrt(1d0 + 4d0*gam*delM))/(2d0*delM)
          if(n==0) then
             C0 = 0.5d0*l_RotSrc(l)
             C1 = - C0* (lonW*exp(r2*dLon) - lonE)/(exp((r2-r1)*dlon) - 1d0)
             C2 = - C0* (lonW - lonE*exp(-r1*dLon))/(1d0 - exp((r2-r1)*dlon)) * exp(-r1*lonE)
          else
             C0 = - 0.5d0/gam*l_RotSrc(l)
             C1 = - C0* (1d0 - exp( r2*dLon))/(1d0 - exp(-(r1-r2)*dlon))
             C2 = - C0* (1d0 - exp(-r1*dLon))/(1d0 - exp(-(r1-r2)*dlon))
          end if
          
          if(n==0) then
             do i=0, im-1
                lon = xy_Lon(i,1)
                xy_Psi(i,j) = xy_Psi(i,j) + (C0*lon + C1*exp(-r1*(lonE-lon)) + C2*exp(r2*(lon-lonW)))*y_Pn(j,l)
                xy_Vcos(i,j) = xy_Vcos(i,j) + (C0 + r1*C1*exp(-r1*(lonE-lon)) + r2*C2*exp(r2*(lon-lonW)))*y_Pn(j,l)
             end do
          else
             do i=0, im-1
                lon = xy_Lon(i,1)
                xy_Psi(i,j) = xy_Psi(i,j) + (C0 + C1*exp(r1*(lon-lonE)) + C2*exp(r2*(lon-lonW)))*y_Pn(j,l)
                xy_Vcos(i,j) = xy_Vcos(i,j) + (r1*C1*exp(r1*(lon-lonE)) + r2*C2*exp(r2*(lon-lonW)))*y_Pn(j,l)
             end do
          end if
          
       end do
    end do

  end subroutine calc_streamFuncField

  function w_AlphaOptr_xy(xy_A, xy_B)
    real(DP), dimension(1:jm), intent(in)   :: xy_A
    real(DP), dimension(1:jm), intent(in)   :: xy_B
    real(DP), dimension(lm)       :: w_AlphaOptr_xy

    w_AlphaOptr_xy = (  l_DivMu_y(xy_B) )
  end function w_AlphaOptr_xy

  subroutine Output_data()
    use gtool_history

    character(TOKEN) :: longnames(3), units(3)

    longnames(1) = 'longitude'; units(1) = 'degrees_east'
    longnames(2) = 'latitude'; units(2) = 'degrees_north'
    longnames(3) = 'time'; units(3) = 'sec'
    
    call HistoryCreate( &                        ! ヒストリー作成
         & file='diffusion_1.nc', title='Diffusion equation', &
         & source='Sample program of gtool_history/gtool5',   &
         & institution='GFD_Dennou Club davis project',       &
         & dims=(/'lon', 'lat', 't  '  /), dimsizes=(/ im, jm, 0 /),               &
         & longnames=longnames,       &
         & units=units,                                 &
         & origin=real(0), interval=real(10) )

    call HistoryAddVariable( &                   ! 変数定義
         & varname='Psi', dims=(/'lon','lat'/), &
         & longname='nondimensional stream function', units='1', xtype='double')

    call HistoryAddVariable( &                   ! 変数定義
         & varname='RotWindStress', dims=(/'lon','lat'/), &
         & longname='nondimensional stream function', units='1', xtype='double')

    call HistoryAddVariable( &                   ! 変数定義
         & varname='V', dims=(/'lon','lat'/), &
         & longname='meriodional velocity', units='m/s', xtype='double')

    call HistoryPut('lon', xy_Lon(:,1)*180d0/PI)
    call HistoryPut('lat', y_Lat*180d0/PI)

    call HistoryPut('Psi', xy_Psi)
    call HistoryPut('RotWindStress', xy_RotWindStress)
    call HistoryPut('V', Radius*Omega*xy_Vcos/xy_CosLat)

    call HistoryClose()

  end subroutine Output_data

end program analysticSol
