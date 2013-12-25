program EOS_mod_test
  use dc_test
  use dc_message
  use dc_string
  use dc_types
  use gtool_history
  use EqState_JM95_mod

  implicit none

  character(*), parameter :: configNmlFile = "defaultConfig.nml"
  integer, parameter :: Npress = 51
  integer, parameter :: NSal = 43 
  integer, parameter :: Ntheta = 43
  real(DP), parameter :: pMin = 0d0, pMax = 1000d0
  real(DP), parameter :: thetaMin = -2d0, thetaMax = 40d0
  real(DP), parameter :: SMin = 0d0, SMax = 42d0
  
  real(DP) :: ptemp(Ntheta)
  real(DP) :: press(Npress)
  real(DP) :: sal(NSal)
  real(DP) :: dens(NSal, Ntheta, Npress)

  !
  call set_ThermoVarRange()
  call prepair_output()

  !
  call eval_EOS()

  !
  call HistoryPut('dens', dens)
  call HistoryClose()

contains
subroutine eval_EOS()
  integer :: i,j,k

  call EqState_JM95_Init()
  
  do i=1,Nsal
     do j=1, Ntheta
        do k=1, Npress
           dens(i,j,k) = EqState_JM95_Eval(theta=ptemp(j), s=sal(i), p=press(k))
        end do
     end do
  end do

  dens = dens - 1000d0

  call EqState_JM95_Final()
end subroutine eval_EOS

subroutine set_ThermoVarRange()
  
  integer :: m

  do m=1, Npress
     ptemp(m) = thetaMin + (m-1)*(thetaMax - thetaMin)/real(Ntheta-1)
  end do

  do m=1, NSal
     sal(m) = SMin + (m-1)*(SMax - SMin)/real(NSal-1)
  end do

  do m=1, NPress
     press(m) = pMin + (m-1)*(pMax - pmin)/real(NPress-1)
  end do
end subroutine set_ThermoVarRange

subroutine prepair_output()

  character(STRING) :: longnames(3), units(3)

  longnames(1) = "salinity"
  longnames(2) = "potential temperature"
  longnames(3) = "pressure"

  units(1) = "psu"
  units(2) = "deg C"
  units(3) = "bar"

  call HistoryCreate( &                        ! ヒストリー作成
    & file='dens.nc', title='Density of sea water', &
    & source='Density anomaly of sea water(sigma)',   &
    & institution='GFD_Dennou Club ogcm project',   &
    & dims=(/'sal  ','ptemp', 'press' /), &
    & dimsizes=(/ NSal,Ntheta,NPress /), &
    & longnames=longnames, units=units &
    & )

  call HistoryPut('sal', sal)
  call HistoryPut('press', press)
  call HistoryPut('ptemp', ptemp)

  call HistoryAddVariable( &
       & 'dens', (/ 'sal  ', 'ptemp', 'press' /), &
       & 'density anomaly of sea water calculated by JM95 EOS', 'kg*m-3', 'double' )
 
end subroutine prepair_output

end program EOS_mod_test
