program EOS_mod_test

  ! モジュール引用; Use statement
  !
  use dc_test
  use dc_message
  use dc_string
  use dc_types
  use gtool_history

  use UnitConversion_mod, only: &
       & bar2Pa, degC2K

  use EOS_Linear_mod, only: &
       & EOSTYPENAME_LINEAR, &
       & EOS_Linear_Init, EOS_Linear_Final, &
       & EOS_Linear_Eval, &
       & EOS_Linear_PTemp2Temp, EOS_Linear_Temp2PTemp, &
       & EOS_Linear_HeatCapacity

  use EOS_SimpleNonLinear_mod, only: &
       & EOSTYPENAME_SIMPLENONLINEAR, &
       & EOS_SimpleNonLinear_Init, EOS_SimpleNonLinear_Final, &
       & EOS_SimpleNonLinear_Eval, &
       & EOS_SimpleNonLinear_PTemp2Temp, EOS_SimpleNonLinear_Temp2PTemp, &
       & EOS_SimpleNonLinear_HeatCapacity

  use EOS_JM95_mod, only: &
       & EOSTYPENAME_JM95, &
       & EOS_JM95_Init, EOS_JM95_Final, &
       & EOS_JM95_Eval, &
       & EOS_JM95_PTemp2Temp, EOS_JM95_Temp2PTemp, &
       & EOS_JM95_HeatCapacity

  implicit none

  ! 宣言文; Declaration statement
  !
  character(*), parameter :: configNmlFile = "defaultConfig.nml"
  integer, parameter :: Npress = 51
  integer, parameter :: NSal = 43 
  integer, parameter :: Ntheta = 43
  real(DP), parameter :: pMin = 0d0, pMax = 1000d0
  real(DP), parameter :: thetaMin = -2d0, thetaMax = 40d0
  real(DP), parameter :: SMin = 0d0, SMax = 42d0
  real(DP), parameter :: RefDens = 1.027d03  
  real(DP) :: ptemp(Ntheta), temp(Ntheta)
  real(DP) :: press(Npress)
  real(DP) :: sal(NSal)

  type :: EOS
     character(TOKEN) :: name
     real(DP), dimension(NSal, Ntheta, Npress) :: &
          & dens, Cp
     real(DP), dimension(NSal, Ntheta, Npress) :: &
          & PTempTempDiff, PTemp
  end type EOS

  type(EOS), save :: EOSs(3) 
  type(gt_history) :: hst_density, hst_PTempTempDiff, hst_HeatCapacity

  integer :: EOSId

  ! 実行文; Executable statement
  !

  !

  ! Initialize each module of EOS for seawater. 
  !
  call EOS_Linear_Init()
  call EOS_SimpleNonLinear_Init()
  call EOS_JM95_Init()

  EOSs(1)%name = EOSTYPENAME_LINEAR
  EOSs(2)%name = EOSTYPENAME_SimpleNonLINEAR
  EOSs(3)%name = EOSTYPENAME_JM95


  ! Set the coordinate using thermodynamic variables. 
  call set_ThermoVarCoord()

  ! Preparation for outputing data of density field. 
  call prepair_output(hst_density, "dens", "density annomaly" )
  call prepair_output(hst_HeatCapacity, "Cp", "heat capacity", isPTemp=.false. )
  call prepair_output(hst_PTempTempDiff, "PTempTempDiff", &
       & "difference between potential temperature and temperature", isPTemp=.false. )


  ! Calculate the density field using each EOS module.
  call eval_EOS()

  ! Output the calculated density field. 
  do EOSId=1, 3
     call HistoryPut(EOSs(EOSId)%name, EOSs(EOSId)%dens, hst_density)
     call HistoryPut(EOSs(EOSId)%name, EOSs(EOSId)%Cp, hst_HeatCapacity)
     call HistoryPut(EOSs(EOSId)%name, EOSs(EOSId)%PTempTempDiff, hst_PTempTempDiff)
  end do

  call HistoryClose(hst_density)
  call HistoryClose(hst_HeatCapacity)
  call HistoryClose(hst_PTempTempDiff)

  ! Finalize each module of EOS for seawater.
  call EOS_Linear_Final()
  call EOS_SimpleNonLinear_Final()
  call EOS_JM95_Final()
  
contains

subroutine eval_EOS()  

  integer :: i,j,k
  
  do i=1,Nsal
     do j=1, Ntheta
        do k=1, Npress

           ! Obtain perturbed component of density.
           EOSs(1)%dens(i,j,k) = &
                & EOS_Linear_Eval( theta=degC2K(ptemp(j)), s=sal(i), p=bar2Pa(press(k)) ) &
                & + RefDens - 1000d0

           EOSs(1)%Cp(i,j,k) = &
                & EOS_Linear_HeatCapacity( T=degC2K(temp(j)), S=sal(i), P=bar2Pa(press(k)) )

           EOSs(1)%PTempTempDiff(i,j,k) = &
                & EOS_Linear_Temp2PTemp( T0=degC2K(temp(j)), S0=sal(i), P0=bar2Pa(press(k)), Pref=0d0 ) - degC2K(temp(j))


           ! Obtain perturbed component of density.
           EOSs(2)%dens(i,j,k) = &
                & EOS_SimpleNonLinear_Eval( theta=degC2K(ptemp(j)), s=sal(i), p=bar2Pa(press(k)) ) &
                & + RefDens -1000d0

           EOSs(2)%Cp(i,j,k) = &
                & EOS_SimpleNonLinear_HeatCapacity( T=degC2K(temp(j)), S=sal(i), P=bar2Pa(press(k)) )

           EOSs(2)%PTempTempDiff(i,j,k) = &
                & EOS_SimpleNonLinear_Temp2PTemp( T0=degC2K(temp(j)), S0=sal(i), P0=bar2Pa(press(k)), Pref=0d0 ) - degC2K(temp(j))

           ! Obtain total density.
           EOSs(3)%dens(i,j,k) = &
                & EOS_JM95_Eval(theta=ptemp(j), s=sal(i), p=press(k)) &
                & - 1000d0

           EOSs(3)%Cp(i,j,k) = &
                & EOS_JM95_HeatCapacity( T=temp(j), S=sal(i), P=press(k) )

           EOSs(3)%PTempTempDiff(i,j,k) = &
                & EOS_JM95_Temp2PTemp( T0=temp(j), S0=sal(i), P0=press(k), Pref=0d0 ) - temp(j) 

        end do
     end do
  end do

end subroutine eval_EOS


subroutine set_ThermoVarCoord()
  
  integer :: m

  do m=1, Npress
     ptemp(m) = thetaMin + (m-1)*(thetaMax - thetaMin)/real(Ntheta-1)
  end do

  temp(:) = ptemp(:)

  do m=1, NSal
     sal(m) = SMin + (m-1)*(SMax - SMin)/real(NSal-1)
  end do

  do m=1, NPress
     press(m) = pMin + (m-1)*(pMax - pmin)/real(NPress-1)
  end do
end subroutine set_ThermoVarCoord

subroutine prepair_output(hst, varName, varLongName, isPTemp)

  ! 宣言文; Declaration statement
  !
  type(gt_history), intent(inout)  :: hst
  character(*), intent(in) :: varName, varLongName
  logical, intent(in), optional :: isPTemp

  ! 作業変数
  ! Work variable
  !
  character(TOKEN) :: dims(3), units(3)
  character(TOKEN) :: longDimNames(3), longName
  integer :: id
  logical :: usePTempAxis

  ! 実行文; Executable statement
  ! 

  usePTempAxis = .true.
  if(present(isPTemp)) usePTempAxis = isPTemp

  dims(1) = 'sal';  
  longDimNames(1) = "salinity"
  if(usePTempAxis) then
     dims(2) = 'ptemp'
     longDimNames(2) = "potential temperature"
  else
     dims(2) = 'temp'
     longDimNames(2) = "temperature"
  end if
  dims(3) = 'press'
  longDimNames(3) = "pressure"

  units(1) = "psu"
  units(2) = "deg C"
  units(3) = "bar"

  call HistoryCreate( &                        ! ヒストリー作成
    & file=trim(varName)//".nc", title=trim(varName), &
    & source='Density anomaly of sea water(sigma)',   &
    & institution='GFD_Dennou Club ogcm project',   &
    & dims=dims, &
    & dimsizes=(/ NSal,Ntheta,NPress /), &
    & longnames=longDimNames, units=units, history=hst &
    & )


  call HistoryPut(dims(1), sal, hst)
  if(usePTempAxis) then
     call HistoryPut(dims(2), ptemp, hst)
  else
     call HistoryPut(dims(2), temp, hst)
  end if
  call HistoryPut(dims(3), press, hst)

  do id=1, 3
     longName = CPrintf("%a calculated by %a", ca=(/ trim(varLongName), trim(EOSs(id)%name) /))
     call HistoryAddVariable( &
          & trim(EOSs(id)%name), dims, trim(longName), &
          & 'kg*m-3', 'double', history=hst )
  end do
end subroutine prepair_output

end program EOS_mod_test
