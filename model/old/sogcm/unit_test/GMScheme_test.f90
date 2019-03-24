program GMScheme_test
  
  ! Use statement
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use dc_string, only: &
       & CPrintf

  use dc_test

  use Constants_mod, only: &
       & Constants_Init, Constants_Final, &
       & PI, RPlanet

  use GridSet_mod, only: &
       & GridSet_Init, GridSet_Final, &
       & GridSet_construct, &
       & iMax, jMax, kMax, nMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon


  use SpmlUtil_mod, only: &
       & SpmlUtil_Init, SpmlUtil_Final, &
       & g_Sig, IntSig_BtmToTop, wz_xyz, xyz_wz

  use VariableSet_mod

  use TemporalIntegUtil_mod2

  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final

  use EOS_Linear_mod, only: &
       & EOSTYPE_LINEAR

  use SGSEddyMixing_mod

  use Exp_APEOGCirc_mod, only: &
       & Exp_Init => Exp_APEOGCirc_Init, &
       & Exp_Final => Exp_APEOGCirc_Final, &
       & Exp_SetInitCond => SetInitCondition

  use InitCond_mod

  use gtool_history

  ! Declaration statement
  implicit none

  character(*),  parameter :: PROGRAM_NAME = "GMScheme_test"
  character(*), parameter :: configNmlFile = "defaultConfig_axisym_GM.nml"


  real(DP), parameter :: outputIntrVal = 1000d0*86400d0
  real(DP), parameter :: EndTime = 2000d0*365d0*86400d0
  real(DP), parameter :: DelTime = 4*24d0*3600d0
  real(DP) :: CurrentTime
  real(DP), allocatable :: xyz_DensPot(:,:,:)
  
  !***********************************************************************
  ! Executable statement
  !

  ! Initialize some modules. 
  call initialize()

  CurrentTime = 0d0
  call set_initialPerturbField()
  call output_data()


  !
  do while(CurrentTime <= EndTime)
     if(mod(CurrentTime, DelTime*500d0)==0) &
          & write(*,*) "Current time=", CurrentTime/86400d0, "[day]"

     call advance_fieldData()
     CurrentTime = CurrentTime + DelTime
     if(mod(CurrentTime, outputIntrVal)==0) call output_data()
     
     xyz_PTempEddN = xyz_PTempEddA; xyz_SaltN = xyz_SaltA
  end do

  ! Finalize
  call SGSEddyMixing_Final()

contains

subroutine advance_fieldData()

  use SGSConvAdjust_mod
  
  real(DP) :: wz_PTempBasic(lMax,0:kMax), xyz_PTemp(0:iMax-1,jMax,0:kMax)
  logical :: xyz_isAdjustOccur(0:iMax-1,jMax,0:kMax)

  real(DP), dimension(lMax,0:kMax) :: &
       & wz_PTempRKTmp, wz_PTempEdd, wz_SaltRKTmp, wz_Salt, &
       & wz_PTempRHS, wz_SaltRHS
  integer :: k, stage
  
  !
  wz_PTempBasic = wz_xyz(spread(spread(z_PTempBasic,1,jMax),1,iMax))
  forAll(k=0:kMax) xyz_PTemp(:,:,k) = z_PTempBasic(k) + xyz_PTempEddN(:,:,k)

  call SGSConvAdjust_perform(xyz_PTemp, xyz_SaltN, xy_totDepthBasic, xyz_isAdjustOccur)

  forAll(k=0:kMax) xyz_PTempEddN(:,:,k) = xyz_PTemp(:,:,k) - z_PTempBasic(k)

  wz_PTempEdd = wz_xyz(xyz_PTempEddN); wz_Salt = wz_xyz(xyz_SaltN)
  do Stage=1, 2
     wz_PTempRHS = 0d0; wz_SaltRHS = 0d0

     call SGSEddyMixing_AddMixingTerm(wz_PTempRHS, wz_SaltRHS, &
       & wz_PTempBasic+wz_PTempEdd, wz_Salt, xy_totDepthBasic)

     wz_PTempEdd = timeIntRK(wz_xyz(xyz_PTempEddN), wz_PTempRHS, 2, Stage, wz_PTempRKTmp)
     wz_Salt = timeIntRK(wz_xyz(xyz_SaltN), wz_SaltRHS, 2, Stage, wz_SaltRKTmp)
  end do

  xyz_PTempEddA = xyz_wz(wz_PTempEdd)
  xyz_SaltA = xyz_wz(wz_Salt)

end subroutine advance_fieldData

subroutine initialize()

  integer :: k

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  
  call GridSet_construct()
  call VariableSet_Init()

  call EOSDriver_Init(EOSTYPE_LINEAR)
  call SGSEddyMixing_Init( SGSEddyMixParamType=SGSEddyMixing_GM90, &
       & KappaRedi=1000d0, KappaGM=0d0, isVarsOutput=.true.) 

  call TemporalIntegUtil_Init(DelTime)

  !
  allocate(xyz_DensPot(0:iMax-1,jMax,0:kMax))

  !
  call Exp_Init(configNmlFile)
  call InitCond_Init()
  call InitCond_Set(Exp_SetInitCond)


  !
  call HistoryCreate( &               
    & file='GMScheme_test.nc', title='GMScheme tests', &
    & source='GMScheme_test',   &
    & institution='dcmodel project',       &
    & dims=(/'lon ', 'lat ', 'sig ', 'time'/), dimsizes=(/ iMax, jMax, kMax+1, 0 /), &
    & longnames=(/ 'longitude', 'latitude ', 'sig      ', 'time     '/),           &
    & units=(/ 'degree_east ', 'degree_north', '(1)         ', 's           '/),   &
    & origin=real(0), interval=real(outputIntrVal) )

  call HistoryPut('lon', xyz_lon(:,1,0)*180d0/PI)
  call HistoryPut('lat', xyz_lat(0,:,0)*180d0/PI)
  call HistoryPut('sig', g_Sig)

  call HistoryAddVariable( &
    & varname='PTemp', dims=(/ 'lon ', 'lat ', 'sig ', 'time' /), &
    & longname='potential temperature', units='K')

  call HistoryAddVariable( &
    & varname='Salt', dims=(/ 'lon ', 'lat ', 'sig ', 'time' /), &
    & longname='salinity', units='psu')
  
  call HistoryAddVariable( &
    & varname='DensPot', dims=(/ 'lon ', 'lat ', 'sig ', 'time' /), &
    & longname='potential density', units='kg/m3')

  call SGSEddyMixing_PrepareOutput(0d0, EndTime, outputIntrVal, "")

end subroutine initialize

subroutine set_initialPerturbField()

  integer :: k
  real(DP) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
  logical :: xyz_isAdjustOccur(0:iMax-1,jMax,0:kMax)

  real(DP), parameter :: PTempEddAmplFactor = 2.5d0

  do k=0, kMax
     xyz_SaltN = 36d0
     xyz_PTempEddN(:,:,k) = PTempEddAmplFactor*sin(2d0*xyz_Lat(:,:,k))**2*sin(-PI*g_Sig(k))
     xyz_PTemp(:,:,k) = z_PTempBasic(k) + xyz_PTempEddN(:,:,k)
  end do


end subroutine set_initialPerturbField

subroutine output_data()
  real(DP) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
  integer :: k

  forAll(k=0:kMax) &
       & xyz_PTemp(:,:,k) = xyz_PTempEddN(:,:,k) + z_PTempBasic(k)

  call HistoryPut('PTemp', xyz_PTemp)
  call HistoryPut('Salt', xyz_SaltN)
!  call HistoryPut('DensPot', xyz_DensPot)

  call SGSEddyMixing_Output()

end subroutine output_data

end program GMScheme_test
