program SGSSlowConvAdjust_test
  
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
       & iMax, jMax, kMax, nMax, lMax, tMax

  use SGSSlowConvAdjust_mod, only: &
       & SGSSlowConvAdjust_Init, SGSSlowConvAdjust_Final, &
       & SGSSlowConvAdjust_perform, &
       & DEFAULT_MIXTIME, DEFAULT_CONVLAYER_DEPTH

  use SpmlUtil_mod, only: &
       & SpmlUtil_Init, SpmlUtil_Final, &
       & g_Sig, IntSig_BtmToTop

  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final

  use EOS_Linear_mod, only: &
       & EOSTYPE_LINEAR

  ! Declaration statement
  implicit none

  character(*),  parameter :: PROGRAM_NAME = "SGSSlowConvAdjust_test"
  character(*), parameter :: configNmlFile = "defaultConfig_z1D.nml"
  real(DP), parameter :: totDepth = 5.2d3
  real(DP), parameter :: GCMTimeStep = 2*3600d0
  real(DP), parameter :: ConvMixTime = DEFAULT_MIXTIME
  real(DP), allocatable :: xy_totDepth(:,:), z(:)


  !***********************************************************************
  ! Executable statement
  !

  ! Initialize some modules. 
  call initialize()
  
  !
  call check_StaticStableColumn()
  call check_StaticUnstableColumn1()
  call check_StaticUnstableColumn2()

  ! Finalize
  call SGSSlowConvAdjust_Final()

contains
subroutine initialize()
  integer :: k

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call EOSDriver_Init(EOSTYPE_LINEAR)
  call SGSSlowConvAdjust_Init(GCMTimeStep, ConvMixTime, DEFAULT_CONVLAYER_DEPTH)

  allocate(z(0:kMax))
  allocate(xy_totDepth(0:iMax-1,jMax))

  xy_totDepth(:,:) = totDepth
  z(:) = g_Sig(:)*totDepth

end subroutine initialize


subroutine check_StaticStableColumn()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt
  logical, dimension(0:iMax-1,jMax) :: xy_isAdjustOccur

  integer :: k
  real(DP), parameter :: TempSurf = 300d0
  real(DP), parameter :: delTheta = 1d0

  xyz_Salt = 35d0
  do k=0, kMax
     xyz_PTemp(:,:,k) = TempSurf + z(k)/totDepth*delTheta
  end do

  call MessageNotify('M', PROGRAM_NAME, &
       & "Check convective adjustment routine for static stable column..")

  write(*,*) "PTemp before the adjustment:", xyz_PTemp(0,1,:)
  call SGSSlowConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xy_isAdjustOccur)
  write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)

  call  AssertEqual( &
       &  message="Check that the convective adjustment has not occured.", &
       &  answer=.false., check=xy_isAdjustOccur(0,1) )

end subroutine check_StaticStableColumn

subroutine check_StaticUnStableColumn1()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt, xyz_PTempOld
  logical, dimension(0:iMax-1,jMax) :: xy_isAdjustOccur
  
  integer :: k, n
  real(DP), parameter :: TempSurf = 300d0
  real(DP), parameter :: delTheta = 1d0
  real(DP) :: Res

  xyz_Salt = 35d0
  do k=0, kMax
     xyz_PTemp(:,:,k) = TempSurf - z(k)/totDepth*delTheta
  end do
  xyz_PTempOld = xyz_PTemp
  write(*,*) "Inital PTemp: ", xyz_PTempOld(0,1,:)
!  call Output(xyz_PTemp(0,1,:), 0)

  do n=1, int(20*ConvMixTime/GCMTimeStep)
     call SGSSlowConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xy_isAdjustOccur)
     Res = IntSig_BtmToTop(xyz_PTemp(0,1,:)-300.5d0)/IntSig_BtmToTop(xyz_PTemp(0,1,:))
     xyz_PTempOld = xyz_PTemp


     if(mod(int(n*GCMTimeStep/3600d0),4)==0) then
        write(*,*) "t=", int(n*GCMTimeStep/3600d0), "[h].., Res=", Res
        write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)
!        call Output(xyz_PTemp(0,1,:), n)
     end if

  end do

  call AssertLessThan( &
       & message="Check the profile of potential temperature after the convective adjustment.", &
       & answer=1d-8, check=Res &
       & )
end subroutine check_StaticUnStableColumn1

subroutine check_StaticUnStableColumn2()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt, xyz_PTempOld
  logical, dimension(0:iMax-1,jMax) :: xy_isAdjustOccur
  
  integer :: k, n
  real(DP), parameter :: TempSurf = 300d0
  real(DP), parameter :: delTheta = 1d0
  real(DP), parameter :: adjustedPTemp = (delTheta/3d0 + TempSurf)
  real(DP) :: Res
  character(STRING) :: message, fileName
  
  xyz_Salt = 35d0
  do k=0, kMax
     xyz_PTemp(:,:,k) = delTheta/totDepth**2 * z(k)**2 + TempSurf
  end do
  xyz_PTempOld = xyz_PTemp
  write(*,*) "Inital PTemp: ", xyz_PTempOld(0,1,:)
  call Output(xyz_PTemp(0,1,:), 0)

  do n=1, int(20*ConvMixTime/GCMTimeStep)
     call SGSSlowConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xy_isAdjustOccur)
     xyz_PTempOld = xyz_PTemp
     Res = IntSig_BtmToTop(xyz_PTemp(0,1,:)-adjustedPTemp)/IntSig_BtmToTop(xyz_PTemp(0,1,:))


     if(mod(int(n*GCMTimeStep/3600d0),4)==0) then
        write(*,*) "t=", int(n*GCMTimeStep/3600d0), "[h].., Res=", Res
        write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)
        call Output(xyz_PTemp(0,1,:), n)
     end if
  end do

  call MessageNotify('M', PROGRAM_NAME, &
       & "The adjusted potential temperature must be  %f[K]..", d=(/ adjustedPTemp /)) 
  call AssertLessThan( &
       & message="Check the profile of potential temperature after the convective adjustment.", &
       & answer=1d-5, check=Res &
       & )
end subroutine check_StaticUnStableColumn2

subroutine output(z_PTemp, n)
  real(DP), intent(in) :: z_PTemp(0:kMax)
  integer, intent(in) :: n

  integer :: k
  character(TOKEN) :: fileName

  fileName = CPrintf("SlowConv_%d.dat", i=(/ n /))
  open(10, file=trim(fileName), status="replace")
  write(10,*) "# Sig, PTemp"
  do k=0, kMax
     write(10,*) g_Sig(k), z_PTemp(k)
  end do
  close(10)

end subroutine output

end program SGSSlowConvAdjust_test
