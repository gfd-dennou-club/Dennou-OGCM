program SGSConvAdjust_test
  
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

  use SGSConvAdjust_mod

  use SpmlUtil_mod, only: &
       & SpmlUtil_Init, SpmlUtil_Final, &
       & g_Sig, IntSig_BtmToTop

  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final

  use EOS_Linear_mod, only: &
       & EOSTYPE_LINEAR

  ! Declaration statement
  implicit none

  character(*),  parameter :: PROGRAM_NAME = "SGSConvAdjust_test"
  character(*), parameter :: configNmlFile = "defaultConfig_z1D.nml"
  real(DP), parameter :: totDepth = 5.2d3

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
  call SGSConvAdjust_Final()

contains
subroutine initialize()

  integer :: k

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call EOSDriver_Init(EOSTYPE_LINEAR)
  call SGSConvAdjust_Init()

  allocate(z(0:kMax))
  allocate(xy_totDepth(0:iMax-1,jMax))

  xy_totDepth(:,:) = totDepth
  z(:) = g_Sig(:)*totDepth

end subroutine initialize


subroutine check_StaticStableColumn()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt
  logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_isAdjustOccur

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
  call SGSConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xyz_isAdjustOccur)
  write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)

  call  AssertEqual( &
       &  message="Check that the convective adjustment has not occured.", &
       &  answer=.false., check=xyz_isAdjustOccur(0,1,kMax/2) )

end subroutine check_StaticStableColumn

subroutine check_StaticUnStableColumn1()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt
  logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_isAdjustOccur
  
  integer :: k
  real(DP), parameter :: TempSurf = 300d0
  real(DP), parameter :: delTheta = 1d0

  xyz_Salt = 35d0
  do k=0, kMax
     xyz_PTemp(:,:,k) = TempSurf - z(k)/totDepth*delTheta
  end do

  write(*,*) "IntPTemp before..", IntSig_BtmToTop(xyz_PTemp(0,1,:))

  write(*,*) "PTemp before the adjustment:", xyz_PTemp(0,1,:)
  call SGSConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xyz_isAdjustOccur)
  write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)

  write(*,*) xyz_isAdjustOccur(:,:,:)
  call  AssertEqual( &
       &  message="Check that the convective adjustment has occured.", &
       &  answer=.true., check=xyz_isAdjustOccur(0,1,kMax/2) )

  call AssertLessThan( &
       & message="Check the profile of potential temperature after the convective adjustment.", &
       & answer=1d-12, check=IntSig_BtmToTop(xyz_PTemp(0,1,:)-300.5d0)/IntSig_BtmToTop(xyz_PTemp(0,1,:)) &
       & )

  write(*,*) "IntPTemp after..", IntSig_BtmToTop(xyz_PTemp(0,1,:))

end subroutine check_StaticUnStableColumn1

subroutine check_StaticUnStableColumn2()

  real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
       & xyz_PTemp, xyz_Salt
  logical, dimension(0:iMax-1,jMax, 0:kMax) :: xyz_isAdjustOccur
  
  integer :: k
  real(DP), parameter :: TempSurf = 300d0
  real(DP), parameter :: delTheta = 1d0
  real(DP), parameter :: adjustedPTemp = (delTheta/3d0 + TempSurf)

  real(DP) :: IntPTempBefore

  xyz_Salt = 35d0
  do k=0, kMax
     xyz_PTemp(:,:,k) = delTheta/totDepth**2 * z(k)**2 + TempSurf
  end do

  IntPTempBefore = IntSig_BtmToTop(xyz_PTemp(0,1,:))
  write(*,*) "IntPTemp before..", IntPTempBefore

  write(*,*) "PTemp before the adjustment:", xyz_PTemp(0,1,:)
  call SGSConvAdjust_perform(xyz_PTemp, xyz_Salt, xy_totDepth, xyz_isAdjustOccur)
  write(*,*) "PTemp after the adjustment:", xyz_PTemp(0,1,:)

  call  AssertEqual( &
       &  message="Check that the convective adjustment has occured.", &
       &  answer=.true., check=xyz_isAdjustOccur(0,1,kMax/2) )

  call MessageNotify('M', PROGRAM_NAME, &
       & "The adjusted potential temperature must be  %f[K]..", d=(/ adjustedPTemp /)) 
  call AssertLessThan( &
       & message="Check the profile of potential temperature after the convective adjustment.", &
       & answer=1d-5, check=IntSig_BtmToTop(xyz_PTemp(0,1,:)-adjustedPTemp)/IntSig_BtmToTop(xyz_PTemp(0,1,:)) &
       & )

  write(*,*) "IntPTemp after..", abs(IntPTempBefore-IntSig_BtmToTop(xyz_PTemp(0,1,:)))/IntPTempBefore

end subroutine check_StaticUnStableColumn2

end program SGSConvAdjust_test
