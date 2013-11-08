program HydroBouEqSolverVDiffProc_mod_test
  use dc_types
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use at_module
  use HydroBouEqSolverVDiffProc_mod, only: &
       & HydroBouEqSolverVDiffProc_Init, HydroBouEqSolverVDiffProc_Final, &
       & construct_vDiffProcMat, vDiffImplicitSolve => Solve
  use dc_test
  use dc_message
  use dc_string

  implicit none

  character(*), parameter :: configNmlFile = "defaultConfig.nml"
  real(DP), parameter :: Kappa = 1d-02
  real(DP), allocatable :: xy_totDepth(:,:)
  real(DP), allocatable ::  xyz_Temp(:,:,:), wt_Temp(:,:)
  integer :: i,j,n
  integer, parameter :: nStep = 1000
  real(DP) :: dt
  real(DP) :: l2Error

  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call GridSet_construct()
  call HydroBouEqSolverVDiffProc_Init()

  !
  allocate(xy_totDepth(0:iMax-1,jMax))
  allocate(xyz_Temp(0:iMax-1,jMax,0:kMax), wt_Temp(lMax,0:tMax))

  xy_totDepth = 1d0
  dt = (2d0/Kappa)/real(nStep)
!!$  forAll(i=0:iMax-1,j=1:jMax) xyz_Temp(i,j,:) = sin(-PI*g_Sig)
  forAll(i=0:iMax-1,j=1:jMax) xyz_Temp(i,j,:) = cos(-PI*g_Sig)

  call MessageNotify('M', 'HydroBouEqSolverVDiffProc_mod_test', &
       & 'construct_vDiffProcMat..')
!!$  call construct_vDiffProcMat(0.5*Kappa, dt, xy_totDepth, &
!!$       & 'D', 'D' )
  call construct_vDiffProcMat(0.5*Kappa, dt, xy_totDepth, &
       & 'N', 'N' )


  wt_Temp = wt_xyz(xyz_Temp)
  do n=1, nStep

     xyz_Temp = xyz_wt( wt_Temp + 0.5d0*dt*Kappa*at_Dx_at(at_Dx_at(wt_Temp)) )
     xyz_Temp(:,:,0) = 0d0
     xyz_Temp(:,:,kMax) = 0d0
     wt_Temp = vDiffImplicitSolve(xyz_Temp)

     if(n==1 .or. mod(n,int(nStep*0.1)) == 0) then
        xyz_Temp = xyz_wt(wt_Temp)
!!$        l2Error = sqrt(sum( (xyz_Temp(1,1,:) - sin(-PI*g_Sig)*exp(-Kappa*PI**2*(n*dt)))**2 ))/(1d0/PI)
        l2Error = sqrt(sum( (xyz_Temp(1,1,:) - cos(-PI*g_Sig)*exp(-Kappa*PI**2*(n*dt)))**2 ))/(1d0/PI)

        call MessageNotify("M", "HydroBouEqSolverVDiffProc_mod_test", &
             & "t=%f, l2ErrorNorm=%f", d=(/ n*dt, l2Error /))
        write(*,*) xyz_Temp(1,1,:)
        write(*,*) cos(-PI*g_Sig)*exp(-Kappa*PI**2*(n*dt))
!!$        write(*,*) sin(-PI*g_Sig)*exp(-Kappa*PI**2*(n*dt))
     end if

  end do



  call HydroBouEqSolverVDiffProc_Final()
  call SpmlUtil_Final()
  call GridSet_Final() 
  call Constants_Final()

contains

end program HydroBouEqSolverVDiffProc_mod_test
