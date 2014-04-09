program HydroBouEqSolverVDiffProc_mod_test
  
  ! Use statement
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use dc_string, only: &
       & CPrintf

  use Constants_mod, only: &
       & Constants_Init, Constants_Final, &
       & PI, RPlanet

  use GridSet_mod, only: &
       & GridSet_Init, GridSet_Final, &
       & GridSet_construct, &
       & iMax, jMax, kMax, nMax, lMax, tMax

  use SpmlUtil_mod, only: &
       & SpmlUtil_Init, SpmlUtil_Final

  ! Declaration statement
  implicit none


  character(*), parameter :: PROGRAM_NAME = "HydroBouEqSolverVDiffProc_mod_test"
#ifdef DSOGCM_MODE_AXISYM
  character(*), parameter :: configNmlFile = "defaultConfig_axisym.nml"
#else
  character(*), parameter :: configNmlFile = "defaultConfig.nml"
#endif

  real(DP), parameter :: EndNonDimTime = 4d0        ! normalized by the corresponding to e-folding time. 
  real(DP), parameter :: defaultDtRatio = 5d-04     ! relative to the corresponding to e-folding time.

  integer, parameter :: TS_FWEuler = 1   ! Forward Euler scheme
  integer, parameter :: TS_BWEuler = 2   ! Backward Euler scheme 
  integer, parameter :: TS_CRANK_NICO = 3 ! Crank Nicolson method
 
  !***********************************************************************
  ! Executable statement
  !

  ! Initialization
  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call GridSet_construct()

  ! Test each case. 
  !

  call test_case1()
  call test_case2()
  call test_case3()
  call test_case4()

  ! Finalization
  call SpmlUtil_Final()
  call GridSet_Final() 
  call Constants_Final()

contains

!*****************************************************

  ! Test case 1
  !
  subroutine test_case1()
    call MessageNotify('M', PROGRAM_NAME, &
         & 'Solve 1D diffusion equation with homogeneous bc (D and D)..')

    call perform_timeIntegration( 1d0/PI**2, 'D', 0d0, 'D', 0d0, get_anaSol_homo_DD1, TS_FWEuler,    1d-04 )
    call perform_timeIntegration( 1d0/PI**2, 'D', 0d0, 'D', 0d0, get_anaSol_homo_DD1, TS_BWEuler,    1d-04 )
    call perform_timeIntegration( 1d0/PI**2, 'D', 0d0, 'D', 0d0, get_anaSol_homo_DD1, TS_CRANK_NICO, 5d-9  )
  end subroutine test_case1

  pure function get_anaSol_homo_DD1(x,t) result(anaSol)
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: t
    real(DP) :: anaSol(size(x))
    
    anaSol = exp(-PI**2*t)*sin(PI*x)
  end function get_anaSol_homo_DD1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Test case 2
  !
  subroutine test_case2()
    call MessageNotify('M', PROGRAM_NAME, &
         & 'Solve 1D diffusion equation with homogeneous bc (N and N)..')
  
    call perform_timeIntegration( 1d0/PI**2, 'N', 0d0, 'N', 0d0, get_anaSol_homo_NN1, TS_FWEuler,    1d-04 )
    call perform_timeIntegration( 1d0/PI**2, 'N', 0d0, 'N', 0d0, get_anaSol_homo_NN1, TS_BWEuler,    1d-04 )
    call perform_timeIntegration( 1d0/PI**2, 'N', 0d0, 'N', 0d0, get_anaSol_homo_NN1, TS_CRANK_NICO, 5d-9  )
  end subroutine test_case2


  pure function get_anaSol_homo_NN1(x,t) result(anaSol)
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: t
    real(DP) :: anaSol(size(x))
    
    anaSol = exp(-PI**2*t)*cos(PI*x)
  end function get_anaSol_homo_NN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Test case 3
  !
  subroutine test_case3()
    call MessageNotify('M', PROGRAM_NAME, &
         & 'Solve 1D diffusion equation with homogeneous bc (D and N)..')
  
    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'D', 0d0, 'N', 0d0, get_anaSol_homo_ND1, TS_FWEuler,    1d-04, &
         & dtRatio_=0.25*defaultDtRatio )
    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'D', 0d0, 'N', 0d0, get_anaSol_homo_ND1, TS_BWEuler,    1d-04 )
    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'D', 0d0, 'N', 0d0, get_anaSol_homo_ND1, TS_CRANK_NICO, 5d-9 )

  end subroutine test_case3


  pure function get_anaSol_homo_ND1(x,t) result(anaSol)
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: t
    real(DP) :: anaSol(size(x))
    
    anaSol = exp(-(0.5d0*PI)**2*t)*cos(0.5d0*PI*x)
  end function get_anaSol_homo_ND1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Test case 4
  !
  subroutine test_case4()

    call MessageNotify('M', PROGRAM_NAME, &
         & 'Solve 1D diffusion equation with inhomogeneous bc (D and N)..')

    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'N', -1d0, 'D', 0d0, get_anaSol_inhomo_DN1, TS_FWEuler,    1d-04, &
         & dtRatio_=0.25*defaultDtRatio )
    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'N', -1d0, 'D', 0d0, get_anaSol_inhomo_DN1, TS_BWEuler,    1d-04 )

    ! Although The Crank-Nocolson scheme has second-order accuracy for time discritization, 
    ! the same amount of threshold as low-order scheme is set beacause the effective digits 
    ! of reference solution obtained by get_anaSol_inhomo_DN1 is only 6. 
    call perform_timeIntegration( 1d0/(0.5d0*PI)**2, 'N', -1d0, 'D', 0d0, get_anaSol_inhomo_DN1, TS_CRANK_NICO, 1d-4 )

  end subroutine test_case4

  !> This function returns analystic solution calculated with the sum of series. 
  !> Note that if kMax=2, the maximum number of significant digits for returned solution is *6*.  
  pure function get_anaSol_inhomo_DN1(x,t) result(anaSol)
    real(DP), intent(in) :: x(:)
    real(DP), intent(in) :: t
    real(DP) :: anaSol(size(x))

    integer, parameter :: nMax = 500000
    integer :: n
    real(DP) :: sig
    real(DP) :: coef

    sig = -1d0
    anaSol = 0d0
    do n=0, nMax, 1
       coef = (dble(n) + 0.5d0)*PI
       sig = sign(1d0,-sig)
       anaSol = anaSol + &
            & 2d0*sig/(coef**2) * exp(-(coef**2)*t) * sin(coef*x)
    end do
    anaSol =  anaSol - x
  end function get_anaSol_inhomo_DN1


!****************************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solver for 1D diffusion equation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine perform_timeIntegration( &
       & eFoldTime, & 
       & lBCType, lBCVal, uBCType, uBCVal, &
       & get_anaSolution, &
       & tschemeId, errorLimit, dtRatio_ )

    ! Use statement
    !
    use HydroBouEqSolverVDiffProc_mod, only: &
       & HydroBouEqSolverVDiffProc_Init, HydroBouEqSolverVDiffProc_Final, &
       & construct_vDiffProcMat, vDiffImplicitSolve => Solve

    use SpmlUtil_mod, only: &
         & xyz_wt, wt_xyz, wt_DSig_wt, g_Sig

    use dc_test

    ! Declaration statement
    !
    interface
       pure function get_anaSolution(x,t) result(anaSol)
         use dc_types, only: DP
         real(DP), intent(in) :: x(:)
         real(DP), intent(in) :: t
         real(DP) :: anaSol(size(x))
       end function get_anaSolution
    end interface
    real(DP), intent(in) :: eFoldTime
    character, intent(in) :: lBCType
    real(DP), intent(in) :: lBCVal
    character, intent(in) :: uBCType
    real(DP), intent(in) :: uBCVal
    integer, intent(in) :: tschemeId
    real(DP), intent(in) :: errorLimit
    real(DP), intent(in), optional :: dtRatio_

    ! Work variables
    !
    real(DP) :: dt
    real(DP) :: t
    integer :: nStep
    real(DP), allocatable :: xy_totDepth(:,:)
    real(DP), allocatable ::  xyz_Temp(:,:,:), wt_Temp(:,:), xyz_Work(:,:,:)
    integer :: i, j, n
    real(DP) :: l2Error
    real(DP), allocatable :: vDiffProcMat(:,:,:)
    integer, allocatable :: vDiffProcMatKp(:,:)
    character(STRING) :: message
    real(DP) :: alpha, beta
    real(DP) :: dtRatio

    ! Executable statement

    ! Preparation
    allocate(xy_totDepth(0:iMax-1,jMax))
    allocate(xyz_Temp(0:iMax-1,jMax,0:kMax), xyz_Work(0:iMax-1,jMax,0:kMax), wt_Temp(lMax,0:tMax))

    xy_totDepth = 1d0

    dtRatio = defaultDtRatio
    if(present(dtRatio_)) dtRatio = dtRatio_

    dt = eFoldTime*dtRatio
    nStep = EndNonDimTime*eFoldTime/dt
    call MessageNotify("M", PROGRAM_NAME, "* NStep=%d, e-foldTime=%f, dt=%f, tscheme=%d", &
         & i=(/nStep, tschemeId /), d=(/eFoldTime, dt/) )

    ! Determine coffiecients which depend on a scheme for temporal integration. 
    select case(tschemeId)
    case (TS_FWEuler) 
       alpha = 1d0; beta = 0d0;
    case (TS_BWEuler) 
       alpha = 0d0; beta = 2d0;
    case (TS_CRANK_NICO)
       alpha = 0.5d0; beta = 1d0;
    end select


    ! Construct a matrix used in convert (I - dt*D) operator.
    call HydroBouEqSolverVDiffProc_Init()
    call construct_vDiffProcMat(vDiffProcMat, vDiffProcMatKp, & !(out)
         & 1d0, beta*dt, xy_totDepth, uBCType, lBCType )

    ! Set initial condition
    forAll(i=0:iMax-1,j=1:jMax) xyz_Temp(i,j,:) = get_anaSolution(g_Sig,0d0)
    wt_Temp = wt_xyz(xyz_Temp)


    ! Loop for temporal integration
    !
    t = 0d0
    do n=1, nStep

       if(n==1 .or. mod(n,int(nStep*0.1)) == 0) then
          xyz_Temp = xyz_wt(wt_Temp)
          l2Error = sqrt(sum( (xyz_Temp(1,1,:) - get_anaSolution(g_Sig, t))**2 ))/real(kMax+1)

!!$          write(*,'(21f12.6)') xyz_Temp(1,1,0:kMax)
!!$          write(*,'(21f12.6)') get_anaSolution(g_Sig,t)

          message = CPrintf( &
               & "nondim time=%f, l2ErrorNorm=%f", d=(/ t/eFoldTime, l2Error /))
          call AssertLessThan(message=message, &
               & answer=errorLimit, check=l2Error )
       end if

       t = n*dt
       if(alpha == 0d0) then
          xyz_Work = xyz_wt(wt_Temp)
       else
          xyz_Work = xyz_wt( wt_Temp + alpha*dt*wt_DSig_wt(wt_DSig_wt(wt_Temp)) )
       end if

       ! If explicit temporal integration scheme is used(i.e. beta=0), 
       ! we call vDiffImplicitSolve to statisfy the boundary conditions.  
       xyz_Work(:,:,0) = uBCVal; xyz_Work(:,:,kMax) = lBCVal
       wt_Temp = vDiffImplicitSolve(vDiffProcMat, vDiffProcMatKp, xyz_Work)

    end do

    ! Finalization
    call HydroBouEqSolverVDiffProc_Final()

  end subroutine perform_timeIntegration


end program HydroBouEqSolverVDiffProc_mod_test
