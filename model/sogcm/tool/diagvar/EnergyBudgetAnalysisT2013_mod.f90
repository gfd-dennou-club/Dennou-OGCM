!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module EnergyBudgetAnalysisT2013_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  !
  
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify

  use gtool_history

  !* Dennou-OGCM
  !
  
  use Constants_mod, only: &
       & RPlanet, Grav, RefDens, &
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef, Cp0

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use EOSDriver_mod, only: &
       & EOSDriver_Eval
  
  use DiagVarFileSet_mod, only: &
       & gtool_historyauto_info

  use BoundaryCondO_mod, only: &
       & BoundaryCondO_Update

  !* diagvar
  !

  use DiagnoseUtil_mod
  use DiagVarEval_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: EnergyBudgetAnalysisT2013_Init, EnergyBudgetAnalysisT2013_Final
  public :: EnergyBudgetAnalysisT2013_Perform

  
  ! 公開変数
  ! Public variable
  !  
  character(*), parameter, public :: ENBUDGANAKEY_GLOBALMEANENERGY = 'GlobalMeanEnergy'
  character(*), parameter, public :: ENBUDGANAKEY_ENERGYBUDGET = 'EnergyBudget'

  character(*), parameter, public :: ENBUDGANAKEY_TEAVG = 'TEAvg'
  character(*), parameter, public :: ENBUDGANAKEY_KEAVG = 'KEAvg'
!  character(*), parameter, public :: ENBUDGANAKEY_PEAVG = 'PEAvg'
  character(*), parameter, public :: ENBUDGANAKEY_APEDAVG = 'APEDAvg'
  character(*), parameter, public :: ENBUDGANAKEY_PTEMPAVG = 'PTEMPAvg'
  character(*), parameter, public :: ENBUDGANAKEY_SALTAVG = 'SALTAvg'

  character(*), parameter, public :: ENBUDGANAKEY_DKEDTAVG    = 'dKEdtAvg'  
  character(*), parameter, public :: ENBUDGANAKEY_PE2KEAVG    = 'PE2KEAvg'
  character(*), parameter, public :: ENBUDGANAKEY_WFAVG       = 'WFAvg'
  character(*), parameter, public :: ENBUDGANAKEY_HDISPAVG   = 'HDISPAvg'
  character(*), parameter, public :: ENBUDGANAKEY_VDISPAVG   = 'VDISPAvg'

  character(*), parameter, public :: ENBUDGANAKEY_BFAVG    = 'BFAvg'  

  
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EnergyBudgetAnalysisT2013_mod' !< Module Name

  type(gt_history) :: hst_globalMeanEnergy
  type(gt_history) :: hst_energyBudget

  logical :: globalMeanEnergyFlag
  logical :: energyBudgetAnaFlag

  character(*), parameter :: RefDensProfile_NCFILE = 'RefDens4APED.nc'
  real(DP), allocatable :: z_RefDens(:), t_RefDens(:)


  integer, parameter :: N_GAUSSINTPT = 10
  real(DP), dimension(N_GAUSSINTPT), parameter :: X_GaussINTPt = &
       & (/ -1d0, -0.9195339081664588138289d0, -0.7387738651055050750031d0, &
       &          -0.4779249498104444956612d0, -0.1652789576663870246262d0, &
       &           0.1652789576663870246262d0,  0.4779249498104444956612d0, &
       &    0.7387738651055050750031d0, 0.9195339081664588138289d0, 1d0 /)

  real(DP), dimension(N_GAUSSINTPT), parameter :: Wt_GaussINTPt = &
       & (/ 0.02222222222222222222222d0, 0.1333059908510701111262d0, 0.2248893420631264521195d0, &
       &                                  0.292042683679683757876d0, 0.327539761183897456657d0,  &
       &                                  0.327539761183897456657d0, 0.292042683679683757876d0,  &
       &    0.2248893420631264521195d0,  0.1333059908510701111262d0, 0.02222222222222222222222d0 /)


  interface calc_globalMean
     module procedure calc_globalMean_xy
     module procedure calc_globalMean_xyz
  end interface calc_globalMean

contains

  !>
  !!
  !!
  subroutine EnergyBudgetAnalysisT2013_Init(diagVar_gthsInfo, GlMeanEnFlag, EnBudgeFlag)

    ! モジュール引用; Use statements
    !
    
    ! 宣言文; Declare statements
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    logical, intent(in) :: GlMeanEnFlag, EnBudgeFlag
    

    ! 実行文; Executable statements
    !

    
    globalMeanEnergyFlag = GlMeanEnFlag
    energyBudgetAnaFlag = EnBudgeFlag

    allocate(z_RefDens(0:kMax), t_RefDens(0:tMax))
!    call HistoryGet(RefDensProfile_NCFILE, 'RefDens', z_RefDens)    
    z_RefDens(:) = RefDens
    t_RefDens(:) = t_g(z_RefDens)
    
    call prepair_Output(diagVar_gthsInfo)

  end subroutine EnergyBudgetAnalysisT2013_Init

  !>
  !!
  !!
  subroutine EnergyBudgetAnalysisT2013_Final()


    ! 実行文; Executable statements
    !
    
    if( globalMeanEnergyFlag ) then
       call HistoryClose(hst_globalMeanEnergy)
    end if

    if( energyBudgetAnaFlag ) then
       call HistoryClose(hst_energyBudget)
    end if

    deallocate(z_RefDens, t_RefDens)
    
  end subroutine EnergyBudgetAnalysisT2013_Final


  !> @brief 
  !!
  !!
  subroutine EnergyBudgetAnalysisT2013_Perform()

    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_SaltN, xy_SurfHeightN, &
         & xy_totDepthBasic, z_PTempBasic

    use BoundaryCondO_mod, only: &
         & xy_WindStressU, xy_WindStressV
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: KE, APE, TE, PTEMP, SALT
    real(DP) :: WF, KE2PE, HDISP, VDISP, dKEdt

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PTemp
    real(DP), dimension(0:iMax-1,jMax) :: xy_totDepth
    
    integer :: k
    
    ! 実行文; Executable statement
    !

    if(.not. (globalMeanEnergyFlag.or.energyBudgetAnaFlag)) return

    forAll(k=0:kMax) &
         & xyz_PTemp(:,:,k) = z_PTempBasic(k) + xyz_PTempEddN(:,:,k)

    xy_totDepth(:,:) = xy_totDepthBasic + xy_SurfHeightN
    
    !
    !

    if (globalMeanEnergyFlag) then       
       call calc_globalMeanEngy(KE, APE, TE, &
            & xyz_UN, xyz_VN, xyz_PTemp, xyz_SaltN, xy_totDepth)

       PTEMP = calc_globalMean(xyz_PTemp, xy_totDepth)
       SALT = calc_globalMean(xyz_SaltN, xy_totDepth)
       
       call HistoryPut(ENBUDGANAKEY_KEAVG,    KE,    hst_globalMeanEnergy)
       call HistoryPut(ENBUDGANAKEY_APEDAVG,  APE,   hst_globalMeanEnergy)
       call HistoryPut(ENBUDGANAKEY_TEAVG,    TE,    hst_globalMeanEnergy)
       call HistoryPut(ENBUDGANAKEY_PTEMPAVG, PTEMP, hst_globalMeanEnergy)
       call HistoryPut(ENBUDGANAKEY_SALTAVG,  SALT,  hst_globalMeanEnergy)
    end if
    
    !
    !
    if (energyBudgetAnaFlag) then
       call analyze_EnergyBudget( &
            & WF, KE2PE, HDISP, VDISP, &
            & xyz_UN, xyz_VN, xyz_PTemp, xyz_SaltN, &
            & xy_WindStressU, xy_WindStressV, &
            & xy_totDepth )

       dKEdt = WF - KE2PE + HDISP + VDISP       

       ! Output
       !
       call HistoryPut(ENBUDGANAKEY_WFAVG,        WF, hst_energyBudget)
       call HistoryPut(ENBUDGANAKEY_PE2KEAVG, -KE2PE, hst_energyBudget)
       call HistoryPut(ENBUDGANAKEY_HDISPAVG,  HDISP, hst_energyBudget)
       call HistoryPut(ENBUDGANAKEY_VDISPAVG,  VDISP, hst_energyBudget)
       call HistoryPut(ENBUDGANAKEY_DKEDTAVG,  dKEdt, hst_energyBudget)
    end if


    
  end subroutine EnergyBudgetAnalysisT2013_Perform

  subroutine calc_globalMeanEngy(KE, APE, TE, &
       & xyz_U, xyz_V, xyz_PTemp, xyz_Salt, xy_totDepth )

    real(DP), intent(out) :: KE, APE, TE
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V, xyz_PTemp, xyz_Salt
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_totDepth

    integer :: i, j, k
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_DPE

    !$omp parallel do private(j,i)
    do k=0,kMax
       do j=1, jMax
          do i=0, iMax-1
             xyz_DPE(i,j,k) = calc_DPE( &
                  & xyz_PTemp(i,j,k), xyz_Salt(i,j,k), &
                  & g_Sig(k)*xy_totDepth(i,j),  xy_totDepth(i,j) &
                  & )
          end do
       end do
    end do
    
    KE  = calc_globalMean( RefDens*(xyz_U**2 + xyz_V**2), xy_totDepth )
    APE = calc_globalMean( RefDens*xyz_DPE, xy_totDepth )
    TE  = KE + APE
    
  end subroutine calc_globalMeanEngy

  function calc_globalInt_xyz(xyz, xy_totDepth) result(glInt)
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: glInt

    glInt = IntLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz)*xy_totDepth) * RPlanet**2
    
  end function calc_globalInt_xyz

  function calc_globalInt_xy(xy) result(glInt)
    real(DP), intent(in) :: xy(0:iMax-1,jMax)
    real(DP) :: glInt

    glInt = IntLonLat_xy(xy) * RPlanet**2
    
  end function calc_globalInt_xy

  function calc_globalMean_xyz(xyz, xy_totDepth) result(glMean)
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: glMean

    real(DP) :: globalVol

    globalVol = calc_globalInt_xy(xy_totDepth)
    glMean = calc_globalInt_xyz(xyz, xy_totDepth)/globalVol
    
  end function calc_globalMean_xyz

  function calc_globalMean_xy(xy) result(glMean)
    real(DP), intent(in) :: xy(0:iMax-1,jMax)
    real(DP) :: glMean
    real(DP) :: xy_Unit(0:iMax-1,jMax)

    xy_Unit = 1d0
    glMean = calc_globalInt_xy(xy)/calc_globalInt_xy(xy_Unit)
  end function calc_globalMean_xy
  
  function calc_DPE(PTemp, Salt, Depth, totDepth) result(DPE)
    real(DP), intent(in) :: PTemp, Salt, Depth, totDepth
    real(DP) :: DPE

    real(DP), parameter :: Zref = 0d0
     
    DPE = numIntegral_gaussLobattoOrd10(eval_Buoyancy, Depth, Zref)
    
  contains
    function eval_Buoyancy(a_z) result(a_Buoyancy)
      real(DP), intent(in) :: a_z(N_GAUSSINTPT)
      real(DP) :: a_Buoyancy(N_GAUSSINTPT)

      real(DP), dimension(N_GAUSSINTPT) :: a_RefDens, a_rhoEdd
      integer :: m
      
      call EOSDriver_Eval(rhoEdd=a_rhoEdd, &
           & theta= spread(PTemp,1,N_GAUSSINTPT),   &
           & S    = spread(Salt,1,N_GAUSSINTPT),    &
           & p    = -RefDens*Grav*a_z(:)            &
           & )

      do m=1, N_GAUSSINTPT 
         a_RefDens(m) = Interpolate_t(t_RefDens, a_z(m)/totDepth)
      end do

      a_buoyancy(:) = - Grav*( RefDens+a_rhoEdd - a_RefDens )/RefDens
    end function eval_Buoyancy

  end function calc_DPE

  function numIntegral_gaussLobattoOrd10(f, a, b) result(numInt)
    interface
       function f(x) result(ret)
         use dc_types, only: DP
         real(DP), intent(in) :: x(10)
         real(DP) :: ret(10)
       end function f
    end interface
    real(DP), intent(in) :: a, b

    real(DP) :: numInt

    numInt = 0.5d0*(b-a)*sum( &
         & Wt_GaussINTPt(:)*f( 0.5d0*((b - a)*X_GaussINTPt(:) + a + b) ) &
         & )
    
  end function numIntegral_gaussLobattoOrd10

  subroutine analyze_EnergyBudget( &
       & WF, KE2PE, HDISP, VDISP,  &
       & xyz_U, xyz_V, xyz_PTemp, xyz_Salt, &
       & xy_WindStressU, xy_WindStressV, &
       & xy_totDepth )

    use SpmlUtil_mod, only: g_Sig
    
    use BoundaryCondO_mod, only: &
         & xy_SurfHFlxO, xy_SurfFwFlxO
    
    real(DP), intent(out) :: WF, KE2PE, HDISP, VDISP
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V, xyz_PTemp, xyz_Salt
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_WindStressU, xy_WindStressV, &
         & xy_totDepth

    real(DP), dimension(lMax, 0:kMax) :: &
         & wz_Div, wz_Vor
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_Urf, xyz_Vrf, xyz_DensEdd, &
         & xyz_Div, xyz_SigDot, xyz_W, xyz_Buoyancy, &
         & xyz_CosLat, xyz_Work

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_HDiffU, xyz_HDiffV
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DzU, xyz_DzV, xyz_DzSalt, xyz_DzPTemp

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_Gptemp, xyz_Gsalt

    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_Wice, xy_SIceCon
    
    integer :: i, j, k
    real(DP) :: globalVol
    real(DP), parameter :: Zref = 0d0
    real(DP), parameter :: RefSalt  = 35d0
    real(DP) :: BF_PTemp, BF_Salt, VMix_PTemp, VMix_Salt
    
    xyz_CosLat(:,:,:) = cos(xyz_Lat)
    xyz_Urf(:,:,:) = xyz_U*xyz_CosLat; xyz_Vrf(:,:,:) = xyz_V*xyz_CosLat

    call wz_VectorCosLat2VorDiv( xyz_Urf, xyz_Vrf, & !(in) 
         & wz_Vor, wz_Div   )                        !(out)
    xyz_Div = xyz_wz(wz_Div)
    xyz_SigDot(:,:,:)  = Diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div )

    xyz_DensEdd(:,:,:) = eval_DensEdd(xyz_PTemp, xyz_Salt, xy_totDepth)
    do k=0, kMax
       xyz_Buoyancy(:,:,k) = - Grav*(RefDens + xyz_DensEdd(:,:,k) - z_RefDens(k))/RefDens
       xyz_W(:,:,k) = xy_totDepth(:,:)*xyz_SigDot(:,:,k)
    end do

    !
    !
    xyz_Work(:,:,:) = 1d0
    globalVol = calc_globalInt_xyz(xyz_Work, xy_totDepth)
    
    KE2PE = calc_globalInt_xyz( &
         &  - RefDens*(xyz_Buoyancy*xyz_W), &
         &  xy_totDepth ) / globalVol
    WF = calc_globalInt_xy( &
         & xyz_U(:,:,0)*xy_WindStressU + xyz_V(:,:,0)*xy_WindStressV &
         & ) / globalVol

    xyz_HDiffU(:,:,:) = hViscCoef*xyz_CosLat*xyz_AlphaOptr_wz(wz_Div, -wz_Vor)
    xyz_HDiffV(:,:,:) = hViscCoef*xyz_CosLat*xyz_AlphaOptr_wz(wz_Vor,  wz_Div)    
    HDISP = calc_globalInt_xyz( &
         & RefDens*(xyz_U*xyz_HDiffU + xyz_V*xyz_HDiffV), &
         & xy_totDepth ) / globalVol

    xyz_DzU(:,:,:) = xyz_DSig_xyz(xyz_U)/spread(xy_totDepth, 3, kMax+1)
    xyz_DzV(:,:,:) = xyz_DSig_xyz(xyz_V)/spread(xy_totDepth, 3, kMax+1)
    VDISP = calc_globalInt_xyz( &
         & -RefDens*(vViscCoef*(xyz_DzU**2 + xyz_DzV**2)), &
         & xy_totDepth ) / globalVol

    write(*,*) "WF     =", WF
    write(*,*) "HDISP  =", HDISP
    write(*,*) "VDISP1 =", VDISP
    write(*,*) "PE2KE  =", -KE2PE
    write(*,*) "----------------------------"
    write(*,*) "net    =", WF + HDISP + VDISP- KE2PE 

    !
    !
    xy_SIceCon(:,:) = 0d0; xy_Wice = 0d0
    call BoundaryCondO_Update( xy_SIceCon, xy_Wice )

    !$omp parallel do private(j,i)
    do k=0, kMax
       Do j=1, jMax
          do i=0, iMax-1
             call calc_EAPE_coef( xyz_Gptemp(i,j,k), xyz_Gsalt(i,j,k), & ! (out)
                  & xyz_PTemp(i,j,k), xyz_Salt(i,j,k),                 & ! (in)
                  & g_Sig(k)*xy_totDepth(i,j), Zref                    & ! (in)
                  & )
          end do
       end do
    end do


    BF_PTemp = calc_globalInt_xy( &
         & RefDens*xyz_Gptemp(:,:,0)*xy_SurfHFlxO/(RefDens*Cp0)  &
         & )/globalVol
         
    BF_Salt = calc_globalInt_xy( &
         & RefDens*xyz_Gsalt(:,:,0)*xy_SurfFwFlxO*RefSalt &
         & )/globalVol

    xyz_DzPTemp(:,:,:) = xyz_DSig_xyz(xyz_PTemp)/spread(xy_totDepth, 3, kMax+1)
    xyz_DzSalt(:,:,:) = xyz_DSig_xyz(xyz_Salt)/spread(xy_totDepth, 3, kMax+1)

    VMix_PTemp = calc_globalInt_xyz( &
         & RefDens*vDiffCoef*xyz_DSig_xyz(xyz_PTemp)*xyz_DSig_xyz(xyz_Gptemp)/spread(xy_totDepth**2, 3, kMax+1), &
         & xy_totDepth )/globalVol
    VMix_Salt = calc_globalInt_xyz( &
         & RefDens*vDiffCoef*xyz_DSig_xyz(xyz_Salt)*xyz_DSig_xyz(xyz_Gsalt)/spread(xy_totDepth**2, 3, kMax+1), &
         & xy_totDepth )/ globalVol
    
    write(*,*) "*****************************"
    write(*,*) "BF_PTemp   =", BF_PTemp
    write(*,*) "BF_Salt    =", BF_Salt
    write(*,*) "VMix_PTemp =", VMix_PTemp
    write(*,*) "VMix_Salt  =", VMix_Salt
    write(*,*) "KE2PE      =", KE2PE
    write(*,*) "----------------------------"
    write(*,*) "net        =", &
         & BF_PTemp + BF_Salt + VMix_PTemp + VMix_Salt + KE2PE
    
  end subroutine analyze_EnergyBudget
       
  subroutine calc_EAPE_coef( G_PTemp, G_Salt, &
       & PTemp, Salt, Depth, Zref )
    real(DP), intent(out) :: G_PTemp, G_Salt
    real(DP), intent(in) :: PTemp, Salt, Depth, Zref

    G_PTemp = numIntegral_gaussLobattoOrd10(eval_dBuoyancy_dPTemp, Depth, Zref)
    G_Salt = numIntegral_gaussLobattoOrd10(eval_dBuoyancy_dSalt, Depth, Zref)
    
  contains
    function eval_dBuoyancy_dPTemp(a_z) result(a_dBuoyancy_dPTemp)
      real(DP), intent(in) :: a_z(N_GAUSSINTPT)
      real(DP) :: a_dBuoyancy_dPTemp(N_GAUSSINTPT)

      real(DP) :: aa_rhoEdd(10,2)
      real(DP), parameter :: dPTemp = 1d-3

      call EOSDriver_Eval(rhoEdd=aa_rhoEdd(:,1),     &
           & theta= spread(PTemp - dPTemp, 1, 10),   &
           & S    = spread(Salt, 1, 10),             &
           & p    = -RefDens*Grav*a_z(:)             &
           & )
      call EOSDriver_Eval(rhoEdd=aa_rhoEdd(:,2),     &
           & theta= spread(PTemp + dPTemp, 1, 10),   &
           & S    = spread(Salt, 1, 10),             &
           & p    = -RefDens*Grav*a_z(:)             &
           & )

      a_dBuoyancy_dPTemp(:) = - Grav*(aa_rhoEdd(:,2) - aa_rhoEdd(:,1)) / RefDens /(2d0*dPTemp) 
      
    end function eval_dBuoyancy_dPTemp

    function eval_dBuoyancy_dSalt(a_z) result(a_dBuoyancy_dSalt)
      real(DP), intent(in) :: a_z(N_GAUSSINTPT)
      real(DP) :: a_dBuoyancy_dSalt(N_GAUSSINTPT)

      real(DP) :: aa_rhoEdd(10,2)
      real(DP), parameter :: dSalt = 1d-4

      call EOSDriver_Eval(rhoEdd=aa_rhoEdd(:,1),     &
           & theta= spread(PTemp, 1, 10),            &
           & S    = spread(Salt - dSalt, 1, 10),     &
           & p    = -RefDens*Grav*a_z(:)             &
           & )
      call EOSDriver_Eval(rhoEdd=aa_rhoEdd(:,2),     &
           & theta= spread(PTemp, 1, 10),            &
           & S    = spread(Salt + dSalt, 1, 10),     &
           & p    = -RefDens*Grav*a_z(:)             &
           & )

      a_dBuoyancy_dSalt(:) = - Grav*(aa_rhoEdd(:,2) - aa_rhoEdd(:,1)) / RefDens /(2d0*dSalt)
      
    end function eval_dBuoyancy_dSalt
    
  end subroutine calc_EAPE_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine prepair_Output(diagVar_gthsInfo)
    
    ! 宣言文; Declaration statement
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    
    ! 局所変数
    ! Local variables
    !

    
    ! 実行文; Executable statement
    !

    if ( globalMeanEnergyFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'GlobalMeanEnergy.nc', title='global mean of each energy', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_globalMeanEnergy )  

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_TEAVG, dims=(/'t'/), &
             & longname='global mean of total energy', units='J/m3', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_KEAVG, dims=(/'t'/), &
             & longname='global mean of kinetic energy', units='J/m3', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_APEDAVG, dims=(/'t'/), &
             & longname='global mean of available potentail energy density', units='J/m3', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_PTEMPAVG, dims=(/'t'/), &
             & longname='global mean of potential temerature', units='K', xtype='double',&
             & history=hst_globalMeanEnergy)

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_SALTAVG, dims=(/'t'/), &
             & longname='global mean of salinity', units='psu', xtype='double',&
             & history=hst_globalMeanEnergy)

     end if

    if ( energyBudgetAnaFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'EnergyBudget.nc', title='energy budget analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=(/'t'/), dimsizes=(/ 0 /), &
            & longnames=(/'time'/),&
            & units=(/ diagVar_gthsInfo%intUnit /), &
            & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_energyBudget )

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_WFAVG, dims=(/'t'/), &
             & longname='global mean of the term which relate to KE input due to wind stress.', &
             & units='W/m3', xtype='double', history=hst_energyBudget )
       
        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_PE2KEAVG, dims=(/'t'/), &
             & longname='global mean of the term which relate to energy conversion between KE and APED', &
             & units='W/m3', xtype='double', history=hst_energyBudget )

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_HDISPAVG, dims=(/'t'/), &
             & longname='global mean of the term associated with (horizonatal) viscous dissipation.', & 
             & units='W/m3', xtype='double', history=hst_energyBudget )

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_VDISPAVG, dims=(/'t'/), &
             & longname='global mean of the term associated with (vertical) viscous dissipation.', & 
             & units='W/m3', xtype='double', history=hst_energyBudget )

        call HistoryAddVariable( & 
             & varname=ENBUDGANAKEY_DKEDTAVG, dims=(/'t'/), &
             & longname='global mean of KE tendency (dKE/dt).', & 
             & units='W/m3', xtype='double', history=hst_energyBudget )
        
    end if

  end subroutine prepair_Output

end module EnergyBudgetAnalysisT2013_mod

