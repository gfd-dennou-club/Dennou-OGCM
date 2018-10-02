!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Exp_APECoupleClimate_mod 
  
  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING
  
  use dc_message, only: &
       & MessageNotify

  use dc_iounit, only: &
       & FileOpen
  
  !* Dennou-OGCM

  use SpmlUtil_mod, only: &
       & AvrLonLat_xy    

  use UnitConversion_mod, only: &
       & degC2K
  
  use DOGCM_Admin_Constants_mod, only: &
       & PI, StB,                      &
       & Grav, RPlanet, Omega,         &
       & RefDens, RefSalt, LatentHeat, &
       & Cp0,                          &
       & UNDEFVAL

  use DSIce_Admin_Constants_mod, only: &
       & DensSnow, DensIce, IceMaskMin, LFreeze
  
  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, kMax, nMax, lMax,         &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,   &
       & xyz_Lat, xyz_Lon, xy_Topo,            &
       & DOGCM_Admin_Grid_UpdateVCoord, xyz_Z

  use CalculusDriver_mod, only: &
       & Calculus_IntBtmToTop
  
  use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar, &
       & DOGCM_IO_History_HistPut
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_Exp_Init, DOGCM_Exp_Final
  public :: DOGCM_Exp_SetInitCond
  public :: DOGCM_Exp_Do
  

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Exp_APECoupleClimate_mod' !< Module Name
  integer :: expCaseNum

  integer :: CurrentRunCycle
  integer :: CurrentRunType
  integer, parameter :: RUN_TYPE_COUPLED    = 1
  integer, parameter :: RUN_TYPE_STANDALONE = 2

  character(STRING) :: SfcBCDataDir
  real(DP) :: SfcBCMeanInitTime
  real(DP) :: SfcBCMeanEndTime

  character(STRING) :: RestartDataDir
  real(DP) :: RestartMeanInitTime
  real(DP) :: RestartMeanEndTime


  real(DP) :: OcnMeanDepth
  real(DP), parameter :: DEF_OCN_MEAN_DEPTH = 5.2d3

  real(DP), save :: OcnInitSalt
  real(DP), save :: OcnInitPTemp
  
  real(DP), parameter :: DensFreshWater = 1d3

  real(DP) :: InitMass

  character(*), parameter :: BudgetMonitorFileName = "budget_monitor"
  integer :: BudgetMonitorOutUnit
  
contains

  !>
  !!
  !!
  subroutine DOGCM_Exp_Init(configNmlFile)
    
    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile

    ! 作業変数
    ! Work variables
    !
    
    ! 実行文; Executable statements
    !

    call read_expConfig( &
         & CurrentRunCycle, CurrentRunType,                               & ! (out)
         & SfcBCDataDir, SfcBCMeanInitTime, SfcBCMeanEndTime,             & ! (out)
         & RestartDataDir, RestartMeanInitTime, RestartMeanEndTime,       & ! (out)
         & configNmlFile ) ! (in)


    if (CurrentRunType == RUN_TYPE_STANDALONE ) then
       call DOGCM_IO_History_RegistVar( 'SIceMass', 'T', 'global mean of sea ice mass', 'kg/m2')
    end if

    call FileOpen( BudgetMonitorOutUnit, BudgetMonitorFileName,  mode='w' )
    write(BudgetMonitorOutUnit,'(a)') &
         & "# time[sec] | SIceEn : DelSIceEn SIceNetHFlx SfcHFlxAI BtmHFlxIO | OcnEn : DelOcnEn SfcHFlxO"
    
  end subroutine DOGCM_Exp_Init

  !>
  !!
  !!
  subroutine DOGCM_Exp_Final

    ! 実行文; Executable statements
    !

    close( BudgetMonitorOutUnit )
    
  end subroutine DOGCM_Exp_Final

  !> @brief 
  !!
  !!
  subroutine DOGCM_Exp_SetInitCond()
    
    ! モジュール引用; Use statements
    !

    !* Dennou-OGCM/SIce
    
    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N, TIMELV_ID_B
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyza_U, xyza_V,              &
         & xyzaa_TRC, xya_SSH, xyza_H,  &
         & xya_SfcPres,                 &
         & TRCID_PTEMP, TRCID_SALT

    use DOGCM_Admin_Grid_mod, only: &
         & xy_Topo, xyz_Z, z_CK,    &
         & z_KAXIS_Weight

    use DOGCM_Boundary_vars_mod, only:     &
         & xy_SeaSfcTemp0, xy_SeaSfcSalt0,      &
         & xy_SfcHFlx0_ns, xy_SfcHFlx0_sr,      &
         & xy_DSfcHFlxDTsOcn => xy_DSfcHFlxDTs, &
         & xy_WindStressUOcn => xy_WindStressU, &
         & xy_WindStressVOcn => xy_WindStressV, &
         & xy_FreshWtFlxS0, xy_FreshWtFlx0
 
    use DSIce_Admin_Grid_mod, only: &
         & IAS => IA, JAS => JA, KAS => KA
    
    use DSIce_Admin_Variable_mod, only: &
         & xya_SIceCon, xya_IceThick, xya_SnowThick,   &
         & xya_SIceSfcTemp, xyza_SIceTemp, xyza_SIceEn
    
    use DSIce_Boundary_vars_mod, only:   &
         & xy_SUwRFlxSIce => xy_SUwRFlx, xy_SDwRFlxSIce => xy_SDwRFlx, &
         & xy_LUwRFlxSIce => xy_LUwRFlx, xy_LDwRFlxSIce => xy_LDwRFlx, &
         & xy_LatHFlxSIce => xy_LatHFlx, xy_SenHFlxSIce => xy_SenHFlx, &
         & xy_DLatSenHFlxDTsSIce => xy_DLatSenHFlxDTs,                 &
         & xy_SfcHFlxAI0_ns, xy_SfcHFlxAI0_sr, xy_SfcHFlxAO0,          &
         & xy_DSfcHFlxAODTs, xy_DSfcHFlxAIDTs,                         &
         & xy_SIceSfcTempIni => xy_SIceSfcTemp0,                       &
         & xy_WindStressUAI, xy_WindStressVAI,                         &
         & xy_RainFallSIce => xy_RainFall,                             &
         & xy_SnowFallSIce => xy_SnowFall,                             &
         & xy_EvapSIce => xy_Evap,                                     &
         & xy_EvapAO,                                                  &
         & xy_SIceSfcTemp4SfcHFlx => xy_SIceSfcTemp0 
    
    use SpmlUtil_mod, only: &
         & AvrLonLat_xy, calc_UVCosLat2VorDiv, calc_VorDiv2UV, xy_CosLat


    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: j
    integer :: k

    real(DP) :: xyz_U0(IA,JA,KA)
    real(DP) :: xyz_V0(IA,JA,KA)
    real(DP) :: xyz_SSH0(IA,JA,KA)
    real(DP) :: xyz_PTemp0(IA,JA,KA)
    real(DP) :: xyz_Salt0(IA,JA,KA)
    real(DP) :: xyz_H0(IA,JA,KA)
    real(DP) :: xy_SfcPres0(IA,JA)
    real(DP) :: xy_OcnMixLyrDepth(IA,JA)
    
    real(DP) :: xy_IceThick0(IAS,JAS)
    real(DP) :: xy_SnowThick0(IAS,JAS)
    real(DP) :: xyz_SIceTemp0(IAS,JAS,KAS)
    real(DP) :: xy_SIceSfcTemp0(IAS,JAS)
    real(DP) :: xy_SIceCon0(IAS,JAS)
    real(DP) :: xyz_SIceEn0(IAS,JAS,KAS)

    real(DP) :: xy_WindStressU(IA,JA)
    real(DP) :: xy_WindStressV(IA,JA)    
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    
    real(DP) :: xy_SDwRFlx(IA,JA)
    real(DP) :: xy_LDwRFlx(IA,JA)
    real(DP) :: xy_SUwRFlx(IA,JA)
    real(DP) :: xy_LUwRFlx(IA,JA)
    real(DP) :: xy_LatHFlx(IA,JA)
    real(DP) :: xy_SenHFlx(IA,JA)    
    real(DP) :: xy_DSfcHFlxDTs(IA,JA)
    
    real(DP) :: xy_RainFall(IA,JA)
    real(DP) :: xy_SnowFall(IA,JA)
    real(DP) :: xy_Evap(IA,JA)
    
    real(DP) :: InputWtMass
    real(DP) :: ModMassFac

    real(DP) :: xy_Tmp(IA,JA)
    
    real(DP), parameter :: Day2Sec = 86400d0

    integer :: unit00
    
    ! 実行文; Executable statement
    !

    !-------------------------------------------------------------
    
    
    call MessageNotify("M", module_name, "Set initial condition..")

    !-------------------------------------------------------------
    
    ! The depth from sea surface (topography)
    xy_Topo = OcnMeanDepth
    
    ! Sea surface height
    xya_SSH(:,:,:) = 0d0

    ! Update vertical coordinate
    call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_N),  & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                                 & ! (in)
         & )
    xyza_H(:,:,:,TIMELV_ID_B) = xyza_H(:,:,:,TIMELV_ID_N)

    !------------------------------------------------------------------------------
    ! Set initial condition 
    !------------------------------------------------------------------------------

    !* Set the inital value of ocean variables
    
    if ( CurrentRunType == RUN_TYPE_COUPLED .and. CurrentRunCycle == 1 ) then
       xyz_U0(:,:,:)     = 0d0
       xyz_V0(:,:,:)     = 0d0
       xyz_PTemp0(:,:,:) = OcnInitPTemp
       xyz_Salt0(:,:,:)  = OcnInitSalt
       xy_SfcPres0(:,:)  = 0d0
       xyz_H0(:,:,:)     = xyza_H(:,:,:,TIMELV_ID_N)
    else

       !
!!$       call get_Field4Standalone( SfcBCDataDir, "U.nc", "U",             &
!!$            & MeanInitTime, MeanEndTime, xyz=xyz_U0 )
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "V.nc", "V",             &
!!$            & MeanInitTime, MeanEndTime, xyz=xyz_V0 )
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "PTemp.nc", "PTemp",     &
!!$            & MeanInitTime, MeanEndTime, xyz=xyz_PTemp0 )
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "Salt.nc", "Salt",       &
!!$            & MeanInitTime, MeanEndTime, xyz=xyz_Salt0 )       
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "H.nc", "H",             &
!!$            & MeanInitTime, MeanEndTime, xyz=xyz_H0 )       
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "SfcPres.nc", "SfcPres", &
!!$            & MeanInitTime, MeanEndTime, xy=xy_SfcPres0 )       

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "U",            &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xyz=xyz_U0 )

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "V",            &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xyz=xyz_V0 )

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "PTemp",        &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xyz=xyz_PTemp0 )

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "Salt",         &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xyz=xyz_Salt0 )                 

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "H",            &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xyz=xyz_H0 )      

       call get_Field4Standalone( RestartDataDir, "RestartOcnData.nc", "SfcPres",      &
            & RestartMeanInitTime*Day2Sec, RestartMeanEndTime*Day2Sec, xy=xy_SfcPres0 )      

    end if

    !$omp parallel
    !$omp workshare
    xyza_U(:,:,:,TIMELV_ID_N) = xyz_U0(:,:,:)
    xyza_V(:,:,:,TIMELV_ID_N) = xyz_V0(:,:,:)
    xyzaa_TRC(:,:,:,TRCID_PTEMP,TIMELV_ID_N) = xyz_PTemp0(:,:,:)
    xyzaa_TRC(:,:,:,TRCID_SALT,TIMELV_ID_N) = xyz_Salt0(:,:,:)
    xya_SfcPres(:,:,TIMELV_ID_N) = xy_SfcPres0(:,:)
    !$omp end workshare
    !$omp end parallel

    xy_OcnMixLyrDepth(:,:) = xyz_H0(:,:,KS) * z_KAXIS_Weight(KS)
    
    !------------------------------------------------------------------

    !* Set the inital value of sea-ice variables
    
    if ( CurrentRunType == RUN_TYPE_COUPLED .and. CurrentRunCycle == 1 ) then
    else
!!$       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "IceThick",  &
!!$            & MeanInitTime, MeanEndTime, xy=xy_IceThick0 )
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "SnowThick", &
!!$            & MeanInitTime, MeanEndTime, xy=xy_SnowThick0 )
!!$
!!$       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "SIceCon",   &
!!$            & MeanInitTime, MeanEndTime, xy=xy_SIceCon0 )
!!$       
!!$       call get_IceFieldStandalone( SfcBCDataDir, "history_sice.nc",             & ! (in)
!!$            & MeanInitTime, MeanEndTime,                                         & ! (in)
!!$            & xy_IceThick0, xyz_SIceTemp0, xyz_SIceEn0, xy_SIceSfcTemp0,         & ! (out)
!!$            & xy_SIceCon0, xy_SnowThick0, xyz_PTemp0(:,:,KS), xyz_Salt0(:,:,KS), & ! (inout)
!!$            & xy_OcnMixLyrDepth                                                  & ! (in)
!!$            & )        
       
       call get_Field4Standalone( RestartDataDir, "RestartSIceData.nc", "IceThick",     &
            & RestartMeanEndTime*Day2Sec, RestartMeanEndTime*Day2Sec, xy=xy_IceThick0 )

       call get_Field4Standalone( RestartDataDir, "RestartSIceData.nc", "SnowThick",    &
            & RestartMeanEndTime*Day2Sec, RestartMeanEndTime*Day2Sec, xy=xy_SnowThick0 )

       call get_Field4Standalone( RestartDataDir, "RestartSIceData.nc", "SIceCon",      &
            & RestartMeanEndTime*Day2Sec, RestartMeanEndTime*Day2Sec, xy=xy_SIceCon0 )

       call get_IceFieldStandalone( RestartDataDir, "RestartSIceData.nc",               & ! (in)
            & RestartMeanEndTime*Day2Sec, RestartMeanEndTime*Day2Sec,                   & ! (in)
            & xy_IceThick0, xyz_SIceTemp0, xyz_SIceEn0, xy_SIceSfcTemp0,                & ! (out)
            & xy_SIceCon0, xy_SnowThick0, xyz_PTemp0(:,:,KS), xyz_Salt0(:,:,KS),        & ! (inout)
            & xy_OcnMixLyrDepth                                                         & ! (in)
            & )        

    
       !$omp parallel
       !$omp workshare
       xya_IceThick(:,:,TIMELV_ID_N)    = xy_IceThick0(:,:)
       xya_SnowThick(:,:,TIMELV_ID_N)   = xy_SnowThick0(:,:)
       xya_SIceCon(:,:,TIMELV_ID_N)     = xy_SIceCon0(:,:)
       xya_SIceSfcTemp(:,:,TIMELV_ID_N) = xy_SIceSfcTemp0(:,:)
       xyza_SIceTemp(:,:,:,TIMELV_ID_N) = xyz_SIceTemp0(:,:,:)
       xyza_SIceEn(:,:,:,TIMELV_ID_N)   = xyz_SIceEn0(:,:,:)
       !$omp end workshare
       !$omp end parallel       
    end if

    !----------------------------------------------------------------------------

    !* Set sea surface flux
    !

    if ( CurrentRunType == RUN_TYPE_COUPLED .and. CurrentRunCycle == 1 ) then
       xy_WindStressU(:,:) = 0d0
       xy_WindStressV(:,:) = 0d0

       xy_SDwRFlx(:,:) = 0d0
       xy_LDwRFlx(:,:) = 0d0
       xy_SUwRFlx(:,:) = 0d0
       xy_LUwRFlx(:,:) = 0d0
       xy_LatHFlx(:,:) = 0d0
       xy_SenHFlx(:,:) = 0d0

       xy_RainFall(:,:) = 0d0
       xy_SnowFall(:,:) = 0d0

       xy_DSfcHFlxDTs(:,:) = 0d0
    else
!!$    else if ( CurrentRunType == RUN_TYPE_STANDALONE ) then
       ! Input sea surface wind Stress
    
       call get_Field4Standalone( SfcBCDataDir, "a2o_WindStressX.nc", "a2o_WindStressX", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_WindStressU )
       call get_Field4Standalone( SfcBCDataDir, "a2o_WindStressY.nc", "a2o_WindStressY", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_WindStressV )


       ! Input sea surface heat flux and fresh water flux

       if ( CurrentRunType == RUN_TYPE_STANDALONE ) then

          call get_Field4Standalone( SfcBCDataDir, "a2o_SDwRFlx.nc", "a2o_SDwRFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SDwRFlx )
          call get_Field4Standalone( SfcBCDataDir, "a2o_LDwRFlx.nc", "a2o_LDwRFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_LDwRFlx )
          call get_Field4Standalone( SfcBCDataDir, "a2o_SUwRFlx.nc", "a2o_SUwRFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SUwRFlx )
          call get_Field4Standalone( SfcBCDataDir, "a2o_LUwRFlx.nc", "a2o_LUwRFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_LUwRFlx )
          call get_Field4Standalone( SfcBCDataDir, "a2o_LatHFlx.nc", "a2o_LatHFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_LatHFlx )
          call get_Field4Standalone( SfcBCDataDir, "a2o_SenHFlx.nc", "a2o_SenHFlx", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SenHFlx )
          
          xy_SfcHFlx0_ns = xy_LUwRFlx - xy_LDwRFlx + xy_LatHFlx + xy_SenHFlx
          xy_SfcHFlx0_sr = xy_SUwRFlx - xy_SDwRFlx               
       else
          call get_Field4Standalone( SfcBCDataDir, "SfcHFlx_ns.nc", "SfcHFlx_ns", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SfcHFlx0_ns )

          call get_Field4Standalone( SfcBCDataDir, "SfcHFlx_sr.nc", "SfcHFlx_sr", &
               & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SfcHFlx0_sr )       
       end if
       
       call get_Field4Standalone( SfcBCDataDir, "DSfcHFlxDTs.nc", "DSfcHFlxDTs", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_DSfcHFlxDTs )
       
       call get_Field4Standalone( SfcBCDataDir, "a2o_RainFall.nc", "a2o_RainFall", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_RainFall )

       call get_Field4Standalone( SfcBCDataDir, "a2o_SnowFall.nc", "a2o_SnowFall", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SnowFall )

       call get_Field4Standalone( SfcBCDataDir, "a2s_Evap.nc", "a2s_Evap", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_EvapSIce )

       call get_Field4Standalone( SfcBCDataDir, "a2o_Evap.nc", "a2o_Evap", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_EvapAO )
       
       call get_Field4Standalone( SfcBCDataDir, "PTemp.nc", "PTemp",   &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xyz=xyz_PTemp0 )

       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "SfcHFlxAI0_ns", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SfcHFlxAI0_ns )

       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "SfcHFlxAI0_sr", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SfcHFlxAI0_sr )
       
!!$       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "SfcHFlxAO", &
!!$            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_SfcHFlxAO0 )
!!$
       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "DSfcHFlxAODTs", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_DSfcHFlxAODTs )

       call get_Field4Standalone( SfcBCDataDir, "history_sice.nc", "DSfcHFlxAIDTs", &
            & SfcBCMeanInitTime, SfcBCMeanEndTime, xy=xy_DSfcHFlxAIDTs )
    end if

    xy_Evap(:,:)         = (1d0 - xy_SIceCon0)*xy_EvapAO + xy_SIceCon0*xy_EvapSIce  !xy_LatHFlx/LatentHeat
    xy_FreshWtFlxS0(:,:) = (xy_RainFall + xy_SnowFall - xy_Evap)/DensFreshWater       
    
    if ( CurrentRunType == RUN_TYPE_STANDALONE ) then

       ! Report the original input water mass from atmosphere to ocean and sea-ice components. 
       !
       call FileOpen( unit00, file='HydroCycInfo.txt', mode = 'w')
       write(unit00,*) "Atm -> Ocn, SIce [kg/(m2.s)]:", AvrLonLat_xy(xy_FreshWtFlxS0(IS:IE,JS:JE))*DensFreshWater
       write(unit00,*) "Rain             [kg/(m2.s)]:", AvrLonLat_xy( xy_RainFall(IS:IE,JS:JE) )
       write(unit00,*) "Snow             [kg/(m2.s)]:", AvrLonLat_xy( xy_SnowFall(IS:IE,JS:JE) )
       write(unit00,*) "Evap             [kg/(m2.s)]:", AvrLonLat_xy( xy_Evap(IS:IE,JS:JE) )
       close( unit00 )

       ! Modify the water budgets
       !

       call MessageNotify('M', module_name, 'Modify surface water budget..')
       InputWtMass = AvrLonLat_xy(xy_FreshWtFlxS0(IS:IE,JS:JE))*DensFreshWater
       if (InputWtMass > 0d0) then
          xy_EvapAO = xy_EvapAO + InputWtMass
          xy_EvapSIce = xy_EvapSIce + InputWtMass
       else
          xy_RainFall = xy_RainFall - InputWtMass
       end if

       InitMass = AvrLonLat_xy( &
            & DensSnow*xy_SnowThick0(IS:IE,JS:JE) +  DensIce*xy_IceThick0(IS:IE,JS:JE) )

       InitMass = InitMass - AvrLonLat_xy( Calculus_IntBtmToTop( &
            & (xyz_Salt0(IS:IE,JS:JE,:) - RefSalt)/RefSalt * 1d3,  xyz_H0(IS:IE,JS:JE,:) ) &
            & )

    end if
    
    !* Store input data 
    
    !$omp parallel
    !$omp workshare
    
    ! Ocean -----
    
    xy_SeaSfcTemp0 = xyz_PTemp0(:,:,KS)
    xy_SeaSfcSalt0 = xyz_Salt0(:,:,KS)

    xy_WindStressUOcn = xy_WindStressU
    xy_WindStressVOcn = xy_WindStressV

!!$    xy_SfcHFlx0_ns    = (xy_LUwRFlx - xy_LDwRFlx) + xy_LatHFlx + xy_SenHFlx
!!$    xy_SfcHFlx0_sr    = (xy_SUwRFlx - xy_SDwRFlx)
    xy_DSfcHFlxDTsOcn = xy_DSfcHFlxDTs !+ 4d0*StB*xy_SeaSfcTemp0**3
    xy_FreshWtFlxS0   = (xy_RainFall + xy_SnowFall - xy_EvapAO)/DensFreshWater
    xy_FreshWtFlx0    = xy_FreshWtFlxS0

    ! Sea ice ----

    xy_WindStressUAI = xy_WindStressU
    xy_WindStressVAI = xy_WindStressV
    xy_SUwRFlxSIce   = xy_SUwRFlx
    xy_SDwRFlxSIce   = xy_SDwRFlx
    xy_LUwRFlxSIce   = xy_LUwRFlx
    xy_LDwRFlxSIce   = xy_LDwRFlx
    xy_LatHFlxSIce   = xy_LatHFlx
    xy_SenHFlxSIce   = xy_SenHFlx
    xy_DLatSenHFlxDTsSIce = xy_DSfcHFlxDTs !+ 4d0*StB*xy_SeaSfcTemp0**3
    xy_RainFallSIce = xy_RainFall 
    xy_SnowFallSIce = xy_SnowFall
    xy_SfcHFlxAO0 = xy_SfcHFlx0_ns + xy_SfcHFlx0_sr
!!$    xy_EvapSIce     = xy_Evap

    where (xy_IceThick0 > 0d0 )
       xy_SIceSfcTemp4SfcHFlx =  degC2K( xy_SIceSfcTemp0 )
!!$            &    xy_SIceCon0*degC2K( xy_SIceSfcTemp0 ) &
!!$            & + (1d0 - xy_SIceCon0)*xy_SeaSfcTemp0
    elsewhere
       xy_SIceSfcTemp4SfcHFlx = xy_SeaSfcTemp0
    end where
     
    !$omp end workshare       
    !$omp end parallel

!!$    where(xy_IceThick0 <= 0d0)
!!$       xy_SnowFall = 0d0
!!$       xy_RainFall = 0d0
!!$       xy_EvapSIce = 0d0
!!$    end where
!!$    write(*,*) "Water budget:", AvrLonLat_xy(xy_FreshWtFlxS0(IS:IE,JS:JE))*DensFreshWater
!!$    write(*,*) "Snow:", AvrLonLat_xy(xy_SnowFall(IS:IE,JS:JE))
!!$    write(*,*) "Rain:", AvrLonLat_xy(xy_RainFall(IS:IE,JS:JE))
!!$    write(*,*) "Evap:", AvrLonLat_xy(xy_EvapSIce(IS:IE,JS:JE))
!!$    stop
!!$        
    
    call MessageNotify('M', module_name, 'The setting of initial condition has been finished.')

  end subroutine DOGCM_Exp_SetInitCond

  subroutine DOGCM_Exp_Do()

    ! モジュール引用; Use statements
    !
    use DOGCM_Admin_TInteg_mod, only: &
         & DelTime,                   &
         & TIMELV_ID_N, TIMELV_ID_B,  &
         & CurrentTime

    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC, TRCID_PTEMP, TRCID_SALT, xyza_H

    use DOGCM_IO_History_mod, only: &
         & DOGCM_IO_History_isOutputTiming
    
    use DOGCM_Boundary_vars_mod, only:  &
         & xy_FreshWtFlxS,                       &
         & xy_SfcHFlx_ns, xy_SfcHFlx_sr,         &
         & xy_SfcHFlxIO_ns, xy_SfcHFlxIO_sr,     &
         & xy_OcnSfcCellMask, OCNCELLMASK_OCEAN

    use DSIce_Boundary_vars_mod, only: &
         & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs, xy_DelSfcHFlxAI, &
         & xy_PenSWRFlxAI => xy_PenSDRFlx,                  &
         & xy_BtmHFlxIO, xy_SfcHFlxAO
    
    use DSIce_Admin_Variable_mod, only: &
         & xya_IceThick, xya_SnowThick,               &
         & xya_SIceCon, xyza_SIceEn, xya_SIceSfcTemp

    use DSIce_Admin_Grid_mod, only: &
         & ISS => IS, IES => IE, JSS => JS, JES => JE, &
         & KSS => KS, KES => KE

    use DSIce_Boundary_vars_mod, only:   &
         & xy_RainFall, xy_SnowFall, xy_Evap
    
    use SpmlUtil_mod, only: AvrLonLat_xy
    
    ! 局所変数
    ! Local variables
    !    

    real(DP) :: xy_Tmp(IA,JA)
    real(DP) :: xyz_Tmp(IA,JA,KA)

    real(DP) :: InputMass1
    real(DP) :: InputMass2
    real(DP) :: DelSIceMass
    real(DP) :: DelOcnMass
    real(DP) :: NowMass
    real(DP) :: avg
    
    integer :: TID_A
    integer :: TID_N

    integer :: i
    integer :: j

    real(DP) :: SIceEnA
    real(DP) :: DelSIceEn
    real(DP) :: OcnEnA
    real(DP) :: DelOcnEn
    real(DP) :: SfcHFlxAI
    real(DP) :: BtmHFlxIO
    real(DP) :: SfcHFlxO
    real(DP) :: SIceNetHFlx
    
    ! 実行文; Executable statement
    !

    !-------------------------------------------------------------
 
!!$    return

    if ( DOGCM_IO_History_isOutputTiming(CurrentTime-DelTime) ) then
       TID_A = TIMELV_ID_N
       TID_N = TIMELV_ID_B
       
!!$    xy_Tmp(:,:) = (xy_RainFall + xy_SnowFall - xy_Evap)*DelTime
!!$    call print_AverageInfo(InputMass1, "Input mass: ATM -> SIce,Ocn", xy=xy_Tmp)
!!$
!!$    xy_Tmp(:,:) = xy_FreshWtFlxS*DelTime*1d3
!!$    call print_AverageInfo(InputMass2, "Input mass: Atm,SIce -> Ocn", xy=xy_Tmp)
!!$
!!$    xy_Tmp(:,:) = &
!!$         &   DensIce*(xya_IceThick(:,:,TID_A) - xya_IceThick(:,:,TID_N))            &
!!$         & + DensSnow*(xya_SnowThick(:,:,TID_A) - xya_SnowThick(:,:,TID_N))
!!$    call print_AverageInfo(DelSIceMass, "SIce mass change:", xy=xy_Tmp)
!!$
!!$    xyz_Tmp(:,:,:) = &
!!$         & - (xyzaa_TRC(:,:,:,TRCID_SALT,TID_A) - xyzaa_TRC(:,:,:,TRCID_SALT,TID_N)) &
!!$         & / RefSalt * 1d3
!!$    call print_AverageInfo(DelOcnMass, "Ocn mass change:", xyz=xyz_Tmp)
!!$
!!$    NowMass = AvrLonLat_xy( &
!!$            & DensSnow*xya_SnowThick(IS:IE,JS:JE,TID_A) +  DensIce*xya_IceThick(IS:IE,JS:JE,TID_A) )
!!$    NowMass = NowMass - AvrLonLat_xy( &
!!$         & Calculus_IntBtmToTop( (xyzaa_TRC(IS:IE,JS:JE,:,TRCID_SALT,TID_A) - 35d0)/35d0*1d3, &
!!$         &                       xyza_H(IS:IE,JS:JE,:,TID_N) ))
!!$
!!$    
!!$    write(*,*) "Check net budget: ", &
!!$         & "(SIce)", (InputMass1-InputMass2) - DelSIceMass,   &
!!$         & "(Ocn)",  InputMass2 - DelOcnMass,                 &
!!$         & "(Tot)",  InputMass1 - (DelSIceMass + DelOcnMass), &
!!$         & "(IniMass, NowMass)", InitMass, NowMass
!!$    write(*,*) "--------------------------------------------------------------"
       

       SIceEnA = AvrLonLat_xy(  xya_SIceCon(ISS:IES,JSS:JES,TID_A)*(                    &
            &                   sum(xyza_SIceEn(ISS:IES,JSS:JES,KSS:KSS+1,TID_A), 3)    &
            &                 - DensSnow*LFreeze*xya_SnowThick(ISS:IES,JSS:JES,TID_A) ) &
            &    )

       DelSIceEn = SIceEnA - &
            &    AvrLonLat_xy(  xya_SIceCon(ISS:IES,JSS:JES,TID_N)*(                    &
            &                   sum(xyza_SIceEn(ISS:IES,JSS:JES,KSS:KSS+1,TID_N), 3)    &
            &                 - DensSnow*LFreeze*xya_SnowThick(ISS:IES,JSS:JES,TID_N) ) &
            &    )
              
!!$       where(xy_OcnSfcCellMask == OCNCELLMASK_OCEAN)
!!$          xy_Tmp = 0d0
!!$       elsewhere
!!$          xy_Tmp = - LFreeze*(xy_SnowFall)*DelTime
!!$       end where
!!$       call print_AverageInfo(InputMass2, "SnowEnOnSIce:", xy=xy_Tmp)
!!$
!!$       where(xy_OcnSfcCellMask == OCNCELLMASK_OCEAN)
!!$          xy_Tmp = 0d0
!!$       elsewhere
!!$          xy_Tmp = - LFreeze*( - xy_Evap)*DelTime
!!$       end where
!!$       call print_AverageInfo(InputMass2, "EvapEnOnSIce:", xy=xy_Tmp)

       where(xy_OcnSfcCellMask == OCNCELLMASK_OCEAN)
          xy_Tmp = 0d0
       elsewhere
          xy_Tmp = (xy_SfcHFlxAI + xy_DelSfcHFlxAI - xy_PenSWRFlxAI)*DelTime
       end where
       call print_AverageInfo(SfcHFlxAI, "", xy=xy_Tmp)

       where(xy_OcnSfcCellMask == OCNCELLMASK_OCEAN)
          xy_Tmp = xy_BtmHFlxIO - LFreeze*xy_SnowFall
       elsewhere
          xy_Tmp = xy_BtmHFlxIO
       end where
       call print_AverageInfo(BtmHFlxIO, "", xy=xy_Tmp*DelTime)

       where(xy_OcnSfcCellMask == OCNCELLMASK_OCEAN)
          xy_Tmp = (xy_BtmHFlxIO - LFreeze*xy_SnowFall)
       elsewhere
          xy_Tmp = - xya_SIceCon(:,:,TID_N)*(xy_SfcHFlxAI + xy_DelSfcHFlxAI - xy_PenSWRFlxAI)   &
               &   - xya_SIceCon(:,:,TID_N)*LFreeze*(xy_SnowFall - xy_Evap)                     &
               &   + xy_BtmHFlxIO                                                               &
               &   - (1d0 - xya_SIceCon(:,:,TID_N))*xy_SfcHFlxAO
       end where
       call print_AverageInfo(SIceNetHFlx, "", xy=xy_Tmp*DelTime)

!!$       xy_Tmp(:,:) = (xy_SfcHFlxIO_ns + xy_SfcHFlxIO_sr)*DelTime
!!$       call print_AverageInfo(InputMass2, "SfcHFlxIO:", xy=xy_Tmp)
       
       xy_Tmp(:,:) = (xy_SfcHFlx_ns + xy_SfcHFlx_sr)*DelTime
       call print_AverageInfo(SfcHFlxO, "", xy=xy_Tmp)

       OcnEnA = RefDens*Cp0*AvrLonLat_xy(Calculus_IntBtmToTop(                            &
            & xyzaa_TRC(IS:IE,JS:JE,:,TRCID_PTEMP,TID_A), xyza_H(IS:IE,JS:JE,:,TID_A) ))

       DelOcnEn = OcnEnA - RefDens*Cp0*AvrLonLat_xy(Calculus_IntBtmToTop(                &
            & xyzaa_TRC(IS:IE,JS:JE,:,TRCID_PTEMP,TID_N), xyza_H(IS:IE,JS:JE,:,TID_N) ))
       
!!$       write(*,*) "DelOcnEn:", DelOcnEn
!!$       write(*,*) "----------------------------------------------------"
       write(BudgetMonitorOutUnit, '(e15.6e2,a, e15.6e2,a,4e15.6e2,a, e15.6e2,a,2e15.6e2)') &
            & CurrentTime, " |",                                                 &
            & SIceEnA, " :", DelSIceEn, SIceNetHFlx, SfcHFlxAI, BtmHFlxIO, " |", &
            & OcnEnA,  " :", DelOcnEn, SfcHFlxO

!!$       write(*, '(e15.6e2,a, e15.6e2,a,4e15.6e2,a, e15.6e2,a,2e15.6e2)') &
!!$            & CurrentTime, " |",                                                 &
!!$            & SIceEnA, " :", DelSIceEn, SIceNetHFlx, SfcHFlxAI, BtmHFlxIO, " |", &
!!$            & OcnEnA,  " :", DelOcnEn, SfcHFlxO

    end if
    
  contains
    
    subroutine print_AverageInfo(avgVal, lbl, xyz, xy)

      real(DP), intent(out) :: avgVal
      character(*), intent(in) :: lbl
      real(DP), intent(in), optional :: xyz(IA,JA,KA)
      real(DP), intent(in), optional :: xy(IA,JA)

      if( present(xy) ) then
         avgVal = AvrLonLat_xy(xy(IS:IE,JS:JE))
      end if
      
      if( present(xyz) ) then
         avgVal = AvrLonLat_xy( Calculus_IntBtmToTop(xyz(IS:IE,JS:JE,:), xyza_H(IS:IE,JS:JE,:,TID_N)) )
      end if

      if (len(lbl)>0) then
         write(*,*) trim(lbl), avgVal
      end if
    end subroutine print_AverageInfo
    
  end subroutine DOGCM_Exp_Do  
  
  !------ private subroutines / functions -------------------------------------------------
  
  subroutine get_Field4Standalone( dir, ncFileName, varName, & ! (in)
       & MeanInitTime, MeanEndTime,                          & ! (in)
       & xy, xyz                                             & ! (out)
       & ) 

    ! モジュール引用; Use statements
    !    
    use dc_string, only: &
         & CPrintf
    use gtool_history, only: &
         & HistoryGet, HistoryGetPointer

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: dir
    character(*), intent(in) :: ncFileName
    character(*), intent(in) :: varName
    real(DP), intent(in) :: MeanInitTime
    real(DP), intent(in) :: MeanEndTime    
    real(DP), intent(out), optional :: xy(IA,JA)
    real(DP), intent(out), optional :: xyz(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !    
    character(STRING) :: ncFilePath
    
    real(DP), pointer :: timeSeries(:) => null()
    character(STRING) :: timeRangeStr
    integer :: n
    integer :: n1, n2
    integer, parameter :: Read2DChunkSize = 100
    integer, parameter :: Read3DChunkSize = 10
    integer :: ReadChunkSize
    
    integer :: TSeriesInitId
    integer :: TSeriesEndId
    integer :: TSeriesCutNum
    real(DP) :: xya_TMeanTmp(IA,JA,Read2DChunkSize)
    real(DP) :: xyza_TMeanTmp(IA,JA,KA,Read3DChunkSize)
    

    ! 実行文; Executable statement
    !
    
    call MessageNotify( 'M', module_name,  "get_Field4Standalone: varName=%a", ca=(/ varName /) )
    
    ncFilePath = trim(dir) // ncFileName

    call HistoryGetPointer(ncFilePath, "time", timeSeries)
!!$    write(*,*) "time: size=", size(timeSeries)
!!$    write(*,*) "      range=", timeSeries


    TSeriesInitId = -1
    TSeriesEndId  = -1   
    do n = 1, size(timeSeries)
       if (TSeriesInitId < 0 .and. MeanInitTime <= timeSeries(n)) then
          TSeriesInitId = n
       end if
       if (TSeriesEndId < 0 .and. MeanEndTime <= timeSeries(n)) then
          TSeriesEndId = n
          exit
       end if
       if (n == size(timeSeries)) then
          call MessageNotify( 'W', module_name, &
               & "The end of averaged time range is overwritten by %f.", d=(/ timeSeries(n) /) )
          TSeriesEndId = n
       end if
    end do

    TSeriesCutNum = TSeriesEndId - TSeriesInitId + 1

    if ( TSeriesInitId < 0 .or. TSeriesEndId < 0 ) then
       call MessageNotify( 'E', module_name, &
            & "Time range [%f:%f] over which the field is temporal averaged is invalid. Check!", &
            & d=(/ MeanInitTime, MeanEndTime /) )
    end if

    !--------------------------
    
    if (present(xy)) then
       xy(:,:) = 0d0
       ReadChunkSize = Read2DChunkSize
    else if (present(xyz)) then
       xyz(:,:,:) = 0d0
       ReadChunkSize = Read3DChunkSize
    end if
    
    n1 = TSeriesInitId
    do while( n1 <= TSeriesEndId )
       n2 = min(n1 + ReadChunkSize - 1, TSeriesEndId)
       
       if ( present(xy) ) then
          call HistoryGet( ncFilePath, varName, xya_TMeanTmp(IS:IE,JS:JE,1:n2-n1+1),   &
               & range=CPrintf( "time=%f:%f", d=(/ timeSeries(n1), timeSeries(n2) /) ) &
               & )

          xy(IS:IE,JS:JE) = xy(IS:IE,JS:JE) &
               & + sum(xya_TMeanTmp(IS:IE,JS:JE,1:n2-n1+1),dim=3)/dble(TSeriesCutNum)
       else if( present(xyz) ) then
          call HistoryGet( ncFilePath, varName, xyza_TMeanTmp(IS:IE,JS:JE,KS:KE,1:n2-n1+1),   &
               & range=CPrintf( "time=%f:%f", d=(/ timeSeries(n1), timeSeries(n2) /) )        &
               & )

          xyz(IS:IE,JS:JE,KS:KE) = xyz(IS:IE,JS:JE,KS:KE) &
               & + sum(xyza_TMeanTmp(IS:IE,JS:JE,KS:KE,1:n2-n1+1),dim=4) / dble(TSeriesCutNum)
       end if

       n1 = n2 + 1
    end do

    !--------------------------
    
    if ( associated(timeSeries) ) then
       deallocate( timeSeries )
    end if
    
  end subroutine get_Field4Standalone

  
  !-----------------------------------------------------------------------------------------------------------

  subroutine get_IceFieldStandalone( dir, ncFileName,               & ! (in)
       & MeanInitTime, MeanEndTime,                                 & ! (in)
       & xy_IceThick, xyz_SIceTemp, xyz_SIceEn, xy_SIceSfcTemp,     & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_SeaSfcTemp, xy_SeaSfcSalt,    & ! (inout)
       & xy_OcnMixLyrDepth                                          & ! (in)
       & ) 

    ! モジュール引用; Use statements
    !    
    use dc_string, only: &
         & CPrintf
    use gtool_history, only: &
         & HistoryGet, HistoryGetPointer

    use DSIce_Admin_Constants_mod, only: &
         & LFreeze, DensIce, SaltSeaIce, Mu, &
         & IceThickMin
    use DSIce_Admin_Grid_mod, only: &
         & IS, IE, IA, JS, JE, JA, KS, KE, KA         
    use DSIce_ThermoDyn_Winton2000_mod, only: &
         & calc_E_IceLyr1, calc_E_IceLyr2,      &
         & calc_Temp_IceLyr1, calc_Temp_IceLyr2

    use DOGCM_Admin_Constants_mod, only: &
         & RefDensOcn => RefDens, CpOcn => Cp0
    use DOGCM_Admin_Grid_mod, only: &
         & ISO=>IS, IEO=>IE, IAO=>IA, &
         & JSO=>JS, JEO=>JE, JAO=>JA

    use UnitConversion_mod, only: &
         & K2DegC

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: dir
    character(*), intent(in) :: ncFileName
    real(DP), intent(in) :: MeanInitTime
    real(DP), intent(in) :: MeanEndTime
    real(DP), intent(out) :: xy_IceThick(IA,JA)
    real(DP), intent(out) :: xyz_SIceTemp(IA,JA,KA)
    real(DP), intent(out) :: xyz_SIceEn(IA,JA,KA)
    real(DP), intent(out) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(inout) :: xy_SIceCon(IA,JA)
    real(DP), intent(inout) :: xy_SnowThick(IA,JA)
    real(DP), intent(inout) :: xy_SeaSfcTemp(IAO,JAO)
    real(DP), intent(inout) :: xy_SeaSfcSalt(IAO,JAO)
    real(DP), intent(in) :: xy_OcnMixLyrDepth(IAO,JAO)

    ! 局所変数
    ! Local variables
    !    
    character(STRING) :: ncFilePath
    
    real(DP), pointer :: timeSeries(:) => null()
    character(STRING) :: timeRangeStr
    integer :: n
    integer :: n1, n2
    integer, parameter :: ReadChunkSize = 100
    
    integer :: TSeriesInitId
    integer :: TSeriesEndId
    integer :: TSeriesCutNum
    real(DP) :: xya_IceThick(IA,JA,ReadChunkSize)
    real(DP) :: xya_SIceSfcTemp(IA,JA,ReadChunkSize)
    real(DP) :: xyza_SIceTemp(IA,JA,KA,ReadChunkSize)
    real(DP) :: xyza_SIceEn(IA,JA,KA,ReadChunkSize)

    integer :: i
    integer :: j
    real(DP) :: SIceFrzTemp
    real(DP) :: DelOcnPTemp
    real(DP) :: DelOcnSalt
    
    ! 実行文; Executable statement
    !
    
    call MessageNotify( 'M', module_name,  "get_Field4Standalone: varName=%a", ca=(/ "SIceTemp,SIceSfcTemp" /) )
    
    ncFilePath = trim(dir) // ncFileName
    
    call HistoryGetPointer(ncFilePath, "time", timeSeries)

    TSeriesInitId = -1
    TSeriesEndId  = -1    
    do n = 1, size(timeSeries)
       if (TSeriesInitId < 0 .and. MeanInitTime <= timeSeries(n)) then
          TSeriesInitId = n
       end if
       if (TSeriesEndId < 0 .and. MeanEndTime <= timeSeries(n)) then
          TSeriesEndId = n
          exit
       end if
       if (n == size(timeSeries)) then
          call MessageNotify( 'W', module_name, &
               & "The end of averaged time range is overwritten by %f.", d=(/ timeSeries(n) /) )
          TSeriesEndId = n
       end if
    end do

    TSeriesCutNum = TSeriesEndId - TSeriesInitId + 1

    if ( TSeriesInitId < 0 .or. TSeriesEndId < 0 ) then
       call MessageNotify( 'E', module_name, &
            & "Time range [%f:%f] over which the field is temporal averaged is invalid. Check!", &
            & d=(/ MeanInitTime, MeanEndTime /) )
    end if

    !--------------------------

    xy_IceThick(:,:)     = 0d0
    xyz_SIceEn(:,:,:)    = 0d0
    xy_SIceSfcTemp(:,:)  = 0d0
    
    n1 = TSeriesInitId
    do while( n1 <= TSeriesEndId )
       n2 = min(n1 + ReadChunkSize - 1, TSeriesEndId)

       call HistoryGet( ncFilePath, "IceThick", xya_IceThick(IS:IE,JS:JE,1:n2-n1+1),         &
            & range=CPrintf( "time=%f:%f", d=(/ timeSeries(n1), timeSeries(n2) /) ) )
       xy_IceThick(:,:) = xy_IceThick(:,:) + sum(xya_IceThick(:,:,1:n2-n1+1),dim=3)/dble(TSeriesCutNum)
       
       call HistoryGet( ncFilePath, "SIceSfcTemp", xya_SIceSfcTemp(IS:IE,JS:JE,1:n2-n1+1),   &
            & range=CPrintf( "time=%f:%f", d=(/ timeSeries(n1), timeSeries(n2) /) ) )

       call HistoryGet( ncFilePath, "SIceTemp", xyza_SIceTemp(IS:IE,JS:JE,KS:KE,1:n2-n1+1),  &
            & range=CPrintf( "time=%f:%f", d=(/ timeSeries(n1), timeSeries(n2) /) ) )

       do n = 1, n2-n1+1
          where ( xya_IceThick(:,:,n) > 0d0 ) 
             xyza_SIceEn(:,:,KS,n)   = 0.5d0*DensIce*xya_IceThick(:,:,n)  &
                  &                    * calc_E_IceLyr1(xyza_SIceTemp(:,:,KS,n), SaltSeaIce)
             xyza_SIceEn(:,:,KS+1,n) = 0.5d0*DensIce*xya_IceThick(:,:,n)  &
                  &                    * calc_E_IceLyr2(xyza_SIceTemp(:,:,KS+1,n), SaltSeaIce)
          elsewhere
             xyza_SIceEn(:,:,KS,n)   = 0d0
             xyza_SIceEn(:,:,KS+1,n) = 0d0
             xya_SIceSfcTemp(:,:,n)  = K2DegC(xy_SeaSfcTemp(:,:))
          end where
       end do
       
       xyz_SIceEn(:,:,:)    = xyz_SIceEn(:,:,:) + sum(xyza_SIceEn(:,:,:,1:n2-n1+1),dim=4)/dble(TSeriesCutNum)
       xy_SIceSfcTemp(:,:)  = xy_SIceSfcTemp(:,:) + sum(xya_SIceSfcTemp(:,:,1:n2-n1+1),dim=3)/dble(TSeriesCutNum)
       
       n1 = n2 + 1
    end do

    !--------------------------

    where ( xy_IceThick(:,:) > 0d0 ) 
       xyz_SIceTemp(:,:,KS)   = calc_Temp_IceLyr1( xyz_SIceEn(:,:,KS)/(0.5d0*DensIce*xy_IceThick(:,:)), SaltSeaIce )
       xyz_SIceTemp(:,:,KS+1) = calc_Temp_IceLyr2( xyz_SIceEn(:,:,KS+1)/(0.5d0*DensIce*xy_IceThick(:,:)), SaltSeaIce )
!:old       xy_SIceCon(:,:)        = 1d0
    elsewhere
       xyz_SIceTemp(:,:,KS  ) = UNDEFVAL
       xyz_SIceTemp(:,:,KS+1) = UNDEFVAL
       xyz_SIceEn(:,:,KS  )   = 0d0
       xyz_SIceEn(:,:,KS+1)   = 0d0
       xy_SIceSfcTemp(:,:)    = UNDEFVAL
    end where
    
    do j = JS, JE
       do i = IS, IE

          SIceFrzTemp = - Mu*SaltSeaIce
      
          if (    (  xy_IceThick(i,j) >= IceThickMin .and.                                           &
               &     (  xyz_SIceTemp(i,j,KS) >  0d0 .or. xyz_SIceTemp(i,j,KS+1) >  SIceFrzTemp       &
               &        .or. xy_SIceSfcTemp(i,j) > 0d0                                        )      &
               &  )                                                                                  &
               & .or. (xy_IceThick(i,j) > 0d0 .and. xy_IceThick(i,j) < IceThickMin)                  &
               & ) then

             call MessageNotify( 'M', module_name, "Melt ice and ice layer at grid point(%d,%d)..", &
                  & i=(/ i,j /) )
             write(*,*) "i,j=", i,j, ": SIceTemp=", xyz_SIceTemp(i,j,KS:KE), "SIceSfcTemp=", xy_SIceSfcTemp(i,j), &
                  & "SnowThick=", xy_SnowThick(i,j), &
                  & "IceThick=", xy_IceThick(i,j), "SIceEn=", xyz_SIceEn(i,j,KS:KE),                              &
                  & "SST=", xy_SeaSfcTemp(i,j), "SSS=", xy_SeaSfcSalt(i,j),                                       &
                  & "LFreeze", LFreeze*0.5d0*DensIce*xy_IceThick(i,j)
             
             DelOcnPTemp = - (   &
                  &              DensSnow * xy_SnowThick(i,j) * LFreeze       & ! > 0
                  &           - (xyz_SIceEn(i,j,KS) + xyz_SIceEn(i,j,KS+1))   & ! > 0
                  &           ) / ( CpOcn * RefDensOcn * xy_OcnMixLyrDepth(i,j) )

             DelOcnSalt = - RefSalt/xy_OcnMixLyrDepth(i,j)                                          &
                  &       * (DensSnow*xy_SnowThick(i,j) + DensIce*xy_IceThick(i,j))/DensFreshWater
             
             xy_SeaSfcTemp(i,j) = xy_SeaSfcTemp(i,j) + DelOcnPTemp
             xy_SeaSfcSalt(i,j) = xy_SeaSfcSalt(i,j) + DelOcnSalt
             write(*,*) "==> SST=", xy_SeaSfcTemp(i,j), ", SSS=", xy_SeaSfcSalt(i,j)

             xy_SIceSfcTemp(i,j)       = UNDEFVAL             
             xyz_SIceTemp(i,j,KS:KS+1) = UNDEFVAL
             xyz_SIceEn(i,j,KS:KS+1)   = 0d0
             xy_IceThick(i,j)          = 0d0
             xy_SIceCon(i,j)           = 0d0
             xy_SnowThick(i,j)         = 0d0

          end if
       end do
    end do

    !--------------------------
    
    if ( associated(timeSeries) ) then
       deallocate( timeSeries )
    end if
    
  end subroutine get_IceFieldStandalone

  !-----------------------------------------------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine read_expConfig( &
       & RunCycle, RunType,                                        & ! (out)
       & SfcBCDataDir, SfcBCMeanInitTime, SfcBCMeanEndTime,        & ! (out)
       & RestartDataDir, RestartMeanInitTime, RestartMeanEndTime,  & ! (out)
       & configNmlFileName )

    ! モジュール引用; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! 宣言文; Declaration statement
    !
    integer, intent(out) :: RunCycle
    integer, intent(out) :: RunType
    character(STRING), intent(out) :: SfcBCDataDir
    real(DP), intent(out) :: SfcBCMeanInitTime
    real(DP), intent(out) :: SfcBCMeanEndTime
    character(STRING), intent(out) :: RestartDataDir
    real(DP), intent(out) :: RestartMeanInitTime
    real(DP), intent(out) :: RestartMeanEndTime
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open
    character(TOKEN) :: RunTypeName
    
    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /Exp_APECoupleClimate_nml/ &
         & OcnMeanDepth,                                           &
         & OcnInitPTemp,                                           &
         & OcnInitSalt,                                            &
         & RunCycle, RunTypeName,                                  &
         & SfcBCDataDir, SfcBCMeanInitTime, SfcBCMeanEndTime,      &
         & RestartDataDir, RestartMeanInitTime, RestartMeanEndTime

    ! 実行文; Executable statement
    !

    OcnMeanDepth = DEF_OCN_MEAN_DEPTH
    OcnInitSalt  = 35d0
    OcnInitPTemp  = 280d0
    
    RunCycle     = 1
    RunTypeName  = "Coupled"
    SfcBCDataDir = "./SfcBC"
    SfcBCMeanInitTime = 0d0
    SfcBCMeanEndTime  = 0d0
    RestartDataDir = "./SfcBC"
    RestartMeanInitTime = 0d0
    RestartMeanEndTime  = 0d0

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &                   ! (out)
            & configNmlFileName, mode = 'r' )       ! (in)

       rewind( unit_nml )
       read( unit_nml, &                            ! (in)
            & nml = Exp_APECoupleClimate_nml, &     ! (out)
            & iostat = iostat_nml )                 ! (out)
       close( unit_nml )
    end if

    select case( RunTypeName )
    case ("Coupled")
       RunType = RUN_TYPE_COUPLED
    case ("Standalone")       
       RunType = RUN_TYPE_STANDALONE       
    case default
       call MessageNotify( 'E', module_name, "Unexpected run type '%a' is specified. Check!", &
            & ca=(/ RunTypeName /) )
    end select
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, "RunCycle=%d, RunType=%a", i=(/RunCycle/), ca=(/RunTypeName/))
    call MessageNotify( 'M', module_name, "SfcBCDataDir=%a, Averaged time range %f:%f", &
         & ca=(/SfcBCDataDir/), d=(/ SfcBCMeanInitTime, SfcBCMeanEndTime /) )
    call MessageNotify( 'M', module_name, "RestartDataDir=%a, Averaged time range %f:%f", &
         & ca=(/RestartDataDir/), d=(/ RestartMeanInitTime, RestartMeanEndTime /) )
    call MessageNotify( 'M', module_name, "OcnMeanDepth =%f [m]", d=(/ OcnMeanDepth /) )
    call MessageNotify( 'M', module_name, "OcnInitSalt  =%f [psu]", d=(/ OcnInitSalt /) )
    
  end subroutine read_expConfig

end module DOGCM_Exp_APECoupleClimate_mod

