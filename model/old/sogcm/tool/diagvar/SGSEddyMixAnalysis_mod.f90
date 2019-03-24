!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSEddyMixAnalysis_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  use Constants_mod, only: &
       & PI, RPlanet, Grav, Omega, RefDens, Cp0, &
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use DiagVarFileSet_mod
  use DiagnoseUtil_mod
  use DiagVarEval_mod

  use SGSEddyMixingHelper_mod, only: &
       & prepare_SlopeTapering


  use SGSEddyMixing_mod, only: &
       & SGSEddyMixing_Init, SGSEddyMixing_Final, &
       & SGSEddyMixing_GetParameters, &
       & calc_IsoNeutralSlope, calc_BolusVelocity, &
       & calc_IsopycDiffFlux, calc_SkewFlux

  use SGSEddyMixingHelper_mod, only: &
       & xyz_Dz_xyz

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: SGSEddyMixAnalysis_Init, SGSEddyMixAnalysis_Final
  public :: SGSEddyMixAnalysis_perform

  character(*), parameter, public :: SGSEDDANAKEY_BOLUSU = 'BolusU'
  character(*), parameter, public :: SGSEDDANAKEY_BOLUSV = 'BolusV'
  character(*), parameter, public :: SGSEDDANAKEY_BOLUSW = 'BolusW'
  character(*), parameter, public :: SGSEDDANAKEY_BOLUSMSTREAMFUNC = 'BolusMStreamFunc'
  character(*), parameter, public :: SGSEDDANAKEY_RESMSTREAMFUNC = 'ResMStreamFunc'


  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_VDIFF = 'PT_Tend_VDiff'  
  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_HDIFF = 'PT_Tend_HDiff'  
  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_IPDIFF = 'PT_Tend_IPDiff' ! Isopycnal diffusion
  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_GM     = 'PT_Tend_GM'     ! GM
  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_ADV    = 'PT_Tend_ADV'    ! Advection
  character(*), parameter, public :: SGSEDDANAKEY_PT_TEND_LTD    = 'PT_Tend_LTD'    ! Local time derivative

  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_VDIFF = 'PT_TendVINT_VDiff'
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_HDIFF    = 'PT_TendVINT_HDiff'    ! Horizontal (hyper) diffusion
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_IPDIFF = 'PT_TendVINT_IPDiff'     ! Isopycnal diffusion
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_GM     = 'PT_TendVINT_GM'     ! GM 
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_ADV    = 'PT_TendVINT_ADV'    ! Advection
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_LTD    = 'PT_TendVINT_LTD'    ! Local time derivative
  character(*), parameter, public :: SGSEDDANAKEY_PT_TENDVINT_TOT    = 'PT_TendVINT_TOT'    ! Total (ADV + IPDiff + GM + HDIFF + VDIFF)
  
!!$  character(*), parameter, public :: SGSEDDANAKEY_SALTTEND_IPDIFF  = 'SaltTend_IPDiff'
!!$  character(*), parameter, public :: SGSEDDANAKEY_SALTTEND_GM      = 'SaltTend_GM'  

  character(*), parameter, public :: SGSEDDANAKEY_ENGYFLXLAT_TOTAL  = "EngyFlxLatTot"
  character(*), parameter, public :: SGSEDDANAKEY_ENGYFLXLAT_TOTAL_HDIV  = "EngyFlxLatTotHDiv"
  character(*), parameter, public :: SGSEDDANAKEY_ENGYFLXLAT_EULER  = "EngyFlxLatEuler"
  character(*), parameter, public :: SGSEDDANAKEY_ENGYFLXLAT_BOLUS  = "EngyFlxLatBolus"
  character(*), parameter, public :: SGSEDDANAKEY_ENGYFLXLAT_IPDiff = "EngyFlxLatIPDiff"
  

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSEddyMixAnalysis_mod' !< Module Name

  type(gt_history), save :: hst_SGSEddMix
  logical :: SGSEddyMixAnaFlag

  logical :: BolusVelAnaFlag
  logical :: TendencyAnaFlag 
  logical :: TendencyVIntAnaFlag 
  logical :: EngyFlxLatAnaFlag

  real(DP), allocatable :: xyz_PTempB(:,:,:)
  logical :: FirstStepFlag
  
contains

  !>
  !!
  !!
  subroutine SGSEddyMixAnalysis_Init( &
       configNmlFile, diagVar_gthsInfo, isActivated )

    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    logical, intent(in) :: isActivated

    ! 作業変数
    ! Work variables
    !
    integer :: n


    ! 実行文; Executable statements
    !

    SGSEddyMixAnaFlag = isActivated
    
    if(SGSEddyMixAnaFlag) then
       call SGSEddyMixing_Init(configNmlFileName=configNmlFile)

       BolusVelAnaFlag = .true.
       EngyFlxLatAnaFlag = .true.
       TendencyAnaFlag = .false.
       TendencyVIntAnaFlag = .true.
       
       call prepair_Output(diagVar_gthsInfo)

       FirstStepFlag = .true.
    end if
    
  end subroutine SGSEddyMixAnalysis_Init

  !>
  !!
  !!
  subroutine SGSEddyMixAnalysis_Final()

    ! 実行文; Executable statements
    !

    if( SGSEddyMixAnaFlag ) then
       deallocate( xyz_PTempB )
    end if
    
  end subroutine SGSEddyMixAnalysis_Final

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixAnalysis_perform( &
       xyz_U, xyz_V, xyz_SigDot, xyz_PTemp, xyz_Salt, xy_totDepth )
    
    ! モジュール引用 ; Use statements
    !

    use SGSEddyMixingHelper_mod, only: &
         & prepare_SlopeTapering, xyz_Dz_xyz
    use SGSEddyMixing_mod, only: &
         & calc_IsoNeutralSlope, calc_BolusVelocity

    use HydroBoudEq_TimeInteg_v2_mod, only: &
         calc_vViscDiffCoef_MixLyrSimple
    
    use wa_zonal_module_sjpack
    
    ! 宣言文; Declaration statement
    !

    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)    
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)

    ! 局所変数
    ! Local variables
    !
    integer :: k

    real(DP) :: xyz_Depth(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_totDepth(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DensPot(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_SLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SLat(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_T(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_BolusU(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_BolusV(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_BolusW(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_FlxLon(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_FlxLat(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_FlxSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Tendency(0:iMax-1,jMax,0:kMax)
    
    real(DP) :: xy_HFlxTmp(0:iMax-1,jMax)
    real(DP) :: xy_HFlxTot(0:iMax-1,jMax)
    real(DP) :: xy_VFlxTot(0:iMax-1,jMax)
    real(DP) :: xy_TendencyVInt(0:iMax-1,jMax)
    real(DP) :: xy_TendTotVInt(0:iMax-1,jMax)
    
    real(DP) :: xyz_CosLat(0:iMax-1,jMax,0:kMax)

    real(DP), parameter :: Fact_W2PW = 1d-15


    real(DP) :: w_LaplaEigVal(lMax)
    real(DP) :: w_HDifCoefH(lMax)
    integer, parameter :: HDOrder = 8
    real(DP), parameter :: HDEFoldTimeSec = 5*86400d0
    real(DP) :: VisCoef

    real(DP) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_VViscCoef(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_VDiffCoef(0:iMax-1,jMax,0:kMax)
    
    ! 実行文; Executable statement
    !
    
    if(.not. SGSEddyMixAnaFlag) return


    call MessageNotify('M', module_name, 'Perform the analysis of SGS eddy mixing parameterization .')    

    !* Preparation 
    !

    xyz_CosLat = cos(xyz_Lat)
    !$omp parallel do
    do k=0, kMax
       xyz_totDepth(:,:,k) = xy_totDepth
       xyz_Depth(:,:,k) = xy_totDepth*g_Sig(k)
    end do
    

    ! Calculate potential density, and the slope vector. 

    xyz_DensPot(:,:,:) = eval_DensPot(xyz_PTemp, xyz_Salt, PressRef=0d0)
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
           & xyz_DensPot, xyz_totDepth ) !(in)

    ! Calculate tapering function used to avoid numerical instability. 
    call prepare_SlopeTapering(xyz_T, xyz_SLon, xyz_SLat, & ! (inout)
         & xyz_DensPot, xyz_Depth                         & ! (in)
         & )

!!$    call prepare_DFM08Info( &
!!$         & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth, &
!!$         & xyz_G )


    
    !* Calculate bolus velocities and the corresponding stream function, and
    !  output them. 
    !

    if ( BolusVelAnaFlag ) then
       call calc_BolusVelocity(xyz_BolusU, xyz_BolusV, xyz_BolusW, & ! (out)
            & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T )                 ! (in)    
       call HistoryPut(SGSEDDANAKEY_BOLUSU, xyz_BolusU, hst_SGSEddMix)
       call HistoryPut(SGSEDDANAKEY_BOLUSV, xyz_BolusV, hst_SGSEddMix)
       call HistoryPut(SGSEDDANAKEY_BOLUSW, xyz_BolusW, hst_SGSEddMix)

       call HistoryPut(SGSEDDANAKEY_BOLUSMSTREAMFUNC, &
            & eval_MassStreamFunc(xyz_BolusV, xy_totDepth), hst_SGSEddMix )
       call HistoryPut(SGSEDDANAKEY_RESMSTREAMFUNC,   &
            & eval_MassStreamFunc(xyz_V + xyz_BolusV, xy_totDepth), hst_SGSEddMix)

!!$
!!$       write(*,*) "Sum BolusV=", AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyz_BolusV) )
    end if

    !* Calculte tendencies of potential temperature and salinity due to isopycnal diffusion and GM,
    !  and diagnose the meridional energy heat fluxes. Then they are output.

    ! For potential temperature    
    xy_HFlxTot = 0d0
    xy_VFlxTot = 0d0
    xy_TendTotVInt = 0d0

    if ( TendencyAnaFlag .or. TendencyVIntAnaFlag .or. EngyFlxLatAnaFlag ) then

       !* Isopycnal diffusive flux
       !
       
       call calc_IsopycDiffFlux( xyz_FlxLon, xyz_FlxLat, xyz_FlxSig, &
            xyz_PTemp, wz_xyz(xyz_PTemp), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )

       if ( TendencyAnaFlag .or. TendencyVIntAnaFlag ) then
          call calc_DivFlux( xyz_Tendency, xyz_FlxLon, xyz_FlxLat, xyz_FlxSig )
          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_IPDIFF, xyz_Tendency, hst_SGSEddMix )
          if ( TendencyVIntAnaFlag ) then
             xy_TendencyVInt = xy_IntSig_BtmToTop_xyz(xyz_Tendency)
             xy_TendTotVInt = xy_TendTotVInt + xy_TendencyVInt
             call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_IPDIFF, xy_TendencyVInt, hst_SGSEddMix )
          end if
       end if

       if ( EngyFlxLatAnaFlag ) then
          ! Note that minus sign is needed because xyz_FlxLat = K \nabla PTemp. 
          xy_HFlxTmp = - xy_IntSig_BtmToTop_xyz(xyz_FlxLat) * xy_totDepth 
          xy_HFlxTot = xy_HFlxTot + xy_HFlxTmp
          xy_VFlxTot = xy_VFlxTot + xy_IntSig_BtmToTop_xyz(xyz_Dz_xyz(-xyz_FlxSig)) * xy_totDepth
          call HistoryPut( SGSEDDANAKEY_ENGYFLXLAT_IPDiff,                                   &
               Fact_W2PW * RefDens*Cp0 * RPlanet*y_IntLon_xy(xy_HFlxTmp*xyz_CosLat(:,:,0)),  &
               hst_SGSEddMix )
       end if

       !* GM flux
       !
       
       call calc_SkewFlux( xyz_FlxLon, xyz_FlxLat, xyz_FlxSig, &
            xyz_PTemp, wz_xyz(xyz_PTemp), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )

       if ( TendencyAnaFlag .or. TendencyVIntAnaFlag ) then
          call calc_DivFlux( xyz_Tendency, xyz_FlxLon, xyz_FlxLat, xyz_FlxSig )
          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_GM, xyz_Tendency, hst_SGSEddMix )
          if ( TendencyVIntAnaFlag ) then
             xy_TendencyVInt = xy_IntSig_BtmToTop_xyz(xyz_Tendency)
             xy_TendTotVInt = xy_TendTotVInt + xy_TendencyVInt
             call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_GM, xy_TendencyVInt, hst_SGSEddMix )
          end if
       end if

       if ( EngyFlxLatAnaFlag ) then
          ! Note that minus sign is needed because xyz_FlxLat = K \nabla PTemp. 
          xy_HFlxTmp = - xy_IntSig_BtmToTop_xyz(xyz_FlxLat) * xy_totDepth 
          xy_HFlxTot = xy_HFlxTot + xy_HFlxTmp
          xy_VFlxTot = xy_VFlxTot + xy_IntSig_BtmToTop_xyz(xyz_Dz_xyz(-xyz_FlxSig)) * xy_totDepth
          call HistoryPut( SGSEDDANAKEY_ENGYFLXLAT_BOLUS,                                    &
               Fact_W2PW * RefDens*Cp0 * RPlanet*y_IntLon_xy(xy_HFlxTmp*xyz_CosLat(:,:,0)),  &
               hst_SGSEddMix )
       end if

       !* Advective flux

       xyz_FlxLon = xyz_PTemp * xyz_U
       xyz_FlxLat = xyz_PTemp * xyz_V
       xyz_FlxSig = 0d0
       if ( TendencyAnaFlag .or. TendencyVIntAnaFlag ) then
          call calc_DivFlux( xyz_Tendency, xyz_FlxLon, xyz_FlxLat, xyz_FlxSig )
          xyz_Tendency = - xyz_Tendency &
               + xyz_PTemp * eval_Div(xyz_U, xyz_V) - xyz_SigDot*xyz_DSig_xyz(xyz_PTemp)

          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_ADV, xyz_Tendency, hst_SGSEddMix )
          if ( TendencyVIntAnaFlag ) then
             xy_TendencyVInt = xy_IntSig_BtmToTop_xyz(xyz_Tendency)
             xy_TendTotVInt = xy_TendTotVInt + xy_TendencyVInt
             call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_ADV, xy_TendencyVInt, hst_SGSEddMix )
          end if
       end if
       
       if ( EngyFlxLatAnaFlag ) then
          
          xy_HFlxTmp = xy_IntSig_BtmToTop_xyz( xyz_IntSig_SigToTop_xyz(xyz_V)*xyz_DSig_xyz(xyz_PTemp) )*xy_totDepth
          xy_HFlxTot = xy_HFlxTot + xy_HFlxTmp
          
          call HistoryPut( SGSEDDANAKEY_ENGYFLXLAT_EULER,                                    &
               Fact_W2PW * RefDens*Cp0 * RPlanet*y_IntLon_xy(xy_HFlxTmp*xyz_CosLat(:,:,0)),  &
               hst_SGSEddMix )
          
          call HistoryPut( SGSEDDANAKEY_ENGYFLXLAT_TOTAL,                                    &
               Fact_W2PW * RefDens*Cp0 * RPlanet*y_IntLon_xy(xy_HFlxTot*xyz_CosLat(:,:,0)),  &
               hst_SGSEddMix )               

          xyz_FlxLon = 0d0
          call HistoryPut( SGSEDDANAKEY_ENGYFLXLAT_TOTAL_HDIV, &
               - RefDens*Cp0* y_AvrLon_xy( &
                 xy_w(w_AlphaOptr_xy(xyz_FlxLon(:,:,0)*xyz_CosLat(:,:,0), xy_HFlxTot*xyz_CosLat(:,:,0))) &
                 + xy_VFlxTot),                                                                 &
               hst_SGSEddMix )
       end if

       if ( TendencyAnaFlag .or. TendencyVIntAnaFlag ) then

          !* Vertical diffusion
          !
          call calc_vViscDiffCoef_MixLyrSimple( xyz_VViscCoef, xyz_VDiffCoef, &
               & xyz_Depth  )
          xyz_VDiffCoef = xyz_VDiffCoef + vDiffCoef
          xyz_Tendency =  xyz_DSig_xyz(xyz_VDiffCoef*xyz_DSig_xyz(xyz_PTemp))/spread(xy_totDepth**2,3,kMax+1)

          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_VDIFF, xyz_Tendency, hst_SGSEddMix )

          if ( TendencyVIntAnaFlag ) then
             xy_TendencyVInt = xy_IntSig_BtmToTop_xyz(xyz_Tendency)
             xy_TendTotVInt = xy_TendTotVInt + xy_TendencyVInt
             call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_VDIFF, xy_TendencyVInt, hst_SGSEddMix )
          end if

          !* Horizontal diffusion
          !
          w_LaplaEigVal = rn(:,1) / RPlanet**2
          VisCoef = ( nMax*(NMax + 1) / RPlanet**2 )**(-HDOrder/2) / HDEFoldTimeSec
          w_HDifCoefH = - VisCoef * ( ( - W_LaplaEigVal )**(HDOrder/2) )
          do k=0, kMax
             xyz_Tendency(:,:,k) = xy_w( w_HDifCoefH * w_xy(xyz_PTemp(:,:,k)) )
          end do
          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_HDIFF, xyz_Tendency, hst_SGSEddMix )
          if ( TendencyVIntAnaFlag ) then
             xy_TendencyVInt = xy_IntSig_BtmToTop_xyz(xyz_Tendency)
             xy_TendTotVInt = xy_TendTotVInt + xy_TendencyVInt
             call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_HDIFF, xy_TendencyVInt, hst_SGSEddMix )
          end if
          
          !* Local time derivative
          !
          if (FirstStepFlag) then
             allocate( xyz_PTempB(0:iMax-1,jMax,0:kMax) )
             xyz_PTempB = xyz_PTemp
             FirstStepFlag = .false.
          end if

          xyz_Tendency = (xyz_PTemp - xyz_PTempB)/(20d0*24d0*3600d0)
          xyz_PTempB = xyz_PTemp
          
          if ( TendencyAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TEND_LTD, xyz_Tendency, hst_SGSEddMix )
          if ( TendencyVIntAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_LTD, xy_IntSig_BtmToTop_xyz(xyz_Tendency), hst_SGSEddMix )


          !* Total tendency
          !
          if ( TendencyVIntAnaFlag ) &
               call HistoryPut( SGSEDDANAKEY_PT_TENDVINT_TOT, xy_TendTotVInt, hst_SGSEddMix )
          
       end if

       
    end if

    
!!$    ! For salinity
!!$    if ( TendencyAnaFlag ) then
!!$       ! Isopycnal diffusive flux
!!$       call calc_IsopycDiffFlux( xyz_FlxLon, xyz_FlxLat, xyz_FlxSig, &
!!$            xyz_Salt, wz_xyz(xyz_Salt), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )
!!$       call calc_DivFlux( xyz_Tendency, xyz_FlxLon, xyz_FlxLat, xyz_FlxSig )
!!$       call HistoryPut( SGSEDDANAKEY_SALTTEND_IPDIFF, xyz_Tendency, hst_SGSEddMix )
!!$
!!$       ! GM flux
!!$       call calc_SkewFlux( xyz_FlxLon, xyz_FlxLat, xyz_FlxSig, &
!!$            xyz_PTemp, wz_xyz(xyz_PTemp), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )
!!$       call calc_DivFlux( xyz_Tendency, xyz_FlxLon, xyz_FlxLat, xyz_FlxSig )
!!$       call HistoryPut( SGSEDDANAKEY_SALTTEND_GM, xyz_Tendency, hst_SGSEddMix )
!!$    end if
    
  contains

    subroutine calc_DivFlux( xyz_DivFlx,   & ! (out)
      xyz_FLon, xyz_FLat, xyz_FSig )         ! (in)

      real(DP), intent(out) :: xyz_DivFlx(0:iMax-1,jMax,0:kMax)
      real(DP), intent(in) :: xyz_FLon(0:iMax-1,jMax,0:kMax)
      real(DP), intent(in) :: xyz_FLat(0:iMax-1,jMax,0:kMax)
      real(DP), intent(in) :: xyz_FSig(0:iMax-1,jMax,0:kMax)

      xyz_DivFlx(:,:,:) = &
             xyz_wz( wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat) ) &
             + xyz_Dz_xyz(xyz_FSig)
      
    end subroutine calc_DivFlux
    
  end subroutine SGSEddyMixAnalysis_perform

  
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
    integer, parameter :: nAxes = 4
    integer, parameter :: AXES_LON = 1
    integer, parameter :: AXES_LAT = 2
    integer, parameter :: AXES_SIG = 3
    integer, parameter :: AXES_TIME = 2

    character(TOKEN) :: axes_names(nAxes) 
    character(TOKEN) :: axes_long_names(nAxes)
    character(TOKEN) :: axes_units(nAxes)

    
    character(TOKEN) :: axes_YT(2)
    character(TOKEN) :: axes_YZT(3) 
    character(TOKEN) :: axes_XYT(3) 
    character(TOKEN) :: axes_XYZT(4)

    ! 実行文; Executable statement
    !

    axes_names(1) = "lon"
    axes_long_names(1) = "longitude"
    axes_units(1) = "degree_east"

    axes_names(2) = "lat"
    axes_long_names(2) = "latitude"
    axes_units(2) = "degree_north"

    axes_names(3) = "sig"
    axes_long_names(3) = "nondimensional depth"
    axes_units(3) = "1"

    axes_names(4) = "time"
    axes_long_names(4) = "time"
    axes_units(4) = diagVar_gthsInfo%intUnit

    axes_YT(1)  = axes_names(2); axes_YT(2)  = axes_names(4)
    axes_YZT(:) = axes_names(2:4)
    axes_XYT(1:2) = axes_names(1:2); axes_XYT(3) = axes_names(4)
    axes_XYZT(:) = axes_names(:)
    
    call HistoryCreate( & 
         & file= trim(diagVar_gthsInfo%FilePrefix) // 'SGSEddyMixAnalysis.nc', title='analysis of GM scheme', &
         & source='Dennou-OGCM', &
         & institution='Dennou-OGCM project', &
         & dims=axes_names, dimsizes=(/ iMax, jMax, kMax+1, 0 /), &
         & longnames=axes_long_names, units=axes_units, & 
         & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
         & history=hst_SGSEddMix )  

    !
    call HistoryPut(axes_names(1), xyz_Lon*180d0/PI, hst_SGSEddMix)
    call HistoryPut(axes_names(2), xyz_Lat*180d0/PI, hst_SGSEddMix)
    call HistoryPut(axes_names(3), g_Sig, hst_SGSEddMix)

    !
    if ( BolusVelAnaFlag ) then
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_BOLUSU, dims=axes_XYZT, & 
            & longname='bolus velocity(longitude)', units='m/s', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_BOLUSV, dims=axes_XYZT, & 
            & longname='bolus velocity(latitude)', units='m/s', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_BOLUSW, dims=axes_XYZT, & 
            & longname='bolus velocity(vertical)', units='m/s', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_BOLUSMSTREAMFUNC, dims=axes_YZT, & 
            & longname='mass stream function on meriodinal plane associated with bolus velocity', units='Sv', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_RESMSTREAMFUNC, dims=axes_YZT, & 
            & longname='residual mass stream function on meriodinal plane', units='Sv', &
            & history=hst_SGSEddMix )
    end if

    if ( EngyFlxLatAnaFlag ) then
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_ENGYFLXLAT_TOTAL, dims=axes_YT, &
            & longname='total meridional heat flux', units='PW', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_ENGYFLXLAT_TOTAL_HDIV, dims=axes_YT, &
            & longname='divergence of total meridional heat flux', units='W.m-2', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname="VDiff", dims=axes_YT, &
            & longname='divergence of total meridional heat flux', units='W.m-2', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_ENGYFLXLAT_EULER, dims=axes_YT, &
            & longname='meridional heat flux associated with eulerian tranposrt', units='PW',  &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_ENGYFLXLAT_BOLUS, dims=axes_YT, &
            & longname='meridional heat flux associated with bolus tranposrt', units='PW',     &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_ENGYFLXLAT_IPDiff, dims=axes_YT, &
            & longname='meridional heat flux associated with isopycnal diffusion', units='PW', &
            & history=hst_SGSEddMix )
    end if

    if ( TendencyAnaFlag ) then
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TEND_IPDIFF, dims=axes_XYZT, &
            & longname='tendency of potential temperature due to isopycnal diffusion', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TEND_GM, dims=axes_XYZT, &
            & longname='tendency of potential temperature due to GM', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TEND_ADV, dims=axes_XYZT, &
            & longname='tendency of potential temperature due to advection', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TEND_VDIFF, dims=axes_XYZT, &
            & longname='tendency of potential temperature due to vertical diffusion', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TEND_LTD, dims=axes_XYZT, &
            & longname='tendency of potential temperature due to local time derivative', units='K/s', &
            & history=hst_SGSEddMix )
    end if

    if ( TendencyVIntAnaFlag ) then
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_IPDIFF, dims=axes_XYT, &
            & longname='vertically integrated tendency of potential temperature due to isopycnal diffusion', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_GM, dims=axes_XYT, &
            & longname='vertically integrated  tendency of potential temperature due to GM', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_ADV, dims=axes_XYT, &
            & longname='vertically integrated  tendency of potential temperature due to advection', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_HDIFF, dims=axes_XYT, &
            & longname='vertically integrated  tendency of potential temperature due to horizontal (hyper) diffusion', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_VDIFF, dims=axes_XYT, &
            & longname='vertically integrated  tendency of potential temperature due to vertical diffusion', units='K/s', &
            & history=hst_SGSEddMix )
       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_LTD, dims=axes_XYT, &
            & longname='vertically integrated  local time derivative of potential temperature', units='K/s', &
            & history=hst_SGSEddMix )

       call HistoryAddVariable( &
            & varname=SGSEDDANAKEY_PT_TENDVINT_TOT, dims=axes_XYT, &
            & longname='vertically integrated  total tendency', units='K/s', &
            & history=hst_SGSEddMix )
       
       
!!$       call HistoryAddVariable( &
!!$            & varname=SGSEDDANAKEY_SALTTEND_IPDIFF, dims=axes_XYZT, &
!!$            & longname='tendency of salinity  due to isopycnal diffusion', units='psu/s', &
!!$            & history=hst_SGSEddMix )
!!$       call HistoryAddVariable( &
!!$            & varname=SGSEDDANAKEY_SALTTEND_GM, dims=axes_XYZT, &
!!$            & longname='tendency of salinity due to GM', units='psu/s', &
!!$            & history=hst_SGSEddMix )
    end if

  end subroutine prepair_Output

  
  
end module SGSEddyMixAnalysis_mod
