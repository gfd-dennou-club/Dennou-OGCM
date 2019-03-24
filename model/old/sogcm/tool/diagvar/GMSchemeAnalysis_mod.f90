!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module GMSchemeAnalysis_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  use Constants_mod, only: &
       & PI, RPlanet, Grav, Omega, RefDens, &
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
       & calc_IsoNeutralSlope, calc_BolusVelocity

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: GMSchemeAnalysis_Init, GMSchemeAnalysis_Final
  public :: GMSchemeAnalysis_perform


  ! 公開変数
  ! Public variable
  !
  character(*), parameter, public :: GMANAKEY_BOLUSU = 'BolusU'
  character(*), parameter, public :: GMANAKEY_BOLUSV = 'BolusV'
  character(*), parameter, public :: GMANAKEY_BOLUSW = 'BolusW'
  character(*), parameter, public :: GMANAKEY_BOLUSMSTREAMFUNC = 'BolusMStreamFunc'
  character(*), parameter, public :: GMANAKEY_RESMSTREAMFUNC = 'ResMStreamFunc'
  character(*), parameter, public :: GMANAKEY_PTEMP_ISOPYCDIFF_TENDENCY = 'PTempIsopycDiffTend'
  character(*), parameter, public :: GMANAKEY_SALT_ISOPYCDIFF_TENDENCY = 'SaltIsopycDiffTend'
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'GMSchemeAnalysis_mod' !< Module Name

  type(gt_history), save :: hst_GMAna
  logical :: GMSchemeAnaFlag

contains

  !>
  !!
  !!
  subroutine GMSchemeAnalysis_Init(configNmlFile, diagVar_gthsInfo, isActivated)

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

    GMSchemeAnaFlag = isActivated

    if(GMSchemeAnaFlag) then
       call SGSEddyMixing_Init(configNmlFileName=configNmlFile)
       call prepair_Output(diagVar_gthsInfo)
    end if
    
  end subroutine GMSchemeAnalysis_Init

  !>
  !!
  !!
  subroutine GMSchemeAnalysis_Final()

    ! 実行文; Executable statements
    !

  end subroutine GMSchemeAnalysis_Final

  !> @brief 
  !!
  !!
  subroutine GMSchemeAnalysis_perform(xyz_V, xyz_PTemp, xyz_Salt, xy_totDepth)
    
    ! モジュール引用 ; Use statements
    !
    use SGSEddyMixingHelper_mod, only: &
         & prepare_SlopeTapering
    use SGSEddyMixing_mod, only: &
         & calc_IsoNeutralSlope, calc_BolusVelocity

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_V, xyz_PTemp, xyz_Salt
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_totDepth

    ! 局所変数
    ! Local variables
    !
    integer :: k

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_Depth, xyz_totDepth
    
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DensPot, xyz_SLon, xyz_SLat, xyz_T
    
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_IsopycDiffTerm, &
         & xyz_BolusU, xyz_BolusV, xyz_BolusW
    
    ! 実行文; Executable statement
    !

    if(.not. GMSchemeAnaFlag) return


    call MessageNotify('M', module_name, 'Perform the analysis of GM scheme..')    

    ! Preparation 
    !

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


    
    ! Calculate bolus velocities and the corresponding stream function, and
    ! output them. 
    !

    call calc_BolusVelocity(xyz_BolusU, xyz_BolusV, xyz_BolusW, & ! (out)
         & xyz_SLon, xyz_SLat, xyz_Depth, xyz_T )                 ! (in)    
    call HistoryPut(GMANAKEY_BOLUSU, xyz_BolusU, hst_GMAna)
    call HistoryPut(GMANAKEY_BOLUSV, xyz_BolusV, hst_GMAna)
    call HistoryPut(GMANAKEY_BOLUSW, xyz_BolusW, hst_GMAna)

    
    call HistoryPut(GMANAKEY_BOLUSMSTREAMFUNC, &
         & eval_GM_BolusMStreamFunc(xyz_BolusV, xy_totDepth), hst_GMAna )
    call HistoryPut(GMANAKEY_RESMSTREAMFUNC,   &
         & eval_GM_ResMStreamFunc(xyz_V, xyz_BolusV, xy_totDepth), hst_GMAna)

    write(*,*) "Sum BolusV=", AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyz_BolusV) )
    
    ! Calculte tendencies of potential temperature and salinity due to isopycnal diffusion, and
    ! output them.

    ! For potential temperature
    call eval_Redi_IsopycDiffTermTendency( xyz_IsopycDiffTerm,  & ! (out)
       & xyz_PTemp, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )        ! (in)
    call HistoryPut(GMANAKEY_PTEMP_ISOPYCDIFF_TENDENCY, xyz_IsopycDiffTerm, hst_GMAna)

    ! For salinity
    call eval_Redi_IsopycDiffTermTendency( xyz_IsopycDiffTerm,  & ! (out)
       & xyz_Salt, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )        ! (in)
    call HistoryPut(GMANAKEY_SALT_ISOPYCDIFF_TENDENCY, xyz_IsopycDiffTerm, hst_GMAna)
    
  end subroutine GMSchemeAnalysis_perform


  !> @brief 
  !!
  !!
  subroutine eval_Redi_IsopycDiffTermTendency( xyz_IsopycDiffTerm, &
       & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth )

    ! モジュール引用; Use statement
    !
    use SGSEddyMixing_mod, only: &
         & calc_IsopycDiffFlux
    use SGSEddyMixingHelper_mod, only: &
         & xyz_Dz_xyz
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_IsopycDiffTerm
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_Tracer, xyz_SLon, xyz_SLat, xyz_T, xyz_Depth
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_FLon, xyz_FLat, xyz_FSig, xyz_CosLat
    
    ! 実行文; Executable statement
    !

    
    call calc_IsopycDiffFlux(xyz_FLon, xyz_FLat, xyz_FSig, &
         & xyz_Tracer, wz_xyz(xyz_Tracer), xyz_SLon, xyz_SLat, xyz_T, xyz_Depth)

    xyz_CosLat = cos(xyz_Lat)
    xyz_IsopycDiffTerm(:,:,:) = &
         &    xyz_wz(wz_AlphaOptr_xyz(xyz_FLon*xyz_CosLat, xyz_FLat*xyz_CosLat)) &
         & +  xyz_Dz_xyz(xyz_FSig)
    
  end subroutine eval_Redi_IsopycDiffTermTendency


  !> @brief 
  !!
  !!
  function eval_GM_BolusMStreamFunc(xyz_BolusV, xy_totDepth) result(yz_BolusMSF)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_BolusV(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP) :: yz_BolusMSF(jMax,0:kMax)
    
    ! 実行文; Executable statement
    !
    yz_BolusMSF = eval_MassStreamFunc(xyz_BolusV, xy_totDepth)
    
  end function eval_GM_BolusMStreamFunc

  !> @brief 
  !!
  !!
  function eval_GM_ResMStreamFunc(xyz_V, xyz_BolusV, xy_totDepth) result(yz_ResMSF)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_BolusV(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP) :: yz_ResMSF(jMax,0:kMax)
    
    ! 実行文; Executable statement
    !
    yz_ResMSF = &
         &   eval_MassStreamFunc(xyz_V, xy_totDepth)      &
         & + eval_MassStreamFunc(xyz_BolusV, xy_totDepth)
    
  end function eval_GM_ResMStreamFunc


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
    character(TOKEN) :: axes_XYZT(4), axes_YZT(3)

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

    axes_XYZT(:) = axes_names(:)
    axes_YZT(:) = axes_names(2:4)

    call HistoryCreate( & 
         & file= trim(diagVar_gthsInfo%FilePrefix) // 'GMSchemeAnalysis.nc', title='analysis of GM scheme', &
         & source='Dennou-OGCM', &
         & institution='Dennou-OGCM project', &
         & dims=axes_names, dimsizes=(/ iMax, jMax, kMax+1, 0 /), &
         & longnames=axes_long_names, units=axes_units, & 
         & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
         & history=hst_GMAna )  

    !
    call HistoryPut(axes_names(1), xyz_Lon*180d0/PI, hst_GMAna)
    call HistoryPut(axes_names(2), xyz_Lat*180d0/PI, hst_GMAna)
    call HistoryPut(axes_names(3), g_Sig, hst_GMAna)

    !
    call HistoryAddVariable( &
         & varname=GMANAKEY_BOLUSU, dims=axes_XYZT, & 
         & longname='bolus velocity(longitude)', units='m/s', &
         & history=hst_GMAna )

    call HistoryAddVariable( &
         & varname=GMANAKEY_BOLUSV, dims=axes_XYZT, & 
         & longname='bolus velocity(latitude)', units='m/s', &
         & history=hst_GMAna )

    call HistoryAddVariable( &
         & varname=GMANAKEY_BOLUSW, dims=axes_XYZT, & 
         & longname='bolus velocity(vertical)', units='m/s', &
         & history=hst_GMAna )

    call HistoryAddVariable( &
         & varname=GMANAKEY_BOLUSMSTREAMFUNC, dims=axes_YZT, & 
         & longname='mass stream function on meriodinal plane associated with bolus velocity', units='Sv', &
         & history=hst_GMAna )

    call HistoryAddVariable( &
         & varname=GMANAKEY_RESMSTREAMFUNC, dims=axes_YZT, & 
         & longname='residual mass stream function on meriodinal plane', units='Sv', &
         & history=hst_GMAna )
    
    call HistoryAddVariable( &
         & varname=GMANAKEY_PTEMP_ISOPYCDIFF_TENDENCY, dims=axes_XYZT, &
         & longname='tendency of potential temperature due to isopycnal diffusion', units='K/s', &
         & history=hst_GMAna )

    call HistoryAddVariable( &
         & varname=GMANAKEY_SALT_ISOPYCDIFF_TENDENCY, dims=axes_XYZT, &
         & longname='tendency of salinity  due to isopycnal diffusion', units='psu/s', &
         & history=hst_GMAna )
    

  end subroutine prepair_Output

end module GMSchemeAnalysis_mod

