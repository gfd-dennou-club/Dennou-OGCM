!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module SpectralAnalysis_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  use gtool_history

  use Constants_mod, only: &
       & RPlanet, Grav, RefDens, &
       & hViscCoef, vViscCoef, hHyperViscCoef, vHyperViscCoef, &
       & hDiffCoef, vDiffCoef

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, nMax, tMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use DiagVarFileSet_mod
  use DiagnoseUtil_mod
  use DiagVarEval_mod


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SpectralAnalysis_Init, SpectralAnalysis_Final
  public :: SpectralAnalysis_Perform

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SpectralAnalysis_mod' !< Module Name

  character(*), parameter, public :: SPECTRALANAKEY_ENERGY = 'Energy'
  character(*), parameter, public :: SPECTRALANAKEY_ENSTROPHY = 'Enstrophy'

  logical :: EnergySpectralAnaFlag
  logical :: EnstrophySpectralAnaFlag

  type(gt_history), save :: hst_spectralAna

contains

  !>
  !!
  !!
  subroutine SpectralAnalysis_Init(diagVar_gthsInfo, spectralAnaName)

    ! 宣言文; Declare statements
    !
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo
    character(*), intent(in) :: spectralAnaName(:)

    ! 作業変数
    ! Work variables
    !
    integer :: n


    ! 実行文; Executable statements
    !

    EnergySpectralAnaFlag = .false.
    EnstrophySpectralAnaFlag = .false.

    do n=1, size(spectralAnaName)
       select case(spectralAnaName(n))
       case(SPECTRALANAKEY_ENERGY)
          EnergySpectralAnaFlag = .true.
       case(SPECTRALANAKEY_ENSTROPHY)
          EnstrophySpectralAnaFlag = .true.
       case ('')
       case Default
          call MessageNotify('E', module_name, &
               & "The specified type of spectral analysis '%c' is invalid.", c1=trim(spectralAnaName(n)) )
       end select
    end do

    call prepair_Output(diagVar_gthsInfo)

  end subroutine SpectralAnalysis_Init

  !>
  !!
  !!
  subroutine SpectralAnalysis_Final()

    ! 実行文; Executable statements
    !
    
    if(EnergySpectralAnaFlag) then
       call HistoryClose(hst_spectralAna)
    end if

  end subroutine SpectralAnalysis_Final


  !> @brief 
  !!
  !!
  subroutine SpectralAnalysis_Perform()

    ! モジュール引用 ; Use statements
    !
    use VariableSet_mod, only: &
         & xyz_UN, xyz_VN

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !

    real(DP) :: xyz_Psi(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !

    if(.not.(EnergySpectralAnaFlag .or. EnstrophySpectralAnaFlag)) &
         & return

    call eval_StreamPot(xyz_Psi=xyz_Psi, &
         & xyz_u=xyz_UN, xyz_v=xyz_VN)
    
    if(EnergySpectralAnaFlag) call perform_energySpectralAna()

    if(EnstrophySpectralAnaFlag) call perform_enstrophySpectralAna()

  contains
    subroutine perform_energySpectralAna()
      real(DP) :: n_EnergySpectral(nMax+1)

      n_EnergySpectral = w_IntSig_BtmToTop_wz( &
           & na_EnergyFromStreamfunc_wa(wz_xyz(xyz_Psi)) &
           & ) 

      call HistoryPut(SPECTRALANAKEY_ENERGY, n_EnergySpectral, hst_spectralAna)
    end subroutine perform_energySpectralAna

    subroutine perform_enstrophySpectralAna()
      real(DP) :: n_EnstrophySpectral(lMax)

      n_EnstrophySpectral = w_IntSig_BtmToTop_wz( &
           & na_EnstrophyFromStreamfunc_wa(wz_xyz(xyz_Psi)) &
           & )

      call HistoryPut(SPECTRALANAKEY_ENSTROPHY, n_EnstrophySpectral, hst_spectralAna)
    end subroutine perform_enstrophySpectralAna

  end subroutine SpectralAnalysis_Perform


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!

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
    integer, parameter :: nAxes = 2
    integer, parameter :: AXES_TOTWAVENUM = 1
    integer, parameter :: AXES_TIME = 2


    character(TOKEN) :: axes_names(nAxes) 
    character(TOKEN) :: axes_long_names(nAxes)
    character(TOKEN) :: axes_units(nAxes)
    
    integer :: n
    real(DP) :: n_totWaveNum(nMax+1)

    ! 実行文; Executable statement
    !

    axes_names(1) = "n"
    axes_long_names(1) = "total wave number"
    axes_units(1) = "m-1"

    axes_names(2) = "t"
    axes_long_names(2) = "time"
    axes_units(2) = diagVar_gthsInfo%intUnit

    if ( EnergySpectralAnaFlag ) then
       call HistoryCreate( & 
            & file= trim(diagVar_gthsInfo%FilePrefix) // 'SpectralAnalysis.nc', title='spectral analysis', &
            & source='Dennou-OGCM', &
            & institution='Dennou-OGCM project', &
            & dims=axes_names, dimsizes=(/ lMax, 0 /), &
            & longnames=axes_long_names, units=axes_units, & 
            & origin=real(diagVar_gthsInfo%origin), interval=real(diagVar_gthsInfo%intValue), &
            & history=hst_spectralAna )  

       !
       forAll(n=0:nMax) n_totWaveNum(n+1) = n
       call HistoryPut(axes_names(AXES_TOTWAVENUM), n_totWaveNum, hst_spectralAna)

       call HistoryAddVariable( & 
             & varname=SPECTRALANAKEY_ENERGY, dims=(/ axes_names(AXES_TOTWAVENUM), axes_names(AXES_TIME) /), &
             & longname='energy spectral', units='J*m', xtype='double',&
             & history=hst_spectralAna)

        call HistoryAddVariable( & 
             & varname=SPECTRALANAKEY_ENSTROPHY, dims=(/ axes_names(AXES_TOTWAVENUM), axes_names(AXES_TIME) /), &
             & longname='enstrophy spectral', units='s-2*m', xtype='double',&
             & history=hst_spectralAna)
    end if

  end subroutine prepair_Output

end module SpectralAnalysis_mod

