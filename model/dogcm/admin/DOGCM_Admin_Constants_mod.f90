!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_Constants_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_Admin_Constants_Init, DOGCM_Admin_Constants_Final

  ! 公開変数
  ! Public variables
  !
  real(DP), parameter, public :: UNDEFVAL = -9.9999D30

  
  real(DP), parameter, public :: PI = 3.1415926535897932_DP
                              ! $ \pi $ .
                              ! 円周率.  Circular constant
  real(DP), parameter, public:: GasRUniv = 8.314_DP
                              ! $ R^{*} $ [J K-1 mol-1].
                              ! 普遍気体定数.  Universal gas constant
  real(DP), parameter, public:: StB = 5.67e-8_DP
                              ! $ \sigma_{SB} $ . 
                              ! ステファンボルツマン定数. 
                              ! Stefan-Boltzmann constant

  real(DP), save, public:: RPlanet
                              ! $ a $ [m]. 
                              ! 惑星半径. 
                              ! Radius of planet
  real(DP), save, public:: Omega
                              ! $ \Omega $ [s-1]. 
                              ! 回転角速度. 
                              ! Angular velocity
  real(DP), save, public:: Grav
                              ! $ g $ [m s-2]. 
                              ! 重力加速度. 
                              ! Gravitational acceleration

  real(DP), public, save :: hViscCoef  
  real(DP), public, save :: vViscCoef
  real(DP), public, save :: hHyperViscCoef
  real(DP), public, save :: vHyperViscCoef

  real(DP), public, save :: hDiffCoef  
  real(DP), public, save :: vDiffCoef
  real(DP), public, save :: hHyperDiffCoef
  real(DP), public, save :: vHyperDiffCoef

  real(DP), public, save :: RoughnessParamBottom

  real(DP), public, save :: RefDens    
  real(DP), public, save :: RefTemp    
  real(DP), public, save :: Cp0         
  real(DP), public, save :: RefSoundSpeed
  real(DP), public, save :: ThermalExpanCoef 

  !< Albedo of open ocean [(1)]  
  real(DP), public :: albedoOcean

  !< Emissivity of open ocean surface [(1)]  
  real(DP), public :: emissivOcean

  !
  real(DP), public :: LatentHeat
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'DOGCM_Admin_Constants_mod'
                              ! モジュールの名称. 
                              ! Module name  
contains

  !>
  !!
  !!
  subroutine DOGCM_Admin_Constants_Init(configNmlName)

    ! 引用文; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! メッセージ出力
    ! Message output
    !
    use dc_message, only: MessageNotify

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
                              ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
                              ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /constants_nml/ &
      & RefDens, RefTemp, ThermalExpanCoef, RefSoundSpeed, Cp0, &
      & RPlanet, Omega, Grav, &
      & hViscCoef, vViscCoef, &
      & hHyperViscCoef, vHyperViscCoef, &
      & hDiffCoef, vDiffCoef, &
      & hHyperDiffCoef, vHyperDiffCoef, &
      & RoughnessParamBottom, &
      & emissivOcean, albedoOcean, &
      & LatentHeat

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    RPlanet          = 6.371e6_DP
    Omega            = 2.0_DP * PI / ( 60.0_DP * 60.0_DP * 23.9345_DP )
    Grav             = 9.8_DP

    vViscCoef        = 0d0
    hViscCoef        = 0d0
    vHyperViscCoef   = 0d0
    hHyperViscCoef   = 0d0
    vDiffCoef        = 0d0
    hDiffCoef        = 0d0
    vHyperDiffCoef   = 0d0
    hHyperDiffCoef   = 0d0
    
    RefDens          = 1.027d03
    RefTemp          = 283d0
    Cp0              = 3986d0
    RefSoundSpeed    = 1490d0
    ThermalExpanCoef = 1.67d-04
    
    AlbedoOcean = 0.1d0
    EmissivOcean = 0.97d0

    LatentHeat = 2.501d6
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = constants_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
   end if


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  RPlanet          = %f', d = (/ RPlanet          /) )
    call MessageNotify( 'M', module_name, '  Omega            = %f', d = (/ Omega            /) )
    call MessageNotify( 'M', module_name, '  Grav             = %f', d = (/ Grav             /) )
    call MessageNotify( 'M', module_name, '  vViscCoef        = %f', d = (/ vViscCoef        /) )
    call MessageNotify( 'M', module_name, '  hViscCoef        = %f', d = (/ hViscCoef        /) )
    call MessageNotify( 'M', module_name, '  vHyperViscCoef        = %f', d = (/ vHyperViscCoef        /) )
    call MessageNotify( 'M', module_name, '  hHyperViscCoef        = %f', d = (/ hHyperViscCoef        /) )
    call MessageNotify( 'M', module_name, '  vDiffCoef        = %f', d = (/ vDiffCoef        /) )
    call MessageNotify( 'M', module_name, '  hDiffCoef        = %f', d = (/ hDiffCoef        /) )
    call MessageNotify( 'M', module_name, '  vHyperDiffCoef        = %f', d = (/ vHyperDiffCoef        /) )
    call MessageNotify( 'M', module_name, '  hHyperDiffCoef        = %f', d = (/ hHyperDiffCoef        /) )
    call MessageNotify( 'M', module_name, '  RoughnessParamBottom = %f', d = (/ RoughnessParamBottom  /) )
    call MessageNotify( 'M', module_name, '  RefDens           = %f', d=(/ RefDens /) )
    call MessageNotify( 'M', module_name, '  RefTemp           = %f', d=(/ RefTemp /) )
    call MessageNotify( 'M', module_name, '  Cp0               = %f', d=(/ Cp0 /) )
    call MessageNotify( 'M', module_name, '  RefSoundSpeed     = %f', d=(/ RefSoundSpeed /) )
    call MessageNotify( 'M', module_name, '  ThermalExpanCoef  = %f', d=(/ ThermalExpanCoef /) )
    call MessageNotify( 'M', module_name, '  AlbedoOcean       = %f', d=(/ AlbedoOcean /) )
    call MessageNotify( 'M', module_name, '  EmissivOcean      = %f', d=(/ EmissivOcean /) )

end subroutine DOGCM_Admin_Constants_Init

!>
!!
!!
subroutine DOGCM_Admin_Constants_Final()

  ! 実行文; Executable statements

end subroutine DOGCM_Admin_Constants_Final

end module DOGCM_Admin_Constants_mod

