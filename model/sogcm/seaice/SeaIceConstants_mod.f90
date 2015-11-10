!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SeaIceConstants_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SeaIceConstants_Init, SeaIceConstants_Final

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SeaIceConstants_mod' !< Module Name

  !< Stefan=Boltzmann constant [W.m-2.K-4]
  real(DP), public :: SBConst

  !< Density of ice [kg/m3]
  real(DP), public :: DensIce

  !< Density of snow [kg/m3]  
  real(DP), public :: DensSnow

  !< Density of sea water [kg/m3]  
  real(DP), public :: DensSeaWater

  !< Ice heat capacity (excluding internal melt) [J.kg-1.K-1]
  real(DP), public :: CIce

  !< Latent heat of freezeing [J.kg-1]
  real(DP), public :: LFreeze

  !< Thermal conductivity of sea ice [W.m-1.degC-1]
  real(DP), public :: KIce

  !< Thermal conductivity of snow [W.m-1.degC-1]
  real(DP), public :: KSnow

  !< Constant relating freezing temersture to salinity [degC.ppt-1]
  real(DP), public :: Mu

  !< Salinity of sea ice [ppt]
  real(DP), public :: SaltSeaIce

  !< Freezing temperature of (fresh)water [degC]
  real(DP), public :: FreezeTempWater
  
  !< Freezing temperature of seawater [degC]
  real(DP), public :: FreezeTempSW

  !< Albedo of open ocean
  real(DP), public :: AlbedoOcean
  
  !< Albedo of snow [(1)]
  real(DP), public :: AlbedoSnow

  !< Albedo of melting snow [(1)]  
  real(DP), public :: AlbedoMeltSnow
  
  !< Albedo of ice [(1)]    
  real(DP), public :: AlbedoIce

  !< Emissivity of open ocean surface
  real(DP), public :: EmissivOcean
  
  !< Emissivity of snow surface [(1)]
  real(DP), public :: EmissivSnow

  !< Emissivity of sea-ice surface [(1)]
  real(DP), public :: EmissivIce
  
  !> Fraction of the net incoming solar radiation which penetrates into the interior
  !! of snow-free ice [(1)]
  real(DP), public :: I0  

  !> Base-melting heat transfer coeffiecient (between ice and water) [(1)]
  !! This parameter corresponds
  !! \[
  !     \dfrac{1}{P_{\rm rt}} k^{-1} \ln{(-z/z_0)} + B_T}. 
  !! \]
  !!  Please check Eq.18a,b in Mellor and Kantha(1989).
  !!
  real(DP), public :: BaseMeltHeatTransCoef

  
contains

  !>
  !!
  !!
  subroutine SeaIceConstants_Init(ConfigNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName
    
    ! 実行文; Executable statements
    !
    
    call read_namelist(ConfigNmlName)
    
  end subroutine SeaIceConstants_Init

  !>
  !!
  !!
  subroutine SeaIceConstants_Final()

    ! 実行文; Executable statements
    !

  end subroutine SeaIceConstants_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_namelist(ConfigNmlName)

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
    namelist /SeaIceConstants_nml/ &
         & SBConst, &
         & DensIce, DensSnow, DensSeaWater, &
         & CIce, LFreeze, &
         & KIce, KSnow, &
         & Mu, SaltSeaIce, FreezeTempWater, FreezeTempSW, &
         & AlbedoOcean, AlbedoSnow, AlbedoMeltSnow, AlbedoIce, &
         & EmissivOcean, EmissivSnow, EmissivIce, &
         & I0

    ! 実行文; Executable statements
    !
    
    SBConst = 5.67d-8 !5.7948e-8!5.67d-8*1.0219d0
    DensIce = 905d0
    DensSnow = 330d0
    DensSeaWater = 1030d0
    CIce = 2100d0
    LFreeze = 334d3
    KIce = 2.03d0
    KSnow = 0.31d0
    
    Mu = 0.054d0
    SaltSeaIce = 1d0
    FreezeTempWater = 0d0
    FreezeTempSW = -2d0!-1.8d0
    
    AlbedoOcean = 0.1d0
    AlbedoSnow = 0.8d0
    AlbedoMeltSnow = 0.735d0
    AlbedoIce = 0.65d0

    EmissivOcean = 0.97d0
    EmissivSnow      = 0.97d0
    EmissivIce       = 0.97d0
    
    I0 = 0.34d0!*2d0

    baseMeltHeatTransCoef = 6d-3
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = SeaIceConstants_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
   end if


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, 'SBConst      = %f [W.m-2.K-4]', d=(/ SBConst /))
    call MessageNotify( 'M', module_name, 'DensSnow     = %f [kg/m3]',     d=(/ DensSnow /))    
    call MessageNotify( 'M', module_name, 'DensIce      = %f [kg/m3]',     d=(/ DensIce /))
    call MessageNotify( 'M', module_name, 'DensSeaWater = %f [kg/m3]',     d=(/ DensSeaWater /))
    call MessageNotify( 'M', module_name, 'KSnow        = %f [W.m-1.degC-1]',     d=(/ KSnow /))    
    call MessageNotify( 'M', module_name, 'KIce         = %f [W.m-1.degC-1]',     d=(/ KIce /))    
    call MessageNotify( 'M', module_name, 'AlbedoSnow       = %f [(1)]  ',     d=(/ AlbedoSnow   /))
    call MessageNotify( 'M', module_name, 'AlbedoMeltSnow   = %f [(1)]  ',     d=(/ AlbedoMeltSnow   /))
    call MessageNotify( 'M', module_name, 'AlbedoIce        = %f [(1)]  ',     d=(/ AlbedoIce   /))
    call MessageNotify( 'M', module_name, 'AlbedoOcean      = %f [(1)]  ',     d=(/ AlbedoOcean   /))
    call MessageNotify( 'M', module_name, 'EmissivIce       = %f [(1)]  ',     d=(/ EmissivIce   /))
    call MessageNotify( 'M', module_name, 'EmissivSnow      = %f [(1)]  ',     d=(/ EmissivSnow   /))
    call MessageNotify( 'M', module_name, 'EmissivOcean     = %f [(1)]  ',     d=(/ EmissivOcean   /))
    
  end subroutine read_namelist
end module SeaIceConstants_mod

