!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Exp_BarotRossby_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use DOGCM_Admin_Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefTemp

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
  character(*), parameter:: module_name = 'DOGCM_BarotRossy_mod' !< Module Name

  integer, save :: m               !< The number of zonal mode (zonal wave number, m of spherical harmonics Y^m_n)
  integer, save :: n               !< Total wave number (n of spherical harmonics Y^m_n)

  real(DP), parameter :: GEOPOT0             = 2.94d04            !< (Constant) Depth
  real(DP), parameter :: KEDensAvg           = 1.0d-12            !< The average of kinetic energy density
  
contains

  !>
  !!
  !!
  subroutine DOGCM_Exp_Init(configNml)

    ! 実行文; Executable statements
    !
    character(*), intent(in) :: configNml

    call read_expConfig(configNml)

  end subroutine DOGCM_Exp_Init

  !>
  !!
  !!
  subroutine DOGCM_Exp_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_Exp_Final

  !> @brief 
  !!
  !!
  subroutine DOGCM_Exp_SetInitCond
    
    !
    !
    use DOGCM_Admin_Grid_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & IS,IE, JS,JE, KS,KE,          &
         & xyz_Lat, xyz_Lon, xy_Topo

    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyza_U, xyza_V, xya_SSH, xyzaa_TRC, &
         & TRCID_PTEMP

    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0

    real(DP) :: wz_Psi(lMax,0:kMax)
    real(DP) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_V(0:iMax-1,jMax,0:kMax)
    integer :: l
    real(DP) :: KEDensAvg0

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    ! Set the fluid depth
    xy_Topo = GEOPOT0 / Grav 
    xya_SSH(:,:,TIMELV_ID_N) = 0d0


    ! Set the initial perturbation of velociy field
    l = l_nm(n, m)
    wz_Psi = 0d0
    wz_Psi(l,:) = 1d0
    if( m /= 0 ) wz_Psi(l,:) = wz_Psi(l,:)/sqrt(2d0)

    xyz_U = - xya_GradLat_wa(wz_Psi)/RPlanet
    xyz_V =   xya_GradLon_wa(wz_Psi)/RPlanet

    KEDensAvg0 = AvrLonLat_xy(xyz_U(:,:,1)**2 + xyz_V(:,:,1)**2)
    xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) = sqrt(KEDensAvg/KEDensAvg0)*xyz_U
    xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) = sqrt(KEDensAvg/KEDensAvg0)*xyz_V

    ! Set the neutral stratification(i.e, constant potential temperature) 
    xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP, TIMELV_ID_N) = RefTemp

  end subroutine DOGCM_Exp_SetInitCond

  subroutine DOGCM_Exp_Do()

  end subroutine DOGCM_Exp_Do
  
  !-----------------------------------------------------------------------------

  
  !> @brief 
  !!
  !!
  subroutine read_expConfig(configNmlFileName)

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
    character(*), intent(in) :: configNmlFileName

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
    namelist /numexp_nml/ &
         & m, n

    ! 実行文; Executable statement
    !

    m = 0
    n = 1
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                  ! (in)
            & nml = numexp_nml, &         ! (out)
            & iostat = iostat_nml )       ! (out)
       close( unit_nml )
    end if

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, " Parameters of perturbation: m=%d,n=%d", i=(/ m,n /))

    
  end subroutine read_expConfig

end module DOGCM_Exp_BarotRossby_mod

