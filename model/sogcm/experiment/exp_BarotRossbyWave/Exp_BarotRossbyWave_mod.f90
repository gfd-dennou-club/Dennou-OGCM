!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_BarotRossbyWave_mod 

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message

  use Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefTemp

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_BarotRossbyWave_Init, Exp_BarotRossbyWave_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_BarotRossbyWave_mod' !< Module Name

  integer, save :: m               !< The number of zonal mode (zonal wave number, m of spherical harmonics Y^m_n)
  integer, save :: n               !< Total wave number (n of spherical harmonics Y^m_n)

  real(DP), parameter :: GEOPOT0             = 2.94d04            !< (Constant) Depth
  real(DP), parameter :: KEDensAvg           = 1.0d-12            !< The average of kinetic energy density
  
contains

  !>
  !!
  !!
  subroutine Exp_BarotRossbyWave_Init(configNml)

    ! 実行文; Executable statements
    !
    character(*), intent(in) :: configNml

    call read_expConfig(configNml)

  end subroutine Exp_BarotRossbyWave_Init

  !>
  !!
  !!
  subroutine Exp_BarotRossbyWave_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_BarotRossbyWave_Final

  !> @brief 
  !!
  !!
  subroutine setInitCondition()
    
    !
    !
    use GridSet_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & xyz_Lat, xyz_Lon

    use VariableSet_mod

    use SpmlUtil_mod
    use w_module, only : l_nm, AvrLonLat_xy

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0
    real(DP) :: wz_Psi(lMax,0:kMax)
    integer :: l
    real(DP) :: KEDensAvg0

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    ! Set the fluid depth
    xy_totDepthBasic = GEOPOT0 / Grav 
    xy_SurfHeightN = 0d0


    ! Set the initial perturbation of velociy field
    l = l_nm(n, m)
    wz_Psi = 0d0
    wz_Psi(l,:) = 1d0
    if( m /= 0 ) wz_Psi(l,:) = wz_Psi(l,:)/sqrt(2d0)

    xyz_UN = - xya_GradLat_wa(wz_Psi)/RPlanet
    xyz_VN =   xya_GradLon_wa(wz_Psi)/RPlanet
    KEDensAvg0 = AvrLonLat_xy(xyz_UN(:,:,1)**2 + xyz_VN(:,:,1)**2)
    xyz_UN = sqrt(KEDensAvg/KEDensAvg0)*xyz_UN
    xyz_VN = sqrt(KEDensAvg/KEDensAvg0)*xyz_VN


    ! Set the neutral stratification(i.e, constant potential temperature) 
    z_PTempBasic = RefTemp

  end subroutine setInitCondition

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

end module Exp_BarotRossbyWave_mod

