!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_InternalGravWave_mod 

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_InternalGravWave_Init, Exp_InternalGravWave_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_InternalGravWave_mod' !< Module Name

  real(DP), save :: BruntVaisFreq  !< Brunt-Vaisala Frequency
  integer, save :: m               !< The number of zonal mode (zonal wave number, m of spherical harmonics Y^m_n)
  integer, save :: n               !< Total wave number (n of spherical harmonics Y^m_n)
  integer, save :: l               !< The number of vertical mode

  real(DP), parameter :: H0  = 5.0d03         !< (Constant) Depth
  real(DP), parameter :: KEDensAvg           = 1.0d-12            !< The average of kinetic energy density

contains

  !>
  !!
  !!
  subroutine Exp_InternalGravWave_Init(configNml)

    ! 実行文; Executable statements
    !
    character(*), intent(in) :: configNml

    call read_expConfig(configNml)

  end subroutine Exp_InternalGravWave_Init

  !>
  !!
  !!
  subroutine Exp_InternalGravWave_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_InternalGravWave_Final

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

    use Constants_mod, only: &
         & Grav, ThermalExpanCoef, PI, RPlanet, RefTemp

    use gtool_history

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: z
    integer :: k
    real(DP) :: N2 
    real(DP) :: wz_Chi(lMax,0:kMax)
    integer :: lpos
    real(DP) :: KEDensAvg0

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    N2 = BruntVaisFreq**2

    xy_totDepthBasic = H0
    xy_SurfHeightN = 0d0

    ! Set the initial perturbation of velociy field
    lpos = l_nm(n, m)
    wz_Chi = 0d0
    wz_Chi(lpos,:) = 1d0
    wz_Chi(lpos,:) = wz_Chi(lpos,:)/sqrt(2d0)

    xyz_UN = xya_GradLon_wa(wz_Chi)/RPlanet
    xyz_VN = xya_GradLat_wa(wz_Chi)/RPlanet
    
    do k=0, kMax
       z = H0*g_Sig(k)
       xyz_UN(:,:,k) = xyz_UN(:,:,k)*cos(l*PI*z/H0)
       xyz_VN(:,:,k) = xyz_VN(:,:,k)*cos(l*PI*z/H0)
       z_PTempBasic(k) = RefTemp + N2/(ThermalExpanCoef*Grav) * z
    end do

    KEDensAvg0 = AvrLonLat_xy(xyz_UN(:,:,1)**2 + xyz_VN(:,:,1)**2)
    xyz_UN = sqrt(KEDensAvg/KEDensAvg0)*xyz_UN
    xyz_VN = sqrt(KEDensAvg/KEDensAvg0)*xyz_VN

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
         & BruntVaisFreq, m, n, l

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
    call MessageNotify( 'M', module_name, " BruntVaisFreq     = %f [sec]", d=(/ BruntVaisFreq /))
    call MessageNotify( 'M', module_name, " Parameters of perturbation: m=%d,n=%d,l=%d", i=(/ m,n,l /))

    
  end subroutine read_expConfig


end module Exp_InternalGravWave_mod
