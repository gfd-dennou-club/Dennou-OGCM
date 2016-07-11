!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Exp_IGW_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify


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
  character(*), parameter:: module_name = 'DOGCM_Exp_IGW_mod' !< Module Name

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
  subroutine DOGCM_Exp_SetInitCond()
    
    !
    !
    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N
    
    use DOGCM_Admin_Grid_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & IS, IE, JS, JE, KS, KE,       &
         & xyz_Lat, xyz_Lon, xy_Topo,    &
         & DOGCM_Admin_Grid_UpdateVCoord, xyz_Z

    use DOGCM_Admin_Variable_mod, only: &
         & xyza_U, xyza_V, xyza_H, xya_SSH, xyzaa_TRC, &
         & TRCID_PTEMP

    use SpmlUtil_mod

    use DOGCM_Admin_Constants_mod, only: &
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
    real(DP) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
        
    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    N2 = BruntVaisFreq**2

    xy_Topo = H0
    xya_SSH(:,:,:) = 0d0
    call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_N),  & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                        & ! (in)
         & )
    
    ! Set the initial perturbation of velociy field
    lpos = l_nm(n, m)
    wz_Chi = 0d0
    wz_Chi(lpos,:) = 1d0
    wz_Chi(lpos,:) = wz_Chi(lpos,:)/sqrt(2d0)

    xyz_U = xya_GradLon_wa(wz_Chi)/RPlanet
    xyz_V = xya_GradLat_wa(wz_Chi)/RPlanet
    
    do k=0, kMax
       xyz_U(:,:,k) = xyz_U(:,:,k)*cos(l*PI*xyz_Z(IS:IE,JS:JE,k)/H0)
       xyz_V(:,:,k) = xyz_V(:,:,k)*cos(l*PI*xyz_Z(IS:IE,JS:JE,k)/H0)
       xyz_PTemp(:,:,k) = RefTemp + N2/(ThermalExpanCoef*Grav) * xyz_Z(IS:IE,JS:JE,k)
    end do

    KEDensAvg0 = AvrLonLat_xy(xyz_U(:,:,1)**2 + xyz_V(:,:,1)**2)
    xyza_U(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) = sqrt(KEDensAvg/KEDensAvg0)*xyz_U
    xyza_V(IS:IE,JS:JE,KS:KE,TIMELV_ID_N) = sqrt(KEDensAvg/KEDensAvg0)*xyz_V
    xyzaa_TRC(IS:IE,JS:JE,KS:KE,TRCID_PTEMP,TIMELV_ID_N) = xyz_PTemp

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
         & BruntVaisFreq, m, n, l

    ! 実行文; Executable statement
    !

    m = 0
    n = 1
    l = 1
    
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


end module DOGCM_Exp_IGW_mod
