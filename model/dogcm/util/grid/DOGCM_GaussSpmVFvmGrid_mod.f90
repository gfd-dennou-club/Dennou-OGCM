!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_GaussSpmVFvmGrid_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/SeaIce

  use DOGCM_Admin_Constants_mod, only: &
       & RPlanet

  use GridIndex_mod, only: &
       & I_XY, I_UY, I_XV
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_GaussSpmVFvmGrid_Init, DOGCM_GaussSpmVFvmGrid_Final

  ! 公開変数
  ! Public variables
  !
  
  integer, public :: nMax
  integer, public :: lMax
  integer, public :: tMax

  integer, public :: iMax
  integer, public :: jMax
  integer, public :: jMaxGlobe
  integer, public :: kMax
  
  public :: DOGCM_GaussSpmVFvmGrid_ConstructAxisInfo
  public :: DOGCM_GaussSpmVFvmGrid_ConstructGrid
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_GaussSpmVFvmGrid_mod' !< Module Name

  logical :: isInitialzed = .false.
  
contains

  !>
  !!
  !!
  subroutine DOGCM_GaussSpmVFvmGrid_Init(configNmlFileName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! local variable
    !
    
    ! 実行文; Executable statements
    !

    ! Read the path to grid data file from namelist.
    call read_nmlData(configNmlFileName)

    isInitialzed = .true.
    
  end subroutine DOGCM_GaussSpmVFvmGrid_Init

  !>
  !!
  !!
  subroutine DOGCM_GaussSpmVFvmGrid_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_GaussSpmVFvmGrid_Final

  !------------------------------------

  subroutine DOGCM_GaussSpmVFvmGrid_ConstructAxisInfo( &
       & IAXIS_name, IAXIS_long_name, IAXIS_units, IAXIS_Wt_units,    & ! (out)
       & IS, IE, IA, IHALO,                                           & ! (out)
       & JAXIS_name, JAXIS_long_name, JAXIS_units, JAXIS_Wt_units,    & ! (out)
       & JS, JE, JA, JHALO,                                           & ! (out)
       & KAXIS_name, KAXIS_long_name, KAXIS_units, KAXIS_Wt_units,    & ! (out)
       & KS, KE, KA, KHALO,                                           & ! (out)
       & IM, JM, KM                                                   & ! (in)
       & )

    
    ! 宣言文; Declaration statement
    !
    character(TOKEN), intent(out) :: IAXIS_name
    character(STRING), intent(out) :: IAXIS_long_name
    character(TOKEN), intent(out) :: IAXIS_units
    character(TOKEN), intent(out) :: IAXIS_Wt_units
    integer, intent(out) :: IS
    integer, intent(out) :: IE
    integer, intent(out) :: IA
    integer, intent(out) :: IHALO

    character(TOKEN), intent(out) :: JAXIS_name
    character(STRING), intent(out) :: JAXIS_long_name
    character(TOKEN), intent(out) :: JAXIS_units
    character(TOKEN), intent(out) :: JAXIS_Wt_units
    integer, intent(out) :: JS
    integer, intent(out) :: JE
    integer, intent(out) :: JA
    integer, intent(out) :: JHALO

    character(TOKEN), intent(out) :: KAXIS_name
    character(STRING), intent(out) :: KAXIS_long_name
    character(TOKEN), intent(out) :: KAXIS_units    
    character(TOKEN), intent(out) :: KAXIS_Wt_units
    integer, intent(out) :: KS
    integer, intent(out) :: KE
    integer, intent(out) :: KA
    integer, intent(out) :: KHALO
    
    integer, intent(in) :: IM
    integer, intent(in) :: JM
    integer, intent(in) :: KM
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    ! i axis
    
    IHALO = 1
    IA = IM + 2*IHALO
    IS = IHALO + 1
    IE = IS + IM - 1
    IAXIS_name = 'lon'
    IAXIS_long_name = 'longitude'
    IAXIS_units = 'degree_east'
    IAXIS_Wt_units = 'radian'

    ! j axis
    
    JHALO = 1
    JA = JM + 2*JHALO
    JS = JHALO + 1
    JE = JS + JM - 1
    JAXIS_name = 'lat'
    JAXIS_long_name = 'latitude'
    JAXIS_units = 'degree_north'
    JAXIS_Wt_units = 'radian'

    ! k axis
    
    KHALO = 1
    KA = KM + 2*KHALO
    KS = KHALO + 1
    KE = KS + KM - 1
    KAXIS_name = 'sig'
    KAXIS_long_name = 'general vertical coordinate'
    KAXIS_units = '1'
    KAXIS_Wt_units = '1'
    
    ! Set truncated wave number

    iMax = IM
    
    jMax = JM
    jMaxGlobe = jMax

    kMax = KM 
    tMax = kMax

#ifdef DSOGCM_MODE_AXISYM
    if (iMax /= 1) then
       call MessageNotify("E", module_name, &
            & "The number of grid points in longitude must be 1 in DSOGCM_MODE_AXISYM")
    end if
#endif

#ifdef DSOGCM_MODE_AXISYM
    if (nMax == -1) nMax = ( 2*jMaxGlobe - 1 )/ 3
    lMax = nMax + 1 
#else
    if (nMax == -1) nMax = ( iMax - 1 )/ 3
    lMax = ( nMax + 1 )**2
#endif

    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '    nmax        = %d', i = (/   nmax  /) )
    call MessageNotify( 'M', module_name, '    iMax        = %d', i = (/   iMax  /) )
    call MessageNotify( 'M', module_name, '    jMaxGlobe   = %d', i = (/   jMaxGlobe /) )
    call MessageNotify( 'M', module_name, '    jMax        = %d', i = (/   jMax  /) )
    call MessageNotify( 'M', module_name, '    kMax        = %d', i = (/   kmax  /) )
    call MessageNotify( 'M', module_name, '    lMax        = %d', i = (/   lmax  /) )
    call MessageNotify( 'M', module_name, '    tMax        = %d', i = (/   tMax  /) )

    call MessageNotify( 'M', module_name, ' (IM,JM,KM) = (%d,%d,%d) ', i = (/ IM,JM,KM /) )
    
  end subroutine DOGCM_GaussSpmVFvmGrid_ConstructAxisInfo

  !> @brief 
  !!
  !!
  subroutine DOGCM_GaussSpmVFvmGrid_ConstructGrid( &
       & x_CI, x_CDI, x_FI, x_FDI, x_IAXIS_Weight,      & ! (out)
       & y_CJ, y_CDJ, y_FJ, y_FDJ, y_JAXIS_Weight,      & ! (out)
       & z_CK, z_CDK, z_FK, z_FDK, z_KAXIS_Weight,      & ! (out)
       & xy_Lon, xy_Lat,                                & ! (out)
       & SCALEF_E1, SCALEF_E2,                          & ! (out)
       & IS, IE, IA, IM, IHALO,                         & ! (in)
       & JS, JE, JA, JM, JHALO,                         & ! (in)
       & KS, KE, KA, KM, KHALO                          & ! (in)
       & ) 

    ! モジュール引用; Use statements
    !
    use SpmlUtil_mod, only: &
        & isSpmlUtilInitialzed=>isInitialzed, &
        & get_SpmlGridInfo
    
    ! 宣言文; Declaration statement
    !

    integer, intent(in) :: IS
    integer, intent(in) :: IE
    integer, intent(in) :: IA
    integer, intent(in) :: IM
    integer, intent(in) :: IHALO
    
    integer, intent(in) :: JS
    integer, intent(in) :: JE
    integer, intent(in) :: JA
    integer, intent(in) :: JM    
    integer, intent(in) :: JHALO

    integer, intent(in) :: KS
    integer, intent(in) :: KE
    integer, intent(in) :: KA
    integer, intent(in) :: KM
    integer, intent(in) :: KHALO
    
    real(DP), intent(out) :: x_CI(IA)
    real(DP), intent(out) :: x_CDI(IA)
    real(DP), intent(out) :: x_FI(IA)
    real(DP), intent(out) :: x_FDI(IA)
    real(DP), intent(out) :: x_IAXIS_Weight(IA)

    real(DP), intent(out) :: y_CJ(JA)
    real(DP), intent(out) :: y_CDJ(JA)
    real(DP), intent(out) :: y_FJ(JA)
    real(DP), intent(out) :: y_FDJ(JA)
    real(DP), intent(out) :: y_JAXIS_Weight(JA)

    real(DP), intent(out) :: z_CK(KA)
    real(DP), intent(out) :: z_CDK(KA)
    real(DP), intent(out) :: z_FK(KA)
    real(DP), intent(out) :: z_FDK(KA)
    real(DP), intent(out) :: z_KAXIS_Weight(KA)

    real(DP), intent(out) :: xy_Lon(IA,JA)
    real(DP), intent(out) :: xy_Lat(IA,JA)

    real(DP), intent(out) :: SCALEF_E1(IA,JA,4)
    real(DP), intent(out) :: SCALEF_E2(IA,JA,4)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_LonSpml(IS:IE,JS:JE)
    real(DP) :: xy_LatSpml(IS:IE,JS:JE)
    real(DP) :: z_SigSpml(0:kMax)
    real(DP) :: x_Lon_Weight(IS:IE)
    real(DP) :: y_Lat_Weight(JS:JE)
    real(DP) :: z_Sig_Weight(KS:KE)
    real(DP) :: v_Lat
    
    real(DP) :: DSig

    real(DP) :: Mat(JS-1:JE,JS-1:JE)
    real(DP) :: B(JS-1:JE)
    real(DP) :: work(JM + JM+1)
    integer :: info
    
    integer :: i
    integer :: j
    integer :: k
    
    real(DP), parameter :: PI = acos(-1d0)

    ! 実行文; Executable statement
    !


    if( .not. isSpmlUtilInitialzed() ) then
       call MessageNotify('E', module_name, &
            &  "DOGCM_GaussSpmVFvmGrid_ConstructGrid is called before SpmlUtil_mod is initialized.")
    end if


    call get_SpmlGridInfo( &
         & xy_Lon_=xy_LonSpml, xy_Lat_=xy_LatSpml, z_Sig_=z_SigSpml,  &
         & x_Lon_Weight_=x_Lon_Weight, y_Lat_Weight_=y_Lat_Weight, z_Sig_Weight_=z_Sig_Weight &
         & )

    !- Center position and size of each horizontal cells  ---------------------
    
    xy_Lon(IS:IE,JS:JE) = xy_LonSpml(IS:IE,JS:JE)
    xy_Lat(IS:IE,JS:JE) = xy_LatSpml(IS:IE,JS:JE)

    x_CI(IS:IE) = xy_LonSpml(IS:IE,JS) * 180d0 / PI
    y_CJ(JS:JE) = xy_LatSpml(IS,JS:JE) * 180d0 / PI
    
    x_IAXIS_Weight(IS:IE) = x_Lon_Weight(IS:IE)    
    y_JAXIS_Weight(JS:JE) = y_Lat_Weight(JS:JE)
    
    !- Fill coordinate information at horizontal halo region

    x_CI(1:IS-1)  = x_CI(IE-IHALO+1:IE)
    x_CI(IE+1:IA) = x_CI(IS:IS+IHALO-1)

    y_CJ(1:JS-1)  = y_CJ(JE-JHALO+1:JE)
    y_CJ(JE+1:JA) = y_CJ(JS:JS+JHALO-1)
    
    !- Center position and size of each face cell -----------------

    !
    do i = IS-1, IE
       x_FI(i) = 0.5d0 * (x_CI(i) + x_CI(i+1))
    end do

    x_CDI(IS:IE) = x_FI(IS:IE) - x_FI(IS-1:IE-1)
    if (IM == 1) then
       x_CDI(IS:IE) = 360d0
    end if
    x_FDI(IS-1:IE) = x_CI(IS:IE+1) - x_CI(IS-1:IE)
       
    !
    y_FJ(JS-1) = - 90d0
    y_FJ(JE)   =   90d0

    Mat = 0d0; B = 0d0
    do j = JS, JE-1
       Mat(j,j:j+1) = 0.5d0
       B(j) = y_CJ(j+1) - y_CJ(j)
    end do
    Mat(JS-1,JS) = 0.5d0; B(JS-1) = y_CJ(JS) - y_FJ(JS-1)
    Mat(JE  ,JE) = 0.5d0; B(JE  ) = y_FJ(JE) - y_CJ(JE)
    
    call  DGELS('N', JM+1, JM, 1, Mat, JM+1, B, JM+1, work, size(work), info)

    do j = JS, JE-1
       y_FJ(j) = 0.5d0 * (y_CJ(j) + y_CJ(j+1))
    end do

    y_CDJ(JS:JE) = y_FJ(JS:JE) - y_FJ(JS-1:JE-1)

    y_FDJ(JS-1) = 0d0
    y_FDJ(JE)   = 0d0
    y_FDJ(JS:JE-1) = y_CJ(JS+1:JE) - y_CJ(JS:JE-1)

    
    !
    z_FK(KS-1:KE) = z_SigSpml(:)
!    z_FK(:) = (/ ( - 1d0/dble(KE-KS+ 1)*dble(k-(KS-1)), k=KS-1, KE) /)
    do k = KS, KE
       z_CK(k)  = 0.5d0*(z_FK(k-1) + z_FK(k))
    end do
    z_CK(KS-1) = z_FK(KS-1) + 0.5d0*z_CDK(KS)
    z_CK(KE+1) = z_FK(KE) - 0.5d0*z_CDK(KE)

    do k = KS, KE
       z_CDK(k) = z_FK(k-1) - z_FK(k)
    end do
    z_CDK(KS-1) = z_CDK(KS)
    z_CDK(KE+1) = z_CDK(KE)
    z_KAXIS_Weight(KS:KE) = z_CDK(KS:KE)
    
    do k = KS-1, KE
       z_FDK(k) = z_CK(k) - z_CK(k+1)
    end do

!!$    write(*,*) z_CK(KS:KE)
!!$    write(*,*) z_FK(KS-1:KE)
!!$    write(*,*) z_CDK(KS-1:KE+1)
!!$    write(*,*) z_FDK(KS-1:KE)    
!!$    stop

    !-----------------------------------------------------------------------

    !$omp parallel do private(i, v_Lat)
    do j = JS, JE
       do i = IS, IE
          SCALEF_E1(i,j,I_XY) = RPlanet*cos(xy_Lat(i,j)) * PI/180d0
          SCALEF_E1(i,j,I_UY) = RPlanet*cos(xy_Lat(i,j)) * PI/180d0
          v_Lat = y_FJ(j) * PI/180d0          
          SCALEF_E1(i,j,I_XV) = RPlanet*cos(v_Lat) * PI/180d0

          SCALEF_E2(i,j,I_XY) = RPlanet * PI/180d0
          SCALEF_E2(i,j,I_UY) = RPlanet * PI/180d0
          SCALEF_E2(i,j,I_XV) = RPlanet * PI/180d0
       end do
    end do
    
    !----------------------------------------------------------------------
    
!!$    write(*,*) "-------------------"
!!$    write(*,*) "x_CI=", x_CI(1:IA)
!!$    write(*,*) "x_CDI=", x_CDI(1:IA)
!!$    write(*,*) "x_FI=", x_FI(1:IA)
!!$    write(*,*) "x_FDI=", x_FDI(1:IA)
!!$    write(*,*) "-------------------"
!!$    write(*,*) "y_CJ=", y_CJ(1:JA)
!!$    write(*,*) "y_CDJ=", y_CDJ(1:JA), sum(y_CDJ(JS:JE))
!!$    write(*,*) "y_FJ=", y_FJ(1:JA)
!!$    write(*,*) "y_FDJ=", y_FDJ(1:JA)
!!$    stop

  end subroutine DOGCM_GaussSpmVFvmGrid_ConstructGrid
  
  !------------------------------------------------------

!--------------------------------------------------------
    
  subroutine read_nmlData( configNmlFileName )

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

    integer :: NM

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /grid_spm_nml/ &
         & NM

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    
    NM = -1
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &              ! (in)
            & nml = grid_spm_nml, &   ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    nMax = NM
    
    ! 印字 ; Print
    !
    
  end subroutine read_nmlData
  
end module DOGCM_GaussSpmVFvmGrid_mod
