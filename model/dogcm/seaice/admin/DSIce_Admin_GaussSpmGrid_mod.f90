!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Admin_GaussSpmGrid_mod 

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

  use DSIce_Admin_GridDef_mod, only: &
       & I_XY, I_UY, I_XV
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Admin_GaussSpmGrid_Init, DSIce_Admin_GaussSpmGrid_Final

  public :: DSIce_Admin_GaussSpmGrid_ConstructAxisInfo
  public :: DSIce_Admin_GaussSpmGrid_ConstructGrid
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Admin_GaussSpmGrid_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DSIce_Admin_GaussSpmGrid_Init()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Admin_GaussSpmGrid_Init

  !>
  !!
  !!
  subroutine DSIce_Admin_GaussSpmGrid_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Admin_GaussSpmGrid_Final

  !--------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_GaussSpmGrid_ConstructAxisInfo( &
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
    KAXIS_name = 'sig2'
    KAXIS_long_name = 'general vertical coordinate (sea ice model)'
    KAXIS_units = '1'
    KAXIS_Wt_units = '1'
    
    
  end subroutine DSIce_Admin_GaussSpmGrid_ConstructAxisInfo

  !--------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_GaussSpmGrid_ConstructGrid( &
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
    real(DP) :: x_Lon_Weight(IS:IE)
    real(DP) :: y_Lat_Weight(JS:JE)
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
            &  "DSIce_Admin_GaussSpmGrid_ConstructGrid is called before SpmlUtil_mod is initialized.")
    end if

    
    call get_SpmlGridInfo( &
         & xy_Lon_=xy_LonSpml, xy_Lat_=xy_LatSpml,                    & ! (out)
         & x_Lon_Weight_=x_Lon_Weight, y_Lat_Weight_=y_Lat_Weight     & ! (out)
         & )

    !- Center position and size of each cell ---------------------
    
    xy_Lon(IS:IE,JS:JE) = xy_LonSpml(IS:IE,JS:JE)
    xy_Lat(IS:IE,JS:JE) = xy_LatSpml(IS:IE,JS:JE)

    x_CI(IS:IE) = xy_LonSpml(IS:IE,JS) * 180d0 / PI
    y_CJ(JS:JE) = xy_LatSpml(IS,JS:JE) * 180d0 / PI

    DSig = 1d0 / dble( KM )
    do k = KS, KE
       z_CK(k) =  - 0.5d0*DSig - (k - KS)*DSig
    end do

    x_IAXIS_Weight(IS:IE) = x_Lon_Weight(IS:IE)    
    y_JAXIS_Weight(JS:JE) = y_Lat_Weight(JS:JE)

    z_KAXIS_Weight(KS:KE) = DSig
    z_CDK(KS:KE)          = DSig
    
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
    y_FJ(JS-1) = - PI/2d0
    do j = JS, JE-1
       y_FJ(j) = asin( y_JAXIS_Weight(j) + sin(y_FJ(j-1)) )
    end do
    y_FJ(JE) = PI/2d0
    y_FJ(:) = y_FJ*180d0/PI
    
    y_CDJ(JS:JE) = y_FJ(JS:JE) - y_FJ(JS-1:JE-1)

    y_FDJ(JS-1) = 0d0
    y_FDJ(JE)   = 0d0
    y_FDJ(JS:JE-1) = y_CJ(JS+1:JE) - y_CJ(JS:JE-1)
    
    !
    z_FK(KS-1) = 0d0
    do k = KS, KE
       z_FK(k) = z_FK(k-1) + z_CDK(k)
    end do
    do k = KS, KE
       z_FDK(k) = z_FDK(k-1) + z_FDK(k)
    end do

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

  end subroutine DSIce_Admin_GaussSpmGrid_ConstructGrid

  !-----------------------------------------


end module DSIce_Admin_GaussSpmGrid_mod

