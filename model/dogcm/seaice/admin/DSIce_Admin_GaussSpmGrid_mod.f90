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

  !------------------------------------

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

  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_GaussSpmGrid_ConstructGrid( &
       & x_CI, x_FI, x_IAXIS_Weight,                    & ! (out)
       & y_CJ, y_FJ, y_JAXIS_Weight,                    & ! (out)
       & z_CK, z_FK, z_KAXIS_Weight,                    & ! (out)
       & xy_Lon, xy_Lat,                                & ! (out)
       & IS, IE, IA, IM, JS, JE, JA, JM, KS, KE, KA, KM & ! (in)
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
    integer, intent(in) :: JS
    integer, intent(in) :: JE
    integer, intent(in) :: JA
    integer, intent(in) :: JM    
    integer, intent(in) :: KS
    integer, intent(in) :: KE
    integer, intent(in) :: KA
    integer, intent(in) :: KM
    
    real(DP), intent(out) :: x_CI(IA)
    real(DP), intent(out) :: x_FI(IA)
    real(DP), intent(out) :: x_IAXIS_Weight(IA)

    real(DP), intent(out) :: y_CJ(JA)
    real(DP), intent(out) :: y_FJ(JA)
    real(DP), intent(out) :: y_JAXIS_Weight(JA)

    real(DP), intent(out) :: z_CK(KA)
    real(DP), intent(out) :: z_FK(KA)
    real(DP), intent(out) :: z_KAXIS_Weight(KA)

    real(DP), intent(out) :: xy_Lon(IA,JA)
    real(DP), intent(out) :: xy_Lat(IA,JA)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_LonSpml(IS:IE,JS:JE)
    real(DP) :: xy_LatSpml(IS:IE,JS:JE)
    real(DP) :: x_Lon_Weight(IS:IE)
    real(DP) :: y_Lat_Weight(JS:JE)

    real(DP) :: DSig
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

    xy_Lon(IS:IE,JS:JE) = xy_LonSpml(IS:IE,JS:JE)
    xy_Lat(IS:IE,JS:JE) = xy_LatSpml(IS:IE,JS:JE)

    x_CI(IS:IE) = xy_LonSpml(IS:IE,JS) * 180d0 / PI
    x_FI = 0d0    
    y_CJ(JS:JE) = xy_LatSpml(IS,JS:JE) * 180d0 / PI
    y_FJ = 0d0

    DSig = 1d0 / dble( KM )
    do k = KS, KE
       z_CK(k) =  - 0.5d0*DSig - (k - KS)*DSig
    end do
    z_FK = 0d0

    x_IAXIS_Weight(IS:IE) = x_Lon_Weight(IS:IE)
    y_JAXIS_Weight(JS:JE) = y_Lat_Weight(JS:JE)
    z_KAXIS_Weight(KS:KE) = DSig

  end subroutine DSIce_Admin_GaussSpmGrid_ConstructGrid

end module DSIce_Admin_GaussSpmGrid_mod

