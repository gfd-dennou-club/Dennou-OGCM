!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Exp_APEOGCircI98BC_mod 
  
  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING
  
  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use SpmlUtil_mod, only: &
       & AvrLonLat_xy, xy_w, w_xy

  use UnitConversion_mod, only: &
       & degC2K
  
  use DOGCM_Admin_Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefDens

  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, kMax, nMax, lMax,        &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,  &
       & xyz_Lat, xyz_Lon, xy_Topo,           &
       & DOGCM_Admin_Grid_UpdateVCoord, xyz_Z

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
  character(*), parameter:: module_name = 'DOGCM_Exp_APEOGCircI98BC_mod' !< Module Name
  integer :: expCaseNum

  real(DP), parameter :: SeaWaterFusionPt = 271.35d0
  
contains

  !>
  !!
  !!
  subroutine DOGCM_Exp_Init(configNmlFile)
    
    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile

    ! 実行文; Executable statements
    !

    call read_expConfig(configNmlFile)

  end subroutine DOGCM_Exp_Init

  !>
  !!
  !!
  subroutine DOGCM_Exp_Final

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_Exp_Final

  !> @brief 
  !!
  !!
  subroutine DOGCM_Exp_SetInitCond()
    
    ! モジュール引用; Use statements
    !
    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC, xya_SSH, xyza_H, &
         & TRCID_PTEMP, TRCID_SALT
    
    use DOGCM_Admin_Grid_mod, only: &
         & xy_Topo, xyz_Z, z_CK

    use DOGCM_Boundary_vars_mod, only: &
         & xy_SeaSfcTemp0, xy_SeaSfcSalt0,      &
         & xy_SfcHFlx0_ns, xy_SfcHFlx0_sr,      &
         & xy_DSfcHFlxDTsOcn => xy_DSfcHFlxDTs, &
         & xy_WindStressUOcn => xy_WindStressU, &
         & xy_WindStressVOcn => xy_WindStressV, &
         & xy_FreshWtFlxS0, xy_FreshWtFlx0

    
!    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: INIT_PTEMP = 280d0
    real(DP), parameter :: INIT_SALT = 35d0

    integer :: i
    integer :: j
    integer :: k
    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    !-------------------------------------------------------------
    
    ! The depth from sea surface (topography)
    xy_Topo = 5.2d03
    ! Sea surface height
    xya_SSH(:,:,:) = 0d0
    ! Update vertical coordinate
    call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_N),  & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                        & ! (in)
         & )

    !------------------------------------------------------------------------------

    ! Set initial condition for potential temperature
    
    do k=KS, KE       
       xyzaa_TRC(IS:IE,JS:JE,k,TRCID_PTEMP,:) = INIT_PTEMP
       xyzaa_TRC(IS:IE,JS:JE,k,TRCID_SALT,:)  = INIT_SALT
    end do
    
    call read_ncfile( "./SfcTemp0.nc", "SfcTemp", xy_SeaSfcTemp0 )
    do j=JS, JE
       do i=IS, IE
          xy_SeaSfcTemp0(i,j) = max(SeaWaterFusionPt, xy_SeaSfcTemp0(i,j))
       end do
    end do

    xy_SeaSfcSalt0 = INIT_SALT
    
    !-----------------------------------------------------------------------------
    
    ! Surface wind Stress

    call read_ncfile( "./TauX0.nc", "TauX", xy_WindStressUOcn )
    call read_ncfile( "./TauY0.nc", "TauY", xy_WindStressVOcn )

    xy_WindStressUOcn(IS:IE,JS:JE) =   xy_WindStressUOcn(IS:IE,JS:JE) &
         &                           - AvrLonLat_xy(xy_WindStressUOcn(IS:IE,JS:JE))

    xy_WindStressVOcn(IS:IE,JS:JE) =   xy_WindStressVOcn(IS:IE,JS:JE) &
         &                           - AvrLonLat_xy(xy_WindStressUOcn(IS:IE,JS:JE))

    ! Surface heat flux and fresh water flux (But, they wolud not be used.)

    call read_ncfile( "./FreshWtFlxS0.nc", "FreshWtFlxS", xy_FreshWtFlxS0 )

    xy_FreshWtFlxS0(IS:IE,JS:JE) = xy_w(w_xy(xy_FreshWtFlxS0(IS:IE,JS:JE)))
    
    xy_FreshWtFlxS0(IS:IE,JS:JE) =   xy_FreshWtFlxS0(IS:IE,JS:JE) &
         &                         - AvrLonLat_xy(xy_FreshWtFlxS0(IS:IE,JS:JE))
!!$
!!$    write(*,*) 'total angular momentum=', &
!!$         & AvrLonLat_xy( xy_WindStressU(IS:IE,JS:JE)*cos(xyz_Lat(IS:IE,JS:JE,KS)) )
        
!!$    write(*,*) 'total angular momentum=', &
!!$         & AvrLonLat_xy( xy_WindStressU(IS:IE,JS:JE)*cos(xyz_Lat(IS:IE,JS:JE,KS)) )

  end subroutine DOGCM_Exp_SetInitCond

  subroutine DOGCM_Exp_Do()
    
  end subroutine DOGCM_Exp_Do
  
  
  !------ private subroutines / functions -------------------------------------------------

  subroutine read_ncfile( ncfile, varname, xy_var )

    ! モジュール引用; Use statement
    !    
    use gtool_history, only: &
         & HistoryGet, HistoryGetPointer

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: ncfile
    character(*), intent(in) :: varname
    real(DP), intent(out) :: xy_var(IA,JA)
    
    ! 実行文; Executable statement
    !

    call HistoryGet( ncfile, varname, xy_var(IS:IE,JS:JE) )
    
  end subroutine read_ncfile
  
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
    namelist /Exp_APEOGCircI98BC_nml/ &
         & expCaseNum

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
       read( unit_nml, &                            ! (in)
            & nml = Exp_APEOGCircI98BC_nml, &       ! (out)
            & iostat = iostat_nml )                 ! (out)
       close( unit_nml )
    end if

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )


  end subroutine read_expConfig
  
end module DOGCM_Exp_APEOGCircI98BC_mod

