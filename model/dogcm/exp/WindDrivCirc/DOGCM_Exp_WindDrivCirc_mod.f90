!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Exp_WindDrivCirc_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types
  use dc_message

  !* Dennou-OGCM
  
  use DOGCM_Admin_Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefDens

  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, kMax, nMax, lMax, &
       & IS, IE, JS, JE, KS, KE,       &
       & xyz_Lat, xyz_Lon, xy_Topo,    &
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

  public :: construct_WindStressU_Marshall07

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Exp_WindDrivCirc_mod' !< Module Name
  integer :: expCaseNum

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
    
    !
    !
    use DOGCM_Admin_TInteg_mod, only: &
         & TIMELV_ID_N
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC, xya_SSH, xyza_H, &
         & TRCID_PTEMP

    use DOGCM_Boundary_vars_mod, only: &
         & xy_WindStressU, xy_WindStressV, &
         & xy_SeaSfcTemp

    use DOGCM_Admin_Grid_mod, only: &
         & xy_Topo, xyz_Z, z_CK
    
    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0
    real(DP) :: h0, u0, z
    integer :: k
    real(DP), parameter :: Tau0 = 0.1d0

    integer :: i, j, m
    real(DP) :: TempAvg
    real(DP) :: z_PTemp(0:kMax)


    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 5.2d03
    xy_Topo = h0
    xya_SSH(:,:,:) = 0d0
    call DOGCM_Admin_Grid_UpdateVCoord( xyza_H(:,:,:,TIMELV_ID_N),  & ! (out)
         & xya_SSH(:,:,TIMELV_ID_N)                        & ! (in)
         & )

    xy_WindStressU(IS:IE,JS:JE) = construct_WindStressU_Marshall07(xyz_Lat(IS:IE,JS:JE,KS))
!!$    xy_WindStressU = construct_WindStressU_analysticFunc()
    xy_WindStressV(IS:IE,JS:JE) = 0d0

    write(*,*) 'total angular momentum=', &
         & AvrLonLat_xy( xy_WindStressU(IS:IE,JS:JE)*cos(xyz_Lat(IS:IE,JS:JE,KS)) )

!!$    do k=0, kMax
!!$       z_PTemp (k) = eval_PTempBasic(z_CK(k))
!!$    end do
!!$    
!!$    TempAvg = IntSig_BtmToTop(z_PTemp)
!!$    xy_SeaSfcTemp = TempAvg

    do k=0, kMax
       xyzaa_TRC(:,:,k,TRCID_PTEMP,TIMELV_ID_N)  = 280d0 !z_PTemp(k)
    end do

!!$write(*,*) "-- WindStressU --"
!!$w = w_xy(xy_WindStressU*cos(xy_lat))
!!$do m=1, lMax
!!$   write(*,*) m, nm_l(m), "*", w(m)
!!$end do
!!$stop
!!$write(*,*) "-- PTempBasic --"
!!$write(*,*) z_PTempBasic
!!$stop

  end subroutine DOGCM_Exp_SetInitCond
  
  function construct_WindStressU_Marshall07(xy_lat) result(windStressU)
    real(DP), intent(in) :: xy_lat(:,:)
    real(DP) :: windStressU(size(xy_lat,1), size(xy_lat,2))

    real(DP), parameter :: coef(8) = &
         & (/ 0.0265682, -0.0784899, -0.00880389, 0.0343205, 0.0233334, &
         & 0.000641955, -0.00387676, -0.00150998 /)
    integer :: m

    windStressU = 0d0
    do m=1, size(coef)
       windStressU = windStressU + coef(m)*cos((2*m-1)*xy_lat)
    end do

  end function construct_WindStressU_Marshall07

    function construct_WindStressU_analysticFunc() result(xy)
      real(DP) :: xy(0:iMax-1,jMax)

      integer :: j
      real(DP) :: c1, c2, lat

      c2 = 0.12d0
!      c1 = (PI/48d0 + 1089d0*sqrt(3d0)/35840d0)/(20d0*PI/320d0 - 27d0*sqrt(3d0)/320d0) * c2
      c1 = (20d0*PI/320d0 - 27d0*sqrt(3d0)/320d0)/(PI/48d0 + 1089d0*sqrt(3d0)/35840d0) * c2

      xy = 0d0
      write(*,*) 'callll !!'
      do j=1, jMax
         lat = xyz_Lat(1,j,1)
         if ( abs(lat) <= PI/6d0 ) then 
            xy(:,j) = - c1*sin(3d0*lat)**2*sin(6d0*lat)**2
         else
            xy(:,j) = c2*sin(3d0*(lat - PI/6d0))**2*cos(lat)**2
         end if
      end do

write(*,*) xy(1,1:jMax)
    end function construct_WindStressU_analysticFunc


    function eval_PTempBasic(sig) result(z)
      real(DP), intent(in) :: sig
      real(DP) :: z

      real(DP), parameter :: coef(13) = &
         & (/ 277.121, 4.73219, 2.93132, 1.67006, 0.945594, 0.566825, &
         &    0.382828, 0.295956, 0.197681, 0.128343, 0.0627121, 0.0400944, -0.0106413 /)
      integer :: m

      z = 0d0
      do m=1, size(coef)
         z = z + coef(m)*cos((m-1)*acos(1d0+2*sig))
      end do

    end function eval_PTempBasic


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
       read( unit_nml, &                  ! (in)
            & nml = numexp_nml, &         ! (out)
            & iostat = iostat_nml )       ! (out)
       close( unit_nml )
    end if

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )

    
  end subroutine read_expConfig

end module DOGCM_Exp_WindDrivCirc_mod

