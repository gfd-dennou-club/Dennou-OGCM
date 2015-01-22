!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_APEOGCirc_mod

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message

  use Constants_mod, only: &
       & Grav, PI, RPlanet, Omega, RefDens

  use GridSet_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & xyz_Lat, xyz_Lon

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_APEOGCirc_Init, Exp_APEOGCirc_Final
  public :: SetInitCondition

  public :: construct_WindStressU_Marshall07

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_APEOGCirc_mod' !< Module Name
  integer :: expCaseNum

contains

  !>
  !!
  !!
  subroutine Exp_APEOGCirc_Init(configNmlFile)
    
    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile

    ! 実行文; Executable statements
    !

    call read_expConfig(configNmlFile)

  end subroutine Exp_APEOGCirc_Init

  !>
  !!
  !!
  subroutine Exp_APEOGCirc_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_APEOGCirc_Final

  !> @brief 
  !!
  !!
  subroutine setInitCondition()
    
    !
    !

    use VariableSet_mod

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
real(DP) :: w(lMax)

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 5.2d03
    xy_totDepthBasic = h0
    xy_SurfHeightN = 0d0

    xy_WindStressU = construct_WindStressU_Marshall07(xyz_Lat(:,:,0))
!!$    xy_WindStressU = construct_WindStressU_analysticFunc()
    xy_WindStressV = 0d0

    write(*,*) 'total angular momentum=', AvrLonLat_xy( xy_WindStressU*cos(xyz_Lat(:,:,0)) )

    do k=0, kMax
       z_PTempBasic(k) = eval_PTempBasic(g_Sig(k))
    end do
    xy_SeaSurfTemp = eval_SSTref(xyz_Lat(:,:,0))
    xy_SeaSurfSalt = eval_SSSalref(xyz_Lat(:,:,0))

    do k=0, kMax
       xyz_PTempEddN(:,:,k) = 0d0
       xyz_SaltN(:,:,k) = eval_SaltBasic(g_Sig(k))
    end do

!!$write(*,*) "-- WindStressU ------------"
!!$w = w_xy(xy_WindStressU*cos(xy_lat))
!!$do m=1, lMax
!!$   write(*,*) m, nm_l(m), "*", w(m)
!!$end do
!!$stop
!!$write(*,*) "-- PTempBasic ------------"
!!$write(*,*) z_PTempBasic
!!$stop
!!$write(*,*) "-- SaltBasic ------------"
!!$do k=0,kMax
!!$write(*,*) k, g_Sig(k), eval_SaltBasic(g_Sig(k))
!!$end do
!!$stop
!!$write(*,*) "-- SST reference ---------"
!!$do j=1, jMax
!!$   write(*,*) 180d0*xyz_Lat(0,j,0)/PI, xy_SeaSurfTemp(0,j)
!!$end do
!!$stop

  end subroutine setInitCondition
  
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

    function eval_SaltBasic(sig) result(z)
      real(DP), intent(in) :: sig
      real(DP) :: z

      real(DP), parameter :: coef(0:12) = &
         & (/ 34.6744042104524,   -0.16362182596767,   -0.161989120107287,    -0.1149384812899,  -0.0898883131255983,   &
         &   -0.0878779309300037, -0.0835048504873238,  -0.0726172203869132,  -0.0565830160985759, -0.0458673042288947, &
         &   -0.0349035497016124, -0.0264978663170974, -0.0170818216519464 /)

      integer :: m

      z = coef(0)
      do m=1, 12
         z = z + coef(m)*cos(m*acos(1d0+2*sig))
      end do

    end function eval_SaltBasic

    function eval_SSTref(xy_lat) result(xy_SSTref)
      
      use UnitConversion_mod, only: degC2K

      real(DP), intent(in) :: xy_lat(0:iMax-1,jMax)
      real(DP) :: xy_SSTref(0:iMax-1,jMax)

      real(DP), parameter :: coef(0:16) = &
         & (/  10.7493885142121,     17.1560785291747,   3.11530651540079,   -2.87419686348211, -0.937231317140739, &
         &     0.82804929538719, -0.00955630487540587, -0.572708108007084, 0.00806924122762969,  0.234556304923045, &
         &   -0.117106128617759,   -0.186994434602621, 0.0237770284564268,  0.0561637234381364, -0.051869301436776, &
         & -0.00669052658301952, 0.0754931613076543 /)
      integer :: m

      xy_SSTref = 0d0
      do m=0, 16
         xy_SSTref = xy_SSTref + coef(m)*cos(2d0*m*xy_lat)
      end do

      xy_SSTref = degC2K(xy_SSTref)

    end function eval_SSTref

    function eval_SSSalref(xy_lat) result(xy_SSSalref)
      
      real(DP), intent(in) :: xy_lat(0:iMax-1,jMax)
      real(DP) :: xy_SSSalref(0:iMax-1,jMax)

      real(DP), parameter :: coef(0:16) = &
         & (/   34.6759879091764,     1.10532298396405,  0.0763166880123448,   -0.582942261643657, &
         &     -0.288887879152796,   0.162029019072606,   0.196975818352802,  0.00908477735848993, &
         &     -0.0844736709221092, -0.0375313739082941,  0.0231237387527393,  0.0178066954031123, &
         &     -0.0130241816472578, -0.0209568343944606, -0.00517963837567204, 0.0078795138270612, 0.00428900513312646 /)

      integer :: m

      xy_SSSalref = 0d0
      do m=0, 16
         xy_SSSalref = xy_SSSalref + coef(m)*cos(2d0*m*xy_lat)
      end do

    end function eval_SSSalref

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

end module Exp_APEOGCirc_mod

