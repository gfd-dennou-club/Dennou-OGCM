!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_WindDrivenCirculation_mod 

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
  public :: Exp_WindDrivenCirculation_Init, Exp_WindDrivenCirculation_Final
  public :: SetInitCondition

  public :: construct_WindStressU_Marshall07

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_WindDrivenCirculation_mod' !< Module Name
  integer :: expCaseNum

contains

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirculation_Init(configNmlFile)
    
    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile

    ! 実行文; Executable statements
    !

    call read_expConfig(configNmlFile)

  end subroutine Exp_WindDrivenCirculation_Init

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirculation_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_WindDrivenCirculation_Final

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

    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")

    h0 = 5.2d03
    xy_totDepthBasic = h0
    xy_SurfHeightN = 0d0

    xy_WindStressU = construct_WindStressU_Marshall07(xyz_Lat(:,:,0))
!!$    xy_WindStressU = construct_WindStressU_analysticFunc()
    xy_WindStressV = 0d0

    write(*,*) 'total angular momentum=', AvrLonLat_xy( xy_WindStressU*cos(xyz_Lat(:,:,1)) )

    do k=0, kMax
       z_PTemp(k) = eval_PTempBasic(g_Sig(k))
    end do
    
    TempAvg = IntSig_BtmToTop(z_PTemp)
    z_PTempBasic = TempAvg
    do k=0, kMax
       xyz_PTempEddN(:,:,k) = 0d0!z_PTemp(k) - TempAvg
    end do

!!$write(*,*) "-- WindStressU --"
!!$write(*,*) xy_WindStressU(1,:)
!!$stop
!!$write(*,*) "-- PTempBasic --"
!!$write(*,*) z_PTempBasic
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

end module Exp_WindDrivenCirculation_mod

