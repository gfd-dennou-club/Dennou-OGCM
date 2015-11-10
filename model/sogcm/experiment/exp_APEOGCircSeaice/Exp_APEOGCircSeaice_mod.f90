!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module Exp_APEOGcircSeaice_mod

  ! モジュール引用; Use statements
  !
  use dc_types
  use dc_message

  use Constants_mod, only: &
       & StB, PI, Grav, RPlanet, Omega, RefDens

  use GridSet_mod, only: &
         & iMax, jMax, kMax, nMax, lMax, &
         & xyz_Lat, xyz_Lon

  use SpmlUtil_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_APEOGcircSeaice_Init, Exp_APEOGcircSeaice_Final
  public :: SetInitCondition

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_APEOGcircSeaice_mod' !< Module Name
  integer :: expCaseNum

  character(*), parameter :: SurfBC_DATA_DIR = "./sogcm_SurfBC/"
  
contains

  !>
  !!
  !!
  subroutine Exp_APEOGcircSeaice_Init(configNmlFile)
    
    ! 宣言文; Declare statements
    !
    character(*), intent(in) :: configNmlFile

    ! 実行文; Executable statements
    !

    call read_expConfig(configNmlFile)

  end subroutine Exp_APEOGcircSeaice_Init

  !>
  !!
  !!
  subroutine Exp_APEOGcircSeaice_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_APEOGcircSeaice_Final

  !> @brief 
  !!
  !!
  subroutine setInitCondition()

    ! モジュール引用; Use statements
    !

    use VariableSet_mod, only: &
         & z_PTempBasic, xyz_SaltN, xyz_PTempEddN, &
         & xy_totDepthBasic

    use BoundaryCondO_mod, only: &
         & xy_WindStressU, xy_WindStressV, &
         & xy_SeaSurfTemp, xy_SeaSurfSalt, &
         & xy_SurfHFlxO,    &
         & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, &
         & xy_SurfFwFlxO, &
         & xy_Wrain, xy_Wsnow, xy_Wevap
    
    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: alpha = 0d0
    real(DP) :: h0, z
    integer :: k

    integer :: i, j, m
    real(DP) :: z_PTemp(0:kMax),  z_Salt(0:kMax), xy_CosLat(0:iMax-1,jMax)

    real(DP), parameter :: RefSalt = 35d0
    real(DP), parameter :: SeaWaterFreezeTemp = -1.8d0 + 273.15d0
    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", module_name, "Set initial condition..")
    
    h0 = 5.2d03
    xy_totDepthBasic = h0
    xy_CosLat = cos(xyz_Lat(:,:,0))
    
    !
    do k=0, kMax
       z_PTempBasic(k) = eval_PTempBasic(g_Sig(k)) !-2d0
       z_Salt(k) = eval_SaltBasic(g_Sig(k))
    end do

    !
    !
    
    xy_WindStressU = construct_WindStressU_Marshall07_2(xyz_Lat(:,:,0))!- read_SurfField(trim(SurfBC_DATA_DIR)//"TauX.nc", "TauX")
    xy_WindStressV = 0d0
    xy_WindStressU = xy_w(w_xy(xy_WindStressU))
    
    xy_SeaSurfTemp = read_SurfField(trim(SurfBC_DATA_DIR)//"SurfTemp.nc", "SurfTemp")
    where(xy_SeaSurfTemp < SeaWaterFreezeTemp)
       xy_SeaSurfTemp = SeaWaterFreezeTemp
    end where
    xy_SeaSurfTemp = xy_w(w_xy(xy_SeaSurfTemp))
    
    xy_SeaSurfSalt = eval_SSSalref(xyz_Lat(:,:,0)) ! RefSalt
    xy_SeaSurfSalt = xy_w(w_xy(xy_SeaSurfSalt))
    
    xy_SWDWRFlx = read_SurfField(trim(SurfBC_DATA_DIR)//"RadSDWFLXA.nc", "RadSDWFLXA")
    
    xy_LWDWRFlx = read_SurfField(trim(SurfBC_DATA_DIR)//"RadLDWFLXA.nc", "RadLDWFLXA")
    xy_LatentDWHFlx = - read_SurfField(trim(SurfBC_DATA_DIR)//"Latent.nc", "Latent")
    xy_SensDWHFlx = - read_SurfField(trim(SurfBC_DATA_DIR)//"Sens.nc", "Sens")
    
!    xy_SurfFwFlxO = read_SurfField(trim(SurfBC_DATA_DIR)//"Fw.nc", "Fw")
    xy_Wevap = read_SurfField(trim(SurfBC_DATA_DIR)//"Evap.nc", "Evap")
    xy_Wrain = read_SurfField(trim(SurfBC_DATA_DIR)//"Rain.nc", "Rain")
    xy_Wsnow = read_SurfField(trim(SurfBC_DATA_DIR)//"Snow.nc", "Snow")

    do k=0, kMax
       xyz_PTempEddN(:,:,k) = 0d0
       xyz_SaltN(:,:,k) = z_Salt(k)
    end do

!    xy_SeaSurfTemp = 280d0
!    xy_SeaSurfSalt = 35d0
!    z_PTempBasic(:) = 280d0
!    xyz_SaltN(:,:,:) = 35d0
    xyz_SaltN(:,:,0) = xy_SeaSurfSalt
    xyz_PTempEddN(:,:,0) = 0d0

    ! Consider the insulation effect due to sea ice.
!!$    xy_WindStressU = xy_WindStressU &
!!$         & *0.5d0*(1d0 + tanh(  (xyz_Lat(:,:,0) + ICELAT)/TRANSITION_LATWIDTH )) &
!!$         & *0.5d0*(1d0 - tanh(  (xyz_Lat(:,:,0) - ICELAT)/TRANSITION_LATWIDTH ))
         

    write(*,*) 'total angular momentum=', AvrLonLat_xy( xy_WindStressU*cos(xyz_Lat(:,:,0)) )

    write(*,*) 'avg net heat flux', AvrLonLat_xy( xy_SurfHFlxO ), "[W/m2]"    
    write(*,*) 'avg freshWater flux', AvrLonLat_xy( xy_SurfFwFlxO )*86400d0, "[m/day]"
!!$    write(*,*) 'Adjust fluxes such that their global integration equal to zero.'

  end subroutine setInitCondition

  function adjust_GlobalIntegFlx(xy_flux) result(xy_fluxAdjust)
    real(DP), intent(in) :: xy_flux(0:iMax-1,jMax)
    real(DP) :: xy_fluxAdjust(0:iMax-1,jMax)

    xy_fluxAdjust = xy_flux - AvrLonLat_xy(xy_flux)
    
  end function adjust_GlobalIntegFlx
  
  function read_SurfField(ncFileName, varName) result(xy)

    use gtool_history, only: HistoryGet
    
    character(*), intent(in) :: ncFileName
    character(*), intent(in) :: varName
    real(DP) :: xy(0:iMax-1, jMax)

    integer :: i
    
    call HistoryGet(ncFileName, varName, xy(0,:))
    forAll(i=1:iMax-1) xy(i,:) = xy(0,:)
    
  end function read_SurfField

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
  
    function eval_SSSalref(xy_lat) result(xy_SSSalref)
      
      real(DP), intent(in) :: xy_lat(0:iMax-1,jMax)
      real(DP) :: xy_SSSalref(0:iMax-1,jMax)

!!$      real(DP), parameter :: coef(0:16) = &
!!$         & (/   34.6759879091764,     1.10532298396405,  0.0763166880123448,   -0.582942261643657, &
!!$         &     -0.288887879152796,   0.162029019072606,   0.196975818352802,  0.00908477735848993, &
!!$         &     -0.0844736709221092, -0.0375313739082941,  0.0231237387527393,  0.0178066954031123, &
!!$         &     -0.0130241816472578, -0.0209568343944606, -0.00517963837567204, 0.0078795138270612, 0.00428900513312646 /)


      real(DP), parameter :: coef(0:16) = &
           & (/34.7637464489794,1.051234697145,0.14236961978175,-0.619301930415567,-0.269350521029104, &
           & 0.173395861002524,0.203604009069993,0.0695294147072904,-0.0180873718247458,-0.0330484535144368, &
           & -0.00905946324146973,-0.0162817975967151,-0.00750710204111945,0.00760675898215889,0.00916853734527952, &
           & -0.00346164327996474,-0.0100834132695167 /)

!!$           (/34.7825375509604,1.03842433663076,0.0883273377454709,-0.520140557430786,-0.315009600760316, &
!!$           & 0.127365359039471,0.250630657680549,0.0877657409636654,-0.0466327816722313,-0.04174176874417,-&
!!$           & 0.00177200676103769,-0.000567452261260674,-0.0115915645347134,-0.00438358376247672,0.00866997545821757, &
!!$           & 0.00332368662059798,-0.0064656444960895 /)
      
!!$      

      
      
      integer :: m

      xy_SSSalref = 0d0
      do m=0, 16
         xy_SSSalref = xy_SSSalref + coef(m)*cos(2d0*m*xy_lat)
      end do

    end function eval_SSSalref

    function construct_WindStressU_Marshall07_2(xy_lat) result(windStressU)
      real(DP), intent(in) :: xy_lat(:,:)
      real(DP) :: windStressU(size(xy_lat,1), size(xy_lat,2))

      real(DP), parameter :: coef(0:16) = &
           & (/ 0.0306391950853468,-0.053818490401627,-0.0545487529629642,0.0308100670974264,0.0253979991710819, &
           &    0.0133005673008549,-0.00338827649597324,0.00184347119932358,-0.000396782059666333,0.00222160982930569, &
           &    -0.00181121235994517,-0.000693437465271213,-0.00138827649597365,0.000747066068848734,-0.00047153083167714, &
           &    -0.000366847438814569,0.00050062006733265 /)

      integer :: m

      windStressU = 0d0
      do m=0, 16
         windStressU = windStressU + coef(m)*cos(2d0*m*xy_lat)
      end do

    end function construct_WindStressU_Marshall07_2
    
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

end module Exp_APEOGcircSeaice_mod

