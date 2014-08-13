!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module SpmlUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

#ifdef DSOGCM_MODE_AXISYM

#ifdef DSOGCM_USE_SJPACK
  use wa_zonal_module_sjpack
#else
  use wa_zonal_module
#endif

#else

#ifdef DSOGCM_USE_SJPACK
  use wa_module_sjpack
#else
  use wa_module
#endif

#endif

  use at_module_omp, only: &
       & at_Initial, &
       & g_Sig => g_X, g_Sig_WEIGHT => g_X_WEIGHT, &
       & at_az => at_ag, & 
       & az_at => ag_at, &
       & IntSig_BtmToTop => Int_g, &
       & t_DSig_t => t_Dx_t, at_DSig_at => at_Dx_at, &
       & at_BoundariesGrid_NN, at_BoundariesGrid_DD, &
       & at_BoundariesGrid_ND, at_BoundariesGrid_DN

#ifdef _OPENMP
  use omp_lib
#endif  



  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SpmlUtil_Init, SpmlUtil_Final
  public :: isInitialzed
  public :: g_Sig
 
  ! Procedures for the discritization of horizontal differential or integrator with spherical harmonics spectral method
  public :: wz_AlphaOptr_xyz, w_AlphaOptr_xy, xyz_AlphaOptr_wz, xy_AlphaOptr_w
  public :: wz_Lapla2D_wz, wz_InvLapla2D_wz
  public :: w_Lapla2D_w, w_InvLapla2D_w

  ! Procedures for the discritization of vertical differential or integrator with chebyshev spectral method
  public :: xyt_DSig_xyt, wt_DSig_wt, t_DSig_t, xyt_DSigDSig_xyt, wt_DSigDSig_wt
  public :: IntSig_BtmToTop, xy_IntSig_BtmToTop_xyz, xyz_IntSig_SigToTop_xyz, w_IntSig_BtmToTop_wz

  ! Procedures for the data conversion between real and spectral space with the spectral methods.
  public :: xyz_wt, wt_xyz, wt_xyt, xyz_wz, wz_xyz, wz_wt, wt_wz, xyt_xyz, xyz_xyt
  public :: wt_VorDiv2VectorCosLat
  public :: wt_VectorCosLat2VorDiv, wz_VectorCosLat2VorDiv

  ! Procedures to statisfy the vertical boundary conditions.
  public :: apply_ZBoundaryCond

  !* Cascade
  !

  ! basic transformation
  public :: w_xy, xy_w

  ! Derivate operator
  public :: xya_GradLon_wa, xya_GradLambda_wa, xya_GradLat_wa, xya_GradMu_wa
  public :: xy_GradLon_w, xy_GradLambda_w, xy_GradLat_w, xy_GradMu_w
  public :: w_DivLambda_xy, w_DivMu_xy
  public :: l_nm, nm_l
  public :: xy_Lon, xy_Lat, x_Lon, y_Lat
  public :: xya_wa, wa_xya 
  public :: IntLonLat_xy, ya_IntLon_xya
  public :: AvrLonLat_xy, ya_AvrLon_xya
  public :: a_Interpolate_wa, Interpolate_w
  
  ! Operation for pectral analysis
  public :: nma_EnergyFromStreamfunc_wa, na_EnergyFromStreamfunc_wa
  public :: nma_EnstrophyFromStreamfunc_wa, na_EnstrophyFromStreamfunc_wa

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SpmlUtil_mod' !< Module Name

  real(DP), parameter :: SigMin = -1d0
  real(DP), parameter :: SigMax = 0d0
  real(DP), parameter :: PI = 3.1415926535897932385D0

  real(DP) :: Radius
  integer :: im      !< Number of grid points in longitude
  integer :: jm      !< Number of grid points in latitude
  integer :: km      !< Number of vertical layers
  integer :: nm      !< Maximum truncated wave number in horizontal spectral method
  integer :: tm      !< Maximum truncated wave number in vertical spectral method
  integer :: lm      !< 

  real(DP), allocatable :: tr_vIntCoefMat(:,:)
  real(DP), allocatable :: vDiffProcInvMat(:,:)

  logical :: initalizedFlag = .false.

  integer :: nThread
  integer, allocatable :: llMin(:), llMax(:)
  integer, allocatable :: ljMin(:), ljMax(:)
  integer, allocatable :: lkMin(:), lkMax(:)
  integer, allocatable :: ltMin(:), ltMax(:)
  integer, allocatable :: lxyMin(:), lxyMax(:)

  real(DP), allocatable :: xy_CosLat(:,:)

contains

  !>
  !!
  !!
  subroutine SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np)

    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: iMax, jMax, kMax, nMax, tMax
    real(DP), intent(in) :: RPlanet
    integer, intent(in), optional :: np
    
    integer :: tId

    ! 実行文; Executable statements
    !
    im = iMax; jm = jMax; km = kMax;
    nm = nMax; tm = tMax;

#ifdef DSOGCM_MODE_AXISYM
    lm = nm + 1
#else
    lm = (nm + 1)*(nm + 1)
#endif

    nThread = 1

    if( present(np) ) then
       nThread = np
       call wa_Initial(nMax, iMax, jMax, kMax+1, np)
    else
       call wa_Initial(nMax, iMax, jMax, kMax+1)
    end if

    call MessageNotify("M", module_name, "The number of thread is %d", i=(/ nThread /))
    call at_Initial(kMax, tMax, SigMin, SigMax)
    
    Radius = RPlanet

    ! Allocate the memory of work variables in this module. 
    allocate(xy_CosLat(0:im-1,jm))
    xy_CosLat = cos(xy_Lat)
    call construct_tr_vIntCoefMat()
    
    !
    call assign_IDRange(0, kMax+1, lkMin, lkMax)
    call assign_IDRange(0, tMax+1, ltMin, ltMax)
    call assign_IDRange(1, lm, llMin, llMax)
    call assign_IDRange(1, jMax, ljMin, ljMax)
    call assign_IDRange(1, iMax*jMax, lxyMin, lxyMax)

!!$    write(*,*) 'ljMin', ljMin
!!$    write(*,*) 'ljMax', ljMax
!!$
!!$    write(*,*) 'llMin', llMin(:)
!!$    write(*,*) 'llMax', llMax(:)
!!$    write(*,*) 'lxyMin', lxyMin
!!$    write(*,*) 'lxyMax', lxyMax

    !

    call MessageNotify("M", module_name, "SpmlUtil_mod have been initialized.")
    initalizedFlag = .true.

    contains
      subroutine assign_IDRange(dimFirstId, dimSize, lIdMin, lIdMax)
        integer, intent(in) :: dimFirstId
        integer, intent(in) :: dimSize
        integer, allocatable :: lIdMin(:), lIdMax(:)

        integer :: chunkSize

        allocate(lIdMin(nThread), lIdMax(nThread))
        chunkSize = dimSize/nThread
        lIdMin(1) = dimFirstId; lIdMax(1) = chunkSize
        do tId=2,nThread
           lIdMin(tId) = lIdMax(tId-1) + 1
           lIdMax(tId) = lIdMin(tId) + chunkSize - 1
        end do
        lIdMax(nThread) = dimFirstId + dimSize - 1        
      end subroutine assign_IDRange

  end subroutine SpmlUtil_Init

  !>
  !!
  !!
  subroutine SpmlUtil_Final()

    ! 実行文; Executable statements
    !

    deallocate(xy_CosLat)
    deallocate(tr_vIntCoefMat)

  end subroutine SpmlUtil_Final


  !> @brief 
  !!
  !! @return 
  !!
  logical function isInitialzed() 
    ! 実行文; Executable statement
    !
    isInitialzed = initalizedFlag
  end function isInitialzed


  !--------------- 基本変換 -----------------

    function xyz_wt(wt)

      !
      ! スペクトルデータから 3 次元格子点データへ(逆)変換する.
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ
      real(8), dimension(0:im-1,1:jm,0:km)               :: xyz_wt
      !(out) 3 次元経度緯度動径格子点データ

      xyz_wt = xya_wa(wz_wt(wt))

    end function xyz_wt

    function wt_xyz(xyz)
      !
      ! 3 次元格子点データからスペクトルデータへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in) :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)           :: wt_xyz
      !(out) 2 次元球面調和函数チェビシェフスペクトルデータ

      wt_xyz = wt_wz(wa_xya(xyz))

    end function wt_xyz

    function xyz_wz(wz)
      !
      ! 水平スペクトル・動径格子点データから 3 次元格子点データへ(逆)変換する.
      !
      real(8), dimension(lm,0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数スペクトル・動径格子点データ
      real(8), dimension(0:im-1,1:jm,0:km)               :: xyz_wz
      !(out) 3 次元経度緯度動径格子点データ


      xyz_wz = xya_wa(wz)

    end function xyz_wz

    function wz_xyz(xyz)
      !
      ! 3 次元格子データから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:km)             :: wz_xyz
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      wz_xyz = wa_xya(xyz)

    end function wz_xyz

    function wt_xyt(xyt)
      !
      ! 3 次元格子データから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:tm), intent(in)   :: xyt
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)             :: wt_xyt
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      wt_xyt = wa_xya(xyt)

    end function wt_xyt

    function wz_wt(wt)
      !
      ! スペクトルデータから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ
      real(8), dimension(lm,0:km)             :: wz_wt
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      integer :: thId

!!$      !$omp parallel private(thId)
!!$#ifdef _OPENMP
!!$      thId = omp_get_thread_num() + 1
!!$#else
!!$      thId = 1
!!$#endif
!!$      wz_wt(llMin(thId):llMax(thId),:) = az_at(wt(llMin(thId):llMax(thId),:))
!!$      !$omp end parallel 
!!$
      wz_wt = az_at(wt)

    end function wz_wt

    function wt_wz(wz)
      !
      ! 水平スペクトル・動径格子点データからスペクトルデータへ(正)変換する.
      !
      real(8), dimension(lm,0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数スペクトル・動径格子点データ
      real(8), dimension(lm,0:tm)             :: wt_wz
      !(out) 2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: thId

!!$      !$omp parallel private(thId)
!!$#ifdef _OPENMP
!!$      thId = omp_get_thread_num() + 1
!!$#else
!!$      thId = 1
!!$#endif
!!$      wt_wz(llMin(thId):llMax(thId),:) = at_az(wz(llMin(thId):llMax(thId),:))
!!$      !$omp end parallel 

      wt_wz = at_az(wz)

    end function wt_wz

    function xyt_xyz(xyz)
      real(8), dimension(0:im-1,jm,0:km) :: xyz
      real(8), dimension(0:im-1,jm,0:tm) :: xyt_xyz

      integer :: n
      real(DP) :: at_Work(im*jm, km+1)
      
      at_Work = reshape(xyz, (/im*jm, km+1/))
!!$      
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         at_Work(lxyMin(n):lxyMax(n),:) = at_az(at_Work(lxyMin(n):lxyMax(n),:))
!!$      end do
!!$
      at_Work = at_az(at_Work)
      xyt_xyz = reshape(at_Work, shape(xyt_xyz))

      
    end function xyt_xyz

    function xyz_xyt(xyt)
      real(8), dimension(0:im-1,jm,0:tm) :: xyt
      real(8), dimension(0:im-1,jm,0:km) :: xyz_xyt

      integer :: n
      real(DP) :: az_Work(im*jm, km+1)
      
      az_Work = reshape(xyt, (/im*jm, km+1/))
      
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         az_Work(lxyMin(n):lxyMax(n),:) = az_at(az_Work(lxyMin(n):lxyMax(n),:))
!!$      end do

      az_Work = az_at(az_Work)
      xyz_xyt = reshape(az_Work, shape(xyz_xyt))

    end function xyz_xyt

    !> @brief 
    !!
    !!
    subroutine wt_VorDiv2VectorCosLat(  wt_Vor, wt_Div,  & ! (in)
         & xyz_UCosLat, xyz_VCosLat                               & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt_Vor, wt_Div
      real(8), dimension(0:im-1,jm,0:km), intent(out) :: xyz_UCosLat, xyz_VCosLat 
      
      ! 局所変数
      ! Local variables
      !
      real(8), dimension(lm,0:km) :: wz_Psi, wz_Chi
      real(8), dimension(0:im-1,jm,0:km) :: xyz_Cos2Lat

      ! 実行文; Executable statement
      !

      xyz_Cos2Lat = spread(xy_CosLat**2, 3, km+1)

      wz_Psi = wz_InvLapla2D_wz( wz_wt(wt_Vor) )
      wz_Chi = wz_InvLapla2D_wz( wz_wt(wt_Div) )
      xyz_UCosLat = xyz_Cos2Lat * xyz_AlphaOptr_wz(wz_Chi, -wz_Psi)
      xyz_VCosLat = xyz_Cos2Lat * xyz_AlphaOptr_wz(wz_Psi, wz_Chi)

    end subroutine wt_VorDiv2VectorCosLat

    subroutine wt_VectorCosLat2VorDiv(  xyz_UCosLat, xyz_VCosLat,  & ! (in)
         & wt_Vor, wt_Div                                          & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz_UCosLat, xyz_VCosLat
      real(8), dimension(lm,0:tm), intent(out) :: wt_Vor, wt_Div

      ! 実行文; Executable statement
      !

      wt_Vor = wt_wz( wz_AlphaOptr_xyz(xyz_VCosLat, -xyz_UCosLat) )
      wt_Div = wt_wz( wz_AlphaOptr_xyz(xyz_UCosLat,  xyz_VCosLat) )      

    end subroutine wt_VectorCosLat2VorDiv

    subroutine wz_VectorCosLat2VorDiv(  xyz_UCosLat, xyz_VCosLat,  & ! (in)
         & wz_Vor, wz_Div                                          & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz_UCosLat, xyz_VCosLat
      real(8), dimension(lm,0:km), intent(out) :: wz_Vor, wz_Div

      ! 実行文; Executable statement
      !

      wz_Vor = wz_AlphaOptr_xyz(xyz_VCosLat, -xyz_UCosLat)
      wz_Div = wz_AlphaOptr_xyz(xyz_UCosLat,  xyz_VCosLat)

    end subroutine wz_VectorCosLat2VorDiv

  !--------------- 水平微分計算 -----------------

    function wt_DivLon_xyz(xyz)
      ! 
      ! 格子点データに発散型経度微分 1/acosφ・∂/∂λ を作用させた
      ! スペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)       :: wt_DivLon_xyz
      !(out) 発散型経度微分を作用された 2 次元スペクトルデータ

      wt_DivLon_xyz = wt_wz(wa_DivLon_xya(xyz/Radius))

    end function wt_DivLon_xyz

    function wt_DivLat_xyz(xyz)
      !
      ! 格子データに発散型緯度微分 1/acosφ・∂(f cosφ)/∂φ を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)       :: wt_DivLat_xyz
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ

      wt_DivLat_xyz = wt_wz(wa_divlat_xya(xyz/Radius))

    end function wt_DivLat_xyz

    function wz_AlphaOptr_xyz(xyz_A, xyz_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz_A
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:km)       :: wz_AlphaOptr_xyz
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ

      wz_AlphaOptr_xyz = ( wa_DivLambda_xya(xyz_A) + wa_DivMu_xya(xyz_B) )/Radius

    end function wz_AlphaOptr_xyz

    function w_AlphaOptr_xy(xy_A, xy_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm), intent(in)   :: xy_A
      real(8), dimension(0:im-1,1:jm), intent(in)   :: xy_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm)       :: w_AlphaOptr_xy
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ

      w_AlphaOptr_xy = ( w_DivLambda_xy(xy_A) + w_DivMu_xy(xy_B) )/Radius

    end function w_AlphaOptr_xy

    function xyz_AlphaOptr_wz(wz_A, wz_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(lm,0:km), intent(in)   :: wz_A
      real(8), dimension(lm,0:km), intent(in)   :: wz_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(0:im-1,jm,0:km)       :: xyz_AlphaOptr_wz
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ


      
      xyz_AlphaOptr_wz = (xya_GradLambda_wa(wz_A) + xya_GradMu_wa(wz_B))/ (Radius*spread(xy_CosLat**2,3,km+1))

    end function xyz_AlphaOptr_wz

    function xy_AlphaOptr_w(w_A, w_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(lm), intent(in)   :: w_A
      real(8), dimension(lm), intent(in)   :: w_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(0:im-1,jm)       :: xy_AlphaOptr_w
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ


      xy_AlphaOptr_w = (xy_GradLambda_w(w_A) + xy_GradMu_w(w_B))/ (Radius*xy_CosLat**2)

    end function xy_AlphaOptr_w

    function w_Lapla2D_w(w)
      ! 入力スペクトルデータにラプラシアン
      !
      !     ▽^2 =   1/a^2 cos^2φ・∂^2/∂λ^2 
      !            + 1/a^2 cosφ・∂/∂φ(cosφ∂/∂φ) 
      !
      ! を作用する.
      !
      ! スペクトルデータのラプラシアンとは, 対応する格子点データに
      ! ラプラシアンを作用させたデータのスペクトル変換のことである. 
      !
      real(8), dimension(lm), intent(in) :: w
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm)             :: w_Lapla2D_w
      !(out) ラプラシアンを作用された 2 次元スペクトルデータ

      w_Lapla2D_w = w_Lapla_w(w)/Radius**2

    end function w_Lapla2D_w

    function wz_Lapla2D_wz(wz)
      ! 入力スペクトルデータにラプラシアン
      !
      !     ▽^2 =   1/a^2 cos^2φ・∂^2/∂λ^2 
      !            + 1/a^2 cosφ・∂/∂φ(cosφ∂/∂φ) 
      !
      ! を作用する.
      !
      ! スペクトルデータのラプラシアンとは, 対応する格子点データに
      ! ラプラシアンを作用させたデータのスペクトル変換のことである. 
      !
      real(8), dimension(lm,0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm,0:km)             :: wz_Lapla2D_wz
      !(out) ラプラシアンを作用された 2 次元スペクトルデータ

      wz_Lapla2D_wz = wa_Lapla_wa(wz)/Radius**2

    end function wz_Lapla2D_wz

    function wz_InvLapla2D_wz(wz)
      ! 入力スペクトルデータにラプラシアン
      !
      !     ▽^2 =   1/a^2 cos^2φ・∂^2/∂λ^2 
      !            + 1/a^2 cosφ・∂/∂φ(cosφ∂/∂φ) 
      ! 
      ! の逆 ∇^-2 を作用する.
      !
      ! スペクトルデータのラプラシアンとは, 対応する格子点データに
      ! ラプラシアンを作用させたデータのスペクトル変換のことである. 
      !
      real(8), dimension(lm,0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm,0:km)             :: wz_InvLapla2D_wz
      !(out) ラプラシアンを作用された 2 次元スペクトルデータ

      wz_InvLapla2D_wz = wa_LaplaInv_wa(wz) * Radius**2

    end function wz_InvLapla2D_wz


    function w_InvLapla2D_w(w)
      ! 入力スペクトルデータにラプラシアン
      !
      !     ▽^2 =   1/a^2 cos^2φ・∂^2/∂λ^2 
      !            + 1/a^2 cosφ・∂/∂φ(cosφ∂/∂φ) 
      ! 
      ! の逆 ∇^-2 を作用する.
      !
      ! スペクトルデータのラプラシアンとは, 対応する格子点データに
      ! ラプラシアンを作用させたデータのスペクトル変換のことである. 
      !
      real(8), dimension(lm,0:km), intent(in) :: w
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm)             :: w_InvLapla2D_w
      !(out) ラプラシアンを作用された 2 次元スペクトルデータ

      w_InvLapla2D_w = w_LaplaInv_w(w) * Radius**2

    end function w_InvLapla2D_w

    !--------------- 鉛直微分計算 -----------------
 
    function wt_DSig_wt(wt)
      !
      ! 入力スペクトルデータに鉛直微分 d/dsig を作用する.
      !
      ! スペクトルデータの動径微分とは, 対応する格子点データに動径微分を
      ! 作用させたデータのスペクトル変換のことである.
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm,0:tm)             :: wt_DSig_wt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: n

!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         wt_DSig_wt(llMin(n):llMax(n),:) = at_DSig_at(wt(llMin(n):llMax(n),:))
!!$      end do


      wt_DSig_wt = at_DSig_at(wt)
    end function wt_DSig_wt

    function wt_DSigDSig_wt(wt)
      !
      ! 入力スペクトルデータに鉛直微分 d/dsig を作用する.
      !
      ! スペクトルデータの動径微分とは, 対応する格子点データに動径微分を
      ! 作用させたデータのスペクトル変換のことである.
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(lm,0:tm)             :: wt_DSigDSig_wt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: n

!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         wt_DSigDSig_wt(llMin(n):llMax(n),:) = at_DSig_at(at_DSig_at(wt(llMin(n):llMax(n),:)))
!!$      end do

      wt_DSigDSig_wt = at_DSig_at(at_DSig_at(wt))
    end function wt_DSigDSig_wt

    function xyt_DSig_xyt(xyt)
      !
      ! 入力スペクトルデータに鉛直微分 ∂/sig を作用する.
      !
      ! スペクトルデータの動径微分とは, 対応する格子点データに動径微分を
      ! 作用させたデータのスペクトル変換のことである.
      !
      real(8), dimension(0:im,jm,0:tm), intent(in) :: xyt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(0:im-1,jm,0:tm)             :: xyt_DSig_xyt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: n
      real(DP) :: at_Work(im*jm, tm+1)

!!$      at_Work = reshape(xyt, (/ im*jm, tm+1 /))
!!$
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         at_Work(lxyMin(n):lxyMax(n),:) = at_DSig_at(at_Work(lxyMin(n):lxyMax(n),:))
!!$      end do
!!$
!!$      xyt_DSig_xyt = reshape(at_Work, shape(xyt_DSig_xyt))

      xyt_DSig_xyt = reshape( &
           &         at_DSig_at(reshape(xyt, (/ im*jm, tm+1 /)) ), &
           &         shape(xyt_DSig_xyt) )

    end function xyt_DSig_xyt

    function xyt_DSigDSig_xyt(xyt)
      !
      ! 入力スペクトルデータに鉛直微分 ∂/sig を作用する.
      !
      ! スペクトルデータの動径微分とは, 対応する格子点データに動径微分を
      ! 作用させたデータのスペクトル変換のことである.
      !
      real(8), dimension(0:im,jm,0:tm), intent(in) :: xyt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(0:im-1,jm,0:tm)             :: xyt_DSigDSig_xyt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ


      integer :: n
      real(DP) :: at_Work(im*jm, tm+1)

!!$      at_Work = reshape(xyt, (/ im*jm, tm+1 /))
!!$
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         at_Work(lxyMin(n):lxyMax(n),:) = at_DSig_at(at_DSig_at(at_Work(lxyMin(n):lxyMax(n),:)))
!!$      end do
!!$
!!$      xyt_DSigDSig_xyt = reshape(at_Work, shape(xyt_DSigDSig_xyt))


      xyt_DSigDSig_xyt = reshape( &
           &         at_DSig_at(at_DSig_at( reshape(xyt, (/ im*jm, tm+1 /)) )), &
           &         shape(xyt_DSigDSig_xyt) )

    end function xyt_DSigDSig_xyt

    !--------------- 鉛直積分計算 -----------------
    
    function xy_IntSig_BtmToTop_xyz(xyz) result(xy_Int)

      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      real(8), dimension(0:im-1,1:jm) :: xy_Int

      integer :: i,j,k 

      !$omp parallel do private(i)
      do j=1, jm
         do i=0, im-1
            xy_Int(i,j) = sum(xyz(i,j,:)*g_Sig_WEIGHT(:))
         end do
      end do
      !$omp end parallel do

    end function xy_IntSig_BtmToTop_xyz

    function w_IntSig_BtmToTop_wz(wz) result(w_Int)

      real(8), dimension(lm,0:km), intent(in)   :: wz
      real(8), dimension(lm) :: w_Int

      integer :: l

      !$omp parallel do
      do l=1, lm
         w_Int(l) = sum(wz(l,:)*g_Sig_WEIGHT(:))
      end do
      !$omp end parallel do

    end function w_IntSig_BtmToTop_wz

    function xyz_IntSig_SigToTop_xyz(xyz) result(xyz_Int)

      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      real(8), dimension(0:im-1,1:jm,0:km) :: xyz_Int
      
      integer :: i, tr, ljm

      !$omp parallel private(tr,i, ljm)
      tr = omp_get_thread_num() + 1
      ljm = ljMax(tr)-ljMin(tr)+1

      do i=0, im-1
         xyz_Int(i,ljMin(tr):ljMax(tr),:) = matmul(xyz(i,ljMin(tr):ljMax(tr),:),tr_vIntCoefMat)
      end do
      !$omp end parallel

!!$      xyz_Int = 0d0
!!$      !$omp parallel do private(k2)
!!$      do k=0,km
!!$         do k2=0,km
!!$            xyz_Int(:,:,k) = xyz_Int(:,:,k) + tr_vIntCoefMat(k,k2)*xyz(:,:,k2)
!!$         end do
!!$      end do

    end function xyz_IntSig_SigToTop_xyz


    !> @brief 
    !!
    !!
    subroutine apply_ZBoundaryCond( wt_field, SurfBCType, BtmBCType, w_SurfBCWork, w_BtmBCWork )
      
      ! 宣言文; Declaration statement
      !
      real(DP), intent(inout) :: wt_field(lm, 0:tm)
      character, intent(in) :: SurfBCType
      character, intent(in) :: BtmBCType
      real(DP), intent(in), optional :: w_SurfBCWork(lm)
      real(DP), intent(in), optional :: w_BtmBCWork(lm)

      
      ! 局所変数
      ! Local variables
      !
      real(DP) :: bcWork(size(wt_field,1), 2)
      
      ! 実行文; Executable statement
      !

      if( present(w_SurfBCWork) ) then
         bcWork(:,1) = w_SurfBCWork
      else
         bcWork(:,1) = 0d0
      end if

      if( present(w_BtmBCWork) ) then
         bcWork(:,2) = w_BtmBCWork
      else
         bcWork(:,2) = 0d0
      end if
      
      select case(SurfBCType//BtmBCType)
      case('DD')
         call at_BoundariesGrid_DD(wt_field, bcWork)
      case('DN')
         call at_BoundariesGrid_DN(wt_field, bcWork)
      case('ND')
         call at_BoundariesGrid_ND(wt_field, bcWork)
      case('NN')
         call at_BoundariesGrid_NN(wt_field, bcWork)
      case default
         call MessageNotify('E', module_name, &
              & "The specified boundary condition for bottom boundary is invalid.")
      end select

    end subroutine apply_ZBoundaryCond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine construct_tr_vIntCoefMat()

      real(DP) :: tt_data(0:tm,0:tm)
      real(DP) :: TMat(0:tm,0:km)
      real(DP) :: TIntMat(0:km,0:tm)
      real(DP) :: Sigk, theta, tl
      integer :: k, t, k2

      allocate(tr_vIntCoefMat(0:km,0:km))

      tt_data = 0d0
      do t=0, tm
         tt_data(t,t) = 1d0
      end do
      TMat = transpose( at_az(tt_data) )


      do k=0, km
         theta = PI*k/dble(km)
         Sigk = cos(theta)
         TIntMat(k,0) = 1d0 - Sigk
         TIntMat(k,1) = 0.5d0*(1d0 - Sigk**2)

         do t=2, tm
            TIntMat(k,t) = 1d0/(1d0-t**2) &
                 & - 0.5d0*( cos((t+1)*theta)/dble(t+1) - cos((t-1)*theta)/dble(t-1) )
         end do
      end do
      TIntMat(:,0) = 0.5d0*TIntMat(:,0)
      TIntMat(:,tm) = 0.5d0*TIntMat(:,tm)

     TIntMat = 0.5d0*(SigMax - SigMin) * TIntMat 

     tr_vIntCoefMat = transpose( matmul(TIntMat, TMat) )

    end subroutine construct_tr_vIntCoefMat

end module SpmlUtil_mod

