!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief The module which provides some  operators in calculus using spectral method.
!!
!! This module provides some operators(such as gradient, divergence, curl and etc) acting on field data. 
!! Actually, most subroutins in this module call the subroutines in SPMODEL(https://www.gfd-dennou.org/library/spmodel/index.htm). 
!! 
!! @author Yuta Kawai
!! @since 2013
!!
!!
module SpmlUtil_mod 

  ! モジュール引用; Use statements
  !

  !* gtool
  !

  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify


  !* SPMODEL
  !

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
       & t_g, g_t, &
       & IntSig_BtmToTop => Int_g, &
       & t_DSig_t => t_Dx_t, at_DSig_at => at_Dx_at, &
       & Interpolate_t, a_Interpolate_at, &
       & at_BoundariesGrid_NN, at_BoundariesGrid_DD, &
       & at_BoundariesGrid_ND, at_BoundariesGrid_DN

  !* OpenMP
!$  use omp_lib


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
  public :: xyz_GradLon_wz, xyz_GradLat_wz, xyz_GradLambda_wz, xyz_GradMu_wz

  ! Procedures for the discritization of vertical differential or integrator with chebyshev spectral method
  public :: xyt_DSig_xyt, wt_DSig_wt, t_DSig_t, xyt_DSigDSig_xyt, wt_DSigDSig_wt
  public :: z_DSig_z, xyz_DSig_xyz, xyz_DSigDSig_xyz
  public :: IntSig_BtmToTop, xy_IntSig_BtmToTop_xyz, xyz_IntSig_SigToTop_xyz, w_IntSig_BtmToTop_wz

  ! Procedures for the data conversion between real and spectral space with the spectral methods.
  public :: xyz_wt, wt_xyz, wt_xyt, xyz_wz, wz_xyz, wz_wt, wt_wz, xyt_xyz, xyz_xyt
  public :: w_VorDiv2VectorCosLat, wt_VorDiv2VectorCosLat, wz_VorDiv2VectorCosLat
  public :: w_VectorCosLat2VorDiv_2, wt_VectorCosLat2VorDiv, wz_VectorCosLat2VorDiv

  ! Procedures to statisfy the vertical boundary conditions.
  public :: apply_ZBoundaryCond

  public :: get_HorizontalGrid

  !* Cascade
  !


  ! basic transformation
  public :: w_xy, xy_w
  public :: g_t, t_g
  
  ! Derivate operator
  public :: xya_GradLon_wa, xya_GradLambda_wa, xya_GradLat_wa, xya_GradMu_wa
  public :: xy_GradLon_w, xy_GradLambda_w, xy_GradLat_w, xy_GradMu_w
  public :: w_DivLambda_xy, w_DivMu_xy
  public :: l_nm, nm_l
  public :: xya_wa, wa_xya 
  public :: IntLonLat_xy, ya_IntLon_xya
  public :: AvrLonLat_xy, ya_AvrLon_xya
  public :: a_Interpolate_wa, Interpolate_w

  ! Interpolation
  public :: Interpolate_t, a_Interpolate_at
  
  ! Operation for spectral analysis
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
  real(DP), allocatable :: tr_vDeriv1CoefMat(:,:)
  real(DP), allocatable :: tr_vDeriv2CoefMat(:,:)
  real(DP), allocatable :: vDiffProcInvMat(:,:)

  logical, save :: initalizedFlag = .false.

  integer :: nThread
  integer, save, allocatable :: llMin(:), llMax(:)
  integer, save, allocatable :: ljMin(:), ljMax(:)
  integer, save, allocatable :: lkMin(:), lkMax(:)
  integer, save, allocatable :: ltMin(:), ltMax(:)
  integer, save, allocatable :: lxyMin(:), lxyMax(:)

  real(DP), save, allocatable :: xy_CosLat(:,:)

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
       call wa_Initial(nMax, iMax, jMax, kMax+1, 1)
    else
       call wa_Initial(nMax, iMax, jMax, kMax+1)
    end if

    call MessageNotify("M", module_name, "The number of thread is %d", i=(/ nThread /))
    call at_Initial(kMax, tMax, SigMin, SigMax)
    
    Radius = RPlanet

    ! Allocate the memory of work variables in this module. 
    allocate(xy_CosLat(0:im-1,jm))
    xy_CosLat = cos(xy_Lat)
    call construct_VDerivateOptrMat()
    
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

  !> @brief 
  !!
  !!
  subroutine get_HorizontalGrid(x_Lon_, y_Lat_, xy_Lon_, xy_Lat_)
    
    ! 宣言文; Declaration statement
    !
    real(8), dimension(0:im-1), intent(out), optional :: x_Lon_
    real(8), dimension(jm), intent(out), optional :: y_Lat_
    real(8), dimension(0:im-1,jm), intent(out), optional :: xy_Lon_, xy_Lat_
    
    ! 実行文; Executable statement
    !

    if(present(x_Lon_))  x_Lon_ = x_Lon
    if(present(y_Lat_))  y_Lat_ = y_Lat
    if(present(xy_Lon_))  xy_Lon_ = xy_Lon
    if(present(xy_Lat_))  xy_Lat_ = xy_Lat
    
  end subroutine get_HorizontalGrid


  !--------------- 基本変換 -----------------

    function xyz_wt(wt)

      !
      ! スペクトルデータから 3 次元格子点データへ(逆)変換する.
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ
      real(8), dimension(0:im-1,1:jm,0:km)               :: xyz_wt
      !(out) 3 次元経度緯度動径格子点データ

      integer :: k
      real(DP) :: wz(lm, 0:km)
      
      wz(:,:) = az_at(wt)
      !$omp parallel do
      do k=0, km
         xyz_wt(:,:,k) = xy_w(wz(:,k))
      end do
    end function xyz_wt

    function wt_xyz(xyz)
      !
      ! 3 次元格子点データからスペクトルデータへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in) :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)           :: wt_xyz
      !(out) 2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: k
      real(DP) :: wz(lm, 0:km)
      
      !$omp parallel do
      do k=0, km
         wz(:,k) = w_xy(xyz(:,:,k))
      end do
      wt_xyz = at_az(wz)

    end function wt_xyz

    function xyz_wz(wz)
      !
      ! 水平スペクトル・動径格子点データから 3 次元格子点データへ(逆)変換する.
      !
      real(8), dimension(lm,0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数スペクトル・動径格子点データ
      real(8), dimension(0:im-1,1:jm,0:km)               :: xyz_wz
      !(out) 3 次元経度緯度動径格子点データ

      integer :: k
      
      !$omp parallel do
      do k=0, km
         xyz_wz(:,:,k) = xy_w(wz(:,k))
      end do
    end function xyz_wz

    function wz_xyz(xyz)
      !
      ! 3 次元格子データから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:km)             :: wz_xyz
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      integer :: k
      
      !$omp parallel do
      do k=0, km
         wz_xyz(:,k) = w_xy(xyz(:,:,k))
      end do

    end function wz_xyz

    function wt_xyt(xyt)
      !
      ! 3 次元格子データから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension(0:im-1,1:jm,0:tm), intent(in)   :: xyt
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:tm)             :: wt_xyt
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      integer :: t

      !$omp parallel do
      do t=0, tm
         wt_xyt(:,t) = w_xy(xyt(:,:,t))
      end do
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

      wt_wz(:,:) = at_az(wz)

    end function wt_wz

    function xyt_xyz(xyz)
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz
      real(8), dimension(0:im-1,jm,0:tm) :: xyt_xyz

!      integer :: n
!!$      
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         at_Work(lxyMin(n):lxyMax(n),:) = at_az(at_Work(lxyMin(n):lxyMax(n),:))
!!$      end do
!!$

      xyt_xyz(:,:,:) = reshape(at_az(reshape(xyz, (/im*jm, km+1/))), shape(xyt_xyz))

      
    end function xyt_xyz

    function xyz_xyt(xyt)
      real(8), dimension(0:im-1,jm,0:tm), intent(in) :: xyt
      real(8), dimension(0:im-1,jm,0:km) :: xyz_xyt

!      integer :: n
!      real(DP) :: az_Work(im*jm, km+1)
      
!      az_Work(:,:) =
      
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         az_Work(lxyMin(n):lxyMax(n),:) = az_at(az_Work(lxyMin(n):lxyMax(n),:))
!!$      end do

      xyz_xyt(:,:,:) = reshape(az_at(reshape(xyt, (/im*jm, km+1/))), shape(xyz_xyt))

    end function xyz_xyt

    subroutine w_VorDiv2VectorCosLat( w_Vor, w_Div,  & ! (in)
      & xy_UCosLat, xy_VCosLat                       & ! (out)
      & )

      real(8), dimension(lm) :: w_Vor, w_Div
      real(8), dimension(0:im-1,jm) :: xy_UCosLat, xy_VCosLat
      
      real(8), dimension(lm) :: w_Psi, w_Chi

      w_Psi(:) = w_InvLapla2D_w(w_Vor)
      w_Chi(:) = w_InvLapla2D_w(w_Div)

      xy_UCosLat(:,:) = xy_CosLat**2*xy_AlphaOptr_w(w_Chi, -w_Psi)
      xy_VCosLat(:,:) = xy_CosLat**2*xy_AlphaOptr_w(w_Psi,  w_Chi)

    end subroutine w_VorDiv2VectorCosLat


    !> @brief 
    !!
    !!
    subroutine wt_VorDiv2VectorCosLat(  wt_Vor, wt_Div,  & ! (in)
         & xyz_UCosLat, xyz_VCosLat                      & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(lm,0:tm), intent(in) :: wt_Vor, wt_Div
      real(8), dimension(0:im-1,jm,0:km), intent(out) :: xyz_UCosLat, xyz_VCosLat 
      
      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(lm,0:km) :: wz_Vor, wz_Div
      real(8), dimension(lm) :: w_Psi, w_Chi
      integer :: k

      ! 実行文; Executable statement
      !

      wz_Vor(:,:) = wz_wt(wt_Vor); wz_Div(:,:) = wz_wt(wt_Div)
      call wz_VorDiv2VectorCosLat( wz_Vor, wz_Div, &  ! (in)
           & xyz_UCosLat, xyz_VCosLat              &  ! (out)
           & )

    end subroutine wt_VorDiv2VectorCosLat

    !> @brief 
    !!
    !!
    subroutine wz_VorDiv2VectorCosLat(  wz_Vor, wz_Div,  & ! (in)
         & xyz_UCosLat, xyz_VCosLat                      & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(lm,0:km), intent(in) :: wz_Vor, wz_Div
      real(8), dimension(0:im-1,jm,0:km), intent(out) :: xyz_UCosLat, xyz_VCosLat 
      
      ! 局所変数
      ! Local variables
      !
      integer :: k

      ! 実行文; Executable statement
      !

      !$omp parallel do
      do k=0, km
         call w_VorDiv2VectorCosLat( wz_Vor(:,k), wz_Div(:,k), & ! (in)
              & xyz_UCosLat(:,:,k), xyz_VCosLat(:,:,k)         & ! (out)
              & )
      end do

    end subroutine wz_VorDiv2VectorCosLat

    subroutine w_VectorCosLat2VorDiv_2( xy_UCosLat, xy_VCosLat, &
         & w_Vor, w_Div                                         & 
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(0:im-1,jm), intent(in) :: xy_UCosLat, xy_VCosLat
      real(8), dimension(lm), intent(out) :: w_Vor, w_Div

      ! 実行文; Executable statement
      !
      !

      w_Vor(:) = w_AlphaOptr_xy(xy_VCosLat, -xy_UCosLat)
      w_Div(:) = w_AlphaOptr_xy(xy_UCosLat,  xy_VCosLat)
      
    end subroutine w_VectorCosLat2VorDiv_2

    subroutine wt_VectorCosLat2VorDiv(  xyz_UCosLat, xyz_VCosLat,  & ! (in)
         & wt_Vor, wt_Div                                          & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz_UCosLat, xyz_VCosLat
      real(8), dimension(lm,0:tm), intent(out) :: wt_Vor, wt_Div

      ! 局所変数
      ! Local variables
      !
      real(DP), dimension(lm,0:tm) :: wz_Vor, wz_Div

      ! 実行文; Executable statement
      !

      call wz_VectorCosLat2VorDiv( xyz_UCosLat, xyz_VCosLat, & ! (in)
           & wz_Vor, wz_Div                                  & ! (out)
           & )
      wt_Vor(:,:) = wt_wz(wz_Vor); wt_Div(:,:) = wt_wz(wz_Div)

    end subroutine wt_VectorCosLat2VorDiv

    subroutine wz_VectorCosLat2VorDiv(  xyz_UCosLat, xyz_VCosLat,  & ! (in)
         & wz_Vor, wz_Div                                          & ! (out)
         & )
      
      ! 宣言文; Declaration statement
      !
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz_UCosLat, xyz_VCosLat
      real(8), dimension(lm,0:km), intent(out) :: wz_Vor, wz_Div


      ! 局所変数
      ! Local variables
      !
      integer :: k

      ! 実行文; Executable statement
      !

      !$omp parallel do
      do k=0, km 
         call w_VectorCosLat2VorDiv_2( xyz_UCosLat(:,:,k), xyz_VCosLat(:,:,k), & ! (in)
              & wz_Vor(:,k), wz_Div(:,k)                                                   & ! (out)   
              & )
      end do

    end subroutine wz_VectorCosLat2VorDiv

  !--------------- 水平微分計算 -----------------

    function xyz_GradMu_wz(wz) 
      real(8), dimension(lm,0:km), intent(in) :: wz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_GradMu_wz

      integer :: k

      !$omp parallel do
      do k=0, km
         xyz_GradMu_wz(:,:,k) = xy_GradMu_w(wz(:,k))/Radius
      end do

    end function xyz_GradMu_wz

    function xyz_GradLambda_wz(wz) 
      real(8), dimension(lm,0:km), intent(in) :: wz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_GradLambda_wz

      integer :: k

      !$omp parallel do
      do k=0, km
         xyz_GradLambda_wz(:,:,k) = xy_GradLambda_w(wz(:,k))/Radius
      end do
    end function xyz_GradLambda_wz

    function xyz_GradLat_wz(wz) 
      real(8), dimension(lm,0:km), intent(in) :: wz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_GradLat_wz

      integer :: k

      !$omp parallel do
      do k=0, km
         xyz_GradLat_wz(:,:,k) = xy_GradLat_w(wz(:,k))/Radius
      end do

    end function xyz_GradLat_wz

    function xyz_GradLon_wz(wz) 
      real(8), dimension(lm,0:km), intent(in) :: wz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_GradLon_wz

      integer :: k

      !$omp parallel do
      do k=0, km
         xyz_GradLon_wz(:,:,k) = xy_GradLon_w(wz(:,k))/Radius
      end do
    end function xyz_GradLon_wz

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

      wt_DivLat_xyz(:,:) = wt_wz(wa_divlat_xya(xyz/Radius))
    end function wt_DivLat_xyz

    function wz_AlphaOptr_xyz(xyz_A, xyz_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を    ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz_A
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(lm,0:km)       :: wz_AlphaOptr_xyz
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ

      integer :: k
      
      !$omp parallel do
      do k=0, km
         wz_AlphaOptr_xyz(:,k) = w_AlphaOptr_xy(xyz_A(:,:,k), xyz_B(:,:,k))
      end do

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

      w_AlphaOptr_xy(:) = ( w_DivLambda_xy(xy_A) + w_DivMu_xy(xy_B) )/Radius

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

      integer :: k
      
      !$omp parallel do
      do k=0, km
         xyz_AlphaOptr_wz(:,:,k) = xy_AlphaOptr_w(wz_A(:,k), wz_B(:,k))
      end do
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


      xy_AlphaOptr_w(:,:) = (xy_GradLambda_w(w_A) + xy_GradMu_w(w_B))/ (Radius*xy_CosLat**2)

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

      w_Lapla2D_w(:) = w_Lapla_w(w)/Radius**2

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

      integer :: k

      !$omp parallel do
      do k=0, km
         wz_Lapla2D_wz(:,k) = w_Lapla2D_w(wz(:,k))
      end do

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

      integer :: k

      !$omp parallel do
      do k=0, km
         wz_InvLapla2D_wz(:,k) = w_InvLapla2D_w(wz(:,k))
      end do

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

      w_InvLapla2D_w(:) = w_LaplaInv_w(w) * Radius**2

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
      real(8), dimension(0:im-1,jm,0:tm), intent(in) :: xyt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension(0:im-1,jm,0:tm)             :: xyt_DSig_xyt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ

      integer :: n, i, j
      real(DP) :: at_Work(im*jm, 0:tm)

!!$      at_Work = reshape(xyt, (/ im*jm, tm+1 /))
!!$
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         at_Work(lxyMin(n):lxyMax(n),:) = at_DSig_at(at_Work(lxyMin(n):lxyMax(n),:))
!!$      end do
!!$
!!$      xyt_DSig_xyt = reshape(at_Work, shape(xyt_DSig_xyt))

      xyt_DSig_xyt(:,:,:) = reshape( &
           &         at_DSig_at(reshape(xyt, (/ im*jm, tm+1 /)) ), &
           &         shape(xyt_DSig_xyt) )

!!$      do j=1, jm
!!$         do i=0, im-1
!!$            at_Work(im*(j-1) + i+1,:) = xyt(i,j,:)
!!$         end do
!!$      end do
!!$
!!$      at_Work = at_DSig_at(at_Work)
!!$      do j=1, jm
!!$         do i=0, im-1
!!$            xyt_DSig_xyt(i,j,:) = at_Work(im*(j-1)+i+1,:)
!!$         end do
!!$      end do
      
    end function xyt_DSig_xyt

    function xyt_DSigDSig_xyt(xyt)
      
      real(8), dimension(0:im-1,jm,0:tm), intent(in) :: xyt
      integer :: n
      real(8), dimension(0:im-1,jm,0:tm) :: xyt_DSigDSig_xyt
      
      real(DP) :: at_Work(im*jm, tm+1)

      xyt_DSigDSig_xyt(:,:,:) = reshape( &
           &         at_DSig_at(at_DSig_at( reshape(xyt, (/ im*jm, tm+1 /)) )), &
           &         shape(xyt_DSigDSig_xyt) )

    end function xyt_DSigDSig_xyt

    function z_DSig_z(z)
      real(8), dimension(0:km) :: z
      real(8), dimension(0:km) :: z_DSig_z

      integer :: k

      forAll(k=0:km) &
           & z_DSig_z(k) = sum(tr_vDeriv1CoefMat(:,k)*z(:))
      
    end function z_DSig_z
    
    function xyz_DSig_xyz(xyz)
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_DSig_xyz

      integer :: k, k2

      xyz_DSig_xyz(:,:,:) = 0d0

      !$omp parallel do private(k2)
      do k=0, km
         do k2=0, km
            xyz_DSig_xyz(:,:,k) = xyz_DSig_xyz(:,:,k) + tr_vDeriv1CoefMat(k2,k)*xyz(:,:,k2)
         end do
      end do
      
    end function xyz_DSig_xyz
    
    function xyz_DSigDSig_xyz(xyz)
      real(8), dimension(0:im-1,jm,0:km), intent(in) :: xyz
      real(8), dimension(0:im-1,jm,0:km) :: xyz_DSigDSig_xyz

      integer :: k, k2

      xyz_DSigDSig_xyz(:,:,:) = 0d0

      !$omp parallel do private(k2)
      do k=0, km
         do k2=0, km
            xyz_DSigDSig_xyz(:,:,k) = xyz_DSigDSig_xyz(:,:,k) + tr_vDeriv2CoefMat(k2,k)*xyz(:,:,k2)
         end do
      end do
      
    end function xyz_DSigDSig_xyz
    
    !--------------- 鉛直積分計算 -----------------
    
    function xy_IntSig_BtmToTop_xyz(xyz) result(xy_Int)

      ! 宣言文; Declaration statement
      !            
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

      ! 宣言文; Declaration statement
      !            
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

      ! 宣言文; Declaration statement
      !      
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      real(8), dimension(0:im-1,1:jm,0:km) :: xyz_Int

      ! 局所変数
      ! Local variables
      !
      integer :: k, k2 

      ! 実行文; Executable statement
      !      
      xyz_Int = 0d0
      !$omp parallel do private(k2)
      do k=0,km
         do k2=0,km
            xyz_Int(:,:,k) = xyz_Int(:,:,k) + tr_vIntCoefMat(k2,k)*xyz(:,:,k2)
         end do
      end do

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
    
    subroutine construct_VDerivateOptrMat

      !
      real(DP) :: tt_I(0:tm,0:tm)
      real(DP) :: TMat(0:tm,0:km)
      real(DP) :: TIntMat(0:km,0:tm)
      real(DP) :: Sigk, theta, tl
      integer :: k, t, k2


      !
      
      allocate(tr_vIntCoefMat(0:km,0:km))
      allocate(tr_vDeriv1CoefMat(0:km,0:km))
      allocate(tr_vDeriv2CoefMat(0:km,0:km))

      tt_I = 0d0
      do t=0, tm
         tt_I(t,t) = 1d0
      end do
      TMat = at_az(tt_I)

      !
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

     tr_vIntCoefMat = transpose( matmul(TIntMat, transpose(TMat)) )

     !
     tr_vDeriv1CoefMat = ( az_at(at_DSig_at(TMat)) )
     tr_vDeriv2CoefMat = ( az_at(at_DSig_at(at_DSig_at(TMat))) )
     
   end subroutine construct_VDerivateOptrMat

end module SpmlUtil_mod

