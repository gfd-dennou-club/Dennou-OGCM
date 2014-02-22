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

  use wa_module
!!$  use wa_module_sjpack

!!$  use wa_module, only: &
!!$       & wz_DivLambda_xyz => wa_DivLambda_xya, &
!!$       & wz_DivMu_xyz => wa_DivMu_xya

  use at_module, only: &
       & at_Initial, &
       & g_Sig => g_X, g_Sig_WEIGHT => g_X_WEIGHT, &
       & at_az => at_ag, & 
       & az_at => ag_at, &
       & t_DSig_t => t_Dx_t, at_DSig_at => at_Dx_at, &
       & at_BoundariesGrid_NN, at_BoundariesGrid_DD, &
       & at_BoundariesGrid_ND, at_BoundariesGrid_DN


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SpmlUtil_Init, SpmlUtil_Final
  public :: isInitialzed
  public :: w_xy, xy_w
  public :: wt_DSig_wt, t_DSig_t
  public :: g_Sig
  public :: xyz_wt, wt_xyz, xyz_wz, wz_xyz, wz_wt, wt_wz
 

  public :: wz_AlphaOptr_xyz, w_AlphaOptr_xy, xyz_AlphaOptr_wz, xy_AlphaOptr_w
  public :: wz_Lapla2D_wz, wz_InvLapla2D_wz
  public :: w_InvLapla2D_w

  public :: xy_IntSig_BtmToTop_xyz, xyz_IntSig_SigToTop_xyz, w_IntSig_BtmToTop_wz

  public :: apply_ZBoundaryCond

  ! Cascade
  public :: xya_GradLon_wa, xya_GradLambda_wa, xya_GradLat_wa, xya_GradMu_wa
  public :: w_DivLambda_xy, w_DivMu_xy
  public :: xy_Lon, xy_Lat
  public :: xya_wa, wa_xya 
  public :: AvrLonLat_xy, ya_AvrLon_xya

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
  integer :: im, jm, km, nm, lm
  real(DP), allocatable :: vIntCoefMat(:,:)
  real(DP), allocatable :: vDiffProcInvMat(:,:)

  logical :: initalizedFlag = .false.

  integer :: nThread
  integer, allocatable :: lwMin(:), lwMax(:)
  integer, allocatable :: lkMin(:), lkMax(:)
  integer, allocatable :: ltMin(:), ltMax(:)
  
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
    nm = nMax; lm = tMax;
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

    !
    call construct_vIntCoefMat()
    
    !
!!$    call assign_IDRange(0, kMax+1, lkMin, lkMax)
!!$    call assign_IDRange(0, tMax+1, ltMin, ltMax)
!!$    call assign_IDRange(1, (nm+1)**2, lwMin, lwMax)
!!$
!!$    write(*,*) 'lwMin', lwMin
!!$    write(*,*) 'lwMax', lwMax

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
        chunkSize = (dimSize+1)/nThread
        lIdMin(1) = dimFirstId; lIdMax(1) = chunkSize
        do tId=2,nThread
           lIdMin(tId) = lIdMax(tId-1) + 1
           lIdMax(tId) = lIdMin(tId) + chunkSize + 1
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

    deallocate(vIntCoefMat)

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
      real(8), dimension((nm+1)*(nm+1),0:lm), intent(in) :: wt
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
      real(8), dimension((nm+1)*(nm+1),0:lm)           :: wt_xyz
      !(out) 2 次元球面調和函数チェビシェフスペクトルデータ

      wt_xyz = wt_wz(wa_xya(xyz))

    end function wt_xyz

    function xyz_wz(wz)
      !
      ! 水平スペクトル・動径格子点データから 3 次元格子点データへ(逆)変換する.
      !
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in) :: wz
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
      real(8), dimension((nm+1)*(nm+1),0:km)             :: wz_xyz
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

      wz_xyz = wa_xya(xyz)

    end function wz_xyz

    function wz_wt(wt)
      !
      ! スペクトルデータから水平スペクトル・動径格子点データへ(正)変換する.
      !
      real(8), dimension((nm+1)*(nm+1),0:lm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ
      real(8), dimension((nm+1)*(nm+1),0:km)             :: wz_wt
      !(out) 2 次元球面調和函数スペクトル・動径格子点データ

!!$      integer :: n
!!$      !$omp parallel do
!!$      do n=1, nThread
!!$         wz_wt(lwMin(n):lwMax(n),:) = az_at(wt(lwMin(n):lwMax(n),:))
!!$      end do
      wz_wt = az_at(wt)
    end function wz_wt

    function wt_wz(wz)
      !
      ! 水平スペクトル・動径格子点データからスペクトルデータへ(正)変換する.
      !
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数スペクトル・動径格子点データ
      real(8), dimension((nm+1)*(nm+1),0:lm)             :: wt_wz
      !(out) 2 次元球面調和函数チェビシェフスペクトルデータ

      wt_wz = at_az(wz)

    end function wt_wz

  !--------------- 微分計算 -----------------

    function wt_DSig_wt(wt)
      !
      ! 入力スペクトルデータに動径微分 ∂/∂r を作用する.
      !
      ! スペクトルデータの動径微分とは, 対応する格子点データに動径微分を
      ! 作用させたデータのスペクトル変換のことである.
      !
      real(8), dimension((nm+1)*(nm+1),0:lm), intent(in) :: wt
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension((nm+1)*(nm+1),0:lm)             :: wt_DSig_wt
      !(in) 動径微分された2 次元球面調和函数チェビシェフスペクトルデータ

      wt_DSig_wt = at_DSig_at(wt)

    end function wt_DSig_wt


    function wt_DivLon_xyz(xyz)
      ! 
      ! 格子点データに発散型経度微分 1/acosφ・∂/∂λ を作用させた
      ! スペクトルデータを返す.
      !
      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension((nm+1)*(nm+1),0:lm)       :: wt_DivLon_xyz
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
      real(8), dimension((nm+1)*(nm+1),0:lm)       :: wt_DivLat_xyz
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
      real(8), dimension((nm+1)*(nm+1),0:km)       :: wz_AlphaOptr_xyz
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
      real(8), dimension((nm+1)*(nm+1))       :: w_AlphaOptr_xy
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ

      w_AlphaOptr_xy = ( w_DivLambda_xy(xy_A) + w_DivMu_xy(xy_B) )/Radius

    end function w_AlphaOptr_xy

    function xyz_AlphaOptr_wz(wz_A, wz_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in)   :: wz_A
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in)   :: wz_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(0:im-1,jm,0:km)       :: xyz_AlphaOptr_wz
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ


      real(DP) :: xyz_CosLat(0:im-1,jm,0:km)
      
      xyz_CosLat = cos(spread(xy_Lat,3,km+1))
      xyz_AlphaOptr_wz = ( xya_GradLambda_wa(wz_A) + xya_GradMu_wa(wz_B) )/ (Radius*xyz_CosLat**2)

    end function xyz_AlphaOptr_wz

    function xy_AlphaOptr_w(w_A, w_B)
      !
      ! 格子データに 1/a ( 1/(1-μ^2)・∂A/∂λ + ∂B/∂μ ) を
      ! 作用させたスペクトルデータを返す.
      !
      real(8), dimension((nm+1)*(nm+1)), intent(in)   :: w_A
      real(8), dimension((nm+1)*(nm+1)), intent(in)   :: w_B

      !(in) 3 次元経度緯度動径格子点データ
      real(8), dimension(0:im-1,jm)       :: xy_AlphaOptr_w
      !(out) 発散型緯度微分を作用された 2 次元スペクトルデータ


      xy_AlphaOptr_w = ( xy_GradLambda_w(w_A) + xy_GradMu_w(w_B) )/ (Radius*cos(xy_Lat)**2)

    end function xy_AlphaOptr_w


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
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension((nm+1)*(nm+1),0:km)             :: wz_Lapla2D_wz
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
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in) :: wz
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension((nm+1)*(nm+1),0:km)             :: wz_InvLapla2D_wz
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
      real(8), dimension((nm+1)*(nm+1),0:km), intent(in) :: w
      !(in) 2 次元球面調和函数チェビシェフスペクトルデータ

      real(8), dimension((nm+1)*(nm+1))             :: w_InvLapla2D_w
      !(out) ラプラシアンを作用された 2 次元スペクトルデータ

      w_InvLapla2D_w = w_LaplaInv_w(w) * Radius**2

    end function w_InvLapla2D_w

    !--------------- 鉛直積分計算 -----------------

    function xy_IntSig_BtmToTop_xyz(xyz) result(xy_Int)

      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      real(8), dimension(0:im-1,1:jm) :: xy_Int

      integer :: k 

      xy_Int = 0d0
      !$omp parallel do reduction(+: xy_Int)
      do k=0, km
         xy_Int = xy_Int + xyz(:,:,k)*g_Sig_WEIGHT(k)
      end do

    end function xy_IntSig_BtmToTop_xyz

    function w_IntSig_BtmToTop_wz(wz) result(w_Int)

      real(8), dimension((nm+1)*(nm+1),0:km), intent(in)   :: wz
      real(8), dimension((nm+1)*(nm+1)) :: w_Int

      integer :: k 

      w_Int = 0d0
      !$omp parallel do reduction(+: w_Int) 
      do k=0, km
         w_Int = w_Int + wz(:,k)*g_Sig_WEIGHT(k)
      end do

    end function w_IntSig_BtmToTop_wz

    function xyz_IntSig_SigToTop_xyz(xyz) result(xyz_Int)

      real(8), dimension(0:im-1,1:jm,0:km), intent(in)   :: xyz
      real(8), dimension(0:im-1,1:jm,0:km) :: xyz_Int
      
      integer :: k, k2
      !real(DP) :: xy_sumTmp(0:im-1,1:jm)

      xyz_Int = 0d0
      !$omp parallel do private(k2)
      do k=0,km
         do k2=0,km
            xyz_Int(:,:,k) = xyz_Int(:,:,k) + vIntCoefMat(k,k2)*xyz(:,:,k2)
         end do
      end do

    end function xyz_IntSig_SigToTop_xyz


    !> @brief 
    !!
    !!
    subroutine apply_ZBoundaryCond( wt_field, SurfBCType, BtmBCType, w_SurfBCWork, w_BtmBCWork )
      
      ! 宣言文; Declaration statement
      !
      real(DP), intent(inout) :: wt_field((nm+1)*(nm+1), 0:lm)
      character, intent(in) :: SurfBCType
      character, intent(in) :: BtmBCType
      real(DP), intent(in), optional :: w_SurfBCWork((nm+1)*(nm+1))
      real(DP), intent(in), optional :: w_BtmBCWork((nm+1)*(nm+1))

      
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
    
    subroutine construct_vIntCoefMat()

      real(DP) :: tt_data(0:lm,0:lm)
      real(DP) :: TMat(0:lm,0:km)
      real(DP) :: TIntMat(0:km,0:lm)
      real(DP) :: Sigk, theta, tl
      integer :: k, t, k2

      allocate(vIntCoefMat(0:km,0:km))

      tt_data = 0d0
      do t=0, lm
         tt_data(t,t) = 1d0
      end do
      TMat = transpose( at_az(tt_data) )


      do k=0, km
         theta = PI*k/dble(km)
         Sigk = cos(theta)
         TIntMat(k,0) = 1d0 - Sigk
         TIntMat(k,1) = 0.5d0*(1d0 - Sigk**2)

         do t=2, lm
            TIntMat(k,t) = 1d0/(1d0-t**2) &
                 & - 0.5d0*( cos((t+1)*theta)/dble(t+1) - cos((t-1)*theta)/dble(t-1) )
         end do
      end do
      TIntMat(:,0) = 0.5d0*TIntMat(:,0)
      TIntMat(:,lm) = 0.5d0*TIntMat(:,lm)

     TIntMat = 0.5d0*(SigMax - SigMin) * TIntMat 

     vIntCoefMat = matmul(TIntMat, TMat)

    end subroutine construct_vIntCoefMat

end module SpmlUtil_mod

