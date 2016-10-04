!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
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
  public :: get_SpmlGridInfo

  public :: w_AlphaOptr_xy, xy_AlphaOptr_w
  public :: wa_AlphaOptr_xya, xya_AlphaOptr_wa

  interface calc_VorDiv2UV
     module procedure calc_VorDiv2UV_w_xy
     module procedure calc_VorDiv2UV_wa_xya
  end interface calc_VorDiv2UV
  public :: calc_VorDiv2UV

  interface calc_UVCosLat2VorDiv
     module procedure calc_UVCosLat2VorDiv_xy_w
     module procedure calc_UVCosLat2VorDiv_xya_wa
  end interface calc_UVCosLat2VorDiv
  public :: calc_UVCosLat2VorDiv
  
  public :: z_DSig_z, xyz_DSig_xyz
  public :: xyz_DSigDSig_xyz
  public :: IntSig_BtmToTop, xy_IntSig_BtmToTop_xyz
  public :: xyz_IntSig_SigToTop_xyz, calc_IntSig_BtmToTop
  
  ! 公開変数
  ! Public variables
  !
  real(DP), public, save, allocatable :: xy_CosLat(:,:)
  
  !* Cascade ( from Spml )
  !

  ! grid informatiion
  public :: x_Lon_Weight, y_Lat_Weight
  public :: g_Sig_WEIGHT
  
  ! basic transformation
  public :: w_xy, xy_w
  public :: g_t, t_g
  
  ! Derivate operator
  public :: xya_GradLon_wa, xya_GradLambda_wa, xya_GradLat_wa, xya_GradMu_wa
  public :: xy_GradLon_w, xy_GradLambda_w, xy_GradLat_w, xy_GradMu_w
  public :: w_DivLambda_xy, w_DivMu_xy
  public :: w_Lapla_w, w_LaplaInv_w
  public :: l_nm, nm_l
  public :: xya_wa, wa_xya

  ! Integral or average operator
  public :: IntLonLat_xy, ya_IntLon_xya
  public :: AvrLonLat_xy, ya_AvrLon_xya
  
  ! Interpolation
  public :: a_Interpolate_wa, Interpolate_w  
  public :: Interpolate_t, a_Interpolate_at
  
  ! Operation for spectral analysis
  public :: nma_EnergyFromStreamfunc_wa, na_EnergyFromStreamfunc_wa
  public :: nma_EnstrophyFromStreamfunc_wa, na_EnstrophyFromStreamfunc_wa

  ! Array managed by spml
  public :: rn
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SpmlUtil_mod' !< Module Name

  real(DP) :: Radius
  integer :: im      !< Number of grid points in longitude
  integer :: jm      !< Number of grid points in latitude
  integer :: km      !< Number of vertical layers
  integer :: nm      !< Maximum truncated wave number in horizontal spectral method
  integer :: tm      !< Maximum truncated wave number in vertical spectral method
  integer :: lm      !< Size of array saving data in wavenumber domain

  logical, save :: initalizedFlag = .false.

  integer :: nThread

  real(DP), save, allocatable :: tr_vDeriv1CoefMat(:,:)
  real(DP), save, allocatable :: tr_vDeriv2CoefMat(:,:)
  real(DP), save, allocatable :: tr_vIntCoefMat(:,:)

  real(DP), parameter :: SigMin = -1d0
  real(DP), parameter :: SigMax =  0d0
  
contains

  !>
  !!
  !!
  !>
  !!
  !!
  subroutine SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np)

    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: iMax, jMax, kMax, nMax, tMax
    real(DP), intent(in) :: RPlanet
    integer, intent(in), optional :: np

    ! 作業変数
    ! Work variables
    !
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
    call MessageNotify("M", module_name, "SpmlUtil_mod have been initialized.")
    initalizedFlag = .true.


  end subroutine SpmlUtil_Init

  !>
  !!
  !!
  subroutine SpmlUtil_Final()

    ! 実行文; Executable statements
    !

    if (initalizedFlag) then
       deallocate(xy_CosLat)
       deallocate(tr_vIntCoefMat)
    end if
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
  subroutine get_SpmlGridInfo( &
       & x_Lon_, y_Lat_, z_Sig_, xy_Lon_, xy_Lat_,     &  ! (out)
       & x_Lon_Weight_, y_Lat_Weight_, z_Sig_Weight_   &  ! (out)
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:im-1), intent(out), optional :: x_Lon_, x_Lon_Weight_
    real(DP), dimension(jm), intent(out), optional :: y_Lat_, y_Lat_Weight_
    real(DP), dimension(0:km), intent(out), optional :: z_Sig_, z_Sig_Weight_
    real(DP), dimension(0:im-1,jm), intent(out), optional :: xy_Lon_, xy_Lat_
    
    ! 実行文; Executable statement
    !

    if(present(x_Lon_))  x_Lon_ = x_Lon
    if(present(y_Lat_))  y_Lat_ = y_Lat
    if(present(z_Sig_))  z_Sig_ = g_Sig
    if(present(x_Lon_Weight_))  x_Lon_Weight_ = x_Lon_Weight
    if(present(y_Lat_Weight_))  y_Lat_Weight_ = y_Lat_Weight
    if(present(z_Sig_Weight_))  z_Sig_Weight_ = g_Sig_WEIGHT
    if(present(xy_Lon_))  xy_Lon_ = xy_Lon
    if(present(xy_Lat_))  xy_Lat_ = xy_Lat
    
  end subroutine get_SpmlGridInfo

  !------------ Horizontal derivative -------------------------

  subroutine calc_UVCosLat2VorDiv_xy_w( xy_UCosLat, xy_VCosLat, w_Vor, w_Div )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xy_UCosLat(0:im-1,jm)
    real(DP), intent(in) :: xy_VCosLat(0:im-1,jm)
    real(DP), intent(out) :: w_Vor(lm)
    real(DP), intent(out) :: w_Div(lm)

    ! 実行文; Executable statement
    !
    
    w_Vor(:) = w_AlphaOptr_xy( xy_VCosLat, -xy_UCosLat )
    w_Div(:) = w_AlphaOptr_xy( xy_UCosLat,  xy_VCosLat )
    
  end subroutine calc_UVCosLat2VorDiv_xy_w

  subroutine calc_VorDiv2UV_w_xy( w_Vor, w_Div, xy_U, xy_V )

    ! 宣言文; Declaration statement
    !

    real(DP), intent(in) :: w_Vor(lm)
    real(DP), intent(in) :: w_Div(lm)
    real(DP), intent(out) :: xy_U(0:im-1,jm)
    real(DP), intent(out) :: xy_V(0:im-1,jm)

    ! 作業変数
    ! Work variables
    !
    real(DP) :: w_Psi(lm)
    real(DP) :: w_Chi(lm)

    ! 実行文; Executable statement
    !
    
    w_Psi(:) = w_LaplaInv_w(w_Vor) * Radius**2
    w_Chi(:) = w_LaplaInv_w(w_Div) * Radius**2
    xy_U(:,:) = xy_CosLat * xy_AlphaOptr_w(w_Chi, -w_Psi)
    xy_V(:,:) = xy_CosLat * xy_AlphaOptr_w(w_Psi, w_Chi)
    
  end subroutine calc_VorDiv2UV_w_xy
  
  subroutine calc_UVCosLat2VorDiv_xya_wa(xya_UCosLat, xya_VCosLat, wa_Vor, wa_Div)

    ! 宣言文; Declaration statement
    !

    real(DP), intent(in) :: xya_UCosLat(:,:,:)
    real(DP), intent(in) :: xya_VCosLat(:,:,:)
    real(DP), intent(out) :: wa_Vor(lm,size(xya_UCosLat,3))
    real(DP), intent(out) :: wa_Div(lm,size(xya_UCosLat,3))

    ! 作業変数
    ! Work variables
    
    integer :: k

    ! 実行文; Executable statement
    !

    !$omp parallel do
    do k = 1, size(xya_UCosLat,3)
       call calc_UVCosLat2VorDiv_xy_w( &
            & xya_UCosLat(:,:,k), xya_VCosLat(:,:,k), & ! (in)
            & wa_Vor(:,k), wa_Div(:,k)                & ! (out)
            & )
    end do
    
  end subroutine calc_UVCosLat2VorDiv_xya_wa
  
  subroutine calc_VorDiv2UV_wa_xya(wa_Vor, wa_Div, xya_U, xya_V)

    ! 宣言文; Declaration statement
    !

    real(DP), intent(in) :: wa_Vor(:,:)
    real(DP), intent(in) :: wa_Div(:,:)
    real(DP), intent(out) :: xya_U(0:im-1,jm,size(wa_Vor,2))
    real(DP), intent(out) :: xya_V(0:im-1,jm,size(wa_Vor,2))


    ! 作業変数
    ! Work variables

    integer :: k

    ! 実行文; Executable statement
    !

    !$omp parallel do
    do k = 1, size(wa_Vor,2)
       call calc_VorDiv2UV_w_xy( wa_Vor(:,k), wa_Div(:,k),  & ! (in)
            & xya_U(:,:,k), xya_V(:,:,k)                    & ! (out)
            & )
    end do
    
  end subroutine calc_VorDiv2UV_wa_xya
  
  function w_AlphaOptr_xy(xy_A, xy_B)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: xy_A(0:im-1,jm)
    real(DP), intent(in) :: xy_B(0:im-1,jm)
    real(DP) :: w_AlphaOptr_xy(lm)

    ! 実行文; Executable statement
    !
    
    w_AlphaOptr_xy(:) = ( w_DivLambda_xy(xy_A) + w_DivMu_xy(xy_B) )/Radius
    
  end function w_AlphaOptr_xy

  function xy_AlphaOptr_w(w_A, w_B)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: w_A(lm)
    real(DP), intent(in) :: w_B(lm)
    real(DP) :: xy_AlphaOptr_w(0:im-1,jm)

    ! 実行文; Executable statement
    !

    xy_AlphaOptr_w(:,:) = (xy_GradLambda_w(w_A) + xy_GradMu_w(w_B))/ (Radius*xy_CosLat*xy_CosLat)
    
  end function xy_AlphaOptr_w

  function xya_AlphaOptr_wa(wa_A, wa_B)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: wa_A(:,:)
    real(DP), intent(in) :: wa_B(:,:)
    real(DP) :: xya_AlphaOptr_wa(0:im-1,jm,size(wa_A,2))

    ! 作業変数
    ! Work variables
    integer :: k
    
    ! 実行文; Executable statement
    !

    !$omp parallel do 
    do k = 1,size(wa_A,2)
       xya_AlphaOptr_wa(:,:,k) = xy_AlphaOptr_w(wa_A(:,k), wa_B(:,k))
    end do
    
  end function xya_AlphaOptr_wa

  function wa_AlphaOptr_xya(xya_A, xya_B)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: xya_A(:,:,:)
    real(DP), intent(in) :: xya_B(:,:,:)
    real(DP) :: wa_AlphaOptr_xya(lm,size(xya_A,3))

    ! 作業変数
    ! Work variables
    integer :: k
    
    ! 実行文; Executable statement
    !

    do k = 1, size(xya_A, 3)
       wa_AlphaOptr_xya(:,k) = w_AlphaOptr_xy(xya_A(:,:,k), xya_B(:,:,k))
    end do
    
  end function wa_AlphaOptr_xya
  
  !------------ Vertical derivative -------------------------

  function z_DSig_z(z)

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: z(0:km)
    real(DP) :: z_DSig_z(0:km) 

    ! 局所変数
    ! Local variables
    !    
    integer :: k

    forAll(k=0:km) &
         & z_DSig_z(k) = sum(tr_vDeriv1CoefMat(:,k)*z(:))

  end function z_DSig_z

  function xyz_DSig_xyz(xyz)

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz(0:im-1,jm,0:km)
    real(DP) :: xyz_DSig_xyz(0:im-1,jm,0:km)

    ! 局所変数
    ! Local variables
    !    
    integer :: k, k2

    ! 実行文; Executable statement
    !
    
    xyz_DSig_xyz(:,:,:) = 0d0

    !$omp parallel do private(k2)
    do k=0, km
       do k2=0, km
          xyz_DSig_xyz(:,:,k) = xyz_DSig_xyz(:,:,k) + tr_vDeriv1CoefMat(k2,k)*xyz(:,:,k2)
       end do
    end do

  end function xyz_DSig_xyz

  function xyz_DSigDSig_xyz(xyz)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(in) :: xyz(0:im-1,jm,0:km)
    real(DP) :: xyz_DSigDSig_xyz(0:im-1,jm,0:km)

    ! 局所変数
    ! Local variables
    !    
    integer :: k, k2

    ! 実行文; Executable statement
    !
    
    xyz_DSigDSig_xyz(:,:,:) = 0d0

    !$omp parallel do private(k2)
    do k=0, km
       do k2=0, km
          xyz_DSigDSig_xyz(:,:,k) = xyz_DSigDSig_xyz(:,:,k) + tr_vDeriv2CoefMat(k2,k)*xyz(:,:,k2)
       end do
    end do

  end function xyz_DSigDSig_xyz
    
  !--------------- Vertical Integration  -----------------
    
  function xy_IntSig_BtmToTop_xyz(xyz) result(xy_Int)

    ! 宣言文; Declaration statement
    !            
    real(DP), intent(in) :: xyz(0:im-1,1:jm,0:km)
    real(DP) :: xy_Int(0:im-1,1:jm)

    ! 局所変数
    ! Local variables
    !    
    integer :: i,j,k 

    ! 実行文; Executable statement
    !
    
    !$omp parallel do private(i)
    do j=1, jm
       do i=0, im-1
          xy_Int(i,j) = sum(xyz(i,j,:)*g_Sig_WEIGHT(:))
       end do
    end do
    !$omp end parallel do
    
  end function xy_IntSig_BtmToTop_xyz

  subroutine calc_IntSig_BtmToTop(xy_Int, xyz)

    ! 宣言文; Declaration statement
    !            
    real(DP), intent(out)   :: xy_Int(0:im-1,jm)
    real(DP), intent(in)   :: xyz(0:im-1,jm,0:km)

    ! 局所変数
    ! Local variables
    !    
    integer :: i,j,k 

    ! 実行文; Executable statement
    !
    
    !$omp parallel do private(i)
    do j=1, jm
       do i=0, im-1
          xy_Int(i,j) = sum(xyz(i,j,:)*g_Sig_WEIGHT(:))
       end do
    end do


  end subroutine calc_IntSig_BtmToTop

  function w_IntSig_BtmToTop_wz(wz) result(w_Int)

    ! 宣言文; Declaration statement
    !            
    real(DP), intent(in)   :: wz(lm,0:km)
    real(DP) :: w_Int(lm)

    ! 局所変数
    ! Local variables
    !    
    integer :: l

    ! 実行文; Executable statement
    !

    !$omp parallel do
    do l=1, lm
       w_Int(l) = sum(wz(l,:)*g_Sig_WEIGHT(:))
    end do
    !$omp end parallel do

  end function w_IntSig_BtmToTop_wz

  function xyz_IntSig_SigToTop_xyz(xyz) result(xyz_Int)

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(in)   :: xyz(0:im-1,1:jm,0:km)
    real(DP) :: xyz_Int(0:im-1,1:jm,0:km)

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

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine construct_VDerivateOptrMat

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: tt_I(0:tm,0:tm)
    real(DP) :: TMat(0:tm,0:km)
    real(DP) :: TIntMat(0:km,0:tm)
    real(DP) :: Sigk, theta, tl
    integer :: k, t, k2

    real(DP), parameter :: PI = acos(-1d0)

    ! 実行文; Executable statement
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

