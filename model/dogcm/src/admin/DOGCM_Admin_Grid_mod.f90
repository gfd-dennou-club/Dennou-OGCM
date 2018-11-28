!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_Grid_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & STRING, TOKEN, DP

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use GridIndex_mod, only: &
       & IS, IE, IM, IA, IHALO, IBLOCK,    &
       & JS, JE, JM, JA, JHALO, JBLOCK,    &
       & KS, KE, KM, KA, KHALO, KBLOCK,    &
       & iMax, jMax, jMaxGlobal, kMax,     &
       & lMax, nMax, tMax

  use DOGCM_Admin_Constants_mod, only: &
       & PI
  
  use DOGCM_Admin_GovernEq_mod, only: &
       & SolverType,                   &
       & OCNGOVERNEQ_SOLVER_HSPM_VSPM, &
       & OCNGOVERNEQ_SOLVER_HSPM_VFVM       
  
  use DOGCM_GaussSpmGrid_mod, only: &
       & DOGCM_GaussSpmGrid_Init,                    &
       & DOGCM_GaussSpmGrid_Final,                   &
       & DOGCM_GaussSpmGrid_ConstructAxisInfo,       &
       & DOGCM_GaussSpmGrid_ConstructGrid

  use DOGCM_GaussSpmVFvmGrid_mod, only: &
       & DOGCM_GaussSpmVFvmGrid_Init,               &
       & DOGCM_GaussSpmVFvmGrid_Final,              &
       & DOGCM_GaussSpmVFvmGrid_ConstructAxisInfo,  &
       & DOGCM_GaussSpmVFvmGrid_ConstructGrid
  
  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: DOGCM_Admin_Grid_Init, DOGCM_Admin_Grid_Final

  public :: DOGCM_Admin_Grid_construct
  public :: DOGCM_Admin_Grid_UpdateVCoord
  
  ! Cascade


  ! 公開変数
  ! Public variables
  !


  real(DP), public, allocatable :: xy_Lon(:,:)
  real(DP), public, allocatable :: xy_Lat(:,:)
  real(DP), public, allocatable :: xyz_Lon(:,:,:)
  real(DP), public, allocatable :: xyz_Lat(:,:,:)
  real(DP), public, allocatable :: xyz_Z(:,:,:)
  real(DP), public, allocatable :: xy_Topo(:,:)
  
!  real(DP), public, allocatable :: z_Sig(:)  
!  real(DP), public, allocatable :: z_LyrThickSig(:)

  real(DP), public, allocatable :: x_Lon_Weight(:)
  real(DP), public, allocatable :: y_Lat_Weight(:)
  real(DP), public, allocatable :: z_Sig_Weight(:)

  type, public :: AXIS_INFO
     character(TOKEN) :: name
     character(STRING) :: long_name
     character(TOKEN) :: units
     character(TOKEN) :: weight_units
  end type AXIS_INFO

  type(AXIS_INFO), public, save :: IAXIS_info
  type(AXIS_INFO), public, save :: JAXIS_info
  type(AXIS_INFO), public, save :: KAXIS_info
  type(AXIS_INFO), public, save :: TAXIS_info

  real(DP), public, allocatable :: x_CI(:) ! Position of cell center
  real(DP), public, allocatable :: y_CJ(:)
  real(DP), public, allocatable :: z_CK(:)  
  real(DP), public, allocatable :: x_FI(:) ! Position of face
  real(DP), public, allocatable :: y_FJ(:)
  real(DP), public, allocatable :: z_FK(:)  

  real(DP), public, allocatable :: x_CDI(:) 
  real(DP), public, allocatable :: y_CDJ(:)
  real(DP), public, allocatable :: z_CDK(:)  
  real(DP), public, allocatable :: x_FDI(:) ! 
  real(DP), public, allocatable :: y_FDJ(:)
  real(DP), public, allocatable :: z_FDK(:)  
  real(DP), public, allocatable :: x_RCDI(:) 
  real(DP), public, allocatable :: y_RCDJ(:)
  real(DP), public, allocatable :: z_RCDK(:)  
  real(DP), public, allocatable :: x_RFDI(:) ! 
  real(DP), public, allocatable :: y_RFDJ(:)
  real(DP), public, allocatable :: z_RFDK(:)  

  real(DP), public, allocatable :: x_IAXIS_Weight(:)
  real(DP), public, allocatable :: y_JAXIS_Weight(:)
  real(DP), public, allocatable :: z_KAXIS_Weight(:)

  real(DP), public, allocatable :: SCALEF_E1(:,:,:)
  real(DP), public, allocatable :: SCALEF_E2(:,:,:)
  real(DP), public, allocatable :: SCALEF_E3(:,:,:)

  ! Cascade

  public :: IA, IS, IE, IM, IHALO, IBLOCK
  public :: JA, JS, JE, JM, JHALO, JBLOCK
  public :: KA, KS, KE, KM, KHALO, KBLOCK  
  
  public :: iMax, jMax, jMaxGlobal, kMax
  public :: lMax, nMax, tMax
  
  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'DOGCM_Admin_Grid_mod' !< Module Name

  logical :: isInitialzed = .false.
  
contains
  subroutine DOGCM_Admin_Grid_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !
    use DOGCM_Admin_Constants_mod, only: RPlanet

    use DOGCM_GaussSpmGrid_mod, only: &
         & hsvs_im => iMax, hsvs_jm => jMax, hsvs_km => kMax, &
         & hsvs_nm => nMax, hsvs_lm => lMax, hsvs_tm => tMax, &
         & hsvs_jmGl => jMaxGlobe

    use DOGCM_GaussSpmVFvmGrid_mod, only: &
         & hsvf_im => iMax, hsvf_jm => jMax, hsvf_km => kMax, &
         & hsvf_nm => nMax, hsvf_lm => lMax, hsvf_tm => tMax, &
         & hsvf_jmGl => jMaxGlobe
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! local variable
    !

    ! 実行文; Executable statement
    !

    !---------------------------------------

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_GaussSpmGrid_ConstructAxisInfo( &
            & IAXIS_info%name, IAXIS_info%long_name, IAXIS_info%units, IAXIS_info%weight_units,    & ! (out)
            & IS, IE, IA, IHALO,                                                                   & ! (out)
            & JAXIS_info%name, JAXIS_info%long_name, JAXIS_info%units, JAXIS_info%weight_units,    & ! (out)
            & JS, JE, JA, JHALO,                                                                   & ! (out)
            & KAXIS_info%name, KAXIS_info%long_name, KAXIS_info%units, KAXIS_info%weight_units,    & ! (out)
            & KS, KE, KA, KHALO,                                                                   & ! (out)
            & IM, JM, KM                                                                           & ! (in)
            & )
 
       iMax = hsvs_im; jMax = hsvs_jm; kMax = hsvs_km
       nMax = hsvs_nm; lMax = hsvs_lm; tMax = hsvs_tm
       jMaxGlobal = hsvs_jmGl

    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_GaussSpmVFvmGrid_Init( confignmlFileName )
       call DOGCM_GaussSpmVFvmGrid_ConstructAxisInfo( &
            & IAXIS_info%name, IAXIS_info%long_name, IAXIS_info%units, IAXIS_info%weight_units,    & ! (out)
            & IS, IE, IA, IHALO,                                                                   & ! (out)
            & JAXIS_info%name, JAXIS_info%long_name, JAXIS_info%units, JAXIS_info%weight_units,    & ! (out)
            & JS, JE, JA, JHALO,                                                                   & ! (out)
            & KAXIS_info%name, KAXIS_info%long_name, KAXIS_info%units, KAXIS_info%weight_units,    & ! (out)
            & KS, KE, KA, KHALO,                                                                   & ! (out)
            & IM, JM, KM                                                                           & ! (in)
            & )

       iMax = hsvf_im; jMax = hsvf_jm; kMax = hsvf_km
       nMax = hsvf_nm; lMax = hsvf_lm; tMax = hsvf_tm
       jMaxGlobal = hsvf_jmGl

    case default
       call MessageNotify('E', module_name, "Unexcept Solver is specified. Check!")
    end select

    call MessageNotify('M', module_name, "(IS,JS,KS)=(%d,%d,%d)", (/ IS, JS, KS /) )
    call MessageNotify('M', module_name, "(IE,JE,KE)=(%d,%d,%d)", (/ IE, JE, KE /) )
    call MessageNotify('M', module_name, "(IHALO,JHALO,KHALO)=(%d,%d,%d)", (/ IHALO, JHALO, KHALO /) )

    TAXIS_info%name      = "time"
    TAXIS_info%long_name = "time"
    TAXIS_info%units     = "sec"
    
  end subroutine DOGCM_Admin_Grid_Init

  subroutine DOGCM_Admin_Grid_Final()

    ! 実行文; Executable statement
    !

    if(isInitialzed) then
       deallocate( x_CI, x_CDI, x_RCDI, x_FI, x_FDI, x_RFDI )
       deallocate( y_CJ, y_CDJ, y_RCDJ, y_FJ, y_FDJ, y_RFDJ )
       deallocate( z_CK, z_CDK, z_RCDK, z_FK, z_FDK, z_RFDK )
       deallocate( x_IAXIS_Weight, y_JAXIS_Weight, z_KAXIS_Weight )
       
       deallocate( xy_Lon, xy_Lat )
       deallocate( xyz_Lon, xyz_Lat, xyz_Z )
       deallocate( xy_Topo )
       deallocate( SCALEF_E1, SCALEF_E2 )
    end if
    
  end subroutine DOGCM_Admin_Grid_Final


  !> @brief 
  !!
  !!
  subroutine DOGCM_Admin_Grid_construct()

    ! モジュール引用; Use statements
    !

    
    ! 宣言文; Declaration statement
    !
        
    ! 局所変数
    ! Local variables
    !
    integer :: k

    ! 実行文; Executable statement
    !

    !---------------------------------
    
    
    allocate( x_CI(IA), x_CDI(IA), x_RCDI(IA), x_FI(IA), x_FDI(IA), x_RFDI(IA), x_IAXIS_Weight(IA) )
    allocate( y_CJ(JA), y_CDJ(JA), y_RCDJ(JA), y_FJ(JA), y_FDJ(JA), y_RFDJ(JA), y_JAXIS_Weight(JA) )
    allocate( z_CK(KA), z_CDK(KA), z_RCDK(KA), z_FK(KA), z_FDK(KA), z_RFDK(KA), z_KAXIS_Weight(KA) )

    allocate( xy_Lon(IA,JA), xy_Lat(IA,JA) )
    allocate( xyz_Lon(IA,JA,KA), xyz_Lat(IA,JA,KA), xyz_Z(IA,JA,KA) )
    allocate( xy_Topo(IA,JA) )
    allocate( SCALEF_E1(IA,JA,4), SCALEF_E2(IA,JA,4) )

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)    
       call DOGCM_GaussSpmGrid_ConstructGrid( &
            & x_CI, x_CDI, x_FI, x_FDI, x_IAXIS_Weight,      & ! (out)
            & y_CJ, y_CDJ, y_FJ, y_FDJ, y_JAXIS_Weight,      & ! (out)
            & z_CK, z_CDK, z_FK, z_FDK, z_KAXIS_Weight,      & ! (out)
            & xy_Lon, xy_Lat,                                & ! (out)
            & SCALEF_E1, SCALEF_E2,                          & ! (out)
            & IS, IE, IA, IM, IHALO,                         & ! (in)
            & JS, JE, JA, JM, JHALO,                         & ! (in)
            & KS, KE, KA, KM, KHALO                          & ! (in)
            & )
       
    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_GaussSpmVFvmGrid_ConstructGrid( &
            & x_CI, x_CDI, x_FI, x_FDI, x_IAXIS_Weight,      & ! (out)
            & y_CJ, y_CDJ, y_FJ, y_FDJ, y_JAXIS_Weight,      & ! (out)
            & z_CK, z_CDK, z_FK, z_FDK, z_KAXIS_Weight,      & ! (out)
            & xy_Lon, xy_Lat,                                & ! (out)
            & SCALEF_E1, SCALEF_E2,                          & ! (out)
            & IS, IE, IA, IM, IHALO,                         & ! (in)
            & JS, JE, JA, JM, JHALO,                         & ! (in)
            & KS, KE, KA, KM, KHALO                          & ! (in)
            & ) 

    case default
       call MessageNotify('E', module_name, "Unexcept Solver is specified. Check!")
    end select


    !--------------


    x_RCDI(:) = 1d0/x_CDI(:)
    x_RFDI(:) = 1d0/x_FDI(:)

    y_RCDJ(:) = 1d0/y_CDJ(:)
    y_RFDJ(:) = 1d0/y_FDJ(:)

    z_RCDK(:) = 1d0/z_CDK(:)
    z_RFDK(:) = 1d0/z_FDK(:)
    
    do k = 1, KA
       xyz_Lon(:,:,k) = xy_Lon(:,:)
       xyz_Lat(:,:,k) = xy_Lat(:,:)
    end do

    
    !---------------------------------------------

    isInitialzed = .true.
    
  end subroutine DOGCM_Admin_Grid_construct

  subroutine DOGCM_Admin_Grid_UpdateVCoord( xyz_H, &     ! (out)
       & xy_SSH                               )          ! (in)

    ! モジュール引用; Use statement
    !

    use DOGCM_Admin_BC_mod, only: &
       & KinBC_Surface,         &
       & KinBCTYPE_RigidLid,    &
       & KinBCTYPE_LinFreeSurf
    
    use SpmlUtil_mod, only: &
         & xyz_DSig_xyz
    
    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    
    ! 局所変数
    ! Local variables
    !
    integer :: j
    integer :: k
    real(DP) :: xy_TotDep(IA,JA)
    
    ! 実行文; Executable statement
    !


    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)        

       do k = 1, KA
          xyz_Z(:,:,k) = xy_Topo(:,:) * z_CK(k)
       end do
       xyz_H(IS:IE,JS:JE,KS:KE) = xyz_DSig_xyz( xyz_Z(IS:IE,JS:JE,KS:KE) )

    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)

!!$       !$omp parallel do
!!$       do k=KS, KE
!!$          xyz_Z(:,:,k) = + xy_Topo(:,:) * &
!!$               & (tanh(5d0*(z_CK(k)+0.5d0)) - tanh(0.5d0*5d0)) &
!!$               & /(2d0*tanh(0.5d0*5d0))
!!$          xyz_H(:,:,k) = + xy_Topo(:,:) * &
!!$               & 5d0/(2d0*tanh(0.5d0*5d0)*cosh(5d0*(z_CK(k)+0.5d0))**2)
!!$       end do

       ! sigma-coordinate
       !

       !$omp parallel
       !$omp do
       do j=JS, JE
          xy_TotDep(:,j) = xy_Topo(:,j) + xy_SSH(:,j)
       end do
       !$omp do
       do k=KS, KE
          xyz_Z(:,:,k) = xy_TotDep(:,:) * z_CK(k) ! z = H * s
          xyz_H(:,:,k) = xy_TotDep(:,:)           ! dz/ds = H
       end do
       !$omp do
       do j=JS, JE
          xyz_H(:,j,KS-1) = xyz_H(:,j,KS)
          xyz_H(:,j,KE+1) = xyz_H(:,j,KE)
       end do
       !$omp end parallel

    end select

    
  end subroutine DOGCM_Admin_Grid_UpdateVCoord
  
end module DOGCM_Admin_Grid_mod
