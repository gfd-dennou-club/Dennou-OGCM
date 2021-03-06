!-------------------------------------------------------------
! Copyright (c) 2013-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module LPhys_RediGMHelper_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod, only: &
       & PI, RPlanet, Omega

  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, kMax, lMax, tMax,        &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,  &
       & xyz_Lat, xyz_Lon,                    &
       & KS, KE, z_CK,                        &
       & JBLOCK, xy_Topo
#include "../../admin/DOGCM_Admin_GaussSpmGridIndexDef.h"

 
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: LPhys_RediGMHelper_Init, LPhys_RediGMHelper_Final

  public :: init_taperingScheme, prepare_SlopeTapering
  public :: TaperingFunc_GKW91, TaperingFunc_DM95
  public :: TaperingFunc_LDD95, LinearSlopeTapering

  public :: prepare_DFM08Info
  public :: TaperingDFM08_GM, TaperingDFM08_IDIFF

  public :: xyz_Dz_xyz
  
  ! 公開変数
  ! Public variable
  !
  
  character(*), public, parameter :: TAPERINGTYPE_GKW91_NAME = 'GKW91'
  integer, public, parameter :: TAPERINGTYPE_GKW91_ID = 1
  character(*), public, parameter :: TAPERINGTYPE_DM95_NAME = 'DM95'
  integer, public, parameter :: TAPERINGTYPE_DM95_ID = 2
  character(*), public, parameter :: PBLTAPERINGTYPE_LDD95_NAME = 'LDD95'
  integer, public, parameter :: TAPERINGTYPE_LDD95_ID = 3
  character(*), public, parameter :: PBLTAPERINGTYPE_LINEAR_NAME = 'linear'
  integer, public, parameter :: TAPERINGTYPE_LINEAR_ID = 4


  integer, public, parameter :: TAPERINGTYPE_DFM08_ID = 5
  logical :: DFM08Flag
  type :: DiabaticLyrInfo_DFM08
     real(DP), pointer :: xy_BLD(:,:)
     real(DP), pointer :: xy_TLT(:,:)
     real(DP), pointer :: xy_Lamb(:,:)
  end type DiabaticLyrInfo_DFM08

  type(DiabaticLyrInfo_DFM08), save, public :: DFM08Info
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'LPhys_RediGMHelper_mod' !< Module Name

  integer :: InteriorTaperingType
  integer :: PBLTaperingType

  real(DP) :: SlopeMaxVal
  real(DP) :: DM95_Sd
  
contains

  !>
  !!
  !!
  subroutine LPhys_RediGMHelper_Init()

    ! 宣言文; Declaration statement
    !
      
    ! 実行文; Executable statements
    !

  end subroutine LPhys_RediGMHelper_Init

  !>
  !!
  !!
  subroutine LPhys_RediGMHelper_Final()

    ! 実行文; Executable statements
    !

    if(DFM08Flag) then
       deallocate(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
       nullify(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
    end if

  end subroutine LPhys_RediGMHelper_Final


  subroutine init_taperingScheme( &
       & interiorTaperingName, PBLTaperingName, SlopeMax, & ! (in)
       & Sd, isUsedDFM08                                  & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: interiorTaperingName
    character(*), intent(in) :: PBLTaperingName
    real(DP), intent(in) :: SlopeMax
    real(DP), intent(in), optional :: Sd
    logical, intent(in) ::isUsedDFM08
    
    ! 実行文; Executable statements
    !

    select case(InteriorTaperingName)
    case(TAPERINGTYPE_GKW91_NAME)
       InteriorTaperingType = TAPERINGTYPE_GKW91_ID
    case(TAPERINGTYPE_DM95_NAME)
       InteriorTaperingType = TAPERINGTYPE_DM95_ID
    case default
       call MessageNotify('E', module_name, &
            & ' InteriorTaperingName ''%c'' is invalid.', c1=trim(InteriorTaperingName) )
    end select

    select case(PBLTaperingName)
    case(PBLTAPERINGTYPE_LDD95_NAME)
       PBLTaperingType = TAPERINGTYPE_LDD95_ID
    case(PBLTAPERINGTYPE_LINEAR_NAME)
       PBLTaperingType = TAPERINGTYPE_LINEAR_ID
    case default
       call MessageNotify('E', module_name, &
            & ' PBLTaperingName ''%c'' is invalid.', c1=trim(PBLTaperingName) )       
    end select

    SlopeMaxVal = SlopeMax
    DM95_Sd = Sd

    if(isUsedDFM08) then
       call MessageNotify('E', module_name, &
            & "Tapering scheme by DFM08 is been implemented now..")
       
       DFM08Flag = .true.
       allocate( &
                & DFM08Info%xy_BLD(0:iMax-1,jMax), DFM08Info%xy_TLT(0:iMax-1,jMax), &
                & DFM08Info%xy_Lamb(0:iMax-1,jMax) )
    end if
  end subroutine init_taperingScheme

  subroutine prepare_SlopeTapering(xyz_T,        & !(out)
       & xyz_SLon, xyz_SLat,                     & !(inout)
       & xyz_DensPot, xyz_Depth,                 & !(in)
       & xy_BLD                                  & !(out)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_T(IA,JA,KA)
    real(DP), intent(inout) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(inout) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensPot(IA,JA,KA)
    real(DP), intent(in) :: xyz_Depth(IA,JA,KA)
    real(DP), intent(out), optional :: xy_BLD(IA,JA)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: xyz_T1(IA,JA,KA)
    real(DP) :: xyz_T2(IA,JA,KA)

    integer :: i
    integer :: j
    integer :: k
    
    ! 実行文; Executable statement
    !

    select case(PBLTaperingType)
    case(TAPERINGTYPE_LINEAR_ID)
       call LinearSlopeTapering(xyz_SLon, xyz_SLat, & ! (inout)
            & xyz_DensPot, xyz_Depth,               & ! (in)
            & xy_BLD                                & ! (out)
            & )
    end select
    
    select case(InteriorTaperingType)
    case(TAPERINGTYPE_GKW91_ID)
       call  TaperingFunc_GKW91( xyz_T1, & ! (out)
            & xyz_SLon, xyz_SLat )        ! (in)

    case(TAPERINGTYPE_DM95_ID)
       call  TaperingFunc_DM95( xyz_T1, & ! (out)
            & xyz_SLon, xyz_SLat )        ! (in)
    end select

    select case(PBLTaperingType)
    case(TAPERINGTYPE_LDD95_ID)
       call TaperingFunc_LDD95( xyz_T2,                        & ! (out)
            & xyz_SLon, xyz_SLat, xyz_Depth, xyz_Lat(:,:,KS),  & ! (in)
            & xyz_T1                                           & ! (inout)
            & )
    end select
    
!!$    if(DFM08Flag) then
!!$       call prepare_DFM08Info( &
!!$            & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth )
!!$    end if

    xyz_T = xyz_T1 * xyz_T2

  end subroutine prepare_SlopeTapering
  
  subroutine TaperingFunc_GKW91( xyz_f,              & ! (out)
       & xyz_SLon, xyz_SLat                          & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_f(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)

    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp workshare
    xyz_f(:,:,:) = min(1d0, SlopeMaxVal**2/(xyz_SLon**2 + xyz_SLat**2))
    !$omp end workshare
    !$omp end parallel
    
  end subroutine TaperingFunc_GKW91
  
  subroutine TaperingFunc_DM95( xyz_f,              & ! (out)
       & xyz_SLon, xyz_SLat                         & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_f(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    
    ! 実行文; Executable statement
    !

    !$omp parallel do private(i,j,jj,k) collapse(2)
    do k=KS, KE
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xyz_f(i,j,k) = 0.5d0*(1d0 + tanh((SlopeMaxVal - sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2))/DM95_Sd))
    end do
    end do
    end do
    
  end subroutine TaperingFunc_DM95
  
  subroutine TaperingFunc_LDD95( xyz_f,                      & ! (out)
       & xyz_SLon, xyz_SLat, xyz_Depth, xy_Lat, xyz_fDM95    & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_f(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLon(IA,JA,KA)
    real(DP), intent(in) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(in) :: xyz_Depth(IA,JA,KA)
    real(DP), intent(in) :: xy_Lat(IA,JA)
    real(DP), intent(inout) :: xyz_fDM95(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: z_r(KA)
    real(DP) :: z_Depth(KA)
    real(DP) :: r
    real(DP) :: xy_BarocEddyDispZMax(IA,JA)
    real(DP) :: tmpBarocEddyDispZ
    real(DP) :: xy_BarocEddDispH(IA,JA)    
    real(DP) :: SlopeABS
    real(DP) :: z_SlopeABS(KA)
    real(DP) :: SlopeABSMax

    real(DP), parameter :: c   = 2d0
    real(DP), parameter :: SlopeSqrtCutOff = 4d-3
    real(DP), parameter :: EPS = 1d-13

    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: kStart(1)
    real(DP) :: xi
    
    ! 実行文; Executable statement
    !

    !$omp parallel private(i,j,k,SlopeABS,tmpBarocEddyDispZ,r, xi)
    !$omp do 
    do j=JS, JE
    do i=IS, IS+_IM_-1
       xy_BarocEddDispH(i,j) = c/(2d0*Omega*abs(sin(xy_Lat(i,j))))
       xy_BarocEddDispH(i,j) = min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
       xy_BarocEddyDispZMax(i,j) = SlopeMaxVal*xy_BarocEddDispH(i,j)
    end do
    end do
    
!!$    !$omp do collapse(2)
!!$    do k=KS, KE
!!$    do j=JS, JE
!!$    do i=IS, IE+_IM_-1
!!$       xi = min( 0d0 - xyz_Depth(i,j,k), xyz_Depth(i,j,k) - (-xy_Topo(i,j)) )
!!$       xyz_f(i,j,k) = 1d0
!!$       if (xi < xy_BarocEddyDispZMax(i,j)) then
!!$          SlopeABS = sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
!!$          tmpBarocEddyDispZ = min(SlopeABS, SlopeSqrtCutOff)*xy_BarocEddDispH(i,j) + EPS
!!$          r =  xi/tmpBarocEddyDispZ
!!$          if(r < 1d0) then
!!$             xyz_f(i,j,k) = 0.5d0*(1d0 + sin(PI*(r - 0.5d0)))
!!$          end if
!!$       end if
!!$    end do
!!$    end do
!!$    end do
    !$omp end parallel

!!$    return

    !$omp parallel do private(i, j, k, kStart, SlopeABSMax, SlopeABS, tmpBarocEddyDispZ, z_r, z_Depth, z_SlopeABS) collapse(2)
    do j=JS, JE
    do i=IS, IS+_IM_-1
          
!!$          z_r(:) = 1d0
!!$          z_Depth(:) = xyz_Depth(i,j,:)
!!$          z_SlopeABS(:) = sqrt(xyz_SLon(i,j,:)**2 + xyz_SLat(i,j,:)**2)          
!!$          !
!!$          SlopeABSMax = 0d0
!!$          kStart = min(KE, minloc(abs(xy_BarocEddyDispZMax(i,j)-(0d0-z_Depth(KS:KE))))+1)
!!$
!!$          do k=kStart(1), KS, -1
!!$             SlopeABS = xyz_fDM95(i,j,k)*sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
!!$             if(SlopeABS > SlopeABSMax) SlopeABSMax = SlopeABS
!!$          end do
!!$          tmpBarocEddyDispZ = SlopeABSMax*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
!!$          do k=KS, KE
!!$             z_r(k) = (0d0 - z_Depth(k))/(tmpBarocEddyDispZ + EPS)
!!$             if(z_r(k) > 1d0) then
!!$                z_r(k) = 1d0; exit
!!$             end if
!!$          end do
!!$          
!!$          SlopeABSMax = 0d0          
!!$          kStart = minloc(abs(xy_BarocEddyDispZMax(i,j)-(z_Depth(KS:KE)-z_Depth(KE))))
!!$          do k=kStart(1), KE, +1
!!$             SlopeABS = xyz_fDM95(i,j,k)*sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
!!$             if(SlopeABS > SlopeABSMax) SlopeABSMax = SlopeABS
!!$          end do
!!$          tmpBarocEddyDispZ = SlopeABSMax*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
!!$          do k=KE, KS, -1
!!$             z_r(k) = (z_Depth(k) - z_Depth(KE))/(tmpBarocEddyDispZ + EPS)
!!$             if(z_r(k) > 1d0) then
!!$                z_r(k) = 1d0; exit
!!$             end if
!!$          end do

       !----
       
       z_r(:) = 1d0
       z_Depth(:) = xyz_Depth(i,j,:)
       z_SlopeABS(:) = sqrt(xyz_SLon(i,j,:)**2 + xyz_SLat(i,j,:)**2)          

       SlopeABSMax = 0d0          
       do k=KS, KE
          if ( z_SlopeABS(k) > SlopeABSMax ) SlopeABSMax = z_SlopeABS(k)
          
          if ( 0d0 - z_Depth(k) > xy_BarocEddyDispZMax(i,j) ) then
             tmpBarocEddyDispZ = min(SlopeABSMax, SlopeMaxVal)*xy_BarocEddDispH(i,j)
             z_r(KS:k) = min( (0d0 - z_Depth(KS:k))/(tmpBarocEddyDispZ + EPS), 1d0 )
             exit
          end if
       end do
!!$       if(j == JS+20) then
!!$          write(*,*) "tmpBaroceddydispz", tmpBarocEddyDispZ, "SlopeABSMax", SlopeABSMax, "EddDispH", xy_BarocEddDispH(i,j)
!!$          write(*,*) "r", z_r(KS:KE)
!!$       end if
       
       SlopeABSMax = 0d0
       do k=KE, KS, -1
          if ( z_SlopeABS(k) > SlopeABSMax ) SlopeABSMax = z_SlopeABS(k)

          if ( z_Depth(k) - (-xy_Topo(i,j)) > xy_BarocEddyDispZMax(i,j) ) then
             tmpBarocEddyDispZ = min(SlopeABSMax, SlopeMaxVal)*xy_BarocEddDispH(i,j)
             z_r(k:KE) = min( (z_Depth(k:KE) - (-xy_Topo(i,j)))/(tmpBarocEddyDispZ + EPS), 1d0 )
             exit
          end if
       end do

       !
       xyz_f(i,j,:) = 0.5d0*(1d0 + sin(PI*(z_r(:) - 0.5d0)))

    end do
    end do

!    write(*,*) "----------------"
!    write(*,*) "f", xyz_f(IS,JS+20,KS:KE)
!    write(*,*) "Dep", xyz_Depth(IS,JS+20,KS:KE)
    
    !where(xyz_f /= 1d0)
    !   xyz_fDM95 = 1d0
    !end where
    
  end subroutine TaperingFunc_LDD95

  Subroutine LinearSlopeTapering( xyz_SLon, xyz_SLat,  & ! (inout)
       & xyz_DensPot, xyz_Depth, xy_BLD                & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_SLon(0:iMax-1,jMax,KS:KE)
    real(DP), intent(inout) :: xyz_SLat(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_DensPot(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_Depth(0:iMax-1,jMax,KS:KE)
    real(DP), intent(out), optional :: xy_BLD(0:iMax-1,jMax)

    ! 局所変数
    ! Local variables
    !        
    real(DP), parameter :: Rho_c = 1d-2
    integer :: k_MixLyrBtm
    integer :: k_Start
    real(DP) :: DensPot_MixLyr
    real(DP) :: xy_h(0:iMax-1,jMax)

    integer :: i
    integer :: j
    integer :: k
    
    ! 実行文; Executable statement
    !    

    !$omp parallel do collapse(2) private( &
    !$omp  k_MixLyrBtm, k_Start, DensPot_MixLyr )
    do i=0, iMax-1
       do j=1, jMax

          DensPot_MixLyr = 1d100
          !k_MixLyrBtm = 1
          do k=KS, KE
             if( -xyz_Depth(i,j,k) > 10d0) then
                DensPot_MixLyr = xyz_DensPot(i,j,k)
                k_Start = k;
                exit
             end if
          end do
          do k=k_Start, KE
             if(DensPot_MixLyr + Rho_c < xyz_Denspot(i,j,k)) then
                k_MixLyrBtm = k; exit
             end if
          end do

          xy_h(i,j) = -xyz_Depth(i,j,k_MixLyrBtm)
          xyz_SLon(i,j,0:k_MixLyrBtm) = -xyz_Depth(i,j,0:k_MixLyrBtm)/xy_h(i,j)*xyz_SLon(i,j,k_MixLyrBtm)
          xyz_SLat(i,j,0:k_MixLyrBtm) = -xyz_Depth(i,j,0:k_MixLyrBtm)/xy_h(i,j)*xyz_SLat(i,j,k_MixLyrBtm)
       end do
    end do

    if(present(xy_BLD)) xy_BLD = xy_h
    
  end subroutine LinearSlopeTapering
  
  subroutine prepare_DFM08Info( &
       & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_H, xyz_Depth, &  ! (in)
       & xyz_G                                              &  ! (out)
       & )

    ! 宣言文; Declaration statement
    !        

    real(DP), intent(in) :: xyz_SLon(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_SLat(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_DensPot(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_Depth(0:iMax-1,jMax,KS:KE)
    real(DP), intent(out), optional :: xyz_G(0:iMax-1,jMax,KS:KE)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: DLD
    real(DP) :: BLD
    real(DP) :: TLT    
    integer :: MLD_k
    integer :: DLD_k

    real(DP), parameter :: c = 2d0
    real(DP) :: xyz_DzDensPot(0:iMax-1,jMax,KS:KE)
    real(DP) :: xyz_DzzDensPot(0:iMax-1,jMax,KS:KE)
    real(DP) :: z_DzDensPot(KS:KE)
    real(DP) :: z_DzzDensPot(KS:KE)
    real(DP) :: lamb
    real(DP) :: z
    real(DP) :: SlopeABS
    real(DP) :: ddzDensPot
    real(DP) :: ddzDensPotMax
    real(DP) :: wt(KS:KE)
    
    integer :: DzMax_k
    integer :: minLocId(1)

    integer :: i
    integer :: j    
    integer :: k
    
    ! 実行文; Executable statement
    !    
    
    xyz_DzDensPot = xyz_Dz_xyz(xyz_DensPot, xyz_H)
    xyz_DzzDensPot = xyz_Dz_xyz(xyz_DzDensPot, xyz_H)

    !$omp parallel do collapse(2) private ( &
    !$omp DLD, BLD, TLT, MLD_k, DLD_k,                                                 &
    !$omp z_DzDensPot, z_DzzDensPot, lamb, z, SlopeABS, ddzDensPOt, ddzDensPotMax, wt, &
    !$omp DzMax_k, minLocID )
    do j=1, jMax
       do i=0, iMax-1

          z_DzDensPot(:) = xyz_DzDensPot(i,j,:)
          z_DzzDensPot(:) = xyz_DzzDensPot(i,j,:)

          ddzDensPotMax = 0d0
          do k=KS,KE
             ddzDensPot = -(xyz_DensPot(i,j,KS) - xyz_DensPot(i,j,k))/(xyz_Depth(i,j,KS) - xyz_Depth(i,j,k))
             if(ddzDensPotMax < ddzDensPot) then
                ddzDensPotMax = ddzDensPot
             else if(k > KS) then
                DzMax_k = k
                exit
             end if
          end do
          
          !
          do k=KS, KE
             if((z_DzDensPot(k) + ddzDensPotMax)*(z_DzDensPot(k-1) + ddzDensPotMax) <= 0d0) then
                MLD_k = k; exit
             end if
             if(k==KE) then
                minLocId(:) = minloc(abs(z_DzzDensPot(:) + ddzDensPotMax))
                MLD_k = minLocId(1)
!!$                write(*,*) "j=", j, ": ddzDensPotMax=", ddzDensPotMax
!!$                write(*,*) z_DzDensPot
             end if
          end do
          wt(:) = abs(z_DzDensPot(MLD_k:MLD_k-1:-1) + ddzDensPotMax)/sum(abs(z_DzDensPot(MLD_k-1:MLD_k) + ddzDensPotMax)) 

          !
          BLD = abs(sum(wt(:)*xyz_Depth(i,j,MLD_k-1:MLD_k)))

          SlopeABS = sqrt(xyz_SLon(i,j,MLD_k)**2 + xyz_SLat(i,j,MLD_k)**2)
          TLT =  min(SlopeABS,SlopeMaxVal)*min(max(15d3,c/(2d0*Omega*abs(sin(xyz_Lat(i,j,0))))), 100d3) + 1d0

          DLD = min(BLD + TLT, abs(xyz_Depth(i,j,kMax)))

          DLD_k = MLD_k
          do k=MLD_k, KE
             if( xyz_Depth(i,j,k) <= -DLD ) then
                DLD_k = k; exit
             end if
             if(k==kMax) DLD_k = k
          end do
          
          wt(:) = abs(xyz_Depth(i,j,DLD_k:DLD_k-1:-1) + DLD)/sum(abs(xyz_Depth(i,j,DLD_k-1:DLD_k) + DLD))

          !
          lamb = abs( - sum(wt(:)*z_DzDensPot(DLD_k-1:DLD_k))/sum(wt(:)*z_DzzDensPot(DLD_k-1:DLD_k)) )

          !
          DFM08Info%xy_BLD(i,j) = BLD; DFM08Info%xy_TLT(i,j) = 0d0!TLT; 
          DFM08Info%xy_Lamb(i,j) = 1d0!lamb

!!$          if(DLD_k == kMax) then
!!$             write(*,*) "DensPot", xyz_DensPot(i,j,:)
!!$             write(*,*) "Dz", z_DzDensPot
!!$             write(*,*) "Max", DzMax_k, ddzDensPotMax
!!$             stop
!!$          end if
       end do
    end do


    if(present(xyz_G)) then
       xyz_G = 1d0
       do k=KS, KE
          where (-DFM08Info%xy_BLD(:,:) < xyz_Depth(:,:,k))
             xyz_G(:,:,k) = -xyz_Depth(:,:,k)/(2d0*DFM08Info%xy_BLD + DFM08Info%xy_TLT) &
                  &        *(2d0 + DFM08Info%xy_TLT/DFM08Info%xy_Lamb)
          end where
       end do
    end if
!!$
!!$    write(*,*) "-- DFM08INFO ---------"
!!$    do j=1, jMax/2+1
!!$    write(*,*) "j=", j, "BLD=", DFM08Info%xy_BLD(0,j), "TLT=", DFM08Info%xy_TLT(0,j), "lamb=", DFM08Info%xy_Lamb(0,j)
!!$    end do

  end subroutine prepare_DFM08Info

  subroutine TaperingDFM08_GM(xyz_PsiLon, xyz_PsiLat, &
       & xy_BLD, xy_TLT, xy_Lamb, xyz_Depth )

    ! 宣言文; Declaration statement
    !            
    real(DP), intent(inout), dimension(0:iMax-1,jMax,KS:KE) :: &
         & xyz_PsiLon, xyz_PsiLat
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_BLD, xy_TLT, xy_Lamb
    real(DP), intent(in), dimension(0:iMax-1,jMax,KS:KE) :: xyz_Depth

    integer :: k, i, j
    real(DP), dimension(0:iMax-1,jMax) :: xy_DLD, xy_G, xy_PsiILon, xy_PsiILat
    real(DP) :: G
    logical :: terminate

    ! 実行文; Executable statement
    !    
    
    xy_DLD = xy_BLD + xy_TLT

    xy_PsiILon = 0d0; xy_PsiILat = 0d0
    do k=KS, KE
       where(xyz_Depth(:,:,k) < -xy_DLD .and. xyz_Depth(:,:,k-1) >= -xy_DLD)
          xy_PsiILon(:,:) = xyz_PsiLon(:,:,k);  xy_PsiILat(:,:) = xyz_PsiLat(:,:,k)
       end where
    end do

    terminate = .false.
    do k=KS, KE
       do j=1, jMax
          do i=0, iMax-1
             if(-xy_BLD(i,j) < xyz_Depth(i,j,k)) then
                G = - xyz_Depth(i,j,k)/(2d0*xy_BLD(i,j) + xy_TLT(i,j))*(2d0 + xy_TLT(i,j)/xy_Lamb(i,j))
                xyz_PsiLon(i,j,k) = G*xy_PsiILon(i,j); xyz_PsiLat(i,j,k) = G*xy_PsiILat(i,j);
!!$                if(j== 5) then
!!$                   write(*,*) "*", j, k, ":", xy_BLD(i,j), xy_DLD(i,j), G
!!$                end if

             else if(-xy_DLD(i,j) < xyz_Depth(i,j,k)) then
!!$                G =   - (xyz_Depth(i,j,k) + xy_BLD(i,j))**2/(xy_DLD(i,j)**2 - xy_BLD(i,j)**2)*(1d0 + xy_DLD(i,j)/xy_Lamb(i,j)) &
!!$                     & - xyz_Depth(i,j,k)/(2d0*xy_BLD(i,j) + xy_TLT(i,j))*(2d0 + xy_TLT(i,j)/xy_Lamb(i,j))
!!$
!!$                xyz_PsiLon(i,j,k) = G*xy_PsiILon(i,j); xyz_PsiLat(i,j,k) = G*xy_PsiILat(i,j)
             else if(-xy_DLD(i,j) <  xyz_Depth(i,j,kMax))then
!!$                write(*,*) "Error", j, xy_DLD(i,j)
                terminate = .true.
             end if
          end do
       end do
       if(terminate) stop
!!$       where(-xy_BLD < xyz_Depth(:,:,k))
!!$          xy_G = -xyz_Depth(:,:,k)/(2d0*xy_BLD + xy_TLT)*(2d0 + xy_TLT/xy_Lamb)
!!$          xyz_PsiLon(:,:,k) = xy_G*xy_PsiILon; xyz_PsiLat(:,:,k) = xy_G*xy_PsiILat;
!!$       end where
!!$       where(-xy_DLD < xyz_Depth(:,:,k) .and. -xy_BLD > xyz_Depth(:,:,k) )
!!$          xy_G = 1d0!(xyz_Depth(:,:,k) + xy_BLD)**2/(xy_DLD**2 - xy_BLD**2)*(1d0 + xy_DLD/xy_Lamb) &
!!$!               & - xyz_Depth(:,:,k)/(2d0*xy_BLD + xy_TLT)*(2d0 + xy_TLT/xy_Lamb)
!!$          xyz_PsiLon(:,:,k) = xy_G*xy_PsiILon; xyz_PsiLat(:,:,k) = xy_G*xy_PsiILat;          
!!$       end where
    end do

  end subroutine TaperingDFM08_GM

  subroutine TaperingDFM08_IDIFF(xyz_C, & ! (out)
       & xy_BLD, xy_TLT, xyz_Depth      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xyz_C(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xy_BLD(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_TLT(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_Depth(0:iMax-1,jMax,KS:KE)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_DLD(0:iMax-1,jMax)
    integer :: k

    ! 実行文; Executable statement
    !    
    
    xy_DLD(:,:) = xy_BLD + xy_TLT

    xyz_C(:,:,:) = 0d0    
    do k=KS, KE
       where(-xy_BLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = 1d0
       elsewhere(-xy_DLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = (xyz_Depth(:,:,k) + xy_DLD)/xy_TLT
       end where
    end do
    
  end subroutine TaperingDFM08_IDIFF

  function xyz_Dz_xyz(xyz, xyz_H, isUSedDF) 

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KS:KE)
    logical, intent(in), optional :: isUSedDF
    real(DP) :: xyz_Dz_xyz(0:iMax-1,jMax,KS:KE)
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: s(KS:KE)
    real(DP) :: t(KS:KE)
    integer :: k

    ! 実行文; Executable statement
    !

!!$    if(present(isUsedDF) .and. isUsedDF) then
    
    t(KS+1:KE-1) = z_CK(KS:KE-2) - z_CK(KS+1:KE-1)
    s(KS+1:KE-1) = z_CK(KS+1:KE-1) - z_CK(KS+2:KE)
    
    !$omp parallel 
    !$omp do
    do k = KS+1, KE-1
!!$    do k=1,kMax-1
       xyz_Dz_xyz(:,:,k) = &
            & (s(k)**2*xyz(:,:,k-1) - (s(k)**2-t(k)**2)*xyz(:,:,k) - t(k)**2*xyz(:,:,k+1)) &
            & /(s(k)*t(k)*(s(k) + t(k)))/xyz_H(:,:,k)
    end do
    
    !$omp workshare
    xyz_Dz_xyz(:,:,KS) = &
         & (xyz(:,:,KS) - xyz(:,:,KS+1))/(t(KS+1)*xyz_H(:,:,KS)) 
    xyz_Dz_xyz(:,:,KE) = &
         & (xyz(:,:,KE-1) - xyz(:,:,KE))/(s(KE-1)*xyz_H(:,:,KE))
    !$omp end workshare

    !$omp end parallel

!!$    xyz_Dz_xyz = xyz_DSig_xyz(xyz)/xyz_H

  end function xyz_Dz_xyz
  
end module LPhys_RediGMHelper_mod

