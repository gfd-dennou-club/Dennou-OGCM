!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSEddyMixingHelper_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use Constants_mod, only: &
       & PI, RPlanet, Omega

  use SpmlUtil_mod, only: &
       & g_Sig

  use GridSet_mod, only: &
       & iMax, jMax, kMax, lMax, tMax, &
       & xyz_Lat, xyz_Lon

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SGSEddyMixingHelper_Init, SGSEddyMixingHelper_Final

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
  character(*), parameter:: module_name = 'SGSEddyMixingHelper_mod' !< Module Name

  integer :: InteriorTaperingType
  integer :: PBLTaperingType

  real(DP) :: SlopeMaxVal
  real(DP) :: DM95_Sd
  
contains

  !>
  !!
  !!
  subroutine SGSEddyMixingHelper_Init()

    ! 宣言文; Declaration statement
    !
      
    ! 実行文; Executable statements
    !

  end subroutine SGSEddyMixingHelper_Init

  !>
  !!
  !!
  subroutine SGSEddyMixingHelper_Final()

    ! 実行文; Executable statements
    !

    if(DFM08Flag) then
       deallocate(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
       nullify(DFM08Info%xy_BLD, DFM08Info%xy_TLT, DFM08Info%xy_Lamb)
    end if

  end subroutine SGSEddyMixingHelper_Final


  subroutine init_taperingScheme( &
       & interiorTaperingName, PBLTaperingName, SlopeMax, &
       & Sd, isUsedDFM08 )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: interiorTaperingName, PBLTaperingName
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

  subroutine prepare_SlopeTapering(xyz_T, & !(out)
       & xyz_SLon, xyz_SLat,                     & !(inout)
       & xyz_DensPot, xyz_Depth,                 & !(in)
       & xy_BLD                                  & !(out)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1, jMax, 0:kMax), intent(out) :: xyz_T
    real(DP), dimension(0:iMax-1, jMax, 0:kMax), intent(inout) :: &
         & xyz_SLon, xyz_SLat
    real(DP), dimension(0:iMax-1, jMax, 0:kMax), intent(in) :: xyz_DensPot, xyz_Depth
    real(DP), dimension(0:iMax-1, jMax), intent(out), optional :: xy_BLD

    ! 局所変数
    ! Local variables
    !    
    real(DP), dimension(0:iMax-1, jMax, 0:kMax) :: xyz_T1, xyz_T2
    
    ! 実行文; Executable statement
    !

    select case(PBLTaperingType)
    case(TAPERINGTYPE_LINEAR_ID)
       call LinearSlopeTapering(xyz_SLon, xyz_SLat, &
            & xyz_DensPot, xyz_Depth, xy_BLD)
    end select
    
    select case(InteriorTaperingType)
    case(TAPERINGTYPE_GKW91_ID)
       !$omp parallel workshare
       xyz_T1 = TaperingFunc_GKW91(xyz_SLon, xyz_SLat)
       !$omp end parallel workshare       
    case(TAPERINGTYPE_DM95_ID)
       !$omp parallel workshare       
       xyz_T1 = TaperingFunc_DM95(xyz_SLon, xyz_SLat)
       !$omp end parallel workshare

    end select

    select case(PBLTaperingType)
    case(TAPERINGTYPE_LDD95_ID)
       xyz_T2 = TaperingFunc_LDD95(xyz_SLon, xyz_SLat, xyz_Depth, xyz_T1)
    end select
    
!!$    if(DFM08Flag) then
!!$       call prepare_DFM08Info( &
!!$            & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth )
!!$    end if

    xyz_T = xyz_T1*xyz_T2
  end subroutine prepare_SlopeTapering
  
  elemental function TaperingFunc_GKW91(SLon, SLat) result(f)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: SLon, SLat
    real(DP) :: f

    ! 実行文; Executable statement
    !
    
    f = min(1d0, SlopeMaxVal**2/(SLon**2 + SLat**2))
    
  end function TaperingFunc_GKW91

  elemental function TaperingFunc_DM95(SLon, SLat) result(f)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: SLon, SLat
    real(DP) :: f

    ! 実行文; Executable statement
    !
    f = 0.5d0*(1d0 + tanh((SlopeMaxVal - sqrt(SLon**2 + SLat**2))/DM95_Sd))

  end function TaperingFunc_DM95

  function TaperingFunc_LDD95(xyz_SLon, xyz_SLat, xyz_Depth, xyz_fDM95) result(xyz_f)

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_SLon, xyz_SLat, xyz_Depth
    real(DP), intent(inout) :: xyz_fDM95(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_f(0:iMax-1, jMax, 0:kMax)


    ! 局所変数
    ! Local variables
    !    
    real(DP) :: r, xy_BarocEddyDispZ(0:iMax-1,jMax), tmpBarocEddyDispZ, SlopeABS, SlopeABSMax
    real(DP) :: xy_BarocEddDispH(0:iMax-1,jMax)
    real(DP), parameter :: c = 2d0
    real(DP), parameter :: EPS = 1d-13
    real(DP) :: xyz_r(0:iMax-1,jMax,0:kMax)
    integer :: i, j, k, kStart(1)
    
    ! 実行文; Executable statement
    !

    xy_BarocEddDispH = c/(2d0*Omega*abs(sin(xyz_Lat(:,:,0))))
    
!!$    xyz_f = 1d0
!!$    xyz_f(:,:,0) = 0d0; xyz_f(:,:,kMax) = 0d0
!!$    return


!!$    do k=0,kMax
!!$       do j=1,jMax
!!$          do i=0, iMax-1
!!$             SlopeABS = sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
!!$             tmpBarocEddyDispZ = min(max(15d3, xy_BarocEddDispH(i,j)), 100d3)*SlopeABS
!!$             r =  min(0d0 - xyz_Depth(i,j,k), xyz_Depth(i,j,k) - xyz_Depth(i,j,kMax))/tmpBarocEddyDispZ
!!$             if(r < 1d0) then
!!$                xyz_f(i,j,k) = 0.5d0*(1d0 + sin(PI*(r - 0.5d0)))
!!$             else
!!$                xyz_f(i,j,k) = 1d0
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$    return
!!$    
    do j=1,jMax
       do i=0, iMax-1
          xy_BarocEddyDispZ(i,j) = SlopeMaxVal*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
       end do
    end do

    xyz_r(:,:,:) = 1d0
    do j=1,jMax
       do i=0, iMax-1

          !
          SlopeABSMax = 0d0
          kStart = min(kMax, minloc(abs(xy_BarocEddyDispZ(i,j)-(0d0-xyz_Depth(i,j,:))))+1)
          do k=kStart(1), 0, -1
             SlopeABS = xyz_fDM95(i,j,k)*sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
             if(SlopeABS > SlopeABSMax) SlopeABSMax = SlopeABS
          end do
          tmpBarocEddyDispZ = SlopeABSMax*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
          do k=0, kMax
             xyz_r(i,j,k) = (0d0 - xyz_Depth(i,j,k))/(tmpBarocEddyDispZ + EPS)
             if(xyz_r(i,j,k) > 1d0) then
                xyz_r(i,j,k) = 1d0; exit
             end if
          end do
          
          SlopeABSMax = 0d0          
          kStart = minloc(abs(xy_BarocEddyDispZ(i,j)-(xyz_Depth(i,j,:)-xyz_Depth(i,j,kMax))))
          do k=kStart(1), 0, -1
             SlopeABS = xyz_fDM95(i,j,k)*sqrt(xyz_SLon(i,j,k)**2 + xyz_SLat(i,j,k)**2)
             if(SlopeABS > SlopeABSMax) SlopeABSMax = SlopeABS
          end do
          tmpBarocEddyDispZ = SlopeABSMax*min(max(15d3,xy_BarocEddDispH(i,j)), 100d3)
          do k=kMax, 0, -1
             xyz_r(i,j,k) = (xyz_Depth(i,j,k) - xyz_Depth(i,j,kMax))/(tmpBarocEddyDispZ + EPS)
             if(xyz_r(i,j,k) > 1d0) then
                xyz_r(i,j,k) = 1d0; exit
             end if
          end do
       end do
    end do

    !
    xyz_f = 0.5d0*(1d0 + sin(PI*(xyz_r - 0.5d0)))
    !where(xyz_f /= 1d0)
    !   xyz_fDM95 = 1d0
    !end where
    
  end function TaperingFunc_LDD95

  Subroutine LinearSlopeTapering(xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth, &
    & xy_BLD )

    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(inout) :: &
         xyz_SLon, xyz_SLat
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_DensPot, xyz_Depth
    real(DP), dimension(0:iMax-1,jMax), intent(out), optional :: xy_BLD

    real(DP), parameter :: Rho_c = 1d-2

    integer :: i, j, k
    integer :: k_MixLyrBtm, k_Start
    real(DP) :: DensPot_MixLyr
    real(DP) :: xy_h(0:iMax-1,jMax)
    
    do i=0, iMax-1
       do j=1, jMax

          DensPot_MixLyr = 1d100
          !k_MixLyrBtm = 1
          do k=1, kMax
             if( -xyz_Depth(i,j,k) > 10d0) then
                DensPot_MixLyr = xyz_DensPot(i,j,k)
                k_Start = k;
                exit
             end if
          end do
          do k=k_Start, kMax
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
       & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth, &
       & xyz_G )

    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, xyz_SLat, xyz_DensPot, xyz_Depth
    real(DP), intent(out), optional :: xyz_G(0:iMax-1,jMax,0:kMax)
    
    integer :: i, j
    integer :: MLD_k, DLD_k
    real(DP) :: DLD, BLD, TLT

    real(DP), parameter :: c = 2d0
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_DzDensPot, xyz_DzzDensPot
    real(DP), dimension(0:kMax)  :: z_DzDensPot, z_DzzDensPot
    real(DP) :: wt(2)
    real(DP) :: lamb, z, SlopeABS, ddzDensPot, ddzDensPotMax
    integer :: k, DzMax_k, minLocId(1)
    
    xyz_DzDensPot = xyz_Dz_xyz(xyz_DensPot)
    xyz_DzzDensPot = xyz_Dz_xyz(xyz_DzDensPot)

    do j=1, jMax
       do i=0, iMax-1

          z_DzDensPot(:) = xyz_DzDensPot(i,j,:)
          z_DzzDensPot(:) = xyz_DzzDensPot(i,j,:)

          ddzDensPotMax = 0d0
          do k=1, kMax
             ddzDensPot = -(xyz_DensPot(i,j,0) - xyz_DensPot(i,j,k))/(xyz_Depth(i,j,0) - xyz_Depth(i,j,k))
             if(ddzDensPotMax < ddzDensPot) then
                ddzDensPotMax = ddzDensPot
             else if(k > 1) then
                DzMax_k = k
                exit
             end if
          end do
          
          !
          do k=1, kMax
             if((z_DzDensPot(k) + ddzDensPotMax)*(z_DzDensPot(k-1) + ddzDensPotMax) <= 0d0) then
                MLD_k = k; exit
             end if
             if(k==kMax) then
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
          do k=MLD_k, kMax
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
       do k=0, kMax
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

    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PsiLon, xyz_PsiLat
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_BLD, xy_TLT, xy_Lamb
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Depth

    integer :: k, i, j
    real(DP), dimension(0:iMax-1,jMax) :: xy_DLD, xy_G, xy_PsiILon, xy_PsiILat
    real(DP) :: G
    logical :: terminate
    
    xy_DLD = xy_BLD + xy_TLT

    xy_PsiILon = 0d0; xy_PsiILat = 0d0
    do k=1, kMax
       where(xyz_Depth(:,:,k) < -xy_DLD .and. xyz_Depth(:,:,k-1) >= -xy_DLD)
          xy_PsiILon(:,:) = xyz_PsiLon(:,:,k);  xy_PsiILat(:,:) = xyz_PsiLat(:,:,k)
       end where
    end do

    terminate = .false.
    do k=0, kMax
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

  subroutine TaperingDFM08_IDIFF(xyz_C, &
       & xy_BLD, xy_TLT, xyz_Depth )

    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: xyz_C
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: xy_BLD, xy_TLT
    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: xyz_Depth

    integer :: k
    real(DP), dimension(0:iMax-1,jMax) :: xy_DLD, xy_C

    
    xy_DLD(:,:) = xy_BLD + xy_TLT

    xyz_C(:,:,:) = 0d0
    
    do k=0, kMax
       where(-xy_BLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = 1d0
       elsewhere(-xy_DLD < xyz_Depth(:,:,k))
          xyz_C(:,:,k) = (xyz_Depth(:,:,k) + xy_DLD)/xy_TLT
       end where
    end do
    
  end subroutine TaperingDFM08_IDIFF

  function xyz_Dz_xyz(xyz, isUSedDF) 

    use VariableSet_mod, only: &
         & xy_totDepthBasic

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz(0:iMax-1,jMax,0:kMax)
    logical, intent(in), optional :: isUSedDF
    real(DP) :: xyz_Dz_xyz(0:iMax-1,jMax,0:kMax)
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: s(0:kMax), t(0:kMax), xyt(0:iMax-1,jMax,0:tMax)
    integer :: k

    ! 実行文; Executable statement
    !

!!$    if(present(isUsedDF) .and. isUsedDF) then
    
    
    t(1:kMax-1) = g_Sig(0:kMax-2) - g_Sig(1:kMax-1)
    s(1:kMax-1) = g_Sig(1:kMax-1) - g_Sig(2:kMax)

    !$omp parallel do
    do k=1,kMax-1
       xyz_Dz_xyz(:,:,k) = &
            & (s(k)**2*xyz(:,:,k-1) - (s(k)**2-t(k)**2)*xyz(:,:,k) - t(k)**2*xyz(:,:,k+1)) &
            & /(s(k)*t(k)*(s(k) + t(k)))/xy_totDepthBasic
    end do
    xyz_Dz_xyz(:,:,0) = &
         & (xyz(:,:,0) - xyz(:,:,1))/(g_Sig(0) - g_Sig(1))/xy_totDepthBasic
    xyz_Dz_xyz(:,:,kMax) = &
         & (xyz(:,:,kMax-1) - xyz(:,:,kMax))/(g_Sig(kMax-1) - g_Sig(kMax))/xy_totDepthBasic

!!$ else
!!$        xyz_Dz_xyz = xyz_DSig_xyz(xyz)/spread(xy_totDepthBasic,3,kMax+1)
!!$     end if
  end function xyz_Dz_xyz
  
end module SGSEddyMixingHelper_mod

