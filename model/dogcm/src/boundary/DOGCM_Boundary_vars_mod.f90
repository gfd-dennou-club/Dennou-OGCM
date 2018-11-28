!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to manage some variables associated with boundary conditions. 
!! 
!! @author Yuta Kawai
!!
module DOGCM_Boundary_vars_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DOGCM_Admin_TInteg_mod, only: &
       & TIMELV_ID_A, TIMELV_ID_N, TIMELV_ID_B, &
       & CurrentTime
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_Boundary_vars_Init, DOGCM_Boundary_vars_Final

  public :: DOGCM_Boundary_vars_HistPut
  
  ! 公開変数
  ! Public variable
  !

  ! The non solar heat flux at surface
  real(DP), public, allocatable :: xy_SfcHFlx_ns(:,:)  
  real(DP), public, allocatable :: xy_SfcHFlx0_ns(:,:)  
  real(DP), public, allocatable :: xy_SfcHFlxIO_ns(:,:)
  
  ! The incoming solar heat flux at surface
  real(DP), public, allocatable :: xy_SfcHFlx_sr(:,:)  
  real(DP), public, allocatable :: xy_SfcHFlx0_sr(:,:)  
  real(DP), public, allocatable :: xy_SfcHFlxIO_sr(:,:)

  !
  real(DP), public, allocatable :: xy_DSfcHFlxDTs(:,:)
  
  ! The i component of wind stress in (i,j) coordinate system
  real(DP), public, allocatable :: xy_WindStressU(:,:)
  
  ! The j component of wind stress in (i,j) coordinate system
  real(DP), public, allocatable :: xy_WindStressV(:,:)
  
  ! The freshwater flux used in the sea surface height equation as a volume flux
  real(DP), public, allocatable :: xy_FreshWtFlx(:,:)
  real(DP), public, allocatable :: xy_FreshWtFlx0(:,:)
  real(DP), public, allocatable :: xy_FreshWtFlxIO(:,:)

  ! The freshwater flux used in the salinity equation as a concentration/dilution effect
  real(DP), public, allocatable :: xy_FreshWtFlxS(:,:)
  real(DP), public, allocatable :: xy_FreshWtFlxS0(:,:)
  real(DP), public, allocatable :: xy_FreshWtFlxSIO(:,:)

  ! Sea surface temperature
  real(DP), public, allocatable :: xy_SeaSfcTemp(:,:)

  ! Surface air temperature
  real(DP), public, allocatable :: xy_SfcAirTemp(:,:)
  
  ! Sea surface salinity
  real(DP), public, allocatable :: xy_SeaSfcSalt0(:,:)

  ! Sea surface temperature
  real(DP), public, allocatable :: xy_SeaSfcTemp0(:,:)

  ! Sea surface salinity
  real(DP), public, allocatable :: xy_SeaSfcSalt(:,:)
  
  real(DP), public, allocatable :: xy_SeaSfcU(:,:)
  real(DP), public, allocatable :: xy_SeaSfcV(:,:)
  
  !
  integer, public, allocatable :: xy_OcnSfcCellMask(:,:)
  integer, public, allocatable :: xyz_OcnCellMask(:,:,:)

  !
  real(DP), public, allocatable :: xy_SfcAlbedoAO(:,:)
  
  integer, public, parameter :: OCNCELLMASK_OCEAN = 0
  integer, public, parameter :: OCNCELLMASK_SICE  = 1
  integer, public, parameter :: OCNCELLMASK_LAND  = 2
  integer, public, parameter :: OCNCELLMASK_UNDEF = 4
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Boundary_vars_mod' !< Module Name

  logical :: isInitialzed = .false.

contains

  !>
  !!
  !!
  Subroutine DOGCM_Boundary_vars_Init(configNmlName)

    use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)


    allocate( xy_SfcHFlx_ns(IA,JA), xy_SfcHFlx_sr(IA,JA) )
    allocate( xy_SfcHFlxIO_ns(IA,JA), xy_SfcHFlxIO_sr(IA,JA) )
    allocate( xy_SfcHFlx0_ns(IA,JA), xy_SfcHFlx0_sr(IA,JA) )
    allocate( xy_DSfcHFlxDTs(IA,JA) )
    
    allocate( xy_SfcAirTemp(IA,JA) )

    allocate( xy_WindStressU(IA,JA), xy_WindStressV(IA,JA) )

    allocate( xy_FreshWtFlx(IA,JA), xy_FreshWtFlxS(IA,JA) )
    allocate( xy_FreshWtFlx0(IA,JA), xy_FreshWtFlxS0(IA,JA) )
    allocate( xy_FreshWtFlxIO(IA,JA), xy_FreshWtFlxSIO(IA,JA) )

    allocate( xy_SeaSfcTemp(IA,JA), xy_SeaSfcSalt(IA,JA) )
    allocate( xy_SeaSfcTemp0(IA,JA), xy_SeaSfcSalt0(IA,JA) )
    allocate( xy_SeaSfcU(IA,JA), xy_SeaSfcV(IA,JA) )

    allocate( xy_SfcAlbedoAO(IA,JA) )
    allocate( xy_OcnSfcCellMask(IA,JA), xyz_OcnCellMask(IA,JA,KA) )

    
    !-------------------------------------------------------

    call DOGCM_IO_History_RegistVar( 'WindStressU', "IJT", "wind stress (i axis)", "N/m2"   )
    call DOGCM_IO_History_RegistVar( 'WindStressV', "IJT", "wind stress (j axis)", "N/m2"   )
    call DOGCM_IO_History_RegistVar( 'SfcHFlx_ns', "IJT", "non-solar part of surface heat flux", "W/m2" )
    call DOGCM_IO_History_RegistVar( 'SfcHFlx_sr', "IJT", "solar part of surface heat flux", "W/m2" )
    call DOGCM_IO_History_RegistVar( 'SfcHFlxIO_ns', "IJT", "non-solar part of surface heat flux", "W/m2" )
    call DOGCM_IO_History_RegistVar( 'SfcHFlxIO_sr', "IJT", "solar part of surface heat flux", "W/m2" )
    call DOGCM_IO_History_RegistVar( 'SfcHFlxO', "IJT", "(total) surface heat flux", "W/m2" )
    call DOGCM_IO_History_RegistVar( 'DSfcHFlxDTs', "IJT", "The dependence  of surface heat flux on SST", "W/(m2.K)" )
    call DOGCM_IO_History_RegistVar( 'FreshWtFlx', "IJT", "fresh water flux", "m/s" )
    call DOGCM_IO_History_RegistVar( 'FreshWtFlxS', "IJT", "fresh water flux (used in salinity equation", "m/s" )
    call DOGCM_IO_History_RegistVar( 'FreshWtFlxIO', "IJT", "fresh water flux", "m/s" )
    call DOGCM_IO_History_RegistVar( 'FreshWtFlxSIO', "IJT", "fresh water flux (used in salinity equation", "m/s" )
    call DOGCM_IO_History_RegistVar( 'SeaSfcTemp', "IJT", "Sea surface temperature", "K" )
    call DOGCM_IO_History_RegistVar( 'SeaSfcSalt', "IJT", "Sea surface salinity", "psu" )
    call DOGCM_IO_History_RegistVar( 'OcnSfcCellMask', "IJT", "map of cell type managed by ocean model", "1" )
    
    
    !-------------------------------------------------------

    xy_OcnSfcCellMask(:,:) = OCNCELLMASK_UNDEF
    xyz_OcnCellMask(:,:,:) = OCNCELLMASK_UNDEF

    xy_OcnSfcCellMask(IS:IE,JS:JE) = OCNCELLMASK_OCEAN
    xyz_OcnCellMask(IS:IE,JS:JE,KS:KE) = OCNCELLMASK_OCEAN
    
    !-------------------------------------------------------

    
    isInitialzed = .true.
    
  end subroutine DOGCM_Boundary_vars_Init

  !---------------------------------------------------------------
  
  !>
  !!
  !!
  subroutine DOGCM_Boundary_vars_Final()

    ! 実行文; Executable statements
    !

    if( isInitialzed ) then
       deallocate( xy_SfcHFlx_ns, xy_SfcHFlx_sr )
       deallocate( xy_SfcHFlxIO_ns, xy_SfcHFlxIO_sr )
       deallocate( xy_SfcHFlx0_ns, xy_SfcHFlx0_sr )
       deallocate( xy_DSfcHFlxDTs )

       deallocate( xy_SfcAirTemp )

       deallocate( xy_WindStressU, xy_WindStressV )

       deallocate( xy_FreshWtFlx, xy_FreshWtFlxS )
       deallocate( xy_FreshWtFlx0, xy_FreshWtFlxS0 )
       deallocate( xy_FreshWtFlxIO, xy_FreshWtFlxSIO )
       
       deallocate( xy_SeaSfcTemp, xy_SeaSfcSalt )
       deallocate( xy_SeaSfcTemp0, xy_SeaSfcSalt0 )
       deallocate( xy_SeaSfcU, xy_SeaSfcV )

       deallocate( xy_SfcAlbedoAO )
       
       deallocate( xy_OcnSfcCellMask, xyz_OcnCellMask )
    end if

  end subroutine DOGCM_Boundary_vars_Final

  !-----------------------------------------

  subroutine DOGCM_Boundary_vars_HistPut()

    use DOGCM_IO_History_mod, only:      &
       & DOGCM_IO_History_HistPut,       &
       & DOGCM_IO_History_IsOutputTiming

    ! 実行文; Executable statement
    !

!!$    if( .not. DOGCM_IO_History_isOutputTiming(CurrentTime) ) return
    
    call DOGCM_IO_History_HistPut( "WindStressU", xy_WindStressU(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "WindStressV", xy_WindStressV(IS:IE,JS:JE) )
    
    call DOGCM_IO_History_HistPut( "SfcHFlx_ns", xy_SfcHFlx_ns(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "SfcHFlx_sr", xy_SfcHFlx_sr(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "SfcHFlxIO_ns", xy_SfcHFlxIO_ns(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "SfcHFlxIO_sr", xy_SfcHFlxIO_sr(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "SfcHFlxO", xy_SfcHFlx_sr(IS:IE,JS:JE) + xy_SfcHFlx_ns(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "DSfcHFlxDTs", xy_DSfcHFlxDTs(IS:IE,JS:JE) )
    
    call DOGCM_IO_History_HistPut( "FreshWtFlx", xy_FreshWtFlx(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "FreshWtFlxS", xy_FreshWtFlxS(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "FreshWtFlxIO", xy_FreshWtFlxIO(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "FreshWtFlxSIO", xy_FreshWtFlxSIO(IS:IE,JS:JE) )
    
    call DOGCM_IO_History_HistPut( "SeaSfcTemp", xy_SeaSfcTemp(IS:IE,JS:JE) )
    call DOGCM_IO_History_HistPut( "SeaSfcSalt", xy_SeaSfcSalt(IS:IE,JS:JE) )

    call DOGCM_IO_History_HistPut( "OcnSfcCellMask", dble(xy_OcnSfcCellMask(IS:IE,JS:JE)) )
    
  end subroutine DOGCM_Boundary_vars_HistPut
  
end module DOGCM_Boundary_vars_mod
