!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Boundary_vars_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/SIce
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DSIce_Admin_TInteg_mod, only: &
       & TIMELV_ID_A, TIMELV_ID_N, TIMELV_ID_B, &
       & CurrentTime,                           &
       & nLongTimeLevel
  
  use DSIce_IO_History_mod, only: &
       & DSIce_IO_History_RegistVar,      &
       & DSIce_IO_History_HistPut,        &
       & DSIce_IO_History_IsOutputTiming
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_Boundary_vars_Init, DSIce_Boundary_vars_Final
  public :: DSIce_Boundary_vars_HistPut
  
  
  ! 公開変数
  ! Public variable
  !

  !
  real(DP), public, allocatable :: xy_SfcHFlxAI(:,:)  

  !
  real(DP), public, allocatable :: xy_SfcHFlxAI0(:,:)  
  
  !
  real(DP), public, allocatable :: xy_SfcHFlxAO(:,:)  
  
  !
  real(DP), public, allocatable :: xy_PenSDRFlx(:,:)  

  !
  real(DP), public, allocatable :: xy_DSfcHFlxAIDTs(:,:)  

  !
  real(DP), public, allocatable :: xy_BtmHFlxIO(:,:)  

  !
  real(DP), public, allocatable :: xy_SurplusEn(:,:)

  
  ! The i component of wind stress in (i,j) coordinate system
  real(DP), public, allocatable :: xy_WindStressUAI(:,:)
  real(DP), public, allocatable :: xy_WindStressUIO(:,:)
  
  ! The j component of wind stress in (i,j) coordinate system
  real(DP), public, allocatable :: xy_WindStressVAI(:,:)
  real(DP), public, allocatable :: xy_WindStressVIO(:,:)
  
  ! The freshwater flux used in the sea surface height equation as a volume flux
  real(DP), public, allocatable :: xy_FreshWtFlx(:,:)

  ! The freshwater flux used in the salinity equation as a concentration/dilution effect
  real(DP), public, allocatable :: xy_FreshWtFlxS(:,:)

  ! Rain fall [m/s]
  real(DP), public, allocatable :: xy_RainFall(:,:)

  ! Snow fall [m/s]
  real(DP), public, allocatable :: xy_SnowFall(:,:)

  ! Evaporation [m/s]
  real(DP), public, allocatable :: xy_Evap(:,:)
  
  real(DP), public, allocatable :: xy_SfcAlbedoAI(:,:)

  real(DP), public, allocatable :: xy_SeaSfcTemp(:,:)
  
  real(DP), public, allocatable :: xy_SeaSfcSalt(:,:)

  real(DP), public, allocatable :: xy_OcnFrzTemp(:,:)
  
  real(DP), public, allocatable :: xy_OcnMixLyrDepth(:,:)

  real(DP), public, allocatable :: xy_SIceSfcTemp0(:,:)
  real(DP), public, allocatable :: xy_SDwRFlx(:,:)
  real(DP), public, allocatable :: xy_LDwRFlx(:,:)
  real(DP), public, allocatable :: xy_LatHFlx(:,:)
  real(DP), public, allocatable :: xy_SenHFlx(:,:)
  real(DP), public, allocatable :: xy_DLatSenHFlxDTs(:,:)

  real(DP), public, allocatable :: xy_SeaSfcHFlx_ns(:,:)  
  real(DP), public, allocatable :: xy_SeaSfcHFlx_sr(:,:)  
  real(DP), public, allocatable :: xy_SeaSfcU(:,:)
  real(DP), public, allocatable :: xy_SeaSfcV(:,:)
  
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DSIce_Boundary_vars_mod' !< Module Name

  logical :: isInitialzed = .false.

contains

  !>
  !!
  !!
  Subroutine DSIce_Boundary_vars_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    call MessageNotify( 'M', module_name, "Allocate array memory  ..")

    allocate( xy_SfcHFlxAI(IA,JA), xy_DSfcHFlxAIDTs(IA,JA) )
    allocate( xy_SfcHFlxAI0(IA,JA) )
    allocate( xy_SfcHFlxAO(IA,JA) )
    allocate( xy_PenSDRFlx(IA,JA) )
    allocate( xy_BtmHFlxIO(IA,JA) )
    allocate( xy_SurplusEn(IA,JA) )

    allocate( xy_WindStressUAI(IA,JA), xy_WindStressVAI(IA,JA) )
    allocate( xy_SnowFall(IA,JA), xy_RainFall(IA,JA), xy_Evap(IA,JA) )
    allocate( xy_SfcAlbedoAI(IA,JA) )
    allocate( xy_SDwRFlx(IA,JA), xy_LDwRFlx(IA,JA) )
    allocate( xy_LatHFlx(IA,JA), xy_SenHFlx(IA,JA) )
    allocate( xy_DLatSenHFlxDTs(IA,JA) )

    allocate( xy_WindStressUIO(IA,JA), xy_WindStressVIO(IA,JA) )    
    allocate( xy_SeaSfcHFlx_ns(IA,JA), xy_SeaSfcHFlx_sr(IA,JA) )
    allocate( xy_FreshWtFlx(IA,JA), xy_FreshWtFlxS(IA,JA) )    
    allocate( xy_SeaSfcTemp(IA,JA), xy_SeaSfcSalt(IA,JA) )
    allocate( xy_SIceSfcTemp0(IA,JA) )
    
    allocate( xy_SeaSfcU(IA,JA), xy_SeaSfcV(IA,JA) )
    allocate( xy_OcnFrzTemp(IA,JA) )
    allocate( xy_OcnMixLyrDepth(IA,JA) )
    
    !-------------------------------------------------------

    call DSIce_IO_History_RegistVar( 'WindStressUAI', "IJT", "wind stress (i axis)", "N/m2"   )
    call DSIce_IO_History_RegistVar( 'WindStressVAI', "IJT", "wind stress (j axis)", "N/m2"   )
    call DSIce_IO_History_RegistVar( 'SfcHFlxAI', "IJT", "surface heat flux", "W/m2" )
    call DSIce_IO_History_RegistVar( 'DLatSenHFlxAIDTs', "IJT", "heat flux dependency for surface temperature", "W/(m2.K)" )
    call DSIce_IO_History_RegistVar( 'PenSDRFlx', "IJT", "solar part of surface heat flux", "W/m2" )
    call DSIce_IO_History_RegistVar( 'BtmHFlxIO', "IJT", "bottom heat flux between sea-ice and ocean interfaces", "W/m2" )
    call DSIce_IO_History_RegistVar( 'FreshWtFlx', "IJT", "fresh water flux", "m/s" )
    call DSIce_IO_History_RegistVar( 'FreshWtFlxS', "IJT", "fresh water flux (used in salinity equation", "m/s" )
    call DSIce_IO_History_RegistVar( 'RainFall', "IJT", "rain fall", "kg.m-2.s-1" )
    call DSIce_IO_History_RegistVar( 'SnowFall', "IJT", "snow fall", "kg.m-2.s-1" )
    call DSIce_IO_History_RegistVar( 'Evap', "IJT", "evarporation", "kg.m-2.s-1" )
    call DSIce_IO_History_RegistVar( 'SfcAlbedoAI', 'IJT', 'surface albedo', '1' )
    
    
    !-------------------------------------------------------
    
    isInitialzed = .true.
    
  end subroutine DSIce_Boundary_vars_Init

  !>
  !!
  !!
  subroutine DSIce_Boundary_vars_Final()

    ! 実行文; Executable statements
    !

    if( isInitialzed ) then

       call MessageNotify( 'M', module_name, "Deallocate array memory  ..")

       deallocate( xy_SfcHFlxAI, xy_DSfcHFlxAIDTs, xy_SfcHFlxAO )
       deallocate( xy_SfcHFlxAI0 )
       deallocate( xy_PenSDRFlx )
       deallocate( xy_BtmHFlxIO )

       deallocate( xy_WindStressUAI, xy_WindStressVAI )
       deallocate( xy_RainFall, xy_SnowFall, xy_Evap )
       deallocate( xy_SfcAlbedoAI )
       deallocate( xy_SDwRFlx, xy_LDwRFlx )
       deallocate( xy_LatHFlx, xy_SenHFlx )
       deallocate( xy_DLatSenHFlxDTs )
       deallocate( xy_SIceSfcTemp0 )
       
       deallocate( xy_WindStressUIO, xy_WindStressVIO )
       deallocate( xy_FreshWtFlx, xy_FreshWtFlxS )       
       deallocate( xy_SeaSfcU, xy_SeaSfcV )
       deallocate( xy_SeaSfcTemp, xy_SeaSfcSalt )
       deallocate( xy_OcnFrzTemp )
       deallocate( xy_OcnMixLyrDepth )
    end if

  end subroutine DSIce_Boundary_vars_Final

  !-----------------------------------------

  subroutine DSIce_Boundary_vars_HistPut()

    ! 実行文; Executable statement
    !

!!$    if( .not. DSIce_IO_History_isOutputTiming(CurrentTime) ) return

    call DSIce_IO_History_HistPut( "WindStressUAI", xy_WindStressUAI(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "WindStressVAI", xy_WindStressVAI(IS:IE,JS:JE) )
    
    call DSIce_IO_History_HistPut( "SfcHFlxAI", xy_SfcHFlxAI(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "DSfcHFlxAIDTs", xy_DSfcHFlxAIDTs(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "PenSDRFlx", xy_PenSDRFlx(IS:IE,JS:JE) )

    call DSIce_IO_History_HistPut( "BtmHFlxIO", xy_BtmHFlxIO(IS:IE,JS:JE) )

    call DSIce_IO_History_HistPut( "FreshWtFlx", xy_FreshWtFlx(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "FreshWtFlxS", xy_FreshWtFlxS(IS:IE,JS:JE) )
    
    call DSIce_IO_History_HistPut( "RainFall", xy_RainFall(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "SnowFall", xy_SnowFall(IS:IE,JS:JE) )
    call DSIce_IO_History_HistPut( "Evap", xy_Evap(IS:IE,JS:JE) )

    call DSIce_IO_History_HistPut( "SfcAlbedoAI", xy_SfcAlbedoAI(IS:IE,JS:JE) )
    
  end subroutine DSIce_Boundary_vars_HistPut

  
end module DSIce_Boundary_vars_mod
