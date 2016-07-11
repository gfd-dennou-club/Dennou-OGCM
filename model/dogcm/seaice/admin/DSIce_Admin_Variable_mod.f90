!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_Admin_Variable_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  !* Dennou-SIce

  use DSIce_Admin_Constants_mod, only: &
       & Mu, SaltSeaIce
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, JA, KA,              &
       & IS, IE, IM,              &
       & JS, JE, JM,              &
       & KS, KE, KM

  use DSIce_Admin_TInteg_mod, only: &
       & TIMELV_ID_A, TIMELV_ID_N, TIMELV_ID_B, &
       & CurrentTime,                           &
       & nLongTimeLevel

  use DSIce_IO_History_mod, only: &
       & DSIce_IO_History_RegistVar,     &
       & DSIce_IO_History_HistPut,       &
       & DSIce_IO_History_IsOutputTiming

  use DSIce_IO_Restart_mod, only: &
       & DSIce_IO_Restart_RegistVar,     &
       & DSIce_IO_Restart_HistPut,       &
       & DSIce_IO_Restart_HistGet,       &
       & DSIce_IO_Restart_IsOutputTiming
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Admin_Variable_Init, DSIce_Admin_Variable_Final

  public :: DSIce_Admin_Variable_SetDefualtValue

  public :: DSIce_Admin_Variable_AdvanceTStep

  public :: DSIce_Admin_Variable_HistPut
  public :: DSIce_Admin_Variable_HistGet
  public :: DSIce_Admin_Variable_RestartPut

  
  ! 公開変数
  ! Public variable
  !
  
  real(DP), allocatable, public :: xya_SIceCon(:,:,:)
  real(DP), allocatable, public :: xyza_SIceTemp(:,:,:,:)
  real(DP), allocatable, public :: xya_SIceSfcTemp(:,:,:)
  real(DP), allocatable, public :: xyza_SIceEn(:,:,:,:)
  real(DP), allocatable, public :: xya_IceThick(:,:,:)
  real(DP), allocatable, public :: xya_SnowThick(:,:,:)

  real(DP), allocatable, public :: xy_Wice(:,:)
  
  character(*), parameter, public :: VARSET_KEY_SICECON = 'SIceCon'
  character(*), parameter, public :: VARSET_KEY_SNOWTHICK = 'SnowThick'  
  character(*), parameter, public :: VARSET_KEY_ICETHICK = 'IceThick'
  character(*), parameter, public :: VARSET_KEY_SICETEMP = 'SIceTemp'
  character(*), parameter, public :: VARSET_KEY_SICESURFTEMP = 'SIceSurfTemp'  

  character(*), parameter, public :: VARSET_KEY_SICECONB = 'SIceConB'
  character(*), parameter, public :: VARSET_KEY_SNOWTHICKB = 'SnowThickB'  
  character(*), parameter, public :: VARSET_KEY_ICETHICKB = 'IceThickB'
  character(*), parameter, public :: VARSET_KEY_SICETEMPB = 'SIceTempB'
  character(*), parameter, public :: VARSET_KEY_SICESURFTEMPB = 'SIceSurfTempB'  
  
  character(*), parameter, public :: VARSET_KEY_SICEEN = 'SIceEn'
  character(*), parameter, public :: VARSET_KEY_SICEENB = 'SIceEnB'

  ! 非公開手続き
  ! Private procedure
  !
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Admin_Variable_mod' !< Module Name
  logical :: isInitialzed

contains

  !>
  !!
  !!
  subroutine DSIce_Admin_Variable_Init()

    ! 局所変数
    ! Local variable
    !
    integer :: TA
    
    ! 実行文; Executable statements
    !

    call MessageNotify('M', module_name, &
         & "Allocate array to store data of variables..")
    
    TA = nLongTimeLevel

    !----------------------------------------------
    
    allocate( xya_SIceCon(IA,JA,TA) )
    allocate( xya_SIceSfcTemp(IA,JA,TA) )
    allocate( xyza_SIceEn(IA,JA,KA,TA) )
    allocate( xyza_SIceTemp(IA,JA,KA,TA) )
    allocate( xya_IceThick(IA,JA,TA) )
    allocate( xya_SnowThick(IA,JA,TA) )
    allocate( xy_Wice(IA,JA) )
    
    !----------------------------------------------

    ! history

    call DSIce_IO_History_RegistVar('SIceCon', 'IJT', 'concentration of sea ice', '1')
    call DSIce_IO_History_RegistVar('SIceSfcTemp', 'IJT', 'sea-ice surface temperature', 'K')
    call DSIce_IO_History_RegistVar('SIceEn', 'IJKT', &
         & 'sea-ice enthalpy at each layer per unit surface area', 'J.m-2')
    call DSIce_IO_History_RegistVar('SIceTemp', 'IJKT', &
         & 'sea-ice temperature at each layer', 'K')

    call DSIce_IO_History_RegistVar('IceThick', 'IJT', 'ice-layer thickness', 'm')
    call DSIce_IO_History_RegistVar('SnowThick', 'IJT', 'snow-layer thickness', 'm')

    ! restart
    
    call DSIce_IO_Restart_RegistVar('SIceCon', 'IJT', 'concentration of sea ice', '1')
    call DSIce_IO_Restart_RegistVar('SIceConB', 'IJT', 'concentration of sea ice', '1')
    call DSIce_IO_Restart_RegistVar('SIceSfcTemp', 'IJT', 'sea-ice surface temperature', 'K')
    call DSIce_IO_Restart_RegistVar('SIceSfcTempB', 'IJT', 'sea-ice surface temperature', 'K')
    call DSIce_IO_Restart_RegistVar('SIceEn', 'IJKT', &
         & 'sea-ice enthalpy at each layer per unit surface area', 'J.m-2')
    call DSIce_IO_Restart_RegistVar('SIceEnB', 'IJKT', &
         & 'sea-ice enthalpy at each layer per unit surface area', 'J.m-2')
    call DSIce_IO_Restart_RegistVar('SIceTemp', 'IJKT', &
         & 'sea-ice temperature at each layer', 'K')
    call DSIce_IO_Restart_RegistVar('SIceTempB', 'IJKT', &
         & 'sea-ice temperature at each layer', 'K')

    call DSIce_IO_Restart_RegistVar('IceThick', 'IJT', 'ice-layer thickness', 'm')
    call DSIce_IO_Restart_RegistVar('IceThickB', 'IJT', 'ice-layer thickness', 'm')
    call DSIce_IO_Restart_RegistVar('SnowThick', 'IJT', 'snow-layer thickness', 'm')
    call DSIce_IO_Restart_RegistVar('SnowThickB', 'IJT', 'snow-layer thickness', 'm')
    
    !-----------------------------------------------
    
    isInitialzed = .true.
    
  end subroutine DSIce_Admin_Variable_Init

  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Variable_SetDefualtValue()
    
    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statement
    !

    xya_SIceCon = 0d0
    xya_SIceSfcTemp = 0d0
    xyza_SIceTemp = 0d0
    xyza_SIceEn = 0d0
    xya_SnowThick = 0d0
    xya_IceThick = 0d0
    xy_Wice = 0d0
    
  end subroutine DSIce_Admin_Variable_SetDefualtValue

  !>
  !!
  !!
  subroutine DSIce_Admin_Variable_Final()

    ! 実行文; Executable statements
    !


    if ( isInitialzed ) then
       deallocate( xya_SIceCon )
       deallocate( xyza_SIceTemp ) 
       deallocate( xyza_SIceEn ) 
       deallocate( xya_SIceSfcTemp )
       deallocate( xya_IceThick )
       deallocate( xya_SnowThick )
       deallocate( xy_Wice )
    end if
    
  end subroutine DSIce_Admin_Variable_Final


  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Variable_AdvanceTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    ! 局所変数
    ! Local variables
    !
    integer :: n

    ! 実行文; Executable statement
    !

    do n = nLongTimeLevel, 2, -1
       !$omp parallel
       !$omp workshare
       xya_SIceCon(:,:,n) = xya_SIceCon(:,:,n-1)
       xya_SIceSfcTemp(:,:,n) = xya_SIceSfcTemp(:,:,n-1)
       xya_IceThick(:,:,n) = xya_IceThick(:,:,n-1)
       xya_SnowThick(:,:,n) = xya_SnowThick(:,:,n-1)       
       !$omp end workshare

       !$omp workshare
       xyza_SIceTemp(:,:,:,n) = xyza_SIceTemp(:,:,:,n-1)
       !$omp end workshare
       !$omp end parallel
    end do
    
  end subroutine DSIce_Admin_Variable_AdvanceTStep

  !----------------------------------------------------------


  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Variable_HistPut()

    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statement
    !
    
    if( .not. DSIce_IO_History_isOutputTiming(CurrentTime) ) return

    call DSIce_IO_History_HistPut( 'SIceCon', xya_SIceCon(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_History_HistPut( 'SIceSfcTemp', xya_SIceSfcTemp(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_History_HistPut( 'SIceEn', xyza_SIceEn(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_History_HistPut( 'SIceTemp', xyza_SIceTemp(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_History_HistPut( 'IceThick', xya_IceThick(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_History_HistPut( 'SnowThick', xya_SnowThick(IS:IE,JS:JE, TIMELV_ID_N) )

  end subroutine DSIce_Admin_Variable_HistPut

  !----------------------------------------------------------


  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Variable_RestartPut

    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statement
    !
    
    if( .not. DSIce_IO_Restart_isOutputTiming(CurrentTime) ) return

    call DSIce_IO_Restart_HistPut( 'SIceCon', xya_SIceCon(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'SIceConB', xya_SIceCon(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistPut( 'SIceSfcTemp', xya_SIceSfcTemp(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'SIceSfcTempB', xya_SIceSfcTemp(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistPut( 'SIceEn', xyza_SIceEn(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'SIceEnB', xyza_SIceEn(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistPut( 'SIceTemp', xyza_SIceTemp(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'SIceTempB', xyza_SIceTemp(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistPut( 'IceThick', xya_IceThick(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'IceThickB', xya_IceThick(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistPut( 'SnowThick', xya_SnowThick(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistPut( 'SnowThickB', xya_SnowThick(IS:IE,JS:JE, TIMELV_ID_B) )

  end subroutine DSIce_Admin_Variable_RestartPut
  
  !-------------------------------------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Variable_HistGet()

    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statement
    !

    call DSIce_IO_Restart_HistGet( 'SIceCon', xya_SIceCon(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'SIceConB', xya_SIceCon(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistGet( 'SIceSfcTemp', xya_SIceSfcTemp(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'SIceSfcTempB', xya_SIceSfcTemp(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistGet( 'SIceEn', xyza_SIceEn(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'SIceEnB', xyza_SIceEn(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistGet( 'SIceTemp', xyza_SIceTemp(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'SIceTempB', xyza_SIceTemp(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistGet( 'IceThick', xya_IceThick(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'IceThickB', xya_IceThick(IS:IE,JS:JE, TIMELV_ID_B) )
    call DSIce_IO_Restart_HistGet( 'SnowThick', xya_SnowThick(IS:IE,JS:JE, TIMELV_ID_N) )
    call DSIce_IO_Restart_HistGet( 'SnowThickB', xya_SnowThick(IS:IE,JS:JE, TIMELV_ID_B) )

  end subroutine DSIce_Admin_Variable_HistGet
  
end module DSIce_Admin_Variable_mod

