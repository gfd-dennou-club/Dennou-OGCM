!-------------------------------------------------------------
! Copyright (c) 2015-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_Variable_mod 

  ! モジュール引用; Use statements
  !

  !* gtool

  use dc_types, only: &
       & DP, STRING, TOKEN 

  use dc_message, only: &
       & MessageNotify
  
  !* Dennou-OGCM

  use DOGCM_Admin_TInteg_mod, only: &
       & nLongTimeLevel,                        &
       & TIMELV_ID_A, TIMELV_ID_N, TIMELV_ID_B, &
       & CurrentTime

  use DOGCM_Admin_Grid_mod, only: &
       & IA, JA, KA,              &
       & IS, IE, IM,              &
       & JS, JE, JM,              &
       & KS, KE, KM,              &
       & xyz_Z
    

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_Admin_Variable_Init, DOGCM_Admin_Variable_Final

  public :: DOGCM_Admin_Variable_AdvanceTStep

  public :: DOGCM_Admin_Variable_HistPut
  public :: DOGCM_Admin_Variable_RestartPut

  public :: DOGCM_Admin_Variable_HistGet
  
  ! 公開変数
  ! Public variable
  !

  
  real(DP), public, allocatable :: xyza_U(:,:,:,:)    !< The velocity component in i-direction
  real(DP), public, allocatable :: xyza_V(:,:,:,:)    !< The velocity component in j-direction
  real(DP), public, allocatable :: xyza_OMG(:,:,:,:)  !< The velocity component in k-direction (dia-surface velocity component)
  real(DP), public, allocatable :: xya_SSH(:,:,:)     !< The sea surface height
  real(DP), public, allocatable :: xyza_H(:,:,:,:)    !< The vertical scale factor

  real(DP), public, allocatable :: xyza_HydPres(:,:,:,:) 
  real(DP), public, allocatable :: xya_SfcPres(:,:,:) 

  integer, public, parameter :: TRCID_PTEMP = 1
  integer, public, parameter :: TRCID_SALT  = 2
  integer, public, parameter :: TRC_TOT_NUM = 2
  real(DP), public, allocatable :: xyzaa_TRC(:,:,:,:,:)

  real(DP), public, allocatable :: xyz_VViscCoef(:,:,:)
  real(DP), public, allocatable :: xyz_VDiffCoef(:,:,:)

  real(DP), public, allocatable :: xyz_ConvIndex(:,:,:)

  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Admin_Variable_mod' !< Module Name

  logical :: isInitialzed = .false.
  
contains

  !>
  !!
  !!
  subroutine DOGCM_Admin_Variable_Init()

    ! モジュール引用; Use statements
    !

    
    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: TA

    ! 実行文; Executable statements
    !

    
    TA = nLongTimeLevel

    !- Allocate array -------------------------------------------

    call MessageNotify('M', module_name, "Allocate array to store data of variables..")    
    allocate( xyza_U(IA,JA,KA,TA), xyza_V(IA,JA,KA,TA), xyza_OMG(IA,JA,KA,TA) )
    allocate( xya_SSH(IA,JA,TA), xyza_H(IA,JA,KA,TA) )
    allocate( xyza_HydPres(IA,JA,KA,TA), xya_SfcPres(IA,JA,TA) )
    allocate( xyzaa_TRC(IA,JA,KA,TRC_TOT_NUM,TA) )
    allocate( xyz_VViscCoef(IA,JA,KA), xyz_VDiffCoef(IA,JA,KA) )
    allocate( xyz_ConvIndex(IA,JA,KA) )

    !- Regist output variables for history and restart --------------

    call MessageNotify('M', module_name, "Regist output variables  ..")
    call regist_OuputVariable()


    !------------------------------------
    
    isInitialzed = .true.

    
  end subroutine DOGCM_Admin_Variable_Init

  !--------------------------------------------------------------------------
  
  !>
  !!
  !!
  subroutine DOGCM_Admin_Variable_Final()

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variable
    !
    
    ! 実行文; Executable statements
    !

    if( isInitialzed ) then
       
       call MessageNotify('M', module_name, "Deallocate array to store data of variables..")
       deallocate( xyza_U, xyza_V, xyza_OMG )
       deallocate( xya_SSH, xyza_H )
       deallocate( xyza_HydPres, xya_SfcPres )
       deallocate( xyzaa_TRC )
       deallocate( xyz_VViscCoef, xyz_VDiffCoef )
       deallocate( xyz_ConvIndex )
       
    end if

  end subroutine DOGCM_Admin_Variable_Final

  !----------------------------------------------------------


  !> @brief 
  !!
  !!
  subroutine DOGCM_Admin_Variable_AdvanceTStep()

    ! 宣言文; Declaration statement
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
       xyza_U(:,:,:,n) = xyza_U(:,:,:,n-1)
       xyza_V(:,:,:,n) = xyza_V(:,:,:,n-1)
       xyza_OMG(:,:,:,n) = xyza_OMG(:,:,:,n-1)
       xyza_H(:,:,:,n) = xyza_H(:,:,:,n-1)
       xyza_HydPres(:,:,:,n) = xyza_HydPres(:,:,:,n-1)
       !$omp end workshare

       !$omp workshare
       xyzaa_TRC(:,:,:,:,n) = xyzaa_TRC(:,:,:,:,n-1)
       !$omp end workshare
       
       !$omp workshare
       xya_SSH(:,:,n) = xya_SSH(:,:,n-1)
       xya_SfcPres(:,:,n) = xya_SfcPres(:,:,n-1)
       !$omp end workshare

       !$omp end parallel
    end do

    
  end subroutine DOGCM_Admin_Variable_AdvanceTStep

  !----------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DOGCM_Admin_Variable_HistPut()

    use DOGCM_IO_History_mod, only:      &
       & DOGCM_IO_History_HistPut,       &
       & DOGCM_IO_History_IsOutputTiming

    use EOSDriver_mod, only: &
       & EOSDriver_Eval
    
    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statement
    !
    
    if( .not. DOGCM_IO_History_isOutputTiming(CurrentTime) ) return

    call DOGCM_IO_History_HistPut( 'U', xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_History_HistPut( 'V', xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_History_HistPut( 'OMG', xyza_OMG(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_History_HistPut( 'H', xyza_H(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_History_HistPut( 'SSH', xya_SSH(IS:IE,JS:JE, TIMELV_ID_N) )
    call DOGCM_IO_History_HistPut( 'PTemp', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TIMELV_ID_N) )
    call DOGCM_IO_History_HistPut( 'Salt', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TIMELV_ID_N) )

    call DOGCM_IO_History_HistPut( 'Z', xyz_Z(IS:IE,JS:JE,KS:KE) )

    call DOGCM_IO_History_HistPut( 'HydPres', xyza_HydPres(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_History_HistPut( 'SfcPres', xya_SfcPres(IS:IE,JS:JE, TIMELV_ID_B) )


    call DOGCM_IO_History_HistPut( 'ConvIndex', xyz_ConvIndex(IS:IE,JS:JE,KS:KE) )

    call DOGCM_IO_History_HistPut( 'VViscCoef', xyz_VViscCoef(IS:IE,JS:JE,KS:KE) )
    call DOGCM_IO_History_HistPut( 'VDiffCoef', xyz_VDiffCoef(IS:IE,JS:JE,KS:KE) )
    
    
  end subroutine DOGCM_Admin_Variable_HistPut

  !-------------------------------------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_Admin_Variable_RestartPut()

    use DOGCM_IO_Restart_mod, only:      &
       & DOGCM_IO_Restart_HistPut,       &
       & DOGCM_IO_Restart_IsOutputTiming

    ! 宣言文; Declaration statement
    !

    ! 実行文; Executable statement
    !
    
    if( .not. DOGCM_IO_Restart_isOutputTiming(CurrentTime) ) return

    call DOGCM_IO_Restart_HistPut( 'U', xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'UB', xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistPut( 'V', xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'VB', xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistPut( 'OMG', xyza_OMG(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'OMGB', xyza_OMG(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistPut( 'H', xyza_H(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'HB', xyza_H(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistPut( 'SSH', xya_SSH(IS:IE,JS:JE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'SSHB', xya_SSH(IS:IE,JS:JE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistPut( 'PTemp', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'PTempB', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TIMELV_ID_B) )

    call DOGCM_IO_Restart_HistPut( 'Salt', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'SaltB', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TIMELV_ID_B) )

    call DOGCM_IO_Restart_HistPut( 'Z', xyz_Z(IS:IE,JS:JE,KS:KE) )

    call DOGCM_IO_Restart_HistPut( 'HydPres', xyza_HydPres(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'HydPresB', xyza_HydPres(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )    
    call DOGCM_IO_Restart_HistPut( 'SfcPres', xya_SfcPres(IS:IE,JS:JE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistPut( 'SfcPresB', xya_SfcPres(IS:IE,JS:JE, TIMELV_ID_B) )
    
  end subroutine DOGCM_Admin_Variable_RestartPut

  !-------------------------------------------------------------------------------------------
  
  subroutine DOGCM_Admin_Variable_HistGet()

    use DOGCM_IO_Restart_mod, only:      &
       & DOGCM_IO_Restart_HistGet

    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statement
    !
    
    call DOGCM_IO_Restart_HistGet( 'U', xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'UB', xyza_U(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistGet( 'V', xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'VB', xyza_V(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistGet( 'OMG', xyza_OMG(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'OMGB', xyza_OMG(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistGet( 'H', xyza_H(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'HB', xyza_H(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistGet( 'SSH', xya_SSH(IS:IE,JS:JE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'SSHB', xya_SSH(IS:IE,JS:JE, TIMELV_ID_B) )
    call DOGCM_IO_Restart_HistGet( 'PTemp', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'PTempB', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_PTEMP,TIMELV_ID_B) )

    call DOGCM_IO_Restart_HistGet( 'Salt', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'SaltB', xyzaa_TRC(IS:IE,JS:JE,KS:KE, TRCID_SALT,TIMELV_ID_B) )

    call DOGCM_IO_Restart_HistGet( 'Z', xyz_Z(IS:IE,JS:JE,KS:KE) )

    call DOGCM_IO_Restart_HistGet( 'HydPres', xyza_HydPres(IS:IE,JS:JE,KS:KE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'HydPresB', xyza_HydPres(IS:IE,JS:JE,KS:KE, TIMELV_ID_B) )    
    call DOGCM_IO_Restart_HistGet( 'SfcPres', xya_SfcPres(IS:IE,JS:JE, TIMELV_ID_N) )
    call DOGCM_IO_Restart_HistGet( 'SfcPresB', xya_SfcPres(IS:IE,JS:JE, TIMELV_ID_B) )
    
    
  end subroutine DOGCM_Admin_Variable_HistGet

  !-------- Private subroutines -------------------------------------

  subroutine regist_OuputVariable()

    use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar

    use DOGCM_IO_Restart_mod, only: &
       & DOGCM_IO_Restart_RegistVar


    !- Regist variables for history
    
    call DOGCM_IO_History_RegistVar( 'U', "IJKT", "velocity (i axis)", "m/s"   )
    call DOGCM_IO_History_RegistVar( 'V', "IJKT", "velocity (j axis)", "m/s"   )
    call DOGCM_IO_History_RegistVar( 'OMG', "IJKT", "velocity (k axis)", "m/s" )
    call DOGCM_IO_History_RegistVar( 'SSH', "IJT", "sea surface height", "m"   )
    call DOGCM_IO_History_RegistVar( 'H', "IJKT",                              &
         & "vertical scale factor corresponding to layer thickness", "m"       )
    call DOGCM_IO_History_RegistVar( 'PTemp', "IJKT", "potential temperature", "K" )
    call DOGCM_IO_History_RegistVar( 'Salt', "IJKT", "salinity", "psu" )

    call DOGCM_IO_History_RegistVar( 'Z', "IJKT", "depth from sea surface", "m" )

    call DOGCM_IO_History_RegistVar( 'SfcPres', "IJT", "sea surface pressure", "Pa" )
    call DOGCM_IO_History_RegistVar( 'HydPres', "IJKT", "hydrostatic pressure", "Pa" )
    
    call DOGCM_IO_History_RegistVar( 'VViscCoef', "IJKT", "vertical viscosity", "m2/s" )
    call DOGCM_IO_History_RegistVar( 'VDiffCoef', "IJKT", "vertical diffusivity", "m2/s" )

    call DOGCM_IO_History_RegistVar( 'ConvIndex', "IJKT", &
         & "The number of calling convective adjustment per time step", "times per time step" )

    !- Regist variables for restart

    call DOGCM_IO_Restart_RegistVar( 'U', "IJKT", "velocity (i axis)", "m/s"   )
    call DOGCM_IO_Restart_RegistVar( 'UB', "IJKT", "velocity (i axis)", "m/s"   )
    call DOGCM_IO_Restart_RegistVar( 'V', "IJKT", "velocity (j axis)", "m/s"   )
    call DOGCM_IO_Restart_RegistVar( 'VB', "IJKT", "velocity (j axis)", "m/s"   )
    call DOGCM_IO_Restart_RegistVar( 'OMG', "IJKT", "velocity (k axis)", "m/s" )
    call DOGCM_IO_Restart_RegistVar( 'OMGB', "IJKT", "velocity (k axis)", "m/s" )
    call DOGCM_IO_Restart_RegistVar( 'SSH', "IJT", "sea surface height", "m"   )
    call DOGCM_IO_Restart_RegistVar( 'SSHB', "IJT", "sea surface height", "m"   )
    call DOGCM_IO_Restart_RegistVar( 'H', "IJKT",                              &
         & "vertical scale factor corresponding to layer thickness", "m"       )
    call DOGCM_IO_Restart_RegistVar( 'HB', "IJKT",                              &
         & "vertical scale factor corresponding to layer thickness", "m"       )
    call DOGCM_IO_Restart_RegistVar( 'PTemp', "IJKT", "potential temperature", "K" )
    call DOGCM_IO_Restart_RegistVar( 'PTempB', "IJKT", "potential temperature", "K" )
    call DOGCM_IO_Restart_RegistVar( 'Salt', "IJKT", "salinity", "psu" )
    call DOGCM_IO_Restart_RegistVar( 'SaltB', "IJKT", "salinity", "psu" )

    call DOGCM_IO_Restart_RegistVar( 'SfcPres', "IJT", "sea surface pressure", "Pa" )
    call DOGCM_IO_Restart_RegistVar( 'SfcPresB', "IJT", "sea surface pressure", "Pa" )
    call DOGCM_IO_Restart_RegistVar( 'HydPres', "IJKT", "hydrostatic pressure", "Pa" )
    call DOGCM_IO_Restart_RegistVar( 'HydPresB', "IJKT", "hydrostatic pressure", "Pa" )
    
    call DOGCM_IO_Restart_RegistVar( 'Z', "IJKT", "depth from sea surface", "m" )
    
    !-------------------------------------------------------------
    
  end subroutine regist_OuputVariable
  
end module DOGCM_Admin_Variable_mod

