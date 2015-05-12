!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module PhysicsDriver_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP

  use dc_message, only: &
       & MessageNotify

  use DataFileSet_mod, only: &
       & DataFileSet

  use GovernEqSet_mod, only: &
       & GOVERNEQSET_PHYSICS_CONVADJUST_NAME, &
       & GOVERNEQSET_PHYSICS_EDDYMIX_NAME, &
       & isPhysicsCompActivated

  use SGSEddyMixing_mod, only: &
       & SGSEddyMixing_Init, SGSEddyMixing_Final, &
       & SGSEddyMixing_Output, SGSEddyMixing_PrepareOutput

  use SGSConvAdjust_mod, only: &
       & SGSConvAdjust_Init, SGSConvAdjust_Final

  use SGSSlowConvAdjust_mod, only: &
       & SGSSlowConvAdjust_Init, SGSSlowConvAdjust_Final

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: PhysicsDriver_Init, PhysicsDriver_Final
  public :: PhysicsDriver_OutputData

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'PhysicsDriver_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine PhysicsDriver_Init(datFile, configNmlFileName)

    use TemporalIntegSet_mod, only: &
         & RestartTime, EndTime, CurrentTime, DelTime


    type(DataFileSet), intent(in) :: datFile
    character(*), intent(in) :: configNmlFileName
    
    ! 実行文; Executable statements
    !
    
    if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) then
       call SGSEddyMixing_Init(configNmlFileName=configNmlFileName)
       call SGSEddyMixing_PrepareOutput( &
            & RestartTime, EndTime, datFile%outputIntrvalSec, datFile%FilePrefix &
            & )
    end if

    if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_CONVADJUST_NAME)) then
       call SGSConvAdjust_Init()
       call SGSSlowConvAdjust_Init(GCMTimeStep=DelTime)
    end if
  end subroutine PhysicsDriver_Init

  !>
  !!
  !!
  subroutine PhysicsDriver_Final()

    ! 実行文; Executable statements
    !
    
    if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) &
         & call SGSEddyMixing_Final()

    if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_CONVADJUST_NAME)) then
       call SGSConvAdjust_Final()
!       call SGSSlowConvAdjust_Final()
    end if

  end subroutine PhysicsDriver_Final

  !> @brief 
  !!
  !!
  subroutine PhysicsDriver_OutputData(datFile)

    use TemporalIntegSet_mod, only: &
         & CurrentTime

    use DataFileSet_mod, only: &
         & DataFileSet, DataFileSet_isOutputTiming

    ! 宣言文; Declaration statement
    !
    type(DataFileSet), intent(in) :: datFile
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if( .not. DataFileSet_isOutputTiming(datFile, CurrentTime) ) return 

    call MessageNotify('M', module_name, &
         & "Output data of variables in physics packages." )

    if(isPhysicsCompActivated(GOVERNEQSET_PHYSICS_EDDYMIX_NAME)) then
       call SGSEddyMixing_Output()
    end if

  end subroutine PhysicsDriver_OutputData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module PhysicsDriver_mod

