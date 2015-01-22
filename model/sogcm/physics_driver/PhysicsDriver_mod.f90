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
       & SGSEddyMixType, &
       & SGSConvAdjustType, &
       & GOVERNEQSET_SGSCONVADJUST_INSTANT, &
       & GOVERNEQSET_SGSCONVADJUST_SLOW, &
       & isPhysicsCompActived

  use SGSEddyMixing_mod, only: &
       & SGSEddyMixing_Init, SGSEddyMixing_Final, &
       & SGSEddyMixing_PrepareOutput, &
       & SGSEddyMixing_Output

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
  subroutine PhysicsDriver_Init(datFile)

    use TemporalIntegSet_mod, only: &
         & RestartTime, EndTime, CurrentTime, DelTime


    type(DataFileSet), intent(in) :: datFile

    ! 実行文; Executable statements
    !

    if(isPhysicsCompActived(SGSEddyMixType)) then
       call SGSEddyMixing_Init(SGSEddyMixType, 1d3)
       call SGSEddyMixing_PrepareOutput( &
            & RestartTime, EndTime, datFile%outputIntrvalSec, datFile%FilePrefix &
            & )
    end if

    if(isPhysicsCompActived(SGSConvAdjustType)) then
       select case(SGSConvAdjustType)
       case(GOVERNEQSET_SGSCONVADJUST_INSTANT)
          call SGSConvAdjust_Init()
       case(GOVERNEQSET_SGSCONVADJUST_SLOW)
          call SGSSlowConvAdjust_Init(GCMTimeStep=DelTime)
       end select
    end if
  end subroutine PhysicsDriver_Init

  !>
  !!
  !!
  subroutine PhysicsDriver_Final()

    ! 実行文; Executable statements
    !
    
    if(isPhysicsCompActived(SGSEddyMixType)) &
         & call SGSEddyMixing_Final()

    if(isPhysicsCompActived(SGSConvAdjustType)) then
       select case(SGSConvAdjustType)
       case(GOVERNEQSET_SGSCONVADJUST_INSTANT)
          call SGSConvAdjust_Final()
       case(GOVERNEQSET_SGSCONVADJUST_SLOW)
          call SGSSlowConvAdjust_Final()
       end select
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

    if(isPhysicsCompActived(SGSEddyMixType)) then
       call SGSEddyMixing_Output()
    end if

  end subroutine PhysicsDriver_OutputData

end module PhysicsDriver_mod

