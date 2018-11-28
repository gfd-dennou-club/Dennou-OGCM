!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief The module which provides some  operators in calculus using fintie volume method. 
!!
!! 
!! @author Yuta Kawai
!! @since 2017
!!
!!
module CalculusDriver_mod

  ! モジュール引用; Use statements
  !

  !* gtool
  !

  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: CalculusDriver_Init, CalculusDriver_Final

  abstract interface
     function IntBtmToTop( xyz_q, xyz_e3 ) result( xy_VIntVal )
       use dc_types, only: DP
       implicit none
       real(DP), intent(in) :: xyz_q(:,:,:)
       real(DP), intent(in) :: xyz_e3(:,:,:)
       real(DP)  :: xy_VIntVal(size(xyz_q,1),size(xyz_q,2))       
     end function IntBtmToTop
  end interface

  procedure(IntBtmToTop), pointer :: Calculus_IntBtmToTop => NULL()
  public :: Calculus_IntBtmToTop
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'CalculusDriver_mod' !< Module Name

  logical :: isInitialzed = .false.

  
contains

  subroutine CalculusDriver_Init( solver_name )

    ! モジュール引用; Use statements
    !

    use VSpmUtil_mod, only: &
         & VSpm_IntBtmToTop_xyz
    
    use VFvmUtil_mod, only: &
         & VFvm_IntBtmToTop_xyz
    
    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: solver_name

    ! 実行文; Executable statement
    !
    
    select case(solver_name)
    case('hspm_vspm')
       Calculus_IntBtmToTop => VSpm_IntBtmToTop_xyz
    case('hspm_vfvm')
       Calculus_IntBtmToTop => VFvm_IntBtmToTop_xyz
    case default
       call MessageNotify('E', module_name,                       &
            & "Unexcepted solver name '%a' is specified. Check!", &
            & ca = (/ solver_name /) )
    end select
    
    isInitialzed = .true.
    
  end subroutine CalculusDriver_Init

  !------------------------------

  subroutine CalculusDriver_Final()

    isInitialzed = .false.
    
  end subroutine CalculusDriver_Final

  !------------------------------

  
end module CalculusDriver_mod
