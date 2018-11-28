!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module OptionParser_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: STRING 

  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: OptionParser_Init, OptionParser_Final
  public :: OptionParser_GetInfo

  ! 公開変数
  ! Public variables
  !
  character(*), parameter, public :: DEFAULT_CONFIG_NML = "defaultConfig.nml"

  ! 非公開変数
  ! Private variables

  character(*), parameter:: module_name = 'OptionParser_mod'
                              ! モジュールの名称. 
                              ! Module name  

contains

!>
!!
!!
subroutine OptionParser_Init()

  ! 実行文; Executable statements

end subroutine OptionParser_Init

!>
!!
!!
subroutine OptionParser_Final()

  ! 実行文; Executable statements

end subroutine OptionParser_Final


!> @brief 
!!
!!
subroutine OptionParser_getInfo( &
     & configNmlName, defaultConfigNml )
  
  ! モジュール引用; Use statement
  !
  use dc_args
  use dc_string

  !
  ! 宣言文; Declaration statement
  !
  character(*), intent(out) :: configNmlName
  character(*), intent(in), optional :: defaultConfigNml

  ! 局所変数
  ! Local variables
  !
  type(Args) :: arg
  logical :: optNmlName
  character(STRING) :: valNmlName

  ! 実行文; Executable statement
  !

  call DCArgsOpen(arg)
  call DCArgsOption(arg, StoA("--N",  "--namelist"), optNmlName, valNmlName, &
       help="Specify the namelist name")
  
  call DCArgsDebug(arg)
  call DCArgsHelp(arg)
  call DCArgsStrict(arg)

  call DCArgsClose(arg)

  if( optNmlName ) then
     configNmlName = valNmlName
  else

     ! If defaultConfigNml is presented, 
     ! configNmlName is set it. 
     if( present(defaultConfigNml) ) then
        configNmlName = defaultConfigNml
     else
        configNmlName = DEFAULT_CONFIG_NML
     end if

  end if

end subroutine OptionParser_getInfo

end module OptionParser_mod

