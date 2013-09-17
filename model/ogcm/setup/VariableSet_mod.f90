!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module VariableSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types

  use GeometricField_mod, only: &
       & volScalarField, pointScalarField, surfaceScalarField

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: VariableSet_Init, VariableSet_Final

  ! 公開変数
  ! Public variable
  !
  type(volScalarField), public, save :: zc_lyrThick
  type(surfaceScalarField), public, save :: ze_hNormVel
  type(volScalarField), public, save :: rc_vNormVel
  type(volScalarField), public, save :: zc_Salt
  type(volScalarField), public, save :: zc_Zmid


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VariableSet_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine VariableSet_Init()

    ! モジュール引用; Use statements
    !
    use GeometricField_mod, only: &
         & GeometricField_Init

    use GridSet_mod, only: &
         & plMesh

    ! 実行文; Executable statements
    !

    call GeometricField_Init(zc_lyrThick, plMesh, &
         & "zc_lyrThick", "layer thickness", "m")

    call GeometricField_Init(ze_hNormVel, plMesh, "ze_hNormVel", &
         & "horizontal velocity which is normal component to edge", "m s-1")
    
    call GeometricField_Init(rc_vNormVel, plMesh, "rc_vNormVel", &
         & "vertical velocity which is normal component to edge", "m s-1")
    
    call GeometricField_Init(zc_Salt, plMesh, "zc_Salt", &
         & "salinity", "kg m-3")

    call GeometricField_Init(zc_Zmid, plMesh, "zc_Zmid", &
         & "layer mid-depth location", "m")

  end subroutine VariableSet_Init

  !>
  !!
  !!
  subroutine VariableSet_Final()

    ! モジュール引用; Use statements
    !
    use GeometricField_mod, only: &
         & GeometricField_Final

    ! 実行文; Executable statements
    !

    call GeometricField_Final(zc_lyrThick)
    call GeometricField_Final(ze_hNormVel)
    call GeometricField_Final(rc_vNormVel)
    call GeometricField_Final(zc_Salt)
    call GeometricField_Final(zc_Zmid)

  end subroutine VariableSet_Final

end module VariableSet_mod

