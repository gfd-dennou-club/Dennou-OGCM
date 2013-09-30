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

  use dc_types, only: &
       & DP, STRING 

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

  !>
  type(volScalarField), public, save :: zc_lyrThickB
  type(volScalarField), public, save :: zc_lyrThickN
  type(volScalarField), public, save :: zc_lyrThickA
  type(volScalarField), public, save :: zc_lyrThickBasic

  !>
  type(volScalarField), public, save :: c_totDepthBasic

  !>
  type(surfaceScalarField), public, save :: ze_hNormVelB
  type(surfaceScalarField), public, save :: ze_hNormVelN
  type(surfaceScalarField), public, save :: ze_hNormVelA

  !>
  type(volScalarField), public, save, allocatable :: zc_TracersB(:)
  type(volScalarField), public, save, allocatable :: zc_TracersN(:)
  type(volScalarField), public, save, allocatable :: zc_TracersA(:)
  integer, public, save :: TracerNum
  integer, public, save :: SaltTracerID
  integer, public, save :: PTempTracerID


  type(volScalarField), public, save :: rc_vNormVel
  type(volScalarField), public, save :: c_SurfPress
  type(volScalarField), public, save :: zc_Press
  type(volScalarField), public, save :: zc_Dens

  type(volScalarField), public, save :: zc_Zmid
  real(DP), public, save :: refDens


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
         & plMesh, nVzLyr, &
         & nVrLyr, vHaloSize

    ! 局所変数
    ! Local variable
    !
    integer :: n

    ! 実行文; Executable statements
    !


    call GeometricField_Init(zc_lyrThickB, plMesh, "zc_lyrThick", &
         & "layer thickness", "m", vHaloSize=vHaloSize)
    call GeometricField_Init(zc_lyrThickN, plMesh, "zc_lyrThick", &
         & "layer thickness", "m", vHaloSize=vHaloSize)
    call GeometricField_Init(zc_lyrThickA, plMesh, "zc_lyrThick", &
         & "layer thickness", "m", vHaloSize=vHaloSize)

    call GeometricField_Init(zc_lyrThickBasic, plMesh, "zc_lyrThickBasic", &
         & "layer thickness", "m", vHaloSize=vHaloSize)

    
    call GeometricField_Init(c_totDepthBasic, plMesh, "c_totBasicDepth", &
         & "total unperturbed depth", "m", vLayerNum=1)


    call GeometricField_Init(ze_hNormVelB, plMesh, "ze_hNormVel", &
         & "horizontal velocity which is normal component to edge", "m*s-1", vHaloSize=vHaloSize)
    call GeometricField_Init(ze_hNormVelN, plMesh, "ze_hNormVel", &
         & "horizontal velocity which is normal component to edge", "m*s-1", vHaloSize=vHaloSize)
    call GeometricField_Init(ze_hNormVelA, plMesh, "ze_hNormVel", &
         & "horizontal velocity which is normal component to edge", "m*s-1", vHaloSize=vHaloSize)

    
    call GeometricField_Init(rc_vNormVel, plMesh, "rc_vNormVel", &
         & "vertical velocity which is normal component to edge", "m*s-1", vLayerNum=nVrLyr)
    
    call GeometricField_Init(c_SurfPress, plMesh, "c_SurfPress", &
         & "surface pressure", "kg*m-1*s-2", vLayerNum=1)

    call GeometricField_Init(zc_Press, plMesh, "zc_Press", &
         & "pressure", "kg*m-1*s-2")

    call GeometricField_Init(zc_Dens, plMesh, "zc_Dens", &
         & "density", "kg*m-3")


    TracerNum = 2
    allocate(zc_TracersB(TracerNum), zc_TracersN(TracerNum), zc_TracersA(TracerNum))

    SaltTracerID = 1
    call GeometricField_Init(zc_TracersB(1), plMesh, "zc_Salt", &
         & "salinity", "kg m-3")
    call GeometricField_Init(zc_TracersN(1), plMesh, "zc_Salt", &
         & "salinity", "kg m-3")
    call GeometricField_Init(zc_TracersA(1), plMesh, "zc_Salt", &
         & "salinity", "kg m-3")

    PTempTracerID = 2
    call GeometricField_Init(zc_TracersB(2), plMesh, "zc_PTemp", &
         & "potential temperature", "k")
    call GeometricField_Init(zc_TracersN(2), plMesh, "zc_PTemp", &
         & "potential temperature", "k")
    call GeometricField_Init(zc_TracersA(2), plMesh, "zc_PTemp", &
         & "potential temperature", "k")

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

    ! 局所変数
    ! Local variable
    !
    integer :: n

    ! 実行文; Executable statements
    !

    call GeometricField_Final(zc_lyrThickB)
    call GeometricField_Final(zc_lyrThickN)
    call GeometricField_Final(zc_lyrThickA)
    call GeometricField_Final(zc_lyrThickBasic)
    call GeometricField_Final(c_totDepthBasic)
    call GeometricField_Final(ze_hNormVelB)
    call GeometricField_Final(ze_hNormVelN)
    call GeometricField_Final(ze_hNormVelA)
    call GeometricField_Final(rc_vNormVel)
    call GeometricField_Final(c_SurfPress)
    call GeometricField_Final(zc_Press)
    call GeometricField_Final(zc_Dens)
    call GeometricField_Final(zc_Zmid)

    do n=1, TracerNum
       call GeometricField_Final(zc_TracersB(n))
       call GeometricField_Final(zc_TracersN(n))
       call GeometricField_Final(zc_TracersA(n))
    end do
    deallocate(zc_TracersB, zc_TracersN, zc_TracersA)

  end subroutine VariableSet_Final

end module VariableSet_mod

