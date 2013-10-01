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

  use TemporalIntegSet_mod, only: &
       & nLongTimeLevel

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

  type, public :: VariableSet
     !>
     type(volScalarField) :: zc_lyrThick(nLongTimeLevel)
     type(volScalarField) :: zc_lyrThickBasic

     !>
     type(volScalarField) :: c_totDepthBasic

     !>
     type(surfaceScalarField) :: ze_hNormVel(nLongTimeLevel)

     !>
     type(volScalarField), pointer :: zc_Tracers(:,:)


     type(volScalarField) :: rc_vNormVel
     type(volScalarField) :: c_SurfPress
     type(volScalarField) :: zc_Press
     type(volScalarField) :: zc_Dens
     
     type(volScalarField) :: zc_Zmid
  end type VariableSet

  integer, public, save :: TracerNum
  integer, public, save :: SaltTracerID
  integer, public, save :: PTempTracerID
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
  subroutine VariableSet_Init(this, plMesh)

    ! モジュール引用; Use statements
    !
    use PolyMesh_mod, only: &
         & PolyMesh

    use GeometricField_mod, only: &
         & GeometricField_Init

    use GridSet_mod, only: &
         & nVzLyr, nVrLyr, vHaloSize

    ! 宣言文; Declaration statement
    !
    type(VariableSet), intent(inout) :: this
    type(PolyMesh), intent(in) :: plMesh

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !

    do tl=1, nLongTimeLevel
       call GeometricField_Init(this%zc_lyrThick(tl), plMesh, "zc_lyrThick", &
            & "layer thickness", "m", vHaloSize=vHaloSize)
    end do

    call GeometricField_Init(this%zc_lyrThickBasic, plMesh, "zc_lyrThickBasic", &
         & "layer thickness", "m", vHaloSize=vHaloSize)

    call GeometricField_Init(this%c_totDepthBasic, plMesh, "c_totBasicDepth", &
         & "total unperturbed depth", "m", vLayerNum=1)

    do tl=1, nLongTimeLevel
       call GeometricField_Init(this%ze_hNormVel(tl), plMesh, "ze_hNormVel", &
            & "horizontal velocity which is normal component to edge", "m*s-1", vHaloSize=vHaloSize)
    end do

    
    call GeometricField_Init(this%rc_vNormVel, plMesh, "rc_vNormVel", &
         & "vertical velocity which is normal component to edge", "m*s-1", vLayerNum=nVrLyr)
    
    call GeometricField_Init(this%c_SurfPress, plMesh, "c_SurfPress", &
         & "surface pressure", "kg*m-1*s-2", vLayerNum=1)

    call GeometricField_Init(this%zc_Press, plMesh, "zc_Press", &
         & "pressure", "kg*m-1*s-2")

    call GeometricField_Init(this%zc_Dens, plMesh, "zc_Dens", &
         & "density", "kg*m-3")


    !
    TracerNum = 2
    SaltTracerID = 1
    PTempTracerID = 2
    
    allocate( this%zc_Tracers(TracerNum, nLongTimeLevel) )
    do tl=1, nLongTimeLevel
       call GeometricField_Init(this%zc_Tracers(SaltTracerID,tl), plMesh, "zc_Salt", &
            & "salinity", "kg m-3")
    
       call GeometricField_Init(this%zc_Tracers(PTempTracerID,tl), plMesh, "zc_PTemp", &
            & "potential temperature", "k")
    end do

    call GeometricField_Init(this%zc_Zmid, plMesh, "zc_Zmid", &
         & "layer mid-depth location", "m")

  end subroutine VariableSet_Init

  !>
  !!
  !!
  subroutine VariableSet_Final(this)

    ! モジュール引用; Use statements
    !
    use GeometricField_mod, only: &
         & GeometricField_Final

    ! 宣言文; Declaration statement
    !
    type(VariableSet), intent(inout) :: this


    ! 局所変数
    ! Local variable
    !
    integer :: n, tl
    
    ! 実行文; Executable statements
    !

    do tl=1, nLongTimeLevel
       call GeometricField_Final( this%zc_lyrThick(tl) )
       call GeometricField_Final( this%ze_hNormVel(tl) )
       do n=1, TracerNum
          call GeometricField_Final( this%zc_Tracers(tl,n) )
       end do
    end do

    deallocate(this%zc_Tracers)

    call GeometricField_Final(this%c_totDepthBasic)
    call GeometricField_Final(this%rc_vNormVel)
    call GeometricField_Final(this%c_SurfPress)
    call GeometricField_Final(this%zc_Press)
    call GeometricField_Final(this%zc_Dens)
    call GeometricField_Final(this%zc_Zmid)

  end subroutine VariableSet_Final

end module VariableSet_mod

