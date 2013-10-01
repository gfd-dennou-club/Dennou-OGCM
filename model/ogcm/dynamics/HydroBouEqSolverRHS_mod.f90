!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module HydroBouEqSolverRHS_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING 

  use GeometricField_mod, only: &
       & volScalarField, pointScalarField, surfaceScalarField, &
       & GeometricField_Init, GeometricField_Final, &
       & assignment(=), operator(+), operator(-), operator(*), operator(/), &
       & DeepCopy
  
  use fvCalculus_mod, only: &
       & curl, grad, div

  use Constants_mod, only: &
       & Omega, Grav

  use GridSet_mod, only: &
       & GridSet_getLocalMeshInfo, &
       & nVzLyr, nVrLyr, vHaloSize

  use VariableSet_mod, only: &
       & refDens

  use CGridFieldDataUtil_mod, only: &
       & CalcKineticEnergy, CalcPVFlux, &
       & e_c

  use VGridFieldDataUtil_mod, only: &
       & z_r, r_z, delz

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HydroBouEqSolverRHS_Init, HydroBouEqSolverRHS_Final
  public :: calc_RHShNormVel, calc_RHSTracer, calc_RHSLyrThick

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HydroBouEqSolverRHS_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine HydroBouEqSolverRHS_Init()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolverRHS_Init

  !>
  !!
  !!
  subroutine HydroBouEqSolverRHS_Final()

    ! 実行文; Executable statements
    !

  end subroutine HydroBouEqSolverRHS_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  !> @brief 
  !!
  !!
  subroutine calc_RHShNormVel( ze_RHShNormVel, &
       & ze_hNormVel, rc_vNormVel, zc_lyrThick, zc_ZMid, zc_Press, zc_Dens )
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(inout) :: ze_RHShNormVel
    type(surfaceScalarField), intent(in) :: ze_hNormVel
    type(volScalarField), intent(in) :: rc_vNormVel
    type(volScalarField), intent(in) :: zc_lyrThick
    type(volScalarField), intent(in) :: zc_Zmid
    type(volScalarField), intent(in) :: zc_Press
    type(volScalarField), intent(in) :: zc_Dens
    
    ! 実行文; Executable statement
    !

    ! ddt(u) = - (1/rho0)*( grad(P) + g*rho*grad(zmid) )
    !          - grad(K) - PVFlux - w*ddz(u)
    !
    ze_RHShNormVel = &
         & (-1d0/refDens)*( grad(zc_Press)  + Grav*e_c(zc_Dens)*grad(zc_Zmid) ) &
         & - grad( CalcKineticEnergy(ze_hNormVel) ) &
         & - CalcPVFlux(zc_lyrThick, ze_hNormVel)   &
         & - z_r( e_c(rc_vNormVel) * delz(ze_hNormVel, zc_lyrThick) ) &
         & - z_r(delz(ze_hNormVel, zc_lyrThick))

  end subroutine calc_RHShNormVel

  !> @brief 
  !!
  !!
  subroutine calc_RHSLyrThick( zc_RHSlyrThick, &
       & ze_massHFlux, rc_vNormVel )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: zc_RHSlyrThick 
    type(surfaceScalarField), intent(in) :: ze_massHFlux
    type(volScalarField), intent(in) :: rc_vNormVel
    
    ! 局所変数
    ! Local variables
    !
    integer :: i, nCell
    
    ! 実行文; Executable statement
    !

    call GridSet_getLocalMeshInfo(ze_massHFlux%mesh, nCell=nCell)

    zc_RHSlyrThick = div(ze_massHFlux)
    
    !$omp parallel do
    do i=1, nCell
       zc_RHSlyrThick%data%v_(1:nVzLyr,i) = &
            &   rc_vNormVel%data%v_(2:nVzLyr+1,i) - rc_vNormVel%data%v_(1:nVzLyr,i) &
            & - zc_RHSlyrThick%data%v_(1:nVzLyr,i) 
    end do

  end subroutine calc_RHSLyrThick

  !> @brief 
  !!
  !!
  subroutine calc_RHSTracer( zc_RHSTracer, &
       & zc_Tracer, ze_massHFlux, rc_vNormVel )
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(inout) :: zc_RHSTracer 
    type(volScalarField), intent(in) :: zc_Tracer
    type(surfaceScalarField), intent(in) :: ze_massHFlux
    type(volScalarField), intent(in) :: rc_vNormVel
    
    ! 局所変数
    ! Local variables
    !
    integer :: i, nCell
    
    ! 実行文; Executable statement
    !

    call GridSet_getLocalMeshInfo(ze_massHFlux%mesh, nCell=nCell)
    

    zc_RHSTracer = div( e_c(zc_Tracer)*ze_massHFlux )

    !$omp parallel do
    do i=1, nCell
       zc_RHSTracer%data%v_(1:nVzLyr,i) = &
            & rc_vNormVel%data%v_(2:nVzLyr+1,i) - rc_vNormVel%data%v_(1:nVzLyr,i) &
            & - zc_RHSTracer%data%v_(1:nVzLyr,i) 
    end do

  end subroutine calc_RHSTracer

end module HydroBouEqSolverRHS_mod
