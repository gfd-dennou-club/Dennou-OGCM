!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module VGridFieldDataUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

  use GridSet_mod, only: &
       & plMesh, fvmInfo, &
       & nCell, nEdge, nVertex, &
       & nVzLyr, nVrLyr, vHaloSize

  use GeometricField_mod, only: &
       & volScalarField, surfaceScalarField, &
       & GeometricField_Init, GeometricField_Final, Release, &
       & operator(+), operator(-), operator(/), operator(*)

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  interface r_z
     module procedure rc_zc
     module procedure re_ze
  end interface r_z

  interface z_r
     module procedure zc_rc
     module procedure ze_re
  end interface z_r

  interface delz
     module procedure delz_ze
     module procedure delz_zc
  end interface delz

  interface verticalInt
     module procedure verticalInt_e
     module procedure verticalInt_c
  end interface verticalInt

  public :: VGridFieldDataUtil_Init, VGridFieldDataUtil_Final
  public :: r_z, z_r
  public :: delz
  public :: verticalInt

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VGridFieldDataUtil_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine VGridFieldDataUtil_Init()

    ! 実行文; Executable statements
    !

  end subroutine VGridFieldDataUtil_Init

  !>
  !!
  !!
  subroutine VGridFieldDataUtil_Final()

    ! 実行文; Executable statements
    !

  end subroutine VGridFieldDataUtil_Final

  !> @brief 
  !!
  !!
  function re_ze(ze_field) result(re_field)
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: ze_field
    type(surfaceScalarField) :: re_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: e
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(re_field, plMesh, "re_field", vLayerNum=nVrLyr)
    re_field%TempDataFlag = .true.

    !$omp parallel do
    do e=1, nEdge
       re_field%data%v_(1:nVrLyr,e) = 0.5d0*( ze_field%data%v_(0:nVzLyr,e) + ze_field%data%v_(1:nVzLyr+1,e) )
    end do

    if(ze_field%TempDataFlag) call Release(ze_field)

  end function re_ze


  !> @brief 
  !!
  !!
  function rc_zc(zc_field) result(rc_field)
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(in) :: zc_field
    type(volScalarField) :: rc_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(rc_field, plMesh, "rc_field", vLayerNum=nVrLyr)
    rc_field%TempDataFlag = .true.

    !$omp parallel do
    do i=1, nCell
       rc_field%data%v_(1:nVrLyr,i) = 0.5d0*( zc_field%data%v_(0:nVzLyr,i) + zc_field%data%v_(1:nVzLyr+1,i) )
    end do

    if(zc_field%TempDataFlag) call Release(zc_field)

  end function rc_zc

  !> @brief 
  !!
  !!
  function ze_re(re_field) result(ze_field)
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: re_field
    type(surfaceScalarField) :: ze_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: e
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(ze_field, plMesh, "ze_field", vLayerNum=nVzLyr)
    ze_field%TempDataFlag = .true.

    !$omp parallel do
    do e=1, nEdge
       ze_field%data%v_(1:nVzLyr,e) = 0.5d0*( re_field%data%v_(1:nVrLyr-1,e) + re_field%data%v_(2:nVrLyr,e) )
    end do

    if(re_field%TempDataFlag) call Release(re_field)

  end function ze_re

  !> @brief 
  !!
  !!
  function zc_rc(rc_field) result(zc_field)
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(in) :: rc_field
    type(volScalarField) :: zc_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(zc_field, plMesh, "zc_field", vLayerNum=nVzLyr)
    zc_field%TempDataFlag = .true.

    !$omp parallel do
    do i=1, nCell
       zc_field%data%v_(1:nVzLyr,i) = 0.5d0*( rc_field%data%v_(1:nVrLyr-1,i) + rc_field%data%v_(2:nVrLyr,i) )
    end do

    if(rc_field%TempDataFlag) call Release(rc_field)

  end function zc_rc

  !> @brief 
  !!
  !!
  function delz_zc(zc_field, zc_lyrThick) result(rc_field)
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(in) :: zc_field
    type(volScalarField), intent(in) :: zc_lyrThick
    type(volScalarField) :: rc_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: i
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(rc_field, plMesh, "rc_field", vLayerNum=nVrLyr)
    rc_field%TempDataFlag = .true.

    !$omp parallel do
    do i=1, nCell
       rc_field%data%v_(1:nVrLyr,i) = &
            & ( zc_field%data%v_(0:nVzLyr,i) - zc_field%data%v_(1:nVzLyr+1,i) )  &
            & / ( 0.5d0*(zc_lyrThick%data%v_(1:nVzLyr+1,i) + zc_lyrThick%data%v_(0:nVzLyr,i)) )
    end do

    if(zc_field%TempDataFlag) call Release(zc_field)

  end function delz_zc

  !> @brief 
  !!
  !!
  function delz_ze(ze_field, zc_lyrThick) result(re_field)
    
    ! モジュール引用; Use statements
    !
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: ze_field
    type(volScalarField), intent(in) :: zc_lyrThick
    type(surfaceScalarField) :: re_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: e, cellIds(2)
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(re_field, plMesh, "rc_field", vLayerNum=nVrLyr)
    re_field%TempDataFlag = .true.

    !$omp parallel do private(cellIds)
    do e=1, nEdge
       cellIds(:) = fvmInfo%Face_CellId(1:2,e)
       re_field%data%v_(1:nVrLyr,e) = & 
            & ( ze_field%data%v_(0:nVzLyr,e) - ze_field%data%v_(1:nVzLyr+1,e) ) &
            & / sum( 0.5d0*(zc_lyrThick%data%v_(1:nVzLyr+1, cellIds) + zc_lyrThick%data%v_(0:nVzLyr, cellIds)), 2)
    end do

    if(ze_field%TempDataFlag) call Release(ze_field)

  end function delz_ze

  !> @brief 
  !!
  !! @return 
  !!
  function verticalInt_e(ze_field, ze_lyrThick, avgFlag) result(e_field)
    
    ! 宣言文; Declaration statement
    !
    type(surfaceScalarField), intent(in) :: ze_field
    type(surfaceScalarField), intent(in) :: ze_lyrThick
    logical, optional, intent(in) :: avgFlag
    type(surfaceScalarField) :: e_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: e
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(e_field, plMesh, "e_field", vLayerNum=1) 
    e_field%TempDataFlag = .true.

    if( present(avgFlag) .and. avgFlag ) then
       !$omp parallel do
       do e=1, nEdge
          e_field%data%v_(1,e) = sum( ze_field%data%v_(1:nVzLyr,e)*ze_lyrThick%data%v_(1:nVzLyr,e) ) &
               & / sum( ze_lyrThick%data%v_(1:nVzLyr,e) )
       end do
    else
       !$omp parallel do
       do e=1, nEdge
          e_field%data%v_(1,e) = sum( ze_field%data%v_(1:nVzLyr,e)*ze_lyrThick%data%v_(1:nVzLyr,e) )
       end do
    end if

    if( ze_field%TempDataFlag ) call Release(ze_field)
    if( ze_lyrThick%TempDataFlag ) call Release(ze_lyrThick)

  end function verticalInt_e

  !> @brief 
  !!
  !! @return 
  !!
  function verticalInt_c(zc_field, zc_lyrThick, avgFlag) result(c_field)
    
    ! 宣言文; Declaration statement
    !
    type(volScalarField), intent(in) :: zc_field
    type(volScalarField), intent(in) :: zc_lyrThick
    logical, optional, intent(in) :: avgFlag
    type(volScalarField) :: c_field
    
    ! 局所変数
    ! Local variables
    !
    integer :: c
    
    ! 実行文; Executable statement
    !

    call GeometricField_Init(c_field, plMesh, "e_field", vLayerNum=1) 
    c_field%TempDataFlag = .true.

    if( present(avgFlag) .and. avgFlag ) then
       !$omp parallel do
       do c=1, nCell
          c_field%data%v_(1,c) = sum( zc_field%data%v_(1:nVzLyr,c)*zc_lyrThick%data%v_(1:nVzLyr,c) ) &
               & / sum( zc_lyrThick%data%v_(1:nVzLyr,c) )
       end do
    else
       !$omp parallel do
       do c=1, nEdge
          c_field%data%v_(1,c) = sum( zc_field%data%v_(1:nVzLyr,c)*zc_lyrThick%data%v_(1:nVzLyr,c) )
       end do
    end if

    if( zc_field%TempDataFlag ) call Release(zc_field)
    if( zc_lyrThick%TempDataFlag ) call Release(zc_lyrThick)

  end function verticalInt_c



end module VGridFieldDataUtil_mod

