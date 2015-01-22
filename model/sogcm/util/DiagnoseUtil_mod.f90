!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DiagnoseUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING

  use SpmlUtil_mod

  use Constants_mod
  use GridSet_mod, only: &
       & iMax, jMax, kMax
  

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DiagnoseUtil_Init, DiagnoseUtil_Final
  public :: Diagnose_SigDot, Diagnose_HydroPressEdd, Diagnose_GeoPot


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DiagnoseUtil_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DiagnoseUtil_Init()

    ! 実行文; Executable statements
    !

  end subroutine DiagnoseUtil_Init

  !>
  !!
  !!
  subroutine DiagnoseUtil_Final()

    ! 実行文; Executable statements
    !

  end subroutine DiagnoseUtil_Final


  function diagnose_SigDot( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div ) result(xyz_SigDot)


    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    integer :: k
    real(DP) :: xyz_UrfHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_VrfHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DivHatSig(0:iMax-1,jMax,0:kMax)
    real(DP) :: xy_UrfHat(0:iMax-1,jMax)
    real(DP) :: xy_VrfHat(0:iMax-1,jMax)
    real(DP) :: xy_DivHat(0:iMax-1,jMax)
    real(DP) :: sigWeight
    real(DP) :: xy_DtotDepthDmu(0:iMax-1,jMax)
    real(DP) :: xy_DtotDepthDLambda(0:iMax-1,jMax)
    real(DP) :: sig

 real(DP) :: xyz_dwdz(0:iMax-1,jMax,0:kMax)

!!$    xyz_UrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Urf)
!!$    xyz_VrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Vrf)
    xyz_DivHatSig = xyz_IntSig_SigToTop_xyz(xyz_Div)
!!$    xy_UrfHat = xy_IntSig_BtmToTop_xyz(xyz_Urf)
!!$    xy_VrfHat = xy_IntSig_BtmToTop_xyz(xyz_Vrf)
!!$    xy_DivHat = xy_IntSig_BtmToTop_xyz(xyz_Div)

    xy_DtotDepthDLambda = 0d0!xy_w(w_DivLambda_xy(xy_totDepth))
    xy_DtotDepthDmu =0d0! xy_w(w_DivMu_xy(xy_totDepth))
    
    !$omp parallel do private(sig)
    do k=0, kMax
       sig = g_Sig(k)
       xyz_SigDot(:,:,k) = &
!!$            &   sig*( &
!!$            &        (xy_UrfHat*xy_DtotDepthDLambda + xy_VrfHat*xy_DtotDepthDmu)/xy_totDepth &
!!$            &     +  xy_DivHat )  &
            & + xyz_DivHatSig(:,:,k) !&
!!$            & + (xyz_UrfHatSig(:,:,k)*xy_DtotDepthDLambda + xyz_VrfHatSig(:,:,k)*xy_DtotDepthDmu)/xy_totDepth 
    end do

!!$xyz_dwdz = xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot)))
!!$write(*,*) "Div----"
!!$write(*,*) xyz_Div(1,16,:)
!!$write(*,*) "---- diagnose sigdot ---"
!!$write(*,*) (xyz_Div(1,1:16,0) + xyz_dwdz(1,1:16,0))!/maxval(xyz_Div)    
!    xyz_SigDot(:,:,kMax) = 0d0

  end function diagnose_SigDot

  function diagnose_HydroPressEdd(xy_totDepth, xyz_DensEdd) result(xyz_HydroPressEdd)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_HydroPressEdd(0:iMax-1,jMax,0:kMax)

    xyz_HydroPressEdd(:,:,:) =  Grav*spread(xy_totDepth,3,kMax+1)*(xyz_IntSig_SigToTop_xyz(xyz_DensEdd)) 

  end function diagnose_HydroPressEdd

  function diagnose_GeoPot( xy_totDepth ) result(xyz_GeoPot)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)

    integer :: k

    !$omp parallel do
    do k=0, kMax
       xyz_GeoPot(:,:,k) = Grav*( g_Sig(k)*xy_totDepth )
    end do

  end function diagnose_GeoPot

end module DiagnoseUtil_mod
