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
  public :: Diagnose_SigDot, Diagnose_PressBaroc, Diagnose_GeoPot
  public :: Diagnose_DensEdd

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

    xyz_UrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Urf)
    xyz_VrfHatSig = xyz_IntSig_SigToTop_xyz(xyz_Vrf)
    xyz_DivHatSig = xyz_IntSig_SigToTop_xyz(xyz_Div)
    xy_UrfHat = xy_IntSig_BtmToTop_xyz(xyz_Urf)
    xy_VrfHat = xy_IntSig_BtmToTop_xyz(xyz_Vrf)
    xy_DivHat = xy_IntSig_BtmToTop_xyz(xyz_Div)

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
    
!    xyz_SigDot(:,:,kMax) = 0d0

  end function diagnose_SigDot

  function diagnose_PressBaroc(xy_totDepth, xyz_DensEdd) result(xyz_PressBaroc)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PressBaroc(0:iMax-1,jMax,0:kMax)
    xyz_PressBaroc =  Grav*spread(xy_totDepth,3,kMax+1)*(xyz_IntSig_SigToTop_xyz(xyz_DensEdd)) 

  end function diagnose_PressBaroc

  function diagnose_GeoPot( xy_totDepth, xy_SurfHeight ) result(xyz_GeoPot)

    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SurfHeight(0:iMax-1,jMax)
    real(DP) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)

    integer :: k

    !$omp parallel do
    do k=1, kMax
       xyz_GeoPot(:,:,k) = g_Sig(k)*xy_totDepth + xy_SurfHeight
    end do

  end function diagnose_GeoPot

  function diagnose_DensEdd( xyz_PTempEdd, xyz_GeoPot, z_PTempBasic, refDens ) result(xyz_DensEdd)

    real(DP), intent(in) :: xyz_PTempEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: z_PTempBasic(0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: refDens

    real(DP), parameter :: Cp = 3986d0
    real(DP), parameter :: BetaT = 1.67d-04
    real(DP), parameter :: T0 = 283d0
    real(DP) :: H_T
    real(DP) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)

    xyz_PTemp = xyz_PTempEdd + spread(spread(z_PTempBasic,1,jMax), 1, iMax)
    H_T = Cp/(BetaT*Grav)!Grav/(BetaT*Cp)
write(*,*) "H_T=",H_T
    xyz_DensEdd = - refDens*( BetaT*(xyz_PTemp*exp( (xyz_GeoPot/Grav)/H_T  ) - T0 ) )

  end function diagnose_DensEdd

end module DiagnoseUtil_mod
