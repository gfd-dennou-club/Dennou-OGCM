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
  interface Diagnose_SigDot
!     module procedure Diagnose_SigDot_new
     module procedure Diagnose_SigDot_old
  end interface Diagnose_SigDot

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


  function diagnose_SigDot_new( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div ) result(xyz_SigDot)

    use SpmlUtil_mod, only: &
         tr_vDeriv1CoefMat
    use at_module_omp, only: &
         g_t, t_Dx_t, t_g
    
    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP), intent(in) :: xyz_Urf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vrf(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_SigDot(0:iMax-1,jMax,0:kMax)

    real(DP) :: DSigMat(0:kMax,0:kMax)
    real(DP) :: DSigMatSave(0:kMax,0:kMax)
    real(DP) :: IMat(0:kMax,0:kMax)
    integer :: i, j, k, n 
    integer :: IPiv(0:kMax)
    real(DP) :: b(0:kMax,1)
    integer :: info

    IMat(:,:) = 0d0
    forAll(k=0:kMax) IMat(k,k) = 1d0
    do k=0, kMax
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do
!    DSigMat(0,:) = IMat(0,:)
    DSigMat(kMax,:) = IMat(kMax,:)
!    DSigMat = transpose(DSigMat)
    DSigMatSave = DSigMat
    n = size(IMat,1)
    call DGETRF(n, n, DSigMat, n, IPiv, info)

    !$omp parallel do private(i, b, info)
    do j=1, jMax
       do i=0, iMax-1
          b(:,1) = - xyz_Div(i,j,:)
 !         b(0,1) = 0d0;
          b(kMax,1) = 0d0
          
          call DGETRS('N', n, 1, DSigMat, n, IPiv, b, n, info)
!!$          do k=0,kMax
!!$             write(*,*) DSigMat(:,k+1)
!!$          end do
!!$          write(*,*) "-----"
!!$          call DGESV(n, 1, DSigMat, n, IPiv, b, n, info)
          xyz_SigDot(i,j,:) = b(:,1)
!!$          write(*,*) "DGETRS info=", info

!!$          write(*,*) "j=", j , matmul(DSigMatSave,b(:,1))
!!$          write(*,*) "j=", j, z_DSig_z(b(:,1))
!!$          write(*,*) "Div=", -xyz_Div(i,j,:)
!!$          write(*,*) "b=", b(:,1)
!!$          write(*,*) "SigDot=", xyz_SigDot(i,j,:)
!!$          write(*,*) "IPiv=", IPiv
!!$          stop
       end do
    end do
!!$stop
  end function diagnose_SigDot_new
  
  function diagnose_SigDot_old( xy_totDepth, xyz_Urf, xyz_Vrf, xyz_Div ) result(xyz_SigDot)


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

!!$real(DP) :: xyz_dwdz(0:iMax-1,jMax,0:kMax)
!!$real(DP) :: xyz_Tmp(0:iMax-1,jMax,0:kMax)
!!$real(DP), parameter  :: m = 20d0

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

!!$do k=0,kMax
!!$   xyz_Tmp(:,:,k) = cos(2d0*PI*m*g_Sig(k))
!!$end do
!!$xyz_dwdz = xyz_DSig_xyz(xyz_Tmp)
!!$write(*,*) "num=", xyz_dwdz(0,16,:) / (2d0*PI*m)
!!$write(*,*) "exact=", -sin(2d0*PI*m*g_Sig(:)) 
!!$
!!$xyz_dwdz = xyz_IntSig_SigToTop_xyz(xyz_Tmp)
!!$write(*,*) "num=", xyz_dwdz(0,16,:) * (2d0*PI*m)
!!$write(*,*) "exact=", -sin(g_Sig(:)*2d0*PI*m)
!!$stop
!!$xyz_dwdz = xyz_wt(wt_DSig_wt(wt_xyz(xyz_SigDot)))
!!$xyz_dwdz = xyz_DSig_xyz(xyz_SigDot)
!!$write(*,*) "---- diagnose sigdot ---"
!!$write(*,*) "HDiv=", xyz_Div(0,1:16,0)
!!$write(*,*) "DsigSigDot=", xyz_dwdz(0,1:16,0)
!!$write(*,*) "DivHat=", xy_DivHat(0,1:16)
!!$write(*,*) "-------"
!!$write(*,*) "HDiv=", xyz_Div(0,16,:)
!!$write(*,*) "DsigSigDot=", xyz_dwdz(0,16,:)
!!$write(*,*) "-------****************************"
!!$
!!$write(*,*) (xyz_Div(1,1:16,0) + xyz_dwdz(1,1:16,0))!/maxval(xyz_Div)    
!    xyz_SigDot(:,:,kMax) = 0d0

  end function diagnose_SigDot_old

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
