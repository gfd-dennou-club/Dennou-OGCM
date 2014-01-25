!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DiagVarEval_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, STRING, TOKEN

  use dc_message, only: &
       & MessageNotify
  

  use Constants_mod, only: &
       & RPlanet, Grav, RefDens

  use GridSet_mod, only: &
       & iMax, jMax, kMax, &
       & xyz_Lat, xyz_Lon

  use SpmlUtil_mod

  use DiagnoseUtil_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DiagVarEval_Init, DiagVarEval_Final
  public :: eval_Vor, eval_Div, eval_totPress, eval_DensEdd
  public :: eval_potentialEnergyAvg, eval_kineticEnergyAvg

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DiagVarEval_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DiagVarEval_Init()

    ! 実行文; Executable statements
    !

    call DiagnoseUtil_Init()

  end subroutine DiagVarEval_Init

  !>
  !!
  !!
  subroutine DiagVarEval_Final()

    ! 実行文; Executable statements
    !

    call DiagnoseUtil_Final()

  end subroutine DiagVarEval_Final

  !> @brief 
  !!
  !! @return 
  !!
  function eval_Div(xyz_u, xyz_v) result(xyz_Div)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_v(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Div(0:iMax-1,jMax,0:kMax)

    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    xyz_Div = xyz_wz( &
         & wz_AlphaOptr_xyz(xyz_u*cos(xyz_Lat), xyz_v*cos(xyz_Lat)) &
         & )

  end function eval_Div

  !> @brief 
  !!
  !! @return 
  !!
  function eval_Vor(xyz_u, xyz_v) result(xyz_Vor)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_v(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Vor(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !
    xyz_Vor = xyz_wz( &
         & wz_AlphaOptr_xyz(xyz_v*cos(xyz_Lat), -xyz_u*cos(xyz_Lat)) &
         & )
    
  end function eval_Vor

  !> @brief 
  !!
  !! @return 
  !!
  function eval_DensEdd(xyz_PTemp, xyz_Salt, xyz_totPress) result(xyz_DensEdd)
    
    use EOSDriver_mod, only: EOSDriver_Eval

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_totPress(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    call EOSDriver_Eval( rhoEdd=xyz_DensEdd, &           !(out)
         & theta=xyz_PTemp, S=xyz_Salt, p=xyz_totPress )  !(in)

  end function eval_DensEdd

  !> @brief 
  !!
  !! @return 
  !!
  function eval_totPress(xy_surfPress, xyz_barocPress) result(xyz_totPress)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xy_surfPress(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_barocPress(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_totPress(0:iMax-1,jMax,0:kMax)

    !
    !
    !
    integer :: k

    ! 実行文; Executable statement
    !

    do k=0, kMax
       xyz_totPress(:,:,k) = xy_surfPress(:,:) + xyz_barocPress(:,:,k)
    end do
  end function eval_totPress


  !> @brief 
  !!
  !! @return 
  !!
  function eval_kineticEnergyAvg(xyz_u, xyz_v) result(KEAvg)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_v(0:iMax-1,jMax,0:kMax)
    real(DP) :: KEAvg
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    KEAvg = RefDens*AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyz_u**2 + xyz_v**2) )

  end function eval_kineticEnergyAvg

  !> @brief 
  !!
  !! @return 
  !!
  function eval_potentialEnergyAvg(xyz_DensEdd, xy_totDepth) result(PEAvg)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: PEAvg
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !
    
    PEAvg = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
                      xyz_DensEdd * Diagnose_GeoPot(xy_totDepth) &
            & ))


  end function eval_potentialEnergyAvg

end module DiagVarEval_mod
