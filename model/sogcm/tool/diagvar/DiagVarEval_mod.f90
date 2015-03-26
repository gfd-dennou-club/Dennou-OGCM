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
       & RPlanet, Grav, RefDens, &
       & hDiffCoef, vDiffCoef

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

  !
  public :: eval_Vor, eval_Div, eval_StreamPot
  public :: eval_MassStreamFunc
  public :: eval_PressEdd, eval_DensEdd, eval_DensPot, eval_Temp
  public :: eval_StaticStability

  !
  public :: eval_potentialEnergyAvg, eval_kineticEnergyAvg
  public :: eval_angularMomAvg


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
  subroutine eval_StreamPot(xyz_Psi, xyz_Chi, & ! (out)
       & xyz_u, xyz_v) !(in)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), optional :: xyz_Psi(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out), optional :: xyz_Chi(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_v(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !
    
    if (present(xyz_Psi)) then
       xyz_Psi = xyz_wz( wz_InvLapla2D_wz(wz_xyz(eval_Vor(xyz_u, xyz_v))) )
    end if

    if (present(xyz_Chi)) then
       xyz_Chi = xyz_wz( wz_InvLapla2D_wz(wz_xyz(eval_Div(xyz_u, xyz_v))) )
    end if


  end subroutine eval_StreamPot

  !> @brief 
  !!
  !! @return 
  !!
  function eval_DensEdd(xyz_PTemp, xyz_Salt, xy_totDepth) result(xyz_DensEdd)
    
    use EOSDriver_mod, only: EOSDriver_Eval

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    call EOSDriver_Eval( rhoEdd=xyz_DensEdd, &           !(out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*Diagnose_GeoPot(xy_totDepth) )  !(in)

  end function eval_DensEdd

  !> @brief 
  !!
  !! 
  !! @return 
  !!
  function eval_DensPot(xyz_PTemp, xyz_Salt, PressRef) result(xyz_DensEdd)
    
    use EOSDriver_mod, only: EOSDriver_Eval

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in), optional :: PressRef          
    real(DP) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: PressRef_
    real(DP) :: xyz_Press(0:iMax-1,jMax,0:kMax)
    
    ! 実行文; Executable statement
    !

    PressRef_ = 0d0
    if(present(PressRef)) PressRef_ = PressRef

    xyz_Press = PressRef_
    call EOSDriver_Eval( rhoEdd=xyz_DensEdd, &           !(out)
         & theta=xyz_PTemp, S=xyz_Salt, p=xyz_Press )  !(in)

  end function eval_DensPot

  !> @brief 
  !!
  !! @return 
  !!
  function eval_Temp(xyz_PTemp, xyz_Salt, xy_totDepth) result(xyz_Temp)
    
    use EOSDriver_mod, only: EOSDriver_PTemp2Temp

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Salt(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1,jMax)
    real(DP) :: xyz_Temp(0:iMax-1,jMax,0:kMax)

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    call EOSDriver_PTemp2Temp( InSituTemp=xyz_Temp, &           !(out)
         & theta=xyz_PTemp, S=xyz_Salt, p=-RefDens*Diagnose_GeoPot(xy_totDepth) )  !(in)

  end function eval_Temp

  !> @brief 
  !!
  !! @return 
  !!
  function eval_PressEdd(xy_surfPress, xyz_HydroPressEdd) result(xyz_PressEdd)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xy_surfPress(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_HydroPressEdd(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_PressEdd(0:iMax-1,jMax,0:kMax)

    !
    !
    !
    integer :: k

    ! 実行文; Executable statement
    !

    do k=0, kMax
       xyz_PressEdd(:,:,k) = xy_surfPress(:,:) + xyz_HydroPressEdd(:,:,k)
    end do

  end function eval_PressEdd

  !> @brief Calculate mass streamfunction from meridional velocity.
  !! 
  !! In this module, Sv is used as the unit for mass streamfunction. 
  !! The 1 Sv is equivalent to rho0 * 10^6 m3.s-1 = 10^9 kg.s-1.
  !! (We follow the definition of Sv in Marshall et al.(2007).)
  !!
  function eval_MassStreamFunc(xyz_V, xy_totDepth) result(yz_MassStreamFunc)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP) :: yz_MassStreamFunc(jMax,0:kMax)
    
    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: OneSvUnit = 1d09
    
    ! 実行文; Executable statement
    !
    
    yz_MassStreamFunc = RPlanet*ya_IntLon_xya( &
         & xyz_IntSig_SigToTop_xyz( &
         & - RefDens*xyz_V*spread(cos(xyz_Lat(:,:,0))*xy_totDepth, 3, kMax+1) ) &
         & ) / OneSvUnit

  end function eval_MassStreamFunc

  !> @brief 
  !!
  !! @return 
  !!
  function eval_StaticStability(xyz_PTemp, xy_totDepth) result(xyz_Stability)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_PTemp(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_totDepth(0:iMax-1, jMax)
    real(DP) :: xyz_Stability(0:iMax-1,jMax,0:kMax)
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    xyz_Stability = Grav/xyz_PTemp * xyz_xyt( & 
         & xyt_DSig_xyt(xyt_xyz(xyz_PTemp/spread(xy_totDepth, 3, kMax+1))) & 
         & )


  end function eval_StaticStability



!!!!!!!!!!!!!! Subroutines for Energy Calculation

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
    integer :: k
    
    ! 実行文; Executable statement
    !
    
    KEAvg = 0.5d0*AvrLonLat_xy( xy_IntSig_BtmToTop_xyz(xyz_u**2 + xyz_v**2) )

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
    real(DP) :: temp(0:iMax-1,jMax,0:kMax)
    ! 実行文; Executable statement
    !

    PEAvg = AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
                    xyz_DensEdd*Diagnose_GeoPot(xy_totDepth) &
            & ))

  end function eval_potentialEnergyAvg


  !> @brief 
  !!
  !! @return 
  !!
  function eval_angularMomAvg(xyz_u) result(angularMomAvg)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: xyz_u(0:iMax-1,jMax,0:kMax)
    real(DP) :: angularMomAvg
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    angularMomAvg =  AvrLonLat_xy( xy_IntSig_BtmToTop_xyz( &
         & cos(xyz_Lat)*xyz_u &
         & ) )

  end function eval_angularMomAvg


end module DiagVarEval_mod
