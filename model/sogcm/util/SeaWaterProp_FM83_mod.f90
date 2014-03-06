!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module SeaWaterProp_FM83_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SeaWaterProp_FM83_Init, SeaWaterProp_FM83_Final
  public :: SeaWaterProp_FM83_AdLapseRate, SeaWaterProp_FM83_Temp2PTemp, &
       & SeaWaterProp_FM83_HeatCapacity

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SeaWaterProp_FM83_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine SeaWaterProp_FM83_Init()

    ! 実行文; Executable statements
    !

  end subroutine SeaWaterProp_FM83_Init

  !>
  !!
  !!
  subroutine SeaWaterProp_FM83_Final()

    ! 実行文; Executable statements
    !

  end subroutine SeaWaterProp_FM83_Final

  !> @brief Calculate the adibatic lapse rate from practical salinity, in-situ temperature and pressure. 
  !!
  !! @param [in] S salinity [psu]
  !! @param [in] T in-situ temperature [degree C]
  !! @param [in] p pressure [dbar]
  !! @return adiabatic lapse rate[degree C/dbar]
  !!
  function SeaWaterProp_FM83_AdLapseRate(S, t, p) result(adLapseRate)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: S
    real(DP), intent(in) :: t
    real(DP), intent(in) :: p
    real(DP) :: adLapseRate

    ! 局所変数
    ! Local variables
    !
    real(DP), parameter :: a0 =  3.5803d-05, &
         &                 a1 =  8.5258d-06, &
         &                 a2 = -6.8360d-08, &
         &                 a3 =  6.6228d-10, &
         &                 b0 =  1.8932d-06, &
         &                 b1 = -4.2393d-08, &
         &                 c0 =  1.8741d-08, &
         &                 c1 = -6.7795d-10, &
         &                 c2 =  8.7330d-12, &
         &                 c3 = -5.4481d-14, &
         &                 d0 = -1.1351d-10, &
         &                 d1 =  2.7759d-12, &
         &                 e0 = -4.6206d-13, &
         &                 e1 =  1.8676d-14, &
         &                 e2 = -2.1687d-16
    
    
    ! 実行文; Executable statement
    !
    
    adLapseRate = &
         &    a0 + t*(a1 + t*(a2 + t*a3)) &
         & + (b0 + t*b1)*(S-35d0) &
         & + p*( (c0 + t*(c1 + t*(c2 + t*c3)) + (d0 + t*d1)*(S-35d0)) &
         &       + p*(e0 + t*(e1 + t*e2)) &
         &   )

  end function SeaWaterProp_FM83_AdLapseRate

  !> @brief 
  !!
  !! @param [in] P0 pressure [dbar]
  !! @param [in] S0 salinity [psu]
  !! @param [in] T0 in-situ temperature [degree C]
  !! @param [in] Pref reference pressure [dbar]
  !! @return potential tempatature [degree C/dbar]
  !!
  function SeaWaterProp_FM83_Temp2PTemp(P0, T0, S0, Pref) result(ptemp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: P0
    real(DP), intent(in) :: T0
    real(DP), intent(in) :: S0
    real(DP), intent(in) :: Pref

    real(DP) :: ptemp

    ! 局所変数
    ! Local variables
    !
    real(DP) :: DelP
    real(DP) :: ThetaTmp, DelTheta
    real(DP) :: qTmp

    ! 実行文; Executable statement
    !
    
    delP = PRef - P0
    ThetaTmp = T0

    ! Step1
    DelTheta = DelP*SeaWaterProp_FM83_AdLapseRate(S=S0, t=ThetaTmp, p=P0)
    ThetaTmp = ThetaTmp + 0.5d0*DelTheta
    qTmp = DelTheta

    ! Step2
    DelTheta = DelP*SeaWaterProp_FM83_AdLapseRate(S=S0, t=ThetaTmp, p=P0+0.5d0*DelP)
    ThetaTmp = ThetaTmp + (1d0 - 1d0/sqrt(2d0))*(DelTheta - qTmp)
    qTmp = (2d0 - sqrt(2d0))*DelTheta + (-2d0 + 3d0/sqrt(2d0))*qTmp

    ! Step3
    DelTheta = DelP*SeaWaterProp_FM83_AdLapseRate(S=S0, t=ThetaTmp, p=P0+0.5d0*DelP)
    ThetaTmp = ThetaTmp + (1d0 + 1d0/sqrt(2d0))*(DelTheta - qTmp)
    qTmp = (2d0 + sqrt(2d0))*DelTheta + (-2d0 - 3d0/sqrt(2d0))*qTmp

    ! Step4
    DelTheta = DelP*SeaWaterProp_FM83_AdLapseRate(S=S0, t=ThetaTmp, p=P0+DelP)
    ThetaTmp = ThetaTmp + 1d0/6d0*(DelTheta - 2d0*qTmp)

    ptemp = ThetaTmp

  end function SeaWaterProp_FM83_Temp2PTemp

  !> @brief Calculate the heat capacity from practical salinity, in-situ temperature and pressure. 
  !!
  !! @param [in] S salinity [psu]
  !! @param [in] T in-situ temperature [degree C]
  !! @param [in] p pressure [dbar]
  !! @return heat capacity [degree C/dbar]
  !!
  function SeaWaterProp_FM83_HeatCapacity(S, t, p) result(cp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: S
    real(DP), intent(in) :: t
    real(DP), intent(in) :: p
    real(DP) :: cp

    ! 局所変数
    ! Local variables
    !
    real(DP) :: p_bar

    ! 実行文; Executable statement
    !
    
    p_bar = p/10d0
    cp = cp_S_t_0() + del1_cp_0_t_p()

    if (S > 0d0) then
       cp  = cp + del2_cp_S_t_p()
    end if
 
   contains
      real(DP) function cp_S_t_0()

        real(DP), parameter :: c0 =   4217.4d00, &
             &                 c1 = - 3.720283d00, &
             &                 c2 =   1.412855d-1, &
             &                 c3 = - 2.654387d-3, &
             &                 c4 =   2.093236d-5, &
             &                 a0 = - 7.643573d00, &
             &                 a1 =   1.072763d-1, &
             &                 a2 = - 1.38385d-3,  &
             &                 b0 =   1.770383d-1, &
             &                 b1 = - 4.07718d-3, &
             &                 b2 =   5.148e-5

        cp_S_t_0 = &
             &   c0 + t*(c1 + t*(c2 + t*(c3 + t*c4))) &  ! cp_0_t_0
             & + S*(a0 + t*(a1 + t*a2) + sqrt(S)*(b0 + t*(b1 + t*b2)))
        
      end function cp_S_t_0

      real(DP) function del1_cp_0_t_p()
        real(DP), parameter :: &
             & a0 = - 4.9592d-1,  &
             & a1 =   1.45747d-2, &
             & a2 = - 3.13885d-4, &
             & a3 =   2.0357d-6, &
             & a4 =   1.7168d-8, &
             & b0 =   2.4931d-4, &
             & b1 = - 1.08645d-5, &
             & b2 =   2.87533d-7, &
             & b3 = - 4.0027d-9, &
             & b4 =   2.2956d-11, &
             & c0 = - 5.422d-8, &
             & c1 =   2.6380d-9, &
             & c2 = - 6.5637d-11, &
             & c3 =   6.136d-13

        del1_cp_0_t_p = p_bar*( &
             &     (a0 + t*(a1 + t*(a2 + t*(a3 + t*a4))))  &
             & + p_bar*( &
             &        b0 + t*(b1 + t*(b2 + t*(b3 + t*b4))) &
             &      + p_bar*(c0 + t*(c1 + t*(c2 + t*c3))) &
             &  )  &
             & )

      end function del1_cp_0_t_p

      real(DP) function del2_cp_S_t_p()
        real(DP), parameter :: &
             &  d0 =   4.9247d-3, &
             &  d1 = - 1.28315d-4, &
             &  d2 =   9.802d-7, &
             &  d3 =   2.5941d-8, &
             &  d4 = - 2.9179d-10, &
             &  e0 =- 1.2331d-4, &
             &  e1 = - 1.517d-6, &
             &  e2 =   3.122d-8, &
             &  f0 = - 2.9558d-6, &
             &  f1 =  1.17054d-7, &
             &  f2 = - 2.3905d-9, &
             &  f3 =   1.8448d-11, &
             &  g0 =   9.971d-8, &
             &  h0 =   5.540d-10, &
             &  h1 = - 1.7682d-11, &
             &  h2 =   3.513d-13, &
             &  J1 = - 1.4300d-12
        real(DP) :: sqrtS

        sqrtS = sqrt(S)
        del2_cp_S_t_p = p_bar*S*( &
             &           (d0 + t*(d1 + t*(d2 + t*(d3 + t*d4))))  &
             &   + sqrtS*(e0 + t*(e1 + t*e2)) &
             &   + p_bar*( &
             &        f0 + t*(f1 + t*(f2 + t*f3)) + sqrtS*g0 &
             &      + p_bar*(h0 + t*(h1 + t*h2 + sqrtS*J1))       & 
             &   ) &
             & )

      end function del2_cp_S_t_p

  end function SeaWaterProp_FM83_HeatCapacity

end module SeaWaterProp_FM83_mod

