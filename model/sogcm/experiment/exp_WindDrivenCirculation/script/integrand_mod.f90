module integrand_mod
  
  use dc_types, only :DP
  
  implicit none
  private

  public :: integrand_set
  public :: gauss_integrand
  public :: W_UpperEBL_integrand
  public :: U_UpperEBL_integrand
  public :: V_UpperEBL_integrand

  real(DP), parameter :: PI = acos(-1d0)
  real(DP) :: param_y, param_z, param_p
  real(DP), pointer :: p_integGauss(:)

contains
  subroutine integrand_set(param_y_, param_z_, param_p_, p_integ_gauss_)
    real(DP), intent(in) :: param_y_
    real(DP), intent(in) :: param_z_
    real(DP), intent(in) :: param_p_
    real(DP), intent(in), target :: p_integ_gauss_(:)
    
    param_y = param_y_
    param_z = param_z_
    param_p = param_p_
    p_integGauss => p_integ_gauss_

  end subroutine integrand_set

  real(DP) function gauss_integrand(z, idx)
    real(DP), intent(in) :: z
    integer, intent(in) :: idx

    gauss_integrand = exp(- 0.25d0*z**2/param_p)

  end function gauss_integrand

  real(DP) function W_UpperEBL_integrand(p, idx)
    real(DP), intent(in) :: p
    integer, intent(in) :: idx
    
    real(DP) :: integ_Gauss

    W_UpperEBL_integrand = &
         & 1d0/sqrt(PI) *sqrt(p)*exp(- p**3/3d0)*cos(param_y*p) &
         &   * p_integGauss(idx)

  end function W_UpperEBL_integrand

  
  
   real(DP) function U_UpperEBL_integrand(p, idx)
    real(DP), intent(in) :: p
    integer, intent(in) :: idx

    real(DP) :: p2

    p2 = p**2
    U_UpperEBL_integrand = - 2d0/sqrt(PI) * exp(-p2**3/3d0 - 0.25d0*param_z**2/p2) * cos(p2*param_y) 

  end function U_UpperEBL_integrand

   real(DP) function V_UpperEBL_integrand(p, idx)
    real(DP), intent(in) :: p
    integer, intent(in) :: idx

    real(DP) :: p2

    p2 = p**2
    V_UpperEBL_integrand = 2d0/sqrt(PI) * exp(-p2**3/3d0 - 0.25d0*param_z**2/p2) * sin(p2*param_y) 

  end function V_UpperEBL_integrand
end module integrand_mod

!---------------------------------
