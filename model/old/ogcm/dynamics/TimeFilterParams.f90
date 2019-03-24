program TimeFilterParams

  use dc_types
  use dc_message

  implicit none

  real(DP) :: A0, tau0
  real(DP), parameter :: EPS = 1d-14
  integer, parameter :: itrMax = 10
  integer, parameter :: p = 2
  integer, parameter :: q = 4
  real(DP), parameter :: r = 0.28462d0

  integer :: itr
  real(DP) :: tauMax

  tau0 = (p+2)*(p+q+2)/dble((p+1)*(p+q+1))
  tauMax = tau0
  do itr=1, itrMax
     tauMax = calcTauMax(tau0, tauMax)
     tau0 = calcTau0(tauMax, tau0)
  end do

  call MessageNotify("M", "TimeFilterParams", &
       & "p=%d; q=%d; r=%f", i=(/ p, q /), d=(/ r /))

  call MessageNotify("M", "TimeFilterParams", &
       & "A0=%f; Tau0=%f; TauMax=%f", d=(/ A0, tau0, tauMax /))
contains
function A(tau, tau0) result(val)
  real(DP), intent(in) :: tau, tau0
  real(DP) :: val

  real(DP) :: x

  x = tau/tau0
  val = x**p * ( 1d0 - x**q) - r*x
  
end function A

function dAdTau(tau, tau0) result(val)
  real(DP), intent(in) :: tau, tau0
  real(DP) :: val

  real(DP) :: x

  x = tau/tau0
  val = ( p*x**(p-1) - (p+q)*x*(p+q-1) - r )/tau0
  
end function dAdtau

function calcTau0(tauMax, guess) result(tau0)
  real(DP), intent(in) :: tauMax, guess
  real(DP) :: tau0

  real(DP) :: dTau0, res
  integer :: nItr
  real(DP) :: func, dfuncdtau0
  real(DP) :: x, J0, J1, dJ0dTau0, dJ1dTau0

  tau0 = guess
  nItr = 0
  res = 1d0

  do while(.true.)

     if( abs(res) < EPS ) exit

     x = tauMax/tau0

     J0 = x**(p+1)/dble(p+1) - x**(p+q+1)/dble(p+q+1) - r*x**2/2d0
     dJ0dTau0 = 1d0/tauMax *( -x**(p+2) + x**(p+q+2) + r*x**3 )
     J1 = tau0*( x**(p+2)/dble(p+2) - x**(p+q+2)/dble(p+q+2) - r*x**3/3d0 )
     dJ1dTau0 = J1/tau0 + 1d0/x *( -x**(p+3) + x**(p+q+3) + r*x**4 )
     res = J0 - J1
     dTau0 = - res/( dJ0dTau0 - dJ1dTau0 )
     tau0 = tau0 + dTau0

     nItr = nItr + 1
  end do

  call MessageNotify("M","calcTau0", &
       & "Itr=%d, Residual=%f, tau0=%f, J0=%f", i=(/ nItr /), d=(/ res, tau0, J0*Tau0 /))
  A0 = 1d0/(J0*Tau0)

end function calcTau0

function calcTauMax(tau0, guess) result(tauMax)

  real(DP), intent(in) :: tau0, guess
  real(DP) :: tauMax

  real(DP) :: dTauMax, res
  integer :: nItr

  tauMax = guess
  nItr = 0
  res = 1d0

  do while(.true.)

     if( abs(res) < EPS ) exit
     res = A(tauMax, tau0)
     dTauMax = - res / dAdTau(tauMax, tau0)
     tauMax = tauMax + dTauMax


     nItr = nItr + 1
  end do

  call MessageNotify("M","calcTauMax", &
       & "Itr=%d, Residual=%f, tauMax=%f", i=(/ nItr /), d=(/ res, tauMax /))

end function calcTauMax


end program TimeFilterParams
