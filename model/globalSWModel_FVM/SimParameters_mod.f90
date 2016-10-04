module SimParameters_mod

  use dc_types

  implicit none

  real(DP), parameter :: radius  = 6.37122d06
  real(DP), parameter :: Omega = 7.292d-05
  real(DP), parameter :: Grav = 9.80616
  real(DP), parameter :: meanDepth = 2000d0

  integer, parameter :: delTime       = 900
  integer, parameter :: endTime  = 2*24*3600
  integer, parameter :: endTStep = endTime / delTime 
  integer, parameter :: outputIntrVal = 2*24*3600
  real(DP), parameter :: PI = acos(-1d0)

  character(STRING) :: gridFilePath = "../../tool/grid/grid-glevel5.nc"

!!$contains
end module SimParameters_mod
