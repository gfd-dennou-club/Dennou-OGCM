#define POLYORDER_Nk3

module SimParameters_mod

  use dc_types

  implicit none

  real(DP), parameter :: radius  = 6.37122d06
  real(DP), parameter :: Omega = 7.292d-05
  real(DP), parameter :: Grav = 9.80616d0
  real(DP), parameter :: meanDepth = 1000d0!2.94d04 / Grav
  real(DP), parameter :: refDens = 1000d0
  integer, parameter :: delTime       = 90d0
  integer, parameter :: endTime  = 100*86400
  integer, parameter :: endTStep = endTime / delTime 
  integer, parameter :: outputIntrVal = 43200
  real(DP), parameter :: PI = acos(-1d0)
  real(DP), parameter :: LinearDragCoef = 3d-6

  character(STRING) :: gridFilePath = & !"/home/ykawai/workspace/Dennou-OGCM/tool/glevel4_globe.nc"
                & "/home/ykawai/workspace/Dennou-OGCM/tool/glevel4_shoreline3.nc"!grid/grid-glevel5_itr.nc"
  character(STRING) :: gridUsageFilePath = & !""
                & "/home/ykawai/workspace/Dennou-OGCM/tool/glevel4_shoreline3_MeshUsage.nc"

#ifdef POLYORDER_Nk1
  integer, parameter :: PolyDegree  = 1
#elif defined POLYORDER_Nk2
  integer, parameter :: PolyDegree  = 2
#elif defined POLYORDER_Nk3
  integer, parameter :: PolyDegree  = 3
#endif

  real(DP), parameter :: latlonIntrv = 1.0d0* PI/180d0
  

!!$contains

end module SimParameters_mod
