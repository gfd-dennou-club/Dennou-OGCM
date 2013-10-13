program SpmlUtil_mod_test
  
  use dc_types
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod

  implicit none

  character(*), parameter :: configNmlFile = "../defaultConfig.nml"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!

  call Constants_Init(configNmlFile)
  call GridSet_Init(configNmlFile)
  call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
  call GridSet_construct()


  call check_xy_IntSig_BtmToTop_xyz
  call check_xyz_IntSig_SigToTop_xyz

  call SpmlUtil_Final()
  call GridSet_Final() 
  call Constants_Final()

contains

  subroutine check_xy_IntSig_BtmToTop_xyz
    real(DP) :: val(0:iMax-1,jMax,0:kMax)

  end subroutine check_xy_IntSig_BtmToTop_xyz

  subroutine check_xyz_IntSig_SigToTop_xyz
    real(DP) :: val(0:iMax-1,jMax,0:kMax)

  end subroutine check_xyz_IntSig_SigToTop_xyz

end program SpmlUtil_mod_test
