program SpmlUtil_mod_test
  
  use dc_types
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use dc_test
  use dc_string

  implicit none

  character(*), parameter :: configNmlFile = "defaultConfig.nml"

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
    real(DP) :: intVal(0:iMax-1,jMax)

    integer :: k, m
    real(DP) :: theta(0:kMax), l2norm
    character(STRING) :: message
    real(DP) :: checkVal
 
    theta = PI*(1d0 + g_Sig)
    do m=1,5
       do k=0, kMax
          val(:,:,k) = cos(m*theta(k))
       end do
       intVal = xy_IntSig_BtmToTop_xyz(val)
       checkVal = sin(m*PI)
       l2norm = abs(sum(intVal - checkVal))/dble(iMax*jMax)
       message=CPrintf("xy_IntSig_BtmToTop_xyz: mode num=%d: l2norm=%f", i=(/m*2/), d=(/ l2norm /))
       call  AssertLessThan(message=message, &
            & answer=5d-10, check=l2norm)
    end do
    

  end subroutine check_xy_IntSig_BtmToTop_xyz

  subroutine check_xyz_IntSig_SigToTop_xyz
    real(DP) :: val(0:iMax-1,jMax,0:kMax)
    real(DP) :: intVal(0:iMax-1,jMax, 0:kMax)
    real(DP) :: checkVal(0:iMax-1,jMax, 0:kMax)
    integer :: k, m
    real(DP) :: theta(0:kMax), l2norm
    character(STRING) :: message
 
    theta = PI*(1d0 + g_Sig)
    do m=1,5
       do k=0, kMax
          val(:,:,k) = cos(m*theta(k))
          checkVal(:,:,k) = -sin(m*theta(k))/(m*PI)
       end do
       intVal = xyz_IntSig_SigToTop_xyz(val)

       l2norm = abs(sum(intVal - checkVal))/dble(iMax*jMax*(kMax+1))       
       message=CPrintf("xyz_IntSig_SigToTop_xyz: mode num %d: l2norm %f", i=(/m*2/), d=(/l2norm/))

!!$do k=0,kMax
!!$   write(*,*) "k=",k,intVal(10,10,k), checkVal(10,10,k)
!!$end do
       call  AssertLessThan(message=message, answer=1d-10, check=l2norm)
    end do

  end subroutine check_xyz_IntSig_SigToTop_xyz

end program SpmlUtil_mod_test
