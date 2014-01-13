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


!  call check_xy_IntSig_BtmToTop_xyz
!  call check_xyz_IntSig_SigToTop_xyz
  call check_spectral_expand_limit

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
!!$          val(:,:,k) = cos(m*theta(k))
!!$          checkVal(:,:,k) = -sin(m*theta(k))/(m*PI)
          val(:,:,k) = sin(m*theta(k))
          checkVal(:,:,k) = (cos(m*theta(k)) - cos(m*PI))/(m*PI)
       end do
       intVal = xyz_IntSig_SigToTop_xyz(val)

       l2norm = abs(sum(intVal - checkVal))/dble(iMax*jMax*(kMax+1))       
       message=CPrintf("xyz_IntSig_SigToTop_xyz: mode num %d: l2norm %f", i=(/m*2/), d=(/l2norm/))


do k=0,kMax
   write(*,*) "k=",k, ": theta=", theta(k), "val=", val(10,10,k), ": integral ",  intVal(10,10,k), checkVal(10,10,k)
end do
       call  AssertLessThan(message=message, answer=1d-10, check=l2norm)
    end do

  end subroutine check_xyz_IntSig_SigToTop_xyz

  subroutine check_spectral_expand_limit

    use at_module

    real(DP) :: xyz_oriField(0:iMax-1,jMax,0:kMax), xyz_field(0:iMax-1,jMax,0:kMax)    
    real(DP) :: xyz_intSig(0:iMax-1,jMax,0:kMax), z_intSig_ana(0:kmax)
    real(DP) :: lb   !< Nondimensional characteristic length of boundary layer
    integer :: i

    do i=5,1,-1

    
       lb = abs( g_Sig(i) )
       xyz_oriField(1,1,:) = exp( g_Sig/lb ) - exp(- (g_Sig + 1d0)/lb )
!!$       xyz_orifield = 0d0
!!$       xyz_orifield(:,:,0) = 1d0; xyz_orifield(:,:,kMax)=-1d0
       write(*,*) "=========="
       write(*,*) "* blyrLength:", lb
       write(*,*) xyz_orifield(1,1,:)
!!$       write(*,*) "* Difference between z_field and one after g_t(t_g(z_field))"
!!$       xyz_field(1,1,:) = g_t(t_g(xyz_oriField(1,1,:)))
!!$       write(*,*) xyz_field(1,1,:)
       write(*,*) "Integ:", sum(g_x_weight*xyz_orifield(1,1,:))
       xyz_intSig = xyz_IntSig_SigToTop_xyz(xyz_orifield)
       z_intSig_ana = lb*((1d0+exp(-1d0/lb))-(exp(g_Sig/lb)+exp(-(g_Sig+1d0)/lb)))
!!$       write(*,*) "IntSig", xyz_intSig(1,1,:)
       write(*,*) "IntSig RMS Error", &
            & sqrt( sum( ( xyz_intSig(1,1,:) -  z_intSig_ana )**2 ) )/kMax
    end do
  
write(*,*) "g_Sig:", g_Sig

  end subroutine check_spectral_expand_limit

end program SpmlUtil_mod_test
