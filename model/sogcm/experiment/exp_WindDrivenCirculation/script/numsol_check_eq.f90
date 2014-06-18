program numsol_check_eq
  use dc_types
  use dc_message

  implicit none

  real(DP), parameter :: &
       & RPlanet = 6371e03, &
       & Omega = 7.292e-05, &
       & hViscCoef = 8e02,  &
       & vViscCoef = 1e-02, &
       & RefDens = 1.024e03, &
       & TotDepth = 5.2e03, &
       & TAU0     = 0.01e0


  integer, parameter :: NY = 200, &
       &                NZ = 100



  real(DP), parameter :: YCOORD_MAX = 10d0
  real(DP), parameter :: ZCOORD_MIN = - 5d0

  real(DP) :: RefHVel0, RefVVel0
  real(DP) :: Ro0, Eh0, Ev0, Reh0
  real(DP) :: EBL_WIDTH, EBL_THICK

  
  real(DP), dimension(NY,NZ) :: yz_U, yz_V, yz_W
  real(DP), dimension(NY,NZ) :: yz_y, yz_z

  real(DP) :: param_y, param_z
  real(DP), parameter :: PI = acos(-1d0)

  !
  !

  call calculate_Params()

  call prepair_coordinate(yz_y, yz_z)
  call calc_W_interiorTop(yz_W(:,1))
  call calc_hVel_UpperEBL(yz_U, yz_V)

contains
  !> @brief 
  !!
  !!
  subroutine calculate_Params()
    
    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    !
    call MessageNotify("M", "numsol_check_eq", "[Velocity Scale]")
    RefHVel0 = TAU0/(2d0*Omega*TotDepth)
    RefVVel0 = TotDepth/RPlanet*RefHVel0
    call MessageNotify("M", "numsol_check_eq", &
         & "U=%f [m/s], W=%f[m/s]", d=(/ RefHVel0, RefVVel0 /)) 

    !
    call MessageNotify("M", "numsol_check_eq", "[Nondimensional Parameters]")
    Ro0 = RefHVel0 / (2d0*Omega*RPlanet)
    Ev0 = vViscCoef/(2d0*Omega*TotDepth**2)
    Eh0 = hViscCoef/(2d0*Omega*RPlanet**2)
    Reh0 = RefHVel0*RPlanet/hViscCoef

    call MessageNotify("M", "numsol_check_eq", &
         & "Ro=%f, Ev0=%f, Eh0=%f, Reh0=%f", d=(/ Ro0, Ev0, Eh0, Reh0 /)) 

    call MessageNotify("M", "numsol_check_eq", "[Equatiorial Boundary Layer]")
    EBL_WIDTH = (Eh0)**(1/3.0)
    EBL_THICK = sqrt(Ev0/EBL_WIDTH)
    call MessageNotify("M", "numsol_check_eq", &
         & "width=%f [km], thickness=%f [m]", d=(/ EBL_WIDTH*RPlanet/1000.0, EBL_THICK*TotDepth /)) 
    
  end subroutine calculate_Params

  !> @brief 
  !!
  !!
  subroutine calc_W_interiorTop(y_W)
    
    ! 宣言文; Declaration statement
    !
    real(DP),intent(inout) :: y_W(NY)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: nintegRet
    integer, parameter :: nTry = 4
    integer, parameter :: NMaxs(nTry) = (/ 200 , 500, 1000, 5000 /)
    integer :: try, j

    ! 実行文; Executable statement
    !
 
    call MessageNotify("M", "numsol_check_eq", "calculate W at the top of interior")
    do j=1, NY
       param_y = yz_y(j,1)
       do try=1, nTry
          nintegRet = numIntegrate(W_interiorTop_integrand, 0d0, 100d0, NMaxs(try))
          if (j==1) then
             write(*,*) "NMax=", NMaxs(try), ":", nintegRet!, (nintegRet-PI)/PI*100.0
          end if
       end do
       y_W(j) = nintegRet
    end do

!    write(*,*) y_W
  end subroutine calc_W_interiorTop

  !> @brief 
  !!
  !!
  subroutine calc_hVel_UpperEBL(yz_U, yz_V)
    
    ! 宣言文; Declaration statement
    !
    real(DP),intent(inout) :: yz_U(NY,NZ), yz_V(NY,NZ)
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: nintegRetU, nintegRetV
    integer, parameter :: nTry = 4
    integer, parameter :: NMaxs(nTry) = (/ 200 , 500, 1000, 10000 /)
    integer :: try, j, k

    ! 実行文; Executable statement
    !
 
    call MessageNotify("M", "numsol_check_eq", "calculate hVel within equatorial upper boundary layer")
    do k=1, NZ
       do j=1, NY
          param_y = yz_y(j,k)
          param_z = yz_z(j,k)
          do try=1, nTry
             nintegRetU = numIntegrate(U_UpperEBL_integrand, 1d-18, 200d0, NMaxs(try))
             nintegRetV = numIntegrate(V_UpperEBL_integrand, 1d-18, 200d0, NMaxs(try))
             if (j==41 .and. k==1) then
                write(*,*) "y=", param_y, "z=", param_z, "NMax=", NMaxs(try), ":", nintegRetU, nintegRetV
             end if
          end do
          yz_U(j,k) = nintegRetU
          yz_V(j,k) = nintegRetV
       end do
    end do

    write(*,*) yz_V(:,1)
    write(*,*) yz_U(1,:)
  end subroutine calc_hVel_UpperEBL

  real(DP) function W_interiorTop_integrand(p)
    real(DP), intent(in) :: p
    
    W_interiorTop_integrand = p*exp(-p**3/3d0)*cos(param_y*p)
!    W_interiorTop_integrand = 4d0/(1d0 + p**2)
  end function W_interiorTop_integrand

  real(DP) function U_UpperEBL_integrand(p)
    real(DP), intent(in) :: p

    real(DP) :: p2

    p2 = p**2
    U_UpperEBL_integrand = - 2d0/sqrt(PI) * exp(-p2**3/3d0 - 0.25d0*param_z**2/p2) * cos(p2*param_y) 

  end function U_UpperEBL_integrand

  real(DP) function V_UpperEBL_integrand(p)
    real(DP), intent(in) :: p

    real(DP) :: p2

    p2 = p**2
    V_UpperEBL_integrand = 2d0/sqrt(PI) * exp(-p2**3/3d0 - 0.25d0*param_z**2/p2) * sin(p2*param_y) 

  end function V_UpperEBL_integrand


!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief 
  !!
  !!
  subroutine prepair_coordinate(yz_y, yz_z)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: yz_y(NY,NZ)
    real(DP), intent(inout) :: yz_z(NY,NZ)
    
    ! 局所変数
    ! Local variables
    !
    integer :: j, k
    real(DP) :: dy, dz

    ! 実行文; Executable statement
    !

    dy = YCOORD_MAX/dble(NY)
    dz = ZCOORD_MIN/dble(NZ)

    do k=1, NZ
       do j=1, NY
          yz_y(j,k) = (j-1)*dy
          yz_z(j,k) = (k-1)*dz
       end do
    end do


  end subroutine prepair_coordinate

  
  !> @brief 
  !!
  !!
  real(DP) function numIntegrate(fx, x0, x1, Nmax)
    
    ! 宣言文; Declaration statement
    !
    interface
       real(DP) function fx(x) 
         use dc_types, only: DP

         real(DP), intent(in) :: x
       end function fx
    end interface
    
    real(DP), intent(in) :: x0, x1
    integer, intent(in) :: Nmax

    ! 局所変数
    ! Local variables
    !
    real(DP) :: dx
    real(DP) :: tmpInt
    integer :: n

    ! 実行文; Executable statement
    !
    
    dx = (x1 - x0)/dble(Nmax)
    tmpInt = dx*(fx(x0) - fx(x1))/6d0

    !$omp parallel do reduction(+:tmpInt)
    do n=1, Nmax
       tmpInt = tmpInt + (2d0*fx(x0+(n-0.5)*dx) + fx(x0 + n*dx))*dx/3d0
    end do

    numIntegrate = tmpInt

  end function numIntegrate
  
end program numsol_check_eq
