program numsol_check_eq

  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_history, only: &
       & HistoryCreate, HistoryClose, HistoryAddVariable, &
       & HistoryPut

  implicit none

  character(*), parameter :: OutputFileName = "EqRefSol_Ah800.nc"

  real(DP), parameter :: &
       & RPlanet = 6371d03, &
       & Omega = 7.292d-05, &
       & hViscCoef = 1d04,  &
       & vViscCoef = 1d-02, &
       & RefDens = 1.027d03, &
       & TotDepth = 5.2d03, &
       & TAU0     = 0.0078d0


  integer, parameter :: NY = 100, &
       &                NZ = 100



  real(DP), parameter :: YCOORD_MAX = 10d0
  real(DP), parameter :: ZCOORD_MIN = - 5d0

  real(DP) :: RefHVel0, RefVVel0
  real(DP) :: Ro0, Eh0, Ev0, Reh0
  real(DP) :: EBL_WIDTH, EBL_THICK

  
  real(DP), dimension(NY,NZ) :: yz_U, yz_V, yz_W
  real(DP), dimension(NY,NZ) :: yz_y, yz_z

  real(DP) :: param_y, param_z, param_p
  
  real(DP), allocatable :: p_integGauss(:)
  real(DP), parameter :: PI = acos(-1d0)

  !
  !

  call calculate_Params()
  
  call prepair_output()
  call prepair_coordinate(yz_y, yz_z)
  call calc_Vel_UpperEBL(yz_U, yz_V, yz_W)

  call HistoryClose()

contains

  subroutine prepair_output()

    
    call HistoryCreate( & 
         & file= trim(OutputFileName), title='Reference Solution of circulation near equator', &
         & source='Dennou-OGCM', &
         & institution='Dennou-OGCM project', &
         & dims=(/'y  ', 'lat', 'z  ', 'sig' /), dimsizes=(/ NY, NY, NZ, NZ /), &
         & longnames=(/ 'y  ', 'lat', 'z  ', 'sig' /),&
         & units=(/ '1           ', 'degree_north', '1           ', '1           ' /) &
         &  )  


       call HistoryAddVariable( & 
            & varname='U_nondim', dims=(/ 'y  ', 'sig' /), &
            & longname='zonal velocity', units='1', xtype='double' &
            & )

       call HistoryAddVariable( & 
            & varname='U', dims=(/ 'lat', 'sig' /), &
            & longname='zonal velocity', units='m/s', xtype='double' &
            & )

       call HistoryAddVariable( & 
            & varname='V_nondim', dims=(/ 'y  ', 'sig' /), &
            & longname='meridional velocity', units='1', xtype='double' &
            & )

       call HistoryAddVariable( & 
            & varname='V', dims=(/ 'lat', 'sig' /), &
            & longname='meridional velocity', units='m/s', xtype='double' &
            & )

       call HistoryAddVariable( & 
            & varname='W_nondim', dims=(/ 'y  ', 'sig' /), &
            & longname='vertical velocity', units='1', xtype='double' &
            & )

       call HistoryAddVariable( & 
            & varname='W', dims=(/ 'lat', 'sig' /), &
            & longname='vertical velocity', units='m/s', xtype='double' &
            & )


  end subroutine prepair_output

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
    Ev0 = vViscCoef / ( 2d0*Omega*TotDepth**2 )
    Eh0 = hViscCoef / ( 2d0*Omega*RPlanet**2 )

    EBL_WIDTH = (Eh0)**(1.0/3.0)
    EBL_THICK = sqrt(Ev0/EBL_WIDTH)

    RefHVel0 = TAU0 / ( RefDens*2d0*Omega*EBL_THICK*TotDepth ) 
    RefVVel0 = TotDepth*EBL_THICK/RPlanet*RefHVel0

    Ro0 = RefHVel0 / ( 2d0*Omega*RPlanet )
    Reh0 = RefHVel0*RPlanet/hViscCoef


    call MessageNotify("M", "numsol_check_eq", "[Nondimensional Parameters]")

    call MessageNotify("M", "numsol_check_eq", &
         & "Ro=%f, Ev0=%f, Eh0=%f, Reh0=%f", d=(/ Ro0, Ev0, Eh0, Reh0 /)) 

    !
    call MessageNotify("M", "numsol_check_eq", "[Velocity Scale]")
    call MessageNotify("M", "numsol_check_eq", &
         & "U=%f [m/s], W=%f[m/s]", d=(/ RefHVel0, RefVVel0 /)) 
write(*,*) 'U*=', RefHVel0/EBL_WIDTH
write(*,*) 'EBL_THICK0=', sqrt(Ev0)*TotDepth, ',EBL_THICK=', EBL_THICK*TotDepth
write(*,*) 'McKee73 paramters: E=', EBL_THICK**2/EBL_WIDTH, 'R=', Ro0/Eh0

    call MessageNotify("M", "numsol_check_eq", "[Equatiorial Boundary Layer]")
    call MessageNotify("M", "numsol_check_eq", &
         & "width=%f [km], thickness=%f [m]", d=(/ EBL_WIDTH*RPlanet/1000.0, EBL_THICK*TotDepth /)) 
    
  end subroutine calculate_Params



  !> @brief 
  !!
  !!
  subroutine calc_Vel_UpperEBL(yz_U, yz_V, yz_W)
    use integrand_mod, only: &
         & integrand_set, &
         & gauss_integrand, &
         & U_UpperEBL_integrand, &
         & V_UpperEBL_integrand, &
         & W_UpperEBL_integrand

    implicit none

    ! 宣言文; Declaration statement
    !
    real(DP),intent(inout) :: yz_U(NY,NZ), yz_V(NY,NZ), yz_W(NY,NZ)
    
    ! 局所変数
    ! Local variables
    !
    abstract interface 
       real(DP) function fx(x, xidx) 
         use dc_types, only: DP

         real(DP), intent(in) :: x
         integer, intent(in) :: xidx
       end function fx
    end interface
    procedure (fx), pointer :: integrand_ptr => null()

    real(DP) :: nintegRetU, nintegRetV, nintegRetW
    integer, parameter :: nTry = 5
    integer, parameter :: NMaxs(nTry) = (/ 200 , 500, 1000, 5000, 10000 /)
    integer :: try, j, k, m
    real(DP), parameter :: PMAX = 200d0
    real(DP) :: delp

    ! 実行文; Executable statement
    !
 
    call MessageNotify("M", "numsol_check_eq", "calculate hVel within equatorial upper boundary layer")
    do k=1, NZ

       call MessageNotify("M", "numsol_check_eq", "k: %d/%d", &
            & i=(/ k, NZ /) )
       param_z = yz_z(1,k)
       
       do try=nTry, nTry

          if(allocated(p_integGauss)) deallocate(p_integGauss)

          allocate( p_integGauss(0:2*NMaxs(try)) )
          call integrand_set(param_y, param_z, param_p, p_integGauss)

          delp = PMAX/dble(NMaxs(try))
          do m=0, 2*NMaxs(try)
             param_p = 0.5d0*m*dp + 1d-30
             integrand_ptr => gauss_integrand
             p_integGauss(m) = numIntegrate(integrand_ptr, param_z, 0d0, int(1000*abs(param_z)) )

          end do
!!$if( k==2 ) then
!!$   write(*,*) p_integGauss
!!$endif

          do j=1, NY
             param_y = yz_y(j,k)

             integrand_ptr =>  U_UpperEBL_integrand
             nintegRetU = numIntegrate(U_UpperEBL_integrand, 1d-18, PMAX, NMaxs(try))

             integrand_ptr =>  V_UpperEBL_integrand
             nintegRetV = numIntegrate(V_UpperEBL_integrand, 1d-18, PMAX, NMaxs(try))

             integrand_ptr =>  W_UpperEBL_integrand
             nintegRetW = numIntegrate(W_UpperEBL_integrand, 1d-18, PMAX, NMaxs(try))
             
             if (j==41 .and. k==1) then
                write(*,*) "y=", param_y, "z=", param_z, "NMax=", NMaxs(try), ":", &
                     & nintegRetU, nintegRetV, nintegRetW
             end if

             yz_U(j,k) = nintegRetU - 0.0689278d0
write(*,*) "hoge"
             yz_V(j,k) = nintegRetV
             yz_W(j,k) = nintegRetW
          end do

       end do
    end do

    !write(*,*) yz_V(:,1)
    !write(*,*) yz_U(1,:)
    
    call HistoryPut('U_nondim', yz_U)
    call HistoryPut('V_nondim', yz_V)
    call HistoryPut('W_nondim', yz_W)
    call HistoryPut('U', yz_U*RefHVel0/EBL_WIDTH)
    call HistoryPut('V', yz_V*RefHVel0/EBL_WIDTH)
    call HistoryPut('W', yz_W*(RefVVel0/TotDepth)/EBL_WIDTH**2)

  end subroutine calc_Vel_UpperEBL



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



    call HistoryPut('y', yz_y(:,1))
    call HistoryPut('lat', yz_y(:,1)*EBL_WIDTH*180d0/PI)
    call HistoryPut('z', yz_z(1,:))
    call HistoryPut('sig', yz_z(1,:)*EBL_THICK)

  end subroutine prepair_coordinate

  
  !> @brief 
  !!
  !!
  real(DP) function numIntegrate(fx, x0, x1, Nmax)

    implicit none
    
    ! 宣言文; Declaration statement
    !
    interface 
       real(DP) function fx(x, xidx) 
         use dc_types, only: DP

         real(DP), intent(in) :: x
         integer, intent(in) :: xidx
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
    
    if (x1 == x0 ) then
       numIntegrate = 0d0; return
    end if

    dx = (x1 - x0)/dble(Nmax)
    tmpInt = dx*(fx(x0,0) - fx(x1,2*Nmax))/6d0

    !$omp parallel do reduction(+:tmpInt)
    do n=1, Nmax
       tmpInt = tmpInt + (2d0*fx(x0+(n-0.5)*dx, 2*n-1) + fx(x0 + n*dx, 2*n))*dx/3d0
    end do

    numIntegrate = tmpInt

  end function numIntegrate
  
end program numsol_check_eq
