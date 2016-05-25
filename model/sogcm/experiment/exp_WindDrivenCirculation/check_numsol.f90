program check_numsol
  use dc_types
  
  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use VariableSet_mod
  use BoundaryCondO_mod
  use gtool_history

  use Exp_WindDrivenCirculation_mod, only: &
         & Exp_Init => Exp_WindDrivenCirculation_Init, &
         & Exp_Final => Exp_WindDrivenCirculation_Final, &
         & Exp_SetInitCond => SetInitCondition, &
         & construct_WindStressU_Marshall07

  implicit none

  real(DP), parameter :: TAU0 = 0.12d0
  real(DP), parameter :: eps = 4d-3
  real(DP), parameter :: TotDepth = 5.2e03
  real(DP), parameter :: Ah = 1d3, Av = 1d-2

  real(DP) :: Eh0, Ev0, r_hv

  character(*), parameter :: configNmlFile = 'defaultConfig.nml'
  character(*), parameter :: resolSig = "T42"
  character(*), parameter :: U_ncfilePath = '/home/ykawai/exp_backup/current/exp_Ah1e3T682L60/Run1_U.nc'
  character(*), parameter :: V_ncfilePath = '/home/ykawai/exp_backup/current/exp_Ah1e3T682L60/Run1_V.nc'
  integer, parameter :: NLAT2=289

 
  !
  call setup()
  call set_parameters()

  !
  !
  call calc_anaSol()

  !
  call shutdown()

contains

  subroutine set_parameters()

    Eh0 = 2d0*Ah/(2d0*Omega*RPlanet**2)
    Ev0 = 2d0*Av/(2d0*Omega*TotDepth**2)
    r_hv = 2d0*Eh0/sqrt(Ev0)
    write(*,*) "Eh0=", Eh0, ", Ev0=", Ev0, ", r_hv=", r_hv
    
  end subroutine set_parameters

  subroutine calc_anaSol()
    
    real(DP), dimension(lMax) :: w_Forcing, w_Zeta, w_Zero, w_Psi0, w_Psi1
    real(DP), dimension(0:iMax-1,jMax) :: xy_fmod, xy_Zero, xy_u0, xy_u1A, xy_u1N, xy_vDiffZeta, xy_hDiffZeta, &
         & xy_Zeta, xy_Psi0, xy_PlVorAdv, xy_f, xy_SWindTrque, xy_Friction, xy_Adv, xy_tot
    real(DP) :: xy_mu(0:iMax-1,jMax), xy_CosLat(0:iMax-1,jMax)
    real(DP) :: VelScale
    integer :: j, n

    call HistoryGet(U_ncfilePath, 'U', xyz_UN, range='t=2000' )
    call HistoryGet(V_ncfilePath, 'V', xyz_VN, range='t=2000' )
    velScale = 2d0*TAU0/(RefDens*2d0*Omega*sqrt(Ev0)*TotDepth)


    xy_mu = sin(xyz_Lat(:,:,0))
    xy_CosLat = cos(xyz_Lat(:,:,0))
    xy_Zero = 0d0; w_Zero = 0d0

    xy_u0 = ( xy_WindStressU/TAU0 * abs(xy_mu)/( xy_mu**2 + eps**2 ) ) * sqrt(abs(xy_mu))
    w_Psi0 = w_InvLapla2D_w( w_AlphaOptr_xy(xy_Zero, -xy_u0*xy_CosLat) ) / RPlanet
    xy_Psi0 = xy_w(w_Psi0)
    write(*,*) xy_Psi0(0, 128), xy_u0(0, 128)*velScale

    w_Psi0 = w_InvLapla2D_w( w_AlphaOptr_xy(xy_Zero, -xyz_UN(:,:,(kMax+1)/2)*xy_CosLat) )
    xy_Psi0 = xy_w(w_Psi0)
    write(*,*) xy_Psi0(0, 128), xyz_UN(:,128,(kMax+1)/2)
    
    xy_f = 2d0*Omega*xy_mu

    xy_hDiffZeta =   TotDepth * xy_w( Ah*w_Lapla2D_w(2d0*w_Psi0/RPlanet**2 + w_Lapla2D_w(w_Psi0)) )

    xy_PlVorAdv  = - TotDepth * 2d0*Omega/RPlanet*xy_CosLat*xyz_VN(:,:,(kMax+1)/2)

    xy_SWindTrque = xy_f * ( xy_w(w_AlphaOptr_xy(xy_Zero, -xy_WindStressU*xy_CosLat* xy_f/(xy_f**2 + 1e-5**2 )/RefDens )) )

    xy_Friction = - xy_f * ( xy_w(w_AlphaOptr_xy(xy_Zero, -xyz_UN(:,:,(kMax+1)/2)*xy_CosLat &
         & * (sqrt(Ev0* abs(xy_mu)/(xy_mu**2 + 8e-2**2)))*TotDepth/2d0 &
         & ))) 

    xy_Adv = - TotDepth* xyz_VN(:,:,(kMax+1)/2) * xy_GradLat_w(w_Lapla2D_w(w_Psi0))/RPlanet
    xy_tot = xy_PlVorAdv + xy_SWindTrque + xy_Friction + xy_hDiffZeta + xy_Adv


    do j=1, jMax
       if (mod(j,4) == 0) then
          write(*,*) "lat=", xyz_Lat(0,j,0)*180/PI, &
               & xy_SWindTrque(0,j), xy_Friction(0,j), xy_hDiffZeta(0,j), xy_PlVorAdv(0,j), xy_Adv(0,j), &
               & ":", xy_tot(0,j)
       end if
    end do
    xy_hDiffZeta = xy_w( r_hv*w_Lapla2D_w(w_Lapla2D_w(w_Psi0) ) ) * RPlanet**4
    xy_PlVorAdv = xy_CosLat*xyz_VN(:,:,(kMax+1)/2) /velScale * 2d0/sqrt(Ev0) 
stop

    w_Zeta = 0d0
    xy_Zeta = 0d0
    xy_u1N = 0d0
    do n=1, 5000000
       
       xy_vDiffZeta = abs(xy_mu) * xy_AlphaOptr_w( &
            &    w_Zero, w_xy(- xy_u1N*xy_CosLat/(sqrt(abs(xy_mu))*(1d0 + eps**2/xy_mu**2)) ) &
            & ) *RPlanet
       xy_Zeta = xy_Zeta + 1d-3*( ( &
            & - xy_vDiffZeta  &
            & + xy_hDiffZeta &
            & - xy_PlVorAdv &
            & ))

       w_Psi1 = w_InvLapla2D_w(w_xy(xy_Zeta)) /RPlanet**2
       xy_u1A  = xy_CosLat*xy_AlphaOptr_w(w_Zero, -w_Psi1) *RPlanet
       if(mod(n,500) == 0) then 
          write(*,*) "n=", n, ":", AvrLonLat_xy((xy_u1A - xy_u1N)**2)
          write(*,*) "**", xy_u0(0,128), xy_u1A(0,128), xy_hDiffZeta(0,128), -xy_vDiffZeta(0,128), -xy_PlVorAdv(0,128)
          write(*,*) "**", xy_u0(0,250), xy_u1A(0,250), xy_hDiffZeta(0,250), -xy_vDiffZeta(0,250), -xy_PlVorAdv(0,250)
          write(*,*) "**", xy_u0(0,385), xy_u1A(0,385), xy_hDiffZeta(0,385), -xy_vDiffZeta(0,385), -xy_PlVorAdv(0,385)
       end if
       xy_u1N = xy_u1A
    end do


    write(*,*) "velScale=", velScale/sqrt(sin(45d0/180d0*PI))
    do j=1, jMax/2+2
       write(*,*) "lat=", xyz_Lat(0,j,0)*180/PI, xy_u0(:,j) * velScale, &
            & velScale*xy_WindStressU(0,j)/TAU0/sqrt(abs(sin(xyz_Lat(0,j,0)))), &
            & xyz_UN(0,j, (kMax+1)/2), &
            & xy_u0(:,j)*velScale - xyz_UN(0,j,(kMax+1)/2)!/velScale
    end do

  end subroutine calc_anaSol
  
  function dZeta1dt(w_Psi, w_Forcing) result(w_dZetadt)
    real(DP), intent(in) :: w_Psi(lMax)
    real(DP), intent(in) :: w_Forcing(lMax)
    real(DP) :: w_dZetadt(lMax)

    w_dZetadt = w_Forcing + w_xy(sqrt(sin(xyz_Lat(:,:,0))) * xy_w(w_Lapla2D_w(w_Psi))*RPlanet**2 )
  end function dZeta1dt

  subroutine setup()
    use InitCond_mod

    integer :: j

    call Constants_Init(configNmlFile)
    call GridSet_Init(configNmlFile)
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
    call GridSet_construct()
    call VariableSet_Init()
    call BoundaryCondO_Init()
    call Exp_Init(configNmlFile)
    call InitCond_Init()
    call InitCond_Set(Exp_SetInitCond)
    call InitCond_Final()

  end subroutine setup

  subroutine shutdown()
    call Exp_Final()
    call BoundaryCondO_Final()    
    call VariableSet_Final()
    call SpmlUtil_Final()
    call GridSet_Final()
    call Constants_Final()
  end subroutine shutdown
end program check_numsol
