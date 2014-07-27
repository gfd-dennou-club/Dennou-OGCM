#!/usr/bin/env ruby


require "numru/gphys"
include NumRu
include NMath

Av = 1e-02
#Ah = 1e5
#Ah = 1e4
#Ah = 1e3
Ah = 1e2
HAh = 0e14
TAU0 = 0.121
Dens = 1.027e03
NLAT = 721
NZ   = 501
PI = acos(-1.0)
Omega = 7.292115e-5
Depth = 5.2e03
RPlanet = 6371e03
TIME = 5000

SphericalDOptr=true

##
=begin
GENNCFILE_SUFFIX = "_Ah1e3ref"
EXPDATADIR = ""
GRIDNCFILE = "/home/ykawai/workspace/DQGModel/model/sogcm/exp/Run1_U.nc"
=end

#=begin
#EXPSPECIFIC="Ah1e5Pl341L60"
#EXPSPECIFIC="Ah1e4Pl341L90"
#EXPSPECIFIC="Ah1e3Pl682L90"
EXPSPECIFIC="Ah1e2Pl341L90"
RUNSPECIFIC="Run1"
GENNCFILE_SUFFIX = "_#{EXPSPECIFIC}"
EXPDATADIR = "/home/ykawai/exp_backup/current/exp_#{EXPSPECIFIC}"
GRIDNCFILE = "#{EXPDATADIR}/#{RUNSPECIFIC}_U.nc"
#=end

##########################

def set_parameters(lat)
  $L = 0.4*RPlanet
  $Eh0 = 2*Ah/( 2*Omega*$L**2 )
  $Ev0 = 2*Av/( 2*Omega*Depth**2 )
  $UVel = 2*TAU0/( Dens*2*Omega*Depth*Math.sqrt($Ev0) )
  $Ro = $UVel/( 2*Omega*$L )
  $r_hv = $Eh0/sqrt($Ev0)
  $ddlat_scle = $L/RPlanet

  puts "Eh0=#{$Eh0}"
  puts "Ev0=#{$Ev0}"
  puts "UVelScale=#{$UVel}"
  puts "Ro0=#{$Ro}"
  puts "r_hv=#{$r_hv}"
end

def gen_taux(lat)
  coef = NArray.to_na( [
    0.0265682, -0.0784899, -0.00880389, 0.0343205, 0.0233334, 0.000641955, -0.00387676, -0.00150998
  ] )

  taux = 0.0
  for i in 0..coef.length-1
    taux = taux + coef[i]*cos( (2.0*i+1.0)*lat )
  end

#  p "gen_taux:", taux
  return taux
end

def gen_taux_dnLat(lat, n)
  coef = NArray.to_na( [
    0.0265682, -0.0784899, -0.00880389, 0.0343205, 0.0233334, 0.000641955, -0.00387676, -0.00150998
  ] )

  taux_dnLat = 0.0
  if n%2 == 0 then
    sign = (-1.0)**(n/2)
#    puts "n=#{n}, sign=#{sign}"
    for i in 0..coef.length-1
      taux_dnLat = taux_dnLat + sign * (2.0*i+1.0)**n * coef[i]*cos( (2.0*i+1.0)*lat )
    end
  else
    sign = (-1.0)**((n+1)/2)
#    puts "n=#{n}, sign=#{sign}"

    for i in 0..coef.length-1
      taux_dnLat = taux_dnLat + sign * (2.0*i+1.0)**n * coef[i]*sin( (2.0*i+1.0)*lat )
    end
  end
#  p "gen_taux:", taux_dnLat
  return taux_dnLat
end

def calc_u0_interior(lat, isNonDim=false, n_dlat=0)
  
  u = nil
  m = isNonDim ? $ddlat_scle : 1.0

  taux = gen_taux(lat)/TAU0
  taux_dlat1 = gen_taux_dnLat(lat,1)* m/TAU0
  taux_dlat2 = gen_taux_dnLat(lat,2)* m**2/TAU0
  taux_dlat3 = gen_taux_dnLat(lat,3)* m**3/TAU0
  taux_dlat4 = gen_taux_dnLat(lat,4)* m**4/TAU0
  sinLat = sin(lat)
  sqrtSinLat = sqrt(sinLat.abs)

  u00 = taux/sqrtSinLat
  u00_dlat1 = ( taux_dlat1 - m*taux/(2.0*tan(lat)) )/sqrtSinLat

  case n_dlat
  when 0 then
    u =   u00  \
        + 0*$Ro*( \
                0.25*taux*( taux_dlat1 - tan(lat)*taux) \
                + taux*( u00_dlat1 - tan(lat)*u00 ) \
              )/sqrtSinLat \
        - 0*$Ro* 7.0/20.0*u00*( u00_dlat1 - tan(lat)*u00 )
  when 1 then
    u = u00_dlat1
  when 2 then
    u = ( taux_dlat2 - m*taux_dlat1/tan(lat) + m**2*0.25*taux*(-1.0 + 3.0/sinLat**2) )/sqrtSinLat
  when 3 then
    u = ( taux_dlat3 - m*taux_dlat2/tan(lat) + m**2*taux_dlat1/sinLat**2 + m**2*0.25*taux_dlat1*(- 1.0 + 3.0/sinLat**2) - m**3*0.25*taux*6.0*cos(lat)/sinLat**3 )/sqrtSinLat \
        - m*calc_u0_interior(lat, true, 2)/(2.0*tan(lat))
  when 4 then
    u = ( taux_dlat4 - m*taux_dlat3/tan(lat) + m**2*2*taux_dlat2/sinLat**2 + m**2*0.25*taux_dlat2*(- 1.0 + 3.0/sinLat**2) - m**3*0.25*taux_dlat1*12.0*cos(lat)/sinLat**3 + m**4*0.25*taux*6.0*(1.0 + 2.0*cos(lat)**2)/sin(lat)**4 )/sqrtSinLat \
        - m*( 2*calc_u0_interior(lat, true, 3) + m*calc_u0_interior(lat,true,2)/(2*tan(lat)) )/(2.0*tan(lat)) \
        + m**2*calc_u0_interior(lat, true, 2) / (2.0*sin(lat)**2)

  else
      puts "n_dlat must be less than 4."; exit
  end

  return isNonDim ? u : $UVel*u
end

def integerand_advTerm(lat)
  u0 = calc_u0_interior(lat, true, 0)
  u0_lat = calc_u0_interior(lat, true, 1)
  u0_latlat = calc_u0_interior(lat, true, 2)

  zeta0_dlat1 = - (u0_latlat - u0_lat*tan(lat) - u0/cos(lat)**2)
  integrand = ( (zeta0_dlat1 - 2*u0)/sin(lat) )*zeta0_dlat1 *  cos(lat)/sin(lat)
  return  integrand
end

def calc_u1_interior(lat, isNonDim=false, n_dlat=0, isNondim_ndlat=false)
  
  nLat = lat.length
  u = NArray.float(lat.length)
  m = isNonDim ? $ddlat_scle : 1.0

  u0 = calc_u0_interior(lat, true, 0)
  u0_dlat1 = calc_u0_interior(lat, true, 1)
  u0_dlat2 = calc_u0_interior(lat, true, 2)
  u0_dlat3 = calc_u0_interior(lat, true, 3)
  u0_dlat4 = calc_u0_interior(lat, true, 4)

  zeta0_dlat1 = - (u0_dlat2 - m*u0_dlat1*tan(lat) - m**2*u0/cos(lat)**2)
  zeta0_dlat2 = - (u0_dlat3 - m*u0_dlat2*tan(lat) - m**2*2.0*u0_dlat1/cos(lat)**2  - m**3*u0* 2.0*tan(lat)/cos(lat)**2)
  zeta0_dlat3 = - (u0_dlat4 - m*u0_dlat3*tan(lat) - m**2*3.0*u0_dlat2/cos(lat)**2  - m**3*6.0*u0_dlat1*tan(lat)/cos(lat)**2 - m**4*u0* 2.0*(1.0 + 3.0*tan(lat)**2)/cos(lat)**2) 
  taux = gen_taux(lat)/TAU0
  sinLat = sin(lat)

  ekRatio = 1.0 / (sqrt(sinLat.abs))
  ekRatio_dlat1 = - m*ekRatio/(2*tan(lat)) 

  case n_dlat
  when 0 then
    u = - ekRatio * (zeta0_dlat1 - m**2 * 2*u0)
  when 1 then
    u = - ekRatio * ( zeta0_dlat2 - m**2 * 2*u0_dlat1) \
        - ekRatio_dlat1 * ( zeta0_dlat1 - m**2 * 2*u0 )

  when 2 then
    u = - ekRatio*( zeta0_dlat3 - m**2 * 2*u0_dlat2) \
        - 2*ekRatio_dlat1 * (zeta0_dlat2 - m**2 * 2*u0_dlat1 ) \
        - m**2* ekRatio * (zeta0_dlat1 - m**2 * 2*u0) * ( 1/(2*sin(lat)**2) + 1/(4*tan(lat)**2) )

  else
    puts "n_dlat must be less than 3."; exit
  end

  return isNonDim ? u : $r_hv * $UVel * u
end

def calc_u2_interior(lat, isNonDim=false, n_dlat=0, isNondim_ndlat=false)
  
  u2 = NArray.float(lat.length)
  m = isNonDim ? $ddlat_scle : 1.0
  
  u1 = calc_u1_interior(lat, true, 0)
  u1_dlat1 = calc_u1_interior(lat, true, 1)
  u1_dlat2 = calc_u1_interior(lat, true, 2)

  sinLat = sin(lat)
  zeta1_dlat1 = - (u1_dlat2 - m*u1_dlat1*tan(lat) - m**2*u1/cos(lat)**2)
  case n_dlat
  when 0 then
    u2 = -  (zeta1_dlat1 - m**2 * 2*u1) / sqrt(sinLat.abs) 
  else
    puts "n_dlat must be less than 1."; exit
  end

  return isNonDim ? u2 : $r_hv**2 * $UVel * u2

end

def calc_interiorV0(lat, isNonDim=false, n_dlat=0, isNondim_ndlat=false)

  v = NArray.float(lat.length)
  m = isNonDim ? $ddlat_scle : 1.0

  u0 = calc_u0_interior(lat, true, 0)
  u0_dlat1 = calc_u0_interior(lat, true, 1)
  u0_dlat2 = calc_u0_interior(lat, true, 2)
  u0_dlat3 = calc_u0_interior(lat, true, 3)

  case n_dlat
  when 0 then
    zeta0_dlat1 = - ( u0_dlat2 - m*u0_dlat1*tan(lat) - m**2*u0/cos(lat)**2 )
    v = 0.5*$Eh0 * (-zeta0_dlat1 + m**2 * 2*u0) / ( - sin(lat) )
  when 1 then
    zeta0_dlat2 = - ( u0_dlat3 - m*u0_dlat2*tan(lat) - m**2*2*u0_dlat1/cos(lat)**2 - m**3*2*u0*tan(lat)/cos(lat)**2 )
    v = - 0.5*$Eh0 * (-zeta0_dlat2 + m**2 * 2.0*u0_dlat1) / sin(lat) \
        - m*calc_interiorV0(lat, true, 0) / tan(lat)
  else
    puts "n_dlat must be less than 2."; exit
  end

  return isNonDim ? v : $UVel*v

end

def calc_interiorV1(lat, isNonDim=false, n_dlat=0, isNondim_ndlat=false)

  v = NArray.float(lat.length)
  m = isNonDim ? $ddlat_scle : 1.0

  u1 = calc_u1_interior(lat, true, 0)
  u1_dlat1 = calc_u1_interior(lat, true, 1)
  u1_dlat2 = calc_u1_interior(lat, true, 2)

  case n_dlat
  when 0 then
    zeta1_dlat1 = - ( u1_dlat2 - m*u1_dlat1*tan(lat) - m**2*u1/cos(lat)**2 )
    v = - 0.5*$Eh0 * (-zeta1_dlat1 + m**2 * 2*u1) / sin(lat)
  when 1 then
  else
    puts "n_dlat must be less than 2."; exit
  end

  return isNonDim ? v : $UVel * $r_hv * v

end

def calc_u(lat, sig)

  u = NArray.float(lat.length, sig.length)

  deltaE = sqrt( 2.0*Av/($f.abs*Depth**2) )
  taux = gen_taux(lat)/TAU0
  barotU0 = calc_u0_interior(lat, true)
  barotU1 = calc_u1_interior(lat, true)
  barotU2 = calc_u2_interior(lat, true)
  barotU = barotU0 + $r_hv*barotU1 + $r_hv**2*barotU2

for n in 0..(lat.length-1)/2
  puts "lat=#{lat[n]*180/PI}, #{barotU0[n]}, #{barotU1[n]*$r_hv**0}, #{barotU2[n]*$r_hv**0}" 
end


  alpha = 2*TAU0/(Dens*sqrt(2*Av*2*Omega*sin(lat).abs)*$UVel)
  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
    u[true,k] = $UVel*( \
                  barotU*(1.0 - exp(-xiB) * cos(xiB) ) \
                + 0.5*alpha* taux * exp(-xiT) * sqrt(2)*cos(xiT + PI/4.0) \
               )
  end

  return u
end

def calc_v(lat, sig)

  v = NArray.float(lat.length, sig.length)

  deltaE = sqrt( 2.0*Av/($f.abs*Depth**2) )
  taux = gen_taux(lat)/TAU0

  barotU0 = calc_u0_interior(lat, true)
  barotU1 = calc_u1_interior(lat, true)
  barotU2 = calc_u2_interior(lat, true)
  barotU = barotU0 + $r_hv*barotU1 + $r_hv**2*barotU2

  v0_interior = calc_interiorV0(lat, true)
  v1_interior = calc_interiorV1(lat, true)
  v_interior = v0_interior + $r_hv*v1_interior

  alpha = 2*TAU0/(Dens*sqrt(2*Av*2*Omega*sin(lat).abs)*$UVel)
  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
    v[true,k] = $UVel*( \
                  (barotU*sin(xiB) - v_interior*cos(xiB))*exp(-xiB) \
                + v_interior \
                - 0.5*alpha * taux * exp(-xiT) * Math.sqrt(2)*sin(xiT + PI/4.0)  \
               )
  end

  return v
end

def calc_sigDot(lat, sig)

  sigDot = NArray.float(lat.length, sig.length)
  deltaE = sqrt(2.0*Av/($f.abs*Depth**2))

  m = $ddlat_scle
  taux = gen_taux(lat)/TAU0
  taux_dlat1 = gen_taux_dnLat(lat, 1)* m/TAU0

  u = calc_u0_interior(lat, true, 0) + $r_hv*calc_u1_interior(lat, true, 0)
  u_dlat1 = calc_u0_interior(lat, true, 1) + $r_hv*calc_u1_interior(lat, true, 1)

  v = calc_interiorV0(lat, true) + $r_hv*calc_interiorV1(lat, true)
  v_dlat1 = calc_interiorV0(lat, true, 1) + $r_hv*calc_interiorV1(lat, true, 1)

  v_interior_divy = (v_dlat1 - m*v*tan(lat))/m

  uekmanMassFlux_divy = 0.5*sqrt($Ev0)/m *(taux_dlat1 - m*taux/tan(lat) - m*taux*tan(lat))/sin(lat)
  lekmanMassFlux_divy = 0.5*sqrt($Ev0)/m *( \
           ( u_dlat1 - m*u*tan(lat) - m*u/(2*tan(lat)) )/sqrt(sin(lat).abs) \
         + 2*( v_dlat1 - m*v/(2*tan(lat)) - m*v*tan(lat) )/sqrt(sin(lat).abs) \
         )
#u_dlat1 - u/(2*tan(lat)) - u*tan(lat) 

  
  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
#    sigDot[true,k] = - $UVel*( \
#                           ( lekmanMassFlux_divy - v0_interior_divy*(sig[k]))*( 1.0 - Math.sqrt(2)*exp(-xiB)*sin(xiB + PI/4.0) ) \
#                         - uekmanMassFlux_divy * exp(-xiT)*cos(xiT)  \
#                       )/RPlanet

    sigDot[true,k] = - $UVel*( \
                           ( lekmanMassFlux_divy + v_interior_divy*(1.0 + sig[k]))*( 1.0 - Math.sqrt(2)*exp(-xiB)*sin(xiB + PI/4.0) ) \
                           - (uekmanMassFlux_divy)* exp(-xiT)*cos(xiT)  \
                       )/RPlanet

  end

  return sigDot
end


def to_gphys(narray, name, long_name, units, grid)
  dat_a = VArray.new(narray, \
                    {"long_name"=>long_name, "units"=>units}, name)
  return GPhys.new(grid, dat_a)
end

def gen_grid(ncfile="")
  lat = nil;  z = nil
  if ncfile.length == 0 then
    lat = NArray.float(NLAT).indgen
    z   = NArray.float(NZ).indgen
    lat = lat * 90.0*(PI/180.0)/(NLAT-1)
    z   = - z /(NZ-1)
  else
    lat = GPhys::IO.open_gturl("#{ncfile}@lat").val / 180.0*PI
    z   = GPhys::IO.open_gturl("#{ncfile}@sig").val
  end

  lat_a = VArray.new(lat*180/PI, {"long_name"=>"latitude", "units"=>"degrees_east"}, "lat")
  sig_a = VArray.new(z, {"long_name"=>"Sigma", "units"=>"(1)"}, "sig")
  y_grid = Grid.new(Axis.new.set_pos(lat_a))
  yz_grid = Grid.new(Axis.new.set_pos(lat_a), Axis.new.set_pos(sig_a))

  return lat, z, y_grid, yz_grid
end

lat, z, y_grid, yz_grid = gen_grid(GRIDNCFILE)
set_parameters(lat)

$f  = 2.0*Omega*sin(lat)
he = to_gphys(sqrt(2*Av/$f), "he", "thickness of ekmanlayer", "m", y_grid)
spinupTime = to_gphys(Depth/sqrt(2*Av*$f), "spinupTime", "spinup time due to ekman pumping", "s", y_grid)
taux = to_gphys(gen_taux(lat), "taux", "zonal component of surface stress", "Pa", y_grid)
u = to_gphys(calc_u(lat,z), "U", "zonal flow", "m*s-1", yz_grid)
v = to_gphys(calc_v(lat,z), "V", "meridional flow", "m*s-1", yz_grid)
sigDot = to_gphys(calc_sigDot(lat,z), "SigDot", "vertical velocity", "(1)", yz_grid)


file = NetCDF.create("ekmanChecker_y#{GENNCFILE_SUFFIX}.nc")
GPhys::NetCDF_IO.write(file, he)
GPhys::NetCDF_IO.write(file, taux)
file.close

file = NetCDF.create("ekmanChecker_yz#{GENNCFILE_SUFFIX}.nc")
GPhys::NetCDF_IO.write(file, u)
GPhys::NetCDF_IO.write(file, v)
GPhys::NetCDF_IO.write(file, sigDot)

def get_diffWithNumSol(gturl, anasol)
  numsol = GPhys::IO.open_gturl(gturl)
  solDiff = numsol - anasol
  solDiff.name = "#{numsol.name}_diff"

  l2error = (solDiff.cut('lat'=>10..90).cut('sig'=>-0.1..-0.9)**2).mean.sqrt
  puts "#{solDiff.name} L2 error norm:  #{l2error}"
  return solDiff
end

if EXPDATADIR.length > 0 then
  GPhys::NetCDF_IO.write(file, get_diffWithNumSol("#{EXPDATADIR}/#{RUNSPECIFIC}_U.nc@U,lon=0,t=#{TIME}", u))
  GPhys::NetCDF_IO.write(file, get_diffWithNumSol("#{EXPDATADIR}/#{RUNSPECIFIC}_V.nc@V,lon=0,t=#{TIME}", v))
  GPhys::NetCDF_IO.write(file, get_diffWithNumSol("#{EXPDATADIR}/#{RUNSPECIFIC}_SigDot.nc@SigDot,lon=0,t=#{TIME}", sigDot))
end
file.close

