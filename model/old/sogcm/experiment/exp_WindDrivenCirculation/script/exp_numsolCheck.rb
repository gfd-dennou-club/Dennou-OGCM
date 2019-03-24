#!/usr/bin/env ruby

require "numru/gphys"
include NumRu
include NMath

Av = 1e-02
Ah = 1e3
HAh = 0e14
TAU0 = 0.12
Dens = 1.027e03
NLAT = 721
NZ   = 1001
PI = acos(-1.0)
Omega = 7.292115e-5
Depth = 5.2e03
RPlanet = 6371e03
TIME = 108000

SphericalDOptr=true

EXPSPECIFIC="Ah1e3T341L60"
RUNSPECIFIC="Run1"
GENNCFILE_SUFFIX = "_#{EXPSPECIFIC}"
#EXPDATADIR = "/home/ykawai/Dennou-OGCM/model/sogcm/exp/"
EXPDATADIR = "/home/ykawai/exp_backup/current/exp_#{EXPSPECIFIC}"
GRIDNCFILE = "#{EXPDATADIR}/#{RUNSPECIFIC}_U.nc"

##########################

def gen_taux(lat)
  coef = NArray.to_na( [
    0.0265682, -0.0784899, -0.00880389, 0.0343205, 0.0233334, 0.000641955, -0.00387676, -0.00150998
  ] )

  taux = NArray.float(lat.length)
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

  taux_dnLat = NArray.float(lat.length)
  taux_dnLat[true] = 0.0

  if n%2 == 0 then
    sign = (-1.0)**(n/2)

    for i in 0..coef.length-1
      taux_dnLat = taux_dnLat + sign * (2.0*i+1.0)**n * coef[i]*cos( (2.0*i+1.0)*lat )
    end
  else
    sign = (-1.0)**n
    for i in 0..coef.length-1
      taux_dnLat = taux_dnLat + sign * (2.0*i+1.0)**n * coef[i]*sin( (2.0*i+1.0)*lat )
    end
  end
#  p "gen_taux:", taux
  return taux_dnLat
end

def printInfo(str, dat, lat)
  
  latD = 180.0/PI*lat

  puts "== #{str}"
  for i in 0..latD.length-1
#      if ( (latD[i]).to_i%2.5 == 0 ) then
        l = latD[i].round; d = dat.val[i]
        puts "lat=#{l}: #{d}"
#      end
  end

end

def calc_interiorBarotFlow(lat, isNonDim=false, n_dlat=0)

  f = 2*Omega*sin(lat)
  ev = 2.0*Av / (f.abs*Depth**2)
  velScale = 2.0*TAU0 / (Dens*sqrt(2.0*Av*f.abs))
  reInv = Ah / (velScale*RPlanet)
  ro = velScale / (f.abs*RPlanet)
  r = sqrt(ev)/ro

  coef = 2*reInv/r
  u0 = NArray.float(lat.length)
  u0[true] = 0.0
  for n in 0..2
    u0 = u0 + coef**(n) * gen_taux_dnLat(lat,2.0*n+n_dlat)/TAU0
  end
  p coef

  if n_dlat==0 then
    u0 = u0 + coef*( \
                         -tan(lat)*gen_taux_dnLat(lat,1) - gen_taux(lat)/cos(lat)**2  \
                      )/TAU0
    u0 = u0 + coef**2*( \
                         -tan(lat)*gen_taux_dnLat(lat,2) - gen_taux_dnLat(lat,1)/cos(lat)**2  \
                      )/TAU0
  end

  return isNonDim ? u0 : velScale*u0
end

def get_velScale(lat)
  f = 2.0*Omega*sin(lat)
  
  velScale = 2.0*TAU0/(Dens*sqrt(2.0*Av*f.abs))
  velScale_y = - velScale/(2*tan(lat))
  velScale_yy = - 0.5*(velScale_y/(tan(lat)) - velScale/(sin(lat)**2))

end

def calc_interiorV0(lat, isNonDim=false, n_dlat=0)

  f = 2.0*Omega*sin(lat)
  eh0 = Ah/(f*RPlanet**2)
  eh4 = HAh/(f*RPlanet**4)
#  puts "eh0 = #{eh0.mean}"

  velScale = 2.0*TAU0/(Dens*sqrt(2.0*Av*f.abs))
  velScale_y = - velScale/(2*tan(lat))
  velScale_yy = - 0.5*(velScale_y/(tan(lat)) - velScale/(sin(lat)**2))
  velScale_yyy = - ( velScale_yy/tan(lat) - 2*velScale_y/sin(lat)**2 + 2.0*velScale*cos(lat)/sin(lat)**3 )
  
  taux = gen_taux_dnLat(lat, 0)
  taux_y = gen_taux_dnLat(lat, 1)
  taux_yy = gen_taux_dnLat(lat, 2)

  if (n_dlat==1) then
    v0_y = NArray.float(lat.length)
    v0_y[true] = - eh0*( \
                       velScale*gen_taux_dnLat(lat, 3)              \
                       + 0*2*velScale_y*gen_taux_dnLat(lat, 2) \
                       + 0*velScale_yy*taux_y                                    \
                       - tan(lat)*velScale*taux_yy      \
                       - 0*tan(lat)*velScale_y*taux_y      \
                       - 0*velScale*taux_y/cos(lat)**2   \
                       + 0*2*velScale*taux_y \
             )/(TAU0*velScale)
    
    v0_y[true] = v0_y[true] + 0*eh4*( \
                                velScale*gen_taux_dnLat(lat,5) \
                                )/(TAU0*velScale)

    return v0_y
  end

  v0 = NArray.float(lat.length)
  v0[true] = - eh0*( \
                     velScale*gen_taux_dnLat(lat, 2)              \
                     + 2.0*velScale_y*taux_y \
                     + velScale_yy*taux                                    \
                     - tan(lat)*velScale*taux_y      \
                     - tan(lat)*velScale_y*taux      \
                     - velScale*taux/cos(lat)**2   \
                     + 2*velScale*taux \
             )/(TAU0*velScale)

  v0[true] = v0[true] + eh4*( \
                 velScale*gen_taux_dnLat(lat,4) \
               - 0*tan(lat)*velScale*gen_taux_dnLat(lat,3) \
               - 0*velScale*gen_taux_dnLat(lat,2)/cos(lat)**2 \
               + 0*2*velScale*gen_taux_dnLat(lat,2) \
            )/(velScale*TAU0)

  return isNonDim ? v0 : velScale*v0

end

def calc_u(lat, sig)

  f = 2.0*Omega*sin(lat)
  deltaE = sqrt( 2.0*Av/(f.abs*Depth**2) )
  velScale = 2.0*TAU0/(Dens*sqrt(2.0*Av*f.abs))

  taux = gen_taux(lat)/TAU0
  u = NArray.float(lat.length, sig.length)
  barotU = calc_interiorBarotFlow(lat, true)

  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
    u[true,k] = velScale*( \
                  barotU[true]*(1.0 - exp(-xiB) * cos(xiB) ) \
                + taux/Math.sqrt(2) * exp(-xiT) * cos(xiT + PI/4.0) \
               )
  end

  return u
end

def calc_v(lat, sig)

  f = 2*Omega*sin(lat)
  deltaE = sqrt( 2.0*Av/(f.abs*Depth**2) )

  velScale = 2.0*TAU0/(Dens*sqrt(2.0*Av*f.abs))

  ro = velScale/(f.abs*RPlanet)

  taux = gen_taux(lat)/TAU0
  barotU = calc_interiorBarotFlow(lat, true)
  v0_interior = calc_interiorV0(lat, true)

  v = NArray.float(lat.length, sig.length)
  
  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
    v[true,k] = velScale*( \
                  barotU[true]*exp(-xiB) * sin(xiB) \
                + v0_interior \
                - taux/Math.sqrt(2) * exp(-xiT) * sin(xiT + PI/4.0)  \
               )
  end

  return v
end

def calc_sigDot(lat, sig)

  f = 2.0*Omega*sin(lat)
  deltaE = sqrt(2.0*Av/(f.abs*Depth**2))
  velScale = 2.0*TAU0/(Dens*f.abs*deltaE*Depth)
  velScale_y = - velScale/(tan(lat))


  taux = gen_taux(lat)/TAU0
  taux_y = gen_taux_dnLat(lat, 1)/TAU0

  barotU = calc_interiorBarotFlow(lat, true)  
  lekmanMassFlux_divy = 0.5*deltaE*( \
                         velScale*calc_interiorBarotFlow(lat, true, 1) + velScale_y*barotU \
                       - velScale*barotU*tan(lat)  )

  v0_interior = calc_interiorV0(lat,true)
  v0_interior_divy =   velScale*calc_interiorV0(lat, true, 1) \
                     + velScale_y*v0_interior  \
                     - velScale*v0_interior*tan(lat)

  uekmanMassFlux_divy   = 0.5*deltaE*( \
                                     velScale*taux_y + velScale_y*taux \
                                       - velScale*taux*tan(lat) \
                                     )


  sigDot = NArray.float(lat.length, sig.length)
  
  for k in 0..sig.length-1
    xiB = (1.0 + sig[k])/deltaE
    xiT = - sig[k]/deltaE
    sigDot[true,k] = - ( \
                           ( lekmanMassFlux_divy + 0*v0_interior_divy*(-sig[k]) )*( 1.0 - Math.sqrt(2)*exp(-xiB)*sin(xiB + PI/4.0) ) \
                         - uekmanMassFlux_divy * exp(-xiT)*cos(xiT)  \
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

f  = 2.0*Omega*sin(lat)
he = to_gphys(sqrt(2*Av/f), "he", "thickness of ekmanlayer", "m", y_grid)
#printInfo("ekman layer thickness scale", he, lat)

spinupTime = to_gphys(Depth/sqrt(2*Av*f), "spinupTime", "spinup time due to ekman pumping", "s", y_grid)
#printInfo("spin up time scale [day]", spinupTime/86400.0, lat)

taux = to_gphys(gen_taux(lat), "taux", "zonal component of surface stress", "Pa", y_grid)

ssflowU = to_gphys(he.val*taux.val/(2.0*Dens*Av), "ssflowU", "zonal velocity at surface", "m*s-1", y_grid)
ssflowV = to_gphys(-he.val*taux.val/(2.0*Dens*Av), "ssflowV", "meriodinal velocity at surface", "m*s-1", y_grid)
#printInfo("zonal surface flow", ssflowU, lat)

barotflowU = to_gphys(2.0*taux.val/(Dens*f*he.val), "barotflowU", "barotropic zonal flow", "m*s-1", y_grid)
barotflowU_hVisc = to_gphys(calc_interiorBarotFlow(lat), "barotflowU_hVisc", "barotropic zonal flow with hVisc", "m*s-1", y_grid)
#printInfo("zonal barotropic flow", barotflowU, lat)

u = to_gphys(calc_u(lat,z), "U", "zonal flow", "m*s-1", yz_grid)
v = to_gphys(calc_v(lat,z), "V", "meridional flow", "m*s-1", yz_grid)
sigDot = to_gphys(calc_sigDot(lat,z), "SigDot", "vertical velocity", "(1)", yz_grid)

ro = to_gphys( (2*TAU0/(Dens*sqrt(2*Av*f.abs)))/(f.abs*RPlanet), "Ro", "Rossby numnber for barotropic zonal flow", "(1)", y_grid)
re_inv = to_gphys((Ah*Dens*sqrt(2.0*Av*f.abs) / (2.0*TAU0*RPlanet)), "ReInv", "Inverse of Reynolse number for barotropic zonal flow", "(1)", y_grid)


file = NetCDF.create("ekmanChecker_y#{GENNCFILE_SUFFIX}.nc")
GPhys::NetCDF_IO.write(file, he)
GPhys::NetCDF_IO.write(file, taux)
GPhys::NetCDF_IO.write(file, ssflowU)
GPhys::NetCDF_IO.write(file, ssflowV)
GPhys::NetCDF_IO.write(file, barotflowU)
GPhys::NetCDF_IO.write(file, barotflowU_hVisc)
GPhys::NetCDF_IO.write(file, ro)
GPhys::NetCDF_IO.write(file, re_inv)
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

###########################

NLAT2=145
lat = NArray.sfloat(NLAT2).indgen
lat = lat * 180.0/(NLAT2-1) - 90

#
taux = gen_taux(lat*PI/180.0)
#printInfo("taux", taux, lat*PI/180.0)
lat_a = VArray.new(lat, 
   {"long_name"=>"latitude", "units"=>"degrees_east"}, "lat")
taux_a = VArray.new(taux, 
                    {"long_name"=>"zonal component of wind stress", "unit"=>"N m-2"}, "windStressLon")

# 
glat = Axis.new.set_pos(lat_a)
gtaux = GPhys.new(Grid.new(glat), taux_a)
file = NetCDF.create("windStressLon.nc")
GPhys::NetCDF_IO.write(file, gtaux)
file.close
