#!/usr/env/ruby
# -*- coding: utf-8 -*-

require "numru/ggraph"
include NumRu
include NMath
include Math

## Setting #############################

Ah = 8e02
Av = 1e-02
Dens = 1.024e03
PI = acos(-1.0)
Omega = 7.292e-05
RPlanet = 6371e03

#
H = 5.2e03
LATMIN = 20.0
LATMAX = 47.5
LAT0 = 0.5*(LATMAX + LATMIN) * PI/180.0
L = 0.5*(LATMAX-LATMIN)*PI/180.0*RPlanet
#LATMIN = 46.5
#LATMAX = 90.0
#LAT0  = 90.0*PI/180.0 
#L = RPlanet
#L = (LATMAX-LATMIN)*PI/180.0*RPlanet


NY_AsympSol = 201
NZ_AsympSol = 201

# Polynominal coeffecients for x-component of nondimensional wind stress
# tau_x = T0 + T1*y + T2*y^2 + T3*y^3 + T4*y^4

Y0=-0.272; a = 1.0; b = -0.5;
#Y0=-0.0; a = 1.0; b = 1.0;

TCoef1 = 0.75*(a-b)
TCoef3 = - 0.25*(a-b)
TCoef4 = (-Y0*(Y0-3.0)*TCoef3 - 0.5*(a+b))/(Y0**2 - 1)**2
TCoef2 =  -2.0*TCoef4
TCoef0 = TCoef4 + (a + b)/2.0
Tau0 = 0.124

CreateFigFlag=true
FigFilePrefix="AsymSol_HCompari_"

NumSolDataParentDir="/home/ykawai/exp_backup/current/"
#NumSolExpNameList = ["Kh800T85L40","Kh800T85L60","Kh800T85L80"]
NumSolExpNameList = ["Kh800HAhT42L60","Kh800HAhT85L60","Kh800HAhT170L60"]
#NumSolExpNameList = ["Kh800T42L60","Kh800T85L60","Kh800T170L60"]

NumSolDataPrefix = "Run1_"
NumSolDataLastPeriod = 108000 #[day]

##########################################

NCompariNumSol = NumSolExpNameList.length
NumSolDataDir = Array.new(NCompariNumSol)
NumSolExpNameList.each_with_index{|expName,i|
  NumSolDataDir[i] = "#{NumSolDataParentDir}exp_#{expName}/"
}

PI = acos(-1.0)
##########################################

CurlTauCoefs = Struct.new("CurlTauCoefs", :a0, :a1, :a2, :a3)

#
# Psi = C1*sinh(Re^1/2 y) + C2*cosh(Re^1/2 y) + C3 + C4 y
PartSolCoefs_Psi = Struct.new("PsiSolCoefs_Psi", :c1, :c2, :c3, :c4)

# Calculate scaling coeffecients or nondimensional parameters
f0 = 2*Omega*sin(LAT0)
Ev = 2*Av/(f0*H**2)
Eh = 2*Ah/(f0*L**2)

U = 2*Tau0/(Dens*f0*H*Math.sqrt(Ev))
W = U*H/L
Ro = U/(f0*L)
Re = U*L/Ah

def puts_parameterInfo()

  puts "-- Scaling --"
  puts "L=#{L}[m], H=#{H}[m]"
  puts "U=#{U}[m/s]"
  puts "W=#{W}[m/s]"
  puts "LAT0=#{LAT0*180/PI}[deg]"

  puts "-- Nondimensional Parameters--"
  puts "Ro=#{Ro}"
  puts "Eh=#{Eh}"
  puts "Ev=#{Ev}"
  puts "Re=#{Re}"

  puts "-- characteristic quantity  --"
  puts " Ekman later depth = #{H*Math.sqrt(Ev)} [m] "
end

def create_grid()
  y = NArray.float(NY_AsympSol)
  z = NArray.float(NZ_AsympSol)

  for j in 0..NY_AsympSol-1
      y[j] = -1.0 + 2.0/(NY_AsympSol-1)*j 
  end
  for k in 0..NZ_AsympSol-1
      z[k] = 0.0 + 1.0/(NZ_AsympSol-1)*k
  end

  y_axis = Axis.new.set_pos(VArray.new(y, {"long_name"=>"y", "units"=>"1"}, "y"))
  z_axis = Axis.new.set_pos(VArray.new(z, {"long_name"=>"z", "units"=>"1"}, "z"))

  return Grid.new(y_axis, z_axis)
end

def create_nondimVar(grid, name, long_name, units)
  data = VArray.new(NArray.float(NY_AsympSol,NZ_AsympSol), {"long_name"=>long_name,"units"=>units},name)
  return GPhys.new(grid, data)
end

def create_yz_gphys(grid)
  y = create_nondimVar(grid, "y", "nondimensional y coordinate", "1")
  z = create_nondimVar(grid, "z", "nondimensional z coordinate", "1")

  for k in 0..NZ_AsympSol-1
    y[true,k] = grid.axis("y").pos[true]
  end
  for j in 0..NY_AsympSol-1
    z[j,true] = grid.axis("z").pos[true]
  end

  return y, z
end

def create_taux_gphys(y, t0, t1, t2, t3, t4, grid)
  taux = create_nondimVar(grid, "taux", "x-component of nondimensional surface stress", "(1)")
  
  taux[true,true] = t0 + y[true,true]*(t1 + y[true,true]*(t2 + y[true,true]*(t3 + y[true,true]*t4)))

  p taux
  return taux
end

def create_curl_tau_gphys(y, t0, t1, t2, t3, t4, grid)
  curl_tau = create_nondimVar(grid, "curl_tau", "z-component of rotation for nondimensional surface stress", "(1)")
  
  curl_tau[true,true] = t1 + y[true,true]*(2*t2 + y[true,true]*(3*t3 + 4*y[true,true]*t4))

  return curl_tau
end

################################################
def get_PsiSol_InhomoTerm(y, curlTauInfo)
  a0 = curlTauInfo.a0; a1 = curlTauInfo.a1; a2 = curlTauInfo.a2; a3 = curlTauInfo.a3;
  return - y**2*((a0/2 +a2/Re) + y*((a1/6.0 + a3/Re) + y*(a2/12.0 + y*a3/20.0)))
end

def get_PsiSol_InhomoTerm_dy(y, curlTauInfo)
  a0 = curlTauInfo.a0; a1 = curlTauInfo.a1; a2 = curlTauInfo.a2; a3 = curlTauInfo.a3;
  return - y*(2.0*(a0/2 +a2/Re) + y*(3.0*(a1/6.0 + a3/Re) + y*(4.0*a2/12.0 + 5.0*y*a3/20.0)))
end

def get_PsiSol_InhomoTerm_dyy(y, curlTauInfo)
  a0 = curlTauInfo.a0; a1 = curlTauInfo.a1; a2 = curlTauInfo.a2; a3 = curlTauInfo.a3;
  return - (2.0*(a0/2 +a2/Re) + y*(6.0*(a1/6.0 + a3/Re) + y*(12.0*a2/12.0 + 20.0*y*a3/20.0)))
end

def get_PsiSol_InhomoTerm_dyyy(y, curlTauInfo)
  a0 = curlTauInfo.a0; a1 = curlTauInfo.a1; a2 = curlTauInfo.a2; a3 = curlTauInfo.a3;
  return - (6.0*(a1/6.0 + a3/Re) + y*(24.0*a2/12.0 + 60.0*y*a3/20.0))
end

def get_PsiSol_InhomoTerm_dyyyy(y, curlTauInfo)
  a0 = curlTauInfo.a0; a1 = curlTauInfo.a1; a2 = curlTauInfo.a2; a3 = curlTauInfo.a3;
  return - (24.0*a2/12.0 + 120.0*y*a3/20.0)
end

def get_PsiSol_ConstCoefs(curlTauInfo)
  r=sqrt(Re)
  q_y_Y0 = get_PsiSol_InhomoTerm_dy(Y0,curlTauInfo)
  q_yy_p1 = get_PsiSol_InhomoTerm_dyy(1.0,curlTauInfo)
  q_yy_m1 = get_PsiSol_InhomoTerm_dyy(-1.0,curlTauInfo)

  c1 = - (q_yy_p1 - q_yy_m1)/(2*Re*sinh(r))
  c2 = - (q_yy_p1 + q_yy_m1)/(2*Re*cosh(r))
  c3 = - c2
  c4 = - q_y_Y0 - r*cosh(r*Y0)*c1 - r*sinh(r*Y0)*c2
  
  return c1,c2,c3,c4
end

##########################################

def calc_U(y, z, taux, curlTauInfo, grid)
  include Math
  u = create_nondimVar(grid, "U", "nondimensional zonal velocity", "1")

  r = sqrt(Re); sEv = sqrt(Ev)
  c1,c2,c3,c4 = get_PsiSol_ConstCoefs(curlTauInfo)
  q_y = get_PsiSol_InhomoTerm_dy(y, curlTauInfo)

  u0 = - (r*(c1*(r*y).cosh + c2*(r*y).sinh) + c4 + q_y)

  xi1 = z/sEv
  xi2 = (1-z)/sEv
  u[true,true] = u0[true,true] * ( 1.0 - (-xi1).exp * (xi1).cos  ) \
                 + taux/sqrt(2) * ((-xi2).exp * (xi2 + 0.25*PI).cos)

  return u
end

def calc_V(y, z, taux, curlTauInfo, grid)
  include Math
  v = create_nondimVar(grid, "V", "nondimensional meriodinal velocity", "1")

  r = sqrt(Re); sEv = sqrt(Ev)
  c1,c2,c3,c4 = get_PsiSol_ConstCoefs(curlTauInfo)
  q_y = get_PsiSol_InhomoTerm_dy(y, curlTauInfo)
  q_yyy = get_PsiSol_InhomoTerm_dyyy(y, curlTauInfo)

  u0 = - (r*(c1*(r*y).cosh + c2*(r*y).sinh) + c4 + q_y)  

  # v1 = - 1/Re * d^2/dy^2 u0 
  v1 =  (1/Re) * \
             (r**3*(c1*(r*y).cosh + c2*(r*y).sinh) + q_yyy)

  xi1 = z/sEv
  xi2 = (1-z)/sEv
  v[true,true] = - u0[true,true] * (-xi1).exp * (xi1).sin  \
                 - taux/sqrt(2) * ((-xi2).exp * (xi2 + 0.25*PI).sin) \
                 + Ro*v1[true,true] 
  return v
end

def calc_W(y, z, taux, curlTauInfo, grid)
  include Math
  w = create_nondimVar(grid, "W", "nondimensional vertical velocity", "1")


  r = sqrt(Re); sEv = sqrt(Ev)
  c1,c2,c3,c4 = get_PsiSol_ConstCoefs(curlTauInfo)
  q_yy = get_PsiSol_InhomoTerm_dyy(y, curlTauInfo)
  q_yyyy = get_PsiSol_InhomoTerm_dyyyy(y, curlTauInfo)

  # + v1 = 1/Re * u0_yy 
  v1_y =  (1/Re) * (r**4*(c1*(r*y).sinh + c2*(r*y).cosh) + q_yyyy)
  int_v1_y_dz =  z*v1_y
  zeta0 = (r**2*(c1*(r*y).sinh + c2*(r*y).cosh) + q_yy)  

  xi1 = z/sEv
  xi2 = (1-z)/sEv
  w[true,true] = ( 0.5*sEv*zeta0[true,true] * ( 1.0 - sqrt(2.0)*(-xi1).exp * (xi1 + 0.25*PI).sin ) \
                   - Ro*int_v1_y_dz[true,true] \
                 )*( 1.0 - (-xi2).exp * (xi2).cos )

  return w
end

###################

BASE_VIEWPORT= [0.10,0.90,0.15,0.85]

def view_prepair()

  nRow = (NCompariNumSol + 1)/2

  DCL.swpset('iwidth', 1000)
  DCL.swpset('iheight',450*nRow)

  DCL.gropn(4)
  DCL.swpset("LWAIT", true)
  DCL.sgpset('lfprop', true)
  DCL.sgpset('lcntl', false)
  DCL.uzfact(0.6)  
#  DCL.slrat( 1.0, 0.9)

  DCL.sldiv("y", 2, nRow) 
end

def view_finish()
  DCL.grcls
end

def view_tonecont_graph(gphys, newframe=true, strOptHashes=nil)
  GGraph.tone_and_contour(gphys, newframe, strOptHashes)
  GGraph.color_bar if newframe==true

end

def view_line_graph(gphys, newframe=true, strOptHashes=nil)
  GGraph.line(gphys, true, strOptHashes)
end

def view_comparigraph(ary_gphys, graphType, strOptHashes=nil)

  vxmin = BASE_VIEWPORT[0]
  vxmax = BASE_VIEWPORT[1]
  vymin = BASE_VIEWPORT[2]
  vymax = BASE_VIEWPORT[3]

  nRestFrame = 2*((NCompariNumSol+1)/2)
  ary_gphys.each{|sol|
    GGraph.set_fig({"viewport"=> [ vxmin, vxmax, vymin,vymax ]})
    funcCode = "view_#{graphType}_graph(sol,true,strOptHashes)" 
    puts funcCode
    eval(funcCode)
    nRestFrame -= 1
  }
  
  puts "Number of rest frames is #{nRestFrame}."
  for i in 1..nRestFrame
    DCL.grfrm
    puts "Call DCL.grfrm"
  end

  if CreateFigFlag then
    pngName = "asymSol_#{FigFilePrefix}#{ary_gphys[0].name}.png"
    p "`mv dcl_*.png #{pngName}` .."
    `mv dcl_*.png #{pngName}`
  end


end

####################

def get_NumSol_gphys(varName, expID=0, strCutPos=nil, ncpath="#{NumSolDataDir[expID]}/#{NumSolDataPrefix}#{varName}.nc")

  gturl = "#{ncpath}@#{varName}"
  gturl << ",#{strCutPos}" if strCutPos
  puts "Load data from #{gturl}.." 
  return GPhys::IO.open_gturl(gturl)
end

def get_expsNumSol(varName, cutpos, nondimFactor)
  numsols = []
  for expID in 0..NCompariNumSol-1
    numsols[expID] = get_NumSol_gphys(varName, expID, cutpos)/nondimFactor
    numsols[expID].name = "#{varName}(#{NumSolExpNameList[expID]})"
    numsols[expID].long_name = "#{varName}(#{NumSolExpNameList[expID]})"
  end

  return numsols
end

##########################

puts_parameterInfo()

yz_grid = create_grid()
#p yz_grid.axis("y")
#p yz_grid.axis("z")
y, z = create_yz_gphys(yz_grid)

cutpos = "lat=#{LATMIN}:#{LATMAX},lon=0,t=#{NumSolDataLastPeriod}"
taux_nondim = create_taux_gphys(y, TCoef0,TCoef1,TCoef2,TCoef3, TCoef4, yz_grid)
curltau_nondim = create_curl_tau_gphys(y, TCoef0,TCoef1,TCoef2,TCoef3, TCoef4, yz_grid)

taux_numsol = get_NumSol_gphys("windStressLon", 0, strCutPos="lat=#{LATMIN}:#{LATMAX}", ncpath="./windStressLon.nc")/Tau0
# -d(tau_x)/dy = a0 + a1*y + a2*y^2
curlTauInfo = CurlTauCoefs.new(TCoef1, 2*TCoef2, 3*TCoef3, 4*TCoef4)

p taux_nondim


u_nondim = calc_U(y, z, taux_nondim, curlTauInfo, yz_grid)
u_numsols = get_expsNumSol("U", cutpos, U)
p u_numsols

v_nondim = calc_V(y, z, taux_nondim, curlTauInfo, yz_grid)
v_numsols = get_expsNumSol("V", cutpos, U)
p v_numsols

w_nondim = calc_W(y, z, taux_nondim, curlTauInfo, yz_grid)
sigDot_numsols = get_expsNumSol("SigDot", cutpos, U/L)
p sigDot_numsols

#
view_prepair()

view_comparigraph([taux_nondim, taux_numsol], "line")
view_comparigraph([taux_nondim, curltau_nondim], "line")
view_comparigraph([u_nondim, u_numsols].flatten, "tonecont", {'max'=>1.5, 'min'=>-1.5})
view_comparigraph([v_nondim, v_numsols].flatten, "tonecont", {'max'=>0.25, 'min'=>-0.25})
view_comparigraph([w_nondim, sigDot_numsols].flatten, "tonecont", {'max'=>1.5*sqrt(Ev), 'min'=>-1.5*sqrt(Ev)})

view_finish
