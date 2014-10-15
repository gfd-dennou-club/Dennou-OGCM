module Exp_Williamson94_Case5

  use dc_message
  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod
  use PolyMesh_mod
  
  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  use DGElement_mod
  use DGHelper_mod
  use DGCalcusUtil_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case5, callBack_EndCurrentTimeStep

  real(DP), parameter :: SolidRotate_AngVel = 2d0*PI / (3600d0*24d0*12d0)
  real(DP), parameter :: cosBellCenterLon = 270d0 * PI/180d0 
  real(DP), parameter :: cosBellCenterLat =   30d0 * PI/180d0
  real(DP), parameter :: cosBellRadius = radius/3d0
  real(DP), parameter :: h0 = 2000d0

  real(DP) :: initTE, initMass
contains
subroutine setIniCond_ExpWS94Case5()
  use Exp_Williamson94_Case1, only: &
       & create_CosBellField
  
  real(DP), parameter :: alpha = 0*PI/2d0
  type(Vector3d) :: geo_vel, cart_vel, geoPos
  real(DP) :: r
  type(Vector3d) :: centerPos, b_1, b_2, b1, b2
  
  integer :: nc, nk
  real(DP), allocatable :: wc(:,:)
  type(Triangle) :: tri
  real(DP) :: u0

  centerPos = SphToCartPos(cosBellCenterLon, cosBellCenterLat, radius)
  wc_BtmTopl = create_CosBellField(h0, centerPos, cosBellRadius)
  u0 = 20d0!SolidRotate_AngVel*Radius

  do nc=1, nDGElement
     do nk=1, nDGNodePerElem

        b_1 = get_DGElemCovariantBasis(1,DGElemInfo%node(nk),nc)
        b_2 = get_DGElemCovariantBasis(2,DGElemInfo%node(nk),nc)
        call get_DGElemContravariantBasis(b1, b2, b_1, b_2)

        geoPos = CartToSphPos(wc_DGNodePos(nk,nc))
        geo_vel%v_(1) = u0*(cos(geoPos%v_(2))*cos(alpha)+cos(geoPos%v_(1))*sin(geoPos%v_(2))*sin(alpha))
        geo_vel%v_(2) = -u0*sin(geoPos%v_(1))*sin(alpha)
        geo_Vel%v_(3) = 0d0

        
        cart_vel = SphToCartVec(geo_vel, wc_DGNodePos(nk,nc))
        wc_U1(nk,nc) = cart_vel.dot.b1
        wc_U2(nk,nc) = cart_vel.dot.b2
        wc_h(nk,nc) = - (Radius*Omega*u0 + 0.5d0*u0**2)/Grav &
             & *(sin(geoPos%v_(2))*cos(alpha)-cos(geoPos%v_(1))*cos(geoPos%v_(2))*sin(alpha))**2 &
             & -wc_BtmTopl(nk,nc)
 
    end do
  end do

end subroutine setIniCond_ExpWS94Case5

subroutine callBack_EndCurrentTimeStep(tstep, wc_h, wc_hU1, wc_hU2)
  use OutputData_mod, only:lonlat_interpolate_wc
  use gtool_history
  integer, intent(in) :: tstep
  real(DP), dimension(:,:), intent(in) :: wc_h, wc_hU1, wc_hU2

  integer :: nc, nk
  real(DP) :: CurrentTime
  type(vector3d), dimension(nDGNodePerElem) :: b_1, b_2
  real(DP), dimension(:,:), allocatable :: wc_TE, wc_Mass
  real(DP) :: TE, Mass

  
  CurrentTime = tstep*delTime
  if( mod(tStep*delTime, outputIntrVal) == 0 ) then

     call DGHelper_MallocElemNode(wc_Mass)
     call DGHelper_MallocElemNode(wc_TE)

     wc_Mass = wc_h + meanDepth
     
     !$omp parallel do private(b_1, b_2)
     do nc=1, nDGElement
        b_1 = get_DGElemCovariantBasisAtNodes(1,nc)
        b_2 = get_DGElemCovariantBasisAtNodes(2,nc)
        do nk=1, nDGNodePerElem
           wc_TE(nk,nc) = 0.5d0*l2norm(wc_hU1(nk,nc)*b_1(nk) + wc_hU2(nk,nc)*b_2(nk))**2/wc_Mass(nk,nc) &
                & + Grav*wc_Mass(nk,nc)*(0.5d0*wc_Mass(nk,nc) + wc_BtmTopl(nk,nc))
        end do
     end do

     Mass = integrate_over_globalRigion(wc_Mass)
     TE = integrate_over_globalRigion(wc_TE)

     if(tStep==0) then
        initMass = Mass; initTE = TE
        write(*,*) 'init Mass=', initMass, 'init Energy=', initTE
     end if

     write(*,*) 'MassVariation=', (Mass - initMass)/initMass, &
          & ', EnergyVariation=', (TE - initTE)/initTE
     call HistoryPut('MassVariation', abs(Mass - initMass)/initMass)
     call HistoryPut('TEVariation', abs(TE - initTE)/initTE)

  end if

end subroutine callBack_EndCurrentTimeStep

end module Exp_Williamson94_Case5
