module Exp_Williamson94_Case1

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod
  use PolyMesh_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  use DGHelper_mod
  use DGCalcusUtil_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case1, callBack_EndCurrentTimeStep
  public :: create_CosBellField

  real(DP), parameter :: SolidRotate_AngVel = 2d0*PI / (3600d0*24d0*12d0)

  real(DP), parameter :: cosBellCenterLon = 270d0 * PI/180d0 
  real(DP), parameter :: cosBellCenterLat =   0d0 * PI/180d0
  real(DP), parameter :: cosBellRadius = radius/3d0
  real(DP), parameter :: h0 = 1000d0
  
contains
subroutine setIniCond_ExpWS94Case1()
 

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)*0
  type(Vector3d) :: geo_vel, cart_vel, geoPos
  real(DP) :: r
  type(Vector3d) :: centerPos, b_1, b_2, b1, b2
  
  integer :: nc, nk
  real(DP), allocatable :: wc(:,:)
  
  centerPos = SphToCartPos(cosBellCenterLon, cosBellCenterLat, radius)
  wc_h = create_CosBellField(h0, centerPos, cosBellRadius)

  u0 = SolidRotate_AngVel*Radius

  do nc=1, nDGElement
     do nk=1, nDGNodePerElem

        b_1 = get_DGElemCovariantBasis(1,DGElemInfo%node(nk),nc)
        b_2 = get_DGElemCovariantBasis(2,DGElemInfo%node(nk),nc)
        call get_DGElemContravariantBasis(b1, b2, b_1, b_2)

        geoPos = CartToSphPos(wc_DGNodePos(nk,nc))
        geo_vel = (/ &
             & u0*(cos(geoPos%v_(2))*cos(alpha)+cos(geoPos%v_(1))*sin(geoPos%v_(2))*sin(alpha)), &
             & -u0*sin(geoPos%v_(1))*sin(alpha), 0d0 /)
        cart_vel = SphToCartVec(geo_vel, wc_DGNodePos(nk,nc))
        wc_U1(nk,nc) = cart_vel.dot.b1
        wc_U2(nk,nc) = cart_vel.dot.b2
     end do
  end do

  call DGHelper_MallocElemNode(wc)
  wc = 1d0
  write(*,*) "Global Integrate:", integrate_over_globalRigion(wc), 4d0*PI*Radius**2, &
       & (integrate_over_globalRigion(wc) - 4d0*PI*radius**2)/(4d0*PI*radius**2)
  
!stop
end subroutine setIniCond_ExpWS94Case1

function create_CosBellField(h0, cosBellCenterPos, cosBellRadius) result(wc_CosBell)

  real(DP), intent(in) :: h0
  type(vector3d), intent(in) :: cosBellCenterPos
  real(DP), intent(in) :: cosBellRadius
  real(DP) :: wc_CosBell(nDGNodePerElem, nDGElement)

  integer :: nc, nk
  real(DP) :: r

  wc_CosBell = 0d0

  !$omp parallel do private(nk,r)
  do nc=1, nDGElement
     do nk=1, nDGNodePerElem
        
        r = geodesicArcLength(cosBellCenterPos, wc_DGNodePos(nk,nc))
        if(r < cosBellRadius) then
           wc_CosBell(nk,nc) = 0.5d0*h0*(1d0 + cos(PI*r/cosBellRadius))
        end if

     end do
  end do

end function create_CosBellField

subroutine callBack_EndCurrentTimeStep(tstep, wc_h, wc_hU1, wc_hU2)
  use gtool_history, only: HistoryPut
  use OutputData_mod, only: Output_FieldData

  integer, intent(in) :: tstep
  real(DP), intent(in), dimension(:,:) :: wc_h, wc_hU1, wc_hU2

  integer :: nc, nk
  real(DP) :: CurrentTime
  real(DP), allocatable :: wc_hAnalystic(:,:), wc_hError(:,:)
  type(vector3d) :: centerPos, geo_pos
  real(DP) :: dLon
  real(DP) :: l2norm, linfnorm

  CurrentTime = tstep*delTime
  if( mod(tStep*delTime, outputIntrVal) == 0 ) then


     call DGHelper_MallocElemNode(wc_hAnalystic)
     call DGHelper_MallocElemNode(wc_hError)

     dLon = SolidRotate_AngVel*CurrentTime
     geo_pos = (/ mod(cosBellCenterLon + dLon, 2d0*PI), cosBellCenterLat, radius /)
     wc_hAnalystic = create_CosBellField(h0, SphToCartPos(geo_pos), cosBellRadius)

     wc_hError = abs(wc_h - wc_hAnalystic)
     l2norm = sqrt( integrate_over_globalRigion(wc_hError**2)       &
          &         /integrate_over_globalRigion(wc_hAnalystic**2) )

     linfnorm = maxval(abs(wc_hError))/maxval(abs(wc_hAnalystic))

     call Output_FieldData('hError', wc_hError)
     call HistoryPut('l2norm', l2norm)
     call HistoryPut('linfnorm', linfnorm)
     write(*,*) "h Error:", &
          & "l2:", l2norm, ", linf:", linfnorm
  end if

end subroutine callBack_EndCurrentTimeStep


end module Exp_Williamson94_Case1
