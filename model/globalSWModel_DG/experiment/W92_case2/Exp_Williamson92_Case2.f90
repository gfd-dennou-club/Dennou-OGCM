module Exp_Williamson94_Case2

  use dc_message

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod
  use PolyMesh_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  use OutputData_mod, only: &
       & OutputDataAnalysisInfo

  use DGHelper_mod
  use DGCalcusUtil_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case2, callBack_EndCurrentTimeStep

  real(DP), parameter :: SolidRotate_AngVel = 2d0*PI / (3600d0*24d0*12d0)

  real(DP), allocatable :: wc_hAnalystic(:,:)

contains
subroutine setIniCond_ExpWS94Case2()

  use DGElement_mod

  real(DP) :: u0
  real(DP), parameter :: alpha = 0*PI/2d0
  type(Vector3d) :: geo_vel, cart_vel, geoPos
  real(DP) :: r
  type(Vector3d) :: centerPos, b_1, b_2, b1, b2
  
  integer :: nc, nk
  real(DP), allocatable :: wc(:,:)
  type(Triangle) :: tri

  u0 = SolidRotate_AngVel*Radius

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
             & *(sin(geoPos%v_(2))*cos(alpha)-cos(geoPos%v_(1))*cos(geoPos%v_(2))*sin(alpha))**2
 
    end do
  end do

  call DGHelper_MallocElemNode(wc_hAnalystic)
  wc_hAnalystic = wc_h

end subroutine setIniCond_ExpWS94Case2

subroutine callBack_EndCurrentTimeStep(tstep, wc_h, wc_hU1, wc_hU2)

  use OutputData_mod, only: &
       & lonlat_interpolate_wc, Output_FieldData
  use gtool_history

  integer, intent(in) :: tstep
  real(DP), dimension(:,:), intent(in) :: wc_h, wc_hU1, wc_hU2

  integer :: nc, nk
  real(DP) :: CurrentTime
  type(vector3d) :: centerPos, geo_pos
  real(DP), allocatable :: wc_hError(:,:)
  real(DP) :: l2norm, linfnorm

  CurrentTime = tstep*delTime
  if( mod(tStep*delTime, outputIntrVal) == 0 ) then

     call DGHelper_MallocElemNode(wc_hError)

     wc_hError = abs(wc_h - wc_hAnalystic)

     l2norm = integrate_over_globalRigion(wc_hError)/abs(integrate_over_globalRigion(wc_hAnalystic))
     linfnorm = maxval(abs(wc_hError))/maxval(abs(wc_hAnalystic))

     call Output_FieldData('hError', wc_hError)
     call HistoryPut('l2norm', l2norm)
     call HistoryPut('linfnorm', linfnorm)
     write(*,*) "Error:", &
          & "l2:", l2norm, ", linf:", linfnorm

  end if

end subroutine callBack_EndCurrentTimeStep


end module Exp_Williamson94_Case2
