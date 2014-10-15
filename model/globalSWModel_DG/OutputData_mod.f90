module OutputData_mod

  use dc_types
  use dc_message
  use gtool_history

  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use netcdfDataWriter_mod

  use SimParameters_mod

  use GridSet_mod

  use DGHelper_mod
  use VariableSet_mod, only: &
       & wc_h, wc_hU1, wc_hU2, wc_dhdt, wc_U1, wc_U2, wc_etc, &
       & wc_BtmTopl

  implicit none
  
  interface Output_FieldData
     module procedure Output_FieldData_i
     module procedure Output_FieldData_d
  end interface Output_FieldData

  public :: OutputData, OutputDataAnalysisInfo
  public :: Output_FieldData, lonlat_interpolate_wc

  type(netcdfDataWriter) :: ncWriter
  type(PointScalarField) :: p_zeta

  real(DP), parameter :: lon1 = -90d0 * PI/180d0
  real(DP), parameter :: lon2 = 90d0 * PI/180d0
  real(DP), parameter :: lat1 = -89d0 * PI/180d0
  real(DP), parameter :: lat2 = 89d0 * PI/180d0

  real(DP), dimension(:), allocatable :: y_Lat, x_Lon
  integer :: iMax, jMax
  integer, dimension(:,:), allocatable :: lonlatDGElemId
  type(vector2d), dimension(:,:), allocatable :: lonlatPt2LclMatrix(:,:)
contains

subroutine OutputData_Init()

  !
  !
  integer :: i, j, nc
  type(Vector3d) :: cart_Pos, geo_Pos

  !
  !
  iMax = (lon2 - lon1) / latlonIntrv + 1
  jMax = (lat2 - lat1) / latlonIntrv + 1
  allocate(x_Lon(0:iMax-1), y_Lat(jMax))
  allocate(lonlatDGElemId(0:iMax-1,jMax), lonlatPt2LclMatrix(0:iMax-1,jMax))

  do i=1, iMax
     x_Lon(i-1) = lon1 + (i-1)*latlonIntrv! + 1d-2
  end do
  x_Lon(0) = x_Lon(0) + 1d-5
  x_Lon(iMax-1) = x_Lon(iMax-1) - 1d-5
  
  do j=1, jMax
     y_Lat(j) = lat1 + (j-1)*latlonIntrv !+ 1d-2
  end do
  y_Lat(1) = y_Lat(1) + 1d-9
  y_Lat(jMax) = y_Lat(jMax) - 1d-9

write(*,*) "Satrt to search.."  
lonlatDGElemId = -1

  !$omp parallel do private(i,geo_Pos,cart_Pos,nc)
  do j=1, jMax
     do i=0,iMax-1
        geo_Pos = (/ x_Lon(i), y_Lat(j), Radius /)
        cart_Pos = SphToCartPos(geo_Pos)
        do nc=1, nDGElement
           if( isPtContainedTri(cart_Pos, nc) ) then
              lonlatDGElemId(i,j) = nc
               if(.not.c_SeaCellFlag(nc)) then
                  write(*,*) "Not Sea..", x_Lon(i), y_Lat(j)
               end if
               lonlatPt2LclMatrix(i,j) = get_LocalCoordMatrix(cart_Pos, nc)
              exit
           end if
        end do
        if(lonlatDGElemId(i,j)==-1) then
           write(*,*) "Can't find corresponding point..", i, j, x_Lon(i)*180/PI, y_Lat(j)*180/PI
        end if
     end do
  end do
write(*,*) "Finish to search.."

  !
  !
  !
  call HistoryCreate( &               
    & file='dataAnalysis.nc', title='data analysis', &
    & source='Sample program of gtool_history/gtool5',   &
    & institution='GFD_Dennou Club davis project',       &
    & dims=(/'lon ', 'lat ', 'time'/), dimsizes=(/ iMax, jMax, 0 /), &
    & longnames=(/ 'longitude', 'latitude ', 'time     '/),       &
    & units=(/ 'degree_east ', 'degree_north', 's           '/),              &
    & origin=real(0), interval=real(outputIntrVal) )

  call HistoryPut('lon', x_lon*180d0/PI)
  call HistoryPut('lat', y_lat*180d0/PI)

  call HistoryAddVariable( &
    & varname='BtmTopl', dims=(/ 'lon ', 'lat ' /), &
    & longname='topology at bottom', units='m', xtype='double')

  call HistoryAddVariable( &
    & varname='CellArea', dims=(/ 'lon ', 'lat ' /), &
    & longname='usage of mesh', units='m2', xtype='doube')

  call HistoryAddVariable( &
    & varname='MeshUsage', dims=(/ 'lon ', 'lat ' /), &
    & longname='usage of mesh', units='1', xtype='int')
  
  call HistoryAddVariable( &
    & varname='h', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='layer depth', units='m', xtype='double')

  call HistoryAddVariable( &
    & varname='SurfHeight', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='surface elevetion', units='m', xtype='double')

  call HistoryAddVariable( &
    & varname='dhdt', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='surface elevetion', units='m', xtype='double')

  call HistoryAddVariable( &
    & varname='hError', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='surface elevetion', units='m', xtype='double')

  call HistoryAddVariable( &
    & varname='U', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='zonal velocity', units='m/s', xtype='double')

  call HistoryAddVariable( &
    & varname='V', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='meridional velocity', units='m/s', xtype='double')

  call HistoryAddVariable( &
    & varname='hU1', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='momentum', units='m2/s', xtype='double')

  call HistoryAddVariable( &
    & varname='hU2', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='momentum', units='m2/s', xtype='double')

  call HistoryAddVariable( &
    & varname='etc', dims=(/ 'lon ', 'lat ', 'time'/), &
    & longname='momentum', units='m2/s', xtype='double')

  call HistoryAddVariable( &
    & varname='l2norm', dims=(/'time'/), &
    & longname='l2norm for analystic solution', units='1', xtype='double')

  call HistoryAddVariable( &
    & varname='linfnorm', dims=(/'time'/), &
    & longname='linfnorm for analystic solution', units='1', xtype='double')

  call HistoryAddVariable( &
    & varname='KEDblTime', dims=(/'time'/), &
    & longname='KE doubling time scale', units='year', xtype='double')


  call HistoryAddVariable( &
    & varname='TEVariation', dims=(/'time'/), &
    & longname='variation of total energy', units='1', xtype='double')

  call HistoryAddVariable( &
    & varname='MassVariation', dims=(/'time'/), &
    & longname='varition of  mass', units='1', xtype='double')


end subroutine OutputData_Init


subroutine OutputData_Final()

!  call netcdfDataWriter_Final(ncWriter)

!  call GeometricField_Final(p_zeta)

end subroutine OutputData_Final

subroutine Output_FieldData_d(fieldName, wc)
  character(*), intent(in) :: fieldName
  real(DP), intent(in) :: wc(:,:)

  call HistoryPut(trim(fieldName), lonlat_interpolate_wc(wc))

end subroutine Output_FieldData_d

subroutine Output_FieldData_i(fieldName, wc)
  use DGElement_mod
  character(*), intent(in) :: fieldName
  integer, intent(in) :: wc(:,:)

  call HistoryPut(trim(fieldName), nint(lonlat_interpolate_wc(dble(wc))))

end subroutine Output_FieldData_i

subroutine OutputData( tstep )

  use LagrangePolyn_mod

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName


  integer :: nc, nk
  real(DP), allocatable :: wc_U(:,:), wc_V(:,:)
  type(vector3d) :: b_1, b_2, cartVel, geoVel
  real(DP) :: totDepth

  real(DP), allocatable :: wc_TriArea(:,:)
  real(DP) :: Jacobi(nDGNodePerElem)

  call DGHelper_MallocElemNode(wc_U)
  call DGHelper_MallocElemNode(wc_V)

  if(tstep == 0) then
     call Output_FieldData('BtmTopl', wc_BtmTopl)
     call Output_FieldData('MeshUsage', wc_DGNodeUsage)

     call DGHelper_MallocElemNode(wc_TriArea)
     !$omp parallel do private(Jacobi)
     do nc=1, nDGElement
        Jacobi = get_DGElemJacobian(nc)
        wc_TriArea(:,nc) = TriNk_sinteg(Jacobi)
     end do
     call Output_FieldData('CellArea', wc_TriArea)
  end if

  !$omp parallel do private(nk,totDepth,b_1,b_2,geoVel)
  do nc=1,nDGElement
     do nk=1,nDGNodePerElem
        totDepth = wc_h(nk,nc) + meanDepth
        b_1 = get_DGElemCovariantBasis(1,DGElemInfo%node(nk),nc)
        b_2 = get_DGElemCovariantBasis(2,DGElemInfo%node(nk),nc)
        geoVel = CartToSphVec(wc_U1(nk,nc)*b_1 + wc_U2(nk,nc)*b_2, wc_DGNodePos(nk,nc))
        wc_U(nk,nc) = geoVel%v_(1)
        wc_V(nk,nc) = geoVel%v_(2)
     end do
  end do

  call Output_FieldData('h', wc_h)
  call Output_FieldData('SurfHeight', meanDepth+wc_h+wc_BtmTopl)
!  call Output_FieldData('dhdt', wc_dhdt)
  call Output_FieldData('hU1', wc_hU1)
  call Output_FieldData('hU2', wc_hU2)
  call Output_FieldData('U', wc_U)
  call Output_FieldData('V', wc_V)
  call Output_FieldData('etc', wc_etc)

end subroutine OutputData

subroutine OutputDataAnalysisInfo(tstep, & 
  & l2norm, linfnorm )

  integer, intent(in) :: tstep
  real(DP), optional :: l2norm, linfnorm

  if( present(l2norm) ) call HistoryPut("l2norm", l2norm)
  if( present(linfnorm) ) call HistoryPut("linfnorm", linfnorm)
 
end subroutine OutputDataAnalysisInfo

function lonlat_interpolate_wc(wc_Var) result(lonlat_Var)
  use LagrangePolyn_mod

  real(DP), intent(in) :: wc_Var(:,:)
  real(DP) :: lonlat_Var(0:iMax-1,jMax)

  integer :: i,j, nc

  !$omp parallel do private(i, nc)
  do j=1, jMax
     do i=0,iMax-1
        nc = lonlatDGElemId(i,j)
        lonlat_Var(i,j) = TriNk_interpolate( &
             & lonlatPt2LclMatrix(i,j)%v_(1), lonlatPt2LclMatrix(i,j)%v_(2), &
             & wc_Var(:,nc) )
     end do
  end do

end function lonlat_interpolate_wc

logical function  isPtContainedTri(pt, ic)
  use DGElement_mod, only: Triangle
 
  type(Vector3d), intent(in) :: pt
  integer, intent(in) :: ic

  type(Triangle) :: tri
  integer :: i
  type(vector3d) :: triNormal

  tri = get_DGElementTri(ic)
!!$  triNormal = (tri%node(2) - tri%node(1)).cross.(tri%node(3) - tri%node(1))
!!$
!!$  do i=1,3
!!$     if( &
!!$       &      (((tri%node(1+mod(i,3)) - tri%node(i)).cross.(pt - tri%node(i))).dot.triNormal) < 0d0 &
!!$       & .or. &
!!$       & geodesicArcLength(pt, tri%node(i)) > PI*Radius/10d0                                   &
!!$       & )  then
!!$        isPtContainedTri = .false.
!!$        return
!!$     end if
!!$  end do
!!$
!!$  isPtContainedTri = .true.

  isPtContainedTri = isPtInsideSphericalTri(normalizedVec(pt), &
       & normalizedVec(tri%node(1)), normalizedVec(tri%node(2)), normalizedVec(tri%node(3)) &
       & )

end function isPtContainedTri

function get_LocalCoordMatrix(cartPos, nc) result(pt)
  use DGElement_mod, only: Triangle

  type(vector3d), intent(in) :: cartPos
  integer, intent(in) :: nc
  type(vector2d) :: pt

  type(vector2d) :: dy
  type(Triangle) :: tri
  integer :: i
  real(DP), parameter :: EPS = 1d-6

  pt = (/ 1d0/3d0, 1d0/3d0 /)
  tri = get_DGElementTri(nc)

  do i=1, 100
     dy%v_ = matmul(inverseMat(Jacob(pt)), cartPos%v_(1:2)-gammaE(pt))
     if(l2norm(dy) < EPS) exit
     pt = pt + dy
  end do
contains
function Jacob(y) result(ret)
  type(vector2d), intent(in) ::y
  real(DP) :: ret(2,2)

  type(vector3d) :: b_1, b_2
  
  b_1 = get_DGElemCovariantBasis(1,y,tri)
  b_2 = get_DGElemCovariantBasis(2,y,tri)

  ret(1,:) = (/ b_1%v_(1), b_2%v_(1) /) 
  ret(2,:) = (/ b_1%v_(2), b_2%v_(2) /) 
end function Jacob

function gammaE(y) result(ret)
  type(vector2d), intent(in) :: y
  real(DP) :: ret(2)

  type(vector3d) :: x_p
  
  x_p = mapping_local2globalCoord(y,tri)
  ret = x_p%v_(1:2)/l2norm(x_p)*Radius

end function gammaE

end function get_LocalCoordMatrix

end module OutputData_mod
