module GridSet_mod
  use dc_types
  use dc_message

  use VectorSpace_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use HexTriIcMesh_mod
  use GeometricField_mod

  use DGElement_mod, only: &
       & Triangle, &
       & DGTriElement, DGTriElement_Init, DGTriElement_Final

  use LagrangePolyn_mod
  use SimParameters_mod, only: &
       & gridFilePath, gridUsageFilePath, &
       & Radius, &
       & PolyDegree

  implicit none
  private

  interface get_DGElemCovariantBasis
     module procedure get_DGElemCovariantBasis1
     module procedure get_DGElemCovariantBasis2
  end interface get_DGElemCovariantBasis

  interface get_DGElemChristoffelSymbl
     module procedure get_DGElemChristoffelSymbl1
  end interface get_DGElemChristoffelSymbl

  public :: GridSet_Init, GridSet_Final
  public :: GridSet_prepair
  public :: mapping_local2globalCoord

  public :: get_DGElemCovariantBasis, get_DGElemCovariantBasisAtNodes
  public :: get_DGElemContravariantBasis
  public :: get_DGElemChristoffelSymbl, get_DGElemChristoffelSymblAtNodes

  public :: calc_G_ij, calc_Gij
  public :: get_DGElemJacobian, get_DGElemSIntNodeJacobian 

  public :: inverseMat
  public :: get_DGElementTri

  type(PolyMesh), public :: plMesh  
  type(fvMeshInfo), public :: fvmInfo
  type(HexTriIcMesh), public :: htiMesh


  type(DGTriElement), public :: DGElemInfo
  integer, public :: nDGElement
  integer, public :: nDGFace
  integer, public :: nDGNodePerElem
  integer, public :: nDGNodePerFace
  integer, public :: nDGSIntNodePerElem

  type(vector3d), allocatable, public :: wc_DGNodePos(:,:)

  integer, dimension(:,:), allocatable, public :: wc_DGNodeUsage
  integer, parameter, public :: NODEUSAGE_OCEAN = 1
  integer, parameter, public :: NODEUSAGE_OCEANSHORE = 2
  integer, parameter, public :: NODEUSAGE_LAND = 3
  integer, parameter, public :: NODEUSAGE_LANDSHORE = 4
  logical, dimension(:), allocatable, public :: c_SeaCellFlag
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'GridSet_mod' !< Module Name

contains
subroutine GridSet_Init()

  !
  call DGTriElement_Init(DGElemInfo, PolyDegree)
 

end subroutine GridSet_Init

subroutine GridSet_Final()


  if(allocated(wc_DGNodePos)) then
     deallocate(wc_DGNodePos)
     deallocate(wc_DGNodeUsage)
  end if

  call fvMeshInfo_Final(fvmInfo)
  call HexTriIcMesh_Final(htiMesh)
  call PolyMesh_Final(plMesh)
  call DGTriElement_Final(DGElemInfo)
  
end subroutine GridSet_Final

!> @brief 
!!
!!
subroutine GridSet_prepair()
  
  use SphericalCoord_mod

  ! 宣言文; Declaration statement
  !
  
  
  ! 局所変数
  ! Local variables
  !
  integer :: nk, nc
  integer :: triNodeCellIds(3)
  type(Triangle) :: tri
  type(vector3d) :: nodePos
  integer :: boundaryId
  type(DomainBoundary), pointer :: boundary => null()

  ! 実行文; Executable statement
  !

  ! Load grid data of HexTriIcMesh from netcdf file. 
  call load_gridData()
  call HexTriIcMesh_Init(htiMesh, plMesh, Radius)


  ! Initialize some modules for finite volume method
  call MessageNotify( 'M', module_name, &
       & "Initialize some modules for finite volume method..")

  call fvMeshInfo_Init(fvmInfo, plMesh, dualMeshFlag=.true.)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvmInfo)

  !
  nDGElement = GetPointListSize(plMesh)
  nDGFace    = GetFaceListSize(plMesh)
  nDGNodePerElem = DGElemInfo%nNode
  nDGNodePerFace = PolyDegree + 1
  nDGSIntNodePerElem = DGElemInfo%nSIntNode

  call MessageNotify( 'M', module_name, &
       & "nDGElement=%d, nDGNodePerElement=%d, nDGNodePerFace=%d, nDGSIntNodePerElem=%d", &
       & i=(/ nDGElement, nDGNodePerElem, nDGNodePerFace, nDGSIntNodePerElem  /))

  !
  allocate(wc_DGNodePos(nDGNodePerElem,nDGElement))
  allocate(c_SeaCellFlag(nDGElement))

  do nc=1, nDGElement
     tri = get_DGElementTri(nc)
     do nk=1,nDGNodePerElem
        wc_DGNodePos(nk,nc) = Radius * normalizedVec( &
             & mapping_local2globalCoord(DGElemInfo%node(nk), tri) )

     end do
  end do

  !
  !

  call assign_DGElemNodeUsage()

end subroutine GridSet_prepair

subroutine assign_DGElemNodeUsage()
  
use SimParameters_mod
use SphericalCoord_mod
  integer :: triNodeCellIds(3)
  integer :: triVxUsage(3)
  type(volScalarField) :: v_MeshUsageMap
  type(pointScalarField) :: p_MeshUsageMap
  integer :: nc, i, boundaryEdgeJudge
  integer :: restNodeUsage
type(vector3d) :: sphPos

  call Load_gridUsageData(v_MeshUsageMap, p_MeshUsageMap)

  allocate(wc_DGNodeUsage(nDGNodePerElem,nDGElement))

  wc_DGNodeUsage = 100000000

  do nc=1, nDGElement
     triNodeCellIds(1) = fvmInfo%Point_CellId(1,nc)
     triNodeCellIds(2) = fvmInfo%Point_CellId(3,nc)
     triNodeCellIds(3) = fvmInfo%Point_CellId(2,nc)
     triVxUsage = nint(v_MeshUsageMap%data%v_(1, triNodeCellIds))

     !
     c_SeaCellFlag(nc) = (nint(p_MeshUsageMap%data%v_(1,nc)) == 1)
     wc_DGNodeUsage(:,nc) = nint(p_MeshUsageMap%data%v_(1,nc))

     !
     do i=1, 3
        if(triVxUsage(i)==10.and.triVxUsage(1+mod(i,3))==10)then
           boundaryEdgeJudge = check_BoundaryEdge(triNodeCellIds(i), triNodeCellIds(1+mod(i,3)))
           if(boundaryEdgeJudge /= 0) then
              wc_DGNodeUsage((i-1)*(nDGNodePerFace-1)+1:i*(nDGNodePerFace-1),nc) = boundaryEdgeJudge*10
              wc_DGNodeUsage(mod(i*(nDGNodePerFace-1),3*(nDGNodePerFace-1))+1,nc) = boundaryEdgeJudge*10
           end if
        end if
     end do

  end do

  call GeometricField_Final(v_MeshUsageMap)
  call GeometricField_Final(p_MeshUsageMap)

contains
function check_BoundaryEdge(p1Id, p2Id) result(ret)
  use SphericalCoord_mod
  integer, intent(in) :: p1Id, p2Id
  integer:: ret

  type(DomainBoundary), pointer :: boundary=>null()
  integer :: i, j, nBoundaryCell


  do i=1, getBoundaryListSize(plMesh)
     boundary => plMesh%boundaryList(i)
     nBoundaryCell = size(boundary%boundaryElemIdList)
     do j=1, nBoundaryCell
        if(    p1Id==boundary%boundaryElemIdList(j) .and. &
             & p2Id==boundary%boundaryElemIdList(mod(j,nBoundaryCell)+1) ) then
           ret = -1
           return
        end if
        if(    p2Id==boundary%boundaryElemIdList(j) .and. &
             & p1Id==boundary%boundaryElemIdList(mod(j,nBoundaryCell)+1) ) then
           ret = +1
           return
        end if

     end do
  end do

  ret = 0
end function check_BoundaryEdge

end subroutine assign_DGElemNodeUsage

function get_DGElementTri(nc) result(tri)
use SphericalCoord_mod
use SimParameters_mod
  integer, intent(in) :: nc
  type(Triangle) :: tri

  integer :: triNodeCellIds(3)
  integer :: i, maxloc_(1)
  real(DP) :: arcLen(3)
  integer, parameter :: defaultIds(3) = (/ 1,3,2 /)
  integer :: nodeIds(3)


  nodeIds(:) = defaultIds(:)
  triNodeCellIds = fvmInfo%Point_CellId(nodeIds,nc)
  do i=1,3
     tri%node(i) = plMesh%cellPosList(triNodeCellIds(i))
  end do

!!$  arcLen(1) = geodesicArcLength(tri%node(2),tri%node(3))
!!$  arcLen(2) = geodesicArcLength(tri%node(3),tri%node(1))
!!$  arcLen(3) = geodesicArcLength(tri%node(1),tri%node(2))
!!$  maxloc_ = maxloc(arcLen(:))
!!$  if( maxloc_(1)  /= 1 ) then
!!$     nodeIds = cshift(defaultIds, maxloc_(1)-1)
!!$     triNodeCellIds = fvmInfo%Point_CellId(nodeIds,nc)
!!$     do i=1,3
!!$        tri%node(i) = plMesh%cellPosList(triNodeCellIds(i))
!!$     end do
!!$     write(*,*)  arcLen(1)/arcLen(2), arcLen(1)/arcLen(3)!, arcLen
!!$     call print(RadToDegUnit(CartToSphPos(tri%node(1))))
!!$     call print(RadToDegUnit(CartToSphPos(tri%node(2))))
!!$     call print(RadToDegUnit(CartToSphPos(tri%node(3))))
!!$  end if

end function get_DGElementTri

  

  pure function calc_G_ij(b_1, b_2) result(G_ij)
    type(vector3d), intent(in) :: b_1, b_2
    real(DP) :: G_ij(2,2)

    G_ij(1,1) = b_1.dot.b_1; G_ij(1,2) = b_1.dot.b_2; 
    G_ij(2,1) = b_2.dot.b_1; G_ij(2,2) = b_2.dot.b_2; 
  end function calc_G_ij

  pure function calc_Gij(G_ij) result(Gij)
    real(DP), intent(in) :: G_ij(2,2)
    real(DP) :: Gij(2,2)
    
    real(DP) :: detG_ij

    detG_ij = G_ij(1,1)*G_ij(2,2) - G_ij(1,2)**2
    Gij(1,1) =  G_ij(2,2)/detG_ij
    Gij(1,2) = -G_ij(2,1)/detG_ij
    Gij(2,1) = Gij(1,2)
    Gij(2,2) =  G_ij(1,1)/detG_ij

  end function calc_Gij

  function get_DGElemJacobian(nc) result(Jacobi)
    integer, intent(in) :: nc
    real(DP) :: Jacobi(nDGNodePerElem)

    integer :: nk
    real(DP) :: G_ij(2,2)

    do nk=1, nDGNodePerElem
       Jacobi(nk) = get_Jacobian(DGElemInfo%node(nk), nc)
    end do

  end function get_DGElemJacobian

  function get_DGElemSIntNodeJacobian(nc) result(Jacobi)
    integer, intent(in) :: nc
    real(DP) :: Jacobi(nDGSIntNodePerElem)

    integer :: nk

    do nk=1, nDGSIntNodePerElem
       Jacobi(nk) = get_Jacobian(DGElemInfo%sIntNode(nk), nc)
    end do
  end function get_DGElemSIntNodeJacobian

  function get_Jacobian(y, nc) result(Jacobi)
    type(vector2d), intent(in) :: y
    integer, intent(in) :: nc
    real(DP) :: Jacobi

    real(DP) :: G_ij(2,2)
    
    G_ij(:,:) = calc_G_ij( &
         & get_DGElemCovariantBasis(1, y, nc), &
         & get_DGElemCovariantBasis(2, y, nc) )
    Jacobi = sqrt(G_ij(1,1)*G_ij(2,2) - G_ij(1,2)*G_ij(2,1))
  end function get_Jacobian

  function mapping_local2globalCoord(y,tri) result(ret)
    type(vector2d), intent(in) :: y
    type(Triangle), intent(in) :: tri
    type(vector3d) :: ret

    type(vector3d) :: x0, x1, x2

    x0 = tri%node(1); x1 = tri%node(2); x2 = tri%node(3)
    ret = x0 + y%v_(1)*(x1 - x0) + y%v_(2)*(x2 - x0)

  end function mapping_local2globalCoord

  function get_DGElemChristoffelSymbl1(i,j,k,y,tri) result(ret)
    integer, intent(in) :: i, j, k
    type(vector2d), intent(in) :: y
    type(Triangle), intent(in) :: tri
    real(DP) :: ret

    type(vector3d) :: xp, tmp1, tmp2
    
    tmp1 = 0d0; tmp2 = 0d0
    if(i==k) then
       tmp1 = tri%node(j+1) - tri%node(1)
    end if

    if(i==j) then
       tmp2 = tri%node(k+1) - tri%node(1)
    end if

    xp = mapping_local2globalCoord(y, tri)
    ret = - ( (xp.dot.(tmp1+tmp2))/(xp.dot.xp) )

  end function get_DGElemChristoffelSymbl1

  function get_DGElemChristoffelSymblAtNodes(i,j,k,nc) result(ret)

    integer, intent(in) :: i, j, k
    integer, intent(in) :: nc
    real(DP) :: ret(nDGNodePerElem)

    type(Triangle) :: tri
    integer :: nk

    tri = get_DGElementTri(nc)
    do nk=1, nDGNodePerElem
       ret(nk) = get_DGElemChristoffelSymbl1(i,j,k,DGElemInfo%node(nk),tri)
    end do
    
  end function get_DGElemChristoffelSymblAtNodes

  function get_DGElemCovariantBasis1(i,y,tri) result(ret)

    integer, intent(in) :: i
    type(vector2d), intent(in) :: y
    type(Triangle), intent(in) :: tri
    type(vector3d) :: ret

    real(DP) :: xp_abs
    type(vector3d) :: xp, x0, xi, xp_unitvec

    x0 = tri%node(1); xi = tri%node(i+1)
    xp = mapping_local2globalCoord(y, tri)
    xp_abs = l2norm(xp)
    xp_unitvec = xp/xp_abs

    ret = (Radius/xp_abs) * ((xi - x0) - ((xi - x0).dot.xp_unitvec)*xp_unitvec)
    
  end function get_DGElemCovariantBasis1

  function get_DGElemCovariantBasis2(i,y,nc) result(ret)

    integer, intent(in) :: i
    type(vector2d), intent(in) :: y
    integer, intent(in) :: nc
    type(vector3d) :: ret

    ret = get_DGElemCovariantBasis1(i,y,get_DGElementTri(nc))
  end function get_DGElemCovariantBasis2

  function get_DGElemCovariantBasisAtNodes(i,nc) result(ret)

    integer, intent(in) :: i
    integer, intent(in) :: nc
    type(vector3d) :: ret(nDGNodePerElem)

    type(Triangle) :: tri
    integer :: nk

    tri = get_DGElementTri(nc)
    do nk=1, nDGNodePerElem
       ret(nk) = get_DGElemCovariantBasis1(i,DGElemInfo%node(nk),tri)
    end do
    
  end function get_DGElemCovariantBasisAtNodes

  subroutine get_DGElemContravariantBasis(b1, b2, b_1, b_2)
    type(vector3d), intent(out) :: b1, b2
    type(vector3d), intent(in) :: b_1, b_2

    real(DP) :: Gij(2,2)

    Gij = calc_Gij( calc_G_ij(b_1, b_2) )

    b1 = Gij(1,1)*b_1 + Gij(1,2)*b_2
    b2 = Gij(2,1)*b_1 + Gij(2,2)*b_2
  end subroutine get_DGElemContravariantBasis

  function inverseMat(A) result(Ainv)
    real(DP), intent(in) :: A(:,:)
    real(DP) :: Ainv(size(A,1),size(A,2))

    real(DP) :: work(size(A,2)*128)
    integer :: IPIV(size(A,2))
    integer :: INFO
    integer :: m, n, i

    m = size(A,1)
    n = size(A,2)

    Ainv = A
    call DGETRF(m, n, Ainv, m, IPIV(1:min(m,n)), info)
    call DGETRI(n, Ainv, m, IPIV, work, n*64, info)

    if(info /= 0) then
       write(*,*) shape(A)
       do i=1, m
          write(*,*) A(i,:)
       end do
       call MessageNotify('E', module_name, &
            & "Fail to inverse a matrix." )

!!$    else
!!$       call MessageNotify('M', module_name, &
!!$            & "Success to inverse a matrix." )
    end if

  end function inverseMat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine load_gridData()

  use netcdfDataReader_mod

  type(netcdfDataReader) :: ncReader

  call MessageNotify( 'M', module_name, &
       & "Set up grid. Load grid data from '%a' ..", ca=(/ gridFilePath /) ) 

  call netcdfDataReader_Init(ncReader, gridFilePath, plMesh)
  call netcdfDataReader_Final(ncReader)


end subroutine Load_gridData

subroutine load_gridUsageData(v_MeshUsageMap, p_MeshUsageMap)

  use netcdfDataReader_mod

  type(volScalarField), intent(inout) :: v_MeshUsageMap
  type(pointScalarField), intent(inout) :: p_MeshUsageMap
  type(netcdfDataReader) :: ncReader

  if( len(trim(gridUsageFilePath)) == 0 ) then
     call MessageNotify( 'M', module_name, &
          & "A file to set the usage of grid is not specified. The domain is assumed to be global.")
     call GeometricField_Init(v_MeshUsageMap, plMesh, "v_MeshUsageMap")
     call GeometricField_Init(p_MeshUsageMap, plMesh, "p_MeshUsageMap")
     v_MeshUsageMap = 1d0
     p_MeshUsageMap = 1d0

     return
  end if

  ! Load the usage data form a netcdf file.
  !

  call MessageNotify( 'M', module_name, &
       & "Set up the usage of grid. Load the usage  data from '%a' ..", &
       & ca=(/ gridUsageFilePath /) ) 

  call netcdfDataReader_Init(ncReader, gridUsageFilePath, plMesh, &
       & isLoadMeshData=.false. )
  call netcdfDataReader_get(ncReader, 'v_MeshUsageMap', v_MeshUsageMap)
  call netcdfDataReader_get(ncReader, 'p_MeshUsageMap', p_MeshUsageMap)
  call netcdfDataReader_Final(ncReader)


end subroutine Load_gridUsageData

end module GridSet_mod
