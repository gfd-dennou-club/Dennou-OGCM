!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DGHelper_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP
  
  use dc_message, only: MessageNotify

  use VectorSpace_mod

  use LagrangePolyn_mod

  use SimParameters_mod, only: Radius

  use GridSet_mod, only: &
         & plMesh, fvmInfo, &
         & nDGElement, nDGFace, &
         & nDGNodePerElem, nDGSIntNodePerElem, nDGNodePerFace, &
         & DGElemInfo, &
         & get_DGElemCovariantBasis, get_DGElemContravariantBasis, &
         & get_DGElemCovariantBasisAtNodes, &
         & get_DGElemSIntNodeJacobian, calc_G_ij, inverseMat


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  
  interface DGHelper_MallocElemNode
     module procedure DGHelper_MallocElemNodeScalar
     module procedure DGHelper_MallocElemNodeVec3D
  end interface DGHelper_MallocElemNode

  public :: DGHelper_Init, DGHelper_Final
  public :: DGHelper_prepair

  public :: DGHelper_MallocElemNode, DGHelper_MallocFaceNode
  public :: DGHelper_MallocElemMat, DGHelper_MallocElemSIntMat, DGHelper_MallocFaceMat

  real(DP), dimension(:,:,:), allocatable, public :: wc_Dy1, wc_Dy2, wc_S
  real(DP), dimension(:,:,:), allocatable, public :: wc_Ms
  type(vector3d), dimension(:,:), allocatable, public :: &
       & ws_normVec
  integer, dimension(:,:), allocatable, public :: Face_TriEdgeId
  integer, dimension(:,:), allocatable, public :: Face_EdgeId
  integer, dimension(:,:), allocatable, public :: Face_DGElemId


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DGHelper_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DGHelper_Init()

    !
    !
    use SimParameters_mod, only: &
         & PolyDegree

    ! 実行文; Executable statements
    !

    call DGHelper_MallocElemSIntMat(wc_Dy1)
    call DGHelper_MallocElemSIntMat(wc_Dy2)
    call DGHelper_MallocElemSIntMat(wc_S)
    call DGHelper_MallocFaceMat(wc_Ms)
    allocate(Face_TriEdgeId(2, nDGFace))
    allocate(Face_EdgeId(2, nDGFace))
    allocate(Face_DGElemId(2, nDGFace))
    allocate(ws_normVec(nDGNodePerFace, nDGFace))
    
  end subroutine DGHelper_Init

  !>
  !!
  !!
  subroutine DGHelper_Final()

    ! 実行文; Executable statements
    !
    if(allocated(wc_Dy1)) then
       deallocate(wc_Dy1, wc_Dy2, wc_Ms)
       deallocate(Face_TriEdgeId, Face_EdgeId)
       deallocate(Face_DGElemId)
       deallocate(ws_normVec)
    end if

  end subroutine DGHelper_Final


  subroutine DGHelper_MallocElemNodeScalar(wc_Var)

    use PolyMesh_mod, only: &
         & getPointListSize


    real(DP), intent(inout), allocatable :: wc_Var(:,:)

    integer :: Nk

    if(allocated(wc_Var)) deallocate(wc_Var)

    allocate(wc_Var(nDGNodePerElem, getPointListSize(plMesh)))

  end subroutine DGHelper_MallocElemNodeScalar

  subroutine DGHelper_MallocElemNodeVec3D(wc_Var)

    use PolyMesh_mod, only: &
         & getPointListSize


    type(vector3d), intent(inout), allocatable :: wc_Var(:,:)

    integer :: Nk

    if(allocated(wc_Var)) deallocate(wc_Var)

    allocate(wc_Var(nDGNodePerElem, getPointListSize(plMesh)))

  end subroutine DGHelper_MallocElemNodeVec3D

  subroutine DGHelper_MallocFaceNode(ws_Var)

    use PolyMesh_mod, only: &
         & getFaceListSize


    real(DP), intent(inout), allocatable :: ws_Var(:,:)

    integer :: Nk

    if(allocated(ws_Var)) deallocate(ws_Var)

    allocate(ws_Var(nDGNodePerFace, getFaceListSize(plMesh)))

  end subroutine DGHelper_MallocFaceNode

  subroutine DGHelper_MallocElemMat(wc_Mat)

    use PolyMesh_mod, only: &
         & getPointListSize


    real(DP), intent(inout), allocatable :: wc_Mat(:,:,:)

    integer :: Nk

    if(allocated(wc_Mat)) deallocate(wc_Mat)

    allocate(wc_Mat(nDGNodePerElem,nDGNodePerElem, getPointListSize(plMesh)))

  end subroutine DGHelper_MallocElemMat

  subroutine DGHelper_MallocElemSIntMat(wc_Mat)

    use PolyMesh_mod, only: &
         & getPointListSize


    real(DP), intent(inout), allocatable :: wc_Mat(:,:,:)

    integer :: Nk

    if(allocated(wc_Mat)) deallocate(wc_Mat)

    allocate(wc_Mat(nDGSIntNodePerElem, nDGNodePerElem, getPointListSize(plMesh)))

  end subroutine DGHelper_MallocElemSIntMat

  subroutine DGHelper_MallocFaceMat(wc_Mat)

    use PolyMesh_mod, only: &
         & getPointListSize


    real(DP), intent(inout), allocatable :: wc_Mat(:,:,:)

    integer :: Nk

    if(allocated(wc_Mat)) deallocate(wc_Mat)

    allocate(wc_Mat(nDGNodePerElem,nDGNodePerFace*3, getPointListSize(plMesh)))

  end subroutine DGHelper_MallocFaceMat
  
  !> @brief 
  !!
  !!
  subroutine DGHelper_prepair()

    use GridSet_mod, only: &
         & wc_DGNodePos
    use DGElement_mod, only: &
         & Triangle

    use PolyMesh_mod, only: &
         & Face

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(nDGNodePerElem, nDGNodePerElem) :: & 
         & M_k1k2, Minv_k1k2
    real(DP), dimension(nDGNodePerElem, nDGSIntNodePerElem) :: &
         & Dy1_sk1, Dy2_sk1, S_sk1

    real(DP), dimension(nDGNodePerElem, nDGNodePerFace*3) :: Ms_k1l

    type(vector3d) :: b_1, b_2
    type(vector2d), dimension(nDGNodePerElem) :: DGNodePos
    
    real(DP), dimension(2,2,nDGNodePerElem) :: &
         & G_ij, Gij
    real(DP), dimension(nDGNodePerElem) :: det_G_ij
    real(DP), dimension(nDGNodePerElem) :: Jacob
    integer :: nk, nk1, nk2, nc, nlInt, ns, nl


    type(vector2d) :: edgeTd(3)
    integer :: TriEdge, TriEdge_
    type(Face), pointer :: face_
    integer :: ptId, cellId, faceIds(3)
    integer :: ownElemId, nodeId_OwnDGElem, TriEdge_ownDGElem


    ! 実行文; Executable statement
    !
    
    DGNodePos = DGElemInfo%node

    edgeTd(1) = (/ 1d0, 0d0 /); 
    edgeTd(2) = (/ -1d0/sqrt(2d0), 1d0/sqrt(2d0) /); 
    edgeTd(3) = (/ 0d0,  -1d0 /); 


    !$omp parallel do private(face_, ptID, TriEdge, TriEdge_, TriEdge_ownDGElem, ownElemId, nodeId_OwnDGElem, nl, b_1, b_2)
    do ns=1, nDGFace
       face_ => plMesh%faceList(ns) 

       ptID = face_%vertIdList(2)
       Face_DGElemId(:,ns) = (/ ptID, face_%vertIdList(1) /)
!!$       ptID = face_%vertIdList(1)
!!$       Face_DGElemId(:,ns) = (/ ptID, face_%vertIdList(2) /)

       TriEdge_ownDGElem = -1
       do TriEdge=1,3

          TriEdge_ = 4 - TriEdge !!!!< 

          if( fvmInfo%Point_CellId(TriEdge,Face_DGElemId(1,ns)) == face_%neighCellId ) then
             Face_TriEdgeId(1,ns) = (nDGNodePerFace - 1)*(TriEdge_ - 1) + 1
             Face_EdgeId(1,ns) = TriEdge_ 
             TriEdge_ownDGElem = TriEdge_
          end if
          if( fvmInfo%Point_CellId(TriEdge,Face_DGElemId(2,ns)) == face_%ownCellId ) then
             Face_TriEdgeId(2,ns) = (nDGNodePerFace - 1)*(TriEdge_ - 1) + 1             
             Face_EdgeId(2,ns) = TriEdge_
          end if
       end do

       ownElemId = Face_DGElemId(1,ns)
       nodeId_OwnDGElem = Face_TriEdgeId(1,ns)

       do nl=1, nDGNodePerFace
          b_1 = get_DGElemCovariantBasis(1, DGElemInfo%node(nodeId_OwnDGElem), ownElemId)
          b_2 = get_DGElemCovariantBasis(2, DGElemInfo%node(nodeId_OwnDGElem), ownElemId)

          ws_normVec(nl,ns) = normalizedVec( &
               & (   edgeTd(TriEdge_ownDGElem)%v_(1)*b_1   &
               &   + edgeTd(TriEdge_ownDGElem)%v_(2)*b_2   &
               & ).cross.(wc_DGNodePos(nodeId_OwnDGElem, ownElemId)) &
               & )

          nodeId_OwnDGElem = mod(nodeId_OwnDGElem,3*(nDGNodePerFace-1)) + 1          
       end do
    end do ! End loop for DG face

    !$omp parallel do private(M_k1k2, Dy1_sk1, Dy2_sk1, Ms_k1l, Minv_k1k2, S_sk1)
    do nc=1, nDGElement

       call prepair_DGMatPerElement(M_k1k2, Dy1_sk1, Dy2_sk1, Ms_k1l, S_sk1, &
            & nc )

       Minv_k1k2 = inverseMat(M_k1k2)
       wc_Dy1(:,:,nc) = transpose(matmul(Minv_k1k2, Dy1_sk1))
       wc_Dy2(:,:,nc) = transpose(matmul(Minv_k1k2, Dy2_sk1))
       wc_Ms(:,:,nc) = matmul(Minv_k1k2, Ms_k1l)
       wc_S(:,:,nc) = transpose(matmul(Minv_k1k2, S_sk1))

    end do  ! End loop for DG elements


  end subroutine DGHelper_prepair

  subroutine prepair_DGMatPerElement(M_k1k2, Dy1_sk1, Dy2_sk1, Ms_k1l, S_sk1, &
       & nc )

    real(DP), dimension(nDGNodePerElem,nDGNodePerElem), intent(out) :: M_k1k2
    real(DP), dimension(nDGNodePerElem, nDGSIntNodePerElem), intent(out) :: Dy1_sk1, Dy2_sk1, S_sk1
    real(DP), intent(out) :: Ms_k1l(nDGNodePerElem,3*nDGNodePerFace)
    integer, intent(in) :: nc

    integer :: nk, nk1, nk2, ns
    type(vector3d), dimension(nDGNodePerElem) :: b_1, b_2, b1, b2
    real(DP), dimension(nDGNodePerElem,nDGNodePerElem) :: phi, phi_dy1, phi_dy2
    real(DP), dimension(nDGSIntNodePerElem,nDGNodePerElem) :: s_phi, s_phi_dy1, s_phi_dy2
    real(DP) :: l_phi(nDGNodePerFace), l_zero(nDGNodePerFace)
    real(DP), dimension(2,2,nDGNodePerElem) :: G_ij, Gij
    real(DP) :: det_G_ij(nDGNodePerElem), Jacob(nDGNodePerElem)
    real(DP) :: Jacob_(nDGSIntNodePerElem)
    real(DP) :: y1, y2
    integer :: faceIds(3), TriEdge, nlInt
    integer :: DGElemNodeId_TriEdge(nDGNodePerFace,3)
    real(DP) :: l_Jacob(nDGNodePerFace)
    real(DP), dimension(nDGSIntNodePerElem) :: sintTmp_M, sintTmp_D1, sintTmp_D2
 
    type(vector3d) :: edgeNe     ! vector normal to a edge of triangle on global coordinate
    type(vector2d) :: edgeNd(3)  ! vector normal to a edge of triangle on local coordinate

    edgeNd(1) = (/ 0d0, -1d0 /); 
    edgeNd(2) = (/ 1d0/sqrt(2d0), 1d0/sqrt(2d0) /); 
    edgeNd(3) = (/ -1d0,  0d0 /); 
    
    b_1(:) = get_DGElemCovariantBasisAtNodes(1, nc)
    b_2(:) = get_DGElemCovariantBasisAtNodes(2, nc)

    do nk=1, nDGNodePerElem
       G_ij(:,:,nk) = calc_G_ij(b_1(nk), b_2(nk))
       det_G_ij(nk) = G_ij(1,1,nk)*G_ij(2,2,nk) - G_ij(2,1,nk)*G_ij(1,2,nk)
       Jacob(nk) = sqrt(det_G_ij(nk))

       Gij(:,:,nk) = inverseMat(G_ij(:,:,nk))
       call get_DGElemContravariantBasis( b1(nk), b2(nk), &
            & b_1(nk), b_2(nk))

       y1 = DGElemInfo%node(nk)%v_(1); y2 = DGElemInfo%node(nk)%v_(2)
       phi(nk,:) = TriNk_basis(y1, y2)
       phi_dy1(nk,:) = TriNk_basis_dy1(y1, y2)
       phi_dy2(nk,:) = TriNk_basis_dy2(y1, y2)
    end do

    
    ! Set the derivative matrix to calculate surface integral associated with  the product of 
    ! spatial derivative term of test function and flux term over an element.
    !
    Jacob_(:) = get_DGElemSIntNodeJacobian(nc)

    do ns=1, nDGSIntNodePerElem
       y1 = DGElemInfo%sIntNode(ns)%v_(1); y2 = DGElemInfo%sIntNode(ns)%v_(2)
       s_phi(ns,:) = TriNk_basis(y1,y2)
       s_phi_dy1(ns,:) = TriNk_basis_dy1(y1,y2)
       s_phi_dy2(ns,:) = TriNk_basis_dy2(y1,y2)
    end do

    do nk2=1, nDGNodePerElem
       do nk1=1, nDGNodePerElem
          sintTmp_M(:) = Jacob_(:)*s_phi(:,nk1)*s_phi(:,nk2)
          M_k1k2(nk1,nk2) = TriNk_sinteg_dotProdWt(sintTmp_M)
       end do   ! end of loop for nk1       
    end do      ! end of loop for nk2


    do nk1=1, nDGNodePerElem
       do ns=1, nDGSIntNodePerElem
          Dy1_sk1(nk1,ns) = Jacob_(ns)*s_phi_dy1(ns,nk1)
          Dy2_sk1(nk1,ns) = Jacob_(ns)*s_phi_dy2(ns,nk1)
          S_sk1(nk1,ns) = Jacob_(ns)*s_phi(ns,nk1)
       end do
    end do

    ! Set the matrix used to calculate line integral of flux across the boundary of an element.  
    !

    do nlInt=1,nDGNodePerFace
       do TriEdge=1,3
          DGElemNodeId_TriEdge(nlInt,TriEdge) = (nDGNodePerFace-1)*(TriEdge-1) + nlInt
          if(TriEdge==3 .and. nlInt==nDGNodePerFace) DGElemNodeId_TriEdge(nlInt,TriEdge) = 1
       end do
    end do

    faceIds(:) = fvmInfo%Point_FaceId(3:1:-1, nc)
!!$    faceIds(:) = fvmInfo%Point_FaceId(:, nc)

    do nk1=1, nDGNodePerElem
       do TriEdge=1,3

          do nlInt=1,nDGNodePerFace
             nk = DGElemNodeId_TriEdge(nlInt,TriEdge)

             edgeNe = edgeNd(TriEdge)%v_(1)*b1(nk) + edgeNd(TriEdge)%v_(2)*b2(nk)
             l_Jacob(nlInt) = Jacob(nk)*l2norm(edgeNe)
          end do

          do nlInt=1,nDGNodePerFace
             l_phi(:) = 0d0; l_phi(nlInt) = 1d0; l_zero(:) = 0d0
             select case(TriEdge)
             case(1)
                Ms_k1l(nk1,nDGNodePerFace*(TriEdge-1)+nlInt) &
                     & = TriNk_linteg( &
                     & l_Jacob(:),phi(DGElemNodeId_TriEdge(:,TriEdge),nk1),l_phi(:), l_zero,l_zero,l_zero, l_zero,l_zero,l_zero &
                     & )
             case(2)
                Ms_k1l(nk1,nDGNodePerFace*(TriEdge-1)+nlInt) &
                     & = TriNk_linteg( &
                     & l_zero,l_zero,l_zero, l_Jacob(:),phi(DGElemNodeId_TriEdge(:,TriEdge),nk1),l_phi(:), l_zero,l_zero,l_zero &
                     & )
             case(3)
                Ms_k1l(nk1,nDGNodePerFace*(TriEdge-1)+nlInt) &
                     & = TriNk_linteg( &
                     & l_zero,l_zero,l_zero, l_zero,l_zero,l_zero, l_Jacob(:),phi(DGElemNodeId_TriEdge(:,TriEdge),nk1),l_phi(:) &
                     & )                   
             end select
          end do

          if( Face_DGElemId(2,faceIds(TriEdge)) ==nc ) then
             Ms_k1l(nk1,nDGNodePerFace*(TriEdge-1)+1:nDGNodePerFace*TriEdge) &
                  & = - Ms_k1l(nk1,nDGNodePerFace*TriEdge:nDGNodePerFace*(TriEdge-1)+1:-1)
          end if

       end do    ! end do for TriEdge
    end do   ! end do for nk1

  end subroutine prepair_DGMatPerElement

end module DGHelper_mod

