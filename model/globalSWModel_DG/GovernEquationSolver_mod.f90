module GovernEquationSolver_mod

  !
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify

  !
  !

  use VectorSpace_mod

  use SphericalCoord_mod, only: &
       & SphericalTriArea, &
       & CartToSphPos, radToDegUnit

  use PolyMesh_mod

  use HexTriIcMesh_mod

  use GeometricField_mod

  use fvMeshInfo_mod

  use fvCalculus_mod, only: &
       div, curl, grad

  !
  !
  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod, only: &
       & v_h, v_BtmTopl, v_U1, v_U2, wc_etc

  use DGElement_mod
  use DGHelper_mod

  implicit none
  private

  type(vector3d), allocatable, dimension(:,:) :: wc_b_1, wc_b_2

  public :: GovernEquationSolver_Init, GovernEquationSolver_Final
  public :: Solve_GovernEquation
  
contains
subroutine GovernEquationSolver_Init()

  integer :: nc

  call DGHelper_MallocElemNode(wc_b_1)
  call DGHelper_MallocElemNode(wc_b_2)

  !$omp parallel do 
  do nc=1, nDGElement
     wc_b_1(:,nc) = get_DGElemCovariantBasisAtNodes(1,nc)
     wc_b_2(:,nc) = get_DGElemCovariantBasisAtNodes(2,nc)
  end do

  
end subroutine GovernEquationSolver_Init

subroutine GovernEquationSolver_Final()

  if(allocated(wc_b_1))then
     deallocate(wc_b_1, wc_b_2)
  end if

end subroutine GovernEquationSolver_Final

subroutine Solve_GovernEquation()

  use VariableSet_mod, only: &
       & wc_h, wc_hU1, wc_hU2, &
       & wc_U1, wc_U2

  integer :: i , RKStep
  real(DP) :: dt
  real(DP), parameter :: RKCoef(4) = (/ 1d0, 2d0, 2d0, 1d0 /)

  type :: RKTmpDGVar
     real(DP), allocatable :: val(:,:)
  end type RKTmpDGVar
  
  integer, parameter :: nStage = 3
  integer :: s, j
  real(DP), parameter :: RKCoefMat(1:nStage,0:nStage-1) = &
       & reshape((/  1d0, 0.25d0,   1d0/6d0,  &
       &             0d0, 0.25d0, 1d0/6d0,  &
       &             0d0,   0d0,  2d0/3d0 /),    &
       & (/ 3,3 /))

  type(RKTmpDGVar), dimension(0:nStage) :: wc_hTmp, wc_hv1Tmp, wc_hv2Tmp
  type(RKTmpDGVar), dimension(0:nStage) :: wc_dhdt, wc_dhv1dt, wc_dhv2dt

integer :: nc, nk
  dt = dble(delTime)

  !$omp parallel do
  do s=0,nStage
     call DGHelper_MallocElemNode(wc_hTmp(s)%val)
     call DGHelper_MallocElemNode(wc_hv1Tmp(s)%val)
     call DGHelper_MallocElemNode(wc_hv2Tmp(s)%val)

     call DGHelper_MallocElemNode(wc_dhdt(s)%val)
     call DGHelper_MallocElemNode(wc_dhv1dt(s)%val)
     call DGHelper_MallocElemNode(wc_dhv2dt(s)%val)

     wc_hTmp(s)%val = wc_h
     wc_hv1Tmp(s)%val = (meanDepth + wc_h)*wc_U1;   
     wc_hv2Tmp(s)%val = (meanDepth + wc_h)*wc_U2;   
  end do


  do s=1, nStage

     call calcTendency( &
          & wc_dhdt(s-1)%val, wc_dhv1dt(s-1)%val, wc_dhv2dt(s-1)%val, &
          & wc_hTmp(s-1)%val, wc_hv1Tmp(s-1)%val, wc_hv2Tmp(s-1)%val, &
          & wc_U1, wc_U2 )
     
     do j=0, s-1
        !$omp workshare
        wc_hTmp(s)%val = wc_hTmp(s)%val + &
             & dt*RKCoefMat(s,j)*wc_dhdt(j)%val

        wc_hv1Tmp(s)%val = wc_hv1Tmp(s)%val + &
             & dt*RKCoefMat(s,j)*wc_dhv1dt(j)%val

        wc_hv2Tmp(s)%val = wc_hv2Tmp(s)%val + &
             & dt*RKCoefMat(s,j)*wc_dhv2dt(j)%val
        !$omp end workshare
     end do

     !$omp workshare
!!$     wc_hv1Tmp(s)%val = (meanDepth + wc_hTmp(s)%val) * wc_U1
!!$     wc_hv2Tmp(s)%val = (meanDepth + wc_hTmp(s)%val) * wc_U2
     wc_U1 = wc_hv1Tmp(s)%val/(meanDepth + wc_hTmp(s)%val)
     wc_U2 = wc_hv2Tmp(s)%val/(meanDepth + wc_hTmp(s)%val)
     !$omp end workshare
  end do

  !$omp workshare
  wc_h = wc_hTmp(nStage)%val
  wc_hU1 = wc_hv1Tmp(nStage)%val
  wc_hU2 = wc_hv2Tmp(nStage)%val
  wc_etc = wc_dhdt(0)%val
  !$omp end workshare

!!$  do nc=1, nDGElement
!!$     if(.not.c_SeaCellFlag(nc)) then
!!$        wc_h(:,nc) = 10d0
!!$        wc_hU1(:,nc) = 10d0
!!$        wc_hU2(:,nc) = 10d0
!!$     end if
!!$  end do

end subroutine Solve_GovernEquation

subroutine calcTendency(wc_dhdt, wc_dhU1dt, wc_dhU2dt, &
     & wc_h, wc_hU1, wc_hU2, wc_U1, wc_U2 )

  use LagrangePolyn_mod
  use DGHelper_mod
  use SphericalCoord_mod

  real(DP), intent(inout), dimension(:,:) ::  &
       & wc_dhdt, wc_dhU1dt, wc_dhU2dt, wc_U1, wc_U2

  real(DP), intent(in), dimension(:,:) :: wc_h, wc_hU1, wc_hU2


  real(DP), allocatable, dimension(:,:) :: ws_hFlux
  real(DP), allocatable, dimension(:,:,:) :: lsc_hU1Flux, lsc_hU2Flux
  real(DP), dimension(nDGNodePerElem) :: &
       & w_hFlux, w_hU1Flux, w_hU2Flux, w_totDepth, &
       & w_D_hFlux, w_D_hU1Flux, w_D_hU2Flux, w_D1_hU1Flux, w_D2_hU1Flux, w_hU1Src, w_hU2Src
  real(DP), dimension(nDGSIntNodePerElem) :: &
       & s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2

  integer :: nc, nk, ns, nl, ne
  integer :: c_in, c_out, w_in, w_out
  integer :: faceIds(3), node0IDs(3)
  
  type(vector3d) :: geo_pos
  real(DP) :: numFlux(5)
  type(Triangle) :: tri

  call DGHelper_MallocFaceNode(ws_hFlux)
  allocate(lsc_hU1Flux(nDGNodePerFace,3,nDGElement))
  allocate(lsc_hU2Flux(nDGNodePerFace,3,nDGElement))
  

  !$omp parallel do private(c_in, w_in, c_out, w_out, nl, numFlux) schedule(guided)
  do ns=1, nDGFace

     c_in = Face_DGElemId(1,ns); c_out = Face_DGElemId(2,ns)
     w_in = Face_TriEdgeId(1,ns); w_out = Face_TriEdgeId(2,ns) + nDGNodePerFace - 1
     if(w_out > 3*(nDGNodePerFace - 1)) w_out = mod(w_out, 3*(nDGNodePerFace-1))

     if(c_SeaCellFlag(c_in) .or. c_SeaCellFlag(c_out)) then
        do nl=1, nDGNodePerFace

           numFlux(:) = calc_NumericFlux( &
                & wc_h(w_in,c_in), wc_U1(w_in,c_in), wc_U2(w_in,c_in), &
                & wc_h(w_out,c_out), wc_U1(w_out,c_out), wc_U2(w_out,c_out), &
                & wc_b_1(w_in,c_in), wc_b_2(w_in,c_in), wc_b_1(w_out,c_out), wc_b_2(w_out,c_out), &
                & ws_normVec(nl,ns), nl, ns, w_in, c_in, w_out, c_out )

           ws_hFlux(nl,ns) = numFlux(1); 
           lsc_hU1Flux(nl,Face_EdgeId(1,ns),c_in) = numFlux(2); lsc_hU2Flux(nl,Face_EdgeId(1,ns),c_in) = numFlux(3)
           lsc_hU1Flux(nl,Face_EdgeId(2,ns),c_out) = numFlux(4); lsc_hU2Flux(nl,Face_EdgeId(2,ns),c_out) = numFlux(5)

           w_in = mod(w_in, 3*(nDGNodePerFace-1)) + 1
           w_out = w_out + 3*(nDGNodePerFace - 1) - 1
           if(w_out > 3*(nDGNodePerFace - 1)) w_out = mod(w_out, 3*(nDGNodePerFace-1))
        end do
     end if
  end do

  !$omp parallel do private(faceIds, nk, &
  !$omp & w_D_hFlux, w_hFlux, w_hU1Flux, w_D_hU1Flux, w_hU1Src, w_hU2Flux, w_D_hU2Flux, w_hU2Src, &
  !$omp & s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2) schedule(guided)
  do nc=1, nDGElement
     if(c_SeaCellFlag(nc)) then
        faceIds(:) = fvmInfo%Point_FaceId(3:1:-1, nc)
!!$    faceIds(:) = fvmInfo%Point_FaceId(:, nc)

        call SintPt_Flux(s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2, &
             & nc)

        w_hFlux(:) = matmul( wc_Ms(:,:,nc), reshape(ws_hFlux(:,faceIds(:)),  (/ 3*nDGNodePerFace /)) )
        w_hU1Flux(:) = matmul( wc_Ms(:,:,nc), reshape(lsc_hU1Flux(:,1:3,nc), (/ 3*nDGNodePerFace /)) )
        w_hU2Flux(:) = matmul( wc_Ms(:,:,nc), reshape(lsc_hU2Flux(:,1:3,nc), (/ 3*nDGNodePerFace /)) )

        do nk=1, nDGNodePerElem
           w_D_hFlux(nk) = TriNk_sinteg_dotProdWt( &
                & wc_Dy1(:,nk,nc)*s_hFlux1 + wc_Dy2(:,nk,nc)*s_hFlux2)

           w_D_hU1Flux(nk) = TriNk_sinteg_dotProdWt( &
                & wc_Dy1(:,nk,nc)*s_hU1Flux1 + wc_Dy2(:,nk,nc)*s_hU1Flux2)

           w_D_hU2Flux(nk) = TriNk_sinteg_dotProdWt( &
                & wc_Dy1(:,nk,nc)*s_hU2Flux1 + wc_Dy2(:,nk,nc)*s_hU2Flux2)
        end do


        call calc_hU1U2Src(w_hU1Src, w_hU2Src, &
             & s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2, &
             & nc )
     
        wc_dhdt(:,nc) = &
             &  - w_hFlux(:) + w_D_hFlux(:)

        wc_dhU1dt(:,nc) = &
             & - w_hU1Flux(:) + w_D_hU1Flux(:) + w_hU1Src(:)

        wc_dhU2dt(:,nc) = &
             & - w_hU2Flux(:) + w_D_hU2Flux(:) + w_hU2Src(:)

     end if
  end do

contains
subroutine calc_hU1U2Src(w_hU1Src, w_hU2Src, &
     & s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2, &
     & nc )

  use VariableSet_mod, only: &
       & wc_BtmTopl, &
       & wc_WindStress1, wc_WindStress2

  real(DP), dimension(nDGNodePerElem), intent(out) :: w_hU1Src, w_hU2Src
  real(DP), dimension(nDGSIntNodePerElem), intent(in) :: &
       & s_hFlux1, s_hFlux2, s_hU1Flux1, s_hU1Flux2, s_hU2Flux1, s_hU2Flux2
  integer, intent(in) :: nc

  real(DP), dimension(nDGSIntNodePerElem) :: s_y1, s_y2
  real(DP) :: s_ChristSym(nDGSIntNodePerElem, 2,2,2)
  real(DP), dimension(nDGSIntNodePerElem) :: s_Forcing1, s_Forcing2
  real(DP), dimension(nDGNodePerElem) :: w_BtmTopl
  real(DP) :: Gij(2,2)
  real(DP) :: Tau1, Tau2, BtmTopl_dy1, BtmTopl_dy2, toth

  integer :: nk
  type(vector3d) :: CoriolisForce, b1, b2, b_1, b_2, x_p
  integer :: j, k, ns
  real(DP) :: CoriPram
  type(Triangle) :: tri 

  tri = get_DGElementTri(nc)
  s_y1 = DGElemInfo%sIntNode_y1; s_y2 = DGElemInfo%sIntNode_y2

  do k=1, 2
     do j=1,2
        do ns=1, nDGSIntNodePerElem
           s_ChristSym(ns,1,j,k) = get_DGElemChristoffelSymbl(1,j,k,DGElemInfo%sIntNode(ns), tri)
           s_ChristSym(ns,2,j,k) = get_DGElemChristoffelSymbl(2,j,k,DGElemInfo%sIntNode(ns), tri)
        end do
     end do
  end do

  w_BtmTopl = wc_BtmTopl(:,nc)
  do ns=1, nDGSIntNodePerElem
     b_1 = get_DGElemCovariantBasis(1,DGElemInfo%sIntNode(ns),nc)
     b_2 = get_DGElemCovariantBasis(2,DGElemInfo%sIntNode(ns),nc)
     call get_DGElemContravariantBasis(b1, b2, b_1, b_2)
     Gij(:,:) = calc_Gij(calc_G_ij(b_1, b_2))
     x_p = mapping_local2globalCoord(DGElemInfo%sIntNode(ns), tri)

     !
     CoriPram = 2d0*Omega*x_p%v_(3)/l2norm(x_p)
     CoriolisForce =  CoriPram* &
          & (normalizedVec(x_p)).cross.(s_hFlux1(ns)*b_1 + s_hFlux2(ns)*b_2)

     !
     toth = TriNk_interpolate(s_y1(ns), s_y2(ns), wc_h(:,nc) + meanDepth)
     BtmTopl_dy1 = sum(w_BtmTopl*TriNk_basis_dy1(s_y1(ns),s_y2(ns)))
     BtmTopl_dy2 = sum(w_BtmTopl*TriNk_basis_dy2(s_y1(ns),s_y2(ns)))

     !
     Tau1 = TriNk_interpolate(s_y1(ns), s_y2(ns), wc_WindStress1(:,nc))
     Tau2 = TriNk_interpolate(s_y1(ns), s_y2(ns), wc_WindStress2(:,nc))

     !
     s_Forcing1(ns) = &
          &  - (   s_hU1Flux1(ns)*s_ChristSym(ns,1,1,1) + s_hU1Flux2(ns)*s_ChristSym(ns,1,2,1)   &
          &      + s_hU2Flux1(ns)*s_ChristSym(ns,1,1,2)                                      )   &
          &  - (CoriolisForce.dot.b1)                                   &
          &  + toth*Grav*(Gij(1,1)*BtmTopl_dy1 + Gij(1,2)*BtmTopl_dy2)  &
          &  - LinearDragCoef*s_hFlux1(ns)                              &
          &  + Tau1/RefDens*(toth/meanDepth)
     
     s_Forcing2(ns) = &
          &  -(   s_hU2Flux2(ns)*s_ChristSym(ns,2,2,2) + s_hU1Flux2(ns)*s_ChristSym(ns,2,2,1)   &
          &     + s_hU2Flux1(ns)*s_ChristSym(ns,2,1,2)                                      )   &
          &  - (CoriolisForce.dot.b2)                                   &
          &  + toth*Grav*(Gij(2,1)*BtmTopl_dy1 + Gij(2,2)*BtmTopl_dy2)  &
          &  - LinearDragCoef*s_hFlux2(ns)                              &
          &  + Tau2/RefDens*(toth/meanDepth)
  end do

  do nk=1, nDGNodePerElem
     w_hU1Src(nk) = TriNk_sinteg_dotProdWt( wc_S(:,nk,nc)*s_Forcing1(:) )
     w_hU2Src(nk) = TriNk_sinteg_dotProdWt( wc_S(:,nk,nc)*s_Forcing2(:) )
  end do


end subroutine calc_hU1U2Src

subroutine SintPt_Flux( &
     & sintPt_hFlux1, sintPt_hFlux2, sintPt_hU1Flux1, sintPt_hU1Flux2, sintPt_hU2Flux1, sintPt_hU2Flux2, &
     & nc)

  real(DP), dimension(nDGSIntNodePerElem), intent(out) :: &
       & sintPt_hFlux1, sintPt_hFlux2, &
       & sintPt_hU1Flux1, sintPt_hU1Flux2, sintPt_hU2Flux1, sintPt_hU2Flux2
  integer, intent(in) :: nc

  integer :: i, nk
  real(DP), dimension(nDGSIntNodePerElem) :: &
       & s_hU1, s_hU2, s_Phi
  real(DP) :: G_ij(2,2), s_Gij(nDGSIntNodePerElem,2,2)
  real(DP), dimension(nDGSIntNodePerElem) :: s_y1, s_y2
  
  do nk=1, nDGSIntNodePerElem
     s_Gij(nk,:,:) = calc_Gij(calc_G_ij( &
          & get_DGElemCovariantBasis(1,DGElemInfo%sIntNode(nk),nc), &
          & get_DGElemCovariantBasis(2,DGElemInfo%sIntNode(nk),nc) ))
  end do

  s_y1 = DGElemInfo%sIntNode_y1;   s_y2 = DGElemInfo%sIntNode_y2
  s_hU1(:) = TriNk_interpolate(s_y1, s_y2, wc_hU1(:,nc))
  s_hU2(:) = TriNk_interpolate(s_y1, s_y2, wc_hU2(:,nc))
  s_Phi(:) = TriNk_interpolate(s_y1, s_y2, wc_h(:,nc) + meanDepth)

  sintPt_hFlux1 = s_hU1
  sintPt_hFlux2 = s_hU2
  sintPt_hU1Flux1 = s_hU1**2/s_Phi + 0.5d0*Grav*s_Gij(:,1,1)*s_Phi**2
  sintPt_hU1Flux2 = s_hU1*s_hU2/s_Phi + 0.5d0*Grav*s_Gij(:,1,2)*s_Phi**2
  sintPt_hU2Flux1 = sintPt_hU1Flux2
  sintPt_hU2Flux2 = s_hU2**2/s_Phi + 0.5d0*Grav*s_Gij(:,2,2)*s_Phi**2

end subroutine SintPt_Flux

end subroutine calcTendency


function calc_NumericFlux( &
     & h_in, U1_in, U2_in, h_out, U1_out, U2_out, &
     & b_1_in, b_2_in, b_1_out, b_2_out, normVec,      &
     & nl, ns, w_in, c_in, w_out, c_out ) result(numFlux)

use SphericalCoord_mod

  real(DP), intent(in) :: h_in, U1_in, U2_in
  real(DP), intent(in) :: h_out, U1_out, U2_out
  type(vector3d), intent(in) :: b_1_in, b_2_in
  type(vector3d), intent(in) :: b_1_out, b_2_out  
  type(vector3d), intent(in) :: normVec
  integer :: nl, ns, w_in, w_out, c_in, c_out
  real(DP) :: numFlux(5)

  type(vector3d) :: b1_in, b2_in, b1_out, b2_out
  real(DP) :: numFlux_in(3), numFlux_out(3)
  real(DP) :: toth(2), U1_basisIn(2), U1_basisOut(2), U2_basisIn(2), U2_basisOut(2)
  real(DP) :: in2outBasisMat(2,2), out2inBasisMat(2,2)

  real(DP) :: Un, Ut
  type(vector3d) :: tanVec

  !
  call get_DGElemContravariantBasis(b1_in, b2_in, b_1_in, b_2_in)
  call get_DGElemContravariantBasis(b1_out, b2_out, b_1_out, b_2_out)

  in2outBasisMat(1,1) = b1_out.dot.b_1_in; in2outBasisMat(1,2) = b1_out.dot.b_2_in;
  in2outBasisMat(2,1) = b2_out.dot.b_1_in; in2outBasisMat(2,2) = b2_out.dot.b_2_in;

  out2inBasisMat(1,1) = b1_in.dot.b_1_out; out2inBasisMat(1,2) = b1_in.dot.b_2_out;
  out2inBasisMat(2,1) = b2_in.dot.b_1_out; out2inBasisMat(2,2) = b2_in.dot.b_2_out;


  !
  toth(:) = meanDepth + (/ h_in, h_out /); 
  U1_basisIn(:) = (/ U1_in, U1_out*out2inBasisMat(1,1) + U2_out*out2inBasisMat(1,2) /)
  U2_basisIn(:) = (/ U2_in, U1_out*out2inBasisMat(2,1) + U2_out*out2inBasisMat(2,2) /)
  U1_basisOut(:) = (/ U1_in*in2outBasisMat(1,1) + U2_in*in2outBasisMat(1,2), U1_out /)
  U2_basisOut(:) = (/ U1_in*in2outBasisMat(2,1) + U2_in*in2outBasisMat(2,2), U2_out /)

  if( abs(wc_DGNodeUsage(w_in,c_in))==10 .and. abs(wc_DGNodeUsage(w_out, c_out))==10 &
      .and. wc_DGNodeUsage(w_in,c_in)*wc_DGNodeUsage(w_out,c_out) < 0 &
       & )then

     if(wc_DGNodeUsage(w_in,c_in)==-10) toth(1) = toth(2)
     if(wc_DGNodeUsage(w_out,c_out)==-10) toth(2) = toth(1)

!!$write(*,*) wc_DGNodeUsage(w_in,c_in), wc_DGNodeUsage(w_out,c_out), &
!!$     & ":", c_SeaCellFlag(c_in), c_SeaCellFlag(c_out)

!     if(wc_DGNodeUsage(w_in,c_in)==-10) then
     tanVec = normalizedVec(wc_DGNodePos(w_in,c_in).cross.normVec)
     Un =  U1_out*(b_1_out.dot.normVec) + U2_out*(b_2_out.dot.normVec)
     Ut =  U1_out*(b_1_out.dot.tanVec) + U2_out*(b_2_out.dot.tanVec)
     U1_basisOut(1) = -Un*(normVec.dot.b1_out) + Ut*(tanVec.dot.b1_out)
     U2_basisOut(1) = -Un*(normVec.dot.b2_out) + Ut*(tanVec.dot.b2_out)
!     end if
!!$     U1_basisOut(1) = - U1_basisOut(2)
!!$     U2_basisOut(1) = - U2_basisOut(2)
!     if(wc_DGNodeUsage(w_out,c_out)==-10) then
     tanVec = normalizedVec(wc_DGNodePos(w_out,c_out).cross.normVec)
     Un =  U1_in*(b_1_in.dot.normVec) + U2_in*(b_2_in.dot.normVec)
     Ut =  U1_in*(b_1_in.dot.tanVec) + U2_in*(b_2_in.dot.tanVec)
     U1_basisIn(2) = -Un*(normVec.dot.b1_in) + Ut*(tanVec.dot.b1_in)
     U2_basisIn(2) = -Un*(normVec.dot.b2_in) + Ut*(tanVec.dot.b2_in)
!     end if
!!$     U1_basisIn(2) = - U1_basisIn(1)
!!$     U2_basisIn(2) = - U2_basisIn(1)

!!$     if(wc_DGNodeUsage(w_out,c_out)==-10) then
!!$        write(*,*) c_out, w_out, ":", Un, Ut
!!$     end if

  end if

  tanVec = normalizedVec(wc_DGNodePos(w_out,c_out).cross.normVec) 
  numFlux_in(:) = get_RusanovFlux( &
!!$  numFlux_in(:) = get_RoeFlux( &
       & toth, U1_basisIn, U2_basisIn, b_1_in, b_2_in, b1_in, b2_in, normVec, tanVec )
  numFlux_out(:) = get_RusanovFlux( &
!!$  numFlux_out(:) = get_RoeFlux( &
       & toth, U1_basisOut, U2_basisOut, b_1_out, b_2_out, b1_out, b2_out, normVec, tanVec )
  numFlux(:) = (/ numFlux_in(1), numFlux_in(2), numFlux_in(3), numFlux_out(2), numFlux_out(3) /)

!!$  if(abs(wc_DGNodeUsage(w_in,c_in))==10 .and. abs(wc_DGNodeUsage(w_out, c_out))==10)then
!!$     if(wc_DGNodeUsage(w_in,c_in)==-10) numFlux(1) = numFlux_out(1)
!!$  end if
contains 
function get_RusanovFlux(toth, U1, U2, b_1, b_2, b1, b2, edgeNormal, edgeTangen) result(rusanovFlux)
  real(DP), intent(in) :: toth(2), U1(2), U2(2)
  type(vector3d), intent(in) :: b_1, b_2, b1, b2, edgeNormal, edgeTangen
  real(DP) :: rusanovFlux(3)

  real(DP) :: tmpFlux(2,3), lambda
  real(DP) :: Gij(2,2), orthFac(2,2), deorthFac(2,2), Un(2), Ut(2)
  real(DP) :: rusanovFlux_(3)

  orthFac(1,1) = b_1.dot.edgeNormal; orthFac(1,2) = b_2.dot.edgeNormal
  orthFac(2,1) = b_1.dot.edgeTangen; orthFac(2,2) = b_2.dot.edgeTangen

  deorthFac(1,1) = b1.dot.edgeNormal; deorthFac(1,2) = b1.dot.edgeTangen
  deorthFac(2,1) = b2.dot.edgeNormal; deorthFac(2,2) = b2.dot.edgeTangen

!  Gij = calc_Gij(calc_G_ij(b_1, b_2))
  Un = orthFac(1,1)*U1 + orthFac(1,2)*U2
  Ut = orthFac(2,1)*U1 + orthFac(2,2)*U2

!!$  tmpFlux(:,1) = toth*Un
!!$  tmpFlux(:,2) = toth*U1*Un &
!!$       & + 0.5d0*Grav*toth**2*(orthFac(1)*Gij(1,1) + orthFac(2)*Gij(1,2))
!!$  tmpFlux(:,3) = toth*U2*Un &
!!$       & + 0.5d0*Grav*toth**2*(orthFac(1)*Gij(2,1) + orthFac(2)*Gij(2,2))

  tmpFlux(:,1) = toth*Un
  tmpFlux(:,2) = toth*Un*Un + 0.5d0*Grav*toth**2
  tmpFlux(:,3) = toth*Un*Ut

!  lambda = max( abs(Un(1)), abs(Un(2)) )
!!$  lambda = max( abs(Un(1))+sqrt(Grav*toth(1)), abs(Un(2))+sqrt(Grav*toth(2)) )
!!$  rusanovFlux = 0.5d0*( &
!!$       &   (tmpFlux(1,:)+tmpFlux(2,:)) &
!!$       & - lambda*( (/ toth(2),toth(2)*U1(2),toth(2)*U2(2) /) - (/ toth(1),toth(1)*U1(1),toth(1)*U2(1) /) ) &
!!$       & )

  lambda = max( abs(Un(1))+sqrt(Grav*toth(1)), abs(Un(2))+sqrt(Grav*toth(2)) )
  rusanovFlux_ = 0.5d0*( &
       &   (tmpFlux(1,:)+tmpFlux(2,:)) &
       & - lambda*( (/ toth(2),toth(2)*Un(2),toth(2)*Ut(2) /) - (/ toth(1),toth(1)*Un(1),toth(1)*Ut(1) /) ) &
       & )

  rusanovFlux(1) = rusanovFlux_(1)
  rusanovFlux(2) = deorthFac(1,1)*rusanovFlux_(2) + deorthFac(1,2)*rusanovFlux_(3)
  rusanovFlux(3) = deorthFac(2,1)*rusanovFlux_(2) + deorthFac(2,2)*rusanovFlux_(3)

end function get_RusanovFlux

function get_RoeFlux(toth, U1, U2, b_1, b_2, b1, b2, edgeNormal, edgeTangen) result(roeFlux)

  real(DP), intent(in) :: toth(2), U1(2), U2(2)
  type(vector3d), intent(in) :: b_1, b_2, b1, b2, edgeNormal, edgeTangen
  real(DP) :: roeFlux(3)

  real(DP) :: tmpFlux(2,3)
  real(DP) :: Fr, UnA, UtA, cA, w
  real(DP) :: Gij(2,2), orthFac(2,2), deorthFac(2,2), Un(2), Ut(2)
  real(DP) :: roeFlux_(3)
  
  orthFac(1,1) = b_1.dot.edgeNormal; orthFac(1,2) = b_2.dot.edgeNormal
  orthFac(2,1) = b_1.dot.edgeTangen; orthFac(2,2) = b_2.dot.edgeTangen

  deorthFac(1,1) = b1.dot.edgeNormal; deorthFac(1,2) = b1.dot.edgeTangen
  deorthFac(2,1) = b2.dot.edgeNormal; deorthFac(2,2) = b2.dot.edgeTangen

  Un = orthFac(1,1)*U1 + orthFac(1,2)*U2
  Ut = orthFac(2,1)*U1 + orthFac(2,2)*U2

  w = sqrt(toth(2)/toth(1))
  UnA = (Un(1) + w*Un(2))/(1d0 + w)
  UtA = (Ut(1) + w*Ut(2))/(1d0 + w)
  cA = sqrt(Grav*0.5d0*(toth(1) + toth(2)))
  Fr = UnA/cA

  tmpFlux(:,1) = toth*Un
  tmpFlux(:,2) = toth*Un*Un + 0.5d0*Grav*toth**2
  tmpFlux(:,3) = toth*Un*Ut

  cA = max( abs(Un(1))+sqrt(Grav*toth(1)), abs(Un(2))+sqrt(Grav*toth(2)) )
  roeFlux_(:) = 0.5d0*( &
       &       (tmpFlux(1,:)+tmpFlux(2,:)) &
       &  + Fr*(tmpFlux(1,:)-tmpFlux(2,:)) &
       &  + cA*(1d0 - Fr**2)*( &
       &    (/ toth(1),toth(1)*Un(1),toth(1)*Ut(1) /)-(/ toth(2),toth(2)*Un(2),toth(2)*Ut(2) /) &
       &  ) &
       & )

  roeFlux(1) = roeFlux_(1)
  roeFlux(2) = deorthFac(1,1)*roeFlux_(2) + deorthFac(1,2)*roeFlux_(3)
  roeFlux(3) = deorthFac(2,1)*roeFlux_(2) + deorthFac(2,2)*roeFlux_(3)

end function get_RoeFlux

end function calc_NumericFlux


end module GovernEquationSolver_mod
