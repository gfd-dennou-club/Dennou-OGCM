!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Boundary_spm_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod, only: &
       & RefDens, RefSalt, Cp0
    

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax

  use SpmlUtil_mod
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Admin_BC_mod, only: &
       & inquire_VBCSpecType,                                                    &
       & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress,                 &
       & DynBCTYPE_LinFric,                                                      &
       & ThermBCTYPE_PrescFixedFlux, ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, &
       & ThermBCTYPE_PresFlux_Han1984Method,                                     &
       & ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed,                         & 
       & SaltBCTYPE_PrescFixedFlux, SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat,    &
       & SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed,                           &
       & SaltBCTYPE_PresFlux_Han1984Method,                                      &
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface,          &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom,              &
       & SeaSfcTempRelaxedTime, SeaSfcSaltRelaxedTime

  use DOGCM_Boundary_vars_mod, only: &
       & xy_SfcHFlx_ns, xy_SfcHFlx_sr,        &
       & xy_WindStressU, xy_WindStressV,      &
       & xy_FreshWtFlx, xy_FreshWtFlxS,       &
       & xy_SeaSfcTemp, xy_SeaSfcSalt,        &
       & xy_OcnSfcCellMask, xyz_OcnCellMask,  &
       & OCNCELLMASK_SICE, OCNCELLMASK_OCEAN
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_Boundary_spm_Init, DOGCM_Boundary_spm_Final

  public :: DOGCM_Boundary_spm_ApplyBC
  
  public :: DOGCM_Boundary_spm_InqVBCRHS_UV
  public :: DOGCM_Boundary_spm_InqVBCRHS_TRC
  
  ! 公開変数
  ! Public variable
  !

  !< The depth of mixed layer near sea surface to specify Haney-type boundary condition.
  real(DP), parameter :: MixLyrDepthConst = 50d0
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Boundary_spm_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Boundary_spm_Init( &
       & configNmlName )                   ! (in)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    
  end subroutine DOGCM_Boundary_spm_Init

  !>
  !!
  !!
  subroutine DOGCM_Boundary_spm_Final()

    ! 実行文; Executable statements
    !


  end subroutine DOGCM_Boundary_spm_Final

  !-----------------------------------------
  
  subroutine DOGCM_Boundary_spm_ApplyBC(    &
       & xyz_U, xyz_V, xyza_TRC,                    & ! (inout)
       & xyz_H, xyz_VViscCoef, xyz_VDiffCoef        & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_U(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V(IA,JA,KA)
    real(DP), intent(inout) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    
    real(DP) :: xya_U_VBCRHS(IA,JA,2)
    real(DP) :: xya_V_VBCRHS(IA,JA,2)
    real(DP) :: xya_PTemp_VBCRHS(IA,JA,2)
    real(DP) :: xya_Salt_VBCRHS(IA,JA,2)
    real(DP) :: avr_salt
    
    ! 実行文; Executable statements
    !

    call DOGCM_Boundary_spm_InqVBCRHS_UV( &
         & xya_U_VBCRHS, xya_V_VBCRHS,                 & ! (out)
         & xyz_U, xyz_V, xyz_H, xyz_VViscCoef          & ! (in)
         & )
    
    call solve_VBCEq( xyz_U,                           & ! (inout)
         & DynBC_Surface, DynBC_Bottom, xya_U_VBCRHS   & ! (in)
         & )

    call solve_VBCEq( xyz_V,                           & ! (inout)
         & DynBC_Surface, DynBC_Bottom, xya_V_VBCRHS   & ! (in)
         & )

    call DOGCM_Boundary_spm_InqVBCRHS_TRC( &
         & xya_PTemp_VBCRHS, xya_Salt_VBCRHS,          & ! (out)
         & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),  & ! (in)
         & xyz_H, xyz_VDiffCoef                                      & ! (in)
         & )
    
    call solve_VBCEq( xyza_TRC(:,:,:,TRCID_PTEMP),            & ! (inout)
         & ThermBC_Surface, ThermBC_Bottom, xya_PTemp_VBCRHS  & ! (in)
         & )


!!$    avr_salt = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyza_TRC(:,:,:,TRCID_SALT)))
    call solve_VBCEq( xyza_TRC(:,:,:,TRCID_SALT),         & ! (inout)
         & SaltBC_Surface, SaltBC_Bottom, xya_Salt_VBCRHS     & ! (in)
         & )
!!$   write(*,*) "After BC mod:",  z_DSig_z( xyza_TRC(0,jMax/2,:,TRCID_SAlT) )
!!$   write(*,*) ":-> ", (- avr_salt + AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyza_TRC(:,:,:,TRCID_SALT))))

  end subroutine DOGCM_Boundary_spm_ApplyBC
       
  !-----------------------------

  subroutine DOGCM_Boundary_spm_InqVBCRHS_TRC( &
       & xya_PTemp_VBCRHS, xya_Salt_VBCRHS,            & ! (out)
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_VDiffCoef     & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xya_PTemp_VBCRHS(IA,JA,2)
    real(DP), intent(out) :: xya_Salt_VBCRHS(IA,JA,2)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)


    ! 局所変数
    ! Local variables
    !
    integer :: k
    
    ! RHS for thermal boundary condition
    !

    select case(ThermBC_Surface)
    case(  ThermBCTYPE_PrescTemp                              &
         & )
       xya_PTemp_VBCRHS(:,:,1) = xy_SeaSfcTemp
       
    case(  ThermBCTYPE_PrescFixedFlux, ThermBCTYPE_PrescFlux, &
         & ThermBCTYPE_Adiabat, ThermBCTYPE_TempRelaxed,      &
         & ThermBCTYPE_PresFlux_Han1984Method                 &
         & )
       !$omp parallel
       !$omp workshare
       xya_PTemp_VBCRHS(:,:,1) = &
            & - xyz_H(:,:,KS)/(RefDens*Cp0*xyz_vDiffCoef(:,:,KS))*( &
            &      xy_SfcHFlx_sr + xy_SfcHFlx_ns                  &
            & )
       !$omp end workshare
       !$omp end parallel
       
    case default 
       call throw_UnImplementVBCError('ThermBC_Surface', ThermBC_Surface)
    end select

    select case(ThermBC_Bottom)
    case(ThermBCTYPE_Adiabat)
       xya_PTemp_VBCRHS(:,:,2) = 0d0
    case default
       call throw_UnImplementVBCError('ThermBC_Bottom', ThermBC_Bottom)
    end select

    
    ! RHS for boundary condition of Salinity
    !
    
    select case(SaltBC_Surface)
    case(  SaltBCTYPE_PrescSalt                             &
         & )
       xya_Salt_VBCRHS(:,:,1) = xy_SeaSfcSalt

    case( SaltBCTYPE_Adiabat )
       xya_Salt_VBCRHS(:,:,1) = 0d0
    case(  SaltBCTYPE_PrescFixedFlux, SaltBCTYPE_PrescFlux, &
         & SaltBCTYPE_PresFlux_Han1984Method,               &
         & SaltBCTYPE_SaltRelaxed                           &
         & )
       xya_Salt_VBCRHS(:,:,1) = - xyz_H(:,:,KS)/xyz_VDiffCoef(:,:,KS)*( &
            & xy_FreshWtFlxS * RefSalt )
    case default 
       call throw_UnImplementVBCError('SaltBC_Surface', SaltBC_Surface)       
    end select

    !
    select case(SaltBC_Bottom)
    case(SaltBCTYPE_Adiabat)
       xya_Salt_VBCRHS(:,:,2) = 0d0
    case default 
       call throw_UnImplementVBCError('SaltBC_Bottom', SaltBC_Bottom)       
    End select
    
  end subroutine DOGCM_Boundary_spm_InqVBCRHS_TRC
  
  subroutine DOGCM_Boundary_spm_InqVBCRHS_UV( &
       & xya_U_VBCRHS, xya_V_VBCRHS,                   & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_VViscCoef            & ! (in)
       & )

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xya_U_VBCRHS(IA,JA,2)
    real(DP), intent(out) :: xya_V_VBCRHS(IA,JA,2)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    integer :: k
    real(DP) :: xy_Coef(IA,JA)
    
    ! 実行文; Executable statement
    !

    ! -- RHS for dynamical boundary condition

    select case(DynBC_Surface)
    case(  DynBCTYPE_SpecStress,               &
         & DynBCTYPE_Slip                      &
         & )
       !$omp parallel
       !$omp workshare
       xya_U_VBCRHS(:,:,1) = xyz_H(:,:,KS)/(RefDens*xyz_VViscCoef(:,:,KS))*xy_WindStressU
       xya_V_VBCRHS(:,:,1) = xyz_H(:,:,KS)/(RefDens*xyz_VViscCoef(:,:,KS))*xy_WindStressV
       !$omp end workshare
       !$omp end parallel
       
    case( DynBCTYPE_NoSlip )
       xya_U_VBCRHS(:,:,1) = 0d0
       xya_V_VBCRHS(:,:,1) = 0d0
    case default
       call throw_UnImplementVBCError('DynBC_Surface', DynBC_Surface)              
    end select
    
!!$    xy_Coef(:,:) = xyz_H(:,:,KE)/(RefDens*xyz_VViscCoef(:,:,KE))    
    select case(DynBC_Bottom)
    case(DynBCTYPE_NoSlip)
       xya_U_VBCRHS(:,:,2) = 0d0
       xya_V_VBCRHS(:,:,2) = 0d0       
    case(DynBCTYPE_Slip)
       xya_U_VBCRHS(:,:,2) = 0d0
       xya_V_VBCRHS(:,:,2) = 0d0
    case default
       call throw_UnImplementVBCError('DynBC_Bottom', DynBC_Bottom)       
    end select

  end subroutine DOGCM_Boundary_spm_InqVBCRHS_UV

  subroutine throw_UnImplementVBCError(boundaryLabel, BCTypeID)
    character(*), intent(in) :: boundaryLabel
    integer, intent(in) :: BCTypeID

    call MessageNotify( 'E', module_name, &
         & "'%a=%d' has not been implemeneted yet.", &
         & ca=(/ boundaryLabel /), i=(/ BCTypeID /)  &
         & )
  end subroutine throw_UnImplementVBCError
  
  !- Private subroutines -------------------------------------------------

  
  subroutine solve_VBCEq( xyz,                       & ! (inout)
       & SurfBoundaryID, BtmBoundaryID, xya_VBCRHS   & ! (in)
       & )
    real(DP), intent(inout) :: xyz(IA,JA,KA)
    integer, intent(in) :: SurfBoundaryID
    integer, intent(in) :: BtmBoundaryID
    real(DP), intent(in) :: xya_VBCRHS(IA,JA,2)

    character :: SurfVBCType
    character :: BtmVBCType
    real(DP) :: VBCMat(KS:KE,KS:KE)
    real(DP) :: DSigMat(KS:KE,KS:KE)
    real(DP) :: IMat(KS:KE,KS:KE)    
    integer :: i, j, k, n 
    integer :: IPiv(KS:KE)
    real(DP) :: b(KS:KE)
    integer :: info

    real(DP) :: BCMat(2,2)
    real(DP) :: BCInvMat(2,2)
    
    real(DP) :: a1, a2
    real(DP) :: b1, b2
    real(DP) :: RHS1, RHS2
    
    SurfVBCType = inquire_VBCSpecType(SurfBoundaryID)
    BtmVBCType = inquire_VBCSpecType(BtmBoundaryID)
    
    if (SurfVBCType//BtmVBCType == 'DD') then
       xyz(:,:,KS) = xya_VBCRHS(:,:,1)
       xyz(:,:,KE) = xya_VBCRHS(:,:,2)
       return
    end if

    !
    IMat(:,:) = 0d0
    forAll(k=KS:KE) IMat(k,k) = 1d0
    do k = KS, KE
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do

!!$    select case(SurfVBCType)
!!$    case('D')
!!$       a1 = 1; a2 = 0
!!$    case('N')
!!$       a1 = 0; a2 = 1
!!$    end select
!!$
!!$    select case(BtmVBCType)
!!$    case('D')
!!$       b1 = 1; b2 = 0
!!$    case('N')
!!$       b1 = 0; b2 = 1
!!$    end select
!!$
!!$    BCMat(1,:) = (/ a1 + a2*DSigMat(0,0), a2*DSigMat(0,kMax) /)
!!$    BCMat(2,:) = (/ b2*DSigMat(kMax,0), b1 + b2*DSigMat(kMax,kMax) /)
!!$
!!$    BCInvMat(1,:) = (/ BCMat(2,2), -BCMat(1,2) /)
!!$    BCInvMat(2,:) = (/ -BCMat(2,1), BCMat(1,1) /)
!!$    BCInvMat(:,:) = BCInvMat(:,:)/(BCMat(1,1)*BCMat(2,2) - BCMat(1,2)*BCMat(2,1))
!!$
!!$    do j=1, jMax
!!$       do i=0, iMax-1
!!$          RHS1 = xya_VBCRHS(i,j,1) - a2*sum(DSigMat(0,1:kMax-1)*xyz(i,j,1:kMax-1))
!!$          RHS2 = xya_VBCRHS(i,j,2) - b2*sum(DSigMat(kMax,1:kMax-1)*xyz(i,j,1:kMax-1))
!!$          xyz(i,j,0) = BCInvMat(1,1)*RHS1 + BCInvMat(1,2)*RHS2
!!$          xyz(i,j,kMax) = BCInvMat(2,1)*RHS1 + BCInvMat(2,2)*RHS2          
!!$       end do
!!$    end do
!!$    return
    
    !--------------
    
    VBCMat(:,:) = IMat
    if(SurfVBCType == 'N') VBCMat(KS,:) = DSigMat(KS,:)
    if(BtmVBCType == 'N') VBCMat(KE,:) = DSigMat(KE,:)
    n = size(VBCMat, 1)

    call DGETRF(n, n, VBCMat, n, IPiv, info)

    !$omp parallel do private(i, b, info)
    do j=JS, JE
       do i=IS, IE
          b(:) = xyz(i,j,KS:KE);
          b(KS) = xya_VBCRHS(i,j,1); b(KE) = xya_VBCRHS(i,j,2) 

          call DGETRS('N', n, 1, VBCMat, n, IPiv, b, n, info)

          xyz(i,j,KS) = b(KS); xyz(i,j,KE) = b(KE)
       end do
    end do
    
  end subroutine solve_VBCEq
  
end module DOGCM_Boundary_spm_mod
