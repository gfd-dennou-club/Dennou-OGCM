!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module BoundaryCondO_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_historyauto
  
  !* Dennou-OGCM

  use SpmlUtil_mod
  
  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType, &
       & ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed, & 
       & SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed, & 
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom 

  use GridSet_mod, only: &
       & iMax, jMax, kMax


  use VariableSet_mod, only: &
       & xy_totDepthBasic, z_PTempBasic
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BoundaryCondO_Init, BoundaryCondO_Final
  public :: BoundaryCondO_Update

  public :: apply_VBoundaryCondO
  public :: calc_SurfaceHeatFluxO
  
  ! 公開変数
  ! Public variable
  !

  character(*), parameter :: VARSET_KEY_SURFHFLXO   = 'SurfHFlxO'
  character(*), parameter :: VARSET_KEY_SURFFWFLXO  = 'SurfFwFlxO' 
  
  !
  real(DP), public, save, allocatable :: xy_WindStressU(:,:)

  !
  real(DP), public, save, allocatable :: xy_WindStressV(:,:)

  !
  real(DP), public, save, allocatable :: xy_SeaSurfAirTemp(:,:), xy_SeaSurfTemp(:,:)

  !
  real(DP), public, save, allocatable :: xy_SeaSurfSalt(:,:)

  !< Heat flux at ocean surface (*upward positive*) [W/m2]
  real(DP), public, save, allocatable :: xy_SurfHFlxO(:,:)

  real(DP), public, save, allocatable :: xy_SurfHFlxIO(:,:)
  
  !< Freshwater flux at ocean surface (*upward positive*) [W/m2]   
  real(DP), public, save, allocatable :: xy_SurfFwFlxO(:,:)

  real(DP), public, save, allocatable :: xy_SurfFwFlxIO(:,:)

  !
  real(DP), public, save, allocatable :: xy_Wrain(:,:)
  real(DP), public, save, allocatable :: xy_Wsnow(:,:)
  real(DP), public, save, allocatable :: xy_Wevap(:,:)
  real(DP), public, save, allocatable :: xy_LatentDWHFlx(:,:)
  real(DP), public, save, allocatable :: xy_SensDWHFlx(:,:)
  real(DP), public, save, allocatable :: xy_SWDWRFlx(:,:)
  real(DP), public, save, allocatable :: xy_LWDWRFlx(:,:)

  !
  real(DP), parameter, public :: RefSalt_VBC = 35d0
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BoundaryCondO_mod' !< Module Name
  logical :: outputSurfFlxFlag
  
contains

  !>
  !!
  !!
  subroutine BoundaryCondO_Init(isSurfFlxOutput)

    ! 宣言文; Declareration statemenets
    !
    logical, intent(in), optional :: isSurfFlxOutput
    
    ! 実行文; Executable statements
    !

    ! Initialization
    !

    outputSurfFlxFlag = .false.
    
    ! Initialize arrays to store variables associated with  surface fluxes/
    call malloc2DVar(xy_WindStressU); call malloc2DVar(xy_WindStressV)
    call malloc2DVar(xy_SeaSurfTemp); call malloc2DVar(xy_SeaSurfAirTemp); call malloc2DVar(xy_SeaSurfSalt)
    call malloc2DVar(xy_SurfHFlxO); call malloc2DVar(xy_SurfFwFlxO)

    call malloc2DVar(xy_Wrain); call malloc2DVar(xy_Wsnow); call malloc2DVar(xy_Wevap)
    call malloc2DVar(xy_LatentDWHFlx); call malloc2DVar(xy_SensDWHFlx);
    call malloc2DVar(xy_SWDWRFlx); call malloc2DVar(xy_LWDWRFlx)

    call malloc2DVar(xy_SurfHFlxIO); call malloc2DVar(xy_SurfFwFlxIO)

    xy_SurfHFlxIO = 0d0; xy_SurfFwFlxO = 0d0

    ! Preparation of output surface fluxes
    !
    
    if(present(isSurfFlxOutput)) outputSurfFlxFlag = .true.
    if(outputSurfFlxFlag) then
       call prepare_output()
    end if
    
  contains
    subroutine malloc2DVar( var )
      real(DP), intent(inout), allocatable :: var(:,:)
      allocate( var(0:iMax-1, 1:jMax) )
    end subroutine malloc2DVar    
  end subroutine BoundaryCondO_Init

  !>
  !!
  !!
  subroutine BoundaryCondO_Final()

    ! 実行文; Executable statements
    !

    deallocate( xy_WindStressU, xy_WindStressV )
    deallocate( xy_SeaSurfTemp, xy_SeaSurfAirTemp, xy_SeaSurfSalt )
    deallocate( xy_SurfHFlxO, xy_SurfFwFlxO )

    deallocate( xy_Wrain, xy_Wsnow, xy_Wevap )
    deallocate( xy_LatentDWHFlx, xy_SensDWHFlx, xy_SWDWRFlx, xy_LWDWRFlx )

    deallocate( xy_SurfHFlxIO, xy_SurfFwFlxIO )
    
  end subroutine BoundaryCondO_Final

  !> @brief 
  !!
  !!
  subroutine BoundaryCondO_Update( xy_SIceCon, xy_Wice )

    ! モジュール引用; Use statements
    !
    use VariableSet_mod, only: &
         & xyz_PTempEddN, z_PTempBasic
    
    use TemporalIntegSet_mod, only: &
         & CurrentTime
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_SIceCon, xy_Wice
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DzT(0:iMax-1,jMax,0:kMax)
    
    ! 実行文; Executable statement
    !

    ! Update boundary conditions at sea surface and bottom.
    !
    
    select case(ThermBC_Surface)
    case(ThermBCTYPE_PrescFlux)
       call calc_SurfaceHeatFluxO( xy_SurfHFlxO,        & ! (out)
            & z_PTempBasic(0) + xyz_PTempEddN(:,:,0),   & ! (in)
            & xy_SIceCon                                & ! (in)
            & )
    case default
       call MessageNotify('E', module_name, &
            & 'Specified TermBC_Surface ID(=%a) is invalid.', i=(/ThermBC_Surface/))
    end select

    select case(SaltBC_Surface)
    case(SaltBCTYPE_PrescFlux)
       call calc_SurfaceFreshWaterFluxO( xy_SurfFwFlxO,   & ! (out)
            & xy_SIceCon, xy_Wice )
    case default
       call MessageNotify('E', module_name, &
            & 'Specified SaltBC_Surface ID(=%a) is invalid.', i=(/SaltBC_Surface/))       
    end select

    ! Output surface fluxes
    !
    if(outputSurfFlxFlag) then
       call output_surfaceFlux(CurrentTime)
    end if
    
  end subroutine BoundaryCondO_Update

  !> @brief 
  !!
  !!
  subroutine calc_SurfaceHeatFluxO( xy_SurfHFlxO,                & ! (out) 
    & xy_SurfTemp, xy_SIceCon                                    & ! (in)
    & )

    ! モジュール引用; use statements
    !
    use Constants_mod, only: &
         & StB, &
         & AlbedoOcean, EmissivOcean
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: xy_SurfHFlxO
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_SurfTemp, xy_SIceCon
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax) :: xy_SurfDWHFlx_AO

    ! 実行文; Executable statement
    !
    
    !$omp parallel workshare
    xy_SurfDWHFlx_AO(:,:) = &
         &   xy_LatentDWHFlx                        &
         & + xy_SensDWHFlx                          &
         & + emissivOcean*xy_LWDWRFlx               &
         & + (1d0 - albedoOcean)*xy_SWDWRFlx        &
         & - emissivOcean*(StB*xy_SurfTemp**4)

    xy_SurfHFlxO(:,:) = -(1d0 - xy_SIceCon)*xy_SurfDWHFlx_AO + xy_SIceCon*xy_SurfHFlxIO
    !$omp end  parallel workshare
    
  end subroutine calc_SurfaceHeatFluxO

  !> @brief 
  !!
  !!
  subroutine calc_SurfaceFreshWaterFluxO( xy_SurfFwFlxO,   & ! (out)
    & xy_SIceCon, xy_Wice )
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax), intent(out) :: xy_SurfFwFlxO
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: &
         & xy_SIceCon, xy_Wice
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    !$omp parallel workshare
    xy_SurfFwFlxO = &
         &   (1d0 - xy_SIceCon)*((xy_Wrain + xy_Wsnow) - xy_Wevap) &
         & + xy_SIceCon*(xy_Wrain - xy_Wice)
    !$omp end parallel workshare
    
  end subroutine calc_SurfaceFreshWaterFluxO

  !> @brief
  !!
  !!
  subroutine apply_VBoundaryCondO( &
       & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt &
       & )

    ! モジュール引用; Use statements
    !

    use BoundCondSet_mod, only: &
         & SurfTempRelaxedTime, SurfSaltRelaxedTime

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,2) :: &
         & xya_UVelVBCRHS, xya_VVelVBCRHS, xya_PTempEddVBCRHS, xya_SaltVBCRHS

    ! 実行文; Executable statement
    !

    !
    call calc_VBCRHS( &
       & xya_UVelVBCRHS, xya_VVelVBCRHS, xya_PTempEddVBCRHS, xya_SaltVBCRHS, &
       & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt &
       & )

    !
    call solve_VBCEq( xyz_U, &
         & DynBC_Surface, DynBC_Bottom, xya_UVelVBCRHS )

    call solve_VBCEq( xyz_V, &
         & DynBC_Surface, DynBC_Bottom, xya_VVelVBCRHS )
    
    call solve_VBCEq( xyz_PTempEdd, &
         & ThermBC_Surface, ThermBC_Bottom, xya_PTempEddVBCRHS )
    
    call solve_VBCEq( xyz_Salt, &
         & SaltBC_Surface, SaltBC_Bottom, xya_SaltVBCRHS )

  end subroutine apply_VBoundaryCondO

  subroutine calc_VBCRHS( &
       & xya_UVelVBCRHS, xya_VVelVBCRHS, xya_PTempEddVBCRHS, xya_SaltVBCRHS, &
       & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt &
       & )

    ! モジュール引用; Use statements
    !
    use Constants_mod, only: &
         & RefDens, Cp0

    use VariableSet_mod, only: &
         & xyz_VViscCoef, xyz_VDiffCoef
    
    use BoundCondSet_mod, only: &
         & inquire_VBCSpecType, &
         & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress, &
         & ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed, &
         & SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed, &
         & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
         & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax,2) :: &
         & xya_UVelVBCRHS, xya_VVelVBCRHS, xya_PTempEddVBCRHS, xya_SaltVBCRHS

    real(DP), intent(in), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_U, xyz_V, xyz_PTempEdd, xyz_Salt

    ! 局所変数
    ! Local variables
    !
    integer :: k
    real(DP), dimension(0:iMax-1,jMax) :: xy_Coef
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_PTempBasic, xyz_PTempBasicDSig, xyz_Tmp
    
    ! 実行文; Executable statement
    !

    !
    !

    xy_Coef(:,:) = xy_totDepthBasic/(RefDens*xyz_VViscCoef(:,:,0))
    select case(DynBC_Surface)
    case(DynBCTYPE_SpecStress)
       xya_UVelVBCRHS(:,:,1) = xy_Coef*xy_WindStressU;
       xya_VVelVBCRHS(:,:,1) = xy_Coef*xy_WindStressV
    case(DynBCTYPE_NoSlip)
       xya_UVelVBCRHS(:,:,1) = 0d0; xya_VVelVBCRHS(:,:,1) = 0d0
    case(DynBCTYPE_Slip)
       xya_UVelVBCRHS(:,:,1) = 0d0; xya_VVelVBCRHS(:,:,1) = 0d0
    end select

    select case(DynBC_Bottom)
    case(DynBCTYPE_SpecStress)
       call throw_UnImplementVBCError('DynBC_Bottom', DynBCTYPE_SpecStress)
    case(DynBCTYPE_NoSlip)
       xya_UVelVBCRHS(:,:,2) = 0d0; xya_VVelVBCRHS(:,:,2) = 0d0       
    case(DynBCTYPE_Slip)
       xya_UVelVBCRHS(:,:,2) = 0d0; xya_VVelVBCRHS(:,:,2) = 0d0       
    end select


    !
    !
    xy_Coef = xy_totDepthBasic/(RefDens*Cp0*xyz_vDiffCoef(:,:,0))
    forAll(k=0:kMax) &
         & xyz_PTempBasic(:,:,k) = z_PTempBasic(k)
    xyz_PTempBasicDSig = xyz_DSig_xyz(xyz_PTempBasic)
    
    select case(ThermBC_Surface)
    case(ThermBCTYPE_PrescTemp)
       xya_PTempEddVBCRHS(:,:,1) = xy_SeaSurfTemp - z_PTempBasic(0)
    case(ThermBCTYPE_PrescFlux)
       xya_PTempEddVBCRHS(:,:,1) = xy_Coef*( &
            & - xy_SurfHFlxO &
            & ) - xyz_PTempBasicDSig(:,:,0)
    case(ThermBCTYPE_Adiabat)
       xya_PTempEddVBCRHS(:,:,1) = - xyz_PTempBasicDSig(:,:,0)
    case(ThermBCTYPE_TempRelaxed)
       xya_PTempEddVBCRHS(:,:,1) = xy_SeaSurfTemp - z_PTempBasic(0)
    end select

    select case(ThermBC_Bottom)
    case(ThermBCTYPE_PrescTemp)
       call throw_UnImplementVBCError('ThermBC_Bottom', ThermBCTYPE_PrescTemp)
    case(ThermBCTYPE_PrescFlux)
       call throw_UnImplementVBCError('ThermBC_Bottom', ThermBCTYPE_PrescFlux)
    case(ThermBCTYPE_Adiabat)
       xya_PTempEddVBCRHS(:,:,2) = - xyz_PTempBasicDSig(:,:,kMax)       
    end select

    
    !
    xy_Coef = xy_totDepthBasic/xyz_vDiffCoef(:,:,0)
    select case(SaltBC_Surface)
    case(SaltBCTYPE_PrescSalt)
       xya_SaltVBCRHS(:,:,1) = xy_SeaSurfSalt
    case(SaltBCTYPE_PrescFlux)
       xya_SaltVBCRHS(:,:,1) = xy_Coef*( &
            & -xy_SurfFwFlxO*RefSalt_VBC &
            & )
    case(SaltBCTYPE_Adiabat)
       xya_SaltVBCRHS(:,:,1) = 0d0       
    end select

    !
    select case(SaltBC_Bottom)
    case(SaltBCTYPE_PrescSalt)
       call throw_UnImplementVBCError('SaltBC_Bottom', SaltBCTYPE_PrescSalt)
    case(SaltBCTYPE_PrescFlux)
       call throw_UnImplementVBCError('SaltBC_Bottom', SaltBCTYPE_PrescFlux)
    case(SaltBCTYPE_Adiabat)
       xya_SaltVBCRHS(:,:,2) = 0d0
    End select

    
  contains
    subroutine throw_UnImplementVBCError(boundaryLabel, BCTypeID)
      character(*), intent(in) :: boundaryLabel
      integer, intent(in) :: BCTypeID

      call MessageNotify( 'E', module_name, &
           & "'%a=%d' has not been implemeneted yet.", &
           & ca=(/ boundaryLabel /), i=(/ BCTypeID /)  &
           & )
    end subroutine throw_UnImplementVBCError
    
  end subroutine calc_VBCRHS

  subroutine solve_VBCEq(xyz, &
       & SurfBoundaryID, BtmBoundaryID, xya_VBCRHS &
       & )
    real(DP), intent(inout) :: xyz(0:iMax-1,jMax,0:kMax)
    integer, intent(in) :: &
         & SurfBoundaryID, BtmBoundaryID
    real(DP), intent(in) :: xya_VBCRHS(0:iMax-1,jMax,2)

    character :: SurfVBCType, BtmVBCType
    real(DP), dimension(0:kMax,0:kMax) :: &
         & VBCMat, DSigMat, IMat
    integer :: i, j, k, n 
    real(DP) :: IPiv(0:kMax)
    real(DP) :: b(0:kMax,1)
    integer :: info
    
    SurfVBCType = inquire_VBCSpecType(SurfBoundaryID)
    BtmVBCType = inquire_VBCSpecType(BtmBoundaryID)
    
    if (SurfVBCType//BtmVBCType == 'DD') then
       xyz(:,:,0) = xya_VBCRHS(:,:,1)
       xyz(:,:,kMax) = xya_VBCRHS(:,:,2)
       return
    end if

    !
    IMat(:,:) = 0d0
    forAll(k=0:kMax) IMat(k,k) = 1d0
    do k=0, kMax
       DSigMat(:,k) = z_DSig_z(IMat(:,k))
    end do
    
    VBCMat(:,:) = IMat
    if(SurfVBCType == 'N') VBCMat(0,:) = DSigMat(0,:)
    if(BtmVBCType == 'N') VBCMat(kMax,:) = DSigMat(kMax,:)

    n = size(VBCMat, 1)
    call DGETRF(n, n, VBCMat, n, IPiv, info)

    !$omp parallel do private(i, b, info)
    do j=1, jMax
       do i=0, iMax-1
          b(:,1) = xyz(i,j,:);
          b(0,1) = xya_VBCRHS(i,j,1); b(kMax,1) = xya_VBCRHS(i,j,2) 
          
          call DGETRS('N', n, 1, VBCMat, n, IPiv, b, n, info)
          xyz(i,j,0) = b(0,1); xyz(i,j,kMax) = b(kMax,1)
       end do
    end do
    
  end subroutine solve_VBCEq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepare_output()

    character(TOKEN) :: dims_XYT(3)

    dims_XYT = (/ 'lon', 'lat', 't  ' /)
    
    call HistoryAutoAddVariable( varName=VARSET_KEY_SURFHFLXO, &
         & dims=dims_XYT, longname='net heat flux at sea surface', units='W/m2')

    call HistoryAutoAddVariable( varName=VARSET_KEY_SURFFWFLXO, &
         & dims=dims_XYT, longname='freshwater flux at sea surface', units='m/s')
    
  end subroutine prepare_output

  subroutine output_surfaceFlux(CurrentTime)
    real(DP), intent(in) :: CurrentTime

    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFHFLXO, xy_SurfHFlxO)
    call HistoryAutoPut(CurrentTime, VARSET_KEY_SURFFWFLXO, xy_SurfFwFlxO)
    
  end subroutine output_surfaceFlux
  
end module BoundaryCondO_mod

