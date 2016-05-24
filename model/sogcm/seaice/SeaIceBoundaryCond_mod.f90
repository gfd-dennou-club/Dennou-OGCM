!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SeaIceBoundaryCond_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN

  use GridSet_mod, only: &
       & iMax, jMax, kMax, &
       & z_LyrThickSig

  !* Dennou-OGCM

  use UnitConversion_mod, only: &
       & degC2K, K2degC

!  use SpmlUtil_mod
  
  use VariableSet_mod, only: &
       & xy_totDepthBasic

  
  use SeaIceConstants_mod, only: &
       & SBConst, &
       & FreezeTempWater, &
       & Mu, &
       & DensSeaWater, &
       & AlbedoOcean, AlbedoSnow, AlbedoMeltSnow, AlbedoIce, &
       & EmissivOcean, EmissivSnow, EmissivIce, &
       & I0, &
       & BaseMeltHeatTransCoef, &
       & IceMaskMin

  use BoundCondSet_mod, only: &
       & inquire_VBCSpecType, &
       & ThermBCTYPE_PrescFixedFlux, ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, &
       & ThermBCTYPE_PresFlux_Han1984Method, &
       & ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed,                         & 
       & SaltBCTYPE_PrescFixedFlux, SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, &
       & SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed, &
       & SaltBCTYPE_PresFlux_Han1984Method, &
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom 
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SeaIceBoundaryCond_Init, SeaIceBoundaryCond_Final
  public :: calc_SurfaceHeatFluxSIce, calc_BottomHeatFluxSIce

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SeaIceBoundaryCond_mod' !< Module Name

  integer, parameter :: DEBUG_j = -1

  logical :: init_SurfTempAIOSetFlag
  logical :: init_SurfTempOSetFlag
  real(DP), allocatable :: xy_SurfTempAIOIni(:,:)
  real(DP), allocatable :: xy_SurfTempOIni(:,:)

contains

  !>
  !!
  !!
  subroutine SeaIceBoundaryCond_Init()

    ! 実行文; Executable statements
    !

    init_SurfTempAIOSetFlag = .false.
    init_SurfTempOSetFlag = .false.

  end subroutine SeaIceBoundaryCond_Init

  !>
  !!
  !!
  subroutine SeaIceBoundaryCond_Final()

    ! 実行文; Executable statements
    !

  end subroutine SeaIceBoundaryCond_Final

  !> @brief 
  !!
  !! @param xy_SurfTemp  Temperature of snow or sea-ice(if there is no snow) surface [Units: K].
  !! @param xy_SurfTempO Temperature of ocean surface [Units: K]
  !!
  subroutine calc_SurfaceHeatFluxSIce( xy_SurfHFlxAI, xy_SurfHFlxAO, xy_PenSWRFlxSI,       & ! (out)
    & xy_DSurfHFlxDTsAI,                                                                   & ! (out)
    & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, xy_DSurfHFlxDTsTmp,        & ! (in) 
    & xy_SIceCon, xy_SnowThick, xy_SurfTemp, xy_SurfTempO                                  & ! (in)
    & )

    ! モジュール引用; use statements
    !
    use BoundaryCondO_mod, only: xy_SWUWRFLX, xy_LWUWRFlx

    use Constants_mod, only: UNDEFVAL

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHFlxAI, xy_SurfHFlxAO, xy_PenSWRFlxSI, &
         & xy_DSurfHFlxDTsAI
  
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, xy_DSurfHFlxDTsTmp, &
         & xy_SIceCon, xy_SnowThick, xy_SurfTemp, xy_SurfTempO
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax) :: xy_Emissivity, xy_Albedo
    integer :: i, j

    ! 実行文; Executable statement
    !
    

    xy_Albedo = AlbedoOcean; xy_Emissivity = EmissivOcean
    call calcSurfHFlux( xy_SurfHFlxAO, xy_Albedo, xy_Emissivity, xy_SurfTempO, xy_DSurfHFlxDTsTmp, &
         & xy_DSurfHFlxDTsAI )
    select case(ThermBC_Surface)
    case(ThermBCTYPE_PresFlux_Han1984Method)
       if (.not. init_SurfTempOSetFlag) then
          allocate(xy_SurfTempOIni(0:iMax-1,jMax))
          xy_SurfTempOIni = xy_SurfTempO
          init_SurfTempOSetFlag = .true.
       end if
       xy_SurfHFlxAO(:,:) = xy_SurfHFlxAO(:,:)                   &
            + ( 4d0*xy_Emissivity*(SBConst*xy_SurfTempOIni**3) + xy_DSurfHFlxDTsTmp ) &
               *( xy_SurfTempO - xy_SurfTempOIni )
    end select

    where( xy_SIceCon <= IceMaskMin)
    elsewhere( xy_SnowThick <= 0d0)
       xy_Albedo = AlbedoIce; xy_Emissivity = EmissivIce
       xy_PenSWRFlxSI = - I0*(1d0 - xy_Albedo)*xy_SWDWRFlx       
    elsewhere( xy_SurfTemp >= degC2K(FreezeTempWater) )
       xy_Albedo = AlbedoMeltSnow; xy_Emissivity = EmissivSnow
       xy_PenSWRFlxSI = 0d0       
    elsewhere
       xy_Albedo = AlbedoSnow; xy_Emissivity = EmissivSnow
       xy_PenSWRFlxSI = 0d0
    end where
    
!!$    !$omp parallel do private(i)
!!$    do j=1, jMax
!!$       do i=0, iMax-1
!!$          if (xy_SIceCon(i,j) <= IceMaskMin) then
!!$             xy___
!!$          else
!!$          end if
!!$       end do
!!$    end do
    call calcSurfHFlux( xy_SurfHFlxAI, &
         & xy_Albedo, xy_Emissivity, xy_SurfTemp, xy_DSurfHFlxDTsTmp, &
         & xy_DSurfHFlxDTsAI)

!!$    write(*,*) "Albedo=", xy_Albedo(0,:)
!!$    write(*,*) "SWDWRFlx=", xy_SWDWRFlx(0,:)
!!$    write(*,*) "Pen=", xy_PenSWRFlxSI(0,:)
!!$    write(*,*) "SIceCon=", xy_SIceCon(0,:)
!!$    write(*,*) "SurfTemp=", xy_SurfTemp(0,:)

!!$    select case(ThermBC_Surface)
!!$    case(ThermBCTYPE_PresFlux_Han1984Method)
!!$       if (.not. init_SurfTempAIOSetFlag) then
!!$          allocate(xy_SurfTempAIOIni(0:iMax-1,jMax))
!!$          xy_SurfTempAIOIni = xy_SurfTemp
!!$          init_SurfTempAIOSetFlag = .true.
!!$       end if
!!$       xy_SurfHFlxAI(:,:) = xy_SurfHFlxAI(:,:)                   &
!!$            + ( 4d0*xy_Emissivity*(SBConst*xy_SurfTempAIOIni**3) + xy_DSurfHFlxDTsTmp ) &
!!$               *( xy_SurfTemp - xy_SurfTempAIOIni )
!!$    end select

    where( xy_SIceCon <= IceMaskMin ) 
       xy_SurfHFlxAI = 0d0
       xy_PenSWRFlxSI = 0d0
    end where

!    write(*,*) xy_SurfHFlxAO
!!$    write(*,*) "* SIceSurfFlx:", "albedo=", xy_Albedo(0,DEBUG_j), "Emissivty=", xy_Emissivity(0,DEBUG_j), &
!!$         & "LI=", - xy_LatentDWHFlx(0,DEBUG_j), "SI=", - xy_SensDWHFlx(0,DEBUG_j), &
!!$         & "SW=", - (1d0 - xy_Albedo(0,DEBUG_j))*xy_SWDWRFlx(0,DEBUG_j), &
!!$         & "LW=", -xy_Emissivity(0,DEBUG_j)*xy_LWDWRFlx(0,DEBUG_j), xy_Emissivity(0,DEBUG_j)*SBConst*xy_SurfTemp(0,DEBUG_j)**4
    
  contains

    !>
    !! @param xy_SurfHFlx Surface heat flux (positive upward) [W/m2]
    !!
    subroutine calcSurfHFlux( &
         & xy_SurfHFlux,                            & ! (out)
         & xy_Albedo, xy_Emissivity, xy_SurfTemp,   & ! (in)
         & xy_DSurfHFlxDTsTmp,                      & ! (in)
         & xy_DSurfHFlxDTs                          & ! (out)
         & )

      real(DP), dimension(0:iMax-1,jMax), intent(out) :: xy_SurfHFlux
      real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_Albedo, xy_Emissivity, xy_SurfTemp
      real(DP), dimension(0:iMax-1,jMax), intent(in), optional :: xy_DSurfHFlxDTsTmp
      real(DP), dimension(0:iMax-1,jMax), intent(out), optional :: xy_DSurfHFlxDTs

      !$omp parallel
      !$omp workshare
      xy_SurfHFlux(:,:) = &
           & - xy_LatentDWHFlx - xy_SensDWHFlx        &
           & - xy_Emissivity*xy_LWDWRFlx              &
           & - (1d0 - xy_Albedo)*xy_SWDWRFlx          &
           & + xy_Emissivity*(SBConst*xy_SurfTemp**4)
      !$omp end workshare
      !$omp end parallel

!!$      !$omp parallel
!!$      !$omp workshare
!!$      xy_SurfHFlux(:,:) = &
!!$           & - xy_LatentDWHFlx - xy_SensDWHFlx        &
!!$           & - (xy_LWDWRFlx - xy_LWUWRFlx)            &
!!$           & - (1d0 - xy_Albedo)*xy_SWDWRFlx          
!!$      !$omp end workshare
!!$      !$omp end parallel


!      write(*,*) xy_LWUWRFlx(0,:)
!      write(*,*) SBConst*xy_SurfTemp(0,:)**4
!      stop

      if(present(xy_DSurfHFlxDTs)) then
         !$omp parallel
         !$omp workshare
         xy_DSurfHFlxDTs(:,:) = &
              &   4d0*xy_Emissivity*(SBConst*xy_SurfTemp**3) &
              & + xy_DSurfHFlxDTsTmp
         !$omp end workshare
         !$omp end parallel
      end if
    end subroutine calcSurfHFlux
    
  end subroutine calc_SurfaceHeatFluxSIce


  !> @brief 
  !!
  !! @param xy_BtmHFlx bottom heat flux from ocean (oceanic heat flux)  (positive upward) [W/m2]  
  !! @param xy_SurfTempO Temperature of ocea surface [Units: K]  
  !!
  subroutine calc_BottomHeatFluxSIce( xy_BtmHFlxIO,                    &
       & xy_FreezePot, xy_FreezeTempO,                                 &
       & xy_SIceCon, xy_SurfUO, xy_SurfVO, xy_SurfTempO, xy_SurfSaltO, &
       & DelTime )

    ! モジュール引用; Use statements
    !
    use Constants_mod, only: &
         & CpOcn => Cp0

    use Constants_mod, only: &
         & RefDens, Cp0, vDiffCoef

    use SeaIceConstants_mod, only: &
         & FreezeTempSW

    use VariableSet_mod, only: &
         & xyz_PTempEddN, z_PTempBasic, &
         & xy_totDepthBasic
    use SpmlUtil_mod
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax), intent(out) :: xy_BtmHFlxIO
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: &
         & xy_FreezePot, xy_FreezeTempO, xy_SIceCon, &
         & xy_SurfUO, xy_SurfVO, xy_SurfTempO, xy_SurfSaltO
    
    real(DP), intent(in) :: DelTime
    
    ! 局所変数
    ! Local variables
    !

    !
    real(DP) :: xy_FricVel(0:iMax-1,jMax)
    
    integer :: i, j

    real(DP) :: xyz_Tmp(0:iMax-1,jMax,0:kMax)
    integer :: k
    
    ! 実行文; Executable statement
    !

    !$omp parallel
    !$omp workshare
    
    xy_BtmHFlxIO(:,:) = 0d0
    where(xy_FreezePot > 0d0)

       !* Freezing condition
       !
       
       xy_BtmHFlxIO = - xy_FreezePot
       
    elsewhere(xy_SIceCon > 0d-3 .and. xy_FreezePot <= 0d0)

       !* Melting condition
       !
       xy_FricVel = max(sqrt(xy_SurfUO**2 + xy_SurfVO**2), 5d-3)
       xy_BtmHFlxIO = - BaseMeltHeatTransCoef*xy_FricVel* &
            & xy_FreezePot*DelTime/(z_LyrThickSig(0)*xy_totDepthBasic)

       !
       xy_BtmHFlxIO = min(xy_BtmHFlxIO, -xy_FreezePot)
    end where

    !$omp end workshare
    !$omp end parallel
    
!!$    xyz_Tmp(0,DEBUG_j,:) = z_DSig_z(z_PTempBasic(:) + xyz_PTempEddN(0,DEBUG_j,:))/xy_totDepthBasic(0,DEBUG_j)    
!!$    write(*,*) "* BtmHFlxSIce=", xy_BtmHFlxIO(:,DEBUG_j), ", BtmHFlxO", -vDiffCoef*xyz_Tmp(0,DEBUG_j,0)*refDens*Cp0, &
!!$         & "FreezePot=", xy_FreezePot(:,DEBUG_j), "SurfTempO", K2degC(xy_SurfTempO(:,DEBUG_j)), &
!!$         & " BulkTransCoef=", BaseMeltHeatTransCoef*xy_FricVel(:,DEBUG_j)
!!$    write(*,*) "----", K2degC(xyz_PTempEddN(0,DEBUG_j,0:4) + z_PTempBasic(0:4))
  end subroutine calc_BottomHeatFluxSIce

end module SeaIceBoundaryCond_mod

