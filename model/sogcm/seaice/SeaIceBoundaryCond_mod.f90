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
       & BaseMeltHeatTransCoef
  
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
  
contains

  !>
  !!
  !!
  subroutine SeaIceBoundaryCond_Init()

    ! 実行文; Executable statements
    !

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
    & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx,       & ! (in) 
    & xy_SIceCon, xy_SnowThick, xy_SurfTemp, xy_SurfTempO                         & ! (in)
    & )

    ! モジュール引用; use statements
    !
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHFlxAI, xy_SurfHFlxAO, xy_PenSWRFlxSI
  
    real(DP), intent(in), dimension(0:iMax-1,jMax) :: &
         & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, &
         & xy_SIceCon, xy_SnowThick, xy_SurfTemp, xy_SurfTempO
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax) :: xy_Emissivity, xy_Albedo
    
    ! 実行文; Executable statement
    !
    
    where(xy_SIceCon < 1d0)
       xy_Albedo = AlbedoOcean; xy_Emissivity = EmissivOcean
    end where
    call calcSurfHFlux( xy_SurfHFlxAO, xy_Albedo, xy_Emissivity, xy_SurfTempO)

    where( xy_SIceCon <= 0d0)
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
    call calcSurfHFlux(xy_SurfHFlxAI, xy_Albedo, xy_Emissivity, xy_SurfTemp)

!    write(*,*) xy_SurfHFlxAO
!!$    write(*,*) "* SIceSurfFlx:", "albedo=", xy_Albedo(0,DEBUG_j), "Emissivty=", xy_Emissivity(0,DEBUG_j), &
!!$         & "LI=", - xy_LatentDWHFlx(0,DEBUG_j), "SI=", - xy_SensDWHFlx(0,DEBUG_j), &
!!$         & "SW=", - (1d0 - xy_Albedo(0,DEBUG_j))*xy_SWDWRFlx(0,DEBUG_j), &
!!$         & "LW=", -xy_Emissivity(0,DEBUG_j)*xy_LWDWRFlx(0,DEBUG_j), xy_Emissivity(0,DEBUG_j)*SBConst*xy_SurfTemp(0,DEBUG_j)**4
    
  contains

    !>
    !! @param xy_SurfHFlx Surface heat flux (positive upward) [W/m2]
    !!
    subroutine calcSurfHFlux(xy_SurfHFlux, xy_Albedo, xy_Emissivity, xy_SurfTemp)
      real(DP), dimension(0:iMax-1,jMax), intent(out) :: xy_SurfHFlux
      real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_Albedo, xy_Emissivity, xy_SurfTemp

      !$omp parallel workshare
      xy_SurfHFlux(:,:) = &
           & - xy_LatentDWHFlx - xy_SensDWHFlx       &
           & - xy_Emissivity*xy_LWDWRFlx             &
           & - (1d0 - xy_Albedo)*xy_SWDWRFlx         &
           & + xy_Emissivity*(SBConst*xy_SurfTemp**4)
      !$omp end parallel workshare
    end subroutine calcSurfHFlux
    
  end subroutine calc_SurfaceHeatFluxSIce


  !> @brief 
  !!
  !! @param xy_BtmHFlx bottom heat flux from ocean (oceanic heat flux)  (positive upward) [W/m2]  
  !! @param xy_SurfTempO Temperature of ocea surface [Units: K]  
  !!
  subroutine calc_BottomHeatFluxSIce( xy_BtmHFlxIO, &
       & xy_FreezePot, xy_FreezeTempO, &
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

    !$omp parallel workshare
    
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

    !$omp end parallel workshare
    
!!$    xyz_Tmp(0,DEBUG_j,:) = z_DSig_z(z_PTempBasic(:) + xyz_PTempEddN(0,DEBUG_j,:))/xy_totDepthBasic(0,DEBUG_j)    
!!$    write(*,*) "* BtmHFlxSIce=", xy_BtmHFlxIO(:,DEBUG_j), ", BtmHFlxO", -vDiffCoef*xyz_Tmp(0,DEBUG_j,0)*refDens*Cp0, &
!!$         & "FreezePot=", xy_FreezePot(:,DEBUG_j), "SurfTempO", K2degC(xy_SurfTempO(:,DEBUG_j)), &
!!$         & " BulkTransCoef=", BaseMeltHeatTransCoef*xy_FricVel(:,DEBUG_j)
!!$    write(*,*) "----", K2degC(xyz_PTempEddN(0,DEBUG_j,0:4) + z_PTempBasic(0:4))
  end subroutine calc_BottomHeatFluxSIce

end module SeaIceBoundaryCond_mod

