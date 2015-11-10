!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SeaIceEq_TimeInteg_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use UnitConversion_mod, only: &
       & degC2K, K2degC
  
  use GridSet_mod, only: &
       & iMax, jMax, kMax

  use SeaIceConstants_mod, only: &
       & SBConst, &
       & Mu, SaltSeaIce, &
       & CIce, emissivOcean

  use Constants_mod, only: &
       & CpOcn => Cp0, &
       & DensSeaWater => RefDens

  use SpmlUtil_mod, only: &
       & g_Sig
    
  use GridSet_mod, only: &
       & z_LyrThickSig
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SeaIceEq_TimeInteg_Init, SeaIceEq_TimeInteg_Final
  public :: SeaIceEqSolver_AdvanceTStep
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SeaIceEq_TimeInteg_mod' !< Module Name

  integer, parameter :: DEBUG_j = -1

contains

  !>
  !!
  !!
  subroutine SeaIceEq_TimeInteg_Init()

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         & SeaIceThermDyn_Winton2000_Init
    
    ! 実行文; Executable statements
    !

    call SeaIceThermDyn_Winton2000_Init()
    
  end subroutine SeaIceEq_TimeInteg_Init

  !>
  !!
  !!
  subroutine SeaIceEq_TimeInteg_Final()

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         & SeaIceThermDyn_Winton2000_Final    
    
    ! 実行文; Executable statements
    !

    call SeaIceThermDyn_Winton2000_Final()
    
  end subroutine SeaIceEq_TimeInteg_Final

  !> @brief 
  !!
  !!
  subroutine SeaIceEqSolver_AdvanceTStep(DelTime)

    ! モジュール引用; Use statements
    !

    use BoundaryCondO_mod, only: &
         & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, &
         & xy_SurfHFlxIO, xy_Wsnow, &
         & xy_DSurfHFlxDTsTmp => xy_DSurfHFlxDTs
    
    use VariableSet_mod, only: &
         & z_PTempBasic, xyz_PTempEddN, &
         & xyz_UON => xyz_UN, xyz_VON => xyz_VN, &
         & xyz_SaltON => xyz_SaltN, &
         & xy_totDepthBasic
    
    use VarSetSeaice_mod, only: &
         & xy_SIceConN, xy_SnowThickN, xy_IceThickN, xy_SIceSurfTempN, xyz_SIceTempN, &
         & xy_SIceConA, xy_SnowThickA, xy_IceThickA, xy_SIceSurfTempA, xyz_SIceTempA, &
         & xy_Wice
         

    use SeaIceBoundaryCond_mod, only: &
         & calc_SurfaceHeatFluxSIce, calc_BottomHeatFluxSIce

    use GridSet_mod, only: &
         & z_LyrThickSig
    use SpmlUtil_mod, only: &
         & g_Sig
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: DelTime
    
    ! 局所変数
    ! Local variables
    !
    integer :: i, j, k
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfHFlxAIN, xy_SurfHFlxAON, xy_PenSWRFlxSIN, xy_BtmHFlxION, &
         & xy_DSurfHFlxDTsAIN
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTempO
    real(DP), dimension(0:iMax-1,jMax) :: xy_FreezePot, xy_FreezeTempO, xy_SurfMeltEn, xy_BtmMeltEn
    real(DP) :: SurfTempO, FreezePotTmp, FreezePotRes, dPTemp
    
    ! 実行文; Executable statement
    !    
    do k=0, kMax
       xyz_PTempO(:,:,k) = xyz_PTempEddN(:,:,k) + z_PTempBasic(k)
    end do
   
    
    call calculate_FreezeOcnInfo( &
         & xy_FreezePot, xy_FreezeTempO,    &  ! (out)
         & xyz_PTempO, xyz_SaltON, DelTime  &  ! (in)
         & )
    
    call calc_SurfaceHeatFluxSIce( xy_SurfHFlxAIN, xy_SurfHFlxAON, xy_PenSWRFlxSIN,          & ! (out) 
         & xy_DSurfHFlxDTsAIN,                                                               & ! (out)
         & xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, xy_DSurfHFlxDTsTmp,     & ! (in) 
         & xy_SIceConN, xy_SnowThickN, degC2K(xy_SIceSurfTempN), xyz_PTempO(:,:,0)           & ! (in)
         & )

    call calc_BottomHeatFluxSIce( xy_BtmHFlxION,                         & ! (out)
         & xy_FreezePot, xy_FreezeTempO,                                 & ! (in)x
         & xy_SIceConN, xyz_UON(:,:,0), xyz_VON(:,:,0), xyz_PTempO(:,:,0), xyz_SaltON(:,:,0), & ! (in)
         & DelTime )

    !$omp parallel do private(i, SurfTempO, FreezePotTmp, FreezePotRes, k, dPTemp)
    do j=1, jMax
       do i=0, iMax-1

          !
          xy_BtmMeltEn(i,j) = xy_BtmHFlxION(i,j)
          if(j==DEBUG_j) then
             write(*,*) "---------- Start sea ice process ----------------"             
             write(*,*) "* TempO=", K2degC(xyz_PTempO(i,j,0:1)), ", SaltO=", xyz_SaltON(i,j,0)
             write(*,*) "* SurfTempN=", xy_SIceSurfTempN(i,j), ", IceThickN=", xy_IceThickN(i,j), &
                  &     ", BtmMeltEn=", xy_BtmMeltEn(i,j)
             write(*,*) "* FlxInfo(Down +):", "LW=", xy_LWDWRFlx(i,j), "SW=", xy_SWDWRFlx(i,j), "Latent=", xy_LatentDWHFlx(i,j), &
                  &                   "Sens=", xy_SensDWHFlx(i,j), "LWup=", emissivOcean*SBConst*degC2K(xy_SIceSurfTempN(i,j))**4
          end if
          
          !
          if(xy_IceThickN(i,j) > 0d0) then
             call advance_ThermDyncProc_SIceTemp_zColum( &
                  & xy_SIceSurfTempA(i,j), xyz_SIceTempA(i,j,1), xyz_SIceTempA(i,j,2),        &  ! (out)  
                  & xy_SurfMeltEn(i,j),                                                       &  ! (out)
                  & xy_BtmMeltEn(i,j),                                                        &  ! (inout)
                  & xy_SnowThickN(i,j), xy_IceThickN(i,j),                                    &  ! (in) 
                  & xy_SIceSurfTempN(i,j), xyz_SIceTempN(i,j,1), xyz_SIceTempN(i,j,2),        &  ! (in)
                  & xy_SurfHFlxAIN(i,j),   xy_DSurfHFlxDTsAIN(i,j), xy_PenSWRFlxSIN(i,j),     &  ! (in)
                  & DelTime                                                                   &  ! (in)
                  & )
          else
             xy_SIceSurfTempA(i,j) = 0d0; xyz_SIceTempA(i,j,:) = 0d0
          end if

          !
          !
          if(j==DEBUG_j) then
             write(*,*) " - SurfTempA=", xy_SIceSurfTempA(i,j), "SIceTempA=", xyz_SIceTempA(i,j,:)
             write(*,*) " - SurfMeltEn=", xy_SurfMeltEn(i,j), "BtmMeltEn=", xy_BtmMeltEn(i,j)
          end if

          SurfTempO = xyz_PTempO(i,j,0)
          call advance_SeaIceThermDynProc_SIceThick_zColum( &
               & xy_SnowThickA(i,j), xy_IceThickA(i,j), xy_Wice(i,j),  & !(out)
               & xy_SIceSurfTempA(i,j), xyz_SIceTempA(i,j,1), xyz_SIceTempA(i,j,2), SurfTempO,   &  ! (inout)
               & xy_SnowThickN(i,j), xy_IceThickN(i,j),                                &  ! (in) 
               & xy_SIceSurfTempN(i,j), xyz_SIceTempN(i,j,1), xyz_SIceTempN(i,j,2),    &  ! (in)
               & xy_SurfMeltEn(i,j), xy_BtmMeltEn(i,j), xy_Wsnow(i,j), DelTime,        & ! (in)
               & i, j)
          if(SurfTempO /= xyz_PTempO(i,j,0)) then
             !
             FreezePotRes = xy_FreezePot(i,j)*DelTime/(CpOcn*DensSeaWater)

             do k=0, 0
                ! J/m2/s * s * /
                dPTemp = SurfTempO - xyz_PTempO(i,j,k)
                if(dPTemp > 0d0  .and. FreezePotRes > 0d0) then
                   xyz_PTempO(i,j,k) = SurfTempO
                   FreezePotRes = FreezePotRes - dPTemp*z_LyrThickSig(k)*xy_totDepthBasic(i,j)
                end if
                xyz_PTempEddN(i,j,k) = xyz_PTempO(i,j,k) - z_PTempBasic(k)
             end do
          end if
          
          if(xy_IceThickA(i,j) > 0d0) then
             xy_SIceConA(i,j) = 1d0
          else
             xy_SIceConA(i,j) = 0d0
          end if
          if(j==DEBUG_j) then
             write(*,*) " -- TempO=", -273.15d0+(xyz_PTempO(i,j,0:1)), " SaltO=", xyz_SaltON(i,j,0)                          
             write(*,*) " -- IceThickA=", xy_IceThickA(i,j), ", SnowThickA=", xy_SnowThickA(i,j)
          end if
       end do
    end do

    
    
    !
    !write(*,*) "Recalculate BtmHFlxSIce.."
    call calculate_FreezeOcnInfo( &
         & xy_FreezePot, xy_FreezeTempO,    &  ! (out)
         & xyz_PTempO, xyz_SaltON, DelTime  &  ! (in)
         & )
    
    call calc_BottomHeatFluxSIce( xy_SurfHFlxIO,                         & ! (out)
         & xy_FreezePot, xy_FreezeTempO, & !(in)
         & xy_SIceConA, xyz_UON(:,:,0), xyz_VON(:,:,0), xyz_PTempO(:,:,0), xyz_SaltON(:,:,0), & ! (in)
         & DelTime )
    !write(*,*) "***************************************"
    
  end subroutine SeaIceEqSolver_AdvanceTStep

  !> @brief 
  !!
  !!
  subroutine calculate_FreezeOcnInfo( &
       & xy_FreezePot, xy_FreezeTempO, &
       & xyz_PTempO, xyz_SaltO, DelTime )

    ! モジュール引用; Use statements
    !
    
    use VariableSet_mod, only: &
         & xy_totDepthBasic

    use SeaIceConstants_mod, only: &
         & FreezeTempSW, DensSeaWater
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(0:iMax-1,jMax), intent(out) :: &
         & xy_FreezePot, xy_FreezeTempO
    real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: &
         & xyz_PTempO, xyz_SaltO
    real(DP), intent(in) :: DelTime
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: dPTemp, TempOIntSurf
    integer :: i, j, k
    logical :: FreezeCondFlag
    
    ! 実行文; Executable statement
    !

    xy_FreezeTempO(:,:) = degC2K(FreezeTempSW)!- Mu*xy_SurfSaltO

    do j=1, jMax
       do i=0, iMax-1
          dPTemp = xy_FreezeTempO(i,j) - xyz_PTempO(i,j,0)
          TempOIntSurf = dPTemp*z_LyrThickSig(0)
          xy_FreezePot(i,j) = CpOcn*DensSeaWater/DelTime*TempOIntSurf*xy_totDepthBasic(i,j)
       end do
    end do
    
  end subroutine calculate_FreezeOcnInfo

  !> @brief 
  !!
  !!
  subroutine advance_ThermDyncProc_SIceTemp_zColum( &
       & TsA, T1A, T2A, SurfMeltEn,                        & ! (out)
       & BtmMeltEn,                                        & ! (inout)
       & hsN, hiN, TsN, T1N, T2N,                          & ! (in)
       & SurfHFlxAIN, DSurfHFlxDTs, PenSWRFlxSIN,          & ! (in)
       & DelTime                & ! (in)
       & )

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         & update_effConductiveCoupling, calc_layersTemprature

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: TsA, T1A, T2A
    real(DP), intent(out) :: SurfMeltEn, BtmMeltEn
    real(DP), intent(in) :: hsN, hiN, TsN, T1N, T2N
    real(DP), intent(in) :: SurfHFlxAIN, DSurfHFlxDTs, PenSWRFlxSIN
    real(DP), intent(in) :: DelTime
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: K_12, K_23
    real(DP) :: A, B
    
    
    ! 実行文; Executable statement
    !

    call update_effConductiveCoupling(K_12, K_23, & ! (out)
         & hsN, hiN )                               ! (in)

    B = DSurfHFlxDTs
    A = SurfHFlxAIN - TsN*B

    !
    call calc_layersTemprature( &
         & TsA, T1A, T2A, SurfMeltEn,               &  ! (out)
         & BtmMeltEn,                               &  ! (inout)
         & T1N, T2N, hiN, DelTime, A, B,            &  ! (in)
         & -PenSWRFlxSIN, K_12, K_23, hsN > 0d0      &  ! (in)
         & )
    
  end subroutine advance_ThermDyncProc_SIceTemp_zColum

  !>
  !! @param FWFlx freshwater flux generated by 
  !! @param SurfHFxAIN surface heat flux at time level n(**positive upward**). [W/m2]
  !! @param 
  subroutine advance_SeaIceThermDynProc_SIceThick_zColum( &
       & hsA, hiA, Wice,                        & !(out)
       & TsA, T1A, T2A, TsO,                    & !(inout)
       & hsN, hiN, TsN, T1N, T2N,               & !(in)
       & SurfMeltEn, BtmMeltEn, SnowFall,             & !(in)
       & DelTime, i, j                                & !(in)
       & )

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         & calc_SnowIceLyrMassChange, adjust_IceLyrInternal
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: hsA, hiA
    real(DP), intent(out) :: Wice
    real(DP), intent(inout) ::  TsA, T1A, T2A, TsO
    real(DP), intent(in)  :: hsN, hiN, TsN, T1N, T2N
    real(DP), intent(in) :: SurfMeltEn, BtmMeltEn
    real(DP), intent(in) :: SnowFall
    real(DP), intent(in) :: DelTime
    integer :: i, j
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: excessMeltEnergy
    real(DP) :: dhs, dh1, dh2
    
    ! 実行文; Executable statements
    !

    ! 
    !
    

    !
    call calc_SnowIceLyrMassChange( &
         & dhs, dh1, dh2, excessMeltEnergy,    & ! (out)
         & T1A, T2A, TsO,                      & ! (inout)
         & hsN, hiN, TsA, DelTime,             & ! (in)
         & SurfMeltEn, BtmMeltEn               & ! (in)
         & )

    Wice = (dh1 + dh2 + min(0d0, dhs))/DelTime
    
    if(j==DEBUG_j) then
       write(*,*) "T1, T2=", T1A, T2A
       write(*,*) "dhs, dh1, dh2=", dhs, dh1, dh2
       write(*,*) "SnowFall=", SnowFall*DelTime
    end if

    !
    hsA = hsN + dhs
    if(0.5d0*hiN + dh1> 0d0) hsA = hsA + SnowFall*DelTime
    call adjust_IceLyrInternal( &
         & hsA, hiA, T1A, T2A,                & !(inout)
         & 0.5d0*hiN+dh1, 0.5d0*hiN+dh2       & !(in)
         & )
    
    if(j==DEBUG_j) then
       write(*,*) "Snow2Ice=", hsA - hsN - dhs - SnowFall*DelTime
       write(*,*) "reshape_IceLyr => T1, T2=", T1A, T2A
    end if
  end subroutine advance_SeaIceThermDynProc_SIceThick_zColum
  
end module SeaIceEq_TimeInteg_mod

