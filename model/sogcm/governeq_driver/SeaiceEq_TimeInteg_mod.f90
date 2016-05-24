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
       & Mu, SaltSeaIce, FreezeTempSW, &
       & LFreeze, DensSnow, DensIce,   &
       & CIce, emissivOcean,           &
       & IceThickMin, IceThickMax, &
       & IceMaskMin, &
       & SIceHDiffCoef

  use Constants_mod, only: &
       & UNDEFVAL, PI, RPlanet, &
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

  integer, parameter :: DEBUG_j = -1 !64
  
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
         xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, &
         xy_SurfHFlxIO, xy_SurfFwFlxIO, xy_SurfHFlxAI, xy_SurfHFlxAO, &
         xy_Wrain, xy_Wsnow, xy_Wevap, &
         xy_DSurfHFlxDTsTmp => xy_DSurfHFlxDTs
    
    use VariableSet_mod, only: &
         z_PTempBasic, xyz_PTempEddN, &
         xyz_UON => xyz_UN, xyz_VON => xyz_VN, &
         xyz_SaltON => xyz_SaltN, &
         xy_totDepthBasic
    
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

    real(DP) :: xy_SurfHFlxAIN(0:iMax-1,jMax)
    real(DP) :: xy_DSurfHFlxDTsAIN(0:iMax-1,jMax)
    real(DP) :: xy_PenSWRFlxSIN(0:iMax-1,jMax)
    real(DP) :: xy_SurfHFlxAON(0:iMax-1,jMax)

    real(DP) :: xy_BtmHFlxAO(0:iMax-1,jMax)    
    real(DP) :: xy_BtmHFlxIO(0:iMax-1,jMax)
    
    real(DP) :: xy_FreezePot(0:iMax-1,jMax)
    real(DP) :: xy_FreezeTempO(0:iMax-1,jMax)
    real(DP) :: xy_SurfMeltEn(0:iMax-1,jMax)
    real(DP) :: xy_BtmMeltEn(0:iMax-1,jMax)
    real(DP) :: xy_ExcessMeltEn(0:iMax-1,jMax)
    
    real(DP) :: xyz_PTempO(0:iMax-1,jMax,0:kMax)

    real(DP) :: T1_, T2_, Ts_, hs_
    real(DP) :: SurfTempO, dPTemp
    real(DP) :: FreezeEnOpnOcn
    
    ! 実行文; Executable statement
    !    

    !* Initialization
    !
    
    do k=0, kMax
       xyz_PTempO(:,:,k) = xyz_PTempEddN(:,:,k) + z_PTempBasic(k)
    end do

    xy_ExcessMeltEn = 0d0
    xy_Wice         = 0d0
    
    call calculate_FreezeOcnInfo( &
         xy_FreezePot, xy_FreezeTempO,      &  ! (out)
         xyz_PTempO, xyz_SaltON, DelTime  )    ! (in)
         
    
    !* Calculate heat flux at the top and bottom surface
    !
    
    call calc_SurfaceHeatFluxSIce( xy_SurfHFlxAIN, xy_SurfHFlxAON, xy_PenSWRFlxSIN,       & ! (out) 
         xy_DSurfHFlxDTsAIN,                                                              & ! (out)
         xy_SWDWRFlx, xy_LWDWRFlx, xy_LatentDWHFlx, xy_SensDWHFlx, xy_DSurfHFlxDTsTmp,    & ! (in) 
         xy_SIceConN, xy_SnowThickN, degC2K(xy_SIceSurfTempN), xyz_PTempO(:,:,0)        )   ! (in)
    xy_SurfHFlxAI = xy_SurfHFlxAIN
    xy_SurfHFlxAO = xy_SurfHFlxAON

    call calc_BottomHeatFluxSIce( xy_BtmHFlxIO,                                              & ! (out)
         xy_FreezePot, xy_FreezeTempO,                                                       & ! (in)
         xy_SIceConN, xyz_UON(:,:,0), xyz_VON(:,:,0), xyz_PTempO(:,:,0), xyz_SaltON(:,:,0),  & ! (in)
         DelTime )                                                                             ! (in)

    !* Thermodynamics processes
    !

    do j=1, jMax
       do i=0, iMax-1
          if (j==DEBUG_j) then
             write(*,*) "Pen=", xy_PenSWRFlxSIN(i,j), "SnoW=", xy_SnowThickN(i,j), &
                "DSurfHFlxDTs=", xy_DSurfHFlxDTsAIN(i,j)
          end if
       end do
    end do

    !$omp parallel do private(i)
    do j=1, jMax
       do i=0, iMax-1

          xy_BtmMeltEn(i,j) = xy_BtmHFlxIO(i,j)
          
          if(j==DEBUG_j) then
             write(*,*) "---------- Start sea ice process ----------------"             
             write(*,*) "* TempO=", K2degC(xyz_PTempO(i,j,0:1)), ", SaltO=", xyz_SaltON(i,j,0)
             write(*,*) "* SurfTempN=", xy_SIceSurfTempN(i,j), ", IceThickN=", xy_IceThickN(i,j), &
                  &     ", BtmMeltEn=", xy_BtmMeltEn(i,j)
             write(*,*) "* FlxInfo :", "AI=",  xy_SurfHFlxAIN(i,j), " AO=", xy_SurfHFlxAON(i,j), " IO=", xy_BtmHFlxIO(i,j)
             write(*,*) "* SurfFlxBudget :", "LW=", -xy_LWDWRFlx(i,j), "SW=", -xy_SWDWRFlx(i,j), "Latent=", -xy_LatentDWHFlx(i,j), &
                  "Sens=", -xy_SensDWHFlx(i,j), "LWup=", -emissivOcean*SBConst*degC2K(xy_SIceSurfTempN(i,j))**4!, &
!                  "DSurfHFlxDTs=", xy_DSurfHFlxDTsAIN(i,j)
             write(*,*) "* IceFrac=", xy_SIceConN(i,j), "SIceTempN=", xyz_SIceTempN(i,j,:)
          end if
          
          !
          !if (xy_IceThickN(i,j)/max(xy_SIceConN(i,j), IceMaskMin) >= IceThickMin) then
          if ( xy_SIceConN(i,j) >= IceMaskMin ) then
             xy_SIceSurfTempA(i,j) = xy_SIceSurfTempN(i,j)
             xyz_SIceTempA(i,j,:) = xyz_SIceTempN(i,j,:)
             if ( xy_SIceSurfTempN(i,j) == UNDEFVAL ) xy_SIceSurfTempA(i,j) = 0d0
             if ( xy_IceThickN(i,j)/xy_SIceConN(i,j) < IceThickMin ) then
                write(*,*) "i,j=", i,j, xy_IceThickN(i,j), xy_SIceConN(i,j)
                call MessageNotify('E', module_name, "IceThick < IceThickMin")
             end if
             call advance_ThermDyncProc_SIceTemp_zColum( &
                  & xy_SIceSurfTempA(i,j), xyz_SIceTempA(i,j,1), xyz_SIceTempA(i,j,2),        &  ! (out)  
                  & xy_SurfMeltEn(i,j),                                                       &  ! (out)
                  & xy_BtmMeltEn(i,j),                                                        &  ! (inout)
                  & xy_SnowThickN(i,j)/xy_SIceConN(i,j), xy_IceThickN(i,j)/xy_SIceConN(i,j),  &  ! (in)
                  & xy_SurfHFlxAIN(i,j),   xy_DSurfHFlxDTsAIN(i,j), xy_PenSWRFlxSIN(i,j),     &  ! (in)
                  & DelTime                                                                   &  ! (in)
                  & )
          else if ( xy_IceThickN(i,j) > 0d0 ) then
             write(*,*) "SIceCon <= IceMakMin, but IceThick > 0 !!", &
                  "j=", j, "hi=", xy_IceThickN(i,j), "T=", &
                  xy_SIceSurfTempN(i,j), xyz_SIceTempN(i,j,:)
          else
             xy_SIceSurfTempA(i,j) = UNDEFVAL; xyz_SIceTempA(i,j,:) = UNDEFVAL;
             xy_SurfMeltEn(i,j) = 0d0
          end if

          !
          !
          if(j==DEBUG_j) then
             write(*,*) " - SurfTempA=", xy_SIceSurfTempA(i,j), "SIceTempA=", xyz_SIceTempA(i,j,:)
             write(*,*) " - SurfMeltEn=", xy_SurfMeltEn(i,j), "BtmMeltEn=", xy_BtmMeltEn(i,j)
          end if
          
       end do
    end do

    do j=1, jMax
       do i=0, iMax-1
          if (isNan(xy_SIceSurfTempA(i,j)) .or. isNan(xyz_SIceTempA(i,j,1)) .or. isNan(xyz_SIceTempA(i,j,2)) ) then
             write(*,*) "After calc temp."
             write(*,*) "Nan T1: i,j=", i,j
             write(*,*) "SurfTemp=", xy_SIceSurfTempA(i,:)
             write(*,*) "Temp1=", xyz_SIceTempA(i,:,1)
             write(*,*) "Temp2=", xyz_SIceTempA(i,:,2) 
             write(*,*) "IceThick=", xy_IceThickA(i,:)
             write(*,*) "SIceCon=", xy_SIceConA(i,:)
             stop
          end if
       end do
    end do

    !
    !
    
    !$omp parallel do private(i, SurfTempO, FreezeEnOpnOcn) schedule(guided)
    do j=1, jMax
       do i=0, iMax-1

          SurfTempO = xyz_PTempO(i,j,0)
          xy_IceThickA(i,j) = xy_IceThickN(i,j)
          xy_SnowThickA(i,j) = xy_SnowThickN(i,j)
          FreezeEnOpnOcn = xy_FreezePot(i,j) + xy_SurfHFlxAON(i,j)
          
          call advance_SeaIceThermDynProc_SIceThick_zColum( &
               & xy_SnowThickA(i,j), xy_IceThickA(i,j),                                          &  ! (inout)
               & xy_SIceSurfTempA(i,j), xyz_SIceTempA(i,j,1), xyz_SIceTempA(i,j,2), SurfTempO,   &  ! (inout)
               & xy_Wice(i,j), xy_ExcessMeltEn(i,j),                                             &  ! (inout)
               & xy_SurfMeltEn(i,j), xy_BtmMeltEn(i,j), FreezeEnOpnOcn,                          &  ! (in)
               & xy_Wsnow(i,j), xy_Wevap(i,j), DelTime, i, j )                                                     ! (in)


          !* Replace SST below freezing point with Tfreeze
          ! if(SurfTempO /= xyz_PTempO(i,j,0)) then
          ! 
          ! end if

          if(j==DEBUG_j) then
             write(*,*) " -- TempO=", -273.15d0+(xyz_PTempO(i,j,0:1)), " SaltO=", xyz_SaltON(i,j,0)                          
             write(*,*) " -- IceThickA=", xy_IceThickA(i,j), ", SnowThickA=", xy_SnowThickA(i,j)
          end if
          
       end do
    end do

    do j=1, jMax
       do i=0, iMax-1
          if (isNan(xy_SIceSurfTempA(i,j)) .or. isNan(xyz_SIceTempA(i,j,1)) .or. isNan(xyz_SIceTempA(i,j,2)) ) then
             write(*,*) "After calc thickness.."
             write(*,*) "Nan T1: i,j=", i,j  
             write(*,*) "SurfTemp=", xy_SIceSurfTempA(i,:)
             write(*,*) "Temp1=", xyz_SIceTempA(i,:,1)
             write(*,*) "Temp2=", xyz_SIceTempA(i,:,2) 
             write(*,*) "IceThick=", xy_IceThickA(i,:)
             write(*,*) "SIceCon=", xy_SIceConA(i,:)
             stop
          end if
       end do
    end do

!!$    do j=1, jMax
!!$       write(*,*) "j=", j , "IceThick", xy_IceThickA(i,j), "SnowThick", xy_SnowThicka(i,j), &
!!$            "SurfTemp=", xy_SIceSurfTempA(i,j), "IceTep=", xyz_SIceTempA(i,j,:)
!!$    end do
!!$    do  j=1, jMax/2
!!$       write(*,*) "*** j=", j
!!$       write(*,*) "SnowMassDel:", DensSnow*(xy_SnowThickA(0,j) - xy_SnowThickN(0,j))
!!$       write(*,*) "IceMassDel:", DensIce*(xy_IceThickA(0,j) - xy_IceThickN(0,j))
!!$       write(*,*) "Wice:", xy_Wice(0,j)*DelTime*1d3
!!$    end do
!!$stop

    !* Dynamical processes 
    !
    call advance_SeaIceAdvectProc( & 
         xy_SIceConA, xy_IceThickA, xy_SnowThickA, xyz_SIceTempA, xy_SIceSurfTempA, & ! (inout)
         DelTime )                                                                    ! (in)

    
    !* Adjust the concentration and thickness of sea ice
    !
    call adjust_SeaIceField( &
         xy_SIceConA, xy_IceThickA, xy_SnowThickA, xyz_SIceTempA, xy_SIceSurfTempA, & ! (inout)
         xy_ExcessMeltEn, xy_Wice,                                                  & ! (inout)
         DelTime )
    
    !* Update heat flux and freshwater flux at the bottom of sea ice
    !
    
    !$omp parallel do private(i)
    do j=1, jMax
       do i=0, iMax-1
          xy_SurfHFlxIO(i,j)  = xy_ExcessMeltEn(i,j)/DelTime
          xy_SurfFwFlxIO(i,j) = - xy_Wice(i,j)
          
          if( xy_SIceConA(i,j) >= IceMaskMin ) then

             xy_SurfHFlxIO(i,j) = xy_SurfHFlxIO(i,j) &
                  + (1d0 - xy_SIceConA(i,j)) * (   min(xy_SurfHFlxAON(i,j), -xy_FreezePot(i,j)) &
                                                       + DensSnow*LFreeze*xy_Wsnow(i,j) ) &
                  + xy_SIceConA(i,j) * xy_BtmHFlxIO(i,j)

             xy_SurfFwFlxIO(i,j) = xy_SurfFwFlxIO(i,j) &
                  + (1d0 - xy_SIceConA(i,j)) * (xy_Wrain(i,j) + xy_Wsnow(i,j) - xy_Wevap(i,j))        &
                  + xy_SIceConA(i,j) * xy_Wrain(i,j)
          end if
       end do
    end do
    
!!$    write(*,*) "ConN:", xy_SIceConN(0,:)
!!$    write(*,*) "ConA:", xy_SIceConA(0,:)
!!$    write(*,*) "xy_SurfHFlxIO=", xy_SurfHFlxIO(0,:)
!!$    write(*,*) "------------------------------------"
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
       & TsA, T1A, T2A, SurfMeltEn,                        & ! (inout)
       & BtmMeltEn,                                        & ! (inout)
       & hs, hi,                                           & ! (in)
       & SurfHFlxAIN, DSurfHFlxDTs, PenSWRFlxSIN,          & ! (in)
       & DelTime                & ! (in)
       & )

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         & update_effConductiveCoupling, calc_layersTemprature

    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: TsA, T1A, T2A
    real(DP), intent(out) :: SurfMeltEn, BtmMeltEn
    real(DP), intent(in) :: hs, hi
    real(DP), intent(in) :: SurfHFlxAIN, DSurfHFlxDTs, PenSWRFlxSIN
    real(DP), intent(in) :: DelTime
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: Ts_, T1_, T2_
    real(DP) :: K_12, K_23
    real(DP) :: A, B
    
    
    ! 実行文; Executable statement
    !

    Ts_ = TsA; T1_ = T1A; T2_ = T2A
    
    call update_effConductiveCoupling(K_12, K_23, & ! (out)
         & hs, hi )                                 ! (in)

    B = DSurfHFlxDTs
    A = SurfHFlxAIN - Ts_*B

    !
    call calc_layersTemprature( &
         & TsA, T1A, T2A, SurfMeltEn,                &  ! (out)
         & BtmMeltEn,                                &  ! (inout)
         & T1_, T2_, hi, DelTime, A, B,                &  ! (in)
         & -PenSWRFlxSIN, K_12, K_23, hs > 0d0      &  ! (in)
         & )
    
  end subroutine advance_ThermDyncProc_SIceTemp_zColum

  !>
  !! @param FWFlx freshwater flux generated by 
  !! @param SurfHFxAIN surface heat flux at time level n(**positive upward**). [W/m2]
  !! @param 
  subroutine advance_SeaIceThermDynProc_SIceThick_zColum( &
       & hsEffA, hiEffA, TsA, T1A, T2A, TsO,              &  ! (inout)
       & Wice, excessMeltEn,                              &  ! (inout)
       & SurfMeltEn, BtmMeltEn, FreezeEnOpnOcn,           &  ! (in)
       & SnowFall, Sublim, DelTime, i, j                          &  ! (in)
       & )

    ! モジュール引用; Use statements
    !
    use SeaIceThermDyn_Winton2000_mod, only: &
         calc_SnowIceLyrMassChange, adjust_IceLyrInternal, &
         calc_E_IceLyr1, calc_E_IceLyr2
         
    use VarSetSeaice_mod, only: &
         xy_SIceConN, xy_SIceConA
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: hsEffA, hiEffA
    real(DP), intent(inout) :: Wice, excessMeltEn
    real(DP), intent(inout) ::  TsA, T1A, T2A, TsO
    real(DP), intent(in) :: SurfMeltEn
    real(DP), intent(in) :: BtmMeltEn
    real(DP), intent(in) :: FreezeEnOpnOcn    
    real(DP), intent(in) :: SnowFall
    real(DP), intent(in) :: Sublim
    real(DP), intent(in) :: DelTime
    integer :: i, j
    
    ! 局所変数
    ! Local variables
    !
    real(DP) :: dhs, dh1, dh2
    real(DP) :: hs, hi, IceFrac
    real(DP) :: excessMeltEn_
    real(DP) :: GiceEff, GsnowEff, GiceFrac, Gice, Gsnow, Gsnow2ice

    real(DP) :: GiceOpnOcn
    
    real(DP), parameter :: hIceEPS = 1d-12
    
    ! 実行文; Executable statements
    !

    ! 
    !
    IceFrac = xy_SIceConN(i,j)
    hi = hiEffA / max(iceFrac, IceMaskMin)
    hs = hsEffA / max(iceFrac, IceMaskMin)

    Gice = 0d0; Gsnow = 0d0
    Gsnow2ice = 0d0
    GiceOpnOcn = max(0d0, FreezeEnOpnOcn)/(- DensIce * calc_E_IceLyr2(FreezeTempSW, SaltSeaIce))
    

    if (iceFrac >= IceMaskMin) then
       !
       call calc_SnowIceLyrMassChange( &
            & dhs, dh1, dh2, excessMeltEn_,       & ! (out)
            & T1A, T2A, TsO,                      & ! (inout)
            & hs, hi, TsA, DelTime,               & ! (in)
            & SurfMeltEn, BtmMeltEn               & ! (in)
            & )


       Gice = (dh1 + dh2)/DelTime
       Gsnow = dhs/DelTime + SnowFall* 1d3 / DensSnow ! - Sublim    
       excessMeltEn = excessMeltEn &
            + iceFrac * excessMeltEn_    

       Wice = Wice &
            + iceFrac * (DensIce*Gice + DensSnow*min(0d0, dhs/DelTime))/1d3! DensSeaWater
    
       Gsnow = Gsnow - Sublim * 1d3 / DensSnow
       if (hs <  - Gsnow*DelTime) then
          Gice = Gice + (hs/DelTime + Gsnow) * DensSnow/DensIce
          Gsnow = - hs/DelTime
       end if

       hs = hs + Gsnow*DelTime
    end if
 
   if(j==DEBUG_j) then
       write(*,*) "T1, T2=", T1A, T2A
       write(*,*) "dhs, dh1, dh2=", dhs, dh1, dh2
       write(*,*) "SnowFall=", SnowFall*DelTime, "Sublim=", Sublim*DelTime
    end if

    !
    !    
    if ( hi == 0d0 .and. GiceOpnOcn > 0d0 ) then
       T1A = - Mu*SaltSeaIce
       T2A = FreezeTempSW
       hi = GiceOpnOcn*DelTime
       hs = 0d0
       call adjust_IceLyrInternal( &
            & hs, hi, T1A, T2A,                    & !(inout)
            & 0d0, hi                              & !(in)
            & )
       Wice = Wice &
            + (1d0 - iceFrac) * DensIce*GiceOpnOcn/1d3

    else if ( hi + Gice*DelTime > 0d0 ) then
       call adjust_IceLyrInternal( &
            & hs, hi, T1A, T2A,                    & !(inout)
            & 0.5d0*hi + dh1, 0.5d0*hi + dh2       & !(in)
            & )
       
       Gsnow2ice = ((hsEffA/iceFrac + Gsnow*DelTime) - hs)/DelTime
       Gsnow = Gsnow - Gsnow2ice 
       Gice  = Gice  + DensSnow/DensIce*Gsnow2ice

       if (iceFrac < IceMaskMin) then
          call MessageNotify('E', module_name, &
               "Ice fraction must be larger than IceMaskMin. Unexpected error have occured!")
       end if
    end if

    if(j==DEBUG_j) then
       write(*,*) "Gsnow*dt=", Gsnow*DelTime 
       write(*,*) "GSnow2Ice*dt=", Gsnow2ice*DelTime
       write(*,*) "reshape_IceLyr => T1, T2=", T1A, T2A
       write(*,*) "GiceOpnOcn=", GiceOpnOcn
       write(*,*) "hice=", hi, "hiceEff=", hiEffA, "Gice*dt", Gice*DelTime, "iceFrac*Gice*dt", Gice*DelTime*iceFrac
!       write(*,*) "FreezeEnOpnOcn=", FreezeEnOpnOcn
       write(*,*) "Wsnow", SnowFall*DelTime*1d3, "WbtmMel", dh2*DensIce
    end if

    !
    !
    
    GiceEff  = iceFrac*Gice + (1d0 - iceFrac)*GiceOpnOcn
    GsnowEff = iceFrac*Gsnow
    GiceFrac =   (1d0 - iceFrac) * GiceOpnOcn/0.2d0 &
               + iceFrac * 0.5d0 * min(0d0, Gice) / max(hiEffA, hIceEPS)

    Wice = Wice &
         + (1d0 - iceFrac) * (DensIce*GiceOpnOcn) / 1d3 !DensSeaWater
    
!    xy_SIceConA(i,j) = iceFrac + GiceFrac*DelTime
    xy_SIceConA(i,j) = &
         (iceFrac + GiceOpnOcn / 0.2d0 * DelTime) &
         / (1d0 + DelTime*( GiceOpnOcn / 0.2d0 - 0.5d0 * min(0d0, Gice) / max(hiEffA, hIceEPS)) )

    hiEffA = hiEffA + GiceEff*DelTime
    hsEffA = hsEffA + GsnowEff*DelTime

    if (isNan(hiEffA) .or. isNan(TsA) .or. isNan(T1A) .or. isNan(T2A)) then
       write(*,*) "-- hiEff is Nan. j=", j
       write(*,*) " hs=", hs, "hi=", hi, "T=", TsA, T1A, T2A
       write(*,*) " Gice=", Gice, "Gsnow=", Gsnow, "Gsnow2ice=", Gsnow2ice, "iceFrac=", iceFrac
       write(*,*) " dh1=", dh1, "dh2=", dh2
    end if

    
!    if(hiEffA <= 0d0 .and. hiEffA > 0d0) xy_SIceConA(i,j) = 1d0
    
  end subroutine advance_SeaIceThermDynProc_SIceThick_zColum

  subroutine advance_SeaIceAdvectProc(          & 
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xya_SIceTemp, xy_SIceSurfTemp, & ! (inout)
       & DelTime )                                                               ! (in)

    use GridSet_mod, only: &
         & xyz_Lat
    
    use SeaIceThermDyn_Winton2000_mod, only: &
         & calc_E_IceLyr1, calc_E_IceLyr2, &
         & calc_Temp_IceLyr1, calc_Temp_IceLyr2

    use SpmlUtil_mod


    real(DP), intent(inout) :: xy_SIceCon(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_IceThick(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SnowThick(0:iMax-1,jMax)
    real(DP), intent(inout) :: xya_SIceTemp(0:iMax-1,jMax,2)
    real(DP), intent(inout) :: xy_SIceSurfTemp(0:iMax-1,jMax)
    real(DP), intent(in) :: DelTime

    real(DP) :: xy_q1(0:iMax-1,jMax)
    real(DP) :: xy_q2(0:iMax-1,jMax)

    real(DP) :: iceVol
    real(DP) :: Ice0En1, Ice0En2
    real(DP) :: xy_CosLat(0:iMax-1,jMax)

    integer :: i, j

!!$    write(*,*) "FracB:", xy_SIceCon(0,:)
!!$    write(*,*) "q2B:", xya_SIceTemp(0,:,2)
    
    xy_CosLat = cos(xyz_Lat(:,:,0))

    !$omp parallel do private(i)
    do j=1, jMax
       do i=0, iMax-1
          if (xy_SIceCon(i,j) >= IceMaskMin) then
             xy_q1(i,j) = calc_E_IceLyr1(xya_SIceTemp(i,j,1), SaltSeaIce)
             xy_q2(i,j) = calc_E_IceLyr2(xya_SIceTemp(i,j,2), SaltSeaIce)
          else
             xy_q1(i,j) = 0d0
             xy_q2(i,j) = - LFreeze
          end if

          if ( j == DEBUG_j ) then
             write(*,*) "A*=", xy_SIceCon(i,j)
             write(*,*) "q1*, q2*=", xy_q1(i,j), xy_q2(i,j)
          end if
       end do
    end do

!write(*,*) "SiceDiff Before:", AvrLonLat_xy( xy_IceThick) 
    call apply_diff_term( xy_IceThick  )  ! (inout)
    call apply_diff_term( xy_SnowThick )  ! (inout)
    call apply_diff_term( xy_SIceCon   )  ! (inout)
    call apply_diff_term( xy_q1        )  ! (inout)
    call apply_diff_term( xy_q2        )  ! (inout)
!write(*,*) "SiceDiff After:", AvrLonLat_xy( xy_IceThick) 

    
    !$omp parallel do private(i)
    do j=1, jMax
       do i=0, iMax-1 
          if (xy_SIceCon(i,j) > 0d0) then
             xya_SIceTemp(i,j,1) = calc_Temp_IceLyr1(xy_q1(i,j), SaltSeaIce)
             xya_SIceTemp(i,j,2) = calc_Temp_IceLyr2(xy_q2(i,j), SaltSeaIce)
             if (isNan(xya_SIceTemp(i,j,1)) .or. isNan(xya_SIceTemp(i,j,2)) ) then
                write(*,*) "Nan T1: i,j=", i,j  
                write(*,*) "q1=", xy_q1(i,:)
                write(*,*) "q2=", xy_q1(i,:) 
                write(*,*) "IceThick=", xy_IceThick(i,:)
                write(*,*) "SIceCon=", xy_SIceCon(i,:)
                stop
             end if
          else
             xya_SIceTemp(i,j,:) = UNDEFVAL
          end if

          if ( j == DEBUG_j ) then
             write(*,*) "q1**, q2**=", xy_q1(i,j), xy_q2(i,j)
             write(*,*) "T1**, T2**=", xya_SIceTemp(i,j,:)
          end if          
       end do
    end do

!!$    write(*,*) "FracA:", xy_SIceCon(0,1)    
!!$    write(*,*) "q2A:", xya_SIceTemp(0,:,2)
!!$    stop
  contains
    subroutine apply_diff_term(xy_phi)

      use GridSet_mod, only: &
           & x_Lon_Weight, y_Lat_Weight
      
      real(DP), intent(inout) :: xy_phi(0:iMax-1,jMax)
      
      real(DP) :: xy_FlxLat(0:iMax-1,0:jMax)
      real(DP) :: DLon
      integer :: i, j
      real(DP), parameter :: EPSILL = 1d-12
      
      DLon = 2d0 * PI / dble(iMax)

      xy_FlxLat(:,0) = 0d0; xy_FlxLat(:,jMax) = 0d0

      !$omp parallel do private(i)
      do j=1, jMax-1
         do i=0, iMax-1
            xy_FlxLat(i,j) = &
                  RPlanet*cos( 0.5d0 * (xyz_Lat(i,j+1,0) + xyz_Lat(i,j,0)) )*DLon      &
               *  SIceHDiffCoef *   (xy_phi(i,j+1) - xy_phi(i,j))                      & 
                                  / ( RPlanet * (xyz_Lat(i,j+1,0) - xyz_Lat(i,j,0)) )
         end do
      end do

      !$omp parallel do private(i)
      do j=1, jMax
         do i=0, iMax-1
            xy_phi(i,j) = xy_phi(i,j) &
                + DelTime * (xy_FlxLat(i,j) - xy_FlxLat(i,j-1)) &
                          / ( RPlanet**2 * x_Lon_Weight(i) * y_Lat_Weight(j) )
         end do
      end do
      
    end subroutine apply_diff_term
    
  end subroutine advance_SeaIceAdvectProc   

  subroutine adjust_SeaIceField( &
       xy_SIceCon, xy_IceThickEff, xy_SnowThickEff, xya_SIceTemp, xy_SIceSurfTemp, & ! (inout)
       xy_ExcessMeltEn, xy_Wice,                                                   & ! (inout)
       DelTime )
    
    use SeaIceThermDyn_Winton2000_mod, only: &
         & calc_E_IceLyr1, calc_E_IceLyr2, &
         & calc_Temp_IceLyr1, calc_Temp_IceLyr2

    use SpmlUtil_mod

    use VarSetSeaice_mod, only: &
         xya_SIceEnA

    real(DP), intent(inout) :: xy_SIceCon(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_IceThickEff(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SnowThickEff(0:iMax-1,jMax)
    real(DP), intent(inout) :: xya_SIceTemp(0:iMax-1,jMax,2)
    real(DP), intent(inout) :: xy_SIceSurfTemp(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_ExcessMeltEn(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_Wice(0:iMax-1,jMax)
    real(DP), intent(in) :: DelTime
    
    integer :: i, j
    real(DP) :: iceFrac, E1, E2
    real(DP) :: hs, hi
    real(DP) :: xy_Wice2(0:iMax-1,jMax)

!    write(*,*) "Adjust IceThick:", xy_IceThickEff(0,:)
!    write(*,*) "Before Adjust:", AvrLonLat_xy(xy_IceThickEff*DensIce + xy_SnowThickEff*DensSnow)
!    xy_Wice2 = 0d0

    !$omp parallel do private(i, iceFrac, E1, E2, hs, hi) schedule(guided)
    do j=1 , jMax
       do i=0, iMax-1

          if ( xy_IceThickEff(i,j) <= 0d0 ) then
             xy_SIceSurfTemp(i,j) = UNDEFVAL
             xya_SIceTemp(i,j,:) = UNDEFVAL
             xy_IceThickEff(i,j) = 0d0
             xy_SnowThickEff(i,j) = 0d0
             xy_SIceCon(i,j) = 0d0
             cycle
          end if
          
          if ( xy_SIceCon(i,j) < IceMaskMin ) then
             !xy_IceThickEff(i,j) = xy_SIceCon(i,j)/IceMaskMin * xy_IceThickEff(i,j) 
             !xy_SnowThickEff(i,j) = xy_SIceCon(i,j)/IceMaskMin * xy_SnowThickEff(i,j) 
             xy_SIceCon(i,j) = IceMaskMin
          end if

          hi = xy_IceThickEff(i,j) / xy_SIceCon(i,j)
          E1 = calc_E_IceLyr1(xya_SIceTemp(i,j,1), SaltSeaIce)
          E2 = calc_E_IceLyr2(xya_SIceTemp(i,j,2), SaltSeaIce)
          if ( hi < IceThickMin ) then
             if(xy_IceThickEff(i,j) >= IceThickMin*IceMaskMin) then
                xy_SIceCon(i,j) = xy_IceThickEff(i,j) / IceThickMin
!                xy_SnowThickEff(i,j) = xy_SnowThickEff(i,j) * IceThickMin/hi
!                xy_IceThickEff(i,j) = IceThickMin * IceMaskMin
             else                
                xy_ExcessMeltEn(i,j) = xy_ExcessMeltEn(i,j) + (                 &
                      - xy_IceThickEff(i,j) * 0.5d0 * DensIce * (E1 + E2)          &
                      + xy_SnowThickEff(i,j) * DensSnow * LFreeze                  &
                      ) 

                xy_Wice(i,j) = xy_Wice(i,j) - (                                       &
                     DensIce * xy_IceThickEff(i,j) + DensSnow * xy_SnowThickEff(i,j)  &
                     )  / ( 1d3 * DelTime ) !DensSeaWater * DelTime )
!!$                xy_Wice2(i,j) = - (                                       &
!!$                     DensIce * xy_IceThickEff(i,j) + DensSnow * xy_SnowThickEff(i,j)  &
!!$                     )  !*0d0

!!$                write(*,*) "CALL.. ", i, j, xy_ExcessMeltEn(i,j) / DelTime, xy_Wice(i,j), &
!!$                     xy_IceThick(i,j), xy_SnowThick(i,j), xya_SIceTemp(i,j,:), &
!!$                     xy_SIceSurfTemp(i,j)
                
                xy_IceThickEff(i,j) = 0d0
                xy_SnowThickEff(i,j) = 0d0
                xy_SIceCon(i,j) = 0d0
                xya_SIceTemp(i,j,:) = UNDEFVAL
                xy_SIceSurfTemp(i,j) = UNDEFVAL
                
                xya_SIceEnA(i,j,:) = 0d0
                cycle
             end if
          end if

          xya_SIceEnA(i,j,1) = 0.5d0*E1*xy_IceThickEff(i,j)*DensIce
          xya_SIceEnA(i,j,2) = 0.5d0*E2*xy_IceThickEff(i,j)*DensIce
          if ( xy_SIceSurfTemp(i,j) == UNDEFVAL ) then
             xy_SIceSurfTemp(i,j) = 0d0
          end if
!!$          if ( xya_SIceTemp(i,j,1) == UNDEFVAL ) then
!!$             write(*,*) "Warning.. i,j=", i, j
!!$             write(*,*) "------", xy_IceThick(i,j), xy_SnowThick(i,j), xya_SIceTemp(i,j,:), &
!!$                  xy_SIceSurfTemp(i,j)
!!$             stop
!!$          end if
       end do
    end do

!    write(*,*) "After Adjust:", AvrLonLat_xy(xy_IceThickEff*DensIce + xy_SnowThickEff*DensSnow)
!stop        write(*,*) "Wice2", AvrLonLat_xy(xy_Wice2)

  end subroutine adjust_SeaIceField
  
 end module SeaIceEq_TimeInteg_mod

