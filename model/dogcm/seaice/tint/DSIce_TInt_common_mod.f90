!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_TInt_common_mod

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

  use DOGCM_Admin_Constants_mod, only: &
       & UNDEFVAL,                     &
       & CpOcn => Cp0,                 &
       & DensSeaWater => RefDens
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM
        
  use DSIce_Admin_Constants_mod, only: &
       & DensFreshWater,               &
       & DensIce, DensSnow,            &
       & Mu, SaltSeaIce,               &
       & LFreeze,                      &
       & IceMaskMin, IceThickMin
  
  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs,   &
       & xy_PenSDRFlx,                     &
       & xy_SeaSfcTemp, xy_SeaSfcSalt,     &
       & xy_OcnFrzTemp, xy_BtmHFlxIO,      &
       & xy_OcnMixLyrDepth,                &
       & xy_RainFall, xy_SnowFall, xy_Evap
       

  use DSIce_Dyn_driver_mod, only: &
       & DSIce_Dyn_driver_Do
  
  use DSIce_ThermoDyn_Winton2000_mod, only: &
        & DSIce_ThermoDyn_Winton2000_CalcIceTemp,       &
        & DSIce_ThermoDyn_Winton2000_CalcLyrMassChange, &
        & DSIce_ThermoDyn_Winton2000_AdjustLyrInternal, &
        & calc_E_IceLyr1, calc_E_IceLyr2
  

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_TInt_common_Init, DSIce_TInt_common_Final

  public :: DSIce_TInt_common_advance_Dyn
  public :: DSIce_TInt_common_advance_ThermoDyn
  public :: DSIce_TInt_common_advance_Phys
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_TInt_common_mod' !< Module Name


contains

  !>
  !!
  !!
  Subroutine DSIce_TInt_common_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)

    
  end subroutine DSIce_TInt_common_Init

  !>
  !!
  !!
  subroutine DSIce_TInt_common_Final()

    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statements
    !


  end subroutine DSIce_TInt_common_Final

  !-----------------------------------------------------------------------------------

  subroutine DSIce_TInt_common_advance_Phys( &
       & dt                                                      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: dt

    ! 実行文; Executable statements
    !
    
    != Physical process ======================================================

    
  end subroutine DSIce_TInt_common_advance_Phys

  !----------------------------------------------------------------------
  
  subroutine DSIce_TInt_common_advance_Dyn(  & 
       & dt                                                  & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !

    real(DP), intent(in) :: dt
    

    
    ! 作業変数
    ! Work variables
    
    
    ! 実行文; Executable statements
    !

    

  end subroutine DSIce_TInt_common_advance_Dyn

  !----------------------------------------------------------------------
  
  subroutine DSIce_TInt_common_advance_ThermoDyn(  &
       & xy_SIceConA, xy_IceThickA, xy_SnowThickA, xy_SIceSfcTempA, xyz_SIceTempA,  & ! (out)
       & xy_Wice,                                                                   & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0, xy_SIceSfcTemp0, xyz_SIceTemp0,  & ! (in)
       & dt                                                                         & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceConA(IA,JA)
    real(DP), intent(out) :: xy_IceThickA(IA,JA)
    real(DP), intent(out) :: xy_SnowThickA(IA,JA)
    real(DP), intent(out) :: xy_SIceSfcTempA(IA,JA)
    real(DP), intent(out) :: xyz_SIceTempA(IA,JA,KA)
    real(DP), intent(out) :: xy_Wice(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp0(IA,JA)
    real(DP), intent(in) :: xyz_SIceTemp0(IA,JA,KA)
    real(DP), intent(in) :: dt
    
    
    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j
    real(DP) :: SfcResHFlx
    real(DP) :: BtmResHFlx
    real(DP) :: SfcFrzTemp
    real(DP) :: z_SIceTempA(KS:KE)
    real(DP) :: z_SIceTemp0(KS:KE)
    real(DP) :: dSnowThick
    real(DP) :: dIceThick1
    real(DP) :: dIceThick2
    real(DP) :: FrzPot
    real(DP) :: excessMeltEn
    real(DP) :: OcnFrzTemp

    integer :: j_dbg 
    
    ! 実行文; Executable statements
    !

    j_dbg = - 5 !JS + 4
    
    !$omp parallel do private( &
    !$omp   i, SfcResHFlx, BtmResHFlx, z_SIceTempA, z_SIceTemp0, SfcFrzTemp, OcnFrzTemp,   &
    !$omp   dSnowThick, dIceThick1, dIceThick2, excessMeltEn, FrzPot                       &
    !$omp ) schedule(guided)
    do j = JS, JE
       do i = IS, IE

          z_SIceTempA(:) = xyz_SIceTempA(i,j,KS:KE) 
          z_SIceTemp0(:) = xyz_SIceTemp0(i,j,KS:KE) 

          SfcFrzTemp = 0d0
          if (xy_SnowThick0(i,j) == 0d0) SfcFrzTemp = - Mu*SaltSeaIce

          OcnFrzTemp = K2degC( xy_OcnFrzTemp(i,j) )
          FrzPot = UNDEFVAL
          
          !---------------------------

          if ( xy_SIceCon0(i,j) >= IceMaskMin ) then 
             call DSIce_ThermoDyn_Winton2000_CalcIceTemp( &
                  & xy_SIceSfcTempA(i,j), z_SIceTempA(KS:KE),            & ! (out)
                  & SfcResHFlx, BtmResHFlx,                              & ! (out)
                  & xy_SIceSfcTemp0(i,j), z_SIceTemp0(KS:KE),            & ! (in)
                  & xy_IceThick0(i,j), xy_SnowThick0(i,j),               & ! (in)
                  & xy_SfcHFlxAI(i,j), xy_DSfcHFlxAIDTs(i,j),            & ! (in)
                  & xy_PenSDRFlx(i,j), SfcFrzTemp,                       & ! (in)
                  & xy_BtmHFlxIO(i,j), OcnFrzTemp,                       & ! (in)
                  & dt,  (i==IS .and. j==j_dbg) )

             if (i==IS .and. j==j_dbg) then
                write(*,*) "-calcTemp----------------------"
                write(*,*) "SfcResHFlx=", SfcResHFlx, "BtmResHFlx=", BtmResHFlx
                write(*,*) "TempN=", xy_SIceSfcTemp0(i,j), xyz_SIceTemp0(i,j,KS:KE)
                write(*,*) "TempA=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)
              end if
             
             !---------------------------
          
             call DSIce_ThermoDyn_Winton2000_CalcLyrMassChange( &
                  & dSnowThick, dIceThick1, dIceThick2,         & ! (out)
                  & xy_Wice(i,j), excessMeltEn,                 & ! (out)
                  & z_SIceTempA(KS:KE),                         & ! (inout)
                  & xy_SnowThick0(i,j), xy_IceThick0(i,j),      & ! (in)
                  & SfcResHFlx, BtmResHFlx, OcnFrzTemp,         & ! (in)
                  & dt,  (i==IS .and. j==j_dbg) )                 ! (in)

             xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) - excessMeltEn/dt

             dSnowThick = dSnowThick + xy_SnowFall(i,j)
             

             if (i==IS .and. j==j_dbg) then
                write(*,*) "-calcLyrMassChange----------------------"
                write(*,*) "dhs=", dSnowThick, "dhi1=", dIceThick1, "dhi2=", dIceThick2
                write(*,*) "Temp=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)                
                write(*,*) "excessMeltEn=", excessMeltEn
                write(*,*) "-------------------------"
             end if
             
             !---------------------------

             call DSIce_ThermoDyn_Winton2000_AdjustLyrInternal( &
                  & xy_SnowThickA(i,j), xy_IceThickA(i,j),               & ! (out)
                  & z_SIceTempA(KS:KE),                                  & ! (inout)
                  & xy_SnowThick0(i,j), dSnowThick,                      & ! (in)
                  & xy_IceThick0(i,j), dIceThick1, dIceThick2            & ! (in)
                  & )
!!$             if (i==IS .and. j==JS) then
!!$                write(*,*) "-adjustLyrInternal ----------------------"
!!$                write(*,*) "dhs=", dSnowThick, "dhi1=", dIceThick1, "dhi2=", dIceThick2
!!$                write(*,*) "excessMeltEn=", excessMeltEn, "Temp=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)
!!$                write(*,*) "-------------------------"
!!$             end if

             !---------------------------

             if ( xy_IceThickA(i,j) <= 0d0 ) then
                xy_SIceConA(i,j) = 0d0
                xy_SIceSfcTempA(i,j) = UNDEFVAL
                z_SIceTempA(KS:KE) = UNDEFVAL
             end if
             
          else

             !- Initalize sea ice ------------------------------------------------

             FrzPot =    CpOcn*DensSeaWater*xy_OcnMixLyrDepth(i,j)/dt         &
                  &    * ( OcnFrzTemp - K2degC(xy_SeaSfcTemp(i,j)) )
             if ( FrzPot > 0d0 ) then
                xy_IceThickA(i,j) = &
                     & FrzPot * dt &
                     &  /( - DensIce*calc_E_IceLyr2(OcnFrzTemp,SaltSeaIce) )
                xy_Wice(i,j) = DensIce*xy_IceThickA(i,j)/dt
                
                z_SIceTempA(KS:KE)   = OcnFrzTemp
                xy_SIceSfcTempA(i,j) = OcnFrzTemp
                xy_SIceConA(i,j) = 1d0
                xy_BtmHFlxIO(i,j) = - FrzPot

             else
                xy_SIceSfcTempA(i,j) = UNDEFVAL
                z_SIceTempA(KS:KE) = UNDEFVAL
             end if

          end if

!          if( xy_IceThickA(i,j) > 1d3) j_dbg = j

          if (xy_IceThickA(i,j) < IceThickMin) then
                          
             xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) + ( &
                  &         - DensIce*xy_IceThickA(i,j)*0.5d0*(                      &
                  &               calc_E_IceLyr1(xyz_SIceTempA(i,j,1),SaltSeaIce)    &
                  &             + calc_E_IceLyr2(xyz_SIceTempA(i,j,2),SaltSeaIce) )  &
                  &         + DensSnow*xy_SnowThickA(i,j)*LFreeze                    &
                  &         )/dt

             xy_Wice(i,j) =   xy_Wice(i,j) - (                                       &
                  &  DensIce*xy_IceThickA(i,j) + DensSnow*xy_SnowThickA(i,j)         &
                  &  )/dt

             xy_SIceConA(i,j) = 0d0
             xy_IceThickA(i,j) = 0d0
             xy_SnowThickA(i,j) = 0d0
             z_SIceTempA(KS:KE) = UNDEFVAL
          end if
          
          xyz_SIceTempA(i,j,KS:KE) = z_SIceTempA(KS:KE)

          !$omp critical
          if (i==IS .and. j==j_dbg) then
             write(*,*) "(i,j)=", i, j
             write(*,*) &
                  & "hiA=", xy_IceThickA(i,j), "hsA=", xy_SnowThickA(i,j),           &
                  & "TsfcA=", xy_SIceSfcTempA(i,j), "TA=", xyz_SIceTempA(i,j,KS:KE), &
                  & "SIceConA=", xy_SIceConA(i,j)
             write(*,*) &
                  & "hi0=", xy_IceThick0(i,j), "hs0=", xy_SnowThick0(i,j),           &
                  & "Tsfc0=", xy_SIceSfcTemp0(i,j), "T0=", xyz_SIceTemp0(i,j,KS:KE), &
                  & "SIceCon0=", xy_SIceCon0(i,j)
             
             write(*,*) "SST=", xy_SeaSfcTemp(i,j), "OcnFrzTemp=", OcnFrzTemp, &
                  & "FrzPot=", FrzPot/( - DensIce*calc_E_IceLyr2(OcnFrzTemp,SaltSeaIce) )*dt

             if( xy_IceThickA(i,j) > 1d3 ) stop
          end if
          !$omp end critical
          
       end do
    end do

!    write(*,*) "SnowFall=", xy_SnowFall(IS:IE,JS:JE)
!    write(*,*) "IceThickA=", xy_IceThickA(IS:IE,JS:JE)
!    write(*,*) "SnowThickA=", xy_SnowThickA(IS:IE,JS:JE)
    
  end subroutine DSIce_TInt_common_advance_ThermoDyn
  
end module DSIce_TInt_common_mod
