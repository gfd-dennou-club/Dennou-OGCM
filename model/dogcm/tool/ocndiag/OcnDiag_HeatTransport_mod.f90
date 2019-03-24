module OcnDiag_HeatTransport_mod

  ! モジュール引用; Use statement
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  use DOGCM_Admin_Constants_mod
  
  use GridIndex_mod, only: &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA

  use DOGCM_Admin_Grid_mod, only: &
       & xyz_Lat

  use DSIce_Admin_Constants_mod
  
  use DSIce_Admin_Grid_mod, only: &
       & ISS=>IS, IES=>IE, IAS=>IA, &
       & JSS=>JS, JES=>JE, JAS=>JA, &
       & KSS=>KS, KES=>KE, KAS=>KA

  use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar, &
       & DOGCM_IO_History_HistPut
  
  use SpmlUtil_mod, only: &
       & AvrLonLat_xy
  
  use CalculusDriver_mod, only: &
       & Calculus_IntBtmToTop

  ! 宣言文; Declareration statements
  !  
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !  
  public :: OcnDiag_HeatTransport_Init, OcnDiag_HeatTransport_Final

  public :: OcnDiag_HeatTransport_prepair_Output
  public :: OcnDiag_HTDiagnose
  public :: OcnDiag_SIceHTDiagnose
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'OcnDiag_HeatTransport_mod' !< Module Name

  real(DP), parameter :: Unit2PWFac = 1d-15
  
contains

  subroutine OcnDiag_HeatTransport_Init()
  end subroutine OcnDiag_HeatTransport_Init

  subroutine OcnDiag_HeatTransport_Final()
  end subroutine OcnDiag_HeatTransport_Final

  subroutine OcnDiag_HeatTransport_prepair_Output()

    call DOGCM_IO_History_RegistVar( 'EulerHT', 'JT', 'heat transport(eulerian contribution)', 'PW' )
    call DOGCM_IO_History_RegistVar( 'BolusHT', 'JT', 'heat transport(bolus contribution)', 'PW' )
    call DOGCM_IO_History_RegistVar( 'IsoDiffHT', 'JT', 'heat transport(isopycnal diffusion contribution)', 'PW' )
    call DOGCM_IO_History_RegistVar( 'OcnHT', 'JT', 'total heat transport(ocean)', 'PW' )    
    call DOGCM_IO_History_RegistVar( 'OcnHT_conv', 'JT', 'convergence of total heat transport(ocean)', 'W/m2' )    
    call DOGCM_IO_History_RegistVar( 'NumDiffTend', 'JT', &
      & 'tendency of ocean heat energy due to horizontal numerical filter', 'W/m2' )

    call DOGCM_IO_History_RegistVar( 'SIceHT', 'JT', 'heat transport(sea ice)', 'PW' )

    call DOGCM_IO_History_RegistVar( 'MSF_GM', 'JKT', 'meridional stream function due to bolus transport', 'Sv')
    
  end subroutine OcnDiag_HeatTransport_prepair_Output
  
  subroutine OcnDiag_SIceHTDiagnose( &
       & y_SIceHT,                                                  & ! (out)
       & xy_SIceU, xv_SIceV, xy_IceThick, xy_SnowThick, xyz_SIceEn  & ! (in)
       & )

    use DSIce_Admin_Grid_mod, only: xy_Lat

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: y_SIceHT(JAS)
    real(DP), intent(in) :: xy_SIceU(IAS,JAS)
    real(DP), intent(in) :: xv_SIceV(IAS,JAS)
    real(DP), intent(in) :: xy_IceThick(IAS,JAS)
    real(DP), intent(in) :: xy_SnowThick(IAS,JAS)
    real(DP), intent(in) :: xyz_SIceEn(IAS,JAS,KAS)

    real(DP) :: xy_SIceV(IAS,JAS)
    
    ! 実行文; Executable statement
    !
    
    xy_SIceV(:,JSS:JES) = 0.5d0*(xv_SIceV(:,JSS-1:JES-1) + xv_SIceV(:,JSS:JES))
    
    y_SIceHT(:) = IntLon( xy_SIceV*( &
         & + sum(xyz_SIceEn(:,:,KSS:KES),3) - DensSnow*xy_SnowThick*LFreeze &
         & ) )

    call DOGCM_IO_History_HistPut( 'SIceHT', y_SIceHT(JSS:JES)*Unit2PWFac )
    
  contains
    function IntLon(xy) result(y)
      real(DP), intent(in) :: xy(IAS,JAS)
      real(DP) :: y(JAS)

      y(:)   = sum( xy(ISS:IES,:) * 2d0*PI*RPlanet*cos(xy_Lat(ISS:IES,:)), 1)
    end function IntLon
  end subroutine OcnDiag_SIceHTDiagnose
  
  subroutine OcnDiag_HTDiagnose( &
       & y_EulerHT, y_BolusHT, y_IsoDiffHT, y_NumDiffTend,                 & ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo & ! (in)
       & )

    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: y_EulerHT(JA)
    real(DP), intent(out) :: y_BolusHT(JA)
    real(DP), intent(out) :: y_IsoDiffHT(JA)
    real(DP), intent(out) :: y_NumDiffTend(JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)

    ! 局所変数
    ! Local variables
    !
    
    real(DP) :: xyz_SLat(IA,JA,KA)
    real(DP) :: xyz_EulerHT(IA,JA,KA)
    real(DP) :: xyz_BolusHT(IA,JA,KA)
    real(DP) :: xyz_IsoDiffHT(IA,JA,KA)
    real(DP) :: xy_CosLat(IA,JA)
    real(DP) :: xyz_PTemp_numdiff(IA,JA,KA)
    real(DP) :: y_OcnHT(JA)
    real(DP) :: y_OcnHT_conv(JA)

    real(DP) :: xyz_PsiLonGM(IA,JA,KA)
    real(DP) :: yz_MSF_GM(JA,KA)
    integer :: k
    
    ! 実行文; Executable statement
    !
    
    xy_CosLat(:,:) = cos(xyz_Lat(:,:,KS))
    xyz_EulerHT(:,:,:) = RefDens*Cp0*xyz_V*xyz_PTemp
    
    call get_RediGM_HTFlx( &
         & xyz_IsoDiffHT, xyz_BolusHT, xyz_SLat, xyz_PTemp_numdiff, xyz_PsiLonGM, & ! (out)
         & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo )                             ! (in)
 
    y_EulerHT(:)   = IntLonZ(xyz_EulerHT)
    y_BolusHT(:)   = IntLonZ(xyz_BolusHT)
    y_IsoDiffHT(:) = IntLonZ(xyz_IsoDiffHT)

    y_OcnHT(:) = y_EulerHT + y_BolusHT + y_IsoDiffHT
    y_NumDiffTend(:) = IntLonZ(RefDens*Cp0*xyz_PTemp_numdiff)/(2d0*PI*RPlanet*xy_CosLat(IS,:))
    y_OcnHT_conv(:) = - DivLat(y_OcnHT) + y_NumDiffTend

    call DOGCM_IO_History_HistPut( 'EulerHT', y_EulerHT(JS:JE)*Unit2PWFac )
    call DOGCM_IO_History_HistPut( 'BolusHT', y_BolusHT(JS:JE)*Unit2PWFac )
    call DOGCM_IO_History_HistPut( 'IsoDiffHT', y_IsoDiffHT(JS:JE)*Unit2PWFac )
    call DOGCM_IO_History_HistPut( 'OcnHT', y_OcnHT(JS:JE)*Unit2PWFac )
    call DOGCM_IO_History_HistPut( 'OcnHT_conv', y_OcnHT_conv(JS:JE) )
    call DOGCM_IO_History_HistPut( 'NumDiffTend', y_NumDiffTend(JS:JE) )

    do k=KS, KE
       yz_MSF_GM(:,k) = sum( xyz_PsiLonGM(IS:IE,:,k) * RefDens*2d0*PI*RPlanet*xy_CosLat(IS:IE,:), 1)
    end do
    call DOGCM_IO_History_HistPut( 'MSF_GM', yz_MSF_GM(JS:JE,KS:KE)/1d9 )
    
  contains
    function DivLat(y_flxLatIntLon) result(y_div)

      use SpmlUtil_mod, only: &
           & w_DivMu_xy, xy_w
      
      real(DP), intent(in) :: y_flxLatIntLon(JA)
      real(DP) :: y_div(JA)

      real(DP) :: xy_div(IA,JA)
      real(DP) :: xy_flxLatCos(IA,JA)
      integer :: j
      
      do j = JS,JE
         xy_flxLatCos(:,j) = y_flxLatIntLon(j)/(2d0*PI*RPlanet)
      end do
      xy_div(IS:IE,JS:JE) = xy_w(w_DivMu_xy(xy_flxLatCos(IS:IE,JS:JE)))/RPlanet
      y_div(:) = xy_div(IS,:)
      
    end function DivLat

    function IntLonZ(xyz) result(y)
      real(DP), intent(in) :: xyz(IA,JA,KA)
      real(DP) :: y(JA)

      real(DP) :: xy_Tmp(IA,JA)
      
      xy_Tmp(:,:) = Calculus_IntBtmToTop( xyz, xyz_H )
      y(:)   = sum( xy_Tmp(IS:IE,:) * 2d0*PI*RPlanet*xy_CosLat(IS:IE,:), 1)
    end function IntLonZ
    
  end subroutine OcnDiag_HTDiagnose

  subroutine get_RediGM_HTFlx( &
       & xyz_Redi_HTFlx, xyz_GM_HTFlx, xyz_SLat, xyz_PTemp_numdiff,  &  ! (out)
       & xyz_PsiLonGM,                                               &  ! (out)
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo                  &  ! (in)
       & )

    use EOSDriver_mod, only: &
         & EOSDriver_Eval
    use LPhys_RediGM_hspm_vfvm_mod, only: &
         & calc_GradTRC,            &
         & calc_IsoNeutralSlope,    &
         & calc_IsopycDiffFlux,     &
         & calc_SkewFlux,           &
         & LPhys_RediGM_hspm_vfvm_GetParameters
    
    use LPhys_RediGMHelper_mod, only: &
         & prepare_SlopeTapering

    use DOGCM_Admin_TInteg_mod, only: &
         & DelTime
    use LPhys_DIFF_spm_mod, only: &
         & LPhys_DIFF_spm_LMixTRCRHSImpl
    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyz_Redi_HTFlx(IA,JA,KA)
    real(DP), intent(out) :: xyz_GM_HTFlx(IA,JA,KA)
    real(DP), intent(inout) :: xyz_SLat(IA,JA,KA)
    real(DP), intent(out) :: xyz_PTemp_numdiff(IA,JA,KA)
    real(DP), intent(out) :: xyz_PsiLonGM(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: xyz_DensPot(IA,JA,KA)
    real(DP) :: xy_RefPress(IA,JA)
    
    ! Slopes for tracer iso-neutral mixing
    real(DP) :: xyz_SLon(IA,JA,KA)
    real(DP) :: xyr_SLon(IA,JA,KA)
    real(DP) :: xyr_SLat(IA,JA,KA)

    ! Coeffecient for the slope  tapering 
    real(DP) :: xyz_T(IA,JA,KA)

    !
    real(DP) :: xyzaa_HGradTRC(IA,JA,KA,2,2)
    real(DP) :: xyra_DSigTRC(IA,JA,KA,2)

    real(DP) :: xyz_FLonRedi(IA,JA,KA)
    real(DP) :: xyz_FLatRedi(IA,JA,KA)
    real(DP) :: xyr_FSigRedi(IA,JA,KA)

    real(DP) :: xyz_FLonGM(IA,JA,KA)
    real(DP) :: xyz_FLatGM(IA,JA,KA)
    real(DP) :: xyr_FSigGM(IA,JA,KA)
    
    integer :: k
    real(DP) :: xyz_TRC(IA,JA,KA)    
    integer :: TRCID

    real(DP) :: xyza_TRC(IA,JA,KA,2)
    real(DP) :: xyza_TRC_RHS_numdiff(IA,JA,KA,2)

    real(DP) :: KappaGM
    
    ! 実行文; Executable statement
    !

    xyza_TRC(:,:,:,1) = xyz_PTemp
    xyza_TRC(:,:,:,2) = xyz_Salt
    xyza_TRC_RHS_numdiff(:,:,:,:) = 0d0
    call LPhys_DIFF_spm_LMixTRCRHSImpl( xyza_TRC_RHS_numdiff,        & ! (inout)
         & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef, DelTime       &  ! (in)
         & )
    xyz_PTemp_numdiff = xyza_TRC_RHS_numdiff(:,:,:,1)
    
    !$omp parallel do private(xy_RefPress)
    do k=KS,KE
       xy_RefPress(:,:) = 0d0
       call EOSDriver_Eval( xyz_DensPot(IS:IE,JS:JE,k),                                    & ! (out)
            & xyz_PTemp(IS:IE,JS:JE,k), xyz_Salt(IS:IE,JS:JE,k), xy_RefPress(IS:IE,JS:JE) )  ! (in)
    end do

    ! Calculate the components of the isoneutral slope. 
    !

    call calc_GradTRC( &
       & xyzaa_HGradTRC, xyra_DSigTRC,            & ! (out)
       & xyz_PTemp, xyz_Salt                      & ! (in)
       & )    

    call calc_IsoNeutralSlope( &
         & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat,            & ! (out)
         & xyz_PTemp, xyz_Salt, xyzaa_HGradTRC, xyra_DSigTRC, & ! (in)
         & xyz_H, xyz_Z )                                       !(in)

    call prepare_SlopeTapering( xyz_T, xyz_SLon, xyz_SLat, & ! (inout)
         & xyz_DensPot, xyz_Z                              & ! (in)
         & )
       
    TRCID          = 1
    xyz_TRC(:,:,:) = xyz_PTemp(:,:,:)
    
    call calc_IsopycDiffFlux( &
         & xyz_FLonRedi, xyz_FLatRedi, xyr_FSigRedi,                     &  ! (out)      
         & xyz_TRC, xyz_H, xyz_Z,                                        &  ! (in)
         & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat, xyz_T,                &  ! (in)
         & xyzaa_HGradTRC(:,:,:,:,TRCID), xyra_DSigTRC(:,:,:,TRCID)      &  ! (in)
         & )

    call calc_SkewFlux( &
         & xyz_FLonGM, xyz_FLatGM, xyr_FSigGM,                           &  ! (out)      
         & xyz_TRC, xyz_H, xyz_Z,                                        &  ! (in)
         & xyz_SLon, xyz_SLat, xyr_SLon, xyr_SLat, xyz_T,                &  ! (in)
         & xyzaa_HGradTRC(:,:,:,:,TRCID), xyra_DSigTRC(:,:,:,TRCID)      &  ! (in)
         & )

    xyz_Redi_HTFlx(:,:,:) = - RefDens*Cp0*xyz_FLatRedi
    xyz_GM_HTFlx(:,:,:)   = - RefDens*Cp0*xyz_FLatGM

    !
    call LPhys_RediGM_hspm_vfvm_GetParameters( &
         & KappaGM=KappaGM )
    xyz_PsiLonGM(:,:,:) = - KappaGM*xyz_T*xyz_SLat
    
  end subroutine get_RediGM_HTFlx
  
end module OcnDiag_HeatTransport_mod
