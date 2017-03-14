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
  public :: OcnDiag_HTDiagnose

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'OcnDiag_HeatTransport_mod' !< Module Name
  
contains

  subroutine OcnDiag_HeatTransport_Init()
  end subroutine OcnDiag_HeatTransport_Init

  subroutine OcnDiag_HeatTransport_Final()
  end subroutine OcnDiag_HeatTransport_Final

  subroutine OcnDiag_HTDiagnose( &
       & y_EulerHT, y_BolusHT, y_IsoDiffHT,                                & ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: y_EulerHT(JA)
    real(DP), intent(out) :: y_BolusHT(JA)
    real(DP), intent(out) :: y_IsoDiffHT(JA)
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

    ! 実行文; Executable statement
    !
    
    xy_CosLat(:,:) = cos(xyz_Lat(:,:,KS))
    xyz_EulerHT(:,:,:) = RefDens*Cp0*xyz_V*xyz_PTemp
    
    call get_RediGM_HTFlx( xyz_IsoDiffHT, xyz_BolusHT, xyz_SLat, &
         & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo )

    y_EulerHT(:)   = IntLonZ(xyz_EulerHT)
    y_BolusHT(:)   = IntLonZ(xyz_BolusHT)
    y_IsoDiffHT(:) = IntLonZ(xyz_IsoDiffHT)

  contains
    function IntLonZ(xyz) result(y)
      real(DP), intent(in) :: xyz(IA,JA,KA)
      real(DP) :: y(JA)

      real(DP) :: xy_Tmp(IA,JA)
      
      xy_Tmp(:,:) = Calculus_IntBtmToTop( xyz, xyz_H )
      y(:)   = sum( xy_Tmp(IS:IE,:) * 2d0*PI*RPlanet*xy_CosLat(IS:IE,:), 1)
    end function IntLonZ
  end subroutine OcnDiag_HTDiagnose

  subroutine get_RediGM_HTFlx( &
       & xyz_Redi_HTFlx, xyz_GM_HTFlx, xyz_SLat,               &  ! (out)
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z, xy_Topo            &  ! (in)
       & )

    use EOSDriver_mod, only: &
         & EOSDriver_Eval
    use LPhys_RediGM_hspm_vfvm_mod, only: &
         & calc_GradTRC,            &
         & calc_IsoNeutralSlope,    &
         & calc_IsopycDiffFlux,     &
         & calc_SkewFlux
    use LPhys_RediGMHelper_mod, only: &
         & prepare_SlopeTapering

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyz_Redi_HTFlx(IA,JA,KA)
    real(DP), intent(out) :: xyz_GM_HTFlx(IA,JA,KA)
    real(DP), intent(inout) :: xyz_SLat(IA,JA,KA)
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

    ! 実行文; Executable statement
    !
    
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
    
  end subroutine get_RediGM_HTFlx
  
end module OcnDiag_HeatTransport_mod
