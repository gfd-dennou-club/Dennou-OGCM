!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSEddyMixing_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use gtool_history

  !
  use Constants_mod, only: &
       & RPlanet

  use SpmlUtil_mod

  use DiagnoseUtil_mod

  use EOSDriver_mod

  use GridSet_mod

  

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SGSEddyMixing_Init, SGSEddyMixing_Final
  public :: SGSEddyMixing_PrepareOutput, SGSEddyMixing_Output
  public :: SGSEddyMixing_AddMixingTerm

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSEddyMixing_mod' !< Module Name
  integer, parameter :: SGSEddyMixing_Redi = 1
  integer, parameter :: SGSEddyMixing_GM90 = 2

  real(DP) :: SGSEddyMixType
  real(DP) :: Kapper_rho  !< Isopycnal diffusivity
  real(DP) :: SlopeMaxVal

  type(gt_history) :: hst_SGSEddyMix
  logical :: OutputFlag
  
contains

  !>
  !!
  !!
  subroutine SGSEddyMixing_Init(SGSEddyMixParamType, KapperRho)

    ! 実行文; Executable statements
    !
    integer, intent(in) :: SGSEddyMixParamType
    real(DP), intent(in) :: KapperRho

    SGSEddyMixType = SGSEddyMixParamType
    Kapper_rho = 1000d0!KapperRho
    SlopeMaxVal = 2.5d-3
    OutputFlag = .false.

  end subroutine SGSEddyMixing_Init

  !>
  !!
  !!
  subroutine SGSEddyMixing_Final()

    ! 実行文; Executable statements
    !
    if(OutputFlag) call HistoryClose(hst_SGSEddyMix)

  end subroutine SGSEddyMixing_Final

  !> @brief 
  !!
  !!
  subroutine SGSEddyMixing_AddMixingTerm(wz_PTempRHS, wz_SaltRHS, &
       & wz_PTemp, wz_Salt, xy_totDepth )
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(lMax,0:kMax), intent(inout) :: wz_PTempRHS, wz_SaltRHS
    real(DP), dimension(lMax,0:kMax), intent(in)  :: wz_PTemp, wz_Salt
    real(DP), dimension(0:iMax-1,jMax), intent(in) :: xy_totDepth

    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_DensPot, xyz_RefPress, xyz_CosLat, xyz_totDepth

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
         & xyz_SLon, &  !< The longitude componet of th slope of isopycnal surface
         & xyz_SLat     !< The latitude componet of th slope of isopycnal surface

    real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_PTemp, xyz_Salt

    integer :: k

    ! 実行文; Executable statement
    !
    
    xyz_CosLat = cos(xyz_Lat)
    forAll(k=0:kMax) xyz_totDepth(:,:,k) = xy_totDepth

    xyz_PTemp = xyz_wz(wz_PTemp)
    xyz_Salt = xyz_wz(wz_Salt)

    ! Calculate the potential density. 
    xyz_RefPress = 0d0
    call EOSDriver_Eval(xyz_DensPot, &  !(out)
         & xyz_PTemp, xyz_Salt, xyz_RefPress ) !(in)

    ! Calculate the components of the isoneutral slope. 
    call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
         & xyz_DensPot, xyz_totDepth ) !(in)

    ! 
    
    wz_PTempRHS = wz_PTempRHS + calc_IsoNeutralDiffTerm(wz_PTemp, xyz_PTemp)
    wz_SaltRHS = wz_SaltRHS + calc_IsoNeutralDiffTerm(wz_Salt, xyz_Salt)

    if(SGSEddyMixType == SGSEddyMixing_GM90) then
       wz_PTempRHS = wz_PTempRHS + calc_EddyInducedVelAdvTerm(wz_PTemp, xyz_PTemp)
       wz_SaltRHS = wz_SaltRHS + calc_EddyInducedVelAdvTerm(wz_Salt, xyz_Salt)
    end if

!!$    write(*,*) "MeanSlope", AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(sqrt(xyz_SLon**2+xyz_Slat**2)))
!!$    write(*,*) "Max Slope", maxval(sqrt(xyz_SLon**2+xyz_Slat**2)), maxloc(sqrt(xyz_SLon**2+xyz_Slat**2))

    contains


      function calc_IsoNeutralDiffTerm(wz_Tracer, xyz_Tracer) result(wz_DiffTerm)
        real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
        real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
        real(DP) :: wz_DiffTerm(lMax,0:kMax)

        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_GradLonT, xyz_GradLatT, xyz_DzT
        
        xyz_GradLonT = xya_GradLon_wa(wz_Tracer)/RPlanet
        xyz_GradLatT = xya_GradLat_wa(wz_Tracer)/RPlanet
        xyz_DzT = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))/xyz_totDepth

        wz_DiffTerm = &
             &    wz_AlphaOptr_xyz( &
             &       Kapper_rho*(xyz_GradLonT + xyz_SLon*xyz_DzT)*xyz_CosLat, &
             &       Kapper_rho*(xyz_GradLatT + xyz_SLat*xyz_DzT)*xyz_CosLat  &
             &    ) &
             & +  wz_wt(wt_DSig_wt(wt_xyz( &
             &       Kapper_rho/xyz_totDepth*(xyz_SLon*xyz_GradLonT + xyz_SLat*xyz_GradLatT + &
             &           (xyz_SLon**2+xyz_SLat**2)*xyz_DzT) &
             &    )))
        
      end function calc_IsoNeutralDiffTerm
    
      function calc_EddyInducedVelAdvTerm(wz_Tracer, xyz_Tracer) result(wz_AdvTerm)
        real(DP), intent(in) :: wz_Tracer(lMax,0:kMax)
        real(DP), intent(in) :: xyz_Tracer(0:iMax-1,jMax,0:kMax)
        real(DP) :: wz_AdvTerm(lMax,0:kMax)

        real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
             & xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW


!!$        call calc_BolusVelocity(xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, & ! (out)
!!$             & xyz_SLon, xyz_SLat, xyz_totDepth )                                    ! (in)

!!$write(*,*) "U*,V*,W*", &
!!$     & AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_EddInducedU)), &
!!$     & AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_EddInducedV)), &
!!$     & AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyz_EddInducedW))
!!$write(*,*) "U*,V*,W*", &
!!$     & maxval(xyz_EddInducedU), maxval(xyz_EddInducedV), maxval(xyz_EddInducedW)
!!$write(*,*) "U*,V*,W*", &
!!$     & maxloc(xyz_EddInducedU), maxloc(xyz_EddInducedV), maxloc(xyz_EddInducedW)


        wz_AdvTerm = &
             !             & +  wz_AlphaOptr_xyz(xyz_EddInducedU*xyz_Tracer*xyz_CosLat, xyz_EddInducedV*xyz_Tracer*xyz_CosLat) &
             !             & +  wz_wt(wt_DSig_wt(wt_xyz(xyz_EddInducedW*xyz_Tracer/xyz_totDepth)))
             & - wz_AlphaOptr_xyz( &
             &   Kapper_rho*xyz_SLon*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))*xyz_CosLat/xyz_totDepth, &
             &   Kapper_rho*xyz_SLat*xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_Tracer)))*xyz_CosLat/xyz_totDepth  &
             & ) &
             & + wz_wt(wt_DSig_wt(wt_xyz( &
             & Kapper_rho*(xyz_SLon*xya_GradLon_wa(wz_Tracer)/RPlanet + xyz_SLat*xya_GradLat_wa(wz_Tracer)/RPlanet) &
             & /xyz_totDepth )))
      end function calc_EddyInducedVelAdvTerm

    end subroutine SGSEddyMixing_AddMixingTerm


    subroutine calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &
         & xyz_DensPot, xyz_totDepth )
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_SLon, xyz_SLat
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_DensPot, xyz_totDepth

      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: &
           & xyz_GradDensPotLon, xyz_GradDensPotLat, &
           & xyz_DzDensPot, xyz_SLimit

      integer :: i, j, k

      xyz_GradDensPotLon = xya_GradLon_wa(wz_xyz(xyz_DensPot))/RPlanet
      xyz_GradDensPotLat = xya_GradLat_wa(wz_xyz(xyz_DensPot))/RPlanet

      xyz_DzDensPot = xyz_xyt(xyt_DSig_xyt(xyt_xyz(xyz_DensPot)))/xyz_totDepth
      ! Modify the slope with slope clipping.
      xyz_SLimit = -sqrt(xyz_GradDensPotLon**2+xyz_GradDensPotLat**2)/SlopeMaxVal
      where(xyz_DzDensPot >  xyz_SLimit)
         xyz_DzDensPot = xyz_SLimit
      end where

      xyz_SLon = - xyz_GradDensPotLon/xyz_DzDensPot
      xyz_SLat = - xyz_GradDensPotLat/xyz_DzDensPot
      xyz_SLon(:,:,0) = 0d0
      xyz_SLat(:,:,0) = 0d0
      xyz_SLon(:,:,kMax) = 0d0
      xyz_SLat(:,:,kMax) = 0d0

    end subroutine calc_IsoNeutralSlope

    subroutine calc_BolusVelocity(xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW, &
         & xyz_SLon, xyz_SLat, xyz_totDepth )

      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(out) :: xyz_EddInducedU, xyz_EddInducedV, xyz_EddInducedW
      real(DP), dimension(0:iMax-1,jMax,0:kMax), intent(in) :: xyz_SLon, xyz_SLat, xyz_totDepth

      real(DP), dimension(0:iMax-1,jMax,0:kMax) :: xyz_CosLat

      xyz_CosLat = cos(xyz_Lat)
      xyz_EddInducedU = - xyz_xyt(xyt_DSig_xyt(xyt_xyz(Kapper_rho*xyz_SLon)))/xyz_totDepth
      xyz_EddInducedV = - xyz_xyt(xyt_DSig_xyt(xyt_xyz(Kapper_rho*xyz_SLat)))/xyz_totDepth
      xyz_EddInducedW = xyz_wz( &
           & wz_AlphaOptr_xyz(Kapper_rho*xyz_SLon*xyz_CosLat, Kapper_rho*xyz_SLat*xyz_CosLat) &
           & )

    end subroutine calc_BolusVelocity

    !!!!!!!!!!!!!!!!!!!!!

    subroutine SGSEddyMixing_Output()
      
      use VariableSet_mod, only: &
           & z_PTempBasic, xyz_PTempEddN, xyz_SaltN, xy_SurfHeightN, &
           & xy_totDepthBasic

      real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
           & xyz_DensPot, xyz_PTemp, xyz_RefPress, &
           & xyz_SLon, xyz_SLat, xyz_totDepth

      real(DP), dimension(0:iMax-1,jMax,0:kMax)  :: & 
           & xyz_BolusU, xyz_BolusV, xyz_BolusW

      if(.not. OutputFlag) return 

      xyz_PTemp = xyz_PTempEddN + spread(spread(z_PTempBasic,1,jMax), 1, iMax)
      xyz_totDepth = spread(xy_totDepthBasic + xy_SurfHeightN, 3, kMax+1)

      ! Calculate the potential density. 
      xyz_RefPress = 0d0
      call EOSDriver_Eval(xyz_DensPot, &  !(out)
           & xyz_PTemp, xyz_SaltN, xyz_RefPress ) !(in)

      ! Calculate the components of the isoneutral slope. 
      call calc_IsoNeutralSlope(xyz_SLon, xyz_SLat, &  !(out)
           & xyz_DensPot, xyz_totDepth ) !(in)
      
      !
      call calc_BolusVelocity(xyz_BolusU, xyz_BolusV, xyz_BolusW, & !(out)
           & xyz_SLon, xyz_SLat, xyz_totDepth )

      !
      call HistoryPut('SlopeLon', xyz_SLon, hst_SGSEddyMix)
      call HistoryPut('SlopeLat', xyz_SLat, hst_SGSEddyMix)
      call HistoryPut('BolusU', xyz_BolusU, hst_SGSEddyMix)
      call HistoryPut('BolusV', xyz_BolusV, hst_SGSEddyMix)
      call HistoryPut('BolusW', xyz_BolusW, hst_SGSEddyMix)

    end subroutine SGSEddyMixing_Output

    subroutine SGSEddyMixing_PrepareOutput(OriginTime, EndTime, Intrv, FilePrefix)

      use Constants_mod, only: PI

      use GridSet_mod, only: &
           & x_Lon, y_Lat

      use SpmlUtil_mod, only: g_Sig

      !
      !
      real(DP), intent(in) :: OriginTime, EndTime, Intrv
      character(*), intent(in) :: FilePrefix

      !
      !
      character(TOKEN) :: lonName, latName, sigName, timeName
      character(TOKEN) :: dims_XYZT(4)

      !
      OutputFlag = .true.

      !
      lonName = 'lon'; latName='lat'; sigName='sig'; timeName='t'
      dims_XYZT = (/ lonName, latName, sigName, timeName /)

      call HistoryCreate( &                            ! ヒストリー作成
           & file=trim(FilePrefix) // 'SGSEddyMixOutput.nc', title='OGCM Output',             &
           & source='OGCM Output',    &
           & institution='GFD_Dennou Club OGCM project',    &
           & dims=(/lonName,latName,sigName,timeName/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
           & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
           & units=(/'degree_east ','degree_north','(1)         ', 'sec.        '/),  &
           & origin=real(OriginTime), interval=real(Intrv),  &        
           & history=hst_SGSEddyMix  )    

      call HistoryPut(lonName, x_Lon*180d0/PI, hst_SGSEddyMix)
      call HistoryAddAttr(lonName, 'topology', 'circular', hst_SGSEddyMix)
      call HistoryAddAttr(lonName, 'modulo', 360.0, hst_SGSEddyMix)
      call HistoryPut(latName, y_Lat*180d0/PI, hst_SGSEddyMix)
      call HistoryPut(sigName, g_Sig, hst_SGSEddyMix)

      call HistoryAddVariable( varname='BolusU', &
         & dims=dims_XYZT, longname='bolus velocity(longitude) ', units='m/s', &
         & history=hst_SGSEddyMix )

      call HistoryAddVariable( varname='BolusV', &
         & dims=dims_XYZT, longname='bolus velocity(meridional) ', units='m/s', &
         & history=hst_SGSEddyMix )

      call HistoryAddVariable( varname='BolusW', &
         & dims=dims_XYZT, longname='bolus velocity(vertical) ', units='m/s', &
         & history=hst_SGSEddyMix )

      call HistoryAddVariable( varname='SlopeLon', &
         & dims=dims_XYZT, longname='the longitude component of slope of isoneutral ', units='1', &
         & history=hst_SGSEddyMix )

      call HistoryAddVariable( varname='SlopeLat', &
         & dims=dims_XYZT, longname='the latitude component of slope of isoneutral ', units='1', &
         & history=hst_SGSEddyMix )

    end subroutine SGSEddyMixing_PrepareOutput


end module SGSEddyMixing_mod

