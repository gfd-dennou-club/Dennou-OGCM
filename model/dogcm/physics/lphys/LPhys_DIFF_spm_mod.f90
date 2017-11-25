!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module LPhys_DIFF_spm_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax, nMax

  use SpmlUtil_mod, only: &
       & w_xy, xy_w,           &
       & w_Lapla_w,            &
       & calc_VorDiv2UV,       &
       & calc_UVCosLat2VorDiv, &
       & xy_CosLat, rn, nm_l
  
  use DOGCM_Admin_Constants_mod, only: &
       & PI, RPlanet,               &
       & hViscCoef, hHyperViscCoef, &
       & hDiffCoef, hHyperDiffCoef

  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: LPhys_DIFF_spm_Init, LPhys_DIFF_spm_Final
  public :: LPhys_DIFF_spm_PrintParam

  public :: LPhys_DIFF_spm_LMixMomRHS
  public :: LPhys_DIFF_spm_LMixMomRHSImpl
  
  public :: LPhys_DIFF_spm_LMixTRCRHS
  public :: LPhys_DIFF_spm_LMixTRCRHSImpl

  public :: LPhys_DIFF_AdaptiveFilter4SIce
  
  ! 公開変数
  ! Publicvariable
  !  
  real(DP), public :: ViscCoefH
  real(DP), public :: DiffCoefH

  real(DP), public :: NumViscCoefH  
  real(DP), public :: NumDiffCoefH
  integer, public :: NumDiffOrdH
  integer, public :: NumDiffCutWaveNum

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !

  real(DP), allocatable :: w_HDiffCoefH(:)
  real(DP), allocatable :: w_HViscCoefH(:)
  real(DP), allocatable :: w_Filter(:)
  
  character(*), parameter:: module_name = 'LPhys_DIFF_spm_mod' !< Module Name
  logical :: isInitialzed = .false.

  public :: w_Filter, w_HDiffCoefH, w_HViscCoefH
  
contains

  !>
  !!
  !!
  subroutine LPhys_DIFF_spm_Init( &
       & ViscCoef, DiffCoef, NumDiffCoef,  & ! (in)
       & configNmlName                     & ! (in)
       & )       

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in), optional     :: ViscCoef
    real(DP), intent(in), optional     :: DiffCoef
    real(DP), intent(in), optional     :: NumDiffCoef
    character(*), intent(in), optional :: configNmlName   !< Namelist name

    ! 局所変数
    ! Local variables
    !    

    integer :: l
    real(DP) :: w_LaplaEigVal(lMax)
    integer :: nm(2)
    integer :: Nc
    
    ! 実行文; Executable statements
    !

    ! Set default values
    ViscCoefH    = 0d0
    DiffCoefH    = 0d0
    NumDiffCoefH = 0d0

    ! Set some parameters from arguments
    if(present(ViscCoef)) ViscCoefH = ViscCoef
    if(present(DiffCoef)) DiffCoefH = DiffCoef
    if(present(NumDiffCoef)) NumDiffCoefH = NumDiffCoef

    ! If configNmlFileName is specfied, we read namelist file to set the values of parameters.
    if(present(configNmlName)) then
       call read_nmlData(configNmlName)
    end if

    call LPhys_DIFF_spm_PrintParam()

    !----------------------------------------

    allocate(w_HViscCoefH(lMax))    
    allocate(w_HDiffCoefH(lMax))
    allocate(w_Filter(lMax))
    
    w_LaplaEigVal(:) = rn(:,1)/RPlanet**2

!!$    Nc = int( 0.5d0 * nMax )
!!$    write(*,*) "Nc=",  Nc
    
    !    !$omp parallel do private(nm)
    do l = 1, lMax
       w_HDiffCoefH(l) = -(    DiffCoefH*(-w_LaplaEigVal(l))                              &
            &                + NumDiffCoefH*(-w_LaplaEigVal(l))**(NumDiffOrdH/2) )

       w_HViscCoefH(l) = -(     ViscCoefH*(-w_LaplaEigVal(l) - 2d0/RPlanet**2)            &
            &                +  NumDiffCoefH*(  (-w_LaplaEigVal(l))**(NumDiffOrdH/2)      &
            &                +                - (2d0/RPlanet**2)**(NumDiffOrdH/2) )        &
            &             )

!!$       nm(:) = nm_l(l)
!!$       write(*,*) nm(1), ":", (nm(1) - Nc)/dble(nMax - Nc)
!!$       if (nm(1) < Nc) then
!!$          w_Filter(l) = 1d0
!!$       else
!!$          w_Filter(l) = exp(- 0.18d0 * ((nm(1) - Nc)/dble(nMax - Nc))**4)
!!$       end if
       w_Filter(l) = 1d0
    end do
    
    write(*,*) "Filter information (super spectral viscosity & dealiasing)"
    write(*,*) "total wavenumber | SSVFilter(MOM) | SSVFilter(TRC) | Filter for dealiasing"
    do l=1, lMax
       nm(:) = nm_l(l)
       write(*,'(i,3(ES15.7))') nm(1),  &
            & 1d0/(1d0 - 43200d0*w_HViscCoefH(l)), 1d0/(1d0 - 43200d0*w_HDiffCoefH(l)), &
            & w_Filter(l)
    end do
    !----------------------------------------
    
    isInitialzed = .true.
    
  end subroutine LPhys_DIFF_spm_Init

  !>
  !!
  !!
  subroutine LPhys_DIFF_spm_Final()

    ! 実行文; Executable statements
    !

    if( isInitialzed ) then
       deallocate( w_HDiffCoefH, w_HViscCoefH, w_Filter )
    end if
    
  end subroutine LPhys_DIFF_spm_Final


  !> @brief 
  !!
  !! @return 
  !!
  subroutine LPhys_DIFF_spm_PrintParam()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 実行文; Executable statement
    !
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, ' ViscCoefH      = %f [m2/s]', d=(/ ViscCoefH    /) )
    call MessageNotify( 'M', module_name, ' DiffCoefH      = %f [m2/s]', d=(/ DiffCoefH    /) )
    call MessageNotify( 'M', module_name, ' NumDiffOrdH    = %d [1]',    i=(/ NumDiffOrdH  /) )
    call MessageNotify( 'M', module_name, ' NumDiffCoefH   = %f [m%d/s]', d=(/ NumDiffCoefH /), i=(/ NumDiffOrdH /) )

    
  end subroutine LPhys_DIFF_spm_PrintParam

  !--------------------------------------------------------------------------------------------------

  subroutine LPhys_DIFF_spm_LMixMOMRHS(            &
       & xyz_U_RHS, xyz_V_RHS,                                   & ! (inout)
       & xyz_U, xyz_V, xyz_H,                                    & ! (in)
       & hViscCoef, hHyperViscCoef                               & ! (in)
       & )

    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: hViscCoef
    real(DP), intent(in) :: hHyperViscCoef


    ! 局所変数
    ! Local variables
    !    
    integer :: n
    integer :: k
    
    real(DP) :: w_Vor_RHS(lMax)
    real(DP) :: w_Div_RHS(lMax)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    real(DP) :: xy_U_HDiff(IA,JA)
    real(DP) :: xy_V_HDiff(IA,JA)

    ! 実行文; Executable statements
    !

    !$omp parallel do private(w_Vor, w_Div, w_Vor_RHS, w_Div_RHS, xy_U_HDiff, xy_V_HDiff)
    do k=KS, KE
       call calc_UVCosLat2VorDiv( &
            & xyz_U(IS:IE,JS:JE,k)*xy_CosLat, xyz_V(IS:IE,JS:JE,k)*xy_CosLat, &
            & w_Vor, w_Div                                                    &
            & )

       w_Vor_RHS(:) =  &
            &        hViscCoef* 2d0*w_Vor/RPlanet**2                          &
            &      + w_Lapla_w( &
            &           (hViscCoef - 2d0*hHyperViscCoef/RPlanet**2)*w_Vor     &
            &         - hHyperViscCoef*w_Lapla_w(w_Vor)/RPlanet**2            &
            &        )/RPlanet**2

       w_Div_RHS(:) =  &
            &        hViscCoef* 2d0*w_Div/RPlanet**2                          &
            &      + w_Lapla_w( &
            &           2d0*(hViscCoef - 2d0*hHyperViscCoef/RPlanet**2)*w_Div &
            &         - 2d0*hHyperViscCoef*w_Lapla_w(w_Div)/RPlanet**2        &
            &        )/RPlanet**2


       call calc_VorDiv2UV( w_Vor_RHS, w_Div_RHS,                             &
            & xy_U_HDiff(IS:IE,JS:JE), xy_V_HDiff(IS:IE,JS:JE) )

       xyz_U_RHS(:,:,k) = xyz_U_RHS(:,:,k) + xy_U_HDiff
       xyz_V_RHS(:,:,k) = xyz_V_RHS(:,:,k) + xy_V_HDiff
    end do
    
  end subroutine LPhys_DIFF_spm_LMixMOMRHS

  subroutine LPhys_DIFF_spm_LMixMOMRHSImpl(            &
       & xyz_U_RHS, xyz_V_RHS,                            & ! (inout)
       & xyz_U, xyz_V, xyz_H, hViscCoef, hHyperViscCoef,  & ! (in)
       & dt                                               & ! (in)
       & )

    ! 宣言文; Declaration statement
    !          
    real(DP), intent(inout) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: hViscCoef
    real(DP), intent(in) :: hHyperViscCoef
    real(DP), intent(in) :: dt

    ! 局所変数
    ! Local variables
    !    
    integer :: n
    integer :: k

    real(DP) :: w_Fact(lMax)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    real(DP) :: xy_U_HDiff(IA,JA)
    real(DP) :: xy_V_HDiff(IA,JA)

    ! 実行文; Executable statements
    !
    !$omp parallel do private(w_Vor, w_Div, w_Fact, xy_U_HDiff, xy_V_HDiff)
    do k=KS, KE
       call calc_UVCosLat2VorDiv( &
            & xyz_U(IS:IE,JS:JE,k)*xy_CosLat, xyz_V(IS:IE,JS:JE,k)*xy_CosLat, &
            & w_Vor, w_Div                                                    &
            & )

       w_Fact(:) = 1d0/(1d0/w_HViscCoefH(:) - dt)
       call calc_VorDiv2UV( w_Fact*w_Vor, w_Fact*w_Div,                       & ! (in)
            & xy_U_HDiff(IS:IE,JS:JE), xy_V_HDiff(IS:IE,JS:JE)                & ! (out)
            & )

       xyz_U_RHS(:,:,k) = xyz_U_RHS(:,:,k) + xy_U_HDiff(:,:)
       xyz_V_RHS(:,:,k) = xyz_V_RHS(:,:,k) + xy_V_HDiff(:,:)
    end do
    
  end subroutine LPhys_DIFF_spm_LMixMOMRHSImpl
  
  subroutine LPhys_DIFF_spm_LMixTRCRHS(            &
       & xyza_TRC_RHS,                              & ! (inout)
       & xyza_TRC, xyz_H,                           & ! (in)
       & hDiffCoef, hHyperDiffCoef                  & ! (in)
       & )

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(inout) :: xyza_TRC_RHS(0:iMax-1,jMax,KS:KE,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,KS:KE,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KS:KE)
    real(DP), intent(in) :: hDiffCoef
    real(DP), intent(in) :: hHyperDiffCoef

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: w_TRC(lMax)

    integer :: n
    integer :: k
    

    ! 実行文; Executable statements
    !
    
    do n = 1, TRC_TOT_NUM
       !$omp parallel do private(w_TRC)
       do k=KS, KE
          w_TRC(:) = w_xy(xyza_TRC(:,:,k,n))
          
          xyza_TRC_RHS(:,:,k,n) = xyza_TRC_RHS(:,:,k,n) &
               & + xy_w( w_Lapla_w( &
               &     hDiffCoef*w_TRC                             &
               &   - hHyperDiffCoef*w_Lapla_w(w_TRC)/RPlanet**2  &
               & )/RPlanet**2 )
       end do
    end do
    
  end subroutine LPhys_DIFF_spm_LMixTRCRHS

  subroutine LPhys_DIFF_spm_LMixTRCRHSImpl(      &
       & xyza_TRC_RHS,                                 & ! (inout)
       & xyza_TRC, xyz_H,                              & ! (in)
       & hDiffCoef, hHyperDiffCoef,                    & ! (in)
       & dt                                            & ! (in)
       & )

    ! 宣言文; Declaration statement
    !      
    real(DP), intent(inout) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: hDiffCoef
    real(DP), intent(in) :: hHyperDiffCoef
    real(DP), intent(in) :: dt

    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: w_TRC(lMax)
    
    integer :: l
    integer :: n
    integer :: k
    
    
    ! 実行文; Executable statements
    !
    
    do n = 1, TRC_TOT_NUM
       !$omp parallel do private(w_TRC)
       do k=KS, KE
          w_TRC(:) = w_xy(xyza_TRC(IS:IE,JS:JE,k,n))

          xyza_TRC_RHS(IS:IE,JS:JE,k,n) = xyza_TRC_RHS(IS:IE,JS:JE,k,n) +  &
               & (   xy_w( w_TRC(:)/(1d0/w_HDiffCoefH(:) - dt) ) &
               & )
       end do
    end do
    
  end subroutine LPhys_DIFF_spm_LMixTRCRHSImpl

  subroutine LPhys_DIFF_AdaptiveFilter4SIce( &
       & xyza_TRC, xy_OcnSfcCellMask )

    use SpmlUtil_mod
    use DOGCM_Admin_Grid_mod, only: &
         & xy_Lon, xy_Lat, &
         & x_IAXIS_Weight, y_JAXIS_Weight
    use DOGCM_Admin_Variable_mod, only: &
         & TRCID_PTEMP
    use DOGCM_Boundary_Vars_mod, only: &
         & OCNCELLMASK_OCEAN, OCNCELLMASK_SICE
    
    real(DP), intent(inout) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    integer, intent(in) :: xy_OcnSfcCellMask(IA,JA)

    integer :: i
    integer :: j
    integer :: j2
    integer :: k
    integer :: p
    integer :: l
    real(DP) :: w_Filter_Sfc(lMax)
    real(DP) :: xy_Mollifier(IA,JA)
    integer :: nm(2)
    real(DP) :: xi

    real(DP) :: xy_Q(IA,JA)
    real(DP) :: xy_QN(IA,JA)
    real(DP) :: w_Q(lMax)

    real(DP) :: a_SIceEdgeLat(10)
    integer :: NSIceEdge
    real(DP) :: dist
    integer :: xy_filter_order(IA,JA)
    real(DP) :: wy_Pn(lMax,JA)
    real(DP) :: mu
    real(DP) :: cp
    real(DP) :: sigma_c

    real(DP) :: xi_c
    
    xy_QN(:,:) = xyza_TRC(:,:,KS,TRCID_PTEMP)

    p = 4d0 !max( 2d0, (sqrt(nMax*dist/(0.2d0*PI))) )
    cp = 32d0!2**p * 0.75d0 * (9d0*p**2 + 3d0*p + 14d0)/(9d0*p**2 + 12d0*p + 4d0)
    xi_c = 2d0/3d0
    do l=1, lMax
       nm(:) = nm_l(l)
       xi = nm(1)/dble(nMax)
       if (xi < xi_c ) then
          w_Filter_Sfc(l) = 1d0
       else if (xi >= xi_c) then
          w_Filter_Sfc(l) = exp(-cp*((xi - xi_c)**p)) !exp(cp*xi**p/(xi**2 - 1d0))
       end if
       sigma_c = 0.5d0*(1d0 + cos(PI*xi))
       !sigma_c**4*(35d0 - 84d0*sigma_c + 70d0*sigma_c**2 - 20d0*sigma_c**3)
    end do

    write(*,*) w_Filter_Sfc
    xyza_TRC(IS:IE,JS:JE,KS,TRCID_PTEMP) = xy_w(w_xy(xy_QN(IS:IE,JS:JE))*w_Filter_Sfc)
    return
    
    NSIceEdge = 0
    i = IS
    do j = JS, JE-1
       if(     (     (xy_OcnSfcCellMask(i,j  ) == OCNCELLMASK_OCEAN)   &
         &     .and. (xy_OcnSfcCellMask(i,j+1) == OCNCELLMASK_SICE ) ) &
         & .or.(     (xy_OcnSfcCellMask(i,j  ) == OCNCELLMASK_SICE)   &
         &     .and. (xy_OcnSfcCellMask(i,j+1) == OCNCELLMASK_OCEAN ) ) ) then
          NSIceEdge = NSIceEdge + 1
          if (NSIceEdge > size(a_SIceEdgeLat)) then
             write(*,*) "The number of iceline exceed expectations. Check!"
             stop
          end if
          a_SIceEdgeLat(NSIceEdge) = 0.5d0*(xy_Lat(i,j) + xy_Lat(i,j+1))
       end if
    end do
    do j = JS, JE
       mu = sin(xy_Lat(i,j))
       wy_Pn(1,j) = 1d0
       wy_Pn(2,j) = sqrt(3d0)*mu
       do l=3, lMax
          nm(:) = nm_l(l)          
          wy_Pn(l,j) = &
               &   sqrt((2*nm(1) - 1d0)*(2*nm(1) + 1d0))/dble(nm(1)) * mu * wy_Pn(l-1,j)      &
               & - (nm(1) - 1)*sqrt(2*nm(1) + 1d0)/(nm(1)*sqrt(2*nm(1) - 3d0)) * wy_Pn(l-2,j)          
       end do
    end do
    
    do j = JS, JE
    do i = IS, IE
       xy_Q(:,:) = 0d0; xy_Q(i,j) = 1d0
       w_Q(:) = w_xy(xy_Q(IS:IE,JS:JE))

       dist = minval( abs(xy_Lat(i,j) - a_SIceEdgeLat(1:NSIceEdge)) )
       p = 2d0 !max( 2d0, (sqrt(nMax*dist/(0.2d0*PI))) )
       cp = 1d0!2**p * 0.75d0 * (9d0*p**2 + 3d0*p + 14d0)/(9d0*p**2 + 12d0*p + 4d0)
       xy_filter_order(i,j) = p

       do l=1, lMax
          nm(:) = nm_l(l)

          xi = nm(1)/dble(nMax)
          if (xi < 1d0) then
             w_Filter(l) = exp(cp*xi**p/(xi**2 - 1d0))
          else
             w_Filter(l) = 0d0
          end if
       end do

       do j2=JS, JE
          xy_Mollifier(IS,j2) = sum(w_Filter*w_Q*wy_Pn(:,j2)) / (x_IAXIS_Weight(i)*y_JAXIS_Weight(j))
       end do
       xy_QN(i,j) = IntLonLat_xy( xy_Mollifier(IS:IE,JS:JE) * xyza_TRC(IS:IE,JS:JE,KS,TRCID_PTEMP) )
    end do
    end do

    
!!$    write(*,*) "SIceLat", a_SIceEdgeLat(1:NSIceEdge)
!!$    write(*,*) "FilterOrfer:", xy_filter_order(IS,JS:JE)
!!$    write(*,*) "Before:", xyza_TRC(IS,JS:JE,KS,TRCID_PTEMP)
!!$    write(*,*) "After:", xy_QN(IS,JS:JE)
!!$    write(*,*) IntLonLat_xy(xy_QN(IS:IE,JS:JE)), IntLonLat_xy(xyza_TRC(IS:IE,JS:JE,KS,TRCID_PTEMP))
!!$    stop
    xyza_TRC(IS:IE,JS:JE,KS,TRCID_PTEMP) = xy_QN(IS:IE,JS:JE)
    
  end subroutine LPhys_DIFF_AdaptiveFilter4SIce
  
  !--------------------------------------------------------------------------------------------------

  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !
   
    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    use dc_calendar, only: DCCalConvertByUnit

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !

    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 

    character(STRING) :: msg
    
    ! IOSTAT of NAMELIST read

    real(DP) :: NumDiffTimeVal
    real(DP) :: NumDiffTimeValSec
    character(TOKEN) :: NumDiffTimeUnit
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /LPhys_DIFF_nml/ &
         & ViscCoefH,        &
         & DiffCoefH,        &
         & NumDiffOrdH,      &
         & NumDiffCoefH,     &
         & NumDiffTimeVal,   &
         & NumDiffTimeUnit,  &
         & NumDiffCutWaveNum

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    ViscCoefH       = hViscCoef
    DiffCoefH       = hDiffCoef
    NumDiffOrdH     = 4
    NumDiffCoefH    = -1d0
    NumDiffTimeVal  = 1d100
    NumDiffTimeUnit = 'day'
    NumDiffCutWaveNum = 0
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                        ! (in)
            & nml = LPhys_DIFF_nml,          &  ! (out)
            & iostat = iostat_nml, iomsg=msg )  ! (out)
       close( unit_nml )
    end if

    ! Calculate the coefficient of numerical diffusion.

    if ( NumDiffCoefH == -1d0 ) then
       NumDiffTimeValSec = DCCalConvertByUnit(NumDiffTimeVal, NumDiffTimeUnit, 'sec')
       NumDiffCoefH = (nMax*(nMax + 1)/RPlanet**2 )**(-NumDiffOrdH/2) / NumDiffTimeValSec
    end if
    
  end subroutine read_nmlData

  
end module LPhys_DIFF_spm_mod

