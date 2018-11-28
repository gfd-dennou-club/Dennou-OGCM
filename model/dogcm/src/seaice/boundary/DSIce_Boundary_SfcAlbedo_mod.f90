!-------------------------------------------------------------
! Copyright (c) 2015-2017 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to set some variables (e.g., surface and bottom flux) in order to satisfy boundary conditions
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_Boundary_SfcAlbedo_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  use DOGCM_Admin_Constants_mod, only: &
       & UNDEFVAL,                     &
       & CpOcn => Cp0,                 &
       & DensSeaWater => RefDens,      &
       & AlbedoOcean, PI

  use UnitConversion_mod, only: &
       & degC2K, K2degC
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_Boundary_SfcAlbedo_Init, DSIce_Boundary_SfcAlbedo_Final

  public :: DSIce_Boundary_SfcAlbedo_Get
  public :: DSIce_Boundary_SfcAlbedo_Const
  public :: DSIce_Boundary_SfcAlbedo_Simple
  
  ! 公開変数
  ! Public variable
  !

  character(*), public, parameter :: SICEGOVERNEQ_SFCALBEDO_CONST_NAME  = "Const"
  integer, public, parameter :: SICEGOVERNEQ_SFCALBEDO_CONST            = 1

  character(*), public, parameter :: SICEGOVERNEQ_SFCALBEDO_SIMPLE_NAME  = "Simple"
  integer, public, parameter :: SICEGOVERNEQ_SFCALBEDO_SIMPLE            = 2

  character(*), public, parameter :: SICEGOVERNEQ_SFCALBEDO_I07SGS_NAME  = "I07SGS"
  integer, public, parameter :: SICEGOVERNEQ_SFCALBEDO_I07SGS            = 3
  
  integer, public, save :: SIceSfcAlbedoMethod                               !< The method of sea-ice surface albedo

  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DSIce_Boundary_SfcAlbedo_mod' !< Module Name
  
  real(DP) :: AlbedoIceMax
  real(DP) :: AlbedoIceMin
  real(DP) :: AlbedoSnowWarm
  real(DP) :: AlbedoSnowCold
  real(DP) :: AlbedoSnowOld
  real(DP) :: TLVTempAlbedoSnow     !  Threashold limit value of temperature between
                                    !  warm and cold conditions in calculating sea-ice snow surface albedo. 
  real(DP) :: DampTimeAlbedoSnow

  logical :: SFCALBEDO_I07SGS_FlagIceAreaFrac
  
contains

  !>
  !!
  !!
  Subroutine DSIce_Boundary_SfcAlbedo_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call read_nmlData(configNmlName)
    
  end subroutine DSIce_Boundary_SfcAlbedo_Init

  !>
  !!
  !!
  subroutine DSIce_Boundary_SfcAlbedo_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Boundary_SfcAlbedo_Final

  !-----------------------------------------

  subroutine DSIce_Boundary_SfcAlbedo_Get( &
       & xy_SIceAlbedoAI,                                                     & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_IceThick, xy_SIceSfcTemp, xy_SeaSfcTemp & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SIceAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)
    
    select case(SIceSfcAlbedoMethod)
    case( SICEGOVERNEQ_SFCALBEDO_CONST )
       call DSIce_Boundary_SfcAlbedo_Const( &
            & xy_SIceAlbedoAI,                             & ! (out)
            & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp     & ! (in)
            & ) 
    case( SICEGOVERNEQ_SFCALBEDO_SIMPLE )       
       call DSIce_Boundary_SfcAlbedo_Simple( &
            & xy_SIceAlbedoAI,                                           & ! (out)
            & xy_SIceCon, xy_SnowThick, xy_IceThick, xy_SIceSfcTemp,     & ! (in)
            & xy_SeaSfcTemp                                              & ! (in)
            & ) 
    case( SICEGOVERNEQ_SFCALBEDO_I07SGS )       
       call DSIce_Boundary_SfcAlbedo_I07SGSParam( &
            & xy_SIceAlbedoAI,                                          & ! (out)
            & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,  & ! (in)
            & SFCALBEDO_I07SGS_FlagIceAreaFrac                          & ! (in)
            & ) 
    case default
       call MessageNotify( 'E', module_name,                  &
            & 'Unexpected SIceSfcAlbedoMethod is specified.')
    end select
    
  end subroutine DSIce_Boundary_SfcAlbedo_Get

  !-----------------------------------------
  
  subroutine DSIce_Boundary_SfcAlbedo_Const( &
       & xy_SIceAlbedoAI,                         & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp & ! (in)
       & )

    use DSIce_Admin_Constants_mod, only: &
         & AlbedoOcean, AlbedoIce, AlbedoSnow, AlbedoMeltSnow, &
         & FreezeTempWater, IceMaskMin
    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SIceAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)

    ! 局所変数
    ! Local variables
    !        
    integer :: i
    integer :: j

    ! 実行文; Executable statements
    !

    !$omp parallel do collapse(2)
    do j = JS, JE
    do i = IS, IE
       if ( xy_SIceCon(i,j) >= IceMaskMin ) then
          if ( xy_SnowThick(i,j) <= 0d0 ) then
             xy_SIceAlbedoAI(i,j) = AlbedoIce
          elseif( xy_SIceSfcTemp(i,j) >= FreezeTempWater ) then
             xy_SIceAlbedoAI(i,j) = AlbedoMeltSnow
          else
             xy_SIceAlbedoAI(i,j) = AlbedoSnow
          end if
       else ! grid cell not covered by sea ice
          xy_SIceAlbedoAI(i,j) = UNDEFVAL
       end if
    end do
    end do
    
  end subroutine DSIce_Boundary_SfcAlbedo_Const

  subroutine DSIce_Boundary_SfcAlbedo_Simple( &
       & xy_SIceAlbedoAI,                                       & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_IceThick, xy_SIceSfcTemp, & ! (in)
       & xy_SeaSfcTemp                                          &
       & )

    use DSIce_Admin_Constants_mod, only: &
         & AlbedoOcean, AlbedoMeltSnow,  &
         & FreezeTempWater, IceMaskMin
    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SIceAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)

    ! 局所変数
    ! Local variables
    !        
    integer :: i
    integer :: j
    real(DP) :: SIceSfcTempK
    real(DP) :: AlbedoIce
    real(DP) :: AlbedoSnow
    real(DP) :: r

    real(DP) :: AlbSnowCoef1
    real(DP) :: AlbSnowCoef2
    real(DP), parameter :: AlbSnowTransTempWidth = 1d-12 !10d0 !20d0
    real(DP), parameter :: PI = acos(-1d0)
    
    ! 実行文; Executable statements
    !

    AlbSnowCoef1 = (AlbedoSnowWarm - AlbedoSnowCold)/AlbSnowTransTempWidth
    AlbSnowCoef2 = 0.5d0*(AlbedoSnowWarm + AlbedoSnowCold)
     
    !$omp parallel do private(SIceSfcTempK, AlbedoIce, AlbedoSnow, r) collapse(2)
    do j = JS, JE
    do i = IS, IE
       if ( xy_SIceCon(i,j) >= IceMaskMin ) then
          SIceSfcTempK = degC2K( xy_SIceSfcTemp(i,j) )

          r = exp(- xy_IceThick(i,j)/0.25d0)
          AlbedoIce = r*AlbedoIceMin + (1d0 - r)*AlbedoIceMax

!!$          AlbedoSnow = AlbSnowCoef1*(xy_SIceSfcTemp(i,j) - TLVTempAlbedoSnow) + AlbSnowCoef2
!!$          AlbedoSnow = max(min(AlbedoSnow, AlbedoSnowCold), AlbedoSnowWarm)

          AlbedoSnow =    AlbedoSnowWarm                           &
               &        - 0.5d0*(AlbedoSnowWarm - AlbedoSnowCold)  &
               &          * (1d0 - tanh(2d0*(xy_SIceSfcTemp(i,j) - TLVTempAlbedoSnow)/AlbSnowTransTempWidth))
          
          r = exp(- xy_SnowThick(i,j)/0.05d0)
          xy_SIceAlbedoAI(i,j) = r*AlbedoIce + (1d0 - r)*AlbedoSnow
          
       else ! grid cell not covered by sea ice
          if ( xy_SeaSfcTemp(i,j) <= degC2K(TLVTempAlbedoSnow) + 4d0*AlbSnowTransTempWidth) then
!!$             xy_SIceAlbedoAI(i,j) = AlbSnowCoef1*(xy_SeaSfcTemp(i,j) - degC2K(TLVTempAlbedoSnow)) + AlbSnowCoef2
             xy_SIceAlbedoAI(i,j) =    AlbedoSnowWarm                           &
               &        - 0.5d0*(AlbedoSnowWarm - AlbedoSnowCold)  &
               &          * (1d0 - tanh(2d0*(xy_SeaSfcTemp(i,j) - degC2K(TLVTempAlbedoSnow))/AlbSnowTransTempWidth))
          else
             xy_SIceAlbedoAI(i,j) = UNDEFVAL
          end if
       end if
     end do
    end do

  end subroutine DSIce_Boundary_SfcAlbedo_Simple

  subroutine DSIce_Boundary_SfcAlbedo_I07SGSParam( &
       & xy_SIceAlbedoAI,                                           & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,   & ! (in)
       & FlagIceAreaFrac                                            & ! (in)
       & )

    use DSIce_Admin_Constants_mod, only: &
         & IceMaskMin

    use DSIce_Admin_GridDef_mod, only: &
         & I_XY, I_UY, I_XV
    
    use DSIce_Admin_Grid_mod, only: &
         & y_JAXIS_Weight, y_CJ, y_CDJ, y_FJ, y_FDJ,   &
         & SCALEF_E1, SCALEF_E2
    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SIceAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)
    logical, intent(in) :: FlagIceAreaFrac

    ! 局所変数
    ! Local variables
    !        
    integer :: i
    integer :: j
    real(DP) :: xy_SfcTemp(IA,JA)
    real(DP) :: xv_SfcTemp(IA,JA)
    real(DP) :: m1
    real(DP) :: m2
    real(DP) :: xi
    real(DP) :: grad
    real(DP) :: frac
    real(DP) :: frac1
    real(DP) :: frac2
    real(DP) :: TLVAlbedoSfcTemp
    real(DP) :: AlbedoWarm
    real(DP) :: AlbedoCold

    real(DP) :: icelat
    real(DP), parameter :: Deg2Rad = PI/180d0
    
    ! 実行文; Executable statements
    !

    TLVAlbedoSfcTemp = TLVTempAlbedoSnow
    AlbedoWarm = AlbedoOcean
    AlbedoCold = AlbedoSnowCold

    !---------------------------------------------------------
    if (.not. FlagIceAreaFrac) then
       !$omp parallel do collapse(2)
       do j = JS, JE
       do i = IS, IE
          if ( xy_SIceCon(i,j) >= IceMaskMin ) then
             if( xy_SIceSfcTemp(i,j) >= TLVAlbedoSfcTemp ) then
                xy_SIceAlbedoAI(i,j) = AlbedoWarm
             else
                xy_SIceAlbedoAI(i,j) = AlbedoCold
             end if
          else ! grid cell not covered by sea ice
             xy_SIceAlbedoAI(i,j) = UNDEFVAL
          end if
       end do
       end do
    
       return
    end if
    !---------------------------------------------------------
    
    where(xy_SIceCon >= IceMaskMin)
       xy_SfcTemp = xy_SIceSfcTemp
    elsewhere
       xy_SfcTemp = K2degC(xy_SeaSfcTemp)
    end where
    xy_SfcTemp(:,JS-1) = xy_SfcTemp(:,JS)
    xy_SfcTemp(:,JE+1) = xy_SfcTemp(:,JE)
    
    !$omp parallel do private(grad, icelat, frac) collapse(2)
    do j = JS, JE
    do i = IS, IE

       icelat = 1d20
       if ( (xy_SfcTemp(i,j  ) - TLVAlbedoSfcTemp)*(xy_SfcTemp(i,j+1) - TLVAlbedoSfcTemp) < 0d0 ) then
          grad = (xy_SfcTemp(i,j+1) - xy_SfcTemp(i,j))/(y_CJ(j+1) - y_CJ(j))
          icelat = (TLVAlbedoSfcTemp - xy_SfcTemp(i,j))/grad + y_CJ(j)
       end if
       if ( (xy_SfcTemp(i,j-1) - TLVAlbedoSfcTemp)*(xy_SfcTemp(i,j  ) - TLVAlbedoSfcTemp) < 0d0 ) then
          grad = (xy_SfcTemp(i,j) - xy_SfcTemp(i,j-1))/(y_CJ(j) - y_CJ(j-1))
          icelat = (TLVAlbedoSfcTemp - xy_SfcTemp(i,j-1))/grad + y_CJ(j-1)
       end if
       
       if ( y_FJ(j-1) < icelat .and. icelat <= y_FJ(j) ) then
          if (grad > 0d0) then
             frac = (sin(icelat*Deg2Rad) - sin(y_FJ(j-1)*Deg2Rad))/y_JAXIS_Weight(j)
          else             
             frac = (sin(y_FJ(j  )*Deg2Rad) - sin(icelat*Deg2Rad))/y_JAXIS_Weight(j)
          end if
       else if ( xy_SfcTemp(i,j) - TLVAlbedoSfcTemp > 0d0 ) then
          frac = 0d0
       else
          frac = 1d0
       end if

       xy_SIceAlbedoAI(i,j) = frac*AlbedoCold + (1d0 - frac)*AlbedoWarm
    end do
    end do
    
!!$    xv_SfcTemp(:,JS-1) = xy_SfcTemp(:,JS)
!!$    xv_SfcTemp(:,JE  ) = xy_SfcTemp(:,JE)        
!!$
!!$    !$omp parallel private(m1, m2, grad, xi, frac1, frac2, frac)
!!$    !$omp do
!!$    do j = JS, JE-1
!!$    do i = IS, IE
!!$       m1 = y_JAXIS_Weight(j+1)/(y_JAXIS_Weight(j) + y_JAXIS_Weight(j+1))
!!$       m2 = y_JAXIS_Weight(j  )/(y_JAXIS_Weight(j) + y_JAXIS_Weight(j+1))
!!$       xv_SfcTemp(i,j) = m1*xy_SfcTemp(i,j) + m2*xy_SfcTemp(i,j+1)
!!$    end do
!!$    end do
!!$ 
!!$    !$omp do
!!$    do j = JS, JE
!!$    do i = IS, IE
!!$
!!$       !--
!!$       grad = (xv_SfcTemp(i,j) - xy_SfcTemp(i,j))
!!$       xi = (TLVAlbedoSfcTemp - xy_SfcTemp(i,j))/grad
!!$       frac1 = 1d0
!!$       if (abs(grad) > 1d-10 .and. (xi <= 1d0 .and. xi >= 0d0)) then
!!$          if (xy_SfcTemp(i,j) + grad >= TLVAlbedoSfcTemp) then
!!$             frac1 = xi
!!$          else
!!$             frac1 = 1d0 - xi
!!$          end if
!!$       else
!!$          if (xy_SfcTemp(i,j) > TLVAlbedoSfcTemp) frac1 = 0d0
!!$       end if
!!$
!!$       !--
!!$       grad = (xy_SfcTemp(i,j) - xv_SfcTemp(i,j-1))
!!$       xi = (TLVAlbedoSfcTemp - xv_SfcTemp(i,j-1))/grad
!!$       frac2 = 1d0
!!$       if (abs(grad) > 1d-10 .and. (xi <= 1d0 .and. xi >= 0d0)) then
!!$          if (xv_SfcTemp(i,j-1) + grad >= TLVAlbedoSfcTemp) then
!!$             frac2 = xi
!!$          else
!!$             frac2 = 1d0 - xi
!!$          end if
!!$       else
!!$          if (xy_SfcTemp(i,j) > TLVAlbedoSfcTemp) frac2 = 0d0
!!$       end if
!!$
!!$       frac = 0.5d0*(frac1 + frac2)
!!$       xy_SIceAlbedoAI(i,j) = frac*AlbedoCold + (1d0 - frac)*AlbedoWarm
!!$    end do
!!$    end do
!!$    !$omp end parallel
    
  end subroutine DSIce_Boundary_SfcAlbedo_I07SGSParam
  
!-----------------------------------------------------------

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

    !
    use dc_string, only: Split, Replace, StrInclude

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read

    character(TOKEN) :: SIceSfcAlbedoName
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /SeaIce_Boundary_SfcAlbedo_nml/ &
         & SIceSfcAlbedoName,                &
         & AlbedoSnowWarm,                   &
         & AlbedoSnowCold,                   &
         & AlbedoIceMin,                     &
         & AlbedoIceMax,                     &
         & DampTimeAlbedoSnow,               &
         & TLVTempAlbedoSnow,                &
         & DampTimeAlbedoSnow,               &
         & SFCALBEDO_I07SGS_FlagIceAreaFrac
    

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    SIceSfcAlbedoName = SICEGOVERNEQ_SFCALBEDO_CONST_NAME

    AlbedoSnowCold      = 0.85d0
    AlbedoSnowWarm      = 0.50d0
    AlbedoSnowOld       = 0.50d0
    AlbedoIceMin        = 0.70d0
    AlbedoIceMax        = 0.30d0
    TLVTempAlbedoSnow   = - 5d0
    DampTimeAlbedoSnow  = 5d0 * 86400d0

    SFCALBEDO_I07SGS_FlagIceAreaFrac = .false.
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                                              ! (in)
            & nml = SeaIce_Boundary_SfcAlbedo_nml, iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    select case( SIceSfcAlbedoName )
    case( SICEGOVERNEQ_SFCALBEDO_CONST_NAME )
       SIceSfcAlbedoMethod = SICEGOVERNEQ_SFCALBEDO_CONST
    case( SICEGOVERNEQ_SFCALBEDO_SIMPLE_NAME )
       SIceSfcAlbedoMethod = SICEGOVERNEQ_SFCALBEDO_SIMPLE
    case( SICEGOVERNEQ_SFCALBEDO_I07SGS_NAME )
       SIceSfcAlbedoMethod = SICEGOVERNEQ_SFCALBEDO_I07SGS
    case default
       call MessageNotify( 'E', module_name, &
            & "Unexpected method of sea-ice surface albedo calculation, '%c' is specified. Check!", &
            & ca=(/ SIceSfcAlbedoName /) )
    end select

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  - SIceSfcAlbedoName           = %c', c1 = SIceSfcAlbedoName ) 
    if ( SIceSfcAlbedoMethod == SICEGOVERNEQ_SFCALBEDO_SIMPLE ) then
       call MessageNotify( 'M', module_name, ' - AlbedoSnow(warm,cold,old) = (%f,%f,%f)', &
            & d=(/ AlbedoSnowWarm, AlbedoSnowCold, AlbedoSnowOld /) )
       call MessageNotify( 'M', module_name, ' - TLVTempAlbedoSnow         = %f   [K]',   &
            & d=(/ TLVTempAlbedoSnow /) )
       call MessageNotify( 'M', module_name, ' - DampTimeAlbedoSnow        = %f   [sec]', &
            & d=(/ DampTimeAlbedoSnow /) )
       call MessageNotify( 'M', module_name, ' - AlbedoIce(max,min)        = (%f,%f)',    &
            & d=(/ AlbedoIceMax, AlbedoIceMin /) )       
    end if
    if ( SIceSfcAlbedoMethod == SICEGOVERNEQ_SFCALBEDO_I07SGS ) then
       call MessageNotify( 'M', module_name, ' - SFCALBEDO_I07SGS_FlagIceAreaFrac = %b',  &
            & l=(/ SFCALBEDO_I07SGS_FlagIceAreaFrac /) )
    end if
    
  end subroutine read_nmlData
  
end module DSIce_Boundary_SfcAlbedo_mod
