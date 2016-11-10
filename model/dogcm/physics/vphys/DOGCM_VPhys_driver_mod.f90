module DOGCM_VPhys_driver_mod

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
       & KA, KS, KE, KM

  use DOGCM_Admin_Constants_mod, only: &
       & PI,                           &
       & RefDens,                      &
       & Grav,                         &
       & VDiffCoef,                    &
       & VViscCoef
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRCID_PTEMP, TRCID_SALT, &
       & TRC_TOT_NUM
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_VPhys_driver_Init, DOGCM_VPhys_driver_Final

  public :: DOGCM_VPhys_driver_UpdateVViscDiffCoef
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_VPhys_driver_mod' !< Module Name

  real(DP), save :: Av0Mom
  real(DP), save :: Av0TRC

  integer :: VMixCoef_scheme

  ! Vertical eddy viscosity and diffusivity are simply set constant values. 
  !
  character(*), parameter :: VMIXCOEF_SCHEME_CONST_NAME = 'Const'
  integer, parameter :: VMIXCOEF_SCHEME_CONST = 1

  ! Vertical eddy viscosity and diffusivity are simply set in the same way of 'Const'
  ! in ocean interior region, while much larger values are set near the sea surface to
  ! consider sea surface mixing layer. 
  !
  character(*), parameter :: VMIXCOEF_SCHEME_SIMPLE_NAME = 'Simple'
  integer, parameter :: VMIXCOEF_SCHEME_SIMPLE = 2

  ! Vertical eddy viscosity and diffusivity are calculated by TKE turbulent closure scheme
  !
  character(*), parameter :: VMIXCOEF_SCHEME_TKE_NAME = 'TKE'
  integer, parameter :: VMIXCOEF_SCHEME_TKE = 3

  ! Vertical eddy viscosity and diffusivity are calculated by KPP scheme
  !  
  character(*), parameter :: VMIXCOEF_SCHEME_KPP_NAME = 'KPP'
  integer, parameter :: VMIXCOEF_SCHEME_KPP = 4

  ! Vertical eddy viscosity and diffusivity are calculated by PP81 scheme
  !  
  character(*), parameter :: VMIXCOEF_SCHEME_PP81_NAME = 'PP81'
  integer, parameter :: VMIXCOEF_SCHEME_PP81 = 5
  
  
contains

  !>
  !!
  !!
  Subroutine DOGCM_VPhys_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call read_nmlData(configNmlName)

    
  end subroutine DOGCM_VPhys_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_VPhys_driver_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_VPhys_driver_Final

  !-----------------------------------------

  !>
  !!
  !!
  subroutine DOGCM_VPhys_driver_UpdateVViscDiffCoef(   &
       & xyz_VViscCoef, xyz_VDiffCoef,                   & ! (out)
       & xyz_U, xyz_V, xyza_Tracer, xyz_Z, xy_Topo       & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(out) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyza_Tracer(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)
    
    ! 実行文; Executable statements
    !

    select case (VMixCoef_scheme)
    case (VMIXCOEF_SCHEME_CONST)
       xyz_VViscCoef(:,:,:) = Av0Mom
       xyz_VDiffCoef(:,:,:) = Av0TRC
    case (VMIXCOEF_SCHEME_SIMPLE)
       call calc_vViscDiffCoef_MixLyrSimple( &
            & xyz_VViscCoef, xyz_VDiffCoef,    & ! (inout)
            & xyz_Z, xy_Topo                   & ! (in)
            & )
!!$       call add_EnhancedVDiff_coef( &
!!$            & xyz_VDiffCoef,                  & ! (inout)
!!$            & xyza_Tracer(:,:,:,TRCID_PTEMP), & ! (in)
!!$            & xyza_Tracer(:,:,:,TRCID_SALT),  & ! (in)
!!$            & xyz_Z                           & ! (in)
!!$            & )       
    end select
    
  end subroutine DOGCM_VPhys_driver_UpdateVViscDiffCoef

  !---------------------------------------------------

  subroutine calc_vViscDiffCoef_MixLyrSimple( &
       & xyz_VViscCoef, xyz_VDiffCoef,           & ! (inout)
       & xyz_Depth, xy_Topo                      & ! (in)
       & )

    real(DP), intent(inout) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(inout) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in)    :: xyz_Depth(IA,JA,KA)
    real(DP), intent(in)    :: xy_Topo(IA,JA)
    
    real(DP), parameter :: MixLyrDepth  = 40d0
    real(DP), parameter :: LInv         = 1d0/(0.1d0*MixLyrDepth)
    real(DP), parameter :: ViscCoefMax  = 5d-3
    real(DP), parameter :: DiffCoefMax  = 1d-3 

    real(DP) :: xyz_Func(IA,JA,KA)
    integer :: k

    !$omp parallel do
    do k = KS, KE
       xyz_Func(:,:,k) =  &
            &   0.5d0 - atan( ((-xyz_Depth(:,:,k)) - MixLyrDepth)*LInv )/PI             &
            & + 0.5d0 - atan( ((xy_Topo + xyz_Depth(:,:,k)) - MixLyrDepth)*LInv )/PI
       
       xyz_VViscCoef(:,:,k) = Av0Mom + ViscCoefMax*xyz_Func(:,:,k)
       xyz_VDiffCoef(:,:,k) = Av0TRc + DiffCoefMax*xyz_Func(:,:,k)
    end do
!!$    xyz_VDiffCoef(:,:,:) =   ( 0.8d0 + 1.05d0/PI*atan((abs(xyz_Depth) - 2500d0)/222.2d0) ) &
!!$         &                 * 1d-4 + DiffCoefMax*xyz_Func

  end subroutine calc_vViscDiffCoef_MixLyrSimple

  !------------

  subroutine add_EnhancedVDiff_coef( xyz_VDiffCoef, &
       & xyz_PTemp, xyz_Salt, xyz_Depth )

    use EOSDriver_mod, only: &
         & EOSDriver_Eval
    
    real(DP), intent(inout) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_Depth(IA,JA,KA)

    real(DP) :: xya_DensPotPair(IA,JA,2)
    real(DP) :: xya_RefPres(IA,JA,2)

    integer :: i
    integer :: j
    integer :: k
    logical :: xyz_UnstaleFlag(IA,JA,KA)
    real(DP) :: xyz_EVDiffCoef(IA,JA,KA)

    xyz_EVDiffCoef = 0d0
    xyz_UnstaleFlag(:,:,:) = .false.
    do k = KS, KE-1
       xya_RefPres(:,:,1) = - RefDens*Grav*0.5d0*(xyz_Depth(:,:,k) + xyz_Depth(:,:,k+1))
       xya_RefPres(:,:,2) = xya_RefPres(:,:,1)
       
       call EOSDriver_Eval( rhoEdd=xya_DensPotPair(:,:,:), & ! (out)
            & theta=xyz_PTemp(:,:,k:k+1),                & ! (in)
            & S=xyz_Salt(:,:,k:k+1),                     & ! (in)
            & p=xya_RefPres )

       where (xya_DensPotPair(:,:,1) > xya_DensPotPair(:,:,2))
          xyz_UnstaleFlag(:,:,k) = .true.
          xyz_UnstaleFlag(:,:,k+1) = .true.
          xyz_EVDiffCoef(:,:,k) = 5d-3
          xyz_EVDiffCoef(:,:,k+1) = 5d-3
       end where
    end do

    xyz_EVDiffCoef(:,JS-1,:) = xyz_EVDiffCoef(:,JS,:)
    xyz_EVDiffCoef(:,JE+1,:) = xyz_EVDiffCoef(:,JE,:)
    
    do k = KS+1, KE-1
       do j = JS, JE
          do i = IS, IE
!!$             xyz_VDiffCoef(i,j,k) = xyz_VDiffCoef(i,j,k) + &
!!$                  &   0.25d0*0.25d0*( xyz_EVDiffCoef(i,j-1,k-1) + 2d0*xyz_EVDiffCoef(i,j-1,k) + xyz_EVDiffCoef(i,j-1,k+1) ) &
!!$                  & + 0.50d0*0.25d0*( xyz_EVDiffCoef(i,j,k-1) + 2d0*xyz_EVDiffCoef(i,j,k) + xyz_EVDiffCoef(i,j,k+1)       ) &
!!$                  & + 0.25d0*0.25d0*( xyz_EVDiffCoef(i,j+1,k-1) + 2d0*xyz_EVDiffCoef(i,j+1,k) + xyz_EVDiffCoef(i,j+1,k+1) )
             xyz_VDiffCoef(i,j,k) = xyz_VDiffCoef(i,j,k) + &
                  & + 0.25d0*( xyz_EVDiffCoef(i,j,k-1) + 2d0*xyz_EVDiffCoef(i,j,k) + xyz_EVDiffCoef(i,j,k+1) )

          end do
       end do
    end do
    xyz_VDiffCoef(:,:,KS) = xyz_VDiffCoef(:,:,KS+1)
    xyz_VDiffCoef(:,:,KE) = xyz_VDiffCoef(:,:,KE-1) 

!!$    where(xyz_UnstaleFlag)
!!$       xyz_VDiffCoef(:,:,:) = xyz_VDiffCoef(:,:,:) + 3d-3
!!$    end where

  end subroutine add_EnhancedVDiff_coef
    
!!$  subroutine calc_vViscDiffCoef_PP81( &
!!$       & xyz_VViscCoef, xyz_VDiffCoef,                                     &
!!$       & xyz_U, xyz_V, xyz_DensPot, xy_Topo, vViscCoefBG, vDiffCoefBG )
!!$
!!$    real(DP), intent(inout) :: xyz_VViscCoef(IA,JA,KA)
!!$    real(DP), intent(inout) :: xyz_VDiffCoef(IA,JA,KA)
!!$    real(DP), intent(in)    :: xyz_U(IA,JA,KA)
!!$    real(DP), intent(in)    :: xyz_V(IA,JA,KA)
!!$    real(DP), intent(in)    :: xyz_DensPot(IA,JA,KA)    
!!$    real(DP), intent(in)    :: xyz_Depth(IA,JA,KA)
!!$
!!$    real(DP), parameter :: AvRic = 1d-2
!!$    real(DP), parameter :: a = 5d0
!!$    integer, parameter :: n = 2
!!$    
!!$    real(DP) :: xyz_Ri(IA,JA,KA)
!!$    real(DP) :: xyz_VViscCoefPP81(IA,JA,KA)
!!$
!!$    xyz_Ri(:,:,:) = diagnose_RicardsonNumber(xyz_U, xyz_V, xyz_DensPot, xy_totDepth)
!!$    where(xyz_Ri < 0d0)
!!$       xyz_Ri = 0d0
!!$    end where
!!$    
!!$    !$omp parallel
!!$    !$omp workshare
!!$
!!$    xyz_VViscCoefPP81(:,:,:) = AvRic/(1d0 + a*xyz_Ri)**n
!!$    xyz_VViscCoef(:,:,:) = xyz_VViscCoef + xyz_VViscCoefPP81
!!$
!!$    xyz_VDiffCoef(:,:,:) = xyz_VDiffCoef + &
!!$         & xyz_VViscCoefPP81(:,:,:)/(1d0 + a*xyz_Ri) + vDiffCoefBG
!!$
!!$    !$omp end workshare
!!$    !$omp end parallel
!!$
!!$  end subroutine calc_vViscDiffCoef_PP81

  !---------------------------------------

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
    character(TOKEN) :: VMixCoef_scheme_name
    
    ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /SGSPhys_VMixing_nml/ &
         & VMixCoef_scheme_name,   &
         & Av0Mom, Av0TRC

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    VMixCoef_scheme_name = VMIXCOEF_SCHEME_CONST_NAME
    
    Av0Mom = VViscCoef
    Av0TRC = VDiffCoef

    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                        ! (in)
            & nml = SGSPhys_VMixing_nml,     &  ! (out)
            & iostat = iostat_nml, iomsg=msg )  ! (out)
       close( unit_nml )
    end if

    !
    select case( VMixCoef_scheme_name )
    case ( VMIXCOEF_SCHEME_CONST_NAME )
       VMixCoef_scheme = VMIXCOEF_SCHEME_CONST
    case ( VMIXCOEF_SCHEME_SIMPLE_NAME )
       VMixCoef_scheme = VMIXCOEF_SCHEME_SIMPLE
    case default
       call MessageNotify( 'E', module_name, &
            & "The specified scheme for vertical mixing coeffecients (%a) is not supported. Check!", &
            & ca=(/ VMixCoef_scheme_name /) )
    end select
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, ' VMixCoef scheme      = %a', ca=(/ VMixCoef_scheme_name /))
    call MessageNotify( 'M', module_name, ' Av0Mom               = %f', d=(/ Av0Mom /))
    call MessageNotify( 'M', module_name, ' Av0TRC               = %f', d=(/ Av0TRC /))
    
  end subroutine read_nmlData

end module DOGCM_VPhys_driver_mod
