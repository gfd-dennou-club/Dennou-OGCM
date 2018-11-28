!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Dyn_fvm_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/SIce

  use DOGCM_Admin_Constants_mod, only: &
       & RPlanet
  
  use UnitConversion_mod, only: &
       & degC2K, K2degC

  use DSIce_Admin_GridDef_mod, only: &
       & XDIR, YDIR, &
       & I_XY, I_UY, I_XV
  
  use DSIce_Admin_Grid_mod, only: &
       & IS, IE, IA, x_IAXIS_Weight,   &
       & JS, JE, JA, y_JAXIS_Weight,   &
       & KS, KE, KA, z_KAXIS_Weight,   &
       & x_CDI, x_FDI, y_CDJ, y_FDJ,   &
       & SCALEF_E1, SCALEF_E2,         &
       & x_IAXIS_Weight, y_JAXIS_Weight

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Dyn_fvm_Init, DSIce_Dyn_fvm_Final
  public :: DSIce_Dyn_fvm_SIceThickDiffRHS
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Dyn_fvm_mod' !< Module Name


  logical :: initedFlag = .false.
  
contains

  !>
  !!
  !!
  subroutine DSIce_Dyn_fvm_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName
    
    ! 実行文; Executable statements
    !

    initedFlag = .true.
    
  end subroutine DSIce_Dyn_fvm_Init

  !>
  !!
  !!
  subroutine DSIce_Dyn_fvm_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Dyn_fvm_Final
  
  !---------------------------------------------------------------
  
  subroutine DSIce_Dyn_fvm_SIceThickDiffRHS(           & 
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS, & ! (out)
       & xyz_SIceEn_RHS, xy_SIceSfcTemp_RHS,                & ! (out)
       & xy_SIceU, xy_SIceV,                                & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,          & ! (in)
       & xyz_SIceEn0                                        & ! (in)
       & )

    ! モジュール引用; Use statements
    !
    use DOGCM_Admin_Constants_mod, only: &
         & RPlanet
    
    use DSIce_Admin_Constants_mod, only: &
         & DensSnow, DensIce,            &
         & SIceHDiffCoef

    use SpmlUtil_mod
    use DSIce_Admin_TInteg_mod
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xy_SIceCon_RHS(IA,JA)
    real(DP), intent(inout) :: xy_IceThick_RHS(IA,JA)
    real(DP), intent(inout) :: xy_SnowThick_RHS(IA,JA)
    real(DP), intent(inout) :: xyz_SIceEn_RHS(IA,JA,KA)
    real(DP), intent(inout) :: xy_SIceSfcTemp_RHS(IA,JA)
    real(DP), intent(inout) :: xy_SIceU(IA,JA)
    real(DP), intent(inout) :: xy_SIceV(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xyz_SIceEn0(IA,JA,KA)

    ! 作業変数
    ! Work variables

    real(DP) :: xy_OcnMask(IA,JA)
    real(DP) :: xy_SIceMass(IA,JA)
    real(DP) :: SIceMassTmp

    real(DP) :: VolFlx(IA,JA,2)
    real(DP) :: AdvFlx(IA,JA,2)
    integer :: i
    integer :: j

    real(DP) :: gl_mean
    
    ! 実行文; Executable statements
    !

    
    !------------------------------------------

    xy_OcnMask(:,:) = 1d0
    
    !$omp parallel do private(i)
    do j = JS, JE
    do i = IS, IE
       xy_SIceMass(i,j) = &
            & ( DensSnow*xy_SnowThick0(i,j) + DensIce*xy_IceThick0(i,j) )!* xy_SIceCon0(i,j)
       xy_SIceU(i,j) = 0d0
       xy_SIceV(i,j) = 0d0
    end do
    end do

    !------------------------------------------

    !$omp parallel do private(i, SIceMassTmp)
    do j = JS-1, JE
    do i = IS, IE
       SIceMassTmp =  max( xy_SIceMass(i,j), xy_SIceMass(i,j+1) ) &
            &       * xy_OcnMask(i,j)*xy_OcnMask(i,j+1)
          
       if (SIceMassTmp > 1d-10) then
          xy_SIceV(i,j) = xy_SIceV(i,j) - &
               & SIceHDiffCoef*(xy_SIceMass(i,j+1) - xy_SIceMass(i,j)) &
               & / (SIceMassTmp*SCALEF_E2(i,j,I_XV)*y_FDJ(j))
       end if
    end do
    end do
    xy_SIceV(:,JS-1) = 0d0
    xy_SIceV(:,JE)   = 0d0

    !------------------------------------------

    
    call FVM_UD1( xy_SIceCon_RHS, VolFlx,           &
         & xy_SIceCon0, xy_SIceU, xy_SIceV, .false. )
    
    call FVM_UD1( xy_IceThick_RHS, AdvFlx,          &
         & xy_IceThick0, VolFlx(:,:,XDIR), VolFlx(:,:,YDIR), .true. )

    call FVM_UD1( xy_SnowThick_RHS, AdvFlx,         &
         & xy_SnowThick0, VolFlx(:,:,XDIR), VolFlx(:,:,YDIR), .true. )

    call FVM_UD1( xyz_SIceEn_RHS(:,:,KS  ), AdvFlx, &
         & xyz_SIceEn0(:,:,KS  ), VolFlx(:,:,XDIR), VolFlx(:,:,YDIR), .true. )

    call FVM_UD1( xyz_SIceEn_RHS(:,:,KS+1), AdvFlx, &
         & xyz_SIceEn0(:,:,KS+1), VolFlx(:,:,XDIR), VolFlx(:,:,YDIR), .true. )

!!$    gl_mean = (AvrLonLat_xy( xy_IceThick_RHS(IS:IE,JS:JE) )*DelTime) &
!!$         &    /AvrLonLat_xy( xy_IceThick0(IS:IE,JS:JE) )
!!$    write(*,*) "IceThickRHS_glmean=", gl_mean
!!$
!!$    gl_mean = 0d0
!!$    do j = JS, JE
!!$       do i = IS, IE
!!$          gl_mean = gl_mean + &
!!$               &   SCALEF_E1(i,j,I_XY)*SCALEF_E2(i,j,I_XY)*x_CDI(i)*y_CDJ(j) &
!!$               & * xy_IceThick_RHS(i,j)*DelTime
!!$       end do
!!$    end do
!!$    write(*,*) "-IceThickRHS_glmean=", gl_mean/(RPlanet**2*IntLonLat_xy(xy_IceThick0(IS:IE,JS:JE)))
    
  end subroutine DSIce_Dyn_fvm_SIceThickDiffRHS

  subroutine FVM_UD1( xy_q_RHS,  Flx,   & ! (out)
       & xy_q, uy_u, xv_v, FlagVolTrans & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_q_RHS(IA,JA)
    real(DP), intent(out) :: Flx(IA,JA,2)
    real(DP), intent(in) :: xy_q(IA,JA)
    real(DP), intent(in) :: uy_u(IA,JA)
    real(DP), intent(in) :: xv_v(IA,JA)
    logical, intent(in) :: FlagVolTrans
    
    ! 作業変数
    ! Work variables
    !
    integer :: i
    integer :: j

    ! 実行文; Executable statements
    !

    !$omp parallel

    Flx(:,:,XDIR) = 0d0
    
    !$omp do collapse(2)
    do j = JS-1, JE
    do i = IS, IE
       Flx(i,j,YDIR) =  ( &
            &    0.5d0 * (   ( xy_q(i,j+1) + xy_q(i,j) )*    xv_v(i,j)      &
            &              - ( xy_q(i,j+1) - xy_q(i,j) )*abs(xv_v(i,j))  )  &
            & )
       if (.not.FlagVolTrans) Flx(i,j,YDIR) = SCALEF_E1(i,j,I_XV)*x_CDI(i) * Flx(i,j,YDIR)
    end do
    end do


    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
!!$          xy_q_RHS(i,j) = ( &
!!$               &  - ( Flx(i,j,XDIR) - Flx(i-1,j,XDIR) )             &
!!$               &  - ( Flx(i,j,YDIR) - Flx(i,j-1,YDIR) )                &
!!$               & ) / (SCALEF_E1(i,j,I_XY)*SCALEF_E2(i,j,I_XY)*x_CDI(i)*y_CDJ(j))
       xy_q_RHS(i,j) = ( &
!!$               &  - ( Flx(i,j,XDIR) - Flx(i-1,j,XDIR) )             &
            &  - ( Flx(i,j,YDIR) - Flx(i,j-1,YDIR) )                   &
            & ) / (RPlanet**2 * x_IAXIS_Weight(i) * y_JAXIS_Weight(j))

    end do
    end do

    !$omp end parallel
    
  end subroutine FVM_UD1
      
end module DSIce_Dyn_fvm_mod
