!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to set some variables (e.g., surface and bottom flux) in order to satisfy boundary conditions
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Boundary_common_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod, only: &
       & RefDens, Cp0
    

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax

  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Admin_BC_mod, only: &
       & inquire_VBCSpecType,                                                    &
       & DynBCTYPE_NoSlip, DynBCTYPE_Slip, DynBCTYPE_SpecStress,                 &
       & ThermBCTYPE_PrescFixedFlux, ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, &
       & ThermBCTYPE_PresFlux_Han1984Method,                                     &
       & ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed,                         & 
       & SaltBCTYPE_PrescFixedFlux, SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat,    &
       & SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed,                           &
       & SaltBCTYPE_PresFlux_Han1984Method,                                      &
       & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface,          &
       & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom,              &
       & SeaSfcTempRelaxedTime, SeaSfcSaltRelaxedTime

  use DOGCM_Boundary_vars_mod, only: &
       & xy_SfcHFlx_ns, xy_SfcHFlx_sr, xy_DSfcHFlxDTs,   &
       & xy_SfcHFlx0_ns, xy_SfcHFlx0_sr,                 &
       & xy_SfcHFlxIO_ns, xy_SfcHFlxIO_sr,               &
       & xy_WindStressU, xy_WindStressV,      &
       & xy_FreshWtFlx, xy_FreshWtFlxS,       &
       & xy_FreshWtFlx0, xy_FreshWtFlxS0,     &       
       & xy_FreshWtFlxIO, xy_FreshWtFlxSIO,     &       
       & xy_SeaSfcTemp, xy_SeaSfcSalt,        &
       & xy_SeaSfcTemp0, xy_SeaSfcSalt0,      &
       & xy_SeaSfcU, xy_SeaSfcV,              &
       & xy_OcnSfcCellMask, xyz_OcnCellMask,  &
       & OCNCELLMASK_SICE, OCNCELLMASK_OCEAN
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_Boundary_common_Init, DOGCM_Boundary_common_Final

  public :: DOGCM_Boundary_common_UpdateBeforeTstep  
  public :: DOGCM_Boundary_common_UpdateAfterTstep
  
  ! 公開変数
  ! Public variable
  !

  !< The depth of mixed layer near sea surface to specify Haney-type boundary condition.
  real(DP), parameter :: MixLyrDepthConst = 50d0

  !< Reference salinity to calculate virtual salinity flux. 
  real(DP), public, parameter :: RefSalt_VBC = 35d0
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DOGCM_Boundary_common_mod' !< Module Name

  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Boundary_common_Init( &
       & configNmlName )                   ! (in)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    
  end subroutine DOGCM_Boundary_common_Init

  !>
  !!
  !!
  subroutine DOGCM_Boundary_common_Final()

    ! 実行文; Executable statements
    !


  end subroutine DOGCM_Boundary_common_Final

  !-----------------------------------------

  !> Update variables managed by boundary modules before the start of a time step. 
  !!
  subroutine DOGCM_Boundary_common_UpdateBeforeTstep(    &
       & xyz_U, xyz_V, xyza_TRC, xyz_H, xy_SSH           & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    
    ! 実行文; Executable statements
    !

    !-- Update variables associated with dynamical boundary condition

    select case(DynBC_Surface)
    case ( DynBCTYPE_SpecStress )
    case ( DynBCTYPE_NoSlip )
       xy_SeaSfcU(:,:) = 0d0
       xy_SeaSfcV(:,:) = 0d0       
    case ( DynBCTYPE_Slip )
       xy_WindStressU(:,:) = 0d0
       xy_WindStressV(:,:) = 0d0
    end select
    
    !-- Update variables associated with thermal boundary condition

    select case(ThermBC_Surface)
    case ( ThermBCTYPE_PrescTemp )
    case ( ThermBCTYPE_PrescFixedFlux )
       xy_SfcHFlx_ns(:,:) = xy_SfcHFlx0_ns
       xy_SfcHFlx_sr(:,:) = xy_SfcHFlx0_sr
    case ( ThermBCTYPE_PrescFlux      )
       !$omp parallel
       !$omp workshare
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_SfcHFlx_ns(:,:) =   xy_SfcHFlx0_ns + xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) =   xy_SfcHFlx0_sr + xy_SfcHFlxIO_sr
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_SfcHFlx_ns(:,:) = xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) = xy_SfcHFlxIO_sr
       end where
       !$omp end workshare
       !$omp end parallel
    case ( ThermBCTYPE_PresFlux_Han1984Method )
       !$omp parallel
       !$omp workshare
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_SfcHFlx_ns(:,:) =   xy_SfcHFlx0_ns &
               &               + xy_DSfcHFlxDTs(:,:)*(xyza_TRC(:,:,KS,TRCID_PTEMP) - xy_SeaSfcTemp0) &
               &               + xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) =   xy_SfcHFlx0_sr(:,:)   &
               &               + xy_SfcHFlxIO_ns
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_SfcHFlx_ns(:,:) = xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) = xy_SfcHFlxIO_sr
       end where
       !$omp end workshare
       !$omp end parallel
    case ( ThermBCTYPE_Adiabat )
       xy_SfcHFlx_ns(:,:) = 0d0
       xy_SfcHFlx_sr(:,:) = 0d0
    case ( ThermBCTYPE_TempRelaxed )
       !$omp parallel
       !$omp workshare
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_SfcHFlx_ns(:,:) = - (RefDens*Cp0*MixLyrDepthConst/SeaSfcTempRelaxedTime) &
               &                 *(xy_SeaSfcTemp0 - xyza_TRC(:,:,KS,TRCID_PTEMP))     &
               &               + xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) = xy_SfcHFlxIO_sr
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_SfcHFlx_ns(:,:) = xy_SfcHFlxIO_ns
          xy_SfcHFlx_sr(:,:) = xy_SfcHFlxIO_sr
       end where
       !$omp end workshare
       !$omp end parallel
    end select
    
    !-- Update variables associated with boundary condition for salinity 

    select case(SaltBC_Surface)
    case ( SaltBCTYPE_PrescSalt )
    case ( SaltBCTYPE_PrescFixedFlux )
       xy_FreshWtFlxS(:,:) = xy_FreshWtFlxS0
    case ( SaltBCTYPE_PrescFlux      )
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_FreshWtFlxS(:,:) =   xy_FreshWtFlxS0 + xy_FreshWtFlxSIO
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_FreshWtFlxS(:,:) = xy_FreshWtFlxSIO
       end where
    case ( SaltBCTYPE_PresFlux_Han1984Method )
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_FreshWtFlxS(:,:) =   xy_FreshWtFlxS0 + xy_FreshWtFlxSIO
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_FreshWtFlxS(:,:) = xy_FreshWtFlxSIO
       end where
    case ( SaltBCTYPE_Adiabat )
       xy_FreshWtFlxS(:,:) = 0d0
    case ( SaltBCTYPE_SaltRelaxed )
       where( xy_OcnSfcCellMask == OCNCELLMASK_OCEAN )
          xy_FreshWtFlxS(:,:) = - MixLyrDepthConst/(RefSalt_VBC*SeaSfcSaltRelaxedTime) &
               &                  *(xy_SeaSfcSalt0 - xyza_TRC(:,:,KS,TRCID_SALT))      &
               &                + xy_FreshWtFlxSIO(:,:)          
       elsewhere( xy_OcnSfcCellMask == OCNCELLMASK_SICE )
          xy_FreshWtFlxS(:,:) = xy_FreshWtFlxSIO(:,:)
       end where
    end select
    
    
  end subroutine DOGCM_Boundary_common_UpdateBeforeTstep

  !-----------------------------------------
  
  !> Update variables managed by boundary modules at the end of a time step. 
  !!
  subroutine DOGCM_Boundary_common_UpdateAfterTstep(    &
       & xyz_U, xyz_V, xyza_TRC, xyz_H, xy_SSH          & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    
    ! 実行文; Executable statements
    !

    !-- Update variables associated with dynamical boundary condition

    xy_SeaSfcU(:,:) = xyz_U(:,:,KS)
    xy_SeaSfcV(:,:) = xyz_V(:,:,KS)
    
    select case(DynBC_Surface)
    case ( DynBCTYPE_SpecStress )
    case ( DynBCTYPE_NoSlip )
    case ( DynBCTYPE_Slip )
    end select
    
    !-- Update variables associated with thermal boundary condition

    xy_SeaSfcTemp(:,:) = xyza_TRC(:,:,KS,TRCID_PTEMP)
    
    select case(ThermBC_Surface)
    case ( ThermBCTYPE_PrescTemp )
    case ( ThermBCTYPE_PrescFixedFlux )
    case ( ThermBCTYPE_PrescFlux      )
    case ( ThermBCTYPE_PresFlux_Han1984Method )
    case ( ThermBCTYPE_Adiabat )
    case ( ThermBCTYPE_TempRelaxed )
    end select

    !-- Update variables associated with boundary condition for salinity -

    xy_SeaSfcSalt(:,:) = xyza_TRC(:,:,KS,TRCID_SALT)

    select case(SaltBC_Surface)
    case ( SaltBCTYPE_PrescSalt )
    case ( SaltBCTYPE_PrescFixedFlux )
    case ( SaltBCTYPE_PrescFlux      )
    case ( SaltBCTYPE_PresFlux_Han1984Method )
    case ( SaltBCTYPE_Adiabat )
    case ( SaltBCTYPE_SaltRelaxed )
    end select
    
  end subroutine DOGCM_Boundary_common_UpdateAfterTstep
  
end module DOGCM_Boundary_common_mod
