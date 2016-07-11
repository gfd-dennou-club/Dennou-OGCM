!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Dyn_driver_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/SIce

  use UnitConversion_mod, only: &
       & degC2K, K2degC
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, JA, KA, &
       & x_IAXIS_Weight, y_JAXIS_Weight

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Dyn_driver_Init, DSIce_Dyn_driver_Final
  public :: DSIce_Dyn_driver_Do
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Dyn_driver_mod' !< Module Name


  logical :: initedFlag = .false.
  
contains

  !>
  !!
  !!
  subroutine DSIce_Dyn_driver_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName
    
    ! 実行文; Executable statements
    !

    initedFlag = .true.
    
  end subroutine DSIce_Dyn_driver_Init

  !>
  !!
  !!
  subroutine DSIce_Dyn_driver_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Dyn_driver_Final

  subroutine DSIce_Dyn_driver_Do()
  end subroutine DSIce_Dyn_driver_Do
  
!!$  !> @brief 
!!$  !!
!!$  !!
!!$  subroutine DSIce_Dyn_driver_Advance( &
!!$       & xy_SIceCon, xy_IceThick, xy_SnowThick, xya_SIceTemp, xy_SIceSurfTemp, & ! (inout)
!!$       & DelTime )                                                               ! (in)
!!$
!!$    use GridSet_mod, only: &
!!$         & xyz_Lat
!!$    
!!$    use SeaIceThermDyn_Winton2000_mod, only: &
!!$         & calc_E_IceLyr1, calc_E_IceLyr2, &
!!$         & calc_Temp_IceLyr1, calc_Temp_IceLyr2
!!$
!!$    use SpmlUtil_mod
!!$    
!!$    ! 宣言文; Declaration statement
!!$    !
!!$    
!!$    real(DP), intent(inout) :: xy_SIceCon(0:iMax-1,jMax)
!!$    real(DP), intent(inout) :: xy_IceThick(0:iMax-1,jMax)
!!$    real(DP), intent(inout) :: xy_SnowThick(0:iMax-1,jMax)
!!$    real(DP), intent(inout) :: xya_SIceTemp(0:iMax-1,jMax,2)
!!$    real(DP), intent(inout) :: xy_SIceSurfTemp(0:iMax-1,jMax)
!!$    real(DP), intent(in) :: DelTime
!!$    
!!$    ! 局所変数
!!$    ! Local variables
!!$    !
!!$    
!!$    real(DP) :: xy_q1(0:iMax-1,jMax)
!!$    real(DP) :: xy_q2(0:iMax-1,jMax)
!!$    real(DP) :: xy_mflxLon(0:iMax-1,jMax)
!!$    real(DP) :: xy_mflxLat(0:iMax-1,jMax)
!!$
!!$    real(DP) :: iceVol
!!$    real(DP) :: Ice0En1, Ice0En2
!!$    real(DP) :: xy_CosLat(0:iMax-1,jMax)
!!$
!!$    integer :: i, j
!!$
!!$    integer, parameter :: DEBUG_j = -1
!!$    
!!$    ! 実行文; Executable statement
!!$    !
!!$
!!$    xy_CosLat = cos(xyz_Lat(:,:,0))
!!$
!!$    !$omp parallel do private(i)
!!$    do j=1, jMax
!!$       do i=0, iMax-1
!!$          if (xy_SIceCon(i,j) >= IceMaskMin) then
!!$             xy_q1(i,j) = calc_E_IceLyr1(xya_SIceTemp(i,j,1), SaltSeaIce)
!!$             xy_q2(i,j) = calc_E_IceLyr2(xya_SIceTemp(i,j,2), SaltSeaIce)
!!$          else
!!$             xy_q1(i,j) = 0d0
!!$             xy_q2(i,j) = - LFreeze
!!$          end if
!!$
!!$          if ( j == DEBUG_j ) then
!!$             write(*,*) "A*=", xy_SIceCon(i,j)
!!$             write(*,*) "q1*, q2*=", xy_q1(i,j), xy_q2(i,j)
!!$          end if
!!$       end do
!!$    end do
!!$
!!$!write(*,*) "SiceDiff Before:", AvrLonLat_xy( xy_IceThick) 
!!$    call apply_diff_term( xy_IceThick  )  ! (inout)
!!$    call apply_diff_term( xy_SnowThick )  ! (inout)
!!$    call apply_diff_term( xy_SIceCon   )  ! (inout)
!!$    call apply_diff_term( xy_q1        )  ! (inout)
!!$    call apply_diff_term( xy_q2        )  ! (inout)
!!$!write(*,*) "SiceDiff After:", AvrLonLat_xy( xy_IceThick) 
!!$
!!$    
!!$    !$omp parallel do private(i)
!!$    do j=1, jMax
!!$       do i=0, iMax-1 
!!$          if (xy_SIceCon(i,j) > 0d0) then
!!$             xya_SIceTemp(i,j,1) = calc_Temp_IceLyr1(xy_q1(i,j), SaltSeaIce)
!!$             xya_SIceTemp(i,j,2) = calc_Temp_IceLyr2(xy_q2(i,j), SaltSeaIce)
!!$             if (isNan(xya_SIceTemp(i,j,1)) .or. isNan(xya_SIceTemp(i,j,2)) ) then
!!$                write(*,*) "Nan T1: i,j=", i,j  
!!$                write(*,*) "q1=", xy_q1(i,:)
!!$                write(*,*) "q2=", xy_q1(i,:) 
!!$                write(*,*) "IceThick=", xy_IceThick(i,:)
!!$                write(*,*) "SIceCon=", xy_SIceCon(i,:)
!!$                stop
!!$             end if
!!$          else
!!$             xya_SIceTemp(i,j,:) = UNDEFVAL
!!$          end if
!!$
!!$          if ( j == DEBUG_j ) then
!!$             write(*,*) "q1**, q2**=", xy_q1(i,j), xy_q2(i,j)
!!$             write(*,*) "T1**, T2**=", xya_SIceTemp(i,j,:)
!!$          end if          
!!$       end do
!!$    end do
!!$
!!$  contains
!!$    subroutine apply_diff_term(xy_phi)
!!$
!!$      
!!$      real(DP), intent(inout) :: xy_phi(0:iMax-1,jMax)
!!$      
!!$      real(DP) :: xy_FlxLat(0:iMax-1,0:jMax)
!!$      real(DP) :: DLon
!!$      integer :: i, j
!!$      real(DP), parameter :: EPSILL = 1d-12
!!$      
!!$      DLon = 2d0 * PI / dble(iMax)
!!$
!!$      xy_FlxLat(:,0) = 0d0; xy_FlxLat(:,jMax) = 0d0
!!$
!!$      !$omp parallel do private(i)
!!$      do j=1, jMax-1
!!$         do i=0, iMax-1
!!$            xy_FlxLat(i,j) = &
!!$                  RPlanet*cos( 0.5d0 * (xyz_Lat(i,j+1,0) + xyz_Lat(i,j,0)) )*DLon      &
!!$               *  SIceHDiffCoef *   (xy_phi(i,j+1) - xy_phi(i,j))                      & 
!!$                                  / ( RPlanet * (xyz_Lat(i,j+1,0) - xyz_Lat(i,j,0)) )
!!$         end do
!!$      end do
!!$
!!$      !$omp parallel do private(i)
!!$      do j=1, jMax
!!$         do i=0, iMax-1
!!$            xy_phi(i,j) = xy_phi(i,j) &
!!$                + DelTime * (xy_FlxLat(i,j) - xy_FlxLat(i,j-1)) &
!!$                          / ( RPlanet**2 * x_Lon_Weight(i) * y_Lat_Weight(j) )
!!$         end do
!!$      end do
!!$      
!!$    end subroutine apply_diff_term
!!$    
!!$  end subroutine DSIce_Dyn_driver_Advance

end module DSIce_Dyn_driver_mod

