!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module Exp_WindDrivenCirc_mod 

  ! モジュール引用; Use statements
  !
  use dc_message

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod
  use PolyMesh_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  use OutputData_mod, only: &
       & OutputDataAnalysisInfo

  use DGHelper_mod
  use DGCalcusUtil_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: Exp_WindDrivenCirc_Init, Exp_WindDrivenCirc_Final
  public :: setIniCond_ExpWindDrivenCirc, callBack_EndCurrentTimeStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'Exp_WindDrivenCirc_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirc_Init()

    ! 実行文; Executable statements
    !

  end subroutine Exp_WindDrivenCirc_Init

  !>
  !!
  !!
  subroutine Exp_WindDrivenCirc_Final()

    ! 実行文; Executable statements
    !

  end subroutine Exp_WindDrivenCirc_Final
  
  
  subroutine setIniCond_ExpWindDrivenCirc()

    use SphericalCoord_mod

    integer :: nc, nk
    type(vector3d) :: sphPos, sphVel, cartVel
    type(vector3d) :: b1, b2, b_1, b_2

    wc_h = 0d0
    wc_U1 = 0d0
    wc_U2 = 0d0

    do nc=1, nDGElement
       do nk=1, nDGNodePerElem
          if(wc_DGNodeUsage(nk,nc) /= 2) then
             b_1 = get_DGElemCovariantBasis(1,DGElemInfo%node(nk),nc)
             b_2 = get_DGElemCovariantBasis(2,DGElemInfo%node(nk),nc)
             call get_DGElemContravariantBasis(b1, b2, b_1, b_2)

             sphPos = CartToSphPos(wc_DGNodePos(nk,nc))
!!$             sphVel = 0.2d0*(/ -cos( PI*(sphPos%v_(2) - 40d0/180d0*PI)/(9d0/180d0*PI) ), 0d0, 0d0 /)
             sphVel = 0.05d0*(/ -cos( 3d0*sphPos%v_(2) ), 0d0, 0d0 /)
             cartVel = SphToCartVec(sphVel, wc_DGNodePos(nk,nc))
             wc_WindStress1(nk,nc) = (cartVel.dot.b1)
             wc_WindStress2(nk,nc) = (cartVel.dot.b2)
          end if
       end do
    end do

  end subroutine setIniCond_ExpWindDrivenCirc

  subroutine callBack_EndCurrentTimeStep(tstep, wc_h, wc_hU1, wc_hU2)

    use OutputData_mod, only: &
         & lonlat_interpolate_wc, Output_FieldData
    use gtool_history

    integer, intent(in) :: tstep
    real(DP), dimension(:,:), intent(in) :: wc_h, wc_hU1, wc_hU2

  end subroutine callBack_EndCurrentTimeStep

end module Exp_WindDrivenCirc_mod

