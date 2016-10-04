!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module ShoreLine_mod 

  ! モジュール引用; Use statements
  !
  use dc_types,only: DP
  use dc_message, only: MessageNotify

  use VectorSpace_mod
  use SphericalCoord_mod

  use MeshUsageMask_mod, only: &
       & MeshUsageMask, MeshUsageMask_Init, MeshUsageMask_Final, &
       & MaskToOriCoord, OriToMaskCoord, &
       & DEFAULT_MASKNX, DEFAULT_MASKNY

  use SCVoronoiGen_mod, only: &
       & ClosedCurve, ClosedCurve_Init, ClosedCurve_Final

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  type(ClosedCurve), target, public, save :: shoreLines(1)
  type(MeshUsageMask), public, save :: usageMask
  
  public :: ShoreLine_Init, ShoreLine_Final
  public :: set_ShoreLine, get_gridDensityField

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'ShoreLine_mod' !< Module Name

  real(DP), parameter :: PI = acos(-1d0)
  real(DP), parameter :: lon1 = -90d0   *PI/180d0
  real(DP), parameter :: lon2 =  90d0 *PI/180d0
  real(DP), parameter :: lat1 = -90d0 *PI/180d0
  real(DP), parameter :: lat2 =  90d0  *PI/180d0
  real(DP), parameter :: lon1_ = -90d0   *PI/180d0
  real(DP), parameter :: lon2_ =  90d0 *PI/180d0
  real(DP), parameter :: lat1_ = -90d0 *PI/180d0
  real(DP), parameter :: lat2_ =  90d0 *PI/180d0
  integer, parameter :: nLat = 2500
  integer, parameter :: nLon = 2500

  type(vector3d) :: domainCenterPos
  real(DP) :: domainCircRadius

  integer, parameter :: MASKNX = DEFAULT_MASKNX*4
  integer, parameter :: MASKNY = DEFAULT_MASKNY*4

contains

  !>
  !!
  !!
  subroutine ShoreLine_Init()

    real(DP) :: radius(4)
    type(vector3d) :: pts(4)
    integer :: i

    ! 実行文; Executable statements
    !

    call ClosedCurve_Init(shoreLines(1), 2*nLat)
    call MeshUsageMask_Init(usageMask, &
         & -PI, PI, -PI/2d0, PI/2d0, MASKNX, MASKNY)

    pts(1) = (/ lon1_, lat1_, 1d0  /)
    pts(2) = (/ lon1_, lat2_, 1d0  /)
    pts(3) = (/ lon2_, lat2_, 1d0 /)
    pts(4) = (/ lon2_, lat1_, 1d0 /)
    do i=1, 4
       pts(i) = SphToCartPos(pts(i))
    end do

    domainCenterPos = normalizedVec(pts(1) + pts(2) + pts(3) + pts(4))
    do i=1, 4
       radius(i) = geodesicArcLength(domainCenterPos, pts(i))
    end do
    domainCircRadius = maxval(radius)

  end subroutine ShoreLine_Init

  !>
  !!
  !!
  subroutine ShoreLine_Final()

    ! 実行文; Executable statements
    !
    call ClosedCurve_Final(shoreLines(1))
    call MeshUsageMask_Final(usageMask)

  end subroutine ShoreLine_Final

  subroutine set_ShoreLine()

    use SphericalCoord_mod, only: &
         & SphToCartPos

    real(DP) :: dLat, dLon
    integer :: i, j
    type(vector3d) :: sphPos
    type(ClosedCurve), pointer :: shoreLine

    shoreLine => shoreLines(1)

    dLat = (lat2 - lat1)/nLat
    dLon = (lon2 - lon1)/nLon
  
    do j=1,nLat
       shoreLine%pts(j) = (/ lon2, lat2-(j-1)*dLat, 1d0 /)
       shoreLine%pts(nLat+j) = (/ lon1, lat1+(j-1)*dLat, 1d0 /)
    end do

    do i=1,size(shoreLine%pts)
       shoreLine%pts(i) = SphToCartPos(shoreLine%pts(i))
    end do

    do j=1, MASKNY
       do i=1, MASKNX
          call MaskToOriCoord(usageMask, sphPos%v_(1), sphPos%v_(2), i, j)
          if(lon1 <= sphPos%v_(1) .and. lon2 >= sphPos%v_(1) &
               & .and. lat1 <= sphPos%v_(2) .and. lat2 >= sphPos%v_(2) ) then
             usageMask%masks(i,j) = 1  ! Sea 
          else
             usageMask%masks(i,j) = 2  ! Land
          end if
       end do
    end do
  
  end subroutine set_ShoreLine

  function get_gridDensityField(x) result(density)
    type(vector3d), intent(in) :: x
    real(DP) :: density

    real(DP) :: judge
    real(DP), parameter :: gammas = 3d0
    real(DP), parameter :: transitionEPS = 0.24d0
    type(vector3d) :: sphPos
    real(DP) :: dLon1, dLon2
    real(DP), parameter :: boundaryWidth1 = 8d0*PI/180d0
    real(DP), parameter :: boundaryWidth2 = 8d0*PI/180d0

density = 1d0; return;

    judge = (geodesicArcLength(x, domainCenterPos) - domainCircRadius)/transitionEPS
    sphPos = CartToSphPos(x)
    dLon1 = abs(sphPos%v_(1) - lon1)
    dLon2 = abs(sphPos%v_(1) - lon2)
    
    density =   gammas**4 * exp(- dLon1/boundaryWidth1)  &
         &    + gammas**2 * exp(- dLon2/boundaryWidth2)  &
         &    + 1d0

!!$    if(judge < 0d0) then
!!$       density = gammas**4
!!$    else if(judge < 1d0) then
!!$       density = ((1d0 - judge)*gammas + judge)**4
!!$    else
!!$       density = 1d0
!!$    end if

  end function get_gridDensityField

end module ShoreLine_mod

