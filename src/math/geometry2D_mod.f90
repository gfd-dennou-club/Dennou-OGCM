!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module geometry2D_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only:DP
  use VectorSpace_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: getDistance_LineSeg_Pt
  public :: getIntersectPt_lineSegs
  public :: getVerticalPt
  public :: getCircumCircCenter
  public :: isLineSegsIntersect

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'geometry2D_mod' !< Module Name

contains

function getDistance_LineSeg_Pt(lineSegPt1, lineSegPt2, Pt) result(dist)
  type(Vector2d), intent(in) :: lineSegPt1, lineSegPt2, Pt
  real(DP) :: dist

  type(vector2d) :: x, y

  x = Pt - lineSegPt1
  y = lineSegPt2 - lineSegPt1
  
  dist = l2norm( (x.cross.y) )/l2norm( y )

end function getDistance_LineSeg_Pt

function getIntersectPt_lineSegs(lineSeg1Pt1, lineSeg1Pt2, lineSeg2Pt1, lineSeg2Pt2) result(Pt)
  type(Vector2d), intent(in) :: lineSeg1Pt1, lineSeg1Pt2
  type(Vector2d), intent(in) :: lineSeg2Pt1, lineSeg2Pt2
  type(Vector2d) :: Pt

  real(DP) :: dist_lineSeg1_lineSeg2Pt1, dist_lineSeg1_lineSeg2Pt2

  dist_lineSeg1_lineSeg2Pt1 = getDistance_LineSeg_Pt(lineSeg1Pt1, lineSeg1Pt2, lineSeg2Pt1)
  dist_lineSeg1_lineSeg2Pt2 = getDistance_LineSeg_Pt(lineSeg1Pt1, lineSeg1Pt2, lineSeg2Pt2)

  Pt =  lineSeg2Pt1 + &
       &   (dist_lineSeg1_lineSeg2Pt1/(dist_lineSeg1_lineSeg2Pt1+dist_lineSeg1_lineSeg2Pt2)) &
       & * (lineSeg2Pt2 - lineSeg2Pt1)

end function getIntersectPt_lineSegs

function getVerticalPt(lineSegPt1, lineSegPt2, Pt) result(verticalPt)
  type(Vector2d), intent(in) :: lineSegPt1, lineSegPt2, Pt
  type(Vector2d) :: verticalPt

  type(Vector2d) :: x, y

  x = Pt - lineSegPt1
  y = lineSegPt2 - lineSegPt1
  verticalPt = lineSegPt1 + ((x.dot.y)/(y.dot.y))*y

end function getVerticalPt

function getCircumCircCenter(pts) result(centerPt)

  type(Vector2d), intent(in) :: pts(3)
  type(Vector2d) :: centerPt

  real(DP) :: a2, b2 ,c2, d1, d2, d3

  a2 = l2norm(pts(2)-pts(3))**2
  b2 = l2norm(pts(3)-pts(1))**2
  c2 = l2norm(pts(1)-pts(2))**2

  d1 = a2*(b2 + c2 - a2)
  d2 = b2*(c2 + a2 - b2)
  d3 = c2*(a2 + b2 - c2)

  centerPt = (d1*pts(1) + d2*pts(2) + d3*pts(3))/(d1 + d2 + d3)

end function getCircumCircCenter

function isLineSegsIntersect(lineSeg1Pt1, lineSeg1Pt2, lineSeg2Pt1, lineSeg2Pt2) result(judgeFlag)
  type(Vector2d), intent(in) :: lineSeg1Pt1, lineSeg1Pt2
  type(Vector2d), intent(in) :: lineSeg2Pt1, lineSeg2Pt2
  logical :: judgeFlag

  type(vector2d) :: v, v1, v2
  real(DP) :: v1_Cross_v2, t1, t2

  v = lineSeg2Pt1 - lineSeg1Pt1
  v1 = lineSeg1Pt2 - lineSeg1Pt1
  v2 = lineSeg2Pt2 - lineSeg2Pt1

  v1_Cross_v2 = cross(v1, v2)
  t1 = cross(v,v2)/v1_Cross_v2
  t2 = cross(v,v1)/v1_Cross_v2
  if(t1 >= 0d0 .and. t1 <= 1d0 .and. t2 >= 0d0 .and. t2 <= 1d0) then
     judgeFlag = .true.
  else
     judgeFlag = .false.
  end if

contains
  function cross(a,b) result(ret)
    type(vector2d), intent(in) :: a, b
    real(DP) :: ret

    ret = a%v_(1)*b%v_(2) - a%v_(2)*b%v_(1)
  end function cross
end function isLineSegsIntersect
end module geometry2D_mod

