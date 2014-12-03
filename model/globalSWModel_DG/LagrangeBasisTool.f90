program LagrangeBasisTool
  use dc_types
  use dc_message
  implicit none

  integer, parameter :: PolyDegree = 3
  integer, parameter :: Np = (PolyDegree+1)*(PolyDegree+2)/2

  real(DP), parameter :: pts_x(Np) = (/ 0d0,   0.5d0*(1d0-1d0/sqrt(5d0)), 0.5d0*(1d0+1d0/sqrt(5d0)),                 1d0,  &
       &                                       0.5d0*(1d0+1d0/sqrt(5d0)), 0.5d0*(1d0-1d0/sqrt(5d0)),                 0d0,  &
       &                                                                                        0d0,                 0d0,  &
       &                                                                                                         1d0/3d0 /)

  real(DP), parameter :: pts_y(Np) = (/ 0d0,                   0d0,                 0d0,                 0d0,  &
       &                                            0.5d0*(1d0-1d0/sqrt(5d0)), 0.5d0*(1d0+1d0/sqrt(5d0)),                 1d0,  &
       &                                                                 0.5d0*(1d0+1d0/sqrt(5d0)), 0.5d0*(1d0-1d0/sqrt(5d0)),  &
       &                                                                                                         1d0/3d0 /)

!!$  real(DP), parameter :: pts_x(Np) = (/ 0d0,   1d0,   0d0,  &
!!$       &  0.5d0*(1d0+1d0/sqrt(5d0)), 0.5d0*(1d0-1d0/sqrt(5d0)), &
!!$       &                        0d0,                       0d0, &
!!$       &  0.5d0*(1d0-1d0/sqrt(5d0)), 0.5d0*(1d0+1d0/sqrt(5d0)), &
!!$       &                                              1d0/3d0 /)
!!$
!!$  real(DP), parameter :: pts_y(Np) = (/ 0d0,   0d0,   1d0,  &
!!$       &  0.5d0*(1d0-1d0/sqrt(5d0)), 0.5d0*(1d0+1d0/sqrt(5d0)), &
!!$       &  0.5d0*(1d0+1d0/sqrt(5d0)), 0.5d0*(1d0-1d0/sqrt(5d0)), &
!!$       &                        0d0,                       0d0, &
!!$       &                                              1d0/3d0 /)

  real(DP) :: VadermondeMat(Np,Np), IMat(Np,Np)
  real(DP) :: LagrangeBasisCoefMat(Np,Np)
  real(DP) :: x, y
  integer :: i, l, m

  character(*), parameter :: module_name = 'LagrangeBasisTool'

  do i=1, Np
     x = pts_x(i); y = pts_y(i)
     VadermondeMat(i,:) = (/ 1d0, x, y, x**2, x*y, y**2, x**3, x**2*y, x*y**2, y**3 /) 
  end do

  LagrangeBasisCoefMat = inverseMat(VadermondeMat)
  do i=1, Np
     write(*,*) "* i=", i, "--------"
     write(*,*) "---    1:", LagrangeBasisCoefMat(1,i)
     write(*,*) "---    x:", LagrangeBasisCoefMat(2,i)
     write(*,*) "---    y:", LagrangeBasisCoefMat(3,i)
     write(*,*) "---  x^2:", LagrangeBasisCoefMat(4,i)
     write(*,*) "---   xy:", LagrangeBasisCoefMat(5,i)
     write(*,*) "---  y^2:", LagrangeBasisCoefMat(6,i)
     write(*,*) "---  x^3:", LagrangeBasisCoefMat(7,i)
     write(*,*) "--- x^2y:", LagrangeBasisCoefMat(8,i)
     write(*,*) "--- xy^2:", LagrangeBasisCoefMat(9,i)
     write(*,*) "---  y^3:", LagrangeBasisCoefMat(10,i)
  end do

!!$  where(abs(LagrangeBasisCoefMat) < 1d-13)
!!$     LagrangeBasisCoefMat = 0d0
!!$  end where

  open(10, file='LagrangePolyn_TriNk3.dat', status='replace')
  do i=1, Np
     write(10,*) LagrangeBasisCoefMat(:,i)
  end do
  close(10)

contains

  function inverseMat(A) result(Ainv)
    real(DP), intent(in) :: A(:,:)
    real(DP) :: Ainv(size(A,1),size(A,2))

    real(DP) :: work(size(A,2)*128)
    integer :: IPIV(size(A,2))
    integer :: INFO
    integer :: m, n, i

    m = size(A,1)
    n = size(A,2)

    Ainv = A
    call DGETRF(m, n, Ainv, m, IPIV(1:min(m,n)), info)
    call DGETRI(n, Ainv, m, IPIV, work, n*64, info)

    if(info /= 0) then
       write(*,*) shape(A)
       do i=1, m
          write(*,*) A(i,:)
       end do
       call MessageNotify('E', module_name, &
            & "Fail to inverse a matrix." )

!!$    else
!!$       call MessageNotify('M', module_name, &
!!$            & "Success to inverse a matrix." )
    end if

  end function inverseMat

end program LagrangeBasisTool
