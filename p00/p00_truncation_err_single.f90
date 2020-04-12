!--------------------------------------------------------------------------
!-- PHYS704 (2014)
!-- Program 01: Truncation Error of df/dx
!--------------------------------------------------------------------------
  implicit none
  real(kind=4)	:: h,f,f1,err,x,dfdh   ! 8 bytes (double precision)
  integer i

! open an output file with IO=10, name: err.dat
  open(10,file='test.dat')

  x=1.0

  do i=0,20
     h    = 10.d0**(-i)
     dfdh = (f(x+h)-f(x))/h     ! forward derivative
     err  = abs(f1(x)-dfdh)
     write(*,'(i3,e16.8,e16.8)') i,h,err
     write(10,'(i3,e16.8,e16.8)') i,h,err
  enddo

! close the output file
  close(10)

  stop ; end

!--------------------------------------------------------------------------
! function f(x) = x*x
!--------------------------------------------------------------------------
function f(x)
  implicit none
  real(kind=4)	:: f,x
  f=sin(x)
  return
end function f

!--------------------------------------------------------------------------
! function f'(x) = 2*x
!--------------------------------------------------------------------------
function f1(x)
  implicit none
  real(kind=4)	:: f1,x
  f1=cos(x)
  return
end function f1
