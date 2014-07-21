module m_vib
implicit none

contains

subroutine get_gaussian_derivatives(c,dx,zero,first,second)
!the form of the gaussian is exp(-dx^2*c^2/2)
implicit none

real(8),intent(in) :: c, dx

real(8),intent(out) :: zero,first,second

zero   =  exp(-0.5d0*c**2*dx**2)

first  = -zero*c**2*dx

second =  zero*c**2*(c**2*dx**2 - 1d0)

end subroutine get_gaussian_derivatives

subroutine get_hermite_derivatives(degree,alpha,x,zero,first,second)
implicit none

integer,intent(in) :: degree

real(8),intent(in) :: x,alpha

real(8),intent(out) :: zero,first,second

zero = get_hermite_polynomial(degree,x)

first = 2d0*degree*get_hermite_polynomial(degree-1,x)*alpha

second = 4d0*degree*(degree - 1)*get_hermite_polynomial(degree-2,x)*alpha**2

end subroutine get_hermite_derivatives

!not used, but faster than get_hermite_polynomial(
function get_hermite_polynomial_o(degree,x) result(h)
implicit none

integer,intent(in) :: degree

real(8),intent(in) :: x

real(8) :: h

select case(degree)
case(:-1)
   print *, 'error, values in Hermite degree cant be negative'
   stop
case(0)
   h =    1d0
case(1)
   h =    2d0*x
case(2)
   h =    4d0*x**2  -     2d0
case(3)
   h =    8d0*x**3  -    12d0*x
case(4)
   h =   16d0*x**4  -    48d0*x**2 +     12d0
case(5)
   h =   32d0*x**5  -   160d0*x**3 +    120d0*x
case(6)
   h =   64d0*x**6  -   480d0*x**4 +    720d0*x**2 -    120d0
case(7)
   h =  128d0*x**7  -  1344d0*x**5 +   3360d0*x**3 -   1680d0*x
case(8)
   h =  256d0*x**8  -  3584d0*x**6 +  13440d0*x**4 -  13440d0*x**2 +   1680d0
case(9)
   h =  512d0*x**9  -  9216d0*x**7 +  48384d0*x**5 -  80640d0*x**3 +  30240d0*x
case(10)
   h = 1024d0*x**10 - 23040d0*x**8 + 161280d0*x**6 - 403200d0*x**4 + 302400d0*x**2 - 30240d0
case(11:)
   print *, 'error, values in Hermite degree > 10 are not permitted'
   stop
end select

end function get_hermite_polynomial_o

function get_hermite_polynomial(degree,x) result(h)
implicit none

integer :: i
integer,intent(in) :: degree

real(8) :: h,a,c,yn,y1,y0
real(8),intent(in) :: x

a  = 2d0

y0 = 1d0
y1 = 2d0*x

if (degree == 0) then
   h = y0
else if (degree == 1) then
   h = y1
else if (degree >= 2) then
   do i = 2, degree
      c = 2d0*(i - 1d0)

      yn = a*x*y1 - c*y0

      y0 = y1
      y1 = yn
   end do

   h = yn
else if (degree < 0) then
   h = 0d0
end if

end function get_hermite_polynomial

function get_factorial(n) result(f)
implicit none

integer,intent(in) :: n

real(8) :: f
integer :: i

f = 1
do i = 1, n
   f = f*i
end do

end function get_factorial

function get_ho_wavefunc(i,alpha,dx) result(phi)
implicit none

real(8),parameter :: pisqrt=1.772453851d0

integer,intent(in) :: i

real(8),intent(in) :: alpha,dx
 
integer :: degree

real(8) :: phi
real(8) :: hermite_n,factorial,exponentl

degree = i - 1
!there is a dephasing, phi_i corresponds to H_(i-1), phi_1 has H_0, phi_2 has
!H_1, and so on
hermite_n = get_hermite_polynomial(degree,alpha*dx)
factorial = get_factorial(degree)
exponentl = exp(-1d0*(alpha*dx)**2/2d0)

phi = sqrt(alpha/(2d0**degree*factorial*pisqrt))*hermite_n*exponentl

end function get_ho_wavefunc

function get_secder_howavefunc(i,alpha,dx) result(d2p)
implicit none

real(8),parameter :: pisqrt=1.772453851d0

integer,intent(in) :: i

real(8),intent(in) :: alpha,dx

integer :: degree

real(8) :: d2p
real(8) :: a,f,g,df,dg,d2f,d2g,factorial

degree = i - 1
call get_gaussian_derivatives(alpha,dx,g,dg,d2g)
call get_hermite_derivatives(degree,alpha,alpha*dx,f,df,d2f)
factorial = get_factorial(degree)

a = sqrt(alpha/(2d0**degree*factorial*pisqrt))
d2p = a*f*d2g + 2d0*a*dg*df + a*g*d2f

end function get_secder_howavefunc

function integrate_t_phiphi(intpoints,lowlim,upplim,sb1,center1,sb2,center2,alpha) result(pp)
implicit none

integer,intent(in) :: intpoints,sb1,sb2

real(8),intent(in) :: center1,center2,lowlim,upplim,alpha

integer :: i

real(8) :: pp
real(8) :: h,q,phi1,phi2

h = (upplim - lowlim)/intpoints

q = lowlim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_ho_wavefunc(sb2,alpha,q-center2)

pp = phi1*phi2

q = upplim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_ho_wavefunc(sb2,alpha,q-center2)

pp = pp + phi1*phi2

pp = 0.5d0*pp

do i = 1, intpoints - 1
   q = lowlim + i*h
   phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
   phi2 = get_ho_wavefunc(sb2,alpha,q-center2)
   pp = pp + phi1*phi2
end do

pp = pp*h

end function integrate_t_phiphi

function integrate_t_pqcsqp(intpoints,lowlim,upplim,sb1,center1,sb2,center2,alpha,qe) result(pp)
implicit none

integer,intent(in) :: intpoints,sb1,sb2

real(8),intent(in) :: center1,center2,lowlim,upplim,alpha,qe

integer :: i

real(8) :: pp
real(8) :: h,q,phi1,phi2

h = (upplim - lowlim)/intpoints

q = lowlim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_ho_wavefunc(sb2,alpha,q-center2)

pp = phi1*phi2*(q - qe)**2

q = upplim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_ho_wavefunc(sb2,alpha,q-center2)

pp = pp + phi1*phi2*(q - qe)**2

pp = 0.5d0*pp

do i = 1, intpoints - 1
   q = lowlim + i*h
   phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
   phi2 = get_ho_wavefunc(sb2,alpha,q-center2)
   pp = pp + phi1*phi2*(q - qe)**2
end do

pp = pp*h

end function integrate_t_pqcsqp

function integrate_t_phid2p(intpoints,lowlim,upplim,sb1,center1,sb2,center2,alpha) result(pp)
!note no hbar is added here, considered =1 (one)
implicit none

integer,intent(in) :: intpoints,sb1,sb2

real(8),intent(in) :: center1,center2,alpha,lowlim,upplim

integer :: i

real(8) :: pp
real(8) :: h,q,phi1,phi2

h = (upplim - lowlim)/intpoints

q = lowlim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_secder_howavefunc(sb2,alpha,q-center2)
pp = phi1*phi2

q = upplim
phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
phi2 = get_secder_howavefunc(sb2,alpha,q-center2)

pp = pp + (phi1*phi2)

pp = 0.5d0*pp

do i = 1, intpoints - 1
   q = lowlim + i*h
   phi1 = get_ho_wavefunc(sb1,alpha,q-center1)
   phi2 = get_secder_howavefunc(sb2,alpha,q-center2)
   pp = pp + (phi1*phi2)
end do

pp = -pp*h

end function integrate_t_phid2p

end module m_vib
