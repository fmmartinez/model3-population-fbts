module m_map
use ifport
implicit none

real(8),parameter :: pi=3.1415926535d0

contains

subroutine get_force_traceless(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)
implicit none

real(8),dimension(:),intent(in) :: rm,pm,x
real(8),dimension(:),intent(out) :: f

integer :: a,b,i,j,n
integer,intent(in) :: nmap,ng,nb

real(8) :: trace,tn
real(8),dimension(:),intent(in) :: kosc,c2
real(8),dimension(:,:),intent(in) :: lld
real(8),dimension(:,:),allocatable :: dh

allocate(dh(1:nmap,1:nmap))

n = size(x)

f = 0d0
do j = 1, n
   f(j) = -kosc(j)*x(j)
   
   dh = (lld)*c2(j)
   
   trace = 0d0
   do a = 1, nmap
      trace = trace + dh(a,a)
   end do

   tn = trace/nmap
   
   f(j) = f(j) + tn

   do a = 1, nmap
      do b = 1, nmap
         if (a == b) then
            f(j) = f(j) - (dh(a,b) - tn)*(rm(a)*rm(b) + pm(a)*pm(b))
         else
            f(j) = f(j) - dh(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
         end if
      end do
   end do
end do

end subroutine get_force_traceless

subroutine make_hm_traceless(nmap,hm,tn)
implicit none

real(8),dimension(:,:),intent(inout) :: hm

integer :: i
integer,intent(in) :: nmap

real(8) :: trace
real(8),intent(out) :: tn

trace = 0d0
do i = 1, nmap
   trace = trace + hm(i,i)
end do

trace = trace/nmap

do i = 1, nmap
   hm(i,i) = hm(i,i) - trace
end do

tn = trace

end subroutine make_hm_traceless

subroutine get_totalenergy(nmap,hm,pm,rm,x,p,kosc,etotal)
implicit none

integer :: i,j,n
integer,intent(in) :: nmap

real(8),intent(out) :: etotal
real(8),dimension(:),intent(in) :: x,p,kosc,rm,pm
real(8),dimension(:,:),intent(in) :: hm

n = size(x)

etotal = 0d0
!classical part
do i = 1, n
   etotal = etotal + 0.5d0*(p(i)**2 + kosc(i)*x(i)**2)
end do
!mapping part
do i = 1, nmap
   do j = 1, nmap
      if (i == j) then
         etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j) - 1d0))
      else
         etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j)))
      end if
   end do
end do

end subroutine get_totalenergy

subroutine get_totalenergy_traceless(nmap,hm,tn,pm,rm,x,p,kosc,etotal)
implicit none

integer :: i,j,n
integer,intent(in) :: nmap

real(8),intent(in) :: tn
real(8),intent(out) :: etotal
real(8),dimension(:),intent(in) :: x,p,kosc,rm,pm
real(8),dimension(:,:),intent(in) :: hm

n = size(x)

etotal = 0d0
!classical part
do i = 1, n
   etotal = etotal + 0.5d0*(p(i)**2 + kosc(i)*x(i)**2)
end do
!mapping part
etotal = etotal + tn
do i = 1, nmap
   do j = 1, nmap
      etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j)))
   end do
end do

end subroutine get_totalenergy_traceless

subroutine get_coeff(ng,beta,omega,rm,pm,rn,pn,coeff)
implicit none

complex(8),intent(out) :: coeff

integer :: i
integer,intent(in) :: ng

real(8) :: z
real(8),intent(in) :: beta,omega
real(8),dimension(:),intent(in) :: rm,pm,rn,pn
real(8),dimension(:),allocatable :: exp_be,prob

allocate(exp_be(1:ng),prob(1:ng))

exp_be = 0d0
z = 0d0
do i = 1, ng
   !exp_be(i) = exp(-beta*omega*(i - 0.5d0))
   !exp_be(i) = exp(-1.44d0*(i - 0.5d0))
   exp_be(i) = exp(-2.29d-1*(i - 0.5d0))
   !z = z + exp_be(i)
end do

prob = exp_be!/z

coeff = 0d0
!only diagonal because in the ngxng matrix by construction their lambdas are
! (1,0,0...), (0,1,0...), (0,0,1...), etc
do i = 1, ng
   coeff = coeff + cmplx(rm(i),pm(i))*cmplx(rn(i),-pn(i))*prob(i)
end do

deallocate(exp_be)
deallocate(prob)
end subroutine get_coeff

subroutine evolve_rm(nmap,dt,hm,pm,rm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt
real(8),dimension(:),intent(in) :: pm
real(8),dimension(:),intent(inout) :: rm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      rm(i) = rm(i) + dt*hm(i,j)*pm(j)
   end do
end do

end subroutine evolve_rm

subroutine evolve_pm(nmap,dt2,hm,rm,pm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt2
real(8),dimension(:),intent(in) :: rm
real(8),dimension(:),intent(inout) :: pm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      pm(i) = pm(i) - dt2*hm(i,j)*rm(j)
   end do
end do

end subroutine evolve_pm

subroutine get_facts_pop_traceless(nmap,coeff,llg,llb,lld,rm,pm,fact1,fact2,fact3)
implicit none

integer :: a,b,i
integer,intent(in) :: nmap

real(8) :: trt
real(8),intent(in) :: coeff
real(8),intent(out) :: fact1,fact2,fact3
real(8),dimension(:),intent(in) :: rm,pm
real(8),dimension(:,:),intent(in) :: llg,llb,lld

trt = 0d0
do i = 1, nmap
   trt = trt + llg(i,i)
end do
fact1 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact1 = fact1 + 0.5d0*(llg(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact1 = fact1 + 0.5d0*(llg(a,b))*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact1 = coeff*fact1

trt = 0d0
do i = 1, nmap
   trt = trt + llb(i,i)
end do
fact2 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact2 = fact2 + 0.5d0*(llb(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact2 = fact2 + 0.5d0*llb(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact2 = coeff*fact2

trt = 0d0
do i = 1, nmap
   trt = trt + lld(i,i)
end do
fact3 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact3 = fact3 + 0.5d0*(lld(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact3 = fact3 + 0.5d0*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact3 = coeff*fact3
end subroutine get_facts_pop_traceless

subroutine get_facts_pop(nmap,ng,nb,coeff,rm,pm,fact1,fact2,fact3)
implicit none

integer :: a,b
integer,intent(in) :: nmap,ng,nb

real(8),intent(in) :: coeff
real(8),intent(out) :: fact1,fact2,fact3
real(8),dimension(:),intent(in) :: rm,pm

fact1 = 0d0
do a = 1, ng
   fact1 = fact1 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact1 = coeff*fact1

fact2 = 0d0
do a = ng+1, ng+nb
   fact2 = fact2 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact2 = coeff*fact2

fact3 = 0d0
do a = ng+nb+1, nmap
   fact3 = fact3 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact3 = coeff*fact3
end subroutine get_facts_pop

function kronecker_delta(i,j) result (d)
implicit none

integer,intent(in) :: i, j

real(8) :: d

if (i == j) then
   d = 1d0
else
   d = 0d0
end if

end function kronecker_delta

subroutine get_force(nmap,ng,nb,lld,kosc,x,c2,rm,pm,f)
implicit none

integer :: a,b,i,j,n
integer,intent(in) :: nmap,ng,nb

real(8),dimension(:),intent(in) :: kosc,x,c2,rm,pm
real(8),dimension(:),intent(out) :: f
real(8),dimension(:,:),intent(in) :: lld

n = size(x)

f = 0d0
do j = 1, n
   f(j) = -kosc(j)*x(j)
   !exclude lambdas from 1 to ng, because I know those won't contribute
   do a = ng+1, nmap
      do b = ng+1, nmap
         if (a == b) then
            f(j) = f(j) - c2(j)*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b) - 1d0)
         else
            f(j) = f(j) - c2(j)*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
         end if
      end do
   end do
end do

end subroutine get_force

subroutine get_hm2(nmap,ng,nb,mu,et,a1,a2,hs,hm)
implicit none

integer :: i
integer,intent(in) :: nmap,ng,nb

real(8),intent(in) :: mu,et,a1,a2
real(8),dimension(:,:),intent(in) :: hs
real(8),dimension(:,:),intent(out) :: hm

hm = hs
hm(1:ng,ng+1:ng+nb) = hs(1:ng,ng+1:ng+nb)*(-mu*et)
hm(ng+1:ng+nb,1:ng) = hs(ng+1:ng+nb,1:ng)*(-mu*et)
do i = ng+nb+1,nmap
   hm(i,i) = hm(i,i) + a1 + a2
end do
end subroutine get_hm2

subroutine get_lambda_eigenvectors(ng,nb,nd,eg,eb,ed,delta,omega,&
                                    sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,lambda,hsf)
use m_vib
implicit none

integer,parameter :: ip = 10000

character(len=2) :: c_nb,c_nt
character(len=9) :: fmt1,fmt2

integer :: i, j, nt
integer,intent(in) :: nd,ng,nb

real(8) :: omgsq,cg,cb,cd,uint,lint,alpha
real(8),intent(in) :: eg,eb,ed,delta,omega
real(8),dimension(1:3,1:3) :: hs
real(8),dimension(:),allocatable :: eigval
real(8),dimension(:,:),intent(out) :: lambda,sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,hsf
real(8),dimension(:,:),allocatable :: kg,kb,kd,vg,vb,vd,hxp,hxp_temp,sxp

!lapack variables
integer,parameter :: lwork = 30000
real(8),dimension(1:lwork) :: work
integer :: info

!generating some constants
nt = ng + nb + nd

allocate(kg(1:ng,1:ng),kb(1:nb,1:nb),kd(1:nd,1:nd))
allocate(vg(1:ng,1:ng),vb(1:nb,1:nb),vd(1:nd,1:nd))
allocate(hxp(1:nt,1:nt),hxp_temp(1:nt,1:nt),sxp(1:nt,1:nt))
allocate(eigval(1:nt))

omgsq = omega**2

cg = 0d0
cb = 2d0*sqrt(10d0)/omega
cd = cb/2d0

uint = cb + 5.5d0
lint = cg - 5.5d0

alpha = sqrt(omega)

!original subsystem hamiltonian
hs = 0d0
hs(1,1) = eg
hs(2,2) = eb
!hs(2,3) = delta
!hs(3,2) = delta
hs(3,3) = ed

!initialization
kg = 0d0
vg = 0d0
sgg = 0d0

sgb = 0d0

sgd = 0d0
!
sbg = 0d0

kb = 0d0
vb = 0d0
sbb = 0d0

sbd = 0d0
!
sdg = 0d0

sdb = 0d0

kd = 0d0
vd = 0d0
sdd = 0d0

!construction of integrals for the extended hamiltonian
! 1,1
do i = 1, ng
   do j = 1, ng
      kg(i,j) = integrate_t_phid2p(ip,lint,uint,i,cg,j,cg,alpha)
      vg(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cg,j,cg,alpha,cg)
   end do
   sgg(i,i) = 1d0
end do
! 1,2
do i = 1, ng
   do j = 1, nb
      sgb(i,j) = integrate_t_phiphi(ip,lint,uint,i,cg,j,cb,alpha)
   end do
end do
! 1,3
! stays at 0, will not be used
! 2,1
do i = 1, nb
   do j = 1, ng
      sbg(i,j) = integrate_t_phiphi(ip,lint,uint,i,cb,j,cg,alpha)
   end do
end do
! 2,2
do i = 1, nb
   do j = 1, nb
      kb(i,j) = integrate_t_phid2p(ip,lint,uint,i,cb,j,cb,alpha)
      vb(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cb,j,cb,alpha,cb)
   end do
   sbb(i,i) = 1d0
end do
! 2,3
do i = 1, nb
   do j = 1, nd
      sbd(i,j) = integrate_t_phiphi(ip,lint,uint,i,cb,j,cd,alpha)
   end do
end do
! 3,1
! stays at 0, will not be used
! 3,2
do i = 1, nd
   do j = 1, nb
      sdb(i,j) = integrate_t_phiphi(ip,lint,uint,i,cd,j,cb,alpha)
   end do
end do
! 3,3
do i = 1, nd
   do j = 1, nd
      kd(i,j) = integrate_t_phid2p(ip,lint,uint,i,cd,j,cd,alpha)
      vd(i,j) = integrate_t_pqcsqp(ip,lint,uint,i,cd,j,cd,alpha,cd)
   end do
   sdd(i,i) = 1d0
end do

!expansion of subsystem hamiltonian with basis functions
hxp(1:ng,1:ng)       = hs(1,1)*sgg(1:ng,1:ng) + 0.5d0*(kg(1:ng,1:ng) + omgsq*vg(1:ng,1:ng))
hxp(1:ng,ng+1:ng+nb) = hs(1,2)*sgb(1:ng,1:nb)
hxp(1:ng,ng+nb+1:nt) = hs(1,3)*sgd(1:ng,1:nd)

hxp(ng+1:ng+nb,1:ng)       = hs(2,1)*sbg(1:nb,1:ng)
hxp(ng+1:ng+nb,ng+1:ng+nb) = hs(2,2)*sbb(1:nb,1:nb) + 0.5d0*(kb(1:nb,1:nb) + omgsq*vb(1:nb,1:nb))
hxp(ng+1:ng+nb,ng+nb+1:nt) = hs(2,3)*sbd(1:nb,1:nd)

hxp(ng+nb+1:nt,1:ng)       = hs(3,1)*sdg(1:nd,1:ng)
hxp(ng+nb+1:nt,ng+1:ng+nb) = hs(3,2)*sdb(1:nd,1:nb)
hxp(ng+nb+1:nt,ng+nb+1:nt) = hs(3,3)*sdd(1:nd,1:nd) + 0.5d0*(kd(1:nd,1:nd) + omgsq*vd(1:nd,1:nd))

!overlap?
sxp(1:nt,1:nt) = 0d0
do i = 1, nt
   sxp(i,i) = 1d0
end do

!if (nb > 9) then
!   write(c_nb,'(i2)') nb
!else
!   write(c_nb,'(i1)') nb
!end if
!
!if (nt > 9) then
!   write(c_nt,'(i2)') nt
!else
!   write(c_nt,'(i1)') nt
!end if
!
!fmt1 = '('//trim(c_nb)//'f10.5)'
!fmt2 = '('//trim(c_nt)//'f10.5)'

!write(*,*) 'hxp'
!write(*,fmt2) hxp

!write(*,*) 'sbg'
!write(*,fmt1) sbg

!Diagonalization
hxp_temp = hxp
!call dsyev('v','u',nt,hxp_temp,nt,eigval,work,lwork,info)
call dsygv(1,'v','u',nt,hxp_temp,nt,sxp,nt,eigval,work,lwork,info)
if (info /=0) then
   print *, info, 'info from Diagonalization Hc=ec'
   stop
end if

!write(*,*) 'diagonalization result'
!write(*,*) 'eigenvectors'
!do i = 1, nt
!   write(*,*) i, 'E=',eigval(i)
!   write(*,fmt2) hxp_temp(1:nt,i)
!end do

lambda(1:nt,1:nt) = hxp_temp(1:nt,1:nt)

hsf = 0d0

hsf = matmul(matmul(transpose(lambda),hxp),lambda)

!write(*,*) 'hsf'
!write(*,fmt2) hsf

deallocate(hxp,hxp_temp,eigval)

end subroutine get_lambda_eigenvectors

subroutine get_preh(ng,nb,nd,eg,eb,ed,delta,omega,hs)
use m_vib
implicit none

integer,parameter :: ip = 10000

character(len=2) :: c_nt
character(len=9) :: fmt1

integer :: i,j,nt
integer,intent(in) :: ng,nb,nd

real(8) :: cg,cb,cd,uint,lint,alpha
real(8),intent(in) :: eg,eb,ed,delta,omega
real(8),dimension(:,:),intent(out) :: hs

nt = ng + nb + nd

cg = 0d0
cb = 2d0*sqrt(10d0)/omega
cd = cb/2d0

uint = cb + 6d0
lint = cb - 6d0

alpha = sqrt(omega)

hs = 0d0

!fill g|g
do i = 1, ng
   hs(i,i) = eg + (i - 0.5d0)*omega
end do
!fill g|b
do i = 1, ng
   do j = ng+1, ng+nb
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i,cg,j-ng,cb,alpha)
   end do
end do
!fill g|d
!not necessary
!fill b|g
do i = ng+1, ng+nb
   do j = 1, ng
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i-ng,cb,j,cg,alpha)
   end do
end do
!fill b|b
do i = ng+1, ng+nb
   hs(i,i) = eb + (i -ng - 0.5d0)*omega
end do
!fill b|d
do i = ng+1, ng+nb
   do j = ng+nb+1, nt
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng,cb,j-ng-nb,cd,alpha)
   end do
end do
!fill d|g
!not necessary
!fill d|b
do i = ng+nb+1, nt
   do j = ng+1, ng+nb
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng-nb,cd,j-ng,cb,alpha)
   end do
end do
!fill d|d
do i = ng+nb+1, nt
   hs(i,i) = ed + (i -ng -nb - 0.5d0)*omega
end do

if (nt > 9) then
   write(c_nt,'(i2)') nt
else
   write(c_nt,'(i1)') nt
end if

fmt1 = '('//trim(c_nt)//'f10.5)'

!print fmt1, hs

end subroutine get_preh


subroutine iniconq_d(nosc,lumda_d,ome_max,ome,c2,kosc)
implicit none

integer :: i
integer,intent(in) :: nosc

real(8) :: check
real(8),intent(in) :: ome_max,lumda_d
real(8),dimension(:),intent(inout) :: ome,c2,kosc

check=0
do i=1,nosc
!  ome(i)=tan(pi*i/(4*nosc))
   ome(i)=tan(i*datan(ome_max)/nosc)
!  ome(i)=i**2*ome_max/nosc**2
!  rho(i)=sqrt(ome(i)*ome_max)/(1+ome(i)**2)
!  c2(i)=0.5d0*ome(i)*dsqrt(lumda_d/nosc)
   c2(i)=ome(i)*dsqrt(datan(ome_max)*lumda_d/(pi*nosc))
!  c2(i)=ome(i)*sqrt(lumda_d*rho(i)*2/(pi*nosc))
   check=check+c2(i)**2/ome(i)**2
   kosc(i)=ome(i)**2 
end do
write(6,*) check

end subroutine iniconq_d

function gauss_noise2() result(g)
implicit none

real(8),parameter :: pi2=2.0*3.141592654

real(8) :: g,z1,z2

!call random_number(z1)
!call random_number(z2)
z1 = rand()
z2 = rand()
g = sqrt(-2.d0*log(z1))*cos(pi2*z2)

end function gauss_noise2

end module m_map
