program popmodel3
use ifport
use m_map
implicit none

character(len=2) :: c_nb,c_nt
character(len=9) :: fmt1,fmt2

complex(8) :: coeff,fact1,fact2,fact3
complex(8),dimension(:),allocatable :: pop,pop1,pop2,pop3

integer :: a,b,i,j,ng,nb,nd,basispc,nmap,cont
integer :: np,nosc,nmcs,nmds,seed_dimension,bath,init,mcs,it,is,ib
!integer,dimension(:),allocatable :: seed

real(8) :: delta,ome_max,dt,lumda_d,eg,eb,ed,mu,e0,beta,time_j,taw_j,omega_j,check,vomega
real(8) :: dt2,uj,qbeta,lambdacheck,a1,a2,et,gaussian,etotal,tn
real(8),dimension(:),allocatable :: ome,c2,kosc,x,p,fx,rm,pm,rn,pn,facn,popt
real(8),dimension(:,:),allocatable :: hm,lambda,popn,ug,ub,ud,hc
real(8),dimension(:,:),allocatable :: sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,hs,lld
real(8),dimension(:,:),allocatable :: llg,llb,llgb,llbg,lldb,llbd

call iniconc()

call srand(seed_dimension)

nmap = ng + nb + nd

allocate(ome(1:nosc),c2(1:nosc),kosc(1:nosc))
allocate(rm(1:nmap),pm(1:nmap))
allocate(rn(1:nmap),pn(1:nmap))
allocate(sgg(1:ng,1:ng),sgb(1:ng,1:nb),sgd(1:ng,1:nd))
allocate(sbg(1:nb,1:ng),sbb(1:nb,1:nb),sbd(1:nb,1:nd))
allocate(sdg(1:nd,1:ng),sdb(1:nd,1:nb),sdd(1:nd,1:nd))
allocate(hm(1:nmap,1:nmap),hc(1:nmap,1:nmap),hs(1:nmap,1:nmap))
allocate(ug(1:nmap,1:nmap),ub(1:nmap,1:nmap),ud(1:nmap,1:nmap))
allocate(lambda(1:nmap,1:nmap),llg(1:nmap,1:nmap),llb(1:nmap,1:nmap),lld(1:nmap,1:nmap))
allocate(llgb(1:nmap,1:nmap),llbg(1:nmap,1:nmap))
allocate(llbd(1:nmap,1:nmap),lldb(1:nmap,1:nmap))

call iniconq_d(nosc,lumda_d,ome_max,ome,c2,kosc)

allocate(popn(1:nmds+1,1:nmap),facn(1:nmap),popt(1:nmds+1))
allocate(pop(1:nmds+1),pop1(1:nmds+1),pop2(1:nmds+1),pop3(1:nmds+1))
allocate(x(1:nosc),p(1:nosc),fx(1:nosc))

if (ng > 9) then
   write(c_nb,'(i2)') ng
else
   write(c_nb,'(i1)') ng
end if

if (nmap > 9) then
   write(c_nt,'(i2)') nmap
else
   write(c_nt,'(i1)') nmap
end if

fmt1 = '('//trim(c_nb)//'f10.5)'
fmt2 = '('//trim(c_nt)//'f10.5)'

popn = 0d0
pop  = 0d0
pop1 = 0d0
pop2 = 0d0
pop3 = 0d0

dt  = 2d0*pi*dt
dt2 = 0.5d0*dt

!call get_lambda_eigenvectors(ng,nb,nd,eg,eb,ed,delta,vomega,&
!                              sgg,sgb,sgd,sbg,sbb,sbd,sdg,sdb,sdd,lambda,hs)

call get_preh(ng,nb,nd,eg,eb,ed,delta,vomega,hs)

lld = 0d0
do i = ng+nb+1, nmap
   lld(i,i) = 1d0
end do

MC: do mcs = 1, nmcs
   do is=1,nosc
      uj = 0.5d0*beta*dsqrt(kosc(is))
      
      qbeta = beta/(uj/tanh(uj))
      
      p(is) = gauss_noise2()/dsqrt(qbeta)
      x(is) = gauss_noise2()/dsqrt(qbeta*kosc(is))
      
      if(bath == 1) x(is)=x(is)+c2(is)/kosc(is)
   end do
   
   if (init == 3) then
      do i = 1, nmap
         rm(i) = gauss_noise2()/sqrt(2d0)
         rn(i) = gauss_noise2()/sqrt(2d0)
         pm(i) = gauss_noise2()/sqrt(2d0)
         pn(i) = gauss_noise2()/sqrt(2d0)
      end do
   else
      rm = 0d0
      rn = 0d0
      pm = 0d0
      pn = 0d0
   end if
  
   call get_coeff(ng,beta,vomega,rm,pm,rn,pn,coeff)

   call get_force_fb(nmap,ng,nb,lld,kosc,x,c2,rm,pm,rn,pn,fx)

   ib = 1

!   do i = 1, nmap
!      facn(i)    = coeff*(rm(i)**2 + pm(i)**2 - 1)
!      popn(ib,i) = popn(ib,i) + facn(i)
!      popt(ib)   = popt(ib) + facn(i)
!   end do

   call get_facts_pop(nmap,ng,nb,coeff,rm,pm,rn,pn,fact1,fact2,fact3)

   pop(ib)  = pop(ib)  + (fact1+fact2+fact3)
   pop1(ib) = pop1(ib) + (fact1)
   pop2(ib) = pop2(ib) + (fact2)
   pop3(ib) = pop3(ib) + (fact3)
   
   a1 = 0d0
   a2 = 0d0
   do is=1,nosc
      a1 = a1 + 2.d0*c2(is)**2/ome(is)**2
      a2 = a2 + 2.d0*c2(is)*x(is)
   end do
   
   open(747,file='etotal.log')
   
   MD: do it = 1, nmds
      gaussian=sqrt(4.d0*log(2.d0)/(pi*taw_j**2))*exp(-4.d0*log(2.d0)*((it-0.5d0)*dt-time_j)**2/(taw_j**2))
      et = gaussian*e0*cos(omega_j*((it-0.5d0)*dt-time_j))
   
      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      
!      call get_hm(nmap,ng,nb,lmd,basispc,delta,mu,et,a1,a2,kg,kb,kd,vg,vb,vd,hm)
      call get_hm2(nmap,ng,nb,mu,et,a1,a2,hs,hm)
      !call make_hm_traceless(nmap,hm,tn)
!if (it == 500) then
!write(*,*) 'mapping hamiltonian'
!write(*,fmt2) hm
!stop
!end if
      call evolve_pm(nmap,dt2,hm,rm,pm)
      call evolve_pm(nmap,dt2,hm,rn,pn)

      do is = 1, nosc
         x(is) = x(is) + dt*p(is)
      end do

      a2=0.d0
      do is = 1, nosc
          a2 = a2 + 2.d0*c2(is)*x(is)
      end do

!      call update_hm(nmap,ng,nb,lmd,basispc,delta,mu,et,a1,a2,kg,kb,kd,vg,vb,vd,hm)
!      call update_hm2(nmap,ng,nb,delta,mu,et,a1,a2,hc,hm)
      call get_hm2(nmap,ng,nb,mu,et,a1,a2,hs,hm)
      !call make_hm_traceless(nmap,hm,tn)

      call evolve_rm(nmap,dt,hm,pm,rm)
      call evolve_rm(nmap,dt,hm,pn,rn)

      call evolve_pm(nmap,dt2,hm,rm,pm)
      call evolve_pm(nmap,dt2,hm,rn,pn)

      call get_force_fb(nmap,ng,nb,lld,kosc,x,c2,rm,pm,rn,pn,fx)
      
      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      
      ib = it + 1
      
!      do i = 1, nmap
!         facn(i)    = coeff*(rm(i)**2 + pm(i)**2 - 1d0)
!         popn(ib,i) = popn(ib,i) + facn(i)
!         popt(ib)   = popt(ib) + facn(i)
!      end do

      call get_facts_pop(nmap,ng,nb,coeff,rm,pm,rn,pn,fact1,fact2,fact3)
      
      pop(ib)  = pop(ib)  + (fact1+fact2+fact3)
      pop1(ib) = pop1(ib) + (fact1)
      pop2(ib) = pop2(ib) + (fact2)
      pop3(ib) = pop3(ib) + (fact3)
      
      !if (mod(mcs,1000) == 0) then
      !   call get_totalenergy_traceless(nmap,hm,tn,pm,rm,x,p,kosc,etotal)
      !   write(747,*) it, etotal
      !end if
   end do MD

   close(747)

   if (mod(mcs,1000) == 0) then
      open(444,file='temp.out')
!      do i = 1, ng
!         do j = 1, nmds+1
!            pop1(j) = pop1(j) + popn(j,i)
!         end do
!      end do
!
!      do i = ng+1, ng+nb
!         do j = 1, nmds+1
!            pop2(j) = pop2(j) + popn(j,i)
!         end do
!      end do
!
!      do i = ng+nb+1, ng+nb+nd
!         do j = 1, nmds+1
!            pop3(j) = pop3(j) + popn(j,i)
!         end do
!      end do
!      
!      do i = 1, nmds+1
!         pop(i) = pop1(i) + pop2(i) + pop3(i)
!      end do

      do i = 1, nmds+1
         write(444,'(i10,4f20.9)') i-1,dble(pop1(i)/pop(i)),dble(pop2(i)/pop(i)),dble(pop3(i)/pop(i)),dble(pop(i))
      end do
      close(444)

   end if
end do MC

do ib = 1, nmds+1
   pop1(ib) = pop1(ib)/pop(ib)
   pop2(ib) = pop2(ib)/pop(ib)
   pop3(ib) = pop3(ib)/pop(ib)
   write(333,'(i10,4f20.9)') ib-1, dble(pop1(ib)),dble(pop2(ib)),dble(pop3(ib)),dble(pop(ib))!/dnmcs
end do

deallocate(ome,c2,kosc)
deallocate(pop,pop1,pop2,pop3)
deallocate(x,p)

contains

subroutine iniconc()
implicit none

!integer :: i

open (666,file='md.in')
read(666,*)
read(666,*) np,delta,nosc,ome_max
read(666,*)
read(666,*) nmcs,nmds,seed_dimension,dt,lumda_d
read(666,*)
read(666,*) eg,eb,ed,mu,e0,beta,vomega
read(666,*)
read(666,*) time_j,taw_j,omega_j 
read(666,*)
read(666,*) bath,init
read(666,*)
read(666,*) ng,nb,nd
close(666)

!call random_seed(size=seed_dimension)
!allocate (seed(seed_dimension))
!do i=1,seed_dimension
!  seed(i) = 3*2**i-1
!enddo
!call random_seed(put=seed)

end subroutine iniconc

end program popmodel3
