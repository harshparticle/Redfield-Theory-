Module mod_rt
implicit none
real*8, parameter :: clight=2.99792458D10, av=6.0221367D23, hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0, kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
complex*16,parameter :: iota=(0.d0,1.d0)
real*8 pi,wave_to_J
integer,parameter::nw=5000
!!Defining system parameters

integer nquant
real*8, allocatable:: H_diab(:,:),E_exc(:),lambda_diab(:,:),n_be(:,:)
real*8, allocatable:: c_tr(:,:),corr(:,:,:),omg(:,:),spectral(:,:),R(:,:)
complex*16, allocatable::sigma(:,:),sigma_diab(:,:),sigma_exc(:,:),RR(:,:),rate(:,:),tan_calc(:,:)

real*8 gamma,temperature,lambda
integer nsteps
real*8 dt,tot_time,curr_time

integer nold
integer,allocatable :: seed(:)
real*8, allocatable :: work(:)
complex*16, allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)
!Definining robustness parameter
real*8 p2collapse,ita,p3collapse
real*8,dimension(10000)::robustness_parameter

!Defining simulation parameters
real*8 r1,x1,r2,x2,r3,x3,r4,x4,H_12,r5,x5,H_13,r6,x6,H_23,r7,x7,lamb
 
contains
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine setup
   implicit none
   character st_ch
   integer i,size_seed,seed2(2)

   pi=dacos(-1.d0)
   wave_to_J=2*pi*clight*hbar

   open(10,file="rt.inp")
   read(10,*) nquant
   read(10,*) tot_time
   read(10,*) dt
   read(10,*) seed2
   read(10,*) st_ch
   close(10)

   !-------------------------------------------------------------------------------------

   if(st_ch.ne. 'x') then
     write(6,*) "problem in reading input file"
     stop
   end if
   
   !-------------------------------------------------------------------------------------
   
   allocate(H_diab(nquant,nquant),E_exc(nquant),lambda_diab(nquant,nquant),n_be(nquant,nquant))
   allocate(c_tr(nquant,nquant),omg(nquant,nquant),R(nquant,nquant),RR(nquant,nquant),rate(nquant,nquant))
   allocate(corr(nquant,nquant,nw))
   allocate(sigma(nquant,nquant),sigma_diab(nquant,nquant),sigma_exc(nquant,nquant),tan_calc(nquant,nquant))

   nsteps=int(tot_time/dt)

   !--------------------------------------------------------------------------------------

   call random_seed(size=size_seed)
   allocate(seed(size_seed))
   do i=1,size_seed/2
      seed(i)=seed2(1)*(2*i+i*i-7)
   enddo
   do i=size_seed/2+1,size_seed
      seed(i)=seed2(2)*(i/2+34-i**3)
   enddo
   call random_seed(put=seed)
   !-------------------w------------------------------------------------------------------

end subroutine setup
!---------- ----------------------------------------------------------------------------

subroutine main
   implicit none
   integer i,j,k
 
   open(101,file="rob_parameter_10.out")
   open(102,file="eigen_energies_10.out")
   open(103,file="H_12_10.out")
   open(104,file="H_23_10.out")
   open(105,file="H_13_10.out")
   open(106,file="lambda_10.out")
   call setup_parameters
   do k=1,10000
     call changing_parameters
     call init_cond
     open(100,file="sigma_10.out")
     do i=1,nsteps
        call write_sigma
        call evolve
     enddo
     robustness_parameter(k)=(nsteps*real(sigma_diab(3,3)))/(p2collapse)
     write(101,*) robustness_parameter(k),maxval(robustness_parameter),maxloc(robustness_parameter)
     close(100)
   enddo
   close(101)
   close(102)
   close(103)
   close(104)
   close(105)
   close(106)
end subroutine main

!--------------------------------------------------------------------------------------------


subroutine setup_parameters
    implicit none
    integer i

    gamma=1/25.d-15
    temperature=77.d0
    lambda=35.00*wave_to_J/pi

    H_diab(1,1)=510.d0*wave_to_J
    H_diab(2,2)=390.d0*wave_to_J
    H_diab(3,3)=310.d0*wave_to_J
    H_diab(1,2)=50.d0*wave_to_J  ;  H_diab(2,1)=H_diab(1,2)
    H_diab(1,3)=0.d0*wave_to_J    ;  H_diab(3,1)=H_diab(1,3)
    H_diab(2,3)=30.8d0*wave_to_J   ;  H_diab(3,2)=H_diab(3,2)
  

    do i=1,nquant
       lambda_diab(i,i)=lambda*pi
    enddo
end subroutine setup_parameters    
!-----------------------------------------------------------------------------------------------
subroutine changing_parameters
    implicit none
    real*8 en(nquant),vect(nquant,nquant)
    real*8,dimension(3,3)::Hamil_e

    call gaussian_random_number(r1)
    x1=10*r1
    H_diab(1,1)=H_diab(1,1)+x1*wave_to_J

    call gaussian_random_number(r2)
    x2=10*r2
    H_diab(2,2)=H_diab(2,2)+x2*wave_to_J

    call gaussian_random_number(r3)
    x3=10*r3
    H_diab(3,3)=H_diab(3,3)+x3*wave_to_J
    

    call gaussian_random_number(r4)
    x4=5*r4
    H_12=H_diab(1,2)+x4*wave_to_J
    if(H_12/wave_to_J<100 .and. H_12/wave_to_J>-100) then
      H_diab(1,2)=H_12
      H_diab(2,1)=H_12
      write(103,*) H_diab(1,2)/wave_to_J,x4
    endif
    
    call gaussian_random_number(r6)
    x6=5*r6
    H_23=H_diab(2,3)+x6*wave_to_J
    if(H_23/wave_to_J<100 .and. H_23/wave_to_J>-100) then
      H_diab(2,3)=H_23
      H_diab(3,2)=H_23
      write(104,*) H_diab(2,3)/wave_to_J,x6
    endif
   
    call gaussian_random_number(r5)
    x5=5*r5
    H_13=H_diab(1,3)+x5*wave_to_J
    if(H_13/wave_to_J<20 .and. H_13/wave_to_J>-20) then
      H_diab(1,3)=H_13
      H_diab(3,1)=H_13
      write(105,*) H_diab(1,3)/wave_to_J,x5
    endif 

    call gaussian_random_number(r7)
    x7=5*r7
    lamb=lambda+x7*wave_to_J/pi
    if(lamb/wave_to_J*pi>0) then
      lambda=lamb
    endif
    write(106,*) lambda/wave_to_J*pi   
 
    write(102,*) H_diab(1,1)/wave_to_J,H_diab(2,2)/wave_to_J,H_diab(3,3)/wave_to_J
     Hamil_e=H_diab

    call diag(Hamil_e,nquant,en,vect,nquant)
    c_tr=vect
    E_exc=en
    
    


end subroutine changing_parameters
!------------------------------------------------------------------------------------------

subroutine init_cond
   implicit none
  
   p2collapse=0.d0
   sigma_diab=0.d0
   sigma_diab(1,1)=1.d0

   sigma=matmul(transpose(c_tr),matmul(sigma_diab,c_tr))

   curr_time=0.d0
  

end subroutine init_cond
!-----------------------------------------------------------------------------------------


subroutine evolve
    implicit none
    complex*16 sigma_dot(nquant,nquant)

    call calculate_sigmadot(sigma_dot)
   
    p2collapse=p2collapse+(real(sigma_diab(2,2))*dt*1.d15)
    
    
    sigma=sigma+sigma_dot*dt
    curr_time=curr_time+dt
   

end subroutine evolve
!---------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------     


!------------------------------------------------------------------------------------------
subroutine calculate_sigmadot(sigma_dot)
    implicit none
    complex*16,intent(out)::sigma_dot(nquant,nquant)
    integer i,j
    
    call calculate_rate_exc
   
    do i=1,nquant
      do j=1,nquant
         if(i.ne.j) then
           sigma_dot(i,i)=sigma_dot(i,i)+((RR(j,i)*sigma(j,j))-(RR(i,j)*sigma(i,i)))
         endif
      enddo
    enddo

end subroutine calculate_sigmadot
!----------------------------------------------------------------------------------------
subroutine calculate_rate_exc
   implicit none
   integer i,j,n
   
   RR=0.d0
   call calculate_rate
 
   do i=1,nquant
      do j=1,nquant
         if(i.ne.j) then
           do n=1,nquant
              RR(i,j)=RR(i,j)+c_tr(n,i)*c_tr(n,j)*c_tr(n,j)*c_tr(n,i)*R(i,j)
           enddo
         endif
      enddo
   enddo
        
end subroutine calculate_rate_exc
!----------------------------------------------------------------------------------------

subroutine calculate_rate
    implicit none
    integer i,j

    call calculate_freq

    do i=1,nquant
       do j=1,nquant
          R(i,j)=((((2*lambda*gamma*omg(i,j))/((omg(i,j))**2+gamma**2))*((1/(tanh((hbar*omg(i,j))&
          &/(2*kb*temperature))))+1))/hbar)
       enddo
    enddo

end subroutine calculate_rate

!-----------------------------------------------------------------------------------------
subroutine calculate_freq
      implicit none
      integer i,j
      do i=1,nquant
         do j=1,nquant
            omg(i,j)=((E_exc(i)-E_exc(j))/hbar)
         enddo
      enddo

end subroutine calculate_freq   

!------------------------------------------------------------------------------------------
subroutine write_sigma

   implicit none
   integer i,j
   

   
   do i=1,nquant
      do j=1,nquant
         sigma_exc(i,j)=sigma(i,j)
         if(i.ne.j) then
            sigma_exc(i,j)=0.d0
         endif
      enddo
   enddo
   sigma_diab=matmul(c_tr,matmul(sigma_exc,transpose(c_tr)))

   write(100,*) curr_time*1.d15,real(sigma_diab(1,1)),real(sigma_diab(2,2)),real(sigma_diab(3,3))
  
     
end subroutine write_sigma

!--------------------------------------------------------------------------------------

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
    implicit none
    integer,intent(in) :: n,m_values
    real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
    real*8,intent(inout) :: mat(n,n)
    real*8 vl,vu,abstol
    integer il,iu,info,m,AllocateStatus
    integer lwork,liwork

    vl=0.d0;vu=0.d0
    il=1;iu=m_values
    abstol=0.d0
    info=0

    if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or..not.allocated(isuppz)) then
      lwork=-1;liwork=-1
      if(allocated(isuppz))deallocate(isuppz)
      if(allocated(work))deallocate(work)
      if(allocated(iwork))deallocate(iwork)
      allocate(isuppz(2*m_values),work(n),iwork(n))
      call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
      lwork=nint(work(1));liwork=iwork(1)
      deallocate(work,iwork)
      allocate(work(lwork),STAT=AllocateStatus)
      if(allocatestatus.ne.0) write(6,*) "problem in diag, allocation"
      allocate(iwork(liwork),STAT=AllocateStatus)
      if(allocatestatus.ne.0) write(6,*)"problem in diag,allocation"
      nold=n
    endif

    lwork=size(work)
    liwork=size(iwork)

    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    if(info.ne.0) then
      write(6,*) "problem in diagonalisation",info
      stop
    endif

end subroutine diag

!-------------------------------------------------------------------------------------------
subroutine gaussian_random_number(rnd)
  !!generates gaussian distribution with centre 0, sigma 1
  !! q0+sig*rnd gives centre=q0, sigma=sig
  implicit none
  real*8, intent(out)::rnd
  real*8 rnd1,rnd2,pi
  
  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd=dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!--------------------------------------------------------------------------
End Module mod_rt   


 

   


   
