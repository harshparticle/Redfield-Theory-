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
complex*16, allocatable::sigma(:,:),sigma_diab(:,:),RR(:,:),rate(:,:)

real*8 gamma,temperature,lambda
integer nsteps
real*8 dt,tot_time,curr_time

integer nold
integer,allocatable :: seed(:)
real*8, allocatable :: work(:)
complex*16, allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)


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
   allocate(sigma(nquant,nquant),sigma_diab(nquant,nquant))

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
   !-------------------------------------------------------------------------------------

end subroutine setup
!--------------------------------------------------------------------------------------

subroutine main
   implicit none
   integer i,j
  
   open(101,file="rate.out")
   open(102,file="rate_exc.out")
   open(103,file="energies.out")
   open(104,file="total_rate.out")
   open(105,file="c_tr.out")
   open(106,file="R.out")
   open(107,file="population.out")
   open(108,file="sigma_dot.out")
   open(109,file="tensor.out")
   call setup_parameters
   call init_cond
   open(100,file="sigma.out")
   do i=1,nsteps
      call write_sigma
      call evolve
      rate=rate+RR
   enddo
   write(104,*) rate 
   close(100)
   close(101)
   close(102)
   close(103)
   close(104)
   close(105)
   close(106)
   close(107)
   close(108)
   close(109)
end subroutine main

!--------------------------------------------------------------------------------------------


subroutine setup_parameters
    implicit none
    integer i
    real*8 en(nquant),vect(nquant,nquant)

    gamma=1/25.d-15
    temperature=77.d0
    lambda=35.00*wave_to_J/pi

    H_diab(1,1)=107.d0*wave_to_J
    H_diab(2,2)=-107.d0*wave_to_J
    H_diab(1,2)=27.d0*wave_to_J  ;  H_diab(2,1)=H_diab(1,2)
    
  
    
    do i=1,nquant
       lambda_diab(i,i)=lambda*pi
    enddo

    call diag(H_diab,nquant,en,vect,nquant)
    c_tr=vect
    E_exc=en


end subroutine setup_parameters
!------------------------------------------------------------------------------------------

subroutine init_cond
   implicit none
   
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
           write(107,*) sigma(i,i),sigma(j,j),curr_time*1.d15
           write(109,*) RR(i,j),RR(j,i),curr_time*1.d15
           sigma_dot(i,i)=sigma_dot(i,i)+((RR(j,i)*sigma(j,j))-(RR(i,j)*sigma(i,i)))
           write(108,*) sigma_dot, curr_time*1.d15
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
              write(105,*) c_tr(i,n),c_tr(n,j),c_tr(j,n),c_tr(n,i),curr_time*1.d15
              write(106,*) R(i,j) 
              RR(i,j)=RR(i,j)+c_tr(i,n)*c_tr(n,j)*c_tr(j,n)*c_tr(n,i)*R(i,j)
           enddo
           write(102,*) RR(i,j),RR(j,i),curr_time*1.d15
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
          R(i,j)=((((2*lambda*gamma*omg(i,j))/((omg(i,j))**2+gamma**2))*((1/tanh((hbar*omg(i,j))&
          &/(2*kb*temperature)))+1))/hbar)
          write(101,*) R(i,j),i,j,curr_time*1.d15
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
            write(103,*) E_exc(i),E_exc(j),i,j,curr_time*1.d15,H_diab(1,1),H_diab(2,2)&
            &,H_diab(1,2),lambda
         enddo
      enddo

end subroutine calculate_freq   

!------------------------------------------------------------------------------------------
subroutine write_sigma

   implicit none
   
   sigma_diab=matmul(c_tr,matmul(sigma,transpose(c_tr)))
   write(100,*) curr_time*1.d15,real(sigma_diab(1,1)),real(sigma_diab(2,2))
     
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

End Module mod_rt   


 

   


   
