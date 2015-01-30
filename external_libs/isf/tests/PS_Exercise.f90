!> @file
!! Exercise using the Poisson solver
!! @author
!!    Copyright (c) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
program exercise
   use Poisson_Solver
   implicit none
   character(len=1) :: solvertype,afunc
   character(len=64) :: chain
   
   integer :: i1,i2,i3,n1,n2,n3,i1_max,i2_max,i3_max,isf_order
   real :: t0,t1,t2,t3
   real(kind=8) :: pi,max_diff
   real(kind=8) :: sigma,length,hgrid,mu,energy,offset,acell,epot,intrhoS,intrhoF,intpotS,intpotF
   real(kind=8), dimension(:), allocatable :: fake_arr
   real(kind=8), dimension(:,:,:), allocatable :: psi,rhopot,rhoF,rhoS,potF,potS
   type(coulomb_operator) :: kernel
   
   !Use arguments
   call getarg(1,chain)
   if(trim(chain)=='') then
      write(*,'(1x,a)')&
           'Usage: ./PS_Exercise solvertype hgrid isf_order'
      stop
   end if
   read(unit=chain,fmt=*) solvertype
   call getarg(2,chain)
   read(unit=chain,fmt=*) hgrid
   call getarg(3,chain)
   read(unit=chain,fmt=*) isf_order

   if (solvertype == 'P') then
      write(*,'(1x,a)')'You have chosen the periodic Poisson solver: which function do you want to use?'
      write(*,'(1x,a)')'Type "S" for surfaces type and "F" for free type'
      read(5,*)afunc
   elseif (solvertype == 'W') then
      write(*,'(1x,a)')'You have chosen the Poisson solver for wires boundary conditions: which function do you want to use?'
      write(*,'(1x,a)')'Type "S" for surfaces type and "F" for free type'
      read(5,*)afunc
   else
      afunc=solvertype
   end if

!!!  call getarg(4,chain)
!!!  read(unit=chain,fmt=*) ixc
!!!  call getarg(5,chain)
!!!  read(unit=chain,fmt=*) geocode
!!!  call getarg(6,chain)
!!!  read(unit=chain,fmt=*) datacode


   open(unit=12,file='results.out',status='unknown')
   write(12,*)'# isf  hgrid    max_diff              energy                   epot            Time:Kernel  Solver'

   pi=4.d0*datan(1.d0)
   acell=14.d0
   sigma=1.d0

!!!   !!!!!INPUT VARIABLES THAT CAN BE MODIFIED
!!!   !type of the solver:
!!!   !       'F' for Free (Isolated) BC
!!!   !       'S' for Surfaces BC
!!!   !       'P' for Periodic BC
!!!   solvertype='P'
!!!   !degree of interpolating scaling functions
!!!   isf_order=16
!!!   !grid spacing
!!!   hgrid=0.1d0
!!!   !!!!!END OF THE INPUT VARIABLES

   !in the case of a multiple run, consider the possibility to 
   !uncomment the following lines for making a loop cycle on the grid spacing
   !do i=1,10
   !hgrid=0.1d0*real(i,kind=8)

   n1=8*(int(acell/hgrid)/8+1)
   length=n1*hgrid
   n2=n1
   n3=n1
   mu=0.5d0/sigma**2

   write(*,'(1x,a)')''
   write(*,'(1x,a)')'============================================================'
   write(*,'(1x,a,a,a,f5.2,a,i4,a,3(i4))')'solvertype= ',solvertype,' hgrid=',hgrid,&
        ' isf_order=',isf_order,' n1,n2,n3=',n1,n2,n3
   write(*,'(1x,a)')'============================================================'
   write(*,'(1x,a)')''

   !allocate analytic functions for arrays
   allocate(rhoF(n1,n2,n3),rhoS(n1,n2,n3),potF(n1,n2,n3),potS(n1,n2,n3),psi(n1,n2,n3))

   call assign_functions(n1,n2,n3,hgrid,mu,length,rhoF,potF,rhoS,potS,psi,&
        intrhoS,intrhoF,intpotS,intpotF)

   write(*,'(1x,a)')'Integrals of the input functions:'
   if (afunc=='S') then
      write(*,'(2(1x,a,1pe25.16))')'intrhoS=',intrhoS,'intpotS=',intpotS
   else
      write(*,'(2(1x,a,1pe25.16))')'intrhoF=',intrhoF,'intpotF=',intpotF
   end if

   if (solvertype == 'P' .or. solvertype == 'W') then
      write(*,'(1x,a)')'Which value of the zero-th fourier component do you want to put for the solution?'
      read(5,*)offset
   end if

   !after the assignment of the arrays, solve the Poisson equation with the different
   !cases
   allocate(rhopot(n1,n2,n3),fake_arr(1))

   if (afunc=='F') then
   rhopot=rhoF
   else if (afunc=='S') then
   rhopot=rhoS
   end if

   !offset=0.d0!3.053506154731705d0*n1*n2*n3*hgrid**3
   
   call cpu_time(t0)
   kernel=pkernel_init(.true.,0,1,0,&
        solvertype,(/n1,n2,n3/),(/hgrid,hgrid,hgrid/),isf_order)
   call pkernel_set(kernel,.true.)

   !call createKernel(0,1,solvertype,(/n1,n2,n3/),(/hgrid,hgrid,hgrid/),isf_order,kernel,.true.)
   call cpu_time(t1)

   call cpu_time(t2)
!   call PSolver(solvertype,'G',0,1,n1,n2,n3,0,hgrid,hgrid,hgrid,&
!        rhopot,kernel%kernel,fake_arr,energy,zero,zero,offset,.false.,1)
   call H_potential('G',kernel,rhopot,fake_arr,energy,offset,.false.,quiet='yes') !optional argument

   call cpu_time(t3)

   call pkernel_free(kernel)
   !deallocate(kernel)


   !now calculate the maximum difference and compare with the analytic result


   !ANALYTIC POTENTIAL OF REFERENCE
   !                      v
   !                      v
   if (afunc=='F') then
   call compare(n1,n2,n3,potF,rhopot,i1_max,i2_max,i3_max,max_diff)
   else if (afunc=='S') then
   call compare(n1,n2,n3,potS,rhopot,i1_max,i2_max,i3_max,max_diff)
   end if
   !                      ^
   !                      ^
   !BEWARE TO CHANGE THIS TERM IF YOU CHANGE REFERENCE

   call potential_energy(n1,n2,n3,hgrid,psi,rhopot,epot)

   write(*,'(1x,a)')''
   write(*,'(1x,a)')'============================================================'
   write(*,'(1x,3(a,1pe20.12),1x,a,0pf6.3)')'Accuracy=',max_diff,' Hartree energy=',energy,' Potential energy= ',epot
   write(*,'(2(1x,a,0pf6.3))')'time: Kernel=',t1-t0,' Solver=',t3-t2
   write(*,'(1x,a)')'============================================================'
   write(*,'(1x,a)')' The functions are printed in the file "functions.dat" and the results in "results.out"'

   open(unit=11,file='functions.dat',status='unknown')
   if (afunc=='F') then
      !plot the functions
      i3=n3/2
      i1=n1/2
      do i2=1,n2
         write(11,'(i0,3(2x,1pe27.16e3))')i2,rhoF(i1,i2,i3),potF(i1,i2,i3),rhopot(i1,i2,i3)
      end do
   else if (afunc=='S') then
      !plot the functions
      i3=n3/2
      i1=n1/2
      do i2=1,n2
         write(11,'(i0,3(2x,1pe27.16e3))')i2,rhoS(i1,i2,i3),potS(i1,i2,i3),rhopot(i1,i2,i3)
      end do
   end if

   close(unit=11)

   write(12,'(1x,i3,3x,f5.2,2x,1pe12.5,2(2x,1pe25.16),2(3x,0pf6.3))')&
        isf_order,hgrid,max_diff,energy,epot,t1-t0,t3-t2

   deallocate(rhoF,rhoS,potS,potF,rhopot,fake_arr,psi)

   !to be uncommented in the case of a loop cycle
   !end do

   close(unit=12)
      
 end program exercise

 subroutine assign_functions(n1,n2,n3,hgrid,mu,length,rhoF,potF,rhoS,potS,psi,&
      intrhoS,intrhoF,intpotS,intpotF)
   implicit none
   integer, intent(in) :: n1,n2,n3
   real(kind=8), intent(in) :: hgrid,mu,length
   real(kind=8), intent(out) :: intrhoS,intrhoF,intpotS,intpotF
   real(kind=8), dimension(n1,n2,n3), intent(out) :: rhoF,potF,rhoS,potS,psi
   !local variables
   integer :: i1,i2,i3
   real(kind=8) :: a_gauss,a2,pi,factor,offset,x,y,z,r,r2,fx,fy,fz,fx2,fy2,fz2,zero,norm,derf
   !assign the function that must be analytically compared
   pi=4.d0*datan(1.d0)
   a_gauss=1.d0/sqrt(mu)
   a2 = 1.d0/mu
   !Normalization
   factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
   !offset
   !integral of the true periodic solution
   offset=0.d0
   intrhoS=0.d0
   intrhoF=0.d0
   intpotF=0.d0
   intpotS=0.d0
   norm=0.d0
   do i3=1,n3
      z=real(i3,kind=8)*hgrid-0.5d0*length
      call functions(z,length,zero,fz,fz2,5)
      do i2=1,n2
         y=real(i2,kind=8)*hgrid-0.5d0*length
         !call functions(y,1.d0/(2.d0*sigma**2),zero,fy,fy2,2)
         call functions(y,length,zero,fy,fy2,6)
         do i1=1,n1
            x=real(i1,kind=8)*hgrid-0.5d0*length
            call functions(x,length,zero,fx,fx2,5)
            r2=x**2+y**2+z**2
            !rhoF(i1,i2,i3)=exp(-r2/a2)
            !! CHECK WITH KRONECKER DELTA !!
            if(i1 == n1/2 .and. i2 == n2/2 .and. i3 == n3/2) then
               rhoF(i1,i2,i3) = 1.d0
            else
               rhoF(i1,i2,i3) = 0.d0
            end if
            !rhoS(i1,i2,i3)=fx2*fy*fz+fx*fy2*fz+fx*fy*fz2
            rhoS(i1,i2,i3) = rhoF(i1,i2,i3)
            if (r2==0.d0) then
               potF(i1,i2,i3)= 2.d0/(sqrt(pi)*a_gauss)/factor
            else
               r=sqrt(r2)
               potF(i1,i2,i3)= derf(r/a_gauss)/r/factor
            end if
            potS(i1,i2,i3)=-4.d0*pi*fx*fy*fz
            intrhoS=intrhoS+rhoS(i1,i2,i3)
            intrhoF=intrhoF+rhoF(i1,i2,i3)
            intpotF=intpotF+potF(i1,i2,i3)
            intpotS=intpotS+potS(i1,i2,i3)
            psi(i1,i2,i3)=rhoF(i1,i2,i3)*sqrt(2.d0*sqrt(2.d0)*factor)
            norm=norm+potF(i1,i2,i3)*rhoF(i1,i2,i3)
         end do
      end do
   end do
   intrhoS=intrhoS*hgrid**3
   intrhoF=intrhoF*hgrid**3
   intpotF=intpotF*hgrid**3
   intpotS=intpotS*hgrid**3
   norm=norm*hgrid**3

   !print *,'norm',norm
!!!   !plot of the functions used
!!!   do i1=1,n1
!!!      x = hgrid*real(i1,kind=8)-0.5d0*length!valid if hy=hz
!!!      y = hgrid*real(i1,kind=8)-0.5d0*length 
!!!      call functions(x,length,zero,fx,fx2,5)
!!!      call functions(y,length,zero,fy,fy2,6)
!!!      write(20,*)i1,fx,fx2,fy,fy2,rhoS(n1/2,i1,n3/2),potS(n1/2,i1,n3/2),&
!!!           rhoF(i1,n2/2,n3/2),potF(i1,n2/2,n3/2)
!!!   end do

 end subroutine assign_functions

 subroutine potential_energy(n1,n2,n3,hgrid,psi,pot,energy)
   implicit none
   integer, intent(in) :: n1,n2,n3
   real(kind=8), intent(in) :: hgrid
   real(kind=8), dimension(n1*n2*n3), intent(in) :: psi,pot
   real(kind=8), intent(out) :: energy
   !local variables
   integer :: i
   real(kind=8) :: norm
   energy=0.d0
   norm=0.d0
   do i=1,n1*n2*n3
      energy=energy+psi(i)**2*pot(i)
      norm=norm+psi(i)**2
   end do
   norm=norm*hgrid**3
   !print *,norm,'norm'
   energy=energy*hgrid**3
 end subroutine potential_energy


subroutine functions(x,a,b,f,f2,whichone)
  implicit none
  integer, intent(in) :: whichone
  real(kind=8), intent(in) :: x,a,b
  real(kind=8), intent(out) :: f,f2
  !local variables
  real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
  real(kind=8) :: length,frequency,nu,sigma,agauss,derf

  pi = 4.d0*datan(1.d0)
  select case(whichone)
  case(1)
     !constant
     f=1.d0
     f2=0.d0
  case(2)
     !gaussian of sigma s.t. a=1/(2*sigma^2)
     r2=a*x**2
     f=dexp(-r2)
     f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
  case(3)
     !gaussian "shrinked" with a=length of the system
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     f2=factor*dexp(-y**2)
     f=dexp(-y**2)
  case(4)
     !cosine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dcos(r)
     f2=-(frequency*pi/length)**2*dcos(r)
  case(5)
     !exp of a cosine, a=length
     nu=2.d0
     r=pi*nu/a*x
     y=dcos(r)
     yp=dsin(r)
     f=dexp(y)
     factor=(pi*nu/a)**2*(-y+yp**2)
     f2=factor*f
  case(6)
     !gaussian times "shrinked" gaussian, sigma=length/10
     length=a
     r=pi*x/length
     y=dtan(r)
     yp=pi/length*1.d0/(dcos(r))**2
     ys=2.d0*pi/length*y*yp
     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
     g=dexp(-y**2)
     g1=-2.d0*y*yp*g
     g2=factor*dexp(-y**2)
     
     sigma=length/10
     agauss=0.5d0/sigma**2
     r2=agauss*x**2
     h=dexp(-r2)
     h1=-2.d0*agauss*x*h
     h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
     f=g*h
     f2=g2*h+g*h2+2.d0*g1*h1
  case(7)
     !sine with a=length, b=frequency
     length=a
     frequency=b
     r=frequency*pi*x/length
     f=dsin(r)
     f2=-(frequency*pi/length)**2*dsin(r)
  case(8)
     !error function with a=sigma
     factor=sqrt(2.d0/pi)/a
     r=x
     y=x/(sqrt(2.d0)*a)
     if (abs(x)<=1.d-15) then
        f=factor
        f2=-sqrt(2.d0/pi)/(3.d0*a**3)
     else
        f=derf(y)/r
        y=x*x
        y=y/(2.d0*a**2)
        g=dexp(-y)
        h=1.d0/a**2+2.d0/x**2
        f2=-factor*g*h+2.d0*f/x**2
     end if
  end select

end subroutine functions

subroutine compare(n01,n02,n03,potential,density,i1_max,i2_max,i3_max,max_diff)
  implicit none
  integer, intent(in) :: n01,n02,n03
  real(kind=8), dimension(n01,n02,n03), intent(in) :: potential,density
  integer, intent(out) :: i1_max,i2_max,i3_max
  real(kind=8), intent(out) :: max_diff
  !local variables
  integer :: i1,i2,i3
  real(kind=8) :: factor,max_val
  max_diff = 0.d0
  max_val = 0.d0
  i1_max = 1
  i2_max = 1
  i3_max = 1
  do i3=1,n03
     do i2=1,n02 
        do i1=1,n01
           max_val=max(max_val,abs(potential(i1,i2,i3)))
           factor=abs(potential(i1,i2,i3)-density(i1,i2,i3))
           if (max_diff < factor) then
              max_diff = factor
              i1_max = i1
              i2_max = i2
              i3_max = i3
           end if
        end do
     end do
  end do
  max_diff=max_diff/max_val
end subroutine compare
