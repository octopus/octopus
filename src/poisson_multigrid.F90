!! Copyright (C) 2005 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module poisson_multigrid
  use grid
  use output
  use poisson_corrections

  implicit none

  private float_pointer, &
       multigrid_relax

  type float_pointer
     FLOAT, pointer :: p(:)
  end type float_pointer
     

contains 

subroutine poisson_multigrid_solver(gr, pot, rho)
  type(grid_type), intent(inout) :: gr
  FLOAT, intent(out) :: pot(:)
  FLOAT, intent(in)  :: rho(:)
  
  type(mesh_type), pointer :: m_finest, m_coarse
  type(f_der_type), pointer  :: f_der_finest
  FLOAT, allocatable :: test(:), test2(:), ltest(:), ltest2(:)
  integer :: l, cl,t, np
  integer :: ierr

  type(float_pointer), allocatable :: phi(:)
  type(float_pointer), allocatable :: tau(:)
  type(float_pointer), allocatable :: err(:)
  
  integer, parameter :: presteps=2, poststeps=2
  integer, parameter :: maxl=4


  FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

  call push_sub('multigrid_poisson_solver');

  !!correction for treating boundaries 
  allocate(rho_corrected(gr%m%np), vh_correction(gr%m%np))
  call correct_rho(gr%m, maxl, rho, rho_corrected, vh_correction)
  rho_corrected = - M_FOUR*M_PI*rho_corrected
  pot = pot - vh_correction

  cl=gr%mgrid%n_levels

  allocate(phi(0:cl))
  allocate(tau(0:cl))
  allocate(err(0:cl))

  do l=0,cl
     allocate(phi(l)%p(1:gr%mgrid%level(l)%m%np_tot))
     allocate(tau(l)%p(1:gr%mgrid%level(l)%m%np))
     allocate(err(l)%p(1:gr%mgrid%level(l)%m%np))
  end do

  phi(0)%p(1:gr%m%np)=pot(1:gr%m%np)
  phi(0)%p(gr%m%np+1:gr%m%np_tot)=M_ZERO
  tau(0)%p(:)=rho_corrected(:)

  do t=1,10

#if 0
     l=0
     call multigrid_relax(gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, phi(l)%p, tau(l)%p , 1000);
     
#else 

!     cl=1
     do l=0,(cl-1)
!        print*, "level", l
        if(l /= 0) then 
           phi(l)%p=0.0
        end if
        
        !presmoothing 
        call multigrid_relax(gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, phi(l)%p, tau(l)%p , presteps);
        
        !error calcultion
        call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi(l)%p, err(l)%p);
        err(l)%p=tau(l)%p-err(l)%p
        
        !transfer error to coarser grid
        call multigrid_fine2coarse(gr%mgrid, l+1, err(l)%p, tau(l+1)%p)
        
     end do
     
     do l=cl,0,-1
!        print*, "level", l
        !postsmoothing 
        call multigrid_relax(gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, phi(l)%p, tau(l)%p , poststeps);
        
        if(l/= 0) then
           !transfer correction to finer level
           call multigrid_coarse2fine(gr%mgrid, l, phi(l)%p, err(l-1)%p)
           np=gr%mgrid%level(l-1)%m%np;
           phi(l-1)%p(1:np)=phi(l-1)%p(1:np)+err(l-1)%p(1:np);
      end if
      
   end do

#endif    
!   call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
!        "lixo", "corr", gr%mgrid%level(0)%m , gr%sb, err(0)%p, M_ONE, ierr)
   !error calculation 
   call df_laplacian(gr%sb, gr%mgrid%level(0)%f_der, phi(0)%p, err(0)%p);
   err(0)%p=tau(0)%p-err(0)%p
!   call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
!        "lixo", "err", gr%mgrid%level(0)%m , gr%sb, err(0)%p, M_ONE, ierr)
   print*, "mgcycle", t, " res ", dot_product(err(0)%p,err(0)%p)
end do

#if 0
  ! test transfer between coarse and fine grids
!  allocate(test(m_coarse%np), ltest(m_coarse%np), test2(m_finest%np), ltest2(m_finest%np))
!  call multigrid_fine2coarse(gr%mgrid, 1, rho, test)
!  call multigrid_coarse2fine(gr%mgrid, 1, test, test2)

!  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
!     "lixo", "fine", m_finest, m_finest%sb, rho, M_ONE, ierr)
!  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
!     "lixo", "coarse", m_coarse, m_coarse%sb, test, M_ONE, ierr)
!  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
!     "lixo", "fine2", m_finest, m_finest%sb, test2, M_ONE, ierr)

  ! test laplacians
  call df_laplacian(gr%sb, gr%mgrid%level(0)%f_der, rho, ltest2)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "l-fine", m_finest, m_finest%sb, ltest2, M_ONE, ierr)

  call df_laplacian(gr%sb, gr%mgrid%level(1)%f_der, test, ltest)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "l-coarse", m_coarse, m_coarse%sb, ltest, M_ONE, ierr)

#endif

  pot(1:gr%m%np)=phi(0)%p(1:gr%m%np)+vh_correction(1:gr%m%np);

  deallocate(rho_corrected, vh_correction)
  do l=0,cl
     deallocate(phi(l)%p)
     deallocate(tau(l)%p)
     deallocate(err(l)%p)
  end do
  deallocate(phi)
  deallocate(tau)
  deallocate(err)

  call pop_sub();

end subroutine poisson_multigrid_solver

#define LAP f_der%der_discr%lapl

subroutine multigrid_relax(m,f_der,pot,rho,steps) 
  
  type(mesh_type),  intent(in)    :: m
  type(f_der_type), intent(in)    :: f_der
  FLOAT,            intent(inout) :: pot(:)  ! pot(m%np)
  FLOAT,            intent(in)    :: rho(:)  ! rho(m%np)
  integer,          intent(in)    :: steps
  
  integer :: t
  integer :: i,j
  FLOAT :: point_lap, h2
  
  h2=m%h(1)*m%h(1)
  
  do t=0,steps
     do i=1,m%np
        
        point_lap=0.0;
        do j=1,LAP%n
           point_lap=point_lap+LAP%w_re(j,1)*pot(LAP%i(j,i))
        end do
        pot(i)=pot(i)-1.0/(CNST(-6.0))*h2*(point_lap-rho(i)) 
        
!!! here should be the diagonal part of the laplacian operator instead of -6.0
!!! may be we could just use: 
        !pot(i)=pot(i)-1.0/lap%w_re(1,1)*h2*(l-rho(i)) 
!!! but i'm not sure that the diagonal is always in w_re(1,1) 
!!! in any case, -6.0 is a safe bet 
!!! (is diagonal part for the 3 point 1 -2 1 laplacian) 
!!! and i'm not sure which one gives better convergence
        
     end do
  end do
  
end subroutine multigrid_relax
#undef LAP
end module poisson_multigrid
