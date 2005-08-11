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
  use multigrid
  use poisson_corrections
  use math, only : dconjugate_gradients

  implicit none

  integer :: maxmulti,maxcycles
  FLOAT :: threshold
  integer :: presteps, poststeps
  integer :: restriction_method
  integer :: relaxation_method

  FLOAT, pointer :: old_solution(:)
  logical :: initialized=.false.

  integer, parameter :: GAUSS_SEIDEL = 1, &
                        CONJUGATE_GRADIENTS = 5
  
  private

  public  :: poisson_multigrid_solver, &
       poisson_multigrid_init, &
       poisson_multigrid_end


  type mg_float_pointer
     FLOAT, pointer :: p(:)
  end type mg_float_pointer

contains

  subroutine poisson_multigrid_init(m,ml,thr)
    type(mesh_type), intent(in) :: m
    integer, intent(in) :: ml
    FLOAT, intent(in)   :: thr

    call push_sub('poisson_multigrid.poisson_multigrid_init')
    maxmulti = ml
    threshold = thr

    call loct_parse_int(check_inp('PoissonSolverMGPresmoothingSteps'), 2, presteps)

    !%Variable PoissonSolverMGPresmoothingSteps
    !%Type integer
    !%Section 14 Varia
    !%Description
    !% Number of gauss-seidel smoothing steps before coarse level
    !% correction in multigrid Poisson solver. Default is 2.
    !%End

    call loct_parse_int(check_inp('PoissonSolverMGPostsmoothingSteps'), 3, poststeps)

    !%Variable PoissonSolverMGPostsmoothingSteps
    !%Type integer
    !%Section 14 Varia
    !%Description
    !% Number of gauss-seidel smoothing steps after coarse level
    !% correction in multigrid Poisson solver. Default is 3.
    !%End

    call loct_parse_int(check_inp('PoissonSolverMGMaxCycles'), 20, maxcycles)

    !%Variable PoissonSolverMGMaxCycles 
    !%Type integer
    !%Section 14 Varia
    !%Description
    !% Maximum number of multigrid cycles that are performed if
    !% convergence is not achieved. Default is 20.
    !%End


    call loct_parse_int(check_inp('PoissonSolverMGRestrictionMethod'), 2, restriction_method)

    !%Variable PoissonSolverMGRestrictionMethod
    !%Type integer
    !%Section 14 Varia
    !%Description
    !% Method used from fine to coarse grid. Default is fullweight restriction.
    !%Option injection 1
    !% Injection
    !%Option fullweight 2
    !% Fullweight restriction
    !%End

    call loct_parse_int(check_inp('PoissonSolverMGRestrictionMethod'), 2, restriction_method)

    !%Variable PoissonSolverMGRelaxationMethod
    !%Type integer
    !%Section 14 Varia
    !%Description
    !% Method used from fine to relax, i.e. to solve the linear system approximately, in
    !% the multigrid procedure that solve Poisson equation. For the moment, the option
    !% conjugate gradients is experimental...
    !%Option gauss-seidel 1
    !% Gauss-Seidel
    !%Option cg 5
    !% Conjugate-gradients
    !%End
   
    call loct_parse_int(check_inp('PoissonSolverMGRelaxationMethod'), GAUSS_SEIDEL, relaxation_method)
    if(relaxation_method.ne.GAUSS_SEIDEL .and. relaxation_method.ne.CONJUGATE_GRADIENTS) &
      call input_error('PoissonSolverMGRelaxationMethod')

    call build_aux(m)
    call build_phi(m)
    call pop_sub()
    
    allocate(old_solution(1:m%np))

    initialized=.true.

  end subroutine poisson_multigrid_init

  subroutine poisson_multigrid_end()
    deallocate(aux, phi)
    deallocate(old_solution)
    initialized=.false.
  end subroutine poisson_multigrid_end

  subroutine poisson_multigrid_solver(gr, pot, rho)
    type(grid_type), intent(inout) :: gr
    FLOAT, intent(inout) :: pot(:)
    FLOAT, intent(in)  :: rho(:)

    FLOAT :: res
    integer :: l, cl,t, np

    type(mg_float_pointer), pointer :: phi(:)
    type(mg_float_pointer), pointer :: tau(:)
    type(mg_float_pointer), pointer :: err(:)

    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)
    
    call push_sub('poisson_multigrid.multigrid_poisson_solver');
    
    if(.not. initialized) then
       message(1) = 'Multigrid Poisson not initialized.'
       call write_fatal(1)
    end if

    !!correction for treating boundaries 
    allocate(rho_corrected(gr%m%np), vh_correction(gr%m%np))
    call correct_rho(gr%m, maxmulti, rho, rho_corrected, vh_correction)
    rho_corrected = - M_FOUR*M_PI*rho_corrected
    pot = pot - vh_correction


    call gridhier_init(phi, gr%mgrid, add_points_for_boundaries=.true.)
    call gridhier_init(tau, gr%mgrid, add_points_for_boundaries=.false.)
    call gridhier_init(err, gr%mgrid, add_points_for_boundaries=.false.)

!    phi(0)%p(1:gr%m%np)=pot(1:gr%m%np)
    phi(0)%p(1:gr%m%np)=old_solution(1:gr%m%np)

    phi(0)%p(gr%m%np+1:gr%m%np_tot)=M_ZERO
    tau(0)%p(:)=rho_corrected(:)

    cl=gr%mgrid%n_levels
    
    do t=0,maxcycles

       if(t > 0) then !in the first cycle only the residue is calculated
          do l=0,(cl-1)
             
             if(l /= 0) then 
                phi(l)%p=M_ZERO
             end if
          
             !presmoothing 
             call multigrid_relax(gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
                  phi(l)%p, tau(l)%p , presteps)
             
             !error calcultion
             call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi(l)%p, err(l)%p);
             err(l)%p=tau(l)%p-err(l)%p
             
             !transfer error to coarser grid
             call multigrid_fine2coarse(gr%mgrid, l+1, err(l)%p, tau(l+1)%p,&
                  restriction_method)
             
          end do
          
          call multigrid_relax(gr%mgrid%level(cl)%m,gr%mgrid%level(cl)%f_der, &
               phi(cl)%p, tau(cl)%p , presteps)
          
          do l=cl,0,-1
             
             !postsmoothing 
             call multigrid_relax(gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
                  phi(l)%p, tau(l)%p , poststeps)
             
             if(l/= 0) then
                !transfer correction to finer level
                call multigrid_coarse2fine(gr%mgrid, l, phi(l)%p, err(l-1)%p)
                np=gr%mgrid%level(l-1)%m%np;
                phi(l-1)%p(1:np)=phi(l-1)%p(1:np)+err(l-1)%p(1:np);
             end if
             
          end do
       end if
       
       call df_laplacian(gr%sb, gr%mgrid%level(0)%f_der, phi(0)%p, err(0)%p);
       err(0)%p=err(0)%p-tau(0)%p
       res=sqrt(sum((err(0)%p(1:gr%m%np)*gr%m%vol_pp(1:gr%m%np))**2))

!       print*, "mgcycle", t, " res ", res

       if(res < threshold) then; 
          exit
       endif

    end do

    if(res >= threshold) then
       message(1) = 'Multigrid Poisson solver did not converge.'
       write(message(2), '(a,e14.6)') '  Res = ', res
       call write_warning(2)
    endif
    
    old_solution(1:gr%m%np)=phi(0)%p(1:gr%m%np)
    pot(1:gr%m%np)=phi(0)%p(1:gr%m%np)+vh_correction(1:gr%m%np);

    deallocate(rho_corrected, vh_correction)
    
    call gridhier_end(phi,gr%mgrid)
    call gridhier_end(tau,gr%mgrid)
    call gridhier_end(err,gr%mgrid)
    
    call pop_sub();
    
  end subroutine poisson_multigrid_solver

#define LAP f_der%der_discr%lapl


  subroutine multigrid_relax(m,f_der,pot,rho,steps)

    type(mesh_type),  target, intent(in)    :: m
    type(f_der_type), target, intent(in)    :: f_der
    FLOAT,            intent(inout) :: pot(:)  ! pot(m%np)
    FLOAT,            intent(in)    :: rho(:)  ! rho(m%np)
    integer,          intent(in)    :: steps

    integer :: t
    integer :: i,n, iter
    FLOAT :: point_lap, h2, factor
    FLOAT, allocatable :: w(:)

    select case(relaxation_method)

    case(GAUSS_SEIDEL)

      factor=(-M_ONE/(-M_SIX));

      !!search for the diagonal term
      !!DISABLED: currently using the diagonal term prevents convergence
      !    if(LAP%const_w) then 
      !       do i=1,LAP%n
      !          if( 1 == LAP%i(i,1)) then
      !             factor=-1.0/LAP%w_re(i,1);
      !             exit
      !          end if
      !       end do
      !    end if
    
      h2=m%h(1)*m%h(1)
      n=LAP%n
      allocate(w(1:n))
      if(LAP%const_w) then

        w(1:n)=LAP%w_re(1:n,1)
        do t=0,steps
           do i=1,m%np
              point_lap=sum(w(1:n)*pot(LAP%i(1:n,i)))
              pot(i)=pot(i)+factor*h2*(point_lap-rho(i))
           end do
        end do

      else

        do t=0,steps
           do i=1,m%np
              point_lap=sum(LAP%w_re(1:n,i)*pot(LAP%i(1:n,i)))
              pot(i)=pot(i)+factor*h2*(point_lap-rho(i))
           end do
        end do

      end if
      deallocate(w)

    case(CONJUGATE_GRADIENTS)

      iter = steps
      mesh_pointer => m ; der_pointer => f_der%der_discr
      call dconjugate_gradients(m%np, pot(1:m%np), rho(1:m%np), op, dotp, iter, threshold=CNST(1.0e-12))
      nullify(mesh_pointer) ; nullify(der_pointer)

    end select

  end subroutine multigrid_relax

#undef LAP

  subroutine gridhier_init(a, mgrid, add_points_for_boundaries)
    type(mg_float_pointer),pointer :: a(:)
    type(multigrid_type), pointer :: mgrid
    logical, intent(in) :: add_points_for_boundaries
    integer :: cl,l

    cl=mgrid%n_levels

    allocate(a(0:cl))
    do l=0,cl
       if(add_points_for_boundaries) then 
          allocate(a(l)%p(1:mgrid%level(l)%m%np_tot))
       else
          allocate(a(l)%p(1:mgrid%level(l)%m%np))
       end if
    end do
    
  end subroutine gridhier_init
 
  subroutine gridhier_end(a, mgrid)
    type(mg_float_pointer),pointer :: a(:)
    type(multigrid_type), pointer :: mgrid
    integer :: cl,l

    cl=mgrid%n_levels

    do l=0,cl
       deallocate(a(l)%p)
    end do
    deallocate(a)

  end subroutine gridhier_end

end module poisson_multigrid

