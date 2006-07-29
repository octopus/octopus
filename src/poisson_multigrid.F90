!! Copyright (C) 2005-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module poisson_multigrid_m
  use global_m
  use datasets_m
  use lib_oct_parser_m
  use messages_m
  use mesh_m
  use functions_m
  use grid_m
  use output_m
  use multigrid_m
  use poisson_corrections_m
  use math_m, only : dconjugate_gradients
  use varinfo_m

  implicit none

  integer ::                  &
    maxmulti,                 &
    maxcycles,                &
    presteps,                 &
    poststeps,                &
    restriction_method,       &
    relaxation_method

  FLOAT ::                    &
    threshold,                &
    relax_factor

  integer, parameter ::       &
    GAUSS_SEIDEL        = 1,  &
    GAUSS_JACOBI        = 2,  &
    GAUSS_JACOBI2       = 3,  &
    CONJUGATE_GRADIENTS = 5

  logical :: initialized = .false.


  private
  public ::                   &
    poisson_multigrid_solver, &
    poisson_multigrid_init,   &
    poisson_multigrid_end,    &
    mg_float_pointer,         &
    gridhier_init,            &
    gridhier_end


  type mg_float_pointer
    FLOAT, pointer :: p(:)
  end type mg_float_pointer

contains

  ! ---------------------------------------------------------
  subroutine poisson_multigrid_init(m, ml, thr)
    type(mesh_t), intent(in) :: m
    integer, intent(in) :: ml
    FLOAT, intent(in)   :: thr

    call push_sub('poisson_multigrid.poisson_multigrid_init')
    maxmulti  = ml
    threshold = thr

    !%Variable PoissonSolverMGPresmoothingSteps
    !%Type integer
    !%Default 3
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of gauss-seidel smoothing steps before coarse level
    !% correction in multigrid Poisson solver.
    !%End
    call loct_parse_int(check_inp('PoissonSolverMGPresmoothingSteps'), 2, presteps)

    !%Variable PoissonSolverMGPostsmoothingSteps
    !%Type integer
    !%Default 3
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of gauss-seidel smoothing steps after coarse level
    !% correction in multigrid Poisson solver
    !%End
    call loct_parse_int(check_inp('PoissonSolverMGPostsmoothingSteps'), 3, poststeps)

    !%Variable PoissonSolverMGMaxCycles
    !%Type integer
    !%Default 20
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Maximum number of multigrid cycles that are performed if
    !% convergence is not achieved
    !%End
    call loct_parse_int(check_inp('PoissonSolverMGMaxCycles'), 40, maxcycles)

    !%Variable PoissonSolverMGRestrictionMethod
    !%Type integer
    !%Default fullweight
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Method used from fine to coarse grid transfer.
    !%Option injection 1
    !% Injection
    !%Option fullweight 2
    !% Fullweight restriction
    !%End
    call loct_parse_int(check_inp('PoissonSolverMGRestrictionMethod'), 2, restriction_method)
    if(.not.varinfo_valid_option('PoissonSolverMGRestrictionMethod', restriction_method)) &
       call input_error('PoissonSolverMGRestrictionMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRestrictionMethod", restriction_method)

    !%Variable PoissonSolverMGRelaxationMethod
    !%Type integer
    !%Default Gauss-Seidel
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Method used to solve the linear system approximately in each grid for the
    !% multigrid procedure that solve Poisson equation. For the moment, the option
    !% conjugate gradients is experimental. 
    !%Option gauss_seidel 1
    !% Gauss-Seidel
    !%Option gauss_jacobi 2
    !% Gauss-Jacobi
    !%Option gauss_jacobi2 3
    !% Alternative implementation of Gauss-Jacobi.
    !%Option cg 5
    !% Conjugate-gradients
    !%End
    if ( m%use_curvlinear ) then
      call loct_parse_int(check_inp('PoissonSolverMGRelaxationMethod'), GAUSS_JACOBI, relaxation_method)
    else
      call loct_parse_int(check_inp('PoissonSolverMGRelaxationMethod'), GAUSS_SEIDEL, relaxation_method)
    end if

    if(.not.varinfo_valid_option('PoissonSolverMGRelaxationMethod', relaxation_method)) &
      call input_error('PoissonSolverMGRelaxationMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRelaxationMethod", relaxation_method)

    !%Variable PoissonSolverMGRelaxationFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Relaxation factor of the relaxation operator used for the
    !% multigrid method, default is 1.0. This is mainly for debugging,
    !% overrelaxation does nos help in a multigrid scheme.
    !%End
    if ( relaxation_method == GAUSS_JACOBI) then
      call loct_parse_float(check_inp('PoissonSolverMGRelaxationFactor'), CNST(0.5), relax_factor )
    else
      call loct_parse_float(check_inp('PoissonSolverMGRelaxationFactor'), M_ONE, relax_factor)
    end if

    call build_aux(m)
    call build_phi(m)

    initialized = .true.

    call pop_sub()
  end subroutine poisson_multigrid_init


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_end()
    call push_sub('poisson_multigrid.poisson_multigrid_end')

    deallocate(aux, phi)
    initialized = .false.

    call pop_sub()
  end subroutine poisson_multigrid_end


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_solver(gr, pot, rho)
    type(grid_t), intent(inout) :: gr
    FLOAT, intent(inout) :: pot(:)
    FLOAT, intent(in)    :: rho(:)

    integer :: l, cl, t, np, curr_l
    FLOAT :: res, fmg_threshold = CNST(0.1)
    FLOAT, allocatable :: rho_corrected(:), vh_correction(:)

    type(mg_float_pointer), pointer :: phi(:)
    type(mg_float_pointer), pointer :: phi_ini(:)
    type(mg_float_pointer), pointer :: tau(:)
    type(mg_float_pointer), pointer :: err(:)


    call push_sub('poisson_multigrid.poisson_multigrid_solver');

    ASSERT(initialized)

    ! correction for treating boundaries
    ALLOCATE(rho_corrected(gr%m%np), gr%m%np)
    ALLOCATE(vh_correction(gr%m%np), gr%m%np)

    call correct_rho(gr%m, maxmulti, rho, rho_corrected, vh_correction)
    rho_corrected = - M_FOUR*M_PI*rho_corrected
    pot = pot - vh_correction

    call gridhier_init(phi, gr%mgrid, add_points_for_boundaries=.true.)
    call gridhier_init(phi_ini, gr%mgrid, add_points_for_boundaries=.true.)
    call gridhier_init(tau, gr%mgrid, add_points_for_boundaries=.false.)
    call gridhier_init(err, gr%mgrid, add_points_for_boundaries=.false.)

    phi(0)%p(1:gr%m%np) = pot(1:gr%m%np)
    phi(0)%p(gr%m%np+1:gr%m%np_part) = M_ZERO
    tau(0)%p(:) = rho_corrected(:)

    cl = gr%mgrid%n_levels

    curr_l = 0

    do t= 0, maxcycles

      if(t > 0) then ! during the first cycle only the residue is calculated
        call vcycle_fas(curr_l)
      end if

      res = residue(curr_l, phi(curr_l)%p, tau(curr_l)%p, err(curr_l)%p)

      !if ( gr%m%use_curvlinear ) then
      !      print *, "base level", curr_l, "iter", t, " res ", res
      !end if

      if(res < threshold) then
        if(curr_l > 0 ) then
          call multigrid_coarse2fine(gr%mgrid, curr_l, phi(curr_l)%p, phi(curr_l-1)%p)
          curr_l = curr_l-1
        else
          exit
        end if
      end if

      if(0 == t .and. res > fmg_threshold) then
        curr_l = max(cl-1, 0)
        do l= 0, curr_l-1
          call multigrid_fine2coarse(gr%mgrid, l+1, tau(l)%p, tau(l+1)%p,&
            restriction_method)
        end do
      end  if

    end do

    if(res >= threshold) then
      message(1) = 'Multigrid Poisson solver did not converge.'
      write(message(2), '(a,e14.6)') '  Res = ', res
      call write_warning(2)
    end if

    pot(1:gr%m%np) = phi(0)%p(1:gr%m%np) + vh_correction(1:gr%m%np)

    deallocate(rho_corrected, vh_correction)

    call gridhier_end(phi,gr%mgrid)
    call gridhier_end(tau,gr%mgrid)
    call gridhier_end(err,gr%mgrid)
    call gridhier_end(phi_ini,gr%mgrid)

    call pop_sub()

  contains


    ! ---------------------------------------------------------
    FLOAT function residue(l, phi, rho, tmp) result(r)
      integer, intent(in) :: l
      FLOAT, pointer  :: phi(:)
      FLOAT, pointer  :: rho(:)
      FLOAT, pointer  :: tmp(:)
      type(mesh_t), pointer :: m

      m => gr%mgrid%level(l)%m
      call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi, tmp)
      tmp = tmp-rho
      r   = sqrt(sum((tmp(1:m%np)*m%vol_pp(1:m%np))**2))

    end function residue


    ! ---------------------------------------------------------
    subroutine vcycle_fas(fl)
      integer, intent(in) :: fl

      do l = fl, cl
        ! store the initial approximation
        np = gr%mgrid%level(l)%m%np_part
        phi_ini(l)%p(1:np) = phi(l)%p(1:np)

        ! presmoothing
        call multigrid_relax(gr,gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
          phi(l)%p, tau(l)%p , presteps)

        if (l /= cl ) then
          ! transfer of the current solution
          call multigrid_fine2coarse(gr%mgrid, l+1, phi(l)%p, phi(l+1)%p,&
            restriction_method)

          ! error calculation
          call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi(l)%p, err(l)%p)
          err(l)%p = err(l)%p - tau(l)%p

          ! transfer error to coarser grid
          call multigrid_fine2coarse(gr%mgrid, l+1, err(l)%p, tau(l+1)%p,&
            restriction_method)

          ! the other part of the error
          call df_laplacian(gr%sb, gr%mgrid%level(l+1)%f_der, phi(l+1)%p, err(l+1)%p)
          tau(l+1)%p = err(l+1)%p - tau(l+1)%p
        end if
      end do

      do l = cl, fl, -1
        ! postsmoothing
        call multigrid_relax(gr, gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
          phi(l)%p, tau(l)%p , poststeps)

        if(l /= fl) then
          ! calculate correction as the diference with the
          ! original grid in this level
          np = gr%mgrid%level(l)%m%np_part
          phi(l)%p(1:np) = phi(l)%p(1:np) - phi_ini(l)%p(1:np)

          ! transfer correction to finer level
          call multigrid_coarse2fine(gr%mgrid, l, phi(l)%p, err(l-1)%p)
          np = gr%mgrid%level(l-1)%m%np
          phi(l-1)%p(1:np) = phi(l-1)%p(1:np) + err(l-1)%p(1:np)
        end if
      end do

    end subroutine vcycle_fas


    ! ---------------------------------------------------------
    subroutine vcycle_cs(fl)
      integer, intent(in) :: fl

      do l = fl, cl-1
        if(l /= fl) then
          phi(l)%p=M_ZERO
        end if

        ! presmoothing
        call multigrid_relax(gr, gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
          phi(l)%p, tau(l)%p , presteps)

        ! error calcultion
        call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi(l)%p, err(l)%p);
        err(l)%p = tau(l)%p - err(l)%p

        ! transfer error to coarser grid
        call multigrid_fine2coarse(gr%mgrid, l+1, err(l)%p, tau(l+1)%p,&
          restriction_method)
      end do

      call multigrid_relax(gr, gr%mgrid%level(cl)%m,gr%mgrid%level(cl)%f_der, &
        phi(cl)%p, tau(cl)%p , presteps)

      do l = cl, fl, -1
        ! postsmoothing
        call multigrid_relax(gr, gr%mgrid%level(l)%m,gr%mgrid%level(l)%f_der, &
          phi(l)%p, tau(l)%p , poststeps)

        if( l /= fl ) then
          call df_laplacian(gr%sb, gr%mgrid%level(l)%f_der, phi(l)%p, err(l)%p)
          err(l)%p = err(l)%p - tau(l)%p
          np  = gr%mgrid%level(l)%m%np
          res = sqrt(sum((err(l)%p(1:np)*gr%mgrid%level(l)%m%vol_pp(1:np))**2))
        end if

        if(l/= fl) then
          ! transfer correction to finer level
          call multigrid_coarse2fine(gr%mgrid, l, phi(l)%p, err(l-1)%p)
          np = gr%mgrid%level(l-1)%m%np;
          phi(l-1)%p(1:np) = phi(l-1)%p(1:np) + err(l-1)%p(1:np);
        end if
      end do

    end subroutine vcycle_cs

  end subroutine poisson_multigrid_solver


  ! ---------------------------------------------------------
  subroutine multigrid_relax(gr, m,f_der,pot,rho,steps)
    type(grid_t),          intent(in)    :: gr
    type(mesh_t),  target, intent(in)    :: m
    type(f_der_t), target, intent(inout) :: f_der
    FLOAT,            intent(inout) :: pot(:)  ! pot(m%np)
    FLOAT,            intent(in)    :: rho(:)  ! rho(m%np)
    integer,          intent(in)    :: steps

    integer :: t
    integer :: i, n, iter, diag = 1
    FLOAT   :: point_lap, factor
    FLOAT, allocatable :: w(:), lpot(:), ldiag(:), tmp(:)

    select case(relaxation_method)

    case(GAUSS_SEIDEL)
      call search_diagonal()

      n=LAP%n

      ALLOCATE(w(1:n), n)

      if(LAP%const_w) then
        w(1:n) = LAP%w_re(1:n,1)
        do t = 0, steps
          do i = 1, m%np, 1
            point_lap = sum(w(1:n)*pot(LAP%i(1:n,i)))
            pot(i) = pot(i)+factor*(point_lap-rho(i))
          end do
        end do
      else
        do t = 0, steps
          do i = 1, m%np
            point_lap = sum(LAP%w_re(1:n,i)*pot(LAP%i(1:n,i)))
            pot(i) = pot(i) - CNST(0.7)/LAP%w_re(diag,i)*(point_lap-rho(i))
          end do
        end do
      end if

      deallocate(w)

    case(GAUSS_JACOBI)
      call search_diagonal()

      n=LAP%n

      ALLOCATE(w(1:n), n)
      ALLOCATE(tmp(1:m%np), m%np)

      if(LAP%const_w) then

        w(1:n) = LAP%w_re(1:n,1)
        do t = 0, steps
          do i = 1, m%np
            point_lap = sum(w(1:n)*pot(LAP%i(1:n,i)))
!            tmp(i) = -relax_factor/LAP%w_re(diag,1)*(point_lap-rho(i))
            tmp(i) = factor*(point_lap-rho(i))
          end do
          pot(1:m%np) = pot(1:m%np) + tmp(1:m%np)
        end do

      else

        do t = 0, steps
          do i = 1, m%np
            point_lap = sum(LAP%w_re(1:n,i)*pot(LAP%i(1:n,i)))
            tmp(i) = pot(i) + factor*(point_lap-rho(i))
          end do
          pot(1:m%np) = tmp(1:m%np)
        end do

      end if

      deallocate(w)
      deallocate(tmp)


    case(GAUSS_JACOBI2)
      call search_diagonal()

      ALLOCATE(lpot(1:m%np_part), m%np_part)
      ALLOCATE(ldiag(1:m%np), m%np)

      call df_laplacian_diag(gr%sb, f_der, ldiag)

      do t=1,steps
        call df_laplacian(gr%sb, f_der, pot, lpot)
        pot(1:m%np)=pot(1:m%np) - relax_factor/ldiag(1:m%np)*(lpot(1:m%np)-rho(1:m%np)) 
      end do

      deallocate(ldiag)
      deallocate(lpot)
      

    case(CONJUGATE_GRADIENTS)
      iter = steps
      mesh_pointer => m ; der_pointer => f_der%der_discr
      call dconjugate_gradients(m%np, pot(1:m%np), rho(1:m%np), op, dotp, iter, threshold=CNST(1.0e-12))
      nullify(mesh_pointer) ; nullify(der_pointer)

    end select

  contains

    ! ---------------------------------------------------------
    subroutine search_diagonal()
      do i = 1, LAP%n
        if( 1 == LAP%i(i,1)) then
          diag = i
          if(LAP%const_w) then
            factor = CNST(-1.0)/LAP%w_re(i,1)*relax_factor;
          end if
          exit
        end if
      end do
    end subroutine search_diagonal


  end subroutine multigrid_relax


  ! ---------------------------------------------------------
  subroutine gridhier_init(a, mgrid, add_points_for_boundaries)
    type(mg_float_pointer), pointer :: a(:)
    type(multigrid_t),      pointer :: mgrid
    logical,             intent(in) :: add_points_for_boundaries

    integer :: cl, l

    cl = mgrid%n_levels

    ALLOCATE(a(0:cl), cl+1)

    do l = 0, cl
      if(add_points_for_boundaries) then
        ALLOCATE(a(l)%p(1:mgrid%level(l)%m%np_part), mgrid%level(l)%m%np_part)
      else
        ALLOCATE(a(l)%p(1:mgrid%level(l)%m%np), mgrid%level(l)%m%np)
      end if
    end do

  end subroutine gridhier_init


  ! ---------------------------------------------------------
  subroutine gridhier_end(a, mgrid)
    type(mg_float_pointer), pointer :: a(:)
    type(multigrid_t),      pointer :: mgrid

    integer :: cl, l

    cl = mgrid%n_levels

    do l = 0, cl
      deallocate(a(l)%p)
    end do
    deallocate(a)

  end subroutine gridhier_end

end module poisson_multigrid_m
