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
!! $Id$

#include "global.h"

module poisson_multigrid_m
  use datasets_m
  use derivatives_m
  use global_m
  use grid_m
  use gridhier_m
  use lalg_basic_m
  use loct_parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use multigrid_m
  use par_vec_m
  use poisson_corrections_m
  use profiling_m
  use varinfo_m

  implicit none

  integer, parameter ::       &
    GAUSS_SEIDEL        = 1,  &
    GAUSS_JACOBI        = 2,  &
    GAUSS_JACOBI2       = 3

  private
  public ::                   &
    poisson_multigrid_solver, &
    poisson_multigrid_init,   &
    poisson_multigrid_end,    &
    mg_solver_t

  type mg_solver_t

     FLOAT ::                    &
          threshold,                &
          relax_factor
     
     integer ::                     &
          maxcycles,                &
          presteps,                 &
          poststeps,                &
          restriction_method,       &
          relaxation_method
     
     type(poisson_corr_t) :: corrector
     
  end type mg_solver_t

contains

  ! ---------------------------------------------------------
  subroutine poisson_multigrid_init(this, m, ml, thr)
    type(mg_solver_t), intent(inout) :: this
    type(mesh_t),      intent(inout) :: m
    integer,           intent(in)    :: ml
    FLOAT,             intent(in)    :: thr

    call push_sub('poisson_multigrid.poisson_multigrid_init')
    call poisson_corrections_init(this%corrector, ml, m)

    this%threshold = thr

    !%Variable PoissonSolverMGPresmoothingSteps
    !%Type integer
    !%Default 3
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of gauss-seidel smoothing steps before coarse level
    !% correction in the multigrid Poisson solver. By default 1.
    !%End
    call loct_parse_int(datasets_check('PoissonSolverMGPresmoothingSteps'), 1, this%presteps)

    !%Variable PoissonSolverMGPostsmoothingSteps
    !%Type integer
    !%Default 3
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of gauss-seidel smoothing steps after coarse level
    !% correction in the multigrid Poisson solver. By default 4.
    !%End
    call loct_parse_int(datasets_check('PoissonSolverMGPostsmoothingSteps'), 4, this%poststeps)

    !%Variable PoissonSolverMGMaxCycles
    !%Type integer
    !%Default 60
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Maximum number of multigrid cycles that are performed if
    !% convergence is not achieved. By default 60.
    !%End
    call loct_parse_int(datasets_check('PoissonSolverMGMaxCycles'), 60, this%maxcycles)

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
    call loct_parse_int(datasets_check('PoissonSolverMGRestrictionMethod'), 2, this%restriction_method)
    if(.not.varinfo_valid_option('PoissonSolverMGRestrictionMethod', this%restriction_method)) &
       call input_error('PoissonSolverMGRestrictionMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRestrictionMethod", this%restriction_method)

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
    !%End
    if ( m%use_curvilinear ) then
      call loct_parse_int(datasets_check('PoissonSolverMGRelaxationMethod'), GAUSS_JACOBI, this%relaxation_method)
    else
      call loct_parse_int(datasets_check('PoissonSolverMGRelaxationMethod'), GAUSS_SEIDEL, this%relaxation_method)
    end if

    if(.not.varinfo_valid_option('PoissonSolverMGRelaxationMethod', this%relaxation_method)) &
      call input_error('PoissonSolverMGRelaxationMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRelaxationMethod", this%relaxation_method)

    !%Variable PoissonSolverMGRelaxationFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Relaxation factor of the relaxation operator used for the
    !% multigrid method, default is 1.0. This is mainly for debugging,
    !% overrelaxation does nos help in a multigrid scheme.
    !%End
    if ( this%relaxation_method == GAUSS_JACOBI) then
      call loct_parse_float(datasets_check('PoissonSolverMGRelaxationFactor'), CNST(0.6666), this%relax_factor )
    else
      call loct_parse_float(datasets_check('PoissonSolverMGRelaxationFactor'), M_ONE, this%relax_factor)
    end if

    call pop_sub()
  end subroutine poisson_multigrid_init


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_end(this)
    type(mg_solver_t), intent(inout) :: this

    call push_sub('poisson_multigrid.poisson_multigrid_end')

    call poisson_corrections_end(this%corrector)

    call pop_sub()
  end subroutine poisson_multigrid_end


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_solver(this, gr, pot, rho)
    type(mg_solver_t), intent(in)    :: this
    type(grid_t),      intent(inout) :: gr
    FLOAT,             intent(inout) :: pot(:)
    FLOAT,             intent(in)    :: rho(:)

    integer :: lev, cl, t, np, curr_l, ip
    FLOAT :: res, fmg_threshold = CNST(0.1)
    FLOAT, allocatable :: vh_correction(:)

    type(dgridhier_t) :: phi
    type(dgridhier_t) :: phi_ini
    type(dgridhier_t) :: tau
    type(dgridhier_t) :: err

    call push_sub('poisson_multigrid.poisson_multigrid_solver');

    ! correction for treating boundaries
    SAFE_ALLOCATE(vh_correction(1:gr%mesh%np))

    call gridhier_init(phi, gr%mgrid, add_points_for_boundaries=.true.)
    call gridhier_init(phi_ini, gr%mgrid, add_points_for_boundaries=.false.)
    call gridhier_init(tau, gr%mgrid, add_points_for_boundaries=.true.)
    call gridhier_init(err, gr%mgrid, add_points_for_boundaries=.true.)

    call correct_rho(this%corrector, gr%mesh, rho, tau%level(0)%p, vh_correction)
    call lalg_scal(gr%mesh%np, -M_FOUR*M_PI, tau%level(0)%p)

    forall (ip = 1:gr%mesh%np) phi%level(0)%p(ip) = pot(ip) - vh_correction(ip)

    cl = gr%mgrid%n_levels

    curr_l = 0

    do t = 0, this%maxcycles

      if(t > 0) then ! during the first cycle only the residue is calculated
        call vcycle_cs(curr_l)
      end if

      res = residue(curr_l, phi%level(curr_l)%p, tau%level(curr_l)%p, err%level(curr_l)%p)
      
      if(in_debug_mode) then
        write(message(1), *) "Multigrid: base level ", curr_l, " iter ", t, " res ", res
        call write_info(1)
      end if

      if(res < this%threshold) then
        if(curr_l > 0 ) then
          call dmultigrid_coarse2fine(gr%mgrid%level(curr_l)%tt, gr%mgrid%level(curr_l)%der, &
            gr%mgrid%level(curr_l)%mesh, gr%mgrid%level(curr_l - 1)%mesh, &
            phi%level(curr_l)%p, phi%level(curr_l - 1)%p, order = 2)
          curr_l = curr_l - 1
        else
          exit
        end if
      end if

      if(0 == t .and. res > fmg_threshold) then
        curr_l = max(cl - 1, 0)
        do lev = 0, curr_l - 1
          call dmultigrid_fine2coarse(gr%mgrid, lev + 1, tau%level(lev)%p, tau%level(lev + 1)%p, this%restriction_method)
        end do
      end  if

    end do

    if(res >= this%threshold) then
      message(1) = 'Multigrid Poisson solver did not converge.'
      write(message(2), '(a,e14.6)') '  Res = ', res
      call write_warning(2)
    end if

    forall (ip = 1:gr%mesh%np) pot(ip) = phi%level(0)%p(ip) + vh_correction(ip)

    SAFE_DEALLOCATE_A(vh_correction)

    call gridhier_end(phi, gr%mgrid)
    call gridhier_end(tau, gr%mgrid)
    call gridhier_end(err, gr%mgrid)
    call gridhier_end(phi_ini, gr%mgrid)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    FLOAT function residue(level, phi, rho, tmp) result(res)
      integer, intent(in)    :: level
      FLOAT,   intent(inout) :: phi(:)
      FLOAT,   intent(in)    :: rho(:)
      FLOAT,   intent(inout) :: tmp(:)

      type(mesh_t), pointer :: m
      integer :: ip

      call push_sub('poisson_multigrid.residue')

      m => gr%mgrid%level(level)%mesh

      call dderivatives_lapl(gr%mgrid%level(level)%der, phi, tmp)

      forall (ip = 1:m%np) tmp(ip) = tmp(ip) - rho(ip)

      res = dmf_nrm2(m, tmp)

      call pop_sub()
    end function residue

    ! ---------------------------------------------------------
    subroutine vcycle_cs(fl)
      integer, intent(in) :: fl

      integer :: l, ip

      call push_sub('poisson_multigrid.vcycle_fas')

      do l = fl, cl
        np = gr%mgrid%level(l)%mesh%np
        if(l /= fl) phi%level(l)%p(1:np) = M_ZERO

        ! presmoothing
        call multigrid_relax(this, gr%mgrid%level(l)%mesh, gr%mgrid%level(l)%der, &
          phi%level(l)%p, tau%level(l)%p , this%presteps)

        if (l /= cl ) then
          ! error calculation
          call dderivatives_lapl(gr%mgrid%level(l)%der, phi%level(l)%p, err%level(l)%p)
          forall(ip = 1:np) err%level(l)%p(ip) = tau%level(l)%p(ip) - err%level(l)%p(ip)

          ! transfer error as the source in the coarser grid
          call dmultigrid_fine2coarse(gr%mgrid, l+1, err%level(l)%p, tau%level(l+1)%p, this%restriction_method)
        end if
      end do

      do l = cl, fl, -1
       
        ! postsmoothing
        call multigrid_relax(this, gr%mgrid%level(l)%mesh, gr%mgrid%level(l)%der, &
          phi%level(l)%p, tau%level(l)%p, this%poststeps)

        if(l /= fl) then
          ! transfer correction to finer level
          call dmultigrid_coarse2fine(gr%mgrid%level(l)%tt, gr%mgrid%level(l)%der, &
            gr%mgrid%level(l)%mesh, gr%mgrid%level(l - 1)%mesh, phi%level(l)%p, err%level(l-1)%p)

          np = gr%mgrid%level(l-1)%mesh%np
          forall(ip = 1:np) phi%level(l - 1)%p(ip) = phi%level(l - 1)%p(ip) + err%level(l - 1)%p(ip)
        end if
      end do

      call pop_sub()
    end subroutine vcycle_cs

    ! ---------------------------------------------------------
    subroutine vcycle_fas(fl)
      integer, intent(in) :: fl

      integer :: l, ip

      call push_sub('poisson_multigrid.vcycle_fas')

      do l = fl, cl
        ! store the initial approximation
        np = gr%mgrid%level(l)%mesh%np
        call lalg_copy(np, phi%level(l)%p, phi_ini%level(l)%p)

        ! presmoothing
        call multigrid_relax(this, gr%mgrid%level(l)%mesh, gr%mgrid%level(l)%der, &
          phi%level(l)%p, tau%level(l)%p , this%presteps)

        if (l /= cl ) then
          ! transfer of the current solution
          call dmultigrid_fine2coarse(gr%mgrid, l+1, phi%level(l)%p, phi%level(l+1)%p, this%restriction_method)

          ! error calculation
          call dderivatives_lapl(gr%mgrid%level(l)%der, phi%level(l)%p, err%level(l)%p)
          forall(ip = 1:np) err%level(l)%p(ip) = err%level(l)%p(ip) - tau%level(l)%p(ip)

          ! transfer error to coarser grid
          call dmultigrid_fine2coarse(gr%mgrid, l+1, err%level(l)%p, tau%level(l+1)%p, this%restriction_method)

          ! the other part of the error
          call dderivatives_lapl(gr%mgrid%level(l + 1)%der, phi%level(l + 1)%p, err%level(l + 1)%p)
          np = gr%mgrid%level(l + 1)%mesh%np
          forall(ip = 1:np) tau%level(l + 1)%p(ip) = err%level(l + 1)%p(ip) - tau%level(l + 1)%p(ip)
        end if
      end do

      do l = cl, fl, -1
       
        ! postsmoothing
        call multigrid_relax(this, gr%mgrid%level(l)%mesh, gr%mgrid%level(l)%der, &
          phi%level(l)%p, tau%level(l)%p, this%poststeps)

        if(l /= fl) then
          ! calculate correction as the diference with the
          ! original grid in this level
          np = gr%mgrid%level(l)%mesh%np
          forall(ip = 1:np) phi%level(l)%p(ip) = phi%level(l)%p(ip) - phi_ini%level(l)%p(ip)

          ! transfer correction to finer level
          call dmultigrid_coarse2fine(gr%mgrid%level(l)%tt, gr%mgrid%level(l)%der, &
            gr%mgrid%level(l)%mesh, gr%mgrid%level(l - 1)%mesh, phi%level(l)%p, err%level(l - 1)%p)

          np = gr%mgrid%level(l-1)%mesh%np
          forall(ip = 1:np) phi%level(l - 1)%p(ip) = phi%level(l - 1)%p(ip) + err%level(l - 1)%p(ip)
        end if
      end do

      call pop_sub()
    end subroutine vcycle_fas

  end subroutine poisson_multigrid_solver

  ! ---------------------------------------------------------
  subroutine multigrid_relax(this, m, der, pot, rho, steps)
    type(mg_solver_t),   intent(in)    :: this
    type(mesh_t),        intent(in)    :: m
    type(derivatives_t), intent(inout) :: der
    FLOAT,               intent(inout) :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    integer,             intent(in)    :: steps

    integer :: t
    integer :: i, n
    FLOAT   :: point_lap, factor
    FLOAT, allocatable :: lpot(:), ldiag(:)
    type(profile_t), save :: prof

    call push_sub('poisson_multigrid.multigrid_relax')
    call profiling_in(prof, "MG_GAUSS_SEIDEL")

    select case(this%relaxation_method)

    case(GAUSS_SEIDEL)

      factor = CNST(-1.0)/der%lapl%w_re(der%lapl%stencil%center, 1)*this%relax_factor

      do t = 1, steps

        call dderivatives_set_bc(der, pot)

#ifdef HAVE_MPI
        if(m%parallel_in_domains) call dvec_ghost_update(m%vp, pot)
#endif

        n = der%lapl%stencil%size

        if(der%lapl%const_w) then
          call dgauss_seidel(der%lapl%stencil%size, der%lapl%w_re(1, 1), der%lapl%nri, &
            der%lapl%ri(1, 1), der%lapl%rimap_inv(1), der%lapl%rimap_inv(2),        &
            factor, pot(1), rho(1))
        else
          do i = 1, m%np
            point_lap = sum(der%lapl%w_re(1:n,i)*pot(der%lapl%i(1:n,i)))
            pot(i) = pot(i) - CNST(0.7)/der%lapl%w_re(der%lapl%stencil%center,i)*(point_lap-rho(i))
          end do
        end if
      end do
      call profiling_count_operations(m%np*(steps + 1)*(2*n + 3))

    case(GAUSS_JACOBI)

      SAFE_ALLOCATE( lpot(1:m%np))
      SAFE_ALLOCATE(ldiag(1:m%np))

      call derivatives_lapl_diag(der, ldiag)

      do t = 1, steps
        call dderivatives_lapl(der, pot, lpot)
        pot(1:m%np) = pot(1:m%np) - this%relax_factor/ldiag(1:m%np)*(lpot(1:m%np) - rho(1:m%np)) 
      end do

      SAFE_DEALLOCATE_A(ldiag)
      SAFE_DEALLOCATE_A(lpot)

    end select

    call profiling_out(prof)
    call pop_sub()

  end subroutine multigrid_relax

end module poisson_multigrid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
