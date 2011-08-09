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
  use boundaries_m
  use datasets_m
  use derivatives_m
  use global_m
  use gridhier_m
  use lalg_basic_m
  use parser_m
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
  subroutine poisson_multigrid_init(this, mesh, ml, thr)
    type(mg_solver_t), intent(inout) :: this
    type(mesh_t),      intent(inout) :: mesh
    integer,           intent(in)    :: ml
    FLOAT,             intent(in)    :: thr

    PUSH_SUB(poisson_multigrid_init)

    call poisson_corrections_init(this%corrector, ml, mesh)

    this%threshold = thr

    !%Variable PoissonSolverMGPresmoothingSteps
    !%Type integer
    !%Default 1
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of Gauss-Seidel smoothing steps before coarse-level
    !% correction in the multigrid Poisson solver.
    !%End
    call parse_integer(datasets_check('PoissonSolverMGPresmoothingSteps'), 1, this%presteps)

    !%Variable PoissonSolverMGPostsmoothingSteps
    !%Type integer
    !%Default 4
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Number of Gauss-Seidel smoothing steps after coarse-level
    !% correction in the multigrid Poisson solver.
    !%End
    call parse_integer(datasets_check('PoissonSolverMGPostsmoothingSteps'), 4, this%poststeps)

    !%Variable PoissonSolverMGMaxCycles
    !%Type integer
    !%Default 60
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Maximum number of multigrid cycles that are performed if
    !% convergence is not achieved.
    !%End
    call parse_integer(datasets_check('PoissonSolverMGMaxCycles'), 60, this%maxcycles)

    !%Variable PoissonSolverMGRestrictionMethod
    !%Type integer
    !%Default fullweight
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Method used from fine-to-coarse grid transfer.
    !%Option injection 1
    !% Injection
    !%Option fullweight 2
    !% Fullweight restriction
    !%End
    call parse_integer(datasets_check('PoissonSolverMGRestrictionMethod'), 2, this%restriction_method)
    if(.not.varinfo_valid_option('PoissonSolverMGRestrictionMethod', this%restriction_method)) &
       call input_error('PoissonSolverMGRestrictionMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRestrictionMethod", this%restriction_method)

    !%Variable PoissonSolverMGRelaxationMethod
    !%Type integer
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Method used to solve the linear system approximately in each grid for the
    !% multigrid procedure that solves Poisson equation. Default is <tt>gauss_seidel</tt>,
    !% unless curvilinear coordinates are used, in which case the default is <tt>gauss_jacobi</tt>.
    !%Option gauss_seidel 1
    !% Gauss-Seidel.
    !%Option gauss_jacobi 2
    !% Gauss-Jacobi.
    !%Option gauss_jacobi2 3
    !% Alternative implementation of Gauss-Jacobi.
    !%End
    if ( mesh%use_curvilinear ) then
      call parse_integer(datasets_check('PoissonSolverMGRelaxationMethod'), GAUSS_JACOBI, this%relaxation_method)
    else
      call parse_integer(datasets_check('PoissonSolverMGRelaxationMethod'), GAUSS_SEIDEL, this%relaxation_method)
    end if

    if(.not.varinfo_valid_option('PoissonSolverMGRelaxationMethod', this%relaxation_method)) &
      call input_error('PoissonSolverMGRelaxationMethod')
    call messages_print_var_option(stdout, "PoissonSolverMGRelaxationMethod", this%relaxation_method)

    !%Variable PoissonSolverMGRelaxationFactor
    !%Type float
    !%Section Hamiltonian::Poisson::Multigrid
    !%Description
    !% Relaxation factor of the relaxation operator used for the
    !% multigrid method. This is mainly for debugging,
    !% since overrelaxation does not help in a multigrid scheme.
    !% The default is 1.0, except 0.6666 for the <tt>gauss_jacobi</tt> method.
    !%End
    if ( this%relaxation_method == GAUSS_JACOBI) then
      call parse_float(datasets_check('PoissonSolverMGRelaxationFactor'), CNST(0.6666), this%relax_factor )
    else
      call parse_float(datasets_check('PoissonSolverMGRelaxationFactor'), M_ONE, this%relax_factor)
    end if

    POP_SUB(poisson_multigrid_init)
  end subroutine poisson_multigrid_init


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_end(this)
    type(mg_solver_t), intent(inout) :: this

    PUSH_SUB(poisson_multigrid_end)

    call poisson_corrections_end(this%corrector)

    POP_SUB(poisson_multigrid_end)
  end subroutine poisson_multigrid_end


  ! ---------------------------------------------------------
  subroutine poisson_multigrid_solver(this, base_der, pot, rho)
    type(mg_solver_t),           intent(in)    :: this
    type(derivatives_t), target, intent(inout) :: base_der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)

    integer :: lev, cl, t, np, curr_l, ip
    FLOAT :: res, fmg_threshold = CNST(0.1)
    FLOAT, allocatable :: vh_correction(:)

    type(dgridhier_t) :: phi
    type(dgridhier_t) :: phi_ini
    type(dgridhier_t) :: tau
    type(dgridhier_t) :: err

    type(derivatives_t), pointer :: der

    PUSH_SUB(poisson_multigrid_solver)

    ! correction for treating boundaries
    SAFE_ALLOCATE(vh_correction(1:base_der%mesh%np_part))
    
    call gridhier_init(phi, base_der, np_part_size = .true.)
    call gridhier_init(phi_ini, base_der, np_part_size = .false.)
    call gridhier_init(tau, base_der, np_part_size = .true.)
    call gridhier_init(err, base_der, np_part_size = .true.)

    call correct_rho(this%corrector, base_der, rho, tau%level(0)%p, vh_correction)
    call lalg_scal(base_der%mesh%np, -M_FOUR*M_PI, tau%level(0)%p)

    forall (ip = 1:base_der%mesh%np) phi%level(0)%p(ip) = pot(ip) - vh_correction(ip)

    cl = multigrid_number_of_levels(base_der)

    curr_l = 0
    der => base_der

    do t = 0, this%maxcycles

      if(t > 0) then ! during the first cycle only the residue is calculated
        call vcycle_cs(curr_l, der)
      end if

      res = residue(der, phi%level(curr_l)%p, tau%level(curr_l)%p, err%level(curr_l)%p)
      
      if(in_debug_mode) then
        write(message(1), '(a,i5,a,i5,a,e12.6)') "Multigrid: base level ", curr_l, " iter ", t, " res ", res
        call messages_info(1)
      end if

      if(res < this%threshold) then
        if(curr_l > 0 ) then
          call dmultigrid_coarse2fine(der%to_finer, der, der%finer%mesh,phi%level(curr_l)%p, phi%level(curr_l - 1)%p, order = 2)
          curr_l = curr_l - 1
          der => der%finer
        else
          exit
        end if
      end if

      if(0 == t .and. res > fmg_threshold) then
        curr_l = max(cl - 1, 0)

        do lev = 0, curr_l - 1
          call dmultigrid_fine2coarse(der%to_coarser, der, der%coarser%mesh, &
            tau%level(lev)%p, tau%level(lev + 1)%p, this%restriction_method)
          der => der%coarser
        end do
      end  if

    end do

    if(res >= this%threshold) then
      message(1) = 'Multigrid Poisson solver did not converge.'
      write(message(2), '(a,e14.6)') '  Res = ', res
      call messages_warning(2)
    end if

    forall (ip = 1:base_der%mesh%np) pot(ip) = phi%level(0)%p(ip) + vh_correction(ip)

    SAFE_DEALLOCATE_A(vh_correction)

    call gridhier_end(phi)
    call gridhier_end(tau)
    call gridhier_end(err)
    call gridhier_end(phi_ini)

    POP_SUB(poisson_multigrid_solver)

  contains

    ! ---------------------------------------------------------
    FLOAT function residue(der, phi, rho, tmp) result(res)
      type(derivatives_t), intent(in) :: der
      FLOAT,               intent(inout) :: phi(:)
      FLOAT,               intent(in)    :: rho(:)
      FLOAT,               intent(inout) :: tmp(:)

      integer :: ip

      PUSH_SUB(residue)

      call dderivatives_lapl(der, phi, tmp)

      forall (ip = 1:der%mesh%np) tmp(ip) = tmp(ip) - rho(ip)

      res = dmf_nrm2(der%mesh, tmp)

      POP_SUB(residue)
    end function residue

    ! ---------------------------------------------------------
    subroutine vcycle_cs(fl, base_der)
      integer, intent(in) :: fl
      type(derivatives_t), pointer :: base_der

      type(derivatives_t), pointer :: der
      integer :: level, ip

      PUSH_SUB(vcycle_fas)

      der => base_der
      do level = fl, cl
        ! der points to level
        if(level /= fl) phi%level(level)%p(1:der%mesh%np) = M_ZERO

        ! presmoothing
        call multigrid_relax(this, der%mesh, der, phi%level(level)%p, tau%level(level)%p , this%presteps)

        if (level /= cl ) then
          ! error calculation
          call dderivatives_lapl(der, phi%level(level)%p, err%level(level)%p)
          forall(ip = 1:der%mesh%np) err%level(level)%p(ip) = tau%level(level)%p(ip) - err%level(level)%p(ip)

          ! transfer error as the source in the coarser grid
          call dmultigrid_fine2coarse(der%to_coarser, der, der%coarser%mesh, &
            err%level(level)%p, tau%level(level+1)%p, this%restriction_method)

          der => der%coarser
        end if
      end do

      do level = cl, fl, -1
        ! postsmoothing
        call multigrid_relax(this, der%mesh, der, phi%level(level)%p, tau%level(level)%p, this%poststeps)

        if(level /= fl) then
          ! transfer correction to finer level
          call dmultigrid_coarse2fine(der%to_finer, der, der%finer%mesh, phi%level(level)%p, err%level(level-1)%p)

          np = der%finer%mesh%np

          forall(ip = 1:np) phi%level(level - 1)%p(ip) = phi%level(level - 1)%p(ip) + err%level(level - 1)%p(ip)
          der => der%finer
        end if

      end do

      POP_SUB(vcycle_fas)
    end subroutine vcycle_cs

  end subroutine poisson_multigrid_solver

  ! ---------------------------------------------------------
  subroutine multigrid_relax(this, mesh, der, pot, rho, steps)
    type(mg_solver_t),   intent(in)    :: this
    type(mesh_t),        intent(in)    :: mesh
    type(derivatives_t), intent(inout) :: der
    FLOAT,               intent(inout) :: pot(:)
    FLOAT,               intent(in)    :: rho(:)
    integer,             intent(in)    :: steps

    integer :: istep
    integer :: ip, nn
    FLOAT   :: point_lap, factor
    FLOAT, allocatable :: lpot(:), ldiag(:)
    type(profile_t), save :: prof

    PUSH_SUB(multigrid_relax)
    call profiling_in(prof, "MG_GAUSS_SEIDEL")

    select case(this%relaxation_method)

    case(GAUSS_SEIDEL)

      factor = CNST(-1.0)/der%lapl%w_re(der%lapl%stencil%center, 1)*this%relax_factor

      do istep = 1, steps

        call dderivatives_set_bc(der, pot)

#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) call dvec_ghost_update(mesh%vp, pot)
#endif

        nn = der%lapl%stencil%size

        if(der%lapl%const_w) then
          call dgauss_seidel(der%lapl%stencil%size, der%lapl%w_re(1, 1), der%lapl%nri, &
            der%lapl%ri(1, 1), der%lapl%rimap_inv(1), der%lapl%rimap_inv(2),        &
            factor, pot(1), rho(1))
        else
          do ip = 1, mesh%np
            point_lap = sum(der%lapl%w_re(1:nn, ip)*pot(der%lapl%i(1:nn, ip)))
            pot(ip) = pot(ip) - CNST(0.7)/der%lapl%w_re(der%lapl%stencil%center, ip)*(point_lap-rho(ip))
          end do
        end if
      end do
      call profiling_count_operations(mesh%np*(steps + 1)*(2*nn + 3))

    case(GAUSS_JACOBI)

      SAFE_ALLOCATE( lpot(1:mesh%np))
      SAFE_ALLOCATE(ldiag(1:mesh%np))

      call derivatives_lapl_diag(der, ldiag)

      do istep = 1, steps
        call dderivatives_lapl(der, pot, lpot)
        pot(1:mesh%np) = pot(1:mesh%np) - this%relax_factor/ldiag(1:mesh%np)*(lpot(1:mesh%np) - rho(1:mesh%np)) 
      end do

      SAFE_DEALLOCATE_A(ldiag)
      SAFE_DEALLOCATE_A(lpot)

    end select

    call profiling_out(prof)
    POP_SUB(multigrid_relax)

  end subroutine multigrid_relax

end module poisson_multigrid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
