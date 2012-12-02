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
  use lalg_basic_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use multigrid_m
  use operate_f_m
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
    call parse_integer(datasets_check('PoissonSolverMGMaxCycles'), 50, this%maxcycles)

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
  subroutine poisson_multigrid_solver(this, der, pot, rho)
    type(mg_solver_t),           intent(in)    :: this
    type(derivatives_t), target, intent(inout) :: der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)

    integer :: iter, ip
    FLOAT :: resnorm
    FLOAT, allocatable :: vh_correction(:), res(:), cor(:), err(:)

    PUSH_SUB(poisson_multigrid_solver)

    ! correction for treating boundaries
    SAFE_ALLOCATE(vh_correction(1:der%mesh%np_part))
    SAFE_ALLOCATE(res(1:der%mesh%np))
    SAFE_ALLOCATE(cor(1:der%mesh%np_part))
    SAFE_ALLOCATE(err(1:der%mesh%np))
    
    call correct_rho(this%corrector, der, rho, res, vh_correction)
    call lalg_scal(der%mesh%np, -M_FOUR*M_PI, res)

    forall (ip = 1:der%mesh%np) cor(ip) = pot(ip) - vh_correction(ip)

    do iter = 1, this%maxcycles

      call poisson_multigrid_cycle(this, der, cor, res)
      call dderivatives_lapl(der, cor, err)
      forall (ip = 1:der%mesh%np) err(ip) = res(ip) - err(ip)
      resnorm =  dmf_nrm2(der%mesh, err)

      if(resnorm < this%threshold) exit

      if(in_debug_mode) then
        write(message(1), '(a,i5,a,e13.6)') "Multigrid: base level: iter ", iter, " res ", resnorm
        call messages_info(1)
      end if

    end do

    if(resnorm >= this%threshold) then
      message(1) = 'Multigrid Poisson solver did not converge.'
      write(message(2), '(a,e14.6)') '  Res = ', resnorm
      call messages_warning(2)
    end if

    forall (ip = 1:der%mesh%np) pot(ip) = cor(ip) + vh_correction(ip)

    SAFE_DEALLOCATE_A(vh_correction)
    SAFE_DEALLOCATE_A(res)
    SAFE_DEALLOCATE_A(cor)
    SAFE_DEALLOCATE_A(err)

    POP_SUB(poisson_multigrid_solver)
  end subroutine poisson_multigrid_solver

  ! ---------------------------------------------------------

  recursive subroutine poisson_multigrid_cycle(this, der, pot, rho)
    type(mg_solver_t),           intent(in)    :: this
    type(derivatives_t),         intent(inout) :: der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)
    
    integer :: ip, iter
    real(8) :: resnorm
    real(8), allocatable :: res(:), cres(:), cor(:), ccor(:)
    
    PUSH_SUB(poisson_multigrid_cycle)

    SAFE_ALLOCATE(res(1:der%mesh%np_part))

    if(associated(der%coarser)) then
      SAFE_ALLOCATE(cor(1:der%mesh%np_part))
      SAFE_ALLOCATE(cres(1:der%coarser%mesh%np_part))
      SAFE_ALLOCATE(ccor(1:der%coarser%mesh%np_part))
     
      call multigrid_relax(this, der%mesh, der, pot, rho, this%presteps)
      
      call dderivatives_lapl(der, pot, res)
      forall (ip = 1:der%mesh%np) res(ip) = rho(ip) - res(ip)
      
      call dmultigrid_fine2coarse(der%to_coarser, der, der%coarser%mesh, res, cres, this%restriction_method)
      
      ccor = M_ZERO
      call poisson_multigrid_cycle(this, der%coarser, ccor, cres)

      cor = M_ZERO
      call dmultigrid_coarse2fine(der%to_coarser, der%coarser, der%mesh, ccor, cor)

      forall (ip = 1:der%mesh%np) pot(ip) = pot(ip) + cor(ip)

      SAFE_DEALLOCATE_A(cor)
      SAFE_DEALLOCATE_A(cres)
      SAFE_DEALLOCATE_A(ccor)

      call multigrid_relax(this, der%mesh, der, pot, rho, this%poststeps)

    else    

      do iter = 1, this%maxcycles
        call multigrid_relax(this, der%mesh, der, pot, rho, this%presteps + this%poststeps)
        call dderivatives_lapl(der, pot, res)
        forall (ip = 1:der%mesh%np) res(ip) = rho(ip) - res(ip)
        resnorm = dmf_nrm2(der%mesh, res)
        if(resnorm < this%threshold) exit
      end do

    end if

    SAFE_DEALLOCATE_A(res)

    POP_SUB(poisson_multigrid_cycle)

  end subroutine poisson_multigrid_cycle

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
