!! Copyright (C) 2002-2007 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module geom_opt_m
  use datasets_m
  use energy_m
  use external_pot_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use loct_parser_m
  use loct_m
  use loct_math_m
  use mesh_m
  use messages_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use states_calc_m
  use system_m
  use units_m
  use v_ks_m
  use varinfo_m
  use species_pot_m
  use lcao_m

  implicit none

  private
  public :: geom_opt_run

  type geom_opt_t
    integer  :: method
    FLOAT    :: step
    FLOAT    :: tolgrad
    FLOAT    :: toldr
    integer  :: max_iter
    integer  :: what2minimize

    ! shortcuts
    type(scf_t)                  :: scfv
    type(geometry_t),    pointer :: geo
    type(hamiltonian_t), pointer :: hm
    type(system_t),      pointer :: syst
    type(mesh_t),        pointer :: m
    type(states_t),      pointer :: st
  end type geom_opt_t

  type(geom_opt_t) :: g_opt

  integer, parameter :: &
    MINWHAT_ENERGY = 1, &
    MINWHAT_FORCES = 2

contains

  ! ---------------------------------------------------------
  subroutine geom_opt_run(sys, hm, fromscratch)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm
    logical,                     intent(inout) :: fromscratch

    integer :: i, ierr, lcao_start, lcao_start_default
    real(8), allocatable :: x(:)
    type(lcao_t) :: lcao
    real(8) :: energy

    call init_()
    
    ! load wave-functions
    if(.not. fromscratch) then
      call restart_read(trim(restart_dir)//'gs', sys%st, sys%gr, sys%geo, ierr)
      if(ierr .ne. 0) then
        message(1) = "Could not load wave-functions: Starting from scratch"
        call write_warning(1)
        fromscratch = .true.
      end if
    end if

    if(fromscratch) then

      ! Randomly generate the initial wave-functions
      call states_generate_random(sys%st, sys%gr%mesh)
      call states_orthogonalize(sys%st, sys%gr%mesh)

      ! We do not compute the density from the random wave-functions. 
      ! Instead, we try to get a better guess for the density
      call guess_density(sys%gr%mesh, sys%gr%sb, sys%geo, sys%st%qtot, sys%st%d%nspin, &
           sys%st%d%spin_channels, sys%st%rho)

      ! setup Hamiltonian (we do not call system_h_setup here because we do not want to
      ! overwrite the guess density)
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call v_ks_calc(sys%gr, sys%ks, hm, sys%st, calc_eigenval=.true.) ! get potentials
      call states_fermi(sys%st, sys%gr%mesh)                                ! occupations
      call total_energy(hm, sys%gr, sys%st, -1)

      lcao_start_default = LCAO_START_FULL
      if(sys%geo%only_user_def) lcao_start_default = LCAO_START_NONE

      call loct_parse_int(datasets_check('LCAOStart'), lcao_start_default, lcao_start)
      if(.not.varinfo_valid_option('LCAOStart', lcao_start)) call input_error('LCAOStart')
      call messages_print_var_option(stdout, 'LCAOStart', lcao_start)

      if (lcao_start > LCAO_START_NONE) then

        call lcao_init(lcao, sys%gr, sys%geo, sys%st)

        if(lcao_is_available(lcao)) then
          write(message(1),'(a,i4,a)') 'Info: Performing initial LCAO calculation with ', &
               lcao_num_orbitals(lcao), ' orbitals.'
          call write_info(1)
          
          call lcao_wf(lcao, sys%st, sys%gr, sys%geo, hm)
          call lcao_end(lcao)

          !Just populate again the states, so that the eigenvalues are properly written
          call states_fermi(sys%st, sys%gr%mesh)
          call states_write_eigenvalues(stdout, sys%st%nst, sys%st, sys%gr%sb)

          ! Update the density and the Hamiltonian
          if (lcao_start == LCAO_START_FULL) call system_h_setup(sys, hm)

        end if

      end if

    else

      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call write_info(1)
      call system_h_setup(sys, hm)

    end if

    call scf_init(g_opt%scfv, sys%gr, sys%geo, sys%st, hm)

    !Initial point
    ALLOCATE(x(3*g_opt%geo%natoms), 3*g_opt%geo%natoms)
    do i = 0, g_opt%geo%natoms - 1
      x(3*i + 1) = g_opt%geo%atom(i + 1)%x(1)
      x(3*i + 2) = g_opt%geo%atom(i + 1)%x(2)
      x(3*i + 3) = g_opt%geo%atom(i + 1)%x(3)
    end do

    !Minimize
    select case(g_opt%method)
    case(MINMETHOD_NMSIMPLEX)
      ierr = loct_minimize_direct(g_opt%method, 3*g_opt%geo%natoms, x(1), real(g_opt%step, 8),&
           real(g_opt%toldr, 8), g_opt%max_iter, &
           calc_point_ng, write_iter_info_ng, energy)
    case default
      ierr = loct_minimize(g_opt%method, 3*g_opt%geo%natoms, x(1), real(g_opt%step, 8),&
           real(g_opt%tolgrad, 8), real(g_opt%toldr, 8), g_opt%max_iter, &
           calc_point, write_iter_info, energy)
    end select

    if (ierr /= 0) then
      message(1) = "Error occurred during the GSL minimization procedure:"
      call loct_strerror(ierr, message(2))
      call write_fatal(2)
    end if

    ! print out geometry
    do i = 0, g_opt%geo%natoms - 1
      g_opt%geo%atom(i+1)%x(1) = x(3*i + 1)
      g_opt%geo%atom(i+1)%x(2) = x(3*i + 2)
      g_opt%geo%atom(i+1)%x(3) = x(3*i + 3)
    end do
    call atom_write_xyz(".", "min", g_opt%geo, sys%gr%mesh%sb%dim)

    SAFE_DEALLOCATE_A(x)
    call scf_end(g_opt%scfv)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('geom_opt.geom_opt_run')

      call states_allocate_wfns(sys%st, sys%gr%mesh)

      ! shortcuts
      g_opt%m      => sys%gr%mesh
      g_opt%geo    => sys%geo
      g_opt%st     => sys%st
      g_opt%hm     => hm
      g_opt%syst   => sys

      !%Variable GOMethod
      !%Type integer
      !%Default steep
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Method by which the minimization is performed.
      !%Option steep 1
      !% simple steepest descent.
      !%Option cg_fr 2
      !% Fletcher-Reeves conjugate gradient algorithm. The
      !% conjugate gradient algorithm proceeds as a succession of line
      !% minimizations. The sequence of search directions is used to build
      !% up an approximation to the curvature of the function in the
      !% neighborhood of the minimum. 
      !%Option cg_pr 3
      !% Polak-Ribiere conjugate gradient algorithm.
      !%Option cg_bfgs 4
      !% Vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) conjugate gradient algorithm.
      !% It is a quasi-Newton method which builds up an approximation to the second 
      !% derivatives of the function f using the difference between successive gradient
      !% vectors.  By combining the first and second derivatives the algorithm is able 
      !% to take Newton-type steps towards the function minimum, assuming quadratic 
      !% behavior in that region.
      !%Option cg_bfgs2 5
      !% The bfgs2 version of this minimizer is the most efficient version available, 
      !% and is a faithful implementation of the line minimization scheme described in 
      !% Fletcher, _Practical Methods of Optimization_, Algorithms 2.6.2 and 2.6.4.
      !%Option simplex 6
      !% This is experimental, and in fact, *not* recommended unless you just want to
      !% fool around. It is the Nead-Melder simplex algorithm, as implemented in the
      !% GNU Scientific Library (GSL). It does not make use of the gradients (i.e., the
      !% forces) which makes it more inefficient than other schemes. It is included here
      !% for completeness, since it is free.
      !%End
      call loct_parse_int(datasets_check('GOMethod'), MINMETHOD_STEEPEST_DESCENT, g_opt%method)
      if(.not.varinfo_valid_option('GOMethod', g_opt%method)) call input_error('GOMethod')
      call messages_print_var_option(stdout, "GOMethod", g_opt%method)

      !%Variable GOTolerance
      !%Type float
      !%Default 0.001 a.u.
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Convergence criterion to stop the minimization. In units of force; minimization
      !% is stopped when all forces on ions are smaller.
      !% Used in conjunction with GOMinimumMove. If GOTolerance = 0, this criterion is ignored.
      !%End
      call loct_parse_float(datasets_check('GOTolerance'), CNST(0.001)/units_inp%force%factor, g_opt%tolgrad)
      g_opt%tolgrad = g_opt%tolgrad*units_inp%force%factor
      
      !%Variable GOMinimumMove
      !%Type float
      !%Default 0.0 a.u.
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Convergence criterion to stop the minimization. In units of length; minimization
      !% is stopped when all species coordinates change less than GOMinimumMove.
      !% Used in conjunction with GOTolerance. If GOMinimumMove = 0, this criterion is ignored.
      !%
      !% Note that if you use GOMethod = simplex, then you must supply a non-zero GOMinimumMove.
      !%End
      call loct_parse_float(datasets_check('GOMinimumMove'), CNST(0.0)/units_inp%length%factor, g_opt%toldr)
      g_opt%toldr = g_opt%toldr*units_inp%length%factor
      if(g_opt%method == MINMETHOD_NMSIMPLEX .and. g_opt%toldr <= M_ZERO) call input_error('GOMinimumMove')
      
      !%Variable GOStep
      !%Type float
      !%Default 0.5
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Initial step for the geometry optimizer.
      !%End
      call loct_parse_float(datasets_check('GOStep'), M_HALF, g_opt%step) ! WARNING: in some weird units

      !%Variable GOMaxIter
      !%Type integer
      !%Default 200
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Even if previous convergence criterion is not satisfied, minimization will stop
      !% after this number of iterations.
      !%End
      call loct_parse_int(datasets_check('GOMaxIter'), 200, g_opt%max_iter)
      if(g_opt%max_iter <= 0) then
        message(1) = "GoMaxIter has to be larger than 0"
        call write_fatal(1)
      end if

      !%Variable GOWhat2Minimize
      !%Type integer
      !%Default minimize_energy
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% This rather esoteric option allows one to choose which objective function
      !% to minimize during a geometry minimization. The use of this variable may
      !% lead to inconsistencies, so please make sure you know what you are doing!
      !%Option minimize_energy 1
      !% Use the total energy as objective function
      !%Option minimize_forces 2
      !% Use <math>\sqrt{\sum |f_i|^2}</math> as objective function.
      !% Note that in this case one still uses the forces as the gradient of the objective function.
      !% This is, of course, inconsistent, and may lead to very strange behavior.
      !%End
      call loct_parse_int(datasets_check('GOWhat2Minimize'), MINWHAT_ENERGY, g_opt%what2minimize)
      if(.not.varinfo_valid_option('GOWhat2Minimize', g_opt%what2minimize)) call input_error('GOWhat2Minimize')
      call messages_print_var_option(stdout, "GOWhat2Minimize", g_opt%what2minimize)

      call loct_rm("./work-geom.xyz")

      call pop_sub()
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      call states_deallocate_wfns(sys%st)
      nullify(g_opt%m)
      nullify(g_opt%geo)
      nullify(g_opt%st)
      nullify(g_opt%hm)
      nullify(g_opt%syst)
    end subroutine end_

  end subroutine geom_opt_run


  ! ---------------------------------------------------------
  subroutine calc_point(n, x, f, getgrad, df)
    integer,           intent(in)  :: n
    real(8),           intent(in)  :: x(n)
    real(8),           intent(out) :: f
    integer,           intent(in)  :: getgrad
    real(8), optional, intent(out) :: df(n)
    
    integer :: i

    call push_sub("geom_opt.calc_point")

    do i = 0, g_opt%geo%natoms - 1
      g_opt%geo%atom(i+1)%x(1) = x(3*i + 1)
      g_opt%geo%atom(i+1)%x(2) = x(3*i + 2)
      g_opt%geo%atom(i+1)%x(3) = x(3*i + 3)
    end do

    call atom_write_xyz(".", "work-geom", g_opt%geo, g_opt%syst%gr%mesh%sb%dim, append=.true.)

    call epot_generate(g_opt%hm%ep, g_opt%syst%gr, g_opt%syst%geo, g_opt%syst%st)
    call states_calc_dens(g_opt%st, g_opt%m%np)
    call v_ks_calc(g_opt%syst%gr, g_opt%syst%ks, g_opt%hm, g_opt%st, calc_eigenval=.true.)
    call total_energy(g_opt%hm, g_opt%syst%gr, g_opt%st, -1)

    ! do scf calculation
    call scf_run(g_opt%scfv, g_opt%syst%gr, g_opt%geo, g_opt%st, &
      g_opt%syst%ks, g_opt%hm, g_opt%syst%outp, verbosity = VERB_COMPACT)

    ! store results
    if(getgrad .eq. 1) then
      do i = 0, g_opt%geo%natoms - 1
        df(3*i + 1) = -g_opt%geo%atom(i+1)%f(1)
        df(3*i + 2) = -g_opt%geo%atom(i+1)%f(2)
        df(3*i + 3) = -g_opt%geo%atom(i+1)%f(3)
      end do
    end if

    if(g_opt%what2minimize == MINWHAT_FORCES) then
      f = M_ZERO
      do i = 1, g_opt%geo%natoms
        f = f + sum(g_opt%geo%atom(i)%f(:)**2)
      end do
      f = sqrt(f)
    else
      f = g_opt%hm%etot
    end if

    call pop_sub()
  end subroutine calc_point


  ! ---------------------------------------------------------
  ! Same as calc_point, but without the gradients.
  subroutine calc_point_ng(n, x, f)
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    real(8), intent(out) :: f
    
    call calc_point(n, x, f, getgrad=0)

  end subroutine calc_point_ng


  ! ---------------------------------------------------------
  subroutine write_iter_info(geom_iter, n, energy, maxdx, maxdf, x)
    integer,      intent(in) :: geom_iter, n
    REAL_DOUBLE,  intent(in) :: energy, maxdx, maxdf
    REAL_DOUBLE,  intent(in) :: x(n)

    integer :: i
    character(len=256) :: c_geom_iter, title

    call push_sub("geom_opt.write_iter_info")
    
    write(c_geom_iter, '(a,i4.4)') "go.", geom_iter
    write(title, '(f16.10)') energy/units_out%energy%factor
    call atom_write_xyz("geom", trim(c_geom_iter), g_opt%geo, g_opt%syst%gr%mesh%sb%dim, comment=trim(title))

    do i = 0, g_opt%geo%natoms - 1
      g_opt%geo%atom(i+1)%x(1) = x(3*i + 1)
      g_opt%geo%atom(i+1)%x(2) = x(3*i + 2)
      g_opt%geo%atom(i+1)%x(3) = x(3*i + 3)
    end do

    message(1) = ""
    message(2) = ""
    message(3) = "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(message(4),'("+++++++++++++++++++++ MINIMIZATION ITER #: ",I4," ++++++++++++++++++++++")') geom_iter
    write(message(5), '(2x,2(a,f16.10))') "Energy    = ", energy/units_out%energy%factor
    if(maxdf > M_ZERO) then
      write(message(6),'(2X,2(a,f16.10))')  "Max force = ", maxdf/units_out%force%factor
    end if
    write(message(7),'(2X,2(a,f16.10))')  "Max dr    = ", maxdx/units_out%length%factor
    message(8) = message(3)
    message(9) = message(3)
    message(10) = ""
    message(11) = ""
    call write_info(11)

    call pop_sub()
  end subroutine write_iter_info

  ! ---------------------------------------------------------
  ! Same as write_iter_info, but without the gradients.
  subroutine write_iter_info_ng(geom_iter, n, energy, maxdx, x)
    integer, intent(in) :: geom_iter, n
    REAL_DOUBLE, intent(in) :: energy, maxdx
    REAL_DOUBLE, intent(in) :: x(n)

    call write_iter_info(geom_iter, n, energy, maxdx, real(-M_ONE, 8), x)
  end subroutine write_iter_info_ng

end module geom_opt_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
