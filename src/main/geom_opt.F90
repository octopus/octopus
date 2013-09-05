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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module geom_opt_m
  use datasets_m
  use density_m
  use energy_calc_m
  use epot_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use io_m
  use lcao_m
  use loct_m
  use loct_math_m
  use parser_m
  use mesh_m
  use messages_m
  use minimizer_m
  use profiling_m
  use restart_m
  use scf_m
  use simul_box_m
  use species_pot_m
  use states_m
  use states_calc_m
  use system_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use varinfo_m
  use xyz_adjust_m

  implicit none

  private
  public :: geom_opt_run

  type geom_opt_t
    integer  :: method
    FLOAT    :: step
    FLOAT    :: line_tol
    FLOAT    :: tolgrad
    FLOAT    :: toldr
    integer  :: max_iter
    integer  :: what2minimize

    !> shortcuts
    type(scf_t)                  :: scfv
    type(geometry_t),    pointer :: geo
    type(hamiltonian_t), pointer :: hm
    type(system_t),      pointer :: syst
    type(mesh_t),        pointer :: mesh
    type(states_t),      pointer :: st
    integer                      :: dim
    integer                      :: size
    integer                      :: fixed_atom
  end type geom_opt_t

  type(geom_opt_t), save :: g_opt

  integer, parameter :: &
    MINWHAT_ENERGY = 1, &
    MINWHAT_FORCES = 2

contains

  ! ---------------------------------------------------------
  subroutine geom_opt_run(sys, hm, fromscratch)
    type(system_t), target,      intent(inout) :: sys
    type(hamiltonian_t), target, intent(inout) :: hm
    logical,                     intent(inout) :: fromscratch

    integer :: ierr
    REAL_DOUBLE, allocatable :: coords(:)
    REAL_DOUBLE :: energy

    PUSH_SUB(geom_opt_run)

    call init_()
    

    ! load wavefunctions
    if(.not. fromscratch) then
      call restart_read(trim(restart_dir)//GS_DIR, sys%st, sys%gr, ierr)
      if(ierr /= 0) then
        message(1) = "Could not load wavefunctions: Starting from scratch."
        call messages_warning(1)
        fromscratch = .true.
      end if
    end if

    call scf_init(g_opt%scfv, sys%gr, sys%geo, sys%st, hm, conv_force = CNST(1e-8))

    if(fromScratch) then
      call lcao_run(sys, hm, lmm_r = g_opt%scfv%lmm_r)
    else
      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call messages_info(1)
      call system_h_setup(sys, hm)
    end if

    !Initial point
    SAFE_ALLOCATE(coords(1:g_opt%size))
    call to_coords(g_opt, coords)

    if(sys%st%d%pack_states) call states_pack(sys%st)

    !Minimize
    select case(g_opt%method)
    case(MINMETHOD_NMSIMPLEX)
      call minimize_multidim_nograd(g_opt%method, g_opt%size, coords, real(g_opt%step, 8),&
        real(g_opt%toldr, 8), g_opt%max_iter, &
        calc_point_ng, write_iter_info_ng, energy, ierr)
    case default
      call minimize_multidim(g_opt%method, g_opt%size, coords, real(g_opt%step, 8),&
        real(g_opt%line_tol, 8), real(g_opt%tolgrad, 8), real(g_opt%toldr, 8), g_opt%max_iter, &
        calc_point, write_iter_info, energy, ierr)
    end select

    if(ierr == 1025) then
      ! not a GSL error, set by our minimize routines, so we must handle it separately
      message(1) = "Reached maximum number of iterations allowed by GOMaxIter."
      call messages_info(1)
    else if (ierr /= 0) then
      message(1) = "Error occurred during the GSL minimization procedure:"
      call loct_strerror(ierr, message(2))
      call messages_fatal(2)
    end if

    if(sys%st%d%pack_states) call states_unpack(sys%st)
  
    ! print out geometry
    call from_coords(g_opt, coords)
    message(1) = "Writing final coordinates to min.xyz"
    call messages_info(1)
    call geometry_write_xyz(g_opt%geo, './min')

    SAFE_DEALLOCATE_A(coords)
    call scf_end(g_opt%scfv)
    call end_()
    POP_SUB(geom_opt_run)

  contains

    ! ---------------------------------------------------------
    subroutine init_()

      logical :: center, does_exist
      integer :: iter
      character*100 :: filename
      FLOAT :: default_toldr

      PUSH_SUB(geom_opt_run.init_)

      call states_allocate_wfns(sys%st, sys%gr%mesh)

      ! shortcuts
      g_opt%mesh   => sys%gr%mesh
      g_opt%geo    => sys%geo
      g_opt%st     => sys%st
      g_opt%hm     => hm
      g_opt%syst   => sys
      g_opt%dim    =  sys%gr%mesh%sb%dim

      g_opt%size = g_opt%dim*g_opt%geo%natoms

      !%Variable GOCenter
      !%Type logical
      !%Default no
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% (Experimental) If set to yes, Octopus centers the geometry at
      !% every optimization step. It also reduces the degrees of
      !% freedom of the optimization by using the translational
      !% invariance. The default is no.
      !%End
      call parse_logical(datasets_check('GOCenter'), .false.,  center)

      if(center) then
        g_opt%fixed_atom = 1
        g_opt%size = g_opt%size  - g_opt%dim
        call messages_experimental('GOCenter')
      end if
      
      !%Variable GOMethod
      !%Type integer
      !%Default steep
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Method by which the minimization is performed. For more information see the
      !% <a href=http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html>
      !% GSL documentation</a>.
      !%Option steep 1
      !% Simple steepest descent.
      !%Option steep_native -1
      !% (Experimental) Non-gsl implementation of steepest descent.
      !%Option cg_fr 2
      !% Fletcher-Reeves conjugate-gradient algorithm. The
      !% conjugate-gradient algorithm proceeds as a succession of line
      !% minimizations. The sequence of search directions is used to build
      !% up an approximation to the curvature of the function in the
      !% neighborhood of the minimum. 
      !%Option cg_pr 3
      !% Polak-Ribiere conjugate-gradient algorithm.
      !%Option cg_bfgs 4
      !% Vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) conjugate-gradient algorithm.
      !% It is a quasi-Newton method which builds up an approximation to the second 
      !% derivatives of the function <i>f</i> using the difference between successive gradient
      !% vectors.  By combining the first and second derivatives, the algorithm is able 
      !% to take Newton-type steps towards the function minimum, assuming quadratic 
      !% behavior in that region.
      !%Option cg_bfgs2 5
      !% The bfgs2 version of this minimizer is the most efficient version available, 
      !% and is a faithful implementation of the line minimization scheme described in 
      !% Fletcher, <i>Practical Methods of Optimization</i>, Algorithms 2.6.2 and 2.6.4.
      !%Option simplex 6
      !% This is experimental, and in fact, <b>not</b> recommended unless you just want to
      !% fool around. It is the Nead-Melder simplex algorithm, as implemented in the
      !% GNU Scientific Library (GSL). It does not make use of the gradients (<i>i.e.</i>, the
      !% forces) which makes it less efficient than other schemes. It is included here
      !% for completeness, since it is free.
      !%End
      call parse_integer(datasets_check('GOMethod'), MINMETHOD_STEEPEST_DESCENT, g_opt%method)
      if(.not.varinfo_valid_option('GOMethod', g_opt%method)) call input_error('GOMethod')
      call messages_print_var_option(stdout, "GOMethod", g_opt%method)

      !%Variable GOTolerance
      !%Type float
      !%Default 0.001 H/b (0.05 eV/Ang)
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Convergence criterion, for stopping the minimization. In
      !% units of force; minimization is stopped when all forces on
      !% ions are smaller than this criterion, or the
      !% <tt>GOMinimumMove</tt> is satisfied. If <tt>GOTolerance < 0</tt>,
      !% this criterion is ignored.
      !%End
      call parse_float(datasets_check('GOTolerance'), CNST(0.001), g_opt%tolgrad, units_inp%force)
      
      !%Variable GOMinimumMove
      !%Type float
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Convergence criterion, for stopping the minimization. In
      !% units of length; minimization is stopped when the coordinates
      !% of all species change less than <tt>GOMinimumMove</tt>, or the
      !% <tt>GOTolerance</tt> criterion is satisfied.
      !% If <tt>GOMinimumMove < 0</tt>, this criterion is ignored.
      !% Default is -1, except 0.001 b with <tt>GOMethod = simplex</tt>.
      !% Note that if you use <tt>GOMethod = simplex</tt>,
      !% then you must supply a non-zero <tt>GOMinimumMove</tt>.
      !%End
      if(g_opt%method == MINMETHOD_NMSIMPLEX) then
        default_toldr = CNST(0.001)
      else
        default_toldr = -M_ONE
      endif
      call parse_float(datasets_check('GOMinimumMove'), default_toldr, g_opt%toldr, units_inp%length)

      if(g_opt%method == MINMETHOD_NMSIMPLEX .and. g_opt%toldr <= M_ZERO) call input_error('GOMinimumMove')
      
      !%Variable GOStep
      !%Type float
      !%Default 0.5
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Initial step for the geometry optimizer.
      !% WARNING: in some weird units.
      !%End
      call parse_float(datasets_check('GOStep'), M_HALF, g_opt%step)

      !%Variable GOLineTol
      !%Type float
      !%Default 0.1
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Tolerance for line-minimization. Applies only to GSL methods that use the forces.
      !% WARNING: in some weird units.
      !%End
      call parse_float(datasets_check('GOLineTol'), CNST(0.1), g_opt%line_tol)

      !%Variable GOMaxIter
      !%Type integer
      !%Default 200
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Even if the convergence criterion is not satisfied, the minimization will stop
      !% after this number of iterations.
      !%End
      call parse_integer(datasets_check('GOMaxIter'), 200, g_opt%max_iter)
      if(g_opt%max_iter <= 0) then
        message(1) = "GOMaxIter has to be larger than 0"
        call messages_fatal(1)
      end if

      call messages_obsolete_variable('GOWhat2Minimize', 'GOObjective')

      !%Variable GOObjective
      !%Type integer
      !%Default minimize_energy
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% This rather esoteric option allows one to choose which objective function
      !% to minimize during a geometry minimization. The use of this variable may
      !% lead to inconsistencies, so please make sure you know what you are doing!
      !%Option minimize_energy 1
      !% Use the total energy as objective function.
      !%Option minimize_forces 2
      !% Use <math>\sqrt{\sum |f_i|^2}</math> as objective function.
      !% Note that in this case one still uses the forces as the gradient of the objective function.
      !% This is, of course, inconsistent, and may lead to very strange behavior.
      !%End
      call parse_integer(datasets_check('GOObjective'), MINWHAT_ENERGY, g_opt%what2minimize)
      if(.not.varinfo_valid_option('GOObjective', g_opt%what2minimize)) call input_error('GOObjective')
      call messages_print_var_option(stdout, "GOObjective", g_opt%what2minimize)

      call loct_rm("./work-geom.xyz")

      ! clean out old geom/go.XXXX.xyz files. must be consistent with write_iter_info
      iter = 1
      do
        write(filename, '(a,i4.4,a)') "geom/go.", iter, ".xyz"
        inquire(file = trim(filename), exist = does_exist)
        if(does_exist) then
          call loct_rm(trim(filename))
          iter = iter + 1
        else
          exit
        endif
      enddo

      POP_SUB(geom_opt_run.init_)
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(geom_opt_run.end_)

      call states_deallocate_wfns(sys%st)

      nullify(g_opt%mesh)
      nullify(g_opt%geo)
      nullify(g_opt%st)
      nullify(g_opt%hm)
      nullify(g_opt%syst)

      POP_SUB(geom_opt_run.end_)
    end subroutine end_

  end subroutine geom_opt_run


  ! ---------------------------------------------------------
  !> Note: you might think it would be better to change the arguments with '(size)' below to '(:)'.
  !! However, if you do that, then you will (portably) get a segmentation fault.
  subroutine calc_point(size, coords, objective, getgrad, df)
    integer,     intent(in)    :: size
    REAL_DOUBLE, intent(in)    :: coords(size)
    REAL_DOUBLE, intent(inout) :: objective
    integer,     intent(in)    :: getgrad
    REAL_DOUBLE, intent(inout) :: df(size)
    
    integer :: iatom

    PUSH_SUB(calc_point)

    ASSERT(size == g_opt%size)

    call from_coords(g_opt, coords)

    if(g_opt%fixed_atom /= 0) call xyz_adjust_it(g_opt%geo, rotate = .false.)

    call simul_box_atoms_in_box(g_opt%syst%gr%sb, g_opt%geo, warn_if_not = .false., die_if_not = .true.)

    call geometry_write_xyz(g_opt%geo, './work-geom', append = .true.)

    call hamiltonian_epot_generate(g_opt%hm, g_opt%syst%gr, g_opt%geo, g_opt%st)
    call density_calc(g_opt%st, g_opt%syst%gr, g_opt%st%rho)
    call v_ks_calc(g_opt%syst%ks, g_opt%hm, g_opt%st, g_opt%geo, calc_eigenval = .true.)
    call energy_calc_total(g_opt%hm, g_opt%syst%gr, g_opt%st)

    ! do SCF calculation
    call scf_run(g_opt%scfv, g_opt%syst%mc, g_opt%syst%gr, g_opt%geo, g_opt%st, &
      g_opt%syst%ks, g_opt%hm, g_opt%syst%outp, verbosity = VERB_COMPACT)

    ! store results
    if(getgrad  ==  1) call to_grad(g_opt, df)

    if(g_opt%what2minimize == MINWHAT_FORCES) then
      objective = M_ZERO
      do iatom = 1, g_opt%geo%natoms
        objective = objective + sum(g_opt%geo%atom(iatom)%f(1:g_opt%syst%gr%sb%dim)**2)
      end do
      objective = sqrt(objective)
    else
      objective = g_opt%hm%energy%total
    end if

    POP_SUB(calc_point)
  end subroutine calc_point


  ! ---------------------------------------------------------
  !> Same as calc_point, but without the gradients.
  !! No intents here is unfortunately required because the same dummy function will be passed
  !! also to newuoa routines in opt_control, and there the interface has no intents.
  subroutine calc_point_ng(size, coords, objective)
    integer     :: size         !< intent(in)
    REAL_DOUBLE :: coords(size) !< intent(in)
    REAL_DOUBLE :: objective    !< intent(inout)

    integer :: getgrad
    REAL_DOUBLE , allocatable :: df(:)
    
    PUSH_SUB(calc_point_ng)
    
    ASSERT(size == g_opt%size)

    getgrad = 0
    SAFE_ALLOCATE(df(1:size))
    df = M_ZERO

    call calc_point(size, coords, objective, getgrad, df)
    SAFE_DEALLOCATE_A(df)

    POP_SUB(calc_point_ng)
  end subroutine calc_point_ng


  ! ---------------------------------------------------------
  subroutine write_iter_info(geom_iter, size, energy, maxdx, maxdf, coords)
    integer,     intent(in) :: geom_iter
    integer,     intent(in) :: size !< must equal dim * natoms
    REAL_DOUBLE, intent(in) :: energy, maxdx, maxdf
    REAL_DOUBLE, intent(in) :: coords(size)

    character(len=256) :: c_geom_iter, title

    PUSH_SUB(write_iter_info)
    
    write(c_geom_iter, '(a,i4.4)') "go.", geom_iter
    write(title, '(f16.10,2x,a)') units_from_atomic(units_out%energy, energy), trim(units_abbrev(units_out%energy))
    call io_mkdir('geom')
    call geometry_write_xyz(g_opt%geo, 'geom/'//trim(c_geom_iter), comment = trim(title))
    call geometry_write_xyz(g_opt%geo, './last')

    call from_coords(g_opt, coords)

    call messages_new_line()
    call messages_new_line()

    call messages_write("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", new_line = .true.)

    call messages_write("+++++++++++++++++++++ MINIMIZATION ITER #:")
    call messages_write(geom_iter, fmt = "I5")
    call messages_write(" ++++++++++++++++++++++", new_line = .true.)

    call messages_write("  Energy    = ")
    call messages_write(energy, units = units_out%energy, fmt = "f16.10,1x", print_units = .true., new_line = .true.)

    if(maxdf > M_ZERO) then
      call messages_write("  Max force = ")
      call messages_write(maxdf, units = units_out%force, fmt = "f16.10,1x", print_units = .true., new_line = .true.)
    end if

    call messages_write("  Max dr    = ")
    call messages_write(maxdx, units = units_out%length, fmt = "f16.10,1x", print_units = .true., new_line = .true.)

    call messages_write("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", new_line = .true.)
    call messages_write("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", new_line = .true.)
    call messages_new_line()
    call messages_new_line()
    call messages_info()

    POP_SUB(write_iter_info)
  end subroutine write_iter_info

  ! ---------------------------------------------------------

  subroutine to_coords(gopt, coords)
    type(geom_opt_t), intent(in)  :: gopt
    REAL_DOUBLE,      intent(out) :: coords(:)

    integer :: iatom, idir, icoord

    PUSH_SUB(to_coords)

    icoord = 1
    do iatom = 1, gopt%geo%natoms
      if(gopt%fixed_atom == iatom) cycle
      do idir = 1, gopt%dim
        coords(icoord) = gopt%geo%atom(iatom)%x(idir)
        if(gopt%fixed_atom /= 0) coords(icoord) = coords(icoord) - gopt%geo%atom(gopt%fixed_atom)%x(idir)
        icoord = icoord + 1
      end do
    end do

    POP_SUB(to_coords)
  end subroutine to_coords

  ! ---------------------------------------------------------

  subroutine to_grad(gopt, grad)
    type(geom_opt_t), intent(in)  :: gopt
    REAL_DOUBLE,      intent(out) :: grad(:)

    integer :: iatom, idir, icoord

    PUSH_SUB(to_grad)

    icoord = 1
    do iatom = 1, gopt%geo%natoms
      if(gopt%fixed_atom == iatom) cycle
      do idir = 1, gopt%dim
        grad(icoord) = -gopt%geo%atom(iatom)%f(idir)
        if(gopt%fixed_atom /= 0) grad(icoord) = grad(icoord) + gopt%geo%atom(gopt%fixed_atom)%f(idir)
        icoord = icoord + 1
      end do
    end do

    POP_SUB(to_grad)
  end subroutine to_grad

  ! ---------------------------------------------------------

  subroutine from_coords(gopt, coords)
    type(geom_opt_t), intent(inout) :: gopt
    REAL_DOUBLE,      intent(in)    :: coords(:)

    integer :: iatom, idir, icoord

    PUSH_SUB(from_coords)

    icoord = 1
    do iatom = 1, gopt%geo%natoms
      if(gopt%fixed_atom == iatom) cycle      
      do idir = 1, gopt%dim
        gopt%geo%atom(iatom)%x(idir) = coords(icoord)
        if(gopt%fixed_atom /= 0) then
          gopt%geo%atom(iatom)%x(idir) = gopt%geo%atom(iatom)%x(idir) + gopt%geo%atom(gopt%fixed_atom)%x(idir)
        end if
        icoord = icoord + 1
      end do
    end do

    POP_SUB(from_coords)
  end subroutine from_coords
    
  ! ---------------------------------------------------------
  !> Same as write_iter_info, but without the gradients.
  subroutine write_iter_info_ng(geom_iter, size, energy, maxdx, coords)
    integer,     intent(in) :: geom_iter
    integer,     intent(in) :: size !< must equal dim * natoms
    REAL_DOUBLE, intent(in) :: energy, maxdx
    REAL_DOUBLE, intent(in) :: coords(size)

    PUSH_SUB(write_iter_info_ng)
    call write_iter_info(geom_iter, size, energy, maxdx, real(-M_ONE, 8), coords)

    POP_SUB(write_iter_info_ng)
  end subroutine write_iter_info_ng

end module geom_opt_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
