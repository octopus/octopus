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

#include "global.h"

module geom_opt_oct_m
  use density_oct_m
  use energy_calc_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lcao_oct_m
  use loct_oct_m
  use mesh_oct_m
  use messages_oct_m
  use minimizer_oct_m
  use mpi_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pcm_oct_m
  use profiling_oct_m
  use read_coords_oct_m
  use restart_oct_m
  use scf_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use varinfo_oct_m

  implicit none

  private
  public :: geom_opt_run

  type geom_opt_t
    private
    integer  :: method
    FLOAT    :: step
    FLOAT    :: line_tol
    FLOAT    :: fire_mass
    integer  :: fire_integrator
    FLOAT    :: tolgrad
    FLOAT    :: toldr
    integer  :: max_iter
    integer  :: what2minimize

    !> shortcuts
    type(scf_t)                       :: scfv
    type(ions_t),             pointer :: ions
    type(hamiltonian_elec_t), pointer :: hm
    type(electrons_t),        pointer :: syst
    type(mesh_t),             pointer :: mesh
    type(states_elec_t),      pointer :: st
    integer                           :: dim
    integer                           :: size
    integer                           :: fixed_atom
    type(restart_t)                   :: restart_dump
    
  end type geom_opt_t

  type(geom_opt_t), save :: g_opt

  integer, parameter :: &
    MINWHAT_ENERGY = 1, &
    MINWHAT_FORCES = 2

contains

  ! ---------------------------------------------------------
  subroutine geom_opt_run(system, from_scratch)
    class(*),        intent(inout) :: system
    logical,         intent(inout) :: from_scratch

    PUSH_SUB(geom_opt_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = go not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call geom_opt_run_legacy(system, from_scratch)
    end select

    POP_SUB(geom_opt_run)
  end subroutine geom_opt_run

  ! ---------------------------------------------------------
  subroutine geom_opt_run_legacy(sys, fromscratch)
    type(electrons_t),   target, intent(inout) :: sys
    logical,                     intent(inout) :: fromscratch

    integer :: ierr
    REAL_DOUBLE, allocatable :: coords(:)
    REAL_DOUBLE :: energy

    FLOAT, allocatable :: mass(:)
    integer :: iatom, imass
    type(restart_t) :: restart_load

    PUSH_SUB(geom_opt_run_legacy)

    if (sys%space%periodic_dim == 1 .or. sys%space%periodic_dim == 2) then
      message(1) = "Geometry optimization not allowed for systems periodic in 1D or 2D, "
      message(2) = "as in those cases the total energy is not correct."
      call messages_fatal(2, namespace=sys%namespace)
    end if


    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%kpoints%use_symmetries) then
      call messages_experimental("KPoints symmetries with CalculationMode = go")
    end if

    call init_(fromscratch)

    ! load wavefunctions
    if(.not. fromscratch) then
      call restart_init(restart_load, sys%namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)
      if(ierr == 0) then
        call states_elec_load(restart_load, sys%namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, ierr)
      end if
      call restart_end(restart_load)
      if(ierr /= 0) then
        message(1) = "Unable to read wavefunctions: Starting from scratch."
        call messages_warning(1)
        fromscratch = .true.
      end if
    end if

    call scf_init(g_opt%scfv, sys%namespace, sys%gr, sys%ions, sys%st, sys%mc, sys%hm, sys%ks, sys%space)

    if(fromScratch) then
      call lcao_run(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm, lmm_r = g_opt%scfv%lmm_r)
    else
      ! setup Hamiltonian
      message(1) = 'Info: Setting up Hamiltonian.'
      call messages_info(1)
      call v_ks_h_setup(sys%namespace, sys%space, sys%gr, sys%ions, sys%st, sys%ks, sys%hm)
    end if

    !Initial point
    SAFE_ALLOCATE(coords(1:g_opt%size))
    call to_coords(g_opt, coords)

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%pack()

    !Minimize
    select case(g_opt%method)
    case(MINMETHOD_NMSIMPLEX)
      call minimize_multidim_nograd(g_opt%method, g_opt%size, coords, TOFLOAT(g_opt%step),&
        TOFLOAT(g_opt%toldr), g_opt%max_iter, &
        calc_point_ng, write_iter_info_ng, energy, ierr)
    case(MINMETHOD_FIRE)

      SAFE_ALLOCATE(mass(1:g_opt%size))
      imass = 1
      do iatom = 1, sys%ions%natoms
        if(g_opt%fixed_atom == iatom) cycle
        if (g_opt%ions%fixed(iatom)) cycle
        if (g_opt%fire_mass <= M_ZERO) then
          mass(imass:imass + 2) = sys%ions%mass(iatom)
        else
          mass(imass:imass + 2) = g_opt%fire_mass  ! Mass of H
        end if
        imass = imass + 3
      end do

      !TODO: add variable to use Euler integrator
      call minimize_fire(g_opt%size, coords, TOFLOAT(g_opt%step), TOFLOAT(g_opt%tolgrad), &
        g_opt%max_iter, calc_point, write_iter_info, energy, ierr, mass, integrator=g_opt%fire_integrator)
      SAFE_DEALLOCATE_A(mass)

    case default
      call minimize_multidim(g_opt%method, g_opt%size, coords, TOFLOAT(g_opt%step),&
        TOFLOAT(g_opt%line_tol), TOFLOAT(g_opt%tolgrad), TOFLOAT(g_opt%toldr), g_opt%max_iter, &
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

    if(sys%st%d%pack_states .and. hamiltonian_elec_apply_packed(sys%hm)) call sys%st%unpack()
  
    ! print out geometry
    call from_coords(g_opt, coords)
    message(1) = "Writing final coordinates to min.xyz"
    call messages_info(1)
    call g_opt%ions%write_xyz('./min')

    SAFE_DEALLOCATE_A(coords)
    call scf_end(g_opt%scfv)
    ! Because g_opt has the "save" atribute, we need to explicitly empty the criteria list here, or there will be a memory leak.
    call g_opt%scfv%criterion_list%empty()
    call end_()

    POP_SUB(geom_opt_run_legacy)
  contains

    ! ---------------------------------------------------------
    subroutine init_(fromscratch)
      logical,  intent(inout) :: fromscratch

      logical :: center, does_exist
      integer :: iter, iatom
      character(len=100) :: filename
      FLOAT :: default_toldr
      FLOAT :: default_step
      type(read_coords_info) :: xyz

      PUSH_SUB(geom_opt_run_legacy.init_)

      if (sys%space%periodic_dim > 0) then
        call messages_experimental('Geometry optimization for periodic systems')

        message(1) = "Optimization of cell parameters during geometry optimization"
        message(2) = "of periodic systems is currently not implemented."
        call messages_warning(2)
      end if

      call states_elec_allocate_wfns(sys%st, sys%gr%mesh)

      ! shortcuts
      g_opt%mesh   => sys%gr%mesh
      g_opt%ions    => sys%ions
      g_opt%st     => sys%st
      g_opt%hm     => sys%hm
      g_opt%syst   => sys
      g_opt%dim    =  sys%space%dim

      g_opt%size = g_opt%dim*g_opt%ions%natoms

      !%Variable GOCenter
      !%Type logical
      !%Default no
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% (Experimental) If set to yes, Octopus centers the geometry at
      !% every optimization step. It also reduces the degrees of
      !% freedom of the optimization by using the translational
      !% invariance.
      !%End
      call parse_variable(sys%namespace, 'GOCenter', .false.,  center)

      if(center) then
        g_opt%fixed_atom = 1
        g_opt%size = g_opt%size  - g_opt%dim
        call messages_experimental('GOCenter')
      end if

      !Check if atoms are allowed to move and redifine g_opt%size
      do iatom = 1, g_opt%ions%natoms
        if (g_opt%ions%fixed(iatom)) then
          g_opt%size = g_opt%size  - g_opt%dim
        end if
      end do
      
      !%Variable GOMethod
      !%Type integer
      !%Default fire
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
      !%Option fire 8
      !% The FIRE algorithm. See also <tt>GOFireMass</tt> and <tt>GOFireIntegrator</tt>.
      !% Ref: E. Bitzek, P. Koskinen, F. Gahler, M. Moseler, and P. Gumbsch, <i>Phys. Rev. Lett.</i> <b>97</b>, 170201 (2006).
      !%End
      call parse_variable(sys%namespace, 'GOMethod', MINMETHOD_FIRE, g_opt%method)
      if(.not.varinfo_valid_option('GOMethod', g_opt%method)) call messages_input_error(sys%namespace, 'GOMethod')
      call messages_print_var_option(stdout, "GOMethod", g_opt%method)

      !%Variable GOTolerance
      !%Type float
      !%Default 0.001 H/b (0.051 eV/A)
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Convergence criterion, for stopping the minimization. In
      !% units of force; minimization is stopped when all forces on
      !% ions are smaller than this criterion, or the
      !% <tt>GOMinimumMove</tt> is satisfied. If <tt>GOTolerance < 0</tt>,
      !% this criterion is ignored.
      !%End
      call parse_variable(sys%namespace, 'GOTolerance', CNST(0.001), g_opt%tolgrad, units_inp%force)
      
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
      end if
      call parse_variable(sys%namespace, 'GOMinimumMove', default_toldr, g_opt%toldr, units_inp%length)

      if(g_opt%method == MINMETHOD_NMSIMPLEX .and. g_opt%toldr <= M_ZERO) call messages_input_error(sys%namespace, 'GOMinimumMove')
      
      !%Variable GOStep
      !%Type float
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Initial step for the geometry optimizer. The default is 0.5.
      !% WARNING: in some weird units.
      !% For the FIRE minimizer, default value is 0.1 fs,
      !% and corresponds to the initial time-step for the MD.
      !%End
      if(g_opt%method /= MINMETHOD_FIRE ) then
        default_step = M_HALF
        call parse_variable(sys%namespace, 'GOStep', default_step, g_opt%step)
      else
        default_step = CNST(0.1)*unit_femtosecond%factor
        call parse_variable(sys%namespace, 'GOStep', default_step, g_opt%step, unit = units_inp%time)
      end if

      !%Variable GOLineTol
      !%Type float
      !%Default 0.1
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Tolerance for line-minimization. Applies only to GSL methods
      !% that use the forces.
      !% WARNING: in some weird units.
      !%End
      call parse_variable(sys%namespace, 'GOLineTol', CNST(0.1), g_opt%line_tol)

      !%Variable GOMaxIter
      !%Type integer
      !%Default 200
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Even if the convergence criterion is not satisfied, the minimization will stop
      !% after this number of iterations.
      !%End
      call parse_variable(sys%namespace, 'GOMaxIter', 200, g_opt%max_iter)
      if(g_opt%max_iter <= 0) then
        message(1) = "GOMaxIter has to be larger than 0"
        call messages_fatal(1)
      end if

      !%Variable GOFireMass
      !%Type float
      !%Default 1.0 amu
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% The Fire algorithm (<tt>GOMethod = fire</tt>) assumes that all degrees of freedom
      !% are comparable. All the velocities should be on the same
      !% scale,  which  for  heteronuclear  systems  can  be  roughly
      !% achieved by setting all the atom masses equal, to the value
      !% specified by this variable.
      !% By default the mass of a proton is selected (1 amu).
      !% However, a selection of <tt>GOFireMass = 0.01</tt> can, in manys systems, 
      !% speed up the geometry optimization procedure.
      !% If <tt>GOFireMass</tt> <= 0, the masses of each 
      !% species will be used.
      !%End
      call parse_variable(sys%namespace, 'GOFireMass', M_ONE*unit_amu%factor, g_opt%fire_mass, unit = unit_amu)

      !%Variable GOFireIntegrator
      !%Type integer
      !%Default verlet
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% The Fire algorithm (<tt>GOMethod = fire</tt>) uses a molecular dynamics
      !% integrator to compute new geometries and velocities.
      !% Currently, two integrator schemes can be selected 
      !%Option verlet 1
      !% The Velocity Verlet algorithm.
      !%Option euler 0
      !% The Euler method.
      !%End
      call parse_variable(sys%namespace, 'GOFireIntegrator', OPTION__GOFIREINTEGRATOR__VERLET, g_opt%fire_integrator)

      call messages_obsolete_variable(sys%namespace, 'GOWhat2Minimize', 'GOObjective')

      !%Variable GOObjective
      !%Type integer
      !%Default minimize_energy
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% This rather esoteric option allows one to choose which
      !% objective function to minimize during a geometry
      !% minimization. The use of this variable may lead to
      !% inconsistencies, so please make sure you know what you are
      !% doing.
      !%Option minimize_energy 1
      !% Use the total energy as objective function.
      !%Option minimize_forces 2
      !% Use <math>\sqrt{\sum_i \left| f_i \right|^2}</math> as objective function.
      !% Note that in this case one still uses the forces as the gradient of the objective function.
      !% This is, of course, inconsistent, and may lead to very strange behavior.
      !%End
      call parse_variable(sys%namespace, 'GOObjective', MINWHAT_ENERGY, g_opt%what2minimize)
      if(.not.varinfo_valid_option('GOObjective', g_opt%what2minimize)) call messages_input_error(sys%namespace, 'GOObjective')
      call messages_print_var_option(stdout, "GOObjective", g_opt%what2minimize)


      !%Variable XYZGOConstrains
      !%Type string
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% <tt>Octopus</tt> will try to read the coordinate-dependent constrains from the XYZ file 
      !% specified by the variable <tt>XYZGOConstrains</tt>.
      !% Note: It is important for the contrains to maintain the ordering 
      !% in which the atoms were defined in the coordinates specifications.
      !% Moreover, constrains impose fixed absolute coordinates, therefore
      !% constrains are not compatible with GOCenter = yes
      !%End

      !%Variable XSFGOConstrains
      !%Type string
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Like <tt>XYZGOConstrains</tt> but in XCrySDen format, as in <tt>XSFCoordinates</tt>.
      !%End

      !%Variable PDBGOConstrains
      !%Type string
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% Like <tt>XYZGOConstrains</tt> but in PDB format, as in <tt>PDBCoordinates</tt>.
      !%End

      !%Variable GOConstrains
      !%Type block
      !%Section Calculation Modes::Geometry Optimization
      !%Description
      !% If <tt>XYZGOConstrains</tt>, <tt>PDBConstrains</tt>, and <tt>XSFGOConstrains</tt>
      !% are not present, <tt>Octopus</tt> will try to fetch the geometry optimization
      !% contrains from this block. If this block is not present, <tt>Octopus</tt>
      !% will not set any constrains. The format of this block can be
      !% illustrated by this example:
      !%
      !% <tt>%GOConstrains
      !% <br>&nbsp;&nbsp;'C'  |      1 | 0 | 0
      !% <br>&nbsp;&nbsp;'O'  | &nbsp;1 | 0 | 0
      !% <br>%</tt>
      !%
      !% Coordinates with a constrain value of 0 will be optimized, while
      !% coordinates with a constrain different from zero will be kept fixed. So,
      !% in this example the x coordinates of both atoms will remain fixed and the
      !% distance between the two atoms along the x axis will be constant.
      !%
      !% Note: It is important for the constrains to maintain the ordering 
      !% in which the atoms were defined in the coordinates specifications.
      !% Moreover, constrains impose fixed absolute coordinates, therefore
      !% constrains are not compatible with GOCenter = yes
      !%End

      call read_coords_init(xyz)
      call read_coords_read('GOConstrains', xyz, g_opt%ions%space, sys%namespace)
      if(xyz%source /= READ_COORDS_ERR) then
        !Sanity check
        if(g_opt%ions%natoms /= xyz%n) then
          write(message(1), '(a,i4,a,i4)') 'I need exactly ', g_opt%ions%natoms, ' constrains, but I found ', xyz%n
          call messages_fatal(1)
        end if
        ! copy information and adjust units
        do iatom = 1, g_opt%ions%natoms
          where(xyz%atom(iatom)%x == M_ZERO)
            g_opt%ions%atom(iatom)%c = M_ZERO
          elsewhere
            g_opt%ions%atom(iatom)%c = M_ONE
          end where
        end do

        call read_coords_end(xyz)
       
        if(g_opt%fixed_atom /= 0) &
          call messages_not_implemented("GOCenter with constrains") 
      else
        do iatom = 1, g_opt%ions%natoms
          g_opt%ions%atom(iatom)%c = M_ZERO
        end do
      end if

      call io_rm("geom/optimization.log", sys%namespace)

      call io_rm("work-geom.xyz", sys%namespace)

      if(.not. fromScratch) then
        inquire(file = './last.xyz', exist = does_exist)
        if(.not. does_exist) fromScratch = .true.
      end if

      if(.not. fromScratch) call g_opt%ions%read_xyz('./last')
      
      ! clean out old geom/go.XXXX.xyz files. must be consistent with write_iter_info
      iter = 1
      do
        write(filename, '(a,i4.4,a)') "geom/go.", iter, ".xyz"
        inquire(file = trim(filename), exist = does_exist)
        if(does_exist) then
          call io_rm(trim(filename), sys%namespace)
          iter = iter + 1
        else
          exit
        end if
      ! TODO: clean forces directory
      end do

      call restart_init(g_opt%restart_dump, sys%namespace, RESTART_GS, RESTART_TYPE_DUMP, sys%mc, ierr, mesh=sys%gr%mesh)

      POP_SUB(geom_opt_run_legacy.init_)
    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      PUSH_SUB(geom_opt_run_legacy.end_)

      call states_elec_deallocate_wfns(sys%st)

      call restart_end(g_opt%restart_dump)

      nullify(g_opt%mesh)
      nullify(g_opt%ions)
      nullify(g_opt%st)
      nullify(g_opt%hm)
      nullify(g_opt%syst)

      POP_SUB(geom_opt_run_legacy.end_)
    end subroutine end_

  end subroutine geom_opt_run_legacy


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

    if(g_opt%fixed_atom /= 0) then
      call g_opt%ions%translate(g_opt%ions%center())
    end if

    ! When the system is periodic in some directions, the atoms might have moved to a an adjacent cell, so we need to move them back to the original cell
    call g_opt%ions%fold_atoms_into_cell()

    ! Some atoms might have moved outside the simulation box. We stop if this happens.
    do iatom = 1, g_opt%ions%natoms
      if (.not. g_opt%syst%gr%sb%contains_point(g_opt%syst%ions%pos(:, iatom))) then
        if (g_opt%syst%space%periodic_dim /= g_opt%syst%space%dim) then
          ! FIXME: This could fail for partial periodicity systems
          ! because contains_point is too strict with atoms close to
          ! the upper boundary to the cell.
          write(message(1), '(a,i5,a)') "Atom ", iatom, " has moved outside the box during the geometry optimization."
          call messages_fatal(1, namespace=g_opt%syst%namespace)
        end if
      end if
    end do

    call g_opt%ions%write_xyz('./work-geom', append = .true.)

    call scf_mix_clear(g_opt%scfv)

    call hamiltonian_elec_epot_generate(g_opt%hm, g_opt%syst%namespace, g_opt%syst%space, g_opt%syst%gr, g_opt%ions, g_opt%st)
    call density_calc(g_opt%st, g_opt%syst%gr, g_opt%st%rho)
    call v_ks_calc(g_opt%syst%ks, g_opt%syst%namespace, g_opt%syst%space, g_opt%hm, g_opt%st, g_opt%ions, calc_eigenval = .true.)
    call energy_calc_total(g_opt%syst%namespace, g_opt%syst%space, g_opt%hm, g_opt%syst%gr, g_opt%st)

    ! do SCF calculation
    call scf_run(g_opt%scfv, g_opt%syst%namespace, g_opt%syst%space, g_opt%syst%mc, g_opt%syst%gr, g_opt%ions, g_opt%st, &
      g_opt%syst%ks, g_opt%hm, g_opt%syst%outp, verbosity = VERB_COMPACT, restart_dump=g_opt%restart_dump)

    ! store results
    if(getgrad  ==  1) call to_grad(g_opt, df)

    if(g_opt%what2minimize == MINWHAT_FORCES) then
      objective = M_ZERO
      do iatom = 1, g_opt%ions%natoms
        if (g_opt%ions%fixed(iatom)) cycle
        objective = objective + sum(g_opt%ions%tot_force(:, iatom)**2)
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
  !! UPDATE: Because the newuoa routine have disappeared, probably this can be changed.
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

    character(len=256) :: c_geom_iter, title, c_forces_iter
    integer :: iunit

    PUSH_SUB(write_iter_info)
    
    write(c_geom_iter, '(a,i4.4)') "go.", geom_iter
    write(title, '(f16.10,2x,a)') units_from_atomic(units_out%energy, energy), trim(units_abbrev(units_out%energy))
    call io_mkdir('geom', g_opt%ions%namespace)
    call g_opt%ions%write_xyz('geom/'//trim(c_geom_iter), comment = trim(title))
    call g_opt%ions%write_xyz('./last')

    if (g_opt%syst%outp%what(OPTION__OUTPUT__FORCES)) then
    write(c_forces_iter, '(a,i4.4)') "forces.", geom_iter
      if (bitand(g_opt%syst%outp%how(OPTION__OUTPUT__FORCES), OPTION__OUTPUTFORMAT__BILD) /= 0) then
        call write_bild_forces_file('forces', trim(c_forces_iter), g_opt%ions, g_opt%syst%namespace)
      else
        call write_xsf_geometry_file('forces', trim(c_forces_iter), g_opt%ions, g_opt%syst%gr%mesh, &
          g_opt%syst%namespace, write_forces = .true.)
      end if
    end if

    call from_coords(g_opt, coords)

    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open(trim('geom/optimization.log'), g_opt%syst%namespace, &
        action = 'write', position = 'append')

      if(geom_iter == 1) then
        write(iunit, '(a10,3(5x,a20))') '#     iter','energy [' // trim(units_abbrev(units_out%energy)) // ']', & 
                             'max_force [' // trim(units_abbrev(units_out%force)) // ']',&
                             ' max_dr [' // trim(units_abbrev(units_out%length)) // ']'
      end if

      write(iunit, '(i10,3f25.15)')  geom_iter, units_from_atomic(units_out%energy, energy), & 
                                                units_from_atomic(units_out%force,maxdf), &
                                                units_from_atomic(units_out%length,maxdx)

      call io_close(iunit)
    end if

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
    do iatom = 1, gopt%ions%natoms
      if(gopt%fixed_atom == iatom) cycle
      if (gopt%ions%fixed(iatom)) cycle
      do idir = 1, gopt%dim
        coords(icoord) = gopt%ions%pos(idir, iatom)
        if(gopt%fixed_atom /= 0) coords(icoord) = coords(icoord) - gopt%ions%pos(idir, gopt%fixed_atom)
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
    do iatom = 1, gopt%ions%natoms
      if(gopt%fixed_atom == iatom) cycle
      if (gopt%ions%fixed(iatom)) cycle
      do idir = 1, gopt%dim
        if(abs(gopt%ions%atom(iatom)%c(idir)) <= M_EPSILON) then
          grad(icoord) = -gopt%ions%tot_force(idir, iatom)
        else
          grad(icoord) = M_ZERO
        end if
        if (gopt%fixed_atom /= 0) then
          grad(icoord) = grad(icoord) + gopt%ions%tot_force(idir, gopt%fixed_atom)
        end if
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
    do iatom = 1, gopt%ions%natoms
      if(gopt%fixed_atom == iatom) cycle      
      if (gopt%ions%fixed(iatom)) cycle
      do idir = 1, gopt%dim
        if(abs(gopt%ions%atom(iatom)%c(idir)) <= M_EPSILON) &
          gopt%ions%pos(idir, iatom) = coords(icoord)
        if(gopt%fixed_atom /= 0) then
          gopt%ions%pos(idir, iatom) = gopt%ions%pos(idir, iatom) + gopt%ions%pos(idir, gopt%fixed_atom)
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
    call write_iter_info(geom_iter, size, energy, maxdx, -M_ONE, coords)

    POP_SUB(write_iter_info_ng)
  end subroutine write_iter_info_ng

end module geom_opt_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
