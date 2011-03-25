!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module td_write_m
  use c_pointer_m
  use datasets_m
  use excited_states_m
  use forces_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use h_sys_output_m
  use hamiltonian_base_m
  use hamiltonian_m
  use io_m
  use ion_dynamics_m
  use lasers_m
  use loct_m
  use loct_math_m
  use magnetic_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use pert_m
  use profiling_m
  use restart_m
  use spectrum_m
  use states_m
  use states_calc_m
  use states_dim_m
  use types_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use write_iter_m

  implicit none

  private
  public ::         &
    td_write_t,     &
    td_write_init,  &
    td_write_end,   &
    td_write_iter,  &
    td_write_data

  type td_write_prop_t
    private
    type(c_ptr) :: handle
    logical :: write = .false.
  end type td_write_prop_t

  integer, parameter ::   &
    OUT_MULTIPOLES  =  1, &
    OUT_ANGULAR     =  2, &
    OUT_SPIN        =  3, &
    OUT_POPULATIONS =  4, &
    OUT_COORDS      =  5, &
    OUT_ACC         =  6, & 
    OUT_LASER       =  7, &
    OUT_ENERGY      =  8, &
    OUT_PROJ        =  9, &
    OUT_MAGNETS     = 10, &
    OUT_GAUGE_FIELD = 11, &
    OUT_TEMPERATURE = 12, &
    OUT_FTCHD       = 13, &
    OUT_VEL         = 14, &
    OUT_MAX         = 14
  
  type td_write_t
    private
    type(td_write_prop_t) :: out(OUT_MAX)

    integer        :: lmax     ! maximum multipole moment to output
    FLOAT          :: lmm_r    ! radius of the sphere used to compute the local magnetic moments
    type(states_t) :: gs_st    ! The states_type where the ground state is stored, in order to
                                        ! calculate the projections(s) onto it.
    integer        :: n_excited_states  ! number of excited states onto which the projections are calculated.
    type(excited_states_t), pointer :: excited_st(:) ! The excited states.
  end type td_write_t

contains

  ! ---------------------------------------------------------
  subroutine td_write_init(writ, gr, st, hm, geo, ions_move, with_gauge_field, kick, iter, max_iter, dt)
    type(td_write_t),    intent(out)   :: writ
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t),    intent(in)    :: geo
    logical,             intent(in)    :: ions_move
    logical,             intent(in)    :: with_gauge_field
    type(kick_t),        intent(in)    :: kick
    integer,             intent(in)    :: iter
    integer,             intent(in)    :: max_iter
    FLOAT,               intent(in)    :: dt


    FLOAT :: rmin
    integer :: ierr, first, ii, ist, jj, flags, iout, default
    type(block_t) :: blk
    character(len=100) :: filename

    PUSH_SUB(td_write_init)


    !%Variable TDOutput
    !%Type flag
    !%Default multipoles + geometry + temperature + energy
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent simulation.
    !%Option multipoles 1
    !% Outputs the multipole moments of the density to the file <tt>td.general/multipoles</tt>.
    !% This is required to, <i>e.g.</i>, calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>TDDipoleLmax</tt>.
    !%Option angular 2
    !% Outputs the angular momentum of the system that can be used to calculate circular
    !% dichroism (EXPERIMENTAL).
    !%Option spin 4
    !% Outputs the expectation value of the spin, that can be used to calculate magnetic
    !% cicular dichroism (EXPERIMENTAL).
    !%Option populations 8
    !% Outputs the projection of the time-dependent Kohn-Sham Slater determinant
    !% onto the ground-state (or approximations to the excited states) to the file 
    !% <tt>td.general/populations</tt>.
    !%Option geometry 16
    !% If set (and if the atoms are allowed to move), outputs the coordinates, velocities,
    !% and forces of the atoms to the the file <tt>td.general/coordinates</tt>.
    !%Option acceleration 32
    !% When set, outputs the acceleration, calculated from Ehrenfest theorem,
    !% in the file <tt>td.general/acceleration</tt>. This file can then be
    !% processed by the utility <tt>hs-from-acc</tt> in order to obtain the harmonic spectrum.
    !%Option laser 64
    !% If set, and if there are lasers defined in <tt>TDLasers</tt>,
    !% <tt>octopus</tt> outputs the laser field to the file <tt>td.general/laser</tt>.
    !%Option energy 128
    !% If <tt>set</tt>, <tt>octopus</tt> outputs the different components of the energy
    !% to the file <tt>td.general/el_energy</tt>.
    !%Option td_occup 256
    !% If set, outputs the projections of the time-dependent Kohn-Sham
    !% wavefunctions onto the static (zero-time) wavefunctions to the
    !% file <tt>td.general/projections.XXX</tt>.
    !%Option local_mag_moments 512
    !% If set, outputs the local magnetic moments, integrated in sphere centered around each atom.
    !% The radius of the sphere can be set with <tt>LocalMagneticMomentsSphereRadius</tt>.
    !%Option gauge_field 1024
    !% If set, outputs the vector gauge field corresponding to a spatially uniform (but time-dependent) 
    !% external electrical potential. This is only useful in a time-dependent periodic run.
    !%Option temperature 2048
    !% If set, the ionic temperature at each step is printed.
    !%Option ftchd       4096
    !% Write Fourier transform of the electron density to the file <tt>ftchds.X</tt>,
    !% where X depends on the kick (e.g. with sin-shaped perturbation X=qsin).
    !% This is needed for calculating the dynamic structure factor.
    !% In the case that the kick mode is qbessel, the written quantity is integral over
    !% density, multiplied by spherical Bessel function times real spherical harmonic.
    !%Option velocity    8192
    !% When set, outputs the velocity, calculated from Ehrenfest theorem,
    !% in the file <tt>td.general/velocity</tt>. This file can then be
    !% processed by the utility <tt>hs-from-vel</tt> in order to obtain the harmonic spectrum.
    !%End

    ! by default print multipoles, ftchd, coordinates and energy
    default = &
         2**(OUT_MULTIPOLES - 1) +  &
         2**(OUT_FTCHD - 1) +       &
         2**(OUT_COORDS - 1) +      &
         2**(OUT_TEMPERATURE - 1) + &
         2**(OUT_ENERGY - 1) +      &
         2**(OUT_GAUGE_FIELD - 1)

    call parse_integer(datasets_check('TDOutput'), default, flags)

    if(.not.varinfo_valid_option('TDOutput', flags, is_flag = .true.)) call input_error('TDOutput')

    do iout = 1, OUT_MAX
      writ%out(iout)%write = (iand(flags, 2**(iout - 1)) .ne. 0)
    end do

    !special cases
    writ%out(OUT_COORDS)%write = writ%out(OUT_COORDS)%write .and. ions_move
    writ%out(OUT_TEMPERATURE)%write = writ%out(OUT_TEMPERATURE)%write .and. ions_move
    writ%out(OUT_GAUGE_FIELD)%write = writ%out(OUT_GAUGE_FIELD)%write .and. with_gauge_field
    writ%out(OUT_LASER)%write = writ%out(OUT_LASER)%write .and. (hm%ep%no_lasers > 0)

    !%Variable TDDipoleLmax
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::TD Output
    !%Description
    !% Maximum multipole of the density output to the file <tt>td.general/multipoles</tt>
    !% during a time-dependent simulation. Must be 0 &lt; <tt>TDDipoleLmax &lt; 5</tt>.
    !%End
    call parse_integer(datasets_check('TDDipoleLmax'), 1, writ%lmax)
    if (writ%lmax < 0 .or. writ%lmax > 4) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid TDDipoleLmax."
      message(2) = '(0 <= TDDipoleLmax <= 4 )'
      call messages_fatal(2)
    end if

    ! Compatibility test
    if( (writ%out(OUT_ACC)%write) .and. ions_move ) then
      message(1) = 'Error: If harmonic spectrum is to be calculated,'
      message(2) = 'atoms should not be allowed to move.'
      call messages_fatal(2)
    end if

    rmin = geometry_min_distance(geo)

    ! This variable is documented in scf/scf.F90
    call parse_float(datasets_check('LocalMagneticMomentsSphereRadius'), rmin*M_HALF, writ%lmm_r, units_inp%length)

    if( (writ%out(OUT_PROJ)%write)  .or.  (writ%out(OUT_POPULATIONS)%write) ) then
      call states_copy(writ%gs_st, st)

      ! clean up all the stuff we have to reallocate
      SAFE_DEALLOCATE_P(writ%gs_st%zpsi)
      SAFE_DEALLOCATE_P(writ%gs_st%occ)
      SAFE_DEALLOCATE_P(writ%gs_st%eigenval)
      SAFE_DEALLOCATE_P(writ%gs_st%node)
      if(writ%gs_st%d%ispin == SPINORS) then
        SAFE_DEALLOCATE_P(writ%gs_st%spin)
      end if

      call states_look (trim(restart_dir)//'gs', gr%mesh%mpi_grp, ii, jj, writ%gs_st%nst, ierr)

      if(writ%out(OUT_POPULATIONS)%write) then ! do only this when not calculating populations
        ! We will store the ground-state Kohn-Sham system for all processors.
        !%Variable TDProjStateStart
        !%Type integer
        !%Default 1
        !%Section Time-Dependent::TD Output
        !%Description
        !% Only output projections to states above <tt>TDProjStateStart</tt>. Usually one is only interested
        !% in particle-hole projections around the HOMO, so there is no need to calculate (and store)
        !% the projections of all TD states onto all static states. This sets a lower limit. The upper limit
        !% is set by the number of states in the propagation and the number of unoccupied states
        !% available.
        !%End
        call parse_integer(datasets_check('TDProjStateStart'), 1, writ%gs_st%st_start)
      else
        writ%gs_st%st_start = 1
      end if
      writ%gs_st%st_end = writ%gs_st%nst

      ! allocate memory
      SAFE_ALLOCATE(writ%gs_st%occ(1:writ%gs_st%nst, 1:writ%gs_st%d%nik))
      SAFE_ALLOCATE(writ%gs_st%eigenval(1:writ%gs_st%nst, 1:writ%gs_st%d%nik))
      SAFE_ALLOCATE(writ%gs_st%node(1:writ%gs_st%nst))
      writ%gs_st%eigenval = huge(writ%gs_st%eigenval)
      writ%gs_st%occ      = M_ZERO
      if(writ%gs_st%d%ispin == SPINORS) then
        SAFE_ALLOCATE(writ%gs_st%spin(1:3, 1:writ%gs_st%nst, 1:writ%gs_st%d%nik))
      end if
      call states_allocate_wfns(writ%gs_st, gr%mesh, TYPE_CMPLX)
      writ%gs_st%node(:)  = 0
      call restart_read(trim(restart_dir)//'gs', writ%gs_st, gr, geo, ierr)
      if(ierr.ne.0 .and.ierr.ne.(writ%gs_st%st_end-writ%gs_st%st_start+1)*writ%gs_st%d%nik*writ%gs_st%d%dim) then
        message(1) = "Could not load "//trim(restart_dir)//"gs"
        call messages_fatal(1)
      end if
    end if

    ! Build the excited states...
    if(writ%out(OUT_POPULATIONS)%write) then
      !%Variable TDExcitedStatesToProject
      !%Type block
      !%Section Time-Dependent::TD Output
      !%Description
      !% <b>[WARNING: This is a *very* experimental feature]</b> The population of the excited states
      !% (as defined by <Phi_I|Phi(t)> where |Phi(t)> is the many-body time-dependent state at
      !% time <i>t</i>, and |Phi_I> is the excited state of interest) can be approximated -- it is not clear 
      !% how well -- by substituting for those real many-body states the time-dependent Kohn-Sham
      !% determinant and some modification of the Kohn-Sham ground-state determinant (<i>e.g.</i>,
      !% a simple HOMO-LUMO substitution, or the Casida ansatz for excited states in linear-response
      !% theory. If you set <tt>TDOutput</tt> to contain <tt>populations</tt>, you may ask for these approximated
      !% populations for a number of excited states, which will be described in the files specified
      !% in this block: each line should be the name of a file that contains one excited state.
      !%
      !% FIXME: description of the format of the files.
      !%End
      if(parse_block('TDExcitedStatesToProject', blk) == 0) then
        writ%n_excited_states = parse_block_n(blk)
        SAFE_ALLOCATE(writ%excited_st(1:writ%n_excited_states))
        do ist = 1, writ%n_excited_states
          call parse_block_string(blk, ist-1, 0, filename)
          call excited_states_init(writ%excited_st(ist), writ%gs_st, trim(filename)) 
        end do
      else
        writ%n_excited_states = 0
        nullify(writ%excited_st)
      end if
    end if

    if (iter == 0) then
      first = 0
    else
      first = iter + 1
    end if

    call io_mkdir('td.general')

    if(mpi_grp_is_root(mpi_world)) then
      if(writ%out(OUT_MULTIPOLES)%write) &
        call write_iter_init(writ%out(OUT_MULTIPOLES)%handle, &
        first, units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/multipoles")))

      if(writ%out(OUT_FTCHD)%write) then
        select case(kick%qkick_mode)
          case (QKICKMODE_SIN)
            write(filename, '(a)') 'td.general/ftchd.sin'
          case (QKICKMODE_COS)
            write(filename, '(a)') 'td.general/ftchd.cos'
          case (QKICKMODE_BESSEL)
            write(filename, '(a, SP, I0.3, a, I0.3)') 'td.general/ftchd.l', kick%qbessel_l, '_m', kick%qbessel_m
          case default
            write(filename, '(a)') 'td.general/ftchd'
        end select
        call write_iter_init(writ%out(OUT_FTCHD)%handle, &
          first, units_from_atomic(units_out%time, dt), trim(io_workpath(filename)))
      end if  

      if(writ%out(OUT_ANGULAR)%write) &
        call write_iter_init(writ%out(OUT_ANGULAR)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/angular")))

      if(writ%out(OUT_SPIN)%write) &
        call write_iter_init(writ%out(OUT_SPIN)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/spin")))

      if(writ%out(OUT_MAGNETS)%write) &
        call write_iter_init(writ%out(OUT_MAGNETS)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/magnetic_moments")))

      if(writ%out(OUT_COORDS)%write) &
        call write_iter_init(writ%out(OUT_COORDS)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/coordinates")))

      if(writ%out(OUT_TEMPERATURE)%write) &
        call write_iter_init(writ%out(OUT_TEMPERATURE)%handle, first, M_ONE, trim(io_workpath("td.general/temperature")))

      if(writ%out(OUT_POPULATIONS)%write) &
        call write_iter_init(writ%out(OUT_POPULATIONS)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/populations")))

      if(writ%out(OUT_ACC)%write) &
        call write_iter_init(writ%out(OUT_ACC)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/acceleration")))
          
      if(writ%out(OUT_VEL)%write) &
        call write_iter_init(writ%out(OUT_VEL)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/velocity")))

      if(writ%out(OUT_LASER)%write) then
        call write_iter_init(writ%out(OUT_LASER)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/laser")))
        do ii = 0, max_iter
          call td_write_laser(writ%out(OUT_LASER)%handle, gr, hm, dt, ii)
          if(mod(ii, 100).eq.0) call write_iter_flush(writ%out(OUT_LASER)%handle)
        end do
        call write_iter_end(writ%out(OUT_LASER)%handle)
      end if


      if(writ%out(OUT_ENERGY)%write) &
        call write_iter_init(writ%out(OUT_ENERGY)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/energy")))

      if(writ%out(OUT_PROJ)%write) &
        call write_iter_init(writ%out(OUT_PROJ)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/projections")))

      if(writ%out(OUT_GAUGE_FIELD)%write) &
        call write_iter_init(writ%out(OUT_GAUGE_FIELD)%handle, &
        first, units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/gauge_field")))
    end if

    POP_SUB(td_write_init)
  end subroutine td_write_init


  ! ---------------------------------------------------------
  subroutine td_write_end(writ)
    type(td_write_t), intent(inout) :: writ
    integer :: ist, iout
    PUSH_SUB(td_write_end)

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAX
        if(iout.eq.OUT_LASER) cycle
        if(writ%out(iout)%write)  call write_iter_end(writ%out(iout)%handle)
      end do
    end if

    if( writ%out(OUT_POPULATIONS)%write ) then
      do ist = 1, writ%n_excited_states
        call excited_states_kill(writ%excited_st(ist))
      end do
    end if

    if(writ%out(OUT_POPULATIONS)%write .or. writ%out(OUT_PROJ)%write) call states_end(writ%gs_st)

    POP_SUB(td_write_end)
  end subroutine td_write_end


  ! ---------------------------------------------------------
  subroutine td_write_iter(writ, gr, st, hm, geo, kick, dt, iter)
    type(td_write_t),    intent(in)    :: writ
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t),    intent(inout) :: geo
    type(kick_t),        intent(in)    :: kick
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    type(profile_t), save :: prof

    PUSH_SUB(td_write_iter)
    call profiling_in(prof, "TD_WRITE_ITER")

    if(writ%out(OUT_MULTIPOLES)%write) &
      call td_write_multipole(writ%out(OUT_MULTIPOLES)%handle, gr, geo, st, writ%lmax, kick, iter)
    
    if(writ%out(OUT_FTCHD)%write) &
      call td_write_ftchd(writ%out(OUT_FTCHD)%handle, gr, geo, st, kick, iter)

    if(writ%out(OUT_ANGULAR)%write) &
      call td_write_angular(writ%out(OUT_ANGULAR)%handle, gr, geo, hm, st, kick, iter)

    if(writ%out(OUT_SPIN)%write) &
      call td_write_spin(writ%out(OUT_SPIN)%handle, gr, st, iter)

    if(writ%out(OUT_MAGNETS)%write) &
      call td_write_local_magnetic_moments(writ%out(OUT_MAGNETS)%handle, gr, st, geo, writ%lmm_r, iter)

    if(writ%out(OUT_PROJ)%write) &
      call td_write_proj(writ%out(OUT_PROJ)%handle, gr, geo, st, writ%gs_st, kick, iter)

    if(writ%out(OUT_COORDS)%write) &
      call td_write_coordinates(writ%out(OUT_COORDS)%handle, gr, geo, iter)

    if(writ%out(OUT_TEMPERATURE)%write) &
      call td_write_temperature(writ%out(OUT_TEMPERATURE)%handle, geo, iter)

    if(writ%out(OUT_POPULATIONS)%write) &
      call td_write_populations(writ%out(OUT_POPULATIONS)%handle, gr%mesh, st, &
        writ%gs_st, writ%n_excited_states, writ%excited_st, dt, iter)

    if(writ%out(OUT_ACC)%write) &
      call td_write_acc(writ%out(OUT_ACC)%handle, gr, geo, st, hm, dt, iter)
      
    if(writ%out(OUT_VEL)%write) &
      call td_write_vel(writ%out(OUT_VEL)%handle, gr, geo, st, hm, dt, iter)

    ! td_write_laser no longer called here, because the whole laser is printed
    ! out at the beginning.

    if(writ%out(OUT_ENERGY)%write) &
      call td_write_energy(writ%out(OUT_ENERGY)%handle, hm, iter, geo%kinetic_energy)

    if(writ%out(OUT_GAUGE_FIELD)%write) &
      call td_write_gauge_field(writ%out(OUT_GAUGE_FIELD)%handle, hm, gr, iter)

    call profiling_out(prof)
    POP_SUB(td_write_iter)
  end subroutine td_write_iter


  ! ---------------------------------------------------------
  subroutine td_write_data(writ, gr, st, hm, outp, geo, iter)
    type(td_write_t),     intent(in)    :: writ
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    type(h_sys_output_t), intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: iter

    character(len=256) :: filename
    integer :: iout
    type(profile_t), save :: prof

    PUSH_SUB(td_write_data)
    call profiling_in(prof, "TD_WRITE_DATA")

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAX
        if(iout.eq.OUT_LASER) cycle
        if(writ%out(iout)%write)  call write_iter_flush(writ%out(iout)%handle)
      end do
    end if

    ! now write down the rest
    write(filename, '(a,i7.7)') "td.", iter  ! name of directory

    call h_sys_output_all(outp, gr, geo, st, hm, filename)

    call profiling_out(prof)
    POP_SUB(td_write_data)
  end subroutine td_write_data


  ! ---------------------------------------------------------
  subroutine td_write_spin(out_spin, gr, st, iter)
    type(c_ptr), intent(in)    :: out_spin
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(in)    :: st
    integer,           intent(in)    :: iter

    character(len=130) :: aux
    FLOAT :: spin(3)

    PUSH_SUB(td_write_spin)

    ! The expectation value of the spin operator is half the total magnetic moment
    ! This has to be calculated by all nodes
    call magnetic_moment(gr%mesh, st, st%rho, spin)
    spin = M_HALF*spin

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_spin)

        !second line -> columns name
        call write_iter_header_start(out_spin)
        if (st%d%ispin == SPINORS) then
          write(aux, '(a2,18x)') 'Sx'
          call write_iter_header(out_spin, aux)
          write(aux, '(a2,18x)') 'Sy'
          call write_iter_header(out_spin, aux)
        end if
        write(aux, '(a2,18x)') 'Sz'
        call write_iter_header(out_spin, aux)
        call write_iter_nl(out_spin)

        call td_write_print_header_end(out_spin)
      end if

      call write_iter_start(out_spin)
      select case (st%d%ispin)
      case (SPIN_POLARIZED)
         call write_iter_double(out_spin, spin(3), 1)
      case (SPINORS)
        call write_iter_double(out_spin, spin(1:3), 3)
      end select
      call write_iter_nl(out_spin)

    end if

    POP_SUB(td_write_spin)
  end subroutine td_write_spin


  ! ---------------------------------------------------------
  subroutine td_write_local_magnetic_moments(out_magnets, gr, st, geo, lmm_r, iter)
    type(c_ptr),              intent(in)    :: out_magnets
    type(grid_t),             intent(inout) :: gr
    type(states_t),           intent(in)    :: st
    type(geometry_t),         intent(in)    :: geo
    FLOAT,                    intent(in)    :: lmm_r
    integer,                  intent(in)    :: iter

    integer :: ia
    character(len=50) :: aux
    FLOAT, allocatable :: lmm(:,:)

    PUSH_SUB(td_write_local_magnetic_moments)

    !get the atoms` magnetization. This has to be calculated by all nodes
    SAFE_ALLOCATE(lmm(1:3, 1:geo%natoms))
    call magnetic_local_moments(gr%mesh, st, geo, st%rho, lmm_r, lmm)

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_magnets)

        !second line -> columns name
        call write_iter_header_start(out_magnets)
        do ia = 1, geo%natoms
          if (st%d%ispin == SPINORS) then
            write(aux, '(a2,i2.2,16x)') 'mx', ia
            call write_iter_header(out_magnets, aux)
            write(aux, '(a2,i2.2,16x)') 'my', ia
            call write_iter_header(out_magnets, aux)
          end if
          write(aux, '(a2,i2.2,16x)') 'mz', ia
          call write_iter_header(out_magnets, aux)
        end do
        call write_iter_nl(out_magnets)

        call td_write_print_header_end(out_magnets)
      end if

      call write_iter_start(out_magnets)
      do ia = 1, geo%natoms
        select case (st%d%ispin)
        case (SPIN_POLARIZED)
          call write_iter_double(out_magnets, lmm(3, ia), 1)
        case (SPINORS)
          call write_iter_double(out_magnets, lmm(1:3, ia), 3)
        end select
      end do
      call write_iter_nl(out_magnets)
      SAFE_DEALLOCATE_A(lmm)
    end if

    POP_SUB(td_write_local_magnetic_moments)
  end subroutine td_write_local_magnetic_moments


  ! ---------------------------------------------------------
  subroutine td_write_angular(out_angular, gr, geo, hm, st, kick, iter)
    type(c_ptr),            intent(in)    :: out_angular
    type(grid_t),           intent(inout) :: gr
    type(geometry_t),       intent(inout) :: geo
    type(hamiltonian_t),    intent(inout) :: hm
    type(states_t),         intent(inout) :: st
    type(kick_t),           intent(in)    :: kick
    integer,                intent(in)    :: iter

    integer :: ik, ist, idir
    character(len=130) :: aux
    FLOAT :: angular(MAX_DIM)
    type(pert_t)        :: angular_momentum

    PUSH_SUB(td_write_angular)

    call pert_init(angular_momentum, PERTURBATION_MAGNETIC, gr, geo)

    do idir = 1, 3
       call pert_setup_dir(angular_momentum, idir)
       !we have to multiply by 2, because the perturbation returns L/2
       angular(idir) = M_TWO*zpert_expectation_value(angular_momentum, gr, geo, hm, st, st%zpsi, st%zpsi)
    end do

    call pert_end(angular_momentum)

    if(mpi_grp_is_root(mpi_world)) then ! Only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_angular)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        call kick_write(kick, out = out_angular)

        !second line -> columns name
        call write_iter_header_start(out_angular)
        write(aux, '(a4,18x)') '<Lx>'
        call write_iter_header(out_angular, aux)
        write(aux, '(a4,18x)') '<Ly>'
        call write_iter_header(out_angular, aux)
        write(aux, '(a4,18x)') '<Lz>'
        call write_iter_header(out_angular, aux)
        call write_iter_nl(out_angular)

        !third line -> should hold the units.
        call write_iter_string(out_angular, '#[Iter n.]')
        call write_iter_header(out_angular, '[' // trim(units_abbrev(units_out%time)) // ']')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_header(out_angular, '[hbar]')
        call write_iter_nl(out_angular)

        call td_write_print_header_end(out_angular)
      end if

      call write_iter_start(out_angular)
      call write_iter_double(out_angular, angular(1:3), 3)
      call write_iter_nl(out_angular)

    end if

    POP_SUB(td_write_angular)
  end subroutine td_write_angular


  ! ---------------------------------------------------------
  subroutine td_write_multipole(out_multip, gr, geo, st, lmax, kick, iter)
    type(c_ptr),        intent(in) :: out_multip
    type(grid_t),       intent(in) :: gr
    type(geometry_t),   intent(in) :: geo
    type(states_t),     intent(in) :: st
    integer,            intent(in) :: lmax
    type(kick_t),       intent(in) :: kick
    integer,            intent(in) :: iter

    integer :: is, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT, allocatable :: nuclear_dipole(:), multipole(:,:)

    PUSH_SUB(td_write_multipole)

    if(mpi_grp_is_root(mpi_world).and.iter == 0) then
      call td_write_print_header_init(out_multip)

      write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,i2)')      '# lmax         ', lmax
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      call kick_write(kick, out = out_multip)

      call write_iter_header_start(out_multip)

      do is = 1, st%d%nspin
        write(aux,'(a18,i1,a1)') 'Electronic charge(', is,')'; call write_iter_header(out_multip, aux)
        if(lmax>0) then
          write(aux, '(a3,a1,i1,a1)') '<x>', '(', is,')'; call write_iter_header(out_multip, aux)
          write(aux, '(a3,a1,i1,a1)') '<y>', '(', is,')'; call write_iter_header(out_multip, aux)
          write(aux, '(a3,a1,i1,a1)') '<z>', '(', is,')'; call write_iter_header(out_multip, aux)
        end if
        do ll = 2, lmax
          do mm = -ll, ll
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', ll, ', m=', mm, ' (', is,')'
            call write_iter_header(out_multip, aux)
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      ! units
      call write_iter_string(out_multip, '#[Iter n.]')
      call write_iter_header(out_multip, '[' // trim(units_abbrev(units_out%time)) // ']')

      do is = 1, st%d%nspin
        do ll = 0, lmax
          do mm = -ll, ll
            select case(ll)
            case(0)
              call write_iter_header(out_multip, 'Electrons')
            case(1)
              call write_iter_header(out_multip, '[' // trim(units_abbrev(units_out%length)) // ']')
            case default
              write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
              call write_iter_header(out_multip, '[' // trim(aux) // ']')
            end select
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      call td_write_print_header_end(out_multip)
    end if

    SAFE_ALLOCATE(nuclear_dipole(1:MAX_DIM))
    SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin))
    nuclear_dipole(:) = M_ZERO
    multipole   (:,:) = M_ZERO

    do is = 1, st%d%nspin
      call dmf_multipoles(gr%mesh, st%rho(:,is), lmax, multipole(:,is))
    end do
    call geometry_dipole(geo, nuclear_dipole)
    do is = 1, st%d%nspin
      multipole(2:gr%mesh%sb%dim+1, is) = -nuclear_dipole(1:gr%mesh%sb%dim) - multipole(2:gr%mesh%sb%dim+1, is)
    end do

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_multip)
      do is = 1, st%d%nspin
        add_lm = 1
        do ll = 0, lmax
          do mm = -ll, ll
            call write_iter_double(out_multip, units_from_atomic(units_out%length**ll, multipole(add_lm, is)), 1)
            add_lm = add_lm + 1
          end do
        end do
      end do
      call write_iter_nl(out_multip)
    end if

    SAFE_DEALLOCATE_A(nuclear_dipole)
    SAFE_DEALLOCATE_A(multipole)
    POP_SUB(td_write_multipole)
  end subroutine td_write_multipole

  ! ---------------------------------------------------------
  subroutine td_write_ftchd(out_ftchd, gr, geo, st, kick, iter)
    type(c_ptr),        intent(in) :: out_ftchd
    type(grid_t),       intent(in) :: gr
    type(geometry_t),   intent(in) :: geo
    type(states_t),     intent(in) :: st
    type(kick_t),       intent(in) :: kick
    integer,            intent(in) :: iter

    integer :: is, ip, idir
    character(len=120) :: aux, aux2
    FLOAT   :: ftchd_bessel
    CMPLX   :: ftchd
    FLOAT   :: ylm, gylm(1:MAX_DIM)
    FLOAT, allocatable :: integrand_bessel(:)
    CMPLX, allocatable :: integrand(:)

    PUSH_SUB(td_write_ftchd)

    if(mpi_grp_is_root(mpi_world).and.iter == 0) then
      call td_write_print_header_init(out_ftchd)

      write(aux,'(a15, i2)') '# qkickmode    ', kick%qkick_mode
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      if(kick%qkick_mode.eq.QKICKMODE_BESSEL) then
        write(aux,'(a15, i0.3, 1x, i0.3)') '# ll, mm       ', kick%qbessel_l, kick%qbessel_m
        call write_iter_string(out_ftchd, aux)
        call write_iter_nl(out_ftchd)
      end if

      if(kick%qkick_mode.eq.QKICKMODE_BESSEL) then
        write(aux, '(a15, f9.6)') '# qlength      ', kick%qlength
      else ! sin or cos
        write(aux, '(a15)')       '# qvector      '
        do idir = 1, gr%mesh%sb%dim
          write(aux2, '(f9.5)') kick%qvector(idir)
          aux = trim(aux) // trim(aux2)
        enddo
      end if
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      call write_iter_header_start(out_ftchd)
      if(kick%qkick_mode.eq.QKICKMODE_BESSEL) then
        write(aux,'(a17)') 'int(j_l*Y_lm*rho)'
      else
        write(aux,'(a12)') 'Real, Imag'
      end if
      call write_iter_header(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      ! units
      call write_iter_string(out_ftchd, '#[Iter n.]')
      call write_iter_header(out_ftchd, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_nl(out_ftchd)
      call td_write_print_header_end(out_ftchd)

    end if

    ftchd = M_ZERO

    ! If kick mode is exp, sin, or cos, apply the normal Fourier transform
    if(kick%qkick_mode.ne.QKICKMODE_BESSEL) then
      SAFE_ALLOCATE(integrand(1:gr%mesh%np))
      integrand = M_ZERO
      do is = 1, st%d%nspin
        forall(ip = 1:gr%mesh%np)
          integrand(ip) = integrand(ip) + st%rho(ip, is) * exp(-M_zI*sum(gr%mesh%x(ip,:)*kick%qvector(:)))
        end forall
      end do
      ftchd = zmf_integrate(gr%mesh, integrand)
      SAFE_DEALLOCATE_A(integrand)
    else
      ftchd_bessel = M_ZERO
      SAFE_ALLOCATE(integrand_bessel(1:gr%mesh%np))
      integrand_bessel = M_ZERO
      do is = 1, st%d%nspin
        do ip = 1, gr%mesh%np
          call grylmr(gr%mesh%x(ip, 1), gr%mesh%x(ip, 2), gr%mesh%x(ip, 3), kick%qbessel_l, kick%qbessel_m, ylm, gylm)
          integrand_bessel(ip) = integrand_bessel(ip) + st%rho(ip, is) * &
                                 loct_sph_bessel(kick%qbessel_l, kick%qlength*sqrt(sum(gr%mesh%x(ip, :)**2)))*ylm
        end do
      end do
      ftchd_bessel = dmf_integrate(gr%mesh, integrand_bessel)
      SAFE_DEALLOCATE_A(integrand_bessel)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_ftchd)
      if(kick%qkick_mode.eq.QKICKMODE_BESSEL) then
        call write_iter_double(out_ftchd, ftchd_bessel, 1)
      else ! exp, sin, cos
        call write_iter_double(out_ftchd, real(ftchd), 1)
        call write_iter_double(out_ftchd, aimag(ftchd), 1)
      end if
      call write_iter_nl(out_ftchd)
    end if

    POP_SUB(td_write_ftchd)
  end subroutine td_write_ftchd

  ! ---------------------------------------------------------
  subroutine td_write_coordinates(out_coords, gr, geo, iter)
    type(c_ptr),       intent(in) :: out_coords
    type(grid_t),      intent(in) :: gr
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iter

    integer :: iatom, idir
    character(len=50) :: aux
    FLOAT :: tmp(1:MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_coordinates)

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do iatom = 1, geo%natoms
        do idir = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', iatom, ',', idir, ')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do iatom = 1, geo%natoms
        do idir = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', iatom, ',', idir,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do iatom = 1, geo%natoms
        do idir = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', iatom, ',', idir,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(out_coords, &
        'Positions in '   // trim(units_abbrev(units_out%length))   //   &
        ', Velocities in '// trim(units_abbrev(units_out%velocity)) //   &
        ', Forces in '    // trim(units_abbrev(units_out%force)))
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    do iatom = 1, geo%natoms
      tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%length, geo%atom(iatom)%x(1:gr%mesh%sb%dim))
      call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
    end do
    do iatom = 1, geo%natoms
      tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%velocity, geo%atom(iatom)%v(1:gr%mesh%sb%dim))
      call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
    end do
    do iatom = 1, geo%natoms
      tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%force, geo%atom(iatom)%f(1:gr%mesh%sb%dim))
      call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
    end do
    call write_iter_nl(out_coords)

    POP_SUB(td_write_coordinates)
  end subroutine td_write_coordinates


  ! ---------------------------------------------------------
  subroutine td_write_temperature(out_temperature, geo, iter)
    type(c_ptr),       intent(in) :: out_temperature
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iter

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_temperature)

    if(iter == 0) then
      call td_write_print_header_init(out_temperature)

      ! first line: column names
      call write_iter_header_start(out_temperature)
      call write_iter_header(out_temperature, 'Temperature')
      call write_iter_nl(out_temperature)

      ! second line: units
      call write_iter_string(out_temperature, '#[Iter n.]')
      call write_iter_header(out_temperature, '[' // trim(units_abbrev(units_out%time)) // ']')
      call write_iter_string(out_temperature, '        [K]')
      call write_iter_nl(out_temperature)

      call td_write_print_header_end(out_temperature)
    end if

    call write_iter_start(out_temperature)

    call write_iter_double(out_temperature, ion_dynamics_temperature(geo), 1)

    call write_iter_nl(out_temperature)

    POP_SUB(td_write_temperature)
  end subroutine td_write_temperature


  ! ---------------------------------------------------------
  subroutine td_write_populations(out_populations, mesh, st, gs_st, n_excited_states, excited_st, dt, iter)
    type(c_ptr),            intent(in) :: out_populations
    type(mesh_t),           intent(in) :: mesh
    type(states_t),         intent(in) :: st
    type(states_t),         intent(in) :: gs_st
    integer,                intent(in) :: n_excited_states
    type(excited_states_t), intent(in) :: excited_st(:)
    FLOAT,                  intent(in) :: dt
    integer,                intent(in) :: iter
 
    integer :: ist
    character(len=6) :: excited_name
    CMPLX :: gsp
    CMPLX, allocatable :: excited_state_p(:)
    CMPLX, allocatable :: dotprodmatrix(:, :, :)


    PUSH_SUB(td_write_populations)

    SAFE_ALLOCATE(dotprodmatrix(1:gs_st%nst, 1:st%nst, 1:st%d%nik))
    call zstates_matrix(mesh, gs_st, st, dotprodmatrix)

    ! all processors calculate the projection
    gsp = zstates_mpdotp(mesh, gs_st, st, dotprodmatrix)

    if(n_excited_states > 0) then
      SAFE_ALLOCATE(excited_state_p(1:n_excited_states))
      do ist = 1, n_excited_states
        excited_state_p(ist) = zstates_mpdotp(mesh, excited_st(ist), st, dotprodmatrix)
      end do
    end if

    if(mpi_grp_is_root(mpi_world)) then
      if(iter == 0) then
        call td_write_print_header_init(out_populations)

        ! first line -> column names
        call write_iter_header_start(out_populations)
        call write_iter_header(out_populations, 'Re<Phi_gs|Phi(t)>')
        call write_iter_header(out_populations, 'Im<Phi_gs|Phi(t)>')
        do ist = 1, n_excited_states
          write(excited_name,'(a2,i3,a1)') 'P(', ist,')'
          call write_iter_header(out_populations, 'Re<'//excited_name//'|Phi(t)>')
          call write_iter_header(out_populations, 'Im<'//excited_name//'|Phi(t)>')
        end do
        call write_iter_nl(out_populations)

        ! second line -> units
        call write_iter_string(out_populations, '#[Iter n.]')
        call write_iter_header(out_populations, '[' // trim(units_abbrev(units_out%time)) // ']')
        call write_iter_nl(out_populations)

        call td_write_print_header_end(out_populations)
      end if

      ! cannot call write_iter_start, for the step is not 1
      call write_iter_int(out_populations, iter, 1)
      call write_iter_double(out_populations, units_from_atomic(units_out%time, iter*dt),  1)
      call write_iter_double(out_populations, real(gsp),  1)
      call write_iter_double(out_populations, aimag(gsp), 1)
      do ist = 1, n_excited_states
        call write_iter_double(out_populations, real(excited_state_p(ist)),  1)
        call write_iter_double(out_populations, aimag(excited_state_p(ist)), 1)
      end do
      call write_iter_nl(out_populations)
    end if

    if(n_excited_states > 0) then
      SAFE_DEALLOCATE_A(excited_state_p)
    end if
    SAFE_DEALLOCATE_A(dotprodmatrix)
    POP_SUB(td_write_populations)
  end subroutine td_write_populations



  ! ---------------------------------------------------------
  subroutine td_write_acc(out_acc, gr, geo, st, hm, dt, iter)
    type(c_ptr),         intent(in)    :: out_acc
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    integer :: idim
    character(len=7) :: aux
    FLOAT :: acc(MAX_DIM)

    PUSH_SUB(td_write_acc)

    if(iter == 0 .and. mpi_grp_is_root(mpi_world)) then
      call td_write_print_header_init(out_acc)

      ! first line -> column names
      call write_iter_header_start(out_acc)
      do idim = 1, gr%mesh%sb%dim
        write(aux, '(a4,i1,a1)') 'Acc(', idim, ')'
        call write_iter_header(out_acc, aux)
      end do
      call write_iter_nl(out_acc)

      ! second line: units
      call write_iter_string(out_acc, '#[Iter n.]')
      call write_iter_header(out_acc, '[' // trim(units_abbrev(units_out%time)) // ']')
      do idim = 1, gr%mesh%sb%dim
        call write_iter_header(out_acc, '[' // trim(units_abbrev(units_out%acceleration)) // ']')
      end do
      call write_iter_nl(out_acc)
      call td_write_print_header_end(out_acc)
    end if

    call td_calc_tacc(gr, geo, st, hm, acc, dt*iter)

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_acc)
      acc = units_from_atomic(units_out%acceleration, acc)
      call write_iter_double(out_acc, acc, gr%mesh%sb%dim)
      call write_iter_nl(out_acc)
    end if

    POP_SUB(td_write_acc)
  end subroutine td_write_acc
  
  ! ---------------------------------------------------------
  subroutine td_write_vel(out_vel, gr, geo, st, hm, dt, iter)
    type(c_ptr),         intent(in)    :: out_vel
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    integer :: idim
    character(len=7) :: aux
    FLOAT :: vel(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_vel)

    if(iter == 0) then
      call td_write_print_header_init(out_vel)

      ! first line -> column names
      call write_iter_header_start(out_vel)
      do idim = 1, gr%mesh%sb%dim
        write(aux, '(a4,i1,a1)') 'Vel(', idim, ')'
        call write_iter_header(out_vel, aux)
      end do
      call write_iter_nl(out_vel)

      ! second line: units
      call write_iter_string(out_vel, '#[Iter n.]')
      call write_iter_header(out_vel, '[' // trim(units_abbrev(units_out%time)) // ']')
      do idim = 1, gr%mesh%sb%dim
        call write_iter_header(out_vel, '[' // trim(units_abbrev(units_out%velocity)) // ']')
      end do
      call write_iter_nl(out_vel)
      call td_write_print_header_end(out_vel)
    end if

    call td_calc_tvel(gr, geo, st, hm, vel, dt*iter)

    call write_iter_start(out_vel)
    vel = units_from_atomic(units_out%velocity, vel)
    call write_iter_double(out_vel, vel, gr%mesh%sb%dim)
    call write_iter_nl(out_vel)

    POP_SUB(td_write_vel)
  end subroutine td_write_vel


  ! ---------------------------------------------------------
  subroutine td_write_laser(out_laser, gr, hm, dt, iter)
    type(c_ptr),         intent(in) :: out_laser
    type(grid_t),        intent(in) :: gr
    type(hamiltonian_t), intent(in) :: hm
    FLOAT,               intent(in) :: dt
    integer,             intent(in) :: iter

    integer :: il, idir
    FLOAT :: field(MAX_DIM)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    ! no PUSH SUB, called to often

    if(iter == 0) then
      call td_write_print_header_init(out_laser)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', units_from_atomic(units_out%time, dt), &
        " [", trim(units_abbrev(units_out%time)), "]"
      call write_iter_string(out_laser, aux)
      call write_iter_nl(out_laser)

      call write_iter_header_start(out_laser)
      do il = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(il)))
        case(E_FIELD_ELECTRIC)
          do idir = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'E(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_MAGNETIC)
          do idir = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'B(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          do idir = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'A(', idir, ')'
            call write_iter_header(out_laser, aux)
          end do
        end select
      end do
      call write_iter_nl(out_laser)

      call write_iter_string(out_laser, '#[Iter n.]')
      call write_iter_header(out_laser, '[' // trim(units_abbrev(units_out%time)) // ']')

      ! Note that we do not print out units of E, B, or A, but rather units of e*E, e*B, e*A.
      ! (force, force, and energy, respectively). The reason is that the units of E, B or A 
      ! are ugly.
      do il = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(il)))
        case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
          aux = '[' // trim(units_abbrev(units_out%force)) // ']'
          do idir = 1, gr%mesh%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          aux = '[' // trim(units_abbrev(units_out%energy)) // ']'
          do idir = 1, gr%mesh%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        end select
      end do
      call write_iter_nl(out_laser)

      call td_write_print_header_end(out_laser)
    end if

    call write_iter_start(out_laser)

    do il = 1, hm%ep%no_lasers
      field = M_ZERO
      call laser_field(hm%ep%lasers(il), gr%sb, field, iter*dt)
      select case(laser_kind(hm%ep%lasers(il)))
      case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
        field = units_from_atomic(units_out%force, field)
      case(E_FIELD_VECTOR_POTENTIAL)
        field = units_from_atomic(units_out%energy, field)
      end select
      call write_iter_double(out_laser, field, gr%mesh%sb%dim)
    end do

    call write_iter_nl(out_laser)

  end subroutine td_write_laser


  ! ---------------------------------------------------------
  subroutine td_write_energy(out_energy, hm, iter, ke)
    type(c_ptr),         intent(in) :: out_energy
    type(hamiltonian_t), intent(in) :: hm
    integer,             intent(in) :: iter
    FLOAT,               intent(in) :: ke

    integer :: ii

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_energy)

    if(iter == 0) then
      call td_write_print_header_init(out_energy)

      ! first line -> column names
      call write_iter_header_start(out_energy)
      call write_iter_header(out_energy, 'Total')
      call write_iter_header(out_energy, 'Kinetic (ions)')
      call write_iter_header(out_energy, 'Ion-Ion')
      call write_iter_header(out_energy, 'Electronic')
      call write_iter_header(out_energy, 'Eigenvalues')
      call write_iter_header(out_energy, 'Hartree')
      call write_iter_header(out_energy, 'Int[n v_xc]')
      call write_iter_header(out_energy, 'Exchange')
      call write_iter_header(out_energy, 'Correlation')
      call write_iter_nl(out_energy)

      ! second line: units
      call write_iter_string(out_energy, '#[Iter n.]')
      call write_iter_header(out_energy, '[' // trim(units_abbrev(units_out%time)) // ']')
      do ii = 1, 7
        call write_iter_header(out_energy, '[' // trim(units_abbrev(units_out%energy)) // ']')
      end do
      call write_iter_nl(out_energy)
      call td_write_print_header_end(out_energy)
    end if

    call write_iter_start(out_energy)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%total+ke), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, ke), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%ep%eii), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%total-hm%ep%eii), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%eigenvalues), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%hartree), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%intnvxc), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%exchange), 1)
    call write_iter_double(out_energy, units_from_atomic(units_out%energy, hm%energy%correlation), 1)
    call write_iter_nl(out_energy)

    POP_SUB(td_write_energy)
  end subroutine td_write_energy

  ! ---------------------------------------------------------
  subroutine td_write_gauge_field(out_gauge, hm, gr, iter)
    type(c_ptr),         intent(in) :: out_gauge
    type(hamiltonian_t), intent(in) :: hm
    type(grid_t),        intent(in) :: gr
    integer,             intent(in) :: iter
    
    integer :: idir
    character(len=50) :: aux
    FLOAT :: temp(1:MAX_DIM)
    
    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_gauge_field)
    
    if(iter == 0) then
      call td_write_print_header_init(out_gauge)

      ! first line: column names
      call write_iter_header_start(out_gauge)

      do idir = 1, gr%mesh%sb%dim
        write(aux, '(a2,i1,a1)') 'A(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do idir = 1, gr%mesh%sb%dim
        write(aux, '(a6,i1,a1)') 'dA/dt(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do idir = 1, gr%mesh%sb%dim
        write(aux, '(a10,i1,a1)') 'd^2A/dt^2(', idir, ')'
        call write_iter_header(out_gauge, aux)
      end do
      call write_iter_nl(out_gauge)

      ! second line: units
      !call write_iter_string(out_gauge, '#[Iter n.]')
      !call write_iter_header(out_gauge, '[' // trim(units_abbrev(units_out%time)) // ']')
      !call write_iter_string(out_gauge, &
      !  'A Vector potential in '   // trim(units_abbrev(units_out%length)) &
      !  'A dot in '                // trim(units_abbrev(units_out%length)) &
      !  'A dot dot in '            // trim(units_abbrev(units_out%length))
      !call write_iter_nl(out_gauge)

      call td_write_print_header_end(out_gauge)
    end if

    call write_iter_start(out_gauge)

    ! this complicated stuff here is a workaround for a PGI compiler bug in versions before 9.0-3
    temp = gauge_field_get_vec_pot(hm%ep%gfield)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    temp = gauge_field_get_vec_pot_vel(hm%ep%gfield)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy / units_out%time, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    temp = gauge_field_get_vec_pot_acc(hm%ep%gfield)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy / units_out%time**2, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    call write_iter_nl(out_gauge)
    POP_SUB(td_write_gauge_field)
    
  end subroutine td_write_gauge_field

  ! ---------------------------------------------------------
  subroutine td_write_proj(out_proj, gr, geo, st, gs_st, kick, iter)
    type(c_ptr),       intent(in) :: out_proj
    type(grid_t),      intent(in) :: gr
    type(geometry_t),  intent(in) :: geo
    type(states_t),    intent(in) :: st
    type(states_t),    intent(in) :: gs_st
    type(kick_t),      intent(in) :: kick
    integer,           intent(in) :: iter

    CMPLX, allocatable :: projections(:,:,:)
    character(len=80) :: aux
    integer :: ik, ist, uist, idir

    PUSH_SUB(td_write_proj)

    if(iter == 0) then
      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_init(out_proj)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        call kick_write(kick, out = out_proj)

        call write_iter_string(out_proj, "#%")
        call write_iter_nl(out_proj)
 
        write(aux, '(a,i8)') "# nik  ", st%d%nik
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        write(aux, '(a,2i8)') "#  st  ", gs_st%st_start, st%nst
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        write(aux, '(a,2i8)') "# ust  ", gs_st%st_start, gs_st%st_end
        call write_iter_string(out_proj, aux)
        call write_iter_nl(out_proj)

        do ik = 1, st%d%nik
          call write_iter_string(out_proj, "# w(ik)*occ(ist,ik)  ")
          do ist = gs_st%st_start, st%nst
            call write_iter_double(out_proj, st%d%kweights(ik)*st%occ(ist, ik), 1)
          end do
          call write_iter_nl(out_proj)
        end do

        call write_iter_header_start(out_proj)
        do ik = 1, st%d%nik
          do ist = gs_st%st_start, st%nst
            do uist = gs_st%st_start, gs_st%st_end
              write(aux, '(i4,a,i4)') ist, ' -> ', uist
              call write_iter_header(out_proj, 'Re {'//trim(aux)//'}')
              call write_iter_header(out_proj, 'Im {'//trim(aux)//'}')
            end do
          end do
        end do
        call write_iter_nl(out_proj)

      end if

      SAFE_ALLOCATE(projections(gs_st%st_start:st%nst, gs_st%st_start:gs_st%st_end, 1:st%d%nik))
      do idir = 1, 3
        projections(:,:,:) = M_Z0
        call dipole_matrix_elements(idir)

        if(mpi_grp_is_root(mpi_world)) then
          write(aux, '(a,i1,a)') "<i|x_", idir, "|a>"
          call write_iter_string(out_proj, "# ------")
          call write_iter_header(out_proj, aux)
          do ik = 1, st%d%nik
            do ist = gs_st%st_start, st%nst
              do uist = gs_st%st_start, gs_st%st_end
                call write_iter_double(out_proj,  real(projections(ist, uist, ik)), 1)
                call write_iter_double(out_proj, aimag(projections(ist, uist, ik)), 1)
              end do
            end do
          end do
          call write_iter_nl(out_proj)
          
        end if
      end do
      SAFE_DEALLOCATE_A(projections)

      if(mpi_grp_is_root(mpi_world)) then
        call td_write_print_header_end(out_proj)
      end if

    end if

    SAFE_ALLOCATE(projections(gs_st%st_start:st%nst, gs_st%st_start:gs_st%st_end, 1:st%d%nik))
    projections(:,:,:) = M_Z0
    call calc_projections()

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_proj)
      do ik = 1, st%d%nik
        do ist = gs_st%st_start, st%nst
          do uist = gs_st%st_start, gs_st%st_end
            call write_iter_double(out_proj,  real(projections(ist, uist, ik)), 1)
            call write_iter_double(out_proj, aimag(projections(ist, uist, ik)), 1)
          end do
        end do
      end do
      call write_iter_nl(out_proj)
    end if

    SAFE_DEALLOCATE_A(projections)
    POP_SUB(td_write_proj)

  contains
    ! ---------------------------------------------------------
    ! This subroutine calculates:
    ! p(uist, ist, ik) = < phi0(uist, k) | phi(ist, ik) (t) >
    ! ---------------------------------------------------------
    subroutine calc_projections()
      integer :: uist, ist, ik

      PUSH_SUB(td_write_proj.calc_projections)

      do ik = 1, st%d%nik
        do ist = max(gs_st%st_start, st%st_start), st%st_end
          do uist = gs_st%st_start, gs_st%st_end
            projections(ist, uist, ik) = &
              zmf_dotp(gr%mesh, st%d%dim, st%zpsi(:, :, ist, ik), gs_st%zpsi(:, :, uist, ik))
          end do
        end do
      end do
      
      call distribute_projections()

      POP_SUB(td_write_proj.calc_projections)
    end subroutine calc_projections


    ! ---------------------------------------------------------
    subroutine dipole_matrix_elements(dir)
      integer, intent(in) :: dir

      integer :: uist, ist, ik, idim
      FLOAT   :: n_dip(MAX_DIM)
      CMPLX, allocatable :: xpsi(:,:)

      PUSH_SUB(td_write_proj.dipole_matrix_elements)

      call geometry_dipole(geo, n_dip)

      SAFE_ALLOCATE(xpsi(1:gr%mesh%np, 1:st%d%dim))
      
      do ik = 1, st%d%nik
        do ist = max(gs_st%st_start, st%st_start), st%st_end
          do uist = gs_st%st_start, gs_st%st_end
            
            do idim = 1, st%d%dim
              xpsi(1:gr%mesh%np, idim) = gr%mesh%x(1:gr%mesh%np, dir) * gs_st%zpsi(1:gr%mesh%np, idim, uist, ik)
            end do

            projections(ist, uist, ik) = -n_dip(dir) - &
              zmf_dotp(gr%mesh, st%d%dim, st%zpsi(:, :, ist, ik), xpsi(:, :))

          end do
        end do
      end do
      
      SAFE_DEALLOCATE_A(xpsi)

      call distribute_projections()

      POP_SUB(td_write_proj.dipole_matrix_elements)
    end subroutine dipole_matrix_elements

    subroutine distribute_projections
#if defined(HAVE_MPI)
      integer :: k, ik, ist, uist

      if(.not.st%parallel_in_states) return

      PUSH_SUB(td_write_proj.distribute_projections)

      do ik = 1, st%d%nik
        do ist = gs_st%st_start, st%nst
          k = st%node(ist)
          do uist = gs_st%st_start, gs_st%st_end
            call MPI_Bcast(projections(ist, uist, ik), 1, MPI_CMPLX, k, st%mpi_grp%comm, mpi_err)
          end do
        end do
      end do

      POP_SUB(td_write_proj.distribute_projections)
#endif
    end subroutine distribute_projections

  end subroutine td_write_proj


  ! ---------------------------------------------------------
  subroutine td_write_print_header_init(out)
    type(c_ptr), intent(in) :: out

    PUSH_SUB(td_write_print_header_init)

    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_init)
  end subroutine td_write_print_header_init


  ! ---------------------------------------------------------
  subroutine td_write_print_header_end(out)
    type(c_ptr), intent(in) :: out

    PUSH_SUB(td_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_end)
  end subroutine td_write_print_header_end


#include "td_calc_inc.F90"

end module td_write_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
