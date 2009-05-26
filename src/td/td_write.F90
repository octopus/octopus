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
  use external_pot_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use ion_dynamics_m
  use lasers_m
  use loct_m
  use loct_parser_m
  use magnetic_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use h_sys_output_m
  use pert_m
  use restart_m
  use spectrum_m
  use states_m
  use states_dim_m
  use states_calc_m
  use profiling_m
  use units_m
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
    OUT_MAX         = 12
  
  type td_write_t
    private
    type(td_write_prop_t) :: out(OUT_MAX)

    integer        :: lmax     ! maximum multipole moment to output
    FLOAT          :: lmm_r    ! radius of the sphere used to compute the local magnetic moments
    type(states_t) :: gs_st    ! The states_type where the ground state is stored, in order to
                                        ! calculate the projections(s) onto it.
    integer        :: n_excited_states  ! number of excited sates onto which the projections are calculated.
    type(excited_states_t), pointer :: excited_st(:) ! The excited states.
  end type td_write_t

contains

  ! ---------------------------------------------------------
  subroutine td_write_init(w, gr, st, hm, geo, ions_move, with_gauge_field, iter, max_iter, dt)
    type(td_write_t), intent(out) :: w
    type(grid_t),     intent(in)  :: gr
    type(states_t),   intent(in)  :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t), intent(in)  :: geo
    logical,          intent(in)  :: ions_move
    logical,          intent(in)  :: with_gauge_field
    integer,          intent(in)  :: iter
    integer,          intent(in)  :: max_iter
    FLOAT,            intent(in)  :: dt


    FLOAT :: rmin
    integer :: ierr, first, i, j, flags, iout, default
    type(block_t) :: blk
    character(len=100) :: filename

    call push_sub('td_write.td_write_init')


    !%Variable TDOutput
    !%Type flag
    !%Default multipoles + geometry + temperature + energy
    !%Section Time Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent simulation.
    !%Option multipoles 1
    !% Outputs the multipole moments of the density to the file <tt>td.general/multipoles</tt>.
    !% This is required to, e.g., calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>TDDipoleLmax</tt>.
    !%Option angular 2
    !% Outputs the angular momentum of the system that can be used to calculate circular
    !% dichroism (EXPERIMENTAL)
    !%Option spin 4
    !% Outputs the expectation value of the spin, that can be used to calculate magnetic
    !% cicular dichroism (EXPERIMENTAL)
    !%Option populations 8
    !% Outputs the projection of the time-dependent Kohn-Sham Slater determinant
    !% onto the ground-state (or approximations to the excited states) to the file 
    !% <tt>td.general/populations</tt>.
    !%Option geometry 16
    !% If set (and if the atoms are allowed to move), outputs the coordinates, velocities,
    !% and forces of the atoms to the the file <tt>td.general/coordinates</tt>.
    !%Option acceleration 32
    !% When set outputs the acceleration, calculated from Ehrenfest theorem,
    !% in the file <tt>td.general/acceleration</tt>. This file can then be
    !% processed by the utility "hs-from-acc" in order to obtain the harmonic spectrum.
    !%Option laser 64
    !% If set, and if there are lasers defined in <tt>TDLasers</tt>,
    !% <tt>octopus</tt> outputs the laser field to the file <tt>td.general/laser</tt>.
    !%Option energy 128
    !% If <tt>set</tt>, <tt>octopus</tt> outputs the different components of the energy
    !% to the file <tt>td.general/el_energy</tt>.
    !%Option td_occup 256
    !% If set, outputs the projections of the time-dependent Kohn-Sham
    !% wave-functions onto the static (zero time) wave-functions to the
    !% file <tt>td.general/projections.XXX</tt>.
    !%Option local_mag_moments 512
    !% If set, outputs the local magnetic moments, integrated in sphere centered around each atom.
    !% The radius of the sphere can be ser with <tt>LocalMagneticMomentsSphereRadius</tt>
    !%Option gauge_field 1024
    !% If set, outputs the vector gauge field corresponding to a uniform (but time dependent) 
    !% external electrical potential. This is only useful in a time-dependent periodic run
    !%Option temperature 2048
    !% If set, the ionic temperature at each step is printed.
    !%End

    ! by default print multipoles, coordinates and energy
    default = &
         2**(OUT_MULTIPOLES - 1) +  &
         2**(OUT_COORDS - 1) +      &
         2**(OUT_TEMPERATURE - 1) + &
         2**(OUT_ENERGY - 1) +      &
         2**(OUT_GAUGE_FIELD - 1)

    call loct_parse_int(datasets_check('TDOutput'), default, flags)

    if(.not.varinfo_valid_option('TDOutput', flags, is_flag = .true.)) call input_error('TDOutput')

    do iout = 1, OUT_MAX
      w%out(iout)%write = (iand(flags, 2**(iout - 1)) .ne. 0)
    end do

    !special cases
    w%out(OUT_COORDS)%write = w%out(OUT_COORDS)%write .and. ions_move
    w%out(OUT_TEMPERATURE)%write = w%out(OUT_TEMPERATURE)%write .and. ions_move
    w%out(OUT_GAUGE_FIELD)%write = w%out(OUT_GAUGE_FIELD)%write .and. with_gauge_field
    w%out(OUT_LASER)%write = w%out(OUT_LASER)%write .and. (hm%ep%no_lasers > 0)

    !%Variable TDDipoleLmax
    !%Type integer
    !%Default 1
    !%Section Time Dependent::TD Output
    !%Description
    !% Maximum multi-pole of the density output to the file <tt>td.general/multipoles</tt>
    !% during a time-dependent simulation. Must be 0 &lt; <tt>TDDipoleLmax &lt; 5</tt>.
    !%End
    call loct_parse_int(datasets_check('TDDipoleLmax'), 1, w%lmax)
    if (w%lmax < 0 .or. w%lmax > 4) then
      write(message(1), '(a,i6,a)') "Input: '", w%lmax, "' is not a valid TDDipoleLmax"
      message(2) = '(0 <= TDDipoleLmax <= 4 )'
      call write_fatal(2)
    end if

    ! Compatibility test
    if( (w%out(OUT_ACC)%write) .and. ions_move ) then
      message(1) = 'Error. If harmonic spectrum is to be calculated'
      message(2) = 'Atoms should not be allowed to move'
      call write_fatal(2)
    end if

    call geometry_min_distance(geo, rmin)
    call loct_parse_float(datasets_check('LocalMagneticMomentsSphereRadius'), rmin*M_HALF/units_inp%length%factor, w%lmm_r)
    w%lmm_r = w%lmm_r * units_inp%length%factor

    if( (w%out(OUT_PROJ)%write)  .or.  (w%out(OUT_POPULATIONS)%write) ) then
      call states_copy(w%gs_st, st)

      ! clean up all the stuff we have to reallocate
      SAFE_DEALLOCATE_P(w%gs_st%zpsi)
      SAFE_DEALLOCATE_P(w%gs_st%occ)
      SAFE_DEALLOCATE_P(w%gs_st%eigenval)
      SAFE_DEALLOCATE_P(w%gs_st%node)
      if(w%gs_st%d%ispin == SPINORS) then
        SAFE_DEALLOCATE_P(w%gs_st%spin)
      end if

      call states_look (trim(restart_dir)//'gs', gr%mesh%mpi_grp, i, j, w%gs_st%nst, ierr)

      if(w%out(OUT_POPULATIONS)%write) then ! do only this when not calculating populations
        ! We will store the ground-state Kohn-Sham system by all processors.
        !%Variable TDProjStateStart
        !%Type integer
        !%Default 1
        !%Section Time Dependent::TD Output
        !%Description
        !% Only output projections to states above TDProjStateStart. Usually one is only interested
        !% in particle-hole projections around the HOMO, so there is no need to calculate (and store)
        !% the projections of all static onto all TD states. This sets a lower limit. The upper limit
        !% is set by the number of states in the propagation and the number of uncoccupied states
        !% available
        !%End
        call loct_parse_int(datasets_check('TDProjStateStart'), 1, w%gs_st%st_start)
      else
        w%gs_st%st_start = 1
      end if
      w%gs_st%st_end   = w%gs_st%nst

      ! allocate memory
      SAFE_ALLOCATE(w%gs_st%occ(1:w%gs_st%nst, 1:w%gs_st%d%nik))
      SAFE_ALLOCATE(w%gs_st%eigenval(1:w%gs_st%nst, 1:w%gs_st%d%nik))
      SAFE_ALLOCATE(w%gs_st%node(1:w%gs_st%nst))
      if(w%gs_st%d%ispin == SPINORS) then
        SAFE_ALLOCATE(w%gs_st%spin(1:3, 1:w%gs_st%nst, 1:w%gs_st%d%nik))
      end if
      call states_allocate_wfns(w%gs_st, gr%mesh, M_CMPLX)
      w%gs_st%node(:)  = 0
      call restart_read(trim(restart_dir)//'gs', w%gs_st, gr, geo, ierr)
      if(ierr.ne.0 .and.ierr.ne.(w%gs_st%st_end-w%gs_st%st_start+1)*w%gs_st%d%nik*w%gs_st%d%dim) then
        message(1) = "Could not load "//trim(restart_dir)//"gs"
        call write_fatal(1)
      end if
    end if

    ! Build the excited states...
    if(w%out(OUT_POPULATIONS)%write) then
      !%Variable TDExcitedStatesToProject
      !%Type block
      !%Section Time Dependent::TD Output
      !%Description
      !% [WARNING: This is a *very* experimental feature] The population of the excited states
      !% (as defined by <Phi_I|Phi(t)> where |Phi(t)> is the many-body time-dependent state at
      !% time t, and |Phi_I> is the excited state of interest) can be approximated -- it is not clear 
      !% how well--  by substituting those real many-body states by the time-dependent Kohn-Sham
      !% determinant and by some modification of the Kohn-Sham ground state determinant (e.g.,
      !% a simple HOMO-LUMO substitution, or the Casida ansatz for excited states in linear
      !% response theory. If you set TDOutput to contain, you may ask for these approximated
      !% populations for a number of excited states, which will be described in the files specified
      !% in this block: each line should be the name of a file that contains one excited state.
      !%
      !% FIXME: description of the format of the files.
      !%End
      if(loct_parse_block('TDExcitedStatesToProject', blk) == 0) then
        w%n_excited_states = loct_parse_block_n(blk)
        SAFE_ALLOCATE(w%excited_st(1:w%n_excited_states))
        do i = 1, w%n_excited_states
          call loct_parse_block_string(blk, i-1, 0, filename)
          call excited_states_init(w%excited_st(i), w%gs_st, trim(filename)) 
        end do
      else
        w%n_excited_states = 0
        nullify(w%excited_st)
      end if
    end if

    if (iter == 0) then
      first = 0
    else
      first = iter + 1
    end if

    call io_mkdir('td.general')

    if(mpi_grp_is_root(mpi_world)) then
      if(w%out(OUT_MULTIPOLES)%write) &
        call write_iter_init(w%out(OUT_MULTIPOLES)%handle, &
        first, dt/units_out%time%factor, trim(io_workpath("td.general/multipoles")))

      if(w%out(OUT_ANGULAR)%write) &
        call write_iter_init(w%out(OUT_ANGULAR)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/angular")))

      if(w%out(OUT_SPIN)%write) &
        call write_iter_init(w%out(OUT_SPIN)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/spin")))

      if(w%out(OUT_MAGNETS)%write) &
        call write_iter_init(w%out(OUT_MAGNETS)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/magnetic_moments")))

      if(w%out(OUT_COORDS)%write) &
        call write_iter_init(w%out(OUT_COORDS)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/coordinates")))

      if(w%out(OUT_TEMPERATURE)%write) &
        call write_iter_init(w%out(OUT_TEMPERATURE)%handle, first, M_ONE, trim(io_workpath("td.general/temperature")))

      if(w%out(OUT_POPULATIONS)%write) &
        call write_iter_init(w%out(OUT_POPULATIONS)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/populations")))

      if(w%out(OUT_ACC)%write) &
        call write_iter_init(w%out(OUT_ACC)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/acceleration")))

      if(w%out(OUT_LASER)%write) then
        call write_iter_init(w%out(OUT_LASER)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/laser")))
        do i = 0, max_iter
          call td_write_laser(w%out(OUT_LASER)%handle, gr, hm, dt, i)
          if(mod(i, 100).eq.0) call write_iter_flush(w%out(OUT_LASER)%handle)
        end do
        call write_iter_end(w%out(OUT_LASER)%handle)
      end if


      if(w%out(OUT_ENERGY)%write) &
        call write_iter_init(w%out(OUT_ENERGY)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/energy")))

      if(w%out(OUT_PROJ)%write) &
        call write_iter_init(w%out(OUT_PROJ)%handle, first, &
          dt/units_out%time%factor, trim(io_workpath("td.general/projections")))

      if(w%out(OUT_GAUGE_FIELD)%write) &
        call write_iter_init(w%out(OUT_GAUGE_FIELD)%handle, &
        first, dt/units_out%time%factor, trim(io_workpath("td.general/gauge_field")))
    end if

    call pop_sub()
  end subroutine td_write_init


  ! ---------------------------------------------------------
  subroutine td_write_end(w)
    type(td_write_t), intent(inout) :: w
    integer :: i, iout
    call push_sub('td_write.td_write_end')

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAX
        if(iout.eq.OUT_LASER) cycle
        if(w%out(iout)%write)  call write_iter_end(w%out(iout)%handle)
      end do
    end if

    if( w%out(OUT_POPULATIONS)%write ) then
      do i = 1, w%n_excited_states
        call excited_states_kill(w%excited_st(i))
      end do
    end if

    if(w%out(OUT_POPULATIONS)%write .or. w%out(OUT_PROJ)%write) call states_end(w%gs_st)

    call pop_sub()
  end subroutine td_write_end


  ! ---------------------------------------------------------
  subroutine td_write_iter(w, gr, st, hm, geo, kick, dt, i)
    type(td_write_t),       intent(in) :: w
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t),    intent(inout) :: geo
    type(kick_t),           intent(in) :: kick
    FLOAT,                  intent(in) :: dt
    integer,                intent(in) :: i

    type(profile_t), save :: prof

    call push_sub('td_write.td_write_iter')
    call profiling_in(prof, "TD_WRITE_ITER")

    if(w%out(OUT_MULTIPOLES)%write) &
      call td_write_multipole(w%out(OUT_MULTIPOLES)%handle, gr, geo, st, w%lmax, kick, i)
    
    if(w%out(OUT_ANGULAR)%write) &
      call td_write_angular(w%out(OUT_ANGULAR)%handle, gr, geo, hm, st, kick, i)

    if(w%out(OUT_SPIN)%write) &
      call td_write_spin(w%out(OUT_SPIN)%handle, gr, st, i)

    if(w%out(OUT_MAGNETS)%write) &
      call td_write_local_magnetic_moments(w%out(OUT_MAGNETS)%handle, gr, st, geo, w%lmm_r, i)

    if(w%out(OUT_PROJ)%write) &
      call td_write_proj(w%out(OUT_PROJ)%handle, gr, geo, st, w%gs_st, kick, i)

    if(w%out(OUT_COORDS)%write) &
      call td_write_coordinates(w%out(OUT_COORDS)%handle, gr, geo, i)

    if(w%out(OUT_TEMPERATURE)%write) &
      call td_write_temperature(w%out(OUT_TEMPERATURE)%handle, geo, i)

    if(w%out(OUT_POPULATIONS)%write) &
      call td_write_populations(w%out(OUT_POPULATIONS)%handle, gr%mesh, st, &
        w%gs_st, w%n_excited_states, w%excited_st, dt, i)

    if(w%out(OUT_ACC)%write) &
      call td_write_acc(w%out(OUT_ACC)%handle, gr, geo, st, hm, dt, i)

    ! td_write_laser no longer called here, because the whole laser is printed
    ! out at the beginning.

    if(w%out(OUT_ENERGY)%write) &
      call td_write_energy(w%out(OUT_ENERGY)%handle, hm, i, geo%kinetic_energy)

    if(w%out(OUT_GAUGE_FIELD)%write) &
      call td_write_gauge_field(w%out(OUT_GAUGE_FIELD)%handle, hm, gr, i)

    call profiling_out(prof)
    call pop_sub()
  end subroutine td_write_iter


  ! ---------------------------------------------------------
  subroutine td_write_data(w, gr, st, hm, outp, geo, iter)
    type(td_write_t),     intent(in)    :: w
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(in)    :: hm
    type(h_sys_output_t), intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: iter

    character(len=256) :: filename
    integer :: iout
    type(profile_t), save :: prof

    call push_sub('td.td_write_data')
    call profiling_in(prof, "TD_WRITE_DATA")

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAX
        if(iout.eq.OUT_LASER) cycle
        if(w%out(iout)%write)  call write_iter_flush(w%out(iout)%handle)
      end do
    end if

    ! now write down the rest
    write(filename, '(a,i7.7)') "td.", iter  ! name of directory

    call h_sys_output_all(outp, gr, geo, st, hm, filename)

    call profiling_out(prof)
    call pop_sub()
  end subroutine td_write_data


  ! ---------------------------------------------------------
  subroutine td_write_spin(out_spin, gr, st, iter)
    type(c_ptr), intent(in)    :: out_spin
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(in)    :: st
    integer,           intent(in)    :: iter

    character(len=130) :: aux
    FLOAT :: spin(3)

    call push_sub('td_write.td_write_spin')

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

    call pop_sub()
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

    call push_sub('td_write.td_write_local_magnetic_moments')

    !get the atoms magnetization. This has to be calculated by all nodes
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

    call pop_sub()
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
    FLOAT :: angular(MAX_DIM), lsquare
    FLOAT, allocatable :: ang(:, :, :), ang2(:, :)
    type(pert_t)        :: angular_momentum

    call push_sub('td_write.td_write_angular')

    SAFE_ALLOCATE(ang (st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end, 1:3))
    SAFE_ALLOCATE(ang2(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        call zstates_angular_momentum(gr, st%zpsi(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
      end do
    end do
    lsquare    =  states_eigenvalues_sum(st, ang2)

    call pert_init(angular_momentum, PERTURBATION_MAGNETIC, gr, geo)

    do idir = 1, 3
       call pert_setup_dir(angular_momentum, idir)
       !we have to multiply by 2, because is the perturbation returns L/2
       angular(idir) = M_TWO * zpert_expectation_value(angular_momentum, gr, geo, hm, st, st%zpsi, st%zpsi)
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
        write(aux, '(a4,18x)') '<L2>'
        call write_iter_header(out_angular, aux)
        call write_iter_nl(out_angular)

        !third line -> should hold the units. Now unused (assumes atomic units)
        call write_iter_string(out_angular, '#[Iter n.]')
        call write_iter_header(out_angular, '[' // trim(units_out%time%abbrev) // ']')
        call write_iter_nl(out_angular)

        call td_write_print_header_end(out_angular)
      end if

      call write_iter_start(out_angular)
      call write_iter_double(out_angular, angular(1:3), 3)
      call write_iter_double(out_angular, lsquare, 1)
      call write_iter_nl(out_angular)

    end if

    SAFE_DEALLOCATE_A(ang)
    SAFE_DEALLOCATE_A(ang2)
    call pop_sub()
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

    integer :: is, l, m, add_lm
    character(len=120) :: aux
    FLOAT, allocatable :: nuclear_dipole(:), multipole(:,:)

    call push_sub('td_write.td_write_multipole')

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
        do l = 2, lmax
          do m = -l, l
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
            call write_iter_header(out_multip, aux)
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      ! units
      call write_iter_string(out_multip, '#[Iter n.]')
      call write_iter_header(out_multip, '[' // trim(units_out%time%abbrev) // ']')

      do is = 1, st%d%nspin
        do l = 0, lmax
          do m = -l, l
            select case(l)
            case(0)
              call write_iter_header(out_multip, 'Electrons')
            case(1)
              call write_iter_header(out_multip, '[' // trim(units_out%length%abbrev) // ']')
            case default
              write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
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
        do l = 0, lmax
          do m = -l, l
            call write_iter_double(out_multip, multipole(add_lm, is)/units_out%length%factor**l, 1)
            add_lm = add_lm + 1
          end do
        end do
      end do
      call write_iter_nl(out_multip)
    end if

    SAFE_DEALLOCATE_A(nuclear_dipole)
    SAFE_DEALLOCATE_A(multipole)
    call pop_sub()
  end subroutine td_write_multipole


  ! ---------------------------------------------------------
  subroutine td_write_coordinates(out_coords, gr, geo, iter)
    type(c_ptr),       intent(in) :: out_coords
    type(grid_t),      intent(in) :: gr
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iter

    integer :: i, j
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    call push_sub('td_write.td_write_coordinates')

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do i = 1, geo%natoms
        do j = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, geo%natoms
        do j = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, geo%natoms
        do j = 1, gr%mesh%sb%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_out%time%abbrev) // ']')
      call write_iter_string(out_coords, &
        'Positions in '   // trim(units_out%length%abbrev)   //   &
        ', Velocities in '// trim(units_out%velocity%abbrev) //   &
        ', Forces in '    // trim(units_out%force%abbrev))
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    do i = 1, geo%natoms
      call write_iter_double(out_coords, geo%atom(i)%x(1:gr%mesh%sb%dim)/units_out%length%factor,   gr%mesh%sb%dim)
    end do
    do i = 1, geo%natoms
      call write_iter_double(out_coords, geo%atom(i)%v(1:gr%mesh%sb%dim)/units_out%velocity%factor, gr%mesh%sb%dim)
    end do
    do i = 1, geo%natoms
      call write_iter_double(out_coords, geo%atom(i)%f(1:gr%mesh%sb%dim)/units_out%force%factor,    gr%mesh%sb%dim)
    end do
    call write_iter_nl(out_coords)

    call pop_sub()
  end subroutine td_write_coordinates


  ! ---------------------------------------------------------
  subroutine td_write_temperature(out_temperature, geo, iter)
    type(c_ptr),       intent(in) :: out_temperature
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iter

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    call push_sub('td_write.td_write_temperature')

    if(iter == 0) then
      call td_write_print_header_init(out_temperature)

      ! first line: column names
      call write_iter_header_start(out_temperature)
      call write_iter_header(out_temperature, 'Temperature')
      call write_iter_nl(out_temperature)

      ! second line: units
      call write_iter_string(out_temperature, '#[Iter n.]')
      call write_iter_header(out_temperature, '[' // trim(units_out%time%abbrev) // ']')
      call write_iter_string(out_temperature, '        [K]')
      call write_iter_nl(out_temperature)

      call td_write_print_header_end(out_temperature)
    end if

    call write_iter_start(out_temperature)

    call write_iter_double(out_temperature, ion_dynamics_temperature(geo), 1)

    call write_iter_nl(out_temperature)

    call pop_sub()
  end subroutine td_write_temperature


  ! ---------------------------------------------------------
  subroutine td_write_populations(out_populations, m, st, gs_st, n_excited_states, excited_st, dt, iter)
    type(c_ptr),            intent(in) :: out_populations
    type(mesh_t),           intent(in) :: m
    type(states_t),         intent(in) :: st
    type(states_t),         intent(in) :: gs_st
    integer,                intent(in) :: n_excited_states
    type(excited_states_t), intent(in) :: excited_st(:)
    FLOAT,                  intent(in) :: dt
    integer,                intent(in) :: iter
 
    integer :: j
    character(len=6) :: excited_name
    CMPLX :: gsp
    CMPLX, allocatable :: excited_state_p(:)
    CMPLX, allocatable :: dotprodmatrix(:, :, :)


    call push_sub('td_write.td_write_populations')

    SAFE_ALLOCATE(dotprodmatrix(1:gs_st%nst, 1:st%nst, 1:st%d%nik))
    call zstates_matrix(m, gs_st, st, dotprodmatrix)

    ! all processors calculate the projection
    gsp = zstates_mpdotp(m, gs_st, st, dotprodmatrix)

    if(n_excited_states > 0) then
      SAFE_ALLOCATE(excited_state_p(1:n_excited_states))
      do j = 1, n_excited_states
        excited_state_p(j) = zstates_mpdotp(m, excited_st(j), st, dotprodmatrix)
      end do
    end if

    if(mpi_grp_is_root(mpi_world)) then
      if(iter == 0) then
        call td_write_print_header_init(out_populations)

        ! first line -> column names
        call write_iter_header_start(out_populations)
        call write_iter_header(out_populations, 'Re<Phi_gs|Phi(t)>')
        call write_iter_header(out_populations, 'Im<Phi_gs|Phi(t)>')
        do j = 1, n_excited_states
          write(excited_name,'(a2,i3,a1)') 'P(',j,')'
          call write_iter_header(out_populations, 'Re<'//excited_name//'|Phi(t)>')
          call write_iter_header(out_populations, 'Im<'//excited_name//'|Phi(t)>')
        end do
        call write_iter_nl(out_populations)

        ! second line -> units
        call write_iter_string(out_populations, '#[Iter n.]')
        call write_iter_header(out_populations, '[' // trim(units_out%time%abbrev) // ']')
        call write_iter_nl(out_populations)

        call td_write_print_header_end(out_populations)
      end if

      ! can not call write_iter_start, for the step is not 1
      call write_iter_int(out_populations, iter, 1)
      call write_iter_double(out_populations, iter*dt/units_out%time%factor,  1)
      call write_iter_double(out_populations, real(gsp),  1)
      call write_iter_double(out_populations, aimag(gsp), 1)
      do j = 1, n_excited_states
        call write_iter_double(out_populations, real(excited_state_p(j)),  1)
        call write_iter_double(out_populations, aimag(excited_state_p(j)), 1)
      end do
      call write_iter_nl(out_populations)
    end if

    if(n_excited_states > 0) then
      SAFE_DEALLOCATE_A(excited_state_p)
    end if
    SAFE_DEALLOCATE_A(dotprodmatrix)
    call pop_sub()
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

    integer :: i
    character(len=7) :: aux
    FLOAT :: acc(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    if(iter == 0) then
      call td_write_print_header_init(out_acc)

      ! first line -> column names
      call write_iter_header_start(out_acc)
      do i = 1, gr%mesh%sb%dim
        write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
        call write_iter_header(out_acc, aux)
      end do
      call write_iter_nl(out_acc)

      ! second line: units
      call write_iter_string(out_acc, '#[Iter n.]')
      call write_iter_header(out_acc, '[' // trim(units_out%time%abbrev) // ']')
      do i = 1, gr%mesh%sb%dim
        call write_iter_header(out_acc, '[' // trim(units_out%acceleration%abbrev) // ']')
      end do
      call write_iter_nl(out_acc)
      call td_write_print_header_end(out_acc)
    end if

    call td_calc_tacc(gr, geo, st, hm, acc, dt*i)

    call write_iter_start(out_acc)
    call write_iter_double(out_acc, acc/units_out%acceleration%factor, gr%mesh%sb%dim)
    call write_iter_nl(out_acc)

  end subroutine td_write_acc


  ! ---------------------------------------------------------
  subroutine td_write_laser(out_laser, gr, hm, dt, iter)
    type(c_ptr),         intent(in) :: out_laser
    type(grid_t),        intent(in) :: gr
    type(hamiltonian_t), intent(in) :: hm
    FLOAT,               intent(in) :: dt
    integer,             intent(in) :: iter

    integer :: i, j
    FLOAT :: field(MAX_DIM)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    call push_sub('td_write.td_write_laser')

    ! TODO -> confirm these stupid units, especially for the vector field
    if(iter == 0) then
      call td_write_print_header_init(out_laser)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', dt/units_out%time%factor, &
        " [", trim(units_out%time%abbrev), "]"
      call write_iter_string(out_laser, aux)
      call write_iter_nl(out_laser)

      call write_iter_header_start(out_laser)
      do i = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(i)))
        case(E_FIELD_ELECTRIC)
          do j = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'E(', j, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_MAGNETIC)
          do j = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'B(', j, ')'
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          do j = 1, gr%mesh%sb%dim
            write(aux, '(a,i1,a)') 'A(', j, ')'
            call write_iter_header(out_laser, aux)
          end do
        end select
      end do
      call write_iter_nl(out_laser)

      call write_iter_string(out_laser, '#[Iter n.]')
      call write_iter_header(out_laser, '[' // trim(units_out%time%abbrev) // ']')

      ! Note that we do not print out units of E, B, or A, but rather units of e*E, e*B, e*A.
      ! (force, force, and energy, respectively). The reason is that the units of E, B or A 
      ! are ugly.
      do i = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(i)))
        case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
          aux = '[' // trim(units_out%energy%abbrev) // ' / ' // trim(units_inp%length%abbrev) // ']'
          do j = 1, gr%mesh%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        case(E_FIELD_VECTOR_POTENTIAL)
          aux = '[' // trim(units_out%energy%abbrev) // ']'
          do j = 1, gr%mesh%sb%dim
            call write_iter_header(out_laser, aux)
          end do
        end select
      end do
      call write_iter_nl(out_laser)

      call td_write_print_header_end(out_laser)
    end if

    call write_iter_start(out_laser)

    do i = 1, hm%ep%no_lasers
      field = M_ZERO
      call laser_field(gr%sb, hm%ep%lasers(i), field, iter*dt)
      select case(laser_kind(hm%ep%lasers(i)))
      case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
        field = field * units_inp%length%factor / units_inp%energy%factor
      case(E_FIELD_VECTOR_POTENTIAL)
        field = field / units_inp%energy%factor
      end select
      call write_iter_double(out_laser, field, gr%mesh%sb%dim)
    end do

    call write_iter_nl(out_laser)

    call pop_sub()
  end subroutine td_write_laser


  ! ---------------------------------------------------------
  subroutine td_write_energy(out_energy, hm, iter, ke)
    type(c_ptr),         intent(in) :: out_energy
    type(hamiltonian_t), intent(in) :: hm
    integer,             intent(in) :: iter
    FLOAT,               intent(in) :: ke

    integer :: i

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

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
      call write_iter_header(out_energy, '[' // trim(units_out%time%abbrev) // ']')
      do i = 1, 7
        call write_iter_header(out_energy, '[' // trim(units_out%energy%abbrev) // ']')
      end do
      call write_iter_nl(out_energy)
      call td_write_print_header_end(out_energy)
    end if

    call write_iter_start(out_energy)
    call write_iter_double(out_energy, (hm%etot+ke)/units_out%energy%factor, 1)
    call write_iter_double(out_energy, ke/units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%ep%eii /units_out%energy%factor, 1)
    call write_iter_double(out_energy, (hm%etot-hm%ep%eii)/units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%eeigen /units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%ehartree /units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%epot /units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%ex  /units_out%energy%factor, 1)
    call write_iter_double(out_energy, hm%ec  /units_out%energy%factor, 1)
    call write_iter_nl(out_energy)


  end subroutine td_write_energy

  ! ---------------------------------------------------------
  subroutine td_write_gauge_field(out_gauge, hm, gr, iter)
    type(c_ptr),         intent(in) :: out_gauge
    type(hamiltonian_t), intent(in) :: hm
    type(grid_t),        intent(in) :: gr
    integer,             intent(in) :: iter
    
    integer :: j
    character(len=50) :: aux

    
    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    call push_sub('td_write.td_write_out_gauge')
    
    if(iter == 0) then
      call td_write_print_header_init(out_gauge)

      ! first line: column names
      call write_iter_header_start(out_gauge)

      do j = 1, gr%mesh%sb%dim
        write(aux, '(a2,i1,a1)') 'A(', j, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do j = 1, gr%mesh%sb%dim
        write(aux, '(a6,i1,a1)') 'dA/dt(', j, ')'
        call write_iter_header(out_gauge, aux)
      end do
      do j = 1, gr%mesh%sb%dim
        write(aux, '(a10,i1,a1)') 'd^2A/dt^2(', j, ')'
        call write_iter_header(out_gauge, aux)
      end do
      call write_iter_nl(out_gauge)

      ! second line: units
      !call write_iter_string(out_gauge, '#[Iter n.]')
      !call write_iter_header(out_gauge, '[' // trim(units_out%time%abbrev) // ']')
      !call write_iter_string(out_gauge, &
      !  'A Vector potential in '   // trim(units_out%length%abbrev) &
      !  'A dot in '                // trim(units_out%length%abbrev) &
      !  'A dot dot in '            // trim(units_out%length%abbrev)
      !call write_iter_nl(out_gauge)

      call td_write_print_header_end(out_gauge)
    end if

    call write_iter_start(out_gauge)

    ! TODO: put the appropriate units here 
    call write_iter_double(out_gauge, gauge_field_get_vec_pot(hm%ep%gfield), gr%mesh%sb%dim)
    call write_iter_double(out_gauge, gauge_field_get_vec_pot_vel(hm%ep%gfield), gr%mesh%sb%dim)
    call write_iter_double(out_gauge, gauge_field_get_vec_pot_acc(hm%ep%gfield), gr%mesh%sb%dim)
    call write_iter_nl(out_gauge)
    call pop_sub()
    
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

    call push_sub('td_write.td_write_proj')

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
    call pop_sub()

  contains
    ! ---------------------------------------------------------
    ! This subroutine calculates:
    ! p(uist, ist, ik) = < phi0(uist, k) | phi(ist, ik) (t) >
    ! ---------------------------------------------------------
    subroutine calc_projections()
      integer :: uist, ist, ik

      do ik = 1, st%d%nik
        do ist = max(gs_st%st_start, st%st_start), st%st_end
          do uist = gs_st%st_start, gs_st%st_end
            projections(ist, uist, ik) = &
              zmf_dotp(gr%mesh, st%d%dim, st%zpsi(:, :, ist, ik), gs_st%zpsi(:, :, uist, ik))
          end do
        end do
      end do
      
      call distribute_projections()

    end subroutine calc_projections


    ! ---------------------------------------------------------
    subroutine dipole_matrix_elements(dir)
      integer, intent(in) :: dir

      integer :: uist, ist, ik, idim
      FLOAT   :: n_dip(MAX_DIM)
      CMPLX, allocatable :: xpsi(:,:)

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

    end subroutine dipole_matrix_elements

    subroutine distribute_projections
#if defined(HAVE_MPI)
      integer :: k, ik, ist, uist

      if(.not.st%parallel_in_states) return

      do ik = 1, st%d%nik
        do ist = gs_st%st_start, st%nst
          k = st%node(ist)
          do uist = gs_st%st_start, gs_st%st_end
            call MPI_Bcast(projections(ist, uist, ik), 1, MPI_CMPLX, k, st%mpi_grp%comm, mpi_err)
          end do
        end do
      end do

#endif
    end subroutine distribute_projections

  end subroutine td_write_proj


  ! ---------------------------------------------------------
  subroutine td_write_print_header_init(out)
    type(c_ptr), intent(in) :: out

    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

  end subroutine td_write_print_header_init


  ! ---------------------------------------------------------
  subroutine td_write_print_header_end(out)
    type(c_ptr), intent(in) :: out

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

  end subroutine td_write_print_header_end


#include "td_calc.F90"

end module td_write_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
