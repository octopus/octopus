! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module td_write_oct_m
  use io_function_oct_m
  use iso_c_binding
  use comm_oct_m
  use current_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use output_oct_m
  use hamiltonian_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use kick_oct_m
  use lasers_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_dim_oct_m
  use states_restart_oct_m
  use td_calc_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::         &
    td_write_t,     &
    td_write_init,  &
    td_write_end,   &
    td_write_iter,  &
    td_write_data,  &
    td_write_kick,  &
    td_write_output

  type td_write_prop_t
    private
    type(c_ptr) :: handle
    logical :: write = .false.
  end type td_write_prop_t

  integer, parameter ::   &
    OUT_MULTIPOLES  =  1, &
    OUT_COORDS      =  5, &
    OUT_ACC         =  6, & 
    OUT_LASER       =  7, &
    OUT_ENERGY      =  8, &
    OUT_GAUGE_FIELD = 11, &
    OUT_TEMPERATURE = 12, &
    OUT_FTCHD       = 13, &
    OUT_VEL         = 14, &
    OUT_EIGS        = 15, &
    OUT_TOTAL_CURRENT = 17, &
    OUT_SEPARATE_COORDS  = 22, &
    OUT_SEPARATE_VELOCITY= 23, &
    OUT_SEPARATE_FORCES  = 24, &
    OUT_MAX              = 25
  

  type td_write_t
    private
    type(td_write_prop_t) :: out(OUT_MAX)

    integer        :: lmax     !< maximum multipole moment to output
    !> The states_type where the ground state is stored, in order to
    !! calculate the projections(s) onto it.
    type(states_t) :: gs_st    
    integer        :: n_excited_states  !< number of excited states onto which the projections are calculated.
    integer :: compute_interval     !< Compute every compute_interval
  end type td_write_t

contains


  ! ---------------------------------------------------------
  subroutine td_write_kick(mesh, kick, outp, geo, iter)
    type(mesh_t),     intent(in) :: mesh
    type(kick_t),     intent(in) :: kick
    type(output_t),   intent(in) :: outp
    type(geometry_t), intent(in) :: geo
    integer,          intent(in) :: iter

    character(len=256) :: filename
    PUSH_SUB(td_write_kick)

    write(filename, '(a,i7.7)') "td.", iter  ! name of directory
    call output_kick(outp, mesh, geo, kick, filename)

    POP_SUB(td_write_kick)
  end subroutine td_write_kick


  ! ---------------------------------------------------------
  subroutine td_write_init(writ, outp, gr, st, hm, geo, ks, ions_move, with_gauge_field, kick, iter, max_iter, dt, mc)
    type(td_write_t), target, intent(out)   :: writ
    type(output_t),           intent(out)   :: outp
    type(grid_t),             intent(in)    :: gr
    type(states_t),           intent(inout) :: st
    type(hamiltonian_t),      intent(inout) :: hm
    type(geometry_t),         intent(in)    :: geo
    type(v_ks_t),             intent(inout) :: ks
    logical,                  intent(in)    :: ions_move
    logical,                  intent(in)    :: with_gauge_field
    type(kick_t),             intent(in)    :: kick
    integer,                  intent(in)    :: iter
    integer,                  intent(in)    :: max_iter
    FLOAT,                    intent(in)    :: dt
    type(multicomm_t),        intent(in)    :: mc

    FLOAT :: rmin
    integer :: ierr, first, ii, ist, jj, kk, flags, iout, default
    type(block_t) :: blk
    character(len=100) :: filename

    PUSH_SUB(td_write_init)

    !%Variable TDOutput
    !%Type flag
    !%Default multipoles + energy (+ others depending on other options)
    !%Section Time-Dependent::TD Output
    !%Description
    !% Defines what should be output during the time-dependent
    !% simulation. Many of the options can increase the computational
    !% cost of the simulation, so only use the ones that you need. In
    !% most cases the default value is enough, as it is adapted to the
    !% details of the TD run. If the ions are allowed to be moved, additionally
    !% the geometry and the temperature are output. If a laser is
    !% included it will output by default.
    !%
    !% Note: the output files generated by this option are updated
    !% every <tt>RestartWriteInterval</tt> steps.
    !%Option multipoles bit(0)
    !% Outputs the (electric) multipole moments of the density to the file <tt>td.general/multipoles</tt>.
    !% This is required to, <i>e.g.</i>, calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>TDMultipoleLmax</tt>.
    !%Option geometry bit(4)
    !% If set (and if the atoms are allowed to move), outputs the coordinates, velocities,
    !% and forces of the atoms to the the file <tt>td.general/coordinates</tt>. On by default if <tt>MoveIons = yes</tt>.
    !%Option dipole_acceleration bit(5)
    !% When set, outputs the acceleration of the electronic dipole, calculated from the Ehrenfest theorem,
    !% in the file <tt>td.general/acceleration</tt>. This file can then be
    !% processed by the utility <tt>oct-harmonic-spectrum</tt> in order to obtain the harmonic spectrum.
    !%Option laser bit(6)
    !% If set, outputs the laser field to the file <tt>td.general/laser</tt>.
    !% On by default if <tt>TDExternalFields</tt> is set.
    !%Option energy bit(7)
    !% If set, <tt>octopus</tt> outputs the different components of the energy
    !% to the file <tt>td.general/energy</tt>. Will be zero except for every <tt>TDEnergyUpdateIter</tt> iterations.
    !%Option gauge_field bit(10)
    !% If set, outputs the vector gauge field corresponding to a spatially uniform (but time-dependent) 
    !% external electrical potential. This is only useful in a time-dependent periodic run.
    !% On by default if <tt>GaugeVectorField</tt> is set.
    !%Option temperature bit(11)
    !% If set, the ionic temperature at each step is printed. On by default if <tt>MoveIons = yes</tt>.
    !%Option ftchd bit(12)
    !% Write Fourier transform of the electron density to the file <tt>ftchd.X</tt>,
    !% where X depends on the kick (e.g. with sin-shaped perturbation X=sin).
    !% This is needed for calculating the dynamic structure factor.
    !% In the case that the kick mode is qbessel, the written quantity is integral over
    !% density, multiplied by spherical Bessel function times real spherical harmonic.
    !% On by default if <tt>TDMomentumTransfer</tt> is set.
    !%Option dipole_velocity bit(13)
    !% When set, outputs the dipole velocity, calculated from the Ehrenfest theorem,
    !% in the file <tt>td.general/velocity</tt>. This file can then be
    !% processed by the utility <tt>oct-harmonic-spectrum</tt> in order to obtain the harmonic spectrum.
    !%Option eigenvalues bit(14)
    !% Write the KS eigenvalues. 
    !%Option total_current bit(16)
    !% Output the total current (average of the current density over the cell).
    !%Option coordinates_sep bit(21)
    !% Writes geometries in a separate file.
    !%Option velocities_sep bit(22)
    !% Writes velocities in a separate file.
    !%Option forces_sep bit(23)
    !% Writes forces in a separate file.
    !%End

    default = 2**(OUT_MULTIPOLES - 1) +  2**(OUT_ENERGY - 1)
    if(ions_move) default = default + 2**(OUT_COORDS - 1) + 2**(OUT_TEMPERATURE - 1)
    if(with_gauge_field) default = default + 2**(OUT_GAUGE_FIELD - 1)
    if(hm%ep%no_lasers > 0) default = default + 2**(OUT_LASER - 1)
    if(kick%qkick_mode /= QKICKMODE_NONE) default = default + 2**(OUT_FTCHD - 1)

    call parse_variable('TDOutput', default, flags)

    if(.not.varinfo_valid_option('TDOutput', flags, is_flag = .true.)) call messages_input_error('TDOutput')

    do iout = 1, OUT_MAX
      writ%out(iout)%write = (bitand(flags, 2**(iout - 1)) /= 0)
    end do

    !%Variable TDMultipoleLmax
    !%Type integer
    !%Default 1
    !%Section Time-Dependent::TD Output
    !%Description
    !% Maximum electric multipole of the density output to the file <tt>td.general/multipoles</tt>
    !% during a time-dependent simulation. Must be non-negative.
    !%End
    call parse_variable('TDMultipoleLmax', 1, writ%lmax)
    if (writ%lmax < 0) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid TDMultipoleLmax."
      message(2) = '(Must be TDMultipoleLmax >= 0 )'
      call messages_fatal(2)
    end if
    call messages_obsolete_variable('TDDipoleLmax', 'TDMultipoleLmax')

    ! Compatibility test
    if( (writ%out(OUT_ACC)%write) .and. ions_move ) then
      message(1) = 'If harmonic spectrum is to be calculated, atoms should not be allowed to move.'
      call messages_fatal(1)
    end if

    rmin = geometry_min_distance(geo)
    if(geo%natoms == 1) then 
      if(simul_box_is_periodic(gr%sb)) then
        rmin = minval(gr%sb%lsize(1:gr%sb%periodic_dim))
      else
        rmin = CNST(100.0)
      end if
    end if

    !%Variable TDOutputComputeInterval
    !%Type integer
    !%Default 50
    !%Section Time-Dependent::TD Output
    !%Description
    !% The TD output requested are computed
    !% when the iteration number is a multiple of the <tt>TDOutputComputeInterval</tt> variable.
    !% Must be >= 0. If it is 0, then no output is written. 
    !% Implemented only for projections and number of excited electrons for the moment.
    !%End
    call parse_variable('TDOutputComputeInterval', 50, writ%compute_interval)
    if(writ%compute_interval < 0) then
      message(1) = "TDOutputComputeInterval must be >= 0."
      call messages_fatal(1)
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

      if(writ%out(OUT_COORDS)%write) &
        call write_iter_init(writ%out(OUT_COORDS)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/coordinates")))

      if(writ%out(OUT_SEPARATE_COORDS)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_COORDS)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/onlyCoordinates")))

      if(writ%out(OUT_SEPARATE_VELOCITY)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_VELOCITY)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/onlyVelocities")))

      if(writ%out(OUT_SEPARATE_FORCES)%write) &
        call write_iter_init(writ%out(OUT_SEPARATE_FORCES)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/onlyForces")))

      if(writ%out(OUT_TEMPERATURE)%write) &
        call write_iter_init(writ%out(OUT_TEMPERATURE)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/temperature")))

      if(writ%out(OUT_ACC)%write) &
        call write_iter_init(writ%out(OUT_ACC)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/acceleration")))
          
      if(writ%out(OUT_VEL)%write) &
        call write_iter_init(writ%out(OUT_VEL)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/velocity")))

      if(writ%out(OUT_LASER)%write) then
        if(iter .eq. 0) then
          call write_iter_init(writ%out(OUT_LASER)%handle, first, &
            units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/laser")))
          do ii = 0, max_iter
            call td_write_laser(writ%out(OUT_LASER)%handle, gr, hm, dt, ii)
          end do
          call write_iter_end(writ%out(OUT_LASER)%handle)
        end if
      end if

      if(writ%out(OUT_ENERGY)%write) &
        call write_iter_init(writ%out(OUT_ENERGY)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/energy")))
      
      if(writ%out(OUT_GAUGE_FIELD)%write) &
        call write_iter_init(writ%out(OUT_GAUGE_FIELD)%handle, &
        first, units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/gauge_field")))
      
      if(writ%out(OUT_EIGS)%write) &
        call write_iter_init(writ%out(OUT_EIGS)%handle, first, &
        units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/eigenvalues")))
      
      if(writ%out(OUT_TOTAL_CURRENT)%write) then
        call write_iter_init(writ%out(OUT_TOTAL_CURRENT)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("td.general/total_current")))
      end if

    end if
    
    if(writ%out(OUT_TOTAL_CURRENT)%write) then
      call v_ks_calculate_current(ks, .true.)
      call v_ks_calc(ks, hm, st, geo, calc_eigenval=.false., time = iter*dt)
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
        if(iout == OUT_LASER) cycle
        if(writ%out(iout)%write)  call write_iter_end(writ%out(iout)%handle)
      end do
     
    end if

    POP_SUB(td_write_end)
  end subroutine td_write_end


  ! ---------------------------------------------------------
  subroutine td_write_iter(writ, outp, gr, st, hm, geo, kick, dt,ks, iter)
    type(td_write_t),    intent(inout) :: writ !< Write object
    type(output_t),      intent(in)    :: outp
    type(grid_t),        intent(inout) :: gr   !< The grid
    type(states_t),      intent(inout) :: st   !< State object
    type(hamiltonian_t), intent(inout) :: hm   !< Hamiltonian object
    type(geometry_t),    intent(inout) :: geo  !< Geometry object
    type(kick_t),        intent(in)    :: kick !< The kick
    FLOAT,               intent(in)    :: dt   !< Delta T, time interval
    type(v_ks_t),        intent(in)    :: ks  
    integer,             intent(in)    :: iter !< Iteration number
    type(profile_t), save :: prof

    PUSH_SUB(td_write_iter)
    call profiling_in(prof, "TD_WRITE_ITER")

    if(writ%out(OUT_MULTIPOLES)%write) &
      call td_write_multipole(writ%out(OUT_MULTIPOLES)%handle, gr, geo, st, writ%lmax, kick, iter)
    
    if(writ%out(OUT_FTCHD)%write) &
      call td_write_ftchd(writ%out(OUT_FTCHD)%handle, gr, st, kick, iter)

    if(writ%out(OUT_COORDS)%write) &
      call td_write_coordinates(writ%out(OUT_COORDS)%handle, gr, geo, iter)

    if(writ%out(OUT_SEPARATE_COORDS)%write) &
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_COORDS)%handle, gr, geo, iter,1)

    if(writ%out(OUT_SEPARATE_VELOCITY)%write) &
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_VELOCITY)%handle, gr, geo, iter,2)

    if(writ%out(OUT_SEPARATE_FORCES)%write) &
      call td_write_sep_coordinates(writ%out(OUT_SEPARATE_FORCES)%handle, gr, geo, iter,3)

    if(writ%out(OUT_TEMPERATURE)%write) &
      call td_write_temperature(writ%out(OUT_TEMPERATURE)%handle, geo, iter)

    if(writ%out(OUT_ACC)%write) &
      call td_write_acc(writ%out(OUT_ACC)%handle, gr, geo, st, hm, dt, iter)
      
    if(writ%out(OUT_VEL)%write) &
      call td_write_vel(writ%out(OUT_VEL)%handle, gr, st, iter)

    ! td_write_laser no longer called here, because the whole laser is printed
    ! out at the beginning.

    if(writ%out(OUT_ENERGY)%write) &
      call td_write_energy(writ%out(OUT_ENERGY)%handle, hm, iter, geo%kinetic_energy)

    if(writ%out(OUT_GAUGE_FIELD)%write) &
      call td_write_gauge_field(writ%out(OUT_GAUGE_FIELD)%handle, hm, gr, iter)

    if(writ%out(OUT_EIGS)%write) &
      call td_write_eigs(writ%out(OUT_EIGS)%handle, st, iter)

    if(writ%out(OUT_TOTAL_CURRENT)%write) then
      call td_write_total_current(writ%out(OUT_TOTAL_CURRENT)%handle, gr, st, iter)
    end if
    
    call profiling_out(prof)
    POP_SUB(td_write_iter)
  end subroutine td_write_iter


  ! ---------------------------------------------------------
  subroutine td_write_data(writ, gr, st, hm, ks, outp, geo, iter, dt)
    type(td_write_t),     intent(inout) :: writ
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    type(v_ks_t),         intent(in)    :: ks
    type(output_t),       intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: iter
    FLOAT, optional,      intent(in)    :: dt

    integer :: iout
    type(profile_t), save :: prof

    PUSH_SUB(td_write_data)
    call profiling_in(prof, "TD_WRITE_DATA")

    if(mpi_grp_is_root(mpi_world)) then
      do iout = 1, OUT_MAX
        if(iout == OUT_LASER) cycle
        if(writ%out(iout)%write)  call write_iter_flush(writ%out(iout)%handle)
      end do
    end if

    call profiling_out(prof)
    POP_SUB(td_write_data)
  end subroutine td_write_data

  ! ---------------------------------------------------------
  subroutine td_write_output(writ, gr, st, hm, ks, outp, geo, iter, dt)
    type(td_write_t),     intent(inout) :: writ
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    type(v_ks_t),         intent(in)    :: ks
    type(output_t),       intent(in)    :: outp
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: iter
    FLOAT, optional,      intent(in)    :: dt
    
    character(len=256) :: filename
    type(profile_t), save :: prof

    PUSH_SUB(td_write_output)
    call profiling_in(prof, "TD_WRITE_DATA")

    ! now write down the rest
    write(filename, '(a,a,i7.7)') trim(outp%iter_dir),"td.", iter  ! name of directory

    ! this is required if st%X(psi) is used
    call states_sync(st)

    call output_all(outp, gr, geo, st, hm, ks, filename)
    if(present(dt)) then
      call output_scalar_pot(outp, gr, geo, hm, filename, iter*dt)
    else
      if(iter == 0) call output_scalar_pot(outp, gr, geo, hm, filename)
    end if
 
    call profiling_out(prof)
    POP_SUB(td_write_output)
  end subroutine td_write_output

  ! ---------------------------------------------------------
  !> Subroutine to write multipoles to the corresponding file
  subroutine td_write_multipole(out_multip, gr, geo, st, lmax, kick, iter)
    type(c_ptr),        intent(inout) :: out_multip !< C pointer
    type(grid_t),       intent(in) :: gr   !< The grid
    type(geometry_t),   intent(in) :: geo  !< Geometry object
    type(states_t),     intent(in) :: st   !< State object
    integer,            intent(in) :: lmax
    type(kick_t),       intent(in) :: kick !< Kick object
    integer,            intent(in) :: iter !< Iteration number

    integer :: is, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT, allocatable :: ionic_dipole(:), multipole(:,:)

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

    SAFE_ALLOCATE(ionic_dipole(1:gr%mesh%sb%dim))
    SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin))
    ionic_dipole(:) = M_ZERO
    multipole   (:,:) = M_ZERO

    do is = 1, st%d%nspin
      call dmf_multipoles(gr%mesh, st%rho(:,is), lmax, multipole(:,is))
    end do

    if (lmax > 0) then
      call geometry_dipole(geo, ionic_dipole)
      do is = 1, st%d%nspin
        multipole(2:gr%mesh%sb%dim+1, is) = -ionic_dipole(1:gr%mesh%sb%dim)/st%d%nspin - multipole(2:gr%mesh%sb%dim+1, is)
      end do
    end if

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

    SAFE_DEALLOCATE_A(ionic_dipole)
    SAFE_DEALLOCATE_A(multipole)
    POP_SUB(td_write_multipole)
  end subroutine td_write_multipole

  ! ---------------------------------------------------------
  subroutine td_write_ftchd(out_ftchd, gr, st, kick, iter)
    type(c_ptr),        intent(inout) :: out_ftchd
    type(grid_t),       intent(in) :: gr
    type(states_t),     intent(in) :: st
    type(kick_t),       intent(in) :: kick
    integer,            intent(in) :: iter

    integer :: is, ip, idir
    character(len=120) :: aux, aux2
    FLOAT   :: ftchd_bessel
    CMPLX   :: ftchd
    FLOAT   :: ylm
    FLOAT, allocatable :: integrand_bessel(:)
    CMPLX, allocatable :: integrand(:)

    PUSH_SUB(td_write_ftchd)

    if(mpi_grp_is_root(mpi_world).and.iter == 0) then
      call td_write_print_header_init(out_ftchd)

      write(aux,'(a15, i2)') '# qkickmode    ', kick%qkick_mode
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        write(aux,'(a15, i0.3, 1x, i0.3)') '# ll, mm       ', kick%qbessel_l, kick%qbessel_m
        call write_iter_string(out_ftchd, aux)
        call write_iter_nl(out_ftchd)
      end if

      if(kick%qkick_mode == QKICKMODE_BESSEL) then
        write(aux, '(a15, f9.6)') '# qlength      ', kick%qlength
      else ! sin or cos
        write(aux, '(a15)')       '# qvector      '
        do idir = 1, gr%mesh%sb%dim
          write(aux2, '(f9.5)') kick%qvector(idir)
          aux = trim(aux) // trim(aux2)
        end do
      end if
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
      call write_iter_string(out_ftchd, aux)
      call write_iter_nl(out_ftchd)

      call write_iter_header_start(out_ftchd)
      if(kick%qkick_mode == QKICKMODE_BESSEL) then
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
    if(kick%qkick_mode /= QKICKMODE_BESSEL) then
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
          call grylmr(gr%mesh%x(ip, 1), gr%mesh%x(ip, 2), gr%mesh%x(ip, 3), kick%qbessel_l, kick%qbessel_m, ylm)
          integrand_bessel(ip) = integrand_bessel(ip) + st%rho(ip, is) * &
                                 loct_sph_bessel(kick%qbessel_l, kick%qlength*sqrt(sum(gr%mesh%x(ip, :)**2)))*ylm
        end do
      end do
      ftchd_bessel = dmf_integrate(gr%mesh, integrand_bessel)
      SAFE_DEALLOCATE_A(integrand_bessel)
    end if

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_ftchd)
      if(kick%qkick_mode == QKICKMODE_BESSEL) then
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
    type(c_ptr),       intent(inout) :: out_coords
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
  subroutine td_write_sep_coordinates(out_coords, gr, geo, iter,which)
    type(c_ptr),       intent(inout) :: out_coords
    type(grid_t),      intent(in) :: gr
    type(geometry_t),  intent(in) :: geo
    integer,           intent(in) :: iter
    integer,           intent(in) :: which !1=xyz, 2=velocity, 3=force

    integer, parameter :: COORDINATES=1
    integer, parameter :: VELOCITIES=2
    integer, parameter :: FORCES=3
    integer :: iatom, idir
    character(len=50) :: aux
    FLOAT :: tmp(1:MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_sep_coordinates)

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)

      do iatom = 1, geo%natoms
        do idir = 1, gr%mesh%sb%dim
          select case (which)
            case (COORDINATES)
              write(aux, '(a2,i3,a1,i3,a1)') 'x(', iatom, ',', idir, ')'
            case (VELOCITIES)
              write(aux, '(a2,i3,a1,i3,a1)') 'v(', iatom, ',', idir,')'
            case (FORCES)
              write(aux, '(a2,i3,a1,i3,a1)') 'f(', iatom, ',', idir,')'
          end select
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_abbrev(units_out%time)) // ']')
      select case (which)
        case (COORDINATES)
          call write_iter_string(out_coords, &
            'Positions in '   // trim(units_abbrev(units_out%length))) 
        case (VELOCITIES)
          call write_iter_string(out_coords, &
            'Velocities in '  // trim(units_abbrev(units_out%velocity)))
        case (FORCES)
          call write_iter_string(out_coords, &
            'Forces in '    // trim(units_abbrev(units_out%force)))
      end select
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)

    select case (which)
      case (COORDINATES)
        do iatom = 1, geo%natoms
          tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%length, geo%atom(iatom)%x(1:gr%mesh%sb%dim))
          call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
        end do
      case (VELOCITIES)
        do iatom = 1, geo%natoms
           tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%velocity, geo%atom(iatom)%v(1:gr%mesh%sb%dim))
           call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
        end do
      case (FORCES)
        do iatom = 1, geo%natoms
           tmp(1:gr%mesh%sb%dim) = units_from_atomic(units_out%force, geo%atom(iatom)%f(1:gr%mesh%sb%dim))
           call write_iter_double(out_coords, tmp, gr%mesh%sb%dim)
        end do
    end select
       
    call write_iter_nl(out_coords)

    POP_SUB(td_write_sep_coordinates)
  end subroutine td_write_sep_coordinates


  ! ---------------------------------------------------------
  subroutine td_write_temperature(out_temperature, geo, iter)
    type(c_ptr),       intent(inout) :: out_temperature
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

    call write_iter_double(out_temperature, units_from_atomic(unit_kelvin, ion_dynamics_temperature(geo)), 1)

    call write_iter_nl(out_temperature)

    POP_SUB(td_write_temperature)
  end subroutine td_write_temperature


  ! ---------------------------------------------------------
  subroutine td_write_acc(out_acc, gr, geo, st, hm, dt, iter)
    type(c_ptr),         intent(inout) :: out_acc
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

    ! this is required if st%X(psi) is used
    call states_sync(st)

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
  subroutine td_write_vel(out_vel, gr, st, iter)
    type(c_ptr),         intent(inout) :: out_vel
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
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

    call td_calc_tvel(gr, st, vel)

    call write_iter_start(out_vel)
    vel = units_from_atomic(units_out%velocity, vel)
    call write_iter_double(out_vel, vel, gr%mesh%sb%dim)
    call write_iter_nl(out_vel)

    POP_SUB(td_write_vel)
  end subroutine td_write_vel


  ! ---------------------------------------------------------
  subroutine td_write_laser(out_laser, gr, hm, dt, iter)
    type(c_ptr),         intent(inout) :: out_laser
    type(grid_t),        intent(in) :: gr
    type(hamiltonian_t), intent(in) :: hm
    FLOAT,               intent(in) :: dt
    integer,             intent(in) :: iter

    integer :: il, idir
    FLOAT :: field(MAX_DIM)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    ! no PUSH SUB, called too often

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
        case(E_FIELD_SCALAR_POTENTIAL)
          write(aux, '(a,i1,a)') 'e(t)'
          call write_iter_header(out_laser, aux)
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
        case(E_FIELD_SCALAR_POTENTIAL)
          aux = '[adim]'
          call write_iter_header(out_laser, aux)
        end select
      end do
      call write_iter_nl(out_laser)

      call td_write_print_header_end(out_laser)
    end if

    call write_iter_start(out_laser)

    do il = 1, hm%ep%no_lasers
      field = M_ZERO
      call laser_field(hm%ep%lasers(il), field(1:gr%sb%dim), iter*dt)
      select case(laser_kind(hm%ep%lasers(il)))
      case(E_FIELD_ELECTRIC, E_FIELD_MAGNETIC)
        field = units_from_atomic(units_out%force, field)
        call write_iter_double(out_laser, field, gr%mesh%sb%dim)
      case(E_FIELD_VECTOR_POTENTIAL)
        field = units_from_atomic(units_out%energy, field)
        call write_iter_double(out_laser, field, gr%mesh%sb%dim)
      case(E_FIELD_SCALAR_POTENTIAL)
        call write_iter_double(out_laser, field(1), 1)
      end select
    end do

    call write_iter_nl(out_laser)

  end subroutine td_write_laser


  ! ---------------------------------------------------------
  subroutine td_write_energy(out_energy, hm, iter, ke)
    type(c_ptr),         intent(inout) :: out_energy
    type(hamiltonian_t), intent(in) :: hm
    integer,             intent(in) :: iter
    FLOAT,               intent(in) :: ke

    integer :: ii

    integer :: n_columns         

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(td_write_energy)

    n_columns = 9          

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

      ! units

      call write_iter_string(out_energy, '#[Iter n.]')
      call write_iter_header(out_energy, '[' // trim(units_abbrev(units_out%time)) // ']')

      do ii = 1, n_columns
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
  subroutine td_write_eigs(out_eigs, st, iter)
    type(c_ptr),         intent(inout) :: out_eigs
    type(states_t),      intent(in) :: st
    integer,             intent(in) :: iter

    integer             :: ii, is
    character(len=68)   :: buf
    FLOAT, allocatable  :: eigs(:,:)
#if defined(HAVE_MPI) 
    integer :: outcount, ik
#endif

    PUSH_SUB(td_write_eigs)

    SAFE_ALLOCATE(eigs(1:st%nst,1:st%d%kpt%nglobal)) 
           
    eigs(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
      st%eigenval(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end)

#if defined(HAVE_MPI) 

    do ik = st%d%kpt%start, st%d%kpt%end
      call lmpi_gen_allgatherv(st%lnst, st%eigenval(st%st_start:st%st_end,ik), outcount, &
                             eigs(:, ik), st%mpi_grp)
    end do

#endif
  
    if(.not.mpi_grp_is_root(mpi_world)) then 
      SAFE_DEALLOCATE_A(eigs)
      POP_SUB(td_write_eigs)
      return ! only first node outputs        
    end if


    if(iter == 0) then
      call td_write_print_header_init(out_eigs)

      write(buf, '(a15,i2)')      '# nst          ', st%nst
      call write_iter_string(out_eigs, buf)
      call write_iter_nl(out_eigs)

      write(buf, '(a15,i2)')      '# nspin        ', st%d%nspin
      call write_iter_string(out_eigs, buf)
      call write_iter_nl(out_eigs)

      ! first line -> column names
      call write_iter_header_start(out_eigs)
      do is = 1, st%d%kpt%nglobal
        do ii = 1, st%nst
          write(buf, '(a,i4)') 'Eigenvalue ',ii
          call write_iter_header(out_eigs, buf)
        end do
      end do
      call write_iter_nl(out_eigs)

      ! second line: units
      call write_iter_string(out_eigs, '#[Iter n.]')
      call write_iter_header(out_eigs, '[' // trim(units_abbrev(units_out%time)) // ']')
      do is = 1, st%d%kpt%nglobal
        do ii = 1, st%nst
          call write_iter_header(out_eigs, '[' // trim(units_abbrev(units_out%energy)) // ']')
        end do
      end do
      call write_iter_nl(out_eigs)
      call td_write_print_header_end(out_eigs)
    end if

    call write_iter_start(out_eigs)
    do is = 1, st%d%kpt%nglobal
      do ii =1 , st%nst
        call write_iter_double(out_eigs, units_from_atomic(units_out%energy, eigs(ii,is)), 1)
      end do
    end do
    call write_iter_nl(out_eigs)

    SAFE_DEALLOCATE_A(eigs)

    POP_SUB(td_write_eigs)
  end subroutine td_write_eigs

  ! ---------------------------------------------------------
  subroutine td_write_gauge_field(out_gauge, hm, gr, iter)
    type(c_ptr),         intent(inout) :: out_gauge
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
    call gauge_field_get_vec_pot(hm%ep%gfield, temp)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    call gauge_field_get_vec_pot_vel(hm%ep%gfield, temp)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy / units_out%time, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    call gauge_field_get_vec_pot_acc(hm%ep%gfield, temp)
    forall(idir = 1:gr%mesh%sb%dim)
      temp(idir) = units_from_atomic(units_out%energy / units_out%time**2, temp(idir))
    end forall
    call write_iter_double(out_gauge, temp, gr%mesh%sb%dim)

    call write_iter_nl(out_gauge)
    POP_SUB(td_write_gauge_field)
    
  end subroutine td_write_gauge_field

  ! ---------------------------------------------------------
  subroutine td_write_total_current(out_total_current, gr, st, iter)
    type(c_ptr),         intent(inout) :: out_total_current
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: iter

    integer :: idir, ispin
    character(len=50) :: aux
    FLOAT :: total_current(1:MAX_DIM), abs_current(1:MAX_DIM)

    PUSH_SUB(td_write_total_current)

    if(mpi_grp_is_root(mpi_world) .and. iter == 0) then
      call td_write_print_header_init(out_total_current)

      ! first line: column names
      call write_iter_header_start(out_total_current)
      
      do idir = 1, gr%mesh%sb%dim
        write(aux, '(a2,i1,a1)') 'I(', idir, ')'
        call write_iter_header(out_total_current, aux)
      end do

      do idir = 1, gr%mesh%sb%dim
        write(aux, '(a2,i1,a1)') 'IntAbs(j)(', idir, ')'
        call write_iter_header(out_total_current, aux)
      end do
      
      do ispin = 1, st%d%nspin
        do idir = 1, gr%mesh%sb%dim
          write(aux, '(a4,i1,a1,i1,a1)') 'I-sp', ispin, '(', idir, ')'
          call write_iter_header(out_total_current, aux)
        end do
      end do      

      call write_iter_nl(out_total_current)

      call td_write_print_header_end(out_total_current)
    end if
    
    ASSERT(associated(st%current))

    if(mpi_grp_is_root(mpi_world)) &
      call write_iter_start(out_total_current)

    total_current = CNST(0.0)
    do idir = 1, gr%sb%dim
      do ispin = 1, st%d%spin_channels
        total_current(idir) =  total_current(idir) + dmf_integrate(gr%mesh, st%current(:, idir, ispin))
      end do
      total_current(idir) = units_from_atomic(units_out%length/units_out%time, total_current(idir))
    end do

    abs_current = CNST(0.0)
    do idir = 1, gr%sb%dim
      do ispin = 1, st%d%spin_channels
        abs_current(idir) =  abs_current(idir) + dmf_integrate(gr%mesh, abs(st%current(:, idir, ispin)))
      end do
      abs_current(idir) = units_from_atomic(units_out%length/units_out%time, abs_current(idir))
    end do

   if(mpi_grp_is_root(mpi_world)) then
      call write_iter_double(out_total_current, total_current, gr%mesh%sb%dim)
      call write_iter_double(out_total_current, abs_current, gr%mesh%sb%dim)
   end if
  
    do ispin = 1, st%d%nspin
      total_current = CNST(0.0)
      do idir = 1, gr%sb%dim
        total_current(idir) = units_from_atomic(units_out%length/units_out%time, &
                                    dmf_integrate(gr%mesh, st%current(:, idir, ispin)))
      end do
      if(mpi_grp_is_root(mpi_world)) &
        call write_iter_double(out_total_current, total_current, gr%mesh%sb%dim)
    end do

    if(mpi_grp_is_root(mpi_world)) &
      call write_iter_nl(out_total_current)
      
    POP_SUB(td_write_total_current)
  end subroutine td_write_total_current

  ! ---------------------------------------------------------
  

  subroutine td_write_print_header_init(out)
    type(c_ptr), intent(inout) :: out

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
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(td_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(td_write_print_header_end)
  end subroutine td_write_print_header_end

end module td_write_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
