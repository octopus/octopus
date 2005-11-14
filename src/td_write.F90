!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module td_write
  use global
  use syslabels
  use units
  use messages
  use io
  use lib_oct
  use geometry
  use mesh
  use grid
  use mesh_function
  use functions
  use states
  use output
  use restart
  use lasers
  use hamiltonian
  use external_pot
  use grid
  use spectrum
  use mpi_mod
  use varinfo

  implicit none

  private
  public ::         &
    td_write_type,  &
    td_write_init,  &
    td_write_end,   &
    td_write_iter,  &
    td_write_data


  type td_write_type
    integer(POINTER_SIZE) :: &
      out_multip,   &
      out_coords,   &
      out_gsp,      &
      out_acc,      &
      out_laser,    &
      out_energy,   &
      out_proj,     &
      out_angular,  &
      out_spin,     &
      out_magnets

    integer           :: lmax     ! maximum multipole moment to output
    FLOAT             :: lmm_r    ! radius of the sphere used to compute the local magnetic moments
    type(states_type) :: gs_st    ! The states_type where the ground state is stored, in order to
                                  ! calculate the projections(s) onto it.
  end type td_write_type

contains

  ! ---------------------------------------------------------
  subroutine td_write_init(w, gr, st, geo, ions_move, there_are_lasers, iter, dt)
    type(td_write_type)             :: w
    type(grid_type),     intent(in) :: gr
    type(states_type),   intent(in) :: st
    type(geometry_type), intent(in) :: geo
    logical,             intent(in) :: ions_move, there_are_lasers
    integer,             intent(in) :: iter
    FLOAT,               intent(in) :: dt


    FLOAT :: rmin
    integer :: ierr, first, i, flags

    call push_sub('td_write.td_write_handler')


    !%Variable TDOutput
    !%Type flag
    !%Default multipoles + geometry
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
    !%Option gs_proj 8
    !% Outputs the projection of the time-dependent Kohn-Sham Slater determinant
    !% onto the ground-state to the file <tt>td.general/gs_projection</tt>. As the calculation
    !% of the projection is fairly heavy, this is only done every <tt>OutputEvery</tt> 
    !% iterations.
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
    !%Option el_energy 128
    !% If <tt>set</tt>, <tt>octopus</tt> outputs the different components of the electronic energy
    !% to the file <tt>td.general/el_energy</tt>.
    !%Option td_occup 256
    !% If set, outputs the projections of the time-dependent Kohn-Sham
    !% wave-functions onto the static (zero time) wave-functions to the
    !% file <tt>td.general/projections.XXX</tt>.
    !%Option local_mag_moments 512
    !% If set, outputs the local magnetic moments, integrated in sphere centered around each atom.
    !% The radius of the sphere can be ser with <tt>LocalMagneticMomentsSphereRadius</tt>
    !%End
    call loct_parse_int(check_inp('TDOutput'), 1+16, flags)
    if(.not.varinfo_valid_option('TDOutput', flags, is_flag=.true.)) then
      call input_error('TDOutput')
    end if

    w%out_multip  = 0; if(iand(flags,   1).ne.0) w%out_multip  = 1
    w%out_angular = 0; if(iand(flags,   2).ne.0) w%out_angular = 1
    w%out_spin    = 0; if(iand(flags,   4).ne.0) w%out_spin    = 1
    w%out_gsp     = 0; if(iand(flags,   8).ne.0) w%out_gsp     = 1
    w%out_coords  = 0; if(iand(flags,  16).ne.0.and.ions_move) w%out_coords = 1
    w%out_acc     = 0; if(iand(flags,  32).ne.0) w%out_acc     = 1
    w%out_laser   = 0; if(iand(flags,  64).ne.0) w%out_laser   = 1
    w%out_energy  = 0; if(iand(flags, 128).ne.0) w%out_energy  = 1
    w%out_proj    = 0; if(iand(flags, 256).ne.0) w%out_proj    = 1
    w%out_magnets = 0; if(iand(flags, 512).ne.0) w%out_magnets = 1

    !%Variable TDDipoleLmax
    !%Type integer
    !%Default 1
    !%Section Time Dependent::TD Output
    !%Description
    !% Maximum multi-pole of the density output to the file @code{td.general/multipoles} 
    !% during a time-dependent simulation. Must be 0 &lt; <tt>TDDipoleLmax &lt; 5</tt>.
    !%End
    call loct_parse_int(check_inp('TDDipoleLmax'), 1, w%lmax)
    if (w%lmax < 0 .or. w%lmax > 4) then
      write(message(1), '(a,i6,a)') "Input: '", w%lmax, "' is not a valid TDDipoleLmax"
      message(2) = '(0 <= TDDipoleLmax <= 4 )'
      call write_fatal(2)
    end if

    ! Compatibility test
    if( (w%out_acc.ne.0) .and. ions_move ) then
      message(1) = 'Error. If harmonic spectrum is to be calculated'
      message(2) = 'Atoms should not be allowed to move'
      call write_fatal(2)
    end if

    call geometry_min_distance(geo, rmin)
    call loct_parse_float(check_inp('LocalMagneticMomentsSphereRadius'), rmin*M_HALF/units_inp%length%factor, w%lmm_r)
    w%lmm_r = w%lmm_r * units_inp%length%factor

    if( (w%out_proj.ne.0)  .or.  (w%out_gsp.ne.0) ) then
      call states_copy(w%gs_st, st)
      ! WARNING: should be first deallocate, then nullify?
      nullify(w%gs_st%zpsi, w%gs_st%node, w%gs_st%occ, w%gs_st%eigenval, w%gs_st%mag)
      call restart_look (trim(tmpdir)//'restart_gs', gr%m, i, i, w%gs_st%nst, ierr)
      ! We will store the ground-state Kohn-Sham system by all processors.
      w%gs_st%st_start = 1
      w%gs_st%st_end   = w%gs_st%nst
      ! allocate memory
      allocate(w%gs_st%node(w%gs_st%nst))
      allocate(w%gs_st%eigenval(w%gs_st%nst, w%gs_st%d%nik))
      allocate(w%gs_st%occ(w%gs_st%nst, w%gs_st%d%nik))
      if(w%gs_st%d%ispin == SPINORS) then
        allocate(w%gs_st%mag(w%gs_st%nst, w%gs_st%d%nik, 2))
      end if
      allocate(w%gs_st%zpsi(NP, w%gs_st%d%dim, w%gs_st%st_start:w%gs_st%st_end, w%gs_st%d%nik))
      w%gs_st%node(:)  = 0
      call zrestart_read(trim(tmpdir)//'restart_gs', w%gs_st, gr%m, ierr)
      if(ierr.ne.0) then
        message(1) = "Could not load "//trim(tmpdir)//"restart_gs"
        call write_fatal(1)
      end if
    end if

    if (iter == 0) then
      first = 0
    else
      first = iter + 1
    end if

    if(mpi_grp_is_root(mpi_world)) then
      if(w%out_multip.ne.0)  call write_iter_init(w%out_multip,  first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/multipoles")))
      if(w%out_angular.ne.0) call write_iter_init(w%out_angular, first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/angular")))
      if(w%out_spin.ne.0)    call write_iter_init(w%out_spin,    first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/spin")))
      if(w%out_magnets.ne.0) call write_iter_init(w%out_magnets, first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/magnetic_moments")))
      if(w%out_coords.ne.0)  call write_iter_init(w%out_coords,  first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/coordinates")))
      if(w%out_gsp.ne.0)     call write_iter_init(w%out_gsp,     first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/gs_projection")))
      if(w%out_acc.ne.0)     call write_iter_init(w%out_acc,     first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/acceleration")))
      if(w%out_laser.ne.0)   call write_iter_init(w%out_laser,   first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/laser")))
      if(w%out_energy.ne.0)  call write_iter_init(w%out_energy,  first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/el_energy")))
      if(w%out_proj.ne.0)    call write_iter_init(w%out_proj,    first, dt/units_out%time%factor, &
        trim(io_workpath("td.general/projections")))
    end if

    call pop_sub()
  end subroutine td_write_init


  ! ---------------------------------------------------------
  subroutine td_write_end(w)
    type(td_write_type)       :: w
    call push_sub('td_write.td_write_end')

    if(mpi_grp_is_root(mpi_world)) then
      if(w%out_multip.ne.0)  call write_iter_end(w%out_multip)
      if(w%out_angular.ne.0) call write_iter_end(w%out_angular)
      if(w%out_spin.ne.0)    call write_iter_end(w%out_spin)
      if(w%out_magnets.ne.0) call write_iter_end(w%out_magnets)
      if(w%out_coords.ne.0)  call write_iter_end(w%out_coords)
      if(w%out_gsp.ne.0)     call write_iter_end(w%out_gsp)
      if(w%out_acc.ne.0)     call write_iter_end(w%out_acc)
      if(w%out_laser.ne.0)   call write_iter_end(w%out_laser)
      if(w%out_energy.ne.0)  call write_iter_end(w%out_energy)
      if(w%out_proj.ne.0)    call write_iter_end(w%out_proj)
    end if
    if( (w%out_gsp.ne.0) .or. (w%out_proj.ne.0) ) then
      call states_end(w%gs_st)
    end if

    call pop_sub()
  end subroutine td_write_end


  ! ---------------------------------------------------------
  subroutine td_write_iter(w, gr, st, h, geo, kick, dt, i)
    type(td_write_type),    intent(in) :: w
    type(grid_type),        intent(inout) :: gr
    type(states_type),      intent(inout) :: st
    type(hamiltonian_type), intent(inout) :: h
    type(geometry_type),    intent(in)    :: geo
    type(kick_type),        intent(in) :: kick
    FLOAT,                  intent(in) :: dt
    integer,                intent(in) :: i

    call push_sub('td_write.td_write_iter')

    if(w%out_multip.ne.0)   call td_write_multipole(w%out_multip, gr, st, w%lmax, kick, i)
    if(w%out_angular.ne.0)  call td_write_angular(w%out_angular, gr, st, kick, i)
    if(w%out_spin.ne.0)     call td_write_spin(w%out_spin, gr%m, st, i)
    if(w%out_magnets.ne.0)  call td_write_local_magnetic_moments(w%out_magnets, gr%m, st, geo, w%lmm_r, i)
    if(w%out_proj.ne.0)     call td_write_proj(w%out_proj, gr, st, w%gs_st, i)
    if(w%out_coords.ne.0)   call td_write_nbo(w%out_coords, gr, i, geo%kinetic_energy, h%etot)
    if(w%out_gsp.ne.0)      call td_write_gsp(w%out_gsp, gr%m, st, w%gs_st, dt, i)
    if(w%out_acc.ne.0)      call td_write_acc(w%out_acc, gr, st, h, dt, i)
    if(w%out_laser.ne.0)    call td_write_laser(w%out_laser, gr, h, dt, i)
    if(w%out_energy.ne.0)   call td_write_el_energy(w%out_energy, h, i)

    call pop_sub()
  end subroutine td_write_iter


  ! ---------------------------------------------------------
  subroutine td_write_data(w, gr, st, h, outp, geo, dt_, iter)
    type(td_write_type),    intent(in)    :: w
    type(grid_type),        intent(inout) :: gr
    type(states_type),      intent(in)    :: st
    type(hamiltonian_type), intent(in)    :: h
    type(output_type),      intent(in)    :: outp
    type(geometry_type),    intent(in)    :: geo
    FLOAT, intent(in) :: dt_
    integer, intent(in) :: iter

    FLOAT :: dt
    character(len=256) :: filename

    call push_sub('td.td_write_data')

    if(mpi_grp_is_root(mpi_world)) then
      if(w%out_multip.ne.0)  call write_iter_flush(w%out_multip)
      if(w%out_angular.ne.0) call write_iter_flush(w%out_angular)
      if(w%out_spin.ne.0)    call write_iter_flush(w%out_spin)
      if(w%out_magnets.ne.0) call write_iter_flush(w%out_magnets)
      if(w%out_coords.ne.0)  call write_iter_flush(w%out_coords)
      if(w%out_gsp.ne.0)     call write_iter_flush(w%out_gsp)
      if(w%out_acc.ne.0)     call write_iter_flush(w%out_acc)
      if(w%out_laser.ne.0)   call write_iter_flush(w%out_laser)
      if(w%out_energy.ne.0)  call write_iter_flush(w%out_energy)
      if(w%out_proj.ne.0)    call write_iter_flush(w%out_proj)
    end if

    ! now write down the rest
    write(filename, '(a,i7.7)') "td.", iter  ! name of directory

    call zstates_output(st, gr, filename, outp)
    if(iand(outp%what, output_geometry).ne.0) &
      call atom_write_xyz(filename, "geometry", geo)
    call hamiltonian_output(h, gr%m, gr%sb, filename, outp)

    dt = dt_
!!$#if !defined(DISABLE_PES) && defined(HAVE_FFT)
!!$  call PES_output(td%PESv, gr%m, st, iter, outp%iter, dt)
!!$#endif

    call pop_sub()
  end subroutine td_write_data


  ! ---------------------------------------------------------
  subroutine td_write_spin(out_spin, m, st, iter)
    integer(POINTER_SIZE), intent(in) :: out_spin
    type(mesh_type),       intent(in) :: m
    type(states_type),     intent(in) :: st
    integer,               intent(in) :: iter

    character(len=130) :: aux
    FLOAT :: spin(3)

    call push_sub('td_write.td_write_spin')

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      ! The expectation value of the spin operator is half the total magnetic moment
      call states_magnetic_moment(m, st, st%rho, spin)
      spin = M_HALF*spin

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
  subroutine td_write_local_magnetic_moments(out_magnets, m, st, geo, lmm_r, iter)
    integer(POINTER_SIZE), intent(in) :: out_magnets
    type(mesh_type),       intent(in) :: m
    type(states_type),     intent(in) :: st
    type(geometry_type),   intent(in) :: geo
    FLOAT,                 intent(in) :: lmm_r
    integer,               intent(in) :: iter

    integer :: ia
    character(len=50) :: aux
    FLOAT, allocatable :: lmm(:,:)

    call push_sub('td_write.td_write_local_magnetic_moments')

    if(mpi_grp_is_root(mpi_world)) then ! only first node outputs

      !get the atoms magnetization
      allocate(lmm(3, geo%natoms))
      call states_local_magnetic_moments(m, st, geo, st%rho, lmm_r, lmm)

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
      deallocate(lmm)
    end if

    call pop_sub()
  end subroutine td_write_local_magnetic_moments


  ! ---------------------------------------------------------
  subroutine td_write_angular(out_angular, gr, st, kick, iter)
    integer(POINTER_SIZE), intent(in)    :: out_angular
    type(grid_type),       intent(inout) :: gr
    type(states_type),     intent(inout) :: st
    type(kick_type),       intent(in)    :: kick
    integer,               intent(in)    :: iter

    character(len=130) :: aux
    FLOAT :: angular(3)

    call push_sub('td_write.td_write_angular')

    ! The angular momentum has to be calculated by all nodes...
    call zstates_calc_angular(gr, st, angular)

    if(mpi_grp_is_root(mpi_world)) then ! Only first node outputs

      if(iter ==0) then
        call td_write_print_header_init(out_angular)

        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,i2)')      '# direction    ', kick%pol_dir
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,i2)')      '# kick mode    ', kick%delta_strength_mode
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,f18.12)')  '# kick strength', kick%delta_strength
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,i2)')      '# Equiv. axis  ', kick%pol_equiv_axis
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        write(aux, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
        call write_iter_string(out_angular, aux)
        call write_iter_nl(out_angular)

        !second line -> columns name
        call write_iter_header_start(out_angular)
        write(aux, '(a2,18x)') 'Lx'
        call write_iter_header(out_angular, aux)
        write(aux, '(a2,18x)') 'Ly'
        call write_iter_header(out_angular, aux)
        write(aux, '(a2,18x)') 'Lz'
        call write_iter_header(out_angular, aux)
        call write_iter_nl(out_angular)

        !third line -> should hold the units. Now unused (assumes atomic units)
        call write_iter_string(out_angular, '#[Iter n.]')
        call write_iter_header(out_angular, '[a.u]')
        call write_iter_nl(out_angular)

        call td_write_print_header_end(out_angular)
      end if

      call write_iter_start(out_angular)
      call write_iter_double(out_angular, angular(1:3), 3)
      call write_iter_nl(out_angular)

    end if

    call pop_sub()
  end subroutine td_write_angular


  ! ---------------------------------------------------------
  subroutine td_write_multipole(out_multip, gr, st, lmax, kick, iter)
    integer(POINTER_SIZE), intent(in) :: out_multip
    type(grid_type),       intent(in) :: gr
    type(states_type),     intent(in) :: st
    integer,               intent(in) :: lmax
    type(kick_type),       intent(in) :: kick
    integer,               intent(in) :: iter

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

      write(aux, '(a15,3f18.12)') '# pol(1)       ', kick%pol(1:3, 1)
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,3f18.12)') '# pol(2)       ', kick%pol(1:3, 2)
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,3f18.12)') '# pol(3)       ', kick%pol(1:3, 3)
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,i2)')      '# direction    ', kick%pol_dir
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,i2)')      '# kick mode    ', kick%delta_strength_mode
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,f18.12)')  '# kick strengtph', kick%delta_strength
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,i2)')      '# Equiv. axis  ', kick%pol_equiv_axis
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      write(aux, '(a15,3f18.12)') '# wprime       ', kick%wprime(1:3)
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

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

    allocate(nuclear_dipole(1:3), multipole((lmax + 1)**2, st%d%nspin))
    call states_calculate_multipoles(gr, st, lmax, multipole)
    call geometry_dipole(gr%geo, nuclear_dipole)
    do is = 1, st%d%nspin
      multipole(2:4, is) = nuclear_dipole(1:3) - multipole(2:4, is)
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

    deallocate(nuclear_dipole, multipole)
    call pop_sub()
  end subroutine td_write_multipole


  ! ---------------------------------------------------------
  subroutine td_write_nbo(out_coords, gr, iter, ke, pe)
    integer(POINTER_SIZE), intent(in) :: out_coords
    type(grid_type),       intent(in) :: gr
    integer,               intent(in) :: iter
    FLOAT,                 intent(in) :: ke, pe

    integer :: i, j
    character(len=50) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    if(iter == 0) then
      call td_write_print_header_init(out_coords)

      ! first line: column names
      call write_iter_header_start(out_coords)
      call write_iter_header(out_coords, 'Ekin')
      call write_iter_header(out_coords, 'Epot')
      call write_iter_header(out_coords, 'Etot')

      do i = 1, gr%geo%natoms
        do j = 1, NDIM
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, gr%geo%natoms
        do j = 1, NDIM
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, gr%geo%natoms
        do j = 1, NDIM
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)

      ! second line: units
      call write_iter_string(out_coords, '#[Iter n.]')
      call write_iter_header(out_coords, '[' // trim(units_out%time%abbrev) // ']')
      call write_iter_string(out_coords, &
        'Energy in '      // trim(units_out%energy%abbrev)   //   &
        ', Positions in ' // trim(units_out%length%abbrev)   //   &
        ', Velocities in '// trim(units_out%velocity%abbrev) //   &
        ', Forces in '    // trim(units_out%force%abbrev))
      call write_iter_nl(out_coords)

      call td_write_print_header_end(out_coords)
    end if

    call write_iter_start(out_coords)
    call write_iter_double(out_coords,     ke /units_out%energy%factor, 1)
    call write_iter_double(out_coords,     pe /units_out%energy%factor, 1)
    call write_iter_double(out_coords, (ke+pe)/units_out%energy%factor, 1)

    do i = 1, gr%geo%natoms
      call write_iter_double(out_coords, gr%geo%atom(i)%x(1:NDIM)/units_out%length%factor,   NDIM)
    end do
    do i = 1, gr%geo%natoms
      call write_iter_double(out_coords, gr%geo%atom(i)%v(1:NDIM)/units_out%velocity%factor, NDIM)
    end do
    do i = 1, gr%geo%natoms
      call write_iter_double(out_coords, gr%geo%atom(i)%f(1:NDIM)/units_out%force%factor,    NDIM)
    end do
    call write_iter_nl(out_coords)

  end subroutine td_write_nbo


  ! ---------------------------------------------------------
  subroutine td_write_gsp(out_gsp, m, st, gs_st, dt, iter)
    integer(POINTER_SIZE), intent(in) :: out_gsp
    type(mesh_type),       intent(in) :: m
    type(states_type),     intent(in) :: st
    type(states_type),     intent(in) :: gs_st
    FLOAT,                 intent(in) :: dt
    integer,               intent(in) :: iter

    CMPLX :: gsp

    call push_sub('td_write.td_write_gsp')

    ! all processors calculate the projection
    gsp = zstates_mpdotp(m, 1, st, gs_st)

    if(mpi_grp_is_root(mpi_world)) then
      if(iter == 0) then
        call td_write_print_header_init(out_gsp)

        ! first line -> column names
        call write_iter_header_start(out_gsp)
        call write_iter_header(out_gsp, 'Re <Phi_gs|Phi(t)>')
        call write_iter_header(out_gsp, 'Im <Phi_gs|Phi(t)>')
        call write_iter_nl(out_gsp)

        ! second line -> units
        call write_iter_string(out_gsp, '#[Iter n.]')
        call write_iter_header(out_gsp, '[' // trim(units_out%time%abbrev) // ']')
        call write_iter_nl(out_gsp)

        call td_write_print_header_end(out_gsp)
      end if

      ! can not call write_iter_start, for the step is not 1
      call write_iter_int(out_gsp, iter, 1)
      call write_iter_double(out_gsp, iter*dt/units_out%time%factor,  1)
      call write_iter_double(out_gsp, real(gsp),  1)
      call write_iter_double(out_gsp, aimag(gsp), 1)
      call write_iter_nl(out_gsp)
    end if

    call pop_sub()
  end subroutine td_write_gsp


  ! ---------------------------------------------------------
  subroutine td_write_acc(out_acc, gr, st, h, dt, iter)
    integer(POINTER_SIZE),  intent(in)    :: out_acc
    type(grid_type),        intent(inout) :: gr
    type(states_type),      intent(inout) :: st
    type(hamiltonian_type), intent(inout) :: h
    FLOAT,                  intent(in)    :: dt
    integer,                intent(in)    :: iter

    integer :: i
    character(len=7) :: aux
    FLOAT :: acc(3)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    if(iter == 0) then
      call td_write_print_header_init(out_acc)

      ! first line -> column names
      call write_iter_header_start(out_acc)
      do i = 1, NDIM
        write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
        call write_iter_header(out_acc, aux)
      end do
      call write_iter_nl(out_acc)

      ! second line: units
      call write_iter_string(out_acc, '#[Iter n.]')
      call write_iter_header(out_acc, '[' // trim(units_out%time%abbrev) // ']')
      do i = 1, NDIM
        call write_iter_header(out_acc, '[' // trim(units_out%acceleration%abbrev) // ']')
      end do
      call write_iter_nl(out_acc)
      call td_write_print_header_end(out_acc)
    end if

    call td_calc_tacc(gr, st, h, acc, dt*i)

    call write_iter_start(out_acc)
    call write_iter_double(out_acc, acc/units_out%acceleration%factor, NDIM)
    call write_iter_nl(out_acc)

  end subroutine td_write_acc


  ! ---------------------------------------------------------
  subroutine td_write_laser(out_laser, gr, h, dt, iter)
    integer(POINTER_SIZE),  intent(in) :: out_laser
    type(grid_type),        intent(in) :: gr
    type(hamiltonian_type), intent(in) :: h
    FLOAT,                  intent(in) :: dt
    integer,                intent(in) :: iter

    integer :: i
    FLOAT :: field(3)
    character(len=80) :: aux

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    ! TODO -> confirm these stupid units, especially for the vector field
    if(iter == 0) then
      call td_write_print_header_init(out_laser)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', dt/units_out%time%factor, &
        " [", trim(units_out%time%abbrev), "]"
      call write_iter_string(out_laser, aux)
      call write_iter_nl(out_laser)

      ! second line -> column names
      call write_iter_header_start(out_laser)
      do i = 1, NDIM
        write(aux, '(a,i1,a)') 'E(', i, ')'
        call write_iter_header(out_laser, aux)
      end do
      do i = 1, NDIM
        write(aux, '(a,i1,a)') 'A(', i, ')'
        call write_iter_header(out_laser, aux)
      end do
      call write_iter_nl(out_laser)

      ! third line -> units
      call write_iter_string(out_laser, '#[Iter n.]')
      call write_iter_header(out_laser, '[' // trim(units_out%time%abbrev) // ']')

      aux = '[' // trim(units_out%energy%abbrev) // ' / ' // trim(units_inp%length%abbrev) // ']'
      do i = 1, NDIM
        call write_iter_header(out_laser, aux)
      end do

      aux = '[1/' // trim(units_inp%length%abbrev) // ']'
      do i = 1, NDIM
        call write_iter_header(out_laser, aux)
      end do
      call write_iter_nl(out_laser)
      call td_write_print_header_end(out_laser)
    end if

    call write_iter_start(out_laser)

    field = M_ZERO

    call laser_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
    field = field * units_inp%length%factor / units_inp%energy%factor
    call write_iter_double(out_laser, field, NDIM)

    call laser_vector_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
    field = field  * units_inp%length%factor
    call write_iter_double(out_laser, field, NDIM)

    call write_iter_nl(out_laser)

  end subroutine td_write_laser


  ! ---------------------------------------------------------
  subroutine td_write_el_energy(out_energy, h, iter)
    integer(POINTER_SIZE),  intent(in) :: out_energy
    type(hamiltonian_type), intent(in) :: h
    integer,                intent(in) :: iter

    integer :: i

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    if(iter == 0) then
      call td_write_print_header_init(out_energy)

      ! first line -> column names
      call write_iter_header_start(out_energy)
      call write_iter_header(out_energy, 'Total')
      call write_iter_header(out_energy, 'Ion-Ion')
      call write_iter_header(out_energy, 'Exchange')
      call write_iter_header(out_energy, 'Correlation')
      call write_iter_header(out_energy, 'Potentials')
      call write_iter_nl(out_energy)

      ! second line: units
      call write_iter_string(out_energy, '#[Iter n.]')
      call write_iter_header(out_energy, '[' // trim(units_out%time%abbrev) // ']')
      do i = 1, 5
        call write_iter_header(out_energy, '[' // trim(units_out%energy%abbrev) // ']')
      end do
      call write_iter_nl(out_energy)
      call td_write_print_header_end(out_energy)
    end if

    call write_iter_start(out_energy)
    call write_iter_double(out_energy, h%etot/units_out%energy%factor, 1)
    call write_iter_double(out_energy, h%eii /units_out%energy%factor, 1)
    call write_iter_double(out_energy, h%ex  /units_out%energy%factor, 1)
    call write_iter_double(out_energy, h%ec  /units_out%energy%factor, 1)
    call write_iter_double(out_energy, h%epot/units_out%energy%factor, 1)
    call write_iter_nl(out_energy)


  end subroutine td_write_el_energy


  ! ---------------------------------------------------------
  subroutine td_write_proj(out_proj, gr, st, gs_st, iter)
    integer(POINTER_SIZE), intent(in) :: out_proj
    type(grid_type),       intent(in) :: gr
    type(states_type),     intent(in) :: st
    type(states_type),     intent(in) :: gs_st
    integer,               intent(in) :: iter

    CMPLX, allocatable :: projections(:,:,:)
    character(len=20) :: aux
    integer :: ik, ist, uist
#if defined(HAVE_MPI)
    integer :: mpi_err, k
#endif

    call push_sub('td_write.td_write_proj')

    if(mpi_grp_is_root(mpi_world)) then
      if(iter == 0) then
        call td_write_print_header_init(out_proj)

        ! first line -> column names
        call write_iter_header_start(out_proj)
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do uist = 1, gs_st%nst
              write(aux, '(i3,a,i3)') ist, ' -> ', uist
              call write_iter_header(out_proj, 'Re {'//trim(aux)//'}')
              call write_iter_header(out_proj, 'Im {'//trim(aux)//'}')
            end do
          end do
        end do
        call write_iter_nl(out_proj)

        call td_write_print_header_end(out_proj)
      end if
    end if

    allocate(projections(st%nst, gs_st%nst, st%d%nik))
    call states_calc_projection(gr%m, st, gs_st, projections)
#if defined(HAVE_MPI)
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        k = st%node(ist)
        do uist = 1, gs_st%nst
          call mpi_bcast(projections(ist, uist, ik), 1, MPI_CMPLX, k, st%mpi_grp%comm, mpi_err)
        end do
      end do
    end do
#endif

    if(mpi_grp_is_root(mpi_world)) then
      call write_iter_start(out_proj)
      do ik = 1, st%d%nik
        do ist = 1, st%nst
          do uist = 1, gs_st%nst
            call write_iter_double(out_proj,  real(projections(ist, uist, ik)), 1)
            call write_iter_double(out_proj, aimag(projections(ist, uist, ik)), 1)
          end do
        end do
      end do
      call write_iter_nl(out_proj)
    end if

    deallocate(projections)
    call pop_sub()
  end subroutine td_write_proj


  ! ---------------------------------------------------------
  subroutine td_write_print_header_init(out)
    integer(POINTER_SIZE), intent(in) :: out

    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

  end subroutine td_write_print_header_init


  ! ---------------------------------------------------------
  subroutine td_write_print_header_end(out)
    integer(POINTER_SIZE), intent(in) :: out

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

  end subroutine td_write_print_header_end


#include "td_calc.F90"

end module td_write
