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

module timedep_m
  use global_m
  use messages_m
  use profiling_m
  use lib_oct_m
  use lib_oct_parser_m
  use geometry_m
  use mesh_m
  use grid_m
  use mesh_function_m
  use functions_m
  use states_m
  use output_m
  use restart_m
  use lasers_m
  use v_ks_m
  use hamiltonian_m
  use external_pot_m
  use system_m
  use td_rti_m
  use td_write_m
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  use PES_m
#endif
  use grid_m
  use spectrum_m
  use mpi_m
  use varinfo_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_init,              &
    td_end

  ! Parameters.
  integer, parameter ::   &
    STATIC_IONS     = 0,  &
    NORMAL_VERLET   = 3,  &
    VELOCITY_VERLET = 4

  type td_t
    type(td_rti_t) :: tr             ! contains the details of the time evolution
    FLOAT             :: dt             ! time step
    integer           :: max_iter       ! maximum number of iterations to perform
    integer           :: iter           ! the actual iteration
    integer           :: epot_regenerate! Every epot_regenerate, the external potential
    ! regenerated *exactly*.
    integer           :: move_ions      ! how do we move the ions?

    ! The *kick* used in "linear response in the time domain" calculations.
    type(kick_t)   :: kick

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    type(PES_t) :: PESv
#endif
  end type td_t


contains

  ! ---------------------------------------------------------
  subroutine td_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t) :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical :: stopping
    integer :: i, ii, j, idim, ist, ik
    FLOAT, allocatable ::  x1(:,:), x2(:,:), f1(:,:) ! stuff for verlet
    FLOAT :: etime

    call push_sub('td.td_run')

    ! some shortcuts
    gr  => sys%gr
    geo => sys%gr%geo
    st  => sys%st

    call td_init(gr, td, st, sys%outp)

    call states_distribute_nodes(st, sys%mc)

    call zstates_allocate_wfns(st, gr%m)

    call init_wfs()

    call td_write_init(write_handler, gr, st, geo, (td%move_ions>0), (h%ep%no_lasers>0), td%iter, td%dt )

    ! Calculate initial forces and kinetic energy
    if(td%move_ions > 0) then
      if(td%iter > 0) then
        call td_read_nbo()
        call epot_generate(h%ep, gr, st, h%reltype)
        geo%eii = ion_ion_energy(geo)
        h%eii = geo%eii
      end if

      call zepot_forces(gr, h%ep, st, td%iter*td%dt)
      geo%kinetic_energy = kinetic_energy(geo)
      select case(td%move_ions)
      case(NORMAL_VERLET)
        ALLOCATE(x1(geo%natoms, NDIM), NDIM*geo%natoms)
        ALLOCATE(x2(geo%natoms, NDIM), NDIM*geo%natoms)
        do j = 1, geo%natoms
          if(geo%atom(j)%move) then
            x1(j, :) = geo%atom(j)%x(:) - td%dt*geo%atom(j)%v(:) + &
              M_HALF * td%dt**2/geo%atom(j)%spec%weight * &
              geo%atom(j)%f(:)
          else
            x1(j, :) = geo%atom(j)%x(:)
          end if
        end do
      case(VELOCITY_VERLET)
        ALLOCATE(f1(geo%natoms, NDIM), NDIM*geo%natoms)
      end select
    end if


    if(td%iter == 0) then
      call apply_delta_field(td%kick)
      call td_run_zero_iter()
    end if
    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
    call write_info(1)

    ii = 1
    stopping = .false.
    etime = loct_clock()
    do i = td%iter, td%max_iter
      if(clean_stop()) stopping = .true.
      call profiling_in(C_PROFILING_TIME_STEP)
      ! Move the ions.
      if( td%move_ions > 0 ) then
        select case(td%move_ions)
        case(NORMAL_VERLET)
          x2 = x1
          do j = 1, geo%natoms
            if(geo%atom(j)%move) then
              x1(j, :) = geo%atom(j)%x(:)
              geo%atom(j)%x(:) = M_TWO*x1(j, :) - x2(j, :) + &
                td%dt**2/geo%atom(j)%spec%weight * geo%atom(j)%f(:)
              geo%atom(j)%v(:) = (geo%atom(j)%x(:) - x2(j, :)) / (M_TWO*td%dt)
            end if
          end do
        case(VELOCITY_VERLET)
          do j=1, geo%natoms
            if(geo%atom(j)%move) then
              geo%atom(j)%x(:) = geo%atom(j)%x(:) +  td%dt*geo%atom(j)%v(:) + &
                M_HALF*td%dt**2/geo%atom(j)%spec%weight * geo%atom(j)%f(:)
            end if
          end do
        end select

        if(mod(i, td%epot_regenerate) == 0) then
          call epot_generate(h%ep, gr, st, h%reltype)
        else
          call epot_generate(h%ep, gr, st, h%reltype, fast_generation = .true.)
        end if
        geo%eii = ion_ion_energy(geo)
        h%eii = geo%eii
      end if

      ! time iterate wavefunctions
      call td_rti_dt(sys%ks, h, gr, st, td%tr, i*td%dt, td%dt)

      ! mask function?
      call zvmask(gr, h, st)

      ! update density
      call zstates_calc_dens(st, NP, st%rho)

      ! update hamiltonian and eigenvalues (fermi is *not* called)
      call zv_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)
      call hamiltonian_energy(h, st, geo%eii, -1)

      ! Recalculate forces, update velocities...
      if(td%move_ions > 0) then
        if(td%move_ions == VELOCITY_VERLET) then
          do j = 1, geo%natoms
            f1(j, :) = geo%atom(j)%f(:)
          end do
        end if
        call zepot_forces(gr, h%ep, st, i*td%dt)
        if(td%move_ions == VELOCITY_VERLET) then
          do j = 1, geo%natoms
            if(geo%atom(j)%move) then
              geo%atom(j)%v(:) = geo%atom(j)%v(:) + &
                td%dt/(M_TWO*geo%atom(j)%spec%weight) * &
                (f1(j, :) + geo%atom(j)%f(:))
            end if
          end do
        end if
        geo%kinetic_energy = kinetic_energy(geo)
      end if

      call td_write_iter(write_handler, gr, st, h, geo, td%kick, td%dt, i)


#if !defined(DISABLE_PES) && defined(HAVE_FFT)
      call PES_doit(td%PESv, gr%m, st, ii, td%dt, h%ab_pot)
#endif

      ! write info
      write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
        i*td%dt       / units_out%time%factor, &
        (h%etot + geo%kinetic_energy) / units_out%energy%factor, &
        (loct_clock() - etime)/1e6
      call write_info(1)
      etime = loct_clock()

      ! write down data
      ii = ii + 1
      if(ii==sys%outp%iter+1 .or. i == td%max_iter .or. stopping) then ! output
        if(i == td%max_iter) sys%outp%iter = ii - 1
        ii = 1
        call td_save_restart(i)
        call td_write_data(write_handler, gr, st, h, sys%outp, geo, td%dt, i)
      end if

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(C_PROFILING_TIME_STEP)
      if (stopping) exit
    end do

    call td_write_end(write_handler)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine end_()
      ! free memory
      deallocate(st%zpsi)
      call td_end(td)
      call pop_sub()
    end subroutine end_


    ! ---------------------------------------------------------
    subroutine init_wfs()
      integer :: i, is, ierr
      character(len=50) :: filename
      FLOAT :: x

      if(.not.fromScratch) then
        call zrestart_read(trim(tmpdir)//'restart_td', st, gr, ierr, td%iter)
        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(tmpdir)//"restart_td: Starting from scratch"
          call write_warning(1)

          fromScratch = .true.
        else
          ! read potential from previous interactions
          do i = 1, 2
            do is = 1, st%d%nspin
              write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'restart_td/vprev_', i, is
              call dinput_function(filename, gr%m, td%tr%v_old(1:NP, is, i), ierr)
            end do
          end do

        end if
      end if

      if(fromScratch) then
        call zrestart_read(trim(tmpdir)//'restart_gs', st, gr, ierr)
        if(ierr.ne.0) then
          message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
          message(2) = "Please run a ground-state calculation first!"
          call write_fatal(2)
        end if

        ! check if we should deploy user defined wavefunctions. 
        ! according to the settings in the input file the routine 
        ! overwrites orbitals that were read from restart_gs 
        if(loct_parse_isdef(check_inp('UserDefinedStates')).ne.0) then
          call states_read_user_def_orbitals(gr%m, st)
        end if
      end if

      call zstates_calc_dens(st, NP, st%rho)
      call zv_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)
      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_BCAST(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, i)
      end if
#endif
      call hamiltonian_span(h, minval(gr%m%h(1:NDIM)), x)
      call hamiltonian_energy(h, st, geo%eii, -1)

    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      call push_sub('td.td_run_zero_iter')

      call io_mkdir('td.general')
      call td_write_iter(write_handler, gr, st, h, geo, td%kick, td%dt, 0)
      call td_save_restart(0)
      call td_write_data(write_handler, gr, st, h, sys%outp, geo, td%dt, 0)
      call td_rti_run_zero_iter(h, td%tr)

      call pop_sub()
    end subroutine td_run_zero_iter


    ! ---------------------------------------------------------
    ! Applies the delta function electric field E(t) = E_0 delta(t)
    ! where E_0 = - k \hbar / e
    subroutine apply_delta_field(k)
      type(kick_t), intent(in) :: k
      integer :: i
      CMPLX   :: c(2), kick

      call push_sub('td.apply_delta_field')

      ! The wave-functions at time delta t read
      ! psi(delta t) = psi(t) exp(i k x)
      if(k%delta_strength .ne. M_ZERO) then
        write(message(1),'(a,f11.6)')  'Info: Applying delta kick: k = ', k%delta_strength
        select case (k%delta_strength_mode)
        case (KICK_DENSITY_MODE)
          message(2) = "Info: Delta kick mode: Density mode"
        case (KICK_SPIN_MODE)
          message(2) = "Info: Delta kick mode: Spin mode"
        case (KICK_SPIN_DENSITY_MODE)
          message(2) = "Info: Delta kick mode: Density + Spin modes"
        end select
        call write_info(2)
        do i = 1, NP
          kick = M_zI * k%delta_strength * sum(gr%m%x(i, 1:NDIM)*k%pol(1:NDIM, k%pol_dir))

          select case (k%delta_strength_mode)
          case (KICK_DENSITY_MODE)
            c(1) = exp(kick)
            st%zpsi(i,:,:,:) = c(1) * st%zpsi(i,:,:,:)

          case (KICK_SPIN_MODE)
            c(1) = exp(kick)
            c(2) = exp(-kick)
            select case (st%d%ispin)
            case (SPIN_POLARIZED)
              do ik = 1, st%d%nik, 2
                st%zpsi(i,:,:,ik)   = c(1) * st%zpsi(i,:,:,ik)
                st%zpsi(i,:,:,ik+1) = c(2) * st%zpsi(i,:,:,ik+1)
              end do
            case (SPINORS)
              st%zpsi(i,1,:,:) = c(1) * st%zpsi(i,1,:,:)
              st%zpsi(i,2,:,:) = c(2) * st%zpsi(i,2,:,:)
            end select

          case (KICK_SPIN_DENSITY_MODE)
            c(1) = exp(M_TWO*kick)
            select case (st%d%ispin)
            case (SPIN_POLARIZED)
              do ik = 1, st%d%nik, 2
                st%zpsi(i,:,:,ik) = c(1) * st%zpsi(i,:,:,ik)
              end do
            case (SPINORS)
              st%zpsi(i,1,:,:) = c(1) * st%zpsi(i,1,:,:)
            end select

          end select
        end do
      end if

      ! the nuclei velocity will be changed by
      ! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
      ! where M and Z are the ionic mass and charge, respectively.
      if(td%move_ions > 0) then
        do i = 1, geo%natoms
          geo%atom(i)%v(1:NDIM) = geo%atom(i)%v(1:NDIM) - &
            k%delta_strength_mode*k%pol(1:NDIM, k%pol_dir)*geo%atom(i)%spec%z_val / geo%atom(i)%spec%weight
        end do
      end if

      call pop_sub()
    end subroutine apply_delta_field


    ! ---------------------------------------------------------
    subroutine td_read_nbo() ! reads the pos and vel from coordinates file
      integer :: i, iunit, record_length

      record_length = 100 + 3*geo%natoms*3*20
      call io_assign(iunit)
      open(unit = iunit, file = 'td.general/coordinates', &
        action='read', status='old', recl = record_length)
      if(iunit < 0) then
        message(1) = "Could not open file 'td.general/coordinates'"
        message(2) = "Starting simulation from initial geometry"
        call write_warning(2)
        return
      end if

      read(iunit, *); read(iunit, *) ! skip header
      do i = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(88x)', advance='no') ! skip unrelevant information

      do i = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(i)%x(1:NDIM)
        geo%atom(i)%x(:) = geo%atom(i)%x(:) * units_out%length%factor
      end do
      do i = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(i)%v(1:NDIM)
        geo%atom(i)%v(:) = geo%atom(i)%v(:) * units_out%velocity%factor
      end do
      do i = 1, geo%natoms
        read(iunit, '(3es20.12)', advance='no') geo%atom(i)%f(1:NDIM)
        geo%atom(i)%f(:) = geo%atom(i)%f(:) * units_out%force%factor
      end do

      call io_close(iunit)
    end subroutine td_read_nbo


    ! ---------------------------------------------------------
    subroutine td_save_restart(iter)
      integer, intent(in) :: iter

      integer :: i, is, ierr
      character(len=256) :: filename

      call push_sub('td.td_save_restart')

      ! first write resume file
      call zrestart_write(trim(tmpdir)//'restart_td', st, gr, ierr, iter)
      if(ierr.ne.0) then
        message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_td"'
        call write_fatal(1)
      end if

      ! write potential from previous interactions
      if(mpi_grp_is_root(st%mpi_grp)) then
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a6,i2.2,i3.3)') 'vprev_', i, is

            call doutput_function(restart_format, trim(tmpdir)//"restart_td", &
              filename, gr%m, gr%sb, td%tr%v_old(1:NP, is, i), M_ONE, ierr)

            if(ierr.ne.0) then
              write(message(1), '(3a)') 'Unsuccesfull write of "', trim(filename), '"'
              call write_fatal(1)
            end if
          end do
        end do
      end if

      call pop_sub()
    end subroutine td_save_restart

  end subroutine td_run

#include "td_init.F90"

end module timedep_m
