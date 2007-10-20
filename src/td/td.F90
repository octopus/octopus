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

module timedep_m
  use global_m
  use io_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use units_m
  use messages_m
  use mesh_m
  use external_pot_m
  use geometry_m
  use ground_state_m
  use hamiltonian_m
  use lib_oct_m
  use output_m
  use profiling_m
  use scf_m
  use states_m
  use restart_m
  use system_m
  use td_rti_m
  use td_write_m
  use v_ks_m
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  use PES_m
#endif
  use grid_m
  use spectrum_m
  use mpi_m
  use varinfo_m
  use math_m

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
    VELOCITY_VERLET = 1

  integer, parameter :: &
       EHRENFEST = 1,   &
       BO        = 2
  
  type td_t
    type(td_rti_t)    :: tr             ! contains the details of the time evolution
    type(scf_t)       :: scf
    FLOAT             :: dt             ! time step
    integer           :: max_iter       ! maximum number of iterations to perform
    integer           :: iter           ! the actual iteration
    integer           :: move_ions      ! how do we move the ions?
    logical           :: recalculate_gs ! Recalculate ground-state along the evolution.
    
    ! The *kick* used in "linear response in the time domain" calculations.
    type(kick_t)      :: kick

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    type(PES_t)       :: PESv
#endif
    FLOAT             :: mu
    integer           :: dynamics
  end type td_t


contains

  ! ---------------------------------------------------------
  subroutine td_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t)          :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical                   :: stopping
    integer                   :: i, ii, j, ik, ierr
    real(8)                   :: etime
    type(atom_t), allocatable :: atom1(:)
    FLOAT, allocatable        :: A_gauge_tmp(:)
    FLOAT, allocatable        :: A_gauge_dot_tmp(:)
    FLOAT, allocatable        :: A_gauge_ddot_tmp(:)

    call push_sub('td.td_run')

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    call td_init(sys, h, td)

    call states_distribute_nodes(st, sys%mc)
    ! Alocate complex wave-functions during time-propagation
    call states_allocate_wfns(st, gr%m, M_CMPLX)

    call init_wfs()

    call td_write_init(write_handler, gr, st, geo, (td%move_ions>0), h%ep%with_gauge_field, td%iter, td%dt )

    ! Calculate initial forces and kinetic energy
    if(td%move_ions > 0) then
      if(td%iter > 0) then
        call td_read_coordinates()
        call epot_generate(h%ep, gr, geo, sys%mc, st, h%reltype)
        geo%eii = ion_ion_energy(geo)
        h%eii = geo%eii
      end if

      call epot_forces(gr, geo, h%ep, st, td%iter*td%dt)
      geo%kinetic_energy = kinetic_energy(geo)

      call init_verlet()
      
    end if
    
    ! Calculate initial value of the gauge vector field
    if (h%ep%with_gauge_field) then
      if(td%iter > 0) then
        call td_read_gauge_field()
      end if
      call epot_generate_gauge_field(h%ep, gr, st)
      call init_verlet_gauge_field()
    end if

    if(td%iter == 0) call td_run_zero_iter()

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
    call write_info(1)
    call messages_print_stress(stdout)


    ii = 1
    stopping = .false.
    etime = loct_clock()
    ! This is the time-propagation loop. It starts at t=0 and finishes
    ! at td%max_iter*dt. The index i runs from 1 to td%max_iter, and
    ! step "i" means propagation from (i-1)*dt to i*dt.
    propagation: do i = td%iter, td%max_iter

      if(clean_stop()) stopping = .true.
      call profiling_in(C_PROFILING_TIME_STEP)
      

      if( td%move_ions > 0 .or. h%ep%extra_td_pot .ne. '0' .or. h%ep%with_gauge_field) then
        ! Move the ions: only half step, to obtain the external potential 
        ! in the middle of the time slice.
        if( td%move_ions > 0 ) call apply_verlet_1(td%dt*M_HALF)
	if( h%ep%with_gauge_field ) call apply_verlet_gauge_field_1(td%dt*M_HALF)
        call epot_generate(h%ep, gr, sys%geo, sys%mc, st, h%reltype, time = i*td%dt)

        ! Calculate the gauge vector potential at half step
      end if
      
      ! time iterate wavefunctions
      select case(td%dynamics)
      case(EHRENFEST)
        call td_rti_dt(sys%ks, h, gr, st, td%tr, i*td%dt, td%dt / td%mu, td%max_iter)
      case(BO)
        call scf_run(td%scf, sys%gr, geo, st, sys%ks, h, sys%outp, &
          gs_run = .false., verbosity = VERB_NO)
      end select


      ! mask function?
      call zvmask(gr, h, st)

      ! update density
      call states_calc_dens(st, NP, st%rho)

      if(td%move_ions > 0 .or. h%ep%with_gauge_field) then
        ! Now really move the full time step, from the original positions.
        if( td%move_ions > 0 ) then
	  geo%atom(1:geo%natoms) = atom1(1:geo%natoms)
          call apply_verlet_1(td%dt)
	end if
	if( h%ep%with_gauge_field ) then
	  h%ep%A_gauge(1:NDIM) = A_gauge_tmp(1:NDIM)
	  h%ep%A_gauge_dot(1:NDIM) = A_gauge_dot_tmp(1:NDIM)
	  h%ep%A_gauge_ddot(1:NDIM) = A_gauge_ddot_tmp(1:NDIM)
	  call apply_verlet_gauge_field_1(td%dt)
	  call epot_generate_gauge_field(h%ep, gr, st)
	end if
        call epot_generate(h%ep, gr, sys%geo, sys%mc, st, h%reltype, time = i*td%dt)
        
	if ( td%move_ions > 0 ) then
	  geo%eii = ion_ion_energy(geo)
          h%eii = geo%eii
	end if
      end if

      ! update hamiltonian and eigenvalues (fermi is *not* called)
      call v_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)

      ! Get the energies.
      call hamiltonian_energy(h, sys%gr, sys%geo, st, -1)

      ! Recalculate forces, update velocities...
      if(td%move_ions > 0) then
        call epot_forces(gr, sys%geo, h%ep, st, i*td%dt)
        call apply_verlet_2
        geo%kinetic_energy = kinetic_energy(geo)
      end if

      if (h%ep%with_gauge_field) then
        call epot_generate_gauge_field(h%ep, gr, st)
	call apply_verlet_gauge_field_2
      end if
      
      call td_write_iter(write_handler, gr, st, h, geo, td%kick, td%dt, i)

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
      call PES_doit(td%PESv, gr%m, st, ii, td%dt, h%ab_pot)
#endif

      ! write down data
      call check_point()

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(C_PROFILING_TIME_STEP)
      if (stopping) exit

    end do propagation

    call td_write_end(write_handler)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine check_point
      ! write info
      write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
        i*td%dt       / units_out%time%factor, &
        (h%etot + geo%kinetic_energy) / units_out%energy%factor, &
        loct_clock() - etime
      call write_info(1)
      etime = loct_clock()
      ii = ii + 1
      if(ii==sys%outp%iter+1 .or. i == td%max_iter .or. stopping) then ! output
        if(i == td%max_iter) sys%outp%iter = ii - 1
        ii = 1
        call td_save_restart(i)
        call td_write_data(write_handler, gr, st, h, sys%outp, geo, td%dt, i)
        if( (td%move_ions > 0) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.')
          fromScratch = .false.
          call ground_state_run(sys, h, fromScratch)
          call restart_read(trim(tmpdir)//'td', st, gr, geo, ierr, i)
          call messages_print_stress(stdout, "Time-Dependent simulation proceeds")
          write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
          call write_info(1)
          call messages_print_stress(stdout)
        end if
      end if
    end subroutine check_point

    ! ---------------------------------------------------------
    subroutine init_verlet
      select case(td%move_ions)
      case(VELOCITY_VERLET)
        ALLOCATE(atom1(geo%natoms), geo%natoms)
        atom1(1:geo%natoms) = geo%atom(1:geo%natoms)
      end select
    end subroutine init_verlet

    ! ---------------------------------------------------------
    subroutine end_verlet
      select case(td%move_ions)
      case(VELOCITY_VERLET)
        deallocate(atom1)
      end select
    end subroutine end_verlet

    ! ---------------------------------------------------------
    subroutine apply_verlet_1(tau)
      FLOAT, intent(in) :: tau
      select case(td%move_ions)
      case(VELOCITY_VERLET)
        atom1(1:geo%natoms) = geo%atom(1:geo%natoms)
        do j=1, geo%natoms
          if(geo%atom(j)%move) then
            geo%atom(j)%x(:) = geo%atom(j)%x(:) +  tau*geo%atom(j)%v(:) + &
              M_HALF*tau**2/geo%atom(j)%spec%weight * geo%atom(j)%f(:)
          end if
        end do
      end select
    end subroutine apply_verlet_1

    ! ---------------------------------------------------------
    subroutine apply_verlet_2
      select case(td%move_ions)
      case(VELOCITY_VERLET)      
        do j = 1, geo%natoms
          if(geo%atom(j)%move) then
            geo%atom(j)%v(:) = geo%atom(j)%v(:) + &
              td%dt/(M_TWO*geo%atom(j)%spec%weight) * &
              (atom1(j)%f(:) + geo%atom(j)%f(:))
          end if
        end do
      end select
    end subroutine apply_verlet_2

    ! ---------------------------------------------------------
    subroutine init_verlet_gauge_field
      ALLOCATE(A_gauge_tmp(NDIM), NDIM)
      ALLOCATE(A_gauge_dot_tmp(NDIM), NDIM)
      ALLOCATE(A_gauge_ddot_tmp(NDIM), NDIM)
      A_gauge_tmp(1:NDIM) = h%ep%A_gauge(1:NDIM)
      A_gauge_dot_tmp(1:NDIM) = h%ep%A_gauge_dot(1:NDIM)
      A_gauge_ddot_tmp(1:NDIM) = h%ep%A_gauge_ddot(1:NDIM)
    end subroutine init_verlet_gauge_field

    ! ---------------------------------------------------------
    subroutine end_verlet_gauge_field
      deallocate(A_gauge_tmp)
      deallocate(A_gauge_dot_tmp)
      deallocate(A_gauge_ddot_tmp)
    end subroutine end_verlet_gauge_field

   ! ---------------------------------------------------------
    subroutine apply_verlet_gauge_field_1(tau)
      FLOAT, intent(in) :: tau
      A_gauge_tmp(1:NDIM) = h%ep%A_gauge(1:NDIM)
      A_gauge_dot_tmp(1:NDIM) = h%ep%A_gauge_dot(1:NDIM)
      A_gauge_ddot_tmp(1:NDIM) = h%ep%A_gauge_ddot(1:NDIM)
      h%ep%A_gauge(:) = h%ep%A_gauge(:) +  tau*h%ep%A_gauge_dot(:) + M_HALF*tau**2*h%ep%A_gauge_ddot(:)
    end subroutine apply_verlet_gauge_field_1
   
   ! ---------------------------------------------------------
    subroutine apply_verlet_gauge_field_2
      h%ep%A_gauge_dot(:) = h%ep%A_gauge_dot(:)  +  td%dt/M_TWO * (A_gauge_ddot_tmp(:) + h%ep%A_gauge_ddot(:))
    end subroutine apply_verlet_gauge_field_2
   
   ! ---------------------------------------------------------
    subroutine end_()
      ! free memory
      deallocate(st%zpsi)
      if(td%move_ions>0) call end_verlet
      if (h%ep%with_gauge_field) call end_verlet_gauge_field
      call td_end(td)
      call pop_sub()
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()
      integer :: i, is, ierr, ist, jst
      character(len=50) :: filename
      FLOAT :: x
      logical :: only_userdef_istates
      C_POINTER :: blk
      type(states_t) :: stin
      CMPLX, allocatable :: rotation_matrix(:, :)

      if(.not.fromScratch) then
        call restart_read(trim(tmpdir)//'td', st, gr, geo, ierr, td%iter)
        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(tmpdir)//"td: Starting from scratch"
          call write_warning(1)

          fromScratch = .true.
        else
          ! read potential from previous interactions
          do i = 1, 2
            do is = 1, st%d%nspin
              write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'td/vprev_', i, is
              call dinput_function(trim(filename)//'.obf', gr%m, td%tr%v_old(1:NP, is, i), ierr)
              ! If we do not succed, try netcdf
              if(ierr > 0) call dinput_function(trim(filename)//'.ncdf', gr%m, td%tr%v_old(1:NP, is, i), ierr)
              if(ierr > 0) then
                write(message(1), '(3a)') 'Unsuccesfull write of "', trim(filename), '"'
                call write_fatal(1)
              end if
            end do
          end do

        end if
      end if

      if(fromScratch) then

        !%Variable OnlyUserDefinedInitialStates
        !%Type logical
        !%Default no
        !%Section States
        !%Description
        !% If true, then only user defined states from the block UserDefinedStates
        !% will be used as initial states for a time propagation. No attempt is made
        !% to load ground state orbitals from a previous ground state run.
        !%End
        call loct_parse_logical(check_inp('OnlyUserDefinedInitialStates'), .false., only_userdef_istates)

        if(.not. only_userdef_istates) then
          call restart_read(trim(tmpdir)//'gs', st, gr, geo, ierr)
          if(ierr.ne.0) then
            message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"gs'"
            message(2) = "Please run a ground-state calculation first!"
            call write_fatal(2)
          end if
        end if

        ! check if we should deploy user defined wavefunctions. 
        ! according to the settings in the input file the routine 
        ! overwrites orbitals that were read from restart/gs 
        if(loct_parse_isdef(check_inp('UserDefinedStates')).ne.0) then
          call states_read_user_def_orbitals(gr%m, st)
        end if

        !%Variable TransformStates
        !%Type block
        !%Default no
        !%Section States
        !%Description
        !% Before starting the td calculation, the initial states (that are
        !% read from the restart/gs directory, which should have been
        !% generated in a previous ground state calculation) can be "transformed"
        !% among themselves. The block TransformStates gives the transformation matrix
        !% to be used. The number of rows of the matrix should equal the number
        !% of the states present in the time-dependent calculation (the independent
        !% spin and k-point subspaces are all transformed equally); the number of
        !% columns should be equal to the number of states present in the
        !% restart/gs directory. This number may be different: for example,
        !% one could have run previously in "unocc" mode in order to obtain unoccupied
        !% Kohn-Sham states, and therefore restart/gs will contain more states.
        !% These states can be used in the transformation.
        !%
        !% Note that the code will not check the orthormality of the new states!
        !%
        !% Each line provides the coefficients of the new states, in terms of
        !% the old ones.
        !%End
        if(loct_parse_isdef(check_inp('TransformStates')).ne.0) then
          if(loct_parse_block(check_inp('TransformStates'), blk) == 0) then
            stin = st
            deallocate(stin%zpsi)
            call restart_look_and_read("tmp", stin, gr, sys%geo, ierr)
            ALLOCATE(rotation_matrix(st%nst, stin%nst), st%nst*stin%nst)
            do ist = 1, st%nst
              do jst = 1, stin%nst
                call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call rotate_states(gr%m, st, stin, rotation_matrix)
            deallocate(rotation_matrix)
            call states_end(stin)
          else
            message(1) = '"TransformStates" has to be specified as block.'
            call write_info(1)
            call input_error('TransformStates')
          end if
        end if

      end if

      if(fromScratch) call modify_occs()

      call states_calc_dens(st, NP, st%rho)
      call v_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)
      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call hamiltonian_span(h, minval(gr%m%h(1:NDIM)), x)
      call hamiltonian_energy(h, gr, geo, st, -1)

    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      call push_sub('td.td_run_zero_iter')

      call td_write_iter(write_handler, gr, st, h, geo, td%kick, td%dt, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      call apply_delta_field(td%kick)

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
      integer :: i, j
      CMPLX   :: c(2), kick
      FLOAT   :: ylm, r
      FLOAT   :: x(MAX_DIM)
      FLOAT, allocatable :: kick_function(:)

      call push_sub('td.apply_delta_field')

      ! The wave-functions at time delta t read
      ! psi(delta t) = psi(t) exp(i k x)
      delta_strength: if(k%delta_strength .ne. M_ZERO) then

        ALLOCATE(kick_function(NP), NP)
        if(k%n_multipoles > 0) then
          kick_function = M_ZERO
          do j = 1, k%n_multipoles
            do i = 1, NP
              call mesh_r(gr%m, i, r, x = x)
              ylm = loct_ylm(x(1), x(2), x(3), k%l(j), k%m(j))
              kick_function(i) = kick_function(i) + k%weight(j) * (r**k%l(j)) * ylm 
            end do
          end do
        else
          do i = 1, NP
            kick_function(i) = sum(gr%m%x(i, 1:NDIM)*k%pol(1:NDIM, k%pol_dir))
          end do
        end if
        
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
          kick = M_zI * k%delta_strength * kick_function(i)

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

        ! the nuclei velocity will be changed by
        ! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
        ! where M and Z are the ionic mass and charge, respectively.
        if(td%move_ions > 0  .and. k%delta_strength .ne. M_ZERO) then
          do i = 1, geo%natoms
            geo%atom(i)%v(1:NDIM) = geo%atom(i)%v(1:NDIM) - &
              k%delta_strength_mode*k%pol(1:NDIM, k%pol_dir)*geo%atom(i)%spec%z_val / geo%atom(i)%spec%weight
          end do
        end if

        deallocate(kick_function)
      end if delta_strength

      call pop_sub()
    end subroutine apply_delta_field


    ! ---------------------------------------------------------
    subroutine td_read_coordinates() ! reads the pos and vel from coordinates file
      integer :: i, iunit, record_length
      call push_sub('td.td_read_coordinates')

      record_length = 28 + 3*geo%natoms*3*20
      call io_assign(iunit)
      open(unit = iunit, file = 'td.general/coordinates', &
        action='read', status='old', recl = record_length)
      if(iunit < 0) then
        message(1) = "Could not open file 'td.general/coordinates'"
        message(2) = "Starting simulation from initial geometry"
        call write_warning(2)
        return
      end if

      call io_skip_header(iunit)
      do i = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

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
      call pop_sub()
    end subroutine td_read_coordinates

    ! ---------------------------------------------------------
    subroutine td_read_gauge_field()
      
      integer :: i, iunit, record_length
      
      call push_sub('td.td_read_gauge_field')

      record_length = 28 + 3*3*20
      call io_assign(iunit)
      open(unit = iunit, file = 'td.general/A_gauge', &
        action='read', status='old', recl = record_length)
      if(iunit < 0) then
        message(1) = "Could not open file 'td.general/A_gauge'"
        message(2) = "Starting simulation from initial value"
        call write_warning(2)
        return
      end if

      call io_skip_header(iunit)
      do i = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      ! TODO: units are missing
      read(iunit, '(3es20.12)', advance='no') h%ep%A_gauge(1:NDIM)
      read(iunit, '(3es20.12)', advance='no') h%ep%A_gauge_dot(1:NDIM)
      read(iunit, '(3es20.12)', advance='no') h%ep%A_gauge_ddot(1:NDIM)
      !           !* units_out%length%factor
 
      call io_close(iunit)
      call pop_sub()
    end subroutine td_read_gauge_field

    ! ---------------------------------------------------------
    subroutine td_save_restart(iter)
      integer, intent(in) :: iter

      integer :: i, is, ierr
      character(len=256) :: filename

      call push_sub('td.td_save_restart')

      ! first write resume file
      call restart_write(trim(tmpdir)//'td', st, gr, ierr, iter)
      if(ierr.ne.0) then
        message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'td"'
        call write_fatal(1)
      end if

      ! write potential from previous interactions
      if(mpi_grp_is_root(st%mpi_grp)) then
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a6,i2.2,i3.3)') 'vprev_', i, is
            call doutput_function(restart_format, trim(tmpdir)//"td", &
              filename, gr%m, gr%sb, td%tr%v_old(1:NP, is, i), M_ONE, ierr, is_tmp = .true.)
            if(ierr.ne.0) then
              write(message(1), '(3a)') 'Unsuccesfull write of "', trim(filename), '"'
              call write_fatal(1)
            end if
          end do
        end do
      end if

      call pop_sub()
    end subroutine td_save_restart

    ! ---------------------------------------------------------
    subroutine modify_occs()
      C_POINTER :: blk
      integer  :: nrow
      integer  :: spin, state
      FLOAT    :: new_occ
      blk = int(0, POINTER_SIZE)
      
      if(loct_parse_block(check_inp('ModifyOccupations'), blk) == 0) then
        nrow = loct_parse_block_n(blk)
        
        do i=0,(nrow-1)
          call loct_parse_block_int(blk, i, 0, spin)
          call loct_parse_block_int(blk, i, 1, state)
          call loct_parse_block_float(blk, i, 2, new_occ)
          
          
          if ( (state <= sys%st%st_end) .and. (spin <= sys%st%d%nik)) then 
            write(message(1), '(2i2,f12.6,a,f12.6)') spin,state, st%occ(state,spin), ' ->', new_occ
            call write_info(1)
            st%occ(state,spin) = new_occ
          end if

        end do
        call loct_parse_block_end(blk)

      end if

    end subroutine modify_occs

  end subroutine td_run

#include "td_init.F90"

end module timedep_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
