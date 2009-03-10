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
  use cpmd_m
  use energy_m
  use global_m
  use io_m
  use datasets_m
  use io_function_m
  use loct_math_m
  use loct_parser_m
  use units_m
  use messages_m
  use mesh_m
  use external_pot_m
  use gauge_field_m
  use geometry_m
  use ground_state_m
  use h_sys_output_m
  use hamiltonian_m
  use ion_dynamics_m
  use loct_m
  use profiling_m
  use projector_m
  use scf_m
  use scissor_m
  use states_m
  use states_dim_m
  use states_calc_m
  use restart_m
  use system_m
  use td_rti_m
  use td_write_m
  use v_ks_m
#if !defined(DISABLE_PES)
  use PES_m
#endif
  use grid_m
  use spectrum_m
  use mpi_m
  use varinfo_m
  use math_m
  use lasers_m

  implicit none

  private
  public ::               &
    td_t,                 &
    td_run,               &
    td_init,              &
    td_end

  ! Parameters.
  integer, parameter :: &
       EHRENFEST = 1,   &
       BO        = 2,   &
       CP        = 3
  
  type td_t 
    type(td_rti_t)       :: tr             ! contains the details of the time evolution
    type(scf_t)          :: scf
    type(ion_dynamics_t) :: ions
    type(cpmd_t)         :: cp_propagator
    FLOAT                :: dt             ! time step
    integer              :: max_iter       ! maximum number of iterations to perform
    integer              :: iter           ! the actual iteration
    logical              :: recalculate_gs ! Recalculate ground-state along the evolution.
    
    ! The *kick* used in "linear response in the time domain" calculations.
    type(kick_t)         :: kick

#if !defined(DISABLE_PES)
    type(PES_t)          :: PESv
#endif
    FLOAT                :: mu
    integer              :: dynamics
    FLOAT                :: scissor
  end type td_t


contains

  ! ---------------------------------------------------------
  subroutine td_run(sys, hm, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(inout) :: fromScratch

    type(td_t)                :: td
    type(td_write_t)          :: write_handler
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st
    type(geometry_t), pointer :: geo
    logical                   :: stopping
    integer                   :: i, ii, ik, ierr, iatom, ist
    real(8)                   :: etime
    logical                   :: generate
    type(gauge_force_t)       :: gauge_force
    CMPLX, allocatable        :: gspsi(:, :, :, :)

    type(profile_t), save :: prof

    call push_sub('td.td_run')

    ! some shortcuts
    gr  => sys%gr
    geo => sys%geo
    st  => sys%st

    call td_init(sys, hm, td)

    call states_distribute_nodes(st, sys%mc)

    ! Alocate wave-functions during time-propagation
    if(td%dynamics == EHRENFEST) then
      !complex wfs are required for ehrefenst
      call states_allocate_wfns(st, gr%mesh, M_CMPLX)
    else
      call states_allocate_wfns(st, gr%mesh)
    end if

    ! CP has to be initialized after wavefunction type is set
    if(td%dynamics == CP) call cpmd_init(td%cp_propagator, sys%gr, sys%st)

    call init_wfs()

    if(td%scissor > M_EPSILON) then
      ALLOCATE(gspsi(1:gr%mesh%np, 1:st%d%dim, 1:st%nst, 1:st%d%nik), 1)
      gspsi(1:gr%mesh%np, 1:st%d%dim, 1:st%nst, 1:st%d%nik) = st%zpsi(1:gr%mesh%np, 1:st%d%dim, 1:st%nst, 1:st%d%nik)
      call scissor_init(hm%scissor, st, td%scissor, gspsi)
    end if

    call td_write_init(write_handler, gr, st, hm, geo, &
         ion_dynamics_ions_move(td%ions), gauge_field_is_applied(hm%ep%gfield), td%iter, td%max_iter, td%dt)

    ! Calculate initial forces and kinetic energy
    if(ion_dynamics_ions_move(td%ions)) then
      if(td%iter > 0) then
        call td_read_coordinates()
        call hamiltonian_epot_generate(hm, gr, geo, st)
      end if

      call epot_forces(gr, geo, hm%ep, st, td%iter*td%dt)

      geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)

    end if

    ! Calculate initial value of the gauge vector field
    if (gauge_field_is_applied(hm%ep%gfield)) then

      if(td%iter > 0) then
        call td_read_gauge_field()
      else
        call gauge_field_init_vec_pot(hm%ep%gfield, gr%mesh, gr%sb, st, td%dt)
      end if

      call gauge_field_get_force(hm%ep%gfield, gr, geo, hm%ep%proj, hm%phase, st, gauge_force)

      do iatom = 1, geo%natoms
         call projector_init_phases(hm%ep%proj(iatom), gr%sb, st%d%nik, st%d%kpoints, &
              vec_pot = gauge_field_get_vec_pot(hm%ep%gfield)/P_c, vec_pot_var = hm%ep%a_static)
      end do

    end if

    if(td%iter == 0) call td_run_zero_iter()

    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    call messages_print_stress(stdout, "Time-Dependent Simulation")
    if(td%dynamics /= CP) then 
      write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
    else
      write(message(1), '(a7,1x,a14,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'CP Energy ', 'Elapsed Time '
    end if
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
      call profiling_in(prof, "TIME_STEP")

      if(td%scissor > M_EPSILON) then
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            gspsi(1:gr%mesh%np, 1:st%d%dim, ist, ik) = &
                 exp(-M_HALF*M_ZI*td%dt*st%eigenval(ist, ik))*gspsi(1:gr%mesh%np, 1:st%d%dim, ist, ik)
          end do
        end do
      end if

      ! time iterate wavefunctions
      select case(td%dynamics)
      case(EHRENFEST)
        if(ion_dynamics_ions_move(td%ions)) then
          call td_rti_dt(sys%ks, hm, gr, st, td%tr, i*td%dt, td%dt / td%mu, td%max_iter, i, gauge_force, &
            ions = td%ions, geo = sys%geo, ionic_dt = td%dt)
        else
          call td_rti_dt(sys%ks, hm, gr, st, td%tr, i*td%dt, td%dt / td%mu, td%max_iter, i, gauge_force)
        end if
      case(BO)
        call scf_run(td%scf, sys%gr, geo, st, sys%ks, hm, sys%outp, gs_run = .false., verbosity = VERB_NO)
      case(CP)
        if(states_are_real(st)) then
          call dcpmd_propagate(td%cp_propagator, sys%gr, hm, st, i, td%dt)
        else
          call zcpmd_propagate(td%cp_propagator, sys%gr, hm, st, i, td%dt)
        end if
      end select

      ! mask function?
      call zvmask(gr, hm, st)

      ! update density
      call states_calc_dens(st, NP)

      generate = .false.

      if(ion_dynamics_ions_move(td%ions)) then 
        if(td%dynamics /= EHRENFEST .or. .not. td_rti_ions_are_propagated(td%tr)) then
          call ion_dynamics_propagate(td%ions, sys%gr%sb, sys%geo,  i*td%dt, td%dt)
          generate = .true.
        end if
      end if

      if(gauge_field_is_applied(hm%ep%gfield) .and. .not. td_rti_ions_are_propagated(td%tr)) then
        call gauge_field_propagate(hm%ep%gfield, gauge_force, td%dt)
      end if
      
      if(generate) call hamiltonian_epot_generate(hm, gr, sys%geo, st, time = i*td%dt)

      ! update hamiltonian and eigenvalues (fermi is *not* called)
      call v_ks_calc(gr, sys%ks, hm, st, calc_eigenval=.true.)

      ! Get the energies.
      call total_energy(hm, sys%gr, st, -1)

      if (td%dynamics == CP) then
        if(states_are_real(st)) then
          call dcpmd_propagate_vel(td%cp_propagator, sys%gr, hm, st, td%dt)
        else
          call zcpmd_propagate_vel(td%cp_propagator, sys%gr, hm, st, td%dt)
        end if
      end if
      
      ! Recalculate forces, update velocities...
      if(ion_dynamics_ions_move(td%ions)) then
        call epot_forces(gr, sys%geo, hm%ep, st, i*td%dt)

        call ion_dynamics_propagate_vel(td%ions, sys%geo)

        geo%kinetic_energy = ion_dynamics_kinetic_energy(geo)
      end if

      if(gauge_field_is_applied(hm%ep%gfield)) then

        call gauge_field_get_force(hm%ep%gfield, gr, geo, hm%ep%proj, hm%phase, st, gauge_force)

        call gauge_field_propagate_vel(hm%ep%gfield, gauge_force, td%dt)

        do iatom = 1, geo%natoms
          call projector_init_phases(hm%ep%proj(iatom), gr%sb, st%d%nik, st%d%kpoints, &
               vec_pot = gauge_field_get_vec_pot(hm%ep%gfield)/P_c, vec_pot_var = hm%ep%a_static)
        end do

      end if

      call td_write_iter(write_handler, gr, st, hm, geo, td%kick, td%dt, i)

#if !defined(DISABLE_PES)
      call PES_doit(td%PESv, gr%mesh, st, ii, td%dt, hm%ab_pot)
#endif

      if(td%scissor > M_EPSILON) then
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            gspsi(1:gr%mesh%np, 1:st%d%dim, ist, ik) = &
                 exp(-M_HALF*M_ZI*td%dt*st%eigenval(ist, ik))*gspsi(1:gr%mesh%np, 1:st%d%dim, ist, ik)
          end do
        end do
      end if

      ! write down data
      call check_point()

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()

      call profiling_out(prof)
      if (stopping) exit

    end do propagation

    call td_write_end(write_handler)
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine check_point
      ! write info
      if(td%dynamics /= CP) then 
        write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
             i*td%dt       / units_out%time%factor, &
             (hm%etot + geo%kinetic_energy) / units_out%energy%factor, &
             loct_clock() - etime
      else
        write(message(1), '(i7,1x,3f14.6,f14.3, i10)') i, &
             i*td%dt       / units_out%time%factor, &
             (hm%etot + geo%kinetic_energy) / units_out%energy%factor, &
             (hm%etot + geo%kinetic_energy + cpmd_electronic_energy(td%cp_propagator)) / units_out%energy%factor, &
             loct_clock() - etime
      end if
      call write_info(1)
      etime = loct_clock()
      ii = ii + 1
      if(ii==sys%outp%iter+1 .or. i == td%max_iter .or. stopping) then ! output
        if(i == td%max_iter) sys%outp%iter = ii - 1
        ii = 1
        call td_save_restart(i)
        call td_write_data(write_handler, gr, st, hm, sys%outp, geo, i)
        if( (ion_dynamics_ions_move(td%ions)) .and. td%recalculate_gs) then
          call messages_print_stress(stdout, 'Recalculating the ground state.')
          fromScratch = .false.
          call ground_state_run(sys, hm, fromScratch)
          call restart_read(trim(restart_dir)//'td', st, gr, geo, ierr, i)
          call messages_print_stress(stdout, "Time-Dependent simulation proceeds")
          if(td%dynamics /= CP) then 
            write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
          else
            write(message(1), '(a7,1x,a14,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'CP Energy ', 'Elapsed Time '
          end if
          call write_info(1)
          call messages_print_stress(stdout)
        end if
      end if
    end subroutine check_point

   ! ---------------------------------------------------------
    subroutine end_()
      ! free memory
      if(td%dynamics == CP) call cpmd_end(td%cp_propagator)
      call states_deallocate_wfns(st)
      call ion_dynamics_end(td%ions)
      call td_end(td)
      call pop_sub()
    end subroutine end_

    ! ---------------------------------------------------------
    subroutine init_wfs()
      integer :: i, is, ierr, ist, jst, freeze_orbitals
      character(len=50) :: filename
      FLOAT :: x
      type(block_t) :: blk
      type(states_t) :: stin
      CMPLX, allocatable :: rotation_matrix(:, :)

      if(.not.fromscratch) then
        call restart_read(trim(tmpdir)//'td', st, gr, geo, ierr, td%iter)
        
        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(tmpdir)//"td: Starting from scratch"
          call write_warning(1)

          fromScratch = .true.
        end if
      end if

      if(.not.fromscratch .and. st%open_boundaries) then
        call restart_read_ob_intf(trim(restart_dir)//'gs', st, gr, ierr)
      end if

      if(.not. fromscratch .and. td%dynamics == CP) then 
        call cpmd_restart_read(td%cp_propagator, gr, st, ierr)

        if(ierr.ne.0) then
          message(1) = "Could not load "//trim(restart_dir)//"td/cpmd: Starting from scratch"
          call write_warning(1)
          
          fromScratch = .true.
          td%iter = 0
        end if
      end if

      if(.not. fromscratch) then 
        ! read potential from previous interactions
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'td/vprev_', i, is
            call dinput_function(trim(filename)//'.obf', gr%mesh, td%tr%v_old(1:NP, is, i), ierr)
            ! If we do not succeed, try netcdf
            if(ierr > 0) call dinput_function(trim(filename)//'.ncdf', gr%mesh, td%tr%v_old(1:NP, is, i), ierr)
            if(ierr > 0) then
              write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
              call write_fatal(1)
            end if
          end do
        end do
          
      end if

      if(fromScratch) then

        if(.not. st%only_userdef_istates) then
          ! In the open bounary case the ground state wavefunctions are bit too "wide".
          ! Therefore, we need a special routine to extract the middle.
          if(gr%sb%open_boundaries) then
            call restart_read_ob_intf(trim(restart_dir)//'gs', st, gr, ierr)
            if(ierr.ne.0) then
              message(1) = "Could not read interface wave functions from '"//trim(restart_dir)//"gs'"
              message(2) = "Please run an open boundaries ground-state calculation first!"
              call write_fatal(2)
            end if
            call restart_read_ob_central(trim(restart_dir)//'gs', st, gr, ierr)
          else
            call restart_read(trim(restart_dir)//'gs', st, gr, geo, ierr)
          end if
          if(ierr.ne.0) then
            message(1) = "Could not read KS orbitals from '"//trim(restart_dir)//"gs'"
            message(2) = "Please run a ground-state calculation first!"
            call write_fatal(2)
          end if
        end if

        ! check if we should deploy user-defined wavefunctions. 
        ! according to the settings in the input file the routine 
        ! overwrites orbitals that were read from restart/gs 
        if(loct_parse_isdef(datasets_check('UserDefinedStates')).ne.0) then
          call restart_read_user_def_orbitals(gr%mesh, st)
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
        !% Note that the code will not check the orthonormality of the new states!
        !%
        !% Each line provides the coefficients of the new states, in terms of
        !% the old ones.
        !%End
        if(loct_parse_isdef(datasets_check('TransformStates')).ne.0) then
          if(loct_parse_block(datasets_check('TransformStates'), blk) == 0) then
            call states_copy(stin, st)
            deallocate(stin%zpsi)
            call restart_look_and_read(stin, gr, sys%geo)
            ALLOCATE(rotation_matrix(st%nst, stin%nst), st%nst*stin%nst)
            rotation_matrix = M_z0
            do ist = 1, st%nst
              do jst = 1, loct_parse_block_cols(blk, ist-1)
                call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call rotate_states(gr%mesh, st, stin, rotation_matrix)
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

      !%Variable TDFreezeOrbitals
      !%Type integer
      !%Default 0
      !%Section Time Dependent
      !%Description
      !% You have the possibility of "freezing" a number of orbitals during a time-propagation.
      !% The Hartree and exchange-and-correlation potential originated by these orbitals (which
      !% will be the lowest-energy ones) will be added during the propagation, but the orbitals
      !% will not be propagated.
      !% 
      !% WARNING: NOT TESTED YET.
      !%Option sae -1
      !% Single-active-electron approximation. This option is only valid for time-dependent 
      !% calculations ("CalculationMode = td"). Also, the nuclei should not move. 
      !% The idea is that all orbitals except the last one are frozen. The orbitals are to 
      !% be read from a previous ground-state calculation. The active orbital is then treated
      !% as independent (no matter if it contains one electron or two) -- although it will
      !% feel the Hartree and exchange-correlation created by the ground-state electronic
      !% configureation.
      !% 
      !% It is almost equivalent to setting "TDFreezeOrbitals = N-1", where N is the number
      !% of orbitals, but not completely.
      !%End
      call loct_parse_int(datasets_check('TDFreezeOrbitals'), 0, freeze_orbitals)
      if(freeze_orbitals > 0) then
        ! In this case, we first freeze the orbitals, then calculate the Hxc potential.
        call states_freeze_orbitals(st, gr, sys%mc, freeze_orbitals)
        write(message(1),'(a,i4,a,i4,a)') 'Info: The lowest', freeze_orbitals, &
          ' orbitals have been frozen.', st%nst, ' will be propagated.'
        call write_info(1)
        call states_calc_dens(st, NP)
        call v_ks_calc(gr, sys%ks, hm, st, calc_eigenval=.true.)
      elseif(freeze_orbitals < 0) then
        ! This means SAE approximation. We calculate the Hxc first, then freezer all
        ! orbitals minus one.
        write(message(1),'(a)') 'Info: The single-active-electron approximation will be used.'
        call write_info(1)
        call states_calc_dens(st, NP)
        call v_ks_calc(gr, sys%ks, hm, st, calc_eigenval=.true.)
        call states_freeze_orbitals(st, gr, sys%mc, n = st%nst-1)
        call v_ks_freeze_hxc(sys%ks)
      else
        ! Normal run.
        call states_calc_dens(st, NP)
        call v_ks_calc(gr, sys%ks, hm, st, calc_eigenval=.true.)
      end if

      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      if(st%parallel_in_states) then
        call MPI_Bcast(x, 1, MPI_FLOAT, 0, st%mpi_grp%comm, mpi_err)
      end if
#endif
      call hamiltonian_span(hm, minval(gr%mesh%h(1:NDIM)), x)
      call total_energy(hm, gr, st, -1)

    end subroutine init_wfs


    ! ---------------------------------------------------------
    subroutine td_run_zero_iter()
      call push_sub('td.td_run_zero_iter')

      call td_write_iter(write_handler, gr, st, hm, geo, td%kick, td%dt, 0)

      ! I apply the delta electric field *after* td_write_iter, otherwise the
      ! dipole matrix elements in write_proj are wrong
      call apply_delta_field(td%kick)
      call td_rti_run_zero_iter(hm, td%tr)

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
              call mesh_r(gr%mesh, i, r, x = x)
              call loct_ylm(1, x(1), x(2), x(3), k%l(j), k%m(j), ylm)
              kick_function(i) = kick_function(i) + k%weight(j) * (r**k%l(j)) * ylm 
            end do
          end do
        else
          do i = 1, NP
            kick_function(i) = sum(gr%mesh%x(i, 1:NDIM)*k%pol(1:NDIM, k%pol_dir))
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
        if(ion_dynamics_ions_move(td%ions)  .and. k%delta_strength .ne. M_ZERO) then
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
      open(unit = iunit, file = io_workpath('td.general/coordinates'), &
        action='read', status='old', recl = record_length)
      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/coordinates'))//"'."
        message(2) = "Starting simulation from initial geometry."
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
      FLOAT :: vecpot(1:MAX_DIM), vecpot_vel(1:MAX_DIM), dummy(1:MAX_DIM)

      call push_sub('td.td_read_gauge_field')

      record_length = 28 + 3*3*20
      call io_assign(iunit)
      open(unit = iunit, file = io_workpath('td.general/gauge_field'), &
        action='read', status='old', recl = record_length)
      if(iunit < 0) then
        message(1) = "Could not open file '"//trim(io_workpath('td.general/gauge_field'))//"'."
        message(2) = "Starting simulation from initial values."
        call write_warning(2)
        return
      end if

      call io_skip_header(iunit)
      do i = 0, td%iter - 1
        read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(28x)', advance='no') ! skip the time index.

      ! TODO: units are missing
      read(iunit, '(3es20.12)', advance='no') vecpot(1:NDIM)
      read(iunit, '(3es20.12)', advance='no') vecpot_vel(1:NDIM)
      read(iunit, '(3es20.12)', advance='no') dummy(1:NDIM) ! skip the accel field.

      call gauge_field_set_vec_pot(hm%ep%gfield, vecpot)
      call gauge_field_set_vec_pot_vel(hm%ep%gfield, vecpot_vel)
 
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
        message(1) = 'Unsuccessful write of "'//trim(tmpdir)//'td"'
        call write_fatal(1)
      end if

      ! write potential from previous interactions
      if(mpi_grp_is_root(st%mpi_grp)) then
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a6,i2.2,i3.3)') 'vprev_', i, is
            call doutput_function(restart_format, trim(tmpdir)//"td", &
              filename, gr%mesh, gr%sb, td%tr%v_old(1:NP, is, i), M_ONE, ierr, is_tmp = .true.)
            if(ierr.ne.0) then
              write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
              call write_fatal(1)
            end if
          end do
        end do
      end if

      if(td%dynamics == CP) call cpmd_restart_write(td%cp_propagator, gr, st)

      call pop_sub()
    end subroutine td_save_restart

    ! ---------------------------------------------------------
    subroutine modify_occs()
      type(block_t) :: blk
      integer  :: nrow
      integer  :: spin, state
      FLOAT    :: new_occ
      
      if(loct_parse_block(datasets_check('ModifyOccupations'), blk) == 0) then
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
