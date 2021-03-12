!! Copyright (C) 2019-2021 N. Tancogne-Dejean
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

program tdtdm
  use batch_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use electrons_oct_m
  use fft_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use spectrum_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_restart_oct_m
  use symmetrizer_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use xc_oct_m

  implicit none

  integer :: in_file, ii, jj, kk, ierr, ip_h, irow, ifreq, nrow, it
  integer :: ik, ist, uist, istep, ikpoint, irep, out_file, iop, idim
  integer :: time_steps, energy_steps, istart, iend, ntiter, Nreplica, Ntrans
  FLOAT   :: dt, tt, weight, kpoint(MAX_DIM), xx(MAX_DIM)
  integer :: time_steps_ref, ref_file, irep_h
  FLOAT   :: dt_ref, start_time
  FLOAT, allocatable :: Et(:), ftreal(:, :, :), ftimag(:, :, :), tmp(:), omega(:)
  CMPLX, allocatable :: Xiak(:,:,:), Yiak(:,:,:)
  FLOAT, allocatable :: proj_r(:,:,:,:), proj_i(:,:,:,:)
  FLOAT, allocatable :: proj_r_corr(:,:), proj_i_corr(:,:), centers(:,:)
  CMPLX, allocatable :: tdm(:,:), tdm_1D(:,:,:,:)
  CMPLX, allocatable, target :: psi(:,:), upsi(:,:)
  CMPLX, allocatable ::  phase(:,:), ftcmplx(:,:)
  CMPLX, pointer :: psi_sym(:,:), upsi_sym(:,:)
  type(spectrum_t)  :: spectrum
  type(electrons_t), pointer :: sys
  type(batch_t)     :: projb_r, projb_i, ftrealb, ftimagb
  character(len=MAX_PATH_LEN) :: fname, ref_filename
  type(states_elec_t), pointer :: st
  type(states_elec_t) :: gs_st
  type(restart_t) :: restart
  type(unit_t) :: fn_unit
  integer :: kpt_start, kpt_end, supercell(1:MAX_DIM), nomega
  type(block_t) :: blk
  type(symmetrizer_t) :: symmetrizer
  FLOAT :: pos_h(MAX_DIM), ww

  ! Initializion
  call global_init()
  call parser_init()

  call messages_init()
  call io_init()

  call calc_mode_par_init()

  call profiling_init(global_namespace)

  call messages_experimental("oct-tdtdm utility")
  call fft_all_init(global_namespace)
  call unit_system_init(global_namespace)
  call restart_module_init(global_namespace)

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  sys => electrons_t(global_namespace)
  call sys%init_parallelization(mpi_world)

  call spectrum_init(spectrum, global_namespace)

  st => sys%st

  if(sys%st%d%ispin == SPINORS) then
    call messages_not_implemented('oct-tdtdm with spinors')
  end if

  if(st%parallel_in_states) then
   call messages_not_implemented("oct-tdtdm with states parallelization")
  end if

  if(sys%gr%mesh%parallel_in_domains) then
    call messages_not_implemented("oct-tdtdm with domain parallelization")
  end if

  !%Variable TDTDMFrequencies
  !%Type block
  !%Section Utilities::oct-tdtdm
  !%Description
  !% This block defines for which frequencies the analysis is performed.
  !%
  !% Each row of the block indicates a frequency.
  !%End
  if (parse_block(global_namespace, 'TDTDMFrequencies', blk) == 0) then

    nrow = parse_block_n(blk)
    nomega = nrow

    SAFE_ALLOCATE(omega(1:nrow))
    !read frequencies
    do irow = 0, nrow-1
      call parse_block_float(blk, irow, 0, omega(irow+1))
    end do

    call parse_block_end(blk)
  else
    message(1) = "oct-tdtdm: TDTDMFrequencies must be defined."
    call messages_fatal(1)
  end if

  ! We check that the resonant and antiresonant transitions are contained in the
  ! energy range of the Fourier transforms
  if(any(omega > spectrum%max_energy)) then
    message(1) = "One requested frequecy is larger than PropagationSpectrumMaxEnergy."
    message(2) = "Please increase the value of PropagationSpectrumMaxEnergy."
    call messages_fatal(2)
  end if
  if(any(omega > -spectrum%min_energy)) then
    message(1) = "One requested frequency is larger than -PropagationSpectrumMinEnergy."
    message(2) = "Please decrease the value of PropagationSpectrumMinEnergy."
    call messages_fatal(2)
  end if


  call states_elec_copy(gs_st, st, exclude_wfns = .true., exclude_eigenval = .true.)

  SAFE_DEALLOCATE_A(gs_st%node)

  call restart_init(restart, global_namespace, RESTART_PROJ, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh)
  if(ierr == 0) call states_elec_look(restart, ii, jj, gs_st%nst, ierr)
  if(ierr /= 0) then
    message(1) = "oct-tdtdm: Unable to read states information."
    call messages_fatal(1)
  end if

  ! allocate memory
  SAFE_ALLOCATE(gs_st%occ(1:gs_st%nst, 1:gs_st%d%nik))
  SAFE_ALLOCATE(gs_st%eigenval(1:gs_st%nst, 1:gs_st%d%nik))

  !We want all the task to have all the states
  !States can be distibuted for the states we propagate.
  SAFE_ALLOCATE(gs_st%node(1:gs_st%nst))
  gs_st%node(:)  = 0
  call kpoints_distribute(gs_st%d, sys%mc)
  call states_elec_distribute_nodes(gs_st, global_namespace, sys%mc)

  kpt_start = gs_st%d%kpt%start
  kpt_end = gs_st%d%kpt%end

  gs_st%eigenval = huge(gs_st%eigenval)
  gs_st%occ      = M_ZERO
  if(gs_st%d%ispin == SPINORS) then
    SAFE_DEALLOCATE_A(gs_st%spin)
    SAFE_ALLOCATE(gs_st%spin(1:3, 1:gs_st%nst, 1:gs_st%d%nik))
  end if

  call states_elec_allocate_wfns(gs_st, sys%gr%mesh, TYPE_CMPLX)
  call states_elec_load(restart, global_namespace, gs_st, sys%gr, sys%kpoints, ierr)
  if(ierr /= 0 .and. ierr /= (gs_st%st_end-gs_st%st_start+1)*(kpt_end-kpt_start+1)*gs_st%d%dim) then
    message(1) = "oct-tdtdm: Unable to read wavefunctions for TDOutput."
    call messages_fatal(1)
  end if
  call restart_end(restart)


  in_file = io_open('td.general/projections', action='read', status='old', die=.false.)
  if(in_file < 0) then 
    message(1) = "oct-tdtdm: Cannot open file '"//trim(io_workpath('td.general/projections'))//"'"
    call messages_fatal(1)
  end if
  call io_skip_header(in_file)
  call spectrum_count_time_steps(global_namespace, in_file, time_steps, dt)
  dt = units_to_atomic(units_out%time, dt)
  

  SAFE_ALLOCATE(tmp(1:st%nst*gs_st%nst*st%d%nik*2))
  SAFE_ALLOCATE(proj_r(1:time_steps, 1:gs_st%nst, 1:st%nst, 1:st%d%nik))
  SAFE_ALLOCATE(proj_i(1:time_steps, 1:gs_st%nst, 1:st%nst, 1:st%d%nik))

  
  call io_skip_header(in_file)
  
  do ii = 1, time_steps
    read(in_file, *) jj, tt, (tmp(kk), kk = 1, st%nst*gs_st%nst*st%d%nik*2)
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do uist = 1, gs_st%nst
          jj  = (ik-1)*st%nst*gs_st%nst + (ist-1)*gs_st%nst + uist
          proj_r(ii, uist, ist, ik) = tmp((jj-1)*2+1)
          ! Here we add a minus sign, as we want to get <\phi_0 | \psi(t)>
          ! and td_occup computes the complex conjugaute of this
          proj_i(ii, uist, ist, ik) = -tmp((jj-1)*2+2)
        end do
      end do
    end do
  end do
  SAFE_DEALLOCATE_A(tmp)

  call io_close(in_file)

  write(message(1), '(a, i7, a)') "oct-tdtdm: Read ", time_steps, " steps from file '"// &
    trim(io_workpath('td.general/projections'))//"'"
  call messages_info(1)

  start_time = spectrum%start_time

  !%Variable TransientPopulationDelay
  !%Type float
  !%Default 0.
  !%Section Utilities::oct-tdpop2exc
  !%Description
  !% This variable defines the starting time used to compute the Fourier transform of the 
  !% populations. This is needed if TransientPopulationReference is used.
  !%End
  call parse_variable(global_namespace, 'TransientPopulationDelay', start_time, spectrum%start_time )

  if(parse_is_defined(global_namespace, 'TransientPopulationReference') .and. spectrum%start_time > 0) then
    !%Variable TransientPopulationReference
    !%Type string
    !%Default "."
    !%Section Utilities::oct-tdpop2exc
    !%Description
    !% In case of delayed kick, the calculation of the transient populations requires 
    !% to substract a reference calculation, containing the output of td_occup without the kick
    !% This variables defined the directory in which the reference projections file is,
    !% relative to the current folder
    !%End

    call parse_variable(global_namespace, 'TransientPopulationReference', '.', ref_filename)
    ref_file = io_open(trim(ref_filename)//'/projections', action='read', status='old', die=.false.)
    if(ref_file < 0) then
      message(1) = "Cannot open reference file '"//trim(io_workpath(trim(ref_filename)//'/projections'))//"'"
      call messages_fatal(1)
    end if
    call io_skip_header(ref_file)
    call spectrum_count_time_steps(global_namespace, ref_file, time_steps_ref, dt_ref)
    dt_ref = units_to_atomic(units_out%time, dt_ref)
    time_steps_ref = time_steps_ref + 1
    if(time_steps_ref < time_steps) then
      message(1) = "The reference calculation does not contain enought time steps"
      call messages_fatal(1)
    end if

    if(dt_ref /= dt) then
      message(1) = "The time step of the reference calculation is different from the current calculation"
      call messages_fatal(1)
    end if

!    !We remove the reference
!    !TODO: we need to count the column here, to check that enough states are propagated
!    SAFE _ALLOCATE(tmp(1:gs_st%nst*gs_st%nst*st%d%nik*2))
!    call io_skip_header(ref_file)
!
!    SAFE _ALLOCATE(proj_r_ref(1:time_steps, 1:gs_st%nst, 1:gs_st%nst, 1:st%d%nik))
!    SAFE _ALLOCATE(proj_i_ref(1:time_steps, 1:gs_st%nst, 1:gs_st%nst, 1:st%d%nik))
!    do ii = 1, time_steps
!      read(ref_file, *) jj, tt, (tmp(kk), kk = 1, gs_st%nst*gs_st%nst*st%d%nik*2)
!      do ik = 1, st%d%nik
!        do ist = 1, gs_st%nst
!          do uist = 1, gs_st%nst
!            jj  = (ik-1)*gs_st%nst*gs_st%nst + (ist-1)*gs_st%nst + uist
!            proj_r_ref(ii, uist, ist, ik) = tmp((jj-1)*2+1)
!            proj_i_ref(ii, uist, ist, ik) = tmp((jj-1)*2+2)
!         end do
!       end do
!     end do
!    end do
!    SAFE _DEALLOCATE_A(tmp)
!    call io_close(ref_file)
!
!    !Here we use the reference calculation and we insert of closure relation
!    !This is done to avoid propagating unoccupied states
!    do ik = 1, st%d%nik
!      do ist = 1, st%nst
!        do uist = 1, gs_st%nst
!          do uist_ref = 1, gs_st%nst ! Inner sum from the closure relation
!            do ii = 1, time_steps
!              proj_r_corr(ii, uist, ist, ik) = proj_r_corr(ii, uist, ist, ik) + proj_r(ii, uist_ref, ist, ik)*proj_r_ref(ii, uist_ref, uist, ik) + proj_i(ii, uist_ref, ist, ik)*proj_i_ref(ii, uist_ref, uist, ik)
!              proj_i_corr(ii, uist, ist, ik) = proj_i_corr(ii, uist, ist, ik) - proj_r(ii, uist_ref, ist, ik)*proj_i_ref(ii, uist_ref, uist, ik) + proj_i(ii, uist_ref, ist, ik)*proj_r_ref(ii, uist_ref, uist, ik)
!            end do
!          end do
!        end do
!      end do
!    end do
!
!    !We now transform the GS states into the time-evolved states at the time of the kick
!    !Doing so, the code constructing the two-particle wavefunction remains identical to the
!    !equilibrium case
!    time_kick = int(spectrum%start_time / dt)
!    SAFE _ALLOCATE(psi(1:sys%gr%mesh%np, 1:gs_st%d%dim))
!    SAFE _ALLOCATE(upsi(1:sys%gr%mesh%np, 1:gs_st%d%dim))
!
!    call states_copy(ref_st, gs_st)
!
!    do ik = 1, st%d%nik
!      do ist = 1, gs_st%nst
!        psi = M_ZERO
!        do uist = 1, gs_st%nst
!          call states_get_state(ref_st, sys%gr%mesh, uist, ik, upsi)
!          call lalg_axpy(sys%gr%mesh%np, TOCMPLX(proj_r_ref(time_kick, uist, ist, ik), proj_i_ref(time_kick, uist, ist, ik)), upsi(:, 1), psi(:, 1))
!        end do
!        call states_set_state(gs_st, sys%gr%mesh, ist, ik, psi)
!      end do
!    end do
!
!    call states_end(ref_st)
!
!    SAFE _DEALLOCATE_A(psi)
!    SAFE _DEALLOCATE_A(upsi)
!
 !   SAFE _DEALLOCATE_A(proj_r_ref)
 !   SAFE _DEALLOCATE_A(proj_i_ref)

  end if


  ! Phase correction of the projections before doing the Fourier transforms
  ! See Eq. (5) of Williams et al., JCTC 17, 1795 (2021)
  ! We need to multiply C_ik(t)e^{-ie_kt} (the projection of \phi_i(t) on \phi_k^GS) 
  ! by e^{ie_it}, which is obtained by the cc of the projection of \phi_i(t) on \phi_i^GS
  ! Here we only care about optical transitions (so TD occupied to GS unocc)
  SAFE_ALLOCATE(proj_r_corr(1:time_steps, 1:gs_st%nst*st%nst*(kpt_end-kpt_start+1)))
  SAFE_ALLOCATE(proj_i_corr(1:time_steps, 1:gs_st%nst*st%nst*(kpt_end-kpt_start+1)))
  proj_r_corr = M_ZERO
  proj_i_corr = M_ZERO
  do ik = kpt_start, kpt_end
    do ist = 1, st%nst
      do uist = ist+1, gs_st%nst
        jj = (ik-kpt_start)*st%nst*gs_st%nst+(ist-1)*gs_st%nst+uist
        do ii = 1, time_steps
          proj_r_corr(ii, jj) = proj_r(ii, uist, ist, ik) * proj_r(ii, ist, ist, ik) + proj_i(ii, uist, ist, ik) * proj_i(ii, ist, ist, ik)
          proj_i_corr(ii, jj) =-proj_r(ii, uist, ist, ik) * proj_i(ii, ist, ist, ik) + proj_i(ii, uist, ist, ik) * proj_r(ii, ist, ist, ik)
        end do
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(proj_r)
  SAFE_DEALLOCATE_A(proj_i)

  ! Find out the iteration numbers corresponding to the time limits.
  call spectrum_fix_time_limits(spectrum, time_steps, dt, istart, iend, ntiter)
  istart = max(1, istart)
  energy_steps = spectrum_nenergy_steps(spectrum)

  SAFE_ALLOCATE(ftreal(1:energy_steps, 1:st%nst*gs_st%nst*(kpt_end-kpt_start+1), 1:2))
  SAFE_ALLOCATE(ftimag(1:energy_steps, 1:st%nst*gs_st%nst*(kpt_end-kpt_start+1), 1:2))

  call batch_init(projb_r, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), proj_r_corr)
  call batch_init(projb_i, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), proj_i_corr)
  call batch_init(ftrealb, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), ftreal(:,:,1))
  call batch_init(ftimagb, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), ftimag(:,:,1))

  write(message(1), '(a)') "oct-tdtdm: Fourier transforming real part of the projections"
  call messages_info(1)

  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, projb_r)
  call spectrum_signal_damp(spectrum%damp, spectrum%damp_factor, istart, iend, spectrum%start_time, dt, projb_i)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, projb_r, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, projb_r, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftimagb)

  call ftrealb%end()
  call ftimagb%end()

  SAFE_ALLOCATE(ftcmplx(1:energy_steps, 1:st%nst*gs_st%nst*(kpt_end-kpt_start+1)))
  do ii = 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1)
    ftcmplx(1:energy_steps,ii) =  ftreal(1:energy_steps,ii,1) + M_zI*ftimag(1:energy_steps,ii,1) 
  end do
 
  write(message(1), '(a)') "oct-tdtdm: Fourier transforming imaginary part of the projections"
  call messages_info(1)

  call batch_init(ftrealb, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), ftreal(:,:,2))
  call batch_init(ftimagb, 1, 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1), ftimag(:,:,2))

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_COS, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, projb_i, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftrealb)

  call spectrum_fourier_transform(spectrum%method, SPECTRUM_TRANSFORM_SIN, spectrum%noise, &
    istart, iend, spectrum%start_time, dt, projb_i, spectrum%min_energy, spectrum%max_energy, spectrum%energy_step, ftimagb)

  call projb_i%end()
  call projb_r%end()
  call ftrealb%end()
  call ftimagb%end()
  SAFE_DEALLOCATE_A(proj_r_corr)
  SAFE_DEALLOCATE_A(proj_i_corr)

  do ii = 1, st%nst*gs_st%nst*(kpt_end-kpt_start+1)
    ftcmplx(1:energy_steps,ii) = ftcmplx(1:energy_steps,ii) + M_zI*ftreal(1:energy_steps,ii,2) - ftimag(1:energy_steps,ii,2)
  end do

  out_file = io_open('td.general/tdm_ftcmplx', action='write')
  write(out_file, '(a)') '# Energy  Re[TF(Cai(t)] Im[TF(Cai(t)]'
  do it = 1, energy_steps
    ww = (it-1)*spectrum%energy_step + spectrum%min_energy
    write(out_file, '(9999e15.6)') , ww, (TOFLOAT(ftcmplx(it, ii)), aimag(ftcmplx(it, ii)), ii=1, st%nst*gs_st%nst*(kpt_end-kpt_start+1))
  end do
  call io_close(out_file)


  SAFE_DEALLOCATE_A(ftreal)
  SAFE_DEALLOCATE_A(ftimag)

  write(message(1), '(a)') "oct-tdtdm: Constructing the two-particle wavefunctions."
  call messages_info(1)

  !At the moment the supercell is givn by the number of k-points.
  supercell(1:3) = sys%kpoints%nik_axis(1:3)
  Nreplica = product(supercell(1:sys%space%dim))

  ! The center of each replica of the unit cell
  SAFE_ALLOCATE(centers(1:3, 1:Nreplica))
  irep = 1
  do ii = 0, supercell(1)-1
    do jj = 0, supercell(2)-1
      do kk = 0, supercell(3)-1
        centers(1, irep) = -floor((supercell(1)-1)/M_TWO)+ii
        centers(2, irep) = -floor((supercell(2)-1)/M_TWO)+jj
        centers(3, irep) = -floor((supercell(3)-1)/M_TWO)+kk
        centers(1:sys%space%dim, irep) = matmul(sys%gr%sb%latt%rlattice(1:sys%space%dim, 1:sys%space%dim), centers(1:sys%space%dim, irep))
        irep = irep + 1
      end do
    end do
  end do

  ! The phase for each center
  SAFE_ALLOCATE(phase(kpt_start:kpt_end, 1:Nreplica))

  !!In case a vector potential is present, we need to add it.
  !vector_potential(1:sys%gr%mesh%sb%dim) = M_ZERO
  !call laser_init(lasers, parser, n_lasers, sys%gr%mesh)
  !do ilaser = 1, n_lasers
  !  select case(laser_kind(lasers(ilaser)))
  !  case(E_FIELD_VECTOR_POTENTIAL)
  !    aa = M_ZERO
  !    call laser_field(lasers(ilaser), aa(1:sys%gr%mesh%sb%dim), spectrum%start_time)
  !    vector_potential(1:sys%gr%mesh%sb%dim) = vector_potential(1:sys%gr%mesh%sb%dim) - aa(1:sys%gr%mesh%sb%dim)/P_C
  !  end select
  !end do
  !call laser_end(n_lasers, lasers)

  do irep = 1, Nreplica 
    do ik = kpt_start, kpt_end
      ikpoint = gs_st%d%get_kpoint_index(ik)
      kpoint(1:sys%space%dim) = sys%kpoints%get_point(ikpoint)
  !    !We add the vector potential
  !    kpoint(1:sys%space%dim) = kpoint(1:sys%space%dim) + vector_potential(1:sys%space%dim)
      phase(ik, irep) = exp(-M_zI*sum(kpoint(1:sys%space%dim)*centers(1:sys%space%dim, irep)))
    end do
  end do

  !!We need to update the Hamiltonian, such as the phases are correct, in case of a laser field
  !call hamiltonian_update(sys%hm, sys%gr%mesh, sys%gr%der%boundaries, spectrum%start_time)


  ! Position of the hole, here assumed to be on top of the first atom
  ! To be obtained from the input file
  if(sys%gr%mesh%sb%dim > 1) then
    call tdtdm_get_hole_position(pos_h, ip_h)
  end if

  Ntrans = 0
  ! Here we assume that there is a clear gap, so the information at Gamma is enough
  do ist = 1, gs_st%nst
    if(abs(gs_st%occ(ist, 1)) < M_EPSILON) cycle

    do uist = 1, gs_st%nst
      if(abs(gs_st%occ(uist, 1)) > M_EPSILON) cycle
      weight = gs_st%d%kweights(1) * (gs_st%occ(ist, 1)-gs_st%occ(uist, 1))
      if(abs(weight) < M_EPSILON) cycle
      Ntrans = Ntrans + 1
    end do
  end do

  SAFE_ALLOCATE(Xiak(1:st%nst, 1:gs_st%nst, 1:st%d%nik))
  SAFE_ALLOCATE(Yiak(1:st%nst, 1:gs_st%nst, 1:st%d%nik))
  SAFE_ALLOCATE(Et(1:Ntrans*st%d%nik))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:gs_st%d%dim))
  SAFE_ALLOCATE(upsi(1:sys%gr%mesh%np, 1:gs_st%d%dim))

  if(sys%kpoints%use_symmetries) then
    call symmetrizer_init(symmetrizer, sys%gr%mesh, sys%gr%symm)
    SAFE_ALLOCATE(psi_sym(1:sys%gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(upsi_sym(1:sys%gr%mesh%np, 1:st%d%dim))
  end if

  select case(sys%gr%mesh%sb%dim)
  case(2,3)
    SAFE_ALLOCATE(tdm(1:sys%gr%mesh%np, 1:Nreplica))
  case(1)
    SAFE_ALLOCATE(tdm_1D(1:sys%gr%mesh%np, 1:sys%gr%mesh%np, 1:Nreplica, 1:Nreplica))
  end select

  do ifreq = 1, nomega

    write(message(1), '(a, f6.4, a)') "oct-tdtdm: Constructing the two-particle wavefunction at ", omega(ifreq), " Ha."
    call messages_info(1)

    select case(sys%gr%mesh%sb%dim)
    case(2,3)
      tdm = M_z0
    case(1)
      tdm_1D = M_z0
    end select
 
    Et = M_ZERO
    Xiak = M_z0
    Yiak = M_z0

    !Local transition index
    it = (kpt_start-1)*Ntrans + 1

    do ik = kpt_start, kpt_end
      ikpoint = st%d%get_kpoint_index(ik)

      do ist = 1, st%nst
        if(abs(gs_st%occ(ist, ik)) < M_EPSILON) cycle

        call states_elec_get_state(gs_st, sys%gr%mesh, ist, ik, psi)
        if(allocated(sys%hm%hm_base%phase)) then
          call states_elec_set_phase(gs_st%d, psi, sys%hm%hm_base%phase(1:sys%gr%mesh%np, ik), sys%gr%mesh%np, .false.)
        end if

        do uist = 1, gs_st%nst
          if(abs(gs_st%occ(uist, ik)) > M_EPSILON) cycle
  
          weight = gs_st%d%kweights(ik) * (gs_st%occ(ist, ik)-gs_st%occ(uist, ik)) &
                      / kpoints_get_num_symmetry_ops(sys%kpoints, ikpoint)
          if(abs(weight) < M_EPSILON) cycle

          call states_elec_get_state(gs_st, sys%gr%mesh, uist, ik, upsi)
          if(allocated(sys%hm%hm_base%phase)) then
            call states_elec_set_phase(gs_st%d, upsi, sys%hm%hm_base%phase(1:sys%gr%mesh%np, ik), sys%gr%mesh%np, .false.)
          end if

          do ii = 1, kpoints_get_num_symmetry_ops(sys%kpoints, ikpoint)
            iop = kpoints_get_symmetry_ops(sys%kpoints, ikpoint, ii)

            if(sys%kpoints%use_symmetries) then
              do idim = 1, st%d%dim
                call zsymmetrizer_apply_single(symmetrizer, sys%gr%mesh, iop, psi(:,idim), psi_sym(:,idim))
                call zsymmetrizer_apply_single(symmetrizer, sys%gr%mesh, iop, upsi(:,idim), upsi_sym(:,idim))
              end do
            else
              psi_sym => psi
              upsi_sym => upsi
            end if

            ! For a given requested frequency, we get the corresponding values of Xia and Yia
            ! One correspond to the +\Omega frequency, the other one to the -\Omega frequency
            ! For Xiak, we use the fact that TF[f*](\Omega) = (TF[f](-\Omega))^*
            jj = (ik-kpt_start)*st%nst*gs_st%nst+(ist-1)*gs_st%nst+uist
            istep = int((-omega(ifreq)-spectrum%min_energy)/spectrum%energy_step)
            Xiak(ist, uist, ik) = conjg(ftcmplx(istep, jj))
            istep = int((-omega(ifreq)-spectrum%min_energy)/spectrum%energy_step)
            Yiak(ist, uist, ik) = ftcmplx(istep, jj)

            ! We now compute the single mode TDTDM
            ! See Eq. (5) of Williams et al., JCTC 17, 1795 (2021)
            ! We take here the complex conjugate of the 2-body wavefunction
            select case(sys%gr%mesh%sb%dim)
            case(2,3)
              do irep = 1, Nreplica 
                call lalg_axpy(sys%gr%mesh%np, phase(ik, irep)*weight*conjg(Xiak(ist,uist,ik))*conjg(psi(ip_h,1)), upsi(:, 1), tdm(:,irep))
                call lalg_axpy(sys%gr%mesh%np, phase(ik, irep)*weight*conjg(Yiak(ist,uist,ik))*conjg(upsi(ip_h,1)), psi(:, 1), tdm(:,irep))
              end do
            case(1)
              ! In the 1D case, we contruct the full TDTDM of r_e, r_h
              do irep_h = 1, Nreplica
                do irep = 1, Nreplica
                  do ip_h = 1, sys%gr%mesh%np
                    call lalg_axpy(sys%gr%mesh%np, phase(ik, irep)*conjg(phase(ik, irep_h))*weight*conjg(Xiak(ist,uist,ik))*conjg(psi_sym(ip_h,1)), upsi(:, 1), tdm_1D(:, ip_h, irep, irep_h))
                    call lalg_axpy(sys%gr%mesh%np, phase(ik, irep)*conjg(phase(ik, irep_h))*weight*conjg(Yiak(ist,uist,ik))*conjg(upsi_sym(ip_h,1)), psi(:, 1), tdm_1D(:, ip_h, irep, irep_h))
                  end do
                end do
              end do
            end select
 
          end do ! ii

          Et(it) = gs_st%eigenval(uist, ik) - gs_st%eigenval(ist, ik)
          it = it + 1
        end do
      end do
    end do

#if defined(HAVE_MPI)        
    if(gs_st%d%kpt%parallel) then
      if(sys%gr%mesh%sb%dim > 1) then
        call comm_allreduce(gs_st%d%kpt%mpi_grp, tdm)
      else
        call comm_allreduce(gs_st%d%kpt%mpi_grp, tdm_1D)
      end if
      call comm_allreduce(gs_st%d%kpt%mpi_grp, Et)
      call comm_allreduce(gs_st%d%kpt%mpi_grp, Xiak)
      call comm_allreduce(gs_st%d%kpt%mpi_grp, Yiak)
    end if
#endif  

    call tdtdm_output_density()

    call tdtdm_excitonic_weight()

    !call tdtdm_dielectric_function()

    !write(fname, '(a, f0.4)') 'tdm_trans-0', omega(ifreq) 
    !out_file = io_open('td.general/'//trim(fname), action='write')
    !write(out_file, '(a)') '# Energy  |Xai(k)|^2 |Yai(k)|^2'
    !do it = 1, Ntrans*gs_st%d%nik
    !  write(out_file, '(3e15.6)') Et(it), TOFLOAT(Xiak(it)*conjg(Xiak(it))), &
    !                            TOFLOAT(Yiak(it)*conjg(Yiak(it)))
    !end do
    !call io_close(out_file)

  end do ! ifreq

  SAFE_DEALLOCATE_A(Et)
  SAFE_DEALLOCATE_A(Xiak)
  SAFE_DEALLOCATE_A(Yiak)
  SAFE_DEALLOCATE_A(tdm)
  SAFE_DEALLOCATE_A(tdm_1D)

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(upsi)
  if(sys%kpoints%use_symmetries) then 
    SAFE_DEALLOCATE_P(psi_sym)
    SAFE_DEALLOCATE_P(upsi_sym)
  end if
  SAFE_DEALLOCATE_A(ftcmplx)
  SAFE_DEALLOCATE_A(centers)
  SAFE_DEALLOCATE_A(phase)
  SAFE_DEALLOCATE_A(omega)

  SAFE_DEALLOCATE_P(sys) 
  call states_elec_end(gs_st)
  call fft_all_end()
  call io_end()
  call profiling_end(global_namespace)
  call messages_end()
  call parser_end()
  call global_end()

  contains

    ! -----------------------------------------------------------------
    ! Determines the position of the hole, either from the input or using the
    ! first atom in the cell.
    ! This returns the index of the point in the mesh closest to the position.
    subroutine tdtdm_get_hole_position(xx_h, ip_h)
      FLOAT,   intent(out) :: xx_h(1:MAX_DIM)
      integer, intent(out) :: ip_h

      FLOAT :: dmin
      integer :: idir, rankmin

      PUSH_SUB(tdtdm_get_hole_position)

      !%Variable TDTDMHoleCoordinates
      !%Type float
      !%Section Utilities::oct-tdtdm
      !%Description
      !% The position of the hole used to compute the TDTDM,
      !% in Cartesian coordinates.
      !% Note that the code will use the closest grid point.
      !%
      !% The coordinates of the hole are specified in the following way
      !% <tt>%TDTDMHoleCoordinates
      !% <br>&nbsp;&nbsp;hole_x | hole_y | hole_z
      !% <br>%</tt>
      !% 
      !% If not specified, the code will use the coordinate of the first atom in the cell.
      !%End

      if(parse_block(global_namespace, 'TDTDMHoleCoordinates', blk) == 0) then
        if(parse_block_cols(blk,0) < sys%space%dim) then
          call messages_input_error(global_namespace, 'TDTDMHoleCoordinates')
        end if
        do idir = 1, sys%space%dim
          call parse_block_float(blk, 0, idir - 1, xx_h(idir), units_inp%length)
        end do
        call parse_block_end(blk) 
      else
        xx_h(1:sys%space%dim) = sys%geo%atom(1)%x(1:sys%space%dim)
      end if

      write(message(1), '(a, 3(1x,f7.4,a))') "oct-tdtdm: Setting the hole at (", xx_h(1), &
                     ",", xx_h(2), ",", xx_h(3), ")."
      call messages_info(1)
 
      !At the moment, we ignore rankmin
      ASSERT(.not.sys%gr%mesh%parallel_in_domains)
      ip_h = mesh_nearest_point(sys%gr%mesh, xx_h, dmin, rankmin)
     ! xx_h = matmul(xx_h(1:3), sys%gr%sb%latt%klattice)/(M_TWO*M_PI)
     ! xx_h(1:3) = xx_h(1:3) * M_TWO*sys%gr%sb%lsize(1:3) / sys%gr%mesh%spacing(1:3)
     ! ip_h = sys%gr%mesh%idx%lxyz_inv(nint(xx_h(1)), nint(xx_h(2)), nint(xx_h(3)))

      POP_SUB(tdtdm_get_hole_position)
    end subroutine tdtdm_get_hole_position

    subroutine tdtdm_output_density()
      FLOAT, allocatable :: den(:,:), den_1D(:,:,:,:)
      FLOAT :: norm, xx(3), xx_h(3)
      integer :: iunit

      PUSH_SUB(tdtdm_output_density)
  
      ! We compute the TDM density 
      select case(sys%gr%mesh%sb%dim)
      case(2,3)
        SAFE_ALLOCATE(den(1:sys%gr%mesh%np, 1:Nreplica))
        do irep = 1, Nreplica
          do ii = 1, sys%gr%mesh%np
            den(ii, irep) = TOFLOAT(tdm(ii, irep)*conjg(tdm(ii, irep)))
          end do
        end do

        ! Here we renormalize to avoid too small numbers in the outputs
        norm = maxval(den)
        call lalg_scal(sys%gr%mesh%np, Nreplica, M_ONE/norm, den)

      case(1)
        SAFE_ALLOCATE(den_1D(1:sys%gr%mesh%np, 1:sys%gr%mesh%np, 1:Nreplica, 1:Nreplica))
        do irep_h = 1, Nreplica
          do irep = 1, Nreplica
            do ip_h = 1, sys%gr%mesh%np
              do ii = 1, sys%gr%mesh%np
                tdm_1D(ii, ip_h, irep, irep_h) = conjg(tdm_1D(ii, ip_h, irep, irep_h))
                den_1D(ii, ip_h, irep, irep_h) = TOFLOAT(tdm_1D(ii, ip_h, irep, irep_h)*conjg(tdm_1D(ii,ip_h, irep, irep_h)))
              end do
            end do
          end do
        end do
      end select
       
      fn_unit = units_out%length**(-sys%gr%mesh%sb%dim)

      select case(sys%gr%mesh%sb%dim)
      case(2,3)
        write(fname, '(a, f0.4)') 'tdm_density-0', omega(ifreq)
        call dio_function_output_supercell(io_function_fill_how("XCrySDen"), "td.general", fname, sys%gr%mesh, &
          den, centers, supercell, fn_unit, ierr, global_namespace, geo = sys%geo, grp = st%dom_st_kpt_mpi_grp, extra_atom=pos_h)
 
        call dio_function_output_supercell(io_function_fill_how("PlaneZ"), "td.general", fname, sys%gr%mesh, &
          den, centers, supercell, fn_unit, ierr, global_namespace, geo = sys%geo, grp = st%dom_st_kpt_mpi_grp)

        SAFE_DEALLOCATE_A(den)

      case(1)

        call tdtdm_get_hole_position(pos_h, ip_h)

        write(fname, '(a, f0.4)') 'tdm_density-0', omega(ifreq)
        call dio_function_output_supercell(io_function_fill_how("AxisX"), "td.general", fname, sys%gr%mesh, &
            den_1D(:,ip_h,:,(supercell(1)-1)/2), centers, supercell, fn_unit, ierr, global_namespace, geo = sys%geo, grp = st%dom_st_kpt_mpi_grp)

        write(fname, '(a, f0.4)') 'tdm_wfn-0', omega(ifreq)
        call zio_function_output_supercell(io_function_fill_how("AxisX"), "td.general", fname, sys%gr%mesh, &
            tdm_1D(:,ip_h,:,(supercell(1)-1)/2), centers, supercell, fn_unit, ierr, global_namespace, geo = sys%geo, grp = st%dom_st_kpt_mpi_grp)

        ASSERT(.not.sys%gr%mesh%parallel_in_domains)
        write(fname, '(a, f0.4)') 'td.general/tdm_density-0', omega(ifreq)
        iunit = io_open(fname, action='write')
        write(iunit, '(a)', iostat=ierr) '# r_e    r_h    Re(\Psi(r_e,r_h)) Im(\Psi(r_e,r_h)) |\Psi(r_e,r_h)|^2'

        do irep_h = 1, Nreplica
          do ip_h = 1, sys%gr%mesh%np
            xx_h = units_from_atomic(units_out%length, mesh_x_global(sys%gr%mesh, ip_h) + centers(1:sys%gr%mesh%sb%dim, irep_h))

            do irep = 1, Nreplica
              do ii = 1, sys%gr%mesh%np
                xx = units_from_atomic(units_out%length, mesh_x_global(sys%gr%mesh, ii) + centers(1:sys%gr%mesh%sb%dim, irep))
                write(iunit, '(5es23.14E3)', iostat=ierr) xx(1), xx_h(1), &
                              TOFLOAT(units_from_atomic(fn_unit, tdm_1D(ii, ip_h, irep, irep_h))) ,&
                              aimag(units_from_atomic(fn_unit, tdm_1D(ii, ip_h, irep, irep_h))), &
                              units_from_atomic(fn_unit, den_1D(ii, ip_h, irep, irep_h))
              end do
            end do
          end do
        end do

        SAFE_DEALLOCATE_A(den_1D)
      end select
      

      POP_SUB(tdtdm_output_density)
    end subroutine tdtdm_output_density

    subroutine tdtdm_excitonic_weight()
      FLOAT, allocatable :: weight(:,:)

      PUSH_SUB(tdtdm_excitonic_weight)

      SAFE_ALLOCATE(weight(1:st%d%nik, 1:gs_st%nst))
      weight = M_ZERO

      do ik = 1, st%d%nik
        do ist = 1, st%nst
          if(abs(gs_st%occ(ist, ik)) < M_EPSILON) cycle

          do uist = ist+1, gs_st%nst
            if(abs(gs_st%occ(uist, ik)) > M_EPSILON) cycle
            
            weight(ik, ist)  = weight(ik, ist) + abs(Xiak(ist, uist, ik))**2
            weight(ik, uist) = weight(ik,uist) + abs(Yiak(ist, uist, ik))**2
          end do
        end do
      end do

      write(fname, '(a, f0.4)') 'td.general/tdm_weights-0', omega(ifreq)
      out_file = io_open(fname, action='write')
      write(out_file, '(a)') '# ik kx ky kz weights(ist,ik) '
      do ik = 1, st%d%nik 
        ikpoint = st%d%get_kpoint_index(ik)
!        kpoint(1:sys%space%dim) = sys%kpoints%get_point(ikpoint, absolute_coordinates=.false.)
        kpoint(1:sys%space%dim) = sys%kpoints%reduced%point1BZ(1:sys%space%dim,ikpoint)
        write(out_file, '(i4,3e15.6)', advance='no') ik, kpoint(1:3)
        do uist = 1, gs_st%nst-1
          write(out_file, '(e15.6)', advance='no') weight(ik, uist) 
        end do
        write(out_file, '(e15.6)') weight(ik, uist)
      end do
      call io_close(out_file) 

      SAFE_DEALLOCATE_A(weight)

      POP_SUB(tdtdm_excitonic_weight)
    end subroutine tdtdm_excitonic_weight

    ! Because we know the X and Y vectors of Casida, we can get the dielectric function
    ! This can be useful to compare with the real-time TDDFT results
    ! See Eq. 48 of J. Chem. Phys. 146, 064110 (2017)
!    subroutine tdtdm_dielectric_function()
!      CMPLX, allocatable :: eps(:)
!      integer :: iw, nomega
!      FLOAT :: domega
!      CMPLX :: omega
!
!      PUSH_SUB(tdtdm_dielectric_function)
! 
!      nomega = 101
!      domega =CNST(0.1/27.2114)
!
!      SAFE_ALLOCATE(eps(1:Nomega))      
!
!      do iw = 1, nomega
!        eps(iw) = M_Z0
!        omega = (iw-1)*domega + M_zI*CNST(0.1/27.2114)
!        !Loop over optical transitions
!        do it = 1, Ntrans*gs_st%d%nik
!          eps(iw) = (M_ONE/(omega - Et(it)) - M_ONE/(omega + Et(it)))*(Xiak(it) + Yiak(it))*conjg(Xiak(it) + Yiak(it))
!        end do
!      end do
!
!      eps = M_ONE - M_FOUR*M_PI/sys%gr%sb%latt%rcell_volume * eps
!
!      out_file = io_open('td.general/tdm_dielectric_function', action='write')
!      write(out_file, '(a)') '# Energy  Re(eps) Im(eps)'
!      do iw = 1, Nomega
!        write(out_file, '(3e15.6)') (iw-1)*domega, TOFLOAT(eps(iw)), aimag(eps(iw))
!      end do
!      call io_close(out_file)
! 
!
 !     SAFE_DEALLOCATE_A(eps)
 !    
 !     POP_SUB(tdtdm_dielectric_function)
 !   end subroutine tdtdm_dielectric_function

end program tdtdm

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
