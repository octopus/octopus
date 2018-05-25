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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module hamiltonian_oct_m
  use accel_oct_m
  use base_hamiltonian_oct_m
  use base_potential_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use boundaries_oct_m
  use boundary_op_oct_m
  use cmplxscl_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use hamiltonian_base_oct_m
  use epot_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hardware_oct_m
  use io_oct_m
  use io_function_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use oct_exchange_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use pcm_oct_m
  use restart_oct_m
  use scdm_oct_m
  use scissor_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_parallel_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use xc_functl_oct_m
  use XC_F90(lib_m)

  implicit none

  private
  public ::                          &
    hamiltonian_t,                   &
    hamiltonian_init,                &
    hamiltonian_end,                 &
    hamiltonian_span,                &
    dhamiltonian_apply,              &
    zhamiltonian_apply,              &
    dhamiltonian_apply_all,          &
    zhamiltonian_apply_all,          &
    dhamiltonian_apply_batch,        &
    zhamiltonian_apply_batch,        &
    dhamiltonian_diagonal,           &
    zhamiltonian_diagonal,           &
    dmagnus,                         &
    zmagnus,                         &
    dvmask,                          &
    zvmask,                          &
    hamiltonian_inh_term,            &
    hamiltonian_set_inh,             &
    hamiltonian_remove_inh,          &
    hamiltonian_adjoint,             &
    hamiltonian_not_adjoint,         &
    hamiltonian_hermitian,           &
    hamiltonian_epot_generate,       &
    hamiltonian_update,              &
    hamiltonian_update2,             &
    hamiltonian_get_time,            &
    hamiltonian_apply_packed,        &
    dexchange_operator_single,       &
    zexchange_operator_single,       &
    dscdm_exchange_operator,         &
    zscdm_exchange_operator,         &
    zhamiltonian_dervexternal,       &
    zhamiltonian_apply_atom,         &
    hamiltonian_dump_vhxc,           &
    hamiltonian_load_vhxc,           &
    zoct_exchange_operator

  type hamiltonian_t
    !> The Hamiltonian must know what are the "dimensions" of the spaces,
    !! in order to be able to operate on the states.
    type(states_dim_t)       :: d
    type(hamiltonian_base_t) :: hm_base
    type(energy_t), pointer  :: energy
    type(base_hamiltonian_t), pointer :: subsys_hm    !< Subsystems Hamiltonian.
    type(bc_t)               :: bc      !< boundaries
    FLOAT, pointer :: vhartree(:) !< Hartree potential
    FLOAT, pointer :: vxc(:,:)    !< XC potential
    FLOAT, pointer :: vhxc(:,:)   !< XC potential + Hartree potential + Berry potential
    FLOAT, pointer :: axc(:,:,:)  !< XC vector potential divided by c
    FLOAT, pointer :: vtau(:,:)   !< Derivative of e_XC w.r.t. tau
    FLOAT, pointer :: vberry(:,:) !< Berry phase potential from external E_field
    !>cmplxscl: imaginary parts of the potentials
    FLOAT, pointer :: Imvhartree(:) !< Hartree potential
    FLOAT, pointer :: Imvxc(:,:)    !< XC potential
    FLOAT, pointer :: Imvhxc(:,:)   !< XC potential + Hartree potential + Berry potential
    FLOAT, pointer :: Imvtau(:,:)   !< Derivative of e_XC w.r.t. tau

    type(geometry_t), pointer :: geo
    FLOAT :: exx_coef !< how much of EXX to mix

    !> The self-induced vector potential and magnetic field
    logical :: self_induced_magnetic
    FLOAT, pointer :: a_ind(:, :)
    FLOAT, pointer :: b_ind(:, :)

    integer :: theory_level    !< copied from sys%ks
    integer :: xc_family       !< copied from sys%ks
    integer :: xc_flags        !< copied from sys%ks
    logical :: family_is_mgga_with_exc !< obtained from sys%ks

    type(epot_t) :: ep         !< handles the external potential
    type(pcm_t)  :: pcm        !< handles pcm variables
 
    !> absorbing boundaries
    logical :: adjoint

    !> Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span

    !> Mass of the particle (in most cases, mass = 1, electron mass)
    FLOAT :: mass
    !> anisotropic scaling factor for the mass: different along x,y,z etc...
    FLOAT :: mass_scaling(MAX_DIM)

    !> For the Hartree-Fock Hamiltonian, the Fock operator depends on the states.
    type(states_t), pointer :: hf_st
    !> use the SCDM method to compute the action of the Fock operator
    logical :: scdm_EXX

    !> There may be an "inhomogeneous", "source", or "forcing" term (useful for the OCT formalism)
    logical :: inh_term
    type(states_t) :: inh_st

    !> There may also be a exchange-like term, similar to the one necessary for time-dependent
    !! Hartree Fock, also useful only for the OCT equations
    type(oct_exchange_t) :: oct_exchange

    type(scissor_t) :: scissor

    FLOAT :: current_time
    FLOAT :: Imcurrent_time  !< needed when cmplxscl%time = .true.
    logical :: apply_packed  !< This is initialized by the StatesPack variable.
    
    !> If we use a complex-scaled Hamiltonian by complexifying the spatial coordinate with 
    !> the transformation r -> r*exp(i*theta)      
    type(cmplxscl_t) :: cmplxscl  !< complex scaling parameters

    !> For the Rashba spin-orbit coupling
    FLOAT :: rashba_coupling
    type(scdm_t)  :: scdm

    !> For the LDA+U 
    type(lda_u_t) :: lda_u
    integer       :: lda_u_level

    logical :: time_zero
  end type hamiltonian_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 2, &
    HARTREE               = 1, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4, &
    CLASSICAL             = 5, &
    RDMFT                 = 7

  type(profile_t), save :: prof_hamiltonian, prof_kinetic_start, prof_kinetic_finish
  type(profile_t), save :: prof_exx_scdm, prof_exx

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_init(hm, gr, geo, st, theory_level, xc_family, xc_flags, &
        family_is_mgga_with_exc, subsys_hm)
    type(hamiltonian_t),                        intent(out)   :: hm
    type(grid_t),                       target, intent(inout) :: gr
    type(geometry_t),                   target, intent(inout) :: geo
    type(states_t),                     target, intent(inout) :: st
    integer,                                    intent(in)    :: theory_level
    integer,                                    intent(in)    :: xc_family
    integer,                                    intent(in)    :: xc_flags
    logical,                                    intent(in)    :: family_is_mgga_with_exc
    type(base_hamiltonian_t), optional, target, intent(in)    :: subsys_hm

    integer :: iline, icol
    integer :: ncols
    type(block_t) :: blk
    type(profile_t), save :: prof

    logical :: external_potentials_present
    logical :: kick_present

    PUSH_SUB(hamiltonian_init)
    call profiling_in(prof, 'HAMILTONIAN_INIT')
    
    ! make a couple of local copies
    hm%theory_level = theory_level
    hm%xc_family    = xc_family
    hm%xc_flags     = xc_flags
    hm%family_is_mgga_with_exc = family_is_mgga_with_exc
    call states_dim_copy(hm%d, st%d)

    !%Variable ParticleMass
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian
    !%Description
    !% It is possible to make calculations for a particle with a mass
    !% different from one (atomic unit of mass, or mass of the electron).
    !% This is useful to describe non-electronic systems, or for
    !% esoteric purposes.
    !%End
    call parse_variable('ParticleMass', M_ONE, hm%mass)

    !%Variable RashbaSpinOrbitCoupling
    !%Type float
    !%Default 0.0
    !%Section Hamiltonian
    !%Description
    !% (Experimental.) For systems described in 2D (electrons confined to 2D in semiconductor structures), one
    !% may add the Bychkov-Rashba spin-orbit coupling term [Bychkov and Rashba, <i>J. Phys. C: Solid
    !% State Phys.</i> <b>17</b>, 6031 (1984)]. This variable determines the strength
    !% of this perturbation, and has dimensions of energy times length.
    !%End
    call parse_variable('RashbaSpinOrbitCoupling', M_ZERO, hm%rashba_coupling, units_inp%energy * units_inp%length)
    if(parse_is_defined('RashbaSpinOrbitCoupling')) then
      if(gr%sb%dim .ne. 2) then
        write(message(1),'(a)') 'Rashba spin-orbit coupling can only be used for two-dimensional systems.'
        call messages_fatal(1)
      end if
      call messages_experimental('RashbaSpinOrbitCoupling')
    end if

    call hamiltonian_base_init(hm%hm_base, hm%d%nspin, hm%mass, hm%rashba_coupling)

    ASSERT(associated(gr%der%lapl))
    hm%hm_base%kinetic => gr%der%lapl

    SAFE_ALLOCATE(hm%energy)
    call energy_nullify(hm%energy)

    call oct_exchange_nullify(hm%oct_exchange)
    
    !cmplxscl: copy cmplxscl initialized in states.F90
    call cmplxscl_copy(st%cmplxscl, hm%cmplxscl)

    nullify(hm%subsys_hm)
    if(present(subsys_hm))then
      ! Set Subsystems Hamiltonian pointer.
      ASSERT(.not.hm%cmplxscl%space)
      hm%subsys_hm => subsys_hm
    end if

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = hm%vxc(:, 1), vxc_22 = hm%vxc(:, 2), Re(vxc_12) = hm%vxc(:. 3);
    ! Im(vxc_12) = hm%vxc(:, 4)
    SAFE_ALLOCATE(hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin))
    hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin) = M_ZERO

    nullify(hm%vhartree, hm%vxc, hm%vtau, hm%axc)
    if(hm%theory_level /= INDEPENDENT_PARTICLES) then

      SAFE_ALLOCATE(hm%vhartree(1:gr%mesh%np_part))
      hm%vhartree=M_ZERO

      SAFE_ALLOCATE(hm%vxc(1:gr%mesh%np, 1:hm%d%nspin))
      hm%vxc=M_ZERO

      if(hm%family_is_mgga_with_exc) then
        SAFE_ALLOCATE(hm%vtau(1:gr%mesh%np, 1:hm%d%nspin))
        hm%vtau=M_ZERO
      end if

    end if

    nullify(hm%Imvhxc, hm%Imvhartree, hm%Imvxc, hm%Imvtau)

    if(hm%cmplxscl%space) then
      
      SAFE_ALLOCATE(hm%Imvhxc(1:gr%mesh%np, 1:hm%d%nspin))
      hm%Imvhxc(1:gr%mesh%np, 1:hm%d%nspin) = M_ZERO

      if(hm%theory_level /= INDEPENDENT_PARTICLES) then

        SAFE_ALLOCATE(hm%Imvhartree(1:gr%mesh%np))
        hm%Imvhartree=M_ZERO

        SAFE_ALLOCATE(hm%Imvxc(1:gr%mesh%np, 1:hm%d%nspin))
        hm%Imvxc=M_ZERO

        if(hm%family_is_mgga_with_exc) then
          SAFE_ALLOCATE(hm%Imvtau(1:gr%mesh%np, 1:hm%d%nspin))
          hm%Imvtau=M_ZERO
        end if
      end if
      
    end if

    hm%geo => geo
    !Initialize external potential
    call epot_init(hm%ep, gr, hm%geo, hm%d%ispin, hm%d%nik, hm%cmplxscl%space, subsys_hm,hm%xc_family)

    ! Calculate initial value of the gauge vector field
    call gauge_field_init(hm%ep%gfield, gr%sb)

    nullify(hm%vberry)
    if(associated(hm%ep%E_field) .and. simul_box_is_periodic(gr%sb) .and. .not. gauge_field_is_applied(hm%ep%gfield)) then
      ! only need vberry if there is a field in a periodic direction
      ! and we are not setting a gauge field
      if(any(abs(hm%ep%E_field(1:gr%sb%periodic_dim)) > M_EPSILON)) then
        SAFE_ALLOCATE(hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
      end if
    end if

    !Static magnetic field requires complex wavefunctions
    !Static magnetic field or rashba spin-orbit interaction requires complex wavefunctions
    if (associated(hm%ep%B_field) .or. &
      gauge_field_is_applied(hm%ep%gfield) .or. &
      parse_is_defined('RashbaSpinOrbitCoupling')) call states_set_complex(st)

    call parse_variable('CalculateSelfInducedMagneticField', .false., hm%self_induced_magnetic)
    !%Variable CalculateSelfInducedMagneticField
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% The existence of an electronic current implies the creation of a self-induced magnetic
    !% field, which may in turn back-react on the system. Of course, a fully consistent treatment
    !% of this kind of effect should be done in QED theory, but we will attempt a first
    !% approximation to the problem by considering the lowest-order relativistic terms
    !% plugged into the normal Hamiltonian equations (spin-other-orbit coupling terms, etc.).
    !% For the moment being, none of this is done, but a first step is taken by calculating
    !% the induced magnetic field of a system that has a current, by considering the magnetostatic
    !% approximation and Biot-Savart law:
    !%
    !% <math> \nabla^2 \vec{A} + 4\pi\alpha \vec{J} = 0</math>
    !%
    !% <math> \vec{B} = \vec{\nabla} \times \vec{A}</math>
    !%
    !% If <tt>CalculateSelfInducedMagneticField</tt> is set to yes, this <i>B</i> field is
    !% calculated at the end of a <tt>gs</tt> calculation (nothing is done -- yet -- in the <tt>td</tt>case)
    !% and printed out, if the <tt>Output</tt> variable contains the <tt>potential</tt> keyword (the prefix
    !% of the output files is <tt>Bind</tt>).
    !%End
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(hm%a_ind(1:gr%mesh%np_part, 1:gr%sb%dim))
      SAFE_ALLOCATE(hm%b_ind(1:gr%mesh%np_part, 1:gr%sb%dim))

      !(for dim = we could save some memory, but it is better to keep it simple)
    else
      nullify(hm%a_ind, hm%b_ind)
    end if

    ! Boundaries
    call bc_init(hm%bc, gr%mesh, gr%mesh%sb, hm%geo)

    !%Variable MassScaling
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% Scaling factor for anisotropic masses (different masses along each
    !% geometric direction).
    !%
    !% <tt>%MassScaling
    !% <br>&nbsp;&nbsp;1.0 | 1800.0 | 1800.0
    !% <br>%</tt>
    !%
    !% would fix the mass of the particles to be 1800 along the <i>y</i> and <i>z</i>
    !% directions. This can be useful, <i>e.g.</i>, to simulate 3 particles in 1D,
    !% in this case an electron and 2 protons.
    !%
    !%End
    hm%mass_scaling(1:gr%sb%dim) = M_ONE
    if(parse_block('MassScaling', blk) == 0) then
        ncols = parse_block_cols(blk, 0)
        if(ncols > gr%sb%dim) then
          call messages_input_error("MassScaling")
        end if
        iline = 1 ! just deal with 1 line - should be generalized
        do icol = 1, ncols
          call parse_block_float(blk, iline - 1, icol - 1, hm%mass_scaling(icol))
        end do
        call parse_block_end(blk)
    end if

    nullify(hm%hf_st)

    hm%inh_term = .false.
    call oct_exchange_remove(hm%oct_exchange)

    hm%adjoint = .false.

    !%Variable DFTULevel
    !%Type integer
    !%Default no
    !%Section Hamiltonian::XC
    !%Description
    !% (Experimental) This variable selects which DFT+U
    !% expression is added to the Hamiltonian.
    !%Option dft_u_none 0
    !% No +U term is not applied.
    !%Option dft_u_empirical 1
    !% An empiricial Hubbard U is added on the orbitals specified in the block species
    !% with hubbard_l and hubbard_u
    !%Option dft_u_acbn0 2
    !% Octopus determines the effective U term using the 
    !% ACBN0 functional as defined in PRX 5, 011006 (2015)
    !%End
    call parse_variable('DFTULevel', DFT_U_NONE, hm%lda_u_level)
    call messages_print_var_option(stdout,  'DFTULevel', hm%lda_u_level)
    call lda_u_nullify(hm%lda_u)
    if(hm%lda_u_level /= DFT_U_NONE) then
      call messages_experimental('DFT+U')
      call lda_u_init(hm%lda_u, hm%lda_u_level, gr, geo, st)
    end if
 

    nullify(hm%hm_base%phase)
    if (simul_box_is_periodic(gr%sb) .and. &
        .not. (kpoints_number(gr%sb%kpoints) == 1 .and. kpoints_point_is_gamma(gr%sb%kpoints, 1))) &
      call init_phase()
    ! no e^ik phase needed for Gamma-point-only periodic calculations

    !%Variable HamiltonianApplyPacked
    !%Type logical
    !%Default yes
    !%Section Execution::Optimization
    !%Description
    !% If set to yes (the default), Octopus will 'pack' the
    !% wave-functions when operating with them. This might involve some
    !% additional copying but makes operations more efficient.
    !% See also the related <tt>StatesPack</tt> variable.
    !%End
    call parse_variable('HamiltonianApplyPacked', .true., hm%apply_packed)

    external_potentials_present = associated(hm%ep%v_static) .or. &
				  associated(hm%ep%E_field)  .or. &
				  associated(hm%ep%lasers)

    kick_present = hm%ep%kick%delta_strength /= M_ZERO

    call pcm_init(hm%pcm, geo, gr, st%qtot, st%val_charge, external_potentials_present, kick_present )  !< initializes PCM  
    if(hm%pcm%run_pcm .and. hm%theory_level /= KOHN_SHAM_DFT) &
      call messages_not_implemented("PCM for TheoryLevel /= DFT")
    
    !%Variable SCDM_EXX
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% If set to yes, and <tt>TheoryLevel = hartree_fock</tt>,
    !% the Fock operator for exact exchange will be applied with the SCDM method.
    !%End
    call parse_variable('scdm_EXX', .false., hm%scdm_EXX)
    if(hm%scdm_EXX) then
      call messages_experimental("SCDM method for exact exchange")
      if(hm%theory_level /= HARTREE_FOCK) then
        call messages_not_implemented("SCDM for exact exchange in OEP (TheoryLevel = dft)")
      end if
       message(1) = "Info: Using SCDM for exact exchange"
       call messages_info(1)
    end if

    if(hm%theory_level == HARTREE_FOCK .and. st%parallel_in_states) then
#ifdef HAVE_MPI2
      call messages_experimental('Hartree-Fock parallel in states')
#else
      call messages_write('Hartree-Fock parallel in states required MPI 2')
      call messages_fatal()
#endif
    end if

    !%Variable TimeZero
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% (Experimental) If set to yes, the ground state and other time
    !% dependent calculation will assume that they are done at time
    !% zero, so that all time depedent field at that time will be
    !% included.
    !%End
    call parse_variable('TimeZero', .false., hm%time_zero)
    if(hm%time_zero) call messages_experimental('TimeZero')

    call scissor_nullify(hm%scissor)

    call profiling_out(prof)
    POP_SUB(hamiltonian_init)

  contains

    ! ---------------------------------------------------------
    subroutine init_phase
      integer :: ip, ik
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_init.init_phase)

      SAFE_ALLOCATE(hm%hm_base%phase(1:gr%mesh%np_part, hm%d%kpt%start:hm%d%kpt%end))

      kpoint(1:gr%sb%dim) = M_ZERO
      do ik = hm%d%kpt%start, hm%d%kpt%end
        kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(hm%d, ik))
        forall (ip = 1:gr%mesh%np_part)
          hm%hm_base%phase(ip, ik) = exp(-M_zI * sum(gr%mesh%x(ip, 1:gr%sb%dim) * kpoint(1:gr%sb%dim)))
        end forall
      end do

      if(accel_is_enabled()) then
        call accel_create_buffer(hm%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, gr%mesh%np_part*hm%d%kpt%nlocal)
        call accel_write_buffer(hm%hm_base%buff_phase, gr%mesh%np_part*hm%d%kpt%nlocal, hm%hm_base%phase)
        hm%hm_base%buff_phase_qn_start = hm%d%kpt%start
      end if

      ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
      if(hm%lda_u_level /= DFT_U_NONE) then
        call lda_u_build_phase_correction(hm%lda_u, gr%mesh%sb, hm%d )
      end if

      POP_SUB(hamiltonian_init.init_phase)
    end subroutine init_phase

  end subroutine hamiltonian_init

  
  ! ---------------------------------------------------------
  subroutine hamiltonian_end(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_end)

    call hamiltonian_base_end(hm%hm_base)

    nullify(hm%subsys_hm)
    
    if(associated(hm%hm_base%phase) .and. accel_is_enabled()) then
      call accel_release_buffer(hm%hm_base%buff_phase)
    end if

    SAFE_DEALLOCATE_P(hm%hm_base%phase)
    SAFE_DEALLOCATE_P(hm%vhartree)
    SAFE_DEALLOCATE_P(hm%vhxc)
    SAFE_DEALLOCATE_P(hm%vxc)
    SAFE_DEALLOCATE_P(hm%axc)
    SAFE_DEALLOCATE_P(hm%vberry)
    SAFE_DEALLOCATE_P(hm%a_ind)
    SAFE_DEALLOCATE_P(hm%b_ind)
    !cmplxscl
    SAFE_DEALLOCATE_P(hm%Imvhartree)
    SAFE_DEALLOCATE_P(hm%Imvhxc)
    SAFE_DEALLOCATE_P(hm%Imvxc)
    SAFE_DEALLOCATE_P(hm%Imvtau)
    
    if(hm%family_is_mgga_with_exc) then
      SAFE_DEALLOCATE_P(hm%vtau)
    end if

    call epot_end(hm%ep)
    nullify(hm%geo)

    call bc_end(hm%bc)

    call states_dim_end(hm%d) 

    if(hm%scissor%apply) call scissor_end(hm%scissor)

    ! this is a bit ugly, hf_st is initialized in v_ks_calc but deallocated here.
    if(associated(hm%hf_st))  then
      if(hm%hf_st%parallel_in_states) call states_parallel_remote_access_stop(hm%hf_st)
      call states_end(hm%hf_st)
      SAFE_DEALLOCATE_P(hm%hf_st)
    end if

    SAFE_DEALLOCATE_P(hm%energy)
     
    if (hm%pcm%run_pcm) call pcm_end(hm%pcm)

    POP_SUB(hamiltonian_end)
  end subroutine hamiltonian_end


  ! ---------------------------------------------------------
  ! True if the Hamiltonian is Hermitian, false otherwise
  logical function hamiltonian_hermitian(hm)
    type(hamiltonian_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_hermitian)
    hamiltonian_hermitian = .not.((hm%bc%abtype == IMAGINARY_ABSORBING) .or. &
                                  oct_exchange_enabled(hm%oct_exchange)     .or. &
                                  hm%cmplxscl%space)

    POP_SUB(hamiltonian_hermitian)
  end function hamiltonian_hermitian


  ! ---------------------------------------------------------
  subroutine hamiltonian_span(hm, delta, emin)
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT,               intent(in)    :: delta, emin

    PUSH_SUB(hamiltonian_span)

    hm%spectral_middle_point = ((M_PI**2 / (2 * delta**2)) + emin) / M_TWO
    hm%spectral_half_span    = ((M_PI**2 / (2 * delta**2)) - emin) / M_TWO

    POP_SUB(hamiltonian_span)
  end subroutine hamiltonian_span


  ! ---------------------------------------------------------
  pure logical function hamiltonian_inh_term(hm) result(inh)
    type(hamiltonian_t), intent(in) :: hm

    inh = hm%inh_term
  end function hamiltonian_inh_term


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_inh(hm, st)
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(in)    :: st

    PUSH_SUB(hamiltonian_set_inh)

    if(hm%inh_term) call states_end(hm%inh_st)
    call states_copy(hm%inh_st, st)
    hm%inh_term = .true.

    POP_SUB(hamiltonian_set_inh)
  end subroutine hamiltonian_set_inh


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_inh(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_remove_inh)

    if(hm%inh_term) then
      call states_end(hm%inh_st)
      hm%inh_term = .false.
    end if

    POP_SUB(hamiltonian_remove_inh)
  end subroutine hamiltonian_remove_inh

  ! ---------------------------------------------------------
  subroutine hamiltonian_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_adjoint)

    if(.not.hm%adjoint) then
      hm%adjoint = .true.
      if(hm%bc%abtype == IMAGINARY_ABSORBING) then
        hm%bc%mf = -hm%bc%mf
      end if
    end if

    POP_SUB(hamiltonian_adjoint)
  end subroutine hamiltonian_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_not_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_not_adjoint)

    if(hm%adjoint) then
      hm%adjoint = .false.
      if(hm%bc%abtype == IMAGINARY_ABSORBING) then
        hm%bc%mf = -hm%bc%mf
      end if
    end if

    POP_SUB(hamiltonian_not_adjoint)
  end subroutine hamiltonian_not_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_update(this, mesh, time, Imtime)
    type(hamiltonian_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT, optional,     intent(in)    :: time
    FLOAT, optional,     intent(in)    :: Imtime

    integer :: ispin, ip, idir, iatom, ilaser
    type(profile_t), save :: prof, prof_phases
    FLOAT :: aa(1:MAX_DIM), time_
    FLOAT, allocatable :: vp(:,:)

    PUSH_SUB(hamiltonian_update)
    call profiling_in(prof, "HAMILTONIAN_UPDATE")

    this%current_time = M_ZERO
    this%Imcurrent_time = M_ZERO !cmplxscl
    if(present(time)) this%current_time = time
    if(present(Imtime)) this%Imcurrent_time = Imtime !cmplxscl

    time_ = optional_default(time, CNST(0.0))

    ! set everything to zero
    call hamiltonian_base_clear(this%hm_base)

    ! the xc, hartree and external potentials
    call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_POTENTIAL, &
      complex_potential = this%cmplxscl%space .or. this%bc%abtype == IMAGINARY_ABSORBING)


    do ispin = 1, this%d%nspin
      if(ispin <= 2) then
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin) + this%ep%vpsl(ip)
        !> Adds PCM contributions
        if (this%pcm%run_pcm) then
          if (this%pcm%solute) then
            forall (ip = 1:mesh%np)  
              this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + &
                this%pcm%v_e_rs(ip) + this%pcm%v_n_rs(ip)
            end forall
          end if
          if (this%pcm%localf) then
            forall (ip = 1:mesh%np)  
              this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + &
                this%pcm%v_ext_rs(ip)
            end forall
          end if 
        end if

        if(this%cmplxscl%space) then
          forall (ip = 1:mesh%np)
            this%hm_base%Impotential(ip, ispin) = &
              this%hm_base%Impotential(ip, ispin) + this%Imvhxc(ip, ispin) +  this%ep%Imvpsl(ip)
          end forall
        end if
        
        if(this%bc%abtype == IMAGINARY_ABSORBING) then
          forall (ip = 1:mesh%np)
            this%hm_base%Impotential(ip, ispin) = this%hm_base%Impotential(ip, ispin) + this%bc%mf(ip)
          end forall
        end if

      else !Spinors 
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin)
        if(this%cmplxscl%space) then
          forall (ip = 1:mesh%np) this%hm_base%Impotential(ip, ispin) = this%Imvhxc(ip, ispin)
        end if
          
      end if


    end do

    ! the lasers
    if (present(time) .or. this%time_zero) then

      do ilaser = 1, this%ep%no_lasers
        select case(laser_kind(this%ep%lasers(ilaser)))
        case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
          do ispin = 1, this%d%spin_channels
            call laser_potential(this%ep%lasers(ilaser), mesh,  this%hm_base%potential(:, ispin), time_)
          end do
        case(E_FIELD_MAGNETIC)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL + FIELD_UNIFORM_MAGNETIC_FIELD, &
            this%cmplxscl%space)
          ! get the vector potential
          SAFE_ALLOCATE(vp(1:mesh%np, 1:mesh%sb%dim))
          vp(1:mesh%np, 1:mesh%sb%dim) = M_ZERO
          call laser_vector_potential(this%ep%lasers(ilaser), mesh, vp, time_)
          forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
            this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) - vp(ip, idir)/P_C
          end forall
          ! and the magnetic field
          call laser_field(this%ep%lasers(ilaser), this%hm_base%uniform_magnetic_field(1:mesh%sb%dim), time_)
          SAFE_DEALLOCATE_A(vp)
        case(E_FIELD_VECTOR_POTENTIAL)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, this%cmplxscl%space)
          ! get the uniform vector potential associated with a magnetic field
          aa = M_ZERO
          call laser_field(this%ep%lasers(ilaser), aa(1:mesh%sb%dim), time_)
          this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim) &
            - aa(1:mesh%sb%dim)/P_C
        end select
      end do

      ! the gauge field
      if(gauge_field_is_applied(this%ep%gfield)) then
        call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, this%cmplxscl%space)
        call gauge_field_get_vec_pot(this%ep%gfield, aa)
        this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim)  &
          - aa(1:mesh%sb%dim)/P_c
      end if

      ! the electric field for a periodic system through the gauge field
      if(associated(this%ep%e_field) .and. gauge_field_is_applied(this%ep%gfield)) then
        this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) = &
          this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) - time_*this%ep%e_field(1:mesh%sb%periodic_dim)
      end if
      
    end if

    ! the vector potential of a static magnetic field
    if(associated(this%ep%a_static)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL, this%cmplxscl%space)
      forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
        this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) + this%ep%a_static(ip, idir)
      end forall
    end if

    ! and the static magnetic field
    if(associated(this%ep%b_field)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_MAGNETIC_FIELD, this%cmplxscl%space)
      forall (idir = 1:3)
        this%hm_base%uniform_magnetic_field(idir) = this%hm_base%uniform_magnetic_field(idir) + this%ep%b_field(idir)
      end forall
    end if

    call hamiltonian_base_update(this%hm_base, mesh)

    call build_phase()

    call profiling_out(prof)
    POP_SUB(hamiltonian_update)

  contains

    subroutine build_phase()
      integer :: ik, imat, nmat, max_npoints, offset
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_update.build_phase)

      if(simul_box_is_periodic(mesh%sb) .or. allocated(this%hm_base%uniform_vector_potential)) then

        call profiling_in(prof_phases, 'UPDATE_PHASES')
        ! now regenerate the phases for the pseudopotentials
        do iatom = 1, this%ep%natoms
          call projector_init_phases(this%ep%proj(iatom), mesh%sb, this%d, &
            vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end do

        call profiling_out(prof_phases)
      end if

      if(allocated(this%hm_base%uniform_vector_potential)) then
        if(.not. associated(this%hm_base%phase)) then
          SAFE_ALLOCATE(this%hm_base%phase(1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part*this%d%kpt%nlocal)
          end if
        end if

        kpoint(1:mesh%sb%dim) = M_ZERO
        do ik = this%d%kpt%start, this%d%kpt%end
          kpoint(1:mesh%sb%dim) = kpoints_get_point(mesh%sb%kpoints, states_dim_get_kpoint_index(this%d, ik))

          forall (ip = 1:mesh%np_part)
            this%hm_base%phase(ip, ik) = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*(kpoint(1:mesh%sb%dim) &
              + this%hm_base%uniform_vector_potential(1:mesh%sb%dim))))
          end forall
        end do
        if(accel_is_enabled()) then
          call accel_write_buffer(this%hm_base%buff_phase, mesh%np_part*this%d%kpt%nlocal, this%hm_base%phase)
        end if

        ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
        if(this%lda_u_level /= DFT_U_NONE) then
          call lda_u_build_phase_correction(this%lda_u, mesh%sb, this%d, &
               vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end if


      end if

      max_npoints = this%hm_base%max_npoints
      nmat = this%hm_base%nprojector_matrices


      if(associated(this%hm_base%phase) .and. allocated(this%hm_base%projector_matrices)) then

        if(.not. allocated(this%hm_base%projector_phases)) then
          SAFE_ALLOCATE(this%hm_base%projector_phases(1:max_npoints, nmat, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_projector_phases, ACCEL_MEM_READ_ONLY, &
              TYPE_CMPLX, this%hm_base%total_points*this%d%kpt%nlocal)
          end if
        end if

        offset = 0
        do ik = this%d%kpt%start, this%d%kpt%end
          do imat = 1, this%hm_base%nprojector_matrices
            iatom = this%hm_base%projector_to_atom(imat)
            do ip = 1, this%hm_base%projector_matrices(imat)%npoints
              this%hm_base%projector_phases(ip, imat, ik) = this%ep%proj(iatom)%phase(ip, ik)
            end do

            if(accel_is_enabled() .and. this%hm_base%projector_matrices(imat)%npoints > 0) then
              call accel_write_buffer(this%hm_base%buff_projector_phases, &
                this%hm_base%projector_matrices(imat)%npoints, this%hm_base%projector_phases(1:, imat, ik), offset = offset)
            end if
            offset = offset + this%hm_base%projector_matrices(imat)%npoints
          end do
        end do

      end if

      POP_SUB(hamiltonian_update.build_phase)
    end subroutine build_phase

  end subroutine hamiltonian_update


  ! ---------------------------------------------------------
  subroutine hamiltonian_epot_generate(this, gr, geo, st, time)
    type(hamiltonian_t),      intent(inout) :: this
    type(grid_t),             intent(in)    :: gr
    type(geometry_t), target, intent(inout) :: geo
    type(states_t),           intent(inout) :: st
    FLOAT,          optional, intent(in)    :: time

    PUSH_SUB(hamiltonian_epot_generate)

    this%geo => geo
    call epot_generate(this%ep, gr, this%geo, st, this%cmplxscl%space)
    call hamiltonian_base_build_proj(this%hm_base, gr%mesh, this%ep)
    call hamiltonian_update(this, gr%mesh, time)
   
    if (this%pcm%run_pcm) then
     !> Generates the real-space PCM potential due to nuclei which do not change
     !! during the SCF calculation.
     if (this%pcm%solute) &
       call pcm_calc_pot_rs(this%pcm, gr%mesh, geo = geo)

      !> Local field effects due to static electrostatic potentials (if they were).
      !! The laser and the kick are included in subroutine v_ks_hartree (module v_ks).
      !  Interpolation is needed, hence gr%mesh%np_part -> 1:gr%mesh%np
      if( this%pcm%localf .and. associated(this%ep%v_static)) &
        call pcm_calc_pot_rs(this%pcm, gr%mesh, v_ext = this%ep%v_ext(1:gr%mesh%np_part))

    end if

    call lda_u_update_basis(this%lda_u, gr, geo, st, associated(this%hm_base%phase))

    POP_SUB(hamiltonian_epot_generate)
  end subroutine hamiltonian_epot_generate

  ! -----------------------------------------------------------------

  FLOAT function hamiltonian_get_time(this) result(time)
    type(hamiltonian_t),   intent(inout) :: this

    time = this%current_time
  end function hamiltonian_get_time

  ! -----------------------------------------------------------------

  logical pure function hamiltonian_apply_packed(this, mesh) result(apply)
    type(hamiltonian_t),   intent(in) :: this
    type(mesh_t),          intent(in) :: mesh

    apply = this%apply_packed
    if(mesh%use_curvilinear) apply = .false.
    if(hamiltonian_base_has_magnetic(this%hm_base)) apply = .false.
    if(this%rashba_coupling**2 > M_ZERO) apply = .false.
    if(this%ep%non_local .and. .not. this%hm_base%apply_projector_matrices) apply = .false.
    if(this%family_is_mgga_with_exc)  apply = .false. 
    if(this%scissor%apply) apply = .false.
    if(this%bc%abtype == IMAGINARY_ABSORBING .and. accel_is_enabled()) apply = .false.
    if(this%cmplxscl%space .and. accel_is_enabled()) apply = .false.
    if(associated(this%hm_base%phase) .and. accel_is_enabled()) apply = .false.
    
  end function hamiltonian_apply_packed

  ! -----------------------------------------------------------------
  !> This routine computes the action of the derivative of the external potential
  !! with respect to the nuclear positions. It is preliminary, and should be
  !! recoded in a more efficient way.
  subroutine zhamiltonian_dervexternal(hm, geo, gr, ia, dim, psi, dvpsi)
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t),    intent(in)  :: geo
    type(grid_t),        intent(in)  :: gr
    integer,             intent(in)  :: ia
    integer,             intent(in)  :: dim
    CMPLX,               intent(inout)  :: psi(:, :)
    CMPLX,               intent(out) :: dvpsi(:, :, :)

    CMPLX, allocatable :: dpsi(:, :, :), dvlocalpsi(:, :, :), vlocalpsi(:, :)
    integer :: idim, j

    PUSH_SUB(zhamiltonian_dervexternal)

    SAFE_ALLOCATE(vlocalpsi(1:gr%mesh%np_part, 1:dim))
    SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:gr%sb%dim, 1:dim))
    SAFE_ALLOCATE(dvlocalpsi(1:gr%mesh%np_part, 1:gr%sb%dim, 1:dim))

    vlocalpsi = M_ZERO
    dpsi = M_z0
    dvlocalpsi = M_z0

    do idim = 1, dim
      call zderivatives_grad(gr%der, psi(:, idim), dpsi(:, :, idim))
    end do
    call zhamiltonian_apply_atom (hm, geo, gr, ia, psi, vlocalpsi)

    do idim = 1, dim
      call zderivatives_grad(gr%der, vlocalpsi(:, idim), dvlocalpsi(:, :, idim))
    end do
    
    ! Various ways to do the same thing:
    ! (1)
    !    _SAFE_ALLOCATE(dvlocal(1:gr%mesh%np, 1:gr%sb%dim))
    !    call dderivatives_grad(gr%der, vlocal, dvlocal)
    !    do idim = 1, dim
    !      do ip = 1, gr%mesh%np
    !        call mesh_r(gr%mesh, ip, rr, coords = xx, origin = qa)
    !        dvpsi(ip, idim, 1) = (xx(1) / sqrt( (xx(1)**2+M_ONE)**3 ) ) * psi(ip, idim)
    !      end do
    !    end do
    !    _SAFE_DEALLOCATE_A(dvlocal)
    !
    ! (2)
    !    do idim = 1, dim
    !      do ip = 1, gr%mesh%np
    !        dvpsi(ip, idim) = dvlocal(ip, 1) * psi(ip, idim)
    !      end do
    !    end do
    !
    ! (3)

    do j = 1, gr%sb%dim
      call zhamiltonian_apply_atom (hm, geo, gr, ia, dpsi(:, j, :), vlocalpsi)
      dvpsi(:, :, j) = -vlocalpsi(:, :) + dvlocalpsi(:, j, :)
    end do

    SAFE_DEALLOCATE_A(vlocalpsi)
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dvlocalpsi)
    POP_SUB(zhamiltonian_dervexternal)
  end subroutine zhamiltonian_dervexternal


  ! -----------------------------------------------------------------
  subroutine zhamiltonian_apply_atom (hm, geo, gr, ia, psi, vpsi)
    type(hamiltonian_t), intent(inout) :: hm
    type(geometry_t),    intent(in)    :: geo
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: ia
    CMPLX,               intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
    CMPLX,               intent(out)   :: vpsi(:,:) !< (gr%mesh%np, hm%d%dim)

    integer :: idim
    FLOAT, allocatable :: vlocal(:)
    PUSH_SUB(zhamiltonian_apply_atom)

    SAFE_ALLOCATE(vlocal(1:gr%mesh%np_part))
    vlocal = M_ZERO
    call epot_local_potential(hm%ep, gr%der, gr%dgrid, geo, ia, vlocal)

    do idim = 1, hm%d%dim
      vpsi(1:gr%mesh%np, idim)  = vlocal(1:gr%mesh%np) * psi(1:gr%mesh%np, idim)
    end do


    SAFE_DEALLOCATE_A(vlocal)
    POP_SUB(zhamiltonian_apply_atom)
  end subroutine zhamiltonian_apply_atom


  ! -----------------------------------------------------------------
  subroutine hamiltonian_dump_vhxc(restart, hm, mesh, ierr)
    type(restart_t),     intent(in)  :: restart
    type(hamiltonian_t), intent(in)  :: hm
    type(mesh_t),        intent(in)  :: mesh
    integer,             intent(out) :: ierr

    integer :: iunit, err, err2(2), isp
    character(len=12) :: filename
    character(len=100) :: lines(2)

    PUSH_SUB(hamiltonian_dump_vhxc)

    ierr = 0

    if (restart_skip(restart) .or. hm%theory_level == INDEPENDENT_PARTICLES) then
      POP_SUB(hamiltonian_dump_vhxc)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing Vhxc restart."
      call messages_info(1)
    end if

    !write the different components of the Hartree+XC potential
    iunit = restart_open(restart, 'vhxc')
    lines(1) = '#     #spin    #nspin    filename'
    lines(2) = '%vhxc'
    call restart_write(restart, iunit, lines, 2, err)
    if (err /= 0) ierr = ierr + 1

    err2 = 0
    do isp = 1, hm%d%nspin
      if (hm%d%nspin == 1) then
        write(filename, fmt='(a)') 'vhxc'
      else
        write(filename, fmt='(a,i1)') 'vhxc-sp', isp
      end if
      write(lines(1), '(i8,a,i8,a)') isp, ' | ', hm%d%nspin, ' | "'//trim(adjustl(filename))//'"'
      call restart_write(restart, iunit, lines, 1, err)
      if (err /= 0) err2(1) = err2(1) + 1

      if (hm%cmplxscl%space) then
        call zrestart_write_mesh_function(restart, filename, mesh, hm%vhxc(:,isp) + M_zI*hm%imvhxc(:,isp), err)
      else
        call drestart_write_mesh_function(restart, filename, mesh, hm%vhxc(:,isp), err)
      end if
      if (err /= 0) err2(2) = err2(2) + 1

    end do
    if (err2(1) /= 0) ierr = ierr + 2
    if (err2(2) /= 0) ierr = ierr + 4

    lines(1) = '%'
    call restart_write(restart, iunit, lines, 1, err)
    if (err /= 0) ierr = ierr + 4

    ! MGGAs and hybrid MGGAs have an extra term that also needs to be dumped
    if (hm%family_is_mgga_with_exc) then
      lines(1) = '#     #spin    #nspin    filename'
      lines(2) = '%vtau'
      call restart_write(restart, iunit, lines, 2, err)
      if (err /= 0) ierr = ierr + 8

      err2 = 0
      do isp = 1, hm%d%nspin
        if (hm%d%nspin == 1) then
          write(filename, fmt='(a)') 'vtau'
        else
          write(filename, fmt='(a,i1)') 'vtau-sp', isp
        end if
        write(lines(1), '(i8,a,i8,a)') isp, ' | ', hm%d%nspin, ' | "'//trim(adjustl(filename))//'"'
        call restart_write(restart, iunit, lines, 1, err)
        if (err /= 0) err2(1) = err2(1) + 16

        if (hm%cmplxscl%space) then
          call zrestart_write_mesh_function(restart, filename, mesh, hm%vtau(:,isp) + M_zI*hm%imvtau(:,isp), err)
        else
          call drestart_write_mesh_function(restart, filename, mesh, hm%vtau(:,isp), err)
        end if
        if (err /= 0) err2(1) = err2(1) + 1

      end do
      if (err2(1) /= 0) ierr = ierr + 32
      if (err2(2) /= 0) ierr = ierr + 64

      lines(1) = '%'
      call restart_write(restart, iunit, lines, 1, err)
      if (err /= 0) ierr = ierr + 128
    end if

    call restart_close(restart, iunit)

    if (debug%info) then
      message(1) = "Debug: Writing Vhxc restart done."
      call messages_info(1)
    end if

    POP_SUB(hamiltonian_dump_vhxc)
  end subroutine hamiltonian_dump_vhxc


  ! ---------------------------------------------------------
  subroutine hamiltonian_load_vhxc(restart, hm, mesh, ierr)
    type(restart_t),     intent(in)    :: restart
    type(hamiltonian_t), intent(inout) :: hm
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(out)   :: ierr

    integer :: err, err2, isp
    character(len=12) :: filename
    CMPLX, allocatable :: zv(:)

    PUSH_SUB(hamiltonian_load_vhxc)

    ierr = 0

    if (restart_skip(restart) .or. hm%theory_level == INDEPENDENT_PARTICLES) then
      ierr = -1
      POP_SUB(hamiltonian_load_vhxc)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading Vhxc restart."
      call messages_info(1)
    end if

    if (hm%cmplxscl%space) then
      SAFE_ALLOCATE(zv(1:mesh%np))
    end if

    err2 = 0
    do isp = 1, hm%d%nspin
      if (hm%d%nspin==1) then
        write(filename, fmt='(a)') 'vhxc'
      else
        write(filename, fmt='(a,i1)') 'vhxc-sp', isp
      end if

      if (hm%cmplxscl%space) then
        call zrestart_read_mesh_function(restart, filename, mesh, zv, err)
        hm%vhxc(:,isp) =  real(zv, REAL_PRECISION)
        hm%imvhxc(:,isp) = aimag(zv)
      else
        call drestart_read_mesh_function(restart, filename, mesh, hm%vhxc(:,isp), err)
      end if
      if (err /= 0) err2 = err2 + 1

    end do
    if (err2 /= 0) ierr = ierr + 1

    ! MGGAs and hybrid MGGAs have an extra term that also needs to be read
    err2 = 0
    if (hm%family_is_mgga_with_exc) then
      do isp = 1, hm%d%nspin
        if (hm%d%nspin == 1) then
          write(filename, fmt='(a)') 'vtau'
        else
          write(filename, fmt='(a,i1)') 'vtau-sp', isp
        end if

        if (hm%cmplxscl%space) then
          call zrestart_read_mesh_function(restart, filename, mesh, zv, err)
          hm%vtau(:,isp) =  real(zv, REAL_PRECISION)
          hm%imvtau(:,isp) = aimag(zv)
        else
          call drestart_read_mesh_function(restart, filename, mesh, hm%vtau(:,isp), err)
        end if
        if (err /= 0) err2 = err2 + 1

      end do

      if (err2 /= 0) ierr = ierr + 2
    end if

    if (hm%cmplxscl%space) then
      SAFE_DEALLOCATE_A(zv)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading Vhxc restart done."
      call messages_info(1)
    end if

    POP_SUB(hamiltonian_load_vhxc)
  end subroutine hamiltonian_load_vhxc

  ! ---------------------------------------------------------
  ! This is an extension of "hamiltonian_update2" to be used by the
  ! CFM4 propagator. It updates the Hamiltonian by considering a
  ! weighted sum of the external potentials at times time(1) and time(2),
  ! weighted by alpha(1) and alpha(2).
  subroutine hamiltonian_update2(this, mesh, time, mu,  Imtime)
    type(hamiltonian_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: time(1:2)
    FLOAT,               intent(in)    :: mu(1:2)
    FLOAT, optional,     intent(in)    :: Imtime(1:2)

    integer :: ispin, ip, idir, iatom, ilaser, itime
    type(profile_t), save :: prof, prof_phases
    FLOAT :: aa(1:MAX_DIM), time_
    FLOAT, allocatable :: vp(:,:)

    FLOAT, allocatable :: velectric(:)

    PUSH_SUB(hamiltonian_update2)
    call profiling_in(prof, "HAMILTONIAN_UPDATE")

    this%current_time = M_ZERO
    this%Imcurrent_time = M_ZERO !cmplxscl
    this%current_time = time(1)
    if(present(Imtime)) this%Imcurrent_time = Imtime(1) !cmplxscl

    ! set everything to zero
    call hamiltonian_base_clear(this%hm_base)

    ! the xc, hartree and external potentials
    call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_POTENTIAL, &
      complex_potential = this%cmplxscl%space .or. this%bc%abtype == IMAGINARY_ABSORBING)


    do ispin = 1, this%d%nspin
      if(ispin <= 2) then
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin) + this%ep%vpsl(ip)
        !> Adds PCM contributions
        if (this%pcm%run_pcm) then
          forall (ip = 1:mesh%np)
            this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + &
              this%pcm%v_e_rs(ip) + this%pcm%v_n_rs(ip)
          end forall
        end if

        if(this%cmplxscl%space) then
          forall (ip = 1:mesh%np)
            this%hm_base%Impotential(ip, ispin) = &
              this%hm_base%Impotential(ip, ispin) + this%Imvhxc(ip, ispin) +  this%ep%Imvpsl(ip)
          end forall
        end if

        if(this%bc%abtype == IMAGINARY_ABSORBING) then
          forall (ip = 1:mesh%np)
            this%hm_base%Impotential(ip, ispin) = this%hm_base%Impotential(ip, ispin) + this%bc%mf(ip)
          end forall
        end if

      else !Spinors
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin)
        if(this%cmplxscl%space) then
          forall (ip = 1:mesh%np) this%hm_base%Impotential(ip, ispin) = this%Imvhxc(ip, ispin)
        end if

      end if


    end do


    do itime = 1, 2
      time_ = time(itime)

      do ilaser = 1, this%ep%no_lasers
        select case(laser_kind(this%ep%lasers(ilaser)))
        case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
          SAFE_ALLOCATE(velectric(1:mesh%np))
          do ispin = 1, this%d%spin_channels
            velectric = M_ZERO
            call laser_potential(this%ep%lasers(ilaser), mesh,  velectric, time_)
            this%hm_base%potential(:, ispin) = this%hm_base%potential(:, ispin) + mu(itime) * velectric(:)
          end do
          SAFE_DEALLOCATE_A(velectric)
        case(E_FIELD_MAGNETIC)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL + FIELD_UNIFORM_MAGNETIC_FIELD, &
            this%cmplxscl%space)
          ! get the vector potential
          SAFE_ALLOCATE(vp(1:mesh%np, 1:mesh%sb%dim))
          vp(1:mesh%np, 1:mesh%sb%dim) = M_ZERO
          call laser_vector_potential(this%ep%lasers(ilaser), mesh, vp, time_)
          forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
            this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) &
              - mu(itime) * vp(ip, idir)/P_C
          end forall
          ! and the magnetic field
          call laser_field(this%ep%lasers(ilaser), this%hm_base%uniform_magnetic_field(1:mesh%sb%dim), time_)
          SAFE_DEALLOCATE_A(vp)
        case(E_FIELD_VECTOR_POTENTIAL)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, this%cmplxscl%space)
          ! get the uniform vector potential associated with a magnetic field
          aa = M_ZERO
          call laser_field(this%ep%lasers(ilaser), aa(1:mesh%sb%dim), time_)
          this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim) &
            - mu(itime) * aa(1:mesh%sb%dim)/P_C
        end select
      end do

      ! the gauge field
      if(gauge_field_is_applied(this%ep%gfield)) then
        call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, this%cmplxscl%space)
        call gauge_field_get_vec_pot(this%ep%gfield, aa)
        this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim)  &
          - aa(1:mesh%sb%dim)/P_c
      end if

      ! the electric field for a periodic system through the gauge field
      if(associated(this%ep%e_field) .and. gauge_field_is_applied(this%ep%gfield)) then
        this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) = &
          this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) - time_*this%ep%e_field(1:mesh%sb%periodic_dim)
      end if

    end do

    ! the vector potential of a static magnetic field
    if(associated(this%ep%a_static)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL, this%cmplxscl%space)
      forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
        this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) + this%ep%a_static(ip, idir)
      end forall
    end if

    ! and the static magnetic field
    if(associated(this%ep%b_field)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_MAGNETIC_FIELD, this%cmplxscl%space)
      forall (idir = 1:3)
        this%hm_base%uniform_magnetic_field(idir) = this%hm_base%uniform_magnetic_field(idir) + this%ep%b_field(idir)
      end forall
    end if

    call hamiltonian_base_update(this%hm_base, mesh)

    call build_phase()

    call profiling_out(prof)
    POP_SUB(hamiltonian_update2)

  contains

    subroutine build_phase()
      integer :: ik, imat, nmat, max_npoints, offset
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_update2.build_phase)

      if(simul_box_is_periodic(mesh%sb) .or. allocated(this%hm_base%uniform_vector_potential)) then

        call profiling_in(prof_phases, 'UPDATE_PHASES')
        ! now regenerate the phases for the pseudopotentials
        do iatom = 1, this%ep%natoms
          call projector_init_phases(this%ep%proj(iatom), mesh%sb, this%d, &
            vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end do

        call profiling_out(prof_phases)
      end if

      if(allocated(this%hm_base%uniform_vector_potential)) then
        if(.not. associated(this%hm_base%phase)) then
          SAFE_ALLOCATE(this%hm_base%phase(1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part*this%d%kpt%nlocal)
          end if
        end if

        kpoint(1:mesh%sb%dim) = M_ZERO
        do ik = this%d%kpt%start, this%d%kpt%end
          kpoint(1:mesh%sb%dim) = kpoints_get_point(mesh%sb%kpoints, states_dim_get_kpoint_index(this%d, ik))

          forall (ip = 1:mesh%np_part)
            this%hm_base%phase(ip, ik) = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*(kpoint(1:mesh%sb%dim) &
              + this%hm_base%uniform_vector_potential(1:mesh%sb%dim))))
          end forall
        end do
        if(accel_is_enabled()) then
          call accel_write_buffer(this%hm_base%buff_phase, mesh%np_part*this%d%kpt%nlocal, this%hm_base%phase)
        end if
      end if

      max_npoints = this%hm_base%max_npoints
      nmat = this%hm_base%nprojector_matrices


      if(associated(this%hm_base%phase) .and. allocated(this%hm_base%projector_matrices)) then

        if(.not. allocated(this%hm_base%projector_phases)) then
          SAFE_ALLOCATE(this%hm_base%projector_phases(1:max_npoints, nmat, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_projector_phases, ACCEL_MEM_READ_ONLY, &
              TYPE_CMPLX, this%hm_base%total_points*this%d%kpt%nlocal)
          end if
        end if

        offset = 0
        do ik = this%d%kpt%start, this%d%kpt%end
          do imat = 1, this%hm_base%nprojector_matrices
            iatom = this%hm_base%projector_to_atom(imat)
            do ip = 1, this%hm_base%projector_matrices(imat)%npoints
              this%hm_base%projector_phases(ip, imat, ik) = this%ep%proj(iatom)%phase(ip, ik)
            end do

            if(accel_is_enabled() .and. this%hm_base%projector_matrices(imat)%npoints > 0) then
              call accel_write_buffer(this%hm_base%buff_projector_phases, &
                this%hm_base%projector_matrices(imat)%npoints, this%hm_base%projector_phases(1:, imat, ik), offset = offset)
            end if
            offset = offset + this%hm_base%projector_matrices(imat)%npoints
          end do
        end do

      end if

      POP_SUB(hamiltonian_update2.build_phase)
    end subroutine build_phase

  end subroutine hamiltonian_update2


#include "undef.F90"
#include "real.F90"
#include "hamiltonian_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_inc.F90"

end module hamiltonian_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
