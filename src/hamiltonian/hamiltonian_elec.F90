!! Copyright (C) 2002-2020 M. Marques, A. Castro, A. Rubio, G. Bertsch, 
!!                         N. Tancogne-Dejean, M. Lueders
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

module hamiltonian_elec_oct_m
  use accel_oct_m
  use atom_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use boundaries_oct_m
  use boundary_op_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use exchange_operator_oct_m
  use external_potential_oct_m
  use hamiltonian_elec_base_oct_m
  use epot_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_abst_oct_m
  use interaction_oct_m
  use interaction_partner_oct_m
  use ion_electron_local_potential_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use lasers_oct_m
  use lda_u_oct_m
  use linked_list_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use nlcc_oct_m
  use oct_exchange_oct_m
  use parser_oct_m
  use par_vec_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use pcm_oct_m
  use restart_oct_m
  use scissor_oct_m
  use simul_box_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use wfs_elec_oct_m
  use xc_oct_m
  use xc_f03_lib_m

  implicit none

  private
  public ::                          &
    hamiltonian_elec_t,                   &
    hamiltonian_elec_init,                &
    hamiltonian_elec_end,                 &
    dhamiltonian_elec_apply_single,       &
    zhamiltonian_elec_apply_single,       &
    zhamiltonian_elec_apply_all,          &
    dhamiltonian_elec_apply_batch,        &
    zhamiltonian_elec_apply_batch,        &
    dhamiltonian_elec_diagonal,           &
    zhamiltonian_elec_diagonal,           &
    magnus,                               &
    dvmask,                          &
    zvmask,                          &
    hamiltonian_elec_inh_term,            &
    hamiltonian_elec_set_inh,             &
    hamiltonian_elec_remove_inh,          &
    hamiltonian_elec_adjoint,             &
    hamiltonian_elec_not_adjoint,         &
    hamiltonian_elec_epot_generate,       &
    hamiltonian_elec_needs_current,       &
    hamiltonian_elec_update,              &
    hamiltonian_elec_update_pot,          &
    hamiltonian_elec_update2,             &
    hamiltonian_elec_get_time,            &
    hamiltonian_elec_apply_packed,        &
    zhamiltonian_elec_apply_atom,         &
    hamiltonian_elec_dump_vhxc,           &
    hamiltonian_elec_load_vhxc,           &
    hamiltonian_elec_set_vhxc

  type, extends(hamiltonian_abst_t) :: hamiltonian_elec_t
    ! Components are public by default

    !> The Hamiltonian must know what are the "dimensions" of the spaces,
    !! in order to be able to operate on the states.
    type(space_t), private :: space
    type(states_elec_dim_t)  :: d
    type(hamiltonian_elec_base_t) :: hm_base
    type(energy_t), allocatable  :: energy
    type(bc_t)               :: bc      !< boundaries
    FLOAT, allocatable :: vhartree(:) !< Hartree potential
    FLOAT, allocatable :: vxc(:,:)    !< XC potential
    FLOAT, allocatable :: vhxc(:,:)   !< XC potential + Hartree potential + Berry potential
    FLOAT, allocatable :: vtau(:,:)   !< Derivative of e_XC w.r.t. tau
    FLOAT, allocatable :: vberry(:,:) !< Berry phase potential from external E_field

    type(derivatives_t), pointer :: der !< pointer to derivatives
    
    type(geometry_t), pointer :: geo
    FLOAT :: exx_coef !< how much of EXX to mix

    type(poisson_t)          :: psolver      !< Poisson solver
    type(poisson_t), pointer :: psolver_fine !< Poisson solver on the fine grid

    !> The self-induced vector potential and magnetic field
    logical :: self_induced_magnetic
    FLOAT, allocatable :: a_ind(:, :)
    FLOAT, allocatable :: b_ind(:, :)

    integer :: theory_level    !< copied from sys%ks
    type(xc_t), pointer :: xc  !< pointer to xc object

    type(epot_t) :: ep         !< handles the external potential
    type(pcm_t)  :: pcm        !< handles pcm variables
 
    !> absorbing boundaries
    logical, private :: adjoint

    !> Mass of the particle (in most cases, mass = 1, electron mass)
    FLOAT, private :: mass
    !> anisotropic scaling factor for the mass: different along x,y,z etc...
    FLOAT, private :: mass_scaling(MAX_DIM)

    !> There may be an "inhomogeneous", "source", or "forcing" term (useful for the OCT formalism)
    logical, private :: inh_term
    type(states_elec_t) :: inh_st

    !> There may also be a exchange-like term, similar to the one necessary for time-dependent
    !! Hartree Fock, also useful only for the OCT equations
    type(oct_exchange_t) :: oct_exchange

    type(scissor_t) :: scissor

    FLOAT :: current_time
    logical, private :: apply_packed  !< This is initialized by the StatesPack variable.
    
    !> For the LDA+U 
    type(lda_u_t) :: lda_u
    integer       :: lda_u_level

    logical, public :: time_zero

    type(exchange_operator_t), public :: exxop

    type(kpoints_t), pointer, public :: kpoints => null()

    type(partner_list_t) :: external_potentials  !< List with all the external potentials
    FLOAT, allocatable, public  :: v_ext_pot(:)  !< the potential comming from external potentials
    FLOAT, allocatable, public  :: v_static(:)   !< static scalar potential
    
    type(ion_electron_local_potential_t) :: v_ie_loc !< Ion-electron local potential interaction
    type(nlcc_t) :: nlcc !< Ion-electron NLCC interaction

    ! At the moment this is not treated as an external potential
    class(lasers_t), pointer :: ext_lasers => null()      !< lasers 

  contains
    procedure :: update_span => hamiltonian_elec_span
    procedure :: dapply => dhamiltonian_elec_apply
    procedure :: zapply => zhamiltonian_elec_apply
    procedure :: dmagnus_apply => dhamiltonian_elec_magnus_apply
    procedure :: zmagnus_apply => zhamiltonian_elec_magnus_apply
    procedure :: is_hermitian => hamiltonian_elec_hermitian
  end type hamiltonian_elec_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 2, &
    HARTREE               = 1, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4, &
    RDMFT                 = 7

  type(profile_t), save :: prof_hamiltonian, prof_kinetic_start, prof_kinetic_finish
  type(profile_t), save :: prof_exx

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_init(hm, namespace, space, gr, geo, st, theory_level, xc, mc, kpoints, need_exchange)
    type(hamiltonian_elec_t),           target, intent(inout) :: hm
    type(namespace_t),                          intent(in)    :: namespace
    type(space_t),                              intent(in)    :: space
    type(grid_t),                       target, intent(inout) :: gr
    type(geometry_t),                   target, intent(inout) :: geo
    type(states_elec_t),                target, intent(inout) :: st
    integer,                                    intent(in)    :: theory_level
    type(xc_t),                         target, intent(in)    :: xc
    type(multicomm_t),                          intent(in)    :: mc
    type(kpoints_t),                    target, intent(in)    :: kpoints
    logical, optional,                          intent(in)    :: need_exchange

    integer :: iline, icol, il
    integer :: ncols
    type(block_t) :: blk
    type(profile_t), save :: prof

    logical :: need_exchange_
    FLOAT :: rashba_coupling


    PUSH_SUB(hamiltonian_elec_init)
    call profiling_in(prof, 'HAMILTONIAN_ELEC_INIT')
    
    ! make a couple of local copies
    hm%space = space
    hm%theory_level = theory_level
    call states_elec_dim_copy(hm%d, st%d)

    hm%kpoints => kpoints

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
    call parse_variable(namespace, 'ParticleMass', M_ONE, hm%mass)

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
    call parse_variable(namespace, 'RashbaSpinOrbitCoupling', M_ZERO, rashba_coupling, units_inp%energy*units_inp%length)
    if(parse_is_defined(namespace, 'RashbaSpinOrbitCoupling')) then
      if(gr%sb%dim .ne. 2) then
        write(message(1),'(a)') 'Rashba spin-orbit coupling can only be used for two-dimensional systems.'
        call messages_fatal(1, namespace=namespace)
      end if
      call messages_experimental('RashbaSpinOrbitCoupling')
    end if

    call hamiltonian_elec_base_init(hm%hm_base, hm%d%nspin, hm%mass, rashba_coupling)

    ASSERT(associated(gr%der%lapl))
    hm%hm_base%kinetic => gr%der%lapl

    SAFE_ALLOCATE(hm%energy)
    call energy_nullify(hm%energy)

    call oct_exchange_nullify(hm%oct_exchange)

    !Keep pointers to derivatives, geometry and xc
    hm%der => gr%der
    hm%geo => geo
    hm%xc => xc

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = hm%vxc(:, 1), vxc_22 = hm%vxc(:, 2), Re(vxc_12) = hm%vxc(:. 3);
    ! Im(vxc_12) = hm%vxc(:, 4)
    SAFE_ALLOCATE(hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin))
    hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin) = M_ZERO

    if(hm%theory_level /= INDEPENDENT_PARTICLES) then

      SAFE_ALLOCATE(hm%vhartree(1:gr%mesh%np_part))
      hm%vhartree=M_ZERO

      SAFE_ALLOCATE(hm%vxc(1:gr%mesh%np, 1:hm%d%nspin))
      hm%vxc=M_ZERO

      if (family_is_mgga_with_exc(hm%xc)) then
        SAFE_ALLOCATE(hm%vtau(1:gr%mesh%np, 1:hm%d%nspin))
        hm%vtau=M_ZERO
      end if

    end if

    !Initialize Poisson solvers
    call poisson_init(hm%psolver, namespace, space, gr%der, mc, st%qtot)

    if(poisson_is_multigrid(hm%psolver)) then
      SAFE_ALLOCATE(hm%psolver%mgrid)
      call multigrid_init(hm%psolver%mgrid, namespace, space, gr%cv, gr%mesh, gr%der, gr%stencil, mc)
    end if


    nullify(hm%psolver_fine)
    if (gr%have_fine_mesh) then
      SAFE_ALLOCATE(hm%psolver_fine)
      call poisson_init(hm%psolver_fine, namespace, space, gr%fine%der, mc, st%qtot, label = " (fine mesh)")
    else
      hm%psolver_fine => hm%psolver
    end if
  
    ! Initialize external potential
    call epot_init(hm%ep, namespace, gr, hm%geo, hm%psolver, hm%d%ispin, hm%xc%family, mc, hm%kpoints)

    !Temporary construction of the ion-electron interactions
    call hm%v_ie_loc%init(gr%mesh, hm%psolver, hm%geo, namespace)
    if(hm%ep%nlcc) then
      call hm%nlcc%init(gr%mesh, hm%geo)
      SAFE_ALLOCATE(st%rho_core(1:gr%fine%mesh%np))
      st%rho_core(:) = M_ZERO
    end if

    ! Calculate initial value of the gauge vector field
    call gauge_field_init(hm%ep%gfield, namespace, hm%kpoints)

    ! Temporary place for the initialization of the lasers
    hm%ext_lasers => lasers_t(namespace, gr%mesh)
    call lasers_check_symmetries(hm%ext_lasers, kpoints)

    !Static magnetic field or rashba spin-orbit interaction requires complex wavefunctions
    if (parse_is_defined(namespace, 'StaticMagneticField') .or. gauge_field_is_applied(hm%ep%gfield) .or. &
      parse_is_defined(namespace, 'RashbaSpinOrbitCoupling')) then
      call states_set_complex(st)
    end if

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
    call parse_variable(namespace, 'CalculateSelfInducedMagneticField', .false., hm%self_induced_magnetic)
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(hm%a_ind(1:gr%mesh%np_part, 1:gr%sb%dim))
      SAFE_ALLOCATE(hm%b_ind(1:gr%mesh%np_part, 1:gr%sb%dim))

      !(for dim = we could save some memory, but it is better to keep it simple)
    end if

    ! Boundaries
    call bc_init(hm%bc, namespace, gr%mesh, gr%sb, hm%geo)

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
    if(parse_block(namespace, 'MassScaling', blk) == 0) then
      ncols = parse_block_cols(blk, 0)
      if(ncols > gr%sb%dim) then
        call messages_input_error(namespace, "MassScaling")
      end if
      iline = 1 ! just deal with 1 line - should be generalized
      do icol = 1, ncols
        call parse_block_float(blk, iline - 1, icol - 1, hm%mass_scaling(icol))
      end do
      call parse_block_end(blk)
    end if

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
    call parse_variable(namespace, 'DFTULevel', DFT_U_NONE, hm%lda_u_level)
    call messages_print_var_option(stdout,  'DFTULevel', hm%lda_u_level)
    call lda_u_nullify(hm%lda_u)
    if(hm%lda_u_level /= DFT_U_NONE) then
      call messages_experimental('DFT+U')
      call lda_u_init(hm%lda_u, namespace, space, hm%lda_u_level, gr, geo, st, hm%psolver, hm%kpoints)

      !In the present implementation of DFT+U, in case of spinors, we have off-diagonal terms
      !in spin space which break the assumption of the generalized Bloch theorem
      if(kick_get_type(hm%ep%kick) == KICK_MAGNON_MODE .and. gr%der%boundaries%spiral) then
        call messages_not_implemented("DFT+U with generalized Bloch theorem and magnon kick")
      end if 
    end if
 
    if (.not. hm%kpoints%gamma_only()) then
      call init_phase()
    end if
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
    call parse_variable(namespace, 'HamiltonianApplyPacked', .true., hm%apply_packed)

    if(hm%theory_level == HARTREE_FOCK .and. st%parallel_in_states) then
#ifdef HAVE_MPI2
      call messages_experimental('Hartree-Fock parallel in states')
#else
      call messages_write('Hartree-Fock parallel in states required MPI 2')
      call messages_fatal(namespace=namespace)
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
    call parse_variable(namespace, 'TimeZero', .false., hm%time_zero)
    if(hm%time_zero) call messages_experimental('TimeZero')

    call scissor_nullify(hm%scissor)

    !Cam parameters are irrelevant here and are updated later
    call exchange_operator_nullify(hm%exxop)
    need_exchange_ = optional_default(need_exchange, .false.)
    if (hm%theory_level == HARTREE_FOCK .or. hm%theory_level == HARTREE &
          .or. hm%theory_level == RDMFT .or. need_exchange_) then
      call exchange_operator_init(hm%exxop, namespace, space, st, gr%sb, gr%der, mc, hm%kpoints, M_ONE, M_ZERO, M_ZERO)
    end if

    if (hm%apply_packed .and. accel_is_enabled()) then
      ! Check if we can actually apply the hamiltonian packed
      if (gr%mesh%use_curvilinear) then
        if(accel_allow_CPU_only()) then
          hm%apply_packed = .false.
          call messages_write('Cannot use CUDA or OpenCL as curvilinear coordinates are used.')
          call messages_warning(namespace=namespace)
        else
          call messages_write('Cannot use CUDA or OpenCL as curvilinear coordinates are used.', new_line = .true.)
          call messages_write('Calculation will not be continued. To force execution, set AllowCPUonly = yes.' )
          call messages_fatal(namespace=namespace)          
        end if
      end if

      if(hm%bc%abtype == IMAGINARY_ABSORBING) then
        if(accel_allow_CPU_only()) then
          hm%apply_packed = .false.
          call messages_write('Cannot use CUDA or OpenCL as imaginary absorbing boundaries are enabled.')
          call messages_warning(namespace=namespace)
        else
          call messages_write('Cannot use CUDA or OpenCL as imaginary absorbing boundaries are enabled.', new_line = .true.)
          call messages_write('Calculation will not be continued. To force execution, set AllowCPUonly = yes.' )
          call messages_fatal(namespace=namespace)          
        end if
      end if

      if (.not. space%is_periodic()) then
        do il = 1, hm%ext_lasers%no_lasers
          if (laser_kind(hm%ext_lasers%lasers(il)) == E_FIELD_VECTOR_POTENTIAL) then
            if(accel_allow_CPU_only()) then
              hm%apply_packed = .false.
              call messages_write('Cannot use CUDA or OpenCL as a phase is applied to the states.')
              call messages_warning(namespace=namespace)
              exit
            else
              call messages_write('Cannot use CUDA or OpenCL as a phase is applied to the states.', new_line = .true.)
              call messages_write('Calculation will not be continued. To force execution, set AllowCPUonly = yes.' )
              call messages_fatal(namespace=namespace)          
            end if    
          end if
        end do
      end if
    end if

    !We are building the list of external potentials
    !This is done here at the moment, because we pass directly the mesh
    !TODO: Once the abstract Hamiltonian knows about an abstract basis, we might move this to the 
    !      abstract Hamiltonian 
    call load_external_potentials(hm%external_potentials, namespace)

    !Some checks which are electron specific, like k-points
    call external_potentials_checks()

    !At the moment we do only have static external potential, so we never update them
    call build_external_potentials()

    !Build the resulting interactions
    !TODO: This will be moved to the actual interactions
    call build_interactions()

    call profiling_out(prof)
    POP_SUB(hamiltonian_elec_init)

  contains

    ! ---------------------------------------------------------
    subroutine init_phase
      integer :: ip, ik, sp, ip_inner_global
      FLOAT   :: kpoint(1:MAX_DIM), x_global(1:gr%sb%dim)
      

      PUSH_SUB(hamiltonian_elec_init.init_phase)

      SAFE_ALLOCATE(hm%hm_base%phase(1:gr%mesh%np_part, hm%d%kpt%start:hm%d%kpt%end))
      SAFE_ALLOCATE(hm%hm_base%phase_corr(gr%mesh%np+1:gr%mesh%np_part, hm%d%kpt%start:hm%d%kpt%end))
      hm%hm_base%phase_corr = M_ONE

      if(gr%der%boundaries%spiralBC) then
        sp = gr%mesh%np
        if(gr%mesh%parallel_in_domains) sp = gr%mesh%np + gr%mesh%vp%np_ghost

        ! We decided to allocate the array from 1:np_part-sp as this is less error prone when passing 
        ! the array to other routines, or in particular creating a C-style pointer from phase_spiral(1,1).
        ! We will also update phase_corr and possible other similar arrays.

        SAFE_ALLOCATE(hm%hm_base%phase_spiral(1:gr%mesh%np_part-sp, 1:2))

        ! loop over boundary points
        do ip = sp + 1, gr%mesh%np_part
          ! get corresponding inner point
          ip_inner_global = mesh_periodic_point(gr%mesh, space, ip)
          x_global = mesh_x_global(gr%mesh, ip_inner_global)
          hm%hm_base%phase_spiral(ip-sp, 1) = &
            exp(M_zI * sum((gr%mesh%x(ip, 1:gr%sb%dim)-x_global(1:gr%sb%dim)) * gr%der%boundaries%spiral_q(1:gr%sb%dim)))
          hm%hm_base%phase_spiral(ip-sp, 2) = &
            exp(-M_zI * sum((gr%mesh%x(ip, 1:gr%sb%dim)-x_global(1:gr%sb%dim)) * gr%der%boundaries%spiral_q(1:gr%sb%dim)))
        end do

        if(accel_is_enabled()) then
          call accel_create_buffer(hm%hm_base%buff_phase_spiral, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, (gr%mesh%np_part-sp)*2)
          call accel_write_buffer(hm%hm_base%buff_phase_spiral, (gr%mesh%np_part-sp)*2, hm%hm_base%phase_spiral)
        endif
      end if


      kpoint(1:gr%sb%dim) = M_ZERO
      do ik = hm%d%kpt%start, hm%d%kpt%end
        kpoint(1:gr%sb%dim) = hm%kpoints%get_point(hm%d%get_kpoint_index(ik))
        do ip = 1, gr%mesh%np_part
          hm%hm_base%phase(ip, ik) = exp(-M_zI * sum(gr%mesh%x(ip, 1:gr%sb%dim) * kpoint(1:gr%sb%dim)))
        end do

        ! loop over boundary points
        sp = gr%mesh%np
        if(gr%mesh%parallel_in_domains) sp = gr%mesh%np + gr%mesh%vp%np_ghost
        do ip = sp + 1, gr%mesh%np_part
          ! get corresponding inner point
          ip_inner_global = mesh_periodic_point(gr%mesh, space, ip)

          ! compute phase correction from global coordinate (opposite sign!)
          x_global = mesh_x_global(gr%mesh, ip_inner_global)
          hm%hm_base%phase_corr(ip, ik) = hm%hm_base%phase(ip, ik)* &
            exp(M_zI * sum(x_global(1:gr%sb%dim) * kpoint(1:gr%sb%dim)))
        end do
      end do

      if(accel_is_enabled()) then
        call accel_create_buffer(hm%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, gr%mesh%np_part*hm%d%kpt%nlocal)
        call accel_write_buffer(hm%hm_base%buff_phase, gr%mesh%np_part*hm%d%kpt%nlocal, hm%hm_base%phase)
        hm%hm_base%buff_phase_qn_start = hm%d%kpt%start
      end if

      ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
      if(hm%lda_u_level /= DFT_U_NONE) then
        call lda_u_build_phase_correction(hm%lda_u, gr%sb%dim, hm%d, gr%der%boundaries, namespace, hm%kpoints)
      end if

      POP_SUB(hamiltonian_elec_init.init_phase)
    end subroutine init_phase

    ! ---------------------------------------------------------
    subroutine build_external_potentials()
      type(list_iterator_t) :: iter
      class(*), pointer :: potential

      PUSH_SUB(hamiltonian_elec_init.build_external_potentials)

      SAFE_ALLOCATE(hm%v_ext_pot(1:gr%mesh%np)) 
      hm%v_ext_pot(1:gr%mesh%np) = M_ZERO

      call iter%start(hm%external_potentials)
      do while (iter%has_next())
        potential => iter%get_next() 
        select type (potential)
        class is (external_potential_t)

          call potential%allocate_memory(gr%mesh)
          call potential%calculate(namespace, gr%mesh, hm%psolver)
          !To preserve the old behavior, we are adding the various potentials
          !to the corresponding arrays
          select case(potential%type)
          case(EXTERNAL_POT_USDEF, EXTERNAL_POT_FROM_FILE, EXTERNAL_POT_CHARGE_DENSITY)
            call lalg_axpy(gr%mesh%np, M_ONE, potential%pot, hm%v_ext_pot)

          case(EXTERNAL_POT_STATIC_BFIELD)
            if (.not. allocated(hm%ep%B_field)) then
              SAFE_ALLOCATE(hm%ep%B_field(1:3)) !Cannot be gr%sb%dim
              hm%ep%B_field(1:3) = M_ZERO
            end if
            hm%ep%B_field(1:3) = hm%ep%B_field(1:3) + potential%B_field(1:3)
            
            if (.not. allocated(hm%ep%A_static)) then
              SAFE_ALLOCATE(hm%ep%A_static(1:gr%mesh%np, 1:gr%sb%dim))
              hm%ep%A_static(1:gr%mesh%np, 1:gr%sb%dim) = M_ZERO
            end if
            call lalg_axpy(gr%mesh%np, gr%sb%dim, M_ONE, potential%A_static, hm%ep%A_static)

          case(EXTERNAL_POT_STATIC_EFIELD)
            if (.not. allocated(hm%ep%E_field)) then
              SAFE_ALLOCATE(hm%ep%E_field(1:gr%sb%dim))
              hm%ep%E_field(1:gr%sb%dim) = M_ZERO
            end if
            hm%ep%E_field(1:gr%sb%dim) = hm%ep%E_field(1:gr%sb%dim) + potential%E_field(1:gr%sb%dim)

            !In the fully periodic case, we use Berry phases
            if (space%periodic_dim < space%dim) then
              if (.not. allocated(hm%v_static)) then
                SAFE_ALLOCATE(hm%v_static(1:gr%mesh%np))
                hm%v_static(1:gr%mesh%np) = M_ZERO
              end if
              if (.not. allocated(hm%ep%v_ext)) then
                SAFE_ALLOCATE(hm%ep%v_ext(1:gr%mesh%np_part))
                hm%ep%v_ext(1:gr%mesh%np_part) = M_ZERO
              end if     
              call lalg_axpy(gr%mesh%np, M_ONE, potential%pot, hm%v_static)
              call lalg_axpy(gr%mesh%np, M_ONE, potential%v_ext, hm%ep%v_ext)
            end if
          end select
          call potential%deallocate_memory()

        class default
          ASSERT(.false.)
        end select
      end do

      POP_SUB(hamiltonian_elec_init.build_external_potentials)
    end subroutine build_external_potentials

    ! ---------------------------------------------------------
    subroutine external_potentials_checks()
      type(list_iterator_t) :: iter
      class(*), pointer :: potential

      PUSH_SUB(hamiltonian_elec_init.external_potentials_checks)

      call iter%start(hm%external_potentials)
      do while (iter%has_next())
        potential => iter%get_next()
        select type (potential)
        class is (external_potential_t)

          if(potential%type == EXTERNAL_POT_STATIC_EFIELD .and. hm%d%nik > 1) then
            message(1) = "Applying StaticElectricField in a periodic direction is only accurate for large supercells."
            message(2) = "Single-point Berry phase is not appropriate when k-point sampling is needed."
            call messages_warning(2, namespace=namespace)
          end if

        class default
          ASSERT(.false.)
        end select
      end do

      POP_SUB(hamiltonian_elec_init.external_potentials_checks)
    end subroutine external_potentials_checks


    !The code in this routines needs to know about the external potentials.
    !This will be treated in the future by the interactions directly.
    subroutine build_interactions()
      logical :: external_potentials_present
      logical :: kick_present

      PUSH_SUB(hamiltonian_elec_init.build_interactions)      

      if (allocated(hm%ep%E_field) .and. space%is_periodic() .and. .not. gauge_field_is_applied(hm%ep%gfield)) then
        ! only need vberry if there is a field in a periodic direction
        ! and we are not setting a gauge field
        if (any(abs(hm%ep%E_field(1:space%periodic_dim)) > M_EPSILON)) then
          SAFE_ALLOCATE(hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
          hm%vberry = M_ZERO
        end if
      end if

      external_potentials_present = epot_have_external_potentials(hm%ep) .or. (hm%ext_lasers%no_lasers > 0) &
                                    .or. allocated(hm%v_static) 
      

      kick_present = epot_have_kick(hm%ep)

      call pcm_init(hm%pcm, namespace, geo, gr, st%qtot, st%val_charge, external_potentials_present, kick_present )  !< initializes PCM
      if (hm%pcm%run_pcm) then
        if (hm%theory_level /= KOHN_SHAM_DFT) call messages_not_implemented("PCM for TheoryLevel /= DFT", namespace=namespace)
        if (gr%have_fine_mesh) call messages_not_implemented("PCM with UseFineMesh", namespace=namespace)
      end if

      POP_SUB(hamiltonian_elec_init.build_interactions)

    end subroutine build_interactions


  end subroutine hamiltonian_elec_init

  
  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_end(hm)
    type(hamiltonian_elec_t), target, intent(inout) :: hm

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: potential


    PUSH_SUB(hamiltonian_elec_end)

    call hamiltonian_elec_base_end(hm%hm_base)

    if (allocated(hm%hm_base%phase) .and. accel_is_enabled()) then
      call accel_release_buffer(hm%hm_base%buff_phase)
    end if

    if(allocated(hm%hm_base%phase_spiral) .and. accel_is_enabled()) then
      call accel_release_buffer(hm%hm_base%buff_phase_spiral)
    end if

    SAFE_DEALLOCATE_A(hm%hm_base%phase)
    SAFE_DEALLOCATE_A(hm%hm_base%phase_corr)
    SAFE_DEALLOCATE_A(hm%hm_base%phase_spiral)
    SAFE_DEALLOCATE_A(hm%vhartree)
    SAFE_DEALLOCATE_A(hm%vhxc)
    SAFE_DEALLOCATE_A(hm%vxc)
    SAFE_DEALLOCATE_A(hm%vberry)
    SAFE_DEALLOCATE_A(hm%a_ind)
    SAFE_DEALLOCATE_A(hm%b_ind)
    SAFE_DEALLOCATE_A(hm%v_ext_pot)
    
    if (family_is_mgga_with_exc(hm%xc)) then
      SAFE_DEALLOCATE_A(hm%vtau)
    end if

    if (associated(hm%psolver_fine, hm%psolver)) then
      nullify(hm%psolver_fine)
    else
      call poisson_end(hm%psolver_fine)
      SAFE_DEALLOCATE_P(hm%psolver_fine)
    end if
    call poisson_end(hm%psolver)

    nullify(hm%xc)

    ! No call to safe_deallocate macro here, as it gives an ICE with gfortran
    if (associated(hm%ext_lasers)) then
      deallocate(hm%ext_lasers)
    end if


    call epot_end(hm%ep)
    nullify(hm%geo)

    call bc_end(hm%bc)

    call states_elec_dim_end(hm%d) 

    if(hm%scissor%apply) call scissor_end(hm%scissor)

    call exchange_operator_end(hm%exxop)
    call lda_u_end(hm%lda_u)

    SAFE_DEALLOCATE_A(hm%energy)

    if (hm%pcm%run_pcm) call pcm_end(hm%pcm)

    call hm%v_ie_loc%end()   
    call hm%nlcc%end()

    call iter%start(hm%external_potentials)
    do while (iter%has_next())
      potential => iter%get_next()
      SAFE_DEALLOCATE_P(potential)
    end do
    call hm%external_potentials%empty()
    SAFE_DEALLOCATE_A(hm%v_static)

    POP_SUB(hamiltonian_elec_end)
  end subroutine hamiltonian_elec_end


  ! ---------------------------------------------------------
  ! True if the Hamiltonian is Hermitian, false otherwise
  logical function hamiltonian_elec_hermitian(hm)
    class(hamiltonian_elec_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_elec_hermitian)
    hamiltonian_elec_hermitian = .not.((hm%bc%abtype == IMAGINARY_ABSORBING) .or. &
                                  oct_exchange_enabled(hm%oct_exchange))

    POP_SUB(hamiltonian_elec_hermitian)
  end function hamiltonian_elec_hermitian


  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_span(hm, delta, emin)
    class(hamiltonian_elec_t), intent(inout) :: hm
    FLOAT,                     intent(in)    :: delta, emin

    PUSH_SUB(hamiltonian_elec_span)

    hm%spectral_middle_point = ((M_PI**2 / (2 * delta**2)) + emin) / M_TWO
    hm%spectral_half_span    = ((M_PI**2 / (2 * delta**2)) - emin) / M_TWO

    POP_SUB(hamiltonian_elec_span)
  end subroutine hamiltonian_elec_span


  ! ---------------------------------------------------------
  pure logical function hamiltonian_elec_inh_term(hm) result(inh)
    type(hamiltonian_elec_t), intent(in) :: hm

    inh = hm%inh_term
  end function hamiltonian_elec_inh_term


  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_set_inh(hm, st)
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(states_elec_t), intent(in)    :: st

    PUSH_SUB(hamiltonian_elec_set_inh)

    if(hm%inh_term) call states_elec_end(hm%inh_st)
    call states_elec_copy(hm%inh_st, st)
    hm%inh_term = .true.

    POP_SUB(hamiltonian_elec_set_inh)
  end subroutine hamiltonian_elec_set_inh


  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_remove_inh(hm)
    type(hamiltonian_elec_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_elec_remove_inh)

    if(hm%inh_term) then
      call states_elec_end(hm%inh_st)
      hm%inh_term = .false.
    end if

    POP_SUB(hamiltonian_elec_remove_inh)
  end subroutine hamiltonian_elec_remove_inh

  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_adjoint(hm)
    type(hamiltonian_elec_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_elec_adjoint)

    if(.not.hm%adjoint) then
      hm%adjoint = .true.
      if(hm%bc%abtype == IMAGINARY_ABSORBING) then
        hm%bc%mf = -hm%bc%mf
      end if
    end if

    POP_SUB(hamiltonian_elec_adjoint)
  end subroutine hamiltonian_elec_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_not_adjoint(hm)
    type(hamiltonian_elec_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_elec_not_adjoint)

    if(hm%adjoint) then
      hm%adjoint = .false.
      if(hm%bc%abtype == IMAGINARY_ABSORBING) then
        hm%bc%mf = -hm%bc%mf
      end if
    end if

    POP_SUB(hamiltonian_elec_not_adjoint)
  end subroutine hamiltonian_elec_not_adjoint


  ! ---------------------------------------------------------
  !> (re-)build the Hamiltonian for the next application:
  subroutine hamiltonian_elec_update(this, mesh, namespace, time)
    type(hamiltonian_elec_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    type(namespace_t),        intent(in)    :: namespace
    FLOAT, optional,          intent(in)    :: time

    integer :: ispin, ip, idir, iatom, ilaser
    type(profile_t), save :: prof, prof_phases
    FLOAT :: aa(1:MAX_DIM), time_
    FLOAT, allocatable :: vp(:,:)

    PUSH_SUB(hamiltonian_elec_update)
    call profiling_in(prof, "HAMILTONIAN_ELEC_UPDATE")

    this%current_time = M_ZERO
    if(present(time)) this%current_time = time

    time_ = optional_default(time, CNST(0.0))

    ! set everything to zero
    call hamiltonian_elec_base_clear(this%hm_base)

    ! the xc, hartree and external potentials
    call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_POTENTIAL, &
      complex_potential = this%bc%abtype == IMAGINARY_ABSORBING)

    call hamiltonian_elec_update_pot(this, mesh, accel_copy=.false.)

    ! the lasers
    if (present(time) .or. this%time_zero) then

      do ilaser = 1, this%ext_lasers%no_lasers
        select case(laser_kind(this%ext_lasers%lasers(ilaser)))
        case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
          do ispin = 1, this%d%spin_channels
            call laser_potential(this%ext_lasers%lasers(ilaser), mesh,  &
                                    this%hm_base%potential(:, ispin), time_)
          end do
        case(E_FIELD_MAGNETIC)
          call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL + FIELD_UNIFORM_MAGNETIC_FIELD, &
            .false.)
          ! get the vector potential
          SAFE_ALLOCATE(vp(1:mesh%np, 1:mesh%sb%dim))
          vp(1:mesh%np, 1:mesh%sb%dim) = M_ZERO
          call laser_vector_potential(this%ext_lasers%lasers(ilaser), mesh, vp, time_)
          !$omp parallel do schedule(static)
          do ip = 1, mesh%np
            do idir = 1, mesh%sb%dim
              this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) - vp(ip, idir)/P_C
            end do
          end do
          ! and the magnetic field
          call laser_field(this%ext_lasers%lasers(ilaser), this%hm_base%uniform_magnetic_field(1:mesh%sb%dim), time_)
          SAFE_DEALLOCATE_A(vp)
        case(E_FIELD_VECTOR_POTENTIAL)
          call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, .false.)
          ! get the uniform vector potential associated with a magnetic field
          aa = M_ZERO
          call laser_field(this%ext_lasers%lasers(ilaser), aa(1:mesh%sb%dim), time_)
          this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim) &
            - aa(1:mesh%sb%dim)/P_C
        end select
      end do

      ! the gauge field
      if(gauge_field_is_applied(this%ep%gfield)) then
        call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, .false.)
        call gauge_field_get_vec_pot(this%ep%gfield, aa)
        this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim)  &
          - aa(1:mesh%sb%dim)/P_c
      end if

      ! the electric field for a periodic system through the gauge field
      if (allocated(this%ep%e_field) .and. gauge_field_is_applied(this%ep%gfield)) then
        this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) = &
          this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) - time_*this%ep%e_field(1:mesh%sb%periodic_dim)
      end if
      
    end if

    ! the vector potential of a static magnetic field
    if (allocated(this%ep%a_static)) then
      call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL, .false.)
      !$omp parallel do schedule(static)
      do ip = 1, mesh%np
        do idir = 1, mesh%sb%dim
          this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) + this%ep%a_static(ip, idir)
        end do
      end do
    end if

    !The electric field was added to the KS potential
    call hamiltonian_elec_base_accel_copy_pot(this%hm_base, mesh)

    ! and the static magnetic field
    if (allocated(this%ep%b_field)) then
      call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_MAGNETIC_FIELD, .false.)
      do idir = 1, 3
        this%hm_base%uniform_magnetic_field(idir) = this%hm_base%uniform_magnetic_field(idir) + this%ep%b_field(idir)
      end do
    end if

    call hamiltonian_elec_base_update(this%hm_base, mesh)

    call build_phase()

    call profiling_out(prof)
    POP_SUB(hamiltonian_elec_update)

  contains

    subroutine build_phase()
      integer :: ik, imat, nmat, max_npoints, offset
      integer :: ip, ip_inner_global, sp
      FLOAT   :: kpoint(1:MAX_DIM), x_global(1:mesh%sb%dim)
      integer :: iphase, nphase

      PUSH_SUB(hamiltonian_elec_update.build_phase)

      if ((.not. this%kpoints%gamma_only()) .or. allocated(this%hm_base%uniform_vector_potential)) then

        call profiling_in(prof_phases, 'UPDATE_PHASES')
        ! now regenerate the phases for the pseudopotentials
        do iatom = 1, this%ep%natoms
          call projector_init_phases(this%ep%proj(iatom), mesh%sb%dim, this%d, this%der%boundaries, this%kpoints, &
            vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end do

        call profiling_out(prof_phases)
      end if

      if(allocated(this%hm_base%uniform_vector_potential)) then

        if(.not. allocated(this%hm_base%phase)) then
          SAFE_ALLOCATE(this%hm_base%phase(1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part*this%d%kpt%nlocal)
          end if
        end if

        if(.not. allocated(this%hm_base%phase_corr)) then
          SAFE_ALLOCATE(this%hm_base%phase_corr(mesh%np+1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
        end if

        kpoint(1:mesh%sb%dim) = M_ZERO
        do ik = this%d%kpt%start, this%d%kpt%end
          kpoint(1:mesh%sb%dim) = this%kpoints%get_point(this%d%get_kpoint_index(ik))
          !We add the vector potential
          kpoint(1:mesh%sb%dim) = kpoint(1:mesh%sb%dim) + this%hm_base%uniform_vector_potential(1:mesh%sb%dim)

          !$omp parallel do schedule(static)
          do ip = 1, mesh%np_part
            this%hm_base%phase(ip, ik) = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*kpoint(1:mesh%sb%dim)))
          end do

          ! loop over boundary points
          sp = mesh%np
          ! skip ghost points
          if(mesh%parallel_in_domains) sp = mesh%np + mesh%vp%np_ghost
          !$omp parallel do schedule(static) private(ip_inner_global, x_global)
          do ip = sp + 1, mesh%np_part
            ! get corresponding inner point
            ip_inner_global = mesh_periodic_point(mesh, this%space, ip)

            ! compute phase correction from global coordinate (opposite sign!)
            x_global = mesh_x_global(mesh, ip_inner_global)

            this%hm_base%phase_corr(ip, ik) = M_zI * sum(x_global(1:mesh%sb%dim) * kpoint(1:mesh%sb%dim))
            this%hm_base%phase_corr(ip, ik) = exp(this%hm_base%phase_corr(ip, ik))*this%hm_base%phase(ip, ik)
          end do
          
        end do

        if(accel_is_enabled()) then
          call accel_write_buffer(this%hm_base%buff_phase, mesh%np_part*this%d%kpt%nlocal, this%hm_base%phase)
        end if

        ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
        if(this%lda_u_level /= DFT_U_NONE) then
          call lda_u_build_phase_correction(this%lda_u, mesh%sb%dim, this%d, this%der%boundaries, namespace, this%kpoints,  &
               vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end if
      end if

      max_npoints = this%hm_base%max_npoints
      nmat = this%hm_base%nprojector_matrices


      if (allocated(this%hm_base%phase) .and. allocated(this%hm_base%projector_matrices)) then

        nphase = 1
        if(this%der%boundaries%spiralBC) nphase = 3

        if(.not. allocated(this%hm_base%projector_phases)) then
          SAFE_ALLOCATE(this%hm_base%projector_phases(1:max_npoints, 1:nphase, nmat, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_projector_phases, ACCEL_MEM_READ_ONLY, &
              TYPE_CMPLX, this%hm_base%total_points*nphase*this%d%kpt%nlocal)
            ! We need to save nphase, with which the array has been build, 
            ! as the number might change throughout the run
            this%hm_base%nphase = nphase
          end if
        end if

        offset = 0
        do ik = this%d%kpt%start, this%d%kpt%end
          do imat = 1, this%hm_base%nprojector_matrices
            iatom = this%hm_base%projector_to_atom(imat)
            do iphase = 1, nphase
              !$omp parallel do schedule(static)
              do ip = 1, this%hm_base%projector_matrices(imat)%npoints
                this%hm_base%projector_phases(ip, iphase, imat, ik) = this%ep%proj(iatom)%phase(ip, iphase, ik)
              end do

              if(accel_is_enabled() .and. this%hm_base%projector_matrices(imat)%npoints > 0) then
                call accel_write_buffer(this%hm_base%buff_projector_phases, &
                  this%hm_base%projector_matrices(imat)%npoints, this%hm_base%projector_phases(1:, iphase, imat, ik), &
                  offset = offset)
              end if
              offset = offset + this%hm_base%projector_matrices(imat)%npoints
            end do
          end do
        end do

      end if

      POP_SUB(hamiltonian_elec_update.build_phase)
    end subroutine build_phase

  end subroutine hamiltonian_elec_update


  !----------------------------------------------------------------
  ! Update the KS potential of the electronic Hamiltonian
  subroutine hamiltonian_elec_update_pot(this, mesh, accel_copy)
    type(hamiltonian_elec_t), intent(inout) :: this
    type(mesh_t),             intent(in)    :: mesh
    logical,                  intent(in)    :: accel_copy

    integer :: ispin, ip

    PUSH_SUB(hamiltonian_elec_update_pot)

    do ispin = 1, this%d%nspin
      if(ispin <= 2) then
        !$omp parallel do simd schedule(static)
        do ip = 1, mesh%np
          this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin) + this%ep%vpsl(ip) + this%v_ext_pot(ip)
        end do

        !> Adds PCM contributions
        if (this%pcm%run_pcm) then
          if (this%pcm%solute) then
            !$omp parallel do simd schedule(static)
            do ip = 1, mesh%np
              this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + &
                this%pcm%v_e_rs(ip) + this%pcm%v_n_rs(ip)
            end do
          end if
          if (this%pcm%localf) then
            !$omp parallel do simd schedule(static)
            do ip = 1, mesh%np
              this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + &
                this%pcm%v_ext_rs(ip)
            end do
          end if 
        end if

        !> Adds possible absorbing potential
        if(this%bc%abtype == IMAGINARY_ABSORBING) then
          !$omp parallel do simd schedule(static)
          do ip = 1, mesh%np
            this%hm_base%Impotential(ip, ispin) = this%hm_base%Impotential(ip, ispin) + this%bc%mf(ip)
          end do
        end if

      else !Spinors 
        !$omp parallel do simd schedule(static)
        do ip = 1, mesh%np
          this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin)
        end do
          
      end if

    end do

    if(accel_copy) then
      call hamiltonian_elec_base_accel_copy_pot(this%hm_base, mesh)
    end if

    POP_SUB(hamiltonian_elec_update_pot)

  end subroutine hamiltonian_elec_update_pot

  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_epot_generate(this, namespace, gr, geo, st, time)
    type(hamiltonian_elec_t), intent(inout) :: this
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(geometry_t), target, intent(inout) :: geo
    type(states_elec_t),      intent(inout) :: st
    FLOAT,          optional, intent(in)    :: time

    PUSH_SUB(hamiltonian_elec_epot_generate)

    this%geo => geo
    call epot_generate(this%ep, namespace, gr, this%geo, this%d)

    ! Interation terms are treated below

    ! First we add the static electric fields
    if (this%ep%classical_pot > 0) then
      call lalg_axpy(gr%mesh%np, M_ONE, this%ep%Vclassical, this%ep%vpsl)
    end if

    if (allocated(this%ep%e_field) .and. gr%sb%periodic_dim < gr%sb%dim) then
      call lalg_axpy(gr%mesh%np, M_ONE, this%v_static, this%ep%vpsl)
    end if

    ! Here we need to pass this again, else test are failing.
    ! This is not a real problem, as the multisystem framework will indeed to this anyway
    this%v_ie_loc%atoms_dist => geo%atoms_dist
    this%v_ie_loc%atom => geo%atom
    call this%v_ie_loc%calculate()

    ! At the moment we need to add this to ep%vpsl, to keep the behavior of the code
    call lalg_axpy(gr%mesh%np, M_ONE, this%v_ie_loc%potential(:,1), this%ep%vpsl)

    ! Here we need to pass this again, else test are failing.
    ! This is not a real problem, as the multisystem framework will indeed to this anyway
    if(this%ep%nlcc) then
      this%nlcc%atoms_dist => geo%atoms_dist
      this%nlcc%atom => geo%atom
      call this%nlcc%calculate()
      call lalg_copy(gr%mesh%np, this%nlcc%density(:,1), st%rho_core)
    end if

    call hamiltonian_elec_base_build_proj(this%hm_base, gr%mesh, this%ep)
    call hamiltonian_elec_update(this, gr%mesh, namespace, time)

    ! Check if projectors are still compatible with apply_packed on GPU
    if (this%apply_packed .and. accel_is_enabled()) then
      if (this%ep%non_local .and. .not. this%hm_base%apply_projector_matrices) then
        if(accel_allow_CPU_only()) then
          this%apply_packed = .false.
          call messages_write('Cannot use CUDA or OpenCL as relativistic pseudopotentials are used.')
          call messages_warning(namespace=namespace)
        else
          call messages_write('Cannot use CUDA or OpenCL as relativistic pseudopotentials are used.', new_line = .true.)
          call messages_write('Calculation will not be continued. To force execution, set AllowCPUonly = yes.' )
          call messages_fatal(namespace=namespace)          
       end if
      end if

      if (hamiltonian_elec_base_projector_self_overlap(this%hm_base)) then
        if(accel_allow_CPU_only()) then
          this%apply_packed = .false.
          call messages_write('Cannot use CUDA or OpenCL as some pseudopotentials overlap with themselves.')
          call messages_warning(namespace=namespace)
        else
          call messages_write('Cannot use CUDA or OpenCL as some pseudopotentials overlap with themselves.', new_line = .true.)
          call messages_write('Calculation will not be continued. To force execution, set AllowCPUonly = yes.' )
          call messages_fatal(namespace=namespace)          
        end if
      end if
    end if

    if (this%pcm%run_pcm) then
     !> Generates the real-space PCM potential due to nuclei which do not change
     !! during the SCF calculation.
     if (this%pcm%solute) &
       call pcm_calc_pot_rs(this%pcm, gr%mesh, this%psolver, geo = geo)

      !> Local field effects due to static electrostatic potentials (if they were).
      !! The laser and the kick are included in subroutine v_ks_hartree (module v_ks).
      !  Interpolation is needed, hence gr%mesh%np_part -> 1:gr%mesh%np
      if (this%pcm%localf .and. allocated(this%v_static)) then
        call pcm_calc_pot_rs(this%pcm, gr%mesh, this%psolver, v_ext = this%ep%v_ext(1:gr%mesh%np_part))
      end if

    end if

    call lda_u_update_basis(this%lda_u, gr, geo, st, this%psolver, namespace, this%kpoints, &
                                allocated(this%hm_base%phase))

    POP_SUB(hamiltonian_elec_epot_generate)
  end subroutine hamiltonian_elec_epot_generate

  ! -----------------------------------------------------------------

  FLOAT function hamiltonian_elec_get_time(this) result(time)
    type(hamiltonian_elec_t),   intent(inout) :: this

    time = this%current_time
  end function hamiltonian_elec_get_time

  ! -----------------------------------------------------------------

  logical function hamiltonian_elec_apply_packed(this) result(apply)
    type(hamiltonian_elec_t),   intent(in) :: this

    apply = this%apply_packed

  end function hamiltonian_elec_apply_packed


  ! -----------------------------------------------------------------
  subroutine zhamiltonian_elec_apply_atom (hm, namespace, space, atom, mesh, ia, psi, vpsi)
    type(hamiltonian_elec_t), intent(in)  :: hm
    type(namespace_t),        intent(in)  :: namespace
    type(space_t),            intent(in)  :: space
    type(atom_t),             intent(in)  :: atom
    type(mesh_t),             intent(in)  :: mesh
    integer,                  intent(in)  :: ia
    CMPLX,                    intent(in)  :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
    CMPLX,                    intent(out) :: vpsi(:,:) !< (gr%mesh%np, hm%d%dim)

    integer :: idim
    FLOAT, allocatable :: vlocal(:)
    PUSH_SUB(zhamiltonian_elec_apply_atom)

    SAFE_ALLOCATE(vlocal(1:mesh%np_part))
    vlocal = M_ZERO
    call epot_local_potential(hm%ep, namespace, space, mesh, atom, ia, vlocal)

    do idim = 1, hm%d%dim
      vpsi(1:mesh%np, idim) = vlocal(1:mesh%np) * psi(1:mesh%np, idim)
    end do

    SAFE_DEALLOCATE_A(vlocal)
    POP_SUB(zhamiltonian_elec_apply_atom)
  end subroutine zhamiltonian_elec_apply_atom


  ! -----------------------------------------------------------------
  subroutine hamiltonian_elec_dump_vhxc(restart, hm, mesh, ierr)
    type(restart_t),     intent(in)  :: restart
    type(hamiltonian_elec_t), intent(in)  :: hm
    type(mesh_t),        intent(in)  :: mesh
    integer,             intent(out) :: ierr

    integer :: iunit, err, err2(2), isp
    character(len=12) :: filename
    character(len=100) :: lines(2)

    PUSH_SUB(hamiltonian_elec_dump_vhxc)

    ierr = 0

    if (restart_skip(restart) .or. hm%theory_level == INDEPENDENT_PARTICLES) then
      POP_SUB(hamiltonian_elec_dump_vhxc)
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

      call drestart_write_mesh_function(restart, filename, mesh, hm%vhxc(:,isp), err)
      if (err /= 0) err2(2) = err2(2) + 1

    end do
    if (err2(1) /= 0) ierr = ierr + 2
    if (err2(2) /= 0) ierr = ierr + 4

    lines(1) = '%'
    call restart_write(restart, iunit, lines, 1, err)
    if (err /= 0) ierr = ierr + 4

    ! MGGAs and hybrid MGGAs have an extra term that also needs to be dumped
    if (family_is_mgga_with_exc(hm%xc)) then
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

        call drestart_write_mesh_function(restart, filename, mesh, hm%vtau(:,isp), err)
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

    POP_SUB(hamiltonian_elec_dump_vhxc)
  end subroutine hamiltonian_elec_dump_vhxc


  ! ---------------------------------------------------------
  subroutine hamiltonian_elec_load_vhxc(restart, hm, mesh, ierr)
    type(restart_t),     intent(in)    :: restart
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(mesh_t),        intent(in)    :: mesh
    integer,             intent(out)   :: ierr

    integer :: err, err2, isp
    character(len=12) :: filename

    PUSH_SUB(hamiltonian_elec_load_vhxc)

    ierr = 0

    if (restart_skip(restart) .or. hm%theory_level == INDEPENDENT_PARTICLES) then
      ierr = -1
      POP_SUB(hamiltonian_elec_load_vhxc)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading Vhxc restart."
      call messages_info(1)
    end if

    err2 = 0
    do isp = 1, hm%d%nspin
      if (hm%d%nspin==1) then
        write(filename, fmt='(a)') 'vhxc'
      else
        write(filename, fmt='(a,i1)') 'vhxc-sp', isp
      end if

      call drestart_read_mesh_function(restart, filename, mesh, hm%vhxc(:,isp), err)
      if (err /= 0) err2 = err2 + 1

    end do
    if (err2 /= 0) ierr = ierr + 1

    ! MGGAs and hybrid MGGAs have an extra term that also needs to be read
    err2 = 0
    if (family_is_mgga_with_exc(hm%xc)) then
      do isp = 1, hm%d%nspin
        if (hm%d%nspin == 1) then
          write(filename, fmt='(a)') 'vtau'
        else
          write(filename, fmt='(a,i1)') 'vtau-sp', isp
        end if

        call drestart_read_mesh_function(restart, filename, mesh, hm%vtau(:,isp), err)
        if (err /= 0) err2 = err2 + 1

      end do

      if (err2 /= 0) ierr = ierr + 2
    end if

    if (debug%info) then
      message(1) = "Debug: Reading Vhxc restart done."
      call messages_info(1)
    end if

    POP_SUB(hamiltonian_elec_load_vhxc)
  end subroutine hamiltonian_elec_load_vhxc

  ! ---------------------------------------------------------
  ! This is an extension of "hamiltonian_elec_update2" to be used by the
  ! CFM4 propagator. It updates the Hamiltonian by considering a
  ! weighted sum of the external potentials at times time(1) and time(2),
  ! weighted by alpha(1) and alpha(2).
  subroutine hamiltonian_elec_update2(this, mesh, time, mu)
    type(hamiltonian_elec_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: time(1:2)
    FLOAT,               intent(in)    :: mu(1:2)

    integer :: ispin, ip, idir, iatom, ilaser, itime
    type(profile_t), save :: prof, prof_phases
    FLOAT :: aa(1:MAX_DIM), time_
    FLOAT, allocatable :: vp(:,:)

    FLOAT, allocatable :: velectric(:)

    PUSH_SUB(hamiltonian_elec_update2)
    call profiling_in(prof, "HAMILTONIAN_ELEC_UPDATE")

    this%current_time = M_ZERO
    this%current_time = time(1)

    ! set everything to zero
    call hamiltonian_elec_base_clear(this%hm_base)

    ! the xc, hartree and external potentials
    call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_POTENTIAL, &
      complex_potential = this%bc%abtype == IMAGINARY_ABSORBING)

    call hamiltonian_elec_update_pot(this, mesh, accel_copy=.false.)

    do itime = 1, 2
      time_ = time(itime)

      do ilaser = 1, this%ext_lasers%no_lasers
        select case(laser_kind(this%ext_lasers%lasers(ilaser)))
        case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
          SAFE_ALLOCATE(velectric(1:mesh%np))
          do ispin = 1, this%d%spin_channels
            velectric = M_ZERO
            call laser_potential(this%ext_lasers%lasers(ilaser), mesh,  velectric, time_)
            !$omp parallel do simd schedule(static)
            do ip = 1, mesh%np
              this%hm_base%potential(ip, ispin) = this%hm_base%potential(ip, ispin) + mu(itime) * velectric(ip)
            end do
          end do
          SAFE_DEALLOCATE_A(velectric)
        case(E_FIELD_MAGNETIC)
          call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL + FIELD_UNIFORM_MAGNETIC_FIELD, .false.)
          ! get the vector potential
          SAFE_ALLOCATE(vp(1:mesh%np, 1:mesh%sb%dim))
          vp(1:mesh%np, 1:mesh%sb%dim) = M_ZERO
          call laser_vector_potential(this%ext_lasers%lasers(ilaser), mesh, vp, time_)
          do idir = 1, mesh%sb%dim
            !$omp parallel do schedule(static)
            do ip = 1, mesh%np
              this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) &
                - mu(itime) * vp(ip, idir)/P_C
            end do
          end do
          ! and the magnetic field
          call laser_field(this%ext_lasers%lasers(ilaser), this%hm_base%uniform_magnetic_field(1:mesh%sb%dim), time_)
          SAFE_DEALLOCATE_A(vp)
        case(E_FIELD_VECTOR_POTENTIAL)
          call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, .false.)
          ! get the uniform vector potential associated with a magnetic field
          aa = M_ZERO
          call laser_field(this%ext_lasers%lasers(ilaser), aa(1:mesh%sb%dim), time_)
          this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim) &
            - mu(itime) * aa(1:mesh%sb%dim)/P_C
        end select
      end do

      ! the gauge field
      if(gauge_field_is_applied(this%ep%gfield)) then
        call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL, .false.)
        call gauge_field_get_vec_pot(this%ep%gfield, aa)
        this%hm_base%uniform_vector_potential(1:mesh%sb%dim) = this%hm_base%uniform_vector_potential(1:mesh%sb%dim)  &
          - aa(1:mesh%sb%dim)/P_c
      end if

      ! the electric field for a periodic system through the gauge field
      if (allocated(this%ep%e_field) .and. gauge_field_is_applied(this%ep%gfield)) then
        this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) = &
          this%hm_base%uniform_vector_potential(1:mesh%sb%periodic_dim) - time_*this%ep%e_field(1:mesh%sb%periodic_dim)
      end if

    end do

    ! the vector potential of a static magnetic field
    if (allocated(this%ep%a_static)) then
      call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL, .false.)
      do idir = 1, mesh%sb%dim
        !$omp parallel do schedule(static)
        do ip = 1, mesh%np
          this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) + this%ep%a_static(ip, idir)
        end do
      end do
    end if

    !The electric field is added to the KS potential
    call hamiltonian_elec_base_accel_copy_pot(this%hm_base, mesh)

    ! and the static magnetic field
    if (allocated(this%ep%b_field)) then
      call hamiltonian_elec_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_MAGNETIC_FIELD, .false.)
      do idir = 1, 3
        this%hm_base%uniform_magnetic_field(idir) = this%hm_base%uniform_magnetic_field(idir) + this%ep%b_field(idir)
      end do
    end if

    call hamiltonian_elec_base_update(this%hm_base, mesh)

    call build_phase()

    call profiling_out(prof)
    POP_SUB(hamiltonian_elec_update2)

  contains

    subroutine build_phase()
      integer :: ik, imat, nmat, max_npoints, offset, iphase, nphase
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_elec_update2.build_phase)

      if ((.not. this%kpoints%gamma_only()) .or. allocated(this%hm_base%uniform_vector_potential)) then

        call profiling_in(prof_phases, 'UPDATE_PHASES')
        ! now regenerate the phases for the pseudopotentials
        do iatom = 1, this%ep%natoms
          call projector_init_phases(this%ep%proj(iatom), mesh%sb%dim, this%d, this%der%boundaries, this%kpoints,  &
            vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
        end do

        call profiling_out(prof_phases)
      end if

      if(allocated(this%hm_base%uniform_vector_potential)) then
        if (.not. allocated(this%hm_base%phase)) then
          SAFE_ALLOCATE(this%hm_base%phase(1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_phase, ACCEL_MEM_READ_ONLY, TYPE_CMPLX, mesh%np_part*this%d%kpt%nlocal)
          end if
        end if

        kpoint(1:mesh%sb%dim) = M_ZERO
        do ik = this%d%kpt%start, this%d%kpt%end
          kpoint(1:mesh%sb%dim) = this%kpoints%get_point(this%d%get_kpoint_index(ik))

          !$omp parallel do schedule(static)
          do ip = 1, mesh%np_part
            this%hm_base%phase(ip, ik) = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*(kpoint(1:mesh%sb%dim) &
              + this%hm_base%uniform_vector_potential(1:mesh%sb%dim))))
          end do
        end do
        if(accel_is_enabled()) then
          call accel_write_buffer(this%hm_base%buff_phase, mesh%np_part*this%d%kpt%nlocal, this%hm_base%phase)
        end if
      end if

      max_npoints = this%hm_base%max_npoints
      nmat = this%hm_base%nprojector_matrices


      if (allocated(this%hm_base%phase) .and. allocated(this%hm_base%projector_matrices)) then

        nphase = 1
        if(this%der%boundaries%spiralBC) nphase = 3

        if(.not. allocated(this%hm_base%projector_phases)) then
          SAFE_ALLOCATE(this%hm_base%projector_phases(1:max_npoints, nphase, nmat, this%d%kpt%start:this%d%kpt%end))
          if(accel_is_enabled()) then
            call accel_create_buffer(this%hm_base%buff_projector_phases, ACCEL_MEM_READ_ONLY, &
              TYPE_CMPLX, this%hm_base%total_points*nphase*this%d%kpt%nlocal)
          end if
        end if

        offset = 0
        do ik = this%d%kpt%start, this%d%kpt%end
          do imat = 1, this%hm_base%nprojector_matrices
            iatom = this%hm_base%projector_to_atom(imat)
            do iphase = 1, nphase
              !$omp parallel do schedule(static)
              do ip = 1, this%hm_base%projector_matrices(imat)%npoints
                this%hm_base%projector_phases(ip, imat, iphase, ik) = this%ep%proj(iatom)%phase(ip, iphase, ik)
              end do

              if(accel_is_enabled() .and. this%hm_base%projector_matrices(imat)%npoints > 0) then
                call accel_write_buffer(this%hm_base%buff_projector_phases, &
                  this%hm_base%projector_matrices(imat)%npoints, this%hm_base%projector_phases(1:, iphase, imat, ik), & 
                  offset = offset)
              end if
              offset = offset + this%hm_base%projector_matrices(imat)%npoints
            end do
          end do
        end do

      end if

      POP_SUB(hamiltonian_elec_update2.build_phase)
    end subroutine build_phase

  end subroutine hamiltonian_elec_update2

 ! ---------------------------------------------------------
 subroutine hamiltonian_elec_set_vhxc(hm, mesh, vold, vold_tau)
   type(hamiltonian_elec_t), intent(inout)  :: hm
   type(mesh_t),             intent(in)     :: mesh
   FLOAT,                    intent(in)     :: vold(:, :)
   FLOAT, optional,          intent(in)     :: vold_tau(:, :)

   PUSH_SUB(hamiltonian_elec_set_vhxc)

   call lalg_copy(mesh%np, hm%d%nspin, vold, hm%vhxc)
   if(present(vold_tau)) then
     call lalg_copy(mesh%np, hm%d%nspin, vold_tau, hm%vtau)
   end if

   POP_SUB(hamiltonian_elec_set_vhxc)
 end subroutine hamiltonian_elec_set_vhxc

 logical function hamiltonian_elec_needs_current(hm, states_are_real)
    type(hamiltonian_elec_t), intent(in) :: hm
    logical,                  intent(in) :: states_are_real

    hamiltonian_elec_needs_current = .false.

    if( hm%self_induced_magnetic ) then
      if(.not. states_are_real) then
        hamiltonian_elec_needs_current = .true.
      else
        message(1) = 'No current density for real states since it is identically zero.'
        call messages_warning(1)
      end if
    end if

  end function hamiltonian_elec_needs_current

  ! ---------------------------------------------------------
  subroutine zhamiltonian_elec_apply_all(hm, namespace, mesh, st, hst)
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(inout) :: st
    type(states_elec_t),      intent(inout) :: hst

    integer :: ik, ib, ist
    CMPLX, allocatable :: psi(:, :)
    CMPLX, allocatable :: psiall(:, :, :, :)
  
    PUSH_SUB(zhamiltonian_elec_apply_all)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call zhamiltonian_elec_apply_batch(hm, namespace, mesh, st%group%psib(ib, ik), hst%group%psib(ib, ik))
      end do
    end do

    if(oct_exchange_enabled(hm%oct_exchange)) then

      SAFE_ALLOCATE(psiall(mesh%np_part, 1:hst%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

      call states_elec_get_state(st, mesh, psiall)
    
      call oct_exchange_prepare(hm%oct_exchange, mesh, psiall, hm%xc, hm%psolver, namespace)

      SAFE_DEALLOCATE_A(psiall)
    
      SAFE_ALLOCATE(psi(mesh%np_part, 1:hst%d%dim))
    
      do ik = 1, st%d%nik
        do ist = 1, st%nst
          call states_elec_get_state(hst, mesh, ist, ik, psi)
          call oct_exchange_operator(hm%oct_exchange, namespace, mesh, psi, ist, ik)
          call states_elec_set_state(hst, mesh, ist, ik, psi)
        end do
      end do

      SAFE_DEALLOCATE_A(psi)
    
    end if

    POP_SUB(zhamiltonian_elec_apply_all)
  end subroutine zhamiltonian_elec_apply_all


  ! ---------------------------------------------------------

  subroutine magnus(hm, namespace, mesh, psi, hpsi, ik, vmagnus, set_phase)
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(namespace_t),        intent(in)    :: namespace
    type(mesh_t),             intent(in)    :: mesh
    CMPLX,                    intent(inout) :: psi(:,:)
    CMPLX,                    intent(out)   :: hpsi(:,:)
    integer,                  intent(in)    :: ik
    FLOAT,                    intent(in)    :: vmagnus(:, :, :)
    logical, optional,        intent(in)    :: set_phase !< If set to .false. the phase will not be added to the states.

    CMPLX, allocatable :: auxpsi(:, :), aux2psi(:, :)
    integer :: idim, ispin

    PUSH_SUB(magnus)

    ! We will assume, for the moment, no spinors.
    if(hm%d%dim /= 1) &
      call messages_not_implemented("Magnus with spinors", namespace=namespace)

    SAFE_ALLOCATE( auxpsi(1:mesh%np_part, 1:hm%d%dim))
    SAFE_ALLOCATE(aux2psi(1:mesh%np,      1:hm%d%dim))

    ispin = hm%d%get_spin_index(ik)

    ! Compute (T + Vnl)|psi> and store it
    call zhamiltonian_elec_apply_single(hm, namespace, mesh, psi, auxpsi, 1, ik, &
      terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, set_phase = set_phase)

    ! H|psi>  =  (T + Vnl)|psi> + Vpsl|psi> + Vmagnus(t2)|psi> + Vborders
    do idim = 1, hm%d%dim
      call lalg_copy(mesh%np, auxpsi(:, idim), hpsi(:, idim))
      hpsi(1:mesh%np, idim) = hpsi(1:mesh%np, idim) + hm%ep%Vpsl(1:mesh%np)*psi(1:mesh%np,idim)
      call vborders(mesh, hm, psi(:, idim), hpsi(:, idim))
    end do
    hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) + vmagnus(1:mesh%np, ispin, 2)*psi(1:mesh%np, 1)

    ! Add first term of the commutator:  - i Vmagnus(t1) (T + Vnl) |psi>
    hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) - M_zI*vmagnus(1:mesh%np, ispin, 1)*auxpsi(1:mesh%np, 1)

    ! Add second term of commutator:  i (T + Vnl) Vmagnus(t1) |psi>
    auxpsi(1:mesh%np, 1) = vmagnus(1:mesh%np, ispin, 1)*psi(1:mesh%np, 1)
    call zhamiltonian_elec_apply_single(hm, namespace, mesh, auxpsi, aux2psi, 1, ik, &
      terms = TERM_KINETIC + TERM_NON_LOCAL_POTENTIAL, set_phase = set_phase)
    hpsi(1:mesh%np, 1) = hpsi(1:mesh%np, 1) + M_zI*aux2psi(1:mesh%np, 1)

    SAFE_DEALLOCATE_A(auxpsi)
    SAFE_DEALLOCATE_A(aux2psi)
    POP_SUB(magnus)
  end subroutine magnus
  
  ! ---------------------------------------------------------
  subroutine vborders (mesh, hm, psi, hpsi)
    type(mesh_t),             intent(in)    :: mesh
    type(hamiltonian_elec_t), intent(in)    :: hm
    CMPLX,                    intent(in)    :: psi(:)
    CMPLX,                    intent(inout) :: hpsi(:)

    integer :: ip

    PUSH_SUB(vborders)

    if(hm%bc%abtype == IMAGINARY_ABSORBING) then
      do ip = 1, mesh%np
        hpsi(ip) = hpsi(ip) + M_zI*hm%bc%mf(ip)*psi(ip)
      end do
    end if

    POP_SUB(vborders)
  end subroutine vborders

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_elec_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_elec_inc.F90"

end module hamiltonian_elec_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
