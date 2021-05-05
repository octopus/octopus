!! Copyright (C) 2016-2020 N. Tancogne-Dejean
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

module lda_u_oct_m
  use atomic_orbital_oct_m
  use boundaries_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use energy_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_base_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use loewdin_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use orbitalbasis_oct_m
  use orbitalset_oct_m
  use orbitalset_utils_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use unit_system_oct_m
  use wfs_elec_oct_m

  implicit none

  private

  public ::                          &
    lda_u_t,                         &
    lda_u_init,                      &
    dlda_u_apply,                    &
    zlda_u_apply,                    &
    lda_u_update_basis,              &
    lda_u_update_occ_matrices,       &
    lda_u_end,                       &
    lda_u_build_phase_correction,    &
    lda_u_freeze_occ,                &
    lda_u_freeze_u,                  &
    lda_u_periodic_coulomb_integrals,&
    dlda_u_set_occupations,          &
    zlda_u_set_occupations,          &
    dlda_u_get_occupations,          &
    zlda_u_get_occupations,          &
    dlda_u_update_potential,         &
    zlda_u_update_potential,         &
    lda_u_get_effectiveU,            &
    lda_u_set_effectiveU,            &
    lda_u_get_effectiveV,            &
    lda_u_set_effectiveV,            &
    dlda_u_commute_r,                &
    zlda_u_commute_r,                &
    dlda_u_force,                    &
    zlda_u_force,                    &
    lda_u_write_info,                &
    compute_ACBNO_U_kanamori,        &
    dcompute_dftu_energy,            &
    zcompute_dftu_energy

  
  integer, public, parameter ::        &
    DFT_U_NONE                    = 0, &
    DFT_U_EMPIRICAL               = 1, &
    DFT_U_ACBN0                   = 2

  integer, public, parameter ::        &
    DFT_U_FLL                     = 0, &
    DFT_U_AMF                     = 1

  
  type lda_u_t
    private
    integer,            public   :: level = DFT_U_NONE
    FLOAT, allocatable, public   :: dn(:,:,:,:) !> Occupation matrices for the standard scheme
    FLOAT, allocatable           :: dV(:,:,:,:) !> Potentials for the standard scheme

    CMPLX, allocatable, public   :: zn(:,:,:,:)
    CMPLX, allocatable           :: zV(:,:,:,:)
    FLOAT, allocatable, public   :: dn_alt(:,:,:,:) !> Stores the renomalized occ. matrices
    CMPLX, allocatable, public   :: zn_alt(:,:,:,:) !> if the ACBN0 functional is used

    FLOAT, allocatable           :: renorm_occ(:,:,:,:,:) !> On-site occupations (for the ACBN0 functional)

    FLOAT, allocatable           :: coulomb(:,:,:,:,:) !>Coulomb integrals for all the system
                                                       !> (for the ACBN0 functional)
    CMPLX, allocatable           :: zcoulomb(:,:,:,:,:,:,:) !>Coulomb integrals for all the system
                                                            !> (for the ACBN0 functional with spinors)

    type(orbitalbasis_t),        public :: basis                !> The full basis of localized orbitals
    type(orbitalset_t), pointer, public :: orbsets(:) => NULL() !> All the orbital setss of the system
    integer,                     public :: norbsets = 0

    integer,              public :: nspins = 0
    integer,              public :: spin_channels = 0
    integer                      :: nspecies = 0
    integer,              public :: maxnorbs = 0           !> Maximal number of orbitals for all the atoms
    integer                      :: max_np = 0             !> Maximum number of points in all orbitals submesh spheres

    logical                      :: useAllOrbitals = .false.       !> Do we use all atomic orbitals possible
    logical                      :: skipSOrbitals = .true.         !> Not using s orbitals
    logical                      :: freeze_occ = .false.           !> Occupation matrices are not recomputed during TD evolution
    logical                      :: freeze_u = .false.             !> U is not recomputed during TD evolution
    logical,              public :: intersite = .false.            !> intersite V are computed or not
    FLOAT                        :: intersite_radius = M_ZERO      !> Maximal distance for considering neighboring atoms
    logical,              public :: basisfromstates = .false.      !> We can construct the localized basis from user-defined states
    FLOAT                        :: acbn0_screening = M_ONE        !> We use or not the screening in the ACBN0 functional
    integer, allocatable, public :: basisstates(:)
    logical                      :: rot_inv = .false.              !> Use a rotationally invariant formula for U and J (ACBN0 case)
    integer                      :: double_couting = DFT_U_FLL     !> Double-couting term
    integer                      :: sm_poisson = SM_POISSON_DIRECT !> Poisson solver used for computing Coulomb integrals

    type(distributed_t) :: orbs_dist

    integer, public     :: maxneighbors = 0
    FLOAT, allocatable  :: dn_ij(:,:,:,:,:), dn_alt_ij(:,:,:,:,:), dn_alt_ii(:,:,:,:,:)
    CMPLX, allocatable  :: zn_ij(:,:,:,:,:), zn_alt_ij(:,:,:,:,:), zn_alt_ii(:,:,:,:,:)
  end type lda_u_t

contains

  ! ---------------------------------------------------------
  subroutine lda_u_init(this, namespace, space, level, gr, ions, st, psolver, kpoints)
    type(lda_u_t),     target, intent(inout) :: this
    type(namespace_t),         intent(in)    :: namespace
    type(space_t),             intent(in)    :: space
    integer,                   intent(in)    :: level
    type(grid_t),              intent(in)    :: gr
    type(ions_t),      target, intent(in)    :: ions
    type(states_elec_t),       intent(in)    :: st
    type(poisson_t),           intent(in)    :: psolver
    type(kpoints_t),           intent(in)    :: kpoints

    logical :: complex_coulomb_integrals
    integer :: ios, is
    type(block_t) :: blk

    PUSH_SUB(lda_u_init)

    ASSERT(.not. (level == DFT_U_NONE))

    call messages_print_stress(stdout, "DFT+U", namespace=namespace)
    if(gr%mesh%parallel_in_domains) call messages_experimental("dft+u parallel in domains")
    this%level = level

    !%Variable DFTUBasisFromStates
    !%Type logical
    !%Default no
    !%Section Hamiltonian::DFT+U
    !%Description
    !% If set to yes, Octopus will construct the localized basis from
    !% user-defined states. The states are taken at the Gamma point (or the first k-point of the
    !% states in the restart_proj folder.
    !% The states are defined via the block DFTUBasisStates
    !%End
    call parse_variable(namespace, 'DFTUBasisFromStates', .false., this%basisfromstates)
    if(this%basisfromstates) call messages_experimental("DFTUBasisFromStates")

    !%Variable DFTUDoubleCounting
    !%Type integer
    !%Default dft_u_fll
    !%Section Hamiltonian::DFT+U
    !%Description
    !% This variable selects which DFT+U
    !% double counting term is used.
    !%Option dft_u_fll 0
    !% (Default) The Fully Localized Limit (FLL)
    !%Option dft_u_amf 1
    !% (Experimental) Around mean field double counting, as defined in PRB 44, 943 (1991) and PRB 49, 14211 (1994).
    !%End
    call parse_variable(namespace, 'DFTUDoubleCounting', DFT_U_FLL, this%double_couting)
    call messages_print_var_option(stdout,  'DFTUDoubleCounting', this%double_couting)
    if(this%double_couting /= DFT_U_FLL) call messages_experimental("DFTUDoubleCounting = dft_u_amf")
    if(st%d%ispin == SPINORS .and. this%double_couting /= DFT_U_FLL) then
      call messages_not_implemented("AMF double couting with spinors.", namespace=namespace)
    end if

    !%Variable DFTUPoissonSolver
    !%Type integer
    !%Section Hamiltonian::DFT+U
    !%Description
    !% This variable selects which Poisson solver
    !% is used to compute the Coulomb integrals over a submesh.
    !% These are non-periodic Poisson solvers.
    !% If the domain parallelization is activated, the default is the direct sum.
    !% Otherwise, the FFT Poisson solver is used by default.
    !%
    !%Option dft_u_poisson_direct 0
    !% (Default) Direct Poisson solver. Slow.
    !%Option dft_u_poisson_isf 1
    !% (Experimental) ISF Poisson solver on a submesh.
    !% This does not work for non-orthogonal cells nor domain parallelization.
    !%Option dft_u_poisson_psolver 2
    !% (Experimental) PSolver Poisson solver on a submesh.
    !% This does not work for non-orthogonal cells nor domain parallelization.
    !% Requires the PSolver external library.
    !%Option dft_u_poisson_fft 3
    !% FFT Poisson solver on a submesh.
    !% This uses the 0D periodic version of the FFT kernels.
    !% This is not implemented for domain parallelization.
    !%End
    if(gr%mesh%parallel_in_domains) then
      call parse_variable(namespace, 'DFTUPoissonSolver', SM_POISSON_DIRECT, this%sm_poisson)
    else
      call parse_variable(namespace, 'DFTUPoissonSolver', SM_POISSON_FFT, this%sm_poisson)
    end if
    call messages_print_var_option(stdout,  'DFTUPoissonSolver', this%sm_poisson)
    if(this%sm_poisson /= SM_POISSON_DIRECT .and. this%sm_poisson /= SM_POISSON_FFT) then
      call messages_experimental("DFTUPoissonSolver different from dft_u_poisson_direct and dft_u_poisson_fft")
    end if
    if(this%sm_poisson == SM_POISSON_ISF) then
      if(gr%mesh%parallel_in_domains) then
        call messages_not_implemented("DFTUPoissonSolver=dft_u_poisson_isf with domain parallelization")
      end if
      if (ions%latt%nonorthogonal) then
        call messages_not_implemented("DFTUPoissonSolver=dft_u_poisson_isf with non-orthogonal cells")
      end if
    end if
    if(this%sm_poisson == SM_POISSON_PSOLVER) then
#if !((defined HAVE_LIBISF) || (defined HAVE_PSOLVER))
      message(1) = "The PSolver Poisson solver cannot be used since the code was not compiled with the PSolver libary."
      call messages_fatal(1)
#endif
      if(gr%mesh%parallel_in_domains) then
        call messages_not_implemented("DFTUPoissonSolver=dft_u_poisson_psolver with domain parallelization")
      end if
      if (ions%latt%nonorthogonal) then
        call messages_not_implemented("DFTUPoissonSolver=dft_u_poisson_psolver with non-orthogonal cells")
      end if
    end if
    if(this%sm_poisson == SM_POISSON_FFT) then
      if(gr%mesh%parallel_in_domains) then
        call messages_not_implemented("DFTUPoissonSolver=dft_u_poisson_fft with domain parallelization.")
      end if
    end if

    if(this%level == DFT_U_ACBN0 ) then
      !%Variable UseAllAtomicOrbitals
      !%Type logical
      !%Default no
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to yes, Octopus will determine the effective U for all atomic orbitals
      !% from the peusopotential. Only available with ACBN0 functional.
      !% It is strongly recommended to set AOLoewdin=yes when using the option.
      !%End
      call parse_variable(namespace, 'UseAllAtomicOrbitals', .false., this%useAllOrbitals)
      if(this%useAllOrbitals) call messages_experimental("UseAllAtomicOrbitals")

      !%Variable SkipSOrbitals
      !%Type logical
      !%Default no
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to yes, Octopus will determine the effective U for all atomic orbitals
      !% from the peusopotential but s orbitals. Only available with ACBN0 functional.
      !%End
      call parse_variable(namespace, 'SkipSOrbitals', .true., this%skipSOrbitals)
      if(.not.this%SkipSOrbitals) call messages_experimental("SkipSOrbitals")

      !%Variable ACBN0Screening
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to 0, no screening will be included in the ACBN0 functional, and the U
      !% will be estimated from bare Hartree-Fock. If set to 1 (default), the full screening
      !% of the U, as defined in the ACBN0 functional, is used.
      !%End
      call parse_variable(namespace, 'ACBN0Screening', M_ONE, this%acbn0_screening)
      call messages_print_var_value(stdout, 'ACBN0Screening', this%acbn0_screening)

      !%Variable ACBN0RotationallyInvariant
      !%Type logical
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to yes, Octopus will use for U and J a formula which is rotationally invariant.
      !% This is different from the original formula for U and J.
      !% This is activated by default, except in the case of spinors, as this is not yet implemented in this case.
      !%End
      call parse_variable(namespace, 'ACBN0RotationallyInvariant', st%d%ispin /= SPINORS, this%rot_inv)
      call messages_print_var_value(stdout, 'ACBN0RotationallyInvariant', this%rot_inv)
      if(this%rot_inv .and. st%d%ispin == SPINORS ) then
        call messages_not_implemented("Rotationally invariant ACBN0 with spinors.", namespace=namespace)
      end if

      !%Variable ACBN0IntersiteInteraction
      !%Type logical
      !%Default no
      !%Section Hamiltonian::DFT+U
      !%Description
      !% If set to yes, Octopus will determine the effective intersite interaction V
      !% Only available with ACBN0 functional.
      !% It is strongly recommended to set AOLoewdin=yes when using the option.
      !%End
      call parse_variable(namespace, 'ACBN0IntersiteInteraction', .false., this%intersite)
      call messages_print_var_value(stdout, 'ACBN0IntersiteInteraction', this%intersite)
      if(this%intersite) call messages_experimental("ACBN0IntersiteInteraction")

      if(this%intersite) then

        !This is a non local operator. To make this working, one probably needs to apply the
        ! symmetries to the generalized occupation matrices
        if(kpoints%use_symmetries) then
          call messages_not_implemented("Intersite interaction with kpoint symmetries", namespace=namespace)
        end if

        !%Variable ACBN0IntersiteCutoff
        !%Type float
        !%Section Hamiltonian::DFT+U
        !%Description
        !% The cutoff radius defining the maximal intersite distance considered.
        !% Only available with ACBN0 functional with intersite interaction.
        !%End
        call parse_variable(namespace, 'ACBN0IntersiteCutoff', M_ZERO, this%intersite_radius, unit = units_inp%length)
        if(abs(this%intersite_radius) < M_EPSILON) then
          call messages_write("ACBN0IntersiteCutoff must be greater than 0")
          call messages_fatal(1, namespace=namespace)
        end if

      end if

    end if
    
    call lda_u_write_info(this, stdout)

    if(.not.this%basisfromstates) then

      call orbitalbasis_init(this%basis, namespace, space%periodic_dim)

      if (states_are_real(st)) then
        call dorbitalbasis_build(this%basis, ions, gr%mesh, st%d%kpt, st%d%dim, &
          this%skipSOrbitals, this%useAllOrbitals)
      else
        call zorbitalbasis_build(this%basis, ions, gr%mesh, st%d%kpt, st%d%dim, &
          this%skipSOrbitals, this%useAllOrbitals)
      end if
      this%orbsets => this%basis%orbsets
      this%norbsets = this%basis%norbsets
      this%maxnorbs = this%basis%maxnorbs
      this%max_np = this%basis%max_np
      this%nspins = st%d%nspin
      this%spin_channels = st%d%spin_channels
      this%nspecies = ions%nspecies

      !We allocate the necessary ressources
      if (states_are_real(st)) then
        call dlda_u_allocate(this, st)
      else
        call zlda_u_allocate(this, st)
      end if

      call distributed_nullify(this%orbs_dist, this%norbsets)
#ifdef HAVE_MPI
      if(.not. gr%mesh%parallel_in_domains) then
        call distributed_init(this%orbs_dist, this%norbsets, MPI_COMM_WORLD, "orbsets")
      end if
#endif


      if( this%level == DFT_U_ACBN0 ) then

        complex_coulomb_integrals = .false.
        do ios = 1, this%norbsets
          if(this%orbsets(ios)%ndim  > 1) complex_coulomb_integrals = .true.
        end do

        if(.not. complex_coulomb_integrals) then
          write(message(1),'(a)')    'Computing the Coulomb integrals of the localized basis.'
          call messages_info(1)
          if (states_are_real(st)) then
            call dcompute_coulomb_integrals(this, namespace, space, gr%mesh, gr%der, psolver)
          else
            call zcompute_coulomb_integrals(this, namespace, space, gr%mesh, gr%der, psolver)
          end if
        else
          ASSERT(.not.states_are_real(st))
          write(message(1),'(a)')    'Computing complex Coulomb integrals of the localized basis.'
          call messages_info(1)
          call compute_complex_coulomb_integrals(this, gr%mesh, gr%der, st, psolver, namespace, space)
        end if 
      end if

    else

      !%Variable DFTUBasisStates
      !%Type block
      !%Default none
      !%Section Hamiltonian::DFT+U
      !%Description
      !% Each line of this block contains the index of a state to be used to construct the
      !% localized basis. See DFTUBasisFromStates for details.
      !%End
      if(parse_block(namespace, 'DFTUBasisStates', blk) == 0) then
        this%norbsets = 1
        this%maxnorbs = parse_block_n(blk)
        if(this%maxnorbs <1) then
          write(message(1),'(a,i3,a,i3)') 'DFTUBasisStates must contains at least one state.'
          call messages_fatal(1, namespace=namespace)
        end if
        SAFE_ALLOCATE(this%basisstates(1:this%maxnorbs))
        do is = 1, this%maxnorbs
          call parse_block_integer(blk, is-1, 0, this%basisstates(is))
        end do
        call parse_block_end(blk)
      else
        write(message(1),'(a,i3,a,i3)') 'DFTUBasisStates must be specified if DFTUBasisFromStates=yes'
        call messages_fatal(1, namespace=namespace)
      end if

      if (states_are_real(st)) then
        call dorbitalbasis_build_empty(this%basis, gr%mesh, st%d%kpt, st%d%dim, this%maxnorbs)
      else
        call zorbitalbasis_build_empty(this%basis, gr%mesh, st%d%kpt, st%d%dim, this%maxnorbs)
      end if

      this%max_np = gr%mesh%np
      this%nspins = st%d%nspin
      this%spin_channels = st%d%spin_channels
      this%nspecies = 1

      this%orbsets => this%basis%orbsets

      call distributed_nullify(this%orbs_dist, this%norbsets)


      !We allocate the necessary ressources
      if (states_are_real(st)) then
        call dlda_u_allocate(this, st)
      else
        call zlda_u_allocate(this, st)
      end if

    end if

    call messages_print_stress(stdout, namespace=namespace)

    POP_SUB(lda_u_init)
  end subroutine lda_u_init

  ! ---------------------------------------------------------
  subroutine lda_u_end(this)
    implicit none
    type(lda_u_t), intent(inout) :: this

    PUSH_SUB(lda_u_end)

    this%level = DFT_U_NONE

    SAFE_DEALLOCATE_A(this%dn)
    SAFE_DEALLOCATE_A(this%zn)
    SAFE_DEALLOCATE_A(this%dn_alt)
    SAFE_DEALLOCATE_A(this%zn_alt)
    SAFE_DEALLOCATE_A(this%dV)
    SAFE_DEALLOCATE_A(this%zV)
    SAFE_DEALLOCATE_A(this%coulomb)
    SAFE_DEALLOCATE_A(this%zcoulomb)
    SAFE_DEALLOCATE_A(this%renorm_occ)
    SAFE_DEALLOCATE_A(this%dn_ij)
    SAFE_DEALLOCATE_A(this%zn_ij)
    SAFE_DEALLOCATE_A(this%dn_alt_ij)
    SAFE_DEALLOCATE_A(this%zn_alt_ij)
    SAFE_DEALLOCATE_A(this%dn_alt_ii)
    SAFE_DEALLOCATE_A(this%zn_alt_ii)
    SAFE_DEALLOCATE_A(this%basisstates)

    nullify(this%orbsets)
    call orbitalbasis_end(this%basis)

    this%max_np = 0

    if(.not.this%basisfromstates) then
#ifdef HAVE_MPI
      call distributed_end(this%orbs_dist)
#endif
    end if

    POP_SUB(lda_u_end)
  end subroutine lda_u_end

  ! When moving the ions, the basis must be reconstructed
  subroutine lda_u_update_basis(this, space, gr, ions, st, psolver, namespace, kpoints, has_phase)
    type(lda_u_t),     target, intent(inout) :: this
    type(space_t),             intent(in)    :: space
    type(grid_t),              intent(in)    :: gr
    type(ions_t),      target, intent(in)    :: ions
    type(states_elec_t),       intent(in)    :: st
    type(poisson_t),           intent(in)    :: psolver
    type(namespace_t),         intent(in)    :: namespace
    type(kpoints_t),           intent(in)    :: kpoints
    logical,                   intent(in)    :: has_phase

    integer :: ios, maxorbs, nspin

    if(this%level == DFT_U_NONE) return
    !If we use a basis from states, there is nothing to do
    if(this%basisfromstates) return

    PUSH_SUB(lda_u_update_basis)

    !We clean the orbital basis, to be able to reconstruct it
    call orbitalbasis_end(this%basis)
    nullify(this%orbsets)

    !We now reconstruct the basis
    if (states_are_real(st)) then
      call dorbitalbasis_build(this%basis, ions, gr%mesh, st%d%kpt, st%d%dim, &
        this%skipSOrbitals, this%useAllOrbitals, verbose = .false.)
    else
      call zorbitalbasis_build(this%basis, ions, gr%mesh, st%d%kpt, st%d%dim, &
        this%skipSOrbitals, this%useAllOrbitals, verbose = .false.)
    end if
    this%orbsets => this%basis%orbsets

    !In case of intersite interaction we need to reconstruct the basis
    if(this%intersite) then
      this%maxneighbors = 0
      do ios = 1, this%norbsets
        call orbitalset_init_intersite(this%orbsets(ios), namespace, ios, ions, gr%der, psolver, &
          this%orbsets, this%norbsets, this%maxnorbs, this%intersite_radius, st%d%kpt, has_phase, this%sm_poisson)
        this%maxneighbors = max(this%maxneighbors, this%orbsets(ios)%nneighbors)
      end do

      maxorbs = this%maxnorbs
      nspin = this%nspins

      if(states_are_real(st)) then
        SAFE_DEALLOCATE_A(this%dn_ij)
        SAFE_ALLOCATE(this%dn_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%dn_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_ZERO
        SAFE_DEALLOCATE_A(this%dn_alt_ij)
        SAFE_ALLOCATE(this%dn_alt_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%dn_alt_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_ZERO
        SAFE_DEALLOCATE_A(this%dn_alt_ii)
        SAFE_ALLOCATE(this%dn_alt_ii(1:2,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%dn_alt_ii(1:2,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_ZERO
      else
        SAFE_DEALLOCATE_A(this%zn_ij)
        SAFE_ALLOCATE(this%zn_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%zn_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_Z0
        SAFE_DEALLOCATE_A(this%zn_alt_ij)
        SAFE_ALLOCATE(this%zn_alt_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%zn_alt_ij(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_Z0
        SAFE_DEALLOCATE_A(this%zn_alt_ii)
        SAFE_ALLOCATE(this%zn_alt_ii(1:2,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors))
        this%zn_alt_ii(1:2,1:maxorbs,1:nspin,1:this%norbsets,1:this%maxneighbors) = M_Z0
      end if
    end if

    ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
    ! In case of a laser field, the phase is recomputed in hamiltonian_elec_update
    if(has_phase) then
      call lda_u_build_phase_correction(this, space, st%d, gr%der%boundaries, namespace, kpoints)
    else
      !In case there is no phase, we perform the orthogonalization here
      if(this%basis%orthogonalization) then
        call dloewdin_orthogonalize(this%basis, st%d%kpt, namespace)
      else
        if(debug%info .and. space%is_periodic()) then
          call dloewdin_info(this%basis, st%d%kpt, namespace)
        end if
      end if
    end if

    POP_SUB(lda_u_update_basis)

  end subroutine lda_u_update_basis

  ! Interface for the X(update_occ_matrices) routines
  subroutine lda_u_update_occ_matrices(this, namespace, mesh, st, hm_base, energy )
    type(lda_u_t),                 intent(inout) :: this
    type(namespace_t),             intent(in)    :: namespace
    type(mesh_t),                  intent(in)    :: mesh
    type(states_elec_t),           intent(in)    :: st
    type(hamiltonian_elec_base_t), intent(in)    :: hm_base
    type(energy_t),                intent(inout) :: energy

    if(this%level == DFT_U_NONE .or. this%freeze_occ) return
    PUSH_SUB(lda_u_update_occ_matrices)

    if (states_are_real(st)) then
      call dupdate_occ_matrices(this, namespace, mesh, st, energy%dft_u)
    else
      if (allocated(hm_base%phase)) then
        call zupdate_occ_matrices(this, namespace, mesh, st, energy%dft_u, hm_base%phase)
      else
        call zupdate_occ_matrices(this, namespace, mesh, st, energy%dft_u)
      end if
    end if

    POP_SUB(lda_u_update_occ_matrices)
  end subroutine lda_u_update_occ_matrices


  !> Build the phase correction to the global phase for all orbitals
  subroutine lda_u_build_phase_correction(this, space, std, boundaries, namespace, kpoints, vec_pot, vec_pot_var)
    type(lda_u_t),                 intent(inout) :: this
    type(space_t),                 intent(in)    :: space
    type(states_elec_dim_t),       intent(in)    :: std
    type(boundaries_t),            intent(in)    :: boundaries
    type(namespace_t),             intent(in)    :: namespace
    type(kpoints_t),               intent(in)    :: kpoints
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot(:) !< (space%dim)
    FLOAT, optional,  allocatable, intent(in)    :: vec_pot_var(:, :) !< (1:space%dim, 1:ns)

    integer :: ios

    if(boundaries%spiralBC) call messages_not_implemented("DFT+U with spiral boundary conditions.", &
      namespace=namespace)

    !In this case there is no phase difference, as the basis come from states on the full
    !grid and not from spherical meshes around the atoms
    if(this%basisfromstates) return

    PUSH_SUB(lda_u_build_phase_correction)

    do ios = 1, this%norbsets
      call orbitalset_update_phase(this%orbsets(ios), space%dim, std%kpt, kpoints, (std%ispin==SPIN_POLARIZED), &
        vec_pot, vec_pot_var)
    end do

    if(this%basis%orthogonalization) then
      call zloewdin_orthogonalize(this%basis, std%kpt, namespace)
    else
      if(debug%info .and. space%is_periodic()) call zloewdin_info(this%basis, std%kpt, namespace)
    end if

    POP_SUB(lda_u_build_phase_correction)

  end subroutine lda_u_build_phase_correction

  ! ---------------------------------------------------------
  subroutine lda_u_periodic_coulomb_integrals(this, namespace, space, st, der, mc, has_phase)
    type(lda_u_t),                 intent(inout) :: this
    type(namespace_t),             intent(in)    :: namespace
    type(space_t)    ,             intent(in)    :: space
    type(states_elec_t),           intent(in)    :: st
    type(derivatives_t),           intent(in)    :: der
    type(multicomm_t),             intent(in)    :: mc
    logical,                       intent(in)    :: has_phase

    integer :: ik, im, idim

    if(this%level /= DFT_U_ACBN0) return

    ASSERT(this%basisfromstates)

    PUSH_SUB(lda_u_periodic_coulomb_integrals)

    if(states_are_real(st)) then
      call dcompute_periodic_coulomb_integrals(this, namespace, space, der, mc)
    else
      call zcompute_periodic_coulomb_integrals(this, namespace, space, der, mc)
    end if

    ! We rebuild the phase for the orbital projection, similarly to the one of the pseudopotentials
    ! In case of a laser field, the phase is recomputed in hamiltonian_elec_update
    if(has_phase) then
      ASSERT(states_are_complex(st))
      do ik = st%d%kpt%start, st%d%kpt%end
        do im = 1, this%orbsets(1)%norbs
          do idim = 1, st%d%dim
            call lalg_copy(der%mesh%np, this%orbsets(1)%zorb(:,idim, im), &
              this%orbsets(1)%eorb_mesh(:,im,idim,ik))
          end do
        end do
      end do
    end if


    POP_SUB(lda_u_periodic_coulomb_integrals)
  end subroutine lda_u_periodic_coulomb_integrals

  ! ---------------------------------------------------------
  subroutine compute_ACBNO_U_kanamori(this, st, kanamori)
    type(lda_u_t),        intent(in)  :: this
    type(states_elec_t),  intent(in)  :: st
    FLOAT,                intent(out) :: kanamori(:,:)

    if(this%nspins == 1) then
      if(states_are_real(st)) then
        call dcompute_ACBNO_U_kanamori_restricted(this, kanamori)
      else
        call zcompute_ACBNO_U_kanamori_restricted(this, kanamori)
      end if
    else
      if(states_are_real(st)) then
        call dcompute_ACBNO_U_kanamori(this, kanamori)
      else
        call zcompute_ACBNO_U_kanamori(this, kanamori)
      end if
    end if


  end subroutine compute_ACBNO_U_kanamori

  ! ---------------------------------------------------------
  subroutine lda_u_freeze_occ(this)
    type(lda_u_t),     intent(inout) :: this

    this%freeze_occ = .true.
  end subroutine lda_u_freeze_occ

  ! ---------------------------------------------------------
  subroutine lda_u_freeze_u(this)
    type(lda_u_t),     intent(inout) :: this

    this%freeze_u = .true.
  end subroutine lda_u_freeze_u

  ! ---------------------------------------------------------
  subroutine lda_u_set_effectiveU(this, Ueff)
    type(lda_u_t),  intent(inout) :: this
    FLOAT,          intent(in)    :: Ueff(:) !< (this%norbsets)

    integer :: ios

    PUSH_SUB(lda_u_set_effectiveU)

    do ios = 1,this%norbsets
      this%orbsets(ios)%Ueff = Ueff(ios)
    end do

    POP_SUB(lda_u_set_effectiveU)
  end subroutine lda_u_set_effectiveU

  ! ---------------------------------------------------------
  subroutine lda_u_get_effectiveU(this, Ueff)
    type(lda_u_t),  intent(in)    :: this
    FLOAT,          intent(inout) :: Ueff(:) !< (this%norbsets)

    integer :: ios

    PUSH_SUB(lda_u_get_effectiveU)

    do ios = 1,this%norbsets
      Ueff(ios) = this%orbsets(ios)%Ueff
    end do

    POP_SUB(lda_u_get_effectiveU)
  end subroutine lda_u_get_effectiveU

  ! ---------------------------------------------------------
  subroutine lda_u_set_effectiveV(this, Veff)
    type(lda_u_t),  intent(inout) :: this
    FLOAT,          intent(in)    :: Veff(:)

    integer :: ios, ncount

    PUSH_SUB(lda_u_set_effectiveV)

    ncount = 0
    do ios = 1, this%norbsets
      this%orbsets(ios)%V_ij(1:this%orbsets(ios)%nneighbors,0) = Veff(ncount+1:ncount+this%orbsets(ios)%nneighbors)
      ncount = ncount + this%orbsets(ios)%nneighbors
    end do

    POP_SUB(lda_u_set_effectiveV)
  end subroutine lda_u_set_effectiveV

  ! ---------------------------------------------------------
  subroutine lda_u_get_effectiveV(this, Veff)
    type(lda_u_t),  intent(in)    :: this
    FLOAT,          intent(inout) :: Veff(:)

    integer :: ios, ncount

    PUSH_SUB(lda_u_get_effectiveV)

    ncount = 0
    do ios = 1, this%norbsets
      Veff(ncount+1:ncount+this%orbsets(ios)%nneighbors) = this%orbsets(ios)%V_ij(1:this%orbsets(ios)%nneighbors,0)
      ncount = ncount + this%orbsets(ios)%nneighbors
    end do

    POP_SUB(lda_u_get_effectiveV)
  end subroutine lda_u_get_effectiveV

  ! ---------------------------------------------------------
  subroutine lda_u_write_info(this, iunit)
    type(lda_u_t),  intent(in)    :: this
    integer,        intent(in)    :: iunit

    PUSH_SUB(lda_u_write_info)

    write(message(1), '(1x)')
    call messages_info(1, iunit)
    if(this%level == DFT_U_EMPIRICAL) then
      write(message(1), '(a)') "Method:"
      write(message(2), '(a)') "  [1] Dudarev et al., Phys. Rev. B 57, 1505 (1998)"
      call messages_info(2, iunit)
    else
      if(.not. this%intersite) then
        write(message(1), '(a)') "Method:"
        write(message(2), '(a)') "  [1] Agapito et al., Phys. Rev. X 5, 011006 (2015)"
        call messages_info(2, iunit)
      else
        write(message(1), '(a)') "Method:"
        write(message(2), '(a)') "  [1] Tancogne-Dejean, and Rubio, Phys. Rev. B 102, 155117 (2020)"
        call messages_info(2, iunit)    
      end if
    end if
    write(message(1), '(a)') "Implementation:"
    write(message(2), '(a)') "  [1] Tancogne-Dejean, Oliveira, and Rubio, Phys. Rev. B 69, 245133 (2017)"
    write(message(3), '(1x)')
    call messages_info(3, iunit)

    POP_SUB(lda_u_write_info)

  end subroutine lda_u_write_info

#include "dft_u_noncollinear_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "lda_u_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "lda_u_inc.F90"
end module lda_u_oct_m
