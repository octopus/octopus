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

module hamiltonian_m
  use batch_m
  use blas_m
  use datasets_m
  use derivatives_m
  use hamiltonian_base_m
  use epot_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use gridhier_m
  use hardware_m
  use io_m
  use io_function_m
  use kpoints_m
  use lalg_basic_m
  use lasers_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use ob_interface_m
  use ob_lead_m
  use opencl_m
  use parser_m
  use poisson_m
  use profiling_m
  use projector_m
  use simul_box_m
  use smear_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use xc_m
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
    hamiltonian_oct_exchange,        &
    hamiltonian_set_oct_exchange,    &
    hamiltonian_remove_oct_exchange, &
    hamiltonian_adjoint,             &
    hamiltonian_not_adjoint,         &
    hamiltonian_hermitian,           &
    hamiltonian_epot_generate,       &
    hamiltonian_update,              &
    hamiltonian_get_time,            &
    hamiltonian_apply_packed

  public ::                          &
    energy_t,                        &
    energy_copy

  type energy_t
    ! Energies
    FLOAT :: total       !< Total energy E = Eii + Sum[Eigenvalues] - U + Ex + Ec - Int[n v_xc]
    FLOAT :: eigenvalues !< Sum[Eigenvalues]
    FLOAT :: exchange
    FLOAT :: correlation
    FLOAT :: xc_j
    FLOAT :: intnvxc     !< Int[n vxc]
    FLOAT :: hartree     !< Hartree      U = (1/2)*Int [n v_Hartree]
    FLOAT :: kinetic     !< Kinetic energy of the non-interacting (KS) system of electrons
    FLOAT :: extern      !< External     V = <Phi|V|Phi> = Int[n v] (if no non-local pseudos exist)
    FLOAT :: entropy
    FLOAT :: ts          !< TS
    FLOAT :: berry       !< Berry energy correction = -mu.E - <Vberry>
  end type energy_t

  type hamiltonian_t
    ! The Hamiltonian must know what are the "dimensions" of the spaces,
    ! in order to be able to operate on the states.
    type(states_dim_t)       :: d
    type(hamiltonian_base_t) :: hm_base
    type(energy_t), pointer  :: energy
    FLOAT, pointer :: vhartree(:) ! Hartree potential
    FLOAT, pointer :: vxc(:,:)    ! XC potential
    FLOAT, pointer :: vhxc(:,:)   ! XC potential + Hartree potential + Berry potential
    FLOAT, pointer :: axc(:,:,:)  ! XC vector potential divided by c
    FLOAT, pointer :: vtau(:,:)   ! Derivative of e_XC w.r.t. tau
    FLOAT, pointer :: vberry(:,:) ! Berry phase potential from external E_field

    FLOAT :: exx_coef ! how much of EXX to mix

    ! The self-induced vector potential and magnetic field
    logical :: self_induced_magnetic
    FLOAT, pointer :: a_ind(:, :)
    FLOAT, pointer :: b_ind(:, :)

    integer :: theory_level    ! copied from sys%ks
    integer :: xc_family       ! copied from sys%ks

    type(epot_t) :: ep         ! handles the external potential

    ! absorbing boundaries
    logical :: adjoint
    integer  :: ab                ! do we have absorbing boundaries?
    FLOAT :: ab_width             ! width of the absorbing boundary
    FLOAT :: ab_height            ! height of the absorbing boundary
    FLOAT, pointer :: ab_pot(:)   ! where we store the ab potential

    ! Open boundaries.
    type(lead_t) :: lead(2*MAX_DIM)

    ! Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span

    ! Mass of the particle (in most cases, mass = 1, electron mass)
    FLOAT :: mass
    ! anisotropic scaling factor for the mass: different along x,y,z etc...
    FLOAT :: mass_scaling(MAX_DIM)

    ! For the Hartree-Fock Hamiltonian, the Fock operator depends on the states.
    type(states_t), pointer :: hf_st

    ! There may be an "inhomogeneous", "source", or "forcing" term (useful for the OCT formalism)
    logical :: inh_term
    type(states_t) :: inh_st

    ! There may also be a exchange-like term, similar to the one necessary for time-dependent
    ! Hartree Fock, also useful only for the OCT equations
    logical :: oct_exchange
    type(states_t), pointer :: oct_st
    FLOAT, pointer :: oct_fxc(:, :, :)

    CMPLX, pointer :: phase(:, :)

    FLOAT :: current_time
    logical :: apply_packed  !< This is initialized by the StatesPack variable.
  end type hamiltonian_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter :: &
    NOT_ABSORBING       = 0,    &
    IMAGINARY_ABSORBING = 1,    &
    MASK_ABSORBING      = 2,    &
    EXACT_ABSORBING     = 3

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 2, &
    HARTREE               = 1, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4, &
    CLASSICAL             = 5

  type(profile_t), save :: prof_hamiltonian, prof_vlpsi, prof_kinetic

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_init(hm, gr, geo, st, theory_level, xc_family)
    type(hamiltonian_t),    intent(out)   :: hm
    type(grid_t),           intent(inout) :: gr
    type(geometry_t),       intent(inout) :: geo
    type(states_t), target, intent(inout) :: st
    integer,                intent(in)    :: theory_level
    integer,                intent(in)    :: xc_family

    integer :: iline, icol, ispin
    type(states_dim_t), pointer :: states_dim
    integer :: ncols
    type(block_t) :: blk

    PUSH_SUB(hamiltonian_init)

    states_dim => st%d

    ! make a couple of local copies
    hm%theory_level = theory_level
    hm%xc_family    = xc_family
    call states_dim_copy(hm%d, states_dim)

    call hamiltonian_base_init(hm%hm_base, gr%mesh, hm%d%nspin)
    ASSERT(associated(gr%der%lapl))
    hm%hm_base%kinetic => gr%der%lapl

    SAFE_ALLOCATE(hm%energy)

    ! initialize variables
    hm%energy%intnvxc = M_ZERO
    hm%energy%exchange = M_ZERO
    hm%energy%correlation = M_ZERO
    hm%energy%total = M_ZERO
    hm%energy%kinetic = M_ZERO

    nullify(hm%oct_fxc)

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = hm%vxc(:, 1), vxc_22 = hm%vxc(:, 2), Re(vxc_12) = hm%vxc(:. 3);
    ! Im(vxc_12) = hm%vxc(:, 4)
    SAFE_ALLOCATE(hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin))
    hm%vhxc(1:gr%mesh%np, 1:hm%d%nspin) = M_ZERO

    if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
      SAFE_ALLOCATE(hm%vhartree(1:gr%mesh%np))
      SAFE_ALLOCATE(hm%vxc(1:gr%mesh%np, 1:hm%d%nspin))
      if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) then
        SAFE_ALLOCATE(hm%vtau(1:gr%mesh%np, 1:hm%d%nspin))
      else
        nullify(hm%vtau)
      end if

      hm%vhartree(1:gr%mesh%np) = M_ZERO

      do ispin = 1, hm%d%nspin
        hm%vxc(1:gr%mesh%np, ispin) = M_ZERO
        if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) hm%vtau(1:gr%mesh%np, ispin) = M_ZERO
      end do

      if (hm%d%cdft) then
        SAFE_ALLOCATE(hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:hm%d%nspin))
        hm%axc = M_ZERO
      else
        nullify(hm%axc)
      end if
    end if

    !Initialize external potential
    call epot_init(hm%ep, gr, geo, hm%d%ispin, hm%d%nik)

    nullify(hm%vberry)
    if(associated(hm%ep%E_field) .and. simul_box_is_periodic(gr%sb)) then
      ! only need vberry if there is a field in a periodic direction
      if(any(abs(hm%ep%E_field(1:gr%sb%periodic_dim)) > M_EPSILON)) then
        SAFE_ALLOCATE(hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
      endif
    endif

    !Static magnetic field requires complex wavefunctions
    if (associated(hm%ep%B_field) .or. gauge_field_is_applied(hm%ep%gfield)) call states_set_complex(st)

    call parse_logical(datasets_check('CalculateSelfInducedMagneticField'), .false., hm%self_induced_magnetic)
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
    !% calculated at the end of a <tt>gs</tt> calculation (nothing is done -- yet -- in the <tt>td </tt>case)
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

    !%Variable AbsorbingBoundaries
    !%Type integer
    !%Default not_absorbing
    !%Section Time-Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of
    !% standing density waves, one can make the boundaries of the simulation
    !% box absorbing.
    !%Option not_absorbing 0
    !% No absorbing boundaries.
    !%Option sin2 1
    !% A <math>\sin^2</math> imaginary potential is added at the boundaries.
    !%Option mask 2
    !% A mask is applied to the wavefunctions at the boundaries.
    !%Option exact 3
    !% NOT WORKING YET!
    !% An exactly absorbing scheme is used for open boundaries. This feature
    !% comes from transport calculation and assumes that on <tt>OpenBoundariesNLeads</tt>
    !% sides there is a lead connected. No outgoing density is reflected within the leads,
    !% but some minor reflection will occur on the corners of the box.
    !% This is due to the setup of semi-infinite finite width leads connected to the sides.
    !% Warning: This scheme works only with the special Cranck-Nicholson propagator and has
    !% quadratic scaling with time. It may be tuned with the parameter <tt>OpenBoundariesMaxMemCoeffs</tt>.
    !%End
    call parse_integer(datasets_check('AbsorbingBoundaries'), NOT_ABSORBING, hm%ab)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', hm%ab)) call input_error('AbsorbingBoundaries')
    call messages_print_var_option(stdout, "AbsorbingBoundaries", hm%ab)

    nullify(hm%ab_pot)

    if(hm%ab .ne. NOT_ABSORBING) call init_abs_boundaries()

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
    call parse_float(datasets_check('ParticleMass'), M_ONE, hm%mass)

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
    if(parse_block(datasets_check('MassScaling'), blk) == 0) then
        ncols = parse_block_cols(blk, 0)
        if(ncols > gr%sb%dim) then
          call input_error("MassScaling")
        end if
        iline = 1 ! just deal with 1 line - should be generalized
        do icol = 1, ncols
          call parse_block_float(blk, iline - 1, icol - 1, hm%mass_scaling(icol))
        end do
        call parse_block_end(blk)
    end if

    SAFE_ALLOCATE(hm%hf_st)
    call states_null(hm%hf_st)

    hm%inh_term = .false.
    call hamiltonian_remove_oct_exchange(hm)

    hm%adjoint = .false.

    nullify(hm%phase)
    if (simul_box_is_periodic(gr%sb) .and. &
        .not. (kpoints_number(gr%sb%kpoints) == 1 .and. kpoints_point_is_gamma(gr%sb%kpoints, 1))) &
      call init_phase()
    ! no e^ik phase needed for Gamma-point-only periodic calculations

    if(gr%ob_grid%open_boundaries) then
      if(hm%theory_level .ne. INDEPENDENT_PARTICLES) then
        message(1) = 'Open-boundary calculations for interacting electrons are'
        message(2) = 'not yet possible.'
        call messages_fatal(2)
      end if
      call init_lead_h
    end if

    !%Variable StatesPack
    !%Type logical
    !%Default yes
    !%Section Execution::Optimization
    !%Description
    !% If set to yes (the default), Octopus will 'pack' the
    !% wave-functions when operating with them. This involves some
    !% additional copying but makes operations more efficient.
    !%End
    call parse_logical(datasets_check('StatesPack'), .true., hm%apply_packed)

    POP_SUB(hamiltonian_init)

  contains

    ! ---------------------------------------------------------
    subroutine init_phase
      integer :: ip, ik
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_init.init_phase)

      SAFE_ALLOCATE(hm%phase(1:gr%mesh%np_part, hm%d%kpt%start:hm%d%kpt%end))

      kpoint(1:gr%sb%dim) = M_ZERO
      do ik = hm%d%kpt%start, hm%d%kpt%end
        kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
        forall (ip = 1:gr%mesh%np_part)
          hm%phase(ip, ik) = exp(-M_zI * sum(gr%mesh%x(ip, 1:gr%sb%dim) * kpoint(1:gr%sb%dim)))
        end forall
      end do

      POP_SUB(hamiltonian_init.init_phase)
    end subroutine init_phase


    ! ---------------------------------------------------------
    subroutine init_abs_boundaries()
      FLOAT :: dd
      integer :: ip

      PUSH_SUB(hamiltonian_init.init_abs_boundaries)

      !%Variable ABWidth
      !%Type float
      !%Default 0.4 a.u.
      !%Section Time-Dependent::Absorbing Boundaries
      !%Description
      !% Width of the region used to apply the absorbing boundaries.
      !%End
      call parse_float(datasets_check('ABWidth'), CNST(0.4), hm%ab_width, units_inp%length)

      if(hm%ab == 1) then
        !%Variable ABHeight
        !%Type float
        !%Default -0.2 a.u.
        !%Section Time-Dependent::Absorbing Boundaries
        !%Description
        !% When <tt>AbsorbingBoundaries = sin2</tt>, this is the height of the imaginary potential.
        !%End
        call parse_float(datasets_check('ABHeight'), -CNST(0.2), hm%ab_height, units_inp%energy)
      else
        hm%ab_height = M_ONE
      end if

      ! generate boundary potential...
      SAFE_ALLOCATE(hm%ab_pot(1:gr%mesh%np))
      hm%ab_pot = M_ZERO
      do ip = 1, gr%mesh%np
        if(mesh_inborder(gr%mesh, geo, ip, dd, hm%ab_width)) then
          hm%ab_pot(ip) = hm%ab_height * sin(dd * M_PI / (M_TWO * hm%ab_width))**2
        end if
      end do

      POP_SUB(hamiltonian_init.init_abs_boundaries)
    end subroutine init_abs_boundaries


    ! ---------------------------------------------------------
    ! Calculate the blocks of the lead Hamiltonian and read the potential
    ! of the lead unit cell.
    subroutine init_lead_h
      integer               :: np, np_part, il, ierr, pot, ix, iy, is
      integer               :: irow, diag, offdiag
      character             :: channel
      character(len=256)    :: fname, fmt, static_dir
      character(len=6)      :: name
      type(mesh_t), pointer :: mesh
      logical               :: t_inv

      PUSH_SUB(hamiltonian_init.init_lead_h)

      ! Read potential of the leads. We try v0 (for non-interacting electrons)
      ! (Octopus binary and NetCDF format). If DFT (without pseudopotentials)
      ! is used we try to read vh and vks. If one is not existing a warning
      ! is written and a zero potential is assumed.
      ! \todo: spinors
      do il = 1, NLEADS
        np      = gr%intf(il)%np_uc
        np_part = gr%intf(il)%np_part_uc
        static_dir = gr%ob_grid%lead(il)%info%static_dir
        mesh => gr%ob_grid%lead(il)%mesh
        call lead_init_pot(hm%lead(il), np, np_part, hm%d%nspin)

        name = 'v0'
        fname = trim(static_dir)//'/'//trim(name)
        call read_potential(fname, mesh, hm%lead(il)%v0(:), trim(name), il)
        if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
          name = 'vh'
          fname = trim(static_dir)//'/'//trim(name)
          call read_potential(fname, mesh, hm%lead(il)%vh(:), trim(name), il)
          ! not sure if this is correct as the potentials are only output up to is=2 in output_h_inc.F90
          do is = 1, hm%d%ispin
            if(hm%d%ispin == 1) then
              write(name, '(a)') 'vks'
            else
              write(name, '(a,i1)') 'vks-sp', is
            endif
            fname = trim(static_dir)//'/'//trim(name)
            call read_potential(fname, mesh, hm%lead(il)%vks(:, is), trim(name), il)
          end do
        else
          do is = 1, hm%d%ispin
            hm%lead(il)%vks(:, is) = hm%lead(il)%v0(:)
          end do
          hm%lead(il)%vh(:)     = M_ZERO
        end if
      end do

      ! Calculate the diagonal and offdiagonal blocks of the lead Hamiltonian.
      ! First check if we can reduce the size of the interface.
      ! If there is no dependence in transport direction in the Kohn-Sham
      ! potential then it is possible to reduce the size of the unit cell
      ! to the size of the interface itself.
      do il = 1, NLEADS
        t_inv = .true.
        do ispin = 1, hm%d%nspin
          ! \todo generalize to find the smallest periodicity of the leads
          t_inv = t_inv.and.is_lead_transl_inv(gr%der%lapl, hm%lead(il)%vks(:, ispin), gr%intf(il))
        end do
        if (t_inv) then
          if(in_debug_mode) then ! write some info
            write(message(1), '(a,2i8)') 'Reallocate lead unit cell from ', gr%intf(il)%np_uc
            call messages_info(1)
          end if
          ! resize array
          ! so delete the old array intf%index
          call interface_end(gr%intf(il))
          ! then re-initialize interface
          ! \todo generalize to find the smallest periodicity of the leads
          call interface_init(gr%der, gr%intf(il), il, gr%ob_grid%lead(il)%sb%lsize, &
                              derivatives_stencil_extent(gr%der, (il+1)/2))
          np = gr%intf(il)%np_uc
          np_part =  gr%intf(il)%np_part_uc
          call lead_init_kin(hm%lead(il), np, np_part, st%d%dim)
          call lead_resize(gr%intf(il), hm%lead(il), st%d%dim, hm%d%nspin)
          if(in_debug_mode) then ! write some info
            write(message(1), '(a,2i8)') 'to ', gr%intf(il)%np_uc
            call messages_info(1)
          end if
        end if
      end do

      do il = 1, NLEADS
        np = gr%intf(il)%np_uc
        do ispin = 1, hm%d%nspin
          call lead_diag(gr%der%lapl, hm%lead(il)%vks(:, ispin), &
            gr%intf(il), hm%lead(il)%h_diag(:, :, ispin))
          ! In debug mode write the diagonal block to a file.
          if(in_debug_mode) then
            call io_mkdir('debug/open_boundaries')
            write(fname, '(3a,i1.1,a)') 'debug/open_boundaries/diag-', &
              trim(LEAD_NAME(il)), '-', ispin, '.real'
            diag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
            write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
            do irow = 1, np
              write(diag, fmt) real(hm%lead(il)%h_diag(:, irow, ispin))
            end do
            call io_close(diag)
            write(fname, '(3a,i1.1,a)') 'debug/open_boundaries/diag-', &
              trim(LEAD_NAME(il)), '-', ispin, '.imag'
            diag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
            write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
            do irow = 1, np
              write(diag, fmt) aimag(hm%lead(il)%h_diag(:, irow, ispin))
            end do
            call io_close(diag)
          end if
        end do
        call lead_offdiag(gr%der%lapl, gr%intf(il), hm%lead(il)%h_offdiag(:, :))
        if(in_debug_mode) then
          write(fname, '(3a)') 'debug/open_boundaries/offdiag-', &
            trim(LEAD_NAME(il)), '.real'
          offdiag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
          write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
          do irow = 1, np
            write(offdiag, fmt) real(hm%lead(il)%h_offdiag(:, irow))
          end do
          call io_close(offdiag)
          write(fname, '(3a)') 'debug/open_boundaries/offdiag-', &
            trim(LEAD_NAME(il)), '.imag'
          offdiag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
          write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
          do irow = 1, np
            write(offdiag, fmt) aimag(hm%lead(il)%h_offdiag(:, irow))
          end do
          call io_close(offdiag)
        end if
      end do

      POP_SUB(hamiltonian_init.init_lead_h)
    end subroutine init_lead_h

    subroutine read_potential(fname, mesh, pot, potname, il)
      character(len=256),    intent(in)  :: fname
      type(mesh_t), pointer, intent(in)  :: mesh
      FLOAT,                 intent(out) :: pot(:)
      character(len=*),      intent(in)  :: potname
      integer,               intent(in)  :: il

      integer :: ierr

      PUSH_SUB(hamiltonian.hamiltonian_init.read_potential)

      ! first try obf format
      call dio_function_input(trim(fname)//'.obf', mesh, pot, ierr)
      if(ierr.eq.0) then
        message(1) = 'Info: Successfully read '//trim(potname)//' potential of the '&
                     //trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.obf'//'.'
        call messages_info(1)
      else  ! try ncdf format
        call dio_function_input(trim(fname)//'.ncdf', mesh, pot, ierr)
        if(ierr.eq.0) then
          message(1) = 'Info: Successfully read '//trim(potname)//' potential of the '&
                     //trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.ncdf'//'.'
          call messages_info(1)
        end if
      end if
      if (ierr .ne. 0) then
        ! Reading potential failed.
        message(1) = 'Could not read '//trim(potname)//' potential from the file'
        message(2) = trim(fname)
        message(3) = 'for the '//trim(LEAD_NAME(il))//' lead.'
        message(4) = 'Please include'
        message(5) = ''
        message(6) = '  Output = potential'
        message(7) = ''
        message(8) = 'in your periodic run. Octopus now assumes zero potential for'
        message(9) = trim(potname)//' in the leads. This is most likely not what you want.'
        call messages_warning(9)
        pot(:) = M_ZERO
      end if
      POP_SUB(hamiltonian.hamiltonian_init.read_potential)
    end subroutine read_potential

  end subroutine hamiltonian_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_end(hm, gr, geo)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(inout) :: geo

    integer :: il

    PUSH_SUB(hamiltonian_end)

    call hamiltonian_base_end(hm%hm_base)

    SAFE_DEALLOCATE_P(hm%phase)
    SAFE_DEALLOCATE_P(hm%vhartree)
    SAFE_DEALLOCATE_P(hm%vhxc)
    SAFE_DEALLOCATE_P(hm%vxc)
    SAFE_DEALLOCATE_P(hm%axc)
    SAFE_DEALLOCATE_P(hm%vberry)
    SAFE_DEALLOCATE_P(hm%a_ind)
    SAFE_DEALLOCATE_P(hm%b_ind)

    if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) then
      SAFE_DEALLOCATE_P(hm%vtau)
    end if

    call epot_end(hm%ep, geo)

    SAFE_DEALLOCATE_P(hm%ab_pot)

    if(gr%ob_grid%open_boundaries) then
      do il = 1, NLEADS
        call lead_end(hm%lead(il))
      end do
    end if

    call states_dim_end(hm%d)
    call states_end(hm%hf_st)
    SAFE_DEALLOCATE_P(hm%hf_st)

    SAFE_DEALLOCATE_P(hm%energy)

    POP_SUB(hamiltonian_end)
  end subroutine hamiltonian_end


  ! ---------------------------------------------------------
  ! True if the Hamiltonian is Hermitian, false otherwise
  logical function hamiltonian_hermitian(hm)
    type(hamiltonian_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_hermitian)
    hamiltonian_hermitian = .not.((hm%ab .eq. IMAGINARY_ABSORBING) .or. hamiltonian_oct_exchange(hm))

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
  logical function hamiltonian_oct_exchange(hm) result(oct_exchange)
    type(hamiltonian_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_oct_exchange)
    oct_exchange = hm%oct_exchange

    POP_SUB(hamiltonian_oct_exchange)
  end function hamiltonian_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_oct_exchange(hm, st, gr, xc)
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t), target, intent(in) :: st
    type(grid_t),        intent(in)    :: gr
    type(xc_t),          intent(in)    :: xc

    integer :: np, nspin

    PUSH_SUB(hamiltonian_set_oct_exchange)

    ! In this release, no non-local part for the QOCT Hamiltonian.
    nullify(hm%oct_st)

    hm%oct_st => st
    hm%oct_exchange = .true.
    np = gr%mesh%np
    nspin = hm%oct_st%d%nspin

    SAFE_ALLOCATE(hm%oct_fxc(1:np, 1:nspin, 1:nspin))
    hm%oct_fxc = M_ZERO
    call xc_get_fxc(xc, gr%mesh, st%rho, st%d%ispin, hm%oct_fxc)

    POP_SUB(hamiltonian_set_oct_exchange)
  end subroutine hamiltonian_set_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_oct_exchange(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_remove_oct_exchange)

    nullify(hm%oct_st)
    hm%oct_exchange = .false.
    SAFE_DEALLOCATE_P(hm%oct_fxc)

    POP_SUB(hamiltonian_remove_oct_exchange)
  end subroutine hamiltonian_remove_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_adjoint)

    if(.not.hm%adjoint) then
      hm%adjoint = .true.
      if(hm%ab .eq. IMAGINARY_ABSORBING) then
        hm%ab_pot = -hm%ab_pot
      end if
    endif

    POP_SUB(hamiltonian_adjoint)
  end subroutine hamiltonian_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_not_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_not_adjoint)

    if(hm%adjoint) then
      hm%adjoint = .false.
      if(hm%ab .eq. IMAGINARY_ABSORBING) then
        hm%ab_pot = -hm%ab_pot
      end if
    endif

    POP_SUB(hamiltonian_not_adjoint)
  end subroutine hamiltonian_not_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_update(this, mesh, time)
    type(hamiltonian_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT, optional,     intent(in)    :: time

    integer :: ispin, ip, idir, iatom, ilaser
    type(profile_t), save :: prof
    FLOAT :: aa(1:MAX_DIM)
    FLOAT, allocatable :: vp(:,:)

    PUSH_SUB(hamiltonian_update)
    call profiling_in(prof, "HAMILTONIAN_UPDATE")

    if(present(time)) this%current_time = time

    ! set everything to zero
    call hamiltonian_base_clear(this%hm_base)

    ! the xc, hartree and external potentials
    call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_POTENTIAL)
    do ispin = 1, this%d%nspin
      if(ispin <= 2) then
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin) + this%ep%vpsl(ip)
      else
        forall (ip = 1:mesh%np) this%hm_base%potential(ip, ispin) = this%vhxc(ip, ispin)
      end if
    end do

    ! the lasers
    if (present(time)) then

      do ilaser = 1, this%ep%no_lasers
        select case(laser_kind(this%ep%lasers(ilaser)))
        case(E_FIELD_SCALAR_POTENTIAL, E_FIELD_ELECTRIC)
          do ispin = 1, this%d%nspin
            call laser_potential(this%ep%lasers(ilaser), mesh,  this%hm_base%potential(:, ispin), time)
          end do
        case(E_FIELD_MAGNETIC)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL + FIELD_UNIFORM_MAGNETIC_FIELD)
          ! get the vector potential
          SAFE_ALLOCATE(vp(1:mesh%np, 1:mesh%sb%dim))
          vp(1:mesh%np, 1:mesh%sb%dim) = M_ZERO
          call laser_vector_potential(this%ep%lasers(ilaser), mesh, vp, time)
          forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
            this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) - vp(ip, idir)/P_C
          end forall
          ! and the magnetic field
          call laser_field(this%ep%lasers(ilaser), this%hm_base%uniform_magnetic_field(1:mesh%sb%dim), time)
          SAFE_DEALLOCATE_A(vp)
        case(E_FIELD_VECTOR_POTENTIAL)
          call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL)
          ! get the uniform vector potential associated to a magnetic field
          aa = M_ZERO
          call laser_field(this%ep%lasers(ilaser), aa(1:mesh%sb%dim), time)
          this%hm_base%uniform_vector_potential = this%hm_base%uniform_vector_potential - aa/P_C
        end select
      end do

    end if

    ! the gauge field
    if(gauge_field_is_applied(this%ep%gfield)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_VECTOR_POTENTIAL)
      this%hm_base%uniform_vector_potential = this%hm_base%uniform_vector_potential + gauge_field_get_vec_pot(this%ep%gfield)/P_c
    end if

    ! the vector potential of a static magnetic field
    if(associated(this%ep%a_static)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_VECTOR_POTENTIAL)
      forall (idir = 1:mesh%sb%dim, ip = 1:mesh%np)
        this%hm_base%vector_potential(idir, ip) = this%hm_base%vector_potential(idir, ip) + this%ep%a_static(ip, idir)
      end forall
    end if

    ! and the static magnetic field
    if(associated(this%ep%b_field)) then
      call hamiltonian_base_allocate(this%hm_base, mesh, FIELD_UNIFORM_MAGNETIC_FIELD)
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
      integer :: ik, imat, nmat, max_npoints
      FLOAT   :: kpoint(1:MAX_DIM)

      PUSH_SUB(hamiltonian_update.build_phase)

      ! now regenerate the phases for the pseudopotentials
      do iatom = 1, this%ep%natoms
        call projector_init_phases(this%ep%proj(iatom), mesh%sb, this%d, &
          vec_pot = this%hm_base%uniform_vector_potential, vec_pot_var = this%hm_base%vector_potential)
      end do

      if(associated(this%hm_base%uniform_vector_potential)) then
        if(.not. associated(this%phase)) then
          SAFE_ALLOCATE(this%phase(1:mesh%np_part, this%d%kpt%start:this%d%kpt%end))
        end if

        kpoint(1:mesh%sb%dim) = M_ZERO
        do ik = this%d%kpt%start, this%d%kpt%end
          kpoint(1:mesh%sb%dim) = kpoints_get_point(mesh%sb%kpoints, states_dim_get_kpoint_index(this%d, ik))

          forall (ip = 1:mesh%np_part)
            this%phase(ip, ik) = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*(kpoint(1:mesh%sb%dim) &
              + this%hm_base%uniform_vector_potential(1:mesh%sb%dim))))
          end forall
        end do
      end if

      max_npoints = this%hm_base%max_npoints
      nmat = this%hm_base%nprojector_matrices


      if(associated(this%phase) .and. associated(this%hm_base%projector_matrices)) then

        if(.not. associated(this%hm_base%projector_phases)) then
          SAFE_ALLOCATE(this%hm_base%projector_phases(1:max_npoints, nmat, this%d%kpt%start:this%d%kpt%end))
        end if

        do ik = this%d%kpt%start, this%d%kpt%end
          do imat = 1, this%hm_base%nprojector_matrices
            iatom = this%hm_base%projector_to_atom(imat)
            do ip = 1, this%hm_base%projector_matrices(imat)%npoints
              this%hm_base%projector_phases(ip, imat, ik) = this%ep%proj(iatom)%phase(ip, ik)
            end do
          end do
        end do
      end if

      POP_SUB(hamiltonian_update.build_phase)
    end subroutine build_phase

  end subroutine hamiltonian_update


  ! ---------------------------------------------------------
  subroutine hamiltonian_epot_generate(this, gr, geo, st, time)
    type(hamiltonian_t),   intent(inout) :: this
    type(grid_t),          intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(states_t),        intent(inout) :: st
    FLOAT,       optional, intent(in)    :: time

    integer :: np

    PUSH_SUB(hamiltonian_epot_generate)

    if (st%open_boundaries) then
      np = gr%ob_grid%lead(LEFT)%mesh%np
      this%ep%vpsl_lead(1:np, LEFT) = this%lead(LEFT)%v0(1:np)
      np = gr%ob_grid%lead(RIGHT)%mesh%np
      this%ep%vpsl_lead(1:np, RIGHT) = this%lead(RIGHT)%v0(1:np)
    end if

    if(present(time)) then
      call epot_generate(this%ep, gr, geo, st, time)
    else
      call epot_generate(this%ep, gr, geo, st)
    end if
    call hamiltonian_base_build_proj(this%hm_base, gr%mesh, this%ep, geo)
    call hamiltonian_update(this, gr%mesh, time)

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
    if(associated(this%phase)) apply = .false.
    if(hamiltonian_base_has_magnetic(this%hm_base)) apply = .false.
    if(this%ab .eq. IMAGINARY_ABSORBING) apply = .false.
    if(this%theory_level == HARTREE .or. this%theory_level == HARTREE_FOCK) apply = .false.
    if(iand(this%xc_family, XC_FAMILY_MGGA).ne.0)  apply = .false.
    if(this%ep%non_local .and. .not. this%hm_base%apply_projector_matrices) apply = .false.

  end function hamiltonian_apply_packed

  ! -----------------------------------------------

  subroutine energy_copy(ein, eout)
    type(energy_t), intent(in)  :: ein
    type(energy_t), intent(out) :: eout

    PUSH_SUB(energy_copy)

    eout%total = ein%total
    eout%eigenvalues = ein%eigenvalues
    eout%exchange = ein%exchange
    eout%correlation = ein%correlation
    eout%xc_j = ein%xc_j
    eout%intnvxc = ein%intnvxc
    eout%hartree = ein%hartree
    eout%kinetic = ein%kinetic
    eout%extern = ein%extern
    eout%entropy = ein%entropy
    eout%ts = ein%ts
    eout%berry = ein%berry

    POP_SUB(energy_copy)
  end subroutine energy_copy

#include "undef.F90"
#include "real.F90"
#include "hamiltonian_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "hamiltonian_inc.F90"

end module hamiltonian_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
