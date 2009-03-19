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
  use calc_mode_m
  use datasets_m
  use em_field_m
  use derivatives_m
  use external_pot_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use gridhier_m
  use hardware_m
  use lalg_basic_m
  use loct_parser_m
  use XC_F90(lib_m)
  use io_m
  use io_function_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use multigrid_m
  use ob_lead_m
  use profiling_m
  use projector_m
  use simul_box_m
  use smear_m
  use scissor_m
  use states_m
  use states_dim_m
  use units_m
  use varinfo_m
  use lasers_m
  use poisson_m
  
  implicit none

  private
  public ::                          &
    hamiltonian_t,                   &
    hamiltonian_init,                &
    hamiltonian_mg_init,             &
    hamiltonian_end,                 &
    hamiltonian_span,                &
    dhamiltonian_apply,              &
    zhamiltonian_apply,              &
    dhamiltonian_apply_batch,        &
    zhamiltonian_apply_batch,        &
    dhamiltonian_diagonal,           &
    zhamiltonian_diagonal,           &
    dvlpsi_batch,                    &
    zvlpsi_batch,                    &
    dvnlpsi,                         &
    zvnlpsi,                         &
    dvnlpsi_batch,                   &
    zvnlpsi_batch,                   &
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
    hamiltonian_hermitean,           &
    hamiltonian_epot_generate,       &
    hamiltonian_update_potential,    &
    dvexternal,                      &
    zvexternal

  type hamiltonian_t
    ! The Hamiltonian must know what are the "dimensions" of the spaces,
    ! in order to be able to operate on the states.
    type(states_dim_t) :: d

    FLOAT, pointer :: vhartree(:) ! hartree potential
    FLOAT, pointer :: vxc(:,:)    ! xc potential
    FLOAT, pointer :: vhxc(:,:)   ! xc potential + hartree potential
    FLOAT, pointer :: axc(:,:,:)  ! xc vector-potential divided by c
    FLOAT, pointer :: vtau(:,:)   ! Derivative of e_xc w.r.t. tau

    FLOAT :: exx_coef ! how much of exx to mix

    ! The self-induced potential vector and magnetic field
    logical :: self_induced_magnetic
    FLOAT, pointer :: a_ind(:, :)
    FLOAT, pointer :: b_ind(:, :)

    ! Energies
    FLOAT :: etot,    &  ! Total energy E = Eii + Sum[Eigenvalues] - U + Ex + Ec - Int[n v_xc]
             eeigen,  &  ! Sum[Eigenvalues]
             ex,      &  ! Exchange     Ex
             ec,      &  ! Correlation  Ec
             exc_j,   &  ! 
             epot,    &  ! Int[n vxc]   
             ehartree,&  ! Hartree      U = (1/2)*Int [n v_Hartree]
             t0,      &  ! Kinetic energy of the non-interacting (KS) system of electrons
             eext,    &  ! External     V = <Phi|V|Phi> = Int[n v] (if no non-local pseudos exist)
             entropy     ! Entropy (-TS)

    integer :: theory_level    ! copied from sys%ks
    integer :: xc_family       ! copied from sys%ks

    type(epot_t) :: ep         ! handles the external potential

    ! gauge
    integer :: gauge              ! in which gauge shall we work in

    ! absorbing boundaries
    logical :: adjoint
    integer  :: ab                ! do we have absorbing boundaries?
    FLOAT :: ab_width             ! width of the absorbing boundary
    FLOAT :: ab_height            ! height of the absorbing boundary
    FLOAT, pointer :: ab_pot(:)   ! where we store the ab potential

    ! Open boundaries.
    CMPLX, pointer :: lead_h_diag(:, :, :, :)      ! Diagonal block of the lead Hamiltonian.
    CMPLX, pointer :: lead_h_offdiag(:, :, :)      ! Offdiagonal block of the lead Hamiltonian.
    FLOAT, pointer :: lead_vks(:, :, :)            ! (np, nspin, nleads) Kohn-Sham potential of the leads.
    FLOAT, pointer :: lead_vhartree(:, :)          ! (np, nleads) Hartree potential of the leads.
    CMPLX, pointer :: lead_green(:, :, :, :, :, :) ! (np, np, nspin, ncs, nik, nleads) Green function of the leads.
    
    ! Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span

    ! Kinetic Cutoff
    FLOAT :: cutoff

    ! Mass of the particle (in most cases, mass = 1, electron mass)
    FLOAT :: mass
    ! anisotropic scaling factor for the mass: different along x,y,z etc...
    FLOAT :: mass_scaling(MAX_DIM)

    ! For the Hartree Fock Hamiltonian, the Fock operator depends on the states.
    type(states_t) :: st

    ! There may be an "inhomogeneous", "source", or "forcing" term (useful for the OCT formalism)
    logical :: inh_term
    type(states_t), pointer :: inh_st

    ! There may also be a exchange-like term, similar to the one necessary for time-dependent
    ! Hartree Fock, also useful only for the OCT equations
    logical :: oct_exchange
    type(states_t), pointer :: oct_st

    CMPLX, pointer :: phase(:, :)

    logical :: multigrid_initialized
    type(dgridhier_t) :: coarse_v

    type(em_field_t), pointer :: total(:) ! one electromagnetic field per spin channel
    type(scissor_t) :: scissor
  end type hamiltonian_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter :: &
    NOT_ABSORBING       = 0,    &
    IMAGINARY_ABSORBING = 1,    &
    MASK_ABSORBING      = 2

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 2, &
    HARTREE               = 1, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_init(hm, gr, geo, st, theory_level, xc_family)
    type(hamiltonian_t),    intent(out)   :: hm
    type(grid_t),           intent(inout) :: gr
    type(geometry_t),       intent(inout) :: geo
    type(states_t), target, intent(inout) :: st
    integer,                intent(in)    :: theory_level
    integer,                intent(in)    :: xc_family

    integer                     :: i, j, n, ispin
    integer, pointer            :: wfs_type
    FLOAT                       :: d(MAX_DIM)
    type(states_dim_t), pointer :: states_dim

    integer :: ncols
    type(block_t) :: blk

    call push_sub('hamiltonian.hamiltonian_init')

    states_dim => st%d
    wfs_type   => st%wfs_type

    ! make a couple of local copies
    hm%theory_level = theory_level
    hm%xc_family    = xc_family
    call states_dim_copy(hm%d, states_dim)

    ! initialize variables
    hm%epot = M_ZERO
    hm%ex = M_ZERO; hm%ec = M_ZERO
    hm%etot = M_ZERO

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = hm%vxc(:, 1), vxc_22 = hm%vxc(:, 2), Re(vxc_12) = hm%vxc(:. 3);
    ! Im(vxc_12) = hm%vxc(:, 4)
    ALLOCATE(hm%vhartree(NP),        NP)
    ALLOCATE(hm%vxc(NP, hm%d%nspin), NP*hm%d%nspin)
    ALLOCATE(hm%vhxc(NP, hm%d%nspin), NP*hm%d%nspin)
    if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) then
      ALLOCATE(hm%vtau(NP, hm%d%nspin), NP*hm%d%nspin)
    else
      nullify(hm%vtau)
    end if

    hm%vhartree(1:NP) = M_ZERO

    do ispin = 1, hm%d%nspin
      hm%vhxc(1:NP, ispin) = M_ZERO
      hm%vxc(1:NP, ispin) = M_ZERO
      if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) hm%vtau(1:NP, ispin) = M_ZERO
    end do

    if (hm%d%cdft) then
      ALLOCATE(hm%axc(NP, gr%mesh%sb%dim, hm%d%nspin), NP*gr%mesh%sb%dim*hm%d%nspin)
      hm%axc = M_ZERO
    else
      nullify(hm%axc)
    end if

    !Initialize external potential
    call epot_init(hm%ep, gr, geo, hm%d%ispin)

    !Static magnetic field requires complex wave-functions
    if (associated(hm%ep%B_field) .or. gauge_field_is_applied(hm%ep%gfield)) wfs_type = M_CMPLX

    call loct_parse_logical(datasets_check('CalculateSelfInducedMagneticField'), .false., hm%self_induced_magnetic)
    !%Variable CalculateSelfInducedMagneticField
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% The existence of an electronic current implies the creation of a self-induced magnetic
    !% field, which may in turn back-react on the system. Of course, a fully consistent treatment
    !% of this kind of effects should be done of QED theory, but we will attempt a first
    !% approximation to the problem by considering the lowest order relativistic terms
    !% plugged in the normal Hamiltonian equations (spin-other-orbit coupling terms, etc). 
    !% For the moment being, none of this is done, but a first step is taken by calculating
    !% the induced magnetic field of a system that has a current, by considering the magnetostatic
    !% approximation and Biot-Savart law:
    !%
    !% <math> \nabla^2 \vec{A} + 4\pi\alpha \vec{J} = 0</math>
    !%
    !% <math> \vec{B} = \vec{\nabla} \times \vec{A}</math>
    !%
    !% If CalculateSelfInducedMagneticField is set to yes, this <math> B </math> field is
    !% calculated at the end of a gs calculation (nothing is done -- yet -- in the td case)
    !% and printed out, if the Output variable contains the "potential" keyword (the prefix
    !% of the output files are "Bind").
    !%End
    if(hm%self_induced_magnetic) then
      select case(gr%mesh%sb%dim)
      case(3)
        ALLOCATE(hm%a_ind(NP_PART, MAX_DIM), NP_PART*MAX_DIM)
        ALLOCATE(hm%b_ind(NP_PART, MAX_DIM), NP_PART*MAX_DIM)
      case(2)
        ALLOCATE(hm%a_ind(NP_PART, 2), NP_PART*2)
        ALLOCATE(hm%b_ind(NP_PART, 1), NP_PART)
      end select
    else
      nullify(hm%a_ind, hm%b_ind)
    end if

    !%Variable TDGauge
    !%Type integer
    !%Default length
    !%Section Time Dependent
    !%Description
    !% In which gauge to treat the laser field.
    !%Option length 1
    !% Length gauge.
    !%Option velocity 2
    !% Velocity gauge.
    !%End
    call loct_parse_int(datasets_check('TDGauge'), LENGTH, hm%gauge)
    if(.not.varinfo_valid_option('TDGauge', hm%gauge)) call input_error('TDGauge')
    call messages_print_var_option(stdout, "TDGauge", hm%gauge)

    !%Variable AbsorbingBoundaries
    !%Type integer
    !%Default not_absorbing
    !%Section Time Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of 
    !% standing density waves, one can make the boundaries of the simulation 
    !% box absorbing.
    !%Option not_absorbing 0
    !% No absorbing boundaries.
    !%Option sin2 1
    !% A <math>\sin^2</math> imaginary potential is added at the boundaries.
    !%Option mask 2
    !% A mask is applied to the wave-functions at the boundaries.
    !%End
    call loct_parse_int(datasets_check('AbsorbingBoundaries'), NOT_ABSORBING, hm%ab)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', hm%ab)) call input_error('AbsorbingBoundaries')
    call messages_print_var_option(stdout, "AbsorbingBoundaries", hm%ab)
    
    nullify(hm%ab_pot)

    absorbing_boundaries: if(hm%ab.ne.NOT_ABSORBING) then
      !%Variable ABWidth
      !%Type float
      !%Default 0.4 a.u.
      !%Section Time Dependent::Absorbing Boundaries
      !%Description
      !% Width of the region used to apply the absorbing boundaries.
      !%End
      call loct_parse_float(datasets_check('ABWidth'), CNST(0.4)/units_inp%length%factor, hm%ab_width)
      hm%ab_width  = hm%ab_width * units_inp%length%factor
      if(hm%ab == 1) then
        !%Variable ABHeight
        !%Type float
        !%Default -0.2 a.u.
        !%Section Time Dependent::Absorbing Boundaries
        !%Description 
        !% When <tt>AbsorbingBoundaries == sin2</tt>, is the height of the imaginary potential.
        !%End
        call loct_parse_float(datasets_check('ABHeight'), -CNST(0.2)/units_inp%energy%factor, hm%ab_height)
        hm%ab_height = hm%ab_height * units_inp%energy%factor
      else
        hm%ab_height = M_ONE
      end if

      ! generate boundary potential...
      ALLOCATE(hm%ab_pot(NP), NP)
      hm%ab_pot = M_ZERO
      do i = 1, NP
        call mesh_inborder(gr%mesh, geo, i, n, d, hm%ab_width)
        if(n>0) then
          do j = 1, n
            hm%ab_pot(i) = hm%ab_pot(i) + hm%ab_height * sin(d(j)*M_PI/(M_TWO*hm%ab_width))**2
          end do
        end if
        if(abs(hm%ab_pot(i)) > abs(hm%ab_height)) hm%ab_pot(i) = hm%ab_height
      end do

    end if absorbing_boundaries

    ! Cutoff applied to the kinetic term.
    ! it is used *both* in the calculation of the split operator method
    call loct_parse_float(datasets_check('KineticCutoff'), -M_ONE, hm%cutoff)
    if(hm%cutoff > M_ZERO) then
      hm%cutoff = hm%cutoff*units_inp%energy%factor
      write(message(1),'(a,f7.2,a)') 'Info: The kinetic operator will have a cutoff of',&
        hm%cutoff/units_out%energy%factor, units_out%energy%abbrev
      call write_info(1)
    end if

    !%Variable ParticleMass
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian
    !%Description 
    !% It is possible to make calculations for a particle with a mass
    !% different from one (atomic unit of mass, or mass of the electron).
    !% This is useful to describe non-electronic systems, of for
    !% esoteric purposes.
    !%End
    call loct_parse_float(datasets_check('ParticleMass'), M_ONE, hm%mass)

    !%Variable MassScaling
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% scaling factor for anisotropic masses (different masses along each
    !% geometric direction
    !%
    !% <tt>%MassScaling
    !% <br>&nbsp;&nbsp;1.0 | 1800.0 | 1800.0
    !% <br>%</tt>
    !%
    !% would fix the mass of the particles to be 1800 along the y and z
    !% directions. This can be useful, e.g., to simulate 3 particles in 1D,
    !% viz in this case an electron and 2 protons.
    !%
    !%End
    hm%mass_scaling = M_ONE
    if(loct_parse_block(datasets_check('MassScaling'), blk)==0) then
        ncols = loct_parse_block_cols(blk, 0)
        if(ncols > MAX_DIM) then
          call input_error("MassScaling")
        end if
        i=1 ! just deal with 1 line - should be generalized
        do j = 1, ncols
          call loct_parse_block_float(blk, i-1, j-1, hm%mass_scaling(j))
        end do
        call loct_parse_block_end(blk)
    end if

    call states_null(hm%st)

    call hamiltonian_remove_inh(hm)
    call hamiltonian_remove_oct_exchange(hm)

    hm%adjoint = .false.

    nullify(hm%phase)
    if (simul_box_is_periodic(gr%sb)) call init_phase()

    if(gr%sb%open_boundaries) then
      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        message(1) = 'Open boundary calculations for interacting electrons are'
        message(2) = 'not yet possible.'
        call write_fatal(2)
      end if
      call init_lead_h
    end if

    hm%multigrid_initialized = .false.

    ALLOCATE(hm%total(1:hm%d%nspin), hm%d%nspin)

    call scissor_nullify(hm%scissor)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine init_phase
      integer :: ip, ik

      ALLOCATE(hm%phase(1:NP_PART, 1:hm%d%nik), NP_PART*hm%d%nik)

      forall (ik = 1:hm%d%nik, ip = 1:NP_PART)
        hm%phase(ip, ik) = exp(-M_zI*sum(gr%mesh%x(ip, 1:gr%mesh%sb%dim)* hm%d%kpoints(1:gr%mesh%sb%dim, ik)))
      end forall
      
    end subroutine init_phase
    

    ! ---------------------------------------------------------
    ! Calculate the blocks of the lead Hamiltonian and read the potential
    ! of the lead unit cell.
    subroutine init_lead_h
      integer               :: np, il, ierr, alloc_size, ik, ist, pot, ix, iy
      integer               :: green_real, green_imag, irow, diag, offdiag
      character             :: channel
      character(len=1)      :: ln(NLEADS)
      character(len=2)      :: spin
      character(len=256)    :: fname, fmt, fname_real, fname_imag
      FLOAT                 :: energy
      type(mesh_t), pointer :: m

      ln(LEFT)  = 'L'; ln(RIGHT) = 'R'

      np = gr%intf(LEFT)%np

      ! Read potential of the leads. We try vks-x (for DFT without
      ! psuedo-potentials) and v0 (for non-interacting electrons) in
      ! that order (Octopus binary and NetCDF format). If none of the
      ! two can be found, a warning is emittedg and zero potential
      ! assumed.
      ALLOCATE(hm%lead_vks(np, hm%d%nspin, NLEADS), np*hm%d%nspin*NLEADS)
      ALLOCATE(hm%lead_vhartree(np, NLEADS), np*NLEADS)

      do il = 1, NLEADS
        do ispin = 1, hm%d%nspin
          write(channel, '(i1)') ispin

          ! Try vks-ispin first.
          ! OBF.
          fname = trim(gr%sb%lead_static_dir(il))//'/vks-'//trim(channel)//'.obf'
          call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vks(:, ispin, il), ierr)
          if(ierr.eq.0) then
            message(1) = 'Info: Successfully read KS potential of the '//trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
            call write_info(1)
          else
            ! NetCDF.
            fname = trim(gr%sb%lead_static_dir(il))//'/vks-'//trim(channel)//'.ncdf'
            call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vks(:, ispin, il), ierr)
            if(ierr.eq.0) then
              message(1) = 'Info: Successfully read KS potential of the '//trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
              call write_info(1)
            else
              ! Now try v0.
              ! OBF.
              fname = trim(gr%sb%lead_static_dir(il))//'/v0.obf'
              call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vks(:, ispin, il), ierr)
              if(ierr.eq.0) then
                message(1) = 'Info: Successfully read external potential of the '// &
                  trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
                call write_info(1)
              else
                ! NetCDF.
                fname = trim(gr%sb%lead_static_dir(il))//'/v0.ncdf'
                call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vks(:, ispin, il), ierr)
                if(ierr.eq.0) then
                  message(1) = 'Info: Successfully read external potential of the '// &
                    trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
                  call write_info(1)
                else
                  ! Reading potential failed.
                  message(1) = 'Could neither read vks-x nor v0 from the directory'
                  message(2) = trim(gr%sb%lead_static_dir(il))//' for the '//trim(LEAD_NAME(il))//' lead.'
                  message(3) = 'Please include'
                  message(4) = ''
                  message(5) = '  Output = potential'
                  message(6) = ''
                  message(7) = 'in your periodic run. Octopus now assumes zero potential'
                  message(8) = 'in the leads. This is most likely not what you want.'
                  call write_warning(8)
                  hm%lead_vks(:, ispin, il) = M_ZERO
                end if
              end if
            end if
          end if

          ! In debug mode, write potential to file in gnuplot format
          ! (only z=0 plane).  We cannot use doutput_function because
          ! the lead mesh is not completely initialized, in particular
          ! the x array is missing.
          if(in_debug_mode) then
            call io_mkdir('debug/open_boundaries')
            fname = 'debug/open_boundaries/v_lead-'//trim(LEAD_NAME(il))//'-'//trim(channel)
            pot = io_open(trim(fname), action='write', is_tmp=.false., grp=gr%mesh%mpi_grp)
            m => gr%mesh%lead_unit_cell(LEFT)
            do ix = m%idx%nr(1, 1)+m%idx%enlarge(1), m%idx%nr(2, 1)-m%idx%enlarge(1)
              do iy = m%idx%nr(1, 2)+m%idx%enlarge(2), m%idx%nr(2, 2)-m%idx%enlarge(2)
                write(pot, '(2i8,f16.8)') ix, iy, hm%lead_vks(m%idx%Lxyz_inv(ix, iy, 0), ispin, il)
              end do
            end do
            call io_close(pot)
          end if
        end do

        ! Read Hartree potential.
        ! OBF.
        hm%lead_vhartree(:, il) = M_ZERO
        if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
          fname = trim(gr%sb%lead_static_dir(il))//'/vh.obf'
          call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vhartree(:, il), ierr)
          if(ierr.eq.0) then
            message(1) = 'Info: Successfully read Hartree potential of the '//trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
            call write_info(1)
          else
            ! NetCDF.
            fname = trim(gr%sb%lead_static_dir(il))//'/vh.ncdf'
            call dinput_function(trim(fname), gr%mesh%lead_unit_cell(il), hm%lead_vhartree(:, il), ierr)
            if(ierr.eq.0) then
              message(1) = 'Info: Successfully read Hartree potential of the '//trim(LEAD_NAME(il))//' lead from '//trim(fname)//'.'
              call write_info(1)
            else
              message(1) = 'Could not read the Hartree potential of the leads.'
              message(2) = 'The Hartree term will not be calculated correctly.'
              message(3) = 'Include'
              message(4) = ''
              message(5) = '  Output = potential'
              message(6) = ''
              message(7) = 'in your periodic run and make sure the Hartree interaction'
              message(8) = 'is switched on.'
              call write_warning(8)
            end if
          end if
        end if
      end do

      ! Calculate the diagonal and offdiagonal blocks of the lead Hamiltonian.
      ALLOCATE(hm%lead_h_diag(np, np, st%d%dim, NLEADS), np**2*st%d%dim*NLEADS)
      ALLOCATE(hm%lead_h_offdiag(np, np, NLEADS), np**2*NLEADS)
      do il = 1, NLEADS
        do ispin = 1, hm%d%nspin
          call lead_diag(gr%der%lapl, hm%lead_vks(:, ispin, il), &
            gr%intf(il), hm%lead_h_diag(:, :, ispin, il))
          ! In debug mode write the diagonal block to a file.
          if(in_debug_mode) then
            call io_mkdir('debug/open_boundaries')
            write(fname, '(3a,i1.1)') 'debug/open_boundaries/diag-', &
              trim(LEAD_NAME(il)), '-', ispin
            diag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
            write(fmt, '(a,i6,a)') '(', gr%intf(il)%np, 'e14.4)'
            do irow = 1, gr%intf(il)%np
              write(diag, fmt) hm%lead_h_diag(:, :, ispin, il)
            end do
            call io_close(diag)
          end if
        end do
        call lead_offdiag(gr%der%lapl, gr%intf(il), il, &
          hm%lead_h_offdiag(:, :, il))
        if(in_debug_mode) then
          write(fname, '(2a)') 'debug/open_boundaries/offdiag-', &
            trim(LEAD_NAME(il))
          offdiag = io_open(fname, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
          write(fmt, '(a,i6,a)') '(', gr%intf(il)%np, 'e14.4)'
          do irow = 1, gr%intf(il)%np
            write(offdiag, fmt) hm%lead_h_offdiag(:, :, il)
          end do
          call io_close(offdiag)
        end if
      end do

      ! Calculate Green function of the leads.
      ! FIXME: For spinors, this calculation is almost certainly wrong.
      alloc_size = np**2*hm%d%nspin*st%ob_ncs*st%d%nik*NLEADS
      ALLOCATE(hm%lead_green(np, np, hm%d%nspin, st%ob_ncs, st%d%nik, NLEADS), alloc_size)

      if(calc_mode_is(CM_GS)) then
        call messages_print_stress(stdout, 'Lead Green functions')
        message(1) = ' st#  Spin  Lead     Energy'
        call write_info(1)
        do ik = 1, st%d%nik
          do ist = 1, st%ob_ncs
            energy = st%ob_eigenval(ist, ik)
            do il = 1, NLEADS
              do ispin = 1, hm%d%nspin
                select case(hm%d%ispin)
                case(UNPOLARIZED)
                  spin = '--'
                case(SPIN_POLARIZED)
                  if(is_spin_up(ik)) then
                    spin = 'up'
                  else
                    spin = 'dn'
                  end if
                  ! This is nonsene, but at least all indices are present.
                case(SPINORS)
                  if(ispin.eq.1) then
                    spin = 'up'
                  else
                    spin = 'dn'
                  end if
                end select
                write(message(1), '(i4,3x,a2,5x,a1,1x,f12.6)') ist, spin, ln(il), energy
                call write_info(1)
                call lead_green(energy, hm%lead_h_diag(:, :, ispin, il), hm%lead_h_offdiag(:, :, il), &
                  np, hm%lead_green(:, :, ispin, ist, ik, il), gr%sb%h(TRANS_DIR))

                ! Write the entire Green function to a file.
                if(in_debug_mode) then
                  call io_mkdir('debug/open_boundaries')
                  write(fname_real, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/green-', &
                    trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.real'
                  write(fname_imag, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/green-', &
                    trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.imag'
                  green_real = io_open(fname_real, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)
                  green_imag = io_open(fname_imag, action='write', grp=gr%mesh%mpi_grp, is_tmp=.false.)

                  write(fmt, '(a,i6,a)') '(', gr%intf(il)%np, 'e14.4)'
                  do irow = 1, gr%intf(il)%np
                    write(green_real, fmt) real(hm%lead_green(irow, :, ispin, ist, ik, il))
                    write(green_imag, fmt) aimag(hm%lead_green(irow, :, ispin, ist, ik, il))
                  end do
                  call io_close(green_real); call io_close(green_imag)
                end if
              end do
            end do
          end do
        end do
        call messages_print_stress(stdout)
      end if
    end subroutine init_lead_h
  end subroutine hamiltonian_init


  ! ---------------------------------------------------------
  subroutine hamiltonian_mg_init(hm, gr)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(inout) :: gr

    integer :: level

    call push_sub('epot.epot_mg_init')
    
    hm%multigrid_initialized = .true.

    call gridhier_init(hm%coarse_v, gr%mgrid, add_points_for_boundaries=.false.)

    hm%coarse_v%level(0)%p(1:NP) = hm%ep%vpsl(1:NP) + hm%vhxc(1:NP, 1)

    do level = 1, gr%mgrid%n_levels
      call dmultigrid_fine2coarse(gr%mgrid, level, hm%coarse_v%level(level - 1)%p, hm%coarse_v%level(level)%p, INJECTION)
    end do

  end subroutine hamiltonian_mg_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_end(hm, gr, geo)
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(inout) :: geo

    integer :: ispin

    call push_sub('hamiltonian.hamiltonian_end')

    do ispin = 1, hm%d%nspin
      call em_field_end(hm%total(ispin))
    end do

    deallocate(hm%total)

    if(hm%multigrid_initialized) then
      call gridhier_end(hm%coarse_v, gr%mgrid)
    end if

    if(associated(hm%phase)) deallocate(hm%phase)

    if(associated(hm%vhartree)) then
      deallocate(hm%vhartree)
      nullify(hm%vhartree)
    end if
    if(associated(hm%vhxc)) then
      deallocate(hm%vhxc)
      nullify(hm%vhxc)
    end if
    if(associated(hm%vxc)) then
      deallocate(hm%vxc)
      nullify(hm%vxc)
    end if
    if(associated(hm%axc)) then
      deallocate(hm%axc)
      nullify(hm%axc)
    end if
    if(associated(hm%a_ind)) then
      deallocate(hm%a_ind)
      nullify(hm%a_ind)
    end if
    if(associated(hm%b_ind)) then
      deallocate(hm%b_ind)
      nullify(hm%b_ind)
    end if

    if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) then
      deallocate(hm%vtau); nullify(hm%vtau)
    end if

    call epot_end(hm%ep, gr, geo)

    if(associated(hm%ab_pot)) then
      deallocate(hm%ab_pot); nullify(hm%ab_pot)
    end if

    DEALLOC(hm%lead_h_diag)
    DEALLOC(hm%lead_h_offdiag)
    DEALLOC(hm%lead_vks)
    DEALLOC(hm%lead_vhartree)
    DEALLOC(hm%lead_green)

    call states_dim_end(hm%d)
    call scissor_end(hm%scissor)
    call pop_sub()
  end subroutine hamiltonian_end


  ! ---------------------------------------------------------
  ! True if the Hamiltonian is Hermitean, false otherwise
  logical function hamiltonian_hermitean(hm)
    type(hamiltonian_t), intent(in) :: hm
    hamiltonian_hermitean = .not.((hm%ab .eq. IMAGINARY_ABSORBING) .or. hamiltonian_oct_exchange(hm))
  end function hamiltonian_hermitean

  ! ---------------------------------------------------------
  subroutine hamiltonian_span(hm, delta, emin)
    type(hamiltonian_t), intent(inout) :: hm
    FLOAT,               intent(in)    :: delta, emin

    call push_sub('hamiltonian.hamiltonian_span')

    hm%spectral_middle_point = ((M_Pi**2/(2*delta**2)) + emin)/M_TWO
    hm%spectral_half_span    = ((M_Pi**2/(2*delta**2)) - emin)/M_TWO

    call pop_sub()
  end subroutine hamiltonian_span


  ! ---------------------------------------------------------
  logical function hamiltonian_inh_term(hm) result(inh)
    type(hamiltonian_t), intent(in) :: hm
    inh = hm%inh_term
  end function hamiltonian_inh_term


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_inh(hm, st)
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t), target, intent(in) :: st
    hm%inh_st => st
    hm%inh_term = .true.
  end subroutine hamiltonian_set_inh


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_inh(hm)
    type(hamiltonian_t), intent(inout) :: hm
    nullify(hm%inh_st)
    hm%inh_term = .false.
  end subroutine hamiltonian_remove_inh


  ! ---------------------------------------------------------
  logical function hamiltonian_oct_exchange(hm) result(oct_exchange)
    type(hamiltonian_t), intent(in) :: hm
    oct_exchange = hm%oct_exchange
  end function hamiltonian_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_oct_exchange(hm)
    type(hamiltonian_t), intent(inout) :: hm
    ! In this release, no non-local part for the QOCT Hamiltonian.
    nullify(hm%oct_st)
    hm%oct_exchange = .false.
    !hm%oct_st => st
    !hm%oct_exchange = .true.
  end subroutine hamiltonian_set_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_oct_exchange(hm)
    type(hamiltonian_t), intent(inout) :: hm
    nullify(hm%oct_st)
    hm%oct_exchange = .false.
  end subroutine hamiltonian_remove_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm
    if(.not.hm%adjoint) then
      hm%adjoint = .true.
      if(hm%ab .eq. IMAGINARY_ABSORBING) then
        hm%ab_pot = -hm%ab_pot
      end if
    endif
  end subroutine hamiltonian_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_not_adjoint(hm)
    type(hamiltonian_t), intent(inout) :: hm
    if(hm%adjoint) then
      hm%adjoint = .false.
      if(hm%ab .eq. IMAGINARY_ABSORBING) then
        hm%ab_pot = -hm%ab_pot
      end if
    endif
  end subroutine hamiltonian_not_adjoint

  subroutine hamiltonian_update_potential(this, mesh)
    type(hamiltonian_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh

    integer :: ispin, ip

    do ispin = 1, this%d%nspin
      
      if(.not. associated(this%total(ispin)%potential)) then
        ALLOCATE(this%total(ispin)%potential(1:mesh%np), mesh%np)
      end if

      forall (ip = 1:mesh%np) this%total(ispin)%potential(ip) = this%vhxc(ip, ispin) + this%ep%vpsl(ip)
      
    end do

  end subroutine hamiltonian_update_potential

  subroutine hamiltonian_epot_generate(this, gr, geo, st, time)
    type(hamiltonian_t), intent(inout) :: this
    type(grid_t), target,  intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(states_t),        intent(inout) :: st
    FLOAT,       optional, intent(in)    :: time

    if(present(time)) then
      call epot_generate(this%ep, gr, geo, st, time)
    else
      call epot_generate(this%ep, gr, geo, st)
    end if

    call hamiltonian_update_potential(this, gr%mesh)

  end subroutine hamiltonian_epot_generate

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
