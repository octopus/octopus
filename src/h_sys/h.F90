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
  use datasets_m
  use derivatives_m
  use functions_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use lalg_basic_m
  use loct_parser_m
  use mesh_m
  use external_pot_m
  use messages_m
  use mpi_m
  use multigrid_m
  use profiling_m
  use projector_m
  use simul_box_m
  use states_m
  use units_m
  use varinfo_m
  use lasers_m
  use poisson_m

  implicit none

  private
  public ::                &
    hamiltonian_t,         &
    hamiltonian_init,      &
    hamiltonian_mg_init,   &
    hamiltonian_end,       &
    hamiltonian_energy,    &
    hamiltonian_span,      &
    dhamiltonian_eigenval, &
    zhamiltonian_eigenval, &
    delectronic_kinetic_energy,  &
    zelectronic_kinetic_energy,  &
    delectronic_external_energy, &
    zelectronic_external_energy, &
    dhpsi,                 &
    dhpsi_diag,            &
    dvlpsi,                &
    dvnlpsi,               &
    dvnlpsi_diag,          &
    dmagnus,               &
    dkinetic,              &
    dvmask,                &
    zhpsi,                 &
    zhpsi_diag,            &
    zvlpsi,                &
    zvnlpsi,               &
    zvnlpsi_diag,          &
    zmagnus,               &
    zkinetic,              &
    zvmask,                &
    zvlasers,              &
    zvlaser_operator_linear,         &
    zvlaser_operator_quadratic,      &
    hamiltonian_inh_term,            &
    hamiltonian_set_inh,             &
    hamiltonian_remove_inh,          &
    hamiltonian_oct_exchange,        &
    hamiltonian_set_oct_exchange,    &
    hamiltonian_remove_oct_exchange, &
    hamiltonian_adjoint,             &
    hamiltonian_not_adjoint,         &
    hamiltonian_hermitean



  type hamiltonian_t
    ! The Hamiltonian must know what are the "dimensions" of the spaces,
    ! in order to be able to operate on the states.
    type(states_dim_t) :: d

    FLOAT, pointer :: vhartree(:) ! hartree potential
    FLOAT, pointer :: vxc(:,:)    ! xc potential
    FLOAT, pointer :: vhxc(:,:)   ! xc potential + hartree potential
    FLOAT, pointer :: axc(:,:,:)  ! xc vector-potential divided by c

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
             eext        ! External     V = <Phi|V|Phi> = Int[n v] (if no non-local pseudos exist)

    integer :: theory_level    ! copied from sys%ks

    type(epot_t) :: ep         ! handles the external potential

    ! gauge
    integer :: gauge              ! in which gauge shall we work in

    ! absorbing boundaries
    logical :: adjoint
    integer  :: ab                ! do we have absorbing boundaries?
    FLOAT :: ab_width             ! width of the absorbing boundary
    FLOAT :: ab_height            ! height of the absorbing boundary
    FLOAT, pointer :: ab_pot(:)   ! where we store the ab potential

    ! Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span

    ! Kinetic Cutoff
    FLOAT :: cutoff

    ! Mass of the particle (in most cases, mass = 1, electron mass)
    FLOAT :: mass

    ! For the Hartree Fock Hamiltonian, the Fock operator depends on the states.
    type(states_t) :: st

    ! There may be an "inhomogeneous", "source", or "forcing" term (useful for the OCT formalism)
    logical :: inh_term
    type(states_t), pointer :: inh_st

    ! There may also be a exchange-like term, similar to the one necessary for time-dependent
    ! Hartree Fock, also useful only for the OCT equations
    logical :: oct_exchange
    type(states_t), pointer :: oct_st

    ! Handles for the calculation of the kinetic energy in parallel.
    type(der_handle_t), pointer :: handles(:)

    CMPLX, pointer :: phase(:, :)

    type(mg_float_pointer), pointer :: coarse_v(:)

  end type hamiltonian_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter :: &
    NO_ABSORBING        = 0,    &
    IMAGINARY_ABSORBING = 1,    &
    MASK_ABSORBING      = 2

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 1, &
    HARTREE               = 2, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_init(h, gr, geo, states_dim, wfs_type, theory_level)
    type(hamiltonian_t), intent(out)   :: h
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_dim_t),  intent(in)    :: states_dim
    integer,             intent(inout) :: wfs_type
    integer,             intent(in)    :: theory_level

    integer :: i, j, n, ispin
    FLOAT :: d(MAX_DIM)

    call push_sub('h.hamiltonian_init')

    ! make a couple of local copies
    h%theory_level = theory_level
    call states_dim_copy(h%d, states_dim)

    ! Allocate handles
    ALLOCATE(h%handles(h%d%dim), h%d%dim)
    do i = 1, h%d%dim
      call der_handle_init(h%handles(i), gr%m)
    end do

    ! initialize variables
    h%epot = M_ZERO
    h%ex = M_ZERO; h%ec = M_ZERO
    h%etot = M_ZERO

    ! allocate potentials and density of the cores
    ! In the case of spinors, vxc_11 = h%vxc(:, 1), vxc_22 = h%vxc(:, 2), Re(vxc_12) = h%vxc(:. 3);
    ! Im(vxc_12) = h%vxc(:, 4)
    ALLOCATE(h%vhartree(NP),        NP)
    ALLOCATE(h%vxc(NP, h%d%nspin), NP*h%d%nspin)
    ALLOCATE(h%vhxc(NP, h%d%nspin), NP*h%d%nspin)

    !$omp parallel workshare
    h%vhartree(1:NP) = M_ZERO
    !$omp end parallel workshare

    do ispin = 1, h%d%nspin
      !$omp parallel workshare
      h%vhxc(1:NP, ispin) = M_ZERO
      h%vxc(1:NP, ispin) = M_ZERO
      !$omp end parallel workshare
    end do


    if (h%d%cdft) then
      ALLOCATE(h%axc(NP, NDIM, h%d%nspin), NP*NDIM*h%d%nspin)
      h%axc = M_ZERO
    else
      nullify(h%axc)
    end if

    !Initialize external potential
    call epot_init(h%ep, gr, geo, h%d%ispin)


    !Static magnetic field requires complex wave-functions
    if (associated(h%ep%B_field) .or. gauge_field_is_applied(h%ep%gfield)) wfs_type = M_CMPLX

    call loct_parse_logical(check_inp('CalculateSelfInducedMagneticField'), .false., h%self_induced_magnetic)
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
    if(h%self_induced_magnetic) then
      select case(NDIM)
      case(3)
        ALLOCATE(h%a_ind(NP_PART, MAX_DIM), NP_PART*MAX_DIM)
        ALLOCATE(h%b_ind(NP_PART, MAX_DIM), NP_PART*MAX_DIM)
      case(2)
        ALLOCATE(h%a_ind(NP_PART, 2), NP_PART*2)
        ALLOCATE(h%b_ind(NP_PART, 1), NP_PART)
      end select
    else
      nullify(h%a_ind, h%b_ind)
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
    call loct_parse_int(check_inp('TDGauge'), LENGTH, h%gauge)
    if(.not.varinfo_valid_option('TDGauge', h%gauge)) call input_error('TDGauge')
    call messages_print_var_option(stdout, "TDGauge", h%gauge)

    !%Variable AbsorbingBoundaries
    !%Type integer
    !%Default no_absorbing
    !%Section Time Dependent::Absorbing Boundaries
    !%Description
    !% To improve the quality of the spectra by avoiding the formation of 
    !% standing density waves, one can make the boundaries of the simulation 
    !% box absorbing.
    !%Option no_absorbing 0
    !% No absorbing boundaries.
    !%Option sin2 1
    !% A <math>\sin^2</math> imaginary potential is added at the boundaries.
    !%Option mask 2
    !% A mask is applied to the wave-functions at the boundaries.
    !%End
    call loct_parse_int(check_inp('AbsorbingBoundaries'), NO_ABSORBING, h%ab)
    if(.not.varinfo_valid_option('AbsorbingBoundaries', h%ab)) call input_error('AbsorbingBoundaries')
    call messages_print_var_option(stdout, "AbsorbingBoundaries", h%ab)
    
    nullify(h%ab_pot)

    absorbing_boundaries: if(h%ab.ne.NO_ABSORBING) then
      !%Variable ABWidth
      !%Type float
      !%Default 0.4 a.u.
      !%Section Time Dependent::Absorbing Boundaries
      !%Description
      !% Width of the region used to apply the absorbing boundaries.
      !%End
      call loct_parse_float(check_inp('ABWidth'), CNST(0.4)/units_inp%length%factor, h%ab_width)
      h%ab_width  = h%ab_width * units_inp%length%factor
      if(h%ab == 1) then
        !%Variable ABHeight
        !%Type float
        !%Default -0.2 a.u.
        !%Section Time Dependent::Absorbing Boundaries
        !%Description 
        !% When <tt>AbsorbingBoundaries == sin2</tt>, is the height of the imaginary potential.
        !%End
        call loct_parse_float(check_inp('ABHeight'), -CNST(0.2)/units_inp%energy%factor, h%ab_height)
        h%ab_height = h%ab_height * units_inp%energy%factor
      else
        h%ab_height = M_ONE
      end if

      ! generate boundary potential...
      ALLOCATE(h%ab_pot(NP), NP)
      h%ab_pot = M_ZERO
      do i = 1, NP
        call mesh_inborder(gr%m, i, n, d, h%ab_width)
        if(n>0) then
          do j = 1, n
            h%ab_pot(i) = h%ab_pot(i) + h%ab_height * sin(d(j)*M_PI/(M_TWO*h%ab_width))**2
          end do
        end if
        if(abs(h%ab_pot(i)) > abs(h%ab_height)) h%ab_pot(i) = h%ab_height
      end do

    end if absorbing_boundaries

    ! Cutoff applied to the kinetic term.
    ! it is used *both* in the calculation of the derivatives and in the split operator method
    call loct_parse_float(check_inp('KineticCutoff'), -M_ONE, h%cutoff)
    if(h%cutoff > M_ZERO) then
      h%cutoff = h%cutoff * units_inp%energy%factor
      write(message(1),'(a,f7.2,a)') 'Info: The kinetic operator will have a cutoff of',&
        h%cutoff/units_out%energy%factor, units_out%energy%abbrev
      write(message(2),'(a)')        '      (only if DerivativesSpace = 1 is set)'
      call write_info(2)
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
    call loct_parse_float(check_inp('ParticleMass'), M_ONE, h%mass)

    call states_null(h%st)

    call hamiltonian_remove_inh(h)
    call hamiltonian_remove_oct_exchange(h)

    h%adjoint = .false.

    nullify(h%phase)
    if (simul_box_is_periodic(gr%sb)) call init_phase

    call pop_sub()

  contains
    
    subroutine init_phase
      integer :: ip, ik

      ALLOCATE(h%phase(1:NP_PART, 1:h%d%nik), NP_PART*h%d%nik)

      do ik = 1, h%d%nik
        !$omp parallel do
        do ip = 1, NP_PART
          h%phase(ip, ik) = exp(-M_zI*sum(gr%m%x(ip, 1:MAX_DIM)* h%d%kpoints(1:MAX_DIM, ik)))
        end do
        !$omp end parallel do
      end do
      
    end subroutine init_phase
    
  end subroutine hamiltonian_init

  subroutine hamiltonian_mg_init(h, gr)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr

    integer :: level

    call push_sub('epot.epot_mg_init')

    call gridhier_init(h%coarse_v, gr%mgrid, add_points_for_boundaries=.false.)

    h%coarse_v(0)%p(1:NP) = h%ep%vpsl(1:NP) + h%vhxc(1:NP, 1)

    do level = 1, gr%mgrid%n_levels
      call multigrid_fine2coarse(gr%mgrid, level, h%coarse_v(level - 1)%p, h%coarse_v(level)%p, INJECTION)
    end do

  end subroutine hamiltonian_mg_init

  ! ---------------------------------------------------------
  subroutine hamiltonian_end(h, gr, geo)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(inout) :: geo

    integer :: ii

    call push_sub('h.hamiltonian_end')

    if(associated(h%coarse_v)) then
      call gridhier_end(h%coarse_v, gr%mgrid)
    end if

    if(associated(h%phase)) deallocate(h%phase)

    if(associated(h%handles)) then
      do ii = 1, h%d%dim
        call der_handle_end(h%handles(ii))
      end do
      deallocate(h%handles)
      nullify(h%handles)
    end if
    if(associated(h%vhartree)) then
      deallocate(h%vhartree)
      nullify(h%vhartree)
    end if
    if(associated(h%vhxc)) then
      deallocate(h%vhxc)
      nullify(h%vhxc)
    end if
    if(associated(h%vxc)) then
      deallocate(h%vxc)
      nullify(h%vxc)
    end if
    if(associated(h%axc)) then
      deallocate(h%axc)
      nullify(h%axc)
    end if
    if(associated(h%a_ind)) then
      deallocate(h%a_ind)
      nullify(h%a_ind)
    end if
    if(associated(h%b_ind)) then
      deallocate(h%b_ind)
      nullify(h%b_ind)
    end if

    call epot_end(h%ep, gr, geo)

    if(associated(h%ab_pot)) then
      deallocate(h%ab_pot); nullify(h%ab_pot)
    end if

    call states_dim_end(h%d)

    call pop_sub()
  end subroutine hamiltonian_end


  ! ---------------------------------------------------------
  ! True if the Hamiltonian is Hermitean, false otherwise
  logical function hamiltonian_hermitean(h)
    type(hamiltonian_t), intent(in) :: h
    hamiltonian_hermitean = .not.((h%ab .eq. IMAGINARY_ABSORBING) .or. hamiltonian_oct_exchange(h))
  end function hamiltonian_hermitean


  ! ---------------------------------------------------------
  ! This subroutine calculates the total energy of the system. Basically, it
  ! adds up the KS eigenvalues, and then it substracts the whatever double
  ! counts exist (see TDDFT theory for details).
  subroutine hamiltonian_energy(h, gr, geo, st, iunit, full)
    type(hamiltonian_t), intent(inout) :: h
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: iunit
    logical, optional,   intent(in)    :: full

    logical :: full_

    call push_sub('h.hamiltonian_energy')

    full_ = .false.
    if(present(full)) full_ = full

    select case(h%theory_level)
    case(INDEPENDENT_PARTICLES)
      h%eeigen = states_eigenvalues_sum(st)
      h%etot   = h%ep%eii + h%eeigen

    case(HARTREE)
      if(st%wfs_type == M_REAL) then
        h%t0     = delectronic_kinetic_energy(h, gr, st)
        h%eext   = delectronic_external_energy(h, gr, st)
      else
        h%t0     = zelectronic_kinetic_energy(h, gr, st)
        h%eext   = zelectronic_external_energy(h, gr, st)
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot = h%ep%eii + M_HALF*(h%eeigen + h%t0 + h%eext)

    case(HARTREE_FOCK)
      if(st%wfs_type == M_REAL) then
        h%t0     = delectronic_kinetic_energy(h, gr, st)
        h%eext   = delectronic_external_energy(h, gr, st)
      else
        h%t0     = zelectronic_kinetic_energy(h, gr, st)
        h%eext   = zelectronic_external_energy(h, gr, st)
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot = h%ep%eii + M_HALF*(h%eeigen + h%t0 + h%eext - h%epot) + h%ec

    case(KOHN_SHAM_DFT)
      if(full_) then
        if(st%wfs_type == M_REAL) then
          h%t0     = delectronic_kinetic_energy(h, gr, st)
          h%eext   = delectronic_external_energy(h, gr, st)
        else
          h%t0     = zelectronic_kinetic_energy(h, gr, st)
          h%eext   = zelectronic_external_energy(h, gr, st)
        end if
      end if
      h%eeigen = states_eigenvalues_sum(st)
      h%etot   = h%ep%eii + h%eeigen - h%ehartree + h%ex + h%ec - h%epot

    end select
    
    if(gauge_field_is_applied(h%ep%gfield)) then
      h%etot = h%etot + gauge_field_get_energy(h%ep%gfield, gr%sb)
    end if

    if (iunit > 0) then
      write(message(1), '(6x,a, f18.8)')'Total       = ', h%etot     / units_out%energy%factor
      write(message(2), '(6x,a, f18.8)')'Ion-ion     = ', h%ep%eii   / units_out%energy%factor
      write(message(3), '(6x,a, f18.8)')'Eigenvalues = ', h%eeigen   / units_out%energy%factor
      write(message(4), '(6x,a, f18.8)')'Hartree     = ', h%ehartree / units_out%energy%factor
      write(message(5), '(6x,a, f18.8)')'Int[n*v_xc] = ', h%epot     / units_out%energy%factor
      write(message(6), '(6x,a, f18.8)')'Exchange    = ', h%ex       / units_out%energy%factor
      write(message(7), '(6x,a, f18.8)')'Correlation = ', h%ec       / units_out%energy%factor
      call write_info(7, iunit)
      if(full_) then
        write(message(1), '(6x,a, f18.8)')'Kinetic     = ', h%t0 / units_out%energy%factor
        write(message(2), '(6x,a, f18.8)')'External    = ', h%eext / units_out%energy%factor
        call write_info(2, iunit)
      end if
    end if

    call pop_sub()
  end subroutine hamiltonian_energy


  ! ---------------------------------------------------------
  subroutine hamiltonian_span(h, delta, emin)
    type(hamiltonian_t), intent(inout) :: h
    FLOAT,               intent(in)    :: delta, emin

    call push_sub('h.hamiltonian_span')

    h%spectral_middle_point = ((M_Pi**2/(2*delta**2)) + emin)/M_TWO
    h%spectral_half_span    = ((M_Pi**2/(2*delta**2)) - emin)/M_TWO

    call pop_sub()
  end subroutine hamiltonian_span


  ! ---------------------------------------------------------
  logical function hamiltonian_inh_term(h) result(inh)
    type(hamiltonian_t), intent(in) :: h
    inh = h%inh_term
  end function hamiltonian_inh_term


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_inh(h, st)
    type(hamiltonian_t), intent(inout) :: h
    type(states_t), target, intent(in) :: st
    h%inh_st => st
    h%inh_term = .true.
  end subroutine hamiltonian_set_inh


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_inh(h)
    type(hamiltonian_t), intent(inout) :: h
    nullify(h%inh_st)
    h%inh_term = .false.
  end subroutine hamiltonian_remove_inh


  ! ---------------------------------------------------------
  logical function hamiltonian_oct_exchange(h) result(oct_exchange)
    type(hamiltonian_t), intent(in) :: h
    oct_exchange = h%oct_exchange
  end function hamiltonian_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_set_oct_exchange(h, st)
    type(hamiltonian_t), intent(inout) :: h
    type(states_t),      intent(in)    :: st
    ! In this release, no non-local part for the QOCT Hamiltonian.
    nullify(h%oct_st)
    h%oct_exchange = .false.
    !h%oct_st => st
    !h%oct_exchange = .true.
  end subroutine hamiltonian_set_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_remove_oct_exchange(h)
    type(hamiltonian_t), intent(inout) :: h
    nullify(h%oct_st)
    h%oct_exchange = .false.
  end subroutine hamiltonian_remove_oct_exchange


  ! ---------------------------------------------------------
  subroutine hamiltonian_adjoint(h)
    type(hamiltonian_t), intent(inout) :: h
    if(.not.h%adjoint) then
      h%adjoint = .true.
      if(h%ab .eq. IMAGINARY_ABSORBING) then
        h%ab_pot = -h%ab_pot
      end if
    endif
  end subroutine hamiltonian_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_not_adjoint(h)
    type(hamiltonian_t), intent(inout) :: h
    if(h%adjoint) then
      h%adjoint = .false.
      if(h%ab .eq. IMAGINARY_ABSORBING) then
        h%ab_pot = -h%ab_pot
      end if
    endif
  end subroutine hamiltonian_not_adjoint


#include "undef.F90"
#include "real.F90"
#include "h_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "h_inc.F90"

end module hamiltonian_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
