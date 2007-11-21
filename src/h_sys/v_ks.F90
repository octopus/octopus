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

module v_ks_m
  use datasets_m
  use global_m
  use grid_m
  use hamiltonian_m
  use lib_oct_parser_m
  use libxc_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use poisson_m
  use profiling_m
  use states_m
  use varinfo_m
  use xc_OEP_m
  use xc_m
  use magnetic_m
  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_calc

  integer, parameter :: &
    sic_none   = 1,     &  ! no self interaction correction
    sic_pz     = 2,     &  ! SIC a la Perdew Zunger (OEP way)
    sic_amaldi = 3         ! Amaldi correction term

  integer, parameter ::        &
    INDEPENDENT_PARTICLES = 1, &
    KOHN_SHAM_DFT         = 2, &
    HARTREE_FOCK          = 3

  type v_ks_t
    logical :: ip_app
    logical :: hartree_fock

    integer :: theory_level

    integer           :: xc_family  ! the xc stuff
    integer           :: sic_type   ! what kind of Self Interaction Correction to apply
    type(xc_t)        :: xc
    type(xc_OEP_t)    :: oep
  end type v_ks_t


contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(gr, ks, d)
    type(v_ks_t),        intent(out)   :: ks
    type(grid_t),        intent(inout) :: gr
    type(states_dim_t),  intent(in)    :: d

    call push_sub('v_ks.v_ks_init');

    !%Variable TheoryLevel
    !%Type integer
    !%Default dft
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with three different "theory levels": Kohn-Sham (TD)DFT
    !% (the usual one and the default), Hartree-Fock, or independent particles.
    !%Option independent_particles 1
    !%Option dft 2
    !%Option hartree_fock 3
    !%End
    call loct_parse_int(check_inp('TheoryLevel'), KOHN_SHAM_DFT, ks%theory_level)
    if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call input_error('TheoryLevel')

    call obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call obsolete_variable('HartreeFock', 'TheoryLevel')

    ks%ip_app = .false.
    ks%hartree_fock = .false.
    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%ip_app = .true.

    case(HARTREE_FOCK)
      ks%hartree_fock = .true.

      ! initilize xc modules
      call xc_init(ks%xc, NDIM, d%spin_channels, d%cdft, hartree_fock=.true.)
      ks%xc_family = ks%xc%family
      ks%sic_type = sic_none

      call v_ks_write_info(ks, stdout)

    case(KOHN_SHAM_DFT)
      ! initilize xc modules
      call xc_init(ks%xc, NDIM, d%spin_channels, d%cdft, hartree_fock=.false.)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA).ne.0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which Self Interaction Correction to use. Note that
        !% this correction will be applyed to the functional chosen by 'XFunctional' and
        !% 'CFunctional'
        !%Option sic_none 1
        !% No Self Interaction Correction
        !%Option sic_pz 2
        !% SIC a Perdew Zunger, hadled by the OEP technique
        !%Option sic_amaldi 3
        !% Amaldi correction term (NOT WORKING)
        !%End
        call loct_parse_int(check_inp('SICCorrection'), sic_none, ks%sic_type)
        if(.not.varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')

        ! Perdew Zunger corrections
        if(ks%sic_type == sic_pz) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
      else
        ks%sic_type = sic_none
      end if

      call xc_oep_init(ks%oep, ks%xc_family, gr, d)

      call v_ks_write_info(ks, stdout)

    end select

    call pop_sub()
  end subroutine v_ks_init


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks)
    type(v_ks_t), intent(inout) :: ks

    call push_sub('v_ks.v_ks_end');

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      call xc_oep_end(ks%oep)
      call xc_end(ks%xc)
    end select

    call pop_sub();
  end subroutine v_ks_end


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('v_ks.v_ks_write_info');

    call messages_print_stress(iunit, "Theory Level")
    call messages_print_var_option(iunit, "TheoryLevel", ks%theory_level)
    write(iunit, '(1x)')

    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      write(iunit, '(a)') 'Independent Particles'

    case(HARTREE_FOCK)
      call xc_write_info(ks%xc, iunit)

    case(KOHN_SHAM_DFT)
      call xc_write_info(ks%xc, iunit)

      write(iunit, '(1x)')
      call messages_print_var_option(iunit, 'SICCorrection', ks%sic_type)

      if(iand(ks%xc_family, XC_FAMILY_OEP).ne.0) then
        call xc_oep_write_info(ks%oep, iunit)
      end if
    end select

    call messages_print_stress(iunit)

    call pop_sub()
  end subroutine v_ks_write_info


  ! ---------------------------------------------------------
  subroutine v_ks_calc(gr, ks, h, st, calc_eigenval)
    type(grid_t),        intent(inout) :: gr
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(states_t),      intent(inout) :: st
    logical,      optional, intent(in) :: calc_eigenval

    FLOAT :: amaldi_factor

    call push_sub('v_ks.v_ks_calc')

    h%epot     = M_ZERO
    h%ehartree = M_ZERO

    ! check if we should introduce the amaldi SIC correction
    amaldi_factor = M_ONE
    if(ks%sic_type == sic_amaldi) amaldi_factor = (st%qtot-1)/st%qtot


    select case(ks%theory_level)
    case(KOHN_SHAM_DFT, HARTREE_FOCK)

      ! No Hartree or xc if independent electrons
      if(amaldi_factor>M_ZERO) then
        call v_hartree()

        !$omp parallel workshare
        h%vxc      = M_ZERO
        !$omp end parallel workshare
        if(h%d%cdft) h%axc = M_ZERO
        call v_a_xc()

        ! Build Hartree + xc potential

        !$omp parallel workshare
        h%vhxc(1:NP, 1) = h%vxc(1:NP, 1) + h%vhartree(1:NP)
        !$omp end parallel workshare

        if(h%d%ispin > UNPOLARIZED) then
          !$omp parallel workshare
          h%vhxc(1:NP, 2) = h%vxc(1:NP, 2) + h%vhartree(1:NP)
          !$omp end parallel workshare
        end if

        if(h%d%ispin == SPINORS) then
          !$omp parallel workshare
          h%vhxc(1:NP, 3:4) = h%vxc(1:NP, 3:4)
          !$omp end parallel workshare        
        end if
      end if

    case(INDEPENDENT_PARTICLES)

      !$omp parallel workshare
      h%vhxc     = M_ZERO
      !$omp end parallel workshare
      h%epot     = M_ZERO
      h%ex       = M_ZERO
      h%ec       = M_ZERO

    end select
    
    if(ks%theory_level == HARTREE_FOCK) then
      call states_end(h%st)
      call states_copy(h%st, st)

      h%exx_coef = ks%xc%exx_coef
    end if

    ! Calculate the potential vector induced by the electronic current
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    if(h%self_induced_magnetic) call magnetic_induced(gr, st, h%a_ind, h%b_ind)

    if(present(calc_eigenval)) then
      if (st%wfs_type == M_REAL) then
        call dhamiltonian_eigenval(h, gr, st)
      else
        call zhamiltonian_eigenval(h, gr, st)
      end if
    end if

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    ! Hartree contribution to the xc potential
    subroutine v_hartree()
      FLOAT, allocatable :: rho(:)
      integer :: is

      ALLOCATE(rho(NP), NP)

      ! calculate the total density
      !$omp parallel workshare
      rho(1:NP) = st%rho(1:NP, 1)
      !$omp end parallel workshare
      do is = 2, h%d%spin_channels
        !$omp parallel workshare
        rho(1:NP) = rho(1:NP) + st%rho(1:NP, is)
        !$omp end parallel workshare
      end do

      ! Amaldi correction
      if(ks%sic_type == sic_amaldi) rho = amaldi_factor*rho

      ! solve the poisson equation
      call dpoisson_solve(gr, h%vhartree, rho)

      ! Get the Hartree energy
      h%ehartree = M_HALF*dmf_dotp(gr%m, rho, h%vhartree)

      deallocate(rho)
    end subroutine v_hartree

    ! ---------------------------------------------------------
    subroutine v_a_xc()
      FLOAT, allocatable :: rho(:, :)
      integer :: is
      call profiling_in(C_PROFILING_XC)

      h%ex = M_ZERO
      h%ec = M_ZERO
      h%exc_j = M_ZERO

      ! get density taking into account non-linear core corrections, and the Amaldi SIC correction
      ALLOCATE(rho(NP, st%d%nspin), NP*st%d%nspin)
      if(associated(st%rho_core)) then
        do is = 1, st%d%spin_channels
          rho(1:NP, is) = st%rho(1:NP, is) + st%rho_core(1:NP)/st%d%spin_channels
        end do
      else
        rho(1:NP, :) = st%rho(1:NP, :)
      end if

      ! Amaldi correction
      if(ks%sic_type == sic_amaldi) rho(1:NP,:) = amaldi_factor*rho(1:NP,:)

      ! Get the *local* xc term
      if(h%d%cdft) then
        call xc_get_vxc_and_axc(gr, ks%xc, rho, st%j, st%d%ispin, h%vxc, h%axc, &
             h%ex, h%ec, h%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
      else
        call xc_get_vxc(gr, ks%xc, rho, st%d%ispin, h%vxc, h%ex, h%ec, &
             -minval(st%eigenval(st%nst, :)), st%qtot)
      end if
      deallocate(rho)

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to handle specially
        if (st%wfs_type == M_REAL) then
          call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
            gr, h, st, h%vxc, h%ex, h%ec)
        else
          call zxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
            gr, h, st, h%vxc, h%ex, h%ec)
        end if
      end if

      ! Now we calculate Int[n vxc] = h%epot
      select case(h%d%ispin)
      case(UNPOLARIZED)
        h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1))
      case(SPIN_POLARIZED)
        h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1)) &
             + dmf_dotp(gr%m, st%rho(:, 2), h%vxc(:, 2))
      case(SPINORS)
        h%epot = h%epot + dmf_dotp(gr%m, st%rho(:, 1), h%vxc(:, 1)) &
             + dmf_dotp(gr%m, st%rho(:, 2), h%vxc(:, 2)) &
             + M_TWO*dmf_dotp(gr%m, st%rho(:, 3), h%vxc(:, 3)) &
             + M_TWO*dmf_dotp(gr%m, st%rho(:, 4), h%vxc(:, 4))

      end select

      call profiling_out(C_PROFILING_XC)
    end subroutine v_a_xc
  end subroutine v_ks_calc

end module v_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
