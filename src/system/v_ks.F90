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
  use energy_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use lalg_basic_m
  use parser_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use poisson_m
  use poisson_sete_m
  use profiling_m
  use states_m
  use states_dim_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)
  use xc_OEP_m

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_calc,          &
    v_ks_hartree,       &
    v_ks_freeze_hxc


  integer, parameter, public :: &
    sic_none   = 1,     &  ! no self-interaction correction
    sic_pz     = 2,     &  ! Perdew-Zunger SIC (OEP way)
    sic_amaldi = 3         ! Amaldi correction term

  type v_ks_t
    integer :: theory_level

    logical :: frozen_hxc ! For RPA and SAE calculations.

    integer           :: xc_family  ! the XC stuff
    integer           :: sic_type   ! what kind of self-interaction correction to apply
    type(xc_t)        :: xc
    type(xc_OEP_t)    :: oep
  end type v_ks_t


contains

  
  ! ---------------------------------------------------------
  subroutine v_ks_init(gr, ks, d, nel)
    type(v_ks_t),        intent(out)   :: ks
    type(grid_t),        intent(inout) :: gr
    type(states_dim_t),  intent(in)    :: d
    FLOAT,               intent(in)    :: nel ! the total number of electrons

    call push_sub('v_ks.v_ks_init')

    !%Variable TheoryLevel
    !%Type integer
    !%Default dft
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with three different "theory levels":
    !%Option independent_particles 2
    !% Particles will be considered as independent, <i>i.e.</i> as non-interacting.
    !% This mode is mainly used for testing purposes, as the code is usually 
    !% much faster with <tt>independent_particles</tt>.
    !%Option hartree 1
    !% Calculation within the Hartree method. Note that, contrary to popular
    !% belief, the Hartree potential is self-interaction-free. Therefore, this run 
    !% mode will not yield the same result as <tt>dft</tt> without exchange-correlation.
    !% Development version only.
    !%Option hartree_fock 3
    !% This is the traditional Hartree-Fock scheme. Like the Hartree scheme, it is fully
    !% self-interaction-free. This mode is extremely slow. It is often more convenient
    !% to use <tt>dft</tt> within the OEP scheme to get similar (but not the same) results.
    !% Note that within this scheme you can use a correlation functional, or a hybrid
    !% functional (see <tt>XCFunctional</tt>). In the latter case, you will be following the
    !% quantum-chemistry recipe to use hybrids.
    !%Option dft 4
    !% This is the default density-functional theory scheme. Note that you can also use 
    !% hybrids in this scheme, but they will be handled the "DFT" way, <i>i.e.</i>, solving the
    !% OEP equation.
    !%End
    call parse_integer(datasets_check('TheoryLevel'), KOHN_SHAM_DFT, ks%theory_level)
    if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call input_error('TheoryLevel')

    call messages_obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable('HartreeFock', 'TheoryLevel')

    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%sic_type = sic_none
    case(HARTREE)
      call messages_devel_version("Hartree theory level")
    case(HARTREE_FOCK)
      ! initilize xc modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, d%spin_channels, d%cdft, hartree_fock=.true.)
      ks%xc_family = ks%xc%family
      ks%sic_type = sic_none

    case(KOHN_SHAM_DFT)
      ! initilize xc modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, d%spin_channels, d%cdft, hartree_fock=.false.)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA).ne.0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which form of self-interaction correction to use. Note that
        !% this correction will be applied to the functional chosen by <tt>XCFunctional</tt>.
        !%Option sic_none 1
        !% No self-interaction correction.
        !%Option sic_pz 2
        !% Perdew-Zunger SIC, handled by the OEP technique.
        !%Option sic_amaldi 3
        !% Amaldi correction term (NOT WORKING).
        !%End
        call parse_integer(datasets_check('SICCorrection'), sic_none, ks%sic_type)
        if(.not.varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')

        ! Perdew-Zunger corrections
        if(ks%sic_type == sic_pz) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
      else
        ks%sic_type = sic_none
      end if

      call xc_oep_init(ks%oep, ks%xc_family, gr, d)
    end select

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, stdout)

    call pop_sub()
  end subroutine v_ks_init
  ! ---------------------------------------------------------


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


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('v_ks.v_ks_write_info')

    call messages_print_stress(iunit, "Theory Level")
    call messages_print_var_option(iunit, "TheoryLevel", ks%theory_level)

    select case(ks%theory_level)
    case(HARTREE_FOCK)
      write(iunit, '(1x)')
      call xc_write_info(ks%xc, iunit)

    case(KOHN_SHAM_DFT)
      write(iunit, '(1x)')
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


  ! ---------------------------------------------------------
  subroutine v_ks_calc(gr, ks, hm, st, calc_eigenval)
    type(grid_t),           intent(inout) :: gr
    type(v_ks_t),           intent(inout) :: ks
    type(hamiltonian_t),    intent(inout) :: hm
    type(states_t),         intent(inout) :: st
    logical,      optional, intent(in)    :: calc_eigenval

    FLOAT :: amaldi_factor
    integer :: ip, ispin
    type(profile_t), save :: prof
    logical :: calc_eigenval_

    ! The next line is a hack to be able to perform an IP/RPA calculation
    !logical, save :: RPA_first = .true.

    call push_sub('v_ks.v_ks_calc')
    call profiling_in(prof, "KOHN_SHAM_CALC")

    if(in_debug_mode) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call write_info(1)
    end if

    calc_eigenval_ = .false.
    if(present(calc_eigenval)) calc_eigenval_ = calc_eigenval

    ! If the Hxc term is frozen, there is nothing to do, except we 
    ! maybe have to calculate the eigenvalues (and WARNING: MISSING
    ! hm%epot)
    if(ks%frozen_hxc) then
      if(calc_eigenval_) then
        if (states_are_real(st)) then
          call dcalculate_eigenvalues(hm, gr, st)
        else
          call zcalculate_eigenvalues(hm, gr, st)
        end if
      end if
      call pop_sub(); return
    end if

    hm%epot     = M_ZERO

    ! check if we should introduce the Amaldi SIC correction
    amaldi_factor = M_ONE
    if(ks%sic_type == sic_amaldi) amaldi_factor = (st%qtot-1)/st%qtot


    if(ks%theory_level==INDEPENDENT_PARTICLES .or. amaldi_factor==M_ZERO) then
      hm%vhxc     = M_ZERO
      hm%epot     = M_ZERO
      hm%ehartree = M_ZERO
      hm%ex       = M_ZERO
      hm%ec       = M_ZERO
    else
      ! The next 2 lines are a hack to be able to perform an IP/RPA calculation
      !if(RPA_first) then
      !  RPA_first = .false.

        hm%vxc      = M_ZERO
        if(iand(hm%xc_family, XC_FAMILY_MGGA).ne.0) hm%vtau = M_ZERO
        if(hm%d%cdft) hm%axc = M_ZERO
        if(ks%theory_level.ne.HARTREE) call v_a_xc()

        hm%ehartree = M_ZERO
        call v_ks_hartree(gr, st, hm, amaldi_factor)

      !end if

      ! Build Hartree + xc potential

       forall(ip = 1:gr%mesh%np) hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)

      if(hm%d%ispin > UNPOLARIZED) then
        forall(ip = 1:gr%mesh%np) hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
      end if

      if(hm%d%ispin == SPINORS) then
        forall(ip = 1:gr%mesh%np, ispin = 3:4) hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
      end if

    end if
    
    call hamiltonian_update_potential(hm, gr%mesh)
    
    if(ks%theory_level==HARTREE.or.ks%theory_level==HARTREE_FOCK) then
      call states_end(hm%st)
      call states_copy(hm%st, st)
    end if
    if(ks%theory_level==HARTREE_FOCK) then
      hm%exx_coef = ks%xc%exx_coef
    else if (ks%theory_level==HARTREE) then
      hm%exx_coef = M_ONE
    end if

    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    if(hm%self_induced_magnetic) call magnetic_induced(gr, st, hm%a_ind, hm%b_ind)

    if(calc_eigenval_) then
      if (states_are_real(st)) then
        call dcalculate_eigenvalues(hm, gr, st)
      else
        call zcalculate_eigenvalues(hm, gr, st)
      end if
    end if

    call profiling_out(prof)
    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine v_a_xc()
      FLOAT, allocatable :: rho(:, :)
      type(profile_t), save :: prof
      call push_sub('v_ks.v_ks_calc.v_a_xc')
      call profiling_in(prof, "XC")

      hm%ex = M_ZERO
      hm%ec = M_ZERO
      hm%exc_j = M_ZERO

      ! get density taking into account non-linear core corrections, and the Amaldi SIC correction
      SAFE_ALLOCATE(rho(1:gr%fine%mesh%np, 1:st%d%nspin))
      call states_total_density(st, gr%fine%mesh, rho)

      ! Amaldi correction
      if(ks%sic_type == sic_amaldi) then
        rho(1:gr%fine%mesh%np, :) = amaldi_factor*rho(1:gr%fine%mesh%np, :)
      end if

      ! Get the *local* xc term
      if(hm%d%cdft) then
        call xc_get_vxc_and_axc(gr, ks%xc, st, rho, st%current, st%d%ispin, hm%vxc, hm%axc, &
             hm%ex, hm%ec, hm%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
      else
        call xc_get_vxc(gr, ks%xc, st, rho, st%d%ispin, hm%ex, hm%ec, &
             -minval(st%eigenval(st%nst, :)), st%qtot, vxc=hm%vxc, vtau=hm%vtau)
      end if
      SAFE_DEALLOCATE_A(rho)

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if (states_are_real(st)) then
          call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
            gr, hm, st, hm%ex, hm%ec, vxc=hm%vxc)
        else
          call zxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
            gr, hm, st, hm%ex, hm%ec, vxc=hm%vxc)
        end if
      end if

      ! Now we calculate Int[n vxc] = hm%epot
      select case(hm%d%ispin)
      case(UNPOLARIZED)
        hm%epot = hm%epot + dmf_dotp(gr%mesh, st%rho(:, 1), hm%vxc(:, 1))
      case(SPIN_POLARIZED)
        hm%epot = hm%epot + dmf_dotp(gr%mesh, st%rho(:, 1), hm%vxc(:, 1)) &
             + dmf_dotp(gr%mesh, st%rho(:, 2), hm%vxc(:, 2))
      case(SPINORS)
        hm%epot = hm%epot + dmf_dotp(gr%mesh, st%rho(:, 1), hm%vxc(:, 1)) &
             + dmf_dotp(gr%mesh, st%rho(:, 2), hm%vxc(:, 2)) &
             + M_TWO*dmf_dotp(gr%mesh, st%rho(:, 3), hm%vxc(:, 3)) &
             + M_TWO*dmf_dotp(gr%mesh, st%rho(:, 4), hm%vxc(:, 4))

      end select

      call profiling_out(prof)
      call pop_sub()
    end subroutine v_a_xc
  end subroutine v_ks_calc
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! Hartree contribution to the xc potential
  subroutine v_ks_hartree(gr, st, hm, amaldi_factor)
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(in)    :: st
    FLOAT, optional,     intent(in)    :: amaldi_factor

    FLOAT, allocatable :: rho(:)
    integer :: is, ip

    call push_sub('v_ks.v_ks_hartree')

    SAFE_ALLOCATE(rho(1:gr%mesh%np))

    ! calculate the total density
    call lalg_copy(gr%mesh%np, st%rho(:, 1), rho)
    do is = 2, hm%d%spin_channels
      forall(ip = 1:gr%mesh%np) rho(ip) = rho(ip) + st%rho(ip, is)
    end do

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      do is = 1, hm%d%spin_channels
        forall(ip = 1:gr%mesh%np) rho(ip) = rho(ip) + st%frozen_rho(ip, is)
      end do
    end if

    ! Amaldi correction
    if(present(amaldi_factor)) rho = amaldi_factor*rho

    if(poisson_solver_is_iterative()) then
      ! provide a better starting point (in the td case vhxc was interpolated)
      forall(ip = 1:gr%mesh%np) hm%vhartree(ip) = hm%vhxc(ip, 1) - hm%vxc(ip, 1)
    end if

    ! solve the Poisson equation
    call dpoisson_solve(gr, hm%vhartree, rho)
    ! Get the Hartree energy
    hm%ehartree = M_HALF*dmf_dotp(gr%mesh, rho, hm%vhartree)

    if (poisson_get_solver() == POISSON_SETE) then !Roberto
      ! How to get the nuclear density here?
      hm%ep%eii = M_HALF*dmf_dotp(gr%mesh, rho_nuc, hm%ep%vpsl) 
      !hm%ep%eii = 0
      ! ESURF is a module called from later.  Need to write it into somewhere...
      !hm%etot=-poisson_sete_energy(sete_solver)
!      write(525,*) hm%ehartree+poisson_sete_energy(sete_solver),hm%ehartree,poisson_sete_energy(sete_solver)
      hm%ehartree=hm%ehartree+poisson_sete_energy(sete_solver)

!      write(89,*) hm%ehartree*CNST(27.2), poisson_sete_energy(sete_solver)*CNST(27.2), hm%ep%eii*CNST(27.2)
    endif

    SAFE_DEALLOCATE_A(rho)
    call pop_sub()
  end subroutine v_ks_hartree
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_freeze_hxc(ks)
    type(v_ks_t), intent(inout) :: ks

    call push_sub('v_ks.v_ks_freeze_hxc')

    ks%frozen_hxc = .true.
    
    call pop_sub()
  end subroutine v_ks_freeze_hxc
  ! ---------------------------------------------------------



end module v_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
