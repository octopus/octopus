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
  use berry_m
  use datasets_m
  use density_m
  use derivatives_m
  use energy_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use lalg_basic_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use multigrid_m
  use parser_m
  use poisson_m
  use poisson_sete_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_dim_m
  use unit_system_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)
  use xc_functl_m
  use xc_ks_inversion_m
  use xc_OEP_m

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_messages_info,    &
    v_ks_calc,          &
    v_ks_calc_t,        &
    v_ks_calc_start,    &
    v_ks_calc_finish,   &
    v_ks_freeze_hxc

  integer, parameter, public :: &
    SIC_NONE   = 1,     &  ! no self-interaction correction
    SIC_PZ     = 2,     &  ! Perdew-Zunger SIC (OEP way)
    SIC_AMALDI = 3         ! Amaldi correction term

  type v_ks_calc_t
    private
    logical                       :: calculating
    logical                       :: time_present
    FLOAT                         :: time
    FLOAT,                pointer :: density(:, :)
    logical                       :: total_density_alloc
    FLOAT,                pointer :: total_density(:)
    FLOAT                         :: amaldi_factor
    type(energy_t),       pointer :: energy
    type(states_t),       pointer :: hf_st
    FLOAT,                pointer :: vxc(:, :)
    FLOAT,                pointer :: vtau(:, :)
    FLOAT,                pointer :: axc(:, :, :)
    FLOAT,                pointer :: vberry(:, :)
    FLOAT,                pointer :: a_ind(:, :)
    FLOAT,                pointer :: b_ind(:, :)
    logical                       :: calc_energy
  end type v_ks_calc_t

  type v_ks_t
    integer :: theory_level

    logical :: frozen_hxc ! For RPA and SAE calculations.

    integer                  :: xc_family  ! the XC stuff
    integer                  :: sic_type   ! what kind of self-interaction correction to apply
    type(xc_t)               :: xc
    type(xc_OEP_t)           :: oep
    type(xc_ks_inversion_t)  :: ks_inversion
    type(poisson_t), pointer :: hartree_solver
    logical                  :: new_hartree
    logical                  :: tail_correction
    FLOAT                    :: tail_correction_tol
    FLOAT                    :: tc_link_factor
    FLOAT                    :: tc_distance 
    integer                  :: tc_delay
    type(grid_t), pointer    :: gr
    type(v_ks_calc_t)        :: calc
  end type v_ks_t

contains

  
  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, gr, dd, geo, mc, nel)
    type(v_ks_t),         intent(out)   :: ks
    type(grid_t), target, intent(inout) :: gr
    type(states_dim_t),   intent(in)    :: dd
    type(geometry_t),     intent(inout) :: geo
    type(multicomm_t),    intent(in)    :: mc  
    FLOAT,                intent(in)    :: nel ! the total number of electrons

    PUSH_SUB(v_ks_init)

    !%Variable TheoryLevel
    !%Type integer
    !%Default dft
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with different "theory levels":
    !%Option independent_particles 2
    !% Particles will be considered as independent, <i>i.e.</i> as non-interacting.
    !% This mode is mainly used for testing purposes, as the code is usually 
    !% much faster with <tt>independent_particles</tt>.
    !%Option hartree 1
    !% Calculation within the Hartree method (experimental). Note that, contrary to popular
    !% belief, the Hartree potential is self-interaction-free. Therefore, this run 
    !% mode will not yield the same result as <tt>dft</tt> without exchange-correlation.
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
    !%Option classical 5
    !% (Experimental) Only the classical interaction between ions is
    !% considered. This is mainly for testing.
    !%End
    call parse_integer(datasets_check('TheoryLevel'), KOHN_SHAM_DFT, ks%theory_level)
    if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call input_error('TheoryLevel')

    call messages_obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable('HartreeFock', 'TheoryLevel')
    if(ks%theory_level == CLASSICAL) call messages_experimental('Classical theory level')
    
    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%sic_type = SIC_NONE
    case(HARTREE)
      call messages_experimental("Hartree theory level")
    case(HARTREE_FOCK)
      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, dd%spin_channels, dd%cdft, hartree_fock=.true.)
      ks%xc_family = ks%xc%family
      ks%sic_type = SIC_NONE

    case(KOHN_SHAM_DFT)
      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, dd%spin_channels, dd%cdft, hartree_fock=.false.)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA) .ne. 0) then

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
        !% Amaldi correction term.
        !%End
        call parse_integer(datasets_check('SICCorrection'), sic_none, ks%sic_type)
        if(.not. varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')

        ! Perdew-Zunger corrections
        if(ks%sic_type == SIC_PZ) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
      
      else
        ks%sic_type = SIC_NONE
      end if

      !%Variable XCTailCorrection
      !%Type logical
      !%Default no
      !%Section Hamiltonian::XC
      !%Description
      !% (Experimental) This variable enables to apply a correction to
      !% the value of the XC functional in near-zero-density regions.
      !% This zone might have numerical noise or it might 
      !% even be set to zero by <tt>libxc</tt>.
      !% The correction is performed by forcing the "-1/r behaviour" of the XC potential
      !% in the zones where the density is lower then XCTailCorrectionTol.
      !%End
      call parse_logical(datasets_check('XCTailCorrection'), .false., ks%tail_correction)
      
      if(ks%tail_correction) then 
        call messages_experimental("XC tail correction")
        write(message(1),'(a)') 'This correction shouldn''t be used with systems having nodal points of the electron density.'
        call messages_info(1)
        
        !%Variable XCTailCorrectionTol
        !%Type float
        !%Default 5-e12
        !%Section Hamiltonian::XC
        !%Description
        !%This variable sets the total electronic density threshold corresponding
        !% to the starting point of the <tt>XCTailCorrection</tt>.
        !%The value is always assumed to be in atomic units.
        !%End
        call parse_float(datasets_check('XCTailCorrectionTol'), CNST(5e-12), ks%tail_correction_tol)
      
        
        !%Variable XCTailCorrectionLinkFactor
        !%Type float
        !%Default 1
        !%Section Hamiltonian::XC
        !%Description
        !% (Experimental) This variable force a smooth transition beetween the region where the values of the XC functional 
        !% have been previously calculated and the region where the -1/r correction has been applied. 
        !% The region of the transition starts where the electronic total density reaches the value of
        !% (XCTailCorrectionLinkFactor * XCTailCorrectionTol) and ends where the density reaches the value of XCTailCorrectionTol
        !%End
        call parse_float(datasets_check('XCTailCorrectionLinkFactor'), CNST(1.0), ks%tc_link_factor)
      
        !%Variable XCTailCorrectionDelay
        !%Type integer 
        !%Default 0
        !%Section Hamiltonian::XC
        !%Description
        !% (Experimental) This variable allows to skip the application of the tail correction during the first calls of the
        !% subroutine that build the exchange-correlation potential (XCTailCorrectionDelay = number of calls skipped):
        !%this can avoid problems caused by eventual initial guess wavefunctions.
        !%End
        call parse_integer(datasets_check('XCTailCorrectionDelay'), 0, ks%tc_delay)

        !%Variable XCTailCorrectionCMDistance
        !%Type integer 
        !%Default 0
        !%Section Hamiltonian::XC
        !%Description
        !% (Experimental) This variable allows the application of the tail correction to the xc potential only where
        !% the distance of the local point from the center of mass of the system is greater than XCTailCorrectionCMDistance
        !%End
        call parse_float(datasets_check('XCTailCorrectionCMDistance'), M_ZERO, ks%tc_distance)
      end if

      
      



      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_init(ks%oep, ks%xc_family, gr, dd)
      endif
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_init(ks%ks_inversion, ks%xc_family, gr, geo, dd, mc)
      endif
    end select

    ks%frozen_hxc = .false.

    call v_ks_messages_info(ks, stdout)

    ks%new_hartree = .false.
    nullify(ks%hartree_solver)
    if(ks%theory_level /= INDEPENDENT_PARTICLES) then
      if(gr%have_fine_mesh) then
        ks%new_hartree = .true.
        SAFE_ALLOCATE(ks%hartree_solver)
        call poisson_init(ks%hartree_solver, gr%fine%der, geo, mc%master_comm)
      else
        ks%hartree_solver => psolver
      end if
     end if
     
     ks%gr => gr
     ks%calc%calculating = .false.

    POP_SUB(v_ks_init)
  end subroutine v_ks_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks, gr, geo)
    type(v_ks_t),     intent(inout) :: ks
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(inout) :: geo

    PUSH_SUB(v_ks_end)

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_end(ks%ks_inversion, gr, geo)
      endif
      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_end(ks%oep)
      endif
      call xc_end(ks%xc)
    end select

    if(ks%new_hartree) then
      SAFE_DEALLOCATE_P(ks%hartree_solver)
    end if

    POP_SUB(v_ks_end)
  end subroutine v_ks_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_messages_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(v_ks_messages_info)

    call messages_print_stress(iunit, "Theory Level")
    call messages_print_var_option(iunit, "TheoryLevel", ks%theory_level)

    select case(ks%theory_level)
    case(HARTREE_FOCK)
      write(iunit, '(1x)')
      call xc_messages_info(ks%xc, iunit)

    case(KOHN_SHAM_DFT)
      write(iunit, '(1x)')
      call xc_messages_info(ks%xc, iunit)

      write(iunit, '(1x)')
      call messages_print_var_option(iunit, 'SICCorrection', ks%sic_type)

      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_messages_info(ks%oep, iunit)
      end if
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_messages_info(ks%ks_inversion, iunit)
      end if
    end select

    call messages_print_stress(iunit)

    POP_SUB(v_ks_messages_info)
  end subroutine v_ks_messages_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, hm, st, geo, calc_eigenval, time, calc_berry, calc_energy)
    type(v_ks_t),           intent(inout) :: ks
    type(hamiltonian_t),    intent(inout) :: hm
    type(states_t),         intent(inout) :: st
    type(geometry_t), optional, intent(in)    :: geo
    logical,      optional, intent(in)    :: calc_eigenval
    FLOAT,        optional, intent(in)    :: time
    logical,      optional, intent(in)    :: calc_berry ! use this before wfns initialized
    logical,      optional, intent(in)    :: calc_energy
    
    call v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy)
    call v_ks_calc_finish(ks, hm)

    if(optional_default(calc_eigenval, .false.)) then

      if(ks%gr%ob_grid%open_boundaries .and. .not. present(time)) then
        ! We know the eigenvalues.
        st%eigenval(1:st%nst, 1:st%d%nik) = st%ob_eigenval(1:st%nst, 1:st%d%nik)
      else
        call energy_calculate_eigenvalues(hm, ks%gr%der, st)
      end if
      
    end if

  end subroutine v_ks_calc

  ! --------------------------------------------------------- 

  !> This routine starts the calculation of the Kohn-Sham
  !! potential. The routine v_ks_calc_finish must be called to finish
  !! the calculation. The argument hm is not modified. The argument st
  !! can be modified after the function have been used.

  subroutine v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy) 
    type(v_ks_t),                      intent(inout) :: ks 
    type(hamiltonian_t),     target,   intent(in)    :: hm !< This MUST be intent(in), changes to hm are done in v_ks_calc_finish.
    type(states_t),                    intent(inout) :: st
    type(geometry_t) ,       optional, intent(in)    :: geo
    FLOAT,                   optional, intent(in)    :: time 
    logical,                 optional, intent(in)    :: calc_berry !< Use this before wfns initialized.
    logical,                 optional, intent(in)    :: calc_energy

    
    type(profile_t), save :: prof
    type(energy_t), pointer :: energy

    PUSH_SUB(v_ks_calc_start)
    call profiling_in(prof, "KOHN_SHAM_CALC")

    ASSERT(.not. ks%calc%calculating)
    ks%calc%calculating = .true.

    if(in_debug_mode) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call messages_info(1)
    end if

    ks%calc%time_present = present(time)
    if(present(time)) ks%calc%time = time
    ks%calc%calc_energy = optional_default(calc_energy, .true.)

    ! If the Hxc term is frozen, there is nothing to do (WARNING: MISSING ks%calc%energy%intnvxc)
    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_start)
      return
    end if

    SAFE_ALLOCATE(ks%calc%energy)
    energy => ks%calc%energy
    
    call energy_copy(hm%energy, ks%calc%energy)
    
    energy%intnvxc = M_ZERO

    ! check whether we should introduce the Amaldi SIC correction
    ks%calc%amaldi_factor = M_ONE
    if(ks%sic_type == SIC_AMALDI) ks%calc%amaldi_factor = (st%qtot - M_ONE)/st%qtot

    if(ks%theory_level /= INDEPENDENT_PARTICLES .and. ks%calc%amaldi_factor /= M_ZERO) then

      call calculate_density()

      if(poisson_is_async(ks%hartree_solver)) then
        call dpoisson_solve_start(ks%hartree_solver, ks%calc%total_density)
      end if

      if(ks%theory_level .ne. HARTREE) call v_a_xc(geo)
    end if

    if(associated(hm%ep%e_field) .and. simul_box_is_periodic(ks%gr%mesh%sb)) then
      SAFE_ALLOCATE(ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin))
      if(optional_default(calc_berry, .true.)) then
        call berry_potential(st, ks%gr%mesh, hm%ep%E_field, ks%calc%vberry)
      else
        ! before wfns are initialized, cannot calculate this term
        ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = M_ZERO
      endif
    endif

    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK) then
      SAFE_ALLOCATE(ks%calc%hf_st)
      call states_null(ks%calc%hf_st)
      call states_copy(ks%calc%hf_st, st)
    end if

    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(ks%calc%a_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      SAFE_ALLOCATE(ks%calc%b_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      call magnetic_induced(ks%gr%der, st, ks%calc%a_ind, ks%calc%b_ind)
    end if

    call profiling_out(prof)
    POP_SUB(v_ks_calc_start)

  contains
    
    subroutine calculate_density()
      integer :: ip

      PUSH_SUB(v_ks_calc_start.calculate_density)

      ! get density taking into account non-linear core corrections
      SAFE_ALLOCATE(ks%calc%density(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      call states_total_density(st, ks%gr%fine%mesh, ks%calc%density)

      ! Amaldi correction
      if(ks%sic_type == SIC_AMALDI) then
        ks%calc%density(1:ks%gr%fine%mesh%np, 1:st%d%nspin) = &
          ks%calc%amaldi_factor*ks%calc%density(1:ks%gr%fine%mesh%np, 1:st%d%nspin)
      end if

      if(associated(st%rho_core) .or. hm%d%spin_channels > 1) then
        ks%calc%total_density_alloc = .true.

        SAFE_ALLOCATE(ks%calc%total_density(1:ks%gr%fine%mesh%np))

        forall(ip = 1:ks%gr%fine%mesh%np)
          ks%calc%total_density(ip) = sum(ks%calc%density(ip, 1:hm%d%spin_channels))
        end forall

        ! remove non-local core corrections
        if(associated(st%rho_core)) then
          forall(ip = 1:ks%gr%fine%mesh%np)
            ks%calc%total_density(ip) = ks%calc%total_density(ip) - st%rho_core(ip)*ks%calc%amaldi_factor
          end forall
        end if
      else
        ks%calc%total_density_alloc = .false.
        ks%calc%total_density => ks%calc%density(:, 1)
      end if

      POP_SUB(v_ks_calc_start.calculate_density)
    end subroutine calculate_density

    ! ---------------------------------------------------------
    subroutine v_a_xc(geo)
      type(geometry_t), optional, intent(in) :: geo
      type(profile_t), save :: prof
      

      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      energy%exchange = M_ZERO
      energy%correlation = M_ZERO
      energy%xc_j = M_ZERO

      SAFE_ALLOCATE(ks%calc%vxc(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      ks%calc%vxc = M_ZERO

      if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) then
        SAFE_ALLOCATE(ks%calc%vtau(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        ks%calc%vtau = M_ZERO
      end if

      ! Get the *local* XC term
      if(hm%d%cdft) then
        SAFE_ALLOCATE(ks%calc%axc(1:ks%gr%mesh%np, 1:ks%gr%sb%dim, 1:hm%d%nspin))

        ks%calc%axc = M_ZERO
        call xc_get_vxc_and_axc(ks%gr%fine%der, ks%xc, st, ks%calc%density, st%current, st%d%ispin, ks%calc%vxc, ks%calc%axc, &
             energy%exchange, energy%correlation, energy%xc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
      else if(ks%calc%calc_energy) then
        call xc_get_vxc(ks%gr%fine%der, ks%xc, st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst, :)), st%qtot, &
          ex = energy%exchange, ec = energy%correlation, vxc = ks%calc%vxc, vtau = ks%calc%vtau)
      else
        call xc_get_vxc(ks%gr%fine%der, ks%xc, st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst, :)), st%qtot, &
          vxc = ks%calc%vxc, vtau = ks%calc%vtau)
      end if

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
          if (states_are_real(st)) then
            call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type == SIC_PZ),  &
              ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
          else
            call zxc_oep_calc(ks%oep, ks%xc, (ks%sic_type == SIC_PZ),  &
              ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
          end if
        endif

        if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
          ! Also treat KS inversion separately (not part of libxc)
          call xc_ks_inversion_calc(ks%ks_inversion, ks%gr, hm, st, energy%exchange, energy%correlation, vxc = ks%calc%vxc)
        endif
      end if

      if(ks%tail_correction) then 
        ASSERT(present(geo))
        call tail_correction(ks%calc%vxc, geo)
      end if
      
      if(ks%calc%calc_energy) then
        ! Now we calculate Int[n vxc] = energy%intnvxc
        select case(hm%d%ispin)
        case(UNPOLARIZED)
          energy%intnvxc = energy%intnvxc + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), ks%calc%vxc(:, 1))
        case(SPIN_POLARIZED)
          energy%intnvxc = energy%intnvxc + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), ks%calc%vxc(:, 1)) &
            + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 2), ks%calc%vxc(:, 2))
        case(SPINORS)
          energy%intnvxc = energy%intnvxc + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), ks%calc%vxc(:, 1)) &
            + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 2), ks%calc%vxc(:, 2)) &
            + M_TWO*dmf_dotp(ks%gr%fine%mesh, st%rho(:, 3), ks%calc%vxc(:, 3)) &
            + M_TWO*dmf_dotp(ks%gr%fine%mesh, st%rho(:, 4), ks%calc%vxc(:, 4))
        end select
      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc

    !---------------------------------------------


    subroutine tail_correction(vxc,geo)

      FLOAT,             intent(inout) :: vxc(:, :) 
      type(geometry_t),  intent(in) :: geo
      
      FLOAT :: pos(MAX_DIM), distance_origin, distance_cm
      FLOAT, allocatable :: vxcc(:)
      FLOAT ::  s_dens,vnew 
      FLOAT  :: smooth_ratio
      integer :: nspin, is, ip, idim, ik, ist, ik_tmp  
      integer :: ierr, itmp
      integer , save :: counter = 0
      character(len=10) :: vxc_name
      logical :: to_calc
      

      PUSH_SUB(v_ks_calc_start.tail_correction)
      
      SAFE_ALLOCATE(vxcc(1:ks%gr%fine%mesh%np))
      
      counter = counter + 1
      !print *, "vxc tail correction call number" , counter
      nspin = 1
      if (st%d%ispin == SPIN_POLARIZED) nspin = 2
      
      
      spin_cycle: do is = 1,nspin
        kpoint_cycle: do ik_tmp = st%d%kpt%start, st%d%kpt%end, nspin
          ik = ik_tmp + is - 1
          state_cycle: do ist = st%st_start, st%st_end
            
            to_calc = .false.
            if ((st%occ(ist,ik) .ne. M_ZERO) .and. (counter .gt. ks%tc_delay) ) to_calc = .true.
            
            ! "If the state is not occupied and the call counter is greater than the desired value
            ! don't apply the correction and cycle
            if (to_calc .eqv. .false.) cycle 
            
            vxcc(1:ks%gr%fine%mesh%np) = vxc(1:ks%gr%fine%mesh%np, is)
            
            ! some lines for debugging:
            ! set the name of the output file and print the XC potential before the tail correction 
            ! write (vxc_name,'(i10)') is
            ! itmp = verify(vxc_name," ")
            ! vxc_name =  "vxc"//trim(vxc_name(itmp:))
            ! call doutput_function(output_axis_x, "./static", trim(vxc_name)//trim("precorr"), &
            ! ks%gr%fine%mesh, vxcc, unit_one, ierr)
            
            do ip = 1, ks%gr%fine%mesh%np
              s_dens = ks%calc%density(ip, is)
              
              distance_cm = M_ZERO
              call cm_pos(geo,pos)
              pos = M_ZERO
              if (ks%gr%mesh%sb%dim .ne. MAX_DIM) then
                write(message(1),'(a)') "Error: simulation box dimension different from the dimension used in subroutine cm_pos"
                call messages_info(1)
                stop
              end if
              
              do idim = 1,ks%gr%mesh%sb%dim 
                distance_cm = distance_cm +  (ks%gr%fine%mesh%x(ip,idim) - pos(idim) )**2
              end do
              distance_cm = sqrt(distance_cm)
              
              !If:
              !- the density is not exactly zero (and)...
              !- at least we have reached the linking region (and)...
              !- the cm distance is greater than the desired value ...
              ! then apply the correction
              if( (s_dens .ne. M_ZERO) .and. (ks%calc%total_density(ip) .lt. ks%tc_link_factor*ks%tail_correction_tol) & 
                .and. (distance_cm .gt. ks%tc_distance) ) then
                
                if (ks%calc%total_density(ip) .gt. ks%tail_correction_tol ) then 
                  smooth_ratio =  ks%calc%total_density(ip)/ks%tail_correction_tol
                  vnew = - 1/distance_cm
                  vxcc(ip) = ( smooth_ratio*vxcc(ip) + (ks%tc_link_factor-smooth_ratio)*vnew ) / ks%tc_link_factor
                else 
                  vxcc(ip) = -1/distance_cm
                end if
              end if
            end do
            
            !Another output call for debugging
            !!print the XC potential after the tail correction
            !call doutput_function(output_axis_x, "./static", trim(vxc_name)//trim("tcorrected") , &
            !ks%gr%fine%mesh, vxcc, unit_one, ierr)
            
            vxc(1:ks%gr%fine%mesh%np, is) = vxcc(1:ks%gr%fine%mesh%np)
            
            
          end do state_cycle
        end do kpoint_cycle
      end do spin_cycle
      
      
      SAFE_DEALLOCATE_A(vxcc)
      
      POP_SUB(v_ks_calc_start.tail_correction)
      
    end subroutine tail_correction
    
    
  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks, hm)
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: hm

    integer :: ip, ispin

    PUSH_SUB(v_ks_calc_finish)

    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_finish)
      return
    end if

    ASSERT(ks%calc%calculating)
    ks%calc%calculating = .false.

    !change the pointer to the energy object
    SAFE_DEALLOCATE_P(hm%energy)
    hm%energy => ks%calc%energy

    if(associated(hm%ep%e_field) .and. simul_box_is_periodic(ks%gr%mesh%sb)) then
      hm%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin)
      SAFE_DEALLOCATE_P(ks%calc%vberry)
    endif

    if(hm%self_induced_magnetic) then
      hm%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)
      hm%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)

      SAFE_DEALLOCATE_P(ks%calc%a_ind)
      SAFE_DEALLOCATE_P(ks%calc%b_ind)
    end if

    if(ks%theory_level == INDEPENDENT_PARTICLES .or. ks%calc%amaldi_factor == M_ZERO) then

      hm%vhxc = M_ZERO
      hm%energy%intnvxc     = M_ZERO
      hm%energy%hartree     = M_ZERO
      hm%energy%exchange    = M_ZERO
      hm%energy%correlation = M_ZERO

    else

      if(ks%theory_level /= HARTREE) then
        if(ks%gr%have_fine_mesh) then
          do ispin = 1, hm%d%nspin
            call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, &
              ks%calc%vxc(:, ispin), hm%vxc(:, ispin), INJECTION)
            ! some debugging output that I will keep here for the moment, XA
            !          call dio_function_output(1, "./", "vxc_fine", ks%gr%fine%mesh, vxc(:, ispin), unit_one, ierr)
            !          call dio_function_output(1, "./", "vxc_coarse", ks%gr%mesh, hm%vxc(:, ispin), unit_one, ierr)
          end do
          SAFE_DEALLOCATE_P(ks%calc%vxc)
        else
          ! just change the pointer to avoid the copy
          SAFE_DEALLOCATE_P(hm%vxc)
          hm%vxc => ks%calc%vxc
        end if

        if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) then
          do ispin = 1, hm%d%nspin
            call lalg_copy(ks%gr%fine%mesh%np, ks%calc%vtau(:, ispin), hm%vtau(:, ispin))
          end do
          SAFE_DEALLOCATE_P(ks%calc%vtau)
        end if

        if(hm%d%cdft) then
          hm%axc(1:ks%gr%mesh%np, 1:ks%gr%sb%dim, 1:hm%d%nspin) = ks%calc%axc(1:ks%gr%mesh%np, 1:ks%gr%sb%dim, 1:hm%d%nspin)
          SAFE_DEALLOCATE_P(ks%calc%axc)
        end if

      else
        hm%vxc = M_ZERO
      end if

      hm%energy%hartree = M_ZERO
      call v_ks_hartree(ks, hm)

      ! Build Hartree + XC potential
      forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)
      
      if(hm%d%ispin > UNPOLARIZED) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
      end if
      
      if(hm%d%ispin == SPINORS) then
        forall(ispin = 3:4, ip = 1:ks%gr%mesh%np) hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
      end if

      if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK) then

        ! swap the states object
        call states_end(hm%hf_st)
        SAFE_DEALLOCATE_P(hm%hf_st)
        hm%hf_st => ks%calc%hf_st

        select  case(ks%theory_level)
        case(HARTREE_FOCK)
          hm%exx_coef = ks%xc%exx_coef
        case(HARTREE)
          hm%exx_coef = M_ONE
        end select
      end if
      
    end if

    if(ks%calc%time_present) then
      call hamiltonian_update(hm, ks%gr%mesh, ks%calc%time)
    else
      call hamiltonian_update(hm, ks%gr%mesh)
    end if

    SAFE_DEALLOCATE_P(ks%calc%density)
    if(ks%calc%total_density_alloc) then
      SAFE_DEALLOCATE_P(ks%calc%total_density)
    end if

    POP_SUB(v_ks_calc_finish)
  end subroutine v_ks_calc_finish

  ! --------------------------------------------------------- 
  !
  ! Hartree contribution to the KS potential. This function is
  ! designed to be used by v_ks_calc_finish and it cannot be called
  ! directly.
  !
  subroutine v_ks_hartree(ks, hm)
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: hm

    FLOAT, pointer :: pot(:)

    PUSH_SUB(v_ks_hartree)

    ASSERT(associated(ks%hartree_solver))

    if(.not. ks%gr%have_fine_mesh) then
      pot => hm%vhartree
    else
      SAFE_ALLOCATE(pot(1:ks%gr%fine%mesh%np_part))
      pot = M_ZERO
    end if

    if(.not. poisson_is_async(ks%hartree_solver)) then
      ! solve the Poisson equation
      call dpoisson_solve(ks%hartree_solver, pot, ks%calc%total_density)
    else
      ! The calculation was started by v_ks_calc_start.
      call dpoisson_solve_finish(ks%hartree_solver, pot)
    end if

    if(ks%calc%calc_energy) then
      ! Get the Hartree energy
      hm%energy%hartree = M_HALF*dmf_dotp(ks%gr%fine%mesh, ks%calc%total_density, pot)
    end if

    if(ks%gr%have_fine_mesh) then
      ! we use injection to transfer to the fine grid, we cannot use
      ! restriction since the boundary conditions are not zero for the
      ! Hartree potential (and for some XC functionals).
      call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, pot, hm%vhartree, INJECTION)
      ! some debugging output that I will keep here for the moment, XA
      !      call dio_function_output(1, "./", "vh_fine", ks%gr%fine%mesh, pot, unit_one, is)
      !      call dio_function_output(1, "./", "vh_coarse", ks%gr%mesh, hm%vhartree, unit_one, is)
      SAFE_DEALLOCATE_P(pot)
    end if

    if (ks%calc%calc_energy .and. poisson_get_solver(ks%hartree_solver) == POISSON_SETE) then !SEC
      hm%energy%hartree = hm%energy%hartree + poisson_energy(ks%hartree_solver)
      write(89,*) hm%energy%hartree*CNST(2.0)*CNST(13.60569193), poisson_energy(ks%hartree_solver)*CNST(2.0)*CNST(13.60569193), &
        hm%ep%eii*CNST(2.0)*CNST(13.60569193)
    endif

    POP_SUB(v_ks_hartree)
  end subroutine v_ks_hartree
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_freeze_hxc(ks)
    type(v_ks_t), intent(inout) :: ks

    PUSH_SUB(v_ks_freeze_hxc)

    ks%frozen_hxc = .true.
    
    POP_SUB(v_ks_freeze_hxc)
  end subroutine v_ks_freeze_hxc
  ! ---------------------------------------------------------

end module v_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
