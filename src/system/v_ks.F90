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
 
module v_ks_oct_m
  use berry_oct_m
  use current_oct_m
  use density_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use energy_calc_oct_m
  use epot_oct_m 
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use hamiltonian_base_oct_m
  use kick_oct_m
  use index_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use lasers_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use pseudo_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_parallel_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use XC_F90(lib_m)
  use xc_functl_oct_m

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_calc,          &
    v_ks_calc_t,        &
    v_ks_calc_start,    &
    v_ks_calc_finish,   &
    v_ks_freeze_hxc,    &
    v_ks_calculate_current 

  integer, parameter, public :: &
    SIC_NONE   = 1,     &  !< no self-interaction correction
    SIC_AMALDI = 3,     &  !< Amaldi correction term
    SIC_ADSIC  = 4         !< Averaged density SIC

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
    logical                       :: calc_energy
    type(geometry_t),     pointer :: geo
  end type v_ks_calc_t

  type v_ks_t
    integer :: theory_level

    logical :: frozen_hxc !< For RPA and SAE calculations.

    integer                  :: xc_family  !< the XC stuff
    integer                  :: xc_flags   !< the XC flags
    integer                  :: sic_type   !< what kind of self-interaction correction to apply
    type(xc_t)               :: xc
    type(poisson_t), pointer :: hartree_solver
    logical                  :: new_hartree
    type(grid_t), pointer    :: gr
    type(v_ks_calc_t)        :: calc
    logical                  :: calculate_current
    type(current_t)          :: current_calculator
    logical                  :: include_td_field
  end type v_ks_t

contains

  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, gr, st, geo, mc)
    type(v_ks_t),         intent(inout)   :: ks
    type(grid_t), target, intent(inout) :: gr
    type(states_t),       intent(in)    :: st
    type(geometry_t),     intent(inout) :: geo
    type(multicomm_t),    intent(in)    :: mc  

    integer :: x_id, c_id, xk_id, ck_id, default, val, iatom
    logical :: parsed_theory_level
    integer :: pseudo_x_functional, pseudo_c_functional
    
    PUSH_SUB(v_ks_init)

    ! We need to parse TheoryLevel and XCFunctional, this is
    ! complicated because they are interdependent.
    
    !%Variable TheoryLevel
    !%Type integer
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with different "theory levels" that
    !% control how electrons are simulated. The default is
    !% <tt>dft</tt>. When hybrid functionals are requested, through
    !% the <tt>XCFunctional</tt> variable, the default is
    !% <tt>hartree_fock</tt>.
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
    
    ks%xc_family = XC_FAMILY_NONE
    ks%sic_type  = SIC_NONE
   
    ks%theory_level = KOHN_SHAM_DFT
    parsed_theory_level = .false.
    
    ! the user knows what he wants, give her that
    if(parse_is_defined('TheoryLevel')) then
      call parse_variable('TheoryLevel', KOHN_SHAM_DFT, ks%theory_level)
      if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call messages_input_error('TheoryLevel')

      parsed_theory_level = .true.
    end if

    ! parse the XC functional
    default = 0

    call get_functional_from_pseudos(pseudo_x_functional, pseudo_c_functional)

    if(ks%theory_level == KOHN_SHAM_DFT) then
      if(pseudo_x_functional /= PSEUDO_EXCHANGE_ANY) then
        default = pseudo_x_functional
      else
        select case(gr%mesh%sb%dim)
        case(3); default = XC_LDA_X   
        case(2); default = XC_LDA_X_2D
        case(1); default = XC_LDA_X_1D
        end select
      end if
    end if
    
    ASSERT(default >= 0)

    if(ks%theory_level == KOHN_SHAM_DFT) then
      if(pseudo_c_functional /= PSEUDO_CORRELATION_ANY) then
        default = default + 1000*pseudo_c_functional
      else
        select case(gr%mesh%sb%dim)
        case(3); default = default + 1000*XC_LDA_C_PZ_MOD
        case(2); default = default + 1000*XC_LDA_C_2D_AMGB
        case(1); default = default + 1000*XC_LDA_C_1D_CSC
        end select
      end if
    end if

    ASSERT(default >= 0)

    if(.not. parse_is_defined('XCFunctional') &
      .and. (pseudo_x_functional /= PSEUDO_EXCHANGE_ANY .or. pseudo_c_functional /= PSEUDO_CORRELATION_ANY)) then
      call messages_write('Info: the XCFunctional has been selected to match the pseudopotentials', new_line = .true.)
      call messages_write('      used in the calculation.')
      call messages_info()
    end if
    
    ! The description of this variable can be found in file src/xc/functionals_list.F90
    call parse_variable('XCFunctional', default, val)

    ! the first 3 digits of the number indicate the X functional and
    ! the next 3 the C functional.
    c_id = val / 1000
    x_id = val - c_id*1000
    
    if( (x_id /= pseudo_x_functional .and. pseudo_x_functional /= PSEUDO_EXCHANGE_ANY) .or. &
      (c_id /= pseudo_c_functional .and. pseudo_c_functional /= PSEUDO_EXCHANGE_ANY)) then
      call messages_write('The XCFunctional that you selected does not match the one used', new_line = .true.)
      call messages_write('to generate the pseudopotentials.')
      call messages_warning()
    end if

    ! FIXME: we rarely need this. We should only parse when necessary.
    
    !%Variable XCKernel
    !%Type integer
    !%Section Hamiltonian::XC
    !%Description
    !% Defines the exchange-correlation kernel. Only LDA kernels are available currently.
    !% The options are the same as <tt>XCFunctional</tt>.
    !% Note: the kernel is only needed for Casida, Sternheimer, or optimal-control calculations.
    !% Defaults:
    !% <br>1D: <tt>lda_x_1d + lda_c_1d_csc</tt>
    !% <br>2D: <tt>lda_x_2d + lda_c_2d_amgb</tt>
    !% <br>3D: <tt>lda_x + lda_c_pz_mod</tt>
    !%Option xc_functional -1
    !% The same functional defined by <tt>XCFunctional</tt>.
    !%End
    
    call parse_variable('XCKernel', default, val)
    
    if( -1 == val ) then
      ck_id = c_id
      xk_id = x_id
    else
      ck_id = val / 1000
      xk_id = val - ck_id*1000  
    end if
    
    call messages_obsolete_variable('XFunctional', 'XCFunctional')
    call messages_obsolete_variable('CFunctional', 'XCFunctional')

    ! initialize XC modules

    ! This is a bit ugly, theory_level might not be Hartree-Fock now
    ! but it might become Hartree-Fock later. This is safe because it
    ! becomes Hartree-Fock in the cases where the functional is hybrid
    ! and the ifs inside check for both conditions.
    call xc_init(ks%xc, gr%mesh%sb%dim, gr%mesh%sb%periodic_dim, st%qtot, &
      x_id, c_id, xk_id, ck_id, hartree_fock = ks%theory_level == HARTREE_FOCK)

    ks%xc_family = ks%xc%family
    ks%xc_flags  = ks%xc%flags 

    if(.not. parsed_theory_level) then
      default = KOHN_SHAM_DFT

      ! the functional is a hybrid, use Hartree-Fock as theory level by default
      if(bitand(ks%xc_family, XC_FAMILY_HYB_GGA + XC_FAMILY_HYB_MGGA) /= 0) then
        default = HARTREE_FOCK
      end if

      ! In principle we do not need to parse. However we do it for consistency
      call parse_variable('TheoryLevel', default, ks%theory_level)
      if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call messages_input_error('TheoryLevel')
     
    end if

    call messages_obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable('HartreeFock', 'TheoryLevel')

    if(ks%theory_level == CLASSICAL) call messages_experimental('Classical theory level')
    
    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%sic_type = SIC_NONE
    case(HARTREE)
      call messages_experimental("Hartree theory level")
      if(gr%mesh%sb%periodic_dim == gr%mesh%sb%dim) &
        call messages_experimental("Hartree in fully periodic system")
      if(gr%mesh%sb%kpoints%full%npoints > 1) &
        call messages_not_implemented("Hartree with k-points")

    case(HARTREE_FOCK)
      if(gr%mesh%sb%kpoints%full%npoints > 1) &
        call messages_not_implemented("Hartree-Fock with k-points")
      
      ks%sic_type = SIC_NONE

    case(KOHN_SHAM_DFT)

      ! check for SIC
      if(bitand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA) /= 0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which form of self-interaction correction to use. Note that
        !% this correction will be applied to the functional chosen by <tt>XCFunctional</tt>.
        !%Option sic_none 1
        !% No self-interaction correction.
        !%Option sic_amaldi 3
        !% Amaldi correction term.
        !%Option sic_adsic 4
        !% Average-density SIC.
        !% C. Legrand <i>et al.</i>, <i>J. Phys. B</i> <b>35</b>, 1115 (2002). 
        !%End
        call parse_variable('SICCorrection', sic_none, ks%sic_type)
        if(.not. varinfo_valid_option('SICCorrection', ks%sic_type)) call messages_input_error('SICCorrection')

      else
        ks%sic_type = SIC_NONE
      end if

    end select

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, stdout)

    ks%new_hartree = .false.
    nullify(ks%hartree_solver)
    if(ks%theory_level /= INDEPENDENT_PARTICLES) then
      if(gr%have_fine_mesh) then
        ks%new_hartree = .true.
        SAFE_ALLOCATE(ks%hartree_solver)
        call poisson_init(ks%hartree_solver, gr%fine%der, mc, label = " (fine mesh)")
      else
        ks%hartree_solver => psolver
      end if
    end if

    ks%gr => gr
    ks%calc%calculating = .false.

    !The value of ks%calculate_current is set to false or true by Output    
    call current_init(ks%current_calculator)
    
    POP_SUB(v_ks_init)

  contains

    subroutine get_functional_from_pseudos(x_functional, c_functional)
      integer, intent(out) :: x_functional
      integer, intent(out) :: c_functional
      
      integer :: xf, cf, ispecies
      logical :: warned_inconsistent
      
      x_functional = PSEUDO_EXCHANGE_ANY
      c_functional = PSEUDO_CORRELATION_ANY
      
      warned_inconsistent = .false.
      do ispecies = 1, geo%nspecies
        xf = species_x_functional(geo%species(ispecies))
        cf = species_c_functional(geo%species(ispecies))

        if(xf == PSEUDO_EXCHANGE_UNKNOWN .or. cf == PSEUDO_CORRELATION_UNKNOWN) then
          call messages_write("Unknown XC functional for species '"//trim(species_label(geo%species(ispecies)))//"'")
          call messages_warning()
          cycle
        end if

        if(x_functional == PSEUDO_EXCHANGE_ANY) then
          x_functional = xf
        else
          if(xf /= x_functional .and. .not. warned_inconsistent) then
            call messages_write('Inconsistent XC functional detected between species');
            call messages_warning()
            warned_inconsistent = .true.
          end if
        end if

        if(c_functional == PSEUDO_CORRELATION_ANY) then
          c_functional = cf
        else
          if(cf /= c_functional .and. .not. warned_inconsistent) then
            call messages_write('Inconsistent XC functional detected between species');
            call messages_warning()
            warned_inconsistent = .true.
          end if
        end if
        
      end do
      
      ASSERT(x_functional /= PSEUDO_EXCHANGE_UNKNOWN)
      ASSERT(c_functional /= PSEUDO_CORRELATION_UNKNOWN)
      
    end subroutine get_functional_from_pseudos
  end subroutine v_ks_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks)
    type(v_ks_t),     intent(inout) :: ks

    PUSH_SUB(v_ks_end)
    
    call current_end(ks%current_calculator)

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      call xc_end(ks%xc)
    case(HARTREE_FOCK)      
      call xc_end(ks%xc)
    end select

    if(ks%new_hartree) then
      call poisson_end(ks%hartree_solver)
      SAFE_DEALLOCATE_P(ks%hartree_solver)
    end if

    POP_SUB(v_ks_end)
  end subroutine v_ks_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(v_ks_write_info)

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

    end select

    call messages_print_stress(iunit)

    POP_SUB(v_ks_write_info)
  end subroutine v_ks_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, hm, st, geo, calc_eigenval, time, calc_berry, calc_energy, calc_current)
    type(v_ks_t),               intent(inout) :: ks
    type(hamiltonian_t),        intent(inout) :: hm
    type(states_t),             intent(inout) :: st
    type(geometry_t),           intent(in)    :: geo
    logical,          optional, intent(in)    :: calc_eigenval
    FLOAT,            optional, intent(in)    :: time
    logical,          optional, intent(in)    :: calc_berry !< use this before wfns initialized
    logical,          optional, intent(in)    :: calc_energy
    logical,          optional, intent(in)    :: calc_current

    logical :: calc_current_

    PUSH_SUB(v_ks_calc)

    calc_current_ = optional_default(calc_current, .true.)

    call v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy, calc_current_)
    call v_ks_calc_finish(ks, hm)

    if(optional_default(calc_eigenval, .false.)) then
      call energy_calc_eigenvalues(hm, ks%gr%der, st)
    end if

    POP_SUB(v_ks_calc)
  end subroutine v_ks_calc

  ! --------------------------------------------------------- 

  !> This routine starts the calculation of the Kohn-Sham
  !! potential. The routine v_ks_calc_finish must be called to finish
  !! the calculation. The argument hm is not modified. The argument st
  !! can be modified after the function have been used.
  subroutine v_ks_calc_start(ks, hm, st, geo, time, calc_berry, calc_energy, calc_current) 
    type(v_ks_t),            target,   intent(inout) :: ks 
    type(hamiltonian_t),     target,   intent(in)    :: hm !< This MUST be intent(in), changes to hm are done in v_ks_calc_finish.
    type(states_t),                    intent(inout) :: st
    type(geometry_t) ,       target,   intent(in)    :: geo
    FLOAT,                   optional, intent(in)    :: time 
    logical,                 optional, intent(in)    :: calc_berry !< Use this before wfns initialized.
    logical,                 optional, intent(in)    :: calc_energy
    logical,                 optional, intent(in)    :: calc_current

    type(profile_t), save :: prof
    logical :: calc_current_

    PUSH_SUB(v_ks_calc_start)

    calc_current_ = optional_default(calc_current, .true.)

    call profiling_in(prof, "KOHN_SHAM_CALC")

    ASSERT(.not. ks%calc%calculating)
    ks%calc%calculating = .true.

    ks%calc%geo => geo
    
    if(debug%info) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call messages_info(1)
    end if

    ks%calc%time_present = present(time) 

    if(present(time)) then
      ks%calc%time = time
    end if

    ks%calc%calc_energy = optional_default(calc_energy, .true.)

    nullify(ks%calc%vberry)
    if(associated(hm%vberry)) then
      SAFE_ALLOCATE(ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin))
      if(optional_default(calc_berry, .true.)) then
        call berry_potential(st, ks%gr%mesh, hm%ep%E_field, ks%calc%vberry)
      else
        ! before wfns are initialized, cannot calculate this term
        ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = M_ZERO
      end if
    end if

    ! If the Hxc term is frozen, there is nothing more to do (WARNING: MISSING ks%calc%energy%intnvxc)
    if(ks%frozen_hxc) then      
      if(ks%calculate_current .and. calc_current_ ) then
        call states_allocate_current(st, ks%gr)
        call current_calculate(ks%current_calculator, ks%gr%der, hm, geo, st, st%current, st%current_kpt)
      end if      

      POP_SUB(v_ks_calc_start)
      return
    end if

    SAFE_ALLOCATE(ks%calc%energy)

    call energy_copy(hm%energy, ks%calc%energy)

    ks%calc%energy%intnvxc = M_ZERO

    ! check whether we should introduce the Amaldi SIC correction
    ks%calc%amaldi_factor = M_ONE
    if(ks%sic_type == SIC_AMALDI) ks%calc%amaldi_factor = (st%qtot - M_ONE)/st%qtot

    nullify(ks%calc%density, ks%calc%total_density)
    nullify(ks%calc%vxc, ks%calc%vtau, ks%calc%axc)

    if(ks%theory_level /= INDEPENDENT_PARTICLES .and. ks%calc%amaldi_factor /= M_ZERO) then

      call calculate_density()

      if(poisson_is_async(ks%hartree_solver)) then
        call dpoisson_solve_start(ks%hartree_solver, ks%calc%total_density)
      end if

      if(ks%theory_level /= HARTREE) call v_a_xc(hm)
    else
      ks%calc%total_density_alloc = .false.
    end if

    if(ks%calculate_current .and. calc_current_ ) then
      call states_allocate_current(st, ks%gr)
      call current_calculate(ks%current_calculator, ks%gr%der, hm, geo, st, st%current, st%current_kpt)
    end if

    nullify(ks%calc%hf_st) 
    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK) then
      SAFE_ALLOCATE(ks%calc%hf_st)
      call states_copy(ks%calc%hf_st, st)
      if(st%parallel_in_states) call states_parallel_remote_access_start(ks%calc%hf_st)
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
      if(ks%sic_type == SIC_AMALDI) &
        ks%calc%density = ks%calc%amaldi_factor*ks%calc%density

      nullify(ks%calc%total_density)
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

    !ADSIC potential is:
    !V_ADSIC[n] = V_ks[n] - (V_h[n/N] - V_xc[n/N])
    subroutine add_adsic(hm)
      type(hamiltonian_t), intent(in)    :: hm

      integer        :: ip, ispin, ist, ik
      FLOAT, pointer :: vxc_sic(:,:),  Imvxc_sic(:, :), vh_sic(:), rho(:, :), Imrho(:, :), qsp(:)
      
      PUSH_SUB(add_adsic)
      
      if(family_is_mgga(hm%xc_family)) then
        call messages_not_implemented('ADSIC with MGGAs')
      end if
      if (st%d%ispin == SPINORS) then
        call messages_not_implemented('ADSIC with non-collinear spin')      
      end if

      SAFE_ALLOCATE(vxc_sic(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(vh_sic(1:ks%gr%mesh%np))
      SAFE_ALLOCATE(rho(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      SAFE_ALLOCATE(qsp(1:st%d%nspin))
      
      vxc_sic = M_ZERO
      vh_sic = M_ZERO
      qsp = M_ZERO
      do ist = 1, st%nst
        do ispin = 1,st%d%nspin
          do ik = ispin, st%d%nik, st%d%nspin
           qsp(ispin) = qsp(ispin)+ st%occ(ist, ik) * st%d%kweights(ik) 
          enddo
        enddo
      end do

      select case (st%d%ispin)
      case (UNPOLARIZED, SPIN_POLARIZED)
        do ispin = 1, st%d%nspin
          if (abs(qsp(ispin)) <= M_EPSILON) cycle

          rho = M_ZERO
          vxc_sic = M_ZERO

          rho(:, ispin) = ks%calc%density(:, ispin) / qsp(ispin)
          ! TODO : check for solid:   -minval(st%eigenval(st%nst,:))
          call xc_get_vxc(ks%gr%fine%der, ks%xc, &
               st, rho, st%d%ispin, -minval(st%eigenval(st%nst,:)), qsp(ispin), &
               vxc_sic)

          ks%calc%vxc = ks%calc%vxc - vxc_sic
        end do

      case (SPINORS)
        !TODO
      end select

      rho(:, 1) = ks%calc%total_density / st%qtot
      call dpoisson_solve(ks%hartree_solver, vh_sic, rho(:,1))
      forall(ip = 1:ks%gr%mesh%np) ks%calc%vxc(ip,:) = ks%calc%vxc(ip,:) - vh_sic(ip)

      SAFE_DEALLOCATE_P(vxc_sic)
      SAFE_DEALLOCATE_P(vh_sic)                                
      SAFE_DEALLOCATE_P(rho)
      SAFE_DEALLOCATE_P(qsp)

      POP_SUB(add_adsic)
    end subroutine add_adsic


    ! ---------------------------------------------------------
    subroutine v_a_xc(hm)
      type(hamiltonian_t),  intent(in) :: hm 

      type(profile_t), save :: prof
      FLOAT :: factor
      CMPLX :: ctmp
      integer :: ispin, iatom
      FLOAT, allocatable :: coords(:, :)
      FLOAT :: latvec(1:3, 1:3)
      integer, allocatable :: atnum(:)

      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      ks%calc%energy%exchange = M_ZERO
      ks%calc%energy%correlation = M_ZERO
      ks%calc%energy%xc_j = M_ZERO

      SAFE_ALLOCATE(ks%calc%vxc(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      ks%calc%vxc = M_ZERO

      nullify(ks%calc%vtau)
      if(hm%family_is_mgga_with_exc) then
        SAFE_ALLOCATE(ks%calc%vtau(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        ks%calc%vtau = M_ZERO
      end if

      ! Get the *local* XC term
      if(ks%calc%calc_energy) then
        if(hm%family_is_mgga_with_exc) then
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, ks%calc%vxc, &
            ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, deltaxc = ks%calc%energy%delta_xc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr%fine%der, ks%xc, &
            st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, ks%calc%vxc, &
            ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, deltaxc = ks%calc%energy%delta_xc)
        end if
      else
        if(hm%family_is_mgga_with_exc) then
          call xc_get_vxc(ks%gr%fine%der, ks%xc, &
            st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, &
            ks%calc%vxc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr%fine%der, ks%xc, &
            st, ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, &
            ks%calc%vxc)
        end if
      end if

      if (ks%sic_type == SIC_ADSIC) then
        call add_adsic(hm)
      end if

      if(ks%calc%calc_energy) then
        ! Now we calculate Int[n vxc] = energy%intnvxc
        ks%calc%energy%intnvxc = M_ZERO

        do ispin = 1, hm%d%nspin
          if(ispin <= 2) then
            factor = M_ONE
          else
            factor = M_TWO
          end if
          ks%calc%energy%intnvxc = ks%calc%energy%intnvxc + &
            factor*dmf_dotp(ks%gr%fine%mesh, st%rho(:, ispin), ks%calc%vxc(:, ispin))
        end do

      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc

  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks, hm)
    type(v_ks_t), target, intent(inout) :: ks
    type(hamiltonian_t),  intent(inout) :: hm

    integer                           :: ip, ispin

    PUSH_SUB(v_ks_calc_finish)

    ASSERT(ks%calc%calculating)
    ks%calc%calculating = .false.

    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_finish)
      return
    end if

    !change the pointer to the energy object
    SAFE_DEALLOCATE_P(hm%energy)
    hm%energy => ks%calc%energy

    if(associated(hm%vberry)) then
      hm%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = ks%calc%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin)
      SAFE_DEALLOCATE_P(ks%calc%vberry)
    end if

    if(ks%theory_level == INDEPENDENT_PARTICLES .or. abs(ks%calc%amaldi_factor) <= M_EPSILON) then

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

        if(hm%family_is_mgga_with_exc) then
          do ispin = 1, hm%d%nspin
            call lalg_copy(ks%gr%fine%mesh%np, ks%calc%vtau(:, ispin), hm%vtau(:, ispin))
          end do
          SAFE_DEALLOCATE_P(ks%calc%vtau)
        end if

      else
        hm%vxc = M_ZERO
      end if

      hm%energy%hartree = M_ZERO
      call v_ks_hartree(ks, hm)


      ! Build Hartree + XC potential
     
      forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)
      if(associated(hm%vberry)) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vhxc(ip, 1) + hm%vberry(ip, 1)
      end if

      if(hm%d%ispin > UNPOLARIZED) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
        if(associated(hm%vberry)) then
          forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vhxc(ip, 2) + hm%vberry(ip, 2)
        end if
      end if
      
      if(hm%d%ispin == SPINORS) then
        forall(ispin = 3:4, ip = 1:ks%gr%mesh%np) hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
      end if

      ! Note: this includes hybrids calculated with the Fock operator instead of OEP 
      if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK) then

        ! swap the states object
        if(associated(hm%hf_st)) then
          if(hm%hf_st%parallel_in_states) call states_parallel_remote_access_stop(hm%hf_st)
          call states_end(hm%hf_st)
          SAFE_DEALLOCATE_P(hm%hf_st)
        end if
        
        hm%hf_st => ks%calc%hf_st

        select case(ks%theory_level)
        case(HARTREE_FOCK)
          hm%exx_coef = ks%xc%exx_coef
        case(HARTREE)
          hm%exx_coef = M_ONE
        end select
      end if
      
    end if

    if(ks%calc%time_present) then
      call hamiltonian_update(hm, ks%gr%mesh, ks%gr%der%boundaries, time = ks%calc%time)
    else
      call hamiltonian_update(hm, ks%gr%mesh, ks%gr%der%boundaries)
    end if


    SAFE_DEALLOCATE_P(ks%calc%density)
    if(ks%calc%total_density_alloc) then
      SAFE_DEALLOCATE_P(ks%calc%total_density)
    end if
    nullify(ks%calc%total_density)

    POP_SUB(v_ks_calc_finish)
  end subroutine v_ks_calc_finish

  ! --------------------------------------------------------- 
  !
  !> Hartree contribution to the KS potential. This function is
  !! designed to be used by v_ks_calc_finish and it cannot be called
  !! directly.
  !
  subroutine v_ks_hartree(ks, hm)
    type(v_ks_t),                intent(inout) :: ks
    type(hamiltonian_t), target, intent(inout) :: hm

    FLOAT, pointer :: pot(:), aux(:)

    FLOAT, allocatable :: potx(:)
    CMPLX, allocatable :: kick(:)
    FLOAT, allocatable :: kick_real(:)
    integer :: ii

    integer :: asc_unit_test

    FLOAT :: dt

    logical :: kick_time

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

  subroutine v_ks_calculate_current(this, calc_cur)
    type(v_ks_t), intent(inout) :: this
    logical,      intent(in)    :: calc_cur

    PUSH_SUB(v_ks_calculate_current)
    
    this%calculate_current = calc_cur

    POP_SUB(v_ks_calculate_current)
  end subroutine v_ks_calculate_current
 
end module v_ks_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
