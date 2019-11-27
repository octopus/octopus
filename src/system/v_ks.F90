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
  use accel_oct_m
  use berry_oct_m
  use comm_oct_m
  use current_oct_m
  use density_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use energy_calc_oct_m
  use epot_oct_m
  use exchange_operator_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use kick_oct_m
  use lalg_basic_oct_m
  use lasers_oct_m
  use libvdwxc_oct_m
  use magnetic_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use multigrid_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use pseudo_oct_m
  use pcm_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_parallel_oct_m
  use varinfo_oct_m
  use vdw_ts_oct_m
  use xc_oct_m
  use XC_F90(lib_m)
  use xc_functl_oct_m
  use xc_ks_inversion_oct_m
  use xc_OEP_oct_m

  ! from the dftd3 library
  use dftd3_api
  
  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_nullify,       &
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
    SIC_PZ     = 2,     &  !< Perdew-Zunger SIC (OEP way)
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
    type(states_elec_t),  pointer :: hf_st
    FLOAT,                pointer :: vxc(:, :)
    FLOAT,                pointer :: vtau(:, :)
    FLOAT,                pointer :: axc(:, :, :)
    FLOAT,                pointer :: vberry(:, :)
    FLOAT,                pointer :: a_ind(:, :)
    FLOAT,                pointer :: b_ind(:, :)
    logical                       :: calc_energy

    FLOAT, allocatable :: vdw_forces(:, :)
    type(geometry_t), pointer :: geo
  end type v_ks_calc_t

  type v_ks_t
    private
    integer,                  public :: theory_level

    logical,                  public :: frozen_hxc !< For RPA and SAE calculations.

    integer,                  public :: xc_family  !< the XC stuff
    integer,                  public :: xc_flags   !< the XC flags
    integer,                  public :: sic_type   !< what kind of self-interaction correction to apply
    type(xc_t),               public :: xc
    type(xc_OEP_t),           public :: oep
    type(xc_ks_inversion_t),  public :: ks_inversion
    type(grid_t), pointer,    public :: gr
    type(v_ks_calc_t)                :: calc
    logical                          :: calculate_current
    type(current_t)                  :: current_calculator
    integer,                  public :: vdw_correction
    logical                          :: vdw_self_consistent
    type(vdw_ts_t),           public :: vdw_ts
    type(dftd3_calc)                 :: vdw_d3
    logical                          :: include_td_field
  end type v_ks_t

contains
 
  ! ---------------------------------------------------------
  subroutine v_ks_nullify(ks)
    type(v_ks_t),            intent(inout) :: ks

    PUSH_SUB(v_ks_nullify)

    ks%theory_level = -1
    ks%frozen_hxc = .false.
    ks%xc_family = 0
    ks%xc_flags = 0
    ks%sic_type = -1
    ks%calculate_current = .false.
    ks%vdw_correction = -1
    ks%vdw_self_consistent = .false.
    ks%include_td_field = .false.

    POP_SUB(v_ks_nullify)
  end subroutine v_ks_nullify
  

  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, namespace, gr, st, geo, mc)
    type(v_ks_t),            intent(inout) :: ks
    type(namespace_t),       intent(in)    :: namespace
    type(grid_t),    target, intent(inout) :: gr
    type(states_elec_t),     intent(in)    :: st
    type(geometry_t),        intent(inout) :: geo
    type(multicomm_t),       intent(in)    :: mc

    integer :: x_id, c_id, xk_id, ck_id, default, val, iatom
    logical :: parsed_theory_level
    type(dftd3_input) :: d3_input
    character(len=20) :: d3func_def, d3func
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
    !%Option rdmft 7 
    !% (Experimental) Reduced Density Matrix functional theory.
    !%End
    
    ks%xc_family = XC_FAMILY_NONE
    ks%sic_type  = SIC_NONE
   
    ks%theory_level = KOHN_SHAM_DFT
    parsed_theory_level = .false.
    call exchange_operator_nullify(exxop)
    
    ! the user knows what he wants, give her that
    if(parse_is_defined(namespace, 'TheoryLevel')) then
      call parse_variable(namespace, 'TheoryLevel', KOHN_SHAM_DFT, ks%theory_level)
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

    if(.not. parse_is_defined(namespace, 'XCFunctional') &
      .and. (pseudo_x_functional /= PSEUDO_EXCHANGE_ANY .or. pseudo_c_functional /= PSEUDO_CORRELATION_ANY)) then
      call messages_write('Info: the XCFunctional has been selected to match the pseudopotentials', new_line = .true.)
      call messages_write('      used in the calculation.')
      call messages_info()
    end if
    
    ! The description of this variable can be found in file src/xc/functionals_list.F90
    call parse_variable(namespace, 'XCFunctional', default, val)

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
    
    call parse_variable(namespace, 'XCKernel', default, val)
    
    if( -1 == val ) then
      ck_id = c_id
      xk_id = x_id
    else
      ck_id = val / 1000
      xk_id = val - ck_id*1000  
    end if
    
    call messages_obsolete_variable(namespace, 'XFunctional', 'XCFunctional')
    call messages_obsolete_variable(namespace, 'CFunctional', 'XCFunctional')

    ! initialize XC modules

    ! This is a bit ugly, theory_level might not be Hartree-Fock now
    ! but it might become Hartree-Fock later. This is safe because it
    ! becomes Hartree-Fock in the cases where the functional is hybrid
    ! and the ifs inside check for both conditions.
    call xc_init(ks%xc, namespace, gr%mesh%sb%dim, gr%mesh%sb%periodic_dim, st%qtot, &
      x_id, c_id, xk_id, ck_id, hartree_fock = ks%theory_level == HARTREE_FOCK)

    if(bitand(ks%xc%family, XC_FAMILY_LIBVDWXC) /= 0) then
      call libvdwxc_set_geometry(ks%xc%functional(FUNC_C,1)%libvdwxc, namespace, gr%mesh)
    end if

    ks%xc_family = ks%xc%family
    ks%xc_flags  = ks%xc%flags 

    if(.not. parsed_theory_level) then
      default = KOHN_SHAM_DFT

      ! the functional is a hybrid, use Hartree-Fock as theory level by default
      if(bitand(ks%xc_family, XC_FAMILY_HYB_GGA + XC_FAMILY_HYB_MGGA) /= 0) then
        default = HARTREE_FOCK
      end if

      ! In principle we do not need to parse. However we do it for consistency
      call parse_variable(namespace, 'TheoryLevel', default, ks%theory_level)
      if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call messages_input_error('TheoryLevel')
     
    end if

    call messages_obsolete_variable(namespace, 'NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable(namespace, 'HartreeFock', 'TheoryLevel')

    if(ks%theory_level == RDMFT ) call messages_experimental('RDMFT theory level')
    
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
        !%Option sic_pz 2
        !% Perdew-Zunger SIC, handled by the OEP technique.
        !%Option sic_amaldi 3
        !% Amaldi correction term.
        !%Option sic_adsic 4
        !% Average-density SIC.
        !% C. Legrand <i>et al.</i>, <i>J. Phys. B</i> <b>35</b>, 1115 (2002). 
        !%End
        call parse_variable(namespace, 'SICCorrection', sic_none, ks%sic_type)
        if(.not. varinfo_valid_option('SICCorrection', ks%sic_type)) call messages_input_error('SICCorrection')

        ! Perdew-Zunger corrections
        if(ks%sic_type == SIC_PZ) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)

      else
        ks%sic_type = SIC_NONE
      end if

      if(bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        if (gr%have_fine_mesh) call messages_not_implemented("OEP functionals with UseFineMesh")

        call xc_oep_init(ks%oep, namespace, ks%xc_family, gr, st)
      end if

      if(bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_init(ks%ks_inversion, namespace, gr, geo, st, ks%xc, mc)
      end if

    end select

    if (st%d%ispin == SPINORS) then
      if(bitand(ks%xc_family, XC_FAMILY_GGA + XC_FAMILY_HYB_GGA) /= 0) call messages_not_implemented("GGA with spinors")
      if(bitand(ks%xc_family, XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA) /= 0) call messages_not_implemented("MGGA with spinors")
    end if

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, stdout)

    ks%gr => gr
    ks%calc%calculating = .false.

    !The value of ks%calculate_current is set to false or true by Output    
    call current_init(ks%current_calculator, namespace)
    
    !%Variable VDWCorrection
    !%Type integer
    !%Default no
    !%Section Hamiltonian::XC
    !%Description
    !% (Experimental) This variable selects which van der Waals
    !% correction to apply to the correlation functional.
    !%Option none 0
    !% No correction is applied.
    !%Option vdw_ts 1
    !% The scheme of Tkatchenko and Scheffler, Phys. Rev. Lett. 102
    !% 073005 (2009).
    !%Option vdw_d3 3
    !% The DFT-D3 scheme of S. Grimme, J. Antony, S. Ehrlich, and
    !% S. Krieg, J. Chem. Phys. 132, 154104 (2010).
    !%End
    call parse_variable(namespace, 'VDWCorrection', OPTION__VDWCORRECTION__NONE, ks%vdw_correction)
    
    if(ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      call messages_experimental('VDWCorrection')

      select case(ks%vdw_correction)
      case(OPTION__VDWCORRECTION__VDW_TS)

        !%Variable VDWSelfConsistent
        !%Type logical
        !%Default yes
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls whether the VDW correction is applied
        !% self-consistently, the default, or just as a correction to
        !% the total energy. This option only works with vdw_ts.
        !%End
        call parse_variable(namespace, 'VDWSelfConsistent', .true., ks%vdw_self_consistent)

        call vdw_ts_init(ks%vdw_ts, namespace, geo, gr%fine%der)

      case(OPTION__VDWCORRECTION__VDW_D3)
        ks%vdw_self_consistent = .false.

        if(ks%gr%sb%dim /= 3) then
          call messages_write('vdw_d3 can only be used in 3-dimensional systems')
          call messages_fatal()
        end if
        
        do iatom = 1, geo%natoms
          if(.not. species_represents_real_atom(geo%atom(iatom)%species)) then
            call messages_write('vdw_d3 is not implemented when non-atomic species are present')
            call messages_fatal()
          end if
        end do
         
        d3func_def = ''

        ! The list of valid values can be found in 'external_libs/dftd3/core.f90'.
        ! For the moment I include the most common ones.
        if(x_id == OPTION__XCFUNCTIONAL__GGA_X_B88 .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_LYP) &
          d3func_def = 'b-lyp'
        if(x_id == OPTION__XCFUNCTIONAL__GGA_X_PBE .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_PBE) &
          d3func_def = 'pbe'
        if(x_id == OPTION__XCFUNCTIONAL__GGA_X_PBE_SOL .and. c_id*1000 == OPTION__XCFUNCTIONAL__GGA_C_PBE_SOL) &
          d3func_def = 'pbesol'
        if(c_id*1000 == OPTION__XCFUNCTIONAL__HYB_GGA_XC_B3LYP) &
          d3func_def = 'b3-lyp'
        if(c_id*1000 == OPTION__XCFUNCTIONAL__HYB_GGA_XC_PBEH) &
          d3func_def = 'pbe0'

        !%Variable VDWD3Functional
        !%Type string
        !%Section Hamiltonian::XC
        !%Description
        !% (Experimental) You can use this variable to override the
        !% parametrization used by the DFT-D3 van deer Waals
        !% correction. Normally you need not set this variable, as the
        !% proper value will be selected by Octopus (if available).
        !%
        !% This variable takes a string value, the valid values can
        !% be found in the source file 'external_libs/dftd3/core.f90'.
        !% For example you can use:
        !%
        !%  VDWD3Functional = 'pbe'
        !%
        !%End
        if(parse_is_defined(namespace, 'VDWD3Functional')) call messages_experimental('VDWD3Functional')
        call parse_variable(namespace, 'VDWD3Functional', d3func_def, d3func)
        
        if(d3func == '') then
          call messages_write('Cannot find  a matching parametrization  of DFT-D3 for the current')
          call messages_new_line()
          call messages_write('XCFunctional.  Please select a different XCFunctional, or select a')
          call messages_new_line()
          call messages_write('functional for DFT-D3 using the <tt>VDWD3Functional</tt> variable.')
          call messages_fatal()
        end if

        if(ks%gr%sb%periodic_dim /= 0 .and. ks%gr%sb%periodic_dim /= 3) then
          call messages_write('For partially periodic systems,  the vdw_d3 interaction is assumed')
          call messages_new_line()
          call messages_write('to be periodic in three dimensions.')
          call messages_warning()
        end if
          
        call dftd3_init(ks%vdw_d3, d3_input, trim(conf%share)//'/dftd3/pars.dat')
        call dftd3_set_functional(ks%vdw_d3, func = d3func, version = 4, tz = .false.)

      case default
        ASSERT(.false.)
      end select
      
    else
      ks%vdw_self_consistent = .false.
    end if
    
    nullify(ks%calc%hf_st)
    
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
  subroutine v_ks_end(ks, gr)
    type(v_ks_t),     intent(inout) :: ks
    type(grid_t),     intent(inout) :: gr

    PUSH_SUB(v_ks_end)
    
    select case(ks%vdw_correction)
    case(OPTION__VDWCORRECTION__VDW_TS)
      call vdw_ts_end(ks%vdw_ts)
    end select

    call current_end(ks%current_calculator)

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      if(bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_end(ks%ks_inversion, gr)
      end if
      if(bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_end(ks%oep)
      end if
      call xc_end(ks%xc)
    case(HARTREE_FOCK)      
      call xc_end(ks%xc)
    end select

    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK .or. ks%theory_level == RDMFT) then
      if(associated(ks%calc%hf_st)) then
        if(ks%calc%hf_st%parallel_in_states) call states_elec_parallel_remote_access_stop(ks%calc%hf_st)
        call states_elec_end(ks%calc%hf_st)
        SAFE_DEALLOCATE_P(ks%calc%hf_st)
        nullify(ks%calc%hf_st)
      end if
    end if


    call exchange_operator_end(exxop)

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

      if(bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
        call xc_oep_write_info(ks%oep, iunit)
      end if
      if(bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
        call xc_ks_inversion_write_info(ks%ks_inversion, iunit)
      end if

    end select

    call messages_print_stress(iunit)

    POP_SUB(v_ks_write_info)
  end subroutine v_ks_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, namespace, hm, st, geo, calc_eigenval, time, calc_berry, calc_energy, calc_current)
    type(v_ks_t),               intent(inout) :: ks
    type(namespace_t),          intent(in)    :: namespace
    type(hamiltonian_elec_t),   intent(inout) :: hm
    type(states_elec_t),        intent(inout) :: st
    type(geometry_t),           intent(in)    :: geo
    logical,          optional, intent(in)    :: calc_eigenval
    FLOAT,            optional, intent(in)    :: time
    logical,          optional, intent(in)    :: calc_berry !< use this before wfns initialized
    logical,          optional, intent(in)    :: calc_energy
    logical,          optional, intent(in)    :: calc_current

    logical :: calc_current_

    PUSH_SUB(v_ks_calc)

    calc_current_ = optional_default(calc_current, .true.)

    call v_ks_calc_start(ks, namespace, hm, st, geo, time, calc_berry, calc_energy, calc_current_)
    call v_ks_calc_finish(ks, hm, namespace)

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
  subroutine v_ks_calc_start(ks, namespace, hm, st, geo, time, calc_berry, calc_energy, calc_current) 
    type(v_ks_t),              target, intent(inout) :: ks
    type(namespace_t),                 intent(in)    :: namespace
    type(hamiltonian_elec_t),  target, intent(in)    :: hm !< This MUST be intent(in), changes to hm are done in v_ks_calc_finish.
    type(states_elec_t),               intent(inout) :: st
    type(geometry_t) ,         target, intent(in)    :: geo
    FLOAT,                   optional, intent(in)    :: time 
    logical,                 optional, intent(in)    :: calc_berry !< Use this before wfns initialized.
    logical,                 optional, intent(in)    :: calc_energy
    logical,                 optional, intent(in)    :: calc_current

    type(profile_t), save :: prof
    logical :: calc_current_

    PUSH_SUB(v_ks_calc_start)

    calc_current_ = optional_default(calc_current, .true.)  &
                   .and. ks%calculate_current &
                   .and. states_are_complex(st) &
                   .or. hamiltonian_elec_needs_current(hm, states_are_real(st))

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
      if(calc_current_ ) then
        call states_elec_allocate_current(st, ks%gr)
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

      if(poisson_is_async(hm%psolver_fine)) then
        call dpoisson_solve_start(hm%psolver_fine, ks%calc%total_density)
      end if

      if(ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT) call v_a_xc(hm)
    else
      ks%calc%total_density_alloc = .false.
    end if

    if(calc_current_ ) then
      call states_elec_allocate_current(st, ks%gr)
      call current_calculate(ks%current_calculator, ks%gr%der, hm, geo, st, st%current, st%current_kpt)
    end if

    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK .or. ks%theory_level == RDMFT) then
      if(.not. associated(ks%calc%hf_st)) then
        SAFE_ALLOCATE(ks%calc%hf_st)
        call states_elec_copy(ks%calc%hf_st, st)

        if(st%parallel_in_states) then
          if(accel_is_enabled()) then
            call messages_write('State parallelization of Hartree-Fock exchange  is not supported')
            call messages_new_line()
            call messages_write('when running with OpenCL/CUDA. Please use domain parallelization')
            call messages_new_line()
            call messages_write("or disable acceleration using 'DisableAccel = yes'.")
            call messages_fatal()
          end if
          call states_elec_parallel_remote_access_start(ks%calc%hf_st)
        end if
      end if

      select case(ks%theory_level)
        case(HARTREE_FOCK)
          call exchange_operator_reinit(exxop, ks%calc%hf_st, ks%xc%cam_omega, ks%xc%cam_alpha, ks%xc%cam_beta)
        case(HARTREE)
          call exchange_operator_reinit(exxop, ks%calc%hf_st, M_ONE, M_ZERO, M_ZERO)
        case(RDMFT)
          call exchange_operator_reinit(exxop, ks%calc%hf_st, M_ONE, M_ZERO, M_ZERO)
      end select
    end if


    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    nullify(ks%calc%a_ind, ks%calc%b_ind)
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(ks%calc%a_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      SAFE_ALLOCATE(ks%calc%b_ind(1:ks%gr%mesh%np_part, 1:ks%gr%sb%dim))
      call magnetic_induced(ks%gr%der, st, hm%psolver, ks%calc%a_ind, ks%calc%b_ind)
    end if
   
    call profiling_out(prof)
    POP_SUB(v_ks_calc_start)

  contains

    subroutine calculate_density()
      integer :: ip

      PUSH_SUB(v_ks_calc_start.calculate_density)

      ! get density taking into account non-linear core corrections
      SAFE_ALLOCATE(ks%calc%density(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      call states_elec_total_density(st, ks%gr%fine%mesh, ks%calc%density)

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
      type(hamiltonian_elec_t), intent(in)    :: hm

      integer        :: ip, ispin, ist, ik
      FLOAT, pointer :: vxc_sic(:,:),  vh_sic(:), rho(:, :), qsp(:)
      
      PUSH_SUB(add_adsic)
      
      if (family_is_mgga(hm%xc%family)) then
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
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, hm%psolver_fine, namespace, rho, st%d%ispin, &
            -minval(st%eigenval(st%nst,:)), qsp(ispin), exxop, vxc_sic)

          ks%calc%vxc = ks%calc%vxc - vxc_sic
        end do

      case (SPINORS)
        !TODO
      end select

      rho(:, 1) = ks%calc%total_density / st%qtot
      call dpoisson_solve(hm%psolver_fine, vh_sic, rho(:,1))
      forall(ip = 1:ks%gr%mesh%np) ks%calc%vxc(ip,:) = ks%calc%vxc(ip,:) - vh_sic(ip)

      SAFE_DEALLOCATE_P(vxc_sic)
      SAFE_DEALLOCATE_P(vh_sic)                                
      SAFE_DEALLOCATE_P(rho)
      SAFE_DEALLOCATE_P(qsp)

      POP_SUB(add_adsic)
    end subroutine add_adsic


    ! ---------------------------------------------------------
    subroutine v_a_xc(hm)
      type(hamiltonian_elec_t),  intent(in) :: hm

      type(profile_t), save :: prof
      FLOAT :: factor
      CMPLX :: ctmp
      integer :: ispin, iatom
      FLOAT, allocatable :: vvdw(:)
      FLOAT, allocatable :: coords(:, :)
      FLOAT :: vdw_stress(1:3, 1:3), latvec(1:3, 1:3)
      integer, allocatable :: atnum(:)

      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      ks%calc%energy%exchange = M_ZERO
      ks%calc%energy%correlation = M_ZERO
      ks%calc%energy%xc_j = M_ZERO

      SAFE_ALLOCATE(ks%calc%vxc(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      ks%calc%vxc = M_ZERO

      nullify(ks%calc%vtau)
      if (family_is_mgga_with_exc(hm%xc)) then
        SAFE_ALLOCATE(ks%calc%vtau(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
        ks%calc%vtau = M_ZERO
      end if

      ! Get the *local* XC term
      if(ks%calc%calc_energy) then
        if (family_is_mgga_with_exc(hm%xc)) then
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, hm%psolver_fine, namespace, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, exxop, ks%calc%vxc, &
            ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, deltaxc = ks%calc%energy%delta_xc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, hm%psolver_fine, namespace, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, exxop, ks%calc%vxc, &
            ex = ks%calc%energy%exchange, ec = ks%calc%energy%correlation, deltaxc = ks%calc%energy%delta_xc)
        end if
      else
        if (family_is_mgga_with_exc(hm%xc)) then
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, hm%psolver_fine, namespace, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, exxop, &
            ks%calc%vxc, vtau = ks%calc%vtau)
        else
          call xc_get_vxc(ks%gr%fine%der, ks%xc, st, hm%psolver_fine, namespace, &
            ks%calc%density, st%d%ispin, -minval(st%eigenval(st%nst,:)), st%qtot, exxop, ks%calc%vxc)
        end if
      end if

      if (ks%sic_type == SIC_ADSIC) then
        call add_adsic(hm)
      end if

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if(bitand(ks%xc_family, XC_FAMILY_OEP) /= 0) then
          if (states_are_real(st)) then
            call dxc_oep_calc(ks%oep, namespace, ks%xc, (ks%sic_type == SIC_PZ), ks%gr, &
              hm, st, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
          else
            call zxc_oep_calc(ks%oep, namespace, ks%xc, (ks%sic_type == SIC_PZ), ks%gr, &
              hm, st, ks%calc%energy%exchange, ks%calc%energy%correlation, vxc = ks%calc%vxc)
          end if
        end if

        if(bitand(ks%xc_family, XC_FAMILY_KS_INVERSION) /= 0) then
          ! Also treat KS inversion separately (not part of libxc)
          call xc_ks_inversion_calc(ks%ks_inversion, namespace, ks%gr, hm, st, vxc = ks%calc%vxc, &
            time = ks%calc%time)
        end if
      end if

      if(ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
        ASSERT(geo%space%dim == 3)
        
        SAFE_ALLOCATE(vvdw(1:ks%gr%fine%mesh%np))
        SAFE_ALLOCATE(ks%calc%vdw_forces(1:geo%space%dim, 1:geo%natoms))

        select case(ks%vdw_correction)

        case(OPTION__VDWCORRECTION__VDW_TS)
          vvdw = CNST(0.0)
          call vdw_ts_calculate(ks%vdw_ts, namespace, geo, ks%gr%der, ks%gr%sb, st, st%rho, &
            ks%calc%energy%vdw, vvdw, ks%calc%vdw_forces)
           
        case(OPTION__VDWCORRECTION__VDW_D3)

          SAFE_ALLOCATE(coords(1:3, geo%natoms))
          SAFE_ALLOCATE(atnum(geo%natoms))

          do iatom = 1, geo%natoms
            atnum(iatom) = nint(species_z(geo%atom(iatom)%species))
            coords(1:3, iatom) = geo%atom(iatom)%x(1:3)
          end do
          
          if(simul_box_is_periodic(ks%gr%sb)) then
            latvec(1:3, 1:3) = ks%gr%sb%rlattice(1:3, 1:3) !make a copy as rlattice goes up to MAX_DIM
            call dftd3_pbc_dispersion(ks%vdw_d3, coords, atnum, latvec, ks%calc%energy%vdw, ks%calc%vdw_forces, vdw_stress)
          else
            call dftd3_dispersion(ks%vdw_d3, coords, atnum, ks%calc%energy%vdw, ks%calc%vdw_forces)
          end if

          SAFE_DEALLOCATE_A(coords)
          SAFE_DEALLOCATE_A(atnum)
          
        case default
          ASSERT(.false.)
          
        end select
        
        if(ks%vdw_self_consistent) then
          do ispin = 1, hm%d%nspin
            ks%calc%vxc(1:ks%gr%fine%mesh%np, ispin) = ks%calc%vxc(1:ks%gr%fine%mesh%np, ispin) + vvdw(1:ks%gr%fine%mesh%np)
          end do
        end if
        
        SAFE_DEALLOCATE_A(vvdw)    
        
      else

        ks%calc%energy%vdw = CNST(0.0)

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
            factor*dmf_dotp(ks%gr%fine%mesh, st%rho(:, ispin), ks%calc%vxc(:, ispin), reduce = .false.)
        end do
        if(ks%gr%der%mesh%parallel_in_domains) call comm_allreduce(ks%gr%der%mesh%mpi_grp%comm,  ks%calc%energy%intnvxc)

        if(states_are_real(st)) then
          ks%calc%energy%int_dft_u = denergy_calc_electronic(hm, ks%gr%der, st, terms = TERM_DFT_U)
        else
          ctmp = zenergy_calc_electronic(hm, ks%gr%der, st, terms = TERM_DFT_U)
          ks%calc%energy%int_dft_u   = real(ctmp)
        end if

      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc

  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks, hm, namespace)
    type(v_ks_t),     target, intent(inout) :: ks
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(namespace_t),        intent(in)    :: namespace

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

    if(hm%self_induced_magnetic) then
      hm%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%a_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)
      hm%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim) = ks%calc%b_ind(1:ks%gr%mesh%np, 1:ks%gr%sb%dim)

      SAFE_DEALLOCATE_P(ks%calc%a_ind)
      SAFE_DEALLOCATE_P(ks%calc%b_ind)
    end if

    if(associated(hm%ep%v_static)) then
      hm%energy%intnvstatic = dmf_dotp(ks%gr%mesh, ks%calc%total_density, hm%ep%v_static) 
    else
      hm%energy%intnvstatic = M_ZERO
    end if

    if(ks%theory_level == INDEPENDENT_PARTICLES .or. abs(ks%calc%amaldi_factor) <= M_EPSILON) then

      hm%vhxc = M_ZERO
      hm%energy%intnvxc     = M_ZERO
      hm%energy%hartree     = M_ZERO
      hm%energy%exchange    = M_ZERO
      hm%energy%correlation = M_ZERO
    else

      if(ks%theory_level /= HARTREE .and. ks%theory_level /= RDMFT ) then 
        if(ks%gr%have_fine_mesh) then
          do ispin = 1, hm%d%nspin
            call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, &
              ks%calc%vxc(:, ispin), hm%vxc(:, ispin), INJECTION)
            ! This output needs a namespace argument to work from now on. It
            ! hasn't been touched for some years now, so I won't adapt it. - SO
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

        if (family_is_mgga_with_exc(hm%xc)) then
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

    end if

    if(ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      hm%ep%vdw_forces(1:ks%gr%sb%dim, 1:ks%calc%geo%natoms) = ks%calc%vdw_forces(1:ks%gr%sb%dim, 1:ks%calc%geo%natoms)
      SAFE_DEALLOCATE_A(ks%calc%vdw_forces)
    else
      hm%ep%vdw_forces(1:ks%gr%sb%dim, 1:ks%calc%geo%natoms) = CNST(0.0)      
    end if

    if(ks%calc%time_present) then
      call hamiltonian_elec_update(hm, ks%gr%mesh, namespace, time = ks%calc%time)
    else
      call hamiltonian_elec_update(hm, ks%gr%mesh, namespace)
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
    type(hamiltonian_elec_t), target, intent(inout) :: hm

    FLOAT, pointer :: pot(:)

    FLOAT, allocatable :: potx(:)
    CMPLX, allocatable :: kick(:)
    FLOAT, allocatable :: kick_real(:)
    integer :: ii

    logical :: kick_time

    logical :: laser_present, kick_present

    PUSH_SUB(v_ks_hartree)

    if(.not. ks%gr%have_fine_mesh) then
      pot => hm%vhartree
    else
      SAFE_ALLOCATE(pot(1:ks%gr%fine%mesh%np_part))
      pot = M_ZERO
    end if

    if(.not. poisson_is_async(hm%psolver_fine)) then
      ! solve the Poisson equation
      call dpoisson_solve(hm%psolver_fine, pot, ks%calc%total_density)
    else
      ! The calculation was started by v_ks_calc_start.
      call dpoisson_solve_finish(hm%psolver_fine, pot)
    end if


    !> PCM reaction field due to the electronic density
    if (hm%pcm%run_pcm .and. pcm_update(hm%pcm,hm%current_time)) then
      ! Currently this PCM section seems to be inconsistent when one has a fine mesh

      !> Generates the real-space PCM potential due to electrons during the SCF calculation.
      if (hm%pcm%solute) then
        call pcm_calc_pot_rs(hm%pcm, ks%gr%mesh, hm%psolver_fine, v_h = pot, time_present = ks%calc%time_present)
      end if
        
      !> Local field effects due to the applied electrostatic potential representing the laser and the kick (if they were).
      !! For the laser, the latter is only valid in the long-wavelength limit.
      !! Static potentials are included in subroutine hamiltonian_elec_epot_generate (module hamiltonian).
      !! The sign convention for typical potentials and kick are different...
      if( hm%pcm%localf .and. ks%calc%time_present ) then
        laser_present = epot_have_lasers( hm%ep )
        kick_present  = epot_have_kick(   hm%ep )
        if ( laser_present .and. kick_present ) then !< external potential and kick
          SAFE_ALLOCATE(potx(1:ks%gr%mesh%np_part))
          SAFE_ALLOCATE(kick(1:ks%gr%mesh%np_part)) 
          SAFE_ALLOCATE(kick_real(1:ks%gr%mesh%np_part))
          potx = M_ZERO
          kick = M_ZERO
          do ii = 1, hm%ep%no_lasers        
            call laser_potential(hm%ep%lasers(ii), ks%gr%mesh, potx, ks%calc%time)
          end do
          kick_real = M_ZERO
          kick_time = ((hm%pcm%iter-1)*hm%pcm%dt <= hm%ep%kick%time) .and. (hm%pcm%iter*hm%pcm%dt > hm%ep%kick%time)
          if ( kick_time ) then
            call kick_function_get(ks%gr%mesh, hm%ep%kick, kick, to_interpolate = .true.)
            kick = hm%ep%kick%delta_strength * kick
            kick_real = DREAL(kick)
          end if
          call pcm_calc_pot_rs(hm%pcm, ks%gr%mesh, hm%psolver_fine, v_ext = potx, kick = -kick_real, &
            time_present = ks%calc%time_present, kick_time = kick_time )
          SAFE_DEALLOCATE_A(potx)
          SAFE_DEALLOCATE_A(kick)
          SAFE_DEALLOCATE_A(kick_real)
        else if ( laser_present .and. .not.kick_present ) then !< just external potential
          SAFE_ALLOCATE(potx(1:ks%gr%mesh%np_part))
          potx = M_ZERO    
          do ii = 1, hm%ep%no_lasers        
            call laser_potential(hm%ep%lasers(ii), ks%gr%mesh, potx, ks%calc%time)
          end do
          call pcm_calc_pot_rs(hm%pcm, ks%gr%mesh, hm%psolver_fine, v_ext = potx, time_present = ks%calc%time_present)
          SAFE_DEALLOCATE_A(potx)
        else if ( .not.laser_present .and. kick_present ) then !< just kick
          SAFE_ALLOCATE(kick(1:ks%gr%mesh%np_part))
          SAFE_ALLOCATE(kick_real(1:ks%gr%mesh%np_part))
          kick = M_ZERO
          kick_real = M_ZERO
          kick_time =((hm%pcm%iter-1)*hm%pcm%dt <= hm%ep%kick%time) .and. (hm%pcm%iter*hm%pcm%dt > hm%ep%kick%time)
          if ( kick_time ) then
            call kick_function_get(ks%gr%mesh, hm%ep%kick, kick, to_interpolate = .true.)
            kick = hm%ep%kick%delta_strength * kick
            kick_real = DREAL(kick)
          end if
          call pcm_calc_pot_rs(hm%pcm, ks%gr%mesh, hm%psolver_fine, kick = -kick_real, &
            time_present = ks%calc%time_present, kick_time = kick_time)
          SAFE_DEALLOCATE_A(kick)
          SAFE_DEALLOCATE_A(kick_real)
        end if

        ! Calculating the PCM term renormalizing the sum of the single-particle energies
        ! to keep the idea of pcm_corr... but it will be added later on
        hm%energy%pcm_corr = dmf_dotp( ks%gr%fine%mesh, ks%calc%total_density, hm%pcm%v_e_rs + hm%pcm%v_n_rs + hm%pcm%v_ext_rs )
      else
        ! Calculating the PCM term renormalizing the sum of the single-particle energies
        hm%energy%pcm_corr = dmf_dotp( ks%gr%fine%mesh, ks%calc%total_density, hm%pcm%v_e_rs + hm%pcm%v_n_rs )
      end if
        
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
      ! This output needs a namespace argument to work from now on. It
      ! hasn't been touched for some years now, so I won't adapt it. - SO
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
