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

module hamiltonian_mxll_oct_m
  use accel_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use boundaries_oct_m
  use boundary_op_oct_m
  use comm_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use energy_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hardware_oct_m
  use io_oct_m
  use io_function_oct_m
  use lalg_basic_oct_m
  use lasers_oct_m
  use math_oct_m
  use maxwell_boundary_op_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use parser_oct_m
  use partition_oct_m
  use par_vec_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_parallel_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use xc_functl_oct_m
  use XC_F90(lib_m)
  use xyz_adjust_oct_m

  implicit none

  private
  public ::                                     &
    hamiltonian_mxll_t,                         &
    hamiltonian_mxll_null,                      &
    hamiltonian_mxll_init,                      &
    hamiltonian_mxll_end,                       &
    hamiltonian_mxll_apply,                     &
    hamiltonian_mxll_apply_all,                 &
    hamiltonian_mxll_apply_batch,               &
    hamiltonian_mxll_adjoint,                   &
    hamiltonian_mxll_not_adjoint,               &
    hamiltonian_mxll_hermitian,                 &
    hamiltonian_mxll_update,                    &
    hamiltonian_mxll_get_time,                  &
    hamiltonian_mxll_apply_packed,              &
    maxwell_fft_hamiltonian,                    &
    get_electric_field_at_points,               &
    get_magnetic_field_at_points,               &
    maxwell_helmholtz_decomposition_trans_field,&
    maxwell_helmholtz_decomposition_long_field, &
    surface_integral_helmholtz_transverse


  type hamiltonian_mxll_t
    !> The Hamiltonian must know what are the "dimensions" of the spaces,
    !! in order to be able to operate on the states.
    type(states_dim_t)       :: d
    type(bc_t)               :: bc      !< boundaries

    !> absorbing boundaries
    logical :: adjoint

    FLOAT :: current_time
    logical :: apply_packed  !< This is initialized by the StatesPack variable.
    
    logical :: time_zero

    logical                        :: maxwell_ks_calculation

    type(poisson_t)                :: maxwell_poisson_solver
    FLOAT, pointer                 :: maxwell_vector_potential(:,:)

    type(maxwell_bc_t)             :: maxwell_bc
    type(grid_t), pointer          :: maxwell_gr
    type(states_t), pointer        :: maxwell_st

    integer                        :: maxwell_rs_sign

    logical                        :: maxwell_propagation_apply

    integer                        :: maxwell_op_method

    logical                        :: maxwell_lorentz_force
    logical                        :: maxwell_lorentz_force_apply

    integer, pointer               :: maxwell_rs_state_fft_map(:,:,:)
    integer, pointer               :: maxwell_rs_state_fft_map_inv(:,:)

    integer                        :: mx_ma_trans_field_calc_method
    logical                        :: mx_ma_trans_field_calc_corr

    logical                        :: maxwell_bc_add_ab_region  = .false.
    logical                        :: maxwell_bc_zero           = .false.
    logical                        :: maxwell_bc_constant       = .false.
    logical                        :: maxwell_bc_mirror_pec     = .false.
    logical                        :: maxwell_bc_mirror_pmc     = .false.
    logical                        :: maxwell_bc_periodic       = .false.
    logical                        :: maxwell_bc_plane_waves    = .false.
    logical                        :: maxwell_bc_medium         = .false.

    logical                        :: maxwell_plane_waves
    logical                        :: maxwell_plane_waves_apply
    logical                        :: maxwell_spatial_constant
    logical                        :: maxwell_spatial_constant_apply

    logical                        :: maxwell_diamagnetic_current
    logical                        :: maxwell_spin_current

    integer                        :: maxwell_medium_calculation

    logical                        :: maxwell_medium_box = .false.
    integer                        :: maxwell_medium_box_number
    integer, pointer               :: maxwell_medium_box_shape(:)
    FLOAT, pointer                 :: maxwell_medium_box_center(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_size(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_ep(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_mu(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_c(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_ep_factor(:)
    FLOAT, pointer                 :: maxwell_medium_box_mu_factor(:)
    FLOAT, pointer                 :: maxwell_medium_box_sigma_e_factor(:)
    FLOAT, pointer                 :: maxwell_medium_box_sigma_m_factor(:)
    FLOAT, pointer                 :: maxwell_medium_box_sigma_e(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_sigma_m(:,:)
    integer, pointer               :: maxwell_medium_box_points_number(:)
    FLOAT, pointer                 :: maxwell_medium_box_points_map(:,:)
    FLOAT, pointer                 :: maxwell_medium_box_aux_ep(:,:,:)
    FLOAT, pointer                 :: maxwell_medium_box_aux_mu(:,:,:)
    integer, pointer               :: maxwell_medium_box_bdry_number(:)
    FLOAT, pointer                 :: maxwell_medium_box_bdry_map(:,:)
  
    !> maxwell hamiltonian_mxll
    integer                        :: maxwell_operator
    logical                        :: maxwell_current_density_ext_flag
    FLOAT                          :: maxwell_energy
    FLOAT                          :: maxwell_energy_boundaries
    FLOAT                          :: maxwell_e_energy
    FLOAT                          :: maxwell_b_energy
    FLOAT                          :: maxwell_energy_plane_waves
    FLOAT                          :: maxwell_e_energy_plane_waves
    FLOAT                          :: maxwell_b_energy_plane_waves

    FLOAT                          :: maxwell_energy_pml
    FLOAT                          :: maxwell_energy_mask
    FLOAT, pointer                 :: maxwell_energy_density(:)
    FLOAT, pointer                 :: maxwell_energy_density_plane_waves(:)
    FLOAT, pointer                 :: maxwell_e_energy_density(:)
    FLOAT, pointer                 :: maxwell_b_energy_density(:)
    FLOAT                          :: maxwell_energy_trans
    FLOAT                          :: maxwell_energy_long
    FLOAT                          :: maxwell_e_energy_trans
    FLOAT                          :: maxwell_b_energy_trans
    FLOAT                          :: maxwell_energy_incident_waves

    logical                        :: maxwell_cpml_hamiltonian_mxll = .false.

    logical                        :: maxwell_diamag_current = .false.

    integer                        :: current_prop_test = 0

    CMPLX, pointer                 :: test_output(:,:)

    type(cube_t)                   :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map

  end type hamiltonian_mxll_t

  integer, public, parameter :: &
    LENGTH     = 1,             &
    VELOCITY   = 2

  integer, public, parameter ::        &
    INDEPENDENT_PARTICLES = 2, &
    HARTREE               = 1, &
    HARTREE_FOCK          = 3, &
    KOHN_SHAM_DFT         = 4, &
    CLASSICAL             = 5, &
    RDMFT                 = 7

  type(profile_t), save :: prof_hamiltonian_mxll, prof_kinetic_start, prof_kinetic_finish

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_null(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_null)
 
    hm%adjoint = .false.
    
    hm%maxwell_hamiltonian_mxll = .false.
    hm%maxwell_current_density_ext_flag = .false.

    nullify(hm%maxwell_energy_density)
    nullify(hm%maxwell_e_energy_density)
    nullify(hm%maxwell_b_energy_density)

    POP_SUB(hamiltonian_mxll_null)
  end subroutine hamiltonian_mxll_null


  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_init(hm, parser, gr, geo, st)
    type(hamiltonian_mxll_t),                   intent(inout) :: hm
    type(parser_t),                             intent(in)    :: parser
    type(grid_t),                       target, intent(inout) :: gr
    type(geometry_t),                   target, intent(inout) :: geo
    type(states_t),                     target, intent(inout) :: st

    integer :: iline, icol, default
    integer :: ncols
    type(block_t) :: blk
    type(profile_t), save :: prof
    character(len=1024) :: string

    PUSH_SUB(hamiltonian_mxll_init)

    call profiling_in(prof, 'HAMILTONIAN_INIT')
       
    call hamiltonian_mxll_null(hm)

    call states_dim_copy(hm%d, st%d)

    ASSERT(associated(gr%der%lapl))

    maxwell_hm%hm_base%maxwell_operator(1:MAX_DIM) => maxwell_gr%der%grad(1:MAX_DIM)     ! cross product for Maxwell calculation needs dimension >= 2

    maxwell_hm%maxwell_rs_sign = maxwell_st%maxwell_rs_sign

    SAFE_ALLOCATE(maxwell_hm%maxwell_vector_potential(1:maxwell_gr%mesh%np_part,1:maxwell_st%d%dim))
    SAFE_ALLOCATE(maxwell_hm%maxwell_energy_density(1:maxwell_gr%mesh%np_part))
    SAFE_ALLOCATE(maxwell_hm%maxwell_energy_density_plane_waves(1:maxwell_gr%mesh%np_part))
    SAFE_ALLOCATE(maxwell_hm%maxwell_e_energy_density(1:maxwell_gr%mesh%np_part))
    SAFE_ALLOCATE(maxwell_hm%maxwell_b_energy_density(1:maxwell_gr%mesh%np_part))

    maxwell_hm%maxwell_vector_potential = M_ZERO
    maxwell_hm%maxwell_energy_density = M_ZERO
    maxwell_hm%maxwell_e_energy_density = M_ZERO
    maxwell_hm%maxwell_b_energy_density = M_ZERO

    SAFE_ALLOCATE(maxwell_hm%cube%fft)

    !%Variable MaxwellHamiltonianOperator
    !%Type integer
    !%Default riemann_silberstein
    !%Section Hamiltonian
    !%Description
    !% With this variable the the Maxwell Hamiltonian operator can be selected
    !%Option faraday_ampere_old 0
    !% old version
    !%Option faraday_ampere 1
    !% The propagation operation in vacuum with Spin 1 matrices without Gauss law condition.
    !%Option faraday_ampere_medium 2
    !% The propagation operation in medium with Spin 1 matrices without Gauss law condition
    !%Option faraday_ampere_gauss 3
    !% The propagation operation is done by 4x4 matrices also with Gauss laws constraint. 
    !%Option faraday_ampere_gauss_medium 4
    !% The propagation operation is done by 4x4 matrices also with Gauss laws constraint in medium
    !%End
    call parse_variable(parser, 'MaxwellHamiltonianOperator', OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE, &
                        maxwell_hm%maxwell_operator)

    if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
      maxwell_hm%d%dim = maxwell_hm%d%dim+1
    else if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      maxwell_hm%d%dim = 2*maxwell_hm%d%dim
    end if

    !%Variable MaxwellExternalCurrent
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Description follows
    !%End
    call parse_variable(parser, 'MaxwellExternalCurrent', .false., maxwell_hm%maxwell_current_density_ext_flag)

    !%Variable MatterToMaxwellCoupling
    !%Type logical
    !%Default yes
    !%Section Hamiltonian
    !%Description
    !% If the MatterToMaxwellCoupling flag is set, the time-depending motion of the matter 
    !% electrons couples to the electromagnetic fields. In the case of switched off,
    !% the initial electromagnetic fields propagate undisturbed.
    !%End
    call parse_variable(parser, 'MatterToMaxwellCoupling', .true., maxwell_hm%ma_mx_coupling)
    maxwell_hm%ma_mx_coupling_apply = .false.

    maxwell_hm%maxwell_plane_waves_apply = .false.
    maxwell_hm%maxwell_spatial_constant_apply = .false.

    maxwell_hm%maxwell_propagation_apply = .false.

    !%Variable MaxwellMediumCalculation
    !%Type integer
    !%Default current
    !%Section Hamiltonian
    !%Description
    !% The Maxwell Operator e.g. the curl operation can be obtained by
    !% two different methods, the finid-difference or the fast fourier
    !% transform.
    !%Option riemann_silberstein 1
    !% Medium calculation directly via Hamiltonian
    !%Option electric_magnetic_fields 2
    !% Medium calculation via curl of electric field and magnetic field
    !%End
    default_propagator = OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN
    call parse_variable('MaxwellMediumCalculation', default_propagator, maxwell_hm%maxwell_medium_calculation)

    maxwell_hm%maxwell_rs_state_fft_map     => maxwell_st%maxwell_rs_state_fft_map
    maxwell_hm%maxwell_rs_state_fft_map_inv => maxwell_st%maxwell_rs_state_fft_map_inv



    if (hm%maxwell_ks_calculation) then

      call epot_init_maxwell_coupling(hm%ep, gr%mesh)

      !%Variable MaxwellTransFieldCalculationMethod
      !%Type integer
      !%Default 1
      !%Section Hamiltonian
      !%Description
      !% Follows 
      !%Option trans_field_matter 1
      !% Follows
      !%Option trans_field_poisson 3
      !% Follows
      !%Option trans_field_poisson_corr 4
      !% Follows
      !%Option trans_field_poisson_long 5
      !% Follows
      !%Option trans_field_poisson_long_corr 6
      !% Follows
      !%Option trans_field_poisson_no_corr 7
      !% Follows
      !%End
      default = OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON
      call parse_variable('MaxwellTransFieldCalculationMethod', default, hm%mx_ma_trans_field_calc_method)

      !%Variable MaxwellTransFieldCalculationCorr
      !%Type logical
      !%Default yes
      !%Section Hamiltonian
      !%Description
      !% Follows
      !%End
      call parse_variable('MaxwellTransFieldCalculationCorr', .true., hm%mx_ma_trans_field_calc_corr)

      call messages_print_stress(stdout, 'Method of transverse electric field calculation')      
      if (hm%mx_ma_trans_field_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_MATTER) then
        string = ' Transverse electric field calculation by subtracting longitudinal'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
        string = ' part of matter electric field.'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
      else if (hm%mx_ma_trans_field_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON) then
        string = ' Direct transverse electric field calculation with poisson solver'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
      else if (hm%mx_ma_trans_field_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_CORR) then
        string = ' Direct transverse electric field calculation with poisson solver'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
        string = ' after subtracting the longitudinal part of the matter electric field.'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
      else if (hm%mx_ma_trans_field_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_LONG) then
        string = ' Transverse electric field calculation by direct calculation'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
        string = ' of the longitudinal electric field and subtraction.'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
      else if (hm%mx_ma_trans_field_calc_method == OPTION__MAXWELLTRANSFIELDCALCULATIONMETHOD__TRANS_FIELD_POISSON_LONG_CORR) then
        string = ' Transverse electric field calculation by direct calculation'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
        string = ' of the longitudinal electric field after subtracting the'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
        string = ' longitudinal part of the matter electric field.'
        write(message(1),'(a)') trim(string)
        call messages_info(1)
      end if
      call messages_print_stress(stdout)

      !%Variable MaxwellLorentzForce
      !%Type logical
      !%Default yes
      !%Section Hamiltonian
      !%Description
      !% Follows...
      !%End
      call parse_variable('MaxwellLorentzForce', .true., hm%maxwell_lorentz_force)
      hm%maxwell_lorentz_force_apply = .false.

      !%Variable MaxwellDiamagneticCurrent
      !%Type logical
      !%Default yes
      !%Section Hamiltonian
      !%Description
      !% Follows...
      !%End
      call parse_variable('MaxwellDiamagneticCurrent', .true., hm%maxwell_diamagnetic_current)

      !%Variable MaxwellSpinCurrent
      !%Type logical
      !%Default yes
      !%Section Hamiltonian
      !%Description
      !% Follows...
      !%End
      call parse_variable('MaxwellSpinCurrent', .true., hm%maxwell_spin_current)

    end if
    
    if(associated(hm%ep%E_field) .and. simul_box_is_periodic(gr%sb) .and. .not. gauge_field_is_applied(hm%ep%gfield)) then
      ! only need vberry if there is a field in a periodic direction
      ! and we are not setting a gauge field
      if(any(abs(hm%ep%E_field(1:gr%sb%periodic_dim)) > M_EPSILON)) then
        SAFE_ALLOCATE(hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
      end if
    end if

    !Static magnetic field requires complex wavefunctions
    !Static magnetic field or rashba spin-orbit interaction requires complex wavefunctions
    if (associated(hm%ep%B_field)) call states_set_complex(st)

    call parse_variable('CalculateSelfInducedMagneticField', .false., hm%self_induced_magnetic)
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
    if(hm%self_induced_magnetic) then
      SAFE_ALLOCATE(hm%a_ind(1:gr%mesh%np_part, 1:gr%sb%dim))
      SAFE_ALLOCATE(hm%b_ind(1:gr%mesh%np_part, 1:gr%sb%dim))

      !(for dim = we could save some memory, but it is better to keep it simple)
    end if

    ! Boundaries
    call bc_init(hm%bc, gr%mesh, gr%mesh%sb, hm%geo)

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
    if(parse_block('MassScaling', blk) == 0) then
        ncols = parse_block_cols(blk, 0)
        if(ncols > gr%sb%dim) then
          call messages_input_error("MassScaling")
        end if
        iline = 1 ! just deal with 1 line - should be generalized
        do icol = 1, ncols
          call parse_block_float(blk, iline - 1, icol - 1, hm%mass_scaling(icol))
        end do
        call parse_block_end(blk)
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
    call parse_variable('StatesPack', .true., hm%apply_packed)

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
    call parse_variable('TimeZero', .false., hm%time_zero)
    if(hm%time_zero) call messages_experimental('TimeZero')

    call profiling_out(prof)
    POP_SUB(hamiltonian_mxll_init)

  end subroutine hamiltonian_mxll_init

  
  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_end(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_end)

    this%hm_base%maxwell_operator => null()

    SAFE_DEALLOCATE_P(this%maxwell_energy_density)
    SAFE_DEALLOCATE_P(this%cube%fft)

    call bc_end(hm%bc)

    call states_dim_end(hm%d) 

    POP_SUB(hamiltonian_mxll_end)
  end subroutine hamiltonian_mxll_end


  ! ---------------------------------------------------------
  logical function hamiltonian_mxll_hermitian(hm)
    type(hamiltonian_mxll_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_mxll_hermitian)

    hamiltonian_mxll_hermitian = .true.

    POP_SUB(hamiltonian_mxll_hermitian)
  end function hamiltonian_mxll_hermitian


  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_adjoint(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_adjoint)

    if(.not.hm%adjoint) then
      hm%adjoint = .true.
    end if

    POP_SUB(hamiltonian_mxll_adjoint)
  end subroutine hamiltonian_mxll_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_not_adjoint(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_not_adjoint)

    if(hm%adjoint) then
      hm%adjoint = .false.
    end if

    POP_SUB(hamiltonian_mxll_not_adjoint)
  end subroutine hamiltonian_mxll_not_adjoint


  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_update(this, mesh, time)
    type(hamiltonian_mxll_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    FLOAT, optional,     intent(in)    :: time

    PUSH_SUB(hamiltonian_mxll_update)

    this%current_time = M_ZERO
    if(present(time)) this%current_time = time

    POP_SUB(hamiltonian_mxll_update)

  end subroutine hamiltonian_mxll_update

  ! -----------------------------------------------------------------

  FLOAT function hamiltonian_mxll_get_time(this) result(time)
    type(hamiltonian_mxll_t),   intent(inout) :: this

    time = this%current_time

  end function hamiltonian_mxll_get_time

  ! -----------------------------------------------------------------

  logical pure function hamiltonian_mxll_apply_packed(this, mesh) result(apply)
    type(hamiltonian_mxll_t),   intent(in) :: this
    type(mesh_t),          intent(in) :: mesh

    apply = this%apply_packed
    if(mesh%use_curvilinear) apply = .false.
    
  end function hamiltonian_mxll_apply_packed

! ---------------------------------------------------------
subroutine hamiltonian_mxll_apply_batch (hm, der, psib, hpsib, time, terms, set_bc)
  type(hamiltonian_mxll_t),      intent(in)    :: hm
  type(derivatives_t),      intent(in)    :: der
  type(batch_t), target,    intent(inout) :: psib
  type(batch_t), target,    intent(inout) :: hpsib
  FLOAT, optional,          intent(in)    :: time
  integer, optional,        intent(in)    :: terms
  logical, optional,        intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.

  logical :: pack
  integer :: ii
  type(derivatives_handle_batch_t) :: handle
  integer :: terms_
  type(projection_t) :: projection

  call profiling_in(prof_hamiltonian_mxll, "MXLL_HAMILTONIAN")
  PUSH_SUB(hamiltonian_mxll_apply_batch)

  ASSERT(batch_status(psib) == batch_status(hpsib))

  ASSERT(batch_is_ok(psib))
  ASSERT(batch_is_ok(hpsib))
  ASSERT(psib%nst == hpsib%nst)

  !Not implemented at the moment
  ASSERT(.not.present(terms))
  ASSERT(.not.present(set_bc))

!  pack = hamiltonian_mxll_apply_packed(hm, der%mesh) &
!    .and. (accel_is_enabled() .or. psib%nst_linear > 1) &
!    .and. terms_ == TERM_ALL

!  if(pack) then
!    call batch_pack(psib)
!    call batch_pack(hpsib, copy = .false.)
!  end if

  if(present(time)) then
      if(abs(time - maxwell_hm%current_time) > CNST(1e-10)) then
        write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
        write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', maxwell_hm%current_time
        call messages_fatal(2)
      endif
    end if

  call zderivatives_curl(der, psib%states(1)%zpsi, hpsib%states(1)%zpsi)
  hpsib%states(1)%zpsi(:,:) = P_c * hpsib%states(1)%zpsi(:,:)
  
!  if(pack) then
!    call batch_unpack(psib, copy = .false.)
!    call batch_unpack(hpsib)
!  end if

  POP_SUB(hamiltonian_mxll_apply_batch)
  call profiling_out(prof_hamiltonian_mxll)

end subroutine hamiltonian_mxll_apply_batch

! ---------------------------------------------------------

subroutine hamiltonian_mxll_apply (hm, der, psi, hpsi, ist, ik, time, terms, set_bc)
  type(hamiltonian_mxll_t), intent(in)    :: hm
  type(derivatives_t), intent(in)    :: der
  integer,             intent(in)    :: ist       !< the index of the state
  integer,             intent(in)    :: ik        !< the index of the k-point
  R_TYPE,   target,    intent(inout) :: psi(:,:)  !< (gr%mesh%np_part, hm%d%dim)
  R_TYPE,   target,    intent(inout) :: hpsi(:,:) !< (gr%mesh%np, hm%d%dim)
  FLOAT,    optional,  intent(in)    :: time
  integer,  optional,  intent(in)    :: terms
  logical,  optional,  intent(in)    :: set_bc

  type(batch_t) :: psib, hpsib

  PUSH_SUB(hamiltonian_mxll_apply)

!  call batch_init(psib, hm%d%dim, 1)
!  call batch_add_state(psib, ist, psi)
!  call batch_init(hpsib, hm%d%dim, 1)
!  call batch_add_state(hpsib, ist, hpsi)
!
!  call hamiltonian_mxll_apply_batch(hm, der, psib, hpsib, ik, time = time, terms = terms, Imtime = Imtime, set_bc = set_bc)
!
!  call batch_end(psib)
!  call batch_end(hpsib)

   call profiling_in(prof_hamiltonian, "MAXWELLHAMILTONIAN")

   if(present(time)) then
      if(abs(time - maxwell_hm%current_time) > CNST(1e-10)) then
        write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
        write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', maxwell_hm%current_time
        call messages_fatal(2)
      endif
    end if

    select case (maxwell_hm%maxwell_op_method)
      case(OPTION__MAXWELLTDOPERATORMETHOD__MAXWELL_OP_FD)
        call maxwell_hamiltonian_apply_fd(maxwell_hm, maxwell_der, psi, oppsi)
     ! case(OPTION__MAXWELLTDOPERATORMETHOD__MAXWELL_OP_FFT)
     !   call maxwell_hamiltonian_apply_fft(maxwell_hm, maxwell_der, psi, oppsi)
    end select


  POP_SUB(hamiltonian_mxll_apply)
end subroutine hamiltonian_mxll_apply


! ---------------------------------------------------------
subroutine hamiltonian_mxll_apply_all (hm, xc, der, st, hst, time, Imtime)
  type(hamiltonian_mxll_t), intent(inout) :: hm
  type(xc_t),          intent(in)    :: xc
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  type(states_t),      intent(inout) :: hst
  FLOAT, optional,     intent(in)    :: time
  FLOAT, optional,     intent(in)    :: Imtime

  integer :: ik, ib, ist
  R_TYPE, allocatable :: psi(:, :)
  CMPLX,  allocatable :: psiall(:, :, :, :)
  
  PUSH_SUB(X(hamiltonian_mxll_apply_all))

  if(present(Imtime)) then
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call X(hamiltonian_mxll_apply_batch)(hm, der, st%group%psib(ib, ik), hst%group%psib(ib, ik), ik, time, Imtime)
      end do
    end do
  else
    do ik = st%d%kpt%start, st%d%kpt%end
      do ib = st%group%block_start, st%group%block_end
        call X(hamiltonian_mxll_apply_batch)(hm, der, st%group%psib(ib, ik), hst%group%psib(ib, ik), ik, time)
      end do
    end do
  end if

  POP_SUB(hamiltonian_mxll_apply_all)
end subroutine hamiltonian_mxll_apply_all

  ! ---------------------------------------------------------
  subroutine maxwell_hamiltonian_apply_fd(maxwell_hm, maxwell_der, psi, oppsi)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    CMPLX,               intent(inout) :: psi(:,:)
    CMPLX,               intent(inout) :: oppsi(:,:)

    FLOAT              :: kk(MAX_DIM), aa_e(MAX_DIM), aa_m(MAX_DIM), bb_e(MAX_DIM), bb_m(MAX_DIM)
    FLOAT              :: aux_ep(3), aux_mu(3), cc, pml_c(3), pml_aux_ep(3,3), pml_aux_mu(3,3), pml_g_e(3), pml_g_m(3)
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_curl_e(:,:), tmp_curl_b(:,:), tmp_partial_e(:,:), tmp_partial_b(:,:)
    FLOAT, pointer     :: mx_rho(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)
    CMPLX              :: ff_plus(MAX_DIM), ff_minus(MAX_DIM), sigma_e, sigma_m
    CMPLX, allocatable :: tmp(:,:)
    CMPLX, pointer     :: kappa_psi(:,:)
    integer            :: np, np_part, ip, ip_in, array_length_1, array_length_2, idim, il, ii, rs_sign

    PUSH_SUB(maxwell_hamiltonian_apply_fd)

    np = maxwell_der%mesh%np
    np_part = maxwell_der%mesh%np_part
    rs_sign = maxwell_hm%maxwell_rs_sign


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum via partial derivatives:

    if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE) then

      SAFE_ALLOCATE(tmp(np_part,2))
      oppsi       = M_z0

      if (maxwell_hm%maxwell_diamag_current) then
        mx_rho    => maxwell_hm%maxwell_st%maxwell_grid_rho
        kappa_psi => maxwell_hm%maxwell_st%maxwell_kappa_psi 
      end if

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,2), 3, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 2, 3, tmp(:,1))
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 3, 2, tmp(:,2))
      oppsi(1:np_part,1) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,2), 1, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 3, 1, tmp(:,1))
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 1, 3, tmp(:,2))
      oppsi(1:np_part,2) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,2), 2, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 1, 2, tmp(:,1))
      call maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, 2, 1, tmp(:,2))
      oppsi(1:np_part,3) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      if (maxwell_hm%maxwell_bc_constant) then
        do ip_in=1, maxwell_hm%maxwell_bc%constant_points_number
          ip = maxwell_hm%maxwell_bc%constant_points_map(ip_in)
          oppsi(ip,:) = maxwell_hm%maxwell_st%maxwell_rs_state_const(:)
        end do
      end if

!write(*,*) '1', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),1) !, oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),4)
!write(*,*) '2', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),2) !, oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),5) 
!write(*,*) '3', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),3) !, oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),6)
!write(*,*)

      SAFE_DEALLOCATE_A(tmp)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium via partial derivatives:

    else if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then

      SAFE_ALLOCATE(tmp(np_part,4))
      oppsi       = M_z0

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,3), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,6), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,5), tmp(:,4), 3, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 2, 3, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 3, 2, tmp(:,3:4))
      oppsi(1:np_part,1) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,4) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,2), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,6), tmp(:,4), 1, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 3, 1, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 1, 3, tmp(:,3:4))
      oppsi(1:np_part,2) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,5) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,3), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,4), 2, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 1, 2, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, 2, 1, tmp(:,3:4))
      oppsi(1:np_part,3) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,6) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )

!write(*,*) '1', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),1), oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),4)
!write(*,*) '2', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),2), oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),5) 
!write(*,*) '3', oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),3), oppsi(maxwell_der%mesh%idx%lxyz_inv(14,0,0),6)
!write(*,*)

      SAFE_DEALLOCATE_A(tmp)

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation if medium boundaries is set:
      call maxwell_medium_boundaries_calculation(maxwell_hm, maxwell_der, psi, oppsi)

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation for medium boxes:
      call maxwell_medium_boxes_calculation(maxwell_hm, maxwell_der, psi, oppsi)

      SAFE_DEALLOCATE_A(tmp)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum with Gauss condition via partial derivatives:

    else if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then

      SAFE_ALLOCATE(tmp(np_part,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,1) = rs_sign * P_c * ( - M_zI*tmp(1:np_part,1) - tmp(1:np_part,2) - M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,2) = rs_sign * P_c * ( - M_zI*tmp(1:np_part,1) - tmp(1:np_part,2) - M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,3) = rs_sign * P_c * ( tmp(1:np_part,1) - M_zI*tmp(1:np_part,2) + M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,4) = rs_sign * P_c * ( tmp(1:np_part,1) - M_zI*tmp(1:np_part,2) + M_zI*tmp(1:np_part,3) )

      SAFE_DEALLOCATE_A(tmp)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium with Gauss condition via partial derivatives:

    else if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS_MEDIUM) then

      SAFE_ALLOCATE(tmp(np_part,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,1) = P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,2) = P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,3) = P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,4) = P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,5), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,7), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,7), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,5) = - P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,6), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,8), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,8), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,6) = - P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,5), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,7) = - P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(maxwell_der, psi(:,6), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,6), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,8) = - P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      SAFE_DEALLOCATE_A(tmp)

    end if

    POP_SUB(maxwell_hamiltonian_apply_fd)
  end subroutine maxwell_hamiltonian_apply_fd





  ! ---------------------------------------------------------
  subroutine maxwell_pml_hamiltonian(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    CMPLX,               intent(inout) :: psi(:,:)
    integer,             intent(in)    :: dir1
    integer,             intent(in)    :: dir2
    CMPLX,               intent(inout) :: tmp(:)

    PUSH_SUB(maxwell_pml_hamiltonian)

    if ( (maxwell_hm%maxwell_bc%bc_ab_type(dir1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) .and. &
          maxwell_hm%maxwell_cpml_hamiltonian ) then
      if (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) then
        call maxwell_pml_calculation_via_riemann_silberstein(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp(:))
      else if (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) then
        call maxwell_pml_calculation_via_e_b_fields(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp(:))
      end if
    end if

    POP_SUB(maxwell_pml_hamiltonian)
  end subroutine maxwell_pml_hamiltonian


  ! ---------------------------------------------------------
  subroutine maxwell_pml_hamiltonian_medium(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    CMPLX,               intent(inout) :: psi(:,:)
    integer,             intent(in)    :: dir1
    integer,             intent(in)    :: dir2
    CMPLX,               intent(inout) :: tmp(:,:)

    PUSH_SUB(maxwell_pml_hamiltonian_medium)

    if ( (maxwell_hm%maxwell_bc%bc_ab_type(dir1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) .and. &
          maxwell_hm%maxwell_cpml_hamiltonian ) then
      if (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) then
        call maxwell_pml_calculation_via_riemann_silberstein_medium(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp(:,:))
      else if (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) then
        call maxwell_pml_calculation_via_e_b_fields_medium(maxwell_hm, maxwell_der, psi, dir1, dir2, tmp(:,:))
      end if
    end if

    POP_SUB(maxwell_pml_hamiltonian_medium)
  end subroutine maxwell_pml_hamiltonian_medium


  ! ---------------------------------------------------------
  subroutine maxwell_pml_calculation_via_riemann_silberstein(maxwell_hm, maxwell_der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    integer,             intent(in)    :: pml_dir
    CMPLX,               intent(inout) :: psi(:,:)
    integer,             intent(in)    :: field_dir
    CMPLX,               intent(inout) :: pml(:)

    integer            :: ip, ip_in, np_part, rs_sign
    FLOAT              :: pml_c(3)
    CMPLX, allocatable :: tmp_partial(:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein)

    if (maxwell_hm%maxwell_cpml_hamiltonian) then

      rs_sign = maxwell_hm%maxwell_rs_sign

      np_part = maxwell_der%mesh%np_part
      SAFE_ALLOCATE(tmp_partial(np_part))

      call zderivatives_partial(maxwell_der, psi(:,field_dir), tmp_partial(:), pml_dir, set_bc = .false.)
      do ip_in=1, maxwell_hm%maxwell_bc%pml_points_number
        ip       = maxwell_hm%maxwell_bc%pml_points_map(ip_in)
        pml_c(:) = maxwell_hm%maxwell_bc%pml_c(ip_in,:)
        pml_a(:) = maxwell_hm%maxwell_bc%pml_a(ip_in,:)
        pml_b(:) = maxwell_hm%maxwell_bc%pml_b(ip_in,:)
        pml_g(:) = maxwell_hm%maxwell_bc%pml_conv_plus(ip_in,pml_dir,:)
        pml(ip)  = rs_sign * pml_c(pml_dir) * tmp_partial(ip) &
                 + rs_sign * pml_c(pml_dir) * real(pml_a(pml_dir)) * real(tmp_partial(ip)) &
                 + rs_sign * M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip)) &
                 + rs_sign * pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g(field_dir)) &
                 + rs_sign * M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g(field_dir))
      end do

      SAFE_DEALLOCATE_A(tmp_partial)
    end if

    POP_SUB(maxwell_pml_calculation_via_riemann_silberstein)
  end subroutine maxwell_pml_calculation_via_riemann_silberstein


  ! ---------------------------------------------------------
  subroutine maxwell_pml_calculation_via_e_b_fields(maxwell_hm, maxwell_der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    integer,             intent(in)    :: pml_dir
    CMPLX,               intent(in)    :: psi(:,:)
    integer,             intent(in)    :: field_dir
    CMPLX,               intent(inout) :: pml(:)

    integer            :: ip, ip_in, np_part
    FLOAT              :: pml_c(3)
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_partial_e(:), tmp_partial_b(:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)

    PUSH_SUB(maxwell_pml_calculation_via_e_b_fields)

    if (maxwell_hm%maxwell_cpml_hamiltonian) then

      np_part = maxwell_der%mesh%np_part
      SAFE_ALLOCATE(tmp_e(np_part,3))
      SAFE_ALLOCATE(tmp_partial_e(np_part))
      SAFE_ALLOCATE(tmp_b(np_part,3))
      SAFE_ALLOCATE(tmp_partial_b(np_part))

      call maxwell_get_electric_field_state(psi, tmp_e)
      call maxwell_get_magnetic_field_state(psi, maxwell_hm%maxwell_rs_sign, tmp_b)
      call dderivatives_partial(maxwell_der, tmp_e(:,field_dir), tmp_partial_e(:), pml_dir, set_bc = .false.)
      call dderivatives_partial(maxwell_der, tmp_b(:,field_dir), tmp_partial_b(:), pml_dir, set_bc = .false.)
      tmp_partial_e(:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:)
      tmp_partial_b(:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:)
      do ip_in=1, maxwell_hm%maxwell_bc%pml_points_number
        ip       = maxwell_hm%maxwell_bc%pml_points_map(ip_in)
        pml_c(:) = maxwell_hm%maxwell_bc%pml_c(ip_in,:)
        pml_a(:) = maxwell_hm%maxwell_bc%pml_a(ip_in,:)
        pml_b(:) = maxwell_hm%maxwell_bc%pml_b(ip_in,:)
        pml_g(:) = maxwell_hm%maxwell_bc%pml_conv_plus(ip_in,pml_dir,:)
        pml(ip)  = pml_c(pml_dir) * tmp_partial_e(ip) + M_zI * pml_c(pml_dir) * tmp_partial_b(ip) &
                 + pml_c(pml_dir) * real(pml_a(pml_dir)) * tmp_partial_e(ip) &
                 + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * tmp_partial_b(ip) &
                 + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g(field_dir)) &
                 + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g(field_dir))
      end do

      SAFE_DEALLOCATE_A(tmp_e)
      SAFE_DEALLOCATE_A(tmp_partial_e)
      SAFE_DEALLOCATE_A(tmp_b)
      SAFE_DEALLOCATE_A(tmp_partial_b)

    end if

    POP_SUB(maxwell_pml_calculation_via_e_b_fields)
  end subroutine maxwell_pml_calculation_via_e_b_fields


  ! ---------------------------------------------------------
  subroutine maxwell_pml_calculation_via_riemann_silberstein_medium(maxwell_hm, maxwell_der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    integer,             intent(in)    :: pml_dir
    CMPLX,               intent(inout) :: psi(:,:)
    integer,             intent(in)    :: field_dir
    CMPLX,               intent(inout) :: pml(:,:)

    integer            :: ip, ip_in, np_part
    FLOAT              :: pml_c(3)
    CMPLX, allocatable :: tmp_partial(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g_p(3), pml_g_m(3)

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein_medium)

    if (maxwell_hm%maxwell_cpml_hamiltonian) then

      np_part = maxwell_der%mesh%np_part
      SAFE_ALLOCATE(tmp_partial(np_part,1:2))

      call zderivatives_partial(maxwell_der, psi(:,field_dir  ), tmp_partial(:,1), pml_dir, set_bc = .false.)
      call zderivatives_partial(maxwell_der, psi(:,field_dir+3), tmp_partial(:,2), pml_dir, set_bc = .false.)
      do ip_in=1, maxwell_hm%maxwell_bc%pml_points_number
        ip         = maxwell_hm%maxwell_bc%pml_points_map(ip_in)
        pml_c(:)   = maxwell_hm%maxwell_bc%pml_c(ip_in,:)
        pml_a(:)   = maxwell_hm%maxwell_bc%pml_a(ip_in,:)
        pml_b(:)   = maxwell_hm%maxwell_bc%pml_b(ip_in,:)
        pml_g_p(:) = maxwell_hm%maxwell_bc%pml_conv_plus(ip_in,pml_dir,:)
        pml_g_m(:) = maxwell_hm%maxwell_bc%pml_conv_minus(ip_in,pml_dir,:)
        pml(ip,1)  = pml_c(pml_dir) * tmp_partial(ip,1) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * real(tmp_partial(ip,1)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip,1)) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_p(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_p(field_dir))
        pml(ip,2)  = pml_c(pml_dir) * tmp_partial(ip,2) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * real(tmp_partial(ip,2)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip,2)) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_m(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_m(field_dir))
      end do

    end if

    SAFE_DEALLOCATE_A(tmp_partial)

    POP_SUB(maxwell_pml_calculation_via_riemann_silberstein_medium)
  end subroutine maxwell_pml_calculation_via_riemann_silberstein_medium


  ! ---------------------------------------------------------
  subroutine maxwell_pml_calculation_via_e_b_fields_medium(maxwell_hm, maxwell_der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    integer,             intent(in)    :: pml_dir
    CMPLX,               intent(in)    :: psi(:,:)
    integer,             intent(in)    :: field_dir
    CMPLX,               intent(inout) :: pml(:,:)

    integer            :: ip, ip_in, np_part
    FLOAT              :: pml_c(3)
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_partial_e(:,:), tmp_partial_b(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g_p(3), pml_g_m(3)

    PUSH_SUB(maxwell_pml_calculation_via_e_b_fields_medium)

    if (maxwell_hm%maxwell_cpml_hamiltonian) then

      np_part = maxwell_der%mesh%np_part
      SAFE_ALLOCATE(tmp_e(np_part,3))
      SAFE_ALLOCATE(tmp_partial_e(np_part,2))
      SAFE_ALLOCATE(tmp_b(np_part,3))
      SAFE_ALLOCATE(tmp_partial_b(np_part,2))

      call maxwell_get_electric_field_state(psi(:,1:3), tmp_e)
      call maxwell_get_magnetic_field_state(psi(:,1:3), 1, tmp_b)
      call dderivatives_partial(maxwell_der, tmp_e(:,field_dir), tmp_partial_e(:,1), pml_dir, set_bc = .false.)
      call dderivatives_partial(maxwell_der, tmp_b(:,field_dir), tmp_partial_b(:,1), pml_dir, set_bc = .false.)
      call maxwell_get_electric_field_state(psi(:,4:6), tmp_e)
      call maxwell_get_magnetic_field_state(psi(:,4:6), -1, tmp_b)
      call dderivatives_partial(maxwell_der, tmp_e(:,field_dir), tmp_partial_e(:,2), pml_dir, set_bc = .false.)
      call dderivatives_partial(maxwell_der, tmp_b(:,field_dir), tmp_partial_b(:,2), pml_dir, set_bc = .false.)
      tmp_partial_e(:,:) = sqrt(P_ep/M_TWO) * tmp_partial_e(:,:)
      tmp_partial_b(:,:) = sqrt(M_ONE/(M_TWO*P_mu)) * tmp_partial_b(:,:)
      do ip_in=1, maxwell_hm%maxwell_bc%pml_points_number
        ip         = maxwell_hm%maxwell_bc%pml_points_map(ip_in)
        pml_c(:)   = maxwell_hm%maxwell_bc%pml_c(ip_in,:)
        pml_a(:)   = maxwell_hm%maxwell_bc%pml_a(ip_in,:)
        pml_b(:)   = maxwell_hm%maxwell_bc%pml_b(ip_in,:)
        pml_g_p(:) = maxwell_hm%maxwell_bc%pml_conv_plus(ip_in,pml_dir,:)
        pml_g_m(:) = maxwell_hm%maxwell_bc%pml_conv_minus(ip_in,pml_dir,:)
        pml(ip,1)  = pml_c(pml_dir) * tmp_partial_e(ip,1) + M_zI * pml_c(pml_dir) * tmp_partial_b(ip,1) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * tmp_partial_e(ip,1) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * tmp_partial_b(ip,1) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_p(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_p(field_dir))
        pml(ip,2)  = pml_c(pml_dir) * tmp_partial_e(ip,2) + M_zI * pml_c(pml_dir) * tmp_partial_b(ip,2) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * tmp_partial_e(ip,2) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * tmp_partial_b(ip,2) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_m(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_m(field_dir))
      end do

      SAFE_DEALLOCATE_A(tmp_e)
      SAFE_DEALLOCATE_A(tmp_partial_e)
      SAFE_DEALLOCATE_A(tmp_b)
      SAFE_DEALLOCATE_A(tmp_partial_b)

    end if

    POP_SUB(maxwell_pml_calculation_via_e_b_fields_medium)
  end subroutine maxwell_pml_calculation_via_e_b_fields_medium


  ! ---------------------------------------------------------
  subroutine maxwell_medium_boundaries_calculation(maxwell_hm, maxwell_der, psi, oppsi)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    CMPLX,               intent(in)    :: psi(:,:)
    CMPLX,               intent(inout) :: oppsi(:,:)

    integer            :: ip, ip_in, idim
    FLOAT              :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m
    CMPLX              :: ff_plus(3), ff_minus(3)

    PUSH_SUB(maxwell_medium_boundaries_calculation)

    do idim=1, 3
      if ( (maxwell_hm%maxwell_bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) .and. &
           (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) ) then
        do ip_in=1, maxwell_hm%maxwell_bc%medium_points_number(idim)
          ip          = maxwell_hm%maxwell_bc%medium_points_map(ip_in,idim)
          cc          = maxwell_hm%maxwell_bc%medium_c(ip_in,idim)/P_c
          aux_ep(:)   = maxwell_hm%maxwell_bc%medium_aux_ep(ip_in,:,idim)
          aux_mu(:)   = maxwell_hm%maxwell_bc%medium_aux_mu(ip_in,:,idim)
          sigma_e     = maxwell_hm%maxwell_bc%medium_sigma_e(ip_in,idim)
          sigma_m     = maxwell_hm%maxwell_bc%medium_sigma_m(ip_in,idim)
          ff_plus(1)  = psi(ip,1)
          ff_plus(2)  = psi(ip,2)
          ff_plus(3)  = psi(ip,3)
          ff_minus(1) = psi(ip,4)
          ff_minus(2) = psi(ip,5)
          ff_minus(3) = psi(ip,6)
          aux_ep      = dcross_product(aux_ep,real(ff_plus+ff_minus))
          aux_mu      = dcross_product(aux_mu,aimag(ff_plus-ff_minus))
          oppsi(ip,1) = oppsi(ip,1)*cc                                        &
                      - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,4) = oppsi(ip,4)*cc                                        &
                      + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,2) = oppsi(ip,2)*cc                                        &
                      - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2))
          oppsi(ip,5) = oppsi(ip,5)*cc                                        &
                      + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2)) 
          oppsi(ip,3) = oppsi(ip,3)*cc                                        &
                      - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
          oppsi(ip,6) = oppsi(ip,6)*cc                                        &
                      + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
        end do
      end if
    end do

    POP_SUB(maxwell_medium_boundaries_calculation)
  end subroutine maxwell_medium_boundaries_calculation


  ! ---------------------------------------------------------
  subroutine maxwell_medium_boxes_calculation(maxwell_hm, maxwell_der, psi, oppsi)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(derivatives_t), intent(in)    :: maxwell_der
    CMPLX,               intent(in)    :: psi(:,:)
    CMPLX,               intent(inout) :: oppsi(:,:)

    integer            :: ip, ip_in, idim, il, np_part
    FLOAT              :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_curl_e(:,:), tmp_curl_b(:,:)
    CMPLX              :: ff_plus(3), ff_minus(3)

    PUSH_SUB(maxwell_medium_boxes_calculation)

    np_part = maxwell_der%mesh%np_part

    if ( maxwell_hm%maxwell_medium_box .and. &
         (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) ) then
      do il=1, maxwell_hm%maxwell_medium_box_number
        do ip_in=1, maxwell_hm%maxwell_medium_box_points_number(il)
          ip           = maxwell_hm%maxwell_medium_box_points_map(ip_in,il)
          cc           = maxwell_hm%maxwell_medium_box_c(ip_in,il)/P_c
          aux_ep(:)    = maxwell_hm%maxwell_medium_box_aux_ep(ip_in,:,il)
          aux_mu(:)    = maxwell_hm%maxwell_medium_box_aux_mu(ip_in,:,il)
          sigma_e      = maxwell_hm%maxwell_medium_box_sigma_e(ip_in,il)
          sigma_m      = maxwell_hm%maxwell_medium_box_sigma_m(ip_in,il)
          ff_plus(1)   = psi(ip,1)
          ff_plus(2)   = psi(ip,2)
          ff_plus(3)   = psi(ip,3)
          ff_minus(1)  = psi(ip,4)
          ff_minus(2)  = psi(ip,5)
          ff_minus(3)  = psi(ip,6)
          aux_ep       = dcross_product(aux_ep,real(ff_plus+ff_minus))
          aux_mu       = dcross_product(aux_mu,aimag(ff_plus-ff_minus))
          oppsi(ip,1) = oppsi(ip,1)*cc                                        &
                      - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,4) = oppsi(ip,4)*cc                                        &
                      + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,2) = oppsi(ip,2)*cc                                        &
                      - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2))
          oppsi(ip,5) = oppsi(ip,5)*cc                                        &
                      + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2)) 
          oppsi(ip,3) = oppsi(ip,3)*cc                                        &
                      - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))         &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
          oppsi(ip,6) = oppsi(ip,6)*cc                                        &
                      + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))         &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
        end do
      end do
    end if

    if ( maxwell_hm%maxwell_medium_box .and. &
         (maxwell_hm%maxwell_medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) ) then
      SAFE_ALLOCATE(tmp_e(np_part,3))
      SAFE_ALLOCATE(tmp_curl_e(np_part,3))
      SAFE_ALLOCATE(tmp_b(np_part,3))
      SAFE_ALLOCATE(tmp_curl_b(np_part,3))
      call maxwell_get_electric_field_state(psi(:,1:3)+psi(:,4:6), tmp_e, maxwell_hm%maxwell_st%maxwell_ep, np_part)
      call maxwell_get_magnetic_field_state(psi(:,1:3)+psi(:,4:6), maxwell_hm%maxwell_rs_sign, tmp_b, &
                                            maxwell_hm%maxwell_st%maxwell_mu, np_part)
      call dderivatives_curl(maxwell_der, tmp_e, tmp_curl_e, set_bc = .false.)
      call dderivatives_curl(maxwell_der, tmp_b, tmp_curl_b, set_bc = .false.)
      SAFE_DEALLOCATE_A(tmp_e)
      SAFE_DEALLOCATE_A(tmp_b)
      tmp_curl_e(:,1) = sqrt(maxwell_hm%maxwell_st%maxwell_ep(:)/M_TWO) * tmp_curl_e(:,1)
      tmp_curl_e(:,2) = sqrt(maxwell_hm%maxwell_st%maxwell_ep(:)/M_TWO) * tmp_curl_e(:,2)
      tmp_curl_e(:,3) = sqrt(maxwell_hm%maxwell_st%maxwell_ep(:)/M_TWO) * tmp_curl_e(:,3)
      tmp_curl_b(:,1) = sqrt(M_ONE/(M_TWO*maxwell_hm%maxwell_st%maxwell_mu(:))) * tmp_curl_b(:,1)
      tmp_curl_b(:,2) = sqrt(M_ONE/(M_TWO*maxwell_hm%maxwell_st%maxwell_mu(:))) * tmp_curl_b(:,2)
      tmp_curl_b(:,3) = sqrt(M_ONE/(M_TWO*maxwell_hm%maxwell_st%maxwell_mu(:))) * tmp_curl_b(:,3)
      do il=1, maxwell_hm%maxwell_medium_box_number
        do ip_in=1, maxwell_hm%maxwell_medium_box_points_number(il)
          ip           = maxwell_hm%maxwell_medium_box_points_map(ip_in,il)
          cc           = maxwell_hm%maxwell_medium_box_c(ip_in,il)
          sigma_e      = maxwell_hm%maxwell_medium_box_sigma_e(ip_in,il)
          sigma_m      = maxwell_hm%maxwell_medium_box_sigma_m(ip_in,il)
          ff_plus(1)   = psi(ip,1)
          ff_plus(2)   = psi(ip,2)
          ff_plus(3)   = psi(ip,3)
          ff_minus(1)  = psi(ip,4)
          ff_minus(2)  = psi(ip,5)
          ff_minus(3)  = psi(ip,6)
          oppsi(ip,1) =   cc * tmp_curl_e(ip,1) + cc * M_zI * tmp_curl_b(ip,1)   &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))            &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,4) = - cc * tmp_curl_e(ip,1) + cc * M_zI * tmp_curl_b(ip,1)   &
                      - M_zI * sigma_e * real(ff_plus(1)+ff_minus(1))            &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(1)-ff_minus(1))
          oppsi(ip,2) =   cc * tmp_curl_e(ip,2) + cc * M_zI * tmp_curl_b(ip,2)   &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))            &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2))
          oppsi(ip,5) = - cc * tmp_curl_e(ip,2) + cc * M_zI * tmp_curl_b(ip,2)   &
                      - M_zI * sigma_e * real(ff_plus(2)+ff_minus(2))            &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(2)-ff_minus(2)) 
          oppsi(ip,3) =   cc * tmp_curl_e(ip,3) + cc * M_zI * tmp_curl_b(ip,3)   &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))            &
                      - M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
          oppsi(ip,6) = - cc * tmp_curl_e(ip,3) + cc * M_zI * tmp_curl_b(ip,3)   &
                      - M_zI * sigma_e * real(ff_plus(3)+ff_minus(3))            &
                      + M_zI * sigma_m * M_zI * aimag(ff_plus(3)-ff_minus(3))
        end do
      end do
    end if

    POP_SUB(maxwell_medium_boxes_calculation)
  end subroutine maxwell_medium_boxes_calculation


  ! ---------------------------------------------------------
  subroutine maxwell_fft_hamiltonian(maxwell_hm, k_vec, ff_dim, k_index_x, k_index_y, k_index_z, sigma, hm_matrix)
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    FLOAT,               intent(in)    :: k_vec(:,:)
    integer,             intent(in)    :: ff_dim
    integer,             intent(in)    :: k_index_x
    integer,             intent(in)    :: k_index_y
    integer,             intent(in)    :: k_index_z
    FLOAT,               intent(in)    :: sigma(:)
    CMPLX,               intent(inout) :: hm_matrix(:,:)

    FLOAT :: omega
    CMPLX :: ss(3)

    PUSH_SUB(maxwell_fft_hamiltonian)

    if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE) then
      hm_matrix(:,:) = M_z0
      if (maxwell_hm%maxwell_bc%bc_ab_type(1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then   ! TODO different directions
        omega = P_c * sqrt(k_vec(k_index_x,1)**2+k_vec(k_index_y,2)**2+k_vec(k_index_z,3)**2)
        if (omega /= M_ZERO) then
          ss(:) = M_ONE + sigma(:)/(M_zI*P_ep*omega)
        else
          ss(:) = M_ONE
        end if
        hm_matrix(1,1)   =   M_z0
        hm_matrix(1,2)   = - M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(1,3)   =   M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(2,1)   =   M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(2,2)   =   M_z0
        hm_matrix(2,3)   = - M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(3,1)   = - M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(3,2)   =   M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(3,3)   =   M_z0
      else
        hm_matrix(1,1)   =   M_z0
        hm_matrix(1,2)   = - M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(1,3)   =   M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(2,1)   =   M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(2,2)   =   M_z0
        hm_matrix(2,3)   = - M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(3,1)   = - M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(3,2)   =   M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(3,3)   =   M_z0
      end if
    else if (maxwell_hm%maxwell_operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      if (maxwell_hm%maxwell_bc%bc_ab_type(1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) then   ! TODO different directions
        omega = P_c * sqrt(k_vec(k_index_x,1)**2+k_vec(k_index_y,2)**2+k_vec(k_index_z,3)**2)
        if (omega /= M_ZERO) then
          ss(:) = M_ONE + sigma(:)/(M_zI*P_ep*omega)
        else
          ss(:) = M_ONE
        end if
        hm_matrix(1,1)   =   M_z0
        hm_matrix(1,2)   = - M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(1,3)   =   M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(2,1)   =   M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(2,2)   =   M_z0
        hm_matrix(2,3)   = - M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(3,1)   = - M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(3,2)   =   M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(3,3)   =   M_z0
        hm_matrix(4,4)   =   M_z0
        hm_matrix(4,5)   =   M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(4,6)   = - M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(5,4)   = - M_zI * P_c * k_vec(k_index_z,3)/ss(3)
        hm_matrix(5,5)   =   M_z0
        hm_matrix(5,6)   =   M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(6,4)   =   M_zI * P_c * k_vec(k_index_y,2)/ss(2)
        hm_matrix(6,5)   = - M_zI * P_c * k_vec(k_index_x,1)/ss(1)
        hm_matrix(6,6)   =   M_z0
      else
        hm_matrix(1,1)   =   M_z0
        hm_matrix(1,2)   = - M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(1,3)   =   M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(2,1)   =   M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(2,2)   =   M_z0
        hm_matrix(2,3)   = - M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(3,1)   = - M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(3,2)   =   M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(3,3)   =   M_z0
        hm_matrix(4,4)   =   M_z0
        hm_matrix(4,5)   =   M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(4,6)   = - M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(5,4)   = - M_zI * P_c * k_vec(k_index_z,3)
        hm_matrix(5,5)   =   M_z0
        hm_matrix(5,6)   =   M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(6,4)   =   M_zI * P_c * k_vec(k_index_y,2)
        hm_matrix(6,5)   = - M_zI * P_c * k_vec(k_index_x,1)
        hm_matrix(6,6)   =   M_z0
      end if
    end if

    POP_SUB(maxwell_fft_hamiltonian)
  end subroutine maxwell_fft_hamiltonian

  !----------------------------------------------------------
  subroutine get_electric_field_at_points(rs_state, maxwell_gr, maxwell_st, points_number, points, e_field_points)
    CMPLX,               intent(inout) :: rs_state(:,:)
    type(grid_t),        intent(in)    :: maxwell_gr
    type(states_t),      intent(in)    :: maxwell_st
    integer,             intent(in)    :: points_number
    FLOAT,               intent(in)    :: points(:,:)
    FLOAT,               intent(inout) :: e_field_points(:,:)

    integer              :: ig, ip, idim, rankmin_mx, partmin_mx
    FLOAT                :: dmin_mx
    integer, allocatable :: ip_mx_gr_local(:), ip_mx_gr_global(:)
    CMPLX,   allocatable :: ztmp_global(:)

    PUSH_SUB(get_electric_field_at_points)

    ! Fist part: get transverse electric field an evaluation point
    SAFE_ALLOCATE(ztmp_global(maxwell_gr%mesh%np_global))
    SAFE_ALLOCATE(ip_mx_gr_local(1:points_number))
    SAFE_ALLOCATE(ip_mx_gr_global(1:points_number))

    do ig=1, points_number
      ip_mx_gr_local(ig)  = mesh_nearest_point(maxwell_gr%mesh, points(:,ig), dmin_mx, rankmin_mx)
      partmin_mx          = rankmin_mx + 1
      ip_mx_gr_global(ig) = vec_local2global(maxwell_gr%mesh%vp, ip_mx_gr_local(ig), partmin_mx+1)
    end do

    do ig=1, points_number
      do idim=1, maxwell_st%d%dim
        if (maxwell_gr%mesh%parallel_in_domains) then
#ifdef HAVE_MPI
          call vec_allgather(maxwell_gr%mesh%vp, ztmp_global, rs_state(:,idim))
     !     call MPI_Barrier(maxwell_st%mpi_grp%comm, mpi_err)
     !     call MPI_Bcast(maxwell_st%maxwell_rs_state_trans(ip_mx_gr_local,idim), 1, MPI_CMPLX, &
     !                    maxwell_gr%mesh%vp%rank, maxwell_st%mpi_grp%comm-1, mpi_err)
#endif
        else
          ztmp_global = rs_state(:,idim)
        end if
        e_field_points(idim,ig) = sqrt(CNST(8.0)*M_PI) * real( ztmp_global(ip_mx_gr_global(ig)) )
      end do
    end do

    SAFE_DEALLOCATE_A(ztmp_global)
    SAFE_DEALLOCATE_A(ip_mx_gr_local)
    SAFE_DEALLOCATE_A(ip_mx_gr_global)

    POP_SUB(get_electric_field_at_points)

  end subroutine get_electric_field_at_points


  !----------------------------------------------------------
  subroutine get_magnetic_field_at_points(rs_state, maxwell_gr, maxwell_st, points_number, points, b_field_points)
    CMPLX,               intent(inout) :: rs_state(:,:)
    type(grid_t),        intent(in)    :: maxwell_gr
    type(states_t),      intent(in)    :: maxwell_st
    integer,             intent(in)    :: points_number
    FLOAT,               intent(in)    :: points(:,:)
    FLOAT,               intent(inout) :: b_field_points(:,:)

    integer :: ig, ip, idim, rankmin_mx, partmin_mx
    FLOAT   :: dmin_mx
    integer, allocatable :: ip_mx_gr_local(:), ip_mx_gr_global(:)
    CMPLX,   allocatable :: ztmp_global(:)

    PUSH_SUB(get_magnetic_field_at_points)

    ! Fist part: get transverse electric field an evaluation point
    SAFE_ALLOCATE(ztmp_global(maxwell_gr%mesh%np_global))
    SAFE_ALLOCATE(ip_mx_gr_local(points_number))
    SAFE_ALLOCATE(ip_mx_gr_global(points_number))

    do ig=1, points_number
      ip_mx_gr_local(ig)    = mesh_nearest_point(maxwell_gr%mesh, points(:,ig), dmin_mx, rankmin_mx)
      partmin_mx            = rankmin_mx + 1
      ip_mx_gr_global(ig)   = vec_local2global(maxwell_gr%mesh%vp, ip_mx_gr_local(ig), partmin_mx+1)
    end do

    do ig=1, points_number
      do idim=1, maxwell_st%d%dim
        if (maxwell_gr%mesh%parallel_in_domains) then
#ifdef HAVE_MPI
          call vec_allgather(maxwell_gr%mesh%vp, ztmp_global, rs_state(:,idim))
     !     call MPI_Barrier(maxwell_st%mpi_grp%comm, mpi_err)
     !     call MPI_Bcast(maxwell_st%maxwell_rs_state_trans(ip_mx_gr_local,idim), 1, MPI_CMPLX, &
     !                    maxwell_gr%mesh%vp%rank, maxwell_st%mpi_grp%comm-1, mpi_err)
#endif
        else
          ztmp_global = rs_state(:,idim)
        end if
        b_field_points(idim,ig) = sqrt(CNST(8.0)*M_PI)/P_c * real( ztmp_global(ip_mx_gr_global(ig)) )
      end do
    end do

    SAFE_DEALLOCATE_A(ztmp_global)
    SAFE_DEALLOCATE_A(ip_mx_gr_local)
    SAFE_DEALLOCATE_A(ip_mx_gr_global)

    POP_SUB(get_magnetic_field_at_points)

  end subroutine get_magnetic_field_at_points


  !----------------------------------------------------------
  subroutine maxwell_helmholtz_decomposition_trans_field(poisson, gr, maxwell_hm, hm, maxwell_st, transverse_field)
    type(poisson_t),     intent(in)    :: poisson
    type(grid_t),        intent(in)    :: gr
    type(hamiltonian_t), intent(in)    :: maxwell_hm
    type(hamiltonian_t), intent(in)    :: hm
    type(states_t),      intent(in)    :: maxwell_st
    CMPLX,               intent(inout) :: transverse_field(:,:)

    integer            :: idim, rankmin, ip_local, ip_global, ii, jj, kk, ip, ip_in
    FLOAT              :: pos(3), dmin
    CMPLX              :: surface_integral(3)
    CMPLX, allocatable :: field_old(:,:), ztmp(:,:), tmp_poisson(:)

    PUSH_SUB(maxwell_helmholtz_decomposition_trans_field)

    SAFE_ALLOCATE(field_old(1:gr%mesh%np_part,1:3))
    SAFE_ALLOCATE(ztmp(1:gr%mesh%np_part,1:3))
    SAFE_ALLOCATE(tmp_poisson(1:gr%mesh%np_part))

    field_old   = M_z0
    ztmp        = M_z0
    tmp_poisson = M_z0

    field_old = transverse_field
    call zderivatives_curl(gr%der, transverse_field(:,:), ztmp(:,:), set_bc = .false.)
    ! Apply poisson equation to solve helmholtz decomposition integral
    do idim=1, maxwell_st%d%dim
      call zpoisson_solve(poisson, tmp_poisson(:), ztmp(:,idim), .true.)
      ztmp(1:gr%mesh%np_part,idim) = M_ONE / (M_FOUR * M_PI) * tmp_poisson(1:gr%mesh%np_part)
    end do
    !do ip_in=1, maxwell_hm%maxwell_bc%der_bndry_mask_points_number
    !  ip = maxwell_hm%maxwell_bc%der_bndry_mask_points_map(ip_in)
    !  ztmp(ip,:) = maxwell_hm%maxwell_bc%der_bndry_mask(ip_in) * ztmp(ip,:)
    !end do
    call zderivatives_curl(gr%der, ztmp, transverse_field, set_bc = .false.)
    do ip_in=1, maxwell_hm%maxwell_bc%der_bndry_mask_points_number
      ip = maxwell_hm%maxwell_bc%der_bndry_mask_points_map(ip_in)
      transverse_field(ip,:) = maxwell_hm%maxwell_bc%der_bndry_mask(ip_in) * transverse_field(ip,:)
    end do

    SAFE_DEALLOCATE_A(ztmp)
    SAFE_DEALLOCATE_A(tmp_poisson)

    ! correction surface integral
    if (hm%mx_ma_trans_field_calc_corr) then
      call mesh_nearest_point_infos(gr%mesh, hm%mx_ma_coupling_points(:,1), dmin, rankmin, &
                                    ip_local, ip_global)
      ip_local  = 1
      ip_global = 1
      do ii = -gr%der%order, gr%der%order
        do jj = -gr%der%order, gr%der%order
          do kk = -gr%der%order, gr%der%order
            pos(1) = gr%mesh%x(ip_local,1) * ii * gr%mesh%spacing(1)
            pos(2) = gr%mesh%x(ip_local,2) * jj * gr%mesh%spacing(2)
            pos(3) = gr%mesh%x(ip_local,3) * kk * gr%mesh%spacing(3)
            call surface_integral_helmholtz_transverse(gr, maxwell_st, pos, field_old, surface_integral)
            transverse_field(ip_local,:) = transverse_field(ip_local,:) - M_ONE / (M_FOUR * M_PI) * surface_integral(:)
          end do
        end do
      end do
      pos(:) = maxwell_gr%mesh%x(ip_local,:)
    end if

    !do ip=1, maxwell_gr%mesh%np
    !  pos(:) = maxwell_gr%mesh%x(ip,:)
    !  call surface_integral_helmholtz_transverse(maxwell_gr, maxwell_st, pos, field_old, surface_integral)
    !  transverse_field(ip,:) = transverse_field(ip,:) - M_ONE / (M_FOUR * M_PI) * surface_integral(:)
    !end do

    SAFE_DEALLOCATE_A(field_old)

    POP_SUB(maxwell_helmholtz_decomposition_trans_field)
  end subroutine maxwell_helmholtz_decomposition_trans_field


  !----------------------------------------------------------
  subroutine maxwell_helmholtz_decomposition_long_field(poisson, gr, longitudinal_field)
    type(poisson_t),     intent(in)    :: poisson
    type(grid_t),        intent(in)    :: gr
    CMPLX,               intent(inout) :: longitudinal_field(:,:)

    CMPLX, allocatable :: ztmp(:), tmp_poisson(:)        

    SAFE_ALLOCATE(ztmp(1:gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_poisson(1:gr%mesh%np_part))

    ztmp        = M_z0
    tmp_poisson = M_z0

    call zderivatives_div(gr%der, longitudinal_field(:,:), ztmp(:), set_bc = .false.)
    ! Apply poisson equation to solve helmholtz decomposition integral
    call zpoisson_solve(poisson, tmp_poisson(:), ztmp(:), .true.)
    ztmp(1:gr%mesh%np_part) = - M_ONE / (M_FOUR * M_PI) * tmp_poisson(1:gr%mesh%np_part)
    !do ip_in=1, maxwell_hm%maxwell_bc%der_bndry_mask_points_number
    !  ip = maxwell_hm%maxwell_bc%der_bndry_mask_points_map(ip_in)
    !  ztmp(ip) = maxwell_hm%maxwell_bc%der_bndry_mask(ip_in) * ztmp(ip)
    !end do
    call zderivatives_grad(gr%der, ztmp(:), longitudinal_field(:,:), set_bc = .false.)

    !do ip_in=1, maxwell_hm%maxwell_bc%der_bndry_mask_points_number
    !  ip = maxwell_hm%maxwell_bc%der_bndry_mask_points_map(ip_in)
    !  longitudinal_field(ip,:) = maxwell_hm%maxwell_bc%der_bndry_mask(ip_in) * field_long(ip,:)
    !end do

    SAFE_DEALLOCATE_A(ztmp)
    SAFE_DEALLOCATE_A(tmp_poisson)

  end subroutine maxwell_helmholtz_decomposition_long_field


  ! ---------------------------------------------------------
  subroutine surface_integral_helmholtz_transverse(gr, st, pos, field, surface_integral)
    type(grid_t),      intent(in)    :: gr
    type(states_t),    intent(in)    :: st
    FLOAT,             intent(in)    :: pos(:) 
    CMPLX,             intent(in)    :: field(:,:)
    CMPLX,             intent(inout) :: surface_integral(:)

    integer             :: idim, idir, ip_surf, ip_global, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: xx(3)
    CMPLX               :: tmp_sum(3), normal(3)
    CMPLX,  allocatable :: tmp_global(:,:), tmp_surf(:,:,:,:,:)

    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_part_global,1:maxwell_st%d%dim))

    if (gr%mesh%parallel_in_domains) then
      do idim=1, maxwell_st%d%dim
#if defined(HAVE_MPI)
        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), field(:,idim))
        call MPI_Barrieri(gr%mesh%vp%comm, mpi_err)
#endif
      end do
    else
      tmp_global(:,:) = field(:,:)
    end if

    ix_max = st%maxwell_surface_grid_rows_number(1)
    iy_max = st%maxwell_surface_grid_rows_number(2)
    iz_max = st%maxwell_surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_z0
    tmp_sum  = M_z0

    do iy=1, iy_max
      do iz=1, iz_max
        do ip_surf=1, st%maxwell_surface_grid_points_number(1,iy,iz)
          normal    =  M_z0
          normal(1) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(1,1,iy,iz,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) & 
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(1,1,iy,iz,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
          normal    =  M_z0
          normal(1) = +M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(2,1,iy,iz,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) &
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(2,1,iy,iz,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
        end do
        tmp_surf(1,1,iy,iz,:) = tmp_surf(1,1,iy,iz,:) / float(st%maxwell_surface_grid_points_number(1,iy,iz))
        tmp_surf(2,1,iy,iz,:) = tmp_surf(2,1,iy,iz,:) / float(st%maxwell_surface_grid_points_number(1,iy,iz))
      end do
    end do
    do iy=1, iy_max
      do iz=1, iz_max
        tmp_sum(:) = tmp_sum(:) + tmp_surf(1,1,iy,iz,:) * st%maxwell_surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2,1,iy,iz,:) * st%maxwell_surface_grid_element(:)
      end do
    end do

    do ix=1, ix_max
      do iz=1, iz_max
        do ip_surf=1, st%maxwell_surface_grid_points_number(2,ix,iz)
          normal    =  M_z0
          normal(2) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(1,2,ix,iz,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) &
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(1,2,ix,iz,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
          normal    =  M_z0
          normal(2) =  M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(2,2,ix,iz,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) &
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(2,2,ix,iz,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
        end do
        tmp_surf(1,2,ix,iz,:) = tmp_surf(1,2,ix,iz,:) / float(st%maxwell_surface_grid_points_number(2,ix,iz))
        tmp_surf(2,2,ix,iz,:) = tmp_surf(2,2,ix,iz,:) / float(st%maxwell_surface_grid_points_number(2,ix,iz))
      end do
    end do
    do ix=1, ix_max
      do iz=1, iz_max
        tmp_sum(:) = tmp_sum(:) + tmp_surf(1,2,ix,iz,:) * st%maxwell_surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2,2,ix,iz,:) * st%maxwell_surface_grid_element(:)
      end do
    end do

    do ix=1, ix_max
      do iy=1, iy_max
        do ip_surf=1, st%maxwell_surface_grid_points_number(3,ix,iy)
          normal    =  M_z0
          normal(3) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(1,3,ix,iy,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) &
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(1,3,ix,iy,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
          normal    =  M_z0
          normal(3) =  M_z1
          xx(:) = gr%mesh%idx%lxyz(st%maxwell_surface_grid_points_map(2,3,ix,iy,ip_surf),:) &
                * gr%mesh%spacing(:)
          tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) &
             + zcross_product(normal(:),tmp_global(st%maxwell_surface_grid_points_map(2,3,ix,iy,ip_surf),:)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
        end do
        tmp_surf(1,3,ix,iy,:) = tmp_surf(1,3,ix,iy,:) / float(st%maxwell_surface_grid_points_number(3,ix,iy))
        tmp_surf(2,3,ix,iy,:) = tmp_surf(2,3,ix,iy,:) / float(st%maxwell_surface_grid_points_number(3,ix,iy))
      end do
    end do
    do ix=1, ix_max
      do iy=1, iy_max
        tmp_sum(:) = tmp_sum(:) - tmp_surf(1,3,ix,iy,:) * st%maxwell_surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2,3,ix,iy,:) * st%maxwell_surface_grid_element(:)
      end do
    end do

  end subroutine surface_integral_helmholtz_transverse

end module hamiltonian_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
