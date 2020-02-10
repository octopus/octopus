!! Copyright (C) 2019 R. Jestaedt, H. Appel, F. Bonafe, M. Oliveira, N. Tancogne-Dejean
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
  use batch_oct_m
  use batch_ops_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_abst_oct_m
  use hamiltonian_elec_oct_m
  use math_oct_m
  use maxwell_boundary_op_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use nl_operator_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_elec_dim_oct_m
  use states_elec_oct_m
  use states_mxll_oct_m

  implicit none

  private
  public ::                                     &
    hamiltonian_mxll_t,                         &
    hamiltonian_mxll_null,                      &
    hamiltonian_mxll_init,                      &
    hamiltonian_mxll_end,                       &
    dhamiltonian_mxll_apply,                    &
    zhamiltonian_mxll_apply,                    &
    dhamiltonian_mxll_magnus_apply,             &
    zhamiltonian_mxll_magnus_apply,             &
!   !    hamiltonian_mxll_apply_all,                 &
    hamiltonian_mxll_apply_batch,               &
    hamiltonian_mxll_span,                      &
    hamiltonian_mxll_adjoint,                   &
    hamiltonian_mxll_not_adjoint,               &
    hamiltonian_mxll_hermitian,                 &
    hamiltonian_mxll_update,                    &
    hamiltonian_mxll_get_time,                  &
    hamiltonian_mxll_apply_packed,              &
 !   maxwell_fft_hamiltonian,                    &
    maxwell_helmholtz_decomposition_trans_field,&
    maxwell_helmholtz_decomposition_long_field, &
    surface_integral_helmholtz_transverse


   type, extends(hamiltonian_abst_t) :: hamiltonian_mxll_t
    !> The Hamiltonian must know what are the "dimensions" of the spaces,
    !! in order to be able to operate on the states.
    type(states_elec_dim_t)       :: d

    !> absorbing boundaries
    logical :: adjoint

    FLOAT :: current_time
    logical :: apply_packed  !< This is initialized by the StatesPack variable.
    
    logical :: time_zero

    type(nl_operator_t), pointer   :: operators(:) 

!    type(poisson_t)                :: poisson
    FLOAT, pointer                 :: vector_potential(:,:)

    type(bc_mxll_t)                :: bc
    type(derivatives_t), pointer   :: der !< pointer to derivatives
    type(states_mxll_t), pointer   :: st

    integer                        :: rs_sign

    logical                        :: propagation_apply

    integer                        :: op_method

!    logical                        :: lorentz_force
!    logical                        :: lorentz_force_apply

    integer, pointer               :: rs_state_fft_map(:,:,:)
    integer, pointer               :: rs_state_fft_map_inv(:,:)

    logical                        :: mx_ma_coupling
    logical                        :: mx_ma_coupling_apply
    integer                        :: mx_ma_coupling_type
    integer                        :: mx_ma_trans_field_calc_method
    logical                        :: mx_ma_trans_field_calc_corr
    integer                        :: mx_ma_coupling_points_number
    FLOAT, pointer                 :: mx_ma_coupling_points(:,:)
    integer, pointer               :: mx_ma_coupling_points_map(:)
    integer                        :: mx_ma_coupling_order
    logical                        :: ma_mx_coupling
    logical                        :: ma_mx_coupling_apply

    logical                        :: bc_add_ab_region  = .false.
    logical                        :: bc_zero           = .false.
    logical                        :: bc_constant       = .false.
    logical                        :: bc_mirror_pec     = .false.
    logical                        :: bc_mirror_pmc     = .false.
    logical                        :: bc_periodic       = .false.
    logical                        :: bc_plane_waves    = .false.
    logical                        :: bc_medium         = .false.

    logical                        :: plane_waves
    logical                        :: plane_waves_apply
    logical                        :: spatial_constant
    logical                        :: spatial_constant_apply

!    logical                        :: diamagnetic_current
!    logical                        :: spin_current

    integer                        :: medium_calculation

    logical                        :: medium_box = .false.
    integer                        :: medium_box_number
    integer, pointer               :: medium_box_shape(:)
    FLOAT, pointer                 :: medium_box_center(:,:)
    FLOAT, pointer                 :: medium_box_size(:,:)
    FLOAT, pointer                 :: medium_box_ep(:,:)
    FLOAT, pointer                 :: medium_box_mu(:,:)
    FLOAT, pointer                 :: medium_box_c(:,:)
    FLOAT, pointer                 :: medium_box_ep_factor(:)
    FLOAT, pointer                 :: medium_box_mu_factor(:)
    FLOAT, pointer                 :: medium_box_sigma_e_factor(:)
    FLOAT, pointer                 :: medium_box_sigma_m_factor(:)
    FLOAT, pointer                 :: medium_box_sigma_e(:,:)
    FLOAT, pointer                 :: medium_box_sigma_m(:,:)
    integer, pointer               :: medium_box_points_number(:)
    FLOAT, pointer                 :: medium_box_points_map(:,:)
    FLOAT, pointer                 :: medium_box_aux_ep(:,:,:)
    FLOAT, pointer                 :: medium_box_aux_mu(:,:,:)
    integer, pointer               :: medium_box_bdry_number(:)
    FLOAT, pointer                 :: medium_box_bdry_map(:,:)
  
    !> maxwell hamiltonian_mxll
    integer                        :: operator
    logical                        :: current_density_ext_flag
    FLOAT                          :: energy
    FLOAT                          :: energy_boundaries
    FLOAT                          :: e_energy
    FLOAT                          :: b_energy
    FLOAT                          :: energy_plane_waves
    FLOAT                          :: e_energy_plane_waves
    FLOAT                          :: b_energy_plane_waves

    FLOAT                          :: energy_pml
    FLOAT                          :: energy_mask
    FLOAT, pointer                 :: energy_density(:)
    FLOAT, pointer                 :: energy_density_plane_waves(:)
    FLOAT, pointer                 :: e_energy_density(:)
    FLOAT, pointer                 :: b_energy_density(:)
    FLOAT                          :: energy_trans
    FLOAT                          :: energy_long
    FLOAT                          :: e_energy_trans
    FLOAT                          :: b_energy_trans
    FLOAT                          :: energy_incident_waves

    logical                        :: cpml_hamiltonian = .false.

    logical                        :: diamag_current = .false.

!    integer                        :: current_prop_test = 0

!    CMPLX, pointer                 :: test_output(:,:)

    type(cube_t)                   :: cube
    type(mesh_cube_parallel_map_t) :: mesh_cube_map

  contains
    procedure :: update_span => hamiltonian_mxll_span
    procedure :: dapply => dhamiltonian_mxll_apply
    procedure :: zapply => zhamiltonian_mxll_apply
    procedure :: dmagnus_apply => dhamiltonian_mxll_magnus_apply
    procedure :: zmagnus_apply => zhamiltonian_mxll_magnus_apply
    procedure :: is_hermitian => hamiltonian_mxll_hermitian
  end type hamiltonian_mxll_t

  type(profile_t), save :: prof_hamiltonian_mxll

  integer, public, parameter ::      &
    FARADAY_AMPERE_OLD          = 0, &
    FARADAY_AMPERE              = 1, &
    FARADAT_AMPERE_MEDIUM       = 2, &
    FARADAT_AMPERE_GAUSS        = 3, &
    FARADAT_AMPERE_GUASS_MEDIUM = 4

contains

  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_null(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_null)
 
    hm%adjoint = .false.
    
    hm%current_density_ext_flag = .false.

    nullify(hm%energy_density)
    nullify(hm%e_energy_density)
    nullify(hm%b_energy_density)

    POP_SUB(hamiltonian_mxll_null)
  end subroutine hamiltonian_mxll_null


  ! ---------------------------------------------------------
  !> Initializing the Maxwell Hamiltonian
  subroutine hamiltonian_mxll_init(hm, namespace, gr, st)
    type(hamiltonian_mxll_t),                   intent(inout) :: hm
    type(namespace_t),                          intent(in)    :: namespace
    type(grid_t),                       target, intent(inout) :: gr
    type(states_mxll_t),                target, intent(inout) :: st

    integer :: default_propagator
    type(profile_t), save :: prof

    PUSH_SUB(hamiltonian_mxll_init)

    call profiling_in(prof, 'HAMILTONIAN_INIT')
       
    call hamiltonian_mxll_null(hm)

    call states_elec_dim_copy(hm%d, st%d)

    ASSERT(associated(gr%der%lapl))

    hm%operators(1:3) => gr%der%grad(1:3) ! cross product for Maxwell calculation needs dimension >= 2
    hm%der => gr%der
    hm%rs_sign = st%rs_sign

    SAFE_ALLOCATE(hm%vector_potential(1:gr%mesh%np_part,1:st%d%dim))
    SAFE_ALLOCATE(hm%energy_density(1:gr%mesh%np_part))
    SAFE_ALLOCATE(hm%energy_density_plane_waves(1:gr%mesh%np_part))
    SAFE_ALLOCATE(hm%e_energy_density(1:gr%mesh%np_part))
    SAFE_ALLOCATE(hm%b_energy_density(1:gr%mesh%np_part))

    hm%vector_potential = M_ZERO
    hm%energy_density = M_ZERO
    hm%e_energy_density = M_ZERO
    hm%b_energy_density = M_ZERO

    SAFE_ALLOCATE(hm%cube%fft)

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
    call parse_variable(namespace, 'MaxwellHamiltonianOperator', OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE, hm%operator)

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then
      hm%d%dim = hm%d%dim+1
    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then
      hm%d%dim = 2*hm%d%dim
    end if

    !%Variable MaxwellExternalCurrent
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Description follows
    !%End
    call parse_variable(namespace, 'MaxwellExternalCurrent', .false., hm%current_density_ext_flag)

    hm%plane_waves_apply = .false.
    hm%spatial_constant_apply = .false.

    hm%propagation_apply = .false.

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
    call parse_variable(namespace, 'MaxwellMediumCalculation', default_propagator, hm%medium_calculation)
    if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) then
      call messages_not_implemented("Calculation from E and B field not implemented yet.", namespace=namespace)
    end if

    hm%rs_state_fft_map     => st%rs_state_fft_map
    hm%rs_state_fft_map_inv => st%rs_state_fft_map_inv

    call parse_variable(namespace, 'StatesPack', .true., hm%apply_packed)

    call parse_variable(namespace, 'TimeZero', .false., hm%time_zero)
    if(hm%time_zero) call messages_experimental('TimeZero')

    call profiling_out(prof)
    POP_SUB(hamiltonian_mxll_init)
  end subroutine hamiltonian_mxll_init

  
  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_end(hm)
    type(hamiltonian_mxll_t), intent(inout) :: hm

    PUSH_SUB(hamiltonian_mxll_end)

    nullify(hm%operators)

    SAFE_DEALLOCATE_P(hm%energy_density)
    SAFE_DEALLOCATE_P(hm%cube%fft)

    call bc_mxll_end(hm%bc)

    call states_elec_dim_end(hm%d) 

    POP_SUB(hamiltonian_mxll_end)
  end subroutine hamiltonian_mxll_end


  ! ---------------------------------------------------------
  logical function hamiltonian_mxll_hermitian(hm)
    class(hamiltonian_mxll_t), intent(in) :: hm

    PUSH_SUB(hamiltonian_mxll_hermitian)

    hamiltonian_mxll_hermitian = .true.

    POP_SUB(hamiltonian_mxll_hermitian)
  end function hamiltonian_mxll_hermitian


  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_span(hm, delta, emin)
    class(hamiltonian_mxll_t), intent(inout) :: hm
    FLOAT,                     intent(in)    :: delta, emin

    write(message(1),'(a)') 'hamiltonian_mxll_span is not needed and hence not implemented.'
    call messages_fatal(1)

  end subroutine hamiltonian_mxll_span


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
  !> Maxwell Hamiltonian update (here only the time is updated, can maybe be added to another routine)
  subroutine hamiltonian_mxll_update(this, time)
    type(hamiltonian_mxll_t), intent(inout) :: this
    FLOAT,          optional, intent(in)    :: time

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
    type(mesh_t),               intent(in) :: mesh

    apply = this%apply_packed
    if(mesh%use_curvilinear) apply = .false.
    
  end function hamiltonian_mxll_apply_packed

  ! ---------------------------------------------------------
  subroutine hamiltonian_mxll_apply_batch(hm, namespace, der, psib, hpsib, time, terms, set_bc)
    type(hamiltonian_mxll_t),  intent(in)    :: hm
    type(namespace_t),         intent(in)    :: namespace
    type(derivatives_t),       intent(in)    :: der
    type(batch_t), target,     intent(inout) :: psib
    type(batch_t), target,     intent(inout) :: hpsib
    FLOAT, optional,           intent(in)    :: time
    integer, optional,         intent(in)    :: terms
    logical, optional,         intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.

    logical :: pack
    integer :: terms_

    PUSH_SUB(hamiltonian_mxll_apply_batch)
    call profiling_in(prof_hamiltonian_mxll, "MXLL_HAMILTONIAN")

    ASSERT(psib%status() == hpsib%status())

    ASSERT(psib%is_ok())
    ASSERT(hpsib%is_ok())
    ASSERT(psib%nst == hpsib%nst)

    !Not implemented at the moment
    ASSERT(.not.present(terms))
    ASSERT(.not.present(set_bc))

    if(present(time)) then
      if(abs(time - hm%current_time) > CNST(1e-10)) then
        write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
        write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', hm%current_time
        call messages_fatal(2, namespace=namespace)
      endif
    end if

    call zderivatives_curl(der, psib%states(1)%zpsi, hpsib%states(1)%zpsi)
    hpsib%states(1)%zpsi(:,:) = P_c * hpsib%states(1)%zpsi(:,:)
  
    call profiling_out(prof_hamiltonian_mxll)
    POP_SUB(hamiltonian_mxll_apply_batch)
  end subroutine hamiltonian_mxll_apply_batch

  
  ! --------------------------------------------------------
  !> Apply hamiltonian to real states (not possible)
  subroutine dhamiltonian_mxll_apply(hm, namespace, mesh, psib, hpsib, terms, set_bc)
    class(hamiltonian_mxll_t),   intent(in)    :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(mesh_t),                intent(in)    :: mesh
    class(batch_t),      target, intent(inout) :: psib
    class(batch_t),      target, intent(inout) :: hpsib
    integer,           optional, intent(in)    :: terms
    logical,           optional, intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.

    write(message(1),'(a)') 'dhamiltonian_mxll_apply not implemented (states are complex).'
    call messages_fatal(1, namespace=namespace)
    
  end subroutine dhamiltonian_mxll_apply

  ! ---------------------------------------------------------
  !> Applying the Maxwell Hamiltonian on Maxwell states
  subroutine zhamiltonian_mxll_apply(hm, namespace, mesh, psib, hpsib, terms, set_bc)
    class(hamiltonian_mxll_t),   intent(in)    :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(mesh_t),                intent(in)    :: mesh
    class(batch_t),      target, intent(inout) :: psib
    class(batch_t),      target, intent(inout) :: hpsib
    integer,           optional, intent(in)    :: terms
    logical,           optional, intent(in)    :: set_bc !< If set to .false. the boundary conditions are assumed to be set previously.
    
    CMPLX, allocatable :: rs_aux_in(:,:), rs_aux_out(:,:)
    integer :: ii

    PUSH_SUB(zhamiltonian_mxll_apply)

!    call batch_init(psib, hm%d%dim, 1)
!    call psib%add_state(ist, psi)
!    call batch_init(hpsib, hm%d%dim, 1)
!    call hpsib%add_state(ist, hpsi)
!
!    call hamiltonian_mxll_apply_batch(hm, der, psib, hpsib, ik, time = time, terms = terms, Imtime = Imtime, set_bc = set_bc)
!
!    call psib%end()
!    call hpsib%end()

    call profiling_in(prof_hamiltonian_mxll, "MAXWELLHAMILTONIAN")

!    if (present(time)) then
!      if (abs(time - hm%current_time) > CNST(1e-10)) then
!        write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
!        write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', hm%current_time
!        call messages_fatal(2, namespace=namespace)
!      end if
!    end if

    SAFE_ALLOCATE(rs_aux_in(1:mesh%np, 1:3))
    SAFE_ALLOCATE(rs_aux_out(1:mesh%np, 1:3))

    do ii = 1, 3
      call batch_get_state(psib, mesh%np, ii, rs_aux_in(:, ii))
    end do

    call maxwell_hamiltonian_apply_fd(hm, hm%der, rs_aux_in, rs_aux_out)
    !    call hamiltonian_mxll_apply_batch(hm, namespace, hm%der, psib, hpsib)

    do ii = 1, 3
      call batch_set_state(hpsib, ii, mesh%np, rs_aux_out)
    end do

    call profiling_out(prof_hamiltonian_mxll)

    POP_SUB(zhamiltonian_mxll_apply)
  end subroutine zhamiltonian_mxll_apply


  ! ---------------------------------------------------------
!  subroutine hamiltonian_mxll_apply_all(hm, namespace, der, st, hst, time) removed
!  end subroutine hamiltonian_mxll_apply_all

  ! ---------------------------------------------------------
  !> Applying the Maxwell Hamiltonian on Maxwell states with finite difference
  subroutine maxwell_hamiltonian_apply_fd(hm, der, psi, oppsi)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    CMPLX,                    intent(inout) :: psi(:,:)
    CMPLX,                    intent(inout) :: oppsi(:,:)

    FLOAT, pointer     :: mx_rho(:,:)
    CMPLX, allocatable :: tmp(:,:)
    CMPLX, pointer     :: kappa_psi(:,:)
    integer            :: np, np_part, ip, ip_in, rs_sign

    PUSH_SUB(maxwell_hamiltonian_apply_fd)

    np = der%mesh%np
    np_part = der%mesh%np_part
    rs_sign = hm%rs_sign


    !===========================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum via partial derivatives:

    if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE) then

      SAFE_ALLOCATE(tmp(np_part,2))
      oppsi       = M_z0

      if (hm%diamag_current) then
        mx_rho    => hm%st%grid_rho
        kappa_psi => hm%st%kappa_psi 
      end if

      !----------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 3, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 2, 3, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 3, 2, tmp(:,2))
      oppsi(1:np_part,1) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 1, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 3, 1, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 1, 3, tmp(:,2))
      oppsi(1:np_part,2) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 2, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 1, 2, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 2, 1, tmp(:,2))
      oppsi(1:np_part,3) = ( tmp(1:np_part,1)-tmp(1:np_part,2) )

      if (hm%bc_constant) then
        do ip_in = 1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          oppsi(ip,:) = hm%st%rs_state_const(:)
        end do
      end if

      SAFE_DEALLOCATE_A(tmp)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium via partial derivatives:

    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_MEDIUM) then

      SAFE_ALLOCATE(tmp(np_part,4))
      oppsi       = M_z0

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,3), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,4), 3, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 2, 3, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 3, 2, tmp(:,3:4))
      oppsi(1:np_part,1) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,4) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,4), 1, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 3, 1, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 1, 3, tmp(:,3:4))
      oppsi(1:np_part,2) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,5) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,3), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,4), 2, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 1, 2, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 2, 1, tmp(:,3:4))
      oppsi(1:np_part,3) =   ( tmp(1:np_part,1)-tmp(1:np_part,3) )
      oppsi(1:np_part,6) = - ( tmp(1:np_part,2)-tmp(1:np_part,4) )


      SAFE_DEALLOCATE_A(tmp)

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation if medium boundaries is set:
      call maxwell_medium_boundaries_calculation(hm, psi, oppsi)

      !----------------------------------------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation for medium boxes:
      call maxwell_medium_boxes_calculation(hm, der, psi, oppsi)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum with Gauss condition via partial derivatives:

    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS) then

      SAFE_ALLOCATE(tmp(np_part,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,1) = rs_sign * P_c * ( - M_zI*tmp(1:np_part,1) - tmp(1:np_part,2) - M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,2) = rs_sign * P_c * ( - M_zI*tmp(1:np_part,1) - tmp(1:np_part,2) - M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,3) = rs_sign * P_c * ( tmp(1:np_part,1) - M_zI*tmp(1:np_part,2) + M_zI*tmp(1:np_part,3) )

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,4) = rs_sign * P_c * ( tmp(1:np_part,1) - M_zI*tmp(1:np_part,2) + M_zI*tmp(1:np_part,3) )

      SAFE_DEALLOCATE_A(tmp)


    !==============================================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium with Gauss condition via partial derivatives:

    else if (hm%operator == OPTION__MAXWELLHAMILTONIANOPERATOR__FARADAY_AMPERE_GAUSS_MEDIUM) then

      SAFE_ALLOCATE(tmp(np_part,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,1) = P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,2) = P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,3) = P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,4) = P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,5), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,5) = - P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,6), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,8), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,8), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np_part,6) = - P_c*(-M_zI*tmp(1:np_part,1)-tmp(1:np_part,2)-M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,5), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,7) = - P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      call zderivatives_partial(der, psi(:,6), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np_part,8) = - P_c*(tmp(1:np_part,1)-M_zI*tmp(1:np_part,2)+M_zI*tmp(1:np_part,3))

      SAFE_DEALLOCATE_A(tmp)

    end if

    POP_SUB(maxwell_hamiltonian_apply_fd)
  end subroutine maxwell_hamiltonian_apply_fd


  !> Maxwell Hamiltonian is updated for the PML calculation
  subroutine maxwell_pml_hamiltonian(hm, der, psi, dir1, dir2, tmp)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    CMPLX,                    intent(inout) :: psi(:,:)
    integer,                  intent(in)    :: dir1
    integer,                  intent(in)    :: dir2
    CMPLX,                    intent(inout) :: tmp(:)

    PUSH_SUB(maxwell_pml_hamiltonian)

    if ( (hm%bc%bc_ab_type(dir1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) .and. &
          hm%cpml_hamiltonian ) then
      if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) then
        call maxwell_pml_calculation_via_riemann_silberstein(hm, der, psi, dir1, dir2, tmp(:))
      end if
    end if

    POP_SUB(maxwell_pml_hamiltonian)
  end subroutine maxwell_pml_hamiltonian

  ! ---------------------------------------------------------
  !> Maxwell Hamiltonian is updated for the PML calculation
  subroutine maxwell_pml_hamiltonian_medium(hm, der, psi, dir1, dir2, tmp)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    CMPLX,                    intent(inout) :: psi(:,:)
    integer,                  intent(in)    :: dir1
    integer,                  intent(in)    :: dir2
    CMPLX,                    intent(inout) :: tmp(:,:)

    PUSH_SUB(maxwell_pml_hamiltonian_medium)

    if ( (hm%bc%bc_ab_type(dir1) == OPTION__MAXWELLABSORBINGBOUNDARIES__CPML) .and. &
          hm%cpml_hamiltonian ) then
      if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) then
        call maxwell_pml_calculation_via_riemann_silberstein_medium(hm, der, psi, dir1, dir2, tmp(:,:))
!      else if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__ELECTRIC_MAGNETIC_FIELDS) then
!        call maxwell_pml_calculation_via_e_b_fields_medium(hm, der, psi, dir1, dir2, tmp(:,:))
      end if
    end if

    POP_SUB(maxwell_pml_hamiltonian_medium)
  end subroutine maxwell_pml_hamiltonian_medium

  ! ---------------------------------------------------------
  !> Maxwell Hamiltonian is updated for the PML calculation via Riemann-Silberstein vector
  subroutine maxwell_pml_calculation_via_riemann_silberstein(hm, der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    integer,                  intent(in)    :: pml_dir
    CMPLX,                    intent(inout) :: psi(:,:)
    integer,                  intent(in)    :: field_dir
    CMPLX,                    intent(inout) :: pml(:)

    integer            :: ip, ip_in, np_part, rs_sign
    FLOAT              :: pml_c(3)
    CMPLX, allocatable :: tmp_partial(:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g(3)

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein)

    if (hm%cpml_hamiltonian) then

      rs_sign = hm%rs_sign

      np_part = der%mesh%np_part
      SAFE_ALLOCATE(tmp_partial(np_part))

      call zderivatives_partial(der, psi(:,field_dir), tmp_partial(:), pml_dir, set_bc = .false.)
      do ip_in = 1, hm%bc%pml_points_number
        ip       = hm%bc%pml_points_map(ip_in)
        pml_c(:) = hm%bc%pml_c(ip_in, :)
        pml_a(:) = hm%bc%pml_a(ip_in, :)
        pml_b(:) = hm%bc%pml_b(ip_in, :)
        pml_g(:) = hm%bc%pml_conv_plus(ip_in, pml_dir,:)
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
  !> Maxwell Hamiltonian is updated for the PML calculation via Riemann-Silberstein 
  !> vector with medium inside the box
  subroutine maxwell_pml_calculation_via_riemann_silberstein_medium(hm, der, psi, pml_dir, field_dir, pml)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    integer,                  intent(in)    :: pml_dir
    CMPLX,                    intent(inout) :: psi(:,:)
    integer,                  intent(in)    :: field_dir
    CMPLX,                    intent(inout) :: pml(:,:)

    integer            :: ip, ip_in, np_part
    FLOAT              :: pml_c(3)
    CMPLX, allocatable :: tmp_partial(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g_p(3), pml_g_m(3)

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein_medium)

    if (hm%cpml_hamiltonian) then

      np_part = der%mesh%np_part
      SAFE_ALLOCATE(tmp_partial(np_part,1:2))

      call zderivatives_partial(der, psi(:,field_dir  ), tmp_partial(:,1), pml_dir, set_bc = .false.)
      call zderivatives_partial(der, psi(:,field_dir+3), tmp_partial(:,2), pml_dir, set_bc = .false.)
      do ip_in = 1, hm%bc%pml_points_number
        ip         = hm%bc%pml_points_map(ip_in)
        pml_c(:)   = hm%bc%pml_c(ip_in, :)
        pml_a(:)   = hm%bc%pml_a(ip_in, :)
        pml_b(:)   = hm%bc%pml_b(ip_in, :)
        pml_g_p(:) = hm%bc%pml_conv_plus(ip_in, pml_dir, :)
        pml_g_m(:) = hm%bc%pml_conv_minus(ip_in, pml_dir, :)
        pml(ip,1)  = pml_c(pml_dir) * tmp_partial(ip, 1) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * real(tmp_partial(ip, 1)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip, 1)) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_p(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_p(field_dir))
        pml(ip,2)  = pml_c(pml_dir) * tmp_partial(ip, 2) &
                   + pml_c(pml_dir) * real(pml_a(pml_dir)) * real(tmp_partial(ip, 2)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip, 2)) &
                   + pml_c(pml_dir) * real(pml_b(pml_dir)) * real(pml_g_m(field_dir)) &
                   + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_m(field_dir))
      end do

    end if

    SAFE_DEALLOCATE_A(tmp_partial)

    POP_SUB(maxwell_pml_calculation_via_riemann_silberstein_medium)
  end subroutine maxwell_pml_calculation_via_riemann_silberstein_medium

  
  ! ---------------------------------------------------------
  !> Maxwell Hamiltonian for medium boundaries
  subroutine maxwell_medium_boundaries_calculation(hm, psi, oppsi)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    CMPLX,                    intent(in)    :: psi(:,:)
    CMPLX,                    intent(inout) :: oppsi(:,:)

    integer            :: ip, ip_in, idim
    FLOAT              :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m
    CMPLX              :: ff_plus(3), ff_minus(3)

    PUSH_SUB(maxwell_medium_boundaries_calculation)

    do idim = 1, 3
      if ( (hm%bc%bc_type(idim) == OPTION__MAXWELLBOUNDARYCONDITIONS__MAXWELL_MEDIUM) .and. &
           (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) ) then
        do ip_in = 1, hm%bc%medium_points_number(idim)
          ip          = hm%bc%medium_points_map(ip_in, idim)
          cc          = hm%bc%medium_c(ip_in, idim)/P_c
          aux_ep(:)   = hm%bc%medium_aux_ep(ip_in, :, idim)
          aux_mu(:)   = hm%bc%medium_aux_mu(ip_in, :, idim)
          sigma_e     = hm%bc%medium_sigma_e(ip_in, idim)
          sigma_m     = hm%bc%medium_sigma_m(ip_in, idim)
          ff_plus(1)  = psi(ip, 1)
          ff_plus(2)  = psi(ip, 2)
          ff_plus(3)  = psi(ip, 3)
          ff_minus(1) = psi(ip, 4)
          ff_minus(2) = psi(ip, 5)
          ff_minus(3) = psi(ip, 6)
          aux_ep      = dcross_product(aux_ep,real(ff_plus+ff_minus))
          aux_mu      = dcross_product(aux_mu,aimag(ff_plus-ff_minus))
          oppsi(ip, 1) = oppsi(ip, 1)*cc                                         &
                       - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * real(ff_plus(1) + ff_minus(1))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 4) = oppsi(ip, 4)*cc                                         &
                       + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * real(ff_plus(1) + ff_minus(1))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 2) = oppsi(ip, 2)*cc                                         &
                       - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * real(ff_plus(2) + ff_minus(2))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2))
          oppsi(ip, 5) = oppsi(ip, 5)*cc                                         &
                       + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * real(ff_plus(2) + ff_minus(2))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2)) 
          oppsi(ip, 3) = oppsi(ip, 3)*cc                                         &
                       - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * real(ff_plus(3) + ff_minus(3))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
          oppsi(ip, 6) = oppsi(ip, 6)*cc                                         &
                       + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * real(ff_plus(3) + ff_minus(3))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
        end do
      end if
    end do

    POP_SUB(maxwell_medium_boundaries_calculation)
  end subroutine maxwell_medium_boundaries_calculation


  ! ---------------------------------------------------------
  ! > Maxwell Hamiltonian including medium boxes
  subroutine maxwell_medium_boxes_calculation(hm, der, psi, oppsi)
    type(hamiltonian_mxll_t), intent(in)    :: hm
     type(derivatives_t),      intent(in)    :: der
    CMPLX,                    intent(in)    :: psi(:,:)
    CMPLX,                    intent(inout) :: oppsi(:,:)

    integer            :: ip, ip_in, il, np_part
    FLOAT              :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m
    FLOAT, allocatable :: tmp_e(:,:), tmp_b(:,:), tmp_curl_e(:,:), tmp_curl_b(:,:)
    CMPLX              :: ff_plus(3), ff_minus(3)

    PUSH_SUB(maxwell_medium_boxes_calculation)

    np_part = der%mesh%np_part

    if (hm%medium_box .and. &
         (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RIEMANN_SILBERSTEIN) ) then
      do il = 1, hm%medium_box_number
        do ip_in = 1, hm%medium_box_points_number(il)
          ip           = hm%medium_box_points_map(ip_in, il)
          cc           = hm%medium_box_c(ip_in,il)/P_c
          aux_ep(:)    = hm%medium_box_aux_ep(ip_in, :, il)
          aux_mu(:)    = hm%medium_box_aux_mu(ip_in, :, il)
          sigma_e      = hm%medium_box_sigma_e(ip_in, il)
          sigma_m      = hm%medium_box_sigma_m(ip_in, il)
          ff_plus(1)   = psi(ip, 1)
          ff_plus(2)   = psi(ip, 2)
          ff_plus(3)   = psi(ip, 3)
          ff_minus(1)  = psi(ip, 4)
          ff_minus(2)  = psi(ip, 5)
          ff_minus(3)  = psi(ip, 6)
          aux_ep       = dcross_product(aux_ep,real(ff_plus+ff_minus))
          aux_mu       = dcross_product(aux_mu,aimag(ff_plus-ff_minus))
          oppsi(ip, 1) = oppsi(ip,1)*cc                                          &
                       - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * real(ff_plus(1) + ff_minus(1))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 4) = oppsi(ip,4)*cc                                          &
                       + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * real(ff_plus(1) + ff_minus(1))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 2) = oppsi(ip,2)*cc                                          &
                       - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * real(ff_plus(2) + ff_minus(2))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2))
          oppsi(ip, 5) = oppsi(ip,5)*cc                                          &
                       + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * real(ff_plus(2) + ff_minus(2))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2)) 
          oppsi(ip, 3) = oppsi(ip,3)*cc                                          &
                       - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * real(ff_plus(3) + ff_minus(3))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
          oppsi(ip, 6) = oppsi(ip,6)*cc                                          &
                       + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * real(ff_plus(3) + ff_minus(3))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
        end do
      end do
    end if

    POP_SUB(maxwell_medium_boxes_calculation)
  end subroutine maxwell_medium_boxes_calculation

  !----------------------------------------------------------
  !>  Helmholtz decomposition to calculate a transverse field (maybe should be a general math function)
  subroutine maxwell_helmholtz_decomposition_trans_field(poisson, gr, hm_mxll, hm_elec, st, transverse_field)
    type(poisson_t),          intent(in)    :: poisson
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_mxll_t), intent(in)    :: hm_mxll
    type(hamiltonian_elec_t), intent(in)    :: hm_elec
    type(states_mxll_t),      intent(in)    :: st
    CMPLX,                    intent(inout) :: transverse_field(:,:)

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
    do idim = 1, st%d%dim
      call zpoisson_solve(poisson, tmp_poisson(:), ztmp(:, idim), .true.)
      ztmp(1:gr%mesh%np_part, idim) = M_ONE / (M_FOUR * M_PI) * tmp_poisson(1:gr%mesh%np_part)
    end do
    call zderivatives_curl(gr%der, ztmp, transverse_field, set_bc = .false.)
    do ip_in = 1, hm_mxll%bc%der_bndry_mask_points_number
      ip = hm_mxll%bc%der_bndry_mask_points_map(ip_in)
      transverse_field(ip, :) = hm_mxll%bc%der_bndry_mask(ip_in) * transverse_field(ip, :)
    end do

    SAFE_DEALLOCATE_A(ztmp)
    SAFE_DEALLOCATE_A(tmp_poisson)

!    ! correction surface integral
!    if (hm_elec%mx_ma_trans_field_calc_corr) then
!      call mesh_nearest_point_infos(gr%mesh, hm_elec%mx_ma_coupling_points(:,1), dmin, rankmin, &
!                                    ip_local, ip_global)
!      ip_local  = 1
!      ip_global = 1
!      do ii = -gr%der%order, gr%der%order
!        do jj = -gr%der%order, gr%der%order
!          do kk = -gr%der%order, gr%der%order
!            pos(1) = gr%mesh%x(ip_local,1) * ii * gr%mesh%spacing(1)
!            pos(2) = gr%mesh%x(ip_local,2) * jj * gr%mesh%spacing(2)
!            pos(3) = gr%mesh%x(ip_local,3) * kk * gr%mesh%spacing(3)
!            call surface_integral_helmholtz_transverse(gr, st, pos, field_old, surface_integral)
!            transverse_field(ip_local,:) = transverse_field(ip_local,:) - M_ONE / (M_FOUR * M_PI) * surface_integral(:)
!          end do
!        end do
!      end do
!      pos(:) = gr%mesh%x(ip_local,:)
!    end if

    SAFE_DEALLOCATE_A(field_old)

    POP_SUB(maxwell_helmholtz_decomposition_trans_field)
  end subroutine maxwell_helmholtz_decomposition_trans_field

  !----------------------------------------------------------
  !> Helmholtz decomposition to calculate a longitudinal field (maybe should be a general math function)
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
    call zderivatives_grad(gr%der, ztmp(:), longitudinal_field(:,:), set_bc = .false.)

    SAFE_DEALLOCATE_A(ztmp)
    SAFE_DEALLOCATE_A(tmp_poisson)

  end subroutine maxwell_helmholtz_decomposition_long_field

  ! ---------------------------------------------------------
  !> Surface integral of the Helmholtz decomposition to calculate the transverse field
  !> (maybe should be a general math function)
  subroutine surface_integral_helmholtz_transverse(gr, st, pos, field, surface_integral)
    type(grid_t),        intent(in)    :: gr
    type(states_mxll_t), intent(in)    :: st
    FLOAT,               intent(in)    :: pos(:) 
    CMPLX,               intent(in)    :: field(:,:)
    CMPLX,               intent(inout) :: surface_integral(:)

    integer             :: idim, ip_surf, ix, ix_max, iy, iy_max, iz, iz_max, ii_max
    FLOAT               :: xx(3)
    CMPLX               :: tmp_sum(3), normal(3)
    CMPLX,  allocatable :: tmp_global(:,:), tmp_surf(:,:,:,:,:)

    SAFE_ALLOCATE(tmp_global(1:gr%mesh%np_part_global,1:st%d%dim))

!    if (gr%mesh%parallel_in_domains) then
!      do idim=1, st%d%dim
!#if defined(HAVE_MPI)
!        call vec_allgather(gr%mesh%vp, tmp_global(:,idim), field(:,idim))
!        call MPI_Barrieri(gr%mesh%vp%comm, mpi_err)
!#endif
!      end do
!    else
      tmp_global(:,:) = field(:,:)
!    end if

    ix_max = st%surface_grid_rows_number(1)
    iy_max = st%surface_grid_rows_number(2)
    iz_max = st%surface_grid_rows_number(3)
    ii_max = max(ix_max,iy_max,iz_max)

    SAFE_ALLOCATE(tmp_surf(1:2,1:st%d%dim,1:ii_max,1:ii_max,1:st%d%dim))

    tmp_surf = M_z0
    tmp_sum  = M_z0

    do iy = 1, iy_max
      do iz = 1, iz_max
        do ip_surf = 1, st%surface_grid_points_number(1,iy,iz)
          normal    =  M_z0
          normal(1) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(1, 1, iy, iz, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(1, 1, iy, iz,:) = tmp_surf(1, 1, iy, iz, :) & 
             + zcross_product(normal(:), tmp_global(st%surface_grid_points_map(1, 1, iy, iz, ip_surf), :)) &
             / sqrt( (xx(1) - pos(1))**2 + (xx(2) - pos(2))**2 + (xx(3) - pos(3))**2 )
          normal    =  M_z0
          normal(1) = +M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(2, 1, iy, iz, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(2, 1, iy, iz, :) = tmp_surf(2, 1, iy, iz, :) &
             + zcross_product(normal(:), tmp_global(st%surface_grid_points_map(2, 1, iy, iz, ip_surf), :)) &
             / sqrt( (xx(1) - pos(1))**2 + (xx(2) - pos(2))**2 + (xx(3) - pos(3))**2 )
        end do
        tmp_surf(1, 1, iy, iz, :) = tmp_surf(1, 1, iy, iz, :) / float(st%surface_grid_points_number(1, iy, iz))
        tmp_surf(2, 1, iy, iz, :) = tmp_surf(2, 1, iy, iz, :) / float(st%surface_grid_points_number(1, iy, iz))
      end do
    end do
    do iy = 1, iy_max
      do iz = 1, iz_max
        tmp_sum(:) = tmp_sum(:) + tmp_surf(1, 1, iy, iz, :) * st%surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2, 1, iy, iz, :) * st%surface_grid_element(:)
      end do
    end do

    do ix = 1, ix_max
      do iz = 1, iz_max
        do ip_surf = 1, st%surface_grid_points_number(2, ix, iz)
          normal    =  M_z0
          normal(2) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(1, 2, ix, iz, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(1, 2, ix, iz, :) = tmp_surf(1, 2, ix, iz, :) &
             + zcross_product(normal(:), tmp_global(st%surface_grid_points_map(1, 2, ix, iz, ip_surf), :)) &
             / sqrt( (xx(1) - pos(1))**2 + (xx(2) - pos(2))**2 + (xx(3) - pos(3))**2 )
          normal    =  M_z0
          normal(2) =  M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(2, 2, ix, iz, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(2, 2, ix, iz, :) = tmp_surf(2, 2, ix, iz, :) &
             + zcross_product(normal(:), tmp_global(st%surface_grid_points_map(2, 2, ix, iz, ip_surf), :)) &
             / sqrt( (xx(1) - pos(1))**2 + (xx(2) - pos(2))**2 + (xx(3) - pos(3))**2 )
        end do
        tmp_surf(1, 2, ix, iz, :) = tmp_surf(1, 2, ix, iz, :) / float(st%surface_grid_points_number(2, ix, iz))
        tmp_surf(2, 2, ix, iz, :) = tmp_surf(2, 2, ix, iz, :) / float(st%surface_grid_points_number(2, ix, iz))
      end do
    end do
    do ix = 1, ix_max
      do iz = 1, iz_max
        tmp_sum(:) = tmp_sum(:) + tmp_surf(1, 2, ix, iz, :) * st%surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2, 2, ix, iz, :) * st%surface_grid_element(:)
      end do
    end do

    do ix = 1, ix_max
      do iy = 1, iy_max
        do ip_surf = 1, st%surface_grid_points_number(3, ix, iy)
          normal    =  M_z0
          normal(3) = -M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(1, 3, ix, iy, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(1, 3, ix, iy, :) = tmp_surf(1, 3, ix, iy, :) &
             + zcross_product(normal(:),tmp_global(st%surface_grid_points_map(1, 3, ix, iy, ip_surf), :)) &
             / sqrt( (xx(1) - pos(1))**2 + (xx(2) - pos(2))**2 + (xx(3) - pos(3))**2 )
          normal    =  M_z0
          normal(3) =  M_z1
          xx(:) = gr%mesh%idx%lxyz(st%surface_grid_points_map(2, 3, ix, iy, ip_surf), :) &
                * gr%mesh%spacing(1:3)
          tmp_surf(2, 3, ix, iy, :) = tmp_surf(2, 3, ix, iy, :) &
             + zcross_product(normal(:), tmp_global(st%surface_grid_points_map(2, 3, ix, iy, ip_surf), :)) &
             / sqrt( (xx(1)-pos(1))**2 + (xx(2)-pos(2))**2 + (xx(3)-pos(3))**2 )
        end do
        tmp_surf(1, 3, ix, iy, :) = tmp_surf(1, 3, ix, iy, :) / float(st%surface_grid_points_number(3, ix, iy))
        tmp_surf(2, 3, ix, iy, :) = tmp_surf(2, 3, ix, iy, :) / float(st%surface_grid_points_number(3, ix, iy))
      end do
    end do
    do ix = 1, ix_max
      do iy = 1, iy_max
        tmp_sum(:) = tmp_sum(:) - tmp_surf(1, 3, ix, iy, :) * st%surface_grid_element(:)
        tmp_sum(:) = tmp_sum(:) + tmp_surf(2, 3, ix, iy, :) * st%surface_grid_element(:)
      end do
    end do

  end subroutine surface_integral_helmholtz_transverse


  ! ---------------------------------------------------------
  !> Maxwell hamiltonian Magnus (not implemeted)
  subroutine dhamiltonian_mxll_magnus_apply(hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
    class(hamiltonian_mxll_t),   intent(in)    :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(mesh_t),                intent(in)    :: mesh
    class(batch_t),              intent(inout) :: psib
    class(batch_t),              intent(inout) :: hpsib
    FLOAT,                       intent(in)    :: vmagnus(:, :, :)
    logical,           optional, intent(in)    :: set_phase
    
    write(message(1),'(a)') 'hamiltonian_mxll_magnus_apply not implemeted'
    call messages_fatal(1, namespace=namespace)
    
  end subroutine dhamiltonian_mxll_magnus_apply

  ! ---------------------------------------------------------
  !> Maxwell hamiltonian Magnus (not implemeted)
  subroutine zhamiltonian_mxll_magnus_apply(hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
    class(hamiltonian_mxll_t),   intent(in)    :: hm
    type(namespace_t),           intent(in)    :: namespace
    type(mesh_t),                intent(in)    :: mesh
    class(batch_t),              intent(inout) :: psib
    class(batch_t),              intent(inout) :: hpsib
    FLOAT,                       intent(in)    :: vmagnus(:, :, :)
    logical,           optional, intent(in)    :: set_phase
    
    write(message(1),'(a)') 'hamiltonian_mxll_magnus_apply not implemeted'
    call messages_fatal(1, namespace=namespace)
    
  end subroutine zhamiltonian_mxll_magnus_apply

end module hamiltonian_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
