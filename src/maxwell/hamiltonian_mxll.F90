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
  use boundaries_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use energy_mxll_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_abst_oct_m
  use hamiltonian_elec_oct_m
  use math_oct_m
  use linear_medium_em_field_oct_m
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
    hamiltonian_mxll_init,                      &
    hamiltonian_mxll_end,                       &
    dhamiltonian_mxll_apply,                    &
    zhamiltonian_mxll_apply,                    &
    dhamiltonian_mxll_magnus_apply,             &
    zhamiltonian_mxll_magnus_apply,             &
    hamiltonian_mxll_apply_batch,               &
    hamiltonian_mxll_span,                      &
    hamiltonian_mxll_adjoint,                   &
    hamiltonian_mxll_not_adjoint,               &
    hamiltonian_mxll_hermitian,                 &
    hamiltonian_mxll_update,                    &
    hamiltonian_mxll_get_time,                  &
    hamiltonian_mxll_apply_packed,              &
    maxwell_helmholtz_decomposition_trans_field,&
    maxwell_helmholtz_decomposition_long_field

   type, extends(hamiltonian_abst_t) :: hamiltonian_mxll_t
    integer                        :: dim
    !> absorbing boundaries
    logical :: adjoint = .false.

    FLOAT :: current_time
    logical :: apply_packed  !< This is initialized by the StatesPack variable.

    logical :: time_zero

    type(nl_operator_t), pointer   :: operators(:)

    FLOAT, allocatable             :: vector_potential(:,:)

    type(bc_mxll_t)                :: bc
    type(derivatives_t), pointer, private :: der !< pointer to derivatives
    type(states_mxll_t), pointer   :: st

    integer                        :: rs_sign

    logical                        :: propagation_apply = .false.

    integer                        :: op_method

    integer, pointer               :: rs_state_fft_map(:,:,:)
    integer, pointer               :: rs_state_fft_map_inv(:,:)

    logical                        :: mx_ma_coupling = .false.
    logical                        :: mx_ma_coupling_apply = .false.
    integer                        :: mx_ma_trans_field_calc_method
    logical                        :: mx_ma_trans_field_calc_corr = .false.
    integer                        :: mx_ma_coupling_points_number
    FLOAT,   allocatable           :: mx_ma_coupling_points(:,:)
    integer, allocatable           :: mx_ma_coupling_points_map(:)
    integer                        :: mx_ma_coupling_order
    logical                        :: ma_mx_coupling       = .false.
    logical                        :: ma_mx_coupling_apply = .false.

    logical                        :: bc_add_ab_region  = .false.
    logical                        :: bc_zero           = .false.
    logical                        :: bc_constant       = .false.
    logical                        :: bc_mirror_pec     = .false.
    logical                        :: bc_mirror_pmc     = .false.
    logical                        :: bc_periodic       = .false.
    logical                        :: bc_plane_waves    = .false.
    logical                        :: bc_medium         = .false.

    logical                        :: plane_waves                = .false.
    logical                        :: plane_waves_apply          = .false.
    logical                        :: spatial_constant           = .false.
    logical                        :: spatial_constant_apply     = .false.
    logical                        :: spatial_constant_propagate = .false.

    ! TODO: add medium object file
    integer                        :: medium_calculation

    logical                        :: calc_medium_box = .false.
    type(single_medium_box_t), allocatable  :: medium_boxes(:)
    logical                         :: medium_boxes_initialized = .false.

    !> maxwell hamiltonian_mxll
    integer                        :: operator
    logical                        :: current_density_ext_flag = .false.

    type(energy_mxll_t)            :: energy

    logical                        :: cpml_hamiltonian = .false.

    logical                        :: diamag_current = .false.

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

  integer, public, parameter ::      &
    FARADAY_AMPERE_OLD          = 0, &
    FARADAY_AMPERE              = 1, &
    FARADAY_AMPERE_MEDIUM       = 2, &
    FARADAY_AMPERE_GAUSS        = 3, &
    FARADAY_AMPERE_GAUSS_MEDIUM = 4

contains

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

    hm%dim = st%dim
    hm%st => st

    ASSERT(associated(gr%der%lapl))

    hm%operators(1:3) => gr%der%grad(1:3) ! cross product for Maxwell calculation needs dimension >= 2
    hm%der => gr%der
    hm%rs_sign = st%rs_sign

    SAFE_ALLOCATE(hm%vector_potential(1:gr%mesh%np_part,1:st%dim))
    call energy_mxll_allocate(hm%energy, gr%mesh%np_part)

    hm%vector_potential = M_ZERO

    hm%mx_ma_coupling_apply = .false.
    hm%mx_ma_coupling  = .false.
    hm%ma_mx_coupling_apply = .false.
    hm%ma_mx_coupling  = .false.

    !%Variable MaxwellHamiltonianOperator
    !%Type integer
    !%Default faraday_ampere
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
    call parse_variable(namespace, 'MaxwellHamiltonianOperator', FARADAY_AMPERE, hm%operator)

    if (hm%operator == FARADAY_AMPERE_GAUSS) then
      hm%dim = hm%dim+1
    else if (hm%operator == FARADAY_AMPERE_MEDIUM) then
      ! TODO: define this operator automatically if there is a linear medium system present 
      hm%dim = 2*hm%dim
      hm%calc_medium_box = .true.
    end if

    !%Variable ExternalCurrent
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% If an external current density will be used.
    !%End
    call parse_variable(namespace, 'ExternalCurrent', .false., hm%current_density_ext_flag)

    hm%plane_waves_apply = .false.
    hm%spatial_constant_apply = .false.
    hm%spatial_constant_propagate = .true. ! only used if spatially constant field is used

    hm%propagation_apply = .false.

    !%Variable MaxwellMediumCalculation
    !%Type integer
    !%Default RS
    !%Section Hamiltonian
    !%Description
    !% For linear media the calculation of the Maxwell Operator acting on the RS state can be done
    !% directly using the Riemann-Silberstein representation or by calculating the curl of the 
    !% electric and magnetic fields.
    !%Option RS 1
    !% Medium calculation directly via Hamiltonian
    !%Option EM 2
    !% Medium calculation via curl of electric field and magnetic field
    !%End
    default_propagator = OPTION__MAXWELLMEDIUMCALCULATION__RS
    call parse_variable(namespace, 'MaxwellMediumCalculation', default_propagator, hm%medium_calculation)
    if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__EM) then
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

    type(profile_t), save :: prof
    integer :: il

    PUSH_SUB(hamiltonian_mxll_end)

    call profiling_in(prof, "HAMILTONIAN_MXLL_END")

    nullify(hm%operators)

    SAFE_DEALLOCATE_A(hm%vector_potential)
    call energy_mxll_end(hm%energy)

    call bc_mxll_end(hm%bc)

    if (allocated(hm%medium_boxes)) then
      do il = 1, size(hm%medium_boxes)
        call single_medium_box_end(hm%medium_boxes(il))
      end do
    end if

    call profiling_out(prof)

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

    type(profile_t), save :: prof
    type(batch_t), allocatable :: gradb(:)
    integer :: idir, ifield, field_dir, pml_dir, rs_sign
    integer :: ip, ip_in, il
    FLOAT :: pml_c, grad_real, grad_imag
    CMPLX :: pml_a, pml_b, pml_g, grad
    integer, parameter :: field_dirs(3, 2) = reshape([2, 3, 1, 3, 1, 2], [3, 2])
    FLOAT :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m, ff_real(3), ff_imag(3), coeff_real, coeff_imag
    CMPLX :: ff_plus(3), ff_minus(3), hpsi(6)
    integer :: sign_medium(6) = [1, 1, 1, -1, -1, -1]
    logical :: with_medium

    PUSH_SUB(hamiltonian_mxll_apply_batch)
    call profiling_in(prof, "H_MXLL_APPLY_BATCH")

    ASSERT(psib%status() == hpsib%status())

    ASSERT(psib%nst == hpsib%nst)
    ASSERT(hm%st%dim == 3)

    !Not implemented at the moment
    ASSERT(.not.present(terms))
    with_medium = hm%operator == FARADAY_AMPERE_MEDIUM

    if(present(time)) then
      if(abs(time - hm%current_time) > CNST(1e-10)) then
        write(message(1),'(a)') 'hamiltonian_apply_batch time assertion failed.'
        write(message(2),'(a,f12.6,a,f12.6)') 'time = ', time, '; hm%current_time = ', hm%current_time
        call messages_fatal(2, namespace=namespace)
      endif
    end if

    SAFE_ALLOCATE(gradb(1:der%dim))
    do idir = 1, der%dim
      call psib%copy_to(gradb(idir))
    end do
    call zderivatives_batch_grad(der, psib, gradb, set_bc=set_bc)

    if (hm%cpml_hamiltonian) then
      call apply_pml_boundary()
    end if

    call zderivatives_batch_curl(der, hpsib, gradb=gradb)

    call scale_after_curl()

    if (hm%bc_constant .and. .not. with_medium) then
      call apply_constant_boundary()
    end if

    if (with_medium) then
      do idir = 1, 3
        if ((hm%bc%bc_type(idir) == MXLL_BC_MEDIUM) .and. &
            (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS)) then
          call apply_medium_box(hm%bc%medium(idir))
        end if
      end do

      if (hm%calc_medium_box .and. &
           (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) ) then
        do il = 1, size(hm%medium_boxes)
          call apply_medium_box(hm%medium_boxes(il))
        end do
      end if
    end if

    do idir = 1, der%dim
      call gradb(idir)%end()
    end do
    SAFE_DEALLOCATE_A(gradb)

    call profiling_out(prof)
    POP_SUB(hamiltonian_mxll_apply_batch)

    contains
      subroutine apply_pml_boundary
        type(profile_t), save :: prof_pml
        PUSH_SUB(hamiltonian_mxll_apply_batch.apply_pml_boundary)
        call profiling_in(prof_pml, "APPLY_PML_BOUNDARY")
        if (with_medium) then
          rs_sign = 1
        else
          rs_sign = hm%rs_sign
        end if
        do idir = 1, der%dim
          call batch_scal(der%mesh%np, rs_sign * P_c, gradb(idir))
        end do

        do pml_dir = 1, hm%st%dim
          if (hm%bc%bc_ab_type(pml_dir) == MXLL_AB_CPML) then
            select case(gradb(pml_dir)%status())
            case(BATCH_NOT_PACKED)
              do ifield = 1, 2
                field_dir = field_dirs(pml_dir, ifield)
                !$omp parallel do private(ip, pml_c, pml_a, pml_b, pml_g, grad, grad_real, grad_imag)
                do ip_in = 1, hm%bc%pml%points_number
                  ip = hm%bc%pml%points_map(ip_in)
                  pml_c = hm%bc%pml%c(ip_in, pml_dir)
                  pml_a = hm%bc%pml%a(ip_in, pml_dir)
                  pml_b = hm%bc%pml%b(ip_in, pml_dir)
                  pml_g = hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir)
                  grad = gradb(pml_dir)%zff_linear(ip, field_dir)
                  grad_real = pml_c * ((M_ONE + TOFLOAT(pml_a))*TOFLOAT(grad)/P_c &
                       + rs_sign * TOFLOAT(pml_b)*TOFLOAT(pml_g))
                  grad_imag = pml_c * ((M_ONE + aimag(pml_a)) * aimag(grad)/P_c &
                       + rs_sign * aimag(pml_b) * aimag(pml_g))
                  gradb(pml_dir)%zff_linear(ip, field_dir) = TOCMPLX(grad_real, grad_imag)
                  if (with_medium) then
                    pml_g = hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir)
                    grad = gradb(pml_dir)%zff_linear(ip, field_dir+3)
                    grad_real = pml_c * ((M_ONE + TOFLOAT(pml_a))*TOFLOAT(grad)/P_c &
                         + rs_sign * TOFLOAT(pml_b)*TOFLOAT(pml_g))
                    grad_imag = pml_c * ((M_ONE + aimag(pml_a)) * aimag(grad)/P_c &
                         + rs_sign * aimag(pml_b) * aimag(pml_g))
                    gradb(pml_dir)%zff_linear(ip, field_dir+3) = TOCMPLX(grad_real, grad_imag)
                  end if
                end do
              end do
            case(BATCH_PACKED)
              !$omp parallel do private(ip, field_dir, pml_c, pml_a, pml_b, pml_g, grad, grad_real, grad_imag)
              do ip_in = 1, hm%bc%pml%points_number
                ip = hm%bc%pml%points_map(ip_in)
                pml_c = hm%bc%pml%c(ip_in, pml_dir)
                pml_a = hm%bc%pml%a(ip_in, pml_dir)
                pml_b = hm%bc%pml%b(ip_in, pml_dir)
                do ifield = 1, 2
                  field_dir = field_dirs(pml_dir, ifield)
                  pml_g = hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir)
                  grad = gradb(pml_dir)%zff_pack(field_dir, ip)
                  grad_real = pml_c * ((M_ONE + TOFLOAT(pml_a))*TOFLOAT(grad)/P_c &
                       + rs_sign * TOFLOAT(pml_b)*TOFLOAT(pml_g))
                  grad_imag = pml_c * ((M_ONE + aimag(pml_a)) * aimag(grad)/P_c &
                       + rs_sign * aimag(pml_b) * aimag(pml_g))
                  gradb(pml_dir)%zff_pack(field_dir, ip) = TOCMPLX(grad_real, grad_imag)
                  if (with_medium) then
                    pml_g = hm%bc%pml%conv_minus(ip_in, pml_dir, field_dir)
                    grad = gradb(pml_dir)%zff_pack(field_dir+3, ip)
                    grad_real = pml_c * ((M_ONE + TOFLOAT(pml_a))*TOFLOAT(grad)/P_c &
                         + rs_sign * TOFLOAT(pml_b)*TOFLOAT(pml_g))
                    grad_imag = pml_c * ((M_ONE + aimag(pml_a)) * aimag(grad)/P_c &
                         + rs_sign * aimag(pml_b) * aimag(pml_g))
                    gradb(pml_dir)%zff_pack(field_dir+3, ip) = TOCMPLX(grad_real, grad_imag)
                  end if
                end do
              end do
            case(BATCH_DEVICE_PACKED)
              call messages_not_implemented("Maxwell PML on GPU")
            end select
          end if
        end do
        call profiling_out(prof_pml)
        POP_SUB(hamiltonian_mxll_apply_batch.apply_pml_boundary)
      end subroutine apply_pml_boundary

      subroutine scale_after_curl
        type(profile_t), save :: prof_scale
        PUSH_SUB(hamiltonian_mxll_apply_batch.scale_after_curl)
        call profiling_in(prof_scale, "SCALE_AFTER_CURL")
        if (.not. hm%cpml_hamiltonian) then
          ! if we do not need pml, scale after the curl because it is cheaper
          if (with_medium) then
            ! in case of a medium, multiply first 3 components with +, others with -
            call batch_scal(der%mesh%np, sign_medium * P_c, hpsib)
          else
            call batch_scal(der%mesh%np, hm%rs_sign * P_c, hpsib)
          end if
        else
          ! this is needed for PML computations with medium
          if (with_medium) then
            ! in case of a medium, multiply first 3 components with +, others with -
            call batch_scal(der%mesh%np, sign_medium * M_ONE, hpsib)
          end if
        end if
        call profiling_out(prof_scale)
        POP_SUB(hamiltonian_mxll_apply_batch.scale_after_curl)
      end subroutine scale_after_curl

      subroutine apply_constant_boundary
        type(profile_t), save :: prof_bc_const
        PUSH_SUB(hamiltonian_mxll_apply_batch.apply_constant_boundary)
        call profiling_in(prof_bc_const, 'APPLY_CONSTANT_BC')
        select case(hpsib%status())
        case(BATCH_NOT_PACKED)
          do field_dir = 1, hm%st%dim
            do ip_in = 1, hm%bc%constant_points_number
              ip = hm%bc%constant_points_map(ip_in)
              hpsib%zff_linear(ip, field_dir) = hm%st%rs_state_const(field_dir)
            end do
          end do
        case(BATCH_PACKED)
          do ip_in = 1, hm%bc%constant_points_number
            ip = hm%bc%constant_points_map(ip_in)
            do field_dir = 1, hm%st%dim
              hpsib%zff_pack(field_dir, ip) = hm%st%rs_state_const(field_dir)
            end do
          end do
        end select
        call profiling_out(prof_bc_const)
        POP_SUB(hamiltonian_mxll_apply_batch.apply_constant_boundary)
      end subroutine apply_constant_boundary

      subroutine apply_medium_box(medium)
        type(single_medium_box_t),  intent(in) :: medium

        integer :: ifield
        type(profile_t), save :: prof_medium_box

        PUSH_SUB(hamiltonian_mxll_apply_batch.apply_medium_box)
        call profiling_in(prof_medium_box, "MEDIUM_BOX")
        !$omp parallel do private(ip, cc, aux_ep, aux_mu, sigma_e, sigma_m, &
        !$omp ff_plus, ff_minus, hpsi, ff_real, ff_imag, ifield, coeff_real, coeff_imag)
        do ip_in = 1, medium%points_number
          ip          = medium%points_map(ip_in)
          cc          = medium%c(ip_in)/P_c
          aux_ep(1:3) = medium%aux_ep(ip_in, 1:3)
          aux_mu(1:3) = medium%aux_mu(ip_in, 1:3)
          sigma_e     = medium%sigma_e(ip_in)
          sigma_m     = medium%sigma_m(ip_in)
          select case(hpsib%status())
          case(BATCH_NOT_PACKED)
            ff_plus(1:3)  = psib%zff_linear(ip, 1:3)
            ff_minus(1:3) = psib%zff_linear(ip, 4:6)
            hpsi(1:6) = hpsib%zff_linear(ip, 1:6)
          case(BATCH_PACKED)
            ff_plus(1:3)  = psib%zff_pack(1:3, ip)
            ff_minus(1:3) = psib%zff_pack(4:6, ip)
            hpsi(1:6) = hpsib%zff_pack(1:6, ip)
          end select
          ff_real = TOFLOAT(ff_plus+ff_minus)
          ff_imag = aimag(ff_plus-ff_minus)
          aux_ep = dcross_product(aux_ep, ff_real)
          aux_mu = dcross_product(aux_mu, ff_imag)
          do ifield = 1, 3
            coeff_real = - cc * aux_ep(ifield) + sigma_m * ff_imag(ifield)
            coeff_imag = - cc * aux_mu(ifield) - sigma_e * ff_real(ifield)
            hpsi(ifield) = cc * hpsi(ifield) + TOCMPLX(coeff_real, coeff_imag)
            hpsi(ifield+3) = cc * hpsi(ifield+3) + TOCMPLX(-coeff_real, coeff_imag)
          end do
          select case(hpsib%status())
          case(BATCH_NOT_PACKED)
            hpsib%zff_linear(ip, 1:6) = hpsi(1:6)
          case(BATCH_PACKED)
            hpsib%zff_pack(1:6, ip) = hpsi(1:6)
          end select
        end do
        call profiling_out(prof_medium_box)
        POP_SUB(hamiltonian_mxll_apply_batch.apply_medium_box)
      end subroutine apply_medium_box

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
    type(profile_t), save :: prof
    logical :: pml_and_gpu, on_gpu

    PUSH_SUB(zhamiltonian_mxll_apply)

    call profiling_in(prof, 'ZHAMILTONIAN_MXLL_APPLY')


    on_gpu = psib%status() == BATCH_DEVICE_PACKED
    pml_and_gpu = hm%cpml_hamiltonian .and. on_gpu
    if ((hm%operator == FARADAY_AMPERE .and. .not. pml_and_gpu) .or. &
      (hm%operator == FARADAY_AMPERE_MEDIUM .and. .not. on_gpu)) then
      call hamiltonian_mxll_apply_batch(hm, namespace, hm%der, psib, hpsib, set_bc=set_bc)
    else
      SAFE_ALLOCATE(rs_aux_in(1:hm%der%mesh%np_part, 1:hm%dim))
      SAFE_ALLOCATE(rs_aux_out(1:hm%der%mesh%np, 1:hm%dim))
      call boundaries_set(hm%der%boundaries, psib)
      do ii = 1, hm%dim
        call batch_get_state(psib, ii, hm%der%mesh%np_part, rs_aux_in(:, ii))
      end do
      ! This uses the old non-batch implementation
      call maxwell_hamiltonian_apply_fd(hm, hm%der, rs_aux_in, rs_aux_out)
      do ii = 1, hm%dim
        call batch_set_state(hpsib, ii, hm%der%mesh%np, rs_aux_out(:, ii))
      end do
      SAFE_DEALLOCATE_A(rs_aux_in)
      SAFE_DEALLOCATE_A(rs_aux_out)
      
    end if

    call profiling_out(prof)

    POP_SUB(zhamiltonian_mxll_apply)
  end subroutine zhamiltonian_mxll_apply


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
    type(profile_t), save :: prof, prof_method

    PUSH_SUB(maxwell_hamiltonian_apply_fd)

    call profiling_in(prof, 'MAXWELL_HAMILTONIAN_APPLY_FD')

    np = der%mesh%np
    np_part = der%mesh%np_part
    rs_sign = hm%rs_sign


    select case (hm%operator)

    !=================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum via partial derivatives:

    case (FARADAY_AMPERE)
      call profiling_in(prof_method, 'MXLL_HAM_APPLY_FD_FARADAY_AMP')

      SAFE_ALLOCATE(tmp(np,2))
      oppsi       = M_z0

      if (hm%diamag_current) then
        mx_rho    => hm%st%grid_rho
        kappa_psi => hm%st%kappa_psi 
      end if

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 3, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 2, 3, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 3, 2, tmp(:,2))
      oppsi(1:np,1) = ( tmp(1:np,1)-tmp(1:np,2) )

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 1, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 3, 1, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 1, 3, tmp(:,2))
      oppsi(1:np,2) = ( tmp(1:np,1)-tmp(1:np,2) )

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 2, set_bc = .false.)
      tmp = rs_sign * P_c * tmp
      call maxwell_pml_hamiltonian(hm, der, psi, 1, 2, tmp(:,1))
      call maxwell_pml_hamiltonian(hm, der, psi, 2, 1, tmp(:,2))
      oppsi(1:np,3) = ( tmp(1:np,1)-tmp(1:np,2) )

      if (hm%bc_constant) then
        do ip_in = 1, hm%bc%constant_points_number
          ip = hm%bc%constant_points_map(ip_in)
          oppsi(ip,:) = hm%st%rs_state_const(:)
        end do
      end if

      SAFE_DEALLOCATE_A(tmp)

      call profiling_out(prof_method)
    !=================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium via partial derivatives:

    case (FARADAY_AMPERE_MEDIUM)
      call profiling_in(prof_method, 'MXLL_HAM_APP_FAR_AMP_MED')

      SAFE_ALLOCATE(tmp(np,4))
      oppsi       = M_z0

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 1 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,3), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,3), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,4), 3, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 2, 3, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 3, 2, tmp(:,3:4))
      oppsi(1:np,1) =   ( tmp(1:np,1)-tmp(1:np,3) )
      oppsi(1:np,4) = - ( tmp(1:np,2)-tmp(1:np,4) )

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 2 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,4), 1, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 3, 1, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 1, 3, tmp(:,3:4))
      oppsi(1:np,2) =   ( tmp(1:np,1)-tmp(1:np,3) )
      oppsi(1:np,5) = - ( tmp(1:np,2)-tmp(1:np,4) )

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector component 3 calculation:
      tmp = M_z0
      call zderivatives_partial(der, psi(:,2), tmp(:,1), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,3), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,4), 2, set_bc = .false.)
      tmp = P_c * tmp
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 1, 2, tmp(:,1:2))
      call maxwell_pml_hamiltonian_medium(hm, der, psi, 2, 1, tmp(:,3:4))
      oppsi(1:np,3) =   ( tmp(1:np,1)-tmp(1:np,3) )
      oppsi(1:np,6) = - ( tmp(1:np,2)-tmp(1:np,4) )


      SAFE_DEALLOCATE_A(tmp)

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation if medium boundaries is set:
      call maxwell_medium_boundaries_calculation(hm, psi, oppsi)

      !-----------------------------------------------------------------------------------------------
      ! Riemann-Silberstein vector calculation for medium boxes:
      call maxwell_medium_boxes_calculation(hm, der, psi, oppsi)

      call profiling_out(prof_method)
    !=================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in vacuum with Gauss condition via partial derivatives:

    case (FARADAY_AMPERE_GAUSS)

      call profiling_in(prof_method, 'MXLL_HAM_APP_FAR_AMP_GAUSS')

      SAFE_ALLOCATE(tmp(np,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,1) = rs_sign * P_c * ( - M_zI*tmp(1:np,1) - tmp(1:np,2) - M_zI*tmp(1:np,3) )

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,2) = rs_sign * P_c * ( - M_zI*tmp(1:np,1) - tmp(1:np,2) - M_zI*tmp(1:np,3) )

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,3) = rs_sign * P_c * ( tmp(1:np,1) - M_zI*tmp(1:np,2) + M_zI*tmp(1:np,3) )

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,4) = rs_sign * P_c * ( tmp(1:np,1) - M_zI*tmp(1:np,2) + M_zI*tmp(1:np,3) )

      SAFE_DEALLOCATE_A(tmp)
      call profiling_out(prof_method)

    !=================================================================================================
    ! Maxwell Hamiltonian - Hamiltonian operation in medium with Gauss condition via partial derivatives:

    case (FARADAY_AMPERE_GAUSS_MEDIUM)

      call profiling_in(prof_method, 'MXLL_HAM_AP_FAR_AMP_GA_MED')
      SAFE_ALLOCATE(tmp(np,3))
      oppsi       = M_z0
      tmp = M_z0

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,1) = P_c*(-M_zI*tmp(1:np,1)-tmp(1:np,2)-M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,2) = P_c*(-M_zI*tmp(1:np,1)-tmp(1:np,2)-M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,1), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,1), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,3), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,3) = P_c*(tmp(1:np,1)-M_zI*tmp(1:np,2)+M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,2), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,2), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,4), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,4) = P_c*(tmp(1:np,1)-M_zI*tmp(1:np,2)+M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,5), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,5) = - P_c*(-M_zI*tmp(1:np,1)-tmp(1:np,2)-M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,6), tmp(:,1), 3, set_bc = .false.)
      call zderivatives_partial(der, psi(:,8), tmp(:,2), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,8), tmp(:,3), 1, set_bc = .false.)
      oppsi(1:np,6) = - P_c*(-M_zI*tmp(1:np,1)-tmp(1:np,2)-M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,5), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,5), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,7) = - P_c*(tmp(1:np,1)-M_zI*tmp(1:np,2)+M_zI*tmp(1:np,3))

      call zderivatives_partial(der, psi(:,6), tmp(:,1), 2, set_bc = .false.)
      call zderivatives_partial(der, psi(:,6), tmp(:,2), 1, set_bc = .false.)
      call zderivatives_partial(der, psi(:,7), tmp(:,3), 3, set_bc = .false.)
      oppsi(1:np,8) = - P_c*(tmp(1:np,1)-M_zI*tmp(1:np,2)+M_zI*tmp(1:np,3))

      SAFE_DEALLOCATE_A(tmp)
      call profiling_out(prof_method)

    end select

    call profiling_out(prof)

    POP_SUB(maxwell_hamiltonian_apply_fd)
  end subroutine maxwell_hamiltonian_apply_fd


  ! ---------------------------------------------------------
  !> Maxwell Hamiltonian is updated for the PML calculation
  subroutine maxwell_pml_hamiltonian(hm, der, psi, dir1, dir2, tmp)
    type(hamiltonian_mxll_t), intent(in)    :: hm
    type(derivatives_t),      intent(in)    :: der
    CMPLX,                    intent(inout) :: psi(:,:)
    integer,                  intent(in)    :: dir1
    integer,                  intent(in)    :: dir2
    CMPLX,                    intent(inout) :: tmp(:)

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_pml_hamiltonian)

    call profiling_in(prof, 'MAXWELL_PML_HAMILTONIAN')

    if ( (hm%bc%bc_ab_type(dir1) == MXLL_AB_CPML) .and. hm%cpml_hamiltonian ) then
      if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) then
        call maxwell_pml_calculation_via_riemann_silberstein(hm, der, psi, dir1, dir2, tmp(:))
      end if
    end if

   call profiling_out(prof)

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

    type(profile_t), save :: prof

    PUSH_SUB(maxwell_pml_hamiltonian_medium)

    call profiling_in(prof, 'MAXWELL_PML_HAMILTONIAN_MEDIUM')

    if ( (hm%bc%bc_ab_type(dir1) == MXLL_AB_CPML) .and. hm%cpml_hamiltonian ) then
      if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) then
        call maxwell_pml_calculation_via_riemann_silberstein_medium(hm, der, psi, dir1, dir2, tmp(:,:))
!      else if (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__EM) then
!        call maxwell_pml_calculation_via_e_b_fields_medium(hm, der, psi, dir1, dir2, tmp(:,:))
      end if
    end if

    call profiling_out(prof)

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

    integer            :: ip, ip_in, rs_sign
    FLOAT              :: pml_c
    CMPLX, allocatable :: tmp_partial(:)
    CMPLX              :: pml_a, pml_b, pml_g

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein)

    if (hm%cpml_hamiltonian) then

      rs_sign = hm%rs_sign

      SAFE_ALLOCATE(tmp_partial(der%mesh%np_part))

      call zderivatives_partial(der, psi(:,field_dir), tmp_partial(:), pml_dir, set_bc = .false.)
      do ip_in = 1, hm%bc%pml%points_number
        ip       = hm%bc%pml%points_map(ip_in)
        pml_c = hm%bc%pml%c(ip_in, pml_dir)
        pml_a = hm%bc%pml%a(ip_in, pml_dir)
        pml_b = hm%bc%pml%b(ip_in, pml_dir)
        pml_g = hm%bc%pml%conv_plus(ip_in, pml_dir, field_dir)
        pml(ip)  = rs_sign * pml_c * ( tmp_partial(ip) &
                 +  TOFLOAT(pml_a) * TOFLOAT(tmp_partial(ip)) &
                 +  M_zI * aimag(pml_a) * aimag(tmp_partial(ip)) &
                 +  TOFLOAT(pml_b) * TOFLOAT(pml_g) &
                 +  M_zI * aimag(pml_b) * aimag(pml_g))
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

    integer            :: ip, ip_in, np
    FLOAT              :: pml_c(3)
    CMPLX, allocatable :: tmp_partial(:,:)
    CMPLX              :: pml_a(3), pml_b(3), pml_g_p(3), pml_g_m(3)

    PUSH_SUB(maxwell_pml_calculation_via_riemann_silberstein_medium)

    if (hm%cpml_hamiltonian) then

      np = der%mesh%np
      SAFE_ALLOCATE(tmp_partial(np,1:2))

      call zderivatives_partial(der, psi(:,field_dir  ), tmp_partial(:,1), pml_dir, set_bc = .false.)
      call zderivatives_partial(der, psi(:,field_dir+3), tmp_partial(:,2), pml_dir, set_bc = .false.)
      do ip_in = 1, hm%bc%pml%points_number
        ip         = hm%bc%pml%points_map(ip_in)
        pml_c(1:3)   = hm%bc%pml%c(ip_in, 1:3)
        pml_a(1:3)   = hm%bc%pml%a(ip_in, 1:3)
        pml_b(1:3)   = hm%bc%pml%b(ip_in, 1:3)
        pml_g_p(1:3) = hm%bc%pml%conv_plus(ip_in, pml_dir, 1:3)
        pml_g_m(1:3) = hm%bc%pml%conv_minus(ip_in, pml_dir, 1:3)
        pml(ip, 1)   = pml_c(pml_dir) * tmp_partial(ip, 1) &
                     + pml_c(pml_dir) * TOFLOAT(pml_a(pml_dir)) * TOFLOAT(tmp_partial(ip, 1)) &
                     + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip, 1)) &
                     + pml_c(pml_dir) * TOFLOAT(pml_b(pml_dir)) * TOFLOAT(pml_g_p(field_dir)) &
                     + M_zI * pml_c(pml_dir) * aimag(pml_b(pml_dir)) * aimag(pml_g_p(field_dir))
        pml(ip, 2)   = pml_c(pml_dir) * tmp_partial(ip, 2) &
                     + pml_c(pml_dir) * TOFLOAT(pml_a(pml_dir)) * TOFLOAT(tmp_partial(ip, 2)) &
                     + M_zI * pml_c(pml_dir) * aimag(pml_a(pml_dir)) * aimag(tmp_partial(ip, 2)) &
                     + pml_c(pml_dir) * TOFLOAT(pml_b(pml_dir)) * TOFLOAT(pml_g_m(field_dir)) &
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
      if ( (hm%bc%bc_type(idim) == MXLL_BC_MEDIUM) .and. &
           (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) ) then
        do ip_in = 1, hm%bc%medium(idim)%points_number
          ip          = hm%bc%medium(idim)%points_map(ip_in)
          cc          = hm%bc%medium(idim)%c(ip_in)/P_c
          aux_ep(:)   = hm%bc%medium(idim)%aux_ep(ip_in, :)
          aux_mu(:)   = hm%bc%medium(idim)%aux_mu(ip_in, :)
          sigma_e     = hm%bc%medium(idim)%sigma_e(ip_in)
          sigma_m     = hm%bc%medium(idim)%sigma_m(ip_in)
          ff_plus(1)  = psi(ip, 1)
          ff_plus(2)  = psi(ip, 2)
          ff_plus(3)  = psi(ip, 3)
          ff_minus(1) = psi(ip, 4)
          ff_minus(2) = psi(ip, 5)
          ff_minus(3) = psi(ip, 6)
          aux_ep      = dcross_product(aux_ep,TOFLOAT(ff_plus+ff_minus))
          aux_mu      = dcross_product(aux_mu,aimag(ff_plus-ff_minus))
          oppsi(ip, 1) = oppsi(ip, 1)*cc                                         &
                       - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(1) + ff_minus(1))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 4) = oppsi(ip, 4)*cc                                         &
                       + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(1) + ff_minus(1))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 2) = oppsi(ip, 2)*cc                                         &
                       - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(2) + ff_minus(2))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2))
          oppsi(ip, 5) = oppsi(ip, 5)*cc                                         &
                       + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(2) + ff_minus(2))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2))
          oppsi(ip, 3) = oppsi(ip, 3)*cc                                         &
                       - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(3) + ff_minus(3))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
          oppsi(ip, 6) = oppsi(ip, 6)*cc                                         &
                       + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(3) + ff_minus(3))         &
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

    integer            :: ip, ip_in, il
    FLOAT              :: cc, aux_ep(3), aux_mu(3), sigma_e, sigma_m
    CMPLX              :: ff_plus(3), ff_minus(3)

    PUSH_SUB(maxwell_medium_boxes_calculation)

    if (hm%calc_medium_box .and. &
         (hm%medium_calculation == OPTION__MAXWELLMEDIUMCALCULATION__RS) ) then
      do il = 1, size(hm%medium_boxes)
        do ip_in = 1, hm%medium_boxes(il)%points_number
          ip           = hm%medium_boxes(il)%points_map(ip_in)
          cc           = hm%medium_boxes(il)%c(ip_in)/P_c
          aux_ep(1:3)  = hm%medium_boxes(il)%aux_ep(ip_in, 1:3)
          aux_mu(1:3)  = hm%medium_boxes(il)%aux_mu(ip_in, 1:3)
          sigma_e      = hm%medium_boxes(il)%sigma_e(ip_in)
          sigma_m      = hm%medium_boxes(il)%sigma_m(ip_in)
          ff_plus(1)   = psi(ip, 1)
          ff_plus(2)   = psi(ip, 2)
          ff_plus(3)   = psi(ip, 3)
          ff_minus(1)  = psi(ip, 4)
          ff_minus(2)  = psi(ip, 5)
          ff_minus(3)  = psi(ip, 6)
          aux_ep       = dcross_product(aux_ep, TOFLOAT(ff_plus+ff_minus))
          aux_mu       = dcross_product(aux_mu, aimag(ff_plus-ff_minus))
          oppsi(ip, 1) = oppsi(ip,1)*cc                                          &
                       - cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(1) + ff_minus(1))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 4) = oppsi(ip,4)*cc                                          &
                       + cc * aux_ep(1) - cc * M_zI * aux_mu(1)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(1) + ff_minus(1))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(1) - ff_minus(1))
          oppsi(ip, 2) = oppsi(ip,2)*cc                                          &
                       - cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(2) + ff_minus(2))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2))
          oppsi(ip, 5) = oppsi(ip,5)*cc                                          &
                       + cc * aux_ep(2) - cc * M_zI * aux_mu(2)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(2) + ff_minus(2))         &
                       + M_zI * sigma_m * M_zI * aimag(ff_plus(2) - ff_minus(2)) 
          oppsi(ip, 3) = oppsi(ip,3)*cc                                          &
                       - cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(3) + ff_minus(3))         &
                       - M_zI * sigma_m * M_zI * aimag(ff_plus(3) - ff_minus(3))
          oppsi(ip, 6) = oppsi(ip,6)*cc                                          &
                       + cc * aux_ep(3) - cc * M_zI * aux_mu(3)                  &
                       - M_zI * sigma_e * TOFLOAT(ff_plus(3) + ff_minus(3))         &
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

    integer            :: idim, ip, ip_in
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
    do idim = 1, st%dim
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

    ! here we could add surface integral correction from Renes code

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

    PUSH_SUB(maxwell_helmholtz_decomposition_long_field)

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

    POP_SUB(maxwell_helmholtz_decomposition_long_field)

  end subroutine maxwell_helmholtz_decomposition_long_field

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

    call messages_not_implemented ('dhamiltonian_mxll_magnus_apply', namespace=namespace)

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

    call messages_not_implemented ('zhamiltonian_mxll_magnus_apply', namespace=namespace)

  end subroutine zhamiltonian_mxll_magnus_apply

end module hamiltonian_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
