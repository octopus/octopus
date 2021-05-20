!! Copyright (C) 2020 F. Bonafé, H. Appel, R. Jestädt
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

module system_linear_medium_oct_m
  use algorithm_oct_m
  use clock_oct_m
  use cgal_polyhedra_oct_m
  use derivatives_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use iso_c_binding
  use linear_medium_em_field_oct_m
  use messages_oct_m
  use comm_oct_m
  use grid_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use propagator_exp_mid_oct_m
  use propagator_oct_m
  use profiling_oct_m
  use quantity_oct_m
  use system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::           &
    system_linear_medium_t,    &
    system_linear_medium_init

  integer, parameter :: &
    MEDIUM_PARALLELEPIPED = 1,         &
    MEDIUM_BOX_FILE       = 2

  type, extends(system_t) :: system_linear_medium_t
     integer            :: box_shape
     integer            :: edge_profile  !< edge shape profile (smooth or steep)
     FLOAT              :: center(3) !< center of a box
     FLOAT              :: lsize(3)  !< length in each direction of a box
     FLOAT              :: ep_factor !< permitivity before applying edge profile
     FLOAT              :: mu_factor !< permeability before applying edge profile
     FLOAT              :: sigma_e_factor !< electric conductivy before applying edge profile
     FLOAT              :: sigma_m_factor !< magnetic conductivity before applying edge4 profile
     character(len=256) :: filename
     logical            :: check_medium_points = .false.

  contains
    procedure :: init_interaction => system_linear_medium_init_interaction
    procedure :: init_interaction_as_partner => system_linear_medium_init_interaction_as_partner
    procedure :: initial_conditions => system_linear_medium_initial_conditions
    procedure :: do_td_operation => system_linear_medium_do_td
    procedure :: iteration_info => system_linear_medium_iteration_info
    procedure :: output_start => system_linear_medium_output_start
    procedure :: output_write => system_linear_medium_output_write
    procedure :: output_finish => system_linear_medium_output_finish
    procedure :: is_tolerance_reached => system_linear_medium_is_tolerance_reached
    procedure :: update_quantity => system_linear_medium_update_quantity
    procedure :: update_exposed_quantity => system_linear_medium_update_exposed_quantity
    procedure :: copy_quantities_to_interaction => system_linear_medium_copy_quantities_to_interaction
    final :: system_linear_medium_finalize
  end type system_linear_medium_t

  interface system_linear_medium_t
    procedure system_linear_medium_constructor
  end interface system_linear_medium_t

contains

  ! ---------------------------------------------------------
  !> The factory routine (or constructor) allocates a pointer of the
  !! corresponding type and then calls the init routine which is a type-bound
  !! procedure of the corresponding type. With this design, also derived
  !! classes can use the init routine of the parent class.
  function system_linear_medium_constructor(namespace) result(sys)
    class(system_linear_medium_t), pointer    :: sys
    type(namespace_t),           intent(in) :: namespace

    PUSH_SUB(system_linear_medium_constructor)

    SAFE_ALLOCATE(sys)

    call system_linear_medium_init(sys, namespace)

    POP_SUB(system_linear_medium_constructor)
  end function system_linear_medium_constructor

  ! ---------------------------------------------------------
  !> The init routine is a module level procedure
  !! This has the advantage that different classes can have different
  !! signatures for the initialization routines because they are not
  !! type-bound and thus also not inherited.
  ! ---------------------------------------------------------
  subroutine system_linear_medium_init(this, namespace)
    class(system_linear_medium_t), target, intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace

    integer :: nlines, ncols, idim
    type(block_t) :: blk
    type(profile_t), save :: prof

    PUSH_SUB(system_linear_medium_init)

    this%namespace = namespace


    call profiling_in(prof, 'MEDIUM_BOX_INIT')

    !%Variable LinearMediumBoxShape
    !%Type integer
    !%Section Maxwell
    !%Description
    !% This variable defines the shape of the linear medium box.
    !% The default is <tt>medium_parallelepiped</tt>.
    !%Option medium_parallelepiped 1
    !% The medium box will be a parallelepiped whose center and dimensions are taken from
    !% the variable <tt>LinearMediumLsize</tt>.
    !%Option medium_box_file 2
    !% The simulation box will be read from an external file in OFF format, defined by the variable <tt>LinearMediumBoxFile</tt>.
    !%End
    call parse_variable(namespace, 'LinearMediumBoxShape', MEDIUM_PARALLELEPIPED, this%box_shape)

    if(.not.varinfo_valid_option('LinearMediumBoxShape', this%box_shape)) &
        & call messages_input_error(namespace, 'LinearMediumBoxShape')

    select case(this%box_shape)
    case (MEDIUM_PARALLELEPIPED)

      !%Variable LinearMediumBoxSize
      !%Type block
      !%Section Maxwell
      !%Description
      !% Defines center and size of a parallelepiped linear medium box.
      !%
      !% Example:
      !%
      !% <tt>%LinearMediumBox
      !% <br>&nbsp;&nbsp;   center_x | center_y | center_z | x_length | y_length | z_length
      !% <br>%</tt>
      !%End

      if(parse_block(namespace, 'LinearMediumBoxSize', blk) == 0) then
        call messages_print_stress(stdout, trim('Linear Medium box center and size:'))

        nlines = parse_block_n(blk)
        if (nlines /=  1) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of one line')
        end if
        ncols = parse_block_cols(blk, 0)
        if (ncols /= 6) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of six columns')
        end if
        do idim = 1, 3
          call parse_block_float(blk, 0, idim-1, this%center(idim))
          call parse_block_float(blk, 0, idim+2, this%lsize(idim))
        end do
        write(message(1),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', this%center(1), ' | ',&
            this%center(2), ' | ', this%center(3)
        write(message(2),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', this%lsize(1), ' | ', &
            this%lsize(2), ' | ', this%lsize(3)
        write(message(3),'(a)') ""
        call messages_info(3)
      call parse_block_end(blk)

      call messages_print_stress(stdout)
    else
      message(1) = "For parallelepiped box shapes, you must provide a LinearMediumBoxSize block."
      call messages_fatal(1, namespace=namespace)
    end if

  case (MEDIUM_BOX_FILE)
    !%Variable LinearMediumBoxFile
    !%Type string
    !%Section Maxwell
    !%Description
    !% File in OFF format with the shape of the linear medium.
    !%End
    if (parse_is_defined(namespace, 'LinearMediumBoxFile')) then
      call parse_variable(namespace, 'LinearMediumBoxFile', 'mediumboxfile', this%filename)
    else
      message(1) = "When using box_file as the box shape, you must provide a filename through the LinearMediumBoxFile variable."
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable CheckPointsMediumFromFile
    !%Type logical
    !%Default no
    !%Section Maxwell
    !%Description
    !% Whether to re-calculate the points map by artificially shrinking the coordinate system by a factor of
    !% 0.99 to check if the points inside the medium surface are properly detected. This works for only one
    !% medium surface which is centered in the origin of the coordinate system.
    !%End
    call parse_variable(namespace, 'CheckPointsMediumFromFile', .false., this%check_medium_points)

    end select

    !%Variable LinearMediumProperties
    !%Type block
    !%Section Maxwell
    !%Description
    !% Defines electromagnetic parameters for a linear medium box.
    !%
    !% Example:
    !%
    !% <tt>%LinearMediumProperties
    !% <br>&nbsp;&nbsp;   epsilon_factor | mu_factor | sigma_e | sigma_m
    !% <br>%</tt>
    !%
    !% Permittivity factor, permeability factor, electric conductivity and magnetic conductivity of the medium box.
    !%End

    if (parse_block(namespace, 'LinearMediumProperties', blk) == 0) then

      nlines = parse_block_n(blk)
      if (nlines /=  1) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of one line', nlines)
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 4) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of four columns')
      end if
      call parse_block_float(blk, 0, 0, this%ep_factor)
      call parse_block_float(blk, 0, 1, this%mu_factor)
      call parse_block_float(blk, 0, 2, this%sigma_e_factor)
      call parse_block_float(blk, 0, 3, this%sigma_m_factor)
      write(message(1),'(a,es9.2)') 'Box epsilon factor: ', this%ep_factor
      write(message(2),'(a,es9.2)') 'Box mu factor:      ', this%mu_factor
      write(message(3),'(a,es9.2)') 'Box electric sigma: ', this%sigma_e_factor
      write(message(4),'(a,es9.2)') 'Box magnetic sigma: ', this%sigma_m_factor
      write(message(5),'(a)') ""
      call messages_info(5)
      call parse_block_end(blk)

      call messages_print_stress(stdout)
    else
      message(1) = 'You must specify the properties of your linear medium through the LinearMediumProperties block.'
      call messages_fatal(1, namespace=namespace)
    end if

    !%Variable LinearMediumEdgeProfile
    !%Type integer
    !%Section Maxwell
    !%Description
    !% Defines the type of numerical approximation used for the derivatives at the edges of the medium box.
    !% Default is edged. When the box shape is read from file, only the edged profile is supported.
    !%
    !%Option edged 1
    !% Medium box edges are considered steep for derivatives.
    !%Option smooth 2
    !% Medium box edged and softened for derivatives.
    !%End
    call parse_variable(namespace, 'LinearMediumEdgeProfile', OPTION__LINEARMEDIUMEDGEPROFILE__EDGED, this%edge_profile)
    if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__EDGED) then
      write(message(1),'(a,a)')   'Box shape:          ', 'edged'
    else if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__SMOOTH) then
      write(message(1),'(a,a)')   'Box shape:          ', 'smooth'
    end if
    call messages_info(1)

    call this%supported_interactions_as_partner%add(LINEAR_MEDIUM_EM_FIELD)
    this%quantities(PERMITTIVITY)%required = .true.
    this%quantities(PERMEABILITY)%required = .true.
    this%quantities(E_CONDUCTIVITY)%required = .true.
    this%quantities(M_CONDUCTIVITY)%required = .true.

    call profiling_out(prof)

    POP_SUB(system_linear_medium_init)
  end subroutine system_linear_medium_init

  ! ---------------------------------------------------------
  subroutine system_linear_medium_init_interaction(this, interaction)
    class(system_linear_medium_t), target, intent(inout) :: this
    class(interaction_t),                intent(inout) :: interaction

    PUSH_SUB(system_linear_medium_init_interaction)

    select type (interaction)
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_init_interaction)
  end subroutine system_linear_medium_init_interaction

  ! ---------------------------------------------------------
  subroutine system_linear_medium_init_interaction_as_partner(partner, interaction)
    class(system_linear_medium_t),       intent(in)    :: partner
    class(interaction_t),                intent(inout) :: interaction

    integer :: n_points, n_global_points
    integer, allocatable :: tmp(:)

    PUSH_SUB(system_linear_medium_init_interaction_as_partner)

    select type (interaction)
    type is (linear_medium_em_field_t)
      call generate_medium_box(partner, interaction%medium_box, interaction%system_gr)
      interaction%allocated_partner_arrays = .true.

      if (partner%check_medium_points) then
        SAFE_ALLOCATE(tmp(interaction%system_gr%mesh%np))
        n_global_points = 0
        write(message(1),'(a, a, a)')   'Check of points inside surface of medium ', trim(partner%filename), ":"
        call messages_info(1)
        call get_points_map_from_file(partner%filename, interaction%system_gr%mesh, n_points, n_global_points, tmp, CNST(0.99))
        write(message(1),'(a, I8)')'Number of points inside medium (normal coordinates):', interaction%medium_box%global_points_number
        write(message(2),'(a, I8)')'Number of points inside medium (rescaled coordinates):', n_global_points
        write(message(3), '(a)') ""
        call messages_info(3)
        SAFE_DEALLOCATE_A(tmp)
      end if
    class default
      message(1) = "Trying to initialize an unsupported interaction by a linear medium."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_init_interaction_as_partner)
  end subroutine system_linear_medium_init_interaction_as_partner

  ! ---------------------------------------------------------
  subroutine system_linear_medium_initial_conditions(this, from_scratch)
    class(system_linear_medium_t), intent(inout) :: this
    logical,                 intent(in)    :: from_scratch

    PUSH_SUB(system_linear_medium_initial_conditions)

    POP_SUB(system_linear_medium_initial_conditions)
  end subroutine system_linear_medium_initial_conditions

  ! ---------------------------------------------------------
  subroutine system_linear_medium_do_td(this, operation)
    class(system_linear_medium_t),    intent(inout) :: this
    class(algorithmic_operation_t), intent(in)    :: operation

    PUSH_SUB(system_linear_medium_do_td)

    select case (operation%id)
    case (SKIP)
      ! Do nothing
    case (STORE_CURRENT_STATUS)
      ! For the moment we do nothing
    case (EXPMID_START)
      ! Empty for the moment
    case (EXPMID_FINISH)
      ! Empty for the moment
    case (EXPMID_PREDICT_DT_2)
      ! Empty for the moment
    case (UPDATE_HAMILTONIAN)
      ! Empty for the moment
    case (EXPMID_PREDICT_DT)
      ! Empty for the moment
    case default
      message(1) = "Unsupported TD operation."
      call messages_fatal(1, namespace=this%namespace)
    end select


    POP_SUB(system_linear_medium_do_td)
  end subroutine system_linear_medium_do_td

  ! ---------------------------------------------------------
  logical function system_linear_medium_is_tolerance_reached(this, tol) result(converged)
    class(system_linear_medium_t),   intent(in)    :: this
    FLOAT,                     intent(in)    :: tol

    PUSH_SUB(system_linear_medium_is_tolerance_reached)

    ! this routine is never called at present, no reason to be here
    ASSERT(.false.)
    converged = .false.

    POP_SUB(system_linear_medium_is_tolerance_reached)
  end function system_linear_medium_is_tolerance_reached

  ! ---------------------------------------------------------
  subroutine system_linear_medium_iteration_info(this)
    class(system_linear_medium_t), intent(in) :: this

    integer :: idir
    character(len=20) :: fmt

    PUSH_SUB(system_linear_medium_iteration_info)

    write(message(1),'(2X,A,1X,A)') "Linear medium system:", trim(this%namespace%get())

    write(message(2),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(3),'(4x,A,e14.6)') 'Simulation time: ', this%clock%time()
    call messages_info(3)

    POP_SUB(system_linear_medium_iteration_info)
  end subroutine system_linear_medium_iteration_info

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_start(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_output_start)

    ! Output info for first iteration
    call this%output_write()

    POP_SUB(system_linear_medium_output_start)
  end subroutine system_linear_medium_output_start

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_finish(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_output_finish)

    POP_SUB(system_linear_medium_output_finish)
  end subroutine system_linear_medium_output_finish

  ! ---------------------------------------------------------
  subroutine system_linear_medium_output_write(this)
    class(system_linear_medium_t), intent(inout) :: this

    integer :: idir, iat, iout
    character(len=50) :: aux
    character(1) :: out_label(2)
    FLOAT :: tmp(MAX_DIM)

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    PUSH_SUB(system_linear_medium_output_write)

    POP_SUB(system_linear_medium_output_write)
  end subroutine system_linear_medium_output_write

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_quantity(this, iq)
    class(system_linear_medium_t), intent(inout) :: this
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_linear_medium_update_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. this%quantities(iq)%protected)

    select case (iq)
!    case (PERMITTIVITY, PERMEABILITY, E_CONDUCTIVITY, M_CONDUCTIVITY)
!      this%quantities(iq)%clock = this%quantities(iq)%clock + CLOCK_TICK
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_update_quantity)
  end subroutine system_linear_medium_update_quantity

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_exposed_quantity(partner, iq)
    class(system_linear_medium_t), intent(inout) :: partner
    integer,                     intent(in)    :: iq

    PUSH_SUB(system_linear_medium_update_exposed_quantity)

    ! We are not allowed to update protected quantities!
    ASSERT(.not. partner%quantities(iq)%protected)

    select case (iq)
    case (PERMITTIVITY, PERMEABILITY, E_CONDUCTIVITY, M_CONDUCTIVITY)
      partner%quantities(iq)%clock = partner%quantities(iq)%clock + CLOCK_TICK
    case default
      message(1) = "Incompatible quantity."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_update_exposed_quantity)
  end subroutine system_linear_medium_update_exposed_quantity

  ! ---------------------------------------------------------
  subroutine system_linear_medium_copy_quantities_to_interaction(partner, interaction)
    class(system_linear_medium_t),          intent(inout) :: partner
    class(interaction_t),                 intent(inout) :: interaction

    PUSH_SUB(system_linear_medium_copy_quantities_to_interaction)

    select type (interaction)
    type is (linear_medium_em_field_t)
      ! No need to copy anything
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_copy_quantities_to_interaction)
  end subroutine system_linear_medium_copy_quantities_to_interaction


  ! ---------------------------------------------------------
  subroutine system_linear_medium_finalize(this)
    type(system_linear_medium_t), intent(inout) :: this

    type(profile_t), save :: prof

    PUSH_SUB(system_linear_medium_finalize)

    call profiling_in(prof, 'MEDIUM_BOX_END')
    call profiling_out(prof)

    call system_end(this)

    POP_SUB(system_linear_medium_finalize)
  end subroutine system_linear_medium_finalize


  ! Specific routines for this system:

  ! ---------------------------------------------------------
  subroutine generate_medium_box(this, medium_box, gr)
    type(system_linear_medium_t),intent(in)    :: this
    type(single_medium_box_t),   intent(inout) :: medium_box
    type(grid_t),                intent(in)    :: gr

    integer :: il, ip, ip_in, ip_bd, ipp, idim
    integer, allocatable :: tmp_points_map(:), tmp_bdry_map(:)
    FLOAT   :: bounds(2,gr%sb%dim), xx(gr%sb%dim), xxp(gr%sb%dim), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    logical :: inside
    type(profile_t), save :: prof

    PUSH_SUB(generate_medium_box)

    call profiling_in(prof, 'GENERATE_MEDIUM_BOX')

    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%sb%dim))
    SAFE_ALLOCATE(tmp_points_map(gr%mesh%np))
    SAFE_ALLOCATE(tmp_bdry_map(gr%mesh%np))
    tmp_points_map = 0
    tmp_bdry_map = 0

    if (this%box_shape == MEDIUM_BOX_FILE) then

      call get_points_map_from_file(this%filename, gr%mesh, medium_box%points_number,&
           medium_box%global_points_number, tmp_points_map)

    else

      do idim = 1, 3
        bounds(1,idim) = this%center(idim) - this%lsize(idim)/M_TWO
        bounds(2,idim) = this%center(idim) + this%lsize(idim)/M_TWO
      end do
      ip_in = 0
      ip_bd = 0
      do ip = 1, gr%mesh%np
        xx(1:3) = gr%mesh%x(ip, 1:3)
        inside = check_point_in_bounds(xx, bounds(:,:))
        if (check_point_in_bounds(xx, bounds(:,:))) then
          ip_in = ip_in + 1
          tmp_points_map(ip_in) = ip
        end if
        if (check_point_on_bounds(xx, bounds(:,:))) then
          ip_bd = ip_bd + 1
          tmp_bdry_map(ip_bd) = ip
        end if
      end do

      medium_box%points_number = ip_in
      medium_box%bdry_number = ip_bd

      SAFE_ALLOCATE(medium_box%bdry_map(medium_box%bdry_number))
      medium_box%bdry_map = 0
      medium_box%bdry_map = tmp_bdry_map(1:medium_box%bdry_number)

    end if
    call single_medium_box_allocate(medium_box, medium_box%points_number)
    medium_box%points_map(:) = tmp_points_map(1:medium_box%points_number)

    ! TODO: add some check that medium boxes do not overlap (now they are different systems), like this:
    !do ip_in = 1, this%points_number - 1
    !  if (any(this%points_map(ip_in+1:) == this%points_map(ip_in)) .or. &
    !      any(this%points_map(:, il+1:) == this%points_map(ip_in, il))) then
    !    message(1) = 'Linear medium boxes overlap.'
    !    call messages_fatal(1, namespace=namespace)
    !  end if
    !end do

    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__SMOOTH) then
        ASSERT(allocated(medium_box%bdry_map))
 
        xx(1:3) = gr%mesh%x(ip,1:3)
        dd_min = M_HUGE

        do ip_bd = 1, medium_box%bdry_number
          ipp = medium_box%bdry_map(ip_bd)
          xxp(1:3) = gr%mesh%x(ipp,1:3)
          dd = sqrt((xx(1) - xxp(1))**2 + (xx(2) - xxp(2))**2 + (xx(3) - xxp(3))**2)
          if (dd < dd_min) dd_min = dd
        end do

        medium_box%ep(ip_in) = P_ep + ((P_ep * this%ep_factor - P_ep)  &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
        medium_box%mu(ip_in) = P_mu + ((P_mu * this%mu_factor - P_mu) &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
        medium_box%c(ip_in) = M_ONE/sqrt(medium_box%ep(ip_in)*medium_box%mu(ip_in))
        medium_box%sigma_e(ip_in) = this%sigma_e_factor &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) )
        medium_box%sigma_m(ip_in) = this%sigma_m_factor &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) )

      else if (this%edge_profile == OPTION__LINEARMEDIUMEDGEPROFILE__EDGED) then

        medium_box%ep(ip_in) = P_ep * this%ep_factor
        medium_box%mu(ip_in) = P_mu * this%mu_factor
        medium_box%c(ip_in) = M_ONE/sqrt(medium_box%ep(ip_in)*medium_box%mu(ip_in))
        medium_box%sigma_e(ip_in) = this%sigma_e_factor
        medium_box%sigma_m(ip_in) = this%sigma_m_factor

      end if
    end do

    tmp(:) = P_ep
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      tmp(ip)= medium_box%ep(ip_in)
    end do
    call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      medium_box%aux_ep(ip_in, :) = tmp_grad(ip, :)/(M_FOUR * medium_box%ep(ip_in))
    end do

    tmp(:) = P_mu
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      tmp(ip) = medium_box%mu(ip_in)
    end do
    call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
    do ip_in = 1, medium_box%points_number
      ip = medium_box%points_map(ip_in)
      medium_box%aux_mu(ip_in, :) = tmp_grad(ip, :)/(M_FOUR * medium_box%mu(ip_in))
    end do

    !TODO: add print information about the medium box

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

    call profiling_out(prof)

    POP_SUB(generate_medium_box)
  contains

    logical pure function check_point_in_bounds(xx, bounds) result (check)
      FLOAT, intent(in) :: xx(:)
      FLOAT, intent(in) :: bounds(:,:)

      check = .false.
      if ((xx(1) >= bounds(1,1)) .and. (xx(1) <= bounds(2,1)) .and. &
          (xx(2) >= bounds(1,2)) .and. (xx(2) <= bounds(2,2)) .and. &
          (xx(3) >= bounds(1,3)) .and. (xx(3) <= bounds(2,3)) ) then
        check = .true.
      end if

    end function check_point_in_bounds

    logical pure function check_point_on_bounds(xx, bounds) result (check)
      FLOAT, intent(in) :: xx(:)
      FLOAT, intent(in) :: bounds(:,:)

      check = .false.
      if (xx(1) == bounds(1,1) .and. (xx(2) >= bounds(1,2) .and. xx(3) >= bounds(1,3)) &
                               .and. (xx(2) <= bounds(2,2) .and. xx(3) <= bounds(2,3)) .or. &
          xx(2) == bounds(1,2) .and. (xx(1) >= bounds(1,1) .and. xx(3) >= bounds(1,3)) &
                               .and. (xx(1) <= bounds(2,1) .and. xx(3) <= bounds(2,3)) .or. &
          xx(3) == bounds(1,3) .and. (xx(1) >= bounds(1,1) .and. xx(2) >= bounds(1,2)) &
                               .and. (xx(1) <= bounds(2,1) .and. xx(2) <= bounds(2,2)) .or. &
          xx(1) == bounds(2,1) .and. (xx(2) >= bounds(1,2) .and. xx(3) >= bounds(1,3)) &
                               .and. (xx(2) <= bounds(2,2) .and. xx(3) <= bounds(2,3)) .or. &
          xx(2) == bounds(2,2) .and. (xx(1) >= bounds(1,1) .and. xx(3) >= bounds(1,3)) &
                               .and. (xx(1) <= bounds(2,1) .and. xx(3) <= bounds(2,3)) .or. &
          xx(3) == bounds(2,3) .and. (xx(1) >= bounds(1,1) .and. xx(2) >= bounds(1,2)) &
                               .and. (xx(1) <= bounds(2,1) .and. xx(2) <= bounds(2,2)) ) then
        check = .true.
      end if

    end function check_point_on_bounds

  end subroutine generate_medium_box

  ! ----------------------------------------------------------
  !> Populate list of point indices for points inside the polyhedron
  subroutine get_points_map_from_file(filename, mesh, n_points, global_points_number, tmp_map, scale_factor)
    character(len=256),       intent(in)    :: filename
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(out)   :: n_points
    integer,                  intent(out)   :: global_points_number
    integer,                  intent(inout) :: tmp_map(:)
    FLOAT, optional,          intent(in)    :: scale_factor

    integer :: ip_in, ip
    FLOAT   :: xx(3)
    type(cgal_polyhedra_t) :: cgal_poly
    type(profile_t), save :: prof

    PUSH_SUB(get_points_map_from_file)

    call profiling_in(prof, 'GET_POINTS_MAP_FROM_FILE')

    call cgal_polyhedron_init(cgal_poly, trim(filename), verbose = .false.)

    ip_in = 0
    do ip = 1, mesh%np
      if (present(scale_factor)) then
        xx(1:3) = scale_factor * mesh%x(ip, 1:3)
      else
        xx(1:3) = mesh%x(ip, 1:3)
      end if
      if (cgal_polyhedron_point_inside(cgal_poly, xx(1), xx(2), xx(3))) then
        ip_in = ip_in + 1
        tmp_map(ip_in) = ip
      end if
    end do
    n_points = ip_in
    call cgal_polyhedron_end(cgal_poly)

#ifdef HAVE_MPI
    call MPI_Allreduce(ip_in, global_points_number, 1, &
        MPI_INT, MPI_SUM, MPI_COMM_WORLD, mpi_err)
#else
    global_points_number = n_points
#endif

    call profiling_out(prof)

    POP_SUB(get_points_map_from_file)

  end subroutine get_points_map_from_file

end module system_linear_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
