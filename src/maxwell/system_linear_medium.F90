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
  use geometry_oct_m
  use global_oct_m
  use interaction_oct_m
  use interactions_factory_oct_m
  use io_oct_m
  use ion_dynamics_oct_m
  use iso_c_binding
  use lasers_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use propagator_oct_m
  use propagator_verlet_oct_m
  use quantity_oct_m
  use space_oct_m
  use species_oct_m
  use system_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use write_iter_oct_m

  implicit none

  private
  public ::           &
    system_linear_medium_t,    &
    system_linear_medium_init

  integer, parameter :: &
    MEDIUM_PARALLELEPIPED = 1,         &
    MEDIUM_BOX_FILE       = 2

  type, extends(system_t) :: system_linear_medium_t
     integer                         :: box_shape
     integer                         :: number   !< number of linear media boxes
     integer, allocatable            :: edge_profile(:)  !< edge shape profile (smooth or steep)
     FLOAT              :: center(3) !< center of a box
     FLOAT              :: lsize(3)  !< length in each direction of a box
     FLOAT, allocatable              :: ep(:) !< permitivity of the linear media
     FLOAT, allocatable              :: mu(:) !< permeability of the linear media
     FLOAT, allocatable              :: c(:) !< speed of light in the linear media
     FLOAT, allocatable              :: ep_factor !< permitivity before applying edge profile
     FLOAT, allocatable              :: mu_factor !< permeability before applying edge profile
     FLOAT, allocatable              :: sigma_e_factor !< electric conductivy before applying edge profile
     FLOAT, allocatable              :: sigma_m_factor !< magnetic conductivity before applying edge4 profile
     FLOAT, allocatable              :: sigma_e(:) !< electric conductivy of (lossy) medium
     FLOAT, allocatable              :: sigma_m(:) !< magnetic conductivy of (lossy) medium
     integer, allocatable            :: points_number(:)
     integer, allocatable            :: global_points_number(:)
     integer, allocatable            :: points_map(:,:)
     FLOAT, allocatable              :: aux_ep(:,:,:) !< auxiliary array for storing the epsilon derivative profile
     FLOAT, allocatable              :: aux_mu(:,:,:) !< auxiliary array for storing the softened mu profile
     integer, allocatable            :: bdry_number(:)
     integer, allocatable            :: bdry_map(:,:)
     character(len=256)  :: filename
     FLOAT                           :: width !< width of medium medium when used as boundary condition

  contains
    procedure :: init_interaction => system_linear_medium_init_interaction
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
    procedure :: update_interactions_start => system_linear_medium_update_interactions_start
    procedure :: update_interactions_finish => system_linear_medium_update_interactions_finish
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


    PUSH_SUB(system_linear_medium_init)

    this%namespace = namespace


    call profiling_in(prof, 'MEDIUM_BOX_INIT')

    !%Variable LinearMediumBoxShape
    !%Type integer
    !%Section
    !%Description
    !% This variable defines the shape of the linear medium box.
    !% The default is <tt>parallelepiped</tt>.
    !%Option parallelepiped 1
    !% The medium box will be a parallelepiped whose center and dimensions are taken from
    !% the variable <tt>LinearMediumLsize</tt>.
    !%Option box_file 2
    !% The simulation box will be read from an external file in OFF format, defined by the variable <tt>LinearMediumBoxFile</tt>.
    !%End
    call parse_variable(namespace, 'LinearMediumBoxShape', MEDIUM_PARALLELEPIPED, this%box_shape)

    if(.not.varinfo_valid_option('LinearMediumBoxShape', this%box_shape)) call messages_input_error(namespace, 'LinearMediumBoxShape')

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
        calc_medium_box = .true.

        nlines = parse_block_n(blk)
        if (nlines /=  1) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of one line')
        end if
        SAFE_ALLOCATE(medium_box%center(1:3))
        SAFE_ALLOCATE(medium_box%lsize(1:3))
        ncols = parse_block_cols(blk, 0)
        if (ncols /= 6) then
          call messages_input_error(namespace, 'LinearMediumBoxSize', 'should consist of six columns')
        end if
        do idim = 1, 3
          call parse_block_float(blk, 0, idim-1, medium_box%center(idim))
          call parse_block_float(blk, 0, idim+2, medium_box%lsize(idim))
        end do
        write(message(1),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', medium_box%center(1), ' | ',&
            medium_box%center(2), ' | ', medium_box%center(3)
        write(message(2),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', medium_box%lsize(1), ' | ', &
            medium_box%lsize(2), ' | ', medium_box%lsize(3)
        write(message(3),'(a)') ""
        call messages_info(3)
      call parse_block_end(blk)

      call generate_medium_boxes(medium_box, gr, nlines, namespace)

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
    if parse_is_defined(namespace, 'LinearMediumBoxFile') then
      call parse_variable(namespace, 'LinearMediumBoxFile', this%filename)
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
    call parse_variable(namespace, 'CheckPointsMediumFromFile', .false., checkmediumpoints)

    if (checkmediumpoints .and. (nlines > 1)) then
      message(1) = 'Check for points only works for one medium surface, centered at the origin.'
      call messages_fatal(1, namespace=namespace)
    end if

    if (checkmediumpoints) then
      allocate(medium_box_aux, source=medium_box)
      SAFE_ALLOCATE(tmp(gr%mesh%np, 1))
      ip_in_max2 = 0
      write(message(1),'(a, a, a)')   'Check of points inside surface of medium ', trim(medium_box%filename(1)), ":"
      call messages_info(1)
      call get_points_map_from_file(medium_box_aux, ip_in_max2, gr%mesh, tmp, CNST(0.99))

      write(message(1),'(a, I8)')'Number of points inside medium (normal coordinates):', medium_box%global_points_number(1)
      write(message(2),'(a, I8)')'Number of points inside medium (rescaled coordinates):', medium_box_aux%global_points_number(1)
      write(message(3), '(a)') ""
      call messages_info(3)

      deallocate(medium_box_aux)
      SAFE_DEALLOCATE_A(tmp)
    end if

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
    !.
    !%End

    if (parse_block(namespace, 'LinearMediumProperties', blk) == 0) then

      calc_medium_box = .true.

      nlines = parse_block_n(blk)
      if (nlines /=  1) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of one line', nlines)
      end if
      ncols = parse_block_cols(blk, 0)
      if (ncols /= 4) then
        call messages_input_error(namespace, 'LinearMediumProperties', 'should consist of four columns')
      end if
      call parse_block_float(blk, 0, 0, medium_box%ep_factor)
      call parse_block_float(blk, 0, 1, medium_box%mu_factor)
      call parse_block_float(blk, 0, 2, medium_box%sigma_e_factor)
      call parse_block_float(blk, 0, 3, medium_box%sigma_m_factor)
      write(message(1),'(a,es9.2)') 'Box epsilon factor: ', medium_box%ep_factor
      write(message(2),'(a,es9.2)') 'Box mu factor:      ', medium_box%mu_factor
      write(message(3),'(a,es9.2)') 'Box electric sigma: ', medium_box%sigma_e_factor
      write(message(4),'(a,es9.2)') 'Box magnetic sigma: ', medium_box%sigma_m_factor
      write(message(5),'(a)') ""
      call messages_info(5)
      call parse_block_end(blk)

      call generate_medium_boxes(medium_box, gr, nlines, namespace)

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
    call parse_variable(namespace, 'LinearMediumEdgeProfile', OPTION__LINEARMEDIUMBOX__EDGED, medium_box%edge_profile)
    if (medium_box%edge_profile == OPTION__LINEARMEDIUMBOX__EDGED) then
      write(message(1),'(a,a)')   'Box shape:          ', 'edged'
    else if (medium_box%edge_profile == OPTION__LINEARMEDIUMBOX__SMOOTH) then
      write(message(1),'(a,a)')   'Box shape:          ', 'smooth'
    end if
    call messages_info(1)


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

    write(message(1),'(2X,A,1X,A)') "Linar medium system:", trim(this%namespace%get())

    write(message(5),'(4x,A,I8.7)')  'Clock tick:      ', this%clock%get_tick()
    write(message(6),'(4x,A,e14.6)') 'Simulation time: ', this%clock%time()
    call messages_info(6)

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

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_end(this%output_handle(1))
      call write_iter_end(this%output_handle(2))
    end if

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
    class default
      message(1) = "Unsupported interaction."
      call messages_fatal(1)
    end select

    POP_SUB(system_linear_medium_copy_quantities_to_interaction)
  end subroutine system_linear_medium_copy_quantities_to_interaction

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_interactions_start(this)
    class(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_update_interactions_start)

    POP_SUB(system_linear_medium_update_interactions_start)
  end subroutine system_linear_medium_update_interactions_start

  ! ---------------------------------------------------------
  subroutine system_linear_medium_update_interactions_finish(this)
    class(system_linear_medium_t), intent(inout) :: this

    type(interaction_iterator_t) :: iter

    PUSH_SUB(system_linear_medium_update_interactions_finish)

    call iter%start(this%interactions)
    do while (iter%has_next())
      select type (interaction => iter%get_next())
      class default
        message(1) = "Interactions not implemented for linear medium systems."
        call messages_fatal(1)
      end select
    end do

    POP_SUB(system_linear_medium_update_interactions_finish)
  end subroutine system_linear_medium_update_interactions_finish

  ! ---------------------------------------------------------
  subroutine system_linear_medium_finalize(this)
    type(system_linear_medium_t), intent(inout) :: this

    PUSH_SUB(system_linear_medium_finalize)

    call profiling_in(prof, 'MEDIUM_BOX_END')

    SAFE_DEALLOCATE_A(medium_box%center)
    SAFE_DEALLOCATE_A(medium_box%lsize)
    SAFE_DEALLOCATE_A(medium_box%ep_factor)
    SAFE_DEALLOCATE_A(medium_box%mu_factor)
    SAFE_DEALLOCATE_A(medium_box%sigma_e_factor)
    SAFE_DEALLOCATE_A(medium_box%sigma_m_factor)
    SAFE_DEALLOCATE_A(medium_box%edge_profile)
    SAFE_DEALLOCATE_A(medium_box%points_number)
    SAFE_DEALLOCATE_A(medium_box%bdry_number)
    SAFE_DEALLOCATE_A(medium_box%points_map)
    SAFE_DEALLOCATE_A(medium_box%bdry_map)
    SAFE_DEALLOCATE_A(medium_box%aux_ep)
    SAFE_DEALLOCATE_A(medium_box%aux_mu)
    SAFE_DEALLOCATE_A(medium_box%c)
    SAFE_DEALLOCATE_A(medium_box%ep)
    SAFE_DEALLOCATE_A(medium_box%mu)
    SAFE_DEALLOCATE_A(medium_box%sigma_e)
    SAFE_DEALLOCATE_A(medium_box%sigma_m)

    call profiling_out(prof)

    call system_end(this)

    POP_SUB(system_linear_medium_finalize)
  end subroutine system_linear_medium_finalize


  ! Specific routines for this type:

  ! ---------------------------------------------------------
  subroutine generate_medium_boxes(medium_box, gr, nr_of_boxes, namespace)
    type(medium_box_t),  intent(inout)      :: medium_box
    type(grid_t),        intent(in)         :: gr
    integer,             intent(in)         :: nr_of_boxes
    type(namespace_t),   intent(in)         :: namespace

    integer :: il, ip, ip_in, ip_in_max, ip_bd, ip_bd_max, ipp, idim
    integer, allocatable :: tmp_points_map(:,:), tmp_bdry_map(:,:)
    FLOAT   :: bounds(nr_of_boxes,2,gr%sb%dim), xx(gr%sb%dim), xxp(gr%sb%dim), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)
    logical :: inside
    type(profile_t), save :: prof

    PUSH_SUB(generate_medium_boxes)

    call profiling_in(prof, 'GENERATE_MEDIUM_BOXES')

    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%sb%dim))
    SAFE_ALLOCATE(tmp_points_map(gr%mesh%np, nr_of_boxes))
    SAFE_ALLOCATE(tmp_bdry_map(gr%mesh%np, nr_of_boxes))
    tmp_points_map = 0
    tmp_bdry_map = 0

    SAFE_ALLOCATE(medium_box%points_number(nr_of_boxes))
    SAFE_ALLOCATE(medium_box%global_points_number(nr_of_boxes))
    SAFE_ALLOCATE(medium_box%bdry_number(nr_of_boxes))
    medium_box%number = nr_of_boxes

    ip_in_max = 0
    ip_bd_max = 0

    if (allocated(medium_box%filename)) then

       call get_points_map_from_file(medium_box, ip_in_max, gr%mesh, tmp_points_map)

       SAFE_ALLOCATE(medium_box%points_map(ip_in_max, nr_of_boxes))
       SAFE_ALLOCATE(medium_box%bdry_map(1, nr_of_boxes))

       medium_box%points_map = 0
       medium_box%bdry_map = 0

       medium_box%points_map(:,:) = tmp_points_map(1:ip_in_max,:)

    else

       do il = 1, nr_of_boxes
          do idim = 1, 3
             bounds(il,1,idim) = medium_box%center(idim,il) - medium_box%lsize(idim,il)/M_TWO
             bounds(il,2,idim) = medium_box%center(idim,il) + medium_box%lsize(idim,il)/M_TWO
          end do
          ip_in = 0
          ip_bd = 0
          do ip = 1, gr%mesh%np
             xx(1:3) = gr%mesh%x(ip, 1:3)
             inside = check_point_in_bounds(xx, bounds(il,:,:))
             if (check_point_in_bounds(xx, bounds(il,:,:))) then
                ip_in = ip_in + 1
                tmp_points_map(ip_in, il) = ip
             end if
             if (check_point_on_bounds(xx, bounds(il,:,:))) then
                ip_bd = ip_bd + 1
                tmp_bdry_map(ip_bd, il) = ip
             end if
          end do
          if (ip_in > ip_in_max) ip_in_max = ip_in
          if (ip_bd > ip_bd_max) ip_bd_max = ip_bd
          medium_box%points_number(il) = ip_in
          medium_box%bdry_number(il) = ip_bd
       end do

    SAFE_ALLOCATE(medium_box%points_map(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%bdry_map(ip_bd_max,nr_of_boxes))

    medium_box%points_map = 0
    medium_box%bdry_map = 0
    medium_box%points_map = tmp_points_map(1:ip_in_max,1:nr_of_boxes)
    medium_box%bdry_map = tmp_bdry_map(1:ip_bd_max,1:nr_of_boxes)

    end if

    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    do il = 1, nr_of_boxes
      do ip_in = 1, medium_box%points_number(il) - 1
        if (any(medium_box%points_map(ip_in+1:, il) == medium_box%points_map(ip_in, il)) .or. &
            any(medium_box%points_map(:, il+1:) == medium_box%points_map(ip_in, il))) then
          message(1) = 'Linear medium boxes overlap.'
          call messages_fatal(1, namespace=namespace)
        end if
      end do
    end do

    SAFE_ALLOCATE(medium_box%aux_ep(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%aux_mu(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%c(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%ep(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%mu(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%sigma_e(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%sigma_m(ip_in_max,nr_of_boxes))

    do il = 1, nr_of_boxes

      do ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in,il)
        if (medium_box%edge_profile(il) == OPTION__LINEARMEDIUMBOX__SMOOTH) then
          xx(1:3) = gr%mesh%x(ip,1:3)
          dd_min = M_HUGE

          do ip_bd = 1, medium_box%bdry_number(il)
            ipp = medium_box%bdry_map(ip_bd, il)
            xxp(1:3) = gr%mesh%x(ipp,1:3)
            dd = sqrt((xx(1) - xxp(1))**2 + (xx(2) - xxp(2))**2 + (xx(3) - xxp(3))**2)
            if (dd < dd_min) dd_min = dd
          end do

          medium_box%ep(ip_in,il) = P_ep + ((P_ep * medium_box%ep_factor(il) - P_ep)  &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
          medium_box%mu(ip_in,il) = P_mu + ((P_mu * medium_box%mu_factor(il) - P_mu) &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max))))
          medium_box%c(ip_in,il) = M_ONE/sqrt(medium_box%ep(ip_in, il)*medium_box%mu(ip_in, il))
          medium_box%sigma_e(ip_in,il) = medium_box%sigma_e_factor(il) &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) )
          medium_box%sigma_m(ip_in,il) = medium_box%sigma_m_factor(il) &
            * M_ONE/(M_ONE + exp(-M_FIVE/dd_max * (dd_min - M_TWO*dd_max)) )

        else if (medium_box%edge_profile(il) == OPTION__LINEARMEDIUMBOX__EDGED) then

          medium_box%ep(ip_in, il) = P_ep * medium_box%ep_factor(il)
          medium_box%mu(ip_in, il) = P_mu * medium_box%mu_factor(il)
          medium_box%c(ip_in, il) = M_ONE/sqrt(medium_box%ep(ip_in, il)*medium_box%mu(ip_in, il))
          medium_box%sigma_e(ip_in, il) = medium_box%sigma_e_factor(il)
          medium_box%sigma_m(ip_in, il) = medium_box%sigma_m_factor(il)

        end if
      end do

      tmp(:) = P_ep
      do  ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in, il)
        tmp(ip)= medium_box%ep(ip_in, il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in, il)
        medium_box%aux_ep(ip_in, :, il) = tmp_grad(ip, :)/(M_FOUR * medium_box%ep(ip_in, il))
      end do

      tmp(:) = P_mu
      do ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in, il)
        tmp(ip) = medium_box%mu(ip_in, il)
      end do
      call dderivatives_grad(gr%der, tmp, tmp_grad, set_bc = .false.)
      do ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in, il)
        medium_box%aux_mu(ip_in, :, il) = tmp_grad(ip, :)/(M_FOUR * medium_box%mu(ip_in, il))
      end do

      !TODO: add print information about the medium box

    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

    call profiling_out(prof)

    POP_SUB(generate_medium_boxes)
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

  end subroutine generate_medium_boxes

  ! ----------------------------------------------------------
  !> Populate list of point indices for points inside the polyhedron
  subroutine get_points_map_from_file(medium_box, ip_in_max, mesh, tmp_map, scale_factor)
    type(medium_box_t),       intent(inout) :: medium_box
    integer,                  intent(inout) :: ip_in_max
    type(mesh_t),             intent(in)    :: mesh
    integer,                  intent(inout) :: tmp_map(:,:)
    FLOAT, optional,          intent(in)    :: scale_factor

    integer :: il, ip_in, ip
    FLOAT   :: xx(3)
    type(cgal_polyhedra_t), allocatable :: cgal_poly(:)
    type(profile_t), save :: prof

    PUSH_SUB(get_points_map_from_file)

    call profiling_in(prof, 'GET_POINTS_MAP_FROM_FILE')

    SAFE_ALLOCATE(cgal_poly(1:medium_box%number))

    do il = 1, medium_box%number
      call cgal_polyhedron_init(cgal_poly(il), trim(medium_box%filename(il)), verbose = .false.)

      ip_in = 0
      do ip = 1, mesh%np
        if (present(scale_factor)) then
          xx(1:3) = scale_factor * mesh%x(ip, 1:3)
        else
          xx(1:3) = mesh%x(ip, 1:3)
        end if
        if (cgal_polyhedron_point_inside(cgal_poly(il), xx(1), xx(2), xx(3))) then
          ip_in = ip_in + 1
          tmp_map(ip_in, il) = ip
        end if
      end do
      if (ip_in > ip_in_max) ip_in_max = ip_in
      medium_box%points_number(il) = ip_in
      call cgal_polyhedron_end(cgal_poly(il))

#ifdef HAVE_MPI
      call MPI_Allreduce(ip_in, medium_box%global_points_number(il), 1, &
          MPI_INT, MPI_SUM, MPI_COMM_WORLD, mpi_err)
#else
      medium_box%global_points_number(il) = medium_box%points_number(il)
#endif
    end do

    SAFE_DEALLOCATE_A(cgal_poly)

    call profiling_out(prof)

    POP_SUB(get_points_map_from_file)

  end subroutine get_points_map_from_file

end module system_linear_medium_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
