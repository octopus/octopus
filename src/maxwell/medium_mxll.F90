!! Copyright (C) 2019 R. Jestaedt, F. Bonafe, H. Appel
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
!! $Id: propagator.F90 13908 2015-05-05 06:02:30Z xavier $

#include "global.h"

module medium_mxll_oct_m
  use cgal_polyhedra_oct_m
  use derivatives_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m
  use comm_oct_m
  use grid_oct_m
  use mesh_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::           &
    medium_box_t,     &
    medium_box_init,  &
    medium_box_end

   type medium_box_t
     integer                         :: number   !< number of linear media boxes
     integer, allocatable            :: shape(:)  !< edge shape profile (smooth or steep)
     FLOAT, allocatable              :: center(:,:) !< center of each box
     FLOAT, allocatable              :: lsize(:,:)  !< length in each direction of each box
     FLOAT, allocatable              :: ep(:,:) !< permitivity of the linear media
     FLOAT, allocatable              :: mu(:,:) !< permeability of the linear media
     FLOAT, allocatable              :: c(:,:) !< speed of light in the linear media
     FLOAT, allocatable              :: ep_factor(:) !< permitivity before applying edge profile
     FLOAT, allocatable              :: mu_factor(:) !< permeability before applying edge profile
     FLOAT, allocatable              :: sigma_e_factor(:) !< electric conductivy before applying edge profile
     FLOAT, allocatable              :: sigma_m_factor(:) !< magnetic conductivity before applying edge4 profile
     FLOAT, allocatable              :: sigma_e(:,:) !< electric conductivy of (lossy) medium
     FLOAT, allocatable              :: sigma_m(:,:) !< magnetic conductivy of (lossy) medium
     integer, allocatable            :: points_number(:)
     integer, allocatable            :: global_points_number(:)
     integer, allocatable            :: points_map(:,:)
     FLOAT, allocatable              :: aux_ep(:,:,:) !< auxiliary array for storing the epsilon derivative profile
     FLOAT, allocatable              :: aux_mu(:,:,:) !< auxiliary array for storing the softened mu profile
     integer, allocatable            :: bdry_number(:)
     integer, allocatable            :: bdry_map(:,:)
     character(len=256), allocatable :: filename(:)
     FLOAT                           :: width !< width of medium medium when used as boundary condition
   end type medium_box_t

 contains

  ! ---------------------------------------------------------
  subroutine medium_box_init(namespace, medium_box, calc_medium_box, gr)
    type(namespace_t),      intent(in)    :: namespace
    type(grid_t),           intent(in)    :: gr
    type(medium_box_t),     intent(inout) :: medium_box
    logical,                intent(out)   :: calc_medium_box

    integer :: nlines, ncols, idim, il, ip_in_max2
    integer, allocatable :: tmp(:,:)
    type(block_t) :: blk
    type(medium_box_t), allocatable :: medium_box_aux
    logical :: checkmediumpoints
    type(profile_t), save :: prof

    PUSH_SUB(medium_box_init)

    call profiling_in(prof, 'MEDIUM_BOX_INIT')

    !%Variable LinearMediumBox
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Defines parameters for a linear medium box.
    !%
    !% Example:
    !%
    !% <tt>%LinearMediumBox
    !% <br>&nbsp;&nbsp;   center_x | center_y | center_z | x_length | y_length | z_length | epsilon_factor | mu_factor | sigma_e | sigma_m | edged/smooth
    !% <br>%</tt>
    !%
    !% Position of center (three components) and length (three components), followed by permittivity
    !% factor, electric conductivity and magnetic conductivity, and finally type of numerical
    !% approximation used for the derivatives at the edges.
    !%
    !%Option edged 1
    !% Medium box edges are considered steep for derivatives.
    !%Option smooth 2
    !% Medium box edged and softened for derivatives.
    !%End
    if(parse_block(namespace, 'LinearMediumBox', blk) == 0) then

      call messages_print_stress(stdout, trim('Maxwell Medium box:'))
      calc_medium_box = .true.

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)
      SAFE_ALLOCATE(medium_box%center(1:3,1:nlines))
      SAFE_ALLOCATE(medium_box%lsize(1:3,1:nlines))
      SAFE_ALLOCATE(medium_box%ep_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%mu_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%sigma_e_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%sigma_m_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%shape(1:nlines))
      do il = 1, nlines
        ncols = parse_block_cols(blk, il-1)
        if (ncols /= 11) then
          call messages_input_error(namespace, 'LinearMediumBox', 'should consist of eleven columns', row=il-1)
        end if
        do idim = 1, 3
          call parse_block_float(blk, il-1, idim-1, medium_box%center(idim,il))
          call parse_block_float(blk, il-1, idim+2, medium_box%lsize(idim,il))
        end do
        call parse_block_float(blk, il-1, 6, medium_box%ep_factor(il))
        call parse_block_float(blk, il-1, 7, medium_box%mu_factor(il))
        call parse_block_float(blk, il-1, 8, medium_box%sigma_e_factor(il))
        call parse_block_float(blk, il-1, 9, medium_box%sigma_m_factor(il))
        call parse_block_integer(blk, il-1, 10, medium_box%shape(il))
        write(message(1),'(a,I1)')    'Medium box number:  ', il
        write(message(2),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', medium_box%center(1,il), ' | ',&
            medium_box%center(2,il), ' | ', medium_box%center(3,il)
        write(message(3),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', medium_box%lsize(1,il), ' | ', &
            medium_box%lsize(2,il), ' | ', medium_box%lsize(3,il)
        write(message(4),'(a,es9.2)') 'Box epsilon factor: ', medium_box%ep_factor(il)
        write(message(5),'(a,es9.2)') 'Box mu factor:      ', medium_box%mu_factor(il)
        write(message(6),'(a,es9.2)') 'Box electric sigma: ', medium_box%sigma_e_factor(il)
        write(message(7),'(a,es9.2)') 'Box magnetic sigma: ', medium_box%sigma_m_factor(il)
        if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__EDGED) then
          write(message(8),'(a,a)')   'Box shape:          ', 'edged'
        else if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__SMOOTH) then
          write(message(8),'(a,a)')   'Box shape:          ', 'smooth'
        end if
        write(message(9),'(a)') ""
        call messages_info(9)
      end do
      call parse_block_end(blk)

      call generate_medium_boxes(medium_box, gr, nlines, namespace)

      call messages_print_stress(stdout)
    end if

    !%Variable LinearMediumFromFile
    !%Type block
    !%Section Time-Dependent::Propagation
    !%Description
    !% Defines parameters and geometry to create a new linear medium box.
    !%
    !% Example:
    !%
    !% <tt>%LinearMediumFromFile
    !% <br>&nbsp;&nbsp;   medium_surface_file | epsilon_factor | mu_factor | sigma_e | sigma_m | edged/smooth
    !% <br>%</tt>
    !%
    !% Medium surface file, followed by permittivity
    !% factor, electric conductivity and magnetic conductivity, and finally type of numerical
    !% approximation used for the derivatives at the edges.
    !%
    !%Option edged 1
    !% Medium box edges are considered steep for derivatives.
    !%Option smooth 2
    !% Medium box edged and softened for derivatives.
    !%End
    if(parse_block(namespace, 'LinearMediumFromFile', blk) == 0) then

      call messages_print_stress(stdout, trim('Maxwell Medium box:'))
      calc_medium_box = .true.

      nlines = parse_block_n(blk)
      SAFE_ALLOCATE(medium_box%filename(1:nlines))
      SAFE_ALLOCATE(medium_box%ep_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%mu_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%sigma_e_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%sigma_m_factor(1:nlines))
      SAFE_ALLOCATE(medium_box%shape(1:nlines))
      do il = 1, nlines
        ncols = parse_block_cols(blk, il-1)
        call parse_block_string(blk, il-1, 0, medium_box%filename(il))
        call parse_block_float(blk, il-1, 1, medium_box%ep_factor(il))
        call parse_block_float(blk, il-1, 2, medium_box%mu_factor(il))
        call parse_block_float(blk, il-1, 3, medium_box%sigma_e_factor(il))
        call parse_block_float(blk, il-1, 4, medium_box%sigma_m_factor(il))
        call parse_block_integer(blk, il-1, 5, medium_box%shape(il))
        if (medium_box%shape(il) /= OPTION__LINEARMEDIUMBOX__EDGED) then
          call messages_not_implemented("Medium box from file only implemented with edged boundaries.", namespace=namespace)
        end if
        write(message(1),'(a,I1)')    'Medium box number:  ', il
        write(message(2),'(a,a)') 'Box surface file: ', trim(medium_box%filename(il))
        write(message(3),'(a,es9.2)') 'Box epsilon factor: ', medium_box%ep_factor(il)
        write(message(4),'(a,es9.2)') 'Box mu factor:      ', medium_box%mu_factor(il)
        write(message(5),'(a,es9.2)') 'Box electric sigma: ', medium_box%sigma_e_factor(il)
        write(message(6),'(a,es9.2)') 'Box magnetic sigma: ', medium_box%sigma_m_factor(il)
        write(message(7),'(a,a)')   'Box shape:          ', 'edged'
        write(message(8),'(a)') ""
        call messages_info(8)
      end do
      call parse_block_end(blk)

      call generate_medium_boxes(medium_box, gr, nlines, namespace)

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
        SAFE_ALLOCATE(tmp(gr%mesh%np, nlines))
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

      call messages_print_stress(stdout)
    end if

    call profiling_out(prof)

    POP_SUB(medium_box_init)

  end subroutine medium_box_init

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
        if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__SMOOTH) then
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

        else if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__EDGED) then

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

  ! ---------------------------------------------------------
  !> Deallocation of medium_box components
  subroutine medium_box_end(medium_box)
    type(medium_box_t),   intent(inout)    :: medium_box

    type(profile_t), save :: prof

    PUSH_SUB(medium_box_end)

    call profiling_in(prof, 'MEDIUM_BOX_END')

    SAFE_DEALLOCATE_A(medium_box%center)
    SAFE_DEALLOCATE_A(medium_box%lsize)
    SAFE_DEALLOCATE_A(medium_box%ep_factor)
    SAFE_DEALLOCATE_A(medium_box%mu_factor)
    SAFE_DEALLOCATE_A(medium_box%sigma_e_factor)
    SAFE_DEALLOCATE_A(medium_box%sigma_m_factor)
    SAFE_DEALLOCATE_A(medium_box%shape)
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

    POP_SUB(medium_box_end)

  end subroutine medium_box_end

end module medium_mxll_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
