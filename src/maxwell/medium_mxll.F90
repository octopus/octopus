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
#include "undef.F90"
#include "complex.F90"

module medium_mxll_oct_m
  use derivatives_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use comm_oct_m
  use grid_oct_m
  use mesh_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m

  private
  public ::             &
      medium_box_t,     &
      medium_box_init,  &
      medium_box_end,   &
      generate_medium_boxes

   type :: medium_box_t
     integer                        :: number   !< number of linear media boxes
     integer, allocatable           :: shape(:)  !< edge shape profile (smooth or steep)
     FLOAT, allocatable             :: center(:,:) !< center of each box
     FLOAT, allocatable             :: lsize(:,:)  !< length in each direction of each box
     FLOAT, allocatable             :: ep(:,:) !< permitivity of the linear media
     FLOAT, allocatable             :: mu(:,:) !< permeability of the linear media
     FLOAT, allocatable             :: c(:,:) !< speed of light in the linear media
     FLOAT, allocatable             :: ep_factor(:) !< permitivity before applying edge profile
     FLOAT, allocatable             :: mu_factor(:) !< permeability before applying edge profile
     FLOAT, allocatable             :: sigma_e_factor(:) !< electric conductivy before applying edge profile
     FLOAT, allocatable             :: sigma_m_factor(:) !< magnetic conductivity before applying edge4 profile
     FLOAT, allocatable             :: sigma_e(:,:) !< electric conductivy of (lossy) medium
     FLOAT, allocatable             :: sigma_m(:,:) !< magnetic conductivy of (lossy) medium
     integer, allocatable           :: points_number(:)
     integer, allocatable           :: points_map(:,:)
     FLOAT, allocatable             :: aux_ep(:,:,:) !< auxiliary array for storing the epsilon derivative profile
     FLOAT, allocatable             :: aux_mu(:,:,:) !< auxiliary array for storing the softened mu profile
     integer, allocatable           :: bdry_number(:)
     FLOAT, allocatable             :: bdry_map(:,:)
   end type medium_box_t

 contains

  ! ---------------------------------------------------------
  subroutine medium_box_init(namespace, medium_box, calc_medium_box, gr)
    type(namespace_t),      intent(in)    :: namespace
    type(grid_t),           intent(in)    :: gr
    type(medium_box_t),     intent(inout) :: medium_box
    logical,                intent(out)   :: calc_medium_box

    integer :: nlines, ncols, icol, idim, il
    type(block_t) :: blk
    
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
        if (il > 1) then
          write(message(1),'(a)') ""
          write(message(2),'(a,I1)')    'Medium box number:  ', il
          write(message(3),'(a,es9.2,a,es9.2,a,es9.2)') 'Box center:         ', medium_box%center(1,il), ' | ',&
                medium_box%center(2,il), ' | ', medium_box%center(3,il)
          write(message(4),'(a,es9.2,a,es9.2,a,es9.2)') 'Box size  :         ', medium_box%lsize(1,il), ' | ', &
                medium_box%lsize(2,il), ' | ', medium_box%lsize(3,il)
          write(message(5),'(a,es9.2)') 'Box epsilon factor: ', medium_box%ep_factor(il)
          write(message(6),'(a,es9.2)') 'Box mu factor:      ', medium_box%mu_factor(il)
          write(message(7),'(a,es9.2)') 'Box electric sigma: ', medium_box%sigma_e_factor(il)
          write(message(8),'(a,es9.2)') 'Box magnetic sigma: ', medium_box%sigma_m_factor(il)
          if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__EDGED) then
            write(message(9),'(a,a)')   'Box shape:          ', 'edged'
          else if (medium_box%shape(il) == OPTION__LINEARMEDIUMBOX__SMOOTH) then
            write(message(9),'(a,a)')   'Box shape:          ', 'smooth'
          end if
          call messages_info(9)
        else
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
          call messages_info(8)
        end if
      end do
      call parse_block_end(blk)

      call generate_medium_boxes(medium_box, gr, nlines, namespace)

      call messages_print_stress(stdout)
    end if

  end subroutine medium_box_init

  ! ---------------------------------------------------------
  subroutine generate_medium_boxes(medium_box, gr, nr_of_boxes, namespace)
    type(medium_box_t),  intent(inout)      :: medium_box
    type(grid_t),        intent(in)         :: gr
    integer,             intent(in)         :: nr_of_boxes
    type(namespace_t),   intent(in)         :: namespace

    integer :: il, ip, ip_in, ip_in_max, ip_bd, ip_bd_max, ipp, idim
    FLOAT   :: bounds(nr_of_boxes,2,gr%sb%dim), xx(gr%sb%dim), xxp(gr%sb%dim), dd, dd_max, dd_min
    FLOAT, allocatable  :: tmp(:), tmp_grad(:,:)

    PUSH_SUB(generate_medium_boxes)

    SAFE_ALLOCATE(tmp(gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_grad(gr%mesh%np_part,1:gr%mesh%sb%dim))

    SAFE_ALLOCATE(medium_box%points_number(nr_of_boxes))
    SAFE_ALLOCATE(medium_box%bdry_number(nr_of_boxes))
    medium_box%number = nr_of_boxes

    ip_in_max = 0
    ip_bd_max = 0
    do il = 1, nr_of_boxes
      do idim = 1, 3
        bounds(il,1,idim) = medium_box%center(idim,il) - medium_box%lsize(idim,il)/M_TWO
        bounds(il,2,idim) = medium_box%center(idim,il) + medium_box%lsize(idim,il)/M_TWO
      end do
      ip_in = 0
      ip_bd = 0
      do ip = 1, gr%mesh%np
        xx(1:3) = gr%mesh%x(ip, 1:3)
        if (check_point_in_bounds(xx, bounds(il,:,:))) then
          ip_in = ip_in + 1
        end if
        if (check_point_on_bounds(xx, bounds(il,:,:))) then
          ip_bd = ip_bd + 1
        end if
      end do
      if (ip_in > ip_in_max) ip_in_max = ip_in
      if (ip_bd > ip_bd_max) ip_bd_max = ip_bd
      medium_box%points_number(il) = ip_in
      medium_box%bdry_number(il) = ip_bd
    end do

    dd_max = max(2*gr%mesh%spacing(1), 2*gr%mesh%spacing(2), 2*gr%mesh%spacing(3))

    SAFE_ALLOCATE(medium_box%points_map(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%bdry_map(ip_bd_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%aux_ep(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%aux_mu(ip_in_max,1:3,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%c(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%ep(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%mu(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%sigma_e(ip_in_max,nr_of_boxes))
    SAFE_ALLOCATE(medium_box%sigma_m(ip_in_max,nr_of_boxes))

    medium_box%points_map = int(M_zero)
    medium_box%bdry_map = int(M_zero)

    do il = 1, nr_of_boxes
      ip_in = 0
      ip_bd = 0
      do ip = 1, gr%mesh%np
        xx(1:3) = gr%mesh%x(ip,1:3)
        if (check_point_in_bounds(xx, bounds(il,:,:))) then
          ip_in = ip_in + 1
          if (any(medium_box%points_map == ip)) then
            message(1) = 'Linear media boxes overlap.'
            call messages_fatal(1, namespace=namespace)
          else
            medium_box%points_map(ip_in,il) = ip
          end if
        end if
        if (check_point_on_bounds(xx, bounds(il,:,:))) then
          ip_bd = ip_bd + 1
          medium_box%bdry_map(ip_bd,il) = ip
        end if
      end do
    end do

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

      ! print information about the medium box -- get from Renes version in maxwell_propagator.F90

    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_grad)

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

    subroutine get_medium_io_function(medium_func, medium_box, mesh, il, io_func)
      FLOAT,                    intent(in)    :: medium_func(:)
      type(medium_box_t),       intent(in)    :: medium_box
      type(mesh_t),             intent(in)    :: mesh
      integer,                  intent(in)    :: il
      FLOAT,                    intent(inout) :: io_func(:)

      integer :: ip, ip_in

      do ip_in = 1, medium_box%points_number(il)
        ip = medium_box%points_map(ip_in, il)
        io_func(ip) = medium_func(ip_in)
      end do

    end subroutine get_medium_io_function

  end subroutine generate_medium_boxes

  ! ---------------------------------------------------------
  !> Deallocation of medium_box components
  subroutine medium_box_end(medium_box)
    class(medium_box_t),   intent(inout)    :: medium_box

    PUSH_SUB(medium_box_end)

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

    POP_SUB(medium_box_end)

  end subroutine medium_box_end

end module medium_mxll_oct_m
