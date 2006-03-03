!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module mesh_m
  use global_m
  use mpi_m
  use messages_m
  use lib_oct_parser_m
  use units_m
  use math_m
  use geometry_m
  use curvlinear_m
  use simul_box_m
  use lib_adv_alg_m
  use io_m
  use par_vec_m
  use multicomm_m

  implicit none

  private
  public ::            &
    mesh_t,            &
    mesh_init_stage_1, &
    mesh_init_stage_2, &
    mesh_init_stage_3, &
    mesh_end,          &
    mesh_double_box,   &
    mesh_inborder,     &
    mesh_r,            &
    mesh_gcutoff,      &
    mesh_write_info


  ! Describes mesh distribution to nodes.

  ! Some general thing:
  ! All members of type(mesh_t) are equal on all
  ! nodes when running parallel except
  ! - np, np_part
  ! - x, vol_pp
  ! These four are defined for all the points the node is responsible for.
  type mesh_t
    type(simul_box_t), pointer :: sb
    logical :: use_curvlinear

    FLOAT :: h(3)              ! the (constant) spacing between the points

    ! When running serially, the local number of points is
    ! equal to the global number of points.
    ! Otherwise, the next two are different on each node.
    integer  :: np             ! Local number of points in mesh
    integer  :: np_part        ! Local points plus ghost points plus
                               ! boundary points.
    integer  :: np_global      ! Global number of points in mesh.
    integer  :: np_part_global ! Global number of inner points and boundary points.

    integer  :: enlarge(3)     ! number of points to add for boundary conditions

    integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
    integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz

    FLOAT,   pointer :: x_tmp(:,:,:,:)  ! temporary arrays that we have to keep between calls to
    integer, pointer :: Lxyz_tmp(:,:,:) ! init_1 and init_2

    logical            :: parallel_in_domains ! will I run parallel in domains?
    type(mpi_grp_t) :: mpi_grp             ! the mpi group describing parallelization in domains
    type(pv_t)      :: vp                  ! describes parallel vectors defined on the mesh.

    ! some other vars
    integer :: nr(2,3)                  ! dimensions of the box where the points are contained
    integer :: l(3)                     ! literally n(2,:) - n(1,:) + 1 - 2*enlarge(:)

    FLOAT, pointer :: x(:,:)            ! The (local) points,
    FLOAT, pointer :: x_global(:,:)     ! The global points, needed for i/o on
                                        ! the root node and for the poisson solver
                                        ! on all nodes.
                                        ! There is a redundancy in these two
                                        ! entries.
                                        ! In serial: x_global => x.
    FLOAT, pointer :: vol_pp(:)         ! Element of volume for integrations
                                        ! for local points.

  end type mesh_t

contains


  ! finds the dimension of a box doubled in the non-periodic dimensions
  subroutine mesh_double_box(sb, m, db)
    type(simul_box_t), intent(in)  :: sb
    type(mesh_t),      intent(in)  :: m
    integer,              intent(out) :: db(3)

    integer :: i

    db = 1

    ! double mesh with 2n points
    do i = 1, sb%periodic_dim
      db(i) = m%l(i)
    end do
    do i = sb%periodic_dim + 1, sb%dim
      db(i) = nint(sb%fft_alpha*(m%l(i)-1)) + 1
    end do

  end subroutine mesh_double_box


  ! ---------------------------------------------------------
  subroutine mesh_write_info(m, unit)
    type(mesh_t), intent(in) :: m
    integer,         intent(in) :: unit

    if(.not.mpi_grp_is_root(mpi_world)) return

    call push_sub('mesh.mesh_write_info')

    write(message(1),'(a)') 'Main mesh:'

    write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)')  &
       '  Spacing [', trim(units_out%length%abbrev), '] = ',              &
       '(', m%h(1)/units_out%length%factor, ',',                          &
            m%h(2)/units_out%length%factor, ',',                          &
            m%h(3)/units_out%length%factor, ')',                          &
       '   volume/point [', trim(units_out%length%abbrev), '^3] = ',      &
       m%vol_pp(1)/units_out%length%factor**calc_dim

    write(message(3),'(a, i8)') '  # inner mesh = ', m%np_global
    write(message(4),'(a, i8)') '  # total mesh = ', m%np_part_global

    write(message(5),'(3a,f9.3,a)') '  Grid Cutoff [',trim(units_out%energy%abbrev),'] = ', &
       (M_PI**2/(M_TWO*maxval(m%h)**2))/units_out%energy%factor
    call write_info(5, unit)

    call pop_sub()
  end subroutine mesh_write_info


  ! ---------------------------------------------------------
  subroutine mesh_r(m, i, r, a, x)
    type(mesh_t),      intent(in)  :: m
    integer,              intent(in)  :: i
    FLOAT,                intent(out) :: r
    FLOAT,                intent(in),  optional :: a(:) ! a(sb%dim)
    FLOAT,                intent(out), optional :: x(:) ! x(sb%dim)

    FLOAT :: xx(m%sb%dim)

    xx(:) = m%x(i, 1:m%sb%dim)
    if(present(a)) xx(:) = xx(:) - a(1:m%sb%dim)
    r = sqrt(dot_product(xx, xx))

    if(present(x)) x(1:m%sb%dim) = xx(:)
  end subroutine mesh_r


  !/*---------------------------------------------------------------------
  ! Finds out if a given point of a mesh belongs to the "border" of the
  ! mesh. A point belongs to the border of the mesh if it is too close
  ! to any of the walls of the mesh. The criterium is set by input
  ! parameter "width".
  !
  ! m     : the mesh.
  ! i     : the point in the mesh.
  ! n     : on output, the number (0<=n<=3) of "walls" of the mesh that
  !         the point is too close to, in order to consider it belonging
  !         to a mesh.
  ! d     : the distances of the point to the walls, for each of the walls
  !         that the point is too close to.
  ! width : the width of the border.
  !
  ! So, if n>0, the point is in the border.
  ! ----------------------------------------------------------------------*/
  subroutine mesh_inborder(m, i, n, d, width)
    type(mesh_t),      intent(in)  :: m
    integer,              intent(in)  :: i
    FLOAT,                intent(in)  :: width
    integer,              intent(out) :: n
    FLOAT,                intent(out) :: d(3)

    integer :: j
    FLOAT :: x(3), r, dd

    call mesh_r(m, i, r, x=x)
    n = 0
    select case(m%sb%box_shape)
    case(SPHERE)
      dd = r - (m%sb%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      end if
    case(CYLINDER)
      dd = sqrt(x(2)**2 + x(3)**2) - (m%sb%rsize - width)
      if(dd.gt.M_ZERO) then
        n = 1; d(1) = dd
      end if
      dd = abs(x(1)) - (m%sb%xsize - width)
      if(dd.gt.M_ZERO) then
        n = n + 1; d(n) = dd
      end if
    case(MINIMUM,BOX_USDEF)
      message(1) = "Absorbing boundaries are not yet implemented for the 'minimum' box"
      call write_fatal(1)
    case(PARALLELEPIPED)
      do j = 1, m%sb%dim
        dd = abs(x(j)) - (m%sb%lsize(j) - width)
        if(dd.gt.M_ZERO) then
           n = n + 1; d(n) = dd
         end if
       end do
     end select

   end subroutine mesh_inborder


   !/*-------------------------------------------------------------
   ! mesgh_gcutoff returns the "natural" band limitation of the
   ! grid m, in terms of the maximum G vector. For a cubic regular
   ! grid of spacing h is M_PI/h.
   ! --------------------------------------------------------------*/
   FLOAT function mesh_gcutoff(m) result(gmax)
     type(mesh_t), intent(in) :: m
     gmax = M_PI/(maxval(m%h))
   end function mesh_gcutoff


   subroutine mesh_end(m)
     type(mesh_t), intent(inout) :: m

     call push_sub('mesh.mesh_end')

     if(associated(m%Lxyz)) then
       deallocate(m%Lxyz, m%Lxyz_inv, m%x, m%vol_pp)
       nullify(m%Lxyz, m%Lxyz_inv, m%x, m%vol_pp)
     end if

     if(m%parallel_in_domains) then
       if(associated(m%x_global).and.m%vp%rank.eq.m%vp%root) then
         deallocate(m%x_global)
         nullify(m%x_global)
       end if
#if defined(HAVE_MPI)
       call vec_end(m%vp)
#endif
     end if

     call pop_sub()
   end subroutine mesh_end

#include "mesh_init.F90"

end module mesh_m
