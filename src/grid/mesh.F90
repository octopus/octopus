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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module mesh_m
  use datasets_m
  use curvlinear_m
  use geometry_m
  use global_m
  use io_m
  use math_m
  use mesh_lib_m
  use messages_m
  use multicomm_m
  use mpi_m
  use loct_parser_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use units_m

  implicit none
  
  private
  public ::                    &
    mesh_t,                    &
    mesh_plane_t,              &
    mesh_line_t,               &
    mesh_init_from_file,       &
    mesh_lxyz_init_from_file,  &
    mesh_dump,                 &
    mesh_lxyz_dump,            &
    mesh_end,                  &
    mesh_double_box,           &
    mesh_inborder,             &
    mesh_r,                    &
    mesh_gcutoff,              &
    mesh_write_info,           &
    mesh_nearest_point,        &
    mesh_subset_indices,       &
    mesh_periodic_point,       &
    translate_point
  
  ! Describes mesh distribution to nodes.
  
  ! Some general things:
  ! All members of type(mesh_t) are equal on all
  ! nodes when running parallel except
  ! - np, np_part
  ! - x, vol_pp
  ! These four are defined for all the points the node is responsible for.
  type mesh_t
    type(simul_box_t), pointer :: sb
    logical :: use_curvlinear
    
    FLOAT :: h(MAX_DIM)         ! the (constant) spacing between the points
    
    ! When running serially, the local number of points is
    ! equal to the global number of points.
    ! Otherwise, the next two are different on each node.
    integer  :: np               ! Local number of points in mesh
    integer  :: np_part          ! Local points plus ghost points plus
    ! boundary points.
    integer  :: np_global        ! Global number of points in mesh.
    integer  :: np_part_global   ! Global number of inner points and boundary points.
    
    integer  :: enlarge(MAX_DIM) ! number of points to add for boundary conditions
    
    integer, pointer :: Lxyz(:,:)       ! return x, y and z for each point
    integer, pointer :: Lxyz_inv(:,:,:) ! return points # for each xyz
    
    FLOAT,   pointer :: x_tmp(:,:,:,:)  ! temporary arrays that we have to keep between calls to
    integer, pointer :: Lxyz_tmp(:,:,:) ! init_1 and init_2
    
    integer, pointer :: boundary_indices(:,:) ! contains the list of mesh indices for boundary points
    integer          :: boundary_np(6)        ! total number of boundary points
    
    logical         :: parallel_in_domains ! will I run parallel in domains?
    type(mpi_grp_t) :: mpi_grp             ! the mpi group describing parallelization in domains
    type(pv_t)      :: vp                  ! describes parallel vectors defined on the mesh.
    
    ! some other vars
    integer :: nr(2, MAX_DIM)              ! dimensions of the box where the points are contained
    integer :: l(MAX_DIM)                  ! literally n(2,:) - n(1,:) + 1 - 2*enlarge(:)
    
    FLOAT, pointer :: x(:,:)            ! The (local) points,
    FLOAT, pointer :: x_global(:,:)     ! The global points, needed for i/o on
    ! the root node and for the poisson solver
    ! on all nodes.
    ! There is a redundancy in these two
    ! entries.
    ! In serial: x_global => x.
    integer, pointer :: resolution(:, :, :)
    FLOAT, pointer :: vol_pp(:)         ! Element of volume for integrations
    ! for local points.
    integer :: nper                     ! the number of points that correpond to pbc
    integer, pointer :: per_points(:)   ! (1:nper) the list of points that correspond to pbc 
    integer, pointer :: per_map(:)      ! (1:nper) the inner point that corresponds to each pbc point
#ifdef HAVE_MPI
    integer, pointer :: nsend(:)
    integer, pointer :: nrecv(:)
    integer, pointer :: dsend_type(:)
    integer, pointer :: zsend_type(:)
    integer, pointer :: drecv_type(:)
    integer, pointer :: zrecv_type(:)
#endif
    
    type(mesh_t), pointer :: lead_unit_cell(:) ! Meshes of the lead unit cells for open boundary.
                                               ! calculations.
  end type mesh_t
  
  ! This data type defines a plane, and a regular grid defined on 
  ! this plane (or, rather, on a portion of this plane)
  ! n should be a unit vector, that determines the normal of the plane.
  ! origin is a point belonging to the plane
  ! u and v are unit orthogonal vectors belonging to the plane
  ! The grid is generated by the vectors u and v:
  !   x_{i,j} = origin + i*spacing*u + j*spacing*v,
  ! for nu <= i <= mu and nv <= j <= mv
  type mesh_plane_t
    FLOAT :: n(MAX_DIM)
    FLOAT :: u(MAX_DIM), v(MAX_DIM)
    FLOAT :: origin(MAX_DIM)
    FLOAT :: spacing
    integer :: nu, mu, nv, mv
  end type mesh_plane_t
  
  ! This data type defines a line, and a regular grid defined on this
  ! line (or rather, on a portion of this line).
  type mesh_line_t
    FLOAT :: n(MAX_DIM)
    FLOAT :: u(MAX_DIM)
    FLOAT :: origin(MAX_DIM)
    FLOAT :: spacing
    integer :: nu, mu
  end type mesh_line_t
  
  integer, parameter, public ::      &
    LEFT_BOUNDARY_X   =  1,          &
    RIGHT_BOUNDARY_X  =  2,          &
    LEFT_BOUNDARY_Y   =  3,          &
    RIGHT_BOUNDARY_Y  =  4,          &
    LEFT_BOUNDARY_Z   =  5,          &
    RIGHT_BOUNDARY_Z  =  6,          &
    MAX_BOUNDARY_DIM  = RIGHT_BOUNDARY_Z
  
  character(len=17), parameter :: dump_tag = '*** mesh_dump ***'
  
contains

  ! ---------------------------------------------------------
  ! finds the dimension of a box doubled in the non-periodic dimensions
  subroutine mesh_double_box(sb, m, db)
    type(simul_box_t), intent(in)  :: sb
    type(mesh_t),      intent(in)  :: m
    integer,           intent(out) :: db(MAX_DIM)
    
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
    integer,      intent(in) :: unit
    
    if(.not.mpi_grp_is_root(mpi_world)) return
    
    call push_sub('mesh.mesh_write_info')
    
    write(message(1),'(a)') 'Main mesh:'
    
    write(message(2),'(3a, a, f6.3, a, f6.3, a, f6.3, a, 1x, 3a, f8.5)')  &
      '  Spacing [', trim(units_out%length%abbrev), '] = ',              &
      '(', m%h(1)/units_out%length%factor, ',',                          &
      m%h(2)/units_out%length%factor, ',',                          &
      m%h(3)/units_out%length%factor, ')',                          &
      '   volume/point [', trim(units_out%length%abbrev), '^3] = ',      &
      m%vol_pp(1)/units_out%length%factor**m%sb%dim
    
    write(message(3),'(a, i8)') '  # inner mesh = ', m%np_global
    write(message(4),'(a, i8)') '  # total mesh = ', m%np_part_global
    
    write(message(5),'(3a,f9.3,a)') '  Grid Cutoff [',trim(units_out%energy%abbrev),'] = ', &
      (M_PI**2/(M_TWO*maxval(m%h)**2))/units_out%energy%factor
    call write_info(5, unit)
    
    call pop_sub()
  end subroutine mesh_write_info
  
  
  ! ---------------------------------------------------------
  subroutine mesh_r(m, i, r, a, x)
    type(mesh_t), intent(in)  :: m
    integer,      intent(in)  :: i
    FLOAT,        intent(out) :: r
    FLOAT,        intent(in),  optional :: a(:) ! a(sb%dim)
    FLOAT,        intent(out), optional :: x(:) ! x(sb%dim)
    
    FLOAT :: xx(MAX_DIM)
    
    xx(1:m%sb%dim) = m%x(i, 1:m%sb%dim)
    if(present(a)) xx(1:m%sb%dim) = xx(1:m%sb%dim) - a(1:m%sb%dim)
    r = sqrt(dot_product(xx(1:m%sb%dim), xx(1:m%sb%dim)))
    
    if(present(x)) then
      x(1:MAX_DIM) = M_ZERO
      x(1:m%sb%dim) = xx(1:m%sb%dim)
    end if
  end subroutine mesh_r
  
  
  !---------------------------------------------------------------------
  ! Finds out if a given point of a mesh belongs to the "border" of the
  ! mesh. A point belongs to the border of the mesh if it is too close
  ! to any of the walls of the mesh. The criterion is set by input
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
  ! ----------------------------------------------------------------------
  subroutine mesh_inborder(m, i, n, d, width)
    type(mesh_t), intent(in)  :: m
    integer,      intent(in)  :: i
    FLOAT,        intent(in)  :: width
    integer,      intent(out) :: n
    FLOAT,        intent(out) :: d(MAX_DIM)
    
    integer :: j
    FLOAT   :: x(MAX_DIM), r, dd
    
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
      if ( m%sb%periodic_dim.eq.0 ) then
        dd = abs(x(1)) - (m%sb%xsize - width)
        if(dd.gt.M_ZERO) then
          n = n + 1; d(n) = dd
        end if
      end if
    case(MINIMUM,BOX_USDEF)
      message(1) = "Absorbing boundaries are not yet implemented for the 'minimum' box"
      call write_fatal(1)
    case(PARALLELEPIPED)
      do j = m%sb%periodic_dim+1, m%sb%dim
        dd = abs(x(j)) - (m%sb%lsize(j) - width)
        if(dd.gt.M_ZERO) then
          n = n + 1; d(n) = dd
        end if
      end do
    end select
    
  end subroutine mesh_inborder
  
  
  !---------------------------------------------------------------------
  ! Returns the index of the point which is nearest to a given vector
  ! position pos. Variable dmin will hold, on exit, the distance between
  ! pos and this nearest mesh point. rankmin will be zero, if the mesh is
  ! not partitioned, and the rank of the processor which holds the point
  ! ind if the mesh is partitioned.
  ! ----------------------------------------------------------------------
  integer function mesh_nearest_point(mesh, pos, dmin, rankmin) result(ind)
    type(mesh_t), intent(in)  :: mesh
    FLOAT,        intent(in)  :: pos(MAX_DIM)
    FLOAT,        intent(out) :: dmin
    integer,      intent(out) :: rankmin
    
    FLOAT :: d
    integer :: imin, i
#if defined(HAVE_MPI)
    FLOAT :: min_loc_in(2), min_loc_out(2)
#endif
    
    call push_sub('mesh.mesh_nearest_point')
    
    !find the point of the grid that is closer to the atom
    dmin = M_ZERO
    do i = 1, mesh%np
      d = sum((pos(1:mesh%sb%dim) - mesh%x(i, 1:mesh%sb%dim))**2)
      if((d < dmin) .or. (i == 1)) then 
        imin = i
        dmin = d 
      end if
    end do
    
    rankmin = 0
#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      min_loc_in(1) = dmin
      min_loc_in(2) = mesh%np_global * mesh%mpi_grp%rank  + real(imin, REAL_PRECISION) 
      call mpi_allreduce(min_loc_in, min_loc_out, 1, MPI_DOUBLE_PRECISION, &
        MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)
      dmin = min_loc_out(1)
      imin = mod(nint(min_loc_out(2)), mesh%np_global)
      rankmin = nint(min_loc_out(2))/mesh%np_global
    end if
#endif
    
    ind = imin
    call pop_sub()
  end function mesh_nearest_point
  
  
  !--------------------------------------------------------------
  integer function translate_point(mesh, index, tdist) result (tindex)
    type(mesh_t), intent(in) :: mesh
    integer,      intent(in) :: index
    integer,      intent(in) :: tdist(:)
    
    integer :: ixyz(MAX_DIM), id
    
    call push_sub('mesh.translate_point')
    
    ASSERT(index >= 1 .and. index <= mesh%np)
    
    ixyz(:) = mesh%Lxyz(index, :)
    
    do id = 1, mesh%sb%dim
      ixyz(id) = ixyz(id) + tdist(id)
    end do
    
    tindex = mesh_index(mesh%sb%dim, mesh%nr, mesh%Lxyz_inv, ixyz) 
    
    ! check if the translated point is still inside the domain and return
    ! a negative index otherwise to indicate that the requested translation
    ! is not valid
    if (tindex < 1 .or. tindex > mesh%np) then
      tindex = -1
    end if
    
    call pop_sub()
  end function translate_point
  
  
  ! --------------------------------------------------------------
  ! mesh_gcutoff returns the "natural" band limitation of the
  ! grid m, in terms of the maximum G vector. For a cubic regular
  ! grid of spacing h is M_PI/h.
  ! --------------------------------------------------------------
  FLOAT function mesh_gcutoff(m) result(gmax)
    type(mesh_t), intent(in) :: m

    gmax = M_PI/(maxval(m%h))
  end function mesh_gcutoff
  
  
  ! -------------------------------------------------------------- 
  subroutine mesh_dump(mesh, iunit)
    type(mesh_t), intent(in) :: mesh
    integer,      intent(in) :: iunit
    
    call push_sub('mesh.mesh_dump')
    
    write(iunit, '(a)')         dump_tag
    write(iunit, '(a20,7i8)')   'nr(1, :)=           ', mesh%nr(1, 1:mesh%sb%dim)
    write(iunit, '(a20,7i8)')   'nr(2, :)=           ', mesh%nr(2, 1:mesh%sb%dim)
    write(iunit, '(a20,7i8)')   'l(:)=               ', mesh%l(1:mesh%sb%dim)
    write(iunit, '(a20,7i8)')   'enlarge(:)=         ', mesh%enlarge(1:mesh%sb%dim)
    write(iunit, '(a20,1i10)')  'np=                 ', mesh%np
    write(iunit, '(a20,1i10)')  'np_part=            ', mesh%np_part
    write(iunit, '(a20,1i10)')  'np_global=          ', mesh%np_global
    write(iunit, '(a20,1i10)')  'np_part_global=     ', mesh%np_part_global
    
    call pop_sub()
  end subroutine mesh_dump
  
  
  ! -------------------------------------------------------------- 
  ! Read the mesh parameters from file that were written by mesh_dump.
  subroutine mesh_init_from_file(mesh, iunit)
    type(mesh_t), intent(inout) :: mesh
    integer,      intent(in)    :: iunit

    character(len=20)  :: str
    character(len=100) :: line

    call push_sub('mesh.mesh_init_from_file')

    ! Find (and throw away) the dump tag.
    do
      read(iunit, '(a)') line
      if(trim(line).eq.dump_tag) exit
    end do

    read(iunit, '(a20,7i8)')  str, mesh%nr(1, 1:mesh%sb%dim)
    read(iunit, '(a20,7i8)')  str, mesh%nr(2, 1:mesh%sb%dim)
    read(iunit, '(a20,7i8)')  str, mesh%l(1:mesh%sb%dim)
    read(iunit, '(a20,7i8)')  str, mesh%enlarge(1:mesh%sb%dim)
    read(iunit, '(a20,1i10)') str, mesh%np
    read(iunit, '(a20,1i10)') str, mesh%np_part
    read(iunit, '(a20,1i10)') str, mesh%np_global
    read(iunit, '(a20,1i10)') str, mesh%np_part_global
    nullify(mesh%lxyz, mesh%lxyz_inv, mesh%x_tmp, mesh%lxyz_tmp, mesh%boundary_indices, mesh%x, &
      mesh%x_global, mesh%vol_pp, mesh%per_points, mesh%per_map, mesh%lead_unit_cell)
    mesh%parallel_in_domains = .false.

    call pop_sub()
  end subroutine mesh_init_from_file


  ! --------------------------------------------------------------
  subroutine mesh_Lxyz_dump(mesh, iunit)
    type(mesh_t), intent(in) :: mesh
    integer,      intent(in) :: iunit

    integer :: ip

    call push_sub('mesh.Lxyz_dump')

    do ip = 1, mesh%np_part
      write(iunit, '(7i8)') mesh%Lxyz(ip, 1:mesh%sb%dim)
    end do

    call pop_sub()
  end subroutine mesh_Lxyz_dump


  ! --------------------------------------------------------------
  ! Fill the lxyz and lxyz_inv arrays from a file
  subroutine mesh_lxyz_init_from_file(mesh, iunit)
    type(mesh_t), intent(inout) :: mesh
    integer,      intent(in)    :: iunit

    integer :: ip

    call push_sub('mesh.mesh_lxyz_init_from_file')

    do ip = 1, mesh%np_part
      read(iunit, '(7i8)') mesh%lxyz(ip, 1:mesh%sb%dim)
      mesh%lxyz_inv(mesh%lxyz(ip, 1), mesh%lxyz(ip, 2), mesh%lxyz(ip, 3)) = ip
    end do

    call pop_sub()
  end subroutine mesh_lxyz_init_from_file


  ! --------------------------------------------------------------
  ! Extracts the point numbers of a rectangular subset spanned
  ! by the two corner points from and to.
  subroutine mesh_subset_indices(mesh, from, to, indices)
    type(mesh_t), intent(in)  :: mesh
    integer,      intent(in)  :: from(MAX_DIM)
    integer,      intent(in)  :: to(MAX_DIM)
    integer,      intent(out) :: indices(:)

    integer :: lb(MAX_DIM) ! Lower bound of indices.
    integer :: ub(MAX_DIM) ! Upper bound of indices.

    integer :: ix, iy, iz, i

    call push_sub('mesh.mesh_subset_indices')

    ! In debug mode, check for valid indices in from, to first.
    if(in_debug_mode) then
      if(.not.index_valid(mesh, from).or..not.index_valid(mesh, to)) then
        message(1) = 'Failed assertion:'
        message(2) = 'mesh.mesh_subset_indices has been passed points outside the box:'
        message(3) = ''
        write(message (4), '(a, i6, a, i6, a, i6, a)') &
          '  from = (', from(1), ', ', from(2), ', ', from(3), ')'
        write(message(5), '(a, i6, a, i6, a, i6, a)') & 
          '  to   = (', to(1), ', ', to(2), ', ', to(3), ')'
        call write_fatal(5)
      end if
    end if

    lb = min(from, to)
    ub = max(from, to)

    i = 1
    do ix = lb(1), ub(1)
      do iy = lb(2), ub(2)
        do iz = lb(3), ub(3)
          indices(i) = mesh%Lxyz_inv(ix, iy, iz)
          i          = i + 1
        end do
      end do
    end do

    call pop_sub()
  end subroutine mesh_subset_indices


  ! --------------------------------------------------------------
  ! Checks if the (x, y, z) indices of point are valid, i. e.
  ! inside the dimensions of the simulation box.
  logical function index_valid(mesh, point)
    type(mesh_t), intent(in) :: mesh
    integer,      intent(in) :: point(MAX_DIM)

    integer :: i
    logical :: valid

    call push_sub('mesh.index_valid')

    valid = .true.
    do i = 1, 3
      if(point(i).lt.mesh%nr(1, i).or.point(i).gt.mesh%nr(2, i)) then
        valid = .false.
      end if
    end do

    index_valid = valid

    call pop_sub()
  end function index_valid


  ! --------------------------------------------------------------
  recursive subroutine mesh_end(m)
    type(mesh_t), intent(inout) :: m

    integer :: il, ipart
    call push_sub('mesh.mesh_end')

    DEALLOC(m%resolution)
    DEALLOC(m%lxyz)
    DEALLOC(m%lxyz_inv)
    DEALLOC(m%x)
    DEALLOC(m%vol_pp)

    if(m%parallel_in_domains) then
      DEALLOC(m%x_global)
#if defined(HAVE_MPI)
      call vec_end(m%vp)
      if(simul_box_is_periodic(m%sb)) then

        do ipart = 1, m%vp%p
          if(m%nsend(ipart) /= 0) then 
            call MPI_Type_free(m%dsend_type(ipart), mpi_err)
            call MPI_Type_free(m%zsend_type(ipart), mpi_err)
          end if
          if(m%nrecv(ipart) /= 0) then 
            call MPI_Type_free(m%drecv_type(ipart), mpi_err)
            call MPI_Type_free(m%zrecv_type(ipart), mpi_err)
          end if
        end do

        deallocate(m%dsend_type, m%zsend_type)
        deallocate(m%drecv_type, m%zrecv_type)
        deallocate(m%nsend, m%nrecv)
      end if
#endif
    end if

    DEALLOC(m%per_points)
    DEALLOC(m%per_map)

    if(associated(m%lead_unit_cell)) then
      do il = 1, NLEADS
        call mesh_end(m%lead_unit_cell(il))
      end do
    end if
    
    call pop_sub()
  end subroutine mesh_end
  
  ! This function returns the point inside the grid corresponding to
  ! a boundary point when PBC are used. In case the point does not
  ! have a correspondance (i.e. other BC are used in that direction),
  ! the same point is returned. Note that this function returns a
  ! global point number when parallelization in domains is used.
  
  integer function mesh_periodic_point(m, ip) result(ipp)
    type(mesh_t), intent(in)    :: m
    integer,      intent(in)    :: ip
    
    integer :: ix(MAX_DIM), nr(2, MAX_DIM), idim
    
    ix = m%Lxyz(ip, :)
    nr(1, :) = m%nr(1, :) + m%enlarge(:)
    nr(2, :) = m%nr(2, :) - m%enlarge(:)
    
    do idim = 1, m%sb%periodic_dim
      if(ix(idim) < nr(1, idim)) ix(idim) = ix(idim) + m%l(idim)
      if(ix(idim) > nr(2, idim)) ix(idim) = ix(idim) - m%l(idim)
    end do
    
    ipp = m%Lxyz_inv(ix(1), ix(2), ix(3))
    
  end function mesh_periodic_point
  
end module mesh_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
