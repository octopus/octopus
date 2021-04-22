!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module mesh_oct_m
  use basis_set_abst_oct_m
  use box_hypercube_oct_m
  use comm_oct_m
  use curvilinear_oct_m
  use global_oct_m
  use hypercube_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use ions_oct_m
  use mesh_cube_map_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use partition_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none
  
  private
  public ::                        &
    mesh_t,                        &
    mesh_plane_t,                  &
    mesh_line_t,                   &
    mesh_check_dump_compatibility, &
    mesh_end,                      &
    mesh_double_box,               &
    mesh_r,                        &
    mesh_gcutoff,                  &
    mesh_write_info,               &
    mesh_nearest_point,            &
    mesh_periodic_point,           &
    mesh_global_memory,            &
    mesh_local_memory,             &
    mesh_x_global,                 &
    mesh_write_fingerprint,        &
    mesh_read_fingerprint,         &
    mesh_compact_boundaries,       &
    mesh_check_symmetries,         &
    mesh_global_index_to_coords,   &
    mesh_global_index_from_coords, &
    mesh_local_index_to_coords,    &
    mesh_local_index_from_coords,  &
    mesh_local2global,             &
    mesh_global2local

  !> Describes mesh distribution to nodes.
  !!
  !! Some general things:
  !! All members of type(mesh_t) are equal on all
  !! nodes when running parallel except
  !! - np, np_part
  !! - x, vol_pp
  !! These four are defined for all the points the node is responsible for.
  type, extends(basis_set_abst_t) :: mesh_t
    ! Components are public by default
    type(simul_box_t),   pointer :: sb  !< simulation box
    type(curvilinear_t), pointer :: cv  
    type(index_t)                :: idx 
    logical :: use_curvilinear
    
    FLOAT :: spacing(MAX_DIM)         !< the (constant) spacing between the points
    
    !> When running serially, the local number of points is
    !! equal to the global number of points.
    !! Otherwise, the next two are different on each node.
    integer  :: np               !< Local number of points in mesh
    integer  :: np_part          !< Local points plus ghost points plus boundary points.
    integer  :: np_global        !< Global number of points in mesh.
    integer  :: np_part_global   !< Global number of inner points and boundary points.
    !> will I run parallel in domains?
    !! yes or no??
    logical         :: parallel_in_domains 
    type(mpi_grp_t) :: mpi_grp             !< the mpi group describing parallelization in domains
    type(pv_t)      :: vp                  !< describes parallel vectors defined on the mesh.
    type(partition_t) :: partition         !< describes how the inner points are assigned to the domains

    FLOAT,   allocatable :: x(:,:)            !< The (local) \b points
    FLOAT                :: volume_element    !< The global volume element.
    FLOAT                :: surface_element(MAX_DIM)
    FLOAT,   allocatable :: vol_pp(:)         !< Element of volume for curvilinear coordinates.

    type(mesh_cube_map_t) :: cube_map

    logical :: masked_periodic_boundaries
    character(len=256) :: periodic_boundary_mask

  contains
    procedure :: end => mesh_end
    procedure :: init => mesh_init
    procedure :: write_info => mesh_write_info
    procedure :: dmesh_allreduce_0, zmesh_allreduce_0, imesh_allreduce_0
    procedure :: dmesh_allreduce_1, zmesh_allreduce_1, imesh_allreduce_1
    procedure :: dmesh_allreduce_2, zmesh_allreduce_2, imesh_allreduce_2
    procedure :: dmesh_allreduce_3, zmesh_allreduce_3, imesh_allreduce_3
    procedure :: dmesh_allreduce_4, zmesh_allreduce_4, imesh_allreduce_4
    procedure :: dmesh_allreduce_5, zmesh_allreduce_5, imesh_allreduce_5
    generic :: allreduce => dmesh_allreduce_0, zmesh_allreduce_0, imesh_allreduce_0
    generic :: allreduce => dmesh_allreduce_1, zmesh_allreduce_1, imesh_allreduce_1
    generic :: allreduce => dmesh_allreduce_2, zmesh_allreduce_2, imesh_allreduce_2
    generic :: allreduce => dmesh_allreduce_3, zmesh_allreduce_3, imesh_allreduce_3
    generic :: allreduce => dmesh_allreduce_4, zmesh_allreduce_4, imesh_allreduce_4
    generic :: allreduce => dmesh_allreduce_5, zmesh_allreduce_5, imesh_allreduce_5
  end type mesh_t
  
  !> This data type defines a plane, and a regular grid defined on 
  !! this plane (or, rather, on a portion of this plane)
  !! n should be a unit vector, that determines the normal of the plane.
  !! Origin is a point belonging to the plane
  !! u and v are unit orthogonal vectors belonging to the plane
  !! The grid is generated by the vectors u and v:
  !!   x_{i,j} = origin + i*spacing*u + j*spacing*v,
  !! for nu <= i <= mu and nv <= j <= mv
  type mesh_plane_t
    ! Components are public by default
    FLOAT :: n(MAX_DIM)
    FLOAT :: u(MAX_DIM), v(MAX_DIM)
    FLOAT :: origin(MAX_DIM)
    FLOAT :: spacing
    integer :: nu, mu, nv, mv
  end type mesh_plane_t
  
  !> This data type defines a line, and a regular grid defined on this
  !! line (or rather, on a portion of this line).
  type mesh_line_t
    ! Components are public by default
    FLOAT :: n(MAX_DIM)
    FLOAT :: u(MAX_DIM)
    FLOAT :: origin(MAX_DIM)
    FLOAT :: spacing
    integer :: nu, mu
  end type mesh_line_t
  
contains

  subroutine mesh_init(this)
    class(mesh_t), intent(inout) :: this

    PUSH_SUB(mesh_init)

    call this%set_time_dependent(.false.)

    POP_SUB(mesh_init)
  end subroutine mesh_init

! ---------------------------------------------------------
  !> finds the dimension of a box doubled in the non-periodic dimensions
  subroutine mesh_double_box(space, mesh, alpha, db)
    type(space_t),     intent(in)  :: space
    type(mesh_t),      intent(in)  :: mesh
    FLOAT,             intent(in)  :: alpha !< enlargement factor for double box
    integer,           intent(out) :: db(MAX_DIM)

    integer :: idir
    
    PUSH_SUB(mesh_double_box)

    db = 1
    
    ! double mesh with 2n points
    do idir = 1, space%periodic_dim
      db(idir) = mesh%idx%ll(idir)
    end do
    do idir = space%periodic_dim + 1, space%dim
      db(idir) = nint(alpha * (mesh%idx%ll(idir) - 1)) + 1
    end do
    
    POP_SUB(mesh_double_box)
  end subroutine mesh_double_box
  
  
  ! ---------------------------------------------------------
  subroutine mesh_write_info(this, unit)
    class(mesh_t), intent(in) :: this
    integer,      intent(in) :: unit
    
    integer :: ii
    FLOAT :: cutoff

    if(.not.mpi_grp_is_root(mpi_world)) return
    
    PUSH_SUB(mesh_write_info)
    
    write(message(1),'(3a)') '  Spacing [', trim(units_abbrev(units_out%length)), '] = ('
    do ii = 1, this%sb%dim
      if(ii > 1) write(message(1), '(2a)') trim(message(1)), ','
      write(message(1), '(a,f6.3)') trim(message(1)), units_from_atomic(units_out%length, this%spacing(ii))
    end do
    write(message(1), '(5a,f12.5)') trim(message(1)), ') ', &
         '   volume/point [', trim(units_abbrev(units_out%length**this%sb%dim)), '] = ',      &
         units_from_atomic(units_out%length**this%sb%dim, this%vol_pp(1))
    
    write(message(2),'(a, i10)') '  # inner mesh = ', this%np_global
    write(message(3),'(a, i10)') '  # total mesh = ', this%np_part_global
    
    cutoff = mesh_gcutoff(this)**2 / M_TWO
    write(message(4),'(3a,f12.6,a,f12.6)') '  Grid Cutoff [', trim(units_abbrev(units_out%energy)),'] = ', &
      units_from_atomic(units_out%energy, cutoff), '    Grid Cutoff [Ry] = ', cutoff * M_TWO
    call messages_info(4, unit)
    
    POP_SUB(mesh_write_info)
  end subroutine mesh_write_info
  
  
  ! ---------------------------------------------------------
  subroutine mesh_r(mesh, ip, rr, origin, coords)
    type(mesh_t), intent(in)  :: mesh
    integer,      intent(in)  :: ip
    FLOAT,        intent(out) :: rr
    FLOAT,        intent(in),  optional :: origin(:) !< origin(sb%dim)
    FLOAT,        intent(out), optional :: coords(:) !< coords(sb%dim)
   
    FLOAT :: xx(1:mesh%sb%dim)

    ! no push_sub because it is called too frequently
    
    xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
    if(present(origin)) xx(1:mesh%sb%dim) = xx(1:mesh%sb%dim) - origin(1:mesh%sb%dim)
    rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
    
    if(present(coords)) then
      coords(1:mesh%sb%dim) = xx(1:mesh%sb%dim)
    end if

  end subroutine mesh_r
  
  !---------------------------------------------------------------------
  !> Returns the index of the point which is nearest to a given vector
  !! position pos. Variable dmin will hold, on exit, the distance between
  !! pos and this nearest mesh point. rankmin will be zero, if the mesh is
  !! not partitioned, and the rank of the processor which holds the point
  !! ind if the mesh is partitioned.
  ! ----------------------------------------------------------------------
  integer function mesh_nearest_point(mesh, pos, dmin, rankmin) result(ind)
    type(mesh_t), intent(in)  :: mesh
    FLOAT,        intent(in)  :: pos(MAX_DIM)
    FLOAT,        intent(out) :: dmin
    integer,      intent(out) :: rankmin
    
    FLOAT :: dd
    integer :: imin, ip
#if defined(HAVE_MPI)
    FLOAT :: min_loc_in(2), min_loc_out(2)
#endif
    
    PUSH_SUB(mesh_nearest_point)
    
    !find the point of the grid that is closer to the atom
    dmin = M_ZERO
    do ip = 1, mesh%np
      dd = sum((pos(1:mesh%sb%dim) - mesh%x(ip, 1:mesh%sb%dim))**2)
      if((dd < dmin) .or. (ip == 1)) then 
        imin = ip
        dmin = dd
      end if
    end do
    
    rankmin = 0
#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      min_loc_in(1) = dmin
      min_loc_in(2) = mesh%np_global * mesh%mpi_grp%rank  + TOFLOAT(imin) 
      call MPI_Allreduce(min_loc_in, min_loc_out, 1, MPI_2FLOAT, &
        MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)
      dmin = min_loc_out(1)
      imin = mod(nint(min_loc_out(2)), mesh%np_global)
      rankmin = nint(min_loc_out(2))/mesh%np_global
    end if
#endif
    
    ind = imin
    POP_SUB(mesh_nearest_point)
  end function mesh_nearest_point


  ! --------------------------------------------------------------
  !> mesh_gcutoff returns the "natural" band limitation of the
  !! grid mesh, in terms of the maximum G vector. For a cubic regular
  !! grid, it is M_PI/spacing.
  ! --------------------------------------------------------------
  FLOAT function mesh_gcutoff(mesh) result(gmax)
    type(mesh_t), intent(in) :: mesh

    PUSH_SUB(mesh_gcutoff)
    gmax = M_PI / (maxval(mesh%spacing))

    POP_SUB(mesh_gcutoff)
  end function mesh_gcutoff

  ! --------------------------------------------------------------
  subroutine mesh_write_fingerprint(mesh, dir, filename, mpi_grp, namespace, ierr)
    type(mesh_t),     intent(in)  :: mesh
    character(len=*), intent(in)  :: dir
    character(len=*), intent(in)  :: filename
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: ierr

    integer :: iunit

    PUSH_SUB(mesh_write_fingerprint)

    ierr = 0

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action='write', &
      die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1)
      ierr = ierr + 1
    else
      if (mpi_grp_is_root(mpi_grp)) then
        select type (box => mesh%sb%box)
        type is (box_hypercube_t)
          write(iunit, '(a20,l)')  'is_hypercube =         ', .true.
        class default
          write(iunit, '(a20,l)')  'is_hypercube =         ', .false.
          write(iunit, '(a20,i21)')  'np_part_global=     ', mesh%np_part_global
          write(iunit, '(a20,i21)')  'np_global=          ', mesh%np_global
          write(iunit, '(a20,i21)')  'algorithm=          ', 1
          write(iunit, '(a20,i21)')  'checksum=           ', mesh%idx%checksum
        end select
      end if
      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(mesh_write_fingerprint)
  end subroutine mesh_write_fingerprint


  ! -----------------------------------------------------------------------
  !> This function reads the fingerprint of a mesh written in
  !! filename. If the meshes are equal (same fingerprint) return values
  !! are 0, otherwise it returns the size of the mesh stored.
  !! fingerprint cannot be read, it returns ierr /= 0.
  subroutine mesh_read_fingerprint(mesh, dir, filename, mpi_grp, namespace, read_np_part, read_np, ierr)
    type(mesh_t),     intent(in)  :: mesh
    character(len=*), intent(in)  :: dir
    character(len=*), intent(in)  :: filename
    type(mpi_grp_t),  intent(in)  :: mpi_grp
    type(namespace_t),intent(in)  :: namespace
    integer,          intent(out) :: read_np_part
    integer,          intent(out) :: read_np
    integer,          intent(out) :: ierr

    character(len=20)  :: str
    logical :: is_hypercube
    character(len=100) :: lines(4)
    integer :: iunit, algorithm, err
    integer(8) :: checksum

    PUSH_SUB(mesh_read_fingerprint)

    ierr = 0

    read_np_part = 0
    read_np = 0

    iunit = io_open(trim(dir)//"/"//trim(filename), namespace, action='read', &
      status='old', die=.false., grp=mpi_grp)
    if (iunit <= 0) then
      ierr = ierr + 1
      message(1) = "Unable to open file '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1)
    else
      call iopar_read(mpi_grp, iunit, lines, 1, err)
      if (err /= 0) then
        ierr = ierr + 2
      else
        read(lines(1), '(a20,l)')  str, is_hypercube
      end if

      if (is_hypercube) then
        ! We have a hypercube: we will assume everything is OK...
        message(1) = "Simulation box is a hypercube: unable to check mesh compatibility."
        call messages_warning(1)
      else
        call iopar_read(mpi_grp, iunit, lines, 4, err)
        if (err /= 0) then
          ierr = ierr + 4
        else
          read(lines(1), '(a20,i21)')  str, read_np_part
          read(lines(2), '(a20,i21)')  str, read_np
          read(lines(3), '(a20,i21)')  str, algorithm
          read(lines(4), '(a20,i21)')  str, checksum

          ASSERT(read_np_part >= read_np)
            
          if (read_np_part == mesh%np_part_global &
               .and. read_np == mesh%np_global &
               .and. algorithm == 1 &
               .and. checksum == mesh%idx%checksum) then
            read_np_part = 0
            read_np = 0
          end if
        end if

      end if

      call io_close(iunit, grp=mpi_grp)
    end if

    POP_SUB(mesh_read_fingerprint)
  end subroutine mesh_read_fingerprint

  ! ---------------------------------------------------------
  subroutine mesh_check_dump_compatibility(mesh, dir, filename, namespace, mpi_grp, grid_changed, grid_reordered, map, ierr)
    type(mesh_t),         intent(in)  :: mesh
    character(len=*),     intent(in)  :: dir
    character(len=*),     intent(in)  :: filename
    type(namespace_t),    intent(in)  :: namespace
    type(mpi_grp_t),      intent(in)  :: mpi_grp
    logical,              intent(out) :: grid_changed
    logical,              intent(out) :: grid_reordered
    integer, allocatable, intent(out) :: map(:)
    integer,              intent(out) :: ierr

    integer :: ip, read_np_part, read_np, xx(MAX_DIM), err
    integer, allocatable :: read_lxyz(:,:)
    
    PUSH_SUB(mesh_check_dump_compatibility)

    ierr = 0

    grid_changed = .false.
    grid_reordered = .false.

    ! Read the mesh fingerprint
    call mesh_read_fingerprint(mesh, dir, filename, mpi_grp, namespace, read_np_part, read_np, err)
    if (err /= 0) then
      ierr = ierr + 1
      message(1) = "Unable to read mesh fingerprint from '"//trim(dir)//"/"//trim(filename)//"'."
      call messages_warning(1)

    else if (read_np > 0) then
      if (.not. associated(mesh%sb)) then
        ! We can only check the compatibility of two meshes that have different fingerprints if we also
        ! have the simulation box. In the case we do not, we will assume that the fingerprint is enough.
        ierr = ierr + 2
      else
        select type (box => mesh%sb%box)
        type is (box_hypercube_t)
          ! We cannot check the compatibility if the box is an hypercube
        class default
          grid_changed = .true.

          ! perhaps only the order of the points changed, this can only
          ! happen if the number of points is the same and no points maps
          ! to zero (this is checked below)
          grid_reordered = (read_np == mesh%np_global)

          ! the grid is different, so we read the coordinates.
          SAFE_ALLOCATE(read_lxyz(1:read_np_part, 1:mesh%sb%dim))
          ASSERT(allocated(mesh%idx%lxyz))
          call io_binary_read(trim(io_workpath(dir, namespace))//'/lxyz.obf', read_np_part*mesh%sb%dim, read_lxyz, err)
          if (err /= 0) then
            ierr = ierr + 4
            message(1) = "Unable to read index map from '"//trim(dir)//"'."
            call messages_warning(1)
          else
            ! generate the map
            SAFE_ALLOCATE(map(1:read_np))

            do ip = 1, read_np
              xx = 0
              xx(1:mesh%sb%dim) = read_lxyz(ip, 1:mesh%sb%dim)
              if (any(xx(1:mesh%sb%dim) < mesh%idx%nr(1, 1:mesh%sb%dim)) .or. &
                any(xx(1:mesh%sb%dim) > mesh%idx%nr(2, 1:mesh%sb%dim))) then
                map(ip) = 0
                grid_reordered = .false.
              else
                map(ip) = mesh_global_index_from_coords(mesh, [xx(1), xx(2), xx(3)])
                if(map(ip) > mesh%np_global) map(ip) = 0
              end if
            end do
          end if

          SAFE_DEALLOCATE_A(read_lxyz)
        end select
      end if
    end if

    POP_SUB(mesh_check_dump_compatibility)
  end subroutine mesh_check_dump_compatibility


  ! --------------------------------------------------------------
  recursive subroutine mesh_end(this)
    class(mesh_t), intent(inout)   :: this

    PUSH_SUB(mesh_end)

    call mesh_cube_map_end(this%cube_map)

    if(this%idx%is_hypercube) call hypercube_end(this%idx%hypercube)

    SAFE_DEALLOCATE_A(this%idx%lxyz)
    SAFE_DEALLOCATE_A(this%idx%lxyz_inv)
    SAFE_DEALLOCATE_A(this%x)
    SAFE_DEALLOCATE_A(this%vol_pp)

    if(this%parallel_in_domains) then
      call vec_end(this%vp)
      call partition_end(this%partition)
    end if

    POP_SUB(mesh_end)
  end subroutine mesh_end

  
  !> This function returns the point inside the grid corresponding to
  !! a boundary point when PBCs are used. In case the point does not
  !! have a correspondence (i.e. other BCs are used in that direction),
  !! the same point is returned. Note that this function returns a
  !! global point number when parallelization in domains is used.
  ! ---------------------------------------------------------  
  integer function mesh_periodic_point(mesh, space, ip) result(ipg)
    type(mesh_t),  intent(in)    :: mesh
    type(space_t), intent(in)    :: space
    integer,       intent(in)    :: ip     !< local point for which periodic copy is searched
    
    integer :: ix(MAX_DIM), nr(2, MAX_DIM), idim
    FLOAT :: xx(MAX_DIM), rr, ufn_re, ufn_im
    
    ! no push_sub, called too frequently

    call mesh_local_index_to_coords(mesh, ip, ix)
    nr(1, :) = mesh%idx%nr(1, :) + mesh%idx%enlarge(:)
    nr(2, :) = mesh%idx%nr(2, :) - mesh%idx%enlarge(:)
    
    do idim = 1, space%periodic_dim
      if(ix(idim) < nr(1, idim)) ix(idim) = ix(idim) + mesh%idx%ll(idim)
      if(ix(idim) > nr(2, idim)) ix(idim) = ix(idim) - mesh%idx%ll(idim)
    end do
    
    ipg = mesh_global_index_from_coords(mesh, ix)
    ASSERT(ipg > 0)

    if(mesh%masked_periodic_boundaries) then
      call mesh_r(mesh, ip, rr, coords = xx)
      call parse_expression(ufn_re, ufn_im, space%dim, xx, rr, M_ZERO, mesh%periodic_boundary_mask)
      if(int(ufn_re) == 0) ipg = mesh_local2global(mesh, ip) ! Nothing will be done
    end if 
    
  end function mesh_periodic_point
  

  ! ---------------------------------------------------------
  FLOAT pure function mesh_global_memory(mesh) result(memory)
    type(mesh_t), intent(in) :: mesh
    
    memory = M_ZERO
    
    ! lxyz_inv
    memory = memory + SIZEOF_UNSIGNED_INT * product(mesh%idx%nr(2, 1:mesh%sb%dim) - mesh%idx%nr(1, 1:mesh%sb%dim) + M_ONE)
    ! lxyz
    memory = memory + SIZEOF_UNSIGNED_INT * TOFLOAT(mesh%np_part_global) * MAX_DIM

  end function mesh_global_memory


  ! ---------------------------------------------------------
  FLOAT pure function mesh_local_memory(mesh) result(memory)
    type(mesh_t), intent(in) :: mesh
    
    memory = M_ZERO
    
    ! x
    memory = memory + REAL_PRECISION * TOFLOAT(mesh%np_part) * MAX_DIM
  end function mesh_local_memory


  ! ---------------------------------------------------------
  function mesh_x_global(mesh, ip, force) result(xx)
    type(mesh_t),       intent(in) :: mesh
    integer,            intent(in) :: ip
    logical, optional,  intent(in) :: force
    FLOAT                          :: xx(1:mesh%sb%dim)

    FLOAT :: chi(1:MAX_DIM)
    integer :: ix(1:MAX_DIM)
    logical :: force_

! no push_sub because function is called too frequently

    force_ = .false.
    if (present(force)) force_ = force
      
    if(mesh%parallel_in_domains .or. force_) then
      call mesh_global_index_to_coords(mesh, ip, ix)
      chi(1:mesh%sb%dim) = ix(1:mesh%sb%dim) * mesh%spacing(1:mesh%sb%dim)
      chi(mesh%sb%dim + 1:MAX_DIM) = M_ZERO
      xx = M_ZERO ! this initialization is required by gfortran 4.4 or we get NaNs
      call curvilinear_chi2x(mesh%sb, mesh%sb%latt, mesh%cv, chi, xx)
    else
      xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim)
    end if

  end function mesh_x_global


  ! ---------------------------------------------------------
  logical pure function mesh_compact_boundaries(mesh) result(cb)
    type(mesh_t),       intent(in) :: mesh
    
    cb = .not. mesh%use_curvilinear .and. &
         .not. mesh%parallel_in_domains

  end function mesh_compact_boundaries


  ! ---------------------------------------------------------
  subroutine mesh_check_symmetries(mesh, symm, periodic_dim)
    type(mesh_t),        intent(in) :: mesh
    type(symmetries_t),  intent(in) :: symm
    integer,             intent(in) :: periodic_dim

    integer :: iop, ip, idim, nops, ix(1:3)
    FLOAT :: destpoint(1:3), srcpoint(1:3), lsize(1:3), offset(1:3)

    !If all the axis have the same spacing and the same length
    !the grid is by obviously symmetric 
    !Indeed, reduced coordinates are proportional to the point index
    !and the reduced rotation are integer matrices
    !The result of the product is also proportional to an integer
    !and therefore belong to the grid.
    if(mesh%idx%ll(1) == mesh%idx%ll(2) .and.     &
        mesh%idx%ll(2) == mesh%idx%ll(3) .and.    &
         mesh%spacing(1) == mesh%spacing(2) .and. &
          mesh%spacing(2) == mesh%spacing(3) ) return 

    PUSH_SUB(mesh_check_symmetries)

    message(1) = "Checking if the real-space grid is symmetric";
    call messages_info(1)

    lsize(1:3) = TOFLOAT(mesh%idx%ll(1:3))
    offset(1:3) = TOFLOAT(mesh%idx%nr(1, 1:3) + mesh%idx%enlarge(1:3))

    nops = symmetries_number(symm)

    do ip = 1, mesh%np
      !We use floating point coordinates to check if the symmetric point 
      !belong to the grid.
      !If yes, it should have integer reduced coordinates 
      call mesh_local_index_to_coords(mesh, ip, ix)
      destpoint(1:3) = TOFLOAT(ix(1:3)) - offset(1:3)
      ! offset moves corner of cell to origin, in integer mesh coordinates
      ASSERT(all(destpoint >= 0))
      ASSERT(all(destpoint < lsize))

      ! move to center of cell in real coordinates
      destpoint = destpoint - TOFLOAT(int(lsize)/2)

      !convert to proper reduced coordinates
      do idim = 1, 3
        destpoint(idim) = destpoint(idim)/lsize(idim)
      end do

      ! iterate over all points that go to this point by a symmetry operation
      do iop = 1, nops
        srcpoint = symm_op_apply_red(symm%ops(iop), destpoint) 

        !We now come back to what should be an integer, if the symmetric point beloings to the grid
        do idim = 1, 3
          srcpoint(idim) = srcpoint(idim)*lsize(idim)
        end do

        ! move back to reference to origin at corner of cell
        srcpoint = srcpoint + TOFLOAT(int(lsize)/2)

        ! apply periodic boundary conditions in periodic directions 
        do idim = 1, periodic_dim
          if(nint(srcpoint(idim)) < 0 .or. nint(srcpoint(idim)) >= lsize(idim)) then
            srcpoint(idim) = modulo(srcpoint(idim)+M_HALF*SYMPREC, lsize(idim))
          end if
        end do
        ASSERT(all(srcpoint >= -SYMPREC))
        ASSERT(all(srcpoint < lsize))

        srcpoint(1:3) = srcpoint(1:3) + offset(1:3)
 
        if(any(srcpoint-anint(srcpoint)> SYMPREC*M_TWO)) then
          message(1) = "The real-space grid breaks at least one of the symmetries of the system."
          message(2) = "Change your spacing or use SymmetrizeDensity=no."
          call messages_fatal(2)
        end if
      end do
    end do

    POP_SUB(mesh_check_symmetries)
  end subroutine

  !> This function returns the true _global_ index of the point for a given
  !! vector of integer coordinates.
  integer function mesh_global_index_from_coords(mesh, ix) result(index)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ix(:)

    index = index_from_coords(mesh%idx, ix)
  end function mesh_global_index_from_coords

  !> Given a _global_ point index, this function returns the set of
  !! integer coordinates of the point.
  pure subroutine mesh_global_index_to_coords(mesh, ipg, ix)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ipg
    integer,       intent(out)   :: ix(:)

    call index_to_coords(mesh%idx, ipg, ix)
  end subroutine mesh_global_index_to_coords

  !> This function returns the _local_ index of the point for a given
  !! vector of integer coordinates.
  integer function mesh_local_index_from_coords(mesh, ix) result(ip)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ix(:)

    integer :: ipg

    ipg = index_from_coords(mesh%idx, ix)
    ip = mesh_global2local(mesh, ipg)
  end function mesh_local_index_from_coords

  !> Given a _local_ point index, this function returns the set of
  !! integer coordinates of the point.
  subroutine mesh_local_index_to_coords(mesh, ip, ix)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ip
    integer,       intent(out)   :: ix(:)

    integer :: ipg

    ipg = mesh_local2global(mesh, ip)
    call index_to_coords(mesh%idx, ipg, ix)
  end subroutine mesh_local_index_to_coords

  !> This function returns the global mesh index for a given local index
  integer function mesh_local2global(mesh, ip) result(ipg)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ip

    if (.not. mesh%parallel_in_domains) then
      ipg = ip
    else
      if (ip <= mesh%np) then
        ipg = mesh%vp%local(mesh%vp%xlocal + ip - 1)
      else if (ip <= mesh%np + mesh%vp%np_ghost) then
        ipg = mesh%vp%ghost(ip - mesh%np)
      else if (ip <= mesh%np + mesh%vp%np_ghost + mesh%vp%np_bndry) then
        ipg = mesh%vp%bndry(ip - mesh%np - mesh%vp%np_ghost)
      else
        ipg = 0
      end if
    end if
  end function mesh_local2global

  !> This function returns the local mesh index for a given global index
  integer function mesh_global2local(mesh, ipg) result(ip)
    type(mesh_t),  intent(in)    :: mesh
    integer,       intent(in)    :: ipg

    if (.not. mesh%parallel_in_domains) then
      ip = ipg
    else
      ip = vec_global2local(mesh%vp, ipg)
    end if
  end function mesh_global2local

#include "undef.F90"
#include "real.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mesh_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mesh_inc.F90"

end module mesh_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
