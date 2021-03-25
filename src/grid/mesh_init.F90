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

module mesh_init_oct_m
  use checksum_interface_oct_m
  use curvilinear_oct_m
  use global_oct_m
  use iihash_oct_m
  use index_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use mesh_partition_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use partition_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use sort_oct_m
  use space_oct_m
  use stencil_oct_m

  implicit none
  
  private
  public ::                    &
    mesh_init_stage_1,         &
    mesh_init_stage_2,         &
    mesh_init_stage_3

  type(profile_t), save :: mesh_init_prof

  integer, parameter :: INNER_POINT = 1
  integer, parameter :: ENLARGEMENT_POINT = 2
  integer, parameter :: BOUNDARY = -1
  
contains

! ---------------------------------------------------------
subroutine mesh_init_stage_1(mesh, namespace, space, sb, cv, spacing, enlarge)
  type(mesh_t),                intent(inout) :: mesh
  type(namespace_t),           intent(in)    :: namespace
  type(space_t),               intent(in)    :: space
  type(simul_box_t),   target, intent(in)    :: sb
  type(curvilinear_t), target, intent(in)    :: cv
  FLOAT,                       intent(in)    :: spacing(1:MAX_DIM)
  integer,                     intent(in)    :: enlarge(MAX_DIM)

  integer :: idir, jj, delta
  FLOAT   :: x(MAX_DIM), chi(MAX_DIM), spacing_new(-1:1)
  logical :: out
  FLOAT, parameter :: DELTA_ = CNST(1e-12)

  PUSH_SUB(mesh_init_stage_1)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  mesh%sb => sb     ! keep an internal pointer
  mesh%spacing = spacing ! this number can change in the following
  mesh%use_curvilinear = cv%method /= CURV_METHOD_UNIFORM
  mesh%cv => cv

  mesh%idx%dim = space%dim
  mesh%idx%enlarge = enlarge

  ! adjust nr
  mesh%idx%nr = 0
  do idir = 1, space%dim
    chi = M_ZERO
    ! the upper border
    jj = 0
    out = .false.
    do while(.not.out)
      jj = jj + 1
      chi(idir) = TOFLOAT(jj)*mesh%spacing(idir)
      if ( mesh%use_curvilinear ) then
        call curvilinear_chi2x(sb, cv, chi(1:space%dim), x(1:space%dim))
        out = x(idir) > sb%lsize(idir) + DELTA_
      else
        ! do the same comparison here as in simul_box_contains_points
        out = chi(idir) > sb%lsize(idir) + DELTA_
      end if
    end do
    mesh%idx%nr(2, idir) = jj - 1
  end do

  ! we have a symmetric mesh (for now)
  mesh%idx%nr(1,:) = -mesh%idx%nr(2,:)

  ! we have to adjust a couple of things for the periodic directions
  do idir = 1, space%periodic_dim
    if(mesh%idx%nr(2, idir) == 0) then
      ! this happens if Spacing > box size
      mesh%idx%nr(2, idir) =  1
      mesh%idx%nr(1, idir) = -1
    end if

    ! We have to adjust the spacing to be commensurate with the box,
    ! for this we scan the possible values of the grid size around the
    ! one we selected. We choose the size that has the spacing closest
    ! to the requested one.
    do delta = -1, 1
      spacing_new(delta) = CNST(2.0)*sb%lsize(idir)/TOFLOAT(2*mesh%idx%nr(2, idir) + 1 - delta)
      spacing_new(delta) = abs(spacing_new(delta) - spacing(idir))
    end do

    delta = minloc(spacing_new, dim = 1) - 2

    ASSERT(delta >= -1) 
    ASSERT(delta <=  1) 

    mesh%spacing(idir) = CNST(2.0)*sb%lsize(idir)/TOFLOAT(2*mesh%idx%nr(2, idir) + 1 - delta)

    ! we need to adjust the grid by adding or removing one point
    if(delta == -1) then
      mesh%idx%nr(1, idir) = mesh%idx%nr(1, idir) - 1
    else if (delta == 1) then
      mesh%idx%nr(2, idir) = mesh%idx%nr(2, idir) - 1
    end if
    
  end do

  if( any(abs(mesh%spacing(1:space%periodic_dim) - spacing(1:space%periodic_dim)) > CNST(1e-6)) ) then
    call messages_write('The spacing has been modified to make it commensurate with the periodicity of the system.')
    call messages_warning()
  end if

  do idir = space%periodic_dim + 1, space%dim
    if(mesh%idx%nr(2, idir) == 0) then
      write(message(1),'(a,i2)') 'Spacing > box size in direction ', idir
      call messages_fatal(1, namespace=namespace)
    end if
  end do

  mesh%idx%ll(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1


  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_1)
end subroutine mesh_init_stage_1

! ---------------------------------------------------------
!> This subroutine creates the global array of Hilbert indices
!! and the inverse mapping.
subroutine mesh_init_stage_2(mesh, namespace, space, sb, cv, stencil)
  type(mesh_t),        intent(inout) :: mesh
  type(namespace_t),   intent(in)    :: namespace
  type(space_t),       intent(in)    :: space
  type(simul_box_t),   intent(in)    :: sb
  type(curvilinear_t), intent(in)    :: cv
  type(stencil_t),     intent(in)    :: stencil

  integer :: is
  FLOAT   :: chi(MAX_DIM)
  integer :: point(1:MAX_DIM), point_stencil(1:MAX_DIM)
  integer(8) :: global_size, local_size, sizes(1:MAX_DIM)
  integer(8) :: ihilbert, ihilbertb, istart, iend, hilbert_size
  integer :: ip, ip2, ib, ib2, np, np_boundary
  FLOAT :: pos(1:MAX_DIM)
  logical :: found
  integer :: irank, nduplicate
  type(lihash_t) :: hilbert_to_grid, hilbert_to_boundary
  integer(8), allocatable :: boundary_to_hilbert(:), boundary_to_hilbert_global(:)
  integer(8), allocatable :: grid_to_hilbert(:), grid_to_hilbert_initial(:)
  integer, allocatable :: scounts(:), sdispls(:), rcounts(:), rdispls(:)
  integer, allocatable :: initial_sizes(:), initial_offsets(:), final_sizes(:), offsets(:)

  PUSH_SUB(mesh_init_stage_2)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  ! enlarge mesh for boundary points
  mesh%idx%nr(1, 1:MAX_DIM) = mesh%idx%nr(1, 1:MAX_DIM) - mesh%idx%enlarge(1:MAX_DIM)
  mesh%idx%nr(2, 1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) + mesh%idx%enlarge(1:MAX_DIM)
  
  sizes(1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) - mesh%idx%nr(1, 1:MAX_DIM) + 1
  mesh%idx%offset(1:MAX_DIM) = sizes(1:MAX_DIM)/2
  if(space%dim > 1 .and. any(sizes > 2**(63/int(space%dim,8)))) then
    write(message(1), '(A, I10, A, I2, A)') "Error: grid too large, more than ", 2**(63/int(space%dim,8)), &
      " points in one direction for ", space%dim, " dimensions. This is not supported."
    call messages_fatal(1)
  end if
  global_size = product(sizes)
  ! compute the bits per dimension: sizes(i) <= 2**bits
  mesh%idx%bits = maxval(ceiling(log(TOFLOAT(sizes))/log(2.)))

  hilbert_size = 2**(space%dim*mesh%idx%bits)

  ! use block data decomposition of hilbert indices
  istart = floor(TOFLOAT(hilbert_size) * mpi_world%rank/mpi_world%size)
  iend = floor(TOFLOAT(hilbert_size) * (mpi_world%rank+1)/mpi_world%size) - 1
  local_size = iend - istart + 1

  call lihash_init(hilbert_to_grid)
  SAFE_ALLOCATE(grid_to_hilbert_initial(1:local_size))

  ! get grid indices
  ip = 1
  do ihilbert = istart, iend
    call index_hilbert_to_point(mesh%idx, space%dim, ihilbert, point)
    ! first check if point is outside bounding box
    if(any(point(1:space%dim) < mesh%idx%nr(1, 1:space%dim) + mesh%idx%enlarge(1:space%dim))) cycle
    if(any(point(1:space%dim) > mesh%idx%nr(2, 1:space%dim) - mesh%idx%enlarge(1:space%dim))) cycle
    ! then check if point is inside simulation box
    chi(1:space%dim) = TOFLOAT(point(1:space%dim)) * mesh%spacing(1:space%dim)
    call curvilinear_chi2x(sb, cv, chi(:), pos(:))
    if(.not. sb%contains_point(pos)) cycle
    grid_to_hilbert_initial(ip) = ihilbert
    call lihash_insert(hilbert_to_grid, ihilbert, ip)
    ip = ip + 1
  end do
  np = ip - 1

  SAFE_ALLOCATE(initial_sizes(0:mpi_world%size-1))
#ifdef HAVE_MPI
  call MPI_Allgather(np, 1, MPI_INTEGER, initial_sizes(0), 1, MPI_INTEGER, mpi_world%comm, mpi_err)
#else
  initial_sizes(0) = np
#endif
  SAFE_ALLOCATE(initial_offsets(0:mpi_world%size))
  initial_offsets(0) = 0
  do irank = 1, mpi_world%size
    initial_offsets(irank) = initial_offsets(irank-1) + initial_sizes(irank-1)
  end do

  ! now distribute ordered by Hilbert index
  ! use block data decomposition of grid indices
  SAFE_ALLOCATE(offsets(0:mpi_world%size))
  SAFE_ALLOCATE(final_sizes(0:mpi_world%size-1))

  do irank = 0, mpi_world%size
    offsets(irank) = floor(TOFLOAT(sum(initial_sizes)) * irank/mpi_world%size)
  end do
  do irank = 0, mpi_world%size - 1
    final_sizes(irank) = offsets(irank + 1) - offsets(irank)
  end do

  SAFE_ALLOCATE(scounts(0:mpi_world%size-1))
  SAFE_ALLOCATE(sdispls(0:mpi_world%size-1))
  SAFE_ALLOCATE(rcounts(0:mpi_world%size-1))
  SAFE_ALLOCATE(rdispls(0:mpi_world%size-1))
  ! determine communication pattern
  scounts = 0
  do irank = 0, mpi_world%size - 1
    ! get overlap of initial and final distribution
    scounts(irank) = max(0, &
      min(offsets(irank+1), initial_offsets(mpi_world%rank+1)) - &
      max(offsets(irank), initial_offsets(mpi_world%rank)))
  end do
  sdispls(0) = 0
  do irank = 1, mpi_world%size - 1
    sdispls(irank) = sdispls(irank - 1) + scounts(irank - 1)
  end do
  ASSERT(sum(scounts) == initial_sizes(mpi_world%rank))

  rcounts = 0
  do irank = 0, mpi_world%size - 1
    ! get overlap of initial and final distribution
    rcounts(irank) = max(0, &
      min(offsets(mpi_world%rank+1), initial_offsets(irank+1)) - &
      max(offsets(mpi_world%rank), initial_offsets(irank)))
  end do
  rdispls(0) = 0
  do irank = 1, mpi_world%size - 1
    rdispls(irank) = rdispls(irank - 1) + rcounts(irank - 1)
  end do
  ASSERT(sum(rcounts) == final_sizes(mpi_world%rank))

  SAFE_ALLOCATE(grid_to_hilbert(1:final_sizes(mpi_world%rank)))
#ifdef HAVE_MPI
  call MPI_Alltoallv(grid_to_hilbert_initial(1), scounts(0), sdispls(0), MPI_LONG_LONG, &
                     grid_to_hilbert(1), rcounts(0), rdispls(0), MPI_LONG_LONG, &
                     mpi_world%comm, mpi_err)
#else
  grid_to_hilbert(1:np) = grid_to_hilbert_initial(1:np)
#endif

  SAFE_DEALLOCATE_A(grid_to_hilbert_initial)

  SAFE_DEALLOCATE_A(scounts)
  SAFE_DEALLOCATE_A(sdispls)
  SAFE_DEALLOCATE_A(rcounts)
  SAFE_DEALLOCATE_A(rdispls)


  mesh%np = final_sizes(mpi_world%rank)
  mesh%np_global = sum(final_sizes)

  ! get local boundary indices
  call lihash_init(hilbert_to_boundary)
  SAFE_ALLOCATE(boundary_to_hilbert(1:mesh%np*(stencil%size - 1)))
  ib = 1
  do ip = 1, mesh%np
    call index_hilbert_to_point(mesh%idx, space%dim, grid_to_hilbert(ip), point)
    do is = 1, stencil%size
      if(stencil%center == is) cycle
      point_stencil(1:space%dim) = point(1:space%dim) + stencil%points(1:space%dim, is)
      ! check if point is in inner part
      call index_point_to_hilbert(mesh%idx, space%dim, ihilbertb, point_stencil)
      ib2 = lihash_lookup(hilbert_to_grid, ihilbertb, found)
      if(found) cycle
      ! then check if point is inside simulation box
      chi(1:space%dim) = TOFLOAT(point_stencil(1:space%dim)) * mesh%spacing(1:space%dim)
      call curvilinear_chi2x(sb, cv, chi(:), pos(:))
      if(sb%contains_point(pos)) cycle
      ! it has to be a boundary point now
      ! check if already counted
      ib2 = lihash_lookup(hilbert_to_boundary, ihilbertb, found)
      if(found) cycle
      boundary_to_hilbert(ib) = ihilbertb
      call lihash_insert(hilbert_to_boundary, ihilbertb, ib)
      ib = ib + 1
    end do
  end do
  call lihash_end(hilbert_to_boundary)
  np_boundary = ib - 1
  mesh%np_part = mesh%np + np_boundary

#ifdef HAVE_MPI
  ! get global boundary indices
  call MPI_Allgather(np_boundary, 1, MPI_INTEGER, initial_sizes(0), 1, MPI_INTEGER, mpi_world%comm, mpi_err)
#else
  initial_sizes(0) = np_boundary
#endif
  initial_offsets(0) = 0
  do irank = 1, mpi_world%size
    initial_offsets(irank) = initial_offsets(irank-1) + initial_sizes(irank-1)
  end do

  ! boundary
  SAFE_ALLOCATE(boundary_to_hilbert_global(1:sum(initial_sizes)))
#ifdef HAVE_MPI
  call MPI_Allgatherv(boundary_to_hilbert(1), np_boundary, MPI_LONG_LONG, &
    boundary_to_hilbert_global(1), initial_sizes(0), initial_offsets(0), MPI_LONG_LONG, &
    mpi_world%comm, mpi_err)
#else
  boundary_to_hilbert_global(1:sum(initial_sizes)) = boundary_to_hilbert(1:sum(initial_sizes))
#endif

  ! sort boundary
  call sort(boundary_to_hilbert_global)

  ! get number of non-unique elements (could be done in parallel)
  nduplicate = 0
  do ip = 1, sum(initial_sizes) - 1
    if (boundary_to_hilbert_global(ip+1) == boundary_to_hilbert_global(ip)) then
      nduplicate = nduplicate + 1
    end if
  end do
  mesh%np_part_global = mesh%np_global + sum(initial_sizes) - nduplicate

  ! get global indices
  ! inner grid
  SAFE_ALLOCATE(mesh%idx%grid_to_hilbert_global(1:mesh%np_part_global))
#ifdef HAVE_MPI
  call MPI_Allgatherv(grid_to_hilbert(1), mesh%np, MPI_LONG_LONG, &
    mesh%idx%grid_to_hilbert_global(1), final_sizes(0), offsets(0), MPI_LONG_LONG, &
    mpi_world%comm, mpi_err)
#else
  mesh%idx%grid_to_hilbert_global(1:mesh%np) = grid_to_hilbert(1:mesh%np)
#endif

  ! add unique boundary indices
  ip2 = mesh%np_global + 1
  mesh%idx%grid_to_hilbert_global(ip2) = boundary_to_hilbert_global(1)
  ip2 = ip2 + 1
  do ip = 2, sum(initial_sizes)
    if (boundary_to_hilbert_global(ip) /= boundary_to_hilbert_global(ip-1)) then
      mesh%idx%grid_to_hilbert_global(ip2) = boundary_to_hilbert_global(ip)
      ip2 = ip2 + 1
    end if
  end do

  SAFE_DEALLOCATE_A(initial_offsets)
  SAFE_DEALLOCATE_A(initial_sizes)

  ! fill global hash map
  call lihash_init(mesh%idx%hilbert_to_grid_global)
  do ip = 1, mesh%np_part_global
    call lihash_insert(mesh%idx%hilbert_to_grid_global, &
      mesh%idx%grid_to_hilbert_global(ip), ip)
  end do

  SAFE_DEALLOCATE_A(offsets)
  SAFE_DEALLOCATE_A(final_sizes)

  SAFE_DEALLOCATE_A(boundary_to_hilbert)
  SAFE_DEALLOCATE_A(boundary_to_hilbert_global)
  call lihash_end(hilbert_to_grid)
  SAFE_DEALLOCATE_A(grid_to_hilbert)

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_2)
end subroutine mesh_init_stage_2

! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------
subroutine mesh_init_stage_3(mesh, namespace, space, stencil, mc, parent)
  type(mesh_t),              intent(inout) :: mesh
  type(namespace_t),         intent(in)    :: namespace
  type(space_t),             intent(in)    :: space
  type(stencil_t),           intent(in)    :: stencil
  type(multicomm_t),         intent(in)    :: mc
  type(mesh_t),    optional, intent(in)    :: parent

  integer :: ip

  PUSH_SUB(mesh_init_stage_3)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  call mpi_grp_init(mesh%mpi_grp, mc%group_comm(P_STRATEGY_DOMAINS))
  
  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = (mesh%mpi_grp%size > 1)

  ! reorder global index
  call reorder_points()
  call checksum_calculate(1, int(mesh%np_part_global, 8), mesh%idx%grid_to_hilbert_global(1), mesh%idx%checksum)

  if(mesh%parallel_in_domains) then
    call do_partition()
  else
    ! When running serially those two are the same.
    mesh%np      = mesh%np_global
    mesh%np_part = mesh%np_part_global

    ! These must be initialized for vec_gather, vec_scatter to work
    ! as copy operations when running without domain parallelization.
    mesh%vp%np_global = mesh%np_global
    mesh%vp%np_ghost = 0
    mesh%vp%np_bndry = mesh%np_part - mesh%np
    mesh%vp%npart = 1
    mesh%vp%xlocal = 1
  end if

  ! Compute mesh%x
  SAFE_ALLOCATE(mesh%x(1:mesh%np_part, 1:space%dim))
  mesh%x(:, :) = M_ZERO
  do ip = 1, mesh%np_part
    mesh%x(ip, 1:space%dim) = mesh_x_global(mesh, mesh_local2global(mesh, ip))
  end do

  call mesh_cube_map_init(mesh%cube_map, mesh%idx, mesh%np_global)

  call mesh_get_vol_pp(mesh%sb)

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_3)

contains
  subroutine reorder_points()
    integer, allocatable :: initial_sizes(:), initial_offsets(:)
    integer(8), allocatable :: reordered(:)
    integer :: ix, iy, iz, ixb, iyb, izb, ip_inner, ip_boundary, ipg, nn, idir
    integer :: boundary_start, istart, iend, local_size, irank
    integer :: bsize(3), order, default
    type(block_t) :: blk
    integer, parameter :: &
      ORDER_BLOCKS     =  1, &
      ORDER_HILBERT    =  2, &
      ORDER_CUBE       =  3

    PUSH_SUB(mesh_init_stage_3.reorder_points)

    !%Variable MeshOrder
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% This variable controls how the grid points are mapped to a
    !% linear array for global arrays. For runs that are parallel
    !% in domains, the local mesh order may be different (see
    !% <tt>MeshLocalOrder</tt>).
    !% The default is blocks when serial in domains and cube when
    !% parallel in domains with the local mesh order set to blocks.
    !%Option blocks 1
    !% The grid is mapped using small parallelepipedic grids. The size
    !% of the blocks is controlled by <tt>MeshBlockSize</tt>.
    !%Option hilbert 2
    !% A Hilbert space-filling curve is used to map the grid.
    !%Option order_cube 3
    !% The grid is mapped using a full cube, i.e. without blocking.
    !%End
    if (.not. mesh%parallel_in_domains) then
      default = ORDER_BLOCKS
    else
      default = ORDER_CUBE
    end if
    call parse_variable(namespace, 'MeshOrder', default, order)

    select case(order)
    case(ORDER_HILBERT)
      ! nothing to do, points are ordered along a Hilbert curve by default
    case(ORDER_BLOCKS, ORDER_CUBE)
      if (order == ORDER_CUBE) then
        bsize = mesh%idx%ll
      else
        !%Variable MeshBlockSize
        !%Type block
        !%Section Execution::Optimization
        !%Description
        !% To improve memory-access locality when calculating derivatives,
        !% <tt>Octopus</tt> arranges mesh points in blocks. This variable
        !% controls the size of this blocks in the different
        !% directions. The default is selected according to the value of
        !% the <tt>StatesBlockSize</tt> variable. (This variable only affects the
        !% performance of <tt>Octopus</tt> and not the results.)
        !%End

        ! blocking done according to layer conditions in 3d for an L2 cache
        ! of 1024 KB for the default 25pt stencil
        ! See J. Hammer, et al. in Tools for High Performance Computing 2016,
        ! ISBN 978-3-319-56702-0, 1-22 (2017). Proceedings of IPTW 2016,
        ! DOI: 10.1007/978-3-319-56702-0_1
        ! tool: https://rrze-hpc.github.io/layer-condition
        if (conf%target_states_block_size < 32) then
          bsize = (/ 56, 56, 1 /) / abs(conf%target_states_block_size)
        else
          bsize = (/ 10, 2, 1 /)
        end if

        ! no blocking in z direction
        bsize(3) = mesh%idx%ll(3)

        if(parse_block(namespace, 'MeshBlockSize', blk) == 0) then
          nn = parse_block_cols(blk, 0)
          do idir = 1, nn
            call parse_block_integer(blk, 0, idir - 1, bsize(idir))
          end do
        end if
      end if

      ! do the global reordering in parallel
      ! use block data decomposition of global indices
      SAFE_ALLOCATE(initial_offsets(0:mpi_world%size))
      SAFE_ALLOCATE(initial_sizes(0:mpi_world%size-1))
      do irank = 0, mpi_world%size
        initial_offsets(irank) = floor(TOFLOAT(mesh%np_part_global) * irank/mpi_world%size)
      end do
      do irank = 0, mpi_world%size - 1
        initial_sizes(irank) = initial_offsets(irank+1) - initial_offsets(irank)
      end do
      istart = initial_offsets(mpi_world%rank) + 1
      iend = initial_offsets(mpi_world%rank + 1)
      local_size = iend - istart + 1
      boundary_start = max(mesh%np_global-istart+1, 0)
      ASSERT(local_size == initial_sizes(mpi_world%rank))

      ! now comes the reordering loop
      SAFE_ALLOCATE(reordered(1:local_size))
      reordered = 0
      ip_inner = 1
      ip_boundary = 1
      do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsize(3)
        do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize(2)
          do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize(1)

            do iz = izb, min(izb + bsize(3) - 1, mesh%idx%nr(2,3))
              do iy = iyb, min(iyb + bsize(2) - 1, mesh%idx%nr(2,2))
                do ix = ixb, min(ixb + bsize(1) - 1, mesh%idx%nr(2,1))
                  ipg = index_from_coords(mesh%idx, [ix, iy, iz])
                  if (ipg < istart .or. ipg > iend) cycle
                  if (ipg <= mesh%np_global) then
                    reordered(ip_inner) = mesh%idx%grid_to_hilbert_global(ipg)
                    ip_inner = ip_inner + 1
                  else
                    reordered(boundary_start+ip_boundary) = mesh%idx%grid_to_hilbert_global(ipg)
                    ip_boundary = ip_boundary + 1
                  end if
                end do
              end do
            end do

          end do
        end do
      end do
      ASSERT(ip_inner + ip_boundary - 2 == local_size)
      ASSERT(all(reordered > 0))
      ! gather the reordered index
#ifdef HAVE_MPI
      call MPI_Allgatherv(reordered(1), local_size, MPI_LONG_LONG, &
        mesh%idx%grid_to_hilbert_global(1), initial_sizes(0), initial_offsets(0), MPI_LONG_LONG, &
        mpi_world%comm, mpi_err)
#else
      do ipg = 1, mesh%np_part_global
        mesh%idx%grid_to_hilbert_global(ipg) = reordered(ipg)
      end do
#endif
      SAFE_DEALLOCATE_A(reordered)
      SAFE_DEALLOCATE_A(initial_offsets)
      SAFE_DEALLOCATE_A(initial_sizes)

      ! Recreate hash table.
      call lihash_end(mesh%idx%hilbert_to_grid_global)
      call lihash_init(mesh%idx%hilbert_to_grid_global)
      ! Insert points
      do ipg = 1, mesh%np_part_global
        call lihash_insert(mesh%idx%hilbert_to_grid_global, &
          mesh%idx%grid_to_hilbert_global(ipg), ipg)
      end do
    end select

    POP_SUB(mesh_init_stage_3.reorder_points)
  end subroutine reorder_points

  ! ---------------------------------------------------------
  subroutine do_partition()
#ifdef HAVE_MPI
    integer :: jj, ipart, jpart, ip
    integer, allocatable :: gindex(:), gedges(:)
    logical, allocatable :: nb(:, :)
    integer              :: idx(1:MAX_DIM), jx(1:MAX_DIM)
    integer              :: graph_comm, iedge, nnb
    logical              :: use_topo, reorder, partition_print
    integer              :: ierr

    logical :: has_virtual_partition = .false.
    integer :: vsize !< 'virtual' partition size
    type(restart_t) :: restart_load, restart_dump
    integer, allocatable :: part_vec(:)
    
    PUSH_SUB(mesh_init_stage_3.do_partition)

    !Try to load the partition from the restart files
    call restart_init(restart_load, namespace, RESTART_PARTITION, RESTART_TYPE_LOAD, mc, ierr, mesh=mesh, exact=.true.)
    if (ierr == 0) call mesh_partition_load(restart_load, mesh, ierr)
    call restart_end(restart_load)

    if(ierr /= 0) then

      !%Variable MeshPartitionVirtualSize
      !%Type integer
      !%Default mesh mpi_grp size
      !%Section Execution::Parallelization
      !%Description
      !% Gives the possibility to change the partition nodes.
      !% Afterward, it crashes.
      !%End
      call parse_variable(namespace, 'MeshPartitionVirtualSize', mesh%mpi_grp%size, vsize)
      
      if (vsize /= mesh%mpi_grp%size) then
        write(message(1),'(a,I7)') "Changing the partition size to", vsize
        write(message(2),'(a)') "The execution will crash."
        call messages_warning(2)
        has_virtual_partition = .true.
      else 
        has_virtual_partition = .false.
      end if
      
      if(.not. present(parent)) then
        call mesh_partition(mesh, namespace, stencil, vsize)
      else
        ! if there is a parent grid, use its partition
        call mesh_partition_from_parent(mesh, parent)
      end if
      
      !Now that we have the partitions, we save them
      call restart_init(restart_dump, namespace, RESTART_PARTITION, RESTART_TYPE_DUMP, mc, ierr, mesh=mesh)
      call mesh_partition_dump(restart_dump, mesh, vsize, ierr)
      call restart_end(restart_dump)
    end if

    if (has_virtual_partition) then
      call profiling_end(namespace)
      call print_date("Calculation ended on ")
      write(message(1),'(a)') "Execution has ended."
      write(message(2),'(a)') "If you want to run your system, do not use MeshPartitionVirtualSize."
      call messages_warning(2)
      call messages_end()
      call global_end()
      stop
    end if

    !%Variable MeshUseTopology
    !%Type logical
    !%Default false
    !%Section Execution::Parallelization
    !%Description
    !% (experimental) If enabled, <tt>Octopus</tt> will use an MPI virtual
    !% topology to map the processors. This can improve performance
    !% for certain interconnection systems.
    !%End
    call parse_variable(namespace, 'MeshUseTopology', .false., use_topo)

    if(use_topo) then
      ! At the moment we still need the global partition. This will be removed in near future.
      SAFE_ALLOCATE(part_vec(1:mesh%np_part_global))
      call partition_get_global(mesh%partition, part_vec(1:mesh%np_global))


      ! generate a table of neighbours

      SAFE_ALLOCATE(nb(1:mesh%mpi_grp%size, 1:mesh%mpi_grp%size))
      nb = .false.

      do ip = 1, mesh%np_global
        ipart = part_vec(ip)
        call mesh_global_index_to_coords(mesh, ip, idx)
        do jj = 1, stencil%size
          jx(1:MAX_DIM) = idx(1:MAX_DIM) + stencil%points(1:MAX_DIM, jj)
          if(all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
            jpart = part_vec(mesh_global_index_from_coords(mesh, jx))
            if(ipart /= jpart ) nb(ipart, jpart) = .true.
          end if
        end do
      end do
      SAFE_DEALLOCATE_A(part_vec)

      ! now generate the information of the graph 

      SAFE_ALLOCATE(gindex(1:mesh%mpi_grp%size))
      SAFE_ALLOCATE(gedges(1:count(nb)))
      
     ! and now generate it
      iedge = 0
      do ipart = 1, mesh%mpi_grp%size
        do jpart = 1, mesh%mpi_grp%size
          if(nb(ipart, jpart)) then
            iedge = iedge + 1
            gedges(iedge) = jpart - 1
          end if
        end do
        gindex(ipart) = iedge
      end do

      ASSERT(iedge == count(nb))

      reorder = .true.
      call MPI_Graph_create(mesh%mpi_grp%comm, mesh%mpi_grp%size, gindex, gedges, reorder, graph_comm, mpi_err)

      ! we have a new communicator
      call mpi_grp_init(mesh%mpi_grp, graph_comm)

      SAFE_DEALLOCATE_A(nb)
      SAFE_DEALLOCATE_A(gindex)
      SAFE_DEALLOCATE_A(gedges)

    end if

    call vec_init(mesh%mpi_grp%comm, mesh%np_global, mesh%np_part_global, mesh%idx, stencil,&
         space, mesh%partition, mesh%vp, namespace)

    ! check the number of ghost neighbours in parallel
    nnb = 0
    jpart =  mesh%vp%partno
    do ipart = 1, mesh%vp%npart
      if (ipart == jpart) cycle
      if (mesh%vp%ghost_scounts(ipart) /= 0) nnb = nnb + 1
    end do
    ASSERT(nnb >= 0 .and. nnb < mesh%vp%npart)

    ! Set local point numbers.
    mesh%np      = mesh%vp%np_local
    mesh%np_part = mesh%np + mesh%vp%np_ghost + mesh%vp%np_bndry

    !%Variable PartitionPrint
    !%Type logical
    !%Default true
    !%Section Execution::Parallelization
    !%Description
    !% (experimental) If disabled, <tt>Octopus</tt> will not compute
    !% nor print the partition information, such as local points,
    !% no. of neighbours, ghost points and boundary points.
    !%End
    call parse_variable(namespace, 'PartitionPrint', .true., partition_print)
    
    if (partition_print) then
      call mesh_partition_write_info(mesh)
      call mesh_partition_messages_debug(mesh, namespace)
    end if   
#endif

    POP_SUB(mesh_init_stage_3.do_partition)
  end subroutine do_partition


  ! ---------------------------------------------------------
  !> calculate the volume of integration
  subroutine mesh_get_vol_pp(sb)
    type(simul_box_t), intent(in) :: sb

    integer :: jj(1:MAX_DIM), ip, np
    FLOAT   :: chi(MAX_DIM)

    PUSH_SUB(mesh_init_stage_3.mesh_get_vol_pp)

    np = 1
    if(mesh%use_curvilinear) np = mesh%np_part

    SAFE_ALLOCATE(mesh%vol_pp(1:np))

    do ip = 1, np
      mesh%vol_pp(ip) = product(mesh%spacing(1:space%dim))
    end do

    jj(space%dim + 1:MAX_DIM) = 0

    do ip = 1, np
      call mesh_local_index_to_coords(mesh, ip, jj)
      chi(1:space%dim) = jj(1:space%dim)*mesh%spacing(1:space%dim)
      mesh%vol_pp(ip) = mesh%vol_pp(ip)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip, 1:space%dim), chi(1:space%dim))
    end do

    if(mesh%use_curvilinear) then
      mesh%volume_element = CNST(1.0)
    else
      mesh%volume_element = mesh%vol_pp(1)
    end if

    if (space%dim == 3) then
      mesh%surface_element(1) = sqrt(abs(sum(dcross_product(sb%latt%rlattice_primitive(1:3, 2), &
                                                            sb%latt%rlattice_primitive(1:3, 3))**2)))
      mesh%surface_element(2) = sqrt(abs(sum(dcross_product(sb%latt%rlattice_primitive(1:3, 3), &
                                                            sb%latt%rlattice_primitive(1:3, 1))**2)))
      mesh%surface_element(3) = sqrt(abs(sum(dcross_product(sb%latt%rlattice_primitive(1:3, 1), &
                                                            sb%latt%rlattice_primitive(1:3, 2))**2)))
    else
      mesh%surface_element(1:space%dim) = M_ZERO
    end if

    POP_SUB(mesh_init_stage_3.mesh_get_vol_pp)
  end subroutine mesh_get_vol_pp
end subroutine mesh_init_stage_3

end module mesh_init_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
