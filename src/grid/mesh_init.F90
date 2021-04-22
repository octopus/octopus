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
  use box_hypercube_oct_m
  use checksum_interface_oct_m
  use curvilinear_oct_m
  use global_oct_m
  use hypercube_oct_m
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
  select type (box => sb%box)
  type is (box_hypercube_t)
    mesh%idx%is_hypercube = .true.
  class default
    mesh%idx%is_hypercube = .false.
  end select
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
        call curvilinear_chi2x(sb, sb%latt, cv, chi(1:space%dim), x(1:space%dim))
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
!> This subroutine checks if every grid point belongs to the internal
!! mesh, based on the global lxyz_inv matrix. Afterwards, it counts
!! how many points has the mesh and the enlargement.
subroutine mesh_init_stage_2(mesh, space, sb, cv, stencil)
  type(mesh_t),        intent(inout) :: mesh
  type(space_t),       intent(in)    :: space
  type(simul_box_t),   intent(in)    :: sb
  type(curvilinear_t), intent(in)    :: cv
  type(stencil_t),     intent(in)    :: stencil

  integer :: il, ik, ix, iy, iz, is
  integer :: ii, jj, kk
  FLOAT   :: chi(MAX_DIM)
  integer :: nr(1:2, 1:MAX_DIM), res
  logical, allocatable :: in_box(:)
  FLOAT,   allocatable :: xx(:, :)
  FLOAT  , parameter :: DELTA = CNST(1e-12)
  integer :: start_z, end_z
  type(profile_t), save :: prof
#if defined(HAVE_MPI) && defined(HAVE_MPI2)
  type(profile_t), save :: prof_reduce
  integer :: npoints
  integer, allocatable :: start(:), end(:)
#endif

  PUSH_SUB(mesh_init_stage_2)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  ! enlarge mesh for boundary points
  mesh%idx%nr(1, 1:MAX_DIM) = mesh%idx%nr(1, 1:MAX_DIM) - mesh%idx%enlarge(1:MAX_DIM)
  mesh%idx%nr(2, 1:MAX_DIM) = mesh%idx%nr(2, 1:MAX_DIM) + mesh%idx%enlarge(1:MAX_DIM)
  
  if(mesh%idx%is_hypercube) then
    call hypercube_init(mesh%idx%hypercube, space%dim, mesh%idx%nr, mesh%idx%enlarge(1))
    mesh%np_part_global = hypercube_number_total_points(mesh%idx%hypercube)
    mesh%np_global      = hypercube_number_inner_points(mesh%idx%hypercube)
    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_2)
    return
  end if

  nr = mesh%idx%nr

  ! allocate the xyz arrays
  SAFE_ALLOCATE(mesh%idx%lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))

  mesh%idx%lxyz_inv(:,:,:) = 0
  res = 1
  SAFE_ALLOCATE(xx(mesh%idx%nr(1,1):mesh%idx%nr(2,1), 1:MAX_DIM))
  SAFE_ALLOCATE(in_box(mesh%idx%nr(1,1):mesh%idx%nr(2,1)))
  chi = M_ZERO

#if defined(HAVE_MPI) && defined(HAVE_MPI2)
  SAFE_ALLOCATE(start(1:mpi_world%size))
  SAFE_ALLOCATE(end(1:mpi_world%size))
  call multicomm_divide_range(mesh%idx%nr(2,3) - mesh%idx%nr(1,3) + 1, mpi_world%size, start, end)

  start_z = start(mpi_world%rank + 1) - 1 + mesh%idx%nr(1, 3)
  end_z = end(mpi_world%rank + 1) - 1 + mesh%idx%nr(1, 3)

  SAFE_DEALLOCATE_A(start)
  SAFE_DEALLOCATE_A(end)
#else
  start_z = mesh%idx%nr(1, 3)
  end_z = mesh%idx%nr(2, 3)
#endif

  call profiling_in(prof, "MESH_LABEL")
  
  ! We label the points inside the mesh

  do iz = start_z, end_z
    chi(3) = TOFLOAT(iz) * mesh%spacing(3)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      chi(2) = TOFLOAT(iy) * mesh%spacing(2)
      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        chi(1) = TOFLOAT(ix) * mesh%spacing(1)
        call curvilinear_chi2x(sb, sb%latt, cv, chi(:), xx(ix, :))
      end do

      in_box = sb%contains_points(mesh%idx%nr(2,1) - mesh%idx%nr(1,1) + 1, xx)

      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        if (.not.in_box(ix)) cycle
        ASSERT(all((/ix, iy, iz/) <=  mesh%idx%nr(2, 1:3) - mesh%idx%enlarge(1:3)))
        ASSERT(all((/ix, iy, iz/) >=  mesh%idx%nr(1, 1:3) + mesh%idx%enlarge(1:3)))

        mesh%idx%lxyz_inv(ix, iy, iz) = ibset(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)
        do is = 1, stencil%size
          if(stencil%center == is) cycle
          ii = ix + stencil%points(1, is)
          jj = iy + stencil%points(2, is)
          kk = iz + stencil%points(3, is)
          if(any((/ii, jj, kk/) < mesh%idx%nr(1, 1:3)) .or. any((/ii, jj, kk/) >  mesh%idx%nr(2, 1:3))) cycle
          mesh%idx%lxyz_inv(ii, jj, kk) = ibset(mesh%idx%lxyz_inv(ii, jj, kk), ENLARGEMENT_POINT)
        end do
      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(HAVE_MPI2)
  call profiling_in(prof_reduce, "MESH_LABEL_REDUCE")

  npoints = product(nr(2, 1:3) - nr(1, 1:3) + 1)
  call MPI_Allreduce(MPI_IN_PLACE, mesh%idx%lxyz_inv(nr(1, 1), nr(1, 2), nr(1, 3)), npoints, &
    MPI_INTEGER, MPI_BOR, mpi_world%comm, mpi_err)
  ! MPI2 not working could also provoke a segmentation fault in the line above, which we cannot catch.
  if(all(mesh%idx%lxyz_inv(:,:,:) == 0)) then
    message(1) = "Failure of MPI_Allreduce in place for lxyz_inv. MPI2 is not working correctly."
    call messages_fatal(1)
  end if

  call profiling_out(prof_reduce)
#endif

  call profiling_out(prof)

  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(in_box)

  ! count the points
  il = 0
  ik = 0
  do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        if(btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) ik = ik + 1
        if(mesh%idx%lxyz_inv(ix, iy, iz) /= 0) il = il + 1
      end do
    end do
  end do
  mesh%np_part_global = il
  mesh%np_global      = ik

  ASSERT(mesh%np_global > 0)
  ASSERT(mesh%np_part_global > 0)

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
  type(mpi_grp_t) :: mpi_grp

  PUSH_SUB(mesh_init_stage_3)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  call mpi_grp_init(mpi_grp, mc%group_comm(P_STRATEGY_DOMAINS))
  
  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = (mpi_grp%size > 1)

  if(.not. mesh%parallel_in_domains) then
    ! When running parallel, x is computed later.
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part_global, 1:space%dim))
  end if
  
  if(.not. mesh%idx%is_hypercube) then
    call create_x_lxyz()
  else if(.not. mesh%parallel_in_domains) then
    do ip = 1, mesh%np_part_global
      mesh%x(ip, 1:space%dim) = mesh_x_global(mesh, ip, force=.true.)
    end do
  end if

  mesh%mpi_grp = mpi_grp 
  
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

  call mesh_cube_map_init(mesh%cube_map, mesh%idx, mesh%np_global)

  call mesh_get_vol_pp(mesh%sb)

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_3)

contains

  ! ---------------------------------------------------------
  subroutine create_x_lxyz()
    integer :: il, iin, ien, ix, iy, iz, point(1:3)
    integer(8) :: ihilbert
    integer :: ixb, iyb, izb, bsize(1:3)
    type(block_t) :: blk
    integer :: idir, nn, order, size, bits
    FLOAT :: chi(1:MAX_DIM), xx(1:MAX_DIM)

    integer, parameter :: &
      ORDER_BLOCKS     =  1, &
      ORDER_HILBERT    =  2, &
      ORDER_HILBERT_2D =  3

   interface
      subroutine hilbert_index_to_point(dim, nbits, index, point)
        implicit none

        integer,    intent(in)       :: dim
        integer,    intent(in)       :: nbits
        integer(8), intent(in)       :: index
        integer,    intent(out)      :: point !< (1:3)
      end subroutine hilbert_index_to_point
    end interface

    PUSH_SUB(mesh_init_stage_3.create_x_lxyz)

    SAFE_ALLOCATE(mesh%idx%lxyz(1:mesh%np_part_global, 1:MAX_DIM))

    !%Variable MeshOrder
    !%Default blocks
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% This variable controls how the grid points are mapped to a
    !% linear array. This influences the performance of the code.
    !%Option blocks 1
    !% The grid is mapped using small parallelepipedic grids. The size
    !% of the blocks is controlled by <tt>MeshBlockSize</tt>.
    !%Option hilbert 2
    !% (experimental) A Hilbert space-filling curve is used to map the
    !% grid.
    !%Option hilbert_2d 3
    !% (experimental) A Hilbert space-filling curve is used to map the
    !% grid in two of the dimensions.
    !%End

    call parse_variable(namespace, 'MeshOrder', ORDER_BLOCKS, order)

    select case(order)
    case(ORDER_BLOCKS)

      call messages_obsolete_variable(namespace, 'MeshBlockSizeXY', 'MeshBlockSize')
      call messages_obsolete_variable(namespace, 'MeshBlockSizeZ', 'MeshBlockSize')

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

      select case(conf%target_states_block_size)
      case(1)
        bsize(1:3) = (/  2,   1, 200/)
      case(2)
        bsize(1:3) = (/ 10,   4, 200/)
      case(4)
        bsize(1:3) = (/ 10,   4,  80/)
      case(8)
        bsize(1:3) = (/ 10,   2,  30/)
      case(16)
        bsize(1:3) = (/ 15,  15,   4/)
      case(32)
        bsize(1:3) = (/ 15,  15,   2/)
      case(64)
        bsize(1:3) = (/  8,  10,   4/)
      case(128)
        bsize(1:3) = (/  8,   6,   2/)
      case(256)
        bsize(1:3) = (/  4,   6,   2/)
      case default
        bsize(1:3) = (/ 15,  15,   4/)
      end select

      if(parse_block(namespace, 'MeshBlockSize', blk) == 0) then
        nn = parse_block_cols(blk, 0)
        do idir = 1, nn
          call parse_block_integer(blk, 0, idir - 1, bsize(idir))
        end do
      end if 
      ! first we fill the points in the inner mesh
      iin = 0
      ien = mesh%np_global
      do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize(1)
        do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize(2)
          do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsize(3)

            do ix = ixb, min(ixb + bsize(1) - 1, mesh%idx%nr(2,1))
              chi(1) = TOFLOAT(ix) * mesh%spacing(1)
              do iy = iyb, min(iyb + bsize(2) - 1, mesh%idx%nr(2,2))
                chi(2) = TOFLOAT(iy) * mesh%spacing(2)
                do iz = izb, min(izb + bsize(3) - 1, mesh%idx%nr(2,3))
                  chi(3) = TOFLOAT(iz) * mesh%spacing(3)

                  if(btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) then
                    iin = iin + 1
                    il = iin
                  else if (btest(mesh%idx%lxyz_inv(ix, iy, iz), ENLARGEMENT_POINT)) then
                    ien = ien + 1
                    il = ien
                  else
                    cycle
                  end if

                  mesh%idx%lxyz(il, 1) = ix
                  mesh%idx%lxyz(il, 2) = iy
                  mesh%idx%lxyz(il, 3) = iz

                  mesh%idx%lxyz_inv(ix, iy, iz) = il

#ifdef HAVE_MPI
                  if(.not. mesh%parallel_in_domains) then
#endif
                    call curvilinear_chi2x(mesh%sb, mesh%sb%latt, mesh%cv, chi, xx)
                    mesh%x(il, 1:space%dim) = xx(1:space%dim)
#ifdef HAVE_MPI
                  end if
#endif                                   
                end do
              end do
            end do

          end do
        end do
      end do

    case(ORDER_HILBERT)

      call messages_experimental('Hilbert grid ordering')

      size = maxval(mesh%idx%nr(2, 1:space%dim) - mesh%idx%nr(1, 1:space%dim) + 1)

      bits = log2(pad_pow2(size))

      ! since we are using a 64 bit integer, the number of bits is
      ! limited to 21 (3d * 21 = 63)
      ASSERT(bits <= 21)

      iin = 0
      ien = mesh%np_global
      do ihilbert = 0, 2**(bits*3) - 1

        call hilbert_index_to_point(3, bits, ihilbert, point(1))

        ix = point(1) + mesh%idx%nr(1, 1)
        iy = point(2) + mesh%idx%nr(1, 2)
        iz = point(3) + mesh%idx%nr(1, 3)

        if(ix > mesh%idx%nr(2, 1)) cycle
        if(iy > mesh%idx%nr(2, 2)) cycle
        if(iz > mesh%idx%nr(2, 3)) cycle

        if(btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) then
          iin = iin + 1
          il = iin
          !debug output
          !write(12, *) point(1:3), '# ', ihilbert 
        else if (btest(mesh%idx%lxyz_inv(ix, iy, iz), ENLARGEMENT_POINT)) then
          ien = ien + 1
          il = ien
        else
          cycle
        end if

        mesh%idx%lxyz(il, 1) = ix
        mesh%idx%lxyz(il, 2) = iy
        mesh%idx%lxyz(il, 3) = iz

        mesh%idx%lxyz_inv(ix, iy, iz) = il

#ifdef HAVE_MPI
        if(.not. mesh%parallel_in_domains) then
#endif
          chi(1) = TOFLOAT(ix)*mesh%spacing(1)
          chi(2) = TOFLOAT(iy)*mesh%spacing(2)
          chi(3) = TOFLOAT(iz)*mesh%spacing(3)

          call curvilinear_chi2x(mesh%sb, mesh%sb%latt, mesh%cv, chi, xx)
          mesh%x(il, 1:space%dim) = xx(1:space%dim)
#ifdef HAVE_MPI
        end if
#endif                 
      end do

    case(ORDER_HILBERT_2D)

      call messages_experimental('Hilbert 2D grid ordering')

      size = maxval(mesh%idx%nr(2, 1:2) - mesh%idx%nr(1, 1:2) + 1)

      bits = log2(pad_pow2(size))

      bsize = 10
      if(parse_block(namespace, 'MeshBlockSize', blk) == 0) then
        nn = parse_block_cols(blk, 0)
        do idir = 1, nn
          call parse_block_integer(blk, 0, idir - 1, bsize(idir))
        end do
      end if

      ASSERT(bsize(3) > 0)

      iin = 0
      ien = mesh%np_global
      do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsize(3)
        do ihilbert = 0, 2**(bits*2) - 1

          call hilbert_index_to_point(2, bits, ihilbert, point(1))

          ix = point(1) + mesh%idx%nr(1, 1)
          iy = point(2) + mesh%idx%nr(1, 2)

          if(ix > mesh%idx%nr(2, 1)) cycle
          if(iy > mesh%idx%nr(2, 2)) cycle

          do iz = izb, min(izb + bsize(3) - 1, mesh%idx%nr(2,3))

            if(btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) then
              iin = iin + 1
              il = iin
              !debug output
              !write(12, *) ix, iy, iz, '# ', ihilbert 
            else if (btest(mesh%idx%lxyz_inv(ix, iy, iz), ENLARGEMENT_POINT)) then
              ien = ien + 1
              il = ien
            else
              cycle
            end if

            mesh%idx%lxyz(il, 1) = ix
            mesh%idx%lxyz(il, 2) = iy
            mesh%idx%lxyz(il, 3) = iz

            mesh%idx%lxyz_inv(ix, iy, iz) = il

#ifdef HAVE_MPI
            if(.not. mesh%parallel_in_domains) then
#endif
              chi(1) = TOFLOAT(ix)*mesh%spacing(1)
              chi(2) = TOFLOAT(iy)*mesh%spacing(2)
              chi(3) = TOFLOAT(iz)*mesh%spacing(3)

              call curvilinear_chi2x(mesh%sb, mesh%sb%latt, mesh%cv, chi, xx)
              mesh%x(il, 1:space%dim) = xx(1:space%dim)
#ifdef HAVE_MPI
            end if
#endif          
          end do
        end do
      end do

      end select

      ! set the rest to zero
      mesh%idx%lxyz(1:mesh%np_part_global, 4:MAX_DIM) = 0

      call checksum_calculate(1, int(mesh%np_part_global, 8)*space%dim, mesh%idx%lxyz(1, 1), mesh%idx%checksum)

      ASSERT(iin == mesh%np_global)

      if(ien /= mesh%np_part_global) then
        write(message(1), '(a,i9,a,i9)') 'Assertion failure from create_x_lxyz: ien = ', ien, ' /= ', mesh%np_part_global
        write(message(2), '(a)') 'Probably MPI2 is not working correctly.'
        call messages_fatal(2)
      end if

      POP_SUB(mesh_init_stage_3.create_x_lxyz)
    end subroutine create_x_lxyz


  ! ---------------------------------------------------------
  subroutine do_partition()
#ifdef HAVE_MPI
    integer :: ii, jj, ipart, jpart, ip
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

      SAFE_ALLOCATE(nb(1:mpi_grp%size, 1:mpi_grp%size))
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

      SAFE_ALLOCATE(gindex(1:mpi_grp%size))
      SAFE_ALLOCATE(gedges(1:count(nb)))
      
     ! and now generate it
      iedge = 0
      do ipart = 1, mpi_grp%size
        do jpart = 1, mpi_grp%size
          if(nb(ipart, jpart)) then
            iedge = iedge + 1
            gedges(iedge) = jpart - 1
          end if
        end do
        gindex(ipart) = iedge
      end do

      ASSERT(iedge == count(nb))

      reorder = .true.
      call MPI_Graph_create(mpi_grp%comm, mpi_grp%size, gindex, gedges, reorder, graph_comm, mpi_err)

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

    ! Compute mesh%x as it is done in the serial case but only for local points.
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part, 1:space%dim))
    mesh%x(:, :) = M_ZERO
    do ii = 1, mesh%np_part
      jj = mesh_local2global(mesh, ii)
      mesh%x(ii, 1:mesh%sb%dim) = mesh_x_global(mesh, jj)
    end do

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
      chi(1:space%dim) = jj(1:sb%dim)*mesh%spacing(1:space%dim)
      mesh%vol_pp(ip) = mesh%vol_pp(ip)*curvilinear_det_Jac(sb, sb%latt, mesh%cv, mesh%x(ip, 1:sb%dim), chi(1:sb%dim))
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
