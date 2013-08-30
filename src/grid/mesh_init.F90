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
!! $Id$

#include "global.h"

module mesh_init_m
  use cube_m
  use curvilinear_m
  use datasets_m
  use geometry_m
  use global_m
  use hypercube_m
  use index_m
  use io_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_cube_map_m
  use mesh_partition_m
  use messages_m
  use mpi_m
  use multicomm_m
  use ob_grid_m
  use par_vec_m
  use parser_m
  use partition_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use subarray_m

  implicit none
  
  private
  public ::                    &
    mesh_init_stage_1,         &
    mesh_init_stage_2,         &
    mesh_init_stage_3,         &
    mesh_read_lead

  type(profile_t), save :: mesh_init_prof
  
contains

#define ENLARGEMENT_POINT 2
#define INNER_POINT 1
! ---------------------------------------------------------
subroutine mesh_init_stage_1(mesh, sb, cv, spacing, enlarge, ob_grid)
  type(mesh_t),                intent(inout) :: mesh
  type(simul_box_t),   target, intent(in)    :: sb
  type(curvilinear_t), target, intent(in)    :: cv
  FLOAT,                       intent(in)    :: spacing(1:MAX_DIM)
  integer,                     intent(in)    :: enlarge(MAX_DIM)
  type(ob_grid_t),             intent(in)    :: ob_grid


  integer :: idir, jj
  FLOAT   :: x(MAX_DIM), chi(MAX_DIM)
  logical :: out

  PUSH_SUB(mesh_init_stage_1)
  call profiling_in(mesh_init_prof, "MESH_INIT")

  mesh%sb => sb     ! keep an internal pointer
  mesh%idx%sb => sb
  mesh%spacing = spacing ! this number can change in the following
  mesh%use_curvilinear = cv%method /= CURV_METHOD_UNIFORM
  mesh%cv => cv

  ! multiresolution requires the curvilinear coordinates machinery
  mesh%use_curvilinear = mesh%use_curvilinear .or. sb%mr_flag

  mesh%idx%enlarge = enlarge

  if(sb%mr_flag) mesh%idx%enlarge = mesh%idx%enlarge*(2**sb%hr_area%num_radii)

  ! adjust nr
  mesh%idx%nr = 0
  do idir = 1, sb%dim
    chi = M_ZERO
    ! the upper border
    jj = 0
    out = .false.
    do while(.not.out)
      jj = jj + 1
      chi(idir) = real(jj, REAL_PRECISION)*mesh%spacing(idir)
      if ( mesh%use_curvilinear ) then
        call curvilinear_chi2x(sb, cv, chi(1:sb%dim), x(1:sb%dim))
        out = (x(idir) > nearest(sb%lsize(idir), M_ONE))
      else
        out = (chi(idir) > nearest(sb%lsize(idir), M_ONE))
      end if
    end do
    mesh%idx%nr(2, idir) = jj - 1
  end do

  ! we have a symmetric mesh (for now)
  mesh%idx%nr(1,:) = -mesh%idx%nr(2,:)
  
  ! we have to adjust a couple of things for the periodic directions
  do idir = 1, sb%periodic_dim
    !the spacing has to be a divisor of the box size
    mesh%spacing(idir) = sb%lsize(idir)/real(mesh%idx%nr(2, idir))
    !the upper boundary does not have to be included (as it is a copy of the lower boundary)
    mesh%idx%nr(2, idir) = mesh%idx%nr(2, idir) - 1
  end do

  if(ob_grid%open_boundaries) then
    ! The upper boundary must be discarded to preserve periodicity in the leads.
    ! Example in 1D, 2 point unit cell, central region of 8 unit cells, 2
    ! additional unit cells at each end in simulation box:
    !
    ! simulation region:        /-------------------------------\
    ! free state:        . .|. .|. .|. .|. .|. .|. .|. .|. .|. .|. .|. .
    ! scattered state:   . .|. .|. .|. .|. . . . . . . .|. .|. .|. .|. .
    !                       L   |               C               |^  R
    !
    ! The indicated point (^) must be omitted to preserve correct periodicity.
    if(ob_grid%transport_mode) mesh%idx%nr(2, TRANS_DIR) = mesh%idx%nr(2, TRANS_DIR) - 1
  end if

  mesh%idx%ll(:) = mesh%idx%nr(2, :) - mesh%idx%nr(1, :) + 1

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_1)
end subroutine mesh_init_stage_1

! ---------------------------------------------------------
subroutine mesh_read_lead(ob_grid, mesh)
  type(ob_grid_t), target, intent(in)    :: ob_grid
  type(mesh_t),    target, intent(in)    :: mesh

  integer :: il, iunit

  type(mesh_t), pointer :: mm
  integer :: nr(1:2, 1:MAX_DIM)

  PUSH_SUB(mesh_read_lead)

  do il = 1, NLEADS
    mm => ob_grid%lead(il)%mesh
    mm%sb => mesh%sb
    mm%cv => mesh%cv
    mm%parallel_in_domains = mesh%parallel_in_domains
    iunit = io_open(trim(ob_grid%lead(il)%info%restart_dir)//'/'//GS_DIR//'mesh', action='read', is_tmp=.true.)
    call mesh_init_from_file(mm, iunit)
    call io_close(iunit)
    ! Read the lxyz maps.
    nr = mm%idx%nr
    SAFE_ALLOCATE(mm%idx%lxyz(1:mm%np_part, 1:3))
    SAFE_ALLOCATE(mm%idx%lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))
    call mesh_lxyz_init_from_file(mm, trim(ob_grid%lead(il)%info%restart_dir)//'/'//GS_DIR//'lxyz')
  end do

  POP_SUB(mesh_read_lead)
end subroutine mesh_read_lead

! ---------------------------------------------------------
!> This subroutine checks if every grid point belongs to the internal
!! mesh, based on the global lxyz_inv matrix. Afterwards, it counts
!! how many points has the mesh and the enlargement.
subroutine mesh_init_stage_2(mesh, sb, geo, cv, stencil)
  type(mesh_t),        intent(inout) :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(geometry_t),    intent(in)    :: geo
  type(curvilinear_t), intent(in)    :: cv
  type(stencil_t),     intent(in)    :: stencil

  integer :: il, ik, ix, iy, iz, is
  integer :: newi, newj, newk, ii, jj, kk, dx, dy, dz, i_lev
  integer :: jx, jy, jz, res_counter, j_counter
  FLOAT   :: chi(MAX_DIM)
  integer :: nr(1:2, 1:MAX_DIM), res, n_mod
  logical, allocatable :: in_box(:)
  FLOAT,   allocatable :: xx(:, :)
  real(8), parameter :: DELTA = CNST(1e-12)
  integer :: start_z, end_z
  type(profile_t), save :: prof
#if defined(HAVE_MPI) && defined(HAVE_MPI2)
  type(profile_t), save :: prof_reduce
  integer :: npoints
  integer, allocatable :: start(:), end(:)
#endif

  PUSH_SUB(mesh_init_stage_2)
  call profiling_in(mesh_init_prof)

  ! enlarge mesh for boundary points
  mesh%idx%nr(1,:) = mesh%idx%nr(1,:) - mesh%idx%enlarge(:)
  mesh%idx%nr(2,:) = mesh%idx%nr(2,:) + mesh%idx%enlarge(:)
  if(mesh%idx%sb%box_shape == HYPERCUBE) then
    call hypercube_init(mesh%idx%hypercube, sb%dim, mesh%idx%nr, mesh%idx%enlarge(1))
    mesh%np_part_global = hypercube_number_total_points(mesh%idx%hypercube)
    mesh%np_global      = hypercube_number_inner_points(mesh%idx%hypercube)
    nullify(mesh%resolution)
    nullify(mesh%idx%lxyz_inv)
    nullify(mesh%idx%lxyz)
    call profiling_out(mesh_init_prof)
    POP_SUB(mesh_init_stage_2)
    return
  end if

  nr = mesh%idx%nr

  ! allocate the xyz arrays
  SAFE_ALLOCATE(mesh%idx%lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))

  if(sb%mr_flag) then 
    SAFE_ALLOCATE(mesh%resolution(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))
    mesh%resolution(:,:,:) = 0
  else
    nullify(mesh%resolution)
  end if

  mesh%idx%lxyz_inv(:,:,:) = 0
  res = 1
  SAFE_ALLOCATE(xx(1:MAX_DIM, mesh%idx%nr(1,1):mesh%idx%nr(2,1)))
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
    chi(3) = real(iz, REAL_PRECISION) * mesh%spacing(3) + sb%box_offset(3)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      chi(2) = real(iy, REAL_PRECISION) * mesh%spacing(2) + sb%box_offset(2)
      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        chi(1) = real(ix, REAL_PRECISION) * mesh%spacing(1) + sb%box_offset(1)
        call curvilinear_chi2x(sb, cv, chi(:), xx(:, ix))
      end do

      call simul_box_in_box_vec(sb, geo, mesh%idx%nr(2,1) - mesh%idx%nr(1,1) + 1, xx, in_box)

      do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
        ! With multiresolution, only inner (not enlargement) points are marked now
        if(sb%mr_flag) then
          if (in_box(ix) ) then
            ! First check: is the point beyond the multiresolution areas
            n_mod = 2**sb%hr_area%num_radii
            if (sum((xx(:,ix)-sb%hr_area%center(:))**2) > sb%hr_area%radius(sb%hr_area%num_radii)**2 .and. &
                 mod(ix, n_mod) == 0 .and. mod(iy, n_mod) == 0 .and. mod(iz,n_mod) == 0) then
              mesh%idx%lxyz_inv(ix, iy, iz) = ibset(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)
            end if
            ! Other option: must be inside the multiresolution area and satisfy coordinate index conditions
            if(.not.btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) then
              do i_lev = 1,sb%hr_area%num_radii
                n_mod = 2**(i_lev-1)
                if( sum((xx(:,ix)-sb%hr_area%center(:))**2) < sb%hr_area%radius(i_lev)**2 + DELTA .and. &
                    mod(ix, n_mod) == 0 .and. mod(iy, n_mod) == 0 .and. mod(iz,n_mod) == 0) then
                  mesh%idx%lxyz_inv(ix, iy, iz) = ibset(mesh%idx%lxyz_inv(ix,iy, iz), INNER_POINT)
                end if
              end do
            end if
          end if
        else ! the usual way: mark both inner and enlargement points
          if (in_box(ix)) then
            mesh%idx%lxyz_inv(ix, iy, iz) = ibset(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)
            do is = 1, stencil%size
              if(stencil%center == is) cycle
              ii = ix + stencil%points(1, is)
              jj = iy + stencil%points(2, is)
              kk = iz + stencil%points(3, is)
              if(any((/ii, jj, kk/) < mesh%idx%nr(1, 1:3)) .or. any((/ii, jj, kk/) >  mesh%idx%nr(2, 1:3))) cycle
              mesh%idx%lxyz_inv(ii, jj, kk) = ibset(mesh%idx%lxyz_inv(ii, jj, kk), ENLARGEMENT_POINT)
            end do
          end if
        end if
      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(HAVE_MPI2)
  call profiling_in(prof_reduce, "MESH_LABEL_REDUCE")

  npoints = product(nr(2, 1:3) - nr(1, 1:3) + 1)
  call MPI_Allreduce(MPI_IN_PLACE, mesh%idx%lxyz_inv(nr(1, 1), nr(1, 2), nr(1, 3)), npoints, &
    MPI_INTEGER, MPI_BOR, mpi_world%comm, mpi_err)

  call profiling_out(prof_reduce)
#endif

  call profiling_out(prof)

  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(in_box)

  if(sb%mr_flag) then
    ! Calculate the resolution for each point and label the enlargement points
    do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
      chi(3) = real(iz, REAL_PRECISION) * mesh%spacing(3) + sb%box_offset(3)
      do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        chi(2) = real(iy, REAL_PRECISION) * mesh%spacing(2) + sb%box_offset(2)
        do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
          chi(1) = real(ix, REAL_PRECISION) * mesh%spacing(1) + sb%box_offset(1)
           ! skip if not inner point
          if(.not.btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) cycle
          res = -1
          res_counter = 0
          do while(.true.)
            res_counter = res_counter + 1
            ! loop through Cartesian axes (both directions)
            do j_counter = 1,6
              select case(j_counter)
                case(1)
                  jx=ix-res_counter; jy=iy; jz=iz
                case(2)
                  jx=ix+res_counter; jy=iy; jz=iz
                case(3)
                  jx=ix; jy=iy-res_counter; jz=iz
                case(4)
                  jx=ix; jy=iy+res_counter; jz=iz
                case(5)
                  jx=ix; jy=iy; jz=iz-res_counter
                case(6)
                  jx=ix; jy=iy; jz=iz+res_counter
              end select
              if(any((/jx, jy, jz/) < mesh%idx%nr(1, 1:3)) .or. &
                 any((/jx, jy, jz/) > mesh%idx%nr(2, 1:3))) cycle
              ! exit after finding neighboring inner point
              if(btest(mesh%idx%lxyz_inv(jx, jy, jz), INNER_POINT)) then
                res = res_counter
                exit
              end if
            end do
            if(res /= -1) exit
          end do
          mesh%resolution(ix, iy, iz) = res
          ! mark the enlargement points
          do is = 1, stencil%size
            if(stencil%center == is) cycle
            ii = ix + res*stencil%points(1, is)
            jj = iy + res*stencil%points(2, is)
            kk = iz + res*stencil%points(3, is)
            if(any((/ii, jj, kk/) < mesh%idx%nr(1, 1:3)) .or. any((/ii, jj, kk/) >  mesh%idx%nr(2, 1:3))) cycle
            mesh%idx%lxyz_inv(ii, jj, kk) = ibset(mesh%idx%lxyz_inv(ii, jj, kk), ENLARGEMENT_POINT)
            ! If the point is not an inner point, and if its resolution can be decreased, do it now
            if(.not. btest(mesh%idx%lxyz_inv(ii, jj, kk), INNER_POINT)) then
              if(mesh%resolution(ii, jj, kk) == 0 .or. mesh%resolution(ii, jj, kk) > res) &
                   mesh%resolution(ii, jj, kk) = res
            end if
          end do
        end do
      end do
    end do
  end if

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

  ! Errors occur during actual calculation if resolution interfaces are too close to each other. The
  ! following routine checks that everything will be ok.
  if(sb%mr_flag) then
    ! loop through all interpolation points and check that all points used for interpolation exist
    do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
      do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
        do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)

          i_lev = mesh%resolution(ix,iy,iz)

          ! include enlargement points that are neither inner points nor outer boundary points.
          if( .not. btest(mesh%idx%lxyz_inv(ix, iy, iz),ENLARGEMENT_POINT)) cycle
          if(  btest(mesh%idx%lxyz_inv(ix, iy, iz), INNER_POINT)) cycle
          if(  i_lev == 2**mesh%sb%hr_area%num_radii ) cycle

          ! the value of point (ix,iy,iz) is going to be interpolated
          dx = abs(mod(ix, 2**(i_lev)))
          dy = abs(mod(iy, 2**(i_lev)))
          dz = abs(mod(iz, 2**(i_lev)))
          do ii = 1, mesh%sb%hr_area%interp%nn
            do jj = 1, mesh%sb%hr_area%interp%nn
              do kk = 1, mesh%sb%hr_area%interp%nn
                newi = ix + mesh%sb%hr_area%interp%posi(ii)*dx
                newj = iy + mesh%sb%hr_area%interp%posi(jj)*dy
                newk = iz + mesh%sb%hr_area%interp%posi(kk)*dz
                if(any((/newi, newj, newk/) <  mesh%idx%nr(1, 1:3)) .or. &
                   any((/newi, newj, newk/) >  mesh%idx%nr(2, 1:3)) .or. &
                   mesh%idx%lxyz_inv(newi,newj,newk) == 0) then

                     message(1) = 'Multiresolution radii are too close to each other (or outer boundary)'
                     write(message(2),'(7I4)') ix,iy,iz,newi,newj,newk,mesh%resolution(ix,iy,iz)
                     call messages_fatal(2)
                end if
              end do
            end do
          end do
 
        end do
      end do
    end do

  end if

  call profiling_out(mesh_init_prof)
  POP_SUB(mesh_init_stage_2)
end subroutine mesh_init_stage_2


! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------
subroutine mesh_init_stage_3(mesh, stencil, mpi_grp, parent)
  type(mesh_t),              intent(inout) :: mesh
  type(stencil_t), optional, intent(in)    :: stencil
  type(mpi_grp_t), optional, intent(in)    :: mpi_grp
  type(mesh_t),    optional, intent(in)    :: parent

  integer :: ip

  PUSH_SUB(mesh_init_stage_3)
  call profiling_in(mesh_init_prof)

  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = .false.
  if(present(mpi_grp)) mesh%parallel_in_domains = .true.


  if(.not. mesh%parallel_in_domains) then
    ! When running parallel, x is computed later.
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part_global, 1:MAX_DIM))
  end if
  
  if(mesh%idx%sb%box_shape /= HYPERCUBE) then
    call create_x_lxyz()
  else if(.not. mesh%parallel_in_domains) then
    do ip = 1, mesh%np_part_global
      mesh%x(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip, force=.true.)
    end do
  end if

  if(mesh%parallel_in_domains) then
    ASSERT(present(stencil))
    
    call do_partition()
  else
    call mpi_grp_init(mesh%mpi_grp, -1)

    ! When running serially those two are the same.
    mesh%np      = mesh%np_global
    mesh%np_part = mesh%np_part_global

    ! These must be initialized for vec_gather, vec_scatter to work
    ! as copy operations when running without domain parallelization.
    mesh%vp%np_global = mesh%np_global
    mesh%vp%npart = 1
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


    PUSH_SUB(mesh_init_stage_3.create_x_lxyz)

    SAFE_ALLOCATE(mesh%idx%lxyz(1:mesh%np_part_global, 1:MAX_DIM))

    !%Variable MeshOrder
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% This variable controls how the grid points are mapped to a
    !% linear array. This influences the performance of the code. The
    !% default is blocks.
    !%Option blocks 1
    !% The grid is mapped using small parallelepipedic grids. The size
    !% of the blocks is controlled by <tt>MeshBlockSize</tt>.
    !%Option hilbert 2
    !% (experimental) A Hilbert space filling curve is used to map the
    !% grid.
    !%Option hilbert_2d 3
    !% (experimental) A Hilbert space filling curve is used to map the
    !% grid in two of the dimensions.
    !%End

    call parse_integer(datasets_check('MeshOrder'), ORDER_BLOCKS, order)

    select case(order)
    case(ORDER_BLOCKS)

      call messages_obsolete_variable('MeshBlockSizeXY', 'MeshBlockSize')
      call messages_obsolete_variable('MeshBlockSizeZ', 'MeshBlockSize')

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

      if(parse_block('MeshBlockSize', blk) == 0) then
        nn = parse_block_cols(blk, 0)
        do idir = 1, nn
          call parse_block_integer(blk, 0, idir - 1, bsize(idir))
        end do
      end if

      ! When using open boundaries we need to have a mesh block-size of 1
      if (parse_block(datasets_check('OpenBoundaries'), blk) == 0) then
        if (any(bsize > 1)) then
          message(1) = 'When in transport mode, the block-ordering'
          message(2) = 'of the mesh points cannot be chosen freely.'
          message(3) = 'All the values of MeshBlockSize have to be'
          message(4) = 'initialized with the value 1.'
          call messages_fatal(4)
        end if
      end if

      ! first we fill the points in the inner mesh
      iin = 0
      ien = mesh%np_global
      do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize(1)
        do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize(2)
          do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsize(3)

            do ix = ixb, min(ixb + bsize(1) - 1, mesh%idx%nr(2,1))
              chi(1) = real(ix, REAL_PRECISION) * mesh%spacing(1) + mesh%sb%box_offset(1)
              do iy = iyb, min(iyb + bsize(2) - 1, mesh%idx%nr(2,2))
                chi(2) = real(iy, REAL_PRECISION) * mesh%spacing(2) + mesh%sb%box_offset(2)
                do iz = izb, min(izb + bsize(3) - 1, mesh%idx%nr(2,3))
                  chi(3) = real(iz, REAL_PRECISION) * mesh%spacing(3) + mesh%sb%box_offset(3)

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
                    call curvilinear_chi2x(mesh%sb, mesh%cv, chi, xx)
                    mesh%x(il, 1:MAX_DIM) = xx(1:MAX_DIM)
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

      size = maxval(mesh%idx%nr(2, 1:mesh%sb%dim) - mesh%idx%nr(1, 1:mesh%sb%dim) + 1)

      bits = log2(pad_pow2(size))

      ! since we are using a 64 bit integer, the number of bits is
      ! limited to 21 (3d * 21 = 63)
      ASSERT(bits <= 21)

      iin = 0
      ien = mesh%np_global
      do ihilbert = 0, 2**(bits*3) - 1

        call hilbert_index_to_point(3, bits, ihilbert, point)

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
          chi(1) = real(ix, REAL_PRECISION)*mesh%spacing(1) + mesh%sb%box_offset(1)
          chi(2) = real(iy, REAL_PRECISION)*mesh%spacing(2) + mesh%sb%box_offset(2)
          chi(3) = real(iz, REAL_PRECISION)*mesh%spacing(3) + mesh%sb%box_offset(3)

          call curvilinear_chi2x(mesh%sb, mesh%cv, chi, xx)
          mesh%x(il, 1:MAX_DIM) = xx(1:MAX_DIM)
#ifdef HAVE_MPI
        end if
#endif                 
      end do

    case(ORDER_HILBERT_2D)

      call messages_experimental('Hilbert 2D grid ordering')

      size = maxval(mesh%idx%nr(2, 1:2) - mesh%idx%nr(1, 1:2) + 1)

      bits = log2(pad_pow2(size))

      bsize = 10
      if(parse_block('MeshBlockSize', blk) == 0) then
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

          call hilbert_index_to_point(2, bits, ihilbert, point)

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
              chi(1) = real(ix, REAL_PRECISION)*mesh%spacing(1) + mesh%sb%box_offset(1)
              chi(2) = real(iy, REAL_PRECISION)*mesh%spacing(2) + mesh%sb%box_offset(2)
              chi(3) = real(iz, REAL_PRECISION)*mesh%spacing(3) + mesh%sb%box_offset(3)

              call curvilinear_chi2x(mesh%sb, mesh%cv, chi, xx)
              mesh%x(il, 1:MAX_DIM) = xx(1:MAX_DIM)
#ifdef HAVE_MPI
            end if
#endif          
          end do
        end do
      end do

      end select

      ! set the rest to zero
      mesh%idx%lxyz(1:mesh%np_part_global, 4:MAX_DIM) = 0

      call checksum_calculate(1, mesh%np_part_global*mesh%sb%dim, mesh%idx%lxyz(1, 1), mesh%idx%checksum)

      ASSERT(iin == mesh%np_global)
      ASSERT(ien == mesh%np_part_global)

      POP_SUB(mesh_init_stage_3.create_x_lxyz)
    end subroutine create_x_lxyz


  ! ---------------------------------------------------------
  subroutine do_partition()
#ifdef HAVE_MPI
    integer :: ii, jj, ipart, jpart, ip, ix, iy, iz
    integer, allocatable :: gindex(:), gedges(:)
    logical, allocatable :: nb(:, :)
    integer              :: idx(1:MAX_DIM), jx(1:MAX_DIM)
    integer              :: graph_comm, iedge, nnb, npart
    logical              :: use_topo, reorder, partition_print
    integer              :: ierr

    logical :: from_scratch, has_virtual_partition = .false.
    integer :: vsize !< 'virtual' partition size

    PUSH_SUB(mesh_init_stage_3.do_partition)

    mesh%mpi_grp = mpi_grp

    !%Variable MeshPartitionFromScratch
    !%Type logical
    !%Default false
    !%Section Execution::Parallelization
    !%Description
    !% If set to no (the default) Octopus will try to use the mesh
    !% partition from a previous run if available.
    !%End
    call parse_logical(datasets_check('MeshPartitionFromScratch'), .false., from_scratch)

    ierr = -1
    if(.not. from_scratch) call mesh_partition_read(mesh, ierr)


    if(ierr /= 0) then

      !%Variable MeshPartitionVirtualSize
      !%Type integer
      !%Default mesh mpi_grp size
      !%Section Execution::Parallelization
      !%Description
      !% Gives the possibility to change the partition nodes
      !% Afterward, it crashes
      !%End
      call parse_integer(datasets_check('MeshPartitionVirtualSize'), mesh%mpi_grp%size, vsize)
      
      if (vsize /= mesh%mpi_grp%size) then
        write(message(1),'(a,I7)') "Changing the partition size to", vsize
        write(message(2),'(a)') "The execution will crash"
        call messages_warning(2)
        has_virtual_partition = .true.
      else 
        has_virtual_partition = .false.
      end if
      
      if(.not. present(parent)) then
        call mesh_partition(mesh, stencil, vsize)
      else
        ! if there is a parent grid, use its partition
        call mesh_partition_from_parent(mesh, parent)
      end if
      
      !Do the partioning of the boundary points
      call mesh_partition_boundaries(mesh, stencil, vsize)

      !Now that we have the partitions, we save them
      call mesh_partition_write(mesh, vsize)
    end if

    !At the moment we still need the global partition. This will be removed in near future.
    SAFE_ALLOCATE(mesh%vp%part_vec(1:mesh%np_part_global))
    call partition_get_global(mesh%inner_partition, mesh%vp%part_vec(1:mesh%np_global))
    call partition_get_global(mesh%bndry_partition, mesh%vp%part_vec(mesh%np_global+1:mesh%np_part_global))      

    if (has_virtual_partition) then  
      write(message(1),'(a)') "Execution has ended"
      write(message(2),'(a)') "If you want to run your system, you don`t have to use MeshPartitionVirtualSize"
      call messages_warning(2)
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
    call parse_logical(datasets_check('MeshUseTopology'), .false., use_topo)

    if(use_topo) then
      ! this should be integrated in vec_init

      ! generate a table of neighbours

      SAFE_ALLOCATE(nb(1:mpi_grp%size, 1:mpi_grp%size))
      nb = .false.

      do ip = 1, mesh%np_global
        ipart = mesh%vp%part_vec(ip)
        call index_to_coords(mesh%idx, mesh%sb%dim, ip, idx)
        do jj = 1, stencil%size
          jx(1:MAX_DIM) = idx(1:MAX_DIM) + stencil%points(1:MAX_DIM, jj)
          if(all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
            jpart = mesh%vp%part_vec(index_from_coords(mesh%idx, mesh%sb%dim, jx))
            if(ipart /= jpart ) nb(ipart, jpart) = .true.
          end if
        end do
      end do

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
      call MPI_Graph_create(mpi_grp%comm, mpi_grp%size, gindex(1), gedges(1), reorder, graph_comm, mpi_err)

      SAFE_DEALLOCATE_A(nb)
      SAFE_DEALLOCATE_A(gindex)
      SAFE_DEALLOCATE_A(gedges)

    else

      call MPI_Comm_dup(mpi_grp%comm, graph_comm, mpi_err)

    end if

    ! we have a new communicator
    call mpi_grp_init(mesh%mpi_grp, graph_comm)

    npart = partition_get_npart(mesh%inner_partition)
    call vec_init(graph_comm, 0, mesh%np_global, mesh%np_part_global, mesh%idx, stencil,&
         mesh%sb%dim, mesh%sb%periodic_dim, mesh%inner_partition, mesh%bndry_partition, mesh%vp)

    ! check the number of ghost neighbours in parallel
    nnb = 0
    jpart =  mesh%vp%partno
    do ipart = 1, mesh%vp%npart
      if (ipart == jpart) cycle
      if (mesh%vp%np_ghost_neigh_partno(ipart) /= 0) nnb = nnb + 1
    end do
    ASSERT(nnb >= 0 .and. nnb < mesh%vp%npart)

    ! Set local point numbers.
    mesh%np      = mesh%vp%np_local
    mesh%np_part = mesh%np + mesh%vp%np_ghost + mesh%vp%np_bndry

    ! Compute mesh%x as it is done in the serial case but only for local points.
    ! x consists of three parts: the local points, the
    ! ghost points, and the boundary points; in this order
    ! (just as for any other vector, which is distributed).
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part, 1:MAX_DIM))
    mesh%x(:, :) = M_ZERO
    ! Do the inner points
    do ii = 1, mesh%np
      jj = mesh%vp%local(mesh%vp%xlocal + ii - 1)
      mesh%x(ii, 1:MAX_DIM) = mesh_x_global(mesh, jj)
    end do
    ! Do the ghost points
    do ii = 1, mesh%vp%np_ghost
      jj = mesh%vp%ghost(mesh%vp%xghost + ii - 1) 
      mesh%x(ii+mesh%np, 1:MAX_DIM) = mesh_x_global(mesh, jj)
    end do
    ! Do the boundary points
    do ii = 1, mesh%vp%np_bndry
      jj = mesh%vp%bndry(mesh%vp%xbndry + ii - 1)
      mesh%x(ii + mesh%np + mesh%vp%np_ghost, 1:MAX_DIM) = mesh_x_global(mesh, jj)
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
    call parse_logical(datasets_check('PartitionPrint'), .true., partition_print)
    
    if (partition_print) then
      call mesh_partition_write_info(mesh, stencil, mesh%vp%part_vec)
      call mesh_partition_messages_debug(mesh)
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

    integer :: ix, iy, iz, dx, dy, dz, newi, newj, newk, ii, lii, ljj, lkk, nn
    FLOAT,   allocatable :: pos(:), ww(:), vol_tmp(:, :, :)
    integer, allocatable :: posi(:)
    integer :: n_mod, i_lev, nr(1:2, 1:MAX_DIM)
    real(8), parameter :: DELTA = CNST(1e-12)
 
#if defined(HAVE_MPI)
    integer :: kk
#endif

    PUSH_SUB(mesh_init_stage_3.mesh_get_vol_pp)

    np = 1
    if(mesh%use_curvilinear) np = mesh%np_part

    SAFE_ALLOCATE(mesh%vol_pp(1:np))

    forall(ip = 1:np) mesh%vol_pp(ip) = product(mesh%spacing(1:sb%dim))
    jj(sb%dim + 1:MAX_DIM) = 0

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ! Do the inner points.
      do ip = 1, min(np, mesh%np)
        kk = mesh%vp%local(mesh%vp%xlocal + ip - 1)
        call index_to_coords(mesh%idx, sb%dim, kk, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
        mesh%vol_pp(ip) = mesh%vol_pp(ip)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip, :), chi(1:sb%dim))
      end do

      if(mesh%use_curvilinear) then
        ! Do the ghost points.
        do ip = 1, mesh%vp%np_ghost
          kk = mesh%vp%ghost(mesh%vp%xghost + ip - 1)
          call index_to_coords(mesh%idx, sb%dim, kk, jj)
          chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
          mesh%vol_pp(ip + mesh%np) = &
            mesh%vol_pp(ip + mesh%np)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip + mesh%np, :), chi(1:sb%dim))
        end do
        ! Do the boundary points.
        do ip = 1, mesh%vp%np_bndry
          kk = mesh%vp%bndry(mesh%vp%xbndry + ip - 1)
          call index_to_coords(mesh%idx, sb%dim, kk, jj)
          chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
          mesh%vol_pp(ip+mesh%np+mesh%vp%np_ghost) = &
            mesh%vol_pp(ip+mesh%np+mesh%vp%np_ghost) &
            *curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip+mesh%np+mesh%vp%np_ghost, :), chi(1:sb%dim))
        end do
      end if
#endif
    else ! serial mode

      if(mesh%sb%mr_flag) then

        message(1) = 'Info: Point volumes are calculated by solving interpolation coefficients for the intermediate points.'
        call messages_info(1)

        ! The following interpolation routine is essentially the same as in the calculation of the Laplacian

        nn = 2*mesh%sb%hr_area%interp%order

        SAFE_ALLOCATE(ww(1:nn))
        SAFE_ALLOCATE(pos(1:nn))
        SAFE_ALLOCATE(posi(1:nn))

        do ii = 1, mesh%sb%hr_area%interp%order
          posi(ii) = 1 + 2*(ii - 1)
          posi(mesh%sb%hr_area%interp%order + ii) = -posi(ii)
          pos(ii) =  posi(ii)
          pos(mesh%sb%hr_area%interp%order + ii) = -pos(ii)
        end do

        call interpolation_coefficients(nn, pos, M_ZERO, ww)

        ! volumes are initialized even for the intermediate points
        nr(:,:) = mesh%idx%nr(:,:)
        SAFE_ALLOCATE(vol_tmp(nr(1,1):nr(2,1),nr(1,2):nr(2,2),nr(1,3):nr(2,3)))
        vol_tmp(:,:,:) = product(mesh%spacing(1:sb%dim))

        ! The idea is that in the first i_lev loop we find intermediate
        ! points that are odd, i.e. at least one of their indices cannot
        ! be divided by 2 (note that then n_mod=2**1=1).
        !
        ! In the second loop we accept only those intermediate points that
        ! have at least one index that cannot be divided by 4. This
        ! rules out some even points that were included in the first loop.
        !
        ! This continues until the last resolution level. In each step
        ! the point volumes of the neighboring points area modified.

        do i_lev = 1, sb%hr_area%num_radii

          write (message(1),'(a,I2,a,I2)') 'Info: Point volume calculation at stage ',i_lev,'/',sb%hr_area%num_radii
          call messages_info(1)

          ! loop through _all_ the points
          do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
            do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
              do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
 
                ! Skip ordinary points
                if(mesh%idx%lxyz_inv(ix,iy,iz) > 0 .and. &
                     mesh%idx%lxyz_inv(ix,iy,iz) <= mesh%np) cycle

                ! Is it the kind of intermediate point we are looking for?
                n_mod = 2**i_lev
                dx = abs(mod(ix, n_mod))
                dy = abs(mod(iy, n_mod))
                dz = abs(mod(iz, n_mod))
                if(dx+dy+dz == M_ZERO) cycle

                if(abs(vol_tmp(ix, iy, iz)) < DELTA) cycle

                ! The present point (ix,iy,iz) is an intermediate one. When
                ! calculating integrals, the value of the integrand is
                ! interpolated from the neighboring ones, i.e. the values of
                ! the neighboring points are added up with different weights.
                ! The following loop goes through the neighboring points and
                ! modifies their weights, i.e. their volumes.

                do lii = 1, nn
                  do ljj = 1, nn
                    do lkk = 1, nn
                      newi = ix + posi(lii)*dx
                      newj = iy + posi(ljj)*dy
                      newk = iz + posi(lkk)*dz
                      if(any((/newi, newj, newk/) <  mesh%idx%nr(1, 1:3)) .or. &
                         any((/newi, newj, newk/) >  mesh%idx%nr(2, 1:3))) cycle
                      vol_tmp(newi, newj, newk) = vol_tmp(newi, newj, newk) + &
                         vol_tmp(ix, iy, iz) * ww(lii)*ww(ljj)*ww(lkk)
                    end do
                  end do
                end do
                vol_tmp(ix, iy, iz) = M_ZERO
              end do
            end do
          end do

        end do

        ! the volumes are now in vol_tmp table. Move them to vol_pp
        do ip = 1, mesh%np
          ix = mesh%idx%lxyz(ip, 1)
          iy = mesh%idx%lxyz(ip, 2)
          iz = mesh%idx%lxyz(ip, 3)
          mesh%vol_pp(ip) = vol_tmp(ix,iy,iz)
        end do
       
        write (message(1),'(a,F26.12)') 'Info: Point volume calculation finished. Total volume :',sum(mesh%vol_pp(1:mesh%np))
        call messages_info(1)

        SAFE_DEALLOCATE_A(ww)
        SAFE_DEALLOCATE_A(pos)
        SAFE_DEALLOCATE_A(posi)
        SAFE_DEALLOCATE_A(vol_tmp)

      else ! no multiresolution

        do ip = 1, np
          call index_to_coords(mesh%idx, sb%dim, ip, jj)
          chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
          mesh%vol_pp(ip) = mesh%vol_pp(ip)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip, 1:sb%dim), chi(1:sb%dim))
        end do

      end if
    end if

    if(mesh%use_curvilinear) then
      mesh%volume_element = CNST(1.0)
    else
      mesh%volume_element = mesh%vol_pp(1)
    end if

    POP_SUB(mesh_init_stage_3.mesh_get_vol_pp)

  end subroutine mesh_get_vol_pp

end subroutine mesh_init_stage_3

end module mesh_init_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
