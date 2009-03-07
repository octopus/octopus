! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module mesh_init_m
  use datasets_m
  use curvlinear_m
  use geometry_m
  use global_m
  use hypercube_m
  use io_m
  use math_m
  use mesh_m
  use index_m
  use messages_m
  use multicomm_m
  use mpi_m
  use loct_parser_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use stencil_star_m
  use units_m
  use zoltan_m

  implicit none
  
  private
  public ::                    &
    mesh_init_stage_1,         &
    mesh_init_stage_2,         &
    mesh_init_stage_3

  type(profile_t), save :: mesh_init_prof
  
contains

#define ENLARGEMENT_POINT 2
#define INNER_POINT 1
! ---------------------------------------------------------
subroutine mesh_init_stage_1(mesh, sb, geo, cv, enlarge)
  type(mesh_t),              intent(inout) :: mesh
  type(simul_box_t), target, intent(in)    :: sb
  type(geometry_t),          intent(in)    :: geo
  type(curvlinear_t),        intent(in)    :: cv
  integer,                   intent(in)    :: enlarge(MAX_DIM)

  integer :: idir, jj
  FLOAT   :: x(MAX_DIM), chi(MAX_DIM)
  logical :: out

  call push_sub('mesh_init.mesh_init_stage_1')
  call profiling_in(mesh_init_prof, "MESH_INIT")

  mesh%sb => sb     ! keep an internal pointer
  mesh%idx%sb => sb
  mesh%h  =  sb%h ! this number can change in the following
  mesh%use_curvlinear = cv%method.ne.CURV_METHOD_UNIFORM

  ! multiresolution requires the curvilinear coordinates machinery
  mesh%use_curvlinear = mesh%use_curvlinear .or. simul_box_multires(sb)

  mesh%idx%enlarge = enlarge

  if(simul_box_multires(sb)) mesh%idx%enlarge = mesh%idx%enlarge*2

  ! adjust nr
  mesh%idx%nr = 0
  do idir = 1, sb%dim
    chi = M_ZERO
    ! the upper border
    jj = 0
    out = .false.
    do while(.not.out)
      jj = jj + 1
      chi(idir) = real(jj, REAL_PRECISION)*mesh%h(idir)
      if ( mesh%use_curvlinear ) then
        call curvlinear_chi2x(sb, geo, cv, chi(1:sb%dim), x(1:sb%dim))
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
    mesh%h(idir) = sb%lsize(idir)/real(mesh%idx%nr(2, idir))
    !the upper boundary does not have to be included (as it is a copy of the lower boundary)
    mesh%idx%nr(2, idir) = mesh%idx%nr(2, idir) - 1
  end do

  if(sb%open_boundaries) then
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
    mesh%idx%nr(2, TRANS_DIR) = mesh%idx%nr(2, TRANS_DIR) - 1

    call mesh_read_lead()
  else
    nullify(mesh%lead_unit_cell)
  end if

  mesh%idx%ll(:) = mesh%idx%nr(2, :) - mesh%idx%nr(1, :) + 1

  call profiling_out(mesh_init_prof)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine mesh_read_lead()
    integer :: il, iunit

    integer               :: alloc_size
    type(mesh_t), pointer :: m
    integer :: nr(1:2, 1:MAX_DIM)

    call push_sub('mesh_init.mesh_init_stage1_mesh_read_lead')

    ALLOCATE(mesh%lead_unit_cell(NLEADS), NLEADS)

    do il = 1, NLEADS
      m => mesh%lead_unit_cell(il)
      m%sb => mesh%sb
      iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/mesh', action='read', is_tmp=.true.)
      call mesh_init_from_file(m, iunit)
      call io_close(iunit)
      ! Read the lxyz maps.
      nr = m%idx%nr
      alloc_size = (nr(2, 1) - nr(1, 1) + 1)*(nr(2, 2) - nr(1, 2) + 1)*(nr(2, 3) - nr(1, 3) + 1)
      ALLOCATE(m%idx%Lxyz(m%np_part, 3), m%np_global*3)
      ALLOCATE(m%idx%Lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)), alloc_size)
      iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/Lxyz', action='read', is_tmp=.true.)
      call mesh_lxyz_init_from_file(m, iunit)
      call io_close(iunit)
    end do

    call pop_sub()
  end subroutine mesh_read_lead
end subroutine mesh_init_stage_1


! ---------------------------------------------------------
subroutine mesh_init_stage_2(mesh, sb, geo, cv, stencil)
  type(mesh_t),       intent(inout) :: mesh
  type(simul_box_t),  intent(in)    :: sb
  type(geometry_t),   intent(in)    :: geo
  type(curvlinear_t), intent(in)    :: cv
  type(stencil_t),    intent(in)    :: stencil

  integer :: i, j, k, il, ik, ix, iy, iz, is
  FLOAT   :: chi(MAX_DIM)
  integer :: nr(1:2, 1:MAX_DIM), res
  logical :: inside

  call push_sub('mesh_init.mesh_init_stage_2')
  call profiling_in(mesh_init_prof)

  ! enlarge mesh for boundary points
  mesh%idx%nr(1,:) = mesh%idx%nr(1,:) - mesh%idx%enlarge(:)
  mesh%idx%nr(2,:) = mesh%idx%nr(2,:) + mesh%idx%enlarge(:)

  if(mesh%idx%sb%box_shape == HYPERCUBE) then
    call hypercube_init(mesh%idx%hypercube, sb%dim, mesh%idx%nr, mesh%idx%enlarge(1))
    mesh%np_part_global = hypercube_number_total_points(mesh%idx%hypercube)
    mesh%np_global      = hypercube_number_inner_points(mesh%idx%hypercube)

    nullify(mesh%idx%Lxyz_inv)
    nullify(mesh%idx%Lxyz_tmp)
    nullify(mesh%idx%Lxyz)
    nullify(mesh%x_tmp)

    call profiling_out(mesh_init_prof)
    call pop_sub()
    return
  end if

  nr = mesh%idx%nr

  ! allocate the xyz arrays
  i = (nr(2, 1) - nr(1, 1) + 1)*(nr(2, 2) - nr(1, 2) + 1)*(nr(2, 3) - nr(1, 3) + 1)
  ALLOCATE(mesh%idx%Lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)), i)
  ALLOCATE(mesh%idx%Lxyz_tmp(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)), i)
  ALLOCATE(mesh%x_tmp(MAX_DIM, nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)), MAX_DIM*i)

  if(simul_box_multires(sb)) then 
    ALLOCATE(mesh%resolution(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)), i)
  else
    nullify(mesh%resolution)
  end if

  mesh%idx%Lxyz_inv(:,:,:) = 0
  mesh%idx%Lxyz_tmp(:,:,:) = 0
  mesh%x_tmp(:,:,:,:)  = M_ZERO

  ! We label the points inside the mesh + enlargement
  do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
    chi(1) = real(ix, REAL_PRECISION) * mesh%h(1) + sb%box_offset(1)

    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      chi(2) = real(iy, REAL_PRECISION) * mesh%h(2) + sb%box_offset(2)

      do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
        chi(3) = real(iz, REAL_PRECISION) * mesh%h(3) + sb%box_offset(3)

        call curvlinear_chi2x(sb, geo, cv, chi(:), mesh%x_tmp(:, ix, iy, iz))

        inside = simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz), inner_box = .true.)
        if(simul_box_multires(sb)) then
          mesh%resolution(ix, iy, iz) = 1
          if(.not. inside .and. simul_box_multires(sb)) mesh%resolution(ix, iy, iz) = 2
          inside = inside .or. (simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz)) .and. &
             mod(ix, 2) == 0 .and. mod(iy, 2) == 0 .and. mod(iz, 2) == 0)
          res = mesh%resolution(ix, iy, iz)
        else
          res = 1
        end if

        if(inside) then
          do is = 1, stencil%size
            i = ix + res*stencil%points(1, is)
            j = iy + res*stencil%points(2, is)
            k = iz + res*stencil%points(3, is)
            if(  &
                 i >= mesh%idx%nr(1,1) .and. i <= mesh%idx%nr(2,1) .and. &
                 j >= mesh%idx%nr(1,2) .and. j <= mesh%idx%nr(2,2) .and. &
                 k >= mesh%idx%nr(1,3) .and. k <= mesh%idx%nr(2,3)) then
              mesh%idx%Lxyz_tmp(i, j, k) = ENLARGEMENT_POINT
            end if
          end do
        end if

      end do
    end do
  end do

  ! we label the points inside the mesh, and we count the points
  il = 0
  ik = 0
  do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
    do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
      do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)

        inside = simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz), inner_box = .true.)
        inside = inside .or. (simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz)) .and. &
             mod(ix, 2) == 0 .and. mod(iy, 2) == 0 .and. mod(iz, 2) == 0)

        if(inside) then
          mesh%idx%Lxyz_tmp(ix, iy, iz) = INNER_POINT
          ik = ik + 1
        end if

        if(mesh%idx%Lxyz_tmp(ix, iy, iz) > 0) il = il + 1

      end do
    end do
  end do
  mesh%np_part_global = il
  mesh%np_global      = ik

  call profiling_out(mesh_init_prof)
  call pop_sub()
end subroutine mesh_init_stage_2


! ---------------------------------------------------------
! When running parallel in domains, stencil and np_stencil
! are needed to compute the ghost points.
! mpi_grp is the communicator group that will be used for
! this mesh.
! ---------------------------------------------------------
subroutine mesh_init_stage_3(mesh, geo, cv, stencil, mpi_grp, parent)
  type(mesh_t),              intent(inout) :: mesh
  type(geometry_t),          intent(in)    :: geo
  type(curvlinear_t),        intent(in)    :: cv
  type(stencil_t), optional, intent(in)    :: stencil
  type(mpi_grp_t), optional, intent(in)    :: mpi_grp
  type(mesh_t), optional,    intent(in)    :: parent

  integer :: ip, ix(1:MAX_DIM)
  FLOAT   :: chi(1:MAX_DIM)

  call push_sub('mesh_init.mesh_init_stage_3')
  call profiling_in(mesh_init_prof)

  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = .false.
  if(present(mpi_grp)) mesh%parallel_in_domains = .true.


  if(mesh%parallel_in_domains) then
    ! Node 0 has to store all entries from x (in x_global)
    ! as well as the local set in x (see below).
    ALLOCATE(mesh%x_global(mesh%np_part_global, MAX_DIM), mesh%np_part_global*MAX_DIM)
  else
    ! When running parallel, x is computed later.
    ALLOCATE(mesh%x(mesh%np_part_global, MAX_DIM), mesh%np_part_global*MAX_DIM)
    ! in the serial case x_global is the same as x
    mesh%x_global => mesh%x
  end if

  if(mesh%idx%sb%box_shape /= HYPERCUBE) then
    call create_x_Lxyz()
  else
    ! generate x_global directly
    do ip = 1, mesh%np_part_global
      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ix)
      chi(1:mesh%sb%dim) = ix(1:mesh%sb%dim)*mesh%h(1:mesh%sb%dim)
      chi(mesh%sb%dim + 1:MAX_DIM) = M_ZERO
      call curvlinear_chi2x(mesh%sb, geo, cv, chi, mesh%x_global(ip, :))
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
    mesh%vp%np = mesh%np_global
    mesh%vp%npart  = 1
  end if

  call mesh_get_vol_pp(mesh%sb)

  ! these large arrays were allocated in mesh_init_1, and are no longer needed
  DEALLOC(mesh%idx%Lxyz_tmp)
  DEALLOC(mesh%x_tmp)

  call mesh_pbc_init()

  call profiling_out(mesh_init_prof)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine create_x_Lxyz()
    integer :: il, ix, iy, iz
    integer :: ixb, iyb, izb, bsize, bsizez
    type(block_t) :: blk

    call push_sub('mesh_init.mesh_init_stage_3.create_x_Lxyz')

    ALLOCATE(mesh%idx%Lxyz(mesh%np_part_global, MAX_DIM), mesh%np_part_global*MAX_DIM)

    !%Type integer
    !%Default 20
    !%Section Execution::Optimization
    !%Description
    !% To improve memory access locality when calculating derivatives,
    !% Octopus orders mesh points in blocks. This variable controls
    !% the size of this blocks in the X and Y directions. The default
    !% is 20. (This variable only affects the performance of octopus
    !% and not the results.)
    !%End
    call loct_parse_int(datasets_check('MeshBlockSizeXY'), 20, bsize)

    !%Type integer
    !%Default 100
    !%Section Execution::Optimization
    !%Description
    !% To improve memory access locality when calculating derivatives,
    !% Octopus orders mesh points in blocks. This variable controls
    !% the size of this blocks in the Z direction. The default is
    !% 100. (This variable only affects the performance of octopus and
    !% not the results.)
    !%End
    call loct_parse_int(datasets_check('MeshBlockSizeZ'), 100, bsizez)

    ! When using open boundaries we need to have a mesh block-size of 1
    if (loct_parse_block(datasets_check('OpenBoundaries'), blk).eq.0) then
      if (bsize.gt.1 .or. bsizez.gt.1) then
        message(1) = 'When chosen the Transport-Mode the block-ordering'
        message(2) = 'of the mesh points cannot be chosen freely.'
        message(3) = 'Both variables, MeshBlockSizeXY and MeshBlockSizeZ'
        message(4) = 'have to be initialized with the value 1.'
        message(5) = 'Also restarting from datasets needs these to be'
        message(6) = 'calculated with the same block-size of 1.'
        call write_fatal(6)
      end if
    end if

    ! first we fill the points in the inner mesh
    il = 0
    do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize
      do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize
        do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsizez

          do ix = ixb, min(ixb + bsize - 1, mesh%idx%nr(2,1))
            do iy = iyb, min(iyb + bsize - 1, mesh%idx%nr(2,2))
              do iz = izb, min(izb + bsizez - 1, mesh%idx%nr(2,3))
                
                if(mesh%idx%Lxyz_tmp(ix, iy, iz) == INNER_POINT) then
                  il = il + 1
                  mesh%idx%Lxyz(il, 1) = ix
                  mesh%idx%Lxyz(il, 2) = iy
                  mesh%idx%Lxyz(il, 3) = iz
                  mesh%idx%Lxyz_inv(ix,iy,iz) = il
                  if(mesh%parallel_in_domains) then
                    mesh%x_global(il, 1:MAX_DIM) = mesh%x_tmp(1:MAX_DIM, ix, iy, iz)
                  else
                    mesh%x(il, 1:MAX_DIM) = mesh%x_tmp(1:MAX_DIM, ix, iy, iz)
                  end if
                end if
                
              end do
            end do
          end do
          
        end do
      end do
    end do

    ! and now the points from the enlargement
    do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize
      do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize
        do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsizez

          do ix = ixb, min(ixb + bsize - 1, mesh%idx%nr(2,1))
            do iy = iyb, min(iyb + bsize - 1, mesh%idx%nr(2,2))
              do iz = izb, min(izb + bsizez - 1, mesh%idx%nr(2,3))
                
                if(mesh%idx%Lxyz_tmp(ix, iy, iz) == ENLARGEMENT_POINT) then
                  il = il + 1
                  mesh%idx%Lxyz(il, 1) = ix
                  mesh%idx%Lxyz(il, 2) = iy
                  mesh%idx%Lxyz(il, 3) = iz
                  mesh%idx%Lxyz_inv(ix,iy,iz) = il
                  if(mesh%parallel_in_domains) then
                    mesh%x_global(il, :) = mesh%x_tmp(:, ix, iy, iz)
                  else
                    mesh%x(il,:) = mesh%x_tmp(:,ix,iy,iz)
                  end if
                end if
                
              end do
            end do
          end do

        end do
      end do
    end do
    
    call pop_sub()
  end subroutine create_x_Lxyz


  ! ---------------------------------------------------------
  subroutine do_partition()
#if defined(HAVE_METIS) && defined(HAVE_MPI)
    integer :: i, j, ipart, jpart, ip, ix, iy, iz
    integer, allocatable :: part(:), nnb(:)

    call push_sub('mesh_init.mesh_init_stage_3.do_partition')

    mesh%mpi_grp = mpi_grp

    ALLOCATE(part(mesh%np_part_global), mesh%np_part_global)

    if(.not. present(parent)) then
      call mesh_partition(mesh, stencil, part)
    else
      ! if there is a parent grid, use its partition
      do ip = 1, mesh%np_global
        ix = 2*mesh%idx%Lxyz(ip, 1)
        iy = 2*mesh%idx%Lxyz(ip, 2)
        iz = 2*mesh%idx%Lxyz(ip, 3)
        i = parent%idx%Lxyz_inv(ix, iy, iz)
        part(ip) = parent%vp%part(i)
      end do
    end if

    call mesh_partition_boundaries(mesh, stencil, part)

    call vec_init(mesh%mpi_grp%comm, 0, part, mesh%np_global, mesh%np_part_global, mesh%idx, stencil, mesh%sb%dim, mesh%vp)
    deallocate(part)

    ALLOCATE(nnb(1:mesh%vp%npart), mesh%vp%npart)
    nnb = 0
    do jpart = 1, mesh%vp%npart
      do ipart = 1, mesh%vp%npart
        if (ipart == jpart) cycle
        if (mesh%vp%np_ghost_neigh(jpart, ipart) /= 0) nnb(jpart) = nnb(jpart) + 1
      end do
      ASSERT(nnb(jpart) >= 0 .and. nnb(jpart) < mesh%vp%npart)
    end do

    ! Write information about partitions.
    message(1) = 'Info: Mesh partition:'
    message(2) = ''
    call write_info(2)

    write(message(1),'(a)') &
         '                 Neighbours         Ghost points'
    write(message(2),'(a,i5,a,i10)') &
         '      Average  :      ', sum(nnb)/mesh%vp%npart, '           ', sum(mesh%vp%np_ghost)/mesh%vp%npart
    write(message(3),'(a,i5,a,i10)') &
         '      Minimum  :      ', minval(nnb),        '           ', minval(mesh%vp%np_ghost)
    write(message(4),'(a,i5,a,i10)') &
         '      Maximum  :      ', maxval(nnb),        '           ', maxval(mesh%vp%np_ghost)
    message(5) = ''
    call write_info(5)

    do ipart = 1, mesh%vp%npart
      write(message(1),'(a,i5)')  &
           '      Nodes in domain-group  ', ipart
      write(message(2),'(a,i10,a,i10)') &
           '        Neighbours     :', nnb(ipart), &
           '        Local points    :', mesh%vp%np_local(ipart)
      write(message(3),'(a,i10,a,i10)') &
           '        Ghost points   :', mesh%vp%np_ghost(ipart), &
           '        Boundary points :', mesh%vp%np_bndry(ipart)
      call write_info(3)
    end do
    deallocate(nnb)

    message(1) = ''
    call write_info(1)

    ! Set local point numbers.
    mesh%np      = mesh%vp%np_local(mesh%vp%partno)
    mesh%np_part = mesh%np + mesh%vp%np_ghost(mesh%vp%partno) + mesh%vp%np_bndry(mesh%vp%partno)

    ! Compute mesh%x as it is done in the serial case but only for local points.
    ! x consists of three parts: the local points, the
    ! ghost points, and the boundary points; in this order
    ! (just as for any other vector, which is distributed).
    ALLOCATE(mesh%x(mesh%np_part, MAX_DIM), mesh%np_part*MAX_DIM)
    ! Do the inner points
    do i = 1, mesh%np
      j = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + i - 1)
      mesh%x(i, 1:MAX_DIM) = mesh%x_global(j, 1:MAX_DIM)
    end do
    ! Do the ghost points
    do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
      j = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno) + i - 1)
      mesh%x(i+mesh%np, 1:MAX_DIM) = mesh%x_global(j, 1:MAX_DIM)
    end do
    ! Do the boundary points
    do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
      j = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno) + i - 1)
      mesh%x(i + mesh%np + mesh%vp%np_ghost(mesh%vp%partno), 1:MAX_DIM) = mesh%x_global(j, 1:MAX_DIM)
    end do
#endif

    call pop_sub()
  end subroutine do_partition


  ! ---------------------------------------------------------
  ! calculate the volume of integration
  subroutine mesh_get_vol_pp(sb)
    type(simul_box_t), intent(in) :: sb

    integer :: i, jj(1:MAX_DIM), ip
    FLOAT   :: chi(MAX_DIM)
#if defined(HAVE_MPI)
    integer :: k
#endif

    call push_sub('mesh_init.mesh_init_stage_3.mesh_get_vol_pp')

    ALLOCATE(mesh%vol_pp(mesh%np_part), mesh%np_part)

    forall(ip = 1:mesh%np_part) mesh%vol_pp(ip) = product(mesh%h(1:sb%dim))
    jj(sb%dim + 1:MAX_DIM) = M_ZERO

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ! Do the inner points.
      do i = 1, mesh%np
        k = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + i - 1)
        call index_to_coords(mesh%idx, sb%dim, k, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%h(1:sb%dim)
        mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi(1:sb%dim))
      end do
      ! Do the ghost points.
      do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
        k = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno) + i - 1)
        call index_to_coords(mesh%idx, sb%dim, k, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%h(1:sb%dim)
        mesh%vol_pp(i + mesh%np) = mesh%vol_pp(i + mesh%np)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i + mesh%np, :), chi(1:sb%dim))
      end do
      ! Do the boundary points.
      do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
        k = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno) + i - 1)
        call index_to_coords(mesh%idx, sb%dim, k, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%h(1:sb%dim)
        mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) = &
          mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) &
          *curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno), :), chi(1:sb%dim))
      end do
#endif
    else ! serial mode
      do i = 1, mesh%np_part
        call index_to_coords(mesh%idx, sb%dim, i, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%h(1:sb%dim)

        mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, 1:sb%dim), chi(1:sb%dim))
        if(simul_box_multires(mesh%sb)) mesh%vol_pp(i) = mesh%resolution(jj(1), jj(2), jj(3))**sb%dim
      end do
    end if

    call pop_sub()

  end subroutine mesh_get_vol_pp
  
  subroutine mesh_pbc_init()
    integer :: sp, ip, ip_inner, iper, ip_global
#ifdef HAVE_MPI
    integer :: ip_inner_global, ipart, nblocks
    integer, allocatable :: recv_rem_points(:, :)
    integer :: nper_recv
    integer :: maxmax
    integer, allocatable :: recv_points(:, :), blocklengths(:), offsets(:)
    integer, allocatable :: send_points(:, :)
    integer, allocatable :: send_buffer(:)
    integer :: bsize, status(MPI_STATUS_SIZE)
#endif

    call push_sub('mesh_init.mesh_init_stage_3.pbc_init')

    nullify(mesh%per_points)
    nullify(mesh%per_map)

    if (simul_box_is_periodic(mesh%sb)) then

      sp = mesh%np
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) sp = mesh%np + mesh%vp%np_ghost(mesh%vp%partno)
#endif

      !count the number of points that are periodic
      mesh%nper = 0
#ifdef HAVE_MPI
      nper_recv = 0
#endif
      do ip = sp + 1, mesh%np_part

        ip_global = ip

#ifdef HAVE_MPI
        !translate to a global point
        if(mesh%parallel_in_domains) ip_global = mesh%vp%bndry(ip - sp - 1 + mesh%vp%xbndry(mesh%vp%partno))
#endif

        ip_inner = mesh_periodic_point(mesh, ip_global)

#ifdef HAVE_MPI
        !translate back to a local point
        if(mesh%parallel_in_domains) ip_inner = vec_global2local(mesh%vp, ip_inner, mesh%vp%partno)
#endif
        
        ! If the point is the periodic of another point, is not zero
        ! (this might happen in the parallel case) and is inside the
        ! grid then we have to copy it from the grid points.  
        !
        ! If the point index is larger than mesh%np then it is the
        ! periodic copy of a point that is zero, so we don`t count it
        ! as it will be initialized to zero anyway. For different
        ! mixed boundary conditions the last check should be removed.
        !
        if(ip /= ip_inner .and. ip_inner /= 0 .and. ip_inner <= mesh%np) then 
          mesh%nper = mesh%nper + 1
#ifdef HAVE_MPI
        else if(ip /= ip_inner) then
          nper_recv = nper_recv + 1
#endif
        end if
      end do

      ALLOCATE(mesh%per_points(1:mesh%nper), mesh%nper)
      ALLOCATE(mesh%per_map(1:mesh%nper), mesh%nper)

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then
        ALLOCATE(recv_points(1:nper_recv, 1:mesh%vp%npart), nper_recv*mesh%vp%npart)
        ALLOCATE(recv_rem_points(1:nper_recv, 1:mesh%vp%npart), nper_recv*mesh%vp%npart)
        ALLOCATE(mesh%nrecv(1:mesh%vp%npart), mesh%vp%npart)
        mesh%nrecv = 0
      end if
#endif

      iper = 0
      do ip = sp + 1, mesh%np_part

        ip_global = ip

        !translate to a global point
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) ip_global = mesh%vp%bndry(ip - sp - 1 + mesh%vp%xbndry(mesh%vp%partno))
#endif

        ip_inner = mesh_periodic_point(mesh, ip_global)
        
        !translate to local (and keep a copy of the global)
#ifdef HAVE_MPI
        if(mesh%parallel_in_domains) then
          ip_inner_global = ip_inner
          ip_inner = vec_global2local(mesh%vp, ip_inner, mesh%vp%partno)
        end if
#endif

        if(ip /= ip_inner .and. ip_inner /= 0 .and. ip_inner <= mesh%np) then
          iper = iper + 1
          mesh%per_points(iper) = ip
          mesh%per_map(iper) = ip_inner

#ifdef HAVE_MPI
        else if(ip /= ip_inner) then ! the point is in another node
          ! find in which paritition it is
          do ipart = 1, mesh%vp%npart
            if(ipart == mesh%vp%partno) cycle

            ip_inner = vec_global2local(mesh%vp, ip_inner_global, ipart)
            
            if(ip_inner /= 0 .and. ip_inner <= mesh%vp%np_local(ipart)) then
              ! count the points to receive from each node
              mesh%nrecv(ipart) = mesh%nrecv(ipart) + 1
              ! and store the number of the point
              recv_points(mesh%nrecv(ipart), ipart) = ip
              ! and where it is in the other partition
              recv_rem_points(mesh%nrecv(ipart), ipart) = ip_inner

              ASSERT(mesh%mpi_grp%rank /= ipart - 1) ! if we are here, the point must be in another node

              exit
            end if
            
          end do
#endif
        end if
      end do

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then

        ! first we allocate the buffer to be able to use MPI_Bsend
        bsize = mesh%vp%npart - 1 + nper_recv + MPI_BSEND_OVERHEAD*2*(mesh%vp%npart - 1)
        ALLOCATE(send_buffer(1:bsize), bsize)
        call MPI_Buffer_attach(send_buffer, bsize*4, mpi_err)

        ! Now we communicate to each node the points they will have to
        ! send us. Probably this could be done without communication,
        ! but this way it seems simpler to implement.
        
        ! We send the number of points we expect to receive.
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Bsend(mesh%nrecv(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, mpi_err)
        end do

        ! And we receive it
        ALLOCATE(mesh%nsend(1:mesh%vp%npart), mesh%vp%npart)
        mesh%nsend = 0
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Recv(mesh%nsend(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, status, mpi_err)
        end do

        ! Now we send the indexes of the points
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. mesh%nrecv(ipart) == 0) cycle
          call MPI_Bsend(recv_rem_points(:, ipart), mesh%nrecv(ipart), MPI_INTEGER, ipart - 1, 1, mesh%vp%comm, mpi_err)
        end do

        ALLOCATE(send_points(1:maxval(mesh%nsend), 1:mesh%vp%npart), maxval(mesh%nsend)*mesh%vp%npart)

        ! And we receive them
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. mesh%nsend(ipart) == 0) cycle
          call MPI_Recv(send_points(:, ipart), mesh%nsend(ipart), MPI_INTEGER, &
               ipart - 1, 1, mesh%vp%comm, status, mpi_err)
        end do

        ! we no longer need this
        deallocate(recv_rem_points)

        ! Now we have all the indexes required locally, so we can
        ! build the mpi datatypes

        ALLOCATE(mesh%dsend_type(1:mesh%vp%npart), mesh%vp%npart)
        ALLOCATE(mesh%zsend_type(1:mesh%vp%npart), mesh%vp%npart)
        ALLOCATE(mesh%drecv_type(1:mesh%vp%npart), mesh%vp%npart)
        ALLOCATE(mesh%zrecv_type(1:mesh%vp%npart), mesh%vp%npart)

        maxmax = max(maxval(mesh%nsend), maxval(mesh%nrecv))

        ALLOCATE(blocklengths(1:maxmax), maxmax)
        ALLOCATE(offsets(1:maxmax), maxmax)

        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle

          if(mesh%nsend(ipart) > 0) then

            ASSERT(all(send_points(1:mesh%nsend(ipart), ipart) <= mesh%np))

            ! MPI indexes start from zero
            send_points(1:mesh%nsend(ipart), ipart) = send_points(1:mesh%nsend(ipart), ipart) - 1

            call get_blocks(mesh%nsend(ipart), send_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_FLOAT, mesh%dsend_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_CMPLX, mesh%zsend_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%dsend_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%zsend_type(ipart), mpi_err)

          end if
          
          if(mesh%nrecv(ipart) > 0) then
            ASSERT(all(recv_points(1:mesh%nrecv(ipart), ipart) <= mesh%np_part))
            ASSERT(all(recv_points(1:mesh%nrecv(ipart), ipart) > mesh%np))

            recv_points(1:mesh%nrecv(ipart), ipart) = recv_points(1:mesh%nrecv(ipart), ipart) - 1

            call get_blocks(mesh%nrecv(ipart), recv_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_FLOAT, mesh%drecv_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths, offsets, MPI_CMPLX, mesh%zrecv_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%drecv_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%zrecv_type(ipart), mpi_err)

          end if

        end do

        call MPI_Buffer_detach(send_buffer, bsize, mpi_err)

      end if
#endif
      ASSERT(iper == mesh%nper)

    end if

    call pop_sub()
  end subroutine mesh_pbc_init

end subroutine mesh_init_stage_3


#if defined(HAVE_METIS) && defined(HAVE_MPI)
! ---------------------------------------------------------------
! Converts the mesh given by grid points into a graph. Each
! point is a vertex in the graph and closest neighbours are
! connected by an edge (at max. 6 in 3D and 4 in 2D, 2 in
! 1D, less at the boundaries).
! Then calls METIS to get p partitions.
! Stored the mapping point no. -> partition no. into part,
! which has to be allocated beforehand.
! (mesh_partition_end should be called later.)
! In Lxyz_tmp has to be stored which points belong to the
! inner mesh and the enlargement. All other entries have to
! be zero. comm is used to get the number of partitions.
! ---------------------------------------------------------------
subroutine mesh_partition(m, lapl_stencil, part)
  type(mesh_t),            intent(in)  :: m
  type(stencil_t), target, intent(in)  :: lapl_stencil
  integer,                 intent(out) :: part(:)

  integer              :: i, j, inb
  integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
  integer              :: ne             ! Number of edges.
  integer              :: nv             ! Number of vertices.
  integer              :: edgecut        ! Number of edges cut by partitioning.
  ! Number of vertices (nv) is equal to number of
  ! points np_global and maximum number of edges (ne) is 2*m%sb%dim*np_global
  ! (there are a little less because points on the border have less
  ! than two neighbours per dimension).
  ! xadj has nv+1 entries because last entry contains the total
  ! number of edges.
  integer              :: p              ! Number of partitions.
  integer              :: ipart          ! number of the current partition
  integer, allocatable :: xadj(:)        ! Indices of adjacency list in adjncy.
  integer, allocatable :: adjncy(:)      ! Adjacency lists.
  integer              :: options(5)     ! Options to METIS.
  integer              :: iunit          ! For debug output to files.

  type(stencil_t) :: stencil
  integer :: ii, stencil_to_use
  integer :: default_method, method
  integer :: library
  integer, parameter :: METIS = 2, ZOLTAN = 3
  integer, parameter :: STAR = 1, LAPLACIAN = 2
  integer, allocatable :: start(:), final(:), lsize(:)
  
  type(profile_t), save :: prof
  
  call profiling_in(prof, "MESH_PARTITION")
  call push_sub('mesh_init.mesh_partition')

  if(m%np_global == 0) then
    message(1) = 'The mesh is empty and cannot be partitioned.'
    call write_fatal(1)
  end if

  !%Variable MeshPartitionPackage
  !%Type integer
  !%Default zoltan
  !%Section Execution::Parallelization
  !%Description
  !% Decides which library to use to perform the mesh partition. By
  !% default metis is used.
  !%Option metis 2
  !% Metis.
  !%Option zoltan 3
  !% Zoltan.
  !%End
  call loct_parse_int(datasets_check('MeshPartitionPackage'), METIS, library)

  !%Variable MeshPartitionStencil
  !%Type integer
  !%Default star
  !%Section Execution::Parallelization
  !%Description
  !% To partition the mesh is it necessary to calculate the connection
  !% graph connecting the points, this variable selects which stencil
  !% is used to do this. The default is the order one star stencil,
  !% alternatively the stencil used for the laplacian could be used.
  !%Option stencil_star 1
  !% A order one stencil_star.
  !%Option laplacian 2
  !% The stencil used for the laplacian is used to calculate the
  !% partition, this in principle should give a better partition but
  !% it is slower and requires more memory.
  !%End
  call loct_parse_int(datasets_check('MeshPartitionStencil'), STAR, stencil_to_use)

  if (stencil_to_use == STAR) then
    call stencil_star_get_lapl(stencil, m%sb%dim, order = 1)
  else if (stencil_to_use == LAPLACIAN) then
    call stencil_copy(lapl_stencil, stencil)
  else
    call input_error('MeshPartitionStencil')
  end if

  ! Get number of partitions.
  p = m%mpi_grp%size
  ipart = m%mpi_grp%rank + 1

  if(p .lt. 8) then
    default_method = RCB
  else
    default_method = GRAPH
  end if
  ! Documentation is in zoltan.F90`
  call loct_parse_int(datasets_check('MeshPartition'), default_method, method)

  ALLOCATE(start(1:p), p)
  ALLOCATE(final(1:p), p)
  ALLOCATE(lsize(1:p), p)

  select case(library)
  case(METIS)

    start(1:p) = 1
    final(1:p) = m%np_global
    lsize(1:p) = m%np_global

  case(ZOLTAN)

    ! If we use zoltan we divide the space in a basic way, to balance
    ! the memory for the graph. 
    call multicomm_divide_range(m%np_global, p, start, final, lsize)

    do ii = 1, p
      part(start(ii):final(ii)) = ii
    end do

  end select

  ! Shortcut (number of vertices).
  nv = lsize(ipart)
  ALLOCATE(xadj(nv + 1), nv + 1)

  if(library == METIS .or. .not. zoltan_method_is_geometric(method)) then !calculate the graphs
    ALLOCATE(adjncy(stencil%size*nv), stencil%size*nv)

    ! Create graph with each point being
    ! represenetd by a vertice and edges between
    ! neighboured points.
    ne = 1
    ! Iterate over number of vertices.
    do i = 1, nv
      ! Get coordinates of point i (vertex i).
      call index_to_coords(m%idx, m%sb%dim, i, ix)
      ! Set entry in index table.
      xadj(i) = ne
      ! Check all possible neighbours.
      do j = 1, stencil%size 
        ! Store coordinates of possible neighbors, they
        ! are needed several times in the check below.
        jx(1:MAX_DIM) = ix(1:MAX_DIM) + stencil%points(1:MAX_DIM, j)

        ! Only if the neighbour is in the surrounding box,
        ! Lxyz_tmp has an entry for this point, otherweise
        ! it is out of bounds.
        if(all(jx(1:MAX_DIM) >= m%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= m%idx%nr(2, 1:MAX_DIM))) then
          ! Only points inside the mesh or its enlargement
          ! are included in the graph.
          inb = index_from_coords(m%idx, m%sb%dim, jx)
          if(inb /= 0 .and. inb <= nv) then
            ! Store a new edge and increment edge counter.
            adjncy(ne) = inb
            ne         = ne + 1
          end if
        end if
      end do
    end do
    ne         = ne - 1 ! We start with ne=1 for simplicity. This is off by one
    ! in the end --> -1.
    xadj(nv + 1) = ne + 1 ! Set number of edges plus 1 as last index.
    ! The reason is: neighbours of node i are stored
    ! in adjncy(xadj(i):xadj(i+1)-1). Setting the last
    ! index as mentioned makes special handling of
    ! last element unnecessary (this indicing is a
    ! METIS requirement).

    if(in_debug_mode) then
      ! DEBUG output. Write graph to file mesh_graph.txt.
      message(1) = 'Info: Adjacency lists of the graph representing the grid'
      message(2) = 'Info: are stored in debug/mesh_partition/mesh_graph.txt.'
      message(3) = 'Info: Compatible with METIS programs pmetis and kmetis.'
      message(4) = 'Info: First line contains number of vertices and edges.'
      message(5) = 'Info: Edges are not directed and appear twice in the lists.'
      call write_info(5)
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/mesh_partition')
        iunit = io_open('debug/mesh_partition/mesh_graph.txt', action='write')
        write(iunit, *) nv, ne/2
        do i = 1, nv
          write(iunit, *) adjncy(xadj(i):xadj(i+1) - 1)
        end do
        call io_close(iunit)
      end if
    end if

  end if

  select case(library)
  case(METIS)

    options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

    ! Partition graph.
    ! Recursive bisection is better for small number of partitions (<8),
    ! multilevel k-way otherwise (cf. METIS manual).
    ! If the graph contains no vertices, METIS cannot be called. This seems
    ! to happen, e. g., when using minimum BoxShape without any atoms in the
    ! input file.

    select case(method)
    case(RCB)
      message(1) = 'Info: Using Metis multilevel recursive bisection to partition the mesh.'
      call write_info(1)
      call oct_metis_part_graph_recursive(nv, xadj, adjncy, &
           0, 0, 0, 1, p, options, edgecut, part)
    case(GRAPH)
      message(1) = 'Info: Using Metis multilevel k-way algorithm to partition the mesh.'
      call write_info(1)
      call oct_metis_part_graph_kway(nv, xadj, adjncy, &
           0, 0, 0, 1, p, options, edgecut, part)
    case default
      message(1) = 'Error: Selected partition method is not available in Metis'
      call write_fatal(1)
    end select

  case(ZOLTAN)

    call zoltan_method_info(method)

    !assign all points to one node
    call zoltan_partition(method, m%sb%dim, m%np_global, m%np_part_global, &
         m%x_global(1, 1),  start(ipart), xadj(1), adjncy(1), ipart, part(1), m%mpi_grp%comm)

    ! we use xadj as a buffer
    xadj(1:lsize(ipart)) = part(start(ipart):final(ipart))

    ASSERT(all(xadj(1:lsize(ipart)) > 0))

    part(1:m%np_global) = 0 ! so we catch non-initialized values

    ! we collect part from all processors
    call MPI_Allgatherv(xadj, lsize(ipart), MPI_INTEGER, part, lsize, start - 1, MPI_INTEGER, m%mpi_grp%comm, mpi_err)

  end select

  ASSERT(all(part(1:m%np_global) > 0))
  ASSERT(all(part(1:m%np_global) <= p))

  call pop_sub()
  call profiling_out(prof)

end subroutine mesh_partition

subroutine mesh_partition_boundaries(mesh, stencil, part)
  type(mesh_t),    intent(in)    :: mesh
  type(stencil_t), intent(in)    :: stencil
  integer,         intent(inout) :: part(:)

  integer              :: i, j           ! Counter.
  integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
  integer              :: npart
  integer              :: iunit          ! For debug output to files.
  character(len=3)     :: filenum
  integer, allocatable :: votes(:, :)
  integer :: ip, rr

  call push_sub('mesh_init.mesh_partition_boundaries')

  npart = mesh%mpi_grp%size

  ALLOCATE(votes(1:npart, mesh%np_global + 1:mesh%np_part_global), npart*(mesh%np_part_global - mesh%np_global))

  !now assign boundary points

  !count the boundary points that each point needs
  votes = 0
  do i = 1, mesh%np_global
    call index_to_coords(mesh%idx, mesh%sb%dim, i, ix)
    do j = 1, stencil%size
      jx(1:MAX_DIM) = ix(1:MAX_DIM) + stencil%points(1:MAX_DIM, j)
      ip = index_from_coords(mesh%idx, mesh%sb%dim, jx)
      if(ip > mesh%np_global) votes(part(i), ip) = votes(part(i), ip) + 1
    end do
  end do

  rr = 1
  do i = mesh%np_global + 1, mesh%np_part_global
    if(maxval(votes(1:npart, i)) /= minval(votes(1:npart, i))) then
      ! we have a winner that takes the point
      part(i:i) = maxloc(votes(1:npart, i))
    else
      ! points without a winner are assigned in a round robin fashion
      part(i) = rr
      rr = rr + 1
      if(rr > npart) rr = 1
    end if
  end do

  if(in_debug_mode.and.mpi_grp_is_root(mpi_world)) then
    ! Debug output. Write points of each partition in a different file.
    do i = 1, npart

      write(filenum, '(i3.3)') i

      iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
           action='write')
      do j = 1, mesh%np_global
        if(part(j).eq.i) write(iunit, '(i8,3f18.8)') j, mesh%x_global(j, :)
      end do
      call io_close(iunit)

      iunit = io_open('debug/mesh_partition/mesh_partition_all.'//filenum, &
           action='write')
      do j = 1, mesh%np_part_global
        if(part(j).eq.i) write(iunit, '(i8,3f18.8)') j, mesh%x_global(j, :)
      end do
      call io_close(iunit)

    end do
    ! Write points from enlargement to file with number p+1.
    write(filenum, '(i3.3)') npart+1
    iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
         action='write')
    do i = mesh%np_global+1, mesh%np_part_global
      write(iunit, '(i8,3f18.8)') i, mesh%x_global(i, :)
    end do
    call io_close(iunit)
  end if

  call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)

  call pop_sub()

end subroutine mesh_partition_boundaries

#endif
end module mesh_init_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
