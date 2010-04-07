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

module mesh_init_m
  use curvilinear_m
  use datasets_m
  use geometry_m
  use global_m
  use hypercube_m
  use index_m
  use io_m
  use loct_m
  use ob_grid_m
  use parser_m
  use math_m
  use mesh_m
  use messages_m
  use mpi_m
  use multicomm_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use stencil_star_m
  use subarray_m
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

  call push_sub('mesh_init.mesh_init_stage_1')
  call profiling_in(mesh_init_prof, "MESH_INIT")

  mesh%sb => sb     ! keep an internal pointer
  mesh%idx%sb => sb
  mesh%spacing = spacing ! this number can change in the following
  mesh%use_curvilinear = cv%method.ne.CURV_METHOD_UNIFORM
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

    call mesh_read_lead()
  else
    nullify(mesh%lead_unit_cell)
  end if

  mesh%idx%ll(:) = mesh%idx%nr(2, :) - mesh%idx%nr(1, :) + 1

  call profiling_out(mesh_init_prof)
  call pop_sub('mesh_init.mesh_init_stage_1')

contains

  ! ---------------------------------------------------------
  subroutine mesh_read_lead()
    integer :: il, iunit

    type(mesh_t), pointer :: mm
    integer :: nr(1:2, 1:MAX_DIM)

    call push_sub('mesh_init.mesh_init_stage1_mesh_read_lead')

    SAFE_ALLOCATE(mesh%lead_unit_cell(1:NLEADS))

    do il = 1, NLEADS
      mm => mesh%lead_unit_cell(il)
      mm%sb => mesh%sb
      mm%cv => mesh%cv
      mm%parallel_in_domains = mesh%parallel_in_domains
      iunit = io_open(trim(ob_grid%lead(il)%info%restart_dir)//'/'//GS_DIR//'mesh', action='read', is_tmp=.true.)
      call mesh_init_from_file(mm, iunit)
      call io_close(iunit)
      ! Read the lxyz maps.
      nr = mm%idx%nr
      SAFE_ALLOCATE(mm%idx%Lxyz(1:mm%np_part, 1:3))
      SAFE_ALLOCATE(mm%idx%Lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))
      call mesh_lxyz_init_from_file(mm, trim(ob_grid%lead(il)%info%restart_dir)//'/'//GS_DIR//'lxyz')
    end do

    call pop_sub('mesh_init.mesh_init_stage1_mesh_read_lead')
  end subroutine mesh_read_lead
end subroutine mesh_init_stage_1


! ---------------------------------------------------------
subroutine mesh_init_stage_2(mesh, sb, geo, cv, stencil)
  type(mesh_t),       intent(inout) :: mesh
  type(simul_box_t),  intent(in)    :: sb
  type(geometry_t),   intent(in)    :: geo
  type(curvilinear_t), intent(in)    :: cv
  type(stencil_t),    intent(in)    :: stencil

  integer :: i, j, k, il, ik, ix, iy, iz, is
  integer :: newi, newj, newk, ii, jj, kk, dx, dy, dz, i_lev
  integer :: jx, jy, jz, res_counter, j_counter
  FLOAT   :: chi(MAX_DIM)
  integer :: nr(1:2, 1:MAX_DIM), res, n_mod
  logical, allocatable :: in_box(:)
  FLOAT,   allocatable :: xx(:, :)
  real(8), parameter :: DELTA = CNST(1e-12)

  call push_sub('mesh_init.mesh_init_stage_2')
  call profiling_in(mesh_init_prof)

  ! enlarge mesh for boundary points
  mesh%idx%nr(1,:) = mesh%idx%nr(1,:) - mesh%idx%enlarge(:)
  mesh%idx%nr(2,:) = mesh%idx%nr(2,:) + mesh%idx%enlarge(:)

  if(mesh%idx%sb%box_shape == HYPERCUBE) then
    call hypercube_init(mesh%idx%hypercube, sb%dim, mesh%idx%nr, mesh%idx%enlarge(1))
    mesh%np_part_global = hypercube_number_total_points(mesh%idx%hypercube)
    mesh%np_global      = hypercube_number_inner_points(mesh%idx%hypercube)

    nullify(mesh%resolution)
    nullify(mesh%idx%Lxyz_inv)
    nullify(mesh%idx%Lxyz)

    go to 999
  end if

  nr = mesh%idx%nr

  ! allocate the xyz arrays
  SAFE_ALLOCATE(mesh%idx%Lxyz_inv(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))

  if(sb%mr_flag) then 
    SAFE_ALLOCATE(mesh%resolution(nr(1, 1):nr(2, 1), nr(1, 2):nr(2, 2), nr(1, 3):nr(2, 3)))
    mesh%resolution(:,:,:) = 0
  else
    nullify(mesh%resolution)
  end if

  mesh%idx%Lxyz_inv(:,:,:) = 0
  res = 1

  SAFE_ALLOCATE(xx(1:MAX_DIM, mesh%idx%nr(1,1):mesh%idx%nr(2,1)))
  SAFE_ALLOCATE(in_box(mesh%idx%nr(1,1):mesh%idx%nr(2,1)))

  chi = M_ZERO

  ! We label the points inside the mesh
  do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
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
            if (sum((xx(:,ix)-sb%hr_area%center(:))**2).gt. sb%hr_area%radius(sb%hr_area%num_radii)**2 .and. &
                 mod(ix, n_mod).eq.0 .and. mod(iy, n_mod).eq.0 .and. mod(iz,n_mod) .eq. 0) then
              mesh%idx%Lxyz_inv(ix, iy, iz) = ibset(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)
            end if

            ! Other option: must be inside the multiresolution area and satisfy coordinate index conditions
            if(.not.btest(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)) then
              do i_lev = 1,sb%hr_area%num_radii
                n_mod = 2**(i_lev-1)
                if( sum((xx(:,ix)-sb%hr_area%center(:))**2) .lt. sb%hr_area%radius(i_lev)**2 + DELTA .and. &
                    mod(ix, n_mod).eq.0 .and. mod(iy, n_mod).eq.0 .and. mod(iz,n_mod) .eq. 0) then
                  mesh%idx%Lxyz_inv(ix, iy, iz) = ibset(mesh%idx%Lxyz_inv(ix,iy, iz), INNER_POINT)
                end if
              end do
            end if

          end if

        else ! the usual way: mark both inner and enlargement points

          if (in_box(ix)) then

            mesh%idx%Lxyz_inv(ix, iy, iz) = ibset(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)

            do is = 1, stencil%size
              if(stencil%center == is) cycle
    
              i = ix + stencil%points(1, is)
              j = iy + stencil%points(2, is)
              k = iz + stencil%points(3, is)
    
              if(any((/i, j, k/) < mesh%idx%nr(1, 1:3)) .or. any((/i, j, k/) >  mesh%idx%nr(2, 1:3))) cycle
    
              mesh%idx%Lxyz_inv(i, j, k) = ibset(mesh%idx%Lxyz_inv(i, j, k), ENLARGEMENT_POINT)

            end do

          end if
  
        end if

      end do
    end do
  end do

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
          if(.not.btest(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)) cycle
        
          res = -1
          res_counter = 0
 
          do while(.true.)
 
            res_counter = res_counter + 1
 
            ! loop through Cartesian axes (both directions)
            do j_counter = 1,6
              select case(j_counter)
                case(1)
                  jx=ix-res_counter; jy=iy; jz=iz;
                case(2)
                  jx=ix+res_counter; jy=iy; jz=iz;
                case(3)
                  jx=ix; jy=iy-res_counter; jz=iz;
                case(4)
                  jx=ix; jy=iy+res_counter; jz=iz;
                case(5)
                  jx=ix; jy=iy; jz=iz-res_counter;
                case(6)
                  jx=ix; jy=iy; jz=iz+res_counter;
              end select
 
              if(any((/jx, jy, jz/) < mesh%idx%nr(1, 1:3)) .or. &
                 any((/jx, jy, jz/) > mesh%idx%nr(2, 1:3))) cycle
              ! exit after finding neighboring inner point
              if(btest(mesh%idx%Lxyz_inv(jx, jy, jz), INNER_POINT)) then
                res = res_counter
                exit
              end if
 
            end do
 
            if(res.ne.-1) exit
 
          end do
 
          mesh%resolution(ix, iy, iz) = res

          ! mark the enlargement points
          do is = 1, stencil%size
            if(stencil%center == is) cycle
 
            i = ix + res*stencil%points(1, is)
            j = iy + res*stencil%points(2, is)
            k = iz + res*stencil%points(3, is)
 
            if(any((/i, j, k/) < mesh%idx%nr(1, 1:3)) .or. any((/i, j, k/) >  mesh%idx%nr(2, 1:3))) cycle

            mesh%idx%Lxyz_inv(i, j, k) = ibset(mesh%idx%Lxyz_inv(i, j, k), ENLARGEMENT_POINT)

            ! If the point is not an inner point, and if its resolution can be decreased, do it now
            if(.not.btest(mesh%idx%Lxyz_inv(i, j, k), INNER_POINT)) then
              if(mesh%resolution(i,j,k).eq.0 .or. mesh%resolution(i,j,k).gt.res) mesh%resolution(i,j,k) = res
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
        if(btest(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)) ik = ik + 1
        if(mesh%idx%Lxyz_inv(ix, iy, iz) /= 0) il = il + 1
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

          ! include enlargement points that are not inner points nor outer boundary points.
          if( .not. btest(mesh%idx%Lxyz_inv(ix, iy, iz),ENLARGEMENT_POINT)) cycle
          if(  btest(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)) cycle
          if(  i_lev.eq.2**mesh%sb%hr_area%num_radii ) cycle

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
                   mesh%idx%Lxyz_inv(newi,newj,newk).eq.0) then

                     message(1) = 'Multiresolution radii are too close to each other (or outer boundary)'
                     write(message(2),'(7I4)') ix,iy,iz,newi,newj,newk,mesh%resolution(ix,iy,iz)
                     call write_fatal(2)
                end if
              end do
            end do
          end do
 
        end do
      end do
    end do

  end if



999 continue

  call profiling_out(mesh_init_prof)
  call pop_sub('mesh_init.mesh_init_stage_2')
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
  type(mesh_t), optional,    intent(in)    :: parent

  integer :: ip

  call push_sub('mesh_init.mesh_init_stage_3')
  call profiling_in(mesh_init_prof)

  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = .false.
  if(present(mpi_grp)) mesh%parallel_in_domains = .true.


  if(.not. mesh%parallel_in_domains) then
    ! When running parallel, x is computed later.
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part_global, 1:MAX_DIM))
  end if
  
  if(mesh%idx%sb%box_shape /= HYPERCUBE) then
    call create_x_Lxyz()
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
    mesh%vp%np = mesh%np_global
    mesh%vp%npart  = 1
  end if

  call mesh_get_vol_pp(mesh%sb)

  call mesh_pbc_init()

  call profiling_out(mesh_init_prof)
  call pop_sub('mesh_init.mesh_init_stage_3')

contains

  ! ---------------------------------------------------------
  subroutine create_x_Lxyz()
    integer :: il, iin, ien, ix, iy, iz
    integer :: ixb, iyb, izb, bsize, bsizez
    type(block_t) :: blk
    FLOAT :: chi(1:MAX_DIM)

    call push_sub('mesh_init.mesh_init_stage_3.create_x_Lxyz')

    SAFE_ALLOCATE(mesh%idx%Lxyz(1:mesh%np_part_global, 1:MAX_DIM))

    !%Variable MeshBlockSizeXY
    !%Type integer
    !%Default 20
    !%Section Execution::Optimization
    !%Description
    !% To improve memory-access locality when calculating derivatives,
    !% <tt>Octopus</tt> arranges mesh points in blocks. This variable controls
    !% the size of this blocks in the <i>x</i>- and <i>y</i>-directions. The default
    !% is 20. (This variable only affects the performance of <tt>Octopus</tt>
    !% and not the results.)
    !%End
    call parse_integer(datasets_check('MeshBlockSizeXY'), 20, bsize)

    !%Variable MeshBlockSizeZ
    !%Type integer
    !%Default 100
    !%Section Execution::Optimization
    !%Description
    !% To improve memory-access locality when calculating derivatives,
    !% <tt>Octopus</tt> arranges mesh points in blocks. This variable controls
    !% the size of this blocks in the <i>z</i>-direction. The default is
    !% 100. (This variable only affects the performance of <tt>Octopus</tt> and
    !% not the results.)
    !%End
    call parse_integer(datasets_check('MeshBlockSizeZ'), 100, bsizez)

    ! When using open boundaries we need to have a mesh block-size of 1
    if (parse_block(datasets_check('OpenBoundaries'), blk).eq.0) then
      if (bsize.gt.1 .or. bsizez.gt.1) then
        message(1) = 'When in transport mode, the block-ordering'
        message(2) = 'of the mesh points cannot be chosen freely.'
        message(3) = 'Both variables, MeshBlockSizeXY and MeshBlockSizeZ'
        message(4) = 'have to be initialized with the value 1.'
        message(5) = 'Also, when restarting from datasets, these need to be'
        message(6) = 'calculated with the same block-size of 1.'
        call write_fatal(6)
      end if
    end if

    ! first we fill the points in the inner mesh
    iin = 0
    ien = mesh%np_global
    do ixb = mesh%idx%nr(1,1), mesh%idx%nr(2,1), bsize
      do iyb = mesh%idx%nr(1,2), mesh%idx%nr(2,2), bsize
        do izb = mesh%idx%nr(1,3), mesh%idx%nr(2,3), bsizez

          do ix = ixb, min(ixb + bsize - 1, mesh%idx%nr(2,1))
            chi(1) = real(ix, REAL_PRECISION) * mesh%spacing(1) + mesh%sb%box_offset(1)
            do iy = iyb, min(iyb + bsize - 1, mesh%idx%nr(2,2))
              chi(2) = real(iy, REAL_PRECISION) * mesh%spacing(2) + mesh%sb%box_offset(2)
              do iz = izb, min(izb + bsizez - 1, mesh%idx%nr(2,3))
                chi(3) = real(iz, REAL_PRECISION) * mesh%spacing(3) + mesh%sb%box_offset(3)
                
                if(btest(mesh%idx%Lxyz_inv(ix, iy, iz), INNER_POINT)) then
                  iin = iin + 1
                  il = iin
                else if (btest(mesh%idx%Lxyz_inv(ix, iy, iz), ENLARGEMENT_POINT)) then
                  ien = ien + 1
                  il = ien
                else
                  cycle
                end if
                  
                mesh%idx%Lxyz(il, 1) = ix
                mesh%idx%Lxyz(il, 2) = iy
                mesh%idx%Lxyz(il, 3) = iz

                mesh%idx%Lxyz_inv(ix, iy, iz) = il

#ifdef HAVE_MPI
                if(.not. mesh%parallel_in_domains) &
#endif                  
                  call curvilinear_chi2x(mesh%sb, mesh%cv, chi(:), mesh%x(il, 1:MAX_DIM))
                
              end do
            end do
          end do
          
        end do
      end do
    end do

    call checksum_calculate(1, mesh%np_part_global*MAX_DIM, mesh%idx%lxyz(1, 1), mesh%idx%checksum)

    ASSERT(iin == mesh%np_global)
    ASSERT(ien == mesh%np_part_global)
    
    call pop_sub('mesh_init.mesh_init_stage_3.create_x_Lxyz')
  end subroutine create_x_Lxyz


  ! ---------------------------------------------------------
  subroutine do_partition()
#ifdef HAVE_MPI
    integer :: i, j, ipart, jpart, ip, ix, iy, iz
    integer, allocatable :: part(:), nnb(:), gindex(:), gedges(:)
    logical, allocatable :: nb(:, :)
    integer              :: idx(1:MAX_DIM), jx(1:MAX_DIM)
    integer              :: graph_comm, iedge, reorder
    logical              :: use_topo

    call push_sub('mesh_init.mesh_init_stage_3.do_partition')

    mesh%mpi_grp = mpi_grp

    SAFE_ALLOCATE(part(1:mesh%np_part_global))

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
        ipart = part(ip)
        call index_to_coords(mesh%idx, mesh%sb%dim, ip, idx)
        do j = 1, stencil%size
          jx(1:MAX_DIM) = idx(1:MAX_DIM) + stencil%points(1:MAX_DIM, j)
          if(all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
            jpart = part(index_from_coords(mesh%idx, mesh%sb%dim, jx))
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

      reorder = 1
      call MPI_Graph_create(mpi_grp%comm, mpi_grp%size, gindex(1), gedges(1), reorder, graph_comm, mpi_err)

      SAFE_DEALLOCATE_A(nb)
      SAFE_DEALLOCATE_A(gindex)
      SAFE_DEALLOCATE_A(gedges)

    else

      call MPI_Comm_dup(mpi_grp%comm, graph_comm, mpi_err)

    end if

    ! we have a new communicator
    call mpi_grp_init(mesh%mpi_grp, graph_comm)

    call vec_init(graph_comm, 0, part, mesh%np_global, mesh%np_part_global, mesh%idx, stencil, mesh%sb%dim, mesh%vp)
    SAFE_DEALLOCATE_A(part)

    SAFE_ALLOCATE(nnb(1:mesh%vp%npart))
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
    SAFE_DEALLOCATE_A(nnb)

    message(1) = ''
    call write_info(1)

    ! Set local point numbers.
    mesh%np      = mesh%vp%np_local(mesh%vp%partno)
    mesh%np_part = mesh%np + mesh%vp%np_ghost(mesh%vp%partno) + mesh%vp%np_bndry(mesh%vp%partno)

    ! Compute mesh%x as it is done in the serial case but only for local points.
    ! x consists of three parts: the local points, the
    ! ghost points, and the boundary points; in this order
    ! (just as for any other vector, which is distributed).
    SAFE_ALLOCATE(mesh%x(1:mesh%np_part, 1:MAX_DIM))
    mesh%x(:, :) = M_ZERO
    ! Do the inner points
    do i = 1, mesh%np
      j = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + i - 1)
      mesh%x(i, 1:MAX_DIM) = mesh_x_global(mesh, j)
    end do
    ! Do the ghost points
    do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
      j = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno) + i - 1)
      mesh%x(i+mesh%np, 1:MAX_DIM) = mesh_x_global(mesh, j)
    end do
    ! Do the boundary points
    do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
      j = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno) + i - 1)
      mesh%x(i + mesh%np + mesh%vp%np_ghost(mesh%vp%partno), 1:MAX_DIM) = mesh_x_global(mesh, j)
    end do
#endif

    call pop_sub('mesh_init.mesh_init_stage_3.do_partition')
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
    integer :: k
#endif

    call push_sub('mesh_init.mesh_init_stage_3.mesh_get_vol_pp')

    np = 1
    if(mesh%use_curvilinear) np = mesh%np_part

    SAFE_ALLOCATE(mesh%vol_pp(1:np))

    forall(ip = 1:np) mesh%vol_pp(ip) = product(mesh%spacing(1:sb%dim))
    jj(sb%dim + 1:MAX_DIM) = M_ZERO

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ! Do the inner points.
      do ip = 1, min(np, mesh%np)
        k = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
        call index_to_coords(mesh%idx, sb%dim, k, jj)
        chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
        mesh%vol_pp(ip) = mesh%vol_pp(ip)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip, :), chi(1:sb%dim))
      end do

      if(mesh%use_curvilinear) then
        ! Do the ghost points.
        do ip = 1, mesh%vp%np_ghost(mesh%vp%partno)
          k = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno) + ip - 1)
          call index_to_coords(mesh%idx, sb%dim, k, jj)
          chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
          mesh%vol_pp(ip + mesh%np) = &
            mesh%vol_pp(ip + mesh%np)*curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip + mesh%np, :), chi(1:sb%dim))
        end do
        ! Do the boundary points.
        do ip = 1, mesh%vp%np_bndry(mesh%vp%partno)
          k = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno) + ip - 1)
          call index_to_coords(mesh%idx, sb%dim, k, jj)
          chi(1:sb%dim) = jj(1:sb%dim)*mesh%spacing(1:sb%dim)
          mesh%vol_pp(ip+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) = &
            mesh%vol_pp(ip+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) &
            *curvilinear_det_Jac(sb, mesh%cv, mesh%x(ip+mesh%np+mesh%vp%np_ghost(mesh%vp%partno), :), chi(1:sb%dim))
        end do
      end if
#endif
    else ! serial mode

      if(mesh%sb%mr_flag) then

        message(1) = 'Info: Point volumes are calculated by solving interpolation coefficients for the intermediate points.'
        call write_info(1)

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
          call write_info(1)

          ! loop through _all_ the points
          do iz = mesh%idx%nr(1,3), mesh%idx%nr(2,3)
            do iy = mesh%idx%nr(1,2), mesh%idx%nr(2,2)
              do ix = mesh%idx%nr(1,1), mesh%idx%nr(2,1)
 
                ! Skip ordinary points
                if(mesh%idx%Lxyz_inv(ix,iy,iz).gt.0 .and. &
                     mesh%idx%Lxyz_inv(ix,iy,iz).le.mesh%np) cycle

                ! Is it the kind of intermediate point we are looking for?
                n_mod = 2**i_lev
                dx = abs(mod(ix, n_mod))
                dy = abs(mod(iy, n_mod))
                dz = abs(mod(iz, n_mod))
                if(dx+dy+dz.eq.M_ZERO) cycle

                if(abs(vol_tmp(ix, iy, iz)).lt.DELTA) cycle

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
          ix = mesh%idx%Lxyz(ip, 1)
          iy = mesh%idx%Lxyz(ip, 2)
          iz = mesh%idx%Lxyz(ip, 3)
          mesh%vol_pp(ip) = vol_tmp(ix,iy,iz)
        end do
       
        write (message(1),'(a,F26.12)') 'Info: Point volume calculation finished. Total volume :',sum(mesh%vol_pp(1:mesh%np))
        call write_info(1)

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

    call pop_sub('mesh_init.mesh_init_stage_3.mesh_get_vol_pp')

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

      SAFE_ALLOCATE(mesh%per_points(1:mesh%nper))
      SAFE_ALLOCATE(mesh%per_map(1:mesh%nper))

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then
        SAFE_ALLOCATE(recv_points(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(recv_rem_points(1:nper_recv, 1:mesh%vp%npart))
        SAFE_ALLOCATE(mesh%nrecv(1:mesh%vp%npart))
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
            
            if(ip_inner /= 0) then
              if(ip_inner <= mesh%vp%np_local(ipart)) then
                ! count the points to receive from each node
                mesh%nrecv(ipart) = mesh%nrecv(ipart) + 1
                ! and store the number of the point
                recv_points(mesh%nrecv(ipart), ipart) = ip
                ! and where it is in the other partition
                recv_rem_points(mesh%nrecv(ipart), ipart) = ip_inner

                ASSERT(mesh%vp%rank /= ipart - 1) ! if we are here, the point must be in another node
              
                exit
              end if
            end if
            
          end do
#endif
        end if
      end do

#ifdef HAVE_MPI
      if(mesh%parallel_in_domains) then

        ! first we allocate the buffer to be able to use MPI_Bsend
        bsize = mesh%vp%npart - 1 + nper_recv + MPI_BSEND_OVERHEAD*2*(mesh%vp%npart - 1)
        SAFE_ALLOCATE(send_buffer(1:bsize))
        call MPI_Buffer_attach(send_buffer(1), bsize*4, mpi_err)

        ! Now we communicate to each node the points they will have to
        ! send us. Probably this could be done without communication,
        ! but this way it seems simpler to implement.
        
        ! We send the number of points we expect to receive.
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Bsend(mesh%nrecv(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, mpi_err)
        end do

        ! And we receive it
        SAFE_ALLOCATE(mesh%nsend(1:mesh%vp%npart))
        mesh%nsend = 0
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle
          call MPI_Recv(mesh%nsend(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, status, mpi_err)
        end do

        ! Now we send the indexes of the points
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. mesh%nrecv(ipart) == 0) cycle
          call MPI_Bsend(recv_rem_points(1, ipart), mesh%nrecv(ipart), MPI_INTEGER, ipart - 1, 1, mesh%vp%comm, mpi_err)
        end do

        SAFE_ALLOCATE(send_points(1:maxval(mesh%nsend), 1:mesh%vp%npart))

        ! And we receive them
        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno .or. mesh%nsend(ipart) == 0) cycle
          call MPI_Recv(send_points(1, ipart), mesh%nsend(ipart), MPI_INTEGER, &
               ipart - 1, 1, mesh%vp%comm, status, mpi_err)
        end do

        ! we no longer need this
        SAFE_DEALLOCATE_A(recv_rem_points)

        ! Now we have all the indexes required locally, so we can
        ! build the mpi datatypes

        SAFE_ALLOCATE(mesh%dsend_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(mesh%zsend_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(mesh%drecv_type(1:mesh%vp%npart))
        SAFE_ALLOCATE(mesh%zrecv_type(1:mesh%vp%npart))

        maxmax = max(maxval(mesh%nsend), maxval(mesh%nrecv))

        SAFE_ALLOCATE(blocklengths(1:maxmax))
        SAFE_ALLOCATE(offsets(1:maxmax))

        do ipart = 1, mesh%vp%npart
          if(ipart == mesh%vp%partno) cycle

          if(mesh%nsend(ipart) > 0) then

            ASSERT(all(send_points(1:mesh%nsend(ipart), ipart) <= mesh%np))

            ! MPI indexes start from zero
            send_points(1:mesh%nsend(ipart), ipart) = send_points(1:mesh%nsend(ipart), ipart) - 1

            call get_blocks(mesh%nsend(ipart), send_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_FLOAT, mesh%dsend_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_CMPLX, mesh%zsend_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%dsend_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%zsend_type(ipart), mpi_err)

          end if
          
          if(mesh%nrecv(ipart) > 0) then
            ASSERT(all(recv_points(1:mesh%nrecv(ipart), ipart) <= mesh%np_part))
            ASSERT(all(recv_points(1:mesh%nrecv(ipart), ipart) > mesh%np))

            recv_points(1:mesh%nrecv(ipart), ipart) = recv_points(1:mesh%nrecv(ipart), ipart) - 1

            call get_blocks(mesh%nrecv(ipart), recv_points(:, ipart), nblocks, blocklengths, offsets)

            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_FLOAT, mesh%drecv_type(ipart), mpi_err)
            call MPI_Type_indexed(nblocks, blocklengths(1), offsets(1), MPI_CMPLX, mesh%zrecv_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%drecv_type(ipart), mpi_err)
            call MPI_Type_commit(mesh%zrecv_type(ipart), mpi_err)

          end if

        end do

        call MPI_Buffer_detach(send_buffer(1), bsize, mpi_err)

      end if
#endif
      ASSERT(iper == mesh%nper)

    end if

    call pop_sub('mesh_init.mesh_init_stage_3.pbc_init')
  end subroutine mesh_pbc_init

end subroutine mesh_init_stage_3


#ifdef HAVE_MPI
! ---------------------------------------------------------------
!> Converts the mesh given by grid points into a graph. Each
!! point is a vertex in the graph and closest neighbours are
!! connected by an edge (at most 6 in 3D and 4 in 2D, 2 in
!! 1D, fewer at the boundaries).
!! Then calls METIS to get npart partitions.
!! Stored the mapping point no. -> partition no. into part,
!! which has to be allocated beforehand.
!! (mesh_partition_end should be called later.)
! ---------------------------------------------------------------
subroutine mesh_partition(mesh, lapl_stencil, part)
  type(mesh_t),            intent(in)  :: mesh
  type(stencil_t), target, intent(in)  :: lapl_stencil
  integer,                 intent(out) :: part(:)

  integer              :: iv, jp, inb
  integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
  integer              :: ne             !< Number of edges.
  integer              :: nv             !< Number of vertices.
  ! Number of vertices (nv) is equal to number of
  ! points np_global and maximum number of edges (ne) is 2*mesh%sb%dim*np_global
  ! (there are a little fewer because points on the border have fewer
  ! than two neighbours per dimension).
  ! xadj has nv+1 entries because last entry contains the total
  ! number of edges.
  integer              :: npart          !< Number of partitions.
  integer              :: ipart          !< number of the current partition
  integer, allocatable :: xadj(:)        !< Indices of adjacency list in adjncy.
  integer, allocatable :: adjncy(:)      !< Adjacency lists.
  integer              :: iunit          !< For debug output to files.
#ifdef HAVE_METIS
  integer              :: options(5)     !< Options to METIS.
  integer              :: edgecut        !< Number of edges cut by partitioning.
#endif

  type(stencil_t) :: stencil
  integer :: ii, stencil_to_use, ip
  integer :: default_method, method
  integer :: library
  integer, parameter :: METIS = 2, ZOLTAN = 3
  integer, parameter :: STAR = 1, LAPLACIAN = 2
  integer, allocatable :: start(:), final(:), lsize(:)
  FLOAT, allocatable :: xglobal(:, :)

  type(profile_t), save :: prof
  integer :: default
 
  call profiling_in(prof, "MESH_PARTITION")
  call push_sub('mesh_init.mesh_partition')

  if(mesh%np_global == 0) then
    message(1) = 'The mesh is empty and cannot be partitioned.'
    call write_fatal(1)
  end if

  !%Variable MeshPartitionPackage
  !%Type integer
  !%Default metis
  !%Section Execution::Parallelization
  !%Description
  !% Decides which library to use to perform the mesh partition. By
  !% default, METIS is used (if available).
  !%Option metis 2
  !% METIS library.
  !%Option zoltan 3
  !% Zoltan library.
  !%End
  default = ZOLTAN
#ifdef HAVE_METIS
  default = METIS
#endif
  call parse_integer(datasets_check('MeshPartitionPackage'), default, library)

#ifndef HAVE_METIS
  if(library == METIS) then
    message(1) = 'Error: METIS was requested, but Octopus was compiled without it.'
    call write_fatal(1)
  end if
#endif

  !%Variable MeshPartitionStencil
  !%Type integer
  !%Default star
  !%Section Execution::Parallelization
  !%Description
  !% To partition the mesh, it is necessary to calculate the connection
  !% graph connecting the points. This variable selects which stencil
  !% is used to do this. The default is the order-one star stencil.
  !% Alternatively, the stencil used for the Laplacian may be used.
  !%Option stencil_star 1
  !% An order-one star stencil.
  !%Option laplacian 2
  !% The stencil used for the Laplacian is used to calculate the
  !% partition. This in principle should give a better partition, but
  !% it is slower and requires more memory.
  !%End
  call parse_integer(datasets_check('MeshPartitionStencil'), STAR, stencil_to_use)

  if (stencil_to_use == STAR) then
    call stencil_star_get_lapl(stencil, mesh%sb%dim, order = 1)
  else if (stencil_to_use == LAPLACIAN) then
    call stencil_copy(lapl_stencil, stencil)
  else
    call input_error('MeshPartitionStencil')
  end if

  ! Get number of partitions.
  npart = mesh%mpi_grp%size
  ipart = mesh%mpi_grp%rank + 1

  if(npart .lt. 8) then
    default_method = RCB
  else
    default_method = GRAPH
  end if
  ! Documentation is in zoltan.F90
  call parse_integer(datasets_check('MeshPartition'), default_method, method)

  SAFE_ALLOCATE(start(1:npart))
  SAFE_ALLOCATE(final(1:npart))
  SAFE_ALLOCATE(lsize(1:npart))

  select case(library)
  case(METIS)

    start(1:npart) = 1
    final(1:npart) = mesh%np_global
    lsize(1:npart) = mesh%np_global

  case(ZOLTAN)

    ! If we use Zoltan, we divide the space in a basic way, to balance
    ! the memory for the graph. 
    call multicomm_divide_range(mesh%np_global, npart, start, final, lsize)

    do ii = 1, npart
      part(start(ii):final(ii)) = ii
    end do

  end select

  ! Shortcut (number of vertices).
  nv = lsize(ipart)
  SAFE_ALLOCATE(xadj(1:nv + 1))

  if(library == METIS .or. .not. zoltan_method_is_geometric(method)) then !calculate the graphs
    SAFE_ALLOCATE(adjncy(1:stencil%size*nv))

    ! Create graph with each point being
    ! represented by a vertex and edges between
    ! neighbouring points.
    ne = 1
    ! Iterate over number of vertices.
    do iv = 1, nv
      ! Get coordinates of point iv (vertex iv).
      call index_to_coords(mesh%idx, mesh%sb%dim, iv, ix)
      ! Set entry in index table.
      xadj(iv) = ne
      ! Check all possible neighbours.
      do jp = 1, stencil%size 
        ! Store coordinates of possible neighbors, they
        ! are needed several times in the check below.
        jx(1:MAX_DIM) = ix(1:MAX_DIM) + stencil%points(1:MAX_DIM, jp)

        if(all(jx(1:MAX_DIM) >= mesh%idx%nr(1, 1:MAX_DIM)) .and. all(jx(1:MAX_DIM) <= mesh%idx%nr(2, 1:MAX_DIM))) then
          ! Only points inside the mesh or its enlargement
          ! are included in the graph.
          inb = index_from_coords(mesh%idx, mesh%sb%dim, jx)
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
    ! last element unnecessary (this indexing is a
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
        do iv = 1, nv
          write(iunit, *) adjncy(xadj(iv):xadj(iv+1) - 1)
        end do
        call io_close(iunit)
      end if
    end if
    
  else
    SAFE_ALLOCATE(adjncy(1:1))
  end if

  select case(library)
  case(METIS)
#ifdef HAVE_METIS
    options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

    ! Partition graph.
    ! Recursive bisection is better for small number of partitions (<8),
    ! multilevel k-way otherwise (cf. METIS manual).
    ! If the graph contains no vertices, METIS cannot be called. This seems
    ! to happen, e.g., when using minimum BoxShape without any atoms in the
    ! input file.

    select case(method)
    case(RCB)
      message(1) = 'Info: Using METIS multilevel recursive bisection to partition the mesh.'
      call write_info(1)
      call oct_metis_part_graph_recursive(nv, xadj, adjncy, &
           0, 0, 0, 1, npart, options, edgecut, part)
    case(GRAPH)
      message(1) = 'Info: Using METIS multilevel k-way algorithm to partition the mesh.'
      call write_info(1)
      call oct_metis_part_graph_kway(nv, xadj, adjncy, &
           0, 0, 0, 1, npart, options, edgecut, part)
    case default
      message(1) = 'Error: Selected partition method is not available in METIS.'
      call write_fatal(1)
    end select
#endif
  case(ZOLTAN)

    call zoltan_method_info(method)

    SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM))

    do ip = 1, mesh%np_part_global
      xglobal(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip)
    end do

    !assign all points to one node
    call zoltan_partition(method, mesh%sb%dim, mesh%np_global, mesh%np_part_global, &
         xglobal(1, 1),  start(ipart), xadj(1), adjncy(1), ipart, part(1), mesh%mpi_grp%comm)

    SAFE_DEALLOCATE_A(xglobal)

    ! we use xadj as a buffer
    xadj(1:lsize(ipart)) = part(start(ipart):final(ipart))

    ASSERT(all(xadj(1:lsize(ipart)) > 0))

    part(1:mesh%np_global) = 0 ! so we catch non-initialized values

    ! convert start to C notation
    start = start - 1

    ! we collect part from all processors
    call MPI_Allgatherv(xadj(1), lsize(ipart), MPI_INTEGER, part(1), lsize(1), start(1), MPI_INTEGER, mesh%mpi_grp%comm, mpi_err)

  end select

  SAFE_DEALLOCATE_A(start)
  SAFE_DEALLOCATE_A(final)
  SAFE_DEALLOCATE_A(lsize)
  SAFE_DEALLOCATE_A(adjncy)

  ASSERT(all(part(1:mesh%np_global) > 0))
  ASSERT(all(part(1:mesh%np_global) <= npart))

  call stencil_end(stencil)
  call pop_sub('mesh_init.mesh_partition')
  call profiling_out(prof)

end subroutine mesh_partition

subroutine mesh_partition_boundaries(mesh, stencil, part)
  type(mesh_t),    intent(in)    :: mesh
  type(stencil_t), intent(in)    :: stencil
  integer,         intent(inout) :: part(:)

  integer              :: ii, jj         ! Counter.
  integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
  integer              :: npart
  integer              :: iunit          ! For debug output to files.
  character(len=3)     :: filenum
  integer, allocatable :: votes(:), bps(:)
  logical, allocatable :: winner(:)
  integer :: ip, maxvotes

  call push_sub('mesh_init.mesh_partition_boundaries')

  npart = mesh%mpi_grp%size

  SAFE_ALLOCATE(votes(1:npart))
  SAFE_ALLOCATE(bps(1:npart))
  SAFE_ALLOCATE(winner(1:npart))

  !assign boundary points
  bps = 0
  do ii = mesh%np_global + 1, mesh%np_part_global
    !get the coordinates of the point
    call index_to_coords(mesh%idx, mesh%sb%dim, ii, ix)
    votes = 0
    ! check the partition of all the points that are connected through the stencil with this one
    do jj = 1, stencil%size
      jx(1:MAX_DIM) = ix(1:MAX_DIM) + stencil%points(1:MAX_DIM, jj)
      if(any(jx < mesh%idx%nr(1, :)) .or. any(jx > mesh%idx%nr(2, :))) cycle
      ip = index_from_coords(mesh%idx, mesh%sb%dim, jx)
      if(ip > 0 .and. ip <= mesh%np_global) votes(part(ip)) = votes(part(ip)) + 1
    end do
    
    ! now count the votes
    maxvotes = maxval(votes)
    ! from all the ones that have the maximum
    winner = (votes == maxvotes)
    ! select the one that has fewer points currently
    part(ii) = minloc(bps, dim = 1,  mask = winner)
    ! and count it
    bps(part(ii)) = bps(part(ii)) + 1

  end do

  SAFE_DEALLOCATE_A(votes)
  SAFE_DEALLOCATE_A(bps)
  SAFE_DEALLOCATE_A(winner)

  if(in_debug_mode.and.mpi_grp_is_root(mpi_world)) then
    ! Debug output. Write points of each partition in a different file.
    do ii = 1, npart

      write(filenum, '(i3.3)') ii

      iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
           action='write')
      do jj = 1, mesh%np_global
        if(part(jj).eq.ii) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
      end do
      call io_close(iunit)

      iunit = io_open('debug/mesh_partition/mesh_partition_all.'//filenum, &
           action='write')
      do jj = 1, mesh%np_part_global
        if(part(jj).eq.ii) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
      end do
      call io_close(iunit)

    end do
    ! Write points from enlargement to file with number p+1.
    write(filenum, '(i3.3)') npart+1
    iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
         action='write')
    do ii = mesh%np_global+1, mesh%np_part_global
      write(iunit, '(i8,3f18.8)') ii, mesh_x_global(mesh, ii)
    end do
    call io_close(iunit)
  end if

  call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)

  call pop_sub('mesh_init.mesh_partition_boundaries')

end subroutine mesh_partition_boundaries

#endif
end module mesh_init_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
