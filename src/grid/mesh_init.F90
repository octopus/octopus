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

#define ENLARGEMENT_POINT 2
#define INNER_POINT 1
! ---------------------------------------------------------
subroutine mesh_init_stage_1(sb, mesh, geo, cv, enlarge)
  type(simul_box_t), target, intent(in)    :: sb
  type(mesh_t),              intent(inout) :: mesh
  type(geometry_t),          intent(in)    :: geo
  type(curvlinear_t),        intent(in)    :: cv
  integer,                   intent(in)    :: enlarge(MAX_DIM)

  integer :: i, j
  FLOAT   :: x(MAX_DIM), chi(MAX_DIM)
  logical :: out

  call push_sub('mesh_init.mesh_init_stage_1')
  call profiling_in(mesh_init_prof, "MESH_INIT")

  mesh%sb => sb   ! keep an internal pointer
  mesh%h  =  sb%h ! this number can change in the following
  mesh%use_curvlinear = cv%method.ne.CURV_METHOD_UNIFORM
  mesh%enlarge = enlarge

  ! adjust nr
  mesh%nr = 0
  do i = 1, sb%dim
    chi = M_ZERO; j = 0
    out = .false.
    do while(.not.out)
      j      = j + 1
      chi(i) = j*mesh%h(i)
      if ( mesh%use_curvlinear ) then
        call curvlinear_chi2x(sb, geo, cv, chi(1:sb%dim), x(1:sb%dim))
        out = (x(i) > sb%lsize(i)+CNST(1.0e-10))
      else
        out = (chi(i) > sb%lsize(i)+CNST(1.0e-10))
      end if
    end do
    mesh%nr(2, i) = j - 1
  end do

  ! we have a symmetric mesh (for now)
  mesh%nr(1,:) = -mesh%nr(2,:)

  ! we have to adjust a couple of things for the periodic directions
  do i = 1, sb%periodic_dim
    !the spacing has to be a divisor of the box size
    mesh%h(i)     = sb%lsize(i)/real(mesh%nr(2, i))
    !the upper boundary does not have to be included (as it is a copy of the lower boundary)
    mesh%nr(2, i) = mesh%nr(2, i) - 1
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
    mesh%nr(2, TRANS_DIR) = mesh%nr(2, TRANS_DIR) - 1

    call mesh_read_lead()
  else
    nullify(mesh%lead_unit_cell)
  end if

  mesh%l(:) = mesh%nr(2, :) - mesh%nr(1, :) + 1

  call profiling_out(mesh_init_prof)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine mesh_read_lead()
    integer :: il, iunit

    integer               :: alloc_size
    type(mesh_t), pointer :: m

    call push_sub('mesh_init.mesh_read_lead')

    ALLOCATE(mesh%lead_unit_cell(NLEADS), NLEADS)

    do il = 1, NLEADS
      m => mesh%lead_unit_cell(il)
      iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/mesh', action='read', is_tmp=.true.)
      call mesh_init_from_file(m, iunit)
      call io_close(iunit)

      ! Read the lxyz maps.
      ALLOCATE(m%lxyz(m%np_part, 3), m%np_global*3)
      alloc_size = (m%nr(2, 1)-m%nr(1, 1)+1) * (m%nr(2, 2)-m%nr(1, 2)+1) * (m%nr(2, 3)-m%nr(1, 3)+1)
      ALLOCATE(m%lxyz_inv(m%nr(1, 1):m%nr(2, 1), m%nr(1, 2):m%nr(2, 2), m%nr(1, 3):m%nr(2, 3)), alloc_size)
      iunit = io_open(trim(sb%lead_restart_dir(il))//'/gs/Lxyz', action='read', is_tmp=.true.)
      call mesh_lxyz_init_from_file(m, iunit)
      call io_close(iunit)
    end do

    call pop_sub()
  end subroutine mesh_read_lead
end subroutine mesh_init_stage_1


! ---------------------------------------------------------
subroutine mesh_init_stage_2(sb, mesh, geo, cv)
  type(simul_box_t),  intent(in)    :: sb
  type(mesh_t),       intent(inout) :: mesh
  type(geometry_t),   intent(in)    :: geo
  type(curvlinear_t), intent(in)    :: cv

  integer :: i, j, k, il, ik, ix, iy, iz
  FLOAT   :: chi(MAX_DIM)

  call push_sub('mesh_init.mesh_init_stage_2')
  call profiling_in(mesh_init_prof)

  ! enlarge mesh for boundary points
  mesh%nr(1,:) = mesh%nr(1,:) - mesh%enlarge(:)
  mesh%nr(2,:) = mesh%nr(2,:) + mesh%enlarge(:)

  ! allocate the xyz arrays
  i = (mesh%nr(2,1)-mesh%nr(1,1)+1) * (mesh%nr(2,2)-mesh%nr(1,2)+1) * (mesh%nr(2,3)-mesh%nr(1,3)+1)
  ALLOCATE(mesh%Lxyz_inv(mesh%nr(1,1):mesh%nr(2,1), mesh%nr(1,2):mesh%nr(2,2), mesh%nr(1,3):mesh%nr(2,3)),   i)
  ALLOCATE(mesh%Lxyz_tmp(mesh%nr(1,1):mesh%nr(2,1), mesh%nr(1,2):mesh%nr(2,2), mesh%nr(1,3):mesh%nr(2,3)),   i)
  ALLOCATE(mesh%x_tmp(3, mesh%nr(1,1):mesh%nr(2,1), mesh%nr(1,2):mesh%nr(2,2), mesh%nr(1,3):mesh%nr(2,3)), 3*i)

  !$omp parallel workshare
  mesh%Lxyz_inv(:,:,:) = 0
  mesh%Lxyz_tmp(:,:,:) = 0
  mesh%x_tmp(:,:,:,:)  = M_ZERO
  !$omp end parallel workshare

  ! We label 2 the points inside the mesh + enlargement
  !$omp parallel do private(iy, iz, chi, i, j, k) if(.not. mesh%use_curvlinear)
  do ix = mesh%nr(1,1), mesh%nr(2,1)
    chi(1) = real(ix, REAL_PRECISION) * mesh%h(1) + sb%box_offset(1)

    do iy = mesh%nr(1,2), mesh%nr(2,2)
      chi(2) = real(iy, REAL_PRECISION) * mesh%h(2) + sb%box_offset(2)

      do iz = mesh%nr(1,3), mesh%nr(2,3)
        chi(3) = real(iz, REAL_PRECISION) * mesh%h(3) + sb%box_offset(3)

        call curvlinear_chi2x(sb, geo, cv, chi(:), mesh%x_tmp(:, ix, iy, iz))

        if(simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz))) then
          do i = -mesh%enlarge(1), mesh%enlarge(1)
            do j = -mesh%enlarge(2), mesh%enlarge(2)
              do k = -mesh%enlarge(3), mesh%enlarge(3)
                if(  &
                  ix+i>=mesh%nr(1,1).and.ix+i<=mesh%nr(2,1).and. &
                  iy+j>=mesh%nr(1,2).and.iy+j<=mesh%nr(2,2).and. &
                  iz+k>=mesh%nr(1,3).and.iz+k<=mesh%nr(2,3)) mesh%Lxyz_tmp(ix+i, iy+j, iz+k) = ENLARGEMENT_POINT
              end do
            end do
          end do
        end if

      end do
    end do
  end do
  !$omp end parallel do

  ! we label 1 the points inside the mesh, and we count the points
  il = 0
  ik = 0
  do ix = mesh%nr(1,1), mesh%nr(2,1)
    do iy = mesh%nr(1,2), mesh%nr(2,2)
      do iz = mesh%nr(1,3), mesh%nr(2,3)
        if(simul_box_in_box(sb, geo, mesh%x_tmp(:, ix, iy, iz))) then
          mesh%Lxyz_tmp(ix, iy, iz) = INNER_POINT
          ik = ik + 1
        end if

        if(mesh%Lxyz_tmp(ix, iy, iz) > 0) il = il + 1
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
subroutine mesh_init_stage_3(mesh, geo, cv, stencil, np_stencil, mpi_grp)
  type(mesh_t),       intent(inout) :: mesh
  type(geometry_t),   intent(in)    :: geo
  type(curvlinear_t), intent(in)    :: cv
  integer, optional,  intent(in)    :: stencil(:, :)
  integer, optional,  intent(in)    :: np_stencil
  type(mpi_grp_t), optional,  intent(in) :: mpi_grp

  call push_sub('mesh_init.mesh_init_stage_3')
  call profiling_in(mesh_init_prof)

  ! check if we are running in parallel in domains
  mesh%parallel_in_domains = .false.
  if(present(mpi_grp)) mesh%parallel_in_domains = .true.

  call create_x_Lxyz()

  if(mesh%parallel_in_domains) then
    ASSERT(present(stencil).and.present(np_stencil))

    call do_partition()
  else
    call mpi_grp_init(mesh%mpi_grp, -1)

    ! When running serially those two are the same.
    mesh%np      = mesh%np_global
    mesh%np_part = mesh%np_part_global

    ! These must be initialized for vec_gather, vec_scatter to work
    ! as copy operations when running without domain parallelization.
    mesh%vp%np = mesh%np_global
    mesh%vp%p  = 1
  end if

  call mesh_get_vol_pp(mesh%sb)

  ! these large arrays were allocated in mesh_init_1, and are no longer needed
  deallocate(mesh%Lxyz_tmp, mesh%x_tmp)

  call mesh_pbc_init()

  call profiling_out(mesh_init_prof)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine create_x_Lxyz()
    integer :: il, ix, iy, iz

#ifdef USE_OMP
    integer :: ip
#endif

    ALLOCATE(mesh%Lxyz(mesh%np_part_global, 3), mesh%np_part_global*3)
    if(mesh%parallel_in_domains) then
      ! Node 0 has to store all entries from x (in x_global)
      ! as well as the local set in x (see below).
      ALLOCATE(mesh%x_global(mesh%np_part_global, 3), mesh%np_part_global*3)
    else
      ! When running parallel, x is computed later.
      ALLOCATE(mesh%x(mesh%np_part_global, 3), mesh%np_part_global*3)

#ifdef USE_OMP
      !$omp parallel 
      !$omp do
      do ip = 1, mesh%np_global
        mesh%x(ip, 1:3) = M_ZERO
      end do
      !$omp end do nowait
      !$omp do
      do ip = mesh%np_global+1, mesh%np_part_global
        mesh%x(ip, 1:3) = M_ZERO
      end do
      !$omp end do
      !$omp end parallel
#endif

      ! This is a bit ugly: x_global is needed in out_in
      ! but in the serial case it is the same as x
      mesh%x_global => mesh%x
    end if


    ! first we fill the points in the inner mesh
    il = 0
    do ix = mesh%nr(1,1), mesh%nr(2,1)
      do iy = mesh%nr(1,2), mesh%nr(2,2)
        do iz = mesh%nr(1,3), mesh%nr(2,3)
          if(mesh%Lxyz_tmp(ix, iy, iz) == INNER_POINT) then
            il = il + 1
            mesh%Lxyz(il, 1) = ix
            mesh%Lxyz(il, 2) = iy
            mesh%Lxyz(il, 3) = iz
            mesh%Lxyz_inv(ix,iy,iz) = il
            if(mesh%parallel_in_domains) then
              mesh%x_global(il, :) = mesh%x_tmp(:, ix, iy, iz)
            else
              mesh%x(il,:) = mesh%x_tmp(:,ix,iy,iz)
            end if
          end if
        end do
      end do
    end do

    ! and now the points from the enlargement
    do ix = mesh%nr(1,1), mesh%nr(2,1)
      do iy = mesh%nr(1,2), mesh%nr(2,2)
        do iz = mesh%nr(1,3), mesh%nr(2,3)

          if(mesh%Lxyz_tmp(ix, iy, iz) == ENLARGEMENT_POINT) then
            il = il + 1
            mesh%Lxyz(il, 1) = ix
            mesh%Lxyz(il, 2) = iy
            mesh%Lxyz(il, 3) = iz
            mesh%Lxyz_inv(ix,iy,iz) = il
            if(mesh%parallel_in_domains) then
              mesh%x_global(il, :) = mesh%x_tmp(:, ix, iy, iz)
            else
              mesh%x(il,:) = mesh%x_tmp(:,ix,iy,iz)
            end if
          end if

        end do
      end do
    end do
  end subroutine create_x_Lxyz


  ! ---------------------------------------------------------
  subroutine do_partition()
#if defined(HAVE_METIS) && defined(HAVE_MPI)
    integer :: i, j, ipart, jpart
    integer, allocatable :: part(:), nnb(:)

    mesh%mpi_grp = mpi_grp

    ALLOCATE(part(mesh%np_part_global), mesh%np_part_global)
    call mesh_partition(mesh, stencil, np_stencil, part)
    call vec_init(mesh%mpi_grp%comm, 0, part, mesh%np_global, mesh%np_part_global,  &
      mesh%nr, mesh%Lxyz_inv, mesh%Lxyz, stencil, np_stencil, mesh%sb%dim, mesh%vp)
    deallocate(part)

    ALLOCATE(nnb(1:mesh%vp%p), mesh%vp%p)
    nnb = 0
    do jpart = 1, mesh%vp%p
      do ipart = 1, mesh%vp%p
        if (ipart == jpart) cycle
        if (mesh%vp%np_ghost_neigh(jpart, ipart) /= 0) nnb(jpart) = nnb(jpart) + 1
      end do
      ASSERT(nnb(jpart) >= 0 .and. nnb(jpart) < mesh%vp%p)
    end do

    ! Write information about partitions.
    message(1) = 'Info: Mesh partition:'
    call write_info(1)
    do ipart = 1, mesh%vp%p
      write(message(1),'(a,i5)')  &
           'Info: Nodes in domain-group ', ipart
      write(message(2),'(a,i10,a,i10)') &
           '      Neighbours   :', nnb(ipart), &
           '      Local points    :', mesh%vp%np_local(ipart)
      write(message(3),'(a,i10,a,i10)') &
           '      Ghost points :', mesh%vp%np_ghost(ipart), &
           '      Boundary points :', mesh%vp%np_bndry(ipart)
      call write_info(3)
    end do
    deallocate(nnb)

    ! Set local point numbers.
    mesh%np      = mesh%vp%np_local(mesh%vp%partno)
    mesh%np_part = mesh%np + mesh%vp%np_ghost(mesh%vp%partno) + mesh%vp%np_bndry(mesh%vp%partno)

    ! Compute mesh%x as it is done in the serial case but only for local points.
    ! x consists of three parts: the local points, the
    ! ghost points, and the boundary points; in this order
    ! (just as for any other vector, which is distributed).
    ALLOCATE(mesh%x(mesh%np_part, 3), mesh%np_part*3)
    ! Do the inner points
    do i = 1, mesh%np
      j  = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+i-1)
      mesh%x(i, :) = mesh%x_tmp(:, mesh%Lxyz(j, 1), mesh%Lxyz(j, 2), mesh%Lxyz(j, 3))
    end do
    ! Do the ghost points
    do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
      j = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno)+i-1)
      mesh%x(i+mesh%np, :) = mesh%x_tmp(:, mesh%Lxyz(j, 1), mesh%Lxyz(j, 2), mesh%Lxyz(j, 3))
    end do
    ! Do the boundary points
    do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
      j = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno)+i-1)
      mesh%x(i + mesh%np + mesh%vp%np_ghost(mesh%vp%partno), :) = &
        mesh%x_tmp(:, mesh%Lxyz(j, 1), mesh%Lxyz(j, 2), mesh%Lxyz(j, 3))
    end do
#endif
  end subroutine do_partition


  ! ---------------------------------------------------------
  ! calculate the volume of integration
  subroutine mesh_get_vol_pp(sb)
    type(simul_box_t), intent(in) :: sb

    integer :: i
    FLOAT   :: chi(MAX_DIM)
#if defined(HAVE_MPI)
    integer :: k
#endif

    call push_sub('mesh_init.mesh_get_vol_pp')

    ALLOCATE(mesh%vol_pp(mesh%np_part), mesh%np_part)

    !$omp parallel workshare
    mesh%vol_pp(:) = product(mesh%h(1:sb%dim))
    !$omp end parallel workshare

    if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
      ! Do the inner points.
      do i = 1, mesh%np
        k = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+i-1)
        chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
        mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi(1:sb%dim))
      end do
      ! Do the ghost points.
      do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
        k = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno)+i-1)
        chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
        mesh%vol_pp(i+mesh%np) = mesh%vol_pp(i+mesh%np)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np, :), chi(1:sb%dim))
      end do
      ! Do the boundary points.
      do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
        k = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno)+i-1)
        chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
        mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) = &
          mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) &
          *curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno), :), chi(1:sb%dim))
      end do
#endif
    else ! serial mode
      do i = 1, mesh%np_part
        chi(1:sb%dim) = mesh%Lxyz(i, 1:sb%dim) * mesh%h(1:sb%dim)
        mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi(1:sb%dim))
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
        ALLOCATE(recv_points(1:nper_recv, 1:mesh%vp%p), nper_recv*mesh%vp%p)
        ALLOCATE(recv_rem_points(1:nper_recv, 1:mesh%vp%p), nper_recv*mesh%vp%p)
        ALLOCATE(mesh%nrecv(1:mesh%vp%p), mesh%vp%p)
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
          do ipart = 1, mesh%vp%p
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
        bsize = mesh%vp%p - 1 + nper_recv + MPI_BSEND_OVERHEAD*2*(mesh%vp%p - 1)
        ALLOCATE(send_buffer(1:bsize), bsize)
        call MPI_Buffer_attach(send_buffer, bsize*4, mpi_err)

        ! Now we communicate to each node the points they will have to
        ! send us. Probably this could be done without communication,
        ! but this way it seems simpler to implement.
        
        ! We send the number of points we expect to receive.
        do ipart = 1, mesh%vp%p
          if(ipart == mesh%vp%partno) cycle
          call MPI_Bsend(mesh%nrecv(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, mpi_err)
        end do

        ! And we receive it
        ALLOCATE(mesh%nsend(1:mesh%vp%p), mesh%vp%p)
        mesh%nsend = 0
        do ipart = 1, mesh%vp%p
          if(ipart == mesh%vp%partno) cycle
          call MPI_Recv(mesh%nsend(ipart), 1, MPI_INTEGER, ipart - 1, 0, mesh%vp%comm, status, mpi_err)
        end do

        ! Now we send the indexes of the points
        do ipart = 1, mesh%vp%p
          if(ipart == mesh%vp%partno .or. mesh%nrecv(ipart) == 0) cycle
          call MPI_Bsend(recv_rem_points(:, ipart), mesh%nrecv(ipart), MPI_INTEGER, ipart - 1, 1, mesh%vp%comm, mpi_err)
        end do

        ALLOCATE(send_points(1:maxval(mesh%nsend), 1:mesh%vp%p), maxval(mesh%nsend)*mesh%vp%p)

        ! And we receive them
        do ipart = 1, mesh%vp%p
          if(ipart == mesh%vp%partno .or. mesh%nsend(ipart) == 0) cycle
          call MPI_Recv(send_points(:, ipart), mesh%nsend(ipart), MPI_INTEGER, &
               ipart - 1, 1, mesh%vp%comm, status, mpi_err)
        end do

        ! we no longer need this
        deallocate(recv_rem_points)

        ! Now we have all the indexes required locally, so we can
        ! build the mpi datatypes

        ALLOCATE(mesh%dsend_type(1:mesh%vp%p), mesh%vp%p)
        ALLOCATE(mesh%zsend_type(1:mesh%vp%p), mesh%vp%p)
        ALLOCATE(mesh%drecv_type(1:mesh%vp%p), mesh%vp%p)
        ALLOCATE(mesh%zrecv_type(1:mesh%vp%p), mesh%vp%p)

        maxmax = max(maxval(mesh%nsend), maxval(mesh%nrecv))

        ALLOCATE(blocklengths(1:maxmax), maxmax)
        ALLOCATE(offsets(1:maxmax), maxmax)

        do ipart = 1, mesh%vp%p
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
subroutine mesh_partition(m, stencil, np_stencil, part)
  type(mesh_t), intent(in)  :: m
  integer,      intent(in)  :: stencil(:, :)
  integer,      intent(in)  :: np_stencil
  integer,      intent(out) :: part(:)

  integer              :: i, j           ! Counter.
  integer              :: ix, iy, iz     ! Counters to iterate over grid.
  integer              :: jx, jy, jz     ! Coordinates of neighbours.
  integer              :: ne             ! Number of edges.
  integer              :: nv             ! Number of vertices.
  integer              :: d(3, 6)        ! Directions of neighbour points.
  integer              :: edgecut        ! Number of edges cut by partitioning.
                                         ! Number of vertices (nv) is equal to number of
                                         ! points np_global and maximum number of edges (ne) is 2*m%sb%dim*np_global
                                         ! (there are a little less because points on the border have less
                                         ! than two neighbours per dimension).
                                         ! xadj has nv+1 entries because last entry contains the total
                                         ! number of edges.
  integer              :: p              ! Number of partitions.
  integer, allocatable :: xadj(:)        ! Indices of adjacency list in adjncy.
  integer, allocatable :: adjncy(:)      ! Adjacency lists.
  integer              :: options(5)     ! Options to METIS.
  integer              :: iunit          ! For debug output to files.
  character(len=3)     :: filenum
  integer, allocatable :: votes(:, :)
  integer :: ip, rr

  call push_sub('mesh_init.mesh_partition')

  options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

  ! Shortcut (number of vertices).
  nv = m%np_global

  ! Get space for partitioning.
  part = 1

  ALLOCATE(xadj(nv + 1), nv + 1)
  ALLOCATE(adjncy(2*m%sb%dim*nv), 2*m%sb%dim*nv)

  ! Get number of partitions.
  call MPI_Comm_Size(m%mpi_grp%comm, p, mpi_err)

  ! Set directions of possible neighbours.
  ! With this ordering of the directions it is possible
  ! to iterate over d(:, i) with i=1, ..., 2*m%sb%dim,
  ! i. e. this works for dim=1, ..., 3 without any special
  ! cases.
  d(:, 1) = (/ 1,  0,  0/)
  d(:, 2) = (/-1,  0,  0/)
  d(:, 3) = (/ 0,  1,  0/)
  d(:, 4) = (/ 0, -1,  0/)
  d(:, 5) = (/ 0,  0,  1/)
  d(:, 6) = (/ 0,  0, -1/)

  ! Create graph with each point being
  ! represenetd by a vertice and edges between
  ! neighboured points.
  ne = 1
  ! Iterate over number of vertices.
  do i = 1, nv
    ! Get coordinates of point i (vertex i).
    ix      = m%Lxyz(i, 1)
    iy      = m%Lxyz(i, 2)
    iz      = m%Lxyz(i, 3)
    ! Set entry in index table.
    xadj(i) = ne
    ! Check all possible neighbours.
    do j = 1, 2*m%sb%dim
      ! Store coordinates of possible neighbors, they
      ! are needed several times in the check below.
      jx = ix + d(1, j)
      jy = iy + d(2, j)
      jz = iz + d(3, j)
      ! Only if the neighbour is in the surrounding box,
      ! Lxyz_tmp has an entry for this point, otherweise
      ! it is out of bounds.
      if(&
           jx >= m%nr(1, 1) .and. jx <= m%nr(2, 1) .and.  &
           jy >= m%nr(1, 2) .and. jy <= m%nr(2, 2) .and.  &
           jz >= m%nr(1, 3) .and. jz <= m%nr(2, 3)        &
        ) then
        ! Only points inside the mesh or its enlargement
        ! are included in the graph.
        if(m%Lxyz_tmp(jx, jy, jz) /= 0 .and. m%Lxyz_inv(jx, jy, jz) <= nv) then
          ! Store a new edge and increment edge counter.
          adjncy(ne) = m%Lxyz_inv(jx, jy, jz)
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

  ! Partition graph.
  ! Recursive bisection is better for small number of partitions (<8),
  ! multilevel k-way otherwise (cf. METIS manual).
  ! If the graph contains no vertices, METIS cannot be called. This seems
  ! to happen, e. g., when using minimum BoxShape without any atoms in the
  ! input file.
  if(nv.eq.0) then
    message(1) = 'The mesh is empty and cannot be partitioned.'
    call write_fatal(1)
  end if
  if(p == 1) then
    part(:) = 1
  else if(p .lt. 8) then
    message(1) = 'Info: Using multilevel recursive bisection to partition mesh.'
    call write_info(1)
    call oct_metis_part_graph_recursive(nv, xadj, adjncy, &
      0, 0, 0, 1, p, options, edgecut, part)
  else
    message(1) = 'Info: Using multilevel k-way algorithm to partition mesh.'
    call write_info(1)
    call oct_metis_part_graph_kway(nv, xadj, adjncy, &
      0, 0, 0, 1, p, options, edgecut, part)
  end if

  if(in_debug_mode.and.mpi_grp_is_root(mpi_world)) then
    ! Debug output. Write points of each partition in a different file.
    do i = 1, p
      write(filenum, '(i3.3)') i
      iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
         action='write')
      do j = 1, m%np_global
        if(part(j).eq.i) write(iunit, '(i8,3f18.8)') j, m%x_global(j, :)
      end do
      call io_close(iunit)
    end do
    ! Write points from enlargement to file with number p+1.
    write(filenum, '(i3.3)') p+1
    iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
       action='write')
    do i = m%np_global+1, m%np_part_global
      write(iunit, '(i8,3f18.8)') i, m%x_global(i, :)
    end do
    call io_close(iunit)
  end if

  ALLOCATE(votes(1:p, m%np_global + 1:m%np_part_global), p*(m%np_part_global - m%np_global))

  !now assign boundary points

  !count the boundary points that each point needs
  votes = 0
  do i = 1, m%np_global
    do j = 1, np_stencil
      jx = m%Lxyz(i, 1) + stencil(1, j)
      jy = m%Lxyz(i, 2) + stencil(2, j)
      jz = m%Lxyz(i, 3) + stencil(3, j)
      ip = m%Lxyz_inv(jx, jy, jz)
      if(ip > m%np_global) votes(part(i), ip) = votes(part(i), ip) + 1
    end do
  end do

  rr = 1
  do i = m%np_global + 1, m%np_part_global
    if(maxval(votes(1:p, i)) /= minval(votes(1:p, i))) then
      ! we have a winner that takes the point
      part(i:i) = maxloc(votes(1:p, i))
    else
      ! points without a winner are assigned in a round robin fashion
      part(i) = rr
      rr = rr + 1
      if(rr > p) rr = 1
    end if
  end do

  call MPI_Barrier(m%mpi_grp%comm, mpi_err)

  deallocate(xadj, adjncy)
  call pop_sub()

end subroutine mesh_partition
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
