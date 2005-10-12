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

#if defined(HAVE_METIS) && defined(HAVE_MPI)
! Converts the mesh given by grid points into a graph. Each point is
! a vertex in the graph and closest neighbours are connected by an
! edge (at max. 6 in 3D and 4 in 2D, 2 in 1D, less at the boundaries).
! Then calls METIS to get p partitions.
! Stored the mapping point no. -> partition no. into
! part, which has to be allocated beforehand.
! (mesh_partition_end should be called later.)
! In Lxyz_tmp has to be stored which points belong to the inner mesh
! and the enlargement. All other entries have to be zero.
! comm is used to get the number of partitions.
subroutine mesh_partition(m, Lxyz_tmp, comm, part)
  type(mesh_type), intent(in)  :: m
  integer,         intent(in)  :: Lxyz_tmp(m%nr(1,1):m%nr(2,1), &
                                          m%nr(1,2):m%nr(2,2), &
                                          m%nr(1,3):m%nr(2,3))
  integer,         intent(in)  :: comm
  integer,         intent(out) :: part(:)

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
  integer              :: p          ! Number of partitions.
  integer              :: xadj(m%np_part_global+1)
                                     ! Indices of adjacency list
                                     ! in adjncy.
  integer              :: adjncy(2*m%sb%dim*m%np_part_global)
                                     ! Adjacency lists.
  integer              :: options(5) ! Options to METIS.
  integer              :: ierr       ! MPI error code.
  integer, allocatable :: ppp(:)     ! Points per partition.
  integer              :: iunit      ! For debug output to files.
  character(len=3)     :: filenum

  call push_sub('mesh_create.mesh_partition')
  
  options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

  ! Shortcut (number of vertices).
  nv = m%np_part_global

  ! Get space for partitioning.
  part = 1

  ! Get number of partitions.
  call MPI_Comm_size(comm, p, ierr)

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
      jx = ix+d(1, j)
      jy = iy+d(2, j)
      jz = iz+d(3, j)
      ! Only if the neighbour is in the surrounding box,
      ! Lxyz_tmp has an entry for this point, otherweise
      ! it is out of bounds.
      if(jx.ge.m%nr(1, 1).and.jx.le.m%nr(2, 1).and.  &
         jy.ge.m%nr(1, 2).and.jy.le.m%nr(2, 2).and.  &
         jz.ge.m%nr(1, 3).and.jz.le.m%nr(2, 3)) then
        ! Only points inside the mesh or its enlargement
        ! are included in the graph.
        if(Lxyz_tmp(jx, jy, jz).ne.0) then
          ! Store a new edge and increment edge counter.
          adjncy(ne) = m%Lxyz_inv(jx, jy, jz)
          ne         = ne+1
        end if
      end if
    end do
  end do
  ne         = ne-1 ! We start with ne=1 for simplicity. This is off by one
                    ! in the end --> -1.
  xadj(nv+1) = ne+1 ! Set number of edges plus 1 as last index.
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
     if(mpiv%node.eq.0) then
        call io_mkdir('debug/mesh_partition')
        iunit = io_open('debug/mesh_partition/mesh_graph.txt', action='write')
        write(iunit, *) nv, ne/2
        do i = 1, nv 
           write(iunit, *) adjncy(xadj(i):xadj(i+1)-1)
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
  else if(p.lt.8) then
    message(1) = 'Info: Using multilevel recursive bisection to partition mesh.'
    call write_info(1)
    call oct_metis_part_graph_recursive(nv, xadj, adjncy, &
       0, 0, 0, 1, p, options, edgecut, part)
  else
    message(1) = 'Info: Using multilevel k-way algorithm to partition mesh.'
    call write_info(1)
    call oct_metis_part_graph_kway(nv, xadj, adjncy, &
      0, 0, 0, 1, p, options, edgecut, part)
  endif

  ! Write information about partitions.
  ! Count points in each partition.
  allocate(ppp(p))
  ppp = 0
  do i = 1, nv
    ppp(part(i)) = ppp(part(i))+1
  end do
  message(1) = 'Info: Number of points in partitions:'
  write(message(2), '(a,100i7)') 'Info: ', ppp
  call write_info(2)
  deallocate(ppp)

  if(in_debug_mode) then
    ! Debug output. Write points of each partition in a different file.
    if(mpiv%node.eq.0) then
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
  end if

  call pop_sub()

end subroutine mesh_partition
#endif


subroutine mesh_create_xyz(sb, m, cv, geo, stencil, np_stencil, comm)
  type(simul_box_type),  intent(in)    :: sb
  type(mesh_type),       intent(inout) :: m
  type(curvlinear_type), intent(in)    :: cv
  type(geometry_type),   intent(in)    :: geo
  ! When running parallel in domains, stencil and np_stencil
  ! are needed to compute the ghost points.
  ! comm the communicator the parallelization takes places in
  ! and is used to determine the number of partitions.
  integer, optional,     intent(in)    :: stencil(:, :)
  integer, optional,     intent(in)    :: np_stencil
  integer, optional,     intent(in)    :: comm

  integer              :: i, j, k, il, ix, iy, iz
  integer, allocatable :: part(:)
  integer, allocatable :: Lxyz_tmp(:,:,:)
  FLOAT,   allocatable :: x(:,:,:,:)
  FLOAT                :: chi(3)

  call push_sub('mesh_create.mesh_create_xyz')

  ! enlarge mesh in the non-periodic dimensions
    m%nr(1,:) = m%nr(1,:) - m%enlarge(:)
    m%nr(2,:) = m%nr(2,:) + m%enlarge(:)

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(  Lxyz_tmp(m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))
  allocate(       x(3,m%nr(1,1):m%nr(2,1),m%nr(1,2):m%nr(2,2),m%nr(1,3):m%nr(2,3)))

  m%Lxyz_inv(:,:,:) = 0
  Lxyz_tmp(:,:,:)   = 0
  x(:,:,:,:)        = M_ZERO

  ! We label 2 the points inside the mesh + enlargement

  do ix = m%nr(1,1), m%nr(2,1)
    chi(1) = real(ix, PRECISION) * m%h(1) + sb%box_offset(1)

    do iy = m%nr(1,2), m%nr(2,2)
      chi(2) = real(iy, PRECISION) * m%h(2) + sb%box_offset(2)

      do iz = m%nr(1,3), m%nr(2,3)
        chi(3) = real(iz, PRECISION) * m%h(3) + sb%box_offset(3)

        call curvlinear_chi2x(sb, geo, cv, chi(:), x(:, ix, iy, iz))

        if(simul_box_in_box(sb, geo, x(:, ix, iy, iz))) then
          do i = -m%enlarge(1), m%enlarge(1)
            do j = -m%enlarge(2), m%enlarge(2)
              do k = -m%enlarge(3), m%enlarge(3)
                if(  &
                   ix+i>=m%nr(1,1).and.ix+i<=m%nr(2,1).and. &
                   iy+j>=m%nr(1,2).and.iy+j<=m%nr(2,2).and. &
                   iz+k>=m%nr(1,3).and.iz+k<=m%nr(2,3)) Lxyz_tmp(ix+i, iy+j, iz+k) = 2
              end do
            end do
          end do
        end if

      end do
    end do
  end do

  ! we label 1 the points inside the mesh, and we count the points
  il = 0
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(simul_box_in_box(sb, geo, x(:, ix, iy, iz))) Lxyz_tmp(ix, iy, iz) = 1
        if(Lxyz_tmp(ix, iy, iz) > 0) il = il + 1
      end do
    end do
  end do
  m%np_part_global = il
  allocate(m%Lxyz(m%np_part_global, 3))
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  ! Node 0 has to store all entries from x (in x_global)
  ! as well as the local set in x (see below).
    allocate(m%x_global(m%np_part_global, 3))
#else
  ! When running parallel, x is computed later.
  allocate(m%x(m%np_part_global, 3))
  ! This is a bit ugly: x_global is needed in out_in
  ! but in the serial case it is the same as x.
  m%x_global => m%x
#endif


  ! first we fill the points in the inner mesh
  il = 0
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 1) then
          il = il + 1
          m%Lxyz(il, 1) = ix
          m%Lxyz(il, 2) = iy
          m%Lxyz(il, 3) = iz
          m%Lxyz_inv(ix,iy,iz) = il
#if defined(HAVE_MPI) && defined(HAVE_METIS)
          m%x_global(il, :) = x(:, ix, iy, iz)
#else
          m%x(il,:) = x(:,ix,iy,iz)
#endif
        end if
      enddo
    enddo
  enddo
  m%np_global = il

  ! and now the points from the enlargement
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), m%nr(2,3)
        if(Lxyz_tmp(ix, iy, iz) == 2) then
          il = il + 1
          m%Lxyz(il, 1) = ix
          m%Lxyz(il, 2) = iy
          m%Lxyz(il, 3) = iz
          m%Lxyz_inv(ix,iy,iz) = il
#if defined(HAVE_MPI) && defined(HAVE_METIS)
          m%x_global(il, :) = x(:, ix, iy, iz)
#else
          m%x(il,:) = x(:,ix,iy,iz)
#endif
        end if
      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(HAVE_METIS)
    if(present(stencil).and.present(np_stencil).and.present(comm)) then
      allocate(part(m%np_part_global))
      call mesh_partition(m, Lxyz_tmp, comm, part)
      call vec_init(comm, 0, part, m%np_global, m%np_part_global,  &
                    m%nr, m%Lxyz_inv, m%Lxyz, stencil, np_stencil, &
                    m%sb%dim, m%sb%periodic_dim, m%vp)
      deallocate(part)

      ! Set local point numbers.
      m%np      = m%vp%np_local(m%vp%partno)
      m%np_part = m%np + m%vp%np_ghost(m%vp%partno) + m%vp%np_bndry(m%vp%partno)
      ! Compute m%x as it is done in the serial case but
      ! only for local points.
      ! x consists of three parts: the local points, the
      ! ghost points, the boundary points; in this order
      ! (just as for any other vector, which is distributed).
      allocate(m%x(m%np_part, 3))
      ! Do the inner points.
      do i = 1, m%np
        ix = m%Lxyz(m%vp%local(m%vp%xlocal(m%vp%partno)+i-1), 1)
        iy = m%Lxyz(m%vp%local(m%vp%xlocal(m%vp%partno)+i-1), 2)
        iz = m%Lxyz(m%vp%local(m%vp%xlocal(m%vp%partno)+i-1), 3)
        m%x(i, :) = x(:, ix, iy, iz)
      end do
      ! Do the ghost points.
      do i = 1, m%vp%np_ghost(m%vp%partno)
        ix = m%Lxyz(m%vp%ghost(m%vp%xghost(m%vp%partno)+i-1), 1)
        iy = m%Lxyz(m%vp%ghost(m%vp%xghost(m%vp%partno)+i-1), 2)
        iz = m%Lxyz(m%vp%ghost(m%vp%xghost(m%vp%partno)+i-1), 3)
        m%x(i+m%np, :) = x(:, ix, iy, iz)
      end do
      ! Do the boundary points.
      do i = 1, m%vp%np_bndry(m%vp%partno)
        ix = m%Lxyz(m%vp%bndry(m%vp%xbndry(m%vp%partno)+i-1), 1)
        iy = m%Lxyz(m%vp%bndry(m%vp%xbndry(m%vp%partno)+i-1), 2)
        iz = m%Lxyz(m%vp%bndry(m%vp%xbndry(m%vp%partno)+i-1), 3)
        m%x(i+m%np+m%vp%np_ghost(m%vp%partno), :) = x(:, ix, iy, iz)
      end do
    else
      message(1) = 'mesh_create.mesh_create_xyz called without stencil'
      message(2) = 'and np_stencil when running in MPI mode.'
      call write_fatal(2)
    end if
#else
  ! When running serially those two are the same.
  m%np      = m%np_global
  m%np_part = m%np_part_global
#endif

  call mesh_get_vol_pp(sb, geo, cv, m)

  deallocate(Lxyz_tmp)
  deallocate(x)

  call pop_sub()

end subroutine mesh_create_xyz


! calculate the volume of integration
subroutine mesh_get_vol_pp(sb, geo, cv, mesh)
  type(simul_box_type),  intent(in)    :: sb
  type(geometry_type),   intent(in)    :: geo
  type(curvlinear_type), intent(in)    :: cv
  type(mesh_type),       intent(inout) :: mesh

  integer :: i
  FLOAT   :: f, chi(sb%dim)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  integer :: k
#endif

  call push_sub('mesh_create.mesh_get_vol_pp')

  f = M_ONE
  do i = 1, sb%dim
    f = f*mesh%h(i)
  end do

  allocate(mesh%vol_pp(mesh%np_part))
  mesh%vol_pp(:) = f

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  ! Do the inner points.
  do i = 1, mesh%np
    k = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi)
  end do
  ! Do the ghost points.
  do i = 1, mesh%vp%np_ghost(mesh%vp%partno)
    k = mesh%vp%ghost(mesh%vp%xghost(mesh%vp%partno)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i+mesh%np) = mesh%vol_pp(i+mesh%np)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np, :), chi)
  end do
  ! Do the boundary points.
  do i = 1, mesh%vp%np_bndry(mesh%vp%partno)
    k = mesh%vp%bndry(mesh%vp%xbndry(mesh%vp%partno)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) = &
      mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno)) &
      *curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np+mesh%vp%np_ghost(mesh%vp%partno), :), chi)
  end do
#else
  do i = 1, mesh%np_part
    chi(1:sb%dim) = mesh%Lxyz(i, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi)
  end do
#endif

  call pop_sub()

end subroutine mesh_get_vol_pp
