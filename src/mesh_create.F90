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
! Calls mesh_partition_init with the number of partitions
! equal to the number of processors stored in the global mpiv.
! This is what is usually needed.
subroutine mesh_partition_init_default(m, Lxyz_tmp)
  type(mesh_type), intent(inout) :: m
  ! Indices have to be specified, because they do not start with 1.
  integer,         intent(in)    :: Lxyz_tmp(m%nr(1,1):m%nr(2,1), &
                                             m%nr(1,2):m%nr(2,2), &
                                             m%nr(1,3):m%nr(2,3))
  
  call push_sub('mesh_create.mesh_partition_init_default')

  call mesh_partition_init(m, Lxyz_tmp, mpiv%numprocs)

  call pop_sub()
  
end subroutine mesh_partition_init_default


! Converts the mesh given by grid points into a graph. Each point is
! a vertex in the graph and closest neighbours are connected by an
! edge (at max. 6 in 3D and 4 in 2D, 2 in 1D, less at the boundaries).
! Then calls METIS to get p partitions.
! Stored the mapping point no. -> partition no. into
! m%part. m%part is allocated along the way.
! (mesh_partition_end should be called later.)
! In Lxyz_tmp has to be stored which points belong to the inner mesh
! and the enlargement. All other entries have to be zero.
subroutine mesh_partition_init(m, Lxyz_tmp, p)
  type(mesh_type), intent(inout) :: m
  integer,         intent(in)    :: Lxyz_tmp(m%nr(1,1):m%nr(2,1), &
                                             m%nr(1,2):m%nr(2,2), &
                                             m%nr(1,3):m%nr(2,3))
  integer,         intent(in)    :: p 

  integer              :: i, j           ! Counter.
  integer              :: ix, iy, iz     ! Counters to iterate over grid.
  integer              :: jx, jy, jz     ! Coordinates of neighbours.
  integer              :: ne             ! Number of edges.
  integer              :: nv             ! Number of vertices.
  integer              :: d(3, 6)        ! Directions of neighbour points.
  integer              :: edgecut        ! Number of edges cut by partitioning.
  ! Number of vertices (nv) is equal to number of
  ! points np_glob and maximum number of edges (ne) is 2*m%sb%dim*np_glob
  ! (there are a little less because points on the border have less
  ! than two neighbours per dimension).
  ! xadj has nv+1 entries because last entry contains the total
  ! number of edges.
  integer              :: xadj(m%np_tot_glob+1) ! Indices of adjacency list
                                                ! in adjncy.
  integer              :: adjncy(2*m%sb%dim*m%np_tot_glob)
                                                ! Adjacency lists.
  integer              :: options(5)         ! Options to METIS.
  integer, allocatable :: ppp(:)             ! Points per partition.
  integer, pointer     :: part(:)            ! Mapping of nodes to partitions.
#ifdef DEBUG
  integer              :: iunit              ! For debug output to files.
  character(len=3)     :: filenum
#endif

  call push_sub('mesh_create.metis_partition_init')
  
  options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

  ! Shortcut (number of vertices).
  nv = m%np_tot_glob

  ! Get space for partitioning.
  allocate(part(nv))
  part = 1

  ! Store partitioning in mesh.
  m%npart =  p
  m%part  => part

  ! If p=1, we can exit. All points are mapped to
  ! partition 1 (s. a.).
  if(p.eq.1) return

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
      ! Store coordinates of possible neighbour, they
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

#ifdef DEBUG
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
#endif

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
  if(p.lt.8) then
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

#ifdef DEBUG
  ! DEBUG output. Write points of each partition in a different file.
  if(mpiv%node == 0) then
    do i = 1, p 
      write(filenum, '(i3.3)') i
      iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
                      action='write')
      do j = 1, m%np_glob
        ! Somehow all points in the enlargement are mapped
        ! to one point in m%x. So, if the enlargement shall be
        ! written too (s. b.), activate the following line
        ! and take out the usual write and comment out the
        ! block below (this works only for equi-distant grids.
        !!if(part(j).eq.i) write(iunit, '(i8,3i8)') j, m%Lxyz(j, :)
        if(part(j).eq.i) write(iunit, '(4i8)') j, m%Lxyz(j, :)
      end do
      call io_close(iunit)
    end do
    ! Write points from enlargement to file with number p+1.
    !!write(filenum, '(i3.3)') p+1
    !!iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
    !!                action='write')
    !!do i = m%np_glob+1, m%np_tot_glob
    !!  write(iunit, '(i8,3i8)') i, m%Lxyz(i,:)
    !!end do
    !!call io_close(iunit)
  end if
#endif

  call pop_sub()

end subroutine mesh_partition_init


subroutine mesh_partition_end(m)
  type(mesh_type), intent(inout) :: m

  call push_sub('mesh_create.mesh_partition_end')

  ! Free memory used for point-to-partition mapping.
  if(associated(m%part)) then
    deallocate(m%part)
    nullify(m%part)
  end if

  call pop_sub()

end subroutine mesh_partition_end
#endif


subroutine mesh_create_xyz(sb, m, cv, geo, stencil, np_stencil)
  type(simul_box_type),  intent(in)    :: sb
  type(mesh_type),       intent(inout) :: m
  type(curvlinear_type), intent(in)    :: cv
  type(geometry_type),   intent(in)    :: geo
  ! When running parallel in domains, stencil and np_stencil
  ! are needed to compute the ghost points.
  integer, optional,     intent(in)    :: stencil(:, :)
  integer, optional,     intent(in)    :: np_stencil

  integer :: i, j, k, il, ix, iy, iz
  integer, allocatable :: Lxyz_tmp(:,:,:)
  FLOAT,   allocatable :: x(:,:,:,:)
  FLOAT :: chi(3)
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  integer :: ierr
  integer :: rank
#endif

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
  m%np_tot_glob = il
  allocate(m%Lxyz(m%np_tot_glob, 3))
#if defined(HAVE_MPI) && defined(HAVE_METIS)
  ! Node 0 has to store all entries from x (in x_global)
  ! as well as the local set in x (see below).
! BUGFIX
!  if(mpiv%node.eq.0) then
    allocate(m%x_global(m%np_tot_glob, 3))
!  end if
#else
  ! When running parallel, x is computed later.
  allocate(m%x(m%np_tot_glob, 3))
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
  m%np_glob = il

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
  call mesh_partition_init_default(m, Lxyz_tmp)
  if(present(stencil).and.present(np_stencil)) then
    call vec_init_default(m%npart, m%part, m%np_glob, m%np_tot_glob, m%nr, &
                          m%Lxyz_inv, m%Lxyz, stencil, np_stencil, m%vp)

    ! Set local point numbers.
    call MPI_Comm_rank(m%vp%comm, rank, ierr)
    m%np     = m%vp%np_local(rank+1)
    m%np_tot = m%np+m%vp%np_ghost(rank+1)+m%vp%np_bndry(rank+1)
    ! Compute m%x as it is done in the serial case but
    ! only for local points.
    ! x consists of three parts: the local points, the
    ! ghost points, the boundary points; in this order
    ! (just as for any other vector, which is distributed).
    allocate(m%x(m%np_tot, 3))
    ! Do the inner points.
    do i = 1, m%np
      ix = m%Lxyz(m%vp%local(m%vp%xlocal(rank+1)+i-1), 1)
      iy = m%Lxyz(m%vp%local(m%vp%xlocal(rank+1)+i-1), 2)
      iz = m%Lxyz(m%vp%local(m%vp%xlocal(rank+1)+i-1), 3)
      m%x(i, :) = x(:, ix, iy, iz)
    end do
    ! Do the ghost points.
    do i = 1, m%vp%np_ghost(rank+1)
      ix = m%Lxyz(m%vp%ghost(m%vp%xghost(rank+1)+i-1), 1)
      iy = m%Lxyz(m%vp%ghost(m%vp%xghost(rank+1)+i-1), 2)
      iz = m%Lxyz(m%vp%ghost(m%vp%xghost(rank+1)+i-1), 3)
      m%x(i+m%np, :) = x(:, ix, iy, iz)
    end do
    ! Do the boundary points.
    do i = 1, m%vp%np_bndry(rank+1)
      ix = m%Lxyz(m%vp%bndry(m%vp%xbndry(rank+1)+i-1), 1)
      iy = m%Lxyz(m%vp%bndry(m%vp%xbndry(rank+1)+i-1), 2)
      iz = m%Lxyz(m%vp%bndry(m%vp%xbndry(rank+1)+i-1), 3)
      m%x(i+m%np+m%vp%np_ghost(rank+1), :) = x(:, ix, iy, iz)
    end do
  else
    message(1) = 'mesh_create.mesh_create_xyz called without stencil'
    message(2) = 'and np_stencil when running in MPI mode.'
    call write_fatal(2)
  end if
#else
  ! When running serially those two are the same.
  m%np     = m%np_glob
  m%np_tot = m%np_tot_glob
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
  integer :: rank
  integer :: ierr
#endif

  call push_sub('mesh_create.mesh_get_vol_pp')

  f = M_ONE
  do i = 1, sb%dim
    f = f*mesh%h(i)
  end do

  allocate(mesh%vol_pp(mesh%np_tot))
  mesh%vol_pp(:) = f

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  call MPI_Comm_rank(mesh%vp%comm, rank, ierr)

  ! Do the inner points.
  do i = 1, mesh%np
    k = mesh%vp%local(mesh%vp%xlocal(rank+1)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi)
  end do
  ! Do the ghost points.
  do i = 1, mesh%vp%np_ghost(rank+1)
    k = mesh%vp%ghost(mesh%vp%xghost(rank+1)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i+mesh%np) = mesh%vol_pp(i+mesh%np)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np, :), chi)
  end do
  ! Do the boundary points.
  do i = 1, mesh%vp%np_bndry(rank+1)
    k = mesh%vp%bndry(mesh%vp%xbndry(rank+1)+i-1)
    chi(1:sb%dim) = mesh%Lxyz(k, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(rank+1)) = &
      mesh%vol_pp(i+mesh%np+mesh%vp%np_ghost(rank+1)) &
      *curvlinear_det_Jac(sb, geo, cv, mesh%x(i+mesh%np+mesh%vp%np_ghost(rank+1), :), chi)
  end do
#else
  do i = 1, mesh%np_tot
    chi(1:sb%dim) = mesh%Lxyz(i, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i, :), chi)
  end do
#endif

  call pop_sub()

end subroutine mesh_get_vol_pp


! this function takes care of the boundary conditions
! for a given x,y,z it returns the true index of the point

! WARNING: have to get rid of dir, otherwise will not work
integer function mesh_index(m, ix_, dir) result(index)
  type(mesh_type),      intent(in) :: m
  integer,              intent(in) :: ix_(:), dir

  integer :: i, ix(3)  ! ix has to go until 3, not sb%dim

  ix = 0
  ix(1:m%sb%dim) = ix_(1:m%sb%dim) ! make a local copy that we can change

  index = 1
  do i = 1, m%sb%dim
    if(ix(i) < m%nr(1, i)) then       ! first look left
      if(i <= m%sb%periodic_dim) then ! fold point
        ix(i) = ix(i) + abs(m%nr(2,i) - m%nr(1,i) + 1)
      else
        ix(i) = m%nr(1, i)
        index = 0
      end if
    else if(ix(i) > m%nr(2, i)) then  ! the same, but on the right
      if(i <= m%sb%periodic_dim) then
        ix(i) = ix(i) - abs(m%nr(2,i) - m%nr(1,i) + 1)
      else
        ix(i) = m%nr(2, i)
        index = 0
      end if
    end if
  end do

  if(index.ne.0) index = m%Lxyz_inv(ix(1), ix(2), ix(3))

  if(index==0.and.conf%boundary_zero_derivative) then
    do
      index = m%Lxyz_inv(ix(1), ix(2), ix(3))
      if(index.ne.0) exit
      ix(abs(dir)) = ix(abs(dir)) - sign(1, dir)
    end do
  end if


end function mesh_index
