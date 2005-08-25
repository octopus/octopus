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

#if defined(HAVE_METIS)
! Calls mesh_partition_init with the number of partitions
! equal to the number of processors stored in the global mpiv.
! This is what is usually needed.
subroutine mesh_partition_init_default(m, Lxyz_tmp)
   type(mesh_type), intent(inout) :: m
  ! FIXME: Works currently only with specified indices, why?
  integer,         intent(in)    :: Lxyz_tmp(m%nr(1,1):m%nr(2,1), &
                                             m%nr(1,2):m%nr(2,2), &
                                             m%nr(1,3):m%nr(2,3))
  
  call push_sub('mesh_create.mesh_partition_init_default')

  call mesh_partition_init(m, Lxyz_tmp, mpiv%numprocs)

  call pop_sub()
  
end subroutine mesh_partition_init_default


! Converts the mesh given by grid points into a graph. Each point is
! a vertex in the graph and closest neighbours are connected by an
! edge (at max. 6 in 3D and 4 in 2D, less at the borders).
! Asserts working in 2D or 3D (but would do in 1D as well).
! Then calls METIS to get p partitions.
! Stored the mapping point no. -> partition no. into
! m%part. m%part is allocated along the way.
! (mesh_partition_end should not be forgotten to call later.)
subroutine mesh_partition_init(m, Lxyz_tmp, p)
  type(mesh_type), intent(inout) :: m
  ! FIXME: Works currently only with specified indices, why?
  integer,         intent(in)    :: Lxyz_tmp(m%nr(1,1):m%nr(2,1), &
                                             m%nr(1,2):m%nr(2,2), &
                                             m%nr(1,3):m%nr(2,3))
  integer,         intent(in)    :: p 

  integer              :: i              ! Counter.
  integer              :: ix, iy, iz     ! Counters to iterate over grid.
  integer              :: ne             ! Number of edges.
  integer              :: nv             ! Number of vertices.
  integer              :: d(3, 6)        ! Directions of neighbour points.
  integer              :: edgecut        ! Number of edges cut by partitioning.
  ! Number of vertices (nv) is equal to number of
  ! points np and maximum number of edges (ne) is 6*np (there
  ! are a little less because points on the border have less
  ! than 6 neighbours).
  ! xadj has nv+1 entries because last entry contains the total
  ! number of edges.
  integer              :: xadj(m%np+1)   ! Indices of adjacency list in adjncy.
  integer              :: adjncy(6*m%np) ! Adjacency lists.
  integer              :: options(5)     ! Options to METIS.
  integer, allocatable :: ppp(:)         ! Points per partition.
  integer, pointer     :: part(:)        ! Mapping of nodes to partitions.
#ifdef DEBUG
  integer              :: j              ! Counter.
  integer              :: iunit          ! For debug output to files.
  character(len=3)     :: filenum
#endif

  call push_sub('mesh_create.metis_partition')
  
  ASSERT(m%sb%dim==2.or.m%sb%dim==3)

  options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

  ! Get space for partitioning.
  allocate(part(m%np))
  part = 1

  ! If p=1, we can exit. All points are mapped to
  ! partition 1 (s. a.).
  if(p.eq.1) return

  ! Set directions of possible neighbours.
  d(:, 1) = (/ 1,  0,  0/)
  d(:, 2) = (/ 0,  1,  0/)
  d(:, 3) = (/ 0,  0,  1/)
  d(:, 4) = (/-1,  0,  0/)
  d(:, 5) = (/ 0, -1,  0/)
  d(:, 6) = (/ 0,  0, -1/)

  ! Create graph with each point being
  ! represenetd by a vertice and edges between
  ! neighboured points.
  nv = 0
  ne = 1
  do ix = m%nr(1,1), m%nr(2,1)
    do iy = m%nr(1,2), m%nr(2,2)
      do iz = m%nr(1,3), max(0, m%nr(2,3)) ! max is necessary when running in 2 dimensions
        ! Only if the current point (ix, iy, iz) is inside the box,
        ! neighbours are relevant.
        if(Lxyz_tmp(ix, iy, iz).eq.1) then
          nv       = nv+1
          xadj(nv) = ne
          ! Check all possible neighbours.
          do i = 1, 6
            if(Lxyz_tmp(ix+d(1, i), iy+d(2, i), iz+d(3, i)).eq.1) then
              adjncy(ne) = m%Lxyz_inv(ix+d(1, i), iy+d(2, i), iz+d(3, i))
              ne         = ne+1
            end if
          end do
        end if
      end do
    end do
  end do
  ne         = ne-1 ! We start with ne=1 for simplicity. This is off by one
                    ! in the end --> -1.
  xadj(nv+1) = ne+1 ! Set number of edges plus 1 as last index.
                    ! The reason is: neighbours of node i are stored
                    ! in adjncy(xadj(i):xadj(i+1)-1). Setting the last
                    ! index as mentioned makes special handling of
                    ! last element unnecessary.


#ifdef DEBUG
  ! DEBUG output. Write graph to file mesh_graph.txt.
  message(1) = 'Adjacency lists of the graph representing the grid'
  message(2) = 'are stored in debug/mesh_partition/mesh_graph.txt.'
  message(3) = 'Compatible with METIS programs pmetis and kmetis.'
  message(4) = 'First line contains number of vertices and edges.'
  message(5) = 'Edges are not directed and appear twice in the lists.'
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
      do j = 1, m%np
        if(part(j).eq.i) write(iunit, '(i8,3f12.8)') j, m%x(j,:)
      end do
      call io_close(iunit)
    end do
  end if
#endif

  ! Store partitioning in mesh.
  m%npart = p
  m%part => part

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


subroutine mesh_create_xyz(sb, m, cv, geo)
  type(simul_box_type),  intent(in)    :: sb
  type(mesh_type),       intent(inout) :: m
  type(curvlinear_type), intent(in)    :: cv
  type(geometry_type),   intent(in)    :: geo

  integer :: i, j, k, il, ix, iy, iz
  integer, allocatable :: Lxyz_tmp(:,:,:)
  FLOAT,   allocatable :: x(:,:,:,:)
  FLOAT :: chi(3)

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
  m%np_tot = il
  allocate(m%Lxyz(m%np_tot, 3), m%x(m%np_tot, 3))

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
          m%x(il,:) = x(:,ix,iy,iz)
        end if
      enddo
    enddo
  enddo
  m%np = il

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
          m%x(il,:) = x(:,ix,iy,iz)
        end if
      end do
    end do
  end do

  call mesh_get_vol_pp(sb, geo, cv, m)

#if defined(HAVE_MPI) && defined(HAVE_METIS)
  call mesh_partition_init_default(m, Lxyz_tmp)
#endif

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
  FLOAT :: f, chi(sb%dim)

  f = M_ONE
  do i = 1, sb%dim
    f = f*mesh%h(i)
  end do

  allocate(mesh%vol_pp(mesh%np_tot))
  mesh%vol_pp(:) = f

  do i = 1, mesh%np_tot
    chi(1:sb%dim) = mesh%Lxyz(i, 1:sb%dim) * mesh%h(1:sb%dim)
    mesh%vol_pp(i) = mesh%vol_pp(i)*curvlinear_det_Jac(sb, geo, cv, mesh%x(i,:), chi)
  end do

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
