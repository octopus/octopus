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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module mesh_partition_m
  use cube_m
  use curvilinear_m
  use datasets_m
  use geometry_m
  use global_m
  use hypercube_m
  use index_m
  use io_m
  use io_binary_m
  use loct_m
  use math_m
  use mesh_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use partition_m
  use partitioner_m
  use pfft_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use stencil_star_m
  use zoltan_m

  implicit none
  
  private
  public ::                      &
    mesh_partition,              &
    mesh_partition_boundaries,   &
    mesh_partition_write,        &
    mesh_partition_read,         &
    mesh_partition_messages_debug

contains
  
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
  subroutine mesh_partition(mesh, lapl_stencil, part, cube)
    type(mesh_t),    intent(inout)  :: mesh
    type(stencil_t), intent(in)  :: lapl_stencil
    integer,         intent(out) :: part(:) ! 1:mesh%np_part_global
    type(cube_t), optional, intent(in)  :: cube

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
    integer :: ip
    integer :: im, ii, nn
    integer :: stencil_to_use, default_method, method
    integer :: library
    integer, parameter   :: METIS = 2, ZOLTAN = 3, GA = 4, PFFT_PART = 5
    integer, parameter   :: STAR = 1, LAPLACIAN = 2
    integer, allocatable :: istart(:), ifinal(:), lsize(:)
    FLOAT, allocatable   :: xglobal(:, :)
    type(partition_t)    :: partition
    type(partitioner_t)  :: partitioner

    type(profile_t), save :: prof
    integer :: default

    call profiling_in(prof, "MESH_PARTITION")
    PUSH_SUB(mesh_partition)

    if(mesh%np_global == 0) then
      message(1) = 'The mesh is empty and cannot be partitioned.'
      call messages_fatal(1)
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
    !%Option ga 4
    !% (Experimental) Genetic-algorithm optimization of the grid partition.
    !%Option pfft_part 5
    !% (Experimental) Use PFFT to perform the mesh partition.
    !%End
    default = ZOLTAN
#ifdef HAVE_METIS
    default = METIS
#endif
    call parse_integer(datasets_check('MeshPartitionPackage'), default, library)

    if(library == GA) call messages_experimental('Genetic algorithm mesh partition')
    if(library == PFFT_PART) call messages_experimental('PFFT mesh partition')

    if (library == PFFT_PART) then
      ASSERT(present(cube))
    end if

#ifndef HAVE_METIS
    if(library == METIS) then
      message(1) = 'Error: METIS was requested, but Octopus was compiled without it.'
      call messages_fatal(1)
    end if
#endif

#ifndef HAVE_PFFT
    if(library == PFFT_PART) then
      message(1) = 'Error: PFFT was requested, but Octopus was compiled without it.'
      call messages_fatal(1)
    end if
#endif
    mesh%partition_library = library

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

    SAFE_ALLOCATE(istart(1:npart))
    SAFE_ALLOCATE(ifinal(1:npart))
    SAFE_ALLOCATE(lsize(1:npart))

    select case(library)
    case(METIS)

      istart(1:npart) = 1
      ifinal(1:npart) = mesh%np_global
      lsize(1:npart) = mesh%np_global

    case(ZOLTAN)

      ! If we use Zoltan, we divide the space in a basic way, to balance
      ! the memory for the graph. 
      call multicomm_divide_range(mesh%np_global, npart, istart, ifinal, lsize)

      do ii = 1, npart
        part(istart(ii):ifinal(ii)) = ii
      end do

    end select

    if(library /= GA .and. library /= PFFT_PART) then
      ! Shortcut (number of vertices).
      nv = lsize(ipart)
      SAFE_ALLOCATE(xadj(1:nv + 1))

      if(library == METIS .or. .not. zoltan_method_is_geometric(method)) then !calculate the graphs
        SAFE_ALLOCATE(adjncy(1:(stencil%size - 1)*nv))

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
            if(jp == stencil%center) cycle

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
          call messages_info(5)
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
        call messages_info(1)
        call oct_metis_part_graph_recursive(nv, xadj, adjncy, &
          0, 0, 0, 1, npart, options, edgecut, part)
      case(GRAPH)
        message(1) = 'Info: Using METIS multilevel k-way algorithm to partition the mesh.'
        call messages_info(1)
        call oct_metis_part_graph_kway(nv, xadj, adjncy, &
          0, 0, 0, 1, npart, options, edgecut, part)
      case default
        message(1) = 'Error: Selected partition method is not available in METIS.'
        call messages_fatal(1)
      end select
#endif
    case(ZOLTAN)

      call zoltan_method_info(method)

      SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM))

      do ip = 1, mesh%np_part_global
        xglobal(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip)
      end do

#ifdef HAVE_MPI
      !assign all points to one node
      call zoltan_partition(method, mesh%sb%dim, mesh%np_global, mesh%np_part_global, &
        xglobal(1, 1), istart(ipart), xadj(1), adjncy(1), ipart, part(1), mesh%mpi_grp%comm)
#endif

      SAFE_DEALLOCATE_A(xglobal)

      ! we use xadj as a buffer
      xadj(1:lsize(ipart)) = part(istart(ipart):ifinal(ipart))

      ASSERT(all(xadj(1:lsize(ipart)) > 0))

      part(1:mesh%np_global) = 0 ! so we catch non-initialized values

      ! convert start to C notation
      istart = istart - 1

#ifdef HAVE_MPI
      ! we collect part from all processors
      call MPI_Allgatherv(xadj(1), lsize(ipart), MPI_INTEGER, part(1), &
        lsize(1), istart(1), MPI_INTEGER, mesh%mpi_grp%comm, mpi_err)
#endif

    case(GA)

      if(mpi_grp_is_root(mesh%mpi_grp)) then

        call partition_init(partition, mesh)

        call partitioner_init(partitioner)
        call partitioner_perform(partitioner, mesh, stencil, partition)
        call partitioner_end(partitioner)

        part = partition%point_to_part
        call partition_end(partition)
      end if

#ifdef HAVE_MPI
      call MPI_Bcast(part(1), mesh%np_part_global, MPI_INTEGER, 0, mesh%mpi_grp%comm, mpi_err)
#endif

    case(PFFT_PART)
#ifdef HAVE_PFFT
      part = 0
      do ip = 1, mesh%np_global
        call index_to_coords(mesh%idx, mesh%sb%dim, ip, ix(1:3))
        ix(1:3) = ix(1:3) + cube%center(1:3)
        part(ip) = cube%part(ix(1), ix(2), ix(3))
      end do
#endif
    end select

    SAFE_DEALLOCATE_A(istart)
    SAFE_DEALLOCATE_A(ifinal)
    SAFE_DEALLOCATE_A(lsize)
    SAFE_DEALLOCATE_A(adjncy)

    ASSERT(all(part(1:mesh%np_global) > 0))
    ASSERT(all(part(1:mesh%np_global) <= npart))

    call stencil_end(stencil)
    POP_SUB(mesh_partition)
    call profiling_out(prof)

  end subroutine mesh_partition

  ! --------------------------------------------------------

  subroutine mesh_partition_boundaries(mesh, stencil, part)
    type(mesh_t),    intent(in)    :: mesh
    type(stencil_t), intent(in)    :: stencil
    integer,         intent(inout) :: part(:)

    integer              :: ii, jj         ! Counter.
    integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
    integer              :: npart
    integer, allocatable :: votes(:), bps(:)
    logical, allocatable :: winner(:)
    integer :: ip, maxvotes

    PUSH_SUB(mesh_partition_boundaries)

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

    POP_SUB(mesh_partition_boundaries)
  end subroutine mesh_partition_boundaries

  ! ----------------------------------------------------

  subroutine mesh_partition_write(mesh, part)
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(in)    :: part(:)
    
    character(len=6) :: numstring
    integer          :: ierr

    PUSH_SUB(mesh_partition_write)

    ! here we assume that all nodes have the same partition
    if(mpi_grp_is_root(mpi_world)) then

      write(numstring, '(i6.6)') mesh%mpi_grp%size

      call io_mkdir('restart', is_tmp = .true.)
      call io_mkdir('restart/partition', is_tmp = .true.)
      call mesh_write_fingerprint(mesh, 'restart/partition/grid_'//trim(numstring))
      call io_binary_write('restart/partition/partition_'//trim(numstring)//'.obf', mesh%np_part_global, part, ierr)
    end if

    POP_SUB(mesh_partition_write)
  end subroutine mesh_partition_write

  ! ----------------------------------------------------

  subroutine mesh_partition_read(mesh, part, ierr)
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(out)   :: part(:)
    integer,         intent(out)   :: ierr
    
    character(len=6) :: numstring
    integer :: read_np_part

    PUSH_SUB(mesh_partition_read)

    if(mpi_grp_is_root(mesh%mpi_grp)) then

      write(numstring, '(i6.6)') mesh%mpi_grp%size

      call mesh_read_fingerprint(mesh, 'restart/partition/grid_'//trim(numstring), read_np_part, ierr)

      if (ierr == 0) then
        call io_binary_read('restart/partition/partition_'//trim(numstring)//'.obf', mesh%np_part_global, part, ierr)
      else
        ierr = -1
      end if

      if(ierr == 0) then
        message(1) = "Info: Found a compatible mesh partition."
      else
        message(1) = "Info: Could not find a compatible mesh partition."
      end if
      call messages_info(1)
    end if
    
    !broadcast the result
#ifdef HAVE_MPI
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, mesh%mpi_grp%comm, mpi_err)
    if (ierr == 0) call MPI_Bcast(part(1), mesh%np_part_global, MPI_INTEGER, 0, mesh%mpi_grp%comm, mpi_err)
#endif

    POP_SUB(mesh_partition_read)

  end subroutine mesh_partition_read

  ! ----------------------------------------------------

  subroutine mesh_partition_messages_debug(mesh, part)
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(inout) :: part(:)

    integer              :: ii, jj         ! Counter.
    integer              :: npart
    integer              :: iunit          ! For debug output to files.
    character(len=3)     :: filenum

    PUSH_SUB(mesh_partition_messages_debug)

    if(in_debug_mode .and. mpi_grp_is_root(mpi_world)) then

      call io_mkdir('debug/mesh_partition')

      npart = mesh%mpi_grp%size

      ! Debug output. Write points of each partition in a different file.
      do ii = 1, npart

        write(filenum, '(i3.3)') ii

        iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
          action='write')
        do jj = 1, mesh%np_global
          if(part(jj) .eq. ii) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
        end do
        call io_close(iunit)

        iunit = io_open('debug/mesh_partition/mesh_partition_all.'//filenum, &
          action='write')
        do jj = 1, mesh%np_part_global
          if(part(jj) .eq. ii) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
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

#ifdef HAVE_MPI
    call MPI_Barrier(mesh%mpi_grp%comm, mpi_err)
#endif

    POP_SUB(mesh_partition_messages_debug)

  end subroutine mesh_partition_messages_debug

end module mesh_partition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
