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
  use curvilinear_m
  use datasets_m
  use global_m
  use hypercube_m
  use index_m
  use io_m
  use io_binary_m
  use loct_m
  use math_m
  use mesh_m
  use metis_m
  use messages_m
  use mpi_m
  use multicomm_m
  use parser_m
  use partition_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use stencil_star_m

  implicit none
  
  private
  public ::                      &
    mesh_partition,              &
    mesh_partition_boundaries,   &
    mesh_partition_write,        &
    mesh_partition_read,         &
    mesh_partition_write_info,   &
    mesh_partition_messages_debug

  integer, parameter ::  &
       METIS    = 1,     &
       PARMETIS = 2

  integer, parameter :: &
       RCB    = 1,      &
       GRAPH  = 2

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
  subroutine mesh_partition(mesh, lapl_stencil, part)
    type(mesh_t),    intent(inout)  :: mesh
    type(stencil_t), intent(in)  :: lapl_stencil
    integer,         intent(out) :: part(:) !< 1:mesh%np_part_global

    integer              :: iv, jp, inb
    integer              :: ix(1:MAX_DIM), jx(1:MAX_DIM)
    integer              :: ne             !< Number of edges.
    integer              :: nv             !< Number of vertices.
    !! Number of vertices (nv) is equal to number of
    !! points np_global and maximum number of edges (ne) is 2*mesh%sb%dim*np_global
    !! (there are a little fewer because points on the border have fewer
    !! than two neighbours per dimension).
    !! xadj has nv+1 entries because last entry contains the total
    !! number of edges.
    integer              :: npart          !< Number of partitions.
    integer              :: ipart          !< number of the current partition
    integer, allocatable :: xadj(:)        !< Indices of adjacency list in adjncy.
                                           !! Differents for each process if ParMETIS is used.
    integer, allocatable :: adjncy(:)      !< Adjacency lists.
                                           !! Only local part if ParMETIS is used.
    integer              :: iunit          !< For debug output to files.
    
#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
    integer, allocatable :: options(:)     !< Options to (Par)METIS.
    integer              :: edgecut        !< Number of edges cut by partitioning.
    REAL_SINGLE, allocatable :: tpwgts(:)  !< The fraction of vertex weight that should be distributed 
                                           !! to each sub-domain for each balance constraint
#endif
#if defined(HAVE_PARMETIS)
    integer :: ii
    integer, allocatable :: adjncy_local(:)!< Local part of adjacency list
    integer, allocatable :: xadj_local(:)  !< Local part of xadj
    integer, allocatable :: xadj_size_all(:)!< All the sizes of local matrices
    integer, allocatable :: part_local(:)  !< Output of the ParMETIS call
    integer, allocatable :: vtxdist(:)     !< Initial distribution of the points
#endif
    type(stencil_t) :: stencil
    integer :: stencil_to_use, default_method, method
    integer :: library
    integer, parameter   :: STAR = 1, LAPLACIAN = 2
    integer, allocatable :: istart(:), ifinal(:), lsize(:)

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
    !%Section Execution::Parallelization
    !%Description
    !% Decides which library to use to perform the mesh partition.
    !% By default ParMETIS is used when available, otherwise METIS is used.
    !%Option metis 1
    !% METIS library.
    !%Option parmetis 2
    !% (Experimental) Use ParMETIS libary to perform the mesh partition. 
    !% Only available if the code was compiled with ParMETIS support.
    !%End
    default = METIS
#ifdef HAVE_PARMETIS
    default = PARMETIS
#endif
    call parse_integer(datasets_check('MeshPartitionPackage'), default, library)

#if !defined(HAVE_METIS)
    if(library == METIS) then
      message(1) = 'METIS was requested, but Octopus was compiled without it.'
      call messages_fatal(1)
    end if
#endif
#if !defined(HAVE_PARMETIS)
    if(library == PARMETIS) then
      message(1) = 'PARMETIS was requested, but Octopus was compiled without it.'
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

    SAFE_ALLOCATE(istart(1:npart))
    SAFE_ALLOCATE(ifinal(1:npart))
    SAFE_ALLOCATE(lsize(1:npart))

    istart(1:npart) = 1
    ifinal(1:npart) = mesh%np_global
    lsize(1:npart) = mesh%np_global

    ! Shortcut (number of vertices).
    nv = lsize(ipart)
    SAFE_ALLOCATE(xadj(1:nv + 1))
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
      message(1) = "Info: Done writing mesh_graph.txt."
      call messages_info(1)
    end if


    select case(library)
    case(METIS)
#ifdef HAVE_METIS
      SAFE_ALLOCATE(options(1:40))
      SAFE_ALLOCATE(tpwgts(1:npart))
      
      ! The sum of all tpwgts elements has to be 1 and 
      ! we don`t care about the weights. So; 1/npart
      tpwgts(1:npart) = real(1.0, 4)/real(npart, 4)

      options = 0
      call oct_metis_setdefaultoptions(options(1)) ! is equal to: options = -1
      options(METIS_OPTION_NUMBERING) = 1 ! Fortran style: start counting from 1


      if(npart  <  8) then
        default_method = RCB
      else
        default_method = GRAPH
      end if

      !%Variable MeshPartition
      !%Type integer
      !%Section Execution::Parallelization
      !%Description
      !% When using METIS to perform the mesh partitioning, decides which
      !% algorithm is used. By default, <tt>graph</tt> partitioning 
      !% is used for 8 or more partitions, and <tt>rcb</tt> for fewer.
      !%Option rcb 1
      !% Recursive coordinate bisection partitioning.
      !%Option graph 2
      !% Graph partitioning (called 'k-way' by METIS).
      !%End
      call parse_integer(datasets_check('MeshPartition'), default_method, method)

      select case(method)
      case(RCB)
        message(1) = 'Info: Using METIS 5 multilevel recursive bisection to partition the mesh.'
        call messages_info(1)
        call oct_metis_partgraphrecursive(nv, 1, xadj(1), adjncy(1), npart, tpwgts(1), 1.01_4, options(1), edgecut, part(1))
      case(GRAPH)
        message(1) = 'Info: Using METIS 5 multilevel k-way algorithm to partition the mesh.'
        call messages_info(1)
        call oct_metis_partgraphkway(nv, 1, xadj(1), adjncy(1), npart, tpwgts(1), 1.01_4, options(1), edgecut, part(1))
      case default
        message(1) = 'Selected partition method is not available in METIS 5.'
        call messages_fatal(1)
      end select
      
      SAFE_DEALLOCATE_A(options)
#endif
    case(PARMETIS)
#ifdef HAVE_PARMETIS
      SAFE_ALLOCATE(options(1:3))
      options = 0 ! For the moment we use default options
      
      write(message(1),'(a)') 'Info: Using ParMETIS multilevel k-way algorithm to partition the mesh.'
      call messages_info(1)

      ! Compute vtxdist with multicomm_divide_range
      ! vtxdist has to be equal in all processes
      call multicomm_divide_range(mesh%np_global, npart, istart, ifinal, lsize)
      SAFE_ALLOCATE(vtxdist(1:npart+1))
      vtxdist(1) = 1
      do ii = 2, npart + 1
        vtxdist(ii) = vtxdist(ii-1) + lsize(ii-1)
      end do

      ! The sum of all tpwgts elements has to be 1 and 
      ! we don`t care about the weights. So; 1/npart
      SAFE_ALLOCATE(tpwgts(1:npart))
      tpwgts = real(1.0, 4)/real(npart, 4)
      
      ! Recalculate local adjncy and xadj from the global ones
      SAFE_ALLOCATE(xadj_local(1:lsize(ipart) + 1))
      do ii = 1, lsize(ipart) + 1
        xadj_local(ii) =  xadj(ii+vtxdist(ipart)-1) - xadj(vtxdist(ipart)) + 1
      end do
      SAFE_ALLOCATE(adjncy_local(1:xadj_local(lsize(ipart) + 1)))
      do ii = 1, xadj_local(lsize(ipart) + 1)
        adjncy_local(ii) = adjncy(xadj(vtxdist(ipart))+ii-1)
      end do
      
      ! Allocate output matrix
      SAFE_ALLOCATE(part_local(1:lsize(ipart) + 1)) ! Work with a local output
      ! Call to ParMETIS with:
      ! Fortran-style numbering. No imbalance tolerance
      call oct_parmetis_v3_partkway(vtxdist(1), xadj_local(1), adjncy_local(1), & 
           1, mesh%mpi_grp%size, tpwgts(1), 1.05, & 
           options(1), edgecut, part_local(1), mesh%mpi_grp%comm) 

      ASSERT(all(part_local(1:lsize(ipart))>0))
      ! Calculate sending sizes of all the processes
      SAFE_ALLOCATE(xadj_size_all(1:npart))
      do ii = 1, npart
        xadj_size_all(ii) = vtxdist(ii+1) - vtxdist(ii) 
      end do
      ! Adapt to MPI requirements of displacements
      vtxdist = vtxdist - 1  
      ! Gather local outputs to a global old-fashion one
      call MPI_Allgatherv(part_local(1), lsize(ipart), MPI_INTEGER, &
           part(1), xadj_size_all(1), vtxdist(1), MPI_INTEGER, &
           mesh%mpi_grp%comm, mpi_err)
        
      ! Deallocate matrices
      SAFE_DEALLOCATE_A(xadj_size_all)
      SAFE_DEALLOCATE_A(part_local)
      SAFE_DEALLOCATE_A(vtxdist)
      SAFE_DEALLOCATE_A(tpwgts)
      SAFE_DEALLOCATE_A(xadj_local)
      SAFE_DEALLOCATE_A(adjncy_local)
      SAFE_DEALLOCATE_A(options)      
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


  ! ----------------------------------------------------------------------
  subroutine mesh_partition_write_info(mesh, stencil, point_to_part)
    type(mesh_t),      intent(in)  :: mesh
    type(stencil_t),   intent(in)  :: stencil
    integer,           intent(in)  :: point_to_part(:)

    integer :: npart, npoints
    integer, allocatable :: nghost(:), nbound(:), nlocal(:), nneigh(:)
    FLOAT :: quality

    integer :: ip, ipcoords(1:MAX_DIM)
    integer, allocatable :: jpcoords(:, :), jp(:)
    integer :: istencil, ipart, jpart
    type(profile_t), save :: prof
    logical, allocatable :: is_a_neigh(:, :), gotit(:)
    FLOAT :: scal

    PUSH_SUB(mesh_partition_write_info)

    call profiling_in(prof, "MESH_PARTITION_WRITE_INFO")

    ! Build information about the mesh partition
    npart = mesh%mpi_grp%size
    npoints = mesh%np_part_global

    SAFE_ALLOCATE(nghost(1:npart))
    SAFE_ALLOCATE(nbound(1:npart))
    SAFE_ALLOCATE(nlocal(1:npart))
    SAFE_ALLOCATE(nneigh(1:npart))
    SAFE_ALLOCATE(is_a_neigh(1:npart, 1:npart))
    SAFE_ALLOCATE(gotit(1:mesh%np_part_global))
    SAFE_ALLOCATE(jpcoords(1:MAX_DIM, 1:stencil%size))
    SAFE_ALLOCATE(jp(1:stencil%size))

    is_a_neigh = .false.
    nghost = 0
    nbound = 0
    nlocal = 0
    nneigh = 0

    do ipart = 1, npart
      gotit = .false.
      do ip = 1, mesh%np_global
        if(ipart /= point_to_part(ip)) cycle

        INCR(nlocal(ipart), 1)
        call index_to_coords(mesh%idx, mesh%sb%dim, ip, ipcoords)
        
        do istencil = 1, stencil%size
          jpcoords(:, istencil) = ipcoords + stencil%points(:, istencil)
        end do
        
        call index_from_coords_vec(mesh%idx, mesh%sb%dim, stencil%size, jpcoords, jp)
        
        do istencil = 1, stencil%size
          if(stencil%center == istencil) cycle

          if(.not. gotit(jp(istencil))) then
            jpart = point_to_part(jp(istencil))
         
            if(jpart /= ipart) then
              INCR(nghost(ipart), 1)
              is_a_neigh(ipart, jpart) = .true.
            else if(jp(istencil) > mesh%np_global) then
              INCR(nbound(ipart), 1)
            end if
            
            gotit(jp(istencil)) = .true.
          end if
          
        end do
        
      end do
    end do

    forall(ipart = 1:npart)
      nneigh(ipart) = count(is_a_neigh(ipart, 1:npart))
    end forall

    SAFE_DEALLOCATE_A(is_a_neigh)
    SAFE_DEALLOCATE_A(gotit)
    SAFE_DEALLOCATE_A(jpcoords)
    SAFE_DEALLOCATE_A(jp)

    ! Calculate partition quality
    scal = real(npart, REAL_PRECISION)/npoints

    quality = M_ZERO

    quality = quality + (maxval(nlocal) - minval(nlocal))**3
    quality = quality + (sum(TOFLOAT(nghost)**2))

    quality = M_ONE/(M_ONE + quality)


    ! Write information about the partition
    message(1) = &
      'Info: Mesh partition:'
    message(2) = ''
    call messages_info(2)

    write(message(1),'(a,e16.6)') &
      '      Partition quality:', quality
    message(2) = ''
    call messages_info(2)

    write(message(1),'(a)') &
      '                 Neighbours         Ghost points'
    write(message(2),'(a,i5,a,i10)') &
      '      Average  :      ', sum(nneigh)/npart, '           ', sum(nghost)/npart
    write(message(3),'(a,i5,a,i10)') &
      '      Minimum  :      ', minval(nneigh),    '           ', minval(nghost)
    write(message(4),'(a,i5,a,i10)') &
      '      Maximum  :      ', maxval(nneigh),    '           ', maxval(nghost)
    message(5) = ''
    call messages_info(5)

    do ipart = 1, npart
      write(message(1),'(a,i5)')  &
        '      Nodes in domain-group  ', ipart
      write(message(2),'(a,i10,a,i10)') &
        '        Neighbours     :', nneigh(ipart), &
        '        Local points    :', nlocal(ipart)
      write(message(3),'(a,i10,a,i10)') &
        '        Ghost points   :', nghost(ipart), &
        '        Boundary points :', nbound(ipart)
      call messages_info(3)
    end do

    message(1) = ''
    call messages_info(1)

    SAFE_DEALLOCATE_A(nghost)
    SAFE_DEALLOCATE_A(nbound)
    SAFE_DEALLOCATE_A(nlocal)
    SAFE_DEALLOCATE_A(nneigh)

    call profiling_out(prof)

    POP_SUB(mesh_partition_write_info)
  end subroutine mesh_partition_write_info

  ! ----------------------------------------------------
  subroutine mesh_partition_messages_debug(mesh)
    type(mesh_t),    intent(in)    :: mesh
    
    integer              :: ii, jj         ! Counter.
    integer              :: iunit          ! For debug output to files.
    character(len=3)     :: filenum

    if(.not. in_debug_mode) return

    PUSH_SUB(mesh_partition_messages_debug)

    call io_mkdir('debug/mesh_partition')

    ! Debug output. Write points of each partition in a different file.
    write(filenum, '(i3.3)') mesh%mpi_grp%rank+1

    iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
      action='write')
    do jj = 1, mesh%np_global
      if(mesh%vp%part_vec(jj)  ==  mesh%mpi_grp%rank+1) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
    end do
    call io_close(iunit)

    iunit = io_open('debug/mesh_partition/mesh_partition_all.'//filenum, &
      action='write')
    do jj = 1, mesh%np_part_global
      if(mesh%vp%part_vec(jj)  ==  mesh%mpi_grp%rank+1) write(iunit, '(i8,3f18.8)') jj, mesh_x_global(mesh, jj)
    end do
    call io_close(iunit)

    if(mpi_grp_is_root(mpi_world)) then
      ! Write points from enlargement to file with number p+1.
      write(filenum, '(i3.3)') mesh%mpi_grp%size+1
      iunit = io_open('debug/mesh_partition/mesh_partition.'//filenum, &
        action='write')
      do ii = mesh%np_global+1, mesh%np_part_global
        write(iunit, '(i8,3f18.8)') ii, mesh_x_global(mesh, ii)
      end do
      call io_close(iunit)
    end if

    POP_SUB(mesh_partition_messages_debug)

  end subroutine mesh_partition_messages_debug

end module mesh_partition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
