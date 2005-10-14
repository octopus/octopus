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


#include "global.h"

! Notes ragarding the multi communicator part.
!
! This part is for combined parallelization in indices and
! domains.
! Given a number n_index of indices and their ranges
! index_range(1:n_index) communicators for domains and
! indices are created.
!
! Example
! Given to indices j, k with ranges
! index_range(j) = 3, index_range(k) = 4
! and n_node = 48 nodes, we would put
! n_domain_node = n_node/n_domain with
! n_domain = prod(a=j, k)[index_range(a)] = 4
! nodes into each domain parallelization.
! To do collective operations (like sums) over j, k respectively
! n_index_comm = N(n_index) = 7 communicators are introduced
! (for the definition of N see footnote (1)).
! Each of these communicators contains only the root nodes of
! all domain parallelizations participating in this particular
! index. The reason for this is that only root nodes know
! the complete functions (after vec_gather).
! 
! The example above can visualized as follows:
! 
!   j
!  --->
!  |   (1)  (2)  (3)
! k|   (4)  (5)  (6)
!  |   (7)  (8)  (9)
!  V  (10) (11) (12)
!
! (n) are domain parallelizations of four nodes.
! The communicators for j are as follows:
! index_comm(k=1, j) = {1, 2, 3}
! index_comm(k=2, j) = {4, 5, 6}
! index_comm(k=3, j) = {7, 8, 9}
! index_comm(k=4, j) = {10, 11, 12}
! 
! For k they look like this:
! index_comm(j=1, k) = {1, 4, 7, 10}
! index_comm(j=2, k) = {2, 5, 8, 11}
! index_comm(j=3, k) = {3, 6, 9, 12}
!
! {p} means that the root node of domain parallelization
! p is member of the denoted communicator.
!
! In general the communicators for index a(i) out of
! a(1), ..., a(n_index) is addressed by specifying all
! other indices x, y, ... except i:
! index_comm(a=(x, y, ...), i)
!
! For generality, all index communicators are stored in a
! vector and the adressing is done by a function index_comm_number(a, i).
! a(1:n_index) specifies the values for all indices except i
! (the value of a(i) is irrelevant) and i is the number of
! the requested index:
! index_comm_number(a, i) == offset + position
! WHERE
! offset   == sum(x=1, ..., i-1; index_range(x)>1)
!             [prod(y=1, ..., n_index; y|=x)[index_range(y)]]
! position == sum(x=1, ..., n_index; x|=i; index_range(x)>1)[(a(x)-1)*
!             prod(y=x+1, ..., n_index; y|=i)[index_range(y)]] + 1
!
! Note that only indices with ranges greater than one are summed up.
! An index with range one is equal to no parallelization in this index.
!
! For each index i the rest of the indices form a
! n_index-1 dimensional array. The offset of this array in
! index_comm is offset in the above function. It is the number of
! communicators belonging to the indices 1, ..., i-1.
! The position in the array is computed as for any other
! n-dimensional array (for generalities sake it has to be done by hand
! and not by the Fortran compiler).
! With this function, the j communicator for k=2 can be accessed with
! index_comm(index_comm_number((/0, 2/), 1))
! (with j being the first and k the second index).
!
! Some more stripped down cases:
! (*) Only index parallelization (e. g. k-points and states j):
!     There are domain communicators with only one node, i, e.
!     no domain parallelization.
! (*) Only domain parallelization: There would not no index
!     communicators and just one domain parallelization:
!     n_domain = 1, domain_comm(1) = MPI_COMM_WORLD.
!
!
! ---------- 
! (1) N(1)   = 1
!     N(i+1) = N(i)*index_range(i+1) + prod(a=1, ..., i)[index_range(a)]
!
! (2) prod(x=1, ..., n)[f(x)] = f(1)* ... *f(x)
!

module mpi_mod
#if defined(HAVE_MPI)
  use varinfo
  use global
  use messages
  use lib_oct
  use lib_oct_parser
  use syslabels
#if !defined(MPI_H)
  use mpi
#endif


  implicit none


#if defined(MPI_H)
# include "mpif.h"
#endif


  private

  public ::                              &
       MPI_Debug_Statistics,             &
       MPI_Debug_IN, MPI_Debug_OUT,      &
       multicomm_type,                   &
       multicomm_init, multicomm_end,    &
       multicomm_strategy

  public ::                                                        &
       TSD_MPI_Barrier,    TSZ_MPI_Barrier,    TSI_MPI_Barrier,    &
       TSD_MPI_Scatterv,   TSZ_MPI_Scatterv,   TSI_MPI_Scatterv,   &
       TSD_MPI_Gatherv,    TSZ_MPI_Gatherv,    TSI_MPI_Gatherv,    &
       TSD_MPI_Alltoallv,  TSZ_MPI_Alltoallv,  TSI_MPI_Alltoallv,  &
       TSD_MPI_Allgatherv, TSZ_MPI_Allgatherv, TSI_MPI_Allgatherv, &
       TSD_MPI_Bcast,      TSZ_MPI_Bcast,      TSI_MPI_Bcast,      &
       TSD_MPI_Allreduce,  TSZ_MPI_Allreduce,  TSI_MPI_Allreduce

  integer, public, parameter ::  &
       C_MPI_BARRIER    = 1,     &
       C_MPI_SCATTERV   = 2,     &
       C_MPI_GATHERV    = 3,     &
       C_MPI_ALLTOALLV  = 4,     &
       C_MPI_ALLGATHERV = 5,     &
       C_MPI_BCAST      = 6,     &
       C_MPI_ALLREDUCE  = 7

  character(len=15), dimension(C_MPI_ALLREDUCE), public :: mpi_rlabel = &
       (/                &
       'MPI_BARRIER   ', &
       'MPI_SCATTERV  ', &
       'MPI_GATHERV   ', & 
       'MPI_ALLTOALLV ', &
       'MPI_ALLGATHERV', &
       'MPI_BCAST     ', &
       'MPI_ALLREDUCE '  &
       /)       

  integer, public :: call_counter(C_MPI_BARRIER:C_MPI_ALLREDUCE) = 0
  integer, public :: sec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)    = 0
  integer, public :: usec_accum(C_MPI_BARRIER:C_MPI_ALLREDUCE)   = 0

  integer, private :: sec_in, usec_in


  ! Stores all communicators and groups.
  type multicomm_type
    integer          :: n_node               ! Total number of nodes.
    integer          :: n_domain_node        ! Number of nodes per domain
    integer          :: n_index              ! Number of parallel indices.
    integer          :: n_domain_comm        ! Number of domain communicators.
    integer          :: n_index_comm         ! Number of index communicators.

    integer          :: par_strategy         ! What kind of parallelization strategy should we use?
    logical          :: use_domain_par       ! Should we use domain parallelization?

    integer, pointer :: index_range(:)       ! Range of index i is
                                             ! 1, ..., index_range(i).
    integer, pointer :: index_weight(:)      ! How communication intensive
                                             ! are the indices (big weight ->
                                             ! lots of communication)?

    integer, pointer :: domain_comm(:)       ! Domain communicators.
    integer, pointer :: domain_root(:)       ! Ranks of roots from every domain.

    integer, pointer :: index_comm(:)        ! Index communicators
    integer, pointer :: index_contraction(:) ! Maps index-combinations to
                                             ! to domain parallelizations.

    ! mapping for comm_world ranks
    integer, pointer :: domain_comm_of_node(:) ! Domain communicator that node belongs to
  end type multicomm_type


  ! possible parallelization strategies
  integer, public, parameter ::          &
       P_STRATEGY_SERIAL          =   0, & ! single domain, all states, kpoints on a single processor
       P_STRATEGY_ONLY_DOMAINS    =   1, & ! only parallelization in domains
       P_STRATEGY_ONLY_KPOINTS    =  10, & ! no parallelization in domains, only in kpoints
       P_STRATEGY_KPOINTS_DOMAINS =  11, & ! combined parallelization in kpoints and domains
       P_STRATEGY_ONLY_STATES     = 100, & ! only parallelization in state indices
       P_STRATEGY_STATES_DOMAINS  = 101, & ! combined parallelization in states and domains
       P_STRATEGY_STATES_KPOINTS  = 110, & ! no parallelization in domains, only in states and kpoints
       P_STRATEGY_FULL            = 111    ! parallelization in states, kpoints and domains


contains

  ! =========================================================
  ! Routines for the multicommunicator part
  ! =========================================================
  

  ! decide which parallelization strategy we should use (this routine is at the moment 
  ! specialized to the case of two indices: states and kpoints
  subroutine multicomm_strategy(index_range, mc)
    integer, intent(in) :: index_range(:)
    type(multicomm_type), intent(inout) :: mc

    integer :: n_kpoints, n_states

    call push_sub('mpi.multicomm_strategy')

    ! hardwire this here
    n_kpoints = index_range(1)
    n_states  = index_range(2)

    ! the user should be able to select a strategy from the following combinations:
    ! 
    !  1)  000 - serial          (single domain, all states, kpoints on a single processor)
    !  2)  001 - only_domains    (only parallelization in domains)
    !  3)  010 - only_kpoints    (no parallelization in domains, only in kpoints)
    !  4)  011 - kpoints_domains (combined parallelization in kpoints and domains)
    !  5)  100 - only_states     (only parallelization in state indices)
    !  6)  101 - states_domains  (combined parallelization in states and domains)
    !  7)  110 - states_kpoints  (no parallelization in domains, only in states and kpoints)
    !  8)  111 - full            (parallelization in states, kpoints and domains)   
    !
    ! if a selected mode is available depends of course on the calc_mode.

    !%Variable ParallelizationStrategy
    !%Type integer
    !%Section 1 Generalities
    !%Description
    !% Specifies what kind of parallelization strategy octopus should use
    !%Option serial 0
    !% Octopus will run in serial.
    !%Option only_domains 1
    !% Octopus will run parallel in domains.
    !%Option only_kpoints 10
    !% Octopus will run parallel in k-points.
    !%Option kpoints_domains 11
    !% Octopus will run parallel in k-points and domains.
    !%Option only_states 100
    !% Octopus will run parallel in states.
    !%Option states_domains 101
    !% Octopus will run parallel in states and domians.
    !%Option states_kpoints 110
    !% Octopus will run parallel in states and k-points.
    !%Option full 111
    !% Octopus will run parallel in states, k-points and domains.
    !%End
    call loct_parse_int(check_inp('ParallelizationStrategy'),  &
         P_STRATEGY_SERIAL, mc%par_strategy)

    mc%index_weight  = 1

    select case(mc%par_strategy)
    case(P_STRATEGY_SERIAL)
       mc%use_domain_par = .false.
       mc%index_range(1) = 1 
       mc%index_range(2) = 1 
       message(1) = 'Info: Octopus will run in serial.'
    case(P_STRATEGY_ONLY_DOMAINS)
       mc%use_domain_par = .true.             
       mc%index_range(1) = 1 
       mc%index_range(2) = 1 
       message(1) = 'Info: Octopus will run parallel in domains.'
    case(P_STRATEGY_ONLY_KPOINTS)
       mc%use_domain_par = .false.           
       mc%index_range(1) = n_kpoints
       mc%index_range(2) = 1
       message(1) = 'Info: Octopus will run parallel in k-points.'
    case(P_STRATEGY_KPOINTS_DOMAINS)
       mc%use_domain_par = .true.
       mc%index_range(1) = n_kpoints
       mc%index_range(2) = 1
       message(1) = 'Info: Octopus will run parallel in k-points and domains.'
    case(P_STRATEGY_ONLY_STATES)
       mc%use_domain_par = .false.
       mc%index_range(1) = 1
       mc%index_range(2) = n_states             
       message(1) = 'Info: Octopus will run parallel in states.'
    case(P_STRATEGY_STATES_DOMAINS)
       mc%use_domain_par = .true.
       mc%index_range(1) = 1
       mc%index_range(2) = n_states
       message(1) = 'Info: Octopus will run parallel in states and domains.'
    case(P_STRATEGY_STATES_KPOINTS)
       mc%use_domain_par  = .false.
       mc%index_range(1)  = n_kpoints
       mc%index_range(2)  = n_states
       mc%index_weight(1) = 1
       mc%index_weight(2) = 100000
       message(1) = 'Info: Octopus will run parallel in states and k-points.'
    case(P_STRATEGY_FULL)
       mc%use_domain_par  = .true.
       mc%index_range(1)  = n_kpoints
       mc%index_range(2)  = n_states
       mc%index_weight(1) = 1
       mc%index_weight(2) = 100000
       message(1) = 'Info: Octopus will run parallel in states, k-points and domains.'
    case default
       call input_error('ParallelizationStrategy')
    end select
    call write_info(1)

    ! currently we do not have block diagonalizers in octopus, i.e. for the ground state
    ! the k-points are proper multicomm indices, but states are not.
    ! all strategies involving states are therefore not possible for calc_mode = M_GS
    if (calc_mode.eq.1) then ! M_GS    
       mc%index_range(2) = 1
       if (mc%par_strategy.ge.P_STRATEGY_ONLY_STATES) then
          message(1) = 'Warning: Parallelization in states currently not available for CalcluationMode = gs.'
          message(2) = '         Disabling parallelization in states.'
          call write_warning(2)  
       endif
    endif

    ! so far we allow multicommunicators only for gs and td
    if ((calc_mode.ne.1).and.(calc_mode.ne.3)) then 
       message(1) = 'Error: No parallelization strategy defined for this runmode.'
       call write_fatal(1)
    endif

    call pop_sub()
  end subroutine multicomm_strategy


  ! check if a balanced distribution of nodes will be used
  subroutine multicomm_sanity_check(mc)
    type(multicomm_type) :: mc

    call push_sub('mpi.multicomm_sanity_check')

    ! the checks below are still very crude sanity checks that have to be improved
    select case(calc_mode)
    case(1) ! M_GS
       ! check if we would have more domain communicators than processors
       if (mc%n_domain_comm.gt.mc%n_node) then
          message(1) = 'Error: Have more domain communicators than processors:'
          write(message(2), '(a,i4)') '       Number of Processors :',mc%n_node
          write(message(3), '(a,i4)') '       Number of Domains    :',mc%n_domain_comm
          write(message(4), '(a,i4,a)') 'Restart octopus with at least', &
               mc%n_domain_comm,' processors.'
          call write_fatal(4)
       endif
       ! check for balanced node distribution
       if (mc%n_domain_node*mc%n_domain_comm.ne.mc%n_node) then
          message(1) = 'Error: Inbalanced distribution of nodes over domains.'
          write(message(2), '(a,i4,a)') 'Restart octopus with multiples of', &
               mc%n_domain_node*mc%n_domain_comm,' processors.'
          call write_fatal(2)
       endif

    case(3) ! M_TD
       ! FIXME: Add checks for TD
       message(1) = 'Error: TD multi-communicators not implemented yet.'
       call write_fatal(1)
    case default
       message(1) = 'Error: No multi-communicator defined for this runmode.'
       call write_fatal(1)
    end select

    call pop_sub()
  end subroutine multicomm_sanity_check


  ! create index and domain communicators
  subroutine multicomm_init(n_node, n_index, index_range, mc)
    integer, intent(in)  :: n_node, n_index
    integer, intent(in)  :: index_range(:)
    type(multicomm_type) :: mc

    integer :: i, j, k, l1, count, mpierr, group
    integer, allocatable :: domain_ranks(:), index_ranks(:), a(:)
    integer :: domain_group, index_group, MPI_COMM_WORLD_GROUP

    call push_sub('mpi.multicomm_init')

    mc%n_index = n_index  ! size(index_range)
    mc%n_node  = n_node

    allocate(mc%index_range(mc%n_index))
    allocate(mc%index_weight(mc%n_index))
    allocate(a(mc%n_index))

    ! query parallelization strategy
    call multicomm_strategy(index_range, mc)

    ! number of domain communicators
    mc%n_domain_comm = calc_domain_comm(mc%n_node, mc%index_range)

    ! calculate number of index communicators
    mc%n_index_comm = calc_index_comm(mc%n_index, mc%index_range)

    ! how many nodes per domain
    mc%n_domain_node = mc%n_node/mc%n_domain_comm

    ! sanity check of input
    call multicomm_sanity_check(mc)

    message(1) = stars
    write(message(2),'(a,i6)') 'Info: Number of domain communicators:', mc%n_domain_comm
    write(message(3),'(a,i6)') 'Info: Number of index  communicators:', mc%n_index_comm
    write(message(4),'(a,i6)') 'Info: Number of nodes per domain    :', mc%n_domain_node
    call write_info(4)

    ! Contract contract indices to number of available domain
    ! parallelizations (or nodes).
    ! If mc%n_index_comm*m%n_domain_node are available then
    ! each index combination will be on a seperate domain parallelization.
    allocate(mc%index_contraction(product(mc%index_range)))
    call contract_indices(mc%n_index, mc%index_weight, mc%index_range, &
                          mc%n_domain_comm, mc%index_contraction)

    ! allocate space to hold the handles of all required communicators 
    ! and groups
    allocate(mc%domain_comm (mc%n_domain_comm))
    allocate(mc%domain_root (mc%n_domain_comm))
    allocate(mc%index_comm  (mc%n_index_comm ))
    allocate(mc%domain_comm_of_node(0:mc%n_node-1))

    allocate(domain_ranks(mc%n_domain_node))

    ! initially set all communicators and groups to invalid handles
    domain_group    = MPI_GROUP_NULL
    index_group     = MPI_GROUP_NULL
    mc%domain_comm  = MPI_COMM_NULL
    mc%index_comm   = MPI_COMM_NULL

    ! get group of MPI_COMM_WORLD communicator
    call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_COMM_WORLD_GROUP, mpierr)

    ! create domain communicators
    message(1) = 'Info: Ranks of domain groups:'
    call write_info(1)
    count = 0
    l1 = 0
    do j = 1, mc%n_domain_comm
       mc%domain_root(j) = count
       do k = 1, mc%n_domain_node
          domain_ranks(k) = count
          count = count + 1
       enddo
       write(message(1),'(a,i4,a,100i10)') 'Info: Group',j,':',domain_ranks
       call write_info(1)
       call MPI_GROUP_INCL (MPI_COMM_WORLD_GROUP, mc%n_domain_node, &
            domain_ranks, domain_group, mpierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, domain_group,      &
            mc%domain_comm(j), mpierr)
       call MPI_GROUP_FREE (domain_group, mpierr)
       do k = 1, mc%n_domain_node
          ! Keep Domain communicator that node belongs to (rank -> communicator mapping)
          mc%domain_comm_of_node(l1) = mc%domain_comm(j)
          l1 = l1 + 1 
       enddo
    enddo
    message(1) = 'Info: Root nodes of domain groups:'
    write(message(2),'(a,100i10)') 'Info: ', mc%domain_root
    call write_info(2)

    if (mc%n_index_comm.gt.1) then
       message(1) = ''
       message(2) = 'Info: Ranks of index groups:'
       call write_info(2)   
    endif

    ! create index communicators
    do i = 1, mc%n_index
      ! Indices equal to one have no parallelism.
      if(mc%index_range(i).gt.1) then
        a = 1;
        ! This is a do-while loop which iterates over all
        ! indices except the i-th index. The counters for each
        ! index are in a.
        do ! Get number of communicator address by (a, i).
          group = index_comm_number(mc%n_index, mc%index_range, a, i)
          ! Ranks in MPI_COM_WORLD of future member nodes of 
          ! the communicator index_comm(group) are collected in this array.
          allocate(index_ranks(mc%index_range(i)))
          ! Get all members.
          do j = 1, mc%index_range(i)
            a(i)           = j
            index_ranks(j) = mc%domain_root(mc%index_contraction(                &
                             domain_comm_number(mc%n_index, mc%index_range, a)))
          end do
          ! Create communicator.
          call MPI_Group_incl(MPI_COMM_WORLD_GROUP, mc%index_range(j), &
               index_ranks, index_group, mpierr)
          call MPI_Comm_create(MPI_COMM_WORLD, index_group, &
               mc%index_comm(group), mpierr)
          call MPI_Group_free(index_group, mpierr)
          write(message(1),'(a,i4,a,100i10)') 'Info: Group', group, ':', index_ranks
          call write_info(1)
          deallocate(index_ranks)

          ! Go on iterating.
          if(.not.iterate_indices(mc%n_index, mc%index_range, a, i)) exit
        end do
      end if
    end do

    message(1) = stars
    call write_info(1)

    call MPI_Group_free(MPI_COMM_WORLD_GROUP, mpierr)

    deallocate(domain_ranks)

    call pop_sub()
  end subroutine multicomm_init

  
  ! ---------------------------------------------------------
  ! Contracts the index combinations to the available domain
  ! communicators.
  ! This is done by building a graph with index combinations
  ! as vertices. This graph is split into n_domain_comm
  ! parts using METIS. Because the communication intensity
  ! for different indices may be different a weight for each
  ! index may be specified. Heavier weights indicate more
  ! communication so that these indices end up on the same
  ! domain communicator.
  ! The result is a mapping m -> n
  ! where m is a number of an index communicator (with
  ! m = index_comm_number(..., a, i) for suitable a, i) and
  ! n is an index into domain_comm.
  !
  ! Example:
  ! n_index = 2
  ! index_range = (/2, 3/)
  ! index_weight = (/1000, 1/)
  ! n_domain_comm = 4
  !
  ! This looks like this
  !
  !        index 1
  !        +-->
  !        |  (1)   (4)
  ! index 2|  (2)   (5)
  !        V  (3)   (6)
  !
  ! Forming the following graph:
  !     (1)-1000-(4)
  !     1|        |1
  !     (2)-1000-(5)
  !     1|        |1
  !     (3)-1000-(6)
  ! which has to be split into four parts.
  !
  ! The result may be something like
  ! index_contraction(1) = 1
  ! index_contraction(2) = 2
  ! index_contraction(3) = 3
  ! index_contraction(4) = 1
  ! index_contraction(5) = 2
  ! index_contraction(6) = 4
  subroutine contract_indices(n_index, index_weight, index_range, &
                              n_domain_comm, index_contraction) 
    integer, intent(in)  :: n_index
    integer, intent(in)  :: index_weight(:)
    integer, intent(in)  :: index_range(:)
    integer, intent(in)  :: n_domain_comm
    integer, intent(out) :: index_contraction(:)

    integer              :: i, j, k, edgecut
    integer              :: nv         ! Number of vertices.
    integer              :: ne         ! Number of edges.
    integer              :: a(n_index) ! Indexing array.
    integer, allocatable :: xadj(:)    ! Adjacency index.
    integer, allocatable :: adjncy(:)  ! Adjacency lists.
    integer, allocatable :: weight(:)  ! Edge weights.
    integer              :: options(5) ! Options to METIS.

    call push_sub('mpi.contract_indices')

    options = (/1, 2, 1, 1, 0/) ! Use heavy edge matching in METIS.

    ! Space for adjacency and index lists.
    nv = product(index_range)
    allocate(xadj(nv+1))
    allocate(adjncy(2*nv*n_index))
    allocate(weight(2*nv*n_index))

    ! Build graph with index combinations as vertices and
    ! edges as index communicators.
    ne = 1
    a  = 1
    ! Iterate over all index combinations.
    ! The following code is very similair to mesh_partition in
    ! mesh_create.F90. But instead of restricting the graph
    ! to a vertex degree of <= 3 this graph has vertex degree n_index.
    j = 1
    do
      xadj(j) = ne
      ! Consider all neighbours.
      do i = 1, n_index
        ! This is a little bit dirty: Set a to next-neighbour.
        a(i) = a(i)+1
        k    = domain_comm_number(n_index, index_range, a)
        ! If k is existant, add edge.
        if(k.ne.0) then
          adjncy(ne) = k
          weight(ne) = index_weight(i)
          ne         = ne+1
        end if
        ! Set a to prev-neighbour.
        a(i) = a(i)-2
        k    = domain_comm_number(n_index, index_range, a)
        if(k.ne.0) then
          adjncy(ne) = k
          weight(ne) = index_weight(i)
          ne         = ne+1
        end if
        ! Restore a to not mix up iterating over a.
        a(i) = a(i)+1
      end do
      j = j+1
      if(.not.iterate_indices(n_index, index_range, a, 0)) exit
    end do
    xadj(nv+1) = ne
    ne         = ne-1

    call oct_metis_part_graph_recursive(nv, xadj, adjncy, 0,            &
      weight, 1, 1, n_domain_comm, options, edgecut, index_contraction)

    call pop_sub()

  end subroutine contract_indices


  ! ---------------------------------------------------------
  subroutine multicomm_end(mc)
    type(multicomm_type) :: mc

    integer :: j, mpierr

    call push_sub('mpi.multicomm_end')

    ! free domain communicators
    do j = 1, mc%n_domain_comm
       if (mc%domain_comm(j).ne.MPI_COMM_NULL) then
          call MPI_COMM_FREE (mc%domain_comm(j), mpierr)
       endif
    enddo

    ! free index communicators
    do j = 1, mc%n_index_comm
       if (mc%index_comm(j).ne.MPI_COMM_NULL) then
          call MPI_COMM_FREE (mc%index_comm(j), mpierr)
       endif
    enddo

    deallocate(mc%domain_comm_of_node)
    deallocate(mc%index_range)
    deallocate(mc%domain_comm)
    deallocate(mc%domain_root)
    deallocate(mc%index_comm)
    deallocate(mc%index_contraction)
    deallocate(mc%index_weight)

    call pop_sub()
  end subroutine multicomm_end


  ! ---------------------------------------------------------
  ! Helper functions for the indexing business.
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  ! Return number of index communicator for index i and all other
  ! indices set to values a(x), x|=i.
  ! Returns 0 for invalid (out of range) indices in a. If the
  ! result is later on used as an index the bounds
  ! checking of the compiler will alert as the index_comm array
  ! starts with 1.
  integer function index_comm_number(n_index, index_range, a, i) result(number)
    integer, intent(in) :: n_index
    integer, intent(in) :: index_range(:) ! index_range(n_index)
    integer, intent(in) :: a(:)           ! index_range(n_index)
                                          ! a(i) is ignored.
    integer, intent(in) :: i

    integer :: offset
    integer :: position
    integer :: offset_product
    integer :: position_product
    integer :: x, y

    call push_sub('mpi.index_comm_number')

    ! Check ranges first.
    if(.not.range_check(n_index, index_range, a, i)) then
      number = 0
      return
    end if

    ! See comment at top of this file for the formula computed
    ! and a description what is happening.

    ! Compute offset.
    offset = 0
    do  x = 1, i-1
      if(index_range(x).gt.1) then
        offset_product = 1
        do y = 1, n_index
          if(y.ne.x) then
            offset_product = offset_product*index_range(y)
          end if
        end do
        offset = offset+offset_product
      end if
    end do

    ! Compute position.
    position = 0
    do x = 1, n_index
      if(index_range(x).gt.1.and.x.ne.i) then
        position_product = 1
        do y = x+1, n_index
          if(y.ne.i) then
            position_product = position_product*index_range(y)
          end if
        end do
        position = position + (a(x)-1)*position_product
      end if
    end do
    position = position+1

    number = offset+position

    call pop_sub()

  end function index_comm_number


  ! ---------------------------------------------------------
  ! Return the number of domain communicator for indices set to
  ! a(x), x=1, ..., n_index.
  ! This maps (i1, i2, ..., i_n_index) -> n, i. e. computes how to
  ! store a multidimensional array in a linear field (i_n_index is
  ! stored consecutively).
  ! Returns 0 for invalid (out of range) indices in a.
  integer function domain_comm_number(n_index, index_range, a) result(number)
    integer, intent(in) :: n_index
    integer, intent(in) :: index_range(:)
    integer, intent(in) :: a(:)

    integer :: x
    
    call push_sub('mpi.domain_comm_number')

    ! Check ranges first.
    if(.not.range_check(n_index, index_range, a, 0)) then
      number = 0
      return
    end if

    number = 0

    do x = 1, n_index
      if(index_range(x).gt.1) then
        number = number + (a(x)-1)*product(index_range(x+1:n_index))
      end if
    end do
    number = number + 1

    call pop_sub()

  end function domain_comm_number 


  ! ---------------------------------------------------------
  ! Checks, if indices a(j), j=1, ..., n_index, j|=i
  ! are inside their ranges index_range(j).
  ! One index being out of range is sufficient for returning false.
  logical function range_check(n_index, index_range, a, i) result(inbound)
    integer, intent(in) :: n_index
    integer, intent(in) :: index_range(:)
    integer, intent(in) :: a(:)
    integer, intent(in) :: i

    integer :: j

    call push_sub('mpi.range_check')

    inbound = .true.
    j       = 1

    do while(inbound.and.j.le.n_index)
      if(j.ne.i) then
        inbound = a(j).le.index_range(j).and.a(j).gt.0
      end if
      j = j+1
    end do

    call pop_sub()

  end function range_check


  ! ---------------------------------------------------------
  ! Iterate over all indices except the i-th.
  ! This is like haveing n_index-1 do-loops with
  ! counters a(1), ..., a(i-1), a(i+1), ..., a(n_index)
  ! which iterate over 1, index_range(1, ..., i-1, i+1, ..., n_index)
  ! (a(1) is the counter of the outermost loop and a(n_index) is
  ! the counter of the innermost loop).
  ! The new values are found in a after the call.
  ! The return value is true if another iteration was possible
  ! (not all counters reached their end of range), false otherwise.
  ! If not 1 <= i <= n_index iterate over all indices.
  !
  ! A 'multiple' do-loop can be written like this:
  ! a = 1
  ! do
  !   ! Something useful.
  !   if(.not.iterate_indices(n_index, index_ranges, a, i)) exit
  ! end do
  !
  ! If Fortran had a non repelling loop this looked like
  ! do
  !   ! Someting useful.
  ! while(iterate_indices(n_index, index_ranges, a, i)
  logical function iterate_indices(n_index, index_range, a, i) result(iter)
    integer, intent(in)    :: n_index
    integer, intent(in)    :: index_range(:) ! index_range(n_index)
    integer, intent(inout) :: a(:)           ! index_range(n_index)
                                             ! a(i) is ignored.
    integer, intent(in)    :: i

    integer :: j

    call push_sub('mpi.iterate_indices')

    iter = .false.
    j    = n_index 

    ! Loop as long as we either could increase one counter
    ! (iter = true) or no counter is left (j <= 0).
    do while(.not.iter.and.j.gt.0)
      ! Leave out the i-th counter if necessary.
      if(j.ne.i) then
        ! Try to increase the j-tgh counter and exit loop.
        if(a(j).lt.index_range(j)) then
          a(j) = a(j)+1
          iter = .true.
        ! If that is not possible, try next one in next iteration.
        else
          a(j) = 1
          j    = j-1
        end if
      else
        j = j-1
      end if
    end do

    call pop_sub()

  end function iterate_indices


  ! ---------------------------------------------------------
  ! compute the number of required domain communicators
  integer function calc_domain_comm(n_node, index_range) result(n_domain_comm)
    integer, intent(in) :: n_node
    integer, intent(in) :: index_range(:)

    call push_sub('mpi.calc_domain_comm')

    ! In this case, ranges of 1 are no problem as they do not
    ! change the product of index_range(1)*...*index_range(n_index).
    n_domain_comm = product(index_range)

    if(n_domain_comm.gt.n_node) then
      n_domain_comm = n_node
    end if

    call pop_sub()

  end function calc_domain_comm


  ! ---------------------------------------------------------
  ! compute the number of required index communicators
  integer function calc_index_comm(n_index, index_range) result(n_index_comm)
    integer, intent(in) :: n_index 
    integer, intent(in) :: index_range(:)

    integer :: i

    call push_sub('mpi.calc_index_comm')

    ! The formula is mentioned in the comment on top of this file.
    ! Note again only index with ranges greater one are considered.
    n_index_comm = 1
    do i = 2, n_index
      if(index_range(i) > 1) then
        n_index_comm = n_index_comm*index_range(i) + product(index_range(1:i-1))
      end if
    end do

    call pop_sub()

  end function calc_index_comm


  ! =========================================================
  ! Routines to support MPI debugging.
  ! =========================================================

  subroutine MPI_Debug_Statistics()

    integer :: j
    integer :: usec_call(C_MPI_BARRIER:C_MPI_ALLREDUCE)

    if(.not.in_debug_mode) return

    message(1) = ''
    message(2) = hyphens
    message(3) = ''
    write(message(4), '(23x,a,4x,a,8x,a)') 'total time', 'calls', 'usec/call'
    do j = 1, C_MPI_ALLREDUCE
       if (sec_accum(j).eq.0.and.usec_accum(j).eq.0) then
          usec_call(j) = 0
       else
          usec_call(j) = (sec_accum(j)*1000000+usec_accum(j))/call_counter(j)
       endif

       write(message(j+4),'(a,i6,a,i6.6,6x,i4,6x,i10)')          &
            mpi_rlabel(j)//' : ',                                &
            sec_accum(j), '.', usec_accum(j), call_counter(j),   &
            usec_call(j)
    enddo
    message(C_MPI_ALLREDUCE+5) = ''    
    message(C_MPI_ALLREDUCE+6) = hyphens    
    call write_debug(C_MPI_ALLREDUCE+6)

  end subroutine MPI_Debug_Statistics


  ! ---------------------------------------------------------
  subroutine MPI_Debug_In(comm, index)
    integer, intent(in) :: comm, index

    if(.not.in_debug_mode) return

    call_counter(index) = call_counter(index) + 1
    call loct_gettimeofday(sec_in, usec_in)
    call epoch_time_diff(sec_in, usec_in)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i6.6,a,i4.4,a,i6.6)') '* I ',       &
         sec_in, '.', usec_in, ' '//trim(mpi_rlabel(index))//' - ', comm,':',   &
         call_counter(index), ' - ', sec_accum(index), '.', usec_accum(index)
    call write_debug(1)

  end subroutine MPI_Debug_IN


  ! ---------------------------------------------------------
  subroutine MPI_Debug_Out(comm, index)
    integer, intent(in) :: comm, index

    integer :: sec, usec, sec_diff, usec_diff

    if(.not.in_debug_mode) return

    call loct_gettimeofday(sec, usec)
    call epoch_time_diff(sec, usec)
    call mpi_time_accum(index, sec, usec, sec_diff, usec_diff)
    write(message(1),'(a,i6,a,i6.6,a,i3.3,a,i6.6,a,i4.4,a,i6.6,a,i4.4,a,i6.6)') &
         '* O ',                                                                &
         sec, '.', usec, ' '//trim(mpi_rlabel(index))//' - ', comm, ':',        &
         call_counter(index), ' - ', sec_accum(index), '.', usec_accum(index),  &
         ' - ', sec_diff, '.', usec_diff
    call write_debug(1)

  end subroutine MPI_Debug_Out


  ! ---------------------------------------------------------
  subroutine mpi_time_accum(index, sec, usec, sec_diff, usec_diff)
    integer, intent(in)  :: index, sec, usec  
    integer, intent(out) :: sec_diff, usec_diff

    integer :: sec_tmp, usec_tmp

    sec_tmp  = sec
    usec_tmp = usec

    if (usec_tmp-usec_in .lt. 0) then
       usec_tmp = usec_tmp + 1000000
       sec_tmp  = sec_tmp  - 1
    endif
    usec_tmp = usec_tmp - usec_in
    sec_tmp  = sec_tmp  - sec_in    

    usec_diff = usec_tmp
    sec_diff  = sec_tmp

    ! accumulate values
    if (usec_tmp+usec_accum(index) .gt. 1000000) then
       usec_tmp = usec_tmp - 1000000 
       sec_tmp  = sec_tmp  + 1
    endif
    sec_accum(index)  = sec_accum(index)  + sec_tmp
    usec_accum(index) = usec_accum(index) + usec_tmp

  end subroutine mpi_time_accum


#include "undef.F90"
#include "real.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mpi_inc.F90"

#include "undef.F90"
#include "integer.F90"
#include "mpi_inc.F90"

#endif
end module mpi_mod
