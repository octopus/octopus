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
! n_node_domain = n_node/n_domain with
! n_domain = prod(a=j, k)[index_range(a)] = 4
! nodes into each domain parallelization.
! To do collective operations (like sums) over j, k respectively
! n_index_comm = N(n_index) = 7 communicators are introduced
! (for the definition of N see below).
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
! vector and the adressing is done by a function get_comm(a, i).
! a(1:n_index) specifies the values for all indices except i
! (the value of a(i) is irrelevant) and i is the number of
! the requested index:
! get_comm(a, i) == offset + position
! WHERE
! offset   == sum(x=1, ..., i-1)[prod(y=1, ..., n_index, y|=x)[index_range(y)]]
! position == sum(x=1, ..., n_index, x!=i)[(a(x)-1)*
!             prod(y=x+1, ..., n_index, y|=i)[index_range(y)]] + 1
!
! For each index i the rest if the indices form a
! n_index-1 dimensional array. The offset of this array in
! index_comm is offset in the above function. It is the number of
! communicators all indices 1, ..., i-1 have.
! The position in the array is computed as for any other
! n-dimensional array (for generalities sake it has to be done by hand
! and not by the Fortran compiler).
! With this function, the j communicator for k=2 can be accessed with
! index_comm(get_comm((/0, 2/), 1))
! (with j being the first and k the second index).
!
! Some more stripped down cases:
! (*) Only index parallelization (e. g. k-points and states j):
!     There are no domain communicators. In the sketch above
!     (1), ..., (12) would directly denote nodes and not domain
!     parallelizations
! (*) Only domain parallelization: There would not be any index
!     communicators and just one domain parallelization:
!     n_domain = 1, domain_comm(1) = MPI_COMM_WORLD.
!
!
! ---------- 
! (*) N(1)   = 1
!     N(i+1) = N(i)*index_range(i+1) + prod(a=1, ..., i)[index_range(a)]
!
! (*) prod(x=1, ..., n)[f(x)] = f(1)* ... *f(x)
!

module mpi_mod
#if defined(HAVE_MPI)
  use varinfo
  use global
  use lib_oct
  use lib_oct_parser
  use syslabels
  use messages

  implicit none

  private

  public ::                              &
       MPI_Debug_Statistics,             &
       MPI_Debug_IN, MPI_Debug_OUT,      &
       multicomm_type,                   &
       multicomm_init, multicomm_end,    &
       multicomm_strategy,               &
       calc_index_comm, calc_domain_comm

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
    integer          :: n_node          ! Total number of nodes.
    integer          :: n_domain_nodes  ! Number of nodes per domain
    integer          :: n_index         ! Number of parallel indices.
    integer          :: n_domain_comm   ! Number of domain communicators.
    integer          :: n_index_comm    ! Number of index communicators.

    integer          :: par_strategy    ! What kind of parallelization strategy should we use?
    logical          :: use_domain_par  ! Should we use domain parallelization?

    integer, pointer :: index_range(:)  ! Range of index i is
                                        ! 1, ..., index_range(i).

    integer, pointer :: domain_comm(:)  ! Domain communicators ...
    integer, pointer :: domain_group(:) ! ... and corresponding groups.
    integer, pointer :: domain_root(:)  ! Ranks of roots from every domain
    integer, pointer :: index_comm(:)   ! Index communicators ...
    integer, pointer :: index_group(:)  ! ... and corresponding groups.

    ! mapping for comm_world ranks
    integer, pointer :: domain_comm_of_node(:)  ! Domain communicator that node belongs to
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


  ! compute the number of required domain communicators
  integer function calc_domain_comm(dim, index_range) result(n_index)
    integer, intent(in) :: dim
    integer, intent(in) :: index_range(:)

    integer :: j, range_product

    call push_sub('mpi.calc_index_comm')

    range_product = 1
    do j = 1, dim
       range_product = range_product*index_range(j)
    enddo

    n_index = range_product

    call pop_sub()
  end function calc_domain_comm


  ! compute the number of required index communicators
  integer function calc_index_comm(dim, index_range) result(n_index)
    integer, intent(in) :: dim
    integer, intent(in) :: index_range(:)

    integer :: j, k, range_product, rdim, count
    integer, allocatable :: n_index_table(:), proper_range(:)

    call push_sub('mpi.calc_index_comm')
    ! count how many real dimensions we have
    rdim = 0
    do k = 1, dim
       if(index_range(k).gt.1) rdim = rdim + 1
    enddo
    allocate(n_index_table(1:rdim), proper_range(1:rdim))

    ! keep only proper ranges
    count = 1
    do k = 1, dim
       if(index_range(k).gt.1) then
          proper_range(count) = index_range(k)
          count = count + 1
       endif
    enddo
    
    if (rdim.ne.0) n_index_table(1) = 1
    do k = 2, rdim
       range_product = 1
       do j = 1, k-1
          range_product = range_product*proper_range(j)
       enddo
       n_index_table(k) = n_index_table(k-1)*proper_range(k) + range_product
    enddo

    if (rdim.ne.0) then
       n_index = n_index_table(rdim)
    else
       n_index = 1
    endif

    deallocate(n_index_table, proper_range)

    call pop_sub()
  end function calc_index_comm


  ! decide which parallelization strategy we should use
  subroutine multicomm_strategy(index_dim, index_range, mc)
    integer, intent(in) :: index_dim, index_range(:)
    type(multicomm_type), intent(out) :: mc

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
       mc%use_domain_par = .false.
       mc%index_range(1) = n_kpoints
       mc%index_range(2) = n_states
       message(1) = 'Info: Octopus will run parallel in states and k-points.'
    case(P_STRATEGY_FULL)
       mc%use_domain_par = .true.
       mc%index_range(1) = n_kpoints
       mc%index_range(2) = n_states
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
       if (mc%n_domain_nodes*mc%n_domain_comm.ne.mc%n_node) then
          message(1) = 'Error: Inbalanced distribution of nodes over domains.'
          write(message(2), '(a,i4,a)') 'Restart octopus with multiples of', &
               mc%n_domain_nodes*mc%n_domain_comm,' processors.'
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

    integer :: j, k, l1, l2, count, mpierr, range_product, group
    integer, allocatable :: domain_ranks(:), index_ranks(:), stride(:)
    integer :: MPI_COMM_WORLD_GROUP

    call push_sub('mpi.multicomm_init')

    mc%n_index = n_index  ! size(index_range)
    mc%n_node  = n_node

    allocate(stride(1:mc%n_index))
    allocate(mc%index_range(1:mc%n_index))

    ! query parallelization strategy
    call multicomm_strategy(n_index, index_range, mc)

    ! number of domain communicators
    if(mc%use_domain_par) then 
       mc%n_domain_comm = calc_domain_comm(mc%n_index, mc%index_range)
    else
       mc%n_domain_comm = 1
    endif
    ! calculate number of index communicators
    mc%n_index_comm  = calc_index_comm (mc%n_index, mc%index_range)

    ! how many nodes per domain
    call loct_parse_int(check_inp('NumberOfProcessorsForDomain'), &
         mc%n_node/mc%n_domain_comm, mc%n_domain_nodes)

    ! sanity check of input
    call multicomm_sanity_check(mc)

    message(1) = stars
    write(message(2),'(a,i6)') 'Info: Number of domain communicators:', mc%n_domain_comm
    write(message(3),'(a,i6)') 'Info: Number of index  communicators:', mc%n_index_comm
    write(message(4),'(a,i6)') 'Info: Number of nodes per domain    :', mc%n_domain_nodes
    call write_info(4)

    ! allocate space to hold the handles of all required communicators 
    ! and groups
    allocate(mc%domain_comm (mc%n_domain_comm))
    allocate(mc%domain_group(mc%n_domain_comm))
    allocate(mc%domain_root (mc%n_domain_comm))
    allocate(mc%index_comm  (mc%n_index_comm ))
    allocate(mc%index_group (mc%n_index_comm ))
    allocate(mc%domain_comm_of_node(0:mc%n_node-1))

    allocate(domain_ranks(mc%n_domain_nodes))

    ! initially set all communicators and groups to invalid handles
    mc%domain_comm  = MPI_COMM_NULL
    mc%domain_group = MPI_GROUP_NULL
    mc%index_comm   = MPI_COMM_NULL
    mc%index_group  = MPI_GROUP_NULL

    ! get group of MPI_COMM_WORLD communicator
    call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_COMM_WORLD_GROUP, mpierr)

    ! create domain communicators
    message(1) = 'Info: Ranks of domain groups:'
    call write_info(1)
    count = 0
    l1 = 0
    do j = 1, mc%n_domain_comm
       mc%domain_root(j) = count
       do k = 1, mc%n_domain_nodes
          domain_ranks(k) = count
          count = count + 1
       enddo
       write(message(1),'(a,i4,a,100i10)') 'Info: Group',j,':',domain_ranks
       call write_info(1)
       call MPI_GROUP_INCL (MPI_COMM_WORLD_GROUP, mc%n_domain_nodes, &
            domain_ranks, mc%domain_group(j), mpierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, mc%domain_group(j),      &
            mc%domain_comm(j), mpierr)
       do k = 1, mc%n_domain_nodes
          ! Keep Domain communicator that node belongs to (rank -> communicator mapping)
          mc%domain_comm_of_node(l1) = mc%domain_comm(j)
          l1 = l1 + 1
       enddo
    enddo
    message(1) = 'Info: Root nodes of domain groups:'
    write(message(2),'(a,100i10)') 'Info: ',mc%domain_root
    call write_info(2)

    if (mc%n_index_comm.gt.1) then
       message(1) = ''
       message(2) = 'Info: Ranks of index groups:'
       call write_info(2)   
    endif

    ! setup strides
    stride = 1
    do j = 2, mc%n_index
       do k = 2, j
          stride(j) = stride(j)*mc%index_range(k-1)
       enddo
    enddo

    ! create index communicators
    group = 1
    do j = 1, mc%n_index
       range_product = 1
       do k = 1, mc%n_index
          if (j.ne.k) then
             range_product = range_product*mc%index_range(k)
          endif
       enddo
       count = 1
       do l1 = 1, range_product
          if (mc%index_range(j).eq.1) cycle
          allocate(index_ranks(mc%index_range(j)))
          do l2 = 1, mc%index_range(j)
             index_ranks(l2) = mc%domain_root(count)
             count = count + stride(j)
          enddo
          if (count.gt.mc%n_domain_comm) count = count - mc%n_domain_comm + 1
          write(message(1),'(a,i4,a,100i10)') 'Info: Group',group,':',index_ranks
          call write_info(1)
          call MPI_GROUP_INCL (MPI_COMM_WORLD_GROUP, mc%index_range(j), &
               index_ranks, mc%index_group(j), mpierr)
          call MPI_COMM_CREATE(MPI_COMM_WORLD, mc%index_group(j),       &
               mc%index_comm(j), mpierr)
          group = group + 1
          deallocate(index_ranks)
       enddo
    enddo

    message(1) = stars
    call write_info(1)

    call MPI_GROUP_FREE(MPI_COMM_WORLD_GROUP, mpierr)

    deallocate(domain_ranks, stride)

    call pop_sub()
  end subroutine multicomm_init


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
       if (mc%domain_comm(j).ne.MPI_GROUP_NULL) then
          call MPI_GROUP_FREE(mc%domain_group(j), mpierr)
       endif
    enddo

    ! free index communicators
    do j = 1, mc%n_index_comm
       if (mc%index_comm(j).ne.MPI_COMM_NULL) then
          call MPI_COMM_FREE (mc%index_comm(j), mpierr)
       endif
       if (mc%index_comm(j).ne.MPI_GROUP_NULL) then
          call MPI_GROUP_FREE(mc%index_group(j), mpierr)
       endif
    enddo

    deallocate(mc%domain_comm_of_node)
    deallocate(mc%index_range )
    deallocate(mc%domain_comm )
    deallocate(mc%domain_group)
    deallocate(mc%domain_root )
    deallocate(mc%index_comm  )
    deallocate(mc%index_group )

    call pop_sub()
  end subroutine multicomm_end


  ! ---------------------------------------------------------
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
