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

! This part regards the combined parallelization in indices
!
! Given a number n_index of indices and their ranges index_range(1:n_index),
! we divide the n_nodes in groups and create communicators for each group.
! Each group is associated with one index. The min_range indicates the minimum
! number of elements in each processor. For example, given min_range(1) = 25000,
! the algorithm will try to put at least 25000 points in each processor
!
! Example:
! Given 3 indices with ranges
!   index_range(1) = 100000, (number of points in the mesh)
!   index_range(2) = 15,     (number of states)
!   index_range(3) = 2,      (number of k-points)
! and 12 processors, we could get
!   mc%group_sizes = (2, 3, 2)
! which means that space is divided in 8 domains, the states are divided among 4 nodes,
! and the k-points in 2. One can view it as a tree (see degug/parallel_tree.dot)
!
!                       ------------------
!                       /                \         <- Division in k-points (1 k-point per node)
!                  ---------          ---------
!                /     |    \       /     |    \   <- Division in states   (5 states per node)
!              /  \  /  \  /  \   /  \  /  \  /  \ <- Division in domains  (50000 points per node)
!
! To perform collective operations (like a reduce), you can use the communicators
! provided in mc%group_comm(:). For example, to sum over states, the communicator
! to use is mc%group_comm(P_STRATEGY_STATES)
!
! You can use the routine multicomm_strategy_is_parallel to know is a certain
! strategy is parallel.

module multicomm_mod
  use varinfo
  use global
  use lib_oct
  use lib_oct_parser
  use syslabels
  use messages
  use mpi_mod
  use math
  use io

  implicit none

  private

  public ::                          &
    multicomm_type,                  &
    multicomm_tree_type,             &
    multicomm_tree_type_pointer,     &
    multicomm_init, multicomm_end,   &
    multicomm_strategy_is_parallel

  ! possible parallelization strategies
  integer, public, parameter ::      &
    P_STRATEGY_SERIAL  = 0,          & ! single domain, all states, kpoints on a single processor
    P_STRATEGY_DOMAINS = 1,          & ! parallelization domains
    P_STRATEGY_STATES  = 2,          & ! parallelization in kpoints
    P_STRATEGY_KPOINTS = 3,          & ! parallelization in states
    P_STRATEGY_OTHER   = 4             ! something else like e-h pairs

  integer,           parameter :: n_par_types = 4
  character(len=11), parameter :: par_types(0:n_par_types) = &
    (/                               &
    "serial     ",                   &
    "par_domains",                   &
    "par_states ",                   &
    "par_kpoints",                   &
    "par_other  "                    &
    /)

  type multicomm_tree_type_pointer
    type(multicomm_tree_type), pointer :: p
  end type multicomm_tree_type_pointer

  type multicomm_tree_type
    type(multicomm_tree_type_pointer), pointer :: up

    integer          :: n_down
    type(multicomm_tree_type_pointer), pointer :: down(:)

    integer          :: label        ! labels this node
    integer          :: level

    integer          :: n_group
    integer, pointer :: group(:)
  end type multicomm_tree_type

  ! Stores all communicators and groups
  type multicomm_type
    integer          :: n_node       ! Total number of nodes.
    integer          :: n_index      ! Number of parallel indices.

    integer          :: par_strategy ! What kind of parallelization strategy should we use?

    integer, pointer :: group_sizes(:) ! Number of processors in each group
    type(multicomm_tree_type), pointer :: group_tree

    integer, pointer :: who_am_i(:)    ! the path to get to my processor in the tree
    integer, pointer :: group_comm(:)  ! communicators I belong to
    integer, pointer :: group_root(:)  ! root node for each communicator
  end type multicomm_type


contains

  ! create index and domain communicators
  subroutine multicomm_init(mc, parallel_mask, n_node, n_index, index_range, min_range)
    type(multicomm_type), intent(out)  :: mc
    integer,              intent(in)   :: parallel_mask, n_node, n_index
    integer,              intent(inout):: index_range(:)
    integer,              intent(in)   :: min_range(:)

    integer :: i
    integer(POINTER_SIZE) :: blk

    call push_sub('multicomm.multicomm_init')

    ASSERT(n_index <= n_par_types)

    mc%n_node  = n_node
    mc%n_index = n_index  ! size(index_range)

    message(1) = stars
    call write_info(1)

    call strategy()
    if(mc%par_strategy.ne.P_STRATEGY_SERIAL) then
      allocate(mc%group_sizes(mc%n_index))
      mc%group_sizes(:) = 1

      !%Variable ParallelizationGroupRanks
      !%Type block
      !%Section Generalities::Parallel
      !%Description
      !% Specifies the size of the groups used for the parellization. For example
      !% (n_d, n_s, n_k) means we have n_p*n_s*n_k processors and that the k-points
      !% should be divided in n_k groups, the states in n_s groups, and each state
      !% in n_d domains.
      !%End
      if(loct_parse_block(check_inp('ParallelizationGroupRanks'), blk) == 0) then
        call read_block(blk)
        call loct_parse_block_end(blk)
      else
        call assign_nodes()
      end if

      ! reset par_strategy
      do i = 1, mc%n_index
        if(mc%group_sizes(i) == 1) mc%par_strategy = ibclr(mc%par_strategy, i-1)
      end do

      ! reset
      call sanity_check()
      call create_group_tree(mc)
      call group_comm_create()
    end if

    message(1) = stars
    call write_info(1)

    call pop_sub()


  contains

    ! ---------------------------------------------------------
    subroutine strategy()
      integer :: i, j, par_all

      !%Variable ParallelizationStrategy
      !%Type integer
      !%Default par_domains + par_states + par_kpoints
      !%Section Generalities::Parallel
      !%Description
      !% Specifies what kind of parallelization strategy octopus should use.
      !% The values can be combined, for example "par_domains + par_states"
      !% means a combined paralellization in domains and states
      !%Option serial 0
      !% Octopus will run in serial.
      !%Option par_domains 1
      !% Octopus will run parallel in domains.
      !%Option par_states  2
      !% Octopus will run parallel in states.
      !%Option par_kpoints 4
      !% Octopus will run parallel in k-points/spin.
      !%End

      if(mpiv%numprocs > 1) then
        par_all = 0
        do i = 1, n_par_types
          par_all = par_all + 2**(i-1)
        end do

        call loct_parse_int(check_inp('ParallelizationStrategy'),  &
          par_all, mc%par_strategy)

        if(mc%par_strategy<P_STRATEGY_SERIAL.or.mc%par_strategy>par_all) then
          call input_error('ParallelizationStrategy')
        end if

        mc%par_strategy = iand(mc%par_strategy, parallel_mask)

        if(mc%par_strategy == P_STRATEGY_SERIAL) then
          message(1) = "More than one node is available, but the this run mode can not run in parallel"
          message(2) = "Please select a ParallelizationStrategy compatible with"
          j = 2
          do i = 1, n_par_types
            if(iand(parallel_mask, 2**(i-1)).ne.0) then
              j = j + 1
              write(message(j), '(2a)') "  - ", par_types(i)
            end if
          end do
          call write_fatal(j)
        end if
      else
        mc%par_strategy = P_STRATEGY_SERIAL
      end if

      if(mc%par_strategy == P_STRATEGY_SERIAL) then
        message(1) = "Octopus will run in *serial*"
      else
        message(1) = "Octopus will run in *parallel*"
      end if
      call write_info(1)

    end subroutine strategy


    ! ---------------------------------------------------------
    subroutine read_block(blk)
      integer(POINTER_SIZE) :: blk

      integer :: i, n

      n = loct_parse_block_cols(blk, 0)

      mc%group_sizes = 1
      do i = 1, min(n, mc%n_index)
        if(multicomm_strategy_is_parallel(mc, i)) then
          call loct_parse_block_int(blk, 0, i-1, mc%group_sizes(i))
        end if
      end do
    end subroutine read_block


    ! ---------------------------------------------------------
    subroutine assign_nodes()
      integer :: i, n, k, n_divisors, divisors(50)
      integer, allocatable :: n_group_max(:)
      FLOAT   :: f

      allocate(n_group_max(mc%n_index))
      n = mc%n_node

      ! this is the maximum number of processors in each group
      n_group_max(1:mc%n_index) = max(index_range(1:mc%n_index)/min_range(1:mc%n_index), 1)
      do k = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, k)) n_group_max(k) = 1
      end do

      ! for each index
      do k = 1, mc%n_index

        ! distibute the nodes so that domains is the last
        f = real(n, PRECISION)
        do i = k+1, mc%n_index
          f = f/real(n_group_max(i), PRECISION)
        end do

        ! get divisors of n
        n_divisors = 50 ! maximum number of divisors
        call math_divisors(n, n_divisors, divisors)

        ! get the divisor of n >= f
        mc%group_sizes(k) = n
        do i = 1, n_divisors
          if(real(divisors(i), PRECISION) >= f) then
            mc%group_sizes(k) = divisors(i)
            exit
          end if
        end do

        n = n/mc%group_sizes(k)
      end do

      deallocate(n_group_max)

    end subroutine assign_nodes


    ! check if a balanced distribution of nodes will be used
    subroutine sanity_check()
      FLOAT :: frac
      integer :: i, k

      call push_sub('multicomm.sanity_check')

      ! print out some info
      i = 0
      do k = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, k)) cycle
        i = i + 1
        write(message(i),'(3a,i6,a,i8,a)') 'Info: Number of nodes in ', &
          par_types(k), ' group:', mc%group_sizes(k), ' (', index_range(k), ')'
      end do
      call write_info(i)

      ! do we have the correct number of processors
      if(product(mc%group_sizes(1:mc%n_index)).ne.mpiv%numprocs) then
        write(message(1),'(a,i4,a,i4,a)') "Inconsistent number of processors (", &
          product(mc%group_sizes(1:mc%n_index)), ".ne.", mpiv%numprocs, ")"
        message(2) = "You probably have an problem in the block 'ParallelizationGroupRanks'"
        call write_fatal(2)
      end if

      if(any(mc%group_sizes(1:mc%n_index) > index_range(1:mc%n_index))) then
        message(1) = "Could not distribute nodes in parallel job. Most likely you are trying to"
        message(2) = "use too many nodes for the job"
        call write_fatal(2)
      end if

      if(any(index_range(1:mc%n_index)/mc%group_sizes(1:mc%n_index) < min_range(1:mc%n_index))) then
        message(1) = "I have less elements in a parallel group than recommended."
        message(2) = "Maybe you should reduce the number of nodes"
        call write_warning(2)
      end if

      ! calculate fraction of idle time
      frac = M_ONE
      do i = 1, mc%n_index
        frac = frac * (M_ONE - real(mod(index_range(i), mc%group_sizes(i)), PRECISION)/ &
          real(index_range(i), PRECISION))
      end do

      write(message(1), '(a,f5.2,a)') "Info: Octopus will waste at least ", &
        (M_ONE - frac)*CNST(100.), "% of computer time"
      if(frac < CNST(0.8)) then
        message(2) = "I decided this is too much. Change the number of processors and try again."
        message(3) = "Usually number of processors multiple of small primes are best"
        call write_fatal(3)
      else
        call write_info(1)
      end if

      call pop_sub()
    end subroutine sanity_check


    ! ---------------------------------------------------------
    subroutine group_comm_create()
#if defined(HAVE_MPI)
      integer :: i, mpi_err
      integer, allocatable :: me(:)

      allocate(mc%group_comm(mc%n_index), mc%group_root(mc%n_index))
      mc%group_comm = -1

      allocate(me(mc%n_index))
      do i = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, i)) cycle

        me(:) = mc%who_am_i(:)
        me(i) = 1
        call multicomm_proc(mc, me, mc%group_root(i))

        call MPI_Comm_split(MPI_COMM_WORLD, mc%group_root(i), mpiv%node, &
          mc%group_comm(i), mpi_err)
      end do
      deallocate(me)
#endif
    end subroutine group_comm_create

  end subroutine multicomm_init


  ! ---------------------------------------------------------
  subroutine multicomm_proc(mc, index, proc)
    type(multicomm_type), intent(in)  :: mc
    integer,              intent(in)  :: index(:) ! (mc%n_index)
    integer,              intent(out) :: proc     ! the processor corresponding to that index

    type(multicomm_tree_type), pointer :: p
    integer :: i

    p => mc%group_tree
    do i = mc%n_index, 1, -1
      if(.not.multicomm_strategy_is_parallel(mc, i)) cycle
      ASSERT(index(i)>0.and.index(i)<=p%n_down)

      p => p%down(index(i))%p
    end do
    proc = p%group(1)

  end subroutine multicomm_proc


  ! ---------------------------------------------------------
  subroutine create_group_tree(mc)
    type(multicomm_type), intent(inout) :: mc

#if defined(HAVE_MPI)
    integer :: nodes_used, last_level, iunit, label
    integer, allocatable :: index_run(:)
#endif

    allocate(mc%who_am_i(mc%n_index))
    mc%who_am_i = 0

#if defined(HAVE_MPI)
    ! allocate master node
    allocate(mc%group_tree)
    nullify(mc%group_tree%up)

    nodes_used = 0

    ! obtain the lower_level
    last_level = 1
    do
      if((last_level>mc%n_index).or.multicomm_strategy_is_parallel(mc, last_level)) exit
      last_level = last_level + 1
    end do
    ASSERT(last_level <= mc%n_index)

    label = 0 ! will label each node of the tree

    allocate(index_run(mc%n_index))
    index_run = 0
    call nodes_create(mc%group_tree, mc%n_index+1)
    deallocate(index_run)
    ASSERT(.not.all(mc%who_am_i==0))

    call group_create(mc%group_tree)

    ! print nodes
    if(in_debug_mode.and.mpiv%node == 0) then
      iunit = io_open('debug/parallel_tree.dot', action='write')
      write(iunit, '(a)') "digraph G {"
      write(iunit, '(a)') 'node [shape=box,style=filled];'
      call print_dot(mc%group_tree)
      write(iunit, *) "}"
      call io_close(iunit)
    end if

  contains

    ! ---------------------------------------------------------
    recursive subroutine nodes_create(this, level)
      type(multicomm_tree_type), pointer :: this
      integer, intent(in)                :: level

      integer :: i, next_level
      type(multicomm_tree_type), pointer :: p

      ! get next parallel level
      next_level = level - 1
      do
        if((next_level==0).or.multicomm_strategy_is_parallel(mc, next_level)) exit
        next_level = next_level - 1
      end do

      ! set up some variables
      this%label  = label
      label = label + 1
      this%level  = level

      nullify(this%down, this%group)
      this%n_down = 0
      this%n_group = 0

      if(level == last_level) then ! last level, let us populate
        this%n_group = 1
        allocate(this%group(1))
        this%group(1) = nodes_used

        if(nodes_used == mpiv%node) mc%who_am_i = index_run

        nodes_used = nodes_used + 1

      else if(next_level .ne. 0) then

        this%n_down = mc%group_sizes(next_level)
        allocate(this%down(this%n_down))
        do i = 1, this%n_down
          index_run(next_level) = i

          allocate(this%down(i)%p)
          p => this%down(i)%p

          ! initialize parent of structure
          allocate(p%up)
          p%up%p => this ! point up to the parent
          call nodes_create(p, next_level)
        end do
      end if

    end subroutine nodes_create


    ! ---------------------------------------------------------
    recursive subroutine group_create(this)
      type(multicomm_tree_type), pointer :: this

      integer :: i

      if(this%n_group.ne.0) return ! last level

      ! first go through all children
      do i = 1, this%n_down
        call group_create(this%down(i)%p)
      end do

      ! now create group
      this%n_group = this%n_down
      allocate(this%group(this%n_group))
      do i = 1, this%n_down
        this%group(i) = this%down(i)%p%group(1) ! group with root nodes
      end do
    end subroutine group_create


    ! ---------------------------------------------------------
    recursive subroutine print_dot(this)
      type(multicomm_tree_type), pointer :: this

      integer :: i

      write(iunit,'(i4,a)',   advance='no') this%label, ' [label="'
      write(iunit,'(a,i4,a)', advance='no') 'level = ', this%level, '\n'
      write(iunit,'(a)', advance='no') 'group = '
      do i = 1, this%n_group
        if(i.ne.1) write(iunit, '(a)', advance='no') ','
        write(iunit, '(i4)', advance='no') this%group(i)
      end do

      write(iunit,'(a)') '"]'
      do i = 1, this%n_down
        write(iunit,*) this%label, "->", this%down(i)%p%label, ";"

        call print_dot(this%down(i)%p)
      end do
    end subroutine print_dot
#endif

  end subroutine create_group_tree


  ! ---------------------------------------------------------
  subroutine multicomm_end(mc)
    type(multicomm_type) :: mc

#if defined(HAVE_MPI)
    integer :: i, mpi_err

    call push_sub('multicomm.multicomm_end')

    if(mc%par_strategy.ne.P_STRATEGY_SERIAL) then
      ! delete communicators
      do i = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, i)) cycle

        call MPI_Comm_free(mc%group_comm(i), mpi_err)
      end do

      ! now delete the tree
      call delete_tree(mc%group_tree)
      deallocate(mc%group_tree)
      nullify(mc%group_tree)

      ! deallocate the rest of the arrays
      deallocate(mc%group_sizes, mc%who_am_i, mc%group_comm, mc%group_root)
      nullify(mc%group_sizes, mc%who_am_i, mc%group_comm, mc%group_root)
    end if

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    recursive subroutine delete_tree(this)
      type(multicomm_tree_type), pointer :: this

      integer :: i

      deallocate(this%group)

      ! delete children
      do i = 1, this%n_down
        call delete_tree(this%down(i)%p)
        deallocate(this%down(i)%p)
      end do

      if(this%n_down.ne.0) then
        deallocate(this%down)
      end if

    end subroutine delete_tree
#else
    deallocate(mc%group_tree, mc%group_sizes, mc%who_am_i, mc%group_comm, mc%group_root)
    nullify(mc%group_tree, mc%group_sizes, mc%who_am_i, mc%group_comm, mc%group_root)
#endif
  end subroutine multicomm_end

  ! ---------------------------------------------------------
  logical function multicomm_strategy_is_parallel(mc, level) result(r)
    type(multicomm_type), intent(in) :: mc
    integer,              intent(in) :: level

    r = iand(mc%par_strategy, 2**(level-1)).ne.0
  end function multicomm_strategy_is_parallel

end module multicomm_mod
