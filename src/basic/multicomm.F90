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
! which means
! * that space is divided in 2 domains per state,
! * the states are divided in 3 groups, i. e. 5 states per processor, and
! * the whole setting is duplicated because of the 2 kpoints.
!
! To perform collective operations (like a reduce), you can use the communicators
! provided in mc%group_comm(:). For example, to sum over states, the communicator
! to use is mc%group_comm(P_STRATEGY_STATES)
!
! You can use the routine multicomm_strategy_is_parallel to know if a certain
! index is parallelized.

  module multicomm_m
    use datasets_m
    use global_m
    use io_m
    use loct_m
    use loct_parser_m
    use messages_m
    use mpi_m
    use utils_m
    use varinfo_m
    use profiling_m
#if defined(USE_OMP)
    use omp_lib
#endif

  implicit none

  private

  public ::                          &
    multicomm_divide_range,          &
#ifdef USE_OMP
    multicomm_divide_range_omp,      &
#endif
#if defined(HAVE_MPI)
    multicomm_create_all_pairs,      &
#endif
    multicomm_t,                     &
    multicomm_all_pairs_t,           &
    multicomm_init, multicomm_end,   &
    multicomm_all_pairs_copy,        &
    multicomm_strategy_is_parallel,  &
    topology_t

  ! possible parallelization strategies
  integer, public, parameter ::      &
    P_STRATEGY_SERIAL  = 0,          & ! single domain, all states, k-points on a single processor
    P_STRATEGY_DOMAINS = 1,          & ! parallelization domains
    P_STRATEGY_STATES  = 2,          & ! parallelization in states
    P_STRATEGY_KPOINTS = 3,          & ! parallelization in k-points
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

  integer, public, parameter :: MAX_OMP_THREADS = 16
  integer, parameter :: MAX_INDEX = 5

  type topology_t
    integer          :: ng
    integer          :: maxgsize
    integer, pointer :: distance(:, :)
    integer, pointer :: groups(:,:)
    integer, pointer :: gsize(:)
  end type topology_t

  ! Stores all communicators and groups
  type multicomm_t
    integer          :: n_node         ! Total number of nodes.
    integer          :: n_index        ! Number of parallel indices.

    integer          :: par_strategy   ! What kind of parallelization strategy should we use?

    integer, pointer :: group_sizes(:) ! Number of processors in each group.
    integer, pointer :: who_am_i(:)    ! Rank in the "line"-communicators.
    integer, pointer :: group_comm(:)  ! "Line"-communicators I belong to.
    integer          :: dom_st_comm    ! States-domain plane communicator.

    integer          :: nthreads
    logical          :: use_topology
    type(topology_t) :: topo
  end type multicomm_t

  ! An all-pairs communication schedule for a given group.
  type multicomm_all_pairs_t
    type(mpi_grp_t)  :: grp            ! Schedule for this group.
    integer          :: rounds         ! This many comm. rounds.
    integer, pointer :: schedule(:, :) ! This is the schedule.
  end type multicomm_all_pairs_t

contains

  ! ---------------------------------------------------------
  subroutine multicomm_all_pairs_copy(apout, apin)
    type(multicomm_all_pairs_t), intent(inout) :: apout
    type(multicomm_all_pairs_t), intent(in)    :: apin
    integer :: i

    call mpi_grp_copy(apout%grp, apin%grp)
    apout%rounds = apin%rounds
    if(associated(apin%schedule)) then
      i = size(apin%schedule, 1)*size(apin%schedule, 2)
      SAFE_ALLOCATE(apout%schedule(1:size(apin%schedule, 1), 1:size(apin%schedule, 2)))
      apout%schedule = apin%schedule
    end if    

  end subroutine multicomm_all_pairs_copy

  ! ---------------------------------------------------------
  ! create index and domain communicators
  subroutine multicomm_init(mc, parallel_mask, default_mask, n_node, n_index, index_range, min_range)
    type(multicomm_t), intent(out)  :: mc
    integer,           intent(in)   :: parallel_mask, default_mask, n_node, n_index
    integer,           intent(inout):: index_range(:)
    integer,           intent(in)   :: min_range(:)

    integer   :: i
    type(block_t) :: blk

    call push_sub('multicomm.multicomm_init')

    ASSERT(n_index <= n_par_types)

    mc%n_node  = n_node
    mc%n_index = n_index  ! size(index_range)

    call messages_print_stress(stdout, "Parallelization")

    call strategy()

    if (mc%nthreads > 1) then
      write(message(1),'(a, i3)') 'Info: Number of threads ', mc%nthreads
      call write_info(1)
    end if

    if(mc%par_strategy.ne.P_STRATEGY_SERIAL) then
      SAFE_ALLOCATE(mc%group_sizes(1:mc%n_index))
      mc%group_sizes(:) = 1

      !%Variable ParallelizationUseTopology
      !%Type block
      !%Section Execution::Parallelization
      !%Description
      !%
      !% When set to yes, Octopus will try to determine the topology
      !% of the network and, based on this, assign the nodes to the
      !% different possible parallelization strategies and distribute
      !% the processes to reduce the communication overhead. Currently
      !% Octopus is only capable of detecting if processes are on the
      !% same machine or connected by a network. By default it is not
      !% enabled.
      !%
      !% This feature is experimental, so use it with care.
      !%
      !% Warning: currently is not possible to use this variable for
      !% states parallelization for ground-state calculations or
      !% Car-Parrinello molecular dynamics.
      !%
      !%End
  
      call loct_parse_logical(datasets_check('ParallelizationUseTopology'), .false., mc%use_topology)

      mc%use_topology = mc%use_topology .and. multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)
      mc%use_topology = mc%use_topology .and. multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)
      
      if(mc%use_topology) then
        call topology_init(mc%topo)
        mc%use_topology = mc%use_topology .and. topology_groups_are_equal(mc%topo)
      end if

      !%Variable ParallelizationGroupRanks
      !%Type block
      !%Section Execution::Parallelization
      !%Description
      !% Specifies the size of the groups used for the
      !% parallelization. For example (n_d, n_s, n_k) means we have
      !% <i>n_p*n_s*n_k</i> processors and that the <i>k</i>-points should be
      !% divided in <i>n_k</i> groups, the states in <i>n_s</i> groups, and each
      !% state in <i>n_d</i> domains. You can pass the value <tt>fill</tt> to one
      !% field: it will be replaced by the value required to complete
      !% the number of processors in the run.
      !%Option fill -1
      !% Replaced by the value required to complete the number of processors.
      !%End
      if(loct_parse_block(datasets_check('ParallelizationGroupRanks'), blk) == 0) then
        call read_block(blk)
        call loct_parse_block_end(blk)
      else
        call assign_nodes()
      end if

      ! clear parallel strategies that were available but will not be used
      do i = 1, mc%n_index
        if(mc%group_sizes(i) == 1) mc%par_strategy = ibclr(mc%par_strategy, i-1)
      end do

      ! reset
      call sanity_check()
      call cart_topology_create()
      call group_comm_create()
    end if

    call messages_print_stress(stdout)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine strategy()
      integer :: i, j,  par_mask

      !%Variable ParallelizationStrategy
      !%Type flag
      !%Default par_domains + par_states + par_kpoints
      !%Section Execution::Parallelization
      !%Description
      !% Specifies what kind of parallelization strategy Octopus should use.
      !% The values can be combined, for example <tt>par_domains + par_states</tt>
      !% means a combined parallelization in domains and states.
      !%Option serial 0
      !% Octopus will run in serial.
      !%Option par_domains 1
      !% Octopus will run parallel in domains.
      !%Option par_states  2
      !% Octopus will run parallel in states.
      !%Option par_kpoints 4
      !% Octopus will run parallel in <i>k</i>-points/spin.
      !%Option par_other   8
      !% Run-mode dependent. For example, in <tt>casida</tt>, it means parallelization in e-h pairs.
      !%End

      if(mpi_world%size > 1) then

        par_mask = parallel_mask
#if !defined(HAVE_METIS)
        ! parallelization in domains is not allowed
        par_mask = ibclr(par_mask, P_STRATEGY_DOMAINS - 1)
#endif

        call loct_parse_int(datasets_check('ParallelizationStrategy'), default_mask, mc%par_strategy)

        if(.not.varinfo_valid_option('ParallelizationStrategy', mc%par_strategy, is_flag=.true.)) then
          call input_error('ParallelizationStrategy')
        end if

        mc%par_strategy = iand(mc%par_strategy, par_mask)

        if(mc%par_strategy == P_STRATEGY_SERIAL) then
          message(1) = "More than one node is available, but this run mode can not run in parallel."
          message(2) = "Please select a ParallelizationStrategy compatible with"
          j = 2
          do i = 1, n_par_types
            if(iand(par_mask, 2**(i-1)).ne.0) then
              j = j + 1
              write(message(j), '(2a)') "  - ", par_types(i)
            end if
          end do
          call write_fatal(j)
        end if
      else
        mc%par_strategy = P_STRATEGY_SERIAL
      end if

      mc%nthreads = 1
#if defined(USE_OMP)
      !$omp parallel
      !$omp master
      mc%nthreads = omp_get_num_threads()
      !$omp end master
      !$omp end parallel
      if(mc%nthreads > MAX_OMP_THREADS) then
        message(1) = "Number of threads requested is larger than MAX_OMP_THREADS"
        call write_fatal(1)
      end if
#endif

      if(mc%par_strategy == P_STRATEGY_SERIAL .and. mc%nthreads == 1) then
        message(1) = "Octopus will run in *serial*"
      else
        message(1) = "Octopus will run in *parallel*"
      end if
      call write_info(1)

    end subroutine strategy


    ! ---------------------------------------------------------
    subroutine read_block(blk)
      type(block_t), intent(inout) :: blk

      integer :: i, n
      logical :: fill_used

      n = loct_parse_block_cols(blk, 0)

      mc%group_sizes = 1
      do i = 1, min(n, mc%n_index)
        if(multicomm_strategy_is_parallel(mc, i)) then
          call loct_parse_block_int(blk, 0, i-1, mc%group_sizes(i))
        end if
      end do

      fill_used = .false.
      do i = 1, mc%n_index
        if(mc%group_sizes(i) == -1) then
          if(fill_used) then
            message = "Error: The 'fill' value can be used only once in ParallelizationGroupRanks."
            call write_fatal(1)
          end if
          mc%group_sizes(i) = -mpi_world%size/product(mc%group_sizes)
          fill_used = .true.
        end if
      end do

    end subroutine read_block

    ! ---------------------------------------------------------
    subroutine assign_nodes()
      integer :: i, n, k, n_divisors, divisors(50)
      integer, allocatable :: n_group_max(:)
      FLOAT   :: f

      if(mc%use_topology) then
        mc%group_sizes = 1
        mc%group_sizes(P_STRATEGY_DOMAINS) = mc%topo%maxgsize
        mc%group_sizes(P_STRATEGY_STATES)  = mc%topo%ng
        return
      end if

      SAFE_ALLOCATE(n_group_max(1:mc%n_index))

      n = mc%n_node

      ! this is the maximum number of processors in each group
      n_group_max(1:mc%n_index) = max(index_range(1:mc%n_index), 1)
      do k = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, k)) n_group_max(k) = 1
      end do

      ! for each index
      do k = 1, mc%n_index
        if(n_group_max(k) == 1) then ! not parallel in this group
          mc%group_sizes(k) = 1
          cycle
        end if

        ! distibute the nodes so that domains is the last
        f = real(n, REAL_PRECISION)
        do i = k+1, mc%n_index
          f = f/real(n_group_max(i), REAL_PRECISION)
        end do

        ! get divisors of n
        n_divisors = 50 ! maximum number of divisors
        call get_divisors(n, n_divisors, divisors)

        ! get the divisor of n >= f
        mc%group_sizes(k) = n
        do i = 1, n_divisors
          if(real(divisors(i), REAL_PRECISION) >= f) then
            mc%group_sizes(k) = divisors(i)
            exit
          end if
        end do

        n = n/mc%group_sizes(k)
      end do

      SAFE_DEALLOCATE_A(n_group_max)

    end subroutine assign_nodes


    ! ---------------------------------------------------------
    ! check if a balanced distribution of nodes will be used
    subroutine sanity_check()
      FLOAT :: frac
      integer :: i, k, n_max

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
      if(product(mc%group_sizes(1:mc%n_index)).ne.mpi_world%size) then
        write(message(1),'(a,i4,a,i4,a)') "Inconsistent number of processors (", &
          product(mc%group_sizes(1:mc%n_index)), ".ne.", mpi_world%size, ")"
        message(2) = "You probably have a problem in the block 'ParallelizationGroupRanks'"
        call write_fatal(2)
      end if

      if(any(mc%group_sizes(1:mc%n_index) > index_range(1:mc%n_index))) then
        message(1) = "Could not distribute nodes in parallel job. Most likely you are trying to"
        message(2) = "use too many nodes for the job"
        call write_fatal(2)
      end if

      if(any(index_range(1:mc%n_index)/mc%group_sizes(1:mc%n_index) < min_range(1:mc%n_index))) then
        message(1) = "I have fewer elements in a parallel group than recommended."
        message(2) = "Maybe you should reduce the number of nodes"
        call write_warning(2)
      end if

      ! calculate fraction of idle time
      frac = M_ONE
      do i = 1, mc%n_index
        n_max = ceiling(real(index_range(i), REAL_PRECISION)/real(mc%group_sizes(i), REAL_PRECISION))
        k = n_max * mc%group_sizes(i)
        frac = frac * (M_ONE - real(k - index_range(i), REAL_PRECISION)/real(k, REAL_PRECISION))
      end do

      write(message(1), '(a,f5.2,a)') "Info: Octopus will waste at least ", &
        (M_ONE - frac)*CNST(100.), "% of computer time"
      if(frac < CNST(0.8)) then
        message(2) = "I decided this is too much. Change the number of processors and try again."
        message(3) = "Usually a number of processors which is a multiple of small primes is best."
        call write_fatal(3)
      else
        call write_info(1)
      end if

      call pop_sub()
    end subroutine sanity_check


    ! ---------------------------------------------------------
    subroutine cart_topology_create()
#if defined(HAVE_MPI)
      integer :: new_comm
      logical :: reorder, periodic_mask(MAX_INDEX)
#endif

      call push_sub('multicomm.cart_topology_create')

#if defined(HAVE_MPI)
      if(.not. mc%use_topology) then

        periodic_mask = .false.
        reorder = .true.

        ! The domain and states dimensions have to be periodic (2D torus)
        ! in order to circulate matrix blocks.
        if(multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)) then
          periodic_mask(P_STRATEGY_DOMAINS) = .true.
        end if
        if(multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)) then
          periodic_mask(P_STRATEGY_STATES) = .true.
        end if
        ! We allow reordering of ranks, as we intent to replace the
        ! world afterwards.
        ! FIXME: make sure this works! World root may be someone else
        ! afterwards!
        call MPI_Cart_create(mpi_world%comm, mc%n_index, mc%group_sizes, periodic_mask, reorder, new_comm, mpi_err)

        ! Re-initialize the world.
        call mpi_grp_init(mpi_world, new_comm)
      end if
#endif

      call pop_sub()
    end subroutine cart_topology_create


    ! ---------------------------------------------------------
    subroutine group_comm_create()
#if defined(HAVE_MPI)
      logical :: dim_mask(MAX_INDEX)
      integer :: world_group, sub_group, dummy_comm
      integer :: ii
      integer :: i_strategy
#endif

      call push_sub('multicomm.group_comm_create')

      SAFE_ALLOCATE(mc%group_comm(1:mc%n_index))
      SAFE_ALLOCATE(mc%who_am_i(1:mc%n_index))

#if defined(HAVE_MPI)
        if(.not. mc%use_topology) then

          ! The "lines" of the Cartesian grid.
          do i_strategy = 1, mc%n_index
            if(multicomm_strategy_is_parallel(mc, i_strategy)) then
              dim_mask             = .false.
              dim_mask(i_strategy) = .true.
              call MPI_Cart_sub(mpi_world%comm, dim_mask, mc%group_comm(i_strategy), mpi_err)
              call MPI_Comm_rank(mc%group_comm(i_strategy), mc%who_am_i(i_strategy), mpi_err)
            else
              mc%group_comm(i_strategy) = MPI_COMM_NULL
              mc%who_am_i(i_strategy)   = 0
            end if
          end do

          ! The domain-state "planes" of the grid (the ones with periodic dimensions).
          dim_mask                     = .false.
          dim_mask(P_STRATEGY_DOMAINS) = .true.
          dim_mask(P_STRATEGY_STATES)  = .true.
          call MPI_Cart_sub(mpi_world%comm, dim_mask, mc%dom_st_comm, mpi_err)
          
        else

          call MPI_Comm_group(mpi_world%comm, world_group, mpi_err)

          ! domains
          do ii = 1, mc%topo%ng
            call MPI_Group_incl(world_group, mc%topo%gsize(1), mc%topo%groups(1:mc%topo%gsize(1), ii), sub_group, mpi_err)
            call MPI_Comm_create(mpi_world%comm, sub_group, dummy_comm, mpi_err)

            if(dummy_comm /= MPI_COMM_NULL) then
              mc%group_comm(P_STRATEGY_DOMAINS) = dummy_comm
              call MPI_Comm_rank(mc%group_comm(P_STRATEGY_DOMAINS), mc%who_am_i(P_STRATEGY_DOMAINS), mpi_err)
            end if
          end do

          ! states
          do ii = 1, mc%topo%gsize(1)
            call MPI_Group_incl(world_group, mc%topo%ng, mc%topo%groups(ii, 1:mc%topo%ng), sub_group, mpi_err)
            call MPI_Comm_create(mpi_world%comm, sub_group, dummy_comm, mpi_err)

            if(dummy_comm /= MPI_COMM_NULL) then
              mc%group_comm(P_STRATEGY_STATES) = dummy_comm
              call MPI_Comm_rank(mc%group_comm(P_STRATEGY_STATES), mc%who_am_i(P_STRATEGY_STATES), mpi_err)
            end if
          end do

          mc%group_comm(P_STRATEGY_KPOINTS) = MPI_COMM_NULL
          mc%who_am_i(P_STRATEGY_KPOINTS)   = 0
          mc%group_comm(P_STRATEGY_OTHER)   = MPI_COMM_NULL
          mc%who_am_i(P_STRATEGY_OTHER)     = 0

          call MPI_Comm_dup(mpi_world%comm, mc%dom_st_comm, mpi_err)

        end if
#else
        mc%group_comm = -1
        mc%who_am_i   = 0
#endif

      call pop_sub()
    end subroutine group_comm_create
  end subroutine multicomm_init
  

  ! ---------------------------------------------------------
    subroutine multicomm_end(mc)
      type(multicomm_t), intent(inout) :: mc

#if defined(HAVE_MPI)
      integer :: i
#endif

      call push_sub('multicomm.multicomm_end')

    if(mc%par_strategy.ne.P_STRATEGY_SERIAL) then
#if defined(HAVE_MPI)
      ! Delete communicators.
      do i = 1, mc%n_index
        if(.not.multicomm_strategy_is_parallel(mc, i)) cycle
        call MPI_Comm_free(mc%group_comm(i), mpi_err)
      end do
      call MPI_Comm_free(mc%dom_st_comm, mpi_err)
#endif
      ! Deallocate the rest of the arrays.
      SAFE_DEALLOCATE_P(mc%group_sizes)
      SAFE_DEALLOCATE_P(mc%group_comm)
      SAFE_DEALLOCATE_P(mc%who_am_i)
    end if

    call pop_sub()
  end subroutine multicomm_end


  ! ---------------------------------------------------------
  logical function multicomm_strategy_is_parallel(mc, level) result(r)
    type(multicomm_t), intent(in) :: mc
    integer,           intent(in) :: level

    r = iand(mc%par_strategy, 2**(level-1)).ne.0
  end function multicomm_strategy_is_parallel


  ! ---------------------------------------------------------
  ! This routine uses the one-factorization (or near-one-factorization
  ! of a complete graph to construct an all-pair communication
  ! schedule (cf. Wang, X., Blum, E. K., Parker, D. S., and Massey,
  ! D. 1997. The dance party problem and its application to collective
  ! communication in computer networks. Parallel Comput. 23, 8
  ! (Aug. 1997), 1141-1156.
#if defined(HAVE_MPI)

  subroutine multicomm_create_all_pairs(mpi_grp, ap)
    type(mpi_grp_t),             intent(in)  :: mpi_grp
    type(multicomm_all_pairs_t), intent(out) :: ap

    integer :: grp_size, rounds, ir, in

    call push_sub('multicomm.create_all_pairs')

    ap%grp = mpi_grp
    grp_size   = mpi_grp%size

    ! Number of rounds.
    if(mod(grp_size, 2).eq.0) then
      rounds = grp_size-1
    else
      rounds = grp_size
    end if
    ap%rounds = rounds

    ! Calculate schedule.
    SAFE_ALLOCATE(ap%schedule(0:grp_size-1, 1:rounds))
    do ir = 1, rounds
      do in = 0, grp_size-1
        ap%schedule(in, ir) = get_partner(in+1, ir)-1
      end do
    end do

    call pop_sub()

  contains

    ! Those are from the paper cited above.
    integer function get_partner(in, ir)
      integer, intent(in) :: in, ir

      if(mod(grp_size, 2).eq.0) then
        get_partner = get_partner_even(grp_size, in-1, ir-1) + 1
      else
        get_partner = get_partner_odd(grp_size, in-1, ir-1) + 1
      end if
    end function get_partner

    integer function get_partner_even(grp_size, i, r) result(p)
      integer, intent(in) :: grp_size, i, r

      integer :: m

      m = grp_size/2

      if(i.eq.0) then
        p = r+1
      elseif(i.eq.r+1) then
        p = 0
      else
        ! I never know when to use which remainder function, but here
        ! it has to be the modulo one. Do not change that!
        p = modulo(2*r-i+1, 2*m-1)+1
      end if
    end function get_partner_even

    integer function get_partner_odd(grp_size, i, r) result(p)
      integer, intent(in) :: grp_size, i, r

      integer :: m

      m = (grp_size+1)/2

      p = get_partner_even(grp_size+1, i, r)

      if(p.eq.2*m-1) then
        p = i
      end if
    end function get_partner_odd
  end subroutine multicomm_create_all_pairs
#endif

    ! ---------------------------------------------------------
    
    ! this routine tries to guess the distribution of the processors
    ! we got, currently only checks processes that are running in the
    ! same node
    
    subroutine topology_init(this)
      type(topology_t), intent(out) :: this

#ifdef HAVE_MPI
      character(len=25) :: my_name, its_name
      integer :: ir,  wsize, ig

      wsize = mpi_world%size

      !get the system name
      call loct_sysname(my_name)

      SAFE_ALLOCATE(this%distance(1:wsize, 1:wsize))

      do ir = 1, wsize

        if(ir - 1 == mpi_world%rank) then

          call MPI_Bcast(my_name, 256, MPI_CHARACTER, ir - 1, mpi_world%comm, mpi_err)

          this%distance(ir, mpi_world%rank + 1) = 0

        else

          call MPI_Bcast(its_name, 256, MPI_CHARACTER, ir - 1, mpi_world%comm, mpi_err)

          if(my_name == its_name) then
            this%distance(ir, mpi_world%rank + 1) = 1
          else
            this%distance(ir, mpi_world%rank + 1) = 2
          end if

        end if

      end do

      do ir = 1, wsize
        call MPI_Bcast(this%distance(1, ir), wsize, MPI_INTEGER, ir - 1, mpi_world%comm, mpi_err)
      end do
      
      !classify processors in groups

      SAFE_ALLOCATE(this%groups(1:wsize, 1:wsize))
      SAFE_ALLOCATE(this%gsize(1:wsize))

      ! put the first node in the first group
      this%ng = 1
      this%gsize = 0
      this%groups = 0
      this%groups(1, 1) = 1
      this%gsize(1) = 1

      ! check the other processors
      do ir = 2, wsize
        do ig = 1, wsize
          ! if the group has elements
          if(this%gsize(ig) > 0) then
            ! check the distance
            if(this%distance(ir, this%groups(1, ig)) == 1) then
              ! if it is close enough
              ! add it to the group
              this%gsize(ig) = this%gsize(ig) + 1
              this%groups(this%gsize(ig), ig) = ir
              exit
            end if
          else
            ! if the group is empty 
            ! set this processor as the head of a group
            this%gsize(ig) = this%gsize(ig) + 1
            this%groups(1, ig) = ir
            this%ng = this%ng + 1
            exit
          end if
        end do
      end do

      ASSERT(sum(this%gsize(1:this%ng)) == mpi_world%size)

      this%maxgsize = maxval(this%gsize(1:this%ng))

      !convert to mpi ranks
      this%groups = this%groups - 1

#endif
    end subroutine topology_init

    logical function topology_groups_are_equal(this) result(are_equal)
      type(topology_t), intent(in) :: this
      
      are_equal = all(this%gsize(2:this%ng) == this%maxgsize)
    end function topology_groups_are_equal

    subroutine topology_end(this)
      type(topology_t), intent(inout) :: this

      SAFE_DEALLOCATE_P(this%groups)
      SAFE_DEALLOCATE_P(this%distance)
      SAFE_DEALLOCATE_P(this%gsize)
      
    end subroutine topology_end

  !---------------------------------------------------
  ! Function to divide the range of numbers from 1 to nn
  ! between size processors.
  subroutine multicomm_divide_range(nn, tsize, start, final, lsize)
    integer, intent(in)    :: nn
    integer, intent(in)    :: tsize
    integer, intent(out)   :: start(:)
    integer, intent(out)   :: final(:)
    integer, intent(out)   :: lsize(:)

    integer :: ii, jj, rank
    
    if(tsize <= nn ) then
      
      do rank = 0, tsize - 1
        jj = nn / tsize
        ii = nn - jj*tsize
        if(ii > 0 .and. rank < ii) then
          jj = jj + 1
          start(rank + 1) = rank*jj + 1
          final(rank + 1) = start(rank + 1) + jj - 1
        else
          final(rank + 1) = nn - (tsize - rank - 1)*jj
          start(rank + 1) = final(rank + 1) - jj + 1
        end if
      end do

    else
      start = 1
      final = 0
      start(1) = 1
      final(1) = nn
    end if

    lsize(1:tsize) = final(1:tsize) - start(1:tsize) + 1
    
    ASSERT(sum(lsize(1:tsize)) == nn)

  end subroutine multicomm_divide_range

#ifdef USE_OMP
  subroutine multicomm_divide_range_omp(nn, ini, nn_loc)
    integer, intent(in)    :: nn
    integer, intent(out)   :: ini
    integer, intent(out)   :: nn_loc
    
    integer :: start(MAX_OMP_THREADS), end(MAX_OMP_THREADS), lsize(MAX_OMP_THREADS), rank

    call multicomm_divide_range(nn, omp_get_num_threads(), start, end, lsize)

    rank   = 1 + omp_get_thread_num()
    ini    = start(rank)
    nn_loc = lsize(rank)

  end subroutine multicomm_divide_range_omp
#endif

end module multicomm_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
