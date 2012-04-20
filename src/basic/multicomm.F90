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

!> Given a number n_index of indices and their ranges index_range(1:n_index),
!! we divide the n_nodes in groups and create communicators for each group.
!! Each group is associated with one index. The min_range indicates the minimum
!! number of elements in each processor. For example, given min_range(1) = 25000,
!! the algorithm will try to put at least 25000 points in each processor
!!
!! Example:
!! \verbatim
!! Given 3 indices with ranges
!!   index_range(1) = 100000, (number of points in the mesh)
!!   index_range(2) = 15,     (number of states)
!!   index_range(3) = 2,      (number of k-points)
!! \endverbatim
!! and 12 processors, we could get
!! \verbatim
!!   mc%group_sizes = (2, 3, 2)
!! \endverbatim
!! which means that
!! - space is divided in 2 domains per state,
!! - the states are divided in 3 groups, i.e. 5 states per processor, and
!! - the whole setting is duplicated because of the 2 k-points.
!!
!! To perform collective operations (like a reduce), you can use the communicators
!! provided in mc%group_comm(:). For example, to sum over states, the communicator
!! to use is mc%group_comm(P_STRATEGY_STATES)
!!
!! You can use the routine multicomm_strategy_is_parallel to know if a certain
!! index is parallelized.

module multicomm_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
#if defined(HAVE_OPENMP)
  use omp_lib
#endif
  use parser_m
  use profiling_m
  use utils_m
  use varinfo_m
  
  implicit none
  
  private
  
  public ::                          &
    multicomm_divide_range,          &
#ifdef HAVE_OPENMP
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
    multicomm_is_slave,              &
    multicomm_have_slaves

  ! possible parallelization strategies
  integer, public, parameter ::      &
    P_STRATEGY_SERIAL  = 0,          & !< single domain, all states, k-points on a single processor
    P_STRATEGY_DOMAINS = 1,          & !< parallelization in domains
    P_STRATEGY_STATES  = 2,          & !< parallelization in states
    P_STRATEGY_KPOINTS = 3,          & !< parallelization in k-points
    P_STRATEGY_OTHER   = 4             !< something else like e-h pairs

  integer, public, parameter ::      &
    P_MASTER           = 1,          &
    P_SLAVE            = 2
  
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

  !> Stores all communicators and groups
  type multicomm_t
    integer          :: n_node           !< Total number of nodes.
    integer          :: n_index          !< Number of parallel indices.

    integer          :: par_strategy     !< What kind of parallelization strategy should we use?

    integer, pointer :: group_sizes(:)   !< Number of processors in each group.
    integer, pointer :: who_am_i(:)      !< Rank in the "line"-communicators.
    integer, pointer :: group_comm(:)    !< "Line"-communicators I belong to.
    integer          :: dom_st_comm      !< States-domain plane communicator.
    integer          :: st_kpt_comm      !< Kpoints-states plane communicator.
    integer          :: dom_st_kpt_comm  !< Kpoints-states-domain cube communicator.

    integer          :: nthreads         !< Number of OMP threads
    integer          :: node_type        !< Is this node a P_MASTER or a P_SLAVE?
    logical          :: have_slaves      !< are slaves available?
    
    integer          :: full_comm        !< The base communicator.
    integer          :: full_comm_rank   !< The rank in the base communicator.
    integer          :: master_comm      !< The communicator without slaves.
    integer          :: master_comm_rank !< The rank in the communicator without slaves.
    integer          :: slave_intercomm  !< the intercomm to communicate with slaves
  end type multicomm_t

  !> An all-pairs communication schedule for a given group.
  type multicomm_all_pairs_t
    type(mpi_grp_t)  :: grp            !< Schedule for this group.
    integer          :: rounds         !< This many comm. rounds.
    integer, pointer :: schedule(:, :) !< This is the schedule.
  end type multicomm_all_pairs_t

contains

  ! ---------------------------------------------------------
  subroutine multicomm_all_pairs_copy(apout, apin)
    type(multicomm_all_pairs_t), intent(inout) :: apout
    type(multicomm_all_pairs_t), intent(in)    :: apin

    PUSH_SUB(multicomm_all_pairs_copy)

    call mpi_grp_copy(apout%grp, apin%grp)
    apout%rounds = apin%rounds
    if(associated(apin%schedule)) then
      SAFE_ALLOCATE(apout%schedule(1:size(apin%schedule, 1), 1:size(apin%schedule, 2)))
      apout%schedule = apin%schedule
    end if    

    POP_SUB(multicomm_all_pairs_copy)
  end subroutine multicomm_all_pairs_copy

  ! ---------------------------------------------------------
  !> create index and domain communicators
  subroutine multicomm_init(mc, base_grp, parallel_mask, default_mask, n_node, n_index, index_range, min_range)
    type(multicomm_t), intent(out)   :: mc
    type(mpi_grp_t),   intent(inout) :: base_grp
    integer,           intent(in)    :: parallel_mask, default_mask, n_node, n_index
    integer,           intent(inout) :: index_range(:)
    integer,           intent(in)    :: min_range(:)

    integer :: ii, num_slaves, slave_level
    type(block_t) :: blk

    PUSH_SUB(multicomm_init)

    ASSERT(n_index <= n_par_types)

    mc%n_node  = n_node
    mc%n_index = n_index  ! size(index_range)

    call messages_print_stress(stdout, "Parallelization")

    call strategy()

    if (mc%nthreads > 1) then
      write(message(1),'(a, i3)') 'Info: Number of threads ', mc%nthreads
      call messages_info(1)
    end if

    nullify(mc%group_sizes)
    mc%have_slaves = .false.
    if(mc%par_strategy.ne.P_STRATEGY_SERIAL) then
      SAFE_ALLOCATE(mc%group_sizes(1:mc%n_index))
      mc%group_sizes(:) = 1

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
      if(parse_block(datasets_check('ParallelizationGroupRanks'), blk) == 0) then
        call read_block(blk)
        call parse_block_end(blk)
      else
        call assign_nodes()
      end if

      !%Variable ParallelizationNumberSlaves
      !%Type integer
      !%Section Execution::Parallelization
      !%Description
      !% Slaves are nodes used for task parallelization. The number of
      !% such nodes is given by this variable multiplied by the number
      !% of domains used in domain parallelization. The default is 0.
      !%End
      call parse_integer(datasets_check('ParallelizationNumberSlaves'), 0, num_slaves)
      
      ! the slaves must be defined at a certain parallelization level, for the moment this is state parallelization.
      slave_level = P_STRATEGY_STATES
      mc%have_slaves = (num_slaves > 0)

      if(mc%have_slaves) then
#ifndef HAVE_MPI2
        message(1) = 'Task parallelization requires MPI 2.'
        call messages_fatal(1)
#endif
        call messages_experimental('Task parallelization')
      end if

      ! clear parallel strategies that were available but will not be used
      do ii = 1, mc%n_index
        if(mc%group_sizes(ii) == 1) mc%par_strategy = ibclr(mc%par_strategy, ii - 1)
      end do

      ! reset
      call sanity_check()
    end if

    call group_comm_create()

    call messages_print_stress(stdout)

    POP_SUB(multicomm_init)

  contains

    ! ---------------------------------------------------------
    subroutine strategy()
      integer :: jj,  par_mask

      PUSH_SUB(multicomm_init.strategy)

      !%Variable ParallelizationStrategy
      !%Type flag
      !%Section Execution::Parallelization
      !%Description
      !% Specifies what kind of parallelization strategy <tt>Octopus</tt> should use.
      !% The values can be combined: for example, <tt>par_domains + par_states</tt>
      !% means a combined parallelization in domains and states.
      !% Default: <tt>par_domains + par_states</tt> for <tt>CalculationMode = td</tt>,
      !% <tt>par_domains + par_other</tt> for <tt>CalculationMode = casida</tt>,
      !% otherwise <tt>par_domains</tt>.
      !%Option serial 0
      !% <tt>Octopus</tt> will run in serial.
      !%Option par_domains 1
      !% <tt>Octopus</tt> will run parallel in domains.
      !%Option par_states  2
      !% <tt>Octopus</tt> will run parallel in states.
      !%Option par_kpoints 4
      !% <tt>Octopus</tt> will run parallel in <i>k</i>-points/spin.
      !%Option par_other   8
      !% Run-mode-dependent. For example, in <tt>casida</tt>, it means parallelization in <i>e-h</i> pairs.
      !%End

      ! default is set in calc_mode_default_parallel_mask()

     if(base_grp%size > 1) then

        par_mask = parallel_mask

        call parse_integer(datasets_check('ParallelizationStrategy'), default_mask, mc%par_strategy)

        if(.not.varinfo_valid_option('ParallelizationStrategy', mc%par_strategy, is_flag=.true.)) then
          call input_error('ParallelizationStrategy')
        end if

        mc%par_strategy = iand(mc%par_strategy, par_mask)

        if(mc%par_strategy == P_STRATEGY_SERIAL) then
          message(1) = "More than one node is available, but this run mode cannot run in parallel."
          message(2) = "Please select a ParallelizationStrategy compatible with"
          jj = 2
          do ii = 1, n_par_types
            if(iand(par_mask, 2**(ii - 1)) .ne. 0) then
              jj = jj + 1
              write(message(jj), '(2a)') "  -> ", par_types(ii)
            end if
          end do
          jj=jj+1
          write(message(jj),'(a,i6)') "mc%par_strategy is : ",mc%par_strategy
          call messages_fatal(jj, only_root_writes = .true.)
        end if
      else
        mc%par_strategy = P_STRATEGY_SERIAL
      end if

       mc%nthreads = 1
#if defined(HAVE_OPENMP)
      !$omp parallel
      !$omp master
      mc%nthreads = omp_get_num_threads()
      !$omp end master
      !$omp end parallel
      if(mc%nthreads > MAX_OMP_THREADS) then
        message(1) = "Number of threads requested is larger than MAX_OMP_THREADS."
        call messages_fatal(1)
      end if
#endif

      if(mc%par_strategy == P_STRATEGY_SERIAL .and. mc%nthreads == 1) then
        message(1) = "Octopus will run in *serial*"
      else
        message(1) = "Octopus will run in *parallel*"
      end if
      call messages_info(1)
      
      POP_SUB(multicomm_init.strategy)
    end subroutine strategy


    ! ---------------------------------------------------------
    subroutine read_block(blk)
      type(block_t), intent(inout) :: blk

      integer :: ii, nn
      logical :: fill_used

      PUSH_SUB(multicomm_init.read_block)

      nn = parse_block_cols(blk, 0)

      mc%group_sizes = 1
      do ii = 1, min(nn, mc%n_index)
        if(multicomm_strategy_is_parallel(mc, ii)) then
          call parse_block_integer(blk, 0, ii - 1, mc%group_sizes(ii))
        else
          message(1) = 'In ParallelizationGroupRanks, ignoring specification for ' // par_types(ii)
          message(2) = 'This parallelization strategy is not available.'
          call messages_warning(2)
        end if
      end do

      fill_used = .false.
      do ii = 1, mc%n_index
        if(mc%group_sizes(ii) == -1) then
          if(fill_used) then
            message(1) = "Error: The 'fill' value can be used only once in ParallelizationGroupRanks."
            call messages_fatal(1, only_root_writes = .true.)
          end if
          mc%group_sizes(ii) = -base_grp%size / product(mc%group_sizes)
          fill_used = .true.
        end if
      end do

      POP_SUB(multicomm_init.read_block)
    end subroutine read_block

    ! ---------------------------------------------------------
    subroutine assign_nodes()
      integer :: ii, nn, kk, n_divisors, divisors(50)
      integer, allocatable :: n_group_max(:)
      FLOAT   :: ff

      PUSH_SUB(multicomm_init.assign_nodes)

      SAFE_ALLOCATE(n_group_max(1:mc%n_index))

      nn = mc%n_node

      ! this is the maximum number of processors in each group
      n_group_max(1:mc%n_index) = max(index_range(1:mc%n_index), 1)
      do kk = 1, mc%n_index
        if(.not. multicomm_strategy_is_parallel(mc, kk)) n_group_max(kk) = 1
      end do

      ! for each index
      do kk = 1, mc%n_index
        if(n_group_max(kk) == 1) then ! not parallel in this group
          mc%group_sizes(kk) = 1
          cycle
        end if

        ! distibute the nodes so that domains is the last
        ff = real(nn, REAL_PRECISION)
        do ii = kk + 1, mc%n_index
          ff = ff / real(n_group_max(ii), REAL_PRECISION)
        end do

        ! get divisors of nn
        n_divisors = 50 ! maximum number of divisors
        call get_divisors(nn, n_divisors, divisors)

        ! get the divisor of nn >= ff
        mc%group_sizes(kk) = nn
        do ii = 1, n_divisors
          if(real(divisors(ii), REAL_PRECISION) >= ff) then
            mc%group_sizes(kk) = divisors(ii)
            exit
          end if
        end do

        nn = nn / mc%group_sizes(kk)
      end do

      SAFE_DEALLOCATE_A(n_group_max)

      POP_SUB(multicomm_init.assign_nodes)
    end subroutine assign_nodes


    ! ---------------------------------------------------------
    !> check if a balanced distribution of nodes will be used
    subroutine sanity_check()
      FLOAT :: frac
      integer :: ii, kk, n_max
      integer :: real_group_sizes(1:MAX_INDEX)

      PUSH_SUB(multicomm_init.sanity_check)

      if(num_slaves > 0) then

        if(mc%group_sizes(slave_level) < num_slaves + 1) then
          message(1) = 'Too many nodes assigned to task parallelization.'
          call messages_fatal(1)
        end if

        write(message(1),'(a,i6)') 'Info: Number of slaves nodes              :', &
          num_slaves*product(mc%group_sizes(1:slave_level - 1))
        call messages_info(1)

      end if

      ! print out some info
      ii = 0
      do kk = mc%n_index, 1, -1
        real_group_sizes(kk) = mc%group_sizes(kk)
        if(.not. multicomm_strategy_is_parallel(mc, kk)) cycle
        ii = ii + 1
        if(kk == slave_level) INCR(real_group_sizes(kk), -num_slaves)
        write(message(ii),'(3a,i6,a,i8,a)') 'Info: Number of nodes in ', &
          par_types(kk), ' group:', real_group_sizes(kk), ' (', index_range(kk), ')'
      end do
      call messages_info(ii)

      ! do we have the correct number of processors
      if(product(mc%group_sizes(1:mc%n_index)) .ne. base_grp%size) then
        write(message(1),'(a)') 'Inconsistent number of processors:'
        write(message(2),'(a,i6)') '  MPI processes      = ', base_grp%size
        write(message(3),'(a,i6)') '  Required processes = ', product(real_group_sizes(1:mc%n_index))
        message(4) = ''
        message(5) = 'You probably have a problem in the ParallelizationGroupRanks block.'
        call messages_fatal(5, only_root_writes = .true.)
      end if

      if(any(real_group_sizes(1:mc%n_index) > index_range(1:mc%n_index))) then
        message(1) = "Could not distribute nodes in parallel job. Most likely you are trying to"
        message(2) = "use too many nodes for the job."
        call messages_fatal(2, only_root_writes = .true.)
      end if

      if(any(index_range(1:mc%n_index) / real_group_sizes(1:mc%n_index) < min_range(1:mc%n_index))) then
        message(1) = "I have fewer elements in a parallel group than recommended."
        message(2) = "Maybe you should reduce the number of nodes."
        call messages_warning(2)
      end if

      ! calculate fraction of idle time
      frac = M_ONE
      do ii = 1, mc%n_index
        n_max = ceiling(real(index_range(ii), REAL_PRECISION) / real(real_group_sizes(ii), REAL_PRECISION))
        kk = n_max*real_group_sizes(ii)
        frac = frac*(M_ONE - real(kk - index_range(ii), REAL_PRECISION) / real(kk, REAL_PRECISION))
      end do

      write(message(1), '(a,f5.2,a)') "Info: Octopus will waste at least ", &
        (M_ONE - frac)*CNST(100.0), "% of computer time."
      if(frac < CNST(0.8)) then
        message(2) = "Usually a number of processors which is a multiple of small primes is best."
        call messages_warning(3)
      else
        call messages_info(1)
      end if

      POP_SUB(multicomm_init.sanity_check)
    end subroutine sanity_check

    ! ---------------------------------------------------------
    subroutine group_comm_create()
#if defined(HAVE_MPI)
      logical :: dim_mask(MAX_INDEX)
      integer :: i_strategy
      logical :: reorder, periodic_mask(MAX_INDEX)
      integer :: coords(MAX_INDEX)
      integer :: new_comm
#endif

      PUSH_SUB(multicomm_init.group_comm_create)

      mc%node_type = P_MASTER

      SAFE_ALLOCATE(mc%group_comm(1:mc%n_index))
      SAFE_ALLOCATE(mc%who_am_i(1:mc%n_index))

#if defined(HAVE_MPI)
      mc%full_comm = MPI_COMM_NULL
      mc%slave_intercomm = MPI_COMM_NULL
      if(mc%par_strategy /= P_STRATEGY_SERIAL) then
        ! Multilevel parallelization is organized in a hypercube. We
        ! use an MPI Cartesian topology to generate the communicators
        ! that correspond to each level.

        ! create the topology
        periodic_mask = .false.
        reorder = .true.

        ! The domain and states dimensions have to be periodic (2D torus)
        ! in order to circulate matrix blocks.
        periodic_mask(P_STRATEGY_DOMAINS) = multicomm_strategy_is_parallel(mc, P_STRATEGY_DOMAINS)
        periodic_mask(P_STRATEGY_STATES) = multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES)

        ! We allow reordering of ranks.
        call MPI_Cart_create(base_grp%comm, mc%n_index, mc%group_sizes, periodic_mask, reorder, mc%full_comm, mpi_err)

        call MPI_Comm_rank(mc%full_comm, mc%full_comm_rank, mpi_err)

        ! get the coordinates of the current processor
        call MPI_Cart_coords(mc%full_comm, mc%full_comm_rank, mc%n_index, coords, mpi_err)

        ! find out what type of node this is
        if(coords(slave_level) >= mc%group_sizes(slave_level) - num_slaves) then
          mc%node_type = P_SLAVE
        end if

        if(mc%node_type == P_MASTER) then
          INCR(mc%group_sizes(slave_level), -num_slaves)
        else
          mc%group_sizes(slave_level) = num_slaves
        end if

        call MPI_Comm_split(mc%full_comm, mc%node_type, mc%full_comm_rank, new_comm, mpi_err)

        reorder = .false.
        call MPI_Cart_create(new_comm, mc%n_index, mc%group_sizes, periodic_mask, reorder, mc%master_comm, mpi_err)

        call MPI_Comm_free(new_comm, mpi_err)

        call MPI_Comm_rank(mc%master_comm, mc%master_comm_rank, mpi_err)

        ! The "lines" of the Cartesian grid.
        do i_strategy = 1, mc%n_index
          if(multicomm_strategy_is_parallel(mc, i_strategy)) then
            dim_mask             = .false.
            dim_mask(i_strategy) = .true.
            call MPI_Cart_sub(mc%master_comm, dim_mask, mc%group_comm(i_strategy), mpi_err)
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
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%dom_st_comm, mpi_err)

        ! The state-kpoints "planes" of the grid
        dim_mask                     = .false.
        dim_mask(P_STRATEGY_STATES)  = .true.
        dim_mask(P_STRATEGY_KPOINTS) = .true.
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%st_kpt_comm, mpi_err)

        ! The domains-states-kpoints "cubes" of the grid
        dim_mask                     = .false.
        dim_mask(P_STRATEGY_DOMAINS) = .true.
        dim_mask(P_STRATEGY_STATES)  = .true.
        dim_mask(P_STRATEGY_KPOINTS) = .true.
        call MPI_Cart_sub(mc%master_comm, dim_mask, mc%dom_st_kpt_comm, mpi_err)

        if(num_slaves > 0) call create_slave_intercommunicators()
      else
        ! we initialize these communicators so we can use them even in serial
        mc%group_comm = base_grp%comm
        mc%who_am_i   = 0
        mc%master_comm = base_grp%comm
        mc%dom_st_comm = base_grp%comm
        mc%st_kpt_comm = base_grp%comm
        mc%dom_st_kpt_comm = base_grp%comm
      end if
#else
      mc%group_comm = -1
      mc%who_am_i   = 0
      mc%full_comm = -1
      mc%master_comm = -1
      mc%dom_st_comm = -1
      mc%st_kpt_comm = -1
      mc%dom_st_kpt_comm = -1
      mc%slave_intercomm = -1
#endif

      POP_SUB(multicomm_init.group_comm_create)
    end subroutine group_comm_create

    ! -----------------------------------------------------

    subroutine create_slave_intercommunicators()
#ifdef HAVE_MPI2
      integer :: remote_leader
      integer :: tag
      integer :: coords(MAX_INDEX)

      PUSH_SUB(multicomm_init.create_slave_intercommunicators)

      ! create the intercommunicators to communicate with slaves

      ! get the coordinates of the current processor
      call MPI_Cart_coords(mc%full_comm, mc%full_comm_rank, mc%n_index, coords, mpi_err)
      
      !now get the rank of the remote_leader
      if(mc%node_type == P_SLAVE) then
        coords(slave_level) = 0
      else
        coords(slave_level) = mc%group_sizes(slave_level)
      end if
      call MPI_Cart_rank(mc%full_comm, coords, remote_leader, mpi_err)

      ! now create the intercommunicator
      tag = coords(P_STRATEGY_DOMAINS)
      call MPI_Intercomm_create(mc%group_comm(slave_level), 0, base_grp%comm, remote_leader, tag, mc%slave_intercomm, mpi_err)

      POP_SUB(multicomm_init.create_slave_intercommunicators)
#endif
    end subroutine create_slave_intercommunicators

  end subroutine multicomm_init
  

  ! ---------------------------------------------------------
    subroutine multicomm_end(mc)
      type(multicomm_t), intent(inout) :: mc

#if defined(HAVE_MPI)
      integer :: ii
#endif

      PUSH_SUB(multicomm_end)

    if(mc%par_strategy .ne. P_STRATEGY_SERIAL) then
#if defined(HAVE_MPI)
      ! Delete communicators.
      do ii = 1, mc%n_index
        if(.not. multicomm_strategy_is_parallel(mc, ii)) cycle
        call MPI_Comm_free(mc%group_comm(ii), mpi_err)
      end do
      call MPI_Comm_free(mc%dom_st_comm, mpi_err)
      call MPI_Comm_free(mc%st_kpt_comm, mpi_err)
      call MPI_Comm_free(mc%dom_st_kpt_comm, mpi_err)
      call MPI_Comm_free(mc%full_comm, mpi_err)
      call MPI_Comm_free(mc%master_comm, mpi_err)

#ifdef HAVE_MPI2
      if(multicomm_have_slaves(mc)) call MPI_Comm_free(mc%slave_intercomm, mpi_err)
#endif

#endif
      ! Deallocate the rest of the arrays.
      SAFE_DEALLOCATE_P(mc%group_sizes)
    end if

    SAFE_DEALLOCATE_P(mc%group_comm)
    SAFE_DEALLOCATE_P(mc%who_am_i)
    
    POP_SUB(multicomm_end)
  end subroutine multicomm_end


  ! ---------------------------------------------------------
  logical pure function multicomm_strategy_is_parallel(mc, level) result(rr)
    type(multicomm_t), intent(in) :: mc
    integer,           intent(in) :: level

    rr = iand(mc%par_strategy, 2**(level - 1)) .ne. 0

  end function multicomm_strategy_is_parallel


  ! ---------------------------------------------------------
  !> This routine uses the one-factorization (or near-one-factorization
  !! of a complete graph to construct an all-pair communication
  !! schedule (cf. Wang, X., Blum, E. K., Parker, D. S., and Massey,
  !! D. 1997. The dance-party problem and its application to collective
  !! communication in computer networks. Parallel Comput. 23, 8
  !! (Aug. 1997), 1141-1156.
#if defined(HAVE_MPI)

  subroutine multicomm_create_all_pairs(mpi_grp, ap)
    type(mpi_grp_t),             intent(in)  :: mpi_grp
    type(multicomm_all_pairs_t), intent(out) :: ap

    integer :: grp_size, rounds, ir, in

    PUSH_SUB(create_all_pairs)

    ap%grp = mpi_grp
    grp_size = mpi_grp%size

    ! Number of rounds.
    if(mod(grp_size, 2) .eq. 0) then
      rounds = grp_size - 1
    else
      rounds = grp_size
    end if
    ap%rounds = rounds

    ! Calculate schedule.
    SAFE_ALLOCATE(ap%schedule(0:grp_size - 1, 1:rounds))
    do ir = 1, rounds
      do in = 0, grp_size - 1
        ap%schedule(in, ir) = get_partner(in + 1, ir) - 1
      end do
    end do

    POP_SUB(create_all_pairs)

  contains

    ! ---------------------------------------------------------
    !> Those are from the paper cited above.
    integer pure function get_partner(in, ir)
      integer, intent(in) :: in, ir

      ! No PUSH SUB, called too often.

      if(mod(grp_size, 2).eq.0) then
        get_partner = get_partner_even(grp_size, in - 1, ir - 1) + 1
      else
        get_partner = get_partner_odd(grp_size, in - 1, ir - 1) + 1
      end if

     end function get_partner

    ! ---------------------------------------------------------
    integer pure function get_partner_even(grp_size, ii, rr) result(pp)
      integer, intent(in) :: grp_size, ii, rr

      integer :: mm

      ! No PUSH SUB, called too often.

      mm = grp_size / 2

      if(ii .eq. 0) then
        pp = rr + 1
      elseif(ii .eq. rr + 1) then
        pp = 0
      else
        ! I never know when to use which remainder function, but here
        ! it has to be the modulo one. Do not change that!
        pp = modulo(2 * rr - ii + 1, 2 * mm - 1) + 1
      end if

    end function get_partner_even

    ! ---------------------------------------------------------
      integer pure function get_partner_odd(grp_size, ii, rr) result(pp)
      integer, intent(in) :: grp_size, ii, rr

      integer :: mm

      ! No PUSH SUB, called too often.

      mm = (grp_size + 1) / 2

      pp = get_partner_even(grp_size + 1, ii, rr)

      if(pp .eq. 2 * mm - 1) then
        pp = ii
      end if

    end function get_partner_odd

  end subroutine multicomm_create_all_pairs
#endif

  !---------------------------------------------------
  !> Function to divide the range of numbers from 1 to nobjs
  !! between nprocs processors.
  !! THREADSAFE
  subroutine multicomm_divide_range(nobjs, nprocs, istart, ifinal, lsize, scalapack_compat)
    integer,           intent(in)    :: nobjs !< Number of points to divide
    integer,           intent(in)    :: nprocs !< Number of processors
    integer,           intent(out)   :: istart(:)
    integer,           intent(out)   :: ifinal(:)
    integer, optional, intent(out)   :: lsize(:) !< Number of objects in each partition
    logical, optional, intent(in)    :: scalapack_compat

    integer :: ii, jj, rank, size
    logical :: scalapack_compat_
#ifdef HAVE_SCALAPACK
    integer :: nbl
#endif
    
    scalapack_compat_ = .false.
#ifdef HAVE_SCALAPACK
    if(present(scalapack_compat)) scalapack_compat_ = scalapack_compat
#endif
    ! no push_sub, threadsafe
    if(scalapack_compat_) then
#ifdef HAVE_SCALAPACK      
      nbl = nobjs/nprocs
      if (mod(nobjs, nprocs) /= 0) INCR(nbl, 1)
      
      istart(1) = 1
      do rank = 1, nprocs
        size = numroc(nobjs, nbl, rank - 1, 0, nprocs)
        if(size > 0) then
          if(rank > 1) istart(rank) = ifinal(rank - 1) + 1
          ifinal(rank) = istart(rank) + size - 1
        else
          istart(rank) = 1
          ifinal(rank) = 0
        endif
      end do
#endif
    else
      
      if(nprocs <= nobjs) then

        do rank = 0, nprocs - 1
          jj = nobjs / nprocs
          ii = nobjs - jj*nprocs
          if(ii > 0 .and. rank < ii) then
            jj = jj + 1
            istart(rank + 1) = rank*jj + 1
            ifinal(rank + 1) = istart(rank + 1) + jj - 1
          else
            ifinal(rank + 1) = nobjs - (nprocs - rank - 1)*jj
            istart(rank + 1) = ifinal(rank + 1) - jj + 1
          end if
        end do

      else
        do ii = 1, nprocs
          if(ii <= nobjs) then
            istart(ii) = ii
            ifinal(ii) = ii
          else
            istart(ii) = 1
            ifinal(ii) = 0
          end if
        end do
      end if
    end if

    if(present(lsize)) then
      lsize(1:nprocs) = ifinal(1:nprocs) - istart(1:nprocs) + 1
      ASSERT(sum(lsize(1:nprocs)) == nobjs)
    end if

  end subroutine multicomm_divide_range

  ! ---------------------------------------------------------
  !> Function to divide the range of numbers from 1 to nobjs
  !! between all available threads with OpenMP.
  ! THREADSAFE
  subroutine multicomm_divide_range_omp(nobjs, ini, nobjs_loc)
    integer, intent(in)    :: nobjs       !< Number of points to divide
    integer, intent(out)   :: ini         !< Start point of the partition
    integer, intent(out)   :: nobjs_loc   !< Number of objects in each partition
    
    integer :: istart(MAX_OMP_THREADS), ifinal(MAX_OMP_THREADS), lsize(MAX_OMP_THREADS), rank

    ! no push_sub, threadsafe
    rank = 1
#ifdef HAVE_OPENMP
    call multicomm_divide_range(nobjs, omp_get_num_threads(), istart, ifinal, lsize)

    rank   = 1 + omp_get_thread_num()
#endif
    ini    = istart(rank)
    nobjs_loc = lsize(rank)

  end subroutine multicomm_divide_range_omp

  ! ---------------------------------------------------------

  logical pure function multicomm_is_slave(this) result(slave)
    type(multicomm_t), intent(in) :: this
    
    slave = this%node_type == P_SLAVE
  end function multicomm_is_slave

  ! ---------------------------------------------------------

  logical pure function multicomm_have_slaves(this) result(have_slaves)
    type(multicomm_t), intent(in) :: this

    have_slaves = this%have_slaves
  end function multicomm_have_slaves

end module multicomm_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
