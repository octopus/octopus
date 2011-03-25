!! Copyright (C) 2005-2009 Heiko Appel, Florian Lorenzen, Xavier Andrade
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

  !>/* To use the profiling module you simply have to define a
  !!profile_t object (with the save attribute). To initialize it you
  !!can use the method profile_init (the second argument is the label
  !!of the profile) or use implicit initialization, by passing the
  !!optional label argument to profiling_in.
  !!
  !!To profile, call profiling_in and profiling_out around the regions
  !!you want to profile. The objects need not be destroyed as this is
  !!done by profiling_end.
  !!
  !!This module works in the following way:
  !!
  !!Each profile_t object measures two times:
  !!
  !! total time: time spent inside the profiling region
  !! self time: the total time minus the time spent in profiling
  !! regions called inside
  !!
  !!To implement this there is a current pointer that holds the
  !!innermost region active. When a new profile_t is activated, it puts
  !!itself as current and stores the previous current pointer as
  !!parent. When it is deactivated it sets back current as its parent
  !!and subtracts the time spent from the parent`s self time.
  !!
  !!This scheme will fail with recursive calls, but we don`t use that
  !!and I don`t think we will need it.
  !!
  !!To be able to print the profiling results, there is an array with
  !!pointers to all initialized profile_t. For simplicity this array is
  !!static with size MAX_PROFILES, eventually it should be
  !!incremented. It could be replaced by a linked list, but I don`t
  !!think this is necessary.
  !*/
module profiling_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use mpi_m
  use parser_m
#ifdef HAVE_PAPI
  use papi_m
#endif
  use string_m
  use varinfo_m

  implicit none
  private

  public ::                             &
    profile_t,                          &
    profile_pointer_t,                  &
    profile_vars_t,                     &
    profile_init,                       &
    profile_is_initialized,             &
    profiling_init,                     &
    profiling_end,                      &
    profiling_in,                       &
    profiling_out,                      &
    profiling_count_operations,         &
    profiling_count_transfers,          &
    profiling_memory_allocate,          &
    profiling_memory_deallocate,        &
    profiling_output

  integer, parameter ::                 & 
       LABEL_LENGTH = 17,               &  !< Max. number of characters of tag label.
       MAX_PROFILES = 200                  !< Max. number of tags.
  
  type profile_t
    private
    character(LABEL_LENGTH)  :: label
    real(8)                  :: entry_time
    real(8)                  :: total_time
    real(8)                  :: self_time
    real(8)                  :: op_count_current
    real(8)                  :: op_count
    real(8)                  :: tr_count_current
    real(8)                  :: tr_count
    type(profile_t), pointer :: parent
    integer                  :: count
    logical                  :: initialized = .false.
    logical                  :: active
  end type profile_t

  type profile_pointer_t
    private
    type(profile_t), pointer :: p
  end type profile_pointer_t
  
  interface profiling_count_transfers
    module procedure &
         profiling_count_tran_int,         &
         profiling_count_tran_real_4,      &
         profiling_count_tran_real_8,      &
         profiling_count_tran_complex_4,   &
         profiling_count_tran_complex_8
  end interface

  interface profiling_count_operations
    module procedure iprofiling_count_operations
    module procedure rprofiling_count_operations
    module procedure dprofiling_count_operations
  end interface

  integer, parameter, public  ::   &
       PROFILING_TIME        = 1, &
       PROFILING_MEMORY      = 2, &
       PROFILING_MEMORY_FULL = 4

#define MAX_MEMORY_VARS 25

  type profile_vars_t
    integer                  :: mode    !< 1=time, 2=memory, 4=memory_full

    type(profile_pointer_t)  :: current !<the currently active profile
    type(profile_pointer_t)  :: profile_list(MAX_PROFILES) !<the list of all profilers
    integer                  :: last_profile

    integer(8)               :: alloc_count
    integer(8)               :: dealloc_count

    integer(8)               :: memory_limit
    integer(8)               :: total_memory
    integer(8)               :: max_memory
    character(len=256)       :: max_memory_location

    integer(8)               :: large_vars_size(MAX_MEMORY_VARS)
    character(len=256)       :: large_vars(MAX_MEMORY_VARS)

    real(8)                  :: start_time
    integer                  :: mem_iunit

    character(len=256)       :: output_dir
    character(len=6)         :: file_number

    logical                  :: all_nodes
  end type profile_vars_t

  type(profile_vars_t), public :: prof_vars

  !For the moment we will have the profiler objects here, but they
  !should be moved to their respective modules.
  !i.e. DO NOT PUT NEW PROFILES HERE

  type(profile_t), save, public :: C_PROFILING_COMPLETE_DATASET

  type(profile_t), save, public :: &
       C_PROFILING_XC_OEP,         &
       C_PROFILING_XC_EXX,         &
       C_PROFILING_XC_SIC,         &
       C_PROFILING_XC_OEP_FULL,    &
       C_PROFILING_XC_KLI

  type(profile_t), save, public ::    &
       C_PROFILING_BLOCKT,            &
       C_PROFILING_BLOCKT_AR,         &
       C_PROFILING_BLOCKT_MM,         &
       C_PROFILING_BLOCKT_CP,         &
       C_PROFILING_BLOCK_MATR,        &
       C_PROFILING_BLOCK_MATR_CP,     &
       C_PROFILING_BLOCK_MATR_MM,     &
       C_PROFILING_LOBPCG_ESOLVE,     &
       C_PROFILING_LOBPCG_CHOL,       &
       C_PROFILING_LOBPCG_INV,        &
       C_PROFILING_LOBPCG_COPY,       &
       C_PROFILING_LOBPCG_LOOP

contains

  ! ---------------------------------------------------------
  !> Create profiling subdirectory.
  subroutine profiling_init()
    integer :: ii

    PUSH_SUB(profiling_init)

    !%Variable ProfilingMode
    !%Default no
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% Use this variable to run <tt>Octopus</tt> in profiling mode. In this mode
    !% <tt>Octopus</tt> records the time spent in certain areas of the code and
    !% the number of times this code is executed. These numbers
    !% are written in <tt>./profiling.NNN/profiling.nnn</tt> with <tt>nnn</tt> being the
    !% node number (<tt>000</tt> in serial) and <tt>NNN</tt> the number of processors.
    !% This is mainly for development purposes. Note, however, that
    !% <tt>Octopus</tt> should be compiled with <tt>--disable-debug</tt> to do proper
    !% profiling.
    !%Option no 0
    !% No profiling information is generated.
    !%Option prof_time 1
    !% Profile the time spent in defined profiling regions.
    !%Option prof_memory 2
    !% As well as the time, memory usage is reported.
    !%Option prof_memory_full 4
    !% As well as the time, full memory usage is reported.
    !%End

    call parse_integer('ProfilingMode', 0, prof_vars%mode)
    if(.not.varinfo_valid_option('ProfilingMode', prof_vars%mode, is_flag=.true.)) then
      call input_error('ProfilingMode')
    end if

    in_profiling_mode = (prof_vars%mode > 0)
    if(.not.in_profiling_mode) then
    POP_SUB(profiling_init)
      return
    end if

    !%Variable ProfilingAllNodes
    !%Default no
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% This variable controls whether all nodes print the time
    !% profiling output. If set to no, the default, only the root node
    !% will write the profile. If set to yes all nodes will print it.
    !%End

    call parse_logical('ProfilingAllNodes', .false., prof_vars%all_nodes)

    call get_output_dir()

#ifdef HAVE_PAPI
    call papi_init()
#endif

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      prof_vars%mode = ior(prof_vars%mode, PROFILING_MEMORY)
    end if

    ! initialize memory profiling
    if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) then
      prof_vars%alloc_count   = 0
      prof_vars%dealloc_count = 0

      prof_vars%total_memory  = 0
      prof_vars%max_memory    = 0
      prof_vars%max_memory_location = ''
      prof_vars%start_time = loct_clock()
      
      prof_vars%large_vars_size(:) = 0
      prof_vars%large_vars(:) = ''

      !%Variable MemoryLimit
      !%Default -1
      !%Type integer
      !%Section Execution::Optimization
      !%Description
      !% If positive, <tt>Octopus</tt> will stop if more memory than <tt>MemoryLimit</tt>
      !% is requested (in kb). Note that this variable only works when 
      !% <tt>ProfilingMode = prof_memory(_full)</tt>.
      !%End
      call parse_integer('MemoryLimit', -1, ii)
      prof_vars%memory_limit = int(ii, 8)*1024
    end if

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      prof_vars%mem_iunit = io_open(trim(prof_vars%output_dir)//'/memory.'//prof_vars%file_number, action='write')
    end if

    ! initialize time profiling
    prof_vars%last_profile = 0
    nullify(prof_vars%current%p)

    call init_profiles()

    call profiling_in(C_PROFILING_COMPLETE_DATASET)

    POP_SUB(profiling_init)

  contains
    ! ---------------------------------------------------------
    subroutine init_profiles
      PUSH_SUB(profiling_init.init_profiles)

      call profile_init(C_PROFILING_COMPLETE_DATASET, 'COMPLETE_DATASET')
      call profile_init(C_PROFILING_XC_OEP,           'XC_OEP')
      call profile_init(C_PROFILING_XC_EXX,           'XC_EXX')
      call profile_init(C_PROFILING_XC_SIC,           'XC_SIC')
      call profile_init(C_PROFILING_XC_KLI,           'XC_KLI')
      call profile_init(C_PROFILING_XC_OEP_FULL,      'XC_OEP_FULL')
      call profile_init(C_PROFILING_BLOCKT,           'BLOCKT')
      call profile_init(C_PROFILING_BLOCKT_AR,        'BLOCKT_AR')
      call profile_init(C_PROFILING_BLOCKT_MM,        'BLOCKT_MM')
      call profile_init(C_PROFILING_BLOCKT_CP,        'BLOCKT_CP')
      call profile_init(C_PROFILING_BLOCK_MATR,       'BLOCK_MATR')
      call profile_init(C_PROFILING_BLOCK_MATR_CP,    'BLOCK_MATR_CP')
      call profile_init(C_PROFILING_BLOCK_MATR_MM,    'BLOCK_MATR_MM')
      call profile_init(C_PROFILING_LOBPCG_ESOLVE,    'LOBPCG_ESOLVE')
      call profile_init(C_PROFILING_LOBPCG_CHOL,      'LOBPCG_CHOL')
      call profile_init(C_PROFILING_LOBPCG_INV,       'LOBPCG_INV')
      call profile_init(C_PROFILING_LOBPCG_COPY,      'LOBPCG_COPY')
      call profile_init(C_PROFILING_LOBPCG_LOOP,      'LOBPCG_LOOP')

      POP_SUB(profiling_init.init_profiles)
    end subroutine init_profiles

    ! ---------------------------------------------------------
    subroutine get_output_dir()
      character(len=6) :: dirnum

      PUSH_SUB(profiling_init.get_output_dir)

      dirnum  = 'ser '
      prof_vars%file_number = '0000'
#if defined(HAVE_MPI)
      if(mpi_world%size > 1) then
        write(dirnum, '(i6.6)') mpi_world%size
        write(prof_vars%file_number, '(i6.6)') mpi_world%rank
      end if
#endif
      prof_vars%output_dir = 'profiling.'//trim(dirnum)

      if(mpi_grp_is_root(mpi_world)) call io_mkdir(trim(prof_vars%output_dir))

      POP_SUB(profiling_init.get_output_dir)
    end subroutine get_output_dir

  end subroutine profiling_init


  ! ---------------------------------------------------------
  subroutine profiling_end
    integer :: ii
    real(8), parameter :: megabyte = 1048576.0_8

    if(.not. in_profiling_mode) return
    PUSH_SUB(profiling_end)

    call profiling_out(C_PROFILING_COMPLETE_DATASET)
    call profiling_output()

#ifdef HAVE_PAPI
    call papi_end()
#endif

    do ii = 1, prof_vars%last_profile
      call profile_end(prof_vars%profile_list(ii)%p)
    end do

    if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) then
      call messages_print_stress(stdout, "Memory profiling information")
      write(message(1), '(a,i10)') 'Number of   allocations = ', prof_vars%alloc_count
      write(message(2), '(a,i10)') 'Number of deallocations = ', prof_vars%dealloc_count
      write(message(3), '(a,f18.3,a)') 'Maximum total memory allocated = ', prof_vars%max_memory/megabyte, ' Mbytes'
      write(message(4), '(2x,a,a)') 'at ', trim(prof_vars%max_memory_location)
      call messages_info(4)

      message(1) = ''
      message(2) = 'Largest variables allocated:'
      call messages_info(2)
      do ii = 1, MAX_MEMORY_VARS
        write(message(1),'(i2,f18.3,2a)') ii, prof_vars%large_vars_size(ii)/megabyte, ' Mbytes ', trim(prof_vars%large_vars(ii))
        call messages_info(1)
      end do

      call messages_print_stress(stdout)

      if(prof_vars%alloc_count.ne.prof_vars%dealloc_count) then
        message(1) = "Not all memory was deallocated!";
        call messages_warning(1)
      end if
    end if

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      call io_close(prof_vars%mem_iunit)
    end if

    POP_SUB(profiling_end)
  end subroutine profiling_end


  ! ---------------------------------------------------------
  !> Initialize a profile object and add it to the list
  subroutine profile_init(this, label)
    type(profile_t), target, intent(out)   :: this
    character(*),            intent(in)    :: label 
    
    this%label = label
    this%total_time = M_ZERO
    this%self_time  = M_ZERO
    this%entry_time = huge(this%entry_time)
    this%count  = 0
    this%op_count_current = M_ZERO
    this%op_count = M_ZERO
    this%tr_count = M_ZERO
    this%active = .false.
    nullify(this%parent)

    if(.not. in_profiling_mode) return

    prof_vars%last_profile = prof_vars%last_profile + 1

    ASSERT(prof_vars%last_profile .le. MAX_PROFILES)

    prof_vars%profile_list(prof_vars%last_profile)%p => this
    this%initialized = .true.

  end subroutine profile_init

  ! ---------------------------------------------------------
  subroutine profile_end(this)
    type(profile_t), intent(inout) :: this

    PUSH_SUB(profile_end)
    this%initialized = .false.

    POP_SUB(profile_end)
  end subroutine profile_end


  ! ---------------------------------------------------------
  logical function profile_is_initialized(this)
    type(profile_t), intent(in)   :: this

    PUSH_SUB(profile_is_initialized)
    profile_is_initialized = this%initialized

    POP_SUB(profile_is_initialized)
  end function profile_is_initialized


  ! ---------------------------------------------------------
  !> Increment in counter and save entry time.
  !!
  subroutine profiling_in(this, label)
    type(profile_t), target, intent(inout) :: this
    character(*), optional,  intent(in)    :: label 

    real(8) :: now
#ifdef HAVE_PAPI
    real(8) :: ops
#endif

    if(.not.in_profiling_mode) return

    if(.not. this%initialized) then 
      ASSERT(present(label))
      call profile_init(this, label)
    end if

    ASSERT(.not. this%active)
    this%active = .true.
#if defined(HAVE_MPI)
    now = MPI_Wtime()
#else
    now = loct_clock()
#endif
    if(associated(prof_vars%current%p)) then
      !keep a pointer to the parent
      this%parent => prof_vars%current%p
    else 
      !we are orphans
      nullify(this%parent)
    end if
#ifdef HAVE_PAPI
    call papi_get_count_and_reset(ops)
    if(associated(this%parent)) this%parent%op_count_current = this%parent%op_count_current + ops
#endif

    this%op_count_current = M_ZERO
    this%tr_count_current = M_ZERO

    prof_vars%current%p => this
    this%entry_time = now

  end subroutine profiling_in


  ! ---------------------------------------------------------
  !> Increment out counter and sum up difference between entry
  !! and exit time.
  !!
  subroutine profiling_out(this)
    type(profile_t),   intent(inout) :: this

    real(8) :: now, time_spent
#ifdef HAVE_PAPI
    real(8) :: ops
#endif
    
    if(.not.in_profiling_mode) return

    ASSERT(this%active)
    this%active = .false.
#if defined(HAVE_MPI)
    now = MPI_Wtime()
#else
    now = loct_clock()
#endif
    time_spent = now - this%entry_time
    this%total_time = this%total_time + time_spent
    this%self_time  = this%self_time + time_spent
    this%count = this%count + 1

#ifdef HAVE_PAPI
    call papi_get_count_and_reset(ops)
    this%op_count_current = this%op_count_current + ops
#endif

    this%op_count = this%op_count + this%op_count_current
    this%tr_count = this%tr_count + this%tr_count_current

    if(associated(this%parent)) then 
      !remove the spent from the self time of our parent
      this%parent%self_time = this%parent%self_time - time_spent

      ! add the operations to the parent
      this%parent%op_count_current = this%parent%op_count_current + this%op_count_current
      this%parent%tr_count_current = this%parent%tr_count_current + this%tr_count_current

      !and set parent as current
      prof_vars%current%p => this%parent

    else
      nullify(prof_vars%current%p)
    end if

  end subroutine profiling_out


  ! ---------------------------------------------------------

  subroutine iprofiling_count_operations(ops)
    integer,         intent(in)    :: ops
#ifndef HAVE_PAPI
    if(.not.in_profiling_mode) return

    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + dble(ops)
#endif
  end subroutine iprofiling_count_operations


  ! ---------------------------------------------------------

  subroutine rprofiling_count_operations(ops)
    real(4),         intent(in)    :: ops

#ifndef HAVE_PAPI
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + dble(ops)
#endif
  end subroutine rprofiling_count_operations


  ! ---------------------------------------------------------
 
  subroutine dprofiling_count_operations(ops)
    real(8),         intent(in)    :: ops

#ifndef HAVE_PAPI
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + ops
#endif
  end subroutine dprofiling_count_operations


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_int(trf, type)
    integer,         intent(in)    :: trf
    integer,         intent(in)    :: type

    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + dble(4*trf)
  end subroutine profiling_count_tran_int


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_real_4(trf, type)
    integer,         intent(in)    :: trf
    real(4),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + dble(4*trf)
  end subroutine profiling_count_tran_real_4


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_real_8(trf, type)
    integer,         intent(in)    :: trf
    real(8),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + dble(8*trf)
  end subroutine profiling_count_tran_real_8


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_complex_4(trf, type)
    integer,         intent(in)    :: trf
    complex(4),      intent(in)    :: type

    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + dble(8*trf)
  end subroutine profiling_count_tran_complex_4


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_complex_8(trf, type)
    integer,         intent(in)    :: trf
    complex(8),      intent(in)    :: type

    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + dble(16*trf)
  end subroutine profiling_count_tran_complex_8


  ! ---------------------------------------------------------
  real(8) function profile_total_time(this)
    type(profile_t), intent(in) :: this
    profile_total_time = this%total_time
  end function profile_total_time


  ! ---------------------------------------------------------
  real(8) function profile_self_time(this)
    type(profile_t), intent(in) :: this
    profile_self_time = this%self_time
  end function profile_self_time


  ! ---------------------------------------------------------
  real(8) function profile_total_time_per_call(this)
    type(profile_t), intent(in) :: this
    profile_total_time_per_call = this%total_time / dble(this%count)
  end function profile_total_time_per_call


  ! ---------------------------------------------------------
  real(8) function profile_self_time_per_call(this)
    type(profile_t), intent(in) :: this
    profile_self_time_per_call = this%self_time / dble(this%count)
  end function profile_self_time_per_call


  ! ---------------------------------------------------------
  real(8) function profile_throughput(this)
    type(profile_t), intent(in) :: this
    profile_throughput = this%op_count / this%total_time / CNST(1.0e6)
  end function profile_throughput


  ! ---------------------------------------------------------
  real(8) function profile_bandwidth(this)
    type(profile_t), intent(in) :: this
    profile_bandwidth = this%tr_count / (this%total_time*CNST(1024.0)**2)
  end function profile_bandwidth


  ! ---------------------------------------------------------
  integer function profile_num_calls(this)
    type(profile_t), intent(in) :: this
    
    profile_num_calls = this%count
  end function profile_num_calls


  ! ---------------------------------------------------------
  character(LABEL_LENGTH) function profile_label(this)
    type(profile_t), intent(in) :: this
    profile_label = this%label
  end function profile_label


  ! ---------------------------------------------------------
  !> Write profiling results of each node to profiling.NNN/profifling.nnn
  !! The format of each line is
  !! tag-label    pass_in    pass_out    time   time/pass_in
  !!
  !! The last column gives the average time consumed between in and out
  !! (only, if pass_in and pass_out are equal).
  subroutine profiling_output
    integer          :: ii
    integer          :: iunit
    real(8)          :: total_time
    type(profile_t), pointer :: prof

    if(.not.in_profiling_mode) return

    PUSH_SUB(profiling_output)

#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    if(.not. prof_vars%all_nodes .and. .not. mpi_grp_is_root(mpi_world)) then
      POP_SUB(profiling_output)
      return
    end if

    iunit = io_open(trim(prof_vars%output_dir)//'/time.'//prof_vars%file_number, action='write')
    if(iunit.lt.0) then
      message(1) = 'Could not write profiling results.'
      call messages_warning(1)
      POP_SUB(profiling_output)
      return
    end if

    write(iunit, '(2a)')                                                                                    &
      '                                                              CUMULATIVE TIME                    |', &
      '                  SELF TIME'
    write(iunit, '(2a)')                                                                                    &
      '                                            -----------------------------------------------------|', &
      '-------------------------------------------'
    write(iunit, '(2a)')                                                                                    &
      'TAG                   NUMBER_OF_CALLS       TOTAL_TIME    TIME_PER_CALL  MFLOPS MBYTES/S   %TIME |', &
      '        TOTAL_TIME    TIME_PER_CALL  %TIME'
    write(iunit, '(2a)')                                                                    &
      '=================================================================================================|', &
      '==========================================='

    total_time = profile_total_time(C_PROFILING_COMPLETE_DATASET)

    do ii = 1, prof_vars%last_profile
      prof =>  prof_vars%profile_list(ii)%p

      if(profile_num_calls(prof) == 0) cycle

      write(iunit, '(a,i20,2f17.7,f8.1,f9.1,f8.1,a,2f17.7,f8.1)')     &
           profile_label(prof),                             & 
           profile_num_calls(prof),                         &
           profile_total_time(prof),                        &
           profile_total_time_per_call(prof),               &
           profile_throughput(prof),                        &
           profile_bandwidth(prof),                         &
           profile_total_time(prof)/total_time*CNST(100.0), &
           ' | ',                                           &
           profile_self_time(prof),                         &
           profile_self_time_per_call(prof),                &
           profile_self_time(prof)/total_time*CNST(100.0)
    end do

    call io_close(iunit)

    POP_SUB(profiling_output)
  end subroutine profiling_output


  ! ---------------------------------------------------------
  subroutine profiling_make_position_str(var, file, line, str)
    character(len=*), intent(in)  :: var
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: line
    character(len=*), intent(out) :: str

    integer            :: ii, jj, nn
    
    PUSH_SUB(profiling_make_position_str)

    jj = len(var)
    if(var(jj:jj) == ')') then
      nn = 1
      do ii = len(var)-1, 1, -1
        jj = ii - 1
        if(var(ii:ii) == ')') nn = nn + 1
        if(var(ii:ii) == '(') nn = nn - 1
        if(nn == 0) exit
      end do
      if(jj == 0) then
        message(1) = "Internal Error in profiling_memory_log"
        call messages_fatal(1)
      end if
    end if

    write(str, '(4a,i5,a)') var(1:jj), "(", trim(file), ":", line, ")"
    call compact(str)

    POP_SUB(profiling_make_position_str)
  end subroutine profiling_make_position_str


  ! ---------------------------------------------------------
  subroutine profiling_memory_log(type, var, file, line, size)
    character(len=*), intent(in) :: type
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    integer(8),       intent(in) :: size

    character(len=256) :: str
    integer(8) :: mem
    
    PUSH_SUB(profiling_memory_log)

    call profiling_make_position_str(var, file, line, str)

    ! get number of pages
    mem = get_memory_usage()

    write(prof_vars%mem_iunit, '(f16.6,1x,a,3i16,1x,a)') loct_clock() - prof_vars%start_time, &
         trim(type), size, prof_vars%total_memory, mem, trim(str)

    POP_SUB(profiling_memory_log)
  end subroutine profiling_memory_log


  !-----------------------------------------------------
  subroutine profiling_memory_allocate(var, file, line, size_)
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    integer(8),       intent(in) :: size_

    integer :: ii, jj
    integer(8) :: size
    character(len=256) :: str

    PUSH_SUB(profiling_memory_allocate)

    size = size_ ! make a copy that we can change

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then 
      call profiling_memory_log('A ', var, file, line, size)
    end if

    prof_vars%alloc_count  = prof_vars%alloc_count + 1
    prof_vars%total_memory = prof_vars%total_memory + size

    if(prof_vars%memory_limit > 0) then
      if(prof_vars%total_memory > prof_vars%memory_limit) then
        message(1) = "Memory limit set in the input file was passed"
        call messages_fatal(1)
      end if
    end if

    if(prof_vars%total_memory > prof_vars%max_memory) then
      prof_vars%max_memory = prof_vars%total_memory
      call profiling_make_position_str(var, file, line, prof_vars%max_memory_location)
    end if

    call profiling_make_position_str(var, file, line, str)

    ! check if variable is already in stack
    do ii = 1, MAX_MEMORY_VARS
      if(str.eq.prof_vars%large_vars(ii)) then
        if(size > prof_vars%large_vars_size(ii)) then
          ! delete variable by moving stack up
          do jj = ii,  MAX_MEMORY_VARS - 1
            prof_vars%large_vars(jj)      = prof_vars%large_vars(jj + 1)
            prof_vars%large_vars_size(jj) = prof_vars%large_vars_size(jj + 1)
          end do
          prof_vars%large_vars(MAX_MEMORY_VARS) = ""
          prof_vars%large_vars_size(MAX_MEMORY_VARS) = 0
        else
          ! do not consider this variable any longer
          size = -1
        end if
        exit
      end if
    end do

    do ii = 1, MAX_MEMORY_VARS
      if(size > prof_vars%large_vars_size(ii)) then
        ! move the stack one position down
        do jj = MAX_MEMORY_VARS, ii + 1,  -1
          prof_vars%large_vars(jj)      = prof_vars%large_vars(jj - 1)
          prof_vars%large_vars_size(jj) = prof_vars%large_vars_size(jj - 1)
        end do
        prof_vars%large_vars_size(ii) = size
        prof_vars%large_vars(ii) = str
        exit
      end if
    end do
    
    POP_SUB(profiling_memory_allocate)
  end subroutine profiling_memory_allocate


  !-----------------------------------------------------
  subroutine profiling_memory_deallocate(var, file, line, size)
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    integer(8),       intent(in) :: size
    
    PUSH_SUB(profiling_memory_deallocate)

    prof_vars%dealloc_count  = prof_vars%dealloc_count + 1
    prof_vars%total_memory   = prof_vars%total_memory - size

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then 
      call profiling_memory_log('D ', var, file, line, -size)
    end if

    POP_SUB(profiling_memory_deallocate)
  end subroutine profiling_memory_deallocate
 

end module profiling_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
