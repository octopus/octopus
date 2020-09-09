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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

  !>/* To use the profiling module you simply have to define a
  !!profile_t object (with the save attribute). To initialize it you
  !!can use the method profile_init (the second argument is the label
  !!of the profile) or use implicit initialization, by passing the
  !!optional label argument to profiling_in.
  !!If not "save", when the object goes out of scope, it breaks the linked list prof_vars%profile_list. 
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
module profiling_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use namespace_oct_m
  use nvtx_oct_m
  use sort_oct_m
  use string_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none
  private

  public ::                             &
    profile_t,                          &
    profile_pointer_t,                  &
    profile_vars_t,                     &
    profile_is_initialized,             &
    profiling_end,                      &
    profiling_init,                     &
    profiling_in,                       &
    profiling_out,                      &
    profiling_count_operations,         &
    profiling_count_transfers,          &
    profiling_memory_allocate,          &
    profiling_memory_deallocate,        &
    profiling_output

  integer, parameter ::                 & 
       LABEL_LENGTH = 25,               &  !< Max. number of characters of tag label.
       MAX_PROFILES = 200                  !< Max. number of tags.
  
  type profile_t
    private
    character(LABEL_LENGTH)  :: label
    FLOAT                    :: entry_time
    FLOAT                    :: total_time
    FLOAT                    :: min_time
    FLOAT                    :: self_time
    FLOAT                    :: op_count_current
    FLOAT                    :: op_count
    FLOAT                    :: op_count_child
    FLOAT                    :: op_count_child_current
    FLOAT                    :: tr_count_current
    FLOAT                    :: tr_count
    FLOAT                    :: tr_count_child
    FLOAT                    :: tr_count_child_current
    type(profile_t), pointer :: parent
    integer                  :: count
    logical                  :: initialized = .false.
    logical                  :: active = .false.
    logical                  :: exclude
    integer                  :: index
    logical                  :: has_child(MAX_PROFILES)
    FLOAT                    :: timings(MAX_PROFILES)
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
         profiling_count_tran_complex_8,   &
         profiling_count_tran_type
  end interface profiling_count_transfers

  interface profiling_count_operations
    module procedure iprofiling_count_operations
    module procedure rprofiling_count_operations
    module procedure dprofiling_count_operations
  end interface profiling_count_operations
 
  integer, parameter, public  ::  &
       PROFILING_TIME        = 1, &
       PROFILING_MEMORY      = 2, &
       PROFILING_MEMORY_FULL = 4, &
       PROFILING_LIKWID      = 8, &
       PROFILING_IO          = 16

  integer, parameter :: MAX_MEMORY_VARS = 25

  type profile_vars_t
    private
    integer, public          :: mode    !< 1=time, 2=memory, 4=memory_full

    type(profile_pointer_t)  :: current !< the currently active profile
    type(profile_pointer_t)  :: profile_list(MAX_PROFILES) !< the list of all profiles
    integer                  :: last_profile

    integer(8)               :: alloc_count
    integer(8)               :: dealloc_count

    integer(8)               :: memory_limit
    integer(8)               :: total_memory
    integer(8)               :: max_memory
    character(len=256)       :: max_memory_location

    integer(8)               :: large_vars_size(MAX_MEMORY_VARS)
    character(len=256)       :: large_vars(MAX_MEMORY_VARS)

    FLOAT                    :: start_time
    integer                  :: mem_iunit

    character(len=256)       :: output_dir
    character(len=6)         :: file_number

    logical                  :: all_nodes

    logical                  :: output_yaml
    logical                  :: output_tree
  end type profile_vars_t

  type(profile_vars_t), target, public :: prof_vars

  !> For the moment we will have the profiler objects here, but they
  !! should be moved to their respective modules.
  !! i.e. DO NOT PUT NEW PROFILES HERE

  type(profile_t), save, public :: C_PROFILING_COMPLETE_RUN

contains

  ! ---------------------------------------------------------
  !> Create profiling subdirectory.
  subroutine profiling_init(namespace)
    type(namespace_t),          intent(in)    :: namespace
    
    integer :: ii

    PUSH_SUB(profiling_init)

    ! FIXME: nothing is thread-safe here!
    
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
    !% profiling. Warning: you may encounter strange results with OpenMP.
    !%Option no 0
    !% No profiling information is generated.
    !%Option prof_time 1
    !% Profile the time spent in defined profiling regions.
    !%Option prof_memory 2
    !% As well as the time, summary information on memory usage and the largest arrays are reported.
    !%Option prof_memory_full 4
    !% As well as the time and summary memory information, a
    !% log is reported of every allocation and deallocation.
    !%Option likwid 8
    !% Enable instrumentation using LIKWID.
    !%Option prof_io 16
    !% Count the number of file open and close.
    !%End

    call parse_variable(namespace, 'ProfilingMode', 0, prof_vars%mode)
    if(.not.varinfo_valid_option('ProfilingMode', prof_vars%mode)) then
      call messages_input_error(namespace, 'ProfilingMode')
    end if

    in_profiling_mode = (prof_vars%mode > 0)
    if(.not.in_profiling_mode) then
      POP_SUB(profiling_init)
      return
    end if

    !%Variable ProfilingAllNodes
    !%Default no
    !%Type logical
    !%Section Execution::Optimization
    !%Description
    !% This variable controls whether all nodes print the time
    !% profiling output. If set to no, the default, only the root node
    !% will write the profile. If set to yes, all nodes will print it.
    !%End

    call parse_variable(namespace, 'ProfilingAllNodes', .false., prof_vars%all_nodes)

    call get_output_dir()

    if(bitand(prof_vars%mode, PROFILING_MEMORY_FULL) /= 0) then
      prof_vars%mode = ior(prof_vars%mode, PROFILING_MEMORY)
    end if

    ! initialize memory profiling
    if(bitand(prof_vars%mode, PROFILING_MEMORY) /= 0) then
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
      call parse_variable(namespace, 'MemoryLimit', -1, ii)
      prof_vars%memory_limit = int(ii, 8)*1024
    end if

    if(bitand(prof_vars%mode, PROFILING_MEMORY_FULL) /= 0) then
      ! make sure output directory is available before other processes try to write there
#ifdef HAVE_MPI
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
      
      prof_vars%mem_iunit = io_open(trim(prof_vars%output_dir)//'/memory.'//prof_vars%file_number, &
        namespace, action='write')
      write(prof_vars%mem_iunit, '(5a16,a70)') 'Elapsed Time', 'Alloc/Dealloc', 'Size (words)', 'Prof Mem', &
        'Sys Mem', 'Variable Name(Filename:Line)'
    end if

    ! initialize time profiling
    prof_vars%last_profile = 0
    nullify(prof_vars%current%p)

    if(bitand(prof_vars%mode, PROFILING_LIKWID) /= 0) then
#ifdef HAVE_LIKWID
      call likwid_markerInit()
#endif
    end if

    !%Variable ProfilingOutputYAML
    !%Default no
    !%Type logical
    !%Section Execution::Optimization
    !%Description
    !% This variable controls whether the profiling output is additionally
    !% written to a YAML file.
    !%End
    call parse_variable(namespace, 'ProfilingOutputYAML', .false., prof_vars%output_yaml)

    !%Variable ProfilingOutputTree
    !%Default no
    !%Type logical
    !%Section Execution::Optimization
    !%Description
    !% This variable controls whether the profiling output is additionally
    !% written as a tree.
    !%End
    call parse_variable(namespace, 'ProfilingOutputTree', .false., prof_vars%output_tree)

    call profiling_in(C_PROFILING_COMPLETE_RUN, 'COMPLETE_RUN')

    POP_SUB(profiling_init)

  contains

    ! ---------------------------------------------------------
    subroutine get_output_dir()

      PUSH_SUB(profiling_init.get_output_dir)

      write(prof_vars%file_number, '(i6.6)') mpi_world%rank

      prof_vars%output_dir = 'profiling'

      if(mpi_grp_is_root(mpi_world)) call io_mkdir(trim(prof_vars%output_dir), namespace)

      POP_SUB(profiling_init.get_output_dir)
    end subroutine get_output_dir

  end subroutine profiling_init


  ! ---------------------------------------------------------
  subroutine profiling_end(namespace)
    type(namespace_t), intent(in) :: namespace
    integer :: ii
    FLOAT, parameter :: megabyte = CNST(1048576.0)
    integer(8) :: io_open_count, io_close_count
#ifdef HAVE_MPI
    integer(8) :: io_open_count_red, io_close_count_red
#endif

    if(.not. in_profiling_mode) return
    PUSH_SUB(profiling_end)

    call profiling_out(C_PROFILING_COMPLETE_RUN)
    call profiling_output(namespace)

    do ii = 1, prof_vars%last_profile
      prof_vars%profile_list(ii)%p%initialized = .false.
    end do

    if(bitand(prof_vars%mode, PROFILING_MEMORY) /= 0) then
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

      if(prof_vars%alloc_count /= prof_vars%dealloc_count) then
        write(message(1),'(a,i10,a,i10,a)') "Not all memory was deallocated: ", prof_vars%alloc_count, &
          ' allocations and ', prof_vars%dealloc_count, ' deallocations'
        call messages_warning(1, all_nodes = .true.)
      end if
      if(prof_vars%total_memory > 0) then
        write(message(1),'(a,f18.3,a,f18.3,a)') "Remaining allocated memory: ", prof_vars%total_memory/megabyte, &
          ' Mbytes (out of maximum ', prof_vars%max_memory/megabyte, ' Mbytes)'
        call messages_warning(1, all_nodes = .true.)
      end if
    end if

    if(bitand(prof_vars%mode, PROFILING_MEMORY_FULL) /= 0) then
      call io_close(prof_vars%mem_iunit)
    end if

    if(bitand(prof_vars%mode, PROFILING_LIKWID) /= 0) then
#ifdef HAVE_LIKWID
      call likwid_markerClose()
#endif
    end if

    if(bitand(prof_vars%mode, PROFILING_IO) /= 0) then
      call messages_print_stress(stdout, "IO profiling information")
      io_open_count = io_get_open_count()
      io_close_count = io_get_close_count()
      write(message(1), '(a,i10)') 'Number of file open  = ', io_open_count
      write(message(2), '(a,i10)') 'Number of file close = ', io_close_count
#ifdef HAVE_MPI
      call MPI_Allreduce(io_open_count, io_open_count_red, 1, MPI_INTEGER8, MPI_SUM, &
                            mpi_world%comm, mpi_err)
      call MPI_Allreduce(io_close_count, io_close_count_red, 1, MPI_INTEGER8, MPI_SUM, &
                            mpi_world%comm, mpi_err)
      write(message(3), '(a,i10)') 'Global number of file open  = ', io_open_count_red
      write(message(4), '(a,i10)') 'Global number of file close = ', io_close_count_red
      call messages_info(4)
#else
      call messages_info(2)
#endif
      call messages_print_stress(stdout)
    end if

    POP_SUB(profiling_end)
  end subroutine profiling_end


  ! ---------------------------------------------------------
  !> Initialize a profile object and add it to the list
  subroutine profile_init(this, label)
    type(profile_t), target, intent(out)   :: this
    character(*),            intent(in)    :: label 
    
    PUSH_SUB(profile_init)

    this%label = label
    this%total_time = M_ZERO
    this%min_time   = M_HUGE
    this%self_time  = M_ZERO
    this%entry_time = huge(this%entry_time)
    this%count  = 0
    this%op_count_current      = M_ZERO
    this%op_count              = M_ZERO
    this%op_count_child        = M_ZERO
    this%tr_count_current      = M_ZERO
    this%tr_count              = M_ZERO
    this%tr_count_child        = M_ZERO
    this%active = .false.
    nullify(this%parent)
    this%has_child = .false.
    this%timings = M_ZERO
    this%index = 0

    if(.not. in_profiling_mode) then
      POP_SUB(profile_init)
      return
    end if

    prof_vars%last_profile = prof_vars%last_profile + 1

    ASSERT(prof_vars%last_profile  <=  MAX_PROFILES)

    prof_vars%profile_list(prof_vars%last_profile)%p => this
    this%index = prof_vars%last_profile
    this%initialized = .true.

    POP_SUB(profile_init)
  end subroutine profile_init


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
  subroutine profiling_in(this, label, exclude)
    type(profile_t), target,    intent(inout) :: this
    character(*),               intent(in)    :: label
    logical,         optional,  intent(in)    :: exclude !< .true. The time spent here is also excluded from the parent total_time.
                                                         !! Only use it for functions that otherwise would spoil statistics.

    FLOAT :: now

    if(.not.in_profiling_mode) return
    if(.not. not_in_openmp()) return

    ! no PUSH_SUB, called too often

    if(.not. this%initialized) then 
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
      this%parent%has_child(this%index) = .true.
    else 
      !we are orphans
      nullify(this%parent)
    end if

    this%op_count_current = M_ZERO
    this%tr_count_current = M_ZERO
    this%op_count_child_current = M_ZERO
    this%tr_count_child_current = M_ZERO

    prof_vars%current%p => this
    this%entry_time = now

    this%exclude = optional_default(exclude, .false.)

    if(bitand(prof_vars%mode, PROFILING_LIKWID) /= 0) then
#ifdef HAVE_LIKWID
      call likwid_markerStartRegion(trim(label))
#endif
    end if

#ifdef HAVE_NVTX
    call nvtx_range_push(trim(label))
#endif

  end subroutine profiling_in


  ! ---------------------------------------------------------
  !> Increment out counter and sum up difference between entry
  !! and exit time.
  !!
  subroutine profiling_out(this)
    type(profile_t),   intent(inout) :: this

    FLOAT :: now, time_spent

    if(.not.in_profiling_mode) return
    if(.not. not_in_openmp()) return

    ! no PUSH_SUB, called too often

    ASSERT(this%initialized)
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
    if (time_spent < this%min_time) then
      this%min_time = time_spent
    end if

    this%op_count = this%op_count + this%op_count_current
    this%tr_count = this%tr_count + this%tr_count_current
    this%op_count_child = this%op_count_child + this%op_count_child_current
    this%tr_count_child = this%tr_count_child + this%tr_count_child_current
    
    if(associated(this%parent)) then 
      !remove the spent from the self time of our parent
      this%parent%self_time = this%parent%self_time - time_spent
      if(this%exclude) this%parent%total_time = this%parent%total_time - time_spent

      ! add the operations to the parent
      this%parent%op_count_child_current = this%parent%op_count_child_current &
        + this%op_count_current + this%op_count_child_current
      this%parent%tr_count_child_current = this%parent%tr_count_child_current &
        + this%tr_count_current + this%tr_count_child_current

      this%parent%timings(this%index) = this%parent%timings(this%index) + time_spent

      !and set parent as current
      prof_vars%current%p => this%parent

    else
      nullify(prof_vars%current%p)
    end if

    if(bitand(prof_vars%mode, PROFILING_LIKWID) /= 0) then
#ifdef HAVE_LIKWID
      call likwid_markerStopRegion(trim(this%label))
#endif
    end if

#ifdef HAVE_NVTX
    call nvtx_range_pop()
#endif

  end subroutine profiling_out


  ! ---------------------------------------------------------

  subroutine iprofiling_count_operations(ops)
    integer,         intent(in)    :: ops

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often

    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + TOFLOAT(ops)
  end subroutine iprofiling_count_operations


  ! ---------------------------------------------------------

  subroutine rprofiling_count_operations(ops)
    real(4),         intent(in)    :: ops

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + TOFLOAT(ops)
  end subroutine rprofiling_count_operations


  ! ---------------------------------------------------------
 
  subroutine dprofiling_count_operations(ops)
    real(8),         intent(in)    :: ops

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%op_count_current = prof_vars%current%p%op_count_current + ops

  end subroutine dprofiling_count_operations


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_int(trf, type)
    integer,         intent(in)    :: trf
    integer,         intent(in)    :: type

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often    

    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(4*trf)
  end subroutine profiling_count_tran_int


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_real_4(trf, type)
    integer,         intent(in)    :: trf
    real(4),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often

    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(4*trf)

  end subroutine profiling_count_tran_real_4


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_real_8(trf, type)
    integer,         intent(in)    :: trf
    real(8),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(8*trf)

  end subroutine profiling_count_tran_real_8


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_complex_4(trf, type)
    integer,         intent(in)    :: trf
    complex(4),      intent(in)    :: type

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(8*trf)

  end subroutine profiling_count_tran_complex_4


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_complex_8(trf, type)
    integer,         intent(in)    :: trf
    complex(8),      intent(in)    :: type

    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(16*trf)

  end subroutine profiling_count_tran_complex_8


  ! ---------------------------------------------------------
  
  subroutine profiling_count_tran_type(trf, type)
    integer,         intent(in)    :: trf
    type(type_t),    intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    ! no PUSH_SUB, called too often
    
    prof_vars%current%p%tr_count_current = prof_vars%current%p%tr_count_current + TOFLOAT(trf)*types_get_size(type)

  end subroutine profiling_count_tran_type
  
  ! ---------------------------------------------------------
  FLOAT function profile_total_time(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_total_time)
    profile_total_time = this%total_time

    POP_SUB(profile_total_time)
  end function profile_total_time


  ! ---------------------------------------------------------
  FLOAT function profile_self_time(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_self_time)
    profile_self_time = this%self_time

    POP_SUB(profile_self_time)
  end function profile_self_time


  ! ---------------------------------------------------------
  FLOAT function profile_total_time_per_call(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_total_time_per_call)
    profile_total_time_per_call = this%total_time / TOFLOAT(this%count)

    POP_SUB(profile_total_time_per_call)
  end function profile_total_time_per_call


  ! ---------------------------------------------------------
  FLOAT function profile_min_time(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_self_time)
    profile_min_time = this%min_time

    POP_SUB(profile_self_time)
  end function profile_min_time


  ! ---------------------------------------------------------
  FLOAT function profile_self_time_per_call(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_self_time_per_call)
    profile_self_time_per_call = this%self_time / TOFLOAT(this%count)

    POP_SUB(profile_self_time_per_call)
  end function profile_self_time_per_call


  ! ---------------------------------------------------------
  real(8) function profile_total_throughput(this) 
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_throughput)

    if(this%total_time > epsilon(this%total_time)) then
      profile_total_throughput = (this%op_count + this%op_count_child)/this%total_time*CNST(1.0e-6)
    else
      profile_total_throughput = CNST(0.0)
    end if
      
    POP_SUB(profile_throughput)
  end function profile_total_throughput


  ! ---------------------------------------------------------
  
  FLOAT function profile_total_bandwidth(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_bandwidth)

    if(this%total_time > epsilon(this%total_time)) then
      profile_total_bandwidth = (this%tr_count + this%tr_count_child)/(this%total_time*CNST(1024.0)**2)
    else
      profile_total_bandwidth = CNST(0.0)
    end if

    POP_SUB(profile_bandwidth)
  end function profile_total_bandwidth
  
  ! ---------------------------------------------------------
  
  FLOAT function profile_self_throughput(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_throughput)

    if(this%self_time > epsilon(this%self_time)) then
      profile_self_throughput = this%op_count/this%self_time*CNST(1.0e-6)
    else
      profile_self_throughput = CNST(0.0)
    end if
      
    POP_SUB(profile_throughput)
  end function profile_self_throughput

  ! ---------------------------------------------------------

  FLOAT function profile_self_bandwidth(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_bandwidth)

    if(this%self_time > epsilon(this%self_time)) then
      profile_self_bandwidth = this%tr_count/(this%self_time*CNST(1024.0)**2)
    else
      profile_self_bandwidth = CNST(0.0)
    end if

    POP_SUB(profile_bandwidth)
  end function profile_self_bandwidth


  ! ---------------------------------------------------------
  integer function profile_num_calls(this)
    type(profile_t), intent(in) :: this
    
    PUSH_SUB(profile_num_calls)
    profile_num_calls = this%count

    POP_SUB(profile_num_calls)
  end function profile_num_calls


  ! ---------------------------------------------------------
  character(LABEL_LENGTH) function profile_label(this)
    type(profile_t), intent(in) :: this

    PUSH_SUB(profile_label)
    profile_label = this%label

    POP_SUB(profile_label)
  end function profile_label


  ! ---------------------------------------------------------
  !> Write profiling results of each node to profiling.NNN/profiling.nnn
  !! The format of each line is
  !! tag-label    pass_in    pass_out    time   time/pass_in
  !!
  !! The last column gives the average time consumed between in and out
  !! (only, if pass_in and pass_out are equal).
  subroutine profiling_output(namespace)
    type(namespace_t), intent(in) :: namespace
    
    integer          :: ii
    integer          :: iunit
    FLOAT            :: total_time
    type(profile_t), pointer :: prof
    character(len=256) :: filename
    FLOAT,   allocatable :: selftime(:)
    integer, allocatable :: position(:)
    
    if(.not.in_profiling_mode) return
    PUSH_SUB(profiling_output)

#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    if(.not. prof_vars%all_nodes .and. .not. mpi_grp_is_root(mpi_world)) then
      POP_SUB(profiling_output)
      return
    end if

    filename = trim(prof_vars%output_dir)//'/time.'//prof_vars%file_number
    iunit = io_open(trim(filename), namespace, action='write')
    if(iunit < 0) then
      message(1) = 'Failed to open file ' // trim(filename) // ' to write profiling results.'
      call messages_warning(1)
      POP_SUB(profiling_output)
      return
    end if

    write(iunit, '(2a)')                                                                                    &
      '                                                                    CUMULATIVE TIME                ', &
      '                 |                         SELF TIME'
    write(iunit, '(2a)')                                                                                    &
      '                                          ----------------------------------------------------------', &
      '----------------|-------------------------------------------------------------'
    write(iunit, '(2a)')                                                                                    &
      'TAG                           NUM_CALLS      TOTAL_TIME   TIME_PER_CALL        MIN_TIME   ', &
      ' MFLOPS  MBYTES/S   %TIME |       TOTAL_TIME   TIME_PER_CALL    MFLOPS  MBYTES/S   %TIME'
    write(iunit, '(2a)')                                                                    &
      '===================================================================================================', &
      '=================|============================================================='

    total_time = profile_total_time(C_PROFILING_COMPLETE_RUN)

    SAFE_ALLOCATE(selftime(1:prof_vars%last_profile))
    SAFE_ALLOCATE(position(1:prof_vars%last_profile))

    do ii = 1, prof_vars%last_profile
      selftime(ii) = -profile_self_time(prof_vars%profile_list(ii)%p)
      position(ii) = ii
    end do

    call sort(selftime, position)    
    
    do ii = 1, prof_vars%last_profile
      prof =>  prof_vars%profile_list(position(ii))%p
      if(.not. prof%initialized) then
        write(message(1),'(a,i6,a)') "Internal error: Profile number ", position(ii), " is not initialized."
        call messages_fatal(1)
      end if
      if(prof%active) then
        write(message(1),'(a)') "Internal error: Profile '" // trim(profile_label(prof)) // &
          "' is active, i.e. profiling_out was not called."
        call messages_warning(1)
      end if
      
      if(profile_num_calls(prof) == 0) cycle

      write(iunit, '(a,i14,3f16.6,2f10.1,f8.1,a,2f16.6,2f10.1,f8.1)')     &
           profile_label(prof),                             & 
           profile_num_calls(prof),                         &
           profile_total_time(prof),                        &
           profile_total_time_per_call(prof),               &
           profile_min_time(prof),                          &
           profile_total_throughput(prof),                  &
           profile_total_bandwidth(prof),                   &
           profile_total_time(prof)/total_time*CNST(100.0), &
           ' | ',                                           &
           profile_self_time(prof),                         &
           profile_self_time_per_call(prof),                &
           profile_self_throughput(prof),                   &
           profile_self_bandwidth(prof),                    &           
           profile_self_time(prof)/total_time*CNST(100.0)
    end do

    call io_close(iunit)

    if(prof_vars%output_yaml) then
      filename = trim(prof_vars%output_dir)//'/time.'//prof_vars%file_number//'.yaml'
      iunit = io_open(trim(filename), namespace, action='write')
      if(iunit < 0) then
        message(1) = 'Failed to open file ' // trim(filename) // ' to write profiling results.'
        call messages_warning(1)
        POP_SUB(profiling_output)
        return
      end if
      write(iunit, '(2a)') 'schema: [num_calls, total_time, total_throughput, ', &
       'total_bandwidth, self_time, self_throughput, self_bandwidth]'
      write(iunit, '(a)') 'data:'

      do ii = 1, prof_vars%last_profile
        prof =>  prof_vars%profile_list(position(ii))%p
        if(profile_num_calls(prof) == 0) cycle
        write(iunit, '(a,a,a,i6,a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,e10.3,a,e10.3,a)')         &
             '  ', profile_label(prof), ': [',     &
             profile_num_calls(prof),        ', ', &
             profile_total_time(prof),       ', ', &
             profile_total_throughput(prof), ', ', &
             profile_total_bandwidth(prof),  ', ', &
             profile_self_time(prof),        ', ', &
             profile_self_throughput(prof),  ', ', &
             profile_self_bandwidth(prof),   ']'
      end do

      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(selftime)
    SAFE_DEALLOCATE_A(position)

    if(prof_vars%output_tree) then
      filename = trim(prof_vars%output_dir)//'/time.'//prof_vars%file_number//'.tree'
      iunit = io_open(trim(filename), namespace, action='write')
      if(iunit < 0) then
        message(1) = 'Failed to open file ' // trim(filename) // ' to write profiling results.'
        call messages_warning(1)
        POP_SUB(profiling_output)
        return
      end if
      write(iunit, '(a39,a25,a12)')         &
        "Tree level, % of total, % of parent    ", &
        "Region                    ", &
        "   Full time"

      ! output of top-level node
      write(iunit, '(a,f7.2,a,f7.2,a,a,a25,f12.4)')         &
           repeat('-', 0) // '| ',  100.0, "%  ", 100.0, "% ", &
           repeat(' ', 18), profile_label(C_PROFILING_COMPLETE_RUN), &
           total_time
      call output_tree_level(C_PROFILING_COMPLETE_RUN, 1, total_time, iunit)
      write(iunit, '(a)') "// modeline for vim to enable folding (put in ~/.vimrc: set modeline modelineexpr)"
      write(iunit, '(a)') "// vim: fdm=expr fde=getline(v\:lnum)=~'.*\|.*'?len(split(getline(v\:lnum))[0])-1\:0"
      call io_close(iunit)
    end if
    
    POP_SUB(profiling_output)
    contains
      ! Traverse the tree depth-first, pre-order
      ! 
      ! Because the children are added at initialization of the profiles,
      ! each profile is the children of the first parent it occured as the
      ! initialization happens at the first call of profile_init. This means
      ! that a profile at deeper levels of the tree can contain data from
      ! different parents. This is not captured by the simple tree data
      ! structure because it is in principle a graph. Displaying that would
      ! be much more difficult, however.
      recursive subroutine output_tree_level(profile, level, total_time, iunit)
        type(profile_t), intent(in) :: profile
        integer,         intent(in) :: level
        FLOAT,           intent(in) :: total_time
        integer,         intent(in) :: iunit

        integer :: ichild, width

        PUSH_SUB(profiling_output.output_tree_level)
        width = 20
        ! loop over children
        do ichild = 1, MAX_PROFILES
          if (profile%has_child(ichild)) then
            ! print out information on current child with the first marker
            ! placed according to the level of the tree
            write(iunit, '(a,f7.2,a,f7.2,a,a,a25,f12.4)')         &
                 repeat('-', level) // '| ', &
                 profile%timings(ichild)/total_time * 100, "%  ", &
                 profile%timings(ichild)/profile%total_time * 100, "% ", &
                 repeat(' ', width-level-2), &
                 profile_label(prof_vars%profile_list(ichild)%p), &
                 profile%timings(ichild)
            call output_tree_level(prof_vars%profile_list(ichild)%p, &
              level+1, total_time, iunit)
          end if
        end do
        POP_SUB(profiling_output.output_tree_level)
      end subroutine output_tree_level
  end subroutine profiling_output


  ! ---------------------------------------------------------
  subroutine profiling_make_position_str(var, file, line, str)
    character(len=*), intent(in)  :: var
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: line
    character(len=*), intent(out) :: str

    integer            :: ii, jj, nn
    
    ! no push_sub, called too many times

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
    ii = 1
    do while ( file(ii:ii+2) == "../" ) 
      ii = ii + 3
    end do
    write(str, '(4a,i5,a)') var(1:jj), "(", trim(file(ii:len(file))), ":", line, ")"
    call compact(str)

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
    
    ! no push_sub, called too many times

    call profiling_make_position_str(var, file, line, str)

    ! get number of pages
    mem = loct_get_memory_usage()

    write(prof_vars%mem_iunit, '(f16.6,a16,3i16,a70)') loct_clock() - prof_vars%start_time, &
         trim(type), size, prof_vars%total_memory, mem, trim(str)

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

    ! no push_sub, called too many times

    size = size_ ! make a copy that we can change

    prof_vars%alloc_count  = prof_vars%alloc_count + 1
    prof_vars%total_memory = prof_vars%total_memory + size

    if(bitand(prof_vars%mode, PROFILING_MEMORY_FULL) /= 0) then 
      call profiling_memory_log('A ', var, file, line, size)
    end if

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
      if(str == prof_vars%large_vars(ii)) then
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
    
  end subroutine profiling_memory_allocate


  !-----------------------------------------------------
  subroutine profiling_memory_deallocate(var, file, line, size)
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: file
    integer,          intent(in) :: line
    integer(8),       intent(in) :: size
    
    ! no push_sub, called too many times

    prof_vars%dealloc_count  = prof_vars%dealloc_count + 1
    prof_vars%total_memory   = prof_vars%total_memory - size

    if(bitand(prof_vars%mode, PROFILING_MEMORY_FULL) /= 0) then 
      call profiling_memory_log('D ', var, file, line, -size)
    end if

  end subroutine profiling_memory_deallocate


end module profiling_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
