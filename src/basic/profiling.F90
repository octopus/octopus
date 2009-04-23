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

module profiling_m
  use global_m
  use io_m
  use loct_m
  use loct_parser_m
  use messages_m
  use mpi_m
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


  !/* To use the profiling module you simply have to define a
  !profile_t object (with the save attribute). To initialize it you
  !can use the method profile_init (the second argument is the label
  !of the profile) or use implicit initialization, by passing the
  !optional label argument to profiling_in. 
  !
  !To profile, call profiling_in and profiling_out around the regions
  !you want to profile. The objects need not be destroyed as this is
  !done by profiling_end.
  !
  !This module work in the following way: 
  !
  !Each profile_t object measures two times: 
  !
  ! total time: time spent inside the profiling region
  ! self time: the total time minus the time spent in profiling
  ! regions called inside
  !
  !To implement this there is a current pointer that holds the
  !innermost region active. When a new profile_t is activated, it puts
  !itself as current and stores the previous current pointer as
  !parent. When it is deactivated it sets back current as its parent
  !and subtracts the time spent from the parent`s self time.
  !
  !This scheme will fail with recursive calls, but we don`t use that
  !and I don`t think we will need it.
  !
  !To be able to print the profiling results, there is an array with
  !pointers to all initialized profile_t. For simplicity this array is
  !static with size MAX_PROFILES, eventually it should be
  !incremented. It could be replaced by a linked list, but I don`t
  !think this is necessary.
  !*/

  integer, parameter ::                 & 
       LABEL_LENGTH = 17,               &  ! Max. number of characters of tag label.
       MAX_PROFILES = 200                  ! Max. number of tags.
  
  type profile_t
    private
    character(LABEL_LENGTH)  :: label
    real(8)                  :: entry_time
    real(8)                  :: total_time
    real(8)                  :: self_time
    real(8)                  :: op_count
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
    module procedure iprofiling_count_operations, dprofiling_count_operations
  end interface

  integer, parameter, public  ::   &
       PROFILING_TIME        = 1, &
       PROFILING_MEMORY      = 2, &
       PROFILING_MEMORY_FULL = 4

#define MAX_MEMORY_VARS 25

  type profile_vars_t
    integer                  :: mode    ! 1=time, 2=memory, 4=memory_full

    type(profile_pointer_t)  :: current !the currently active profile
    type(profile_pointer_t)  :: profile_list(MAX_PROFILES) !the list of all profilers
    integer                  :: last_profile

    integer(8)               :: alloc_count
    integer(8)               :: dealloc_count
    integer(8)               :: total_memory
    integer(8)               :: max_memory
    character(len=256)       :: max_memory_location

    integer(8)               :: large_vars_size(MAX_MEMORY_VARS)
    character(len=256)       :: large_vars(MAX_MEMORY_VARS)

    real(8)                  :: start_time
    integer                  :: mem_iunit
  end type profile_vars_t

  type(profile_vars_t), public :: prof_vars

  !For the moment we will have the profiler objects here, but they
  !should be moved to their respective modules.
  !i.e. DO NOT PUT NEW PROFILES HERE

  type(profile_t), save, public :: C_PROFILING_COMPLETE_DATASET

  type(profile_t), save, public :: &
       C_PROFILING_KINETIC,        &
       C_PROFILING_VLPSI,          &
       C_PROFILING_VNLPSI,         &
       C_PROFILING_XC,             &
       C_PROFILING_XC_LOCAL,       &
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
  ! Create profiling subdirectory.
  subroutine profiling_init()
    
    character(len=4) :: filenum
    character(len=4) :: dirnum
    integer :: ii

    !%Variable ProfilingMode
    !%Default no
    !%Type integer
    !%Section Execution::Optimization
    !%Description
    !% Use this variable to run octopus in profiling mode. In this mode
    !% octopus records time spent in certain areas of the code and
    !% the number of times this code is executed. These numbers
    !% are written in './profiling.NNN/profiling.nnn' with nnn being the
    !% node number (000 in serial) and NNN the number of processors.
    !% This is mainly for development purposes. Note, however, that
    !% octopus should be compiled with --disable-debug to do proper
    !% profiling.
    !%Option no 0
    !% No profiling information is generated.
    !%Option prof_time 1
    !% Profile the time spent in defined profiling regions.
    !%Option prof_memory 2
    !% Additionally to time, memory usage is reported.
    !%Option prof_memory_full 4
    !% Additionally to time, memory usage is reported.
    !%End

    call loct_parse_int('ProfilingMode', 0, prof_vars%mode)
    if(.not.varinfo_valid_option('ProfilingMode', prof_vars%mode, is_flag=.true.)) then
      call input_error('ProfilingMode')
    end if

    in_profiling_mode = (prof_vars%mode > 0)
    if(.not.in_profiling_mode) return

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      prof_vars%mode = ior(prof_vars%mode, PROFILING_MEMORY)
    end if

    call push_sub('profiling.profiling_init')

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
    end if

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      filenum = '0000'
      dirnum  = 'ser '
#if defined(HAVE_MPI)
      if(mpi_world%size > 1) then
        write(filenum, '(i4.4)') mpi_world%rank
        write(dirnum,  '(i4.4)') mpi_world%size
      end if
#endif
      
      call io_mkdir('profiling.'//trim(dirnum))
      prof_vars%mem_iunit = io_open('profiling.'//trim(dirnum)//'/memory.'//trim(filenum), action='write')
    end if

    ! initialize time profiling
    prof_vars%last_profile = 0
    nullify(prof_vars%current%p)

    call init_profiles()

    call pop_sub()

  contains
    subroutine init_profiles
      call profile_init(C_PROFILING_COMPLETE_DATASET, 'COMPLETE_DATASET')
      call profile_init(C_PROFILING_KINETIC,          'KINETIC')
      call profile_init(C_PROFILING_VLPSI,            'VLPSI')
      call profile_init(C_PROFILING_VNLPSI,           'VNLPSI')
      call profile_init(C_PROFILING_XC,               'XC')
      call profile_init(C_PROFILING_XC_LOCAL,         'XC_LOCAL')
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
    end subroutine init_profiles

  end subroutine profiling_init

  subroutine profiling_end
    integer :: ii
    real(8), parameter :: megabyte = 1048576.0_8

    if(.not. in_profiling_mode) return

    do ii = 1, prof_vars%last_profile
      call profile_end(prof_vars%profile_list(ii)%p)
    end do

    if(iand(prof_vars%mode, PROFILING_MEMORY).ne.0) then
      call messages_print_stress(stdout, "Memory profiling information")
      write(message(1), '(a,i10)') 'Number of   allocations = ', prof_vars%alloc_count
      write(message(2), '(a,i10)') 'Number of deallocations = ', prof_vars%dealloc_count
      write(message(3), '(a,f18.3,a)') 'Maximum total memory allocated = ', prof_vars%max_memory/megabyte, ' Mbytes'
      write(message(4), '(2x,a,a)') 'at ', trim(prof_vars%max_memory_location)
      call write_info(4)

      message(1) = ''
      message(2) = 'Largest variables allocated:'
      call write_info(2)
      do ii = 1, MAX_MEMORY_VARS
        write(message(1),'(i2,f18.3,2a)') ii, prof_vars%large_vars_size(ii)/megabyte, ' Mbytes ', trim(prof_vars%large_vars(ii))
        call write_info(1)
      end do

      call messages_print_stress(stdout)

      if(prof_vars%alloc_count.ne.prof_vars%dealloc_count) then
        message(1) = "Not all memory was deallocated!";
        call write_warning(1)
      end if
    end if

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then
      call io_close(prof_vars%mem_iunit)
    end if

  end subroutine profiling_end


  ! ---------------------------------------------------------
  ! Initialize a profile object and add it to the list
  subroutine profile_init(this, label)
    type(profile_t), target, intent(out)   :: this
    character(*),            intent(in)    :: label 
    
    this%label = label
    this%total_time = M_ZERO
    this%self_time  = M_ZERO
    this%entry_time = HUGE(this%entry_time)
    this%count  = 0
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
    this%initialized = .false.
  end subroutine profile_end


  ! ---------------------------------------------------------
  logical function profile_is_initialized(this)
    type(profile_t), intent(in)   :: this
    profile_is_initialized = this%initialized
  end function profile_is_initialized


  ! ---------------------------------------------------------
  ! Increment in counter and save entry time.
  subroutine profiling_in(this, label)
    type(profile_t), target, intent(inout) :: this
    character(*), optional,  intent(in)    :: label 

    real(8) :: now

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

    prof_vars%current%p => this
    this%entry_time = now
    
  end subroutine profiling_in


  ! ---------------------------------------------------------
  ! Increment out counter and sum up difference between entry
  ! and exit time.
  subroutine profiling_out(this)
    type(profile_t),   intent(inout) :: this

    real(8) :: now, time_spent
    
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
    
    if(associated(this%parent)) then 
      !remove the spent from the self time of our parent
      this%parent%self_time = this%parent%self_time - time_spent
      !and set parent as current
      prof_vars%current%p => this%parent
    else
      nullify(prof_vars%current%p)
    end if

  end subroutine profiling_out


  ! ---------------------------------------------------------
  subroutine iprofiling_count_operations(ops)
    integer,         intent(in)    :: ops

    if(.not.in_profiling_mode) return

    prof_vars%current%p%op_count = prof_vars%current%p%op_count + dble(ops)
  end subroutine iprofiling_count_operations


  ! ---------------------------------------------------------
  subroutine dprofiling_count_operations(ops)
    real(8),         intent(in)    :: ops

    if(.not.in_profiling_mode) return

    prof_vars%current%p%op_count = prof_vars%current%p%op_count + ops
  end subroutine dprofiling_count_operations


  ! ---------------------------------------------------------
  subroutine profiling_count_tran_int(trf, type)
    integer,         intent(in)    :: trf
    integer,         intent(in)    :: type

    if(.not.in_profiling_mode) return

    prof_vars%current%p%tr_count = prof_vars%current%p%tr_count + dble(4*trf)
  end subroutine profiling_count_tran_int


  ! ---------------------------------------------------------
  subroutine profiling_count_tran_real_4(trf, type)
    integer,         intent(in)    :: trf
    real(4),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count = prof_vars%current%p%tr_count + dble(4*trf)
  end subroutine profiling_count_tran_real_4


  ! ---------------------------------------------------------
  subroutine profiling_count_tran_real_8(trf, type)
    integer,         intent(in)    :: trf
    real(8),         intent(in)    :: type
    
    if(.not.in_profiling_mode) return
    
    prof_vars%current%p%tr_count = prof_vars%current%p%tr_count + dble(8*trf)
  end subroutine profiling_count_tran_real_8


  ! ---------------------------------------------------------
  subroutine profiling_count_tran_complex_4(trf, type)
    integer,         intent(in)    :: trf
    complex(4),      intent(in)    :: type

    if(.not.in_profiling_mode) return

    prof_vars%current%p%tr_count = prof_vars%current%p%tr_count + dble(8*trf)
  end subroutine profiling_count_tran_complex_4


  ! ---------------------------------------------------------
  subroutine profiling_count_tran_complex_8(trf, type)
    integer,         intent(in)    :: trf
    complex(8),      intent(in)    :: type

    if(.not.in_profiling_mode) return

    prof_vars%current%p%tr_count = prof_vars%current%p%tr_count + dble(16*trf)
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
    profile_throughput = this%op_count / this%self_time / CNST(1.0e6)
  end function profile_throughput


  ! ---------------------------------------------------------
  real(8) function profile_bandwidth(this)
    type(profile_t), intent(in) :: this
    profile_bandwidth = this%tr_count / (this%self_time*CNST(1024.0)**2)
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
  ! Write profiling results of each node to profiling.NNN/profifling.nnn
  ! The format of each line is
  ! tag-label    pass_in    pass_out    time   time/pass_in
  !
  ! The last column gives the average time consumed between in and out
  ! (only, if pass_in and pass_out are equal).
  subroutine profiling_output
    integer          :: ii
    integer          :: iunit
    character(len=4) :: filenum
    character(len=4) :: dirnum
    real(8)          :: total_time
    type(profile_t), pointer :: prof

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_output')

    filenum = '0000'
    dirnum  = 'ser '
#if defined(HAVE_MPI)
    if(mpi_world%size > 1) then
      write(filenum, '(i4.4)') mpi_world%rank
      write(dirnum, '(i4.4)') mpi_world%size
    end if
#endif

    if(mpi_grp_is_root(mpi_world)) call io_mkdir('profiling.'//trim(dirnum))
#ifdef HAVE_MPI
    call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

    iunit = io_open('profiling.'//trim(dirnum)//'/profiling.'//filenum, action='write')
    if(iunit.lt.0) then
      message(1) = 'Could not write profiling results.'
      call write_warning(1)
      call pop_sub(); return
    end if

    write(iunit, '(2a)')                                                                   &
      '                                                      ACCUMULATIVE TIME         ', &
      '|                SELF TIME'
    write(iunit, '(2a)')                                                                    &
      '                                            ------------------------------------|', &
      '------------------------------------------------------------'
    write(iunit, '(2a)')                                                                    &
      'TAG                   NUMBER_OF_CALLS       TOTAL_TIME    TIME_PER_CALL   %TIME |', &
      '        TOTAL_TIME    TIME_PER_CALL  MFLOPS MBYTES/S   %TIME'
    write(iunit, '(2a)')                                                                    &
      '================================================================================|', &
      '============================================================'

    total_time = profile_total_time(C_PROFILING_COMPLETE_DATASET)

    do ii = 1, prof_vars%last_profile
      prof =>  prof_vars%profile_list(ii)%p

      if(profile_num_calls(prof) == 0) cycle

      write(iunit, '(a,i20,2f17.7,f8.1,a,2f17.7,f8.1,f9.1,f8.1)')     &
           profile_label(prof),                             & 
           profile_num_calls(prof),                         &
           profile_total_time(prof),                        &
           profile_total_time_per_call(prof),               &
           profile_total_time(prof)/total_time*CNST(100.0), &
           ' | ',                                           &
           profile_self_time(prof),                         &
           profile_self_time_per_call(prof),                &
           profile_throughput(prof),                        &
           profile_bandwidth(prof),                         &
           profile_self_time(prof)/total_time*CNST(100.0)
    end do

    call io_close(iunit)

    call pop_sub()
  end subroutine profiling_output


  ! ---------------------------------------------------------
  subroutine profiling_make_position_str(var, file, line, str)
    character(len=*), intent(in)  :: var
    character(len=*), intent(in)  :: file
    integer,          intent(in)  :: line
    character(len=*), intent(out) :: str

    integer            :: ii, jj, nn
    
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
        call write_fatal(1)
      end if
    end if

    write(str, '(4a,i5,a)') var(1:jj), "(", trim(file), ":", line, ")"
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
    
    call profiling_make_position_str(var, file, line, str)

    ! get number of pages
    mem = get_memory_usage()

    write(prof_vars%mem_iunit, '(f16.6,1x,a,3i16,1x,a)') loct_clock() - prof_vars%start_time, &
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

    size = size_ ! make a copy that we can change

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then 
      call profiling_memory_log('A ', var, file, line, size)
    end if

    prof_vars%alloc_count  = prof_vars%alloc_count + 1
    prof_vars%total_memory = prof_vars%total_memory + size
    if(prof_vars%total_memory > prof_vars%max_memory) then
      prof_vars%max_memory = prof_vars%total_memory
      call profiling_make_position_str(var, file, line, prof_vars%max_memory_location)
    end if

    call profiling_make_position_str(var, file, line, str)

    ! check if variable is already in stack
    do ii = 1, MAX_MEMORY_VARS
      if(trim(str).eq.trim(prof_vars%large_vars(ii))) then
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
        prof_vars%large_vars(ii) = trim(str)
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
    
    prof_vars%dealloc_count  = prof_vars%dealloc_count + 1
    prof_vars%total_memory   = prof_vars%total_memory - size

    if(iand(prof_vars%mode, PROFILING_MEMORY_FULL).ne.0) then 
      call profiling_memory_log('D ', var, file, line, -size)
    end if

  end subroutine profiling_memory_deallocate
 

end module profiling_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
