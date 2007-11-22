!! Copyright (C) 2005-2006 Heiko Appel, Florian Lorenzen
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
  use messages_m
  use mpi_m

  implicit none
  private

  public ::                             &
    profile_t,                          &
    profile_init,                       &
    profile_is_initialized,             &
    profiling_init,                     &
    profiling_end,                      &
    profiling_in,                       &
    profiling_out,                      &
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
  !and substracts the time spent from the parent`s self time.
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
    integer                  :: count
    type(profile_t), pointer :: parent
    logical                  :: initialized = .false.
    logical                  :: active
  end type profile_t

  type profile_pointer_t
    private
    type(profile_t), pointer :: p
  end type profile_pointer_t
  
  type(profile_pointer_t)  :: current !the currently active profile
  type(profile_pointer_t)  :: profile_list(MAX_PROFILES) !the list of all profilers
  integer                  :: last_profile
  
  !For the moment we will have the profiler object here, but they
  !should be moved to their respective modules.
  !i.e. DO NOT PUT NEW PROFILES HERE

  type(profile_t), save, public ::    &
       C_PROFILING_COMPLETE_DATASET,  &
       C_PROFILING_MF_INTEGRATE,      &
       C_PROFILING_MF_DOTP,           &
       C_PROFILING_MF_NRM2,           &
       C_PROFILING_GHOST_UPDATE,      &
       C_PROFILING_VEC_INTEGRATE,     &
       C_PROFILING_SCF_CYCLE,         &
       C_PROFILING_MF_DOTP_ALLREDUCE, &
       C_PROFILING_HPSI

  type(profile_t), save, public :: &
       C_PROFILING_KINETIC,        &
       C_PROFILING_VLPSI,          &
       C_PROFILING_VNLPSI,         &
       C_PROFILING_POISSON_SOLVE,  &
       C_PROFILING_XC,             &
       C_PROFILING_XC_LOCAL,       &
       C_PROFILING_XC_OEP,         &
       C_PROFILING_XC_EXX,         &
       C_PROFILING_XC_SIC,         &
       C_PROFILING_XC_KLI

  type(profile_t), save, public :: &
       C_PROFILING_XC_OEP_FULL,    &
       C_PROFILING_TIME_STEP,      &
       C_PROFILING_GRAM_SCHMIDT1,  &
       C_PROFILING_GRAM_SCHMIDT2,  &
       C_PROFILING_EIGEN_SOLVER,   &
       C_PROFILING_DISK_ACCESS,    &
       C_PROFILING_LCAO,           &
       C_PROFILING_LCAO_INIT

  type(profile_t), save, public ::    &
       C_PROFILING_ELEC_PROPAGATOR,   &
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
  subroutine profiling_init
#if defined(HAVE_MPI)
    character(len=3) :: dirnum
#endif

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_init')

#if defined(HAVE_MPI)
    if(mpi_grp_is_root(mpi_world)) then
      write(dirnum, '(i3.3)') mpi_world%size
      call io_mkdir(trim('profiling.'//dirnum))
    end if
#else
    call io_mkdir('profiling.ser')
#endif
    
    last_profile = 0

    nullify(current%p)

    call init_profiles

    call pop_sub()

  contains
    subroutine init_profiles
      call profile_init(C_PROFILING_COMPLETE_DATASET, 'COMPLETE_DATASET')
      call profile_init(C_PROFILING_MF_INTEGRATE,     'MF_INTEGRATE')
      call profile_init(C_PROFILING_MF_DOTP,          'MF_DOTP')
      call profile_init(C_PROFILING_MF_NRM2,          'MF_NRM2')
      call profile_init(C_PROFILING_GHOST_UPDATE,     'GHOST_UPDATE')
      call profile_init(C_PROFILING_VEC_INTEGRATE,    'VEC_INTEGRATE')
      call profile_init(C_PROFILING_SCF_CYCLE,        'SCF_CYCLE')
      call profile_init(C_PROFILING_MF_DOTP_ALLREDUCE,'MF_DOTP_ALLREDUCE')
      call profile_init(C_PROFILING_HPSI,             'HPSI')
      call profile_init(C_PROFILING_KINETIC,          'KINETIC')
      call profile_init(C_PROFILING_VLPSI,            'VLPSI')
      call profile_init(C_PROFILING_VNLPSI,           'VNLPSI')
      call profile_init(C_PROFILING_POISSON_SOLVE,    'POISSON_SOLVE')
      call profile_init(C_PROFILING_XC,               'XC')
      call profile_init(C_PROFILING_XC_LOCAL,         'XC_LOCAL')
      call profile_init(C_PROFILING_XC_OEP,           'XC_OEP')
      call profile_init(C_PROFILING_XC_EXX,           'XC_EXX')
      call profile_init(C_PROFILING_XC_SIC,           'XC_SIC')
      call profile_init(C_PROFILING_XC_KLI,           'XC_KLI')
      call profile_init(C_PROFILING_XC_OEP_FULL,      'XC_OEP_FULL')
      call profile_init(C_PROFILING_TIME_STEP	,     'TIME_STEP')
      call profile_init(C_PROFILING_GRAM_SCHMIDT1,    'GRAM_SCHMIDT1')
      call profile_init(C_PROFILING_GRAM_SCHMIDT2,    'GRAM_SCHMIDT2')
      call profile_init(C_PROFILING_EIGEN_SOLVER,     'EIGEN_SOLVER')
      call profile_init(C_PROFILING_DISK_ACCESS,      'DISK_ACCESS')
      call profile_init(C_PROFILING_LCAO,             'LCAO')
      call profile_init(C_PROFILING_LCAO_INIT,        'LCAO_INIT')
      call profile_init(C_PROFILING_ELEC_PROPAGATOR,  'ELEC_PROPAGATOR')
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
    do ii = 1, last_profile
      call profile_end(profile_list(ii)%p)
    end do
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
    this%active = .false.
    nullify(this%parent)

    if(.not. in_profiling_mode) return

    last_profile = last_profile + 1

    ASSERT(last_profile .le. MAX_PROFILES)

    profile_list(last_profile)%p => this
    
    this%initialized = .true.
    
  end subroutine profile_init


  subroutine profile_end(this)
    type(profile_t), intent(inout) :: this
    this%initialized = .false.
  end subroutine profile_end

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
    if(associated(current%p)) then
      !keep a pointer to the parent
      this%parent => current%p
    else 
      !we are orphans
      nullify(this%parent)
    end if

    current%p => this
    this%entry_time = now
    
  end subroutine profiling_in


  ! ---------------------------------------------------------
  ! Increment out counter and sum up difference between entry
  ! and exit time.
  subroutine profiling_out(this)
    type(profile_t), intent(inout)   :: this
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
      current%p => this%parent
    else
      nullify(current%p)
    end if

  end subroutine profiling_out

  real(8) function profile_total_time(this)
    type(profile_t), intent(in) :: this
    profile_total_time = this%total_time
  end function profile_total_time

  real(8) function profile_self_time(this)
    type(profile_t), intent(in) :: this
    profile_self_time = this%self_time
  end function profile_self_time

  real(8) function profile_total_time_per_call(this)
    type(profile_t), intent(in) :: this
    profile_total_time_per_call = this%total_time / dble(this%count)
  end function profile_total_time_per_call

  real(8) function profile_self_time_per_call(this)
    type(profile_t), intent(in) :: this
    profile_self_time_per_call = this%self_time / dble(this%count)
  end function profile_self_time_per_call

  integer function profile_num_calls(this)
    type(profile_t), intent(in) :: this
    
    profile_num_calls = this%count
  end function profile_num_calls

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
    character(len=3) :: filenum
    character(len=3) :: dirnum
    real(8)          :: total_time
    type(profile_t), pointer :: prof

    if(.not.in_profiling_mode) return

    call push_sub('profiling.profiling_output')

#if defined(HAVE_MPI)
    write(filenum, '(i3.3)') mpi_world%rank
    write(dirnum, '(i3.3)') mpi_world%size
#else
    filenum = '000'
    dirnum  = 'ser'
#endif

    total_time = profile_total_time(C_PROFILING_COMPLETE_DATASET)

    iunit = io_open('profiling.'//dirnum//'/profiling.'//filenum, action='write')
    if(iunit.lt.0) then
      message(1) = 'Could not write profiling results.'
      call write_warning(1)
      call pop_sub()
      return
    end if
    write(iunit, '(2a)')                                                                   &
      '                                                      ACCUMULATIVE TIME          ', &
      '|                SELF TIME'
    write(iunit, '(2a)')                                                                    &
      '                                            -------------------------------------|', &
      '---------------------------------------------'
    write(iunit, '(2a)')                                                                    &
      'TAG                   NUMBER_OF_CALLS       TOTAL_TIME    TIME_PER_CALL    %TIME |', &
      '        TOTAL_TIME    TIME_PER_CALL    %TIME'
    write(iunit, '(2a)')                                                                    &
      '=================================================================================|', &
      '============================================='
    do ii = 1, last_profile
      prof => profile_list(ii)%p

      if(profile_num_calls(prof) == 0) cycle

      write(iunit, '(a,i20,2f17.7,f9.1,a,2f17.7,f9.1)')     &
           profile_label(prof),                             & 
           profile_num_calls(prof),                         &
           profile_total_time(prof),                        &
           profile_total_time_per_call(prof),               &
           profile_total_time(prof)/total_time*CNST(100.0), &
           ' | ',                                           &
           profile_self_time(prof),                         &
           profile_self_time_per_call(prof),                &
           profile_self_time(prof)/total_time*CNST(100.0)
    end do

    call io_close(iunit)

    call pop_sub()
  end subroutine profiling_output

end module profiling_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
