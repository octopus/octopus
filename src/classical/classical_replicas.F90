!! Copyright (C) 2020 Heiko Appel, Kevin Lively
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

module classical_replicas_oct_m
  use classical_particle_oct_m
  use global_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use quickrnd_oct_m
  use replicas_oct_m
  use system_oct_m
  use system_factory_abst_oct_m
  implicit none

  private
  public ::               &
    classical_replicas_t

  type, extends(replicas_t) :: classical_replicas_t

  contains
    procedure :: initial_conditions => classical_replicas_initial_conditions

    final :: classical_replicas_finalizer
  end type classical_replicas_t

  interface classical_replicas_t
    procedure classical_replicas_constructor
  end interface classical_replicas_t

contains

  ! ---------------------------------------------------------------------------------------
  recursive function classical_replicas_constructor(namespace, factory) result(system)
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory
    class(classical_replicas_t),   pointer    :: system

    PUSH_SUB(classical_replicas_constructor)

    SAFE_ALLOCATE(system)

    call replicas_init(system, namespace, factory)

    POP_SUB(classical_replicas_constructor)
  end function classical_replicas_constructor

  ! ---------------------------------------------------------------------------------------
  recursive subroutine classical_replicas_init(this, namespace, factory)
    class(replicas_t),      intent(inout) :: this
    type(namespace_t),            intent(in) :: namespace
    class(system_factory_abst_t), intent(in) :: factory

    PUSH_SUB(classical_replicas_init)

    call multisystem_init(this, namespace, factory)

    POP_SUB(classical_replicas_init)
  end subroutine classical_replicas_init

  ! ---------------------------------------------------------
  subroutine classical_replicas_initial_conditions(this, from_scratch)
    class(classical_replicas_t), intent(inout) :: this
    logical,                     intent(in)    :: from_scratch
          
    class(system_t), pointer :: centroid_particle, classical_particle
    type(system_iterator_t) :: iter
    integer :: idir
    integer :: rand_gen_seed = 0
    integer :: seed, str_len, str_index
    character (len=128) :: namespace
    FLOAT :: rand_num

    PUSH_SUB(classical_replicas_initial_conditions)

    ! the first system in the list has the basename without replica number
    ! for this system we read the initial conditions from the input file
    call iter%start(this%list)
    centroid_particle => iter%get_next()
    call centroid_particle%initial_conditions(from_scratch) 

    do while (iter%has_next())
      classical_particle => iter%get_next()
      rand_gen_seed = rand_gen_seed + 1
      seed = rand_gen_seed * 10000

      select case(this%replica_distribution)
      case(UNIFORM_REPLICA)
        do idir = 1, classical_particle%space%dim
          call quickrnd(seed, rand_num)
!          write(*,*) centroid_particle%get_pos()
!          classical_particle%pos(idir) = centroid_particle%pos(idir) + &
!                  this%distribution_width*centroid_particle%pos(idir)*(rand_num - M_HALF)
        end do
      case (GAUSS_REPLICA)
        message(1) = "Gaussian replica distribution is not implemented so far"
        call messages_fatal(1)
      case (INPUT_REPLICA)
        message(1) = "Replica distributions from input files are not implemented so far"
        call messages_fatal(1)
      end select
    end do

    POP_SUB(classical_replicas_initial_conditions)
  end subroutine classical_replicas_initial_conditions
 

  ! ---------------------------------------------------------
  recursive subroutine classical_replicas_finalizer(this)
    type(classical_replicas_t), intent(inout) :: this

    PUSH_SUB(classical_replicas_finalizer)

    call multisystem_end(this)

    POP_SUB(classical_replicas_finalizer)
  end subroutine classical_replicas_finalizer

end module classical_replicas_oct_m
