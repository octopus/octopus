!! Copyright (C) 2019 M. Oliveira
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

module wfs_elec_oct_m
  use batch_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::                         &
    wfs_elec_t,                     &
    wfs_elec_init

  type, extends(batch_t) :: wfs_elec_t
    private
    integer, public :: ik
    logical, public :: has_phase
  contains
    procedure :: clone_to => wfs_elec_clone_to
    procedure :: clone_to_array => wfs_elec_clone_to_array
    procedure :: copy_to => wfs_elec_copy_to
    procedure :: check_compatibility_with => wfs_elec_check_compatibility_with
    procedure :: end => wfs_elec_end
  end type wfs_elec_t

  !--------------------------------------------------------------
  interface wfs_elec_init
    module procedure  wfs_elec_init_empty
    module procedure dwfs_elec_init_contiguous
    module procedure zwfs_elec_init_contiguous
    module procedure dwfs_elec_init_contiguous_2d
    module procedure zwfs_elec_init_contiguous_2d
  end interface wfs_elec_init

contains

  !--------------------------------------------------------------
  subroutine wfs_elec_init_empty(this, dim, nst, ik)
    type(wfs_elec_t), intent(out)   :: this
    integer,          intent(in)    :: dim
    integer,          intent(in)    :: nst
    integer,          intent(in)    :: ik

    PUSH_SUB(wfs_elec_init_empty)

    this%ik = ik
    this%has_phase = .false.
    call batch_init(this%batch_t, dim, nst)

    POP_SUB(wfs_elec_init_empty)
  end subroutine wfs_elec_init_empty

  !--------------------------------------------------------------
  subroutine wfs_elec_clone_to(this, dest, pack, copy_data)
    class(wfs_elec_t),           intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest
    logical,        optional,    intent(in)    :: pack
    logical,        optional,    intent(in)    :: copy_data

    PUSH_SUB(wfs_elec_clone_to)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE(wfs_elec_t, dest)
    else
      message(1) = "Internal error: destination batch in wfs_elec_clone_to has been previously allocated."
      call messages_fatal(1)
    end if

    select type (dest)
    class is (wfs_elec_t)
      call this%copy_to(dest, pack, copy_data)
    class default
      message(1) = "Internal error: imcompatible batches in wfs_elec_clone_to."
      call messages_fatal(1)
    end select

    POP_SUB(wfs_elec_clone_to)
  end subroutine wfs_elec_clone_to

  !--------------------------------------------------------------
  subroutine wfs_elec_clone_to_array(this, dest, n_batches, pack, copy_data)
    class(wfs_elec_t),           intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest(:)
    integer,                     intent(in)    :: n_batches
    logical,        optional,    intent(in)    :: pack
    logical,        optional,    intent(in)    :: copy_data

    integer :: ib

    PUSH_SUB(wfs_elec_clone_to_array)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE_ARRAY(wfs_elec_t, dest, (1:n_batches))
    else
      message(1) = "Internal error: destination batch in wfs_elec_clone_to_array has been previously allocated."
      call messages_fatal(1)
    end if

    select type (dest)
    class is (wfs_elec_t)
      do ib = 1, n_batches
        call this%copy_to(dest(ib), pack, copy_data)
      end do
    class default
      message(1) = "Internal error: imcompatible batches in wfs_elec_clone_to_array."
      call messages_fatal(1) 
    end select

    POP_SUB(wfs_elec_clone_to_array)
  end subroutine wfs_elec_clone_to_array

  !--------------------------------------------------------------
  subroutine wfs_elec_copy_to(this, dest, pack, copy_data)
    class(wfs_elec_t),  intent(in)    :: this
    class(batch_t),     intent(out)   :: dest
    logical,  optional, intent(in)    :: pack
    logical, optional,  intent(in)    :: copy_data

    PUSH_SUB(wfs_elec_copy_to)

    select type (dest)
    class is (wfs_elec_t)
      dest%ik = this%ik
      dest%has_phase = this%has_phase
      call this%batch_t%copy_to(dest%batch_t, pack, copy_data)
    class default
      message(1) = "Internal error: imcompatible batches in wfs_elec_copy_to."
      call messages_fatal(1)
    end select

    POP_SUB(wfs_elec_copy_to)
  end subroutine wfs_elec_copy_to

  !--------------------------------------------------------------
  subroutine wfs_elec_check_compatibility_with(this, target, only_check_dim)
    class(wfs_elec_t),  intent(in) :: this
    class(batch_t),     intent(in) :: target
    logical,  optional, intent(in) :: only_check_dim

    PUSH_SUB(wfs_elec_check_compatibility_with)

    select type (target)
    class is (wfs_elec_t)
      ASSERT(this%ik == target%ik)
      ASSERT(this%has_phase .eqv. target%has_phase)
    class default
      message(1) = "Internal error: imcompatible batches in wfs_elec_check_compatibility_with."
      call messages_fatal(1)
    end select
    call this%batch_t%check_compatibility_with(target, only_check_dim)

    POP_SUB(wfs_elec_check_compatibility_with)
  end subroutine wfs_elec_check_compatibility_with

  !--------------------------------------------------------------
  subroutine wfs_elec_end(this, copy)
    class(wfs_elec_t),       intent(inout) :: this
    logical,       optional, intent(in)    :: copy

    PUSH_SUB(wfs_elec_end)

    this%ik = -1
    this%has_phase = .false.
    call this%batch_t%end(copy)

    POP_SUB(wfs_elec_end)
  end subroutine wfs_elec_end


#include "real.F90"
#include "wfs_elec_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "wfs_elec_inc.F90"
#include "undef.F90"

end module wfs_elec_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
