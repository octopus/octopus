!! Copyright (C) 2008 X. Andrade
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
!! $Id: batch.F90 4298 2008-06-18 15:03:00Z dstrubbe $

#include "global.h"

module batch_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use messages_m
  use profiling_m
  use varinfo_m

  implicit none

  private
  public ::                     &
       state_t,                 &
       batch_t,                 &
       batch_init,              &
       batch_copy,              &
       batch_end,               &
       batch_add_state,         &
       batch_is_ok

  type state_t
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    integer        :: ist
  end type state_t

  type batch_t
    type(state_t), pointer :: states(:)
    integer                :: nst
    integer                :: current
    integer                :: dim
    FLOAT, pointer :: dpsicont(:, :, :)
    CMPLX, pointer :: zpsicont(:, :, :)
  end type batch_t

  interface batch_init
    module procedure batch_init_empty, dbatch_init_contiguous, zbatch_init_contiguous
  end interface

  interface batch_add_state
    module procedure dbatch_add_state, zbatch_add_state
  end interface

contains

  elemental subroutine state_nullify(this)
    type(state_t), intent(out) :: this
    nullify(this%dpsi, this%zpsi)
  end subroutine state_nullify

  logical elemental function state_is_associated(this) result(is_associated)
    type(state_t), intent(in) :: this

    is_associated = associated(this%dpsi) .or. associated(this%zpsi)
  end function state_is_associated

  subroutine batch_end(this)
    type(batch_t), intent(inout) :: this

    nullify(this%dpsicont)
    nullify(this%zpsicont)
    SAFE_DEALLOCATE_P(this%states)
  end subroutine batch_end

  subroutine batch_init_empty(this, dim, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nst
    
    call push_sub('batch.batch_init_empty')

    this%nst = nst
    this%dim = dim
    this%current = 1
    nullify(this%dpsicont, this%zpsicont)
    
    ALLOCATE(this%states(1:nst), nst)
    call state_nullify(this%states(1:nst))

    call pop_sub()

  end subroutine batch_init_empty

  logical function batch_is_ok(this) result(ok)
    type(batch_t), intent(in)   :: this
    
    ok = (this%nst >= 1) .and. all(state_is_associated(this%states(1:this%nst)))
  end function batch_is_ok

  subroutine batch_copy(bin, bout)
    type(batch_t), intent(in)    :: bin
    type(batch_t), intent(out)   :: bout

    integer :: ii

    call push_sub('batch.batch_copy')

    call batch_init_empty(bout, bin%dim, bin%nst)
    bout%current = bin%current

    do ii = 1, bout%nst
      bout%states(ii)%ist = bin%states(ii)%ist
      if(associated(bin%states(ii)%dpsi)) bout%states(ii)%dpsi => bin%states(ii)%dpsi
      if(associated(bin%states(ii)%zpsi)) bout%states(ii)%zpsi => bin%states(ii)%zpsi
    end do

    call pop_sub()

  end subroutine batch_copy

#include "undef.F90"
#include "real.F90"
#include "batch_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_inc.F90"
#include "undef.F90"


end module batch_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
