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
  use parser_m
  use messages_m
#ifdef HAVE_OPENCL
  use opencl_m
#endif
  use profiling_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::                  &
    batch_state_t,           &
    batch_state_l_t,         &
    batch_t,                 &
    batch_init,              &
    batch_copy,              &
    batch_end,               &
    batch_add_state,         &
    dbatch_new,              &
    zbatch_new,              &
    dbatch_delete,           &
    zbatch_delete,           &
    batch_set,               &
#ifdef HAVE_OPENCL
    batch_to_opencl_buffer,  &
#endif
    batch_is_ok

  !--------------------------------------------------------------
  type batch_state_t
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    integer        :: ist
  end type batch_state_t

  type batch_state_l_t
    FLOAT, pointer :: dpsi(:)
    CMPLX, pointer :: zpsi(:)
  end type batch_state_l_t
  
  type batch_t
    type(batch_state_t), pointer   :: states(:)
    integer                        :: nst
    integer                        :: current
    integer                        :: dim

    !> We also need a linear array with the states in order to calculate derivatives, etc.
    integer                        :: nst_linear
    type(batch_state_l_t), pointer :: states_linear(:)

    !> If the memory is contiguous, we can perform some operations faster.
    FLOAT,               pointer   :: dpsicont(:, :, :)
    CMPLX,               pointer   :: zpsicont(:, :, :)
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure  batch_init_empty
    module procedure  batch_init_empty_linear
    module procedure dbatch_init_contiguous
    module procedure zbatch_init_contiguous
  end interface

  interface batch_add_state
    module procedure dbatch_add_state
    module procedure zbatch_add_state
    module procedure dbatch_add_state_linear
    module procedure zbatch_add_state_linear
  end interface

  interface batch_set
    module procedure dbatch_set
    module procedure zbatch_set
  end interface

contains

  !--------------------------------------------------------------
  subroutine batch_end(this)
    type(batch_t), intent(inout) :: this

    call push_sub('batch.batch_end')

    nullify(this%dpsicont)
    nullify(this%zpsicont)

    SAFE_DEALLOCATE_P(this%states)
    SAFE_DEALLOCATE_P(this%states_linear)

    call pop_sub('batch.batch_end')
  end subroutine batch_end


  !--------------------------------------------------------------
  subroutine batch_init_empty (this, dim, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nst
    
    integer :: ist

    call push_sub('batch.batch_init_empty')
    
    this%nst = nst
    this%dim = dim
    this%current = 1
    nullify(this%dpsicont, this%zpsicont)
    
    SAFE_ALLOCATE(this%states(1:nst))
    do ist = 1, nst
      nullify(this%states(ist)%dpsi)
      nullify(this%states(ist)%zpsi)
    end do
    
    this%nst_linear = nst*dim
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
    end do
    
    call pop_sub('batch.batch_init_empty')
    
  end subroutine batch_init_empty


  !--------------------------------------------------------------
  !> When we are interested in batches of 1D functions
  subroutine batch_init_empty_linear(this, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: nst
    
    integer :: ist

    call push_sub('batch.batch_init_empty_linear')
    
    this%nst = 0
    this%dim = 0
    this%current = 1
    nullify(this%dpsicont, this%zpsicont)
    nullify(this%states)

    this%nst_linear = nst
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
    end do
    
    call pop_sub('batch.batch_init_empty_linear')
    
  end subroutine batch_init_empty_linear


  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    type(batch_t), intent(in)   :: this

    integer :: ist
    
    call push_sub('batch.batch_is_ok')

    ok = (this%nst >= 1)
    if(ok) then
      do ist = 1, this%nst_linear
        ok = ok.and. &
          (associated(this%states_linear(ist)%dpsi) .or. associated(this%states_linear(ist)%zpsi))
      end do
    end if

    call pop_sub('batch.batch_is_ok')
  end function batch_is_ok


  !--------------------------------------------------------------
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

    do ii = 1, bout%nst_linear
      if(associated(bin%states_linear(ii)%dpsi)) bout%states_linear(ii)%dpsi => bin%states_linear(ii)%dpsi
      if(associated(bin%states_linear(ii)%zpsi)) bout%states_linear(ii)%zpsi => bin%states_linear(ii)%zpsi
    end do

    call pop_sub('batch.batch_copy')

  end subroutine batch_copy

  integer pure function batch_type(this) result(btype)
    type(batch_t),      intent(in)    :: this

    if(associated(this%states_linear(1)%dpsi)) btype = TYPE_FLOAT
    if(associated(this%states_linear(1)%zpsi)) btype = TYPE_CMPLX
  end function batch_type

  ! ----------------------------------------------------
#ifdef HAVE_OPENCL
  subroutine batch_to_opencl_buffer(this, opencl, np, flags, buffer)
    type(batch_t),      intent(in)    :: this
    type(opencl_t),     intent(inout) :: opencl
    integer,            intent(in)    :: np
    integer,            intent(in)    :: flags
    type(opencl_mem_t), intent(out)   :: buffer

    integer(SIZEOF_SIZE_T) :: size, pnp
    integer :: ist

    pnp = opencl_padded_size(opencl, np, batch_type(this))
    size = pnp*this%nst_linear

    ASSERT(batch_is_ok(this))

    call opencl_create_buffer(buffer, opencl, flags, batch_type(this), size)

    do ist = 1, this%nst_linear
      if(batch_type(this) == TYPE_FLOAT) then
        call opencl_write_buffer(buffer, opencl, np, this%states_linear(ist)%dpsi, offset = (ist - 1)*pnp)
      else
        call opencl_write_buffer(buffer, opencl, np, this%states_linear(ist)%zpsi, offset = (ist - 1)*pnp)
      end if
    end do

  end subroutine batch_to_opencl_buffer
#endif

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
