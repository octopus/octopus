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
  use blas_m
  use c_pointer_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use parser_m
  use math_m
  use messages_m
  use opencl_m
  use profiling_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::                         &
    batch_state_t,                  &
    batch_state_l_t,                &
    batch_pack_t,                   &
    batch_t,                        &
    batch_init,                     &
    batch_copy,                     &
    batch_end,                      &
    batch_add_state,                &
    dbatch_new,                     &
    zbatch_new,                     &
    dbatch_delete,                  &
    zbatch_delete,                  &
    batch_set,                      &
    batch_set_zero,                 &
    batch_is_packed,                &
    batch_pack_was_modified,        &
    batch_is_ok,                    &
    batch_axpy,                     &
    batch_copy_data

#ifdef HAVE_OPENCL
  public ::                         &
    batch_pack,                     &
    batch_unpack
#endif


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

  type batch_pack_t
    integer                        :: size(1:2)
    integer                        :: size_real(1:2)
#ifdef HAVE_OPENCL
    type(opencl_mem_t)             :: buffer
#endif
  end type batch_pack_t
  
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

    integer                        :: in_buffer_count ! whether there is a copy in the opencl buffer
    logical                        :: dirty     ! if this is true, the buffer has different data
    type(batch_pack_t)             :: pack
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

  interface batch_axpy
    module procedure dbatch_axpy
    module procedure zbatch_axpy
  end interface

  type(profile_t), save :: axpy_prof

contains

  !--------------------------------------------------------------
  subroutine batch_end(this)
    type(batch_t), intent(inout) :: this

    call push_sub('batch.batch_end')

    ASSERT(.not. batch_is_packed(this))

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
    
    this%in_buffer_count = 0

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

    this%in_buffer_count = 0

    call pop_sub('batch.batch_init_empty_linear')
    
  end subroutine batch_init_empty_linear


  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    type(batch_t), intent(in)   :: this

    integer :: ist
    
    call push_sub('batch.batch_is_ok')
    
    ok = (this%nst_linear >= 1)
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

    bout%in_buffer_count     = bin%in_buffer_count
    bout%dirty               = bin%dirty
    bout%pack%size(1:2)      = bin%pack%size(1:2)
    bout%pack%size_real(1:2) = bin%pack%size_real(1:2)
      
#ifdef HAVE_OPENCL
    if(opencl_is_enabled()) then
      bout%pack%buffer = bin%pack%buffer
    end if
#endif    

    call pop_sub('batch.batch_copy')

  end subroutine batch_copy

  ! ----------------------------------------------------
  ! THREADSAFE
  integer pure function batch_type(this) result(btype)
    type(batch_t),      intent(in)    :: this

    if(associated(this%states_linear(1)%dpsi)) btype = TYPE_FLOAT
    if(associated(this%states_linear(1)%zpsi)) btype = TYPE_CMPLX
  end function batch_type

  ! ----------------------------------------------------

  logical pure function batch_is_packed(this) result(in_buffer)
    type(batch_t),      intent(in)    :: this

    in_buffer = this%in_buffer_count > 0
  end function batch_is_packed

  ! ----------------------------------------------------

  integer pure function batch_max_size(this) result(size)
    type(batch_t),      intent(in)    :: this

    integer :: ist

    size = 0
    do ist = 1, this%nst_linear
      if(associated(this%states_linear(ist)%dpsi)) then
        size = max(size, ubound(this%states_linear(ist)%dpsi, dim = 1))
      else
        size = max(size, ubound(this%states_linear(ist)%zpsi, dim = 1))
      end if
    end do

  end function batch_max_size

  ! ----------------------------------------------------

  subroutine batch_pack_was_modified(this)
    type(batch_t),      intent(inout) :: this

    this%dirty = .true.
  end subroutine batch_pack_was_modified

#ifdef HAVE_OPENCL

  ! ----------------------------------------------------

  subroutine batch_pack(this, copy)
    type(batch_t),      intent(inout) :: this
    logical, optional,  intent(in)    :: copy

    logical :: copy_

    call push_sub('batch.batch_pack')

    ASSERT(batch_is_ok(this))

    copy_ = .true.
    if(present(copy)) copy_ = copy

    if(.not. batch_is_packed(this)) then
      this%pack%size(1) = pad_pow2(this%nst_linear)
      this%pack%size(2) = opencl_padded_size(batch_max_size(this))

      this%pack%size_real = this%pack%size
      if(batch_type(this) == TYPE_CMPLX) this%pack%size_real(1) = 2*this%pack%size_real(1)

      call opencl_create_buffer(this%pack%buffer, CL_MEM_READ_WRITE, batch_type(this), product(this%pack%size))

      if(copy_) then
        this%dirty = .false.
        call batch_write_to_opencl_buffer(this)
      else
        this%dirty = .true.
      end if
    end if

    INCR(this%in_buffer_count, 1)

    call pop_sub('batch.batch_pack')
  end subroutine batch_pack

  ! ----------------------------------------------------

  subroutine batch_unpack(this, copy)
    type(batch_t),      intent(inout) :: this
    logical, optional,  intent(in)    :: copy

    logical :: copy_

    call push_sub('batch.batch_unpack')


    if(batch_is_packed(this)) then
      INCR(this%in_buffer_count, -1)

      if(this%in_buffer_count == 0) then
        copy_ = .true.
        if(present(copy)) copy_ = copy
        
        if(copy_ .and. this%dirty) then
          call batch_read_from_opencl_buffer(this)
        end if
        
        call opencl_release_buffer(this%pack%buffer)
      end if
    end if

    call pop_sub('batch.batch_unpack')
  end subroutine batch_unpack

  ! ----------------------------------------------------

  subroutine batch_write_to_opencl_buffer(this)
    type(batch_t),      intent(inout)  :: this

    integer :: ist
    type(opencl_mem_t) :: tmp
    type(c_ptr) :: kernel
    type(profile_t), save :: prof, prof_pack

    call push_sub('batch.batch_write_to_opencl_buffer')
    call profiling_in(prof, "BATCH_TO_BUFFER")

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call opencl_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call opencl_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if

    else
      ! we copy to a temporary array and then we re-arrange data
      call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))

      do ist = 1, this%nst_linear
        if(batch_type(this) == TYPE_FLOAT) then
          kernel = dpack
          call opencl_write_buffer(tmp, ubound(this%states_linear(ist)%dpsi, dim = 1), this%states_linear(ist)%dpsi)
        else
          kernel = zpack
          call opencl_write_buffer(tmp, ubound(this%states_linear(ist)%zpsi, dim = 1), this%states_linear(ist)%zpsi)
        end if

        call opencl_set_kernel_arg(kernel, 0, this%pack%size(1))
        call opencl_set_kernel_arg(kernel, 1, ist - 1)
        call opencl_set_kernel_arg(kernel, 2, tmp)
        call opencl_set_kernel_arg(kernel, 3, this%pack%buffer)

        call profiling_in(prof_pack, "CL_PACK")
        call opencl_kernel_run(kernel, (/this%pack%size(2)/), (/opencl_max_workgroup_size()/))
        call opencl_finish()
        call profiling_out(prof_pack)

      end do

      call opencl_release_buffer(tmp)

    end if

    call profiling_out(prof)
    call pop_sub('batch.batch_write_to_opencl_buffer')
  end subroutine batch_write_to_opencl_buffer

  ! ------------------------------------------------------------------

  subroutine batch_read_from_opencl_buffer(this)
    type(batch_t),      intent(inout) :: this

    integer :: ist
    type(opencl_mem_t) :: tmp
    type(c_ptr) :: kernel
    type(profile_t), save :: prof, prof_unpack

    call push_sub('batch.batch_read_from_opencl_buffer')
    call profiling_in(prof, "BATCH_FROM_BUFFER")

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call opencl_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call opencl_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if
    else
      ! we use a kernel to move to a temporary array and then we read
      call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))

      if(batch_type(this) == TYPE_FLOAT) then
        kernel = dunpack
      else
        kernel = zunpack
      end if

      do ist = 1, this%nst_linear
        call opencl_set_kernel_arg(kernel, 0, this%pack%size(1))
        call opencl_set_kernel_arg(kernel, 1, ist - 1)
        call opencl_set_kernel_arg(kernel, 2, this%pack%buffer)
        call opencl_set_kernel_arg(kernel, 3, tmp)

        call profiling_in(prof_unpack, "CL_UNPACK")
        call opencl_kernel_run(kernel, (/this%pack%size(2)/), (/opencl_max_workgroup_size()/))
        call opencl_finish()
        call profiling_out(prof_unpack)

        if(batch_type(this) == TYPE_FLOAT) then
          call opencl_read_buffer(tmp, ubound(this%states_linear(ist)%dpsi, dim = 1), this%states_linear(ist)%dpsi)
        else
          call opencl_read_buffer(tmp, ubound(this%states_linear(ist)%zpsi, dim = 1), this%states_linear(ist)%zpsi)
        end if
      end do

      call opencl_release_buffer(tmp)
    end if

    call profiling_out(prof)
    call pop_sub('batch.batch_read_from_opencl_buffer')
  end subroutine batch_read_from_opencl_buffer

#endif

!--------------------------------------------------------------

  subroutine batch_set_zero(this)
    type(batch_t),     intent(inout) :: this

    integer :: ist_linear
#ifdef HAVE_OPENCL
    integer :: bsize
#endif

    call push_sub('batch.batch_set_zero')

    if(batch_is_packed(this)) then

#ifdef HAVE_OPENCL
      bsize = product(this%pack%size)
      if(batch_type(this) == TYPE_CMPLX) bsize = bsize*2

      call opencl_set_kernel_arg(set_zero, 0, this%pack%buffer)
      call opencl_kernel_run(set_zero, (/bsize/), (/opencl_max_workgroup_size()/))
      call batch_pack_was_modified(this)
      call opencl_finish()
#endif
      
    else if (associated(this%dpsicont)) then
      this%dpsicont = M_ZERO
    else if (associated(this%zpsicont)) then
      this%zpsicont = M_ZERO
    else

      do ist_linear = 1, this%nst_linear
        if(associated(this%states_linear(ist_linear)%dpsi)) then
          this%states_linear(ist_linear)%dpsi = M_ZERO
        else
          this%states_linear(ist_linear)%zpsi = M_ZERO
        end if
      end do

    end if

    call pop_sub('batch.batch_set_zero')
  end subroutine batch_set_zero

! --------------------------------------------------------------

subroutine batch_copy_data(np, xx, yy)
  integer,           intent(in)    :: np
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist, localsize
  type(profile_t), save :: prof

  call push_sub('batch.batch_copy_data')
  call profiling_in(prof, "BATCH_COPY_DATA")

  ASSERT(batch_type(yy) == batch_type(xx))
  ASSERT(xx%nst_linear == yy%nst_linear)

#ifdef HAVE_OPENCL
  if(batch_is_packed(yy) .or. batch_is_packed(xx)) then
    ASSERT(batch_is_packed(xx))
    ASSERT(batch_is_packed(yy))

    call opencl_set_kernel_arg(kernel_copy, 0, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 1, log2(xx%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_copy, 2, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 3, log2(yy%pack%size_real(1)))

    localsize = opencl_max_workgroup_size()/yy%pack%size_real(1)
    call opencl_kernel_run(kernel_copy, (/yy%pack%size_real(1), pad(np, localsize)/), (/yy%pack%size_real(1), localsize/))

    call batch_pack_was_modified(yy)
    
    call opencl_finish()

  else
#endif

    !$omp parallel do private(ist)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call blas_copy(np, xx%states_linear(ist)%zpsi(1), 1, yy%states_linear(ist)%zpsi(1), 1)
      else
        call blas_copy(np, xx%states_linear(ist)%dpsi(1), 1, yy%states_linear(ist)%dpsi(1), 1)
      end if
    end do

#ifdef HAVE_OPENCL
  end if
#endif

  call profiling_out(prof)
  call pop_sub('batch.batch_copy_data')
end subroutine batch_copy_data

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
