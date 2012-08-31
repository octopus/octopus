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

module batch_ops_m
  use batch_m
  use blas_m
  use c_pointer_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use octcl_kernel_m
  use datasets_m
  use global_m
  use hardware_m
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
    batch_set,                      &
    batch_set_zero,                 &
    batch_axpy,                     &
    batch_copy_data,                &
    batch_set_state,                &
    batch_get_state,                &
    batch_get_points,               &
    batch_set_points

  interface batch_set
    module procedure dbatch_set
    module procedure zbatch_set
  end interface

  interface batch_axpy
    module procedure dbatch_axpy_const
    module procedure zbatch_axpy_const
    module procedure dbatch_axpy_vec
    module procedure zbatch_axpy_vec
  end interface

  interface batch_set_state
    module procedure dbatch_set_state1
    module procedure zbatch_set_state1
    module procedure dbatch_set_state2
    module procedure zbatch_set_state2
  end interface

  interface batch_get_state
    module procedure dbatch_get_state1
    module procedure zbatch_get_state1
    module procedure dbatch_get_state2
    module procedure zbatch_get_state2
  end interface

  interface batch_get_points
    module procedure dbatch_get_points
    module procedure zbatch_get_points
    module procedure batch_get_points_cl
  end interface batch_get_points

  interface batch_set_points
    module procedure dbatch_set_points
    module procedure zbatch_set_points
    module procedure batch_set_points_cl
  end interface batch_set_points

  type(profile_t), save :: axpy_prof, get_points_prof, set_points_prof

contains

  !--------------------------------------------------------------

  subroutine batch_set_zero(this)
    type(batch_t),     intent(inout) :: this

    integer :: ist_linear

    PUSH_SUB(batch_set_zero)

    call batch_pack_was_modified(this)

    if(batch_is_packed(this) .and. opencl_is_enabled()) then
#ifdef HAVE_OPENCL
      call opencl_set_buffer_to_zero(this%pack%buffer, batch_type(this), product(this%pack%size))
#endif
    else if(batch_is_packed(this) .and. batch_type(this) == TYPE_FLOAT) then
      this%pack%dpsi = M_ZERO
    else if(batch_is_packed(this) .and. batch_type(this) == TYPE_CMPLX) then
      this%pack%zpsi = M_ZERO
    else
      do ist_linear = 1, this%nst_linear
        if(associated(this%states_linear(ist_linear)%dpsi)) then
          this%states_linear(ist_linear)%dpsi = M_ZERO
        else
          this%states_linear(ist_linear)%zpsi = M_ZERO
        end if
      end do

    end if

    POP_SUB(batch_set_zero)
  end subroutine batch_set_zero

! --------------------------------------------------------------

subroutine batch_copy_data(np, xx, yy)
  integer,           intent(in)    :: np
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist
  type(profile_t), save :: prof
  integer :: localsize

  PUSH_SUB(batch_copy_data)
  call profiling_in(prof, "BATCH_COPY_DATA")

  ASSERT(batch_type(yy) == batch_type(xx))
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)

#ifdef HAVE_OPENCL
    call opencl_set_kernel_arg(kernel_copy, 0, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 1, log2(xx%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_copy, 2, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 3, log2(yy%pack%size_real(1)))
    
    localsize = opencl_max_workgroup_size()/yy%pack%size_real(1)
    call opencl_kernel_run(kernel_copy, (/yy%pack%size_real(1), pad(np, localsize)/), (/yy%pack%size_real(1), localsize/))
    
    call opencl_finish()
#endif

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_FLOAT) then
      call blas_copy(np*xx%pack%size(1), xx%pack%dpsi(1, 1), 1, yy%pack%dpsi(1, 1), 1)
    else
      call blas_copy(np*xx%pack%size(1), xx%pack%zpsi(1, 1), 1, yy%pack%zpsi(1, 1), 1)
    end if

  case(BATCH_NOT_PACKED)
    !$omp parallel do private(ist)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call blas_copy(np, xx%states_linear(ist)%zpsi(1), 1, yy%states_linear(ist)%zpsi(1), 1)
      else
        call blas_copy(np, xx%states_linear(ist)%dpsi(1), 1, yy%states_linear(ist)%dpsi(1), 1)
      end if
    end do

  end select

  call batch_pack_was_modified(yy)

  call profiling_out(prof)
  POP_SUB(batch_copy_data)
end subroutine batch_copy_data

! --------------------------------------------------------------

subroutine batch_get_points_cl(this, sp, ep, psi, ldpsi)
  type(batch_t),       intent(in)    :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(opencl_mem_t),  intent(inout) :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif

  PUSH_SUB(batch_get_points_cl)
  call profiling_in(get_points_prof, "GET_POINTS")

  select case(batch_status(this))
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_CL_PACKED)

    tsize = types_get_size(batch_type(this))/types_get_size(TYPE_FLOAT)
    offset = this%index(1, 1) - 1

#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(kernel, 'points.cl', 'get_points')

    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, sp)
    call opencl_set_kernel_arg(kernel_ref, 1, ep)
    call opencl_set_kernel_arg(kernel_ref, 2, offset*tsize)
    call opencl_set_kernel_arg(kernel_ref, 3, this%nst_linear*tsize)
    call opencl_set_kernel_arg(kernel_ref, 4, this%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 5, this%pack%size_real(1))
    call opencl_set_kernel_arg(kernel_ref, 6, psi)
    call opencl_set_kernel_arg(kernel_ref, 7, ldpsi*tsize)

    call opencl_kernel_run(kernel_ref, (/this%pack%size_real(1), ep -&
      & sp + 1/), (/this%pack%size_real(1), 1/))
#endif
  end select

  call profiling_out(get_points_prof)

  POP_SUB(batch_get_points_cl)
end subroutine batch_get_points_cl

! --------------------------------------------------------------

subroutine batch_set_points_cl(this, sp, ep, psi, ldpsi)
  type(batch_t),       intent(inout) :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(opencl_mem_t),  intent(in)    :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif

  PUSH_SUB(batch_set_points_cl)
  call profiling_in(set_points_prof, "SET_POINTS")

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_CL_PACKED)

    tsize = types_get_size(batch_type(this))&
      &/types_get_size(TYPE_FLOAT)
    offset = this%index(1, 1) - 1
#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(kernel, 'points.cl', 'set_points')
    
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, sp)
    call opencl_set_kernel_arg(kernel_ref, 1, ep)
    call opencl_set_kernel_arg(kernel_ref, 2, offset*tsize)
    call opencl_set_kernel_arg(kernel_ref, 3, this%nst_linear*tsize)
    call opencl_set_kernel_arg(kernel_ref, 4, psi)
    call opencl_set_kernel_arg(kernel_ref, 5, ldpsi*tsize)
    call opencl_set_kernel_arg(kernel_ref, 6, this%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 7, this%pack%size_real(1))

    call opencl_kernel_run(kernel_ref, (/this%pack%size_real(1), ep -&
      & sp + 1/), (/this%pack%size_real(1), 1/))
#endif
  end select

  call profiling_out(set_points_prof)

  POP_SUB(batch_set_points_cl)
end subroutine batch_set_points_cl

#include "real.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"


end module batch_ops_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
