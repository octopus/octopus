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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module batch_ops_oct_m
  use accel_oct_m
  use batch_oct_m
  use blas_oct_m
  use iso_c_binding
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use messages_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                         &
    batch_set_zero,                 &
    batch_axpy,                     &
    batch_scal,                     &
    batch_xpay,                     &
    batch_set_state,                &
    batch_get_state,                &
    batch_get_points,               &
    batch_set_points,               &
    batch_points_block_size,        &
    batch_mul,                      &
    dbatch_axpy_function,                 &
    zbatch_axpy_function

  interface batch_axpy
    module procedure dbatch_axpy_const
    module procedure zbatch_axpy_const
    module procedure dbatch_axpy_vec
    module procedure zbatch_axpy_vec
  end interface batch_axpy

  interface batch_scal
    module procedure dbatch_scal_const
    module procedure zbatch_scal_const
    module procedure dbatch_scal_vec
    module procedure zbatch_scal_vec
  end interface batch_scal

  interface batch_xpay
    module procedure dbatch_xpay_vec
    module procedure zbatch_xpay_vec
    module procedure dbatch_xpay_const
    module procedure zbatch_xpay_const
  end interface batch_xpay

  interface batch_set_state
    module procedure dbatch_set_state1
    module procedure zbatch_set_state1
    module procedure dbatch_set_state2
    module procedure zbatch_set_state2
    module procedure dbatch_set_state3
    module procedure zbatch_set_state3
  end interface batch_set_state

  interface batch_get_state
    module procedure dbatch_get_state1
    module procedure zbatch_get_state1
    module procedure dbatch_get_state2
    module procedure zbatch_get_state2
    module procedure dbatch_get_state3
    module procedure zbatch_get_state3
  end interface batch_get_state

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

  interface batch_mul
    module procedure dbatch_mul
    module procedure zbatch_mul
  end interface batch_mul

  type(profile_t), save :: scal_prof, xpay_prof, axpy_const_prof, axpy_vec_prof, get_points_prof, set_points_prof
  type(profile_t), save :: mul_prof

contains

  !--------------------------------------------------------------

  subroutine batch_set_zero(this)
    class(batch_t),     intent(inout) :: this

    type(profile_t), save :: prof
    integer :: ist_linear, ist, ip

    PUSH_SUB(batch_set_zero)

    call profiling_in(prof, "BATCH_SET_ZERO")

    select case(this%status())
    case(BATCH_DEVICE_PACKED)
      call accel_set_buffer_to_zero(this%ff_device, this%type(), product(this%pack_size))

    case(BATCH_PACKED)
      if(this%type() == TYPE_FLOAT) then
        !$omp parallel do schedule(static)
        do ip = 1, this%pack_size(2)
          do ist = 1, this%pack_size(1)
            this%dff_pack(ist, ip) = M_ZERO
          end do
        end do
      else
        !$omp parallel do schedule(static)
        do ip = 1, this%pack_size(2)
          do ist = 1, this%pack_size(1)
            this%zff_pack(ist, ip) = M_z0
          end do
        end do
      end if

    case(BATCH_NOT_PACKED)
      if(this%type() == TYPE_FLOAT) then
        do ist_linear = 1, this%nst_linear
          !$omp parallel do schedule(static)
          do ip = 1, ubound(this%dff_linear, dim=1)
            this%dff_linear(ip, ist_linear) = M_ZERO
          end do
        end do
      else
        do ist_linear = 1, this%nst_linear
          !$omp parallel do schedule(static)
          do ip = 1, ubound(this%zff_linear, dim=1)
            this%zff_linear(ip, ist_linear) = M_z0
          end do
        end do
      end if

    case default
      ASSERT(.false.)
      
    end select

    call profiling_out(prof)

    POP_SUB(batch_set_zero)
  end subroutine batch_set_zero

! --------------------------------------------------------------

subroutine batch_get_points_cl(this, sp, ep, psi, ldpsi)
  class(batch_t),      intent(in)    :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(accel_mem_t),   intent(inout) :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
  type(accel_kernel_t), save :: kernel

  PUSH_SUB(batch_get_points_cl)
  call profiling_in(get_points_prof, "GET_POINTS")

  select case(this%status())
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_DEVICE_PACKED)

    tsize = types_get_size(this%type())/types_get_size(TYPE_FLOAT)
    offset = this%linear_to_ist(1) - 1

    call accel_kernel_start_call(kernel, 'points.cl', 'get_points')

    call accel_set_kernel_arg(kernel, 0, sp)
    call accel_set_kernel_arg(kernel, 1, ep)
    call accel_set_kernel_arg(kernel, 2, offset*tsize)
    call accel_set_kernel_arg(kernel, 3, this%nst_linear*tsize)
    call accel_set_kernel_arg(kernel, 4, this%ff_device)
    call accel_set_kernel_arg(kernel, 5, this%pack_size_real(1))
    call accel_set_kernel_arg(kernel, 6, psi)
    call accel_set_kernel_arg(kernel, 7, ldpsi*tsize)

    call accel_kernel_run(kernel, (/this%pack_size_real(1), ep - sp + 1/), (/this%pack_size_real(1), 1/))

  end select

  call profiling_out(get_points_prof)

  POP_SUB(batch_get_points_cl)
end subroutine batch_get_points_cl

! --------------------------------------------------------------

subroutine batch_set_points_cl(this, sp, ep, psi, ldpsi)
  class(batch_t),      intent(inout) :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(accel_mem_t),   intent(in)    :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
  type(accel_kernel_t), save :: kernel

  PUSH_SUB(batch_set_points_cl)
  call profiling_in(set_points_prof, "SET_POINTS")

  select case(this%status())
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_DEVICE_PACKED)

    tsize = types_get_size(this%type())&
      /types_get_size(TYPE_FLOAT)
    offset = this%linear_to_ist(1) - 1

    call accel_kernel_start_call(kernel, 'points.cl', 'set_points')
    
    call accel_set_kernel_arg(kernel, 0, sp)
    call accel_set_kernel_arg(kernel, 1, ep)
    call accel_set_kernel_arg(kernel, 2, offset*tsize)
    call accel_set_kernel_arg(kernel, 3, this%nst_linear*tsize)
    call accel_set_kernel_arg(kernel, 4, psi)
    call accel_set_kernel_arg(kernel, 5, ldpsi*tsize)
    call accel_set_kernel_arg(kernel, 6, this%ff_device)
    call accel_set_kernel_arg(kernel, 7, this%pack_size_real(1))

    call accel_kernel_run(kernel, (/this%pack_size_real(1), ep - sp + 1/), (/this%pack_size_real(1), 1/))

  end select

  call profiling_out(set_points_prof)

  POP_SUB(batch_set_points_cl)
end subroutine batch_set_points_cl

! -------------------------

integer pure function batch_points_block_size(this) result(block_size)
  class(batch_t),       intent(in)    :: this
  
  block_size = 61440

end function batch_points_block_size

#include "real.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_ops_inc.F90"
#include "undef.F90"

end module batch_ops_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
