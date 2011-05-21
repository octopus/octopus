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
!! $Id: batch_inc.F90 4298 2008-06-18 15:03:00Z dstrubbe $

!--------------------------------------------------------------
subroutine X(batch_init_contiguous)(this, dim, st_start, st_end, psi)
  type(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: dim
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, intent(in)    :: psi(:, :, st_start:)

  integer :: ist

  PUSH_SUB(X(batch_init_contiguous))

  ASSERT(st_end >= st_start)

  call batch_init_empty(this, dim, st_end - st_start + 1)

  this%X(psicont) => psi(:, :, st_start:)

  ASSERT(ubound(psi, dim = 3) >= st_end)

  do ist = st_start, st_end
    call X(batch_add_state)(this, ist, psi(:, :, ist))
  end do

  POP_SUB(X(batch_init_contiguous))
end subroutine X(batch_init_contiguous)

!--------------------------------------------------------------
subroutine X(batch_add_state)(this, ist, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ist
  R_TYPE, target, intent(in)    :: psi(:, :)

  integer :: idim, ii

  PUSH_SUB(X(batch_add_state))

  ASSERT(this%current <= this%nst)

  this%states(this%current)%ist    =  ist
  this%states(this%current)%X(psi) => psi

  ! now we also populate the linear array
  do idim = 1, this%dim
    ii = this%dim*(this%current - 1) + idim
    this%states_linear(ii)%X(psi) => psi(:, idim)
    this%index(ii, 1) = ist
    this%index(ii, 2) = idim
  end do

  this%current = this%current + 1

  POP_SUB(X(batch_add_state))
end subroutine X(batch_add_state)

!--------------------------------------------------------------

subroutine X(batch_add_state_linear)(this, psi)
  type(batch_t),  intent(inout) :: this
  R_TYPE, target, intent(in)    :: psi(:)

  PUSH_SUB(X(batch_add_state_linear))

  ASSERT(this%current <= this%nst_linear)
  this%states_linear(this%current)%X(psi) => psi
  this%index(this%current, 1) = this%current

  this%current = this%current + 1

  POP_SUB(X(batch_add_state_linear))
end subroutine X(batch_add_state_linear)


!--------------------------------------------------------------
subroutine X(batch_new)(this, st_start, st_end, np)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  integer,        intent(in)    :: np

  integer :: ist

  PUSH_SUB(X(batch_new))

  SAFE_ALLOCATE(this%X(psicont)(1:np, 1:this%dim, 1:st_end - st_start + 1))

  do ist = st_start, st_end
    call X(batch_add_state)(this, ist,this%X(psicont)(:, :, ist - st_start + 1))
  end do

  POP_SUB(X(batch_new))
end subroutine X(batch_new)


!--------------------------------------------------------------
subroutine X(batch_delete)(this)
  type(batch_t),  intent(inout) :: this

  integer :: ii

  PUSH_SUB(X(batch_delete))

  do ii = 1, this%nst
    nullify(this%states(ii)%X(psi))
  end do

  do ii = 1, this%nst_linear
    nullify(this%states_linear(ii)%X(psi))
  end do

  this%current = 1

  SAFE_DEALLOCATE_P(this%X(psicont))

  POP_SUB(X(batch_delete))
end subroutine X(batch_delete)


!--------------------------------------------------------------

subroutine X(batch_set)(this, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:, :, :)

  integer :: ist, idim

  PUSH_SUB(X(batch_set))

  do ist = 1, this%nst
    do idim = 1, this%dim
      call lalg_copy(np, psi(:, idim, ist), this%states(ist)%X(psi)(:, idim))
    end do
  end do

  POP_SUB(X(batch_set))
end subroutine X(batch_set)

! --------------------------------------------------------------

subroutine X(batch_axpy)(np, aa, xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist
#ifdef HAVE_OPENCL
  integer :: localsize
#endif

  PUSH_SUB(X(batch_axpy))
  call profiling_in(axpy_prof, "BATCH_AXPY")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
#ifdef R_TREAL

    call opencl_set_kernel_arg(kernel_daxpy, 0, aa)
    call opencl_set_kernel_arg(kernel_daxpy, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_daxpy, 2, log2(xx%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_daxpy, 3, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_daxpy, 4, log2(yy%pack%size_real(1)))

    localsize = opencl_max_workgroup_size()/yy%pack%size_real(1)
    call opencl_kernel_run(kernel_daxpy, (/yy%pack%size_real(1), pad(np, localsize)/), (/yy%pack%size_real(1), localsize/))

#else
    call opencl_set_kernel_arg(kernel_zaxpy, 0, real(aa, REAL_PRECISION))
    call opencl_set_kernel_arg(kernel_zaxpy, 1, aimag(aa))
    call opencl_set_kernel_arg(kernel_zaxpy, 2, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_zaxpy, 3, xx%pack%size(1))
    call opencl_set_kernel_arg(kernel_zaxpy, 4, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_zaxpy, 5, yy%pack%size(1))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_zaxpy, (/yy%pack%size(1), pad(np, localsize)/), (/yy%pack%size(1), localsize/yy%pack%size(1)/))

#endif

    call opencl_finish()
#endif

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      call lalg_axpy(xx%pack%size(1), np, aa, xx%pack%zpsi, yy%pack%zpsi)
    else
#ifdef R_TREAL
      call lalg_axpy(xx%pack%size(1), np, aa, xx%pack%dpsi, yy%pack%dpsi)
#endif
    end if

  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call lalg_axpy(np, aa, xx%states_linear(ist)%zpsi, yy%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_axpy(np, aa, xx%states_linear(ist)%dpsi, yy%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(yy)

  call profiling_out(axpy_prof)
  POP_SUB(X(batch_axpy))
end subroutine X(batch_axpy)

! --------------------------------------------------------------


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
