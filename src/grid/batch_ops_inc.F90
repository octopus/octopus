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

subroutine X(batch_axpy_const)(np, aa, xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist
  integer :: localsize

  PUSH_SUB(X(batch_axpy_const))
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
    call opencl_set_kernel_arg(kernel_zaxpy, 0, aa)
    call opencl_set_kernel_arg(kernel_zaxpy, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_zaxpy, 2, xx%pack%size(1))
    call opencl_set_kernel_arg(kernel_zaxpy, 3, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_zaxpy, 4, yy%pack%size(1))

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
  POP_SUB(X(batch_axpy_const))
end subroutine X(batch_axpy_const)

! --------------------------------------------------------------

subroutine X(batch_axpy_vec)(np, aa, xx, yy, a_start)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa(:)
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy
  integer, optional, intent(in)    :: a_start

  integer :: ist, ip, localsize, effsize
  R_TYPE, allocatable     :: aa_linear(:)
  CMPLX,  allocatable     :: zaa_linear(:)
  type(opencl_mem_t)      :: aa_buffer
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif
  
  PUSH_SUB(X(batch_axpy_vec))
  call profiling_in(axpy_prof, "BATCH_AXPY")

  ASSERT(batch_type(yy) == batch_type(xx))
#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  effsize = yy%nst_linear
  if(batch_is_packed(yy)) effsize = yy%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, yy%nst_linear
    aa_linear(ist) = aa(yy%index(ist, 1) - (optional_default(a_start, 1) - 1))
  end do

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, batch_type(yy), yy%pack%size(1))

    if(batch_type(yy) == TYPE_CMPLX) then
      ! convert aa_linear to complex
      SAFE_ALLOCATE(zaa_linear(1:yy%pack%size(1)))
      zaa_linear(1:yy%pack%size(1)) = aa_linear(1:yy%pack%size(1))
      call opencl_write_buffer(aa_buffer, yy%pack%size(1), zaa_linear)
      SAFE_DEALLOCATE_A(zaa_linear)
    else
      call opencl_write_buffer(aa_buffer, yy%pack%size(1), aa_linear)
    end if

    call octcl_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(axpy_vec)), flags = '-D'//R_TYPE_CL)
  
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, aa_buffer)
    call opencl_set_kernel_arg(kernel_ref, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(xx%pack%size(1)))
    call opencl_set_kernel_arg(kernel_ref, 3, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 4, log2(yy%pack%size(1)))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_ref, (/yy%pack%size(1), pad(np, localsize)/), (/yy%pack%size(1), localsize/yy%pack%size(1)/))

    call opencl_finish()
#endif
  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      do ist = 1, yy%pack%size(1)
        do ip = 1, np
          yy%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip) + yy%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      do ist = 1, yy%pack%size(1)
        do ip = 1, np
          yy%pack%dpsi(ist, ip) = aa_linear(ist)*xx%pack%dpsi(ist, ip) + yy%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call lalg_axpy(np, aa_linear(ist), xx%states_linear(ist)%zpsi, yy%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_axpy(np, aa_linear(ist), xx%states_linear(ist)%dpsi, yy%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(yy)

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(axpy_prof)
  POP_SUB(X(batch_axpy_vec))
end subroutine X(batch_axpy_vec)


! --------------------------------------------------------------

subroutine X(batch_scal_vec)(np, aa, xx, a_start)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa(:)
  type(batch_t),     intent(inout) :: xx
  integer, optional, intent(in)    :: a_start

  integer :: ist, ip, localsize, effsize
  R_TYPE, allocatable     :: aa_linear(:)
  CMPLX,  allocatable     :: zaa_linear(:)
  type(opencl_mem_t)      :: aa_buffer
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif
  
  PUSH_SUB(X(batch_axpy_vec))
  call profiling_in(axpy_prof, "BATCH_AXPY")

#ifdef R_TCMPLX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(xx) == TYPE_CMPLX)
#endif

  effsize = xx%nst_linear
  if(batch_is_packed(xx)) effsize = xx%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, xx%nst_linear
    aa_linear(ist) = aa(xx%index(ist, 1) - (optional_default(a_start, 1) - 1))
  end do
  
  select case(batch_status(xx))
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(aa_buffer, CL_MEM_READ_ONLY, batch_type(xx), xx%pack%size(1))

    if(batch_type(xx) == TYPE_CMPLX) then
      ! convert aa_linear to complex
      SAFE_ALLOCATE(zaa_linear(1:xx%pack%size(1)))
      zaa_linear(1:xx%pack%size(1)) = aa_linear(1:xx%pack%size(1))
      call opencl_write_buffer(aa_buffer, xx%pack%size(1), zaa_linear)
      SAFE_DEALLOCATE_A(zaa_linear)
    else
      call opencl_write_buffer(aa_buffer, xx%pack%size(1), aa_linear)
    end if

    call octcl_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(scal_vec)), flags = '-D'//R_TYPE_CL)
  
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, aa_buffer)
    call opencl_set_kernel_arg(kernel_ref, 1, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 2, log2(xx%pack%size(1)))

    localsize = opencl_max_workgroup_size()
    call opencl_kernel_run(kernel_ref, (/xx%pack%size(1), pad(np, localsize)/), (/xx%pack%size(1), localsize/xx%pack%size(1)/))

    call opencl_finish()
#endif
  case(BATCH_PACKED)
    if(batch_type(xx) == TYPE_CMPLX) then
      do ist = 1, xx%pack%size(1)
        do ip = 1, np
          xx%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      do ist = 1, xx%pack%size(1)
        do ip = 1, np
          xx%pack%dpsi(ist, ip) = aa_linear(ist)*xx%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, xx%nst_linear
      if(batch_type(xx) == TYPE_CMPLX) then
        call lalg_scal(np, aa_linear(ist), xx%states_linear(ist)%zpsi)
      else
#ifdef R_TREAL
        call lalg_scal(np, aa_linear(ist), xx%states_linear(ist)%dpsi)
#endif
      end if
    end do
  end select

  call batch_pack_was_modified(xx)

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(axpy_prof)
  POP_SUB(X(batch_axpy_vec))
end subroutine X(batch_scal_vec)


! --------------------------------------------------------------

subroutine X(batch_set_state1)(this, ist, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  integer :: ip
  type(profile_t), save :: prof
  type(opencl_mem_t) :: tmp

  call profiling_in(prof, "BATCH_SET_STATE")

  PUSH_SUB(X(batch_set_state1))

  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) this%states_linear(ist)%dpsi(ip) = psi(ip)
    else
      forall(ip = 1:np) this%states_linear(ist)%zpsi(ip) = psi(ip)
    end if
  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) this%pack%dpsi(ist, ip) = psi(ip)
    else
      forall(ip = 1:np) this%pack%zpsi(ist, ip) = psi(ip)
    end if
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))

    call opencl_write_buffer(tmp, np, psi)

    ! now call an opencl kernel to rearrange the data
    call opencl_set_kernel_arg(X(pack), 0, this%pack%size(1))
    call opencl_set_kernel_arg(X(pack), 1, ist - 1)
    call opencl_set_kernel_arg(X(pack), 2, tmp)
    call opencl_set_kernel_arg(X(pack), 3, this%pack%buffer)
    
    call opencl_kernel_run(X(pack), (/this%pack%size(2), 1/), (/opencl_max_workgroup_size(), 1/))
    
    call opencl_finish()

    call opencl_release_buffer(tmp)
#endif
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_set_state1))
end subroutine X(batch_set_state1)

! --------------------------------------------------------------

subroutine X(batch_set_state2)(this, index, np, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  PUSH_SUB(X(batch_set_state2))

  call X(batch_set_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_set_state2))
end subroutine X(batch_set_state2)

! --------------------------------------------------------------

subroutine X(batch_get_state1)(this, ist, np, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  integer :: ip
  type(profile_t), save :: prof 
  type(opencl_mem_t) :: tmp

  PUSH_SUB(X(batch_get_state1))

  call profiling_in(prof, "BATCH_GET_STATE")

  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) psi(ip) = this%states_linear(ist)%dpsi(ip)
    else
      forall(ip = 1:np) psi(ip) = this%states_linear(ist)%zpsi(ip)
    end if
  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      forall(ip = 1:np) psi(ip) = this%pack%dpsi(ist, ip)
    else
      forall(ip = 1:np) psi(ip) = this%pack%zpsi(ist, ip)
    end if
  case(BATCH_CL_PACKED)
#ifdef HAVE_OPENCL
    call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))

    call opencl_set_kernel_arg(X(unpack), 0, this%pack%size(1))
    call opencl_set_kernel_arg(X(unpack), 1, ist - 1)
    call opencl_set_kernel_arg(X(unpack), 2, this%pack%buffer)
    call opencl_set_kernel_arg(X(unpack), 3, tmp)

    call opencl_kernel_run(X(unpack), (/1, this%pack%size(2)/), (/1, opencl_max_workgroup_size()/))

    call opencl_finish()

    call opencl_read_buffer(tmp, np, psi)

    call opencl_release_buffer(tmp)
#endif
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_get_state1))
end subroutine X(batch_get_state1)

! --------------------------------------------------------------

subroutine X(batch_get_state2)(this, index, np, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  PUSH_SUB(X(batch_get_state2))

  call X(batch_get_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_get_state2))
end subroutine X(batch_get_state2)


! --------------------------------------------------------------

subroutine X(batch_get_points)(this, sp, ep, psi)
  type(batch_t),  intent(in)    :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(inout) :: psi(:, :, sp:)

  integer :: idim, ist, ii
  
  PUSH_SUB(X(batch_get_points))

#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then
      
      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%dpsi(sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%zpsi(sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then
      
      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        psi(ist, idim, sp:ep) = this%pack%dpsi(ii, sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        psi(ist, idim, sp:ep) = this%pack%zpsi(ii, sp:ep)
      end do

    end if

  case(BATCH_CL_PACKED)
    call messages_not_implemented('batch_get_points for CL packed batches')
  end select

  POP_SUB(X(batch_get_points))
end subroutine X(batch_get_points)

! --------------------------------------------------------------

subroutine X(batch_set_points)(this, sp, ep, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(in)    :: psi(:, :, sp:)

  integer :: idim, ist, ii

  PUSH_SUB(X(batch_set_points))

#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        this%states_linear(ii)%dpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        this%states_linear(ii)%zpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        this%pack%dpsi(ii, sp:ep) = psi(ist, idim, sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = this%index(ii, 1)
        idim = this%index(ii, 2)
        this%pack%zpsi(ii, sp:ep) = psi(ist, idim, sp:ep)
      end do

    end if

  case(BATCH_CL_PACKED)
    call messages_not_implemented('batch_set_points for CL packed batches')
  end select

  POP_SUB(X(batch_set_points))
end subroutine X(batch_set_points)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
