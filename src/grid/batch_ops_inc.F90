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

! --------------------------------------------------------------

subroutine X(batch_axpy_const)(np, aa, xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa
  class(batch_t),    intent(in)    :: xx
  class(batch_t),    intent(inout) :: yy

  integer :: ist, dim2, dim3
  integer :: localsize
  CMPLX :: zaa

  PUSH_SUB(X(batch_axpy_const))
  call profiling_in(axpy_const_prof, "BATCH_AXPY_CONST")

  call xx%check_compatibility_with(yy)
#ifdef R_TCOMPLEX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif

  select case(batch_status(xx))
  case(BATCH_DEVICE_PACKED)
    if(batch_type(yy) == TYPE_FLOAT) then

      call accel_set_kernel_arg(kernel_daxpy, 0, np)
      call accel_set_kernel_arg(kernel_daxpy, 1, aa)
      call accel_set_kernel_arg(kernel_daxpy, 2, xx%pack%buffer)
      call accel_set_kernel_arg(kernel_daxpy, 3, log2(xx%pack%size(1)))
      call accel_set_kernel_arg(kernel_daxpy, 4, yy%pack%buffer)
      call accel_set_kernel_arg(kernel_daxpy, 5, log2(yy%pack%size(1)))

      localsize = accel_kernel_workgroup_size(kernel_daxpy)/yy%pack%size(1)
      
      dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))

      call accel_kernel_run(kernel_daxpy, (/yy%pack%size(1), dim2, dim3/), (/yy%pack%size(1), localsize, 1/))
      
    else
      zaa = aa
      call accel_set_kernel_arg(kernel_zaxpy, 0, np)
      call accel_set_kernel_arg(kernel_zaxpy, 1, zaa)
      call accel_set_kernel_arg(kernel_zaxpy, 2, xx%pack%buffer)
      call accel_set_kernel_arg(kernel_zaxpy, 3, log2(xx%pack%size(1)))
      call accel_set_kernel_arg(kernel_zaxpy, 4, yy%pack%buffer)
      call accel_set_kernel_arg(kernel_zaxpy, 5, log2(yy%pack%size(1)))

      localsize = accel_kernel_workgroup_size(kernel_zaxpy)/yy%pack%size(1)

      dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
      dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
      
      call accel_kernel_run(kernel_zaxpy, (/yy%pack%size(1), dim2, dim3/), (/yy%pack%size(1), localsize, 1/))

    end if

    call accel_finish()
    
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

  call profiling_count_operations(xx%nst*np*(R_ADD + R_MUL)*types_get_size(batch_type(xx))/types_get_size(TYPE_FLOAT))

  call profiling_out(axpy_const_prof)
  POP_SUB(X(batch_axpy_const))
end subroutine X(batch_axpy_const)

! --------------------------------------------------------------

subroutine X(batch_axpy_vec)(np, aa, xx, yy, a_start, a_full)
  integer,            intent(in)    :: np
  R_TYPE,             intent(in)    :: aa(:)
  class(batch_t),     intent(in)    :: xx
  class(batch_t),     intent(inout) :: yy
  integer,  optional, intent(in)    :: a_start
  logical,  optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa, dim2, dim3
  R_TYPE, allocatable     :: aa_linear(:)
  integer :: localsize
  integer :: size_factor
  type(accel_mem_t)      :: aa_buffer
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(accel_kernel_t), save :: kernel
  
  PUSH_SUB(X(batch_axpy_vec))
  call profiling_in(axpy_vec_prof, "BATCH_AXPY_VEC")

  call xx%check_compatibility_with(yy)
#ifdef R_TCOMPLEX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif

  effsize = yy%nst_linear
  if(yy%is_packed()) effsize = yy%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, yy%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do

  select case(batch_status(xx))
  case(BATCH_DEVICE_PACKED)
    call accel_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(axpy_vec)), &
      flags = '-D' + R_TYPE_CL)

    if(batch_type(yy) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*yy%pack%size(1)))
      do ist = 1, yy%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 2*yy%pack%size(1))
      call accel_write_buffer(aa_buffer, 2*yy%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, yy%pack%size(1))
      call accel_write_buffer(aa_buffer, yy%pack%size(1), aa_linear)
    end if

    call accel_set_kernel_arg(kernel, 0, np)
    call accel_set_kernel_arg(kernel, 1, aa_buffer)
    call accel_set_kernel_arg(kernel, 2, xx%pack%buffer)
    call accel_set_kernel_arg(kernel, 3, log2(xx%pack%size(1)*size_factor))
    call accel_set_kernel_arg(kernel, 4, yy%pack%buffer)
    call accel_set_kernel_arg(kernel, 5, log2(yy%pack%size(1)*size_factor))

    localsize = accel_kernel_workgroup_size(kernel)/(yy%pack%size(1)*size_factor)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
    
    call accel_kernel_run(kernel, (/yy%pack%size(1)*size_factor, dim2, dim3/), (/yy%pack%size(1)*size_factor, localsize, 1/))

    call accel_finish()

    call accel_release_buffer(aa_buffer)

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
          yy%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip) + yy%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
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

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_count_operations(xx%nst*np*(R_ADD + R_MUL)*types_get_size(batch_type(xx))/types_get_size(TYPE_FLOAT))

  call profiling_out(axpy_vec_prof)
  POP_SUB(X(batch_axpy_vec))
end subroutine X(batch_axpy_vec)

! --------------------------------------------------------------

subroutine X(batch_scal_const)(np, aa, xx)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa
  class(batch_t),    intent(inout) :: xx

  R_TYPE, allocatable :: aavec(:)
  
  PUSH_SUB(X(batch_scal_const))

  SAFE_ALLOCATE(aavec(1:xx%nst))

  aavec(1:xx%nst) = aa

  call X(batch_scal_vec)(np, aavec, xx, a_full = .false.)
  
  SAFE_DEALLOCATE_A(aavec)
  
  POP_SUB(X(batch_scal_const))
end subroutine X(batch_scal_const)

! --------------------------------------------------------------

subroutine X(batch_scal_vec)(np, aa, xx, a_start, a_full)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: aa(:)
  class(batch_t),    intent(inout) :: xx
  integer, optional, intent(in)    :: a_start
  logical, optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa, dim2, dim3
  R_TYPE, allocatable     :: aa_linear(:)
  integer :: localsize
  integer :: size_factor
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(accel_mem_t)      :: aa_buffer
  type(accel_kernel_t), save :: kernel
  
  PUSH_SUB(X(batch_scal_vec))
  call profiling_in(scal_prof, "BATCH_SCAL")

#ifdef R_TCOMPLEX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(xx) == TYPE_CMPLX)
#endif

  effsize = xx%nst_linear
  if(xx%is_packed()) effsize = xx%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, xx%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do
  
  select case(batch_status(xx))
  case(BATCH_DEVICE_PACKED)
    if(batch_type(xx) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*xx%pack%size(1)))
      do ist = 1, xx%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 2*xx%pack%size(1))
      call accel_write_buffer(aa_buffer, 2*xx%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, xx%pack%size(1))
      call accel_write_buffer(aa_buffer, xx%pack%size(1), aa_linear)
    end if

    call accel_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(scal_vec)), flags = '-D' + R_TYPE_CL)

    call accel_set_kernel_arg(kernel, 0, np)
    call accel_set_kernel_arg(kernel, 1, aa_buffer)
    call accel_set_kernel_arg(kernel, 2, xx%pack%buffer)
    call accel_set_kernel_arg(kernel, 3, log2(xx%pack%size(1)*size_factor))

    localsize = accel_kernel_workgroup_size(kernel)/(xx%pack%size(1)*size_factor)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))

    call accel_kernel_run(kernel, (/xx%pack%size(1)*size_factor, dim2, dim3/), (/xx%pack%size(1)*size_factor, localsize, 1/))

    call accel_finish()

    call accel_release_buffer(aa_buffer)
    
  case(BATCH_PACKED)
    if(batch_type(xx) == TYPE_CMPLX) then
      !$omp parallel do
      do ip = 1, np
        do ist = 1, xx%pack%size(1)
          xx%pack%zpsi(ist, ip) = aa_linear(ist)*xx%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
     !$omp parallel do
     do ip = 1, np
        do ist = 1, xx%pack%size(1)
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

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(scal_prof)
  POP_SUB(X(batch_scal_vec))
end subroutine X(batch_scal_vec)

! --------------------------------------------------------------

subroutine X(batch_xpay_vec)(np, xx, aa, yy, a_start, a_full)
  integer,           intent(in)    :: np
  class(batch_t),    intent(in)    :: xx
  R_TYPE,            intent(in)    :: aa(:)
  class(batch_t),    intent(inout) :: yy
  integer, optional, intent(in)    :: a_start
  logical, optional, intent(in)    :: a_full

  integer :: ist, ip, effsize, iaa, dim2, dim3
  R_TYPE, allocatable     :: aa_linear(:)
  integer :: size_factor, localsize
  FLOAT,  allocatable     :: aa_linear_double(:)
  type(accel_mem_t)      :: aa_buffer
  type(accel_kernel_t), save :: kernel
  
  PUSH_SUB(X(batch_xpay_vec))
  call profiling_in(xpay_prof, "BATCH_XPAY")

  call xx%check_compatibility_with(yy)
#ifdef R_TCOMPLEX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif

  effsize = yy%nst_linear
  if(yy%is_packed()) effsize = yy%pack%size(1)
  SAFE_ALLOCATE(aa_linear(1:effsize))

  aa_linear = M_ZERO
  do ist = 1, yy%nst_linear
    iaa = batch_linear_to_ist(xx, ist) - (optional_default(a_start, 1) - 1)
    if(.not. optional_default(a_full, .true.)) iaa = iaa - (batch_linear_to_ist(xx, 1) - 1)
    aa_linear(ist) = aa(iaa)
  end do

  select case(batch_status(xx))
  case(BATCH_DEVICE_PACKED)
    if(batch_type(yy) == TYPE_CMPLX .and. R_TYPE_VAL == TYPE_FLOAT) then
      size_factor = 2
      SAFE_ALLOCATE(aa_linear_double(1:2*yy%pack%size(1)))
      do ist = 1, yy%pack%size(1)
        aa_linear_double(2*ist - 1) = aa_linear(ist)
        aa_linear_double(2*ist) = aa_linear(ist)
      end do
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, 2*yy%pack%size(1))
      call accel_write_buffer(aa_buffer, 2*yy%pack%size(1), aa_linear_double)
      SAFE_DEALLOCATE_A(aa_linear_double)
    else
      size_factor = 1
      call accel_create_buffer(aa_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, yy%pack%size(1))
      call accel_write_buffer(aa_buffer, yy%pack%size(1), aa_linear)
    end if

    call accel_kernel_start_call(kernel, 'axpy.cl', TOSTRING(X(xpay_vec)), flags = '-D' + R_TYPE_CL)

    call accel_set_kernel_arg(kernel, 0, np)
    call accel_set_kernel_arg(kernel, 1, aa_buffer)
    call accel_set_kernel_arg(kernel, 2, xx%pack%buffer)
    call accel_set_kernel_arg(kernel, 3, log2(xx%pack%size(1)*size_factor))
    call accel_set_kernel_arg(kernel, 4, yy%pack%buffer)
    call accel_set_kernel_arg(kernel, 5, log2(yy%pack%size(1)*size_factor))

    localsize = accel_kernel_workgroup_size(kernel)/(yy%pack%size(1)*size_factor)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
    
    call accel_kernel_run(kernel, (/yy%pack%size(1)*size_factor, dim2, dim3/), (/yy%pack%size(1)*size_factor, localsize, 1/))

    call accel_finish()

    call accel_release_buffer(aa_buffer)
    
  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
          yy%pack%zpsi(ist, ip) = xx%pack%zpsi(ist, ip) + aa_linear(ist)*yy%pack%zpsi(ist, ip)
        end do
      end do
    else
#ifdef R_TREAL
      !$omp parallel do private(ip, ist)
      do ip = 1, np
        do ist = 1, yy%pack%size(1)
          yy%pack%dpsi(ist, ip) = xx%pack%dpsi(ist, ip) + aa_linear(ist)*yy%pack%dpsi(ist, ip)
        end do
      end do
#endif
    end if
    
  case(BATCH_NOT_PACKED)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        !omp parallel do 
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = xx%states_linear(ist)%zpsi(ip) + aa_linear(ist)*yy%states_linear(ist)%zpsi(ip)
        end do
      else
#ifdef R_TREAL
        !omp parallel do
        do ip = 1, np
          yy%states_linear(ist)%dpsi(ip) = xx%states_linear(ist)%dpsi(ip) + aa_linear(ist)*yy%states_linear(ist)%dpsi(ip)
        end do
#endif
      end if
    end do
  end select

  call profiling_count_operations(xx%nst_linear*np*(R_ADD + R_MUL))

  SAFE_DEALLOCATE_A(aa_linear)

  call profiling_out(xpay_prof)
  POP_SUB(X(batch_xpay_vec))
end subroutine X(batch_xpay_vec)

! --------------------------------------------------------------

subroutine X(batch_xpay_const)(np, xx, aa, yy)
  integer,           intent(in)    :: np
  class(batch_t),    intent(in)    :: xx
  R_TYPE,            intent(in)    :: aa
  class(batch_t),    intent(inout) :: yy

  integer :: minst, maxst, ii, ist
  R_TYPE, allocatable :: aavec(:)
  
  minst = HUGE(minst)
  maxst = -HUGE(maxst)
  
  do ii = 1, xx%nst_linear
    ist = batch_linear_to_ist(xx, ii)
    minst = min(minst, ist)
    maxst = max(maxst, ist)
  end do


  SAFE_ALLOCATE(aavec(minst:maxst))

  aavec = aa

  call X(batch_xpay_vec)(np, xx, aavec, yy, a_start = minst)

  SAFE_DEALLOCATE_A(aavec)
  
end subroutine X(batch_xpay_const)

! --------------------------------------------------------------

subroutine X(batch_set_state1)(this, ist, np, psi)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  integer :: ip
  type(profile_t), save :: prof
  type(accel_mem_t) :: tmp
  FLOAT, allocatable :: zpsi(:)
  
  call profiling_in(prof, "BATCH_SET_STATE")

  PUSH_SUB(X(batch_set_state1))

  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
#ifdef R_TREAL
      call lalg_copy(np, psi, this%states_linear(ist)%dpsi)
#endif
    else
#ifdef R_TCOMPLEX
      call lalg_copy(np, psi, this%states_linear(ist)%zpsi)
#endif
    end if
  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      !omp parallel do 
      do ip = 1, np
        this%pack%dpsi(ist, ip) = psi(ip)
      end do
    else
      !omp parallel do
      do ip = 1, np
        this%pack%zpsi(ist, ip) = psi(ip)
      end do
    end if
  case(BATCH_DEVICE_PACKED)

    call accel_create_buffer(tmp, ACCEL_MEM_READ_ONLY, batch_type(this), this%pack%size(2))
    
    if(batch_type(this) /= R_TYPE_VAL) then

      ! this is not ideal, we should do the conversion on the GPU, so
      ! that we copy half of the data there
      
      SAFE_ALLOCATE(zpsi(1:np))
      !omp parallel do
      do ip = 1, np
        zpsi(ip) = psi(ip)
      end do
      
      call accel_write_buffer(tmp, np, zpsi)

      SAFE_DEALLOCATE_A(zpsi)
      
    else
      call accel_write_buffer(tmp, np, psi)
    end if


    ! now call an opencl kernel to rearrange the data
    call accel_set_kernel_arg(X(pack), 0, this%pack%size(1))
    call accel_set_kernel_arg(X(pack), 1, np)
    call accel_set_kernel_arg(X(pack), 2, ist - 1)
    call accel_set_kernel_arg(X(pack), 3, tmp)
    call accel_set_kernel_arg(X(pack), 4, this%pack%buffer)
    
    call accel_kernel_run(X(pack), (/this%pack%size(2), 1/), (/accel_max_workgroup_size(), 1/))
    
    call accel_finish()

    call accel_release_buffer(tmp)
    
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_set_state1))
end subroutine X(batch_set_state1)

! --------------------------------------------------------------

subroutine X(batch_set_state2)(this, index, np, psi)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:)

  PUSH_SUB(X(batch_set_state2))

  ASSERT(this%nst_linear > 0)
  call X(batch_set_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_set_state2))
end subroutine X(batch_set_state2)

! --------------------------------------------------------------

subroutine X(batch_set_state3)(this, ii, np, psi)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: ii
  integer,        intent(in)    :: np
  R_TYPE,         intent(in)    :: psi(:, :)

  integer :: i2

  PUSH_SUB(X(batch_set_state3))

  do i2 = 1, this%dim
    call X(batch_set_state1)(this, (ii - 1)*this%dim + i2, np, psi(:, i2))
  end do

  POP_SUB(X(batch_set_state3))
end subroutine X(batch_set_state3)

! --------------------------------------------------------------

subroutine X(batch_get_state1)(this, ist, np, psi)
  class(batch_t), intent(in)    :: this
  integer,        intent(in)    :: ist
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  integer :: ip
  type(profile_t), save :: prof 
  type(accel_mem_t) :: tmp
  FLOAT, allocatable :: dpsi(:)

  PUSH_SUB(X(batch_get_state1))

  call profiling_in(prof, "BATCH_GET_STATE")

  ASSERT(ubound(psi, dim = 1) >= np)
  ASSERT(ist >= 1 .and. ist <= this%nst_linear)
#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%states_linear(ist)%dpsi(ip)
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%states_linear(ist)%zpsi(ip)
      end do
      !$omp end parallel do
    end if

  case(BATCH_PACKED)
    if(batch_type(this) == TYPE_FLOAT) then
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%pack%dpsi(ist, ip)
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do ip = 1, np
        psi(ip) = this%pack%zpsi(ist, ip)
      end do
      !$omp end parallel do
    end if

  case(BATCH_DEVICE_PACKED)

    ASSERT(np <= this%pack%size(2))

    call accel_create_buffer(tmp, ACCEL_MEM_WRITE_ONLY, batch_type(this), this%pack%size(2))
    
    if(batch_type(this) == R_TYPE_VAL) then

      call accel_set_kernel_arg(X(unpack), 0, this%pack%size(1))
      call accel_set_kernel_arg(X(unpack), 1, np)
      call accel_set_kernel_arg(X(unpack), 2, ist - 1)
      call accel_set_kernel_arg(X(unpack), 3, this%pack%buffer)
      call accel_set_kernel_arg(X(unpack), 4, tmp)
      
      call accel_kernel_run(X(unpack), (/1, this%pack%size(2)/), (/1, accel_max_workgroup_size()/))
      
      call accel_finish()
      
      call accel_read_buffer(tmp, np, psi)

    else

      ! the output buffer is complex, we get it as real
      
      ASSERT(batch_type(this) == TYPE_FLOAT)
      ASSERT(R_TYPE_VAL == TYPE_CMPLX)
      
      call accel_set_kernel_arg(dunpack, 0, this%pack%size(1))
      call accel_set_kernel_arg(dunpack, 1, np)
      call accel_set_kernel_arg(dunpack, 2, ist - 1)
      call accel_set_kernel_arg(dunpack, 3, this%pack%buffer)
      call accel_set_kernel_arg(dunpack, 4, tmp)
      
      call accel_kernel_run(dunpack, (/1, this%pack%size(2)/), (/1, accel_max_workgroup_size()/))
      
      SAFE_ALLOCATE(dpsi(1:np))
      
      call accel_finish()

      call accel_read_buffer(tmp, np, dpsi)

      ! and convert to complex on the cpu
      
      !omp parallel do
      do ip = 1, np
        psi(ip) = dpsi(ip)
      end do

      SAFE_DEALLOCATE_A(dpsi)

    end if

    call accel_release_buffer(tmp)
    
  end select

  call profiling_out(prof)

  POP_SUB(X(batch_get_state1))
end subroutine X(batch_get_state1)

! --------------------------------------------------------------

subroutine X(batch_get_state2)(this, index, np, psi)
  class(batch_t), intent(in)    :: this
  integer,        intent(in)    :: index(:)
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:)

  PUSH_SUB(X(batch_get_state2))

  ASSERT(this%nst_linear > 0)
  call X(batch_get_state1)(this, batch_inv_index(this, index), np, psi)

  POP_SUB(X(batch_get_state2))
end subroutine X(batch_get_state2)


! --------------------------------------------------------------

subroutine X(batch_get_state3)(this, ii, np, psi)
  class(batch_t), intent(in)    :: this
  integer,        intent(in)    :: ii
  integer,        intent(in)    :: np
  R_TYPE,         intent(inout) :: psi(:, :)

  integer :: i2

  PUSH_SUB(X(batch_get_state3))

  do i2 = 1, this%dim
    call X(batch_get_state1)(this, (ii - 1)*this%dim + i2, np, psi(:, i2))
  end do

  POP_SUB(X(batch_get_state3))
end subroutine X(batch_get_state3)

! --------------------------------------------------------------

subroutine X(batch_get_points)(this, sp, ep, psi)
  class(batch_t), intent(in)    :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(inout) :: psi(:, :, sp:)

  integer :: idim, ist, ii, ip
  
  PUSH_SUB(X(batch_get_points))
  call profiling_in(get_points_prof, 'GET_POINTS')

#ifdef R_TREAL
  ! cannot get a real value from a complex batch
  ASSERT(batch_type(this) /= TYPE_CMPLX)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then
      
      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%dpsi(sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        psi(ist, idim, sp:ep) = this%states_linear(ii)%zpsi(sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          psi(ist, idim, ip) = this%pack%dpsi(ii, ip)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          psi(ist, idim, ip) = this%pack%zpsi(ii, ip)
        end do
      end do
      !$omp end parallel do

    end if

  case(BATCH_DEVICE_PACKED)
    call messages_not_implemented('batch_get_points for CL packed batches')
  end select

  call profiling_out(get_points_prof)

  POP_SUB(X(batch_get_points))
end subroutine X(batch_get_points)

! --------------------------------------------------------------

subroutine X(batch_set_points)(this, sp, ep, psi)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: sp  
  integer,        intent(in)    :: ep
  R_TYPE,         intent(in)    :: psi(:, :, sp:)

  integer :: idim, ist, ii, ip

  PUSH_SUB(X(batch_set_points))

  call profiling_in(set_points_prof, 'SET_POINTS')

#ifdef R_TCOMPLEX
  ! cannot set a real batch with complex values
  ASSERT(batch_type(this) /= TYPE_FLOAT)
#endif

  select case(batch_status(this))
  case(BATCH_NOT_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        this%states_linear(ii)%dpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    else

      do ii = 1, this%nst_linear
        ist = batch_linear_to_ist(this, ii)
        idim = batch_linear_to_idim(this, ii)
        this%states_linear(ii)%zpsi(sp:ep) = psi(ist, idim, sp:ep)
      end do

    end if

  case(BATCH_PACKED)

    if(batch_type(this) == TYPE_FLOAT) then

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          this%pack%dpsi(ii, ip) = psi(ist, idim, ip)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(ip, ii, ist, idim)
      do ip = sp, ep
        do ii = 1, this%nst_linear
          ist = batch_linear_to_ist(this, ii)
          idim = batch_linear_to_idim(this, ii)
          this%pack%zpsi(ii, ip) = psi(ist, idim, ip)
        end do
      end do
      !$omp end parallel do

    end if

  case(BATCH_DEVICE_PACKED)
    call messages_not_implemented('batch_set_points for CL packed batches')
  end select

  call profiling_out(set_points_prof)

  POP_SUB(X(batch_set_points))
end subroutine X(batch_set_points)

! --------------------------------------------------------------

subroutine X(batch_mul)(np, ff,  xx, yy)
  integer,           intent(in)    :: np
  R_TYPE,            intent(in)    :: ff(:)
  class(batch_t),    intent(in)    :: xx
  class(batch_t),    intent(inout) :: yy

  integer :: ist, ip
  R_TYPE :: mul
#if defined(R_TREAL)
  integer :: iprange
  type(accel_mem_t) :: ff_buffer
#endif
  
  PUSH_SUB(X(batch_mul))
  call profiling_in(mul_prof, "BATCH_MUL")

  call xx%check_compatibility_with(yy)
#ifdef R_TCOMPLEX
  !if aa is complex, the functions must be complex
  ASSERT(batch_type(yy) == TYPE_CMPLX)
#endif

  select case(batch_status(yy))
  case(BATCH_DEVICE_PACKED)

#if defined(R_TREAL)

    ! We reuse here the routine to apply the local potential
    call batch_set_zero(yy)
    
    call accel_create_buffer(ff_buffer, ACCEL_MEM_READ_ONLY, R_TYPE_VAL, np)
    call accel_write_buffer(ff_buffer, np, ff)

    call accel_set_kernel_arg(kernel_vpsi, 0, 0)
    call accel_set_kernel_arg(kernel_vpsi, 1, np)
    call accel_set_kernel_arg(kernel_vpsi, 2, ff_buffer)
    call accel_set_kernel_arg(kernel_vpsi, 3, xx%pack%buffer)
    call accel_set_kernel_arg(kernel_vpsi, 4, log2(xx%pack%size_real(1)))
    call accel_set_kernel_arg(kernel_vpsi, 5, yy%pack%buffer)
    call accel_set_kernel_arg(kernel_vpsi, 6, log2(yy%pack%size_real(1)))
    
    iprange = accel_kernel_workgroup_size(kernel_vpsi)/xx%pack%size_real(1)
    
    call accel_kernel_run(kernel_vpsi, (/xx%pack%size_real(1), pad(np, iprange)/), (/xx%pack%size_real(1), iprange/))

    call accel_release_buffer(ff_buffer)
#else
    call messages_not_implemented("OpenCL batch_mul")
#endif

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      !$omp parallel do private(ip, ist, mul)
      do ip = 1, np
        mul = ff(ip)
        do ist = 1, yy%nst_linear
          yy%pack%zpsi(ist, ip) = mul*xx%pack%zpsi(ist, ip)
        end do
      end do
      !$omp end parallel do
    else
#ifdef R_TREAL
      !$omp parallel do private(ip, ist, mul)
      do ip = 1, np
        mul = ff(ip)
        do ist = 1, yy%nst_linear
          yy%pack%dpsi(ist, ip) = mul*xx%pack%dpsi(ist, ip)
        end do
      end do
      !$omp end parallel do
#endif
    end if

  case(BATCH_NOT_PACKED)
    if(batch_type(yy) == TYPE_CMPLX) then
      do ist = 1, yy%nst_linear
        !$omp parallel do
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = ff(ip)*xx%states_linear(ist)%zpsi(ip)
        end do
        !$omp end parallel do
      end do
    else
#ifdef R_TREAL
      do ist = 1, yy%nst_linear
        !$omp parallel do
        do ip = 1, np
          yy%states_linear(ist)%zpsi(ip) = ff(ip)*xx%states_linear(ist)%zpsi(ip)
        end do
        !$omp end parallel do
      end do

#endif
    end if
  end select

  call profiling_out(mul_prof)
  POP_SUB(X(batch_mul))

end subroutine X(batch_mul)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
