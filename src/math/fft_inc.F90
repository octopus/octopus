!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011 J. Alberdi-Rodriguez, P. Garcia Risueño, M. Oliveira
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


! ---------------------------------------------------------
subroutine X(fft_forward_3d)(fft, in, out, norm)
    type(fft_t),     intent(in)    :: fft
    R_TYPE,          intent(inout) :: in(:,:,:)
    CMPLX,           intent(out)   :: out(:,:,:)
    FLOAT, optional, intent(out)   :: norm 

    integer :: ii, jj, kk, slot, n1, n2, n3
    type(profile_t), save :: prof_fw
#ifdef HAVE_CLFFT
    CMPLX, allocatable :: cin(:, :, :)
    type(accel_mem_t) :: rsbuffer, fsbuffer
#endif

    PUSH_SUB(X(fft_forward_3d))

    call profiling_in(prof_fw, TOSTRING(X(FFT_FORWARD)))

    
    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
      if(all(fft_array(slot)%rs_n(1:3) >= 1)) then
#ifdef R_TREAL
        call fftw_execute_dft_r2c(fft_array(slot)%planf, in(:,:,:), out(:,:,:))
#else
        call fftw_execute_dft(fft_array(slot)%planf, in(:,:,:), out(:,:,:))
#endif
      else
        ii = min(1, fft_array(slot)%rs_n(1))
        jj = min(1, fft_array(slot)%rs_n(2))
        kk = min(1, fft_array(slot)%rs_n(3))
#ifdef R_TREAL
        call fftw_execute_dft_r2c(fft_array(slot)%planf, in(ii:,jj:,kk:), out(ii:,jj:,kk:))
#else
        call fftw_execute_dft(fft_array(slot)%planf, in(ii:,jj:,kk:), out(ii:,jj:,kk:))
#endif
      end if
    case (FFTLIB_NFFT)
      call X(nfft_forward)(fft_array(slot)%nfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%nfft%norm

    case (FFTLIB_PNFFT)
      call X(pnfft_forward)(fft_array(slot)%pnfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%pnfft%norm

    case (FFTLIB_PFFT)
      if (all(fft_array(slot)%rs_n /= 0)) then
        ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == in(1,1,1))
      end if
      if (all(fft_array(slot)%fs_n /= 0)) then
        ASSERT(fft_array(slot)%fs_data(1,1,1) == out(1,1,1))
      end if
      n1 = max(1, fft_array(slot)%rs_n(1))
      n2 = max(1, fft_array(slot)%rs_n(2))
      n3 = max(1, fft_array(slot)%rs_n(3))
#ifdef R_TREAL
      if (fft_array(slot)%rs_n(1) > fft_array(slot)%rs_n_global(1)) then
        do kk = 1, n3
          do jj = 1, n2
            fft_array(slot)%drs_data(n1, jj, kk) = M_ZERO
          end do
        end do
      end if
#endif 
#ifdef HAVE_PFFT
      call pfft_execute(fft_array(slot)%planf)
#endif
    case(FFTLIB_ACCEL)
#ifdef HAVE_CLFFT

      SAFE_ALLOCATE(cin(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)))

      cin(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)) = &
        in(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3))

      call accel_create_buffer(rsbuffer, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%rs_n(1:3)))
      call accel_create_buffer(fsbuffer, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%fs_n(1:3)))

      call accel_write_buffer(rsbuffer, product(fft_array(slot)%rs_n(1:3)), cin)

      call accel_finish()

      call clfftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, accel%command_queue, &
        rsbuffer%mem, fsbuffer%mem, cl_status)
      if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftEnqueueTransform')

      call accel_finish()

      call accel_read_buffer(fsbuffer, product(fft_array(slot)%fs_n(1:3)), out)

      call accel_finish()

      call accel_release_buffer(rsbuffer)
      call accel_release_buffer(fsbuffer)
#endif
    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select

    call fft_operation_count(fft)

    call profiling_out(prof_fw)

    POP_SUB(X(fft_forward_3d))
  end subroutine X(fft_forward_3d)

! ---------------------------------------------------------
  subroutine X(fft_forward_accel)(fft, in, out)
    type(fft_t),       intent(in)    :: fft
    type(accel_mem_t), intent(in)    :: in
    type(accel_mem_t), intent(inout) :: out

    integer :: slot
    type(profile_t), save :: prof_fw
#ifdef HAVE_CLFFT
    type(accel_mem_t)         :: tmp_buf
    integer                    :: bsize
    integer(8)                 :: tmp_buf_size
#endif

    PUSH_SUB(X(fft_forward_accel))

    call profiling_in(prof_fw, TOSTRING(X(FFT_FORWARD_ACCEL)))

    slot = fft%slot
    ASSERT(fft_array(slot)%library == FFTLIB_ACCEL)

#ifdef HAVE_CUDA

#ifdef R_TREAL
    call cuda_fft_execute_d2z(fft_array(slot)%cuda_plan_fw, in%mem, out%mem)
#else
    call cuda_fft_execute_z2z_forward(fft_array(slot)%cuda_plan_fw, in%mem, out%mem)
#endif
    call fft_operation_count(fft)
    call accel_finish()

#endif
    
#ifdef HAVE_CLFFT

    call clfftGetTmpBufSize(fft_array(slot)%cl_plan_bw, tmp_buf_size, cl_status)
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftGetTmpBufSize')

    if(tmp_buf_size > 0) then
      call accel_create_buffer(tmp_buf, ACCEL_MEM_READ_WRITE, TYPE_BYTE, int(tmp_buf_size, 4))
    end if

    if(tmp_buf_size > 0) then
      call clfftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, accel%command_queue, &
        in%mem, out%mem, tmp_buf%mem, cl_status)
    else
      call clfftEnqueueTransform(fft_array(slot)%cl_plan_fw, CLFFT_FORWARD, accel%command_queue, &
        in%mem, out%mem, cl_status)
    end if
    
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftEnqueueTransform')
    
    call fft_operation_count(fft)

    call accel_finish()

    if(tmp_buf_size > 0) call accel_release_buffer(tmp_buf)

#endif

    call profiling_out(prof_fw)

    POP_SUB(X(fft_forward_accel))
  end subroutine X(fft_forward_accel)

  ! ---------------------------------------------------------

  subroutine X(fft_forward_1d)(fft, in, out)
    type(fft_t), intent(in)    :: fft
    R_TYPE,      intent(inout) :: in(:)
    CMPLX,       intent(out)   :: out(:)

    PUSH_SUB(X(fft_forward_1d))

#ifdef R_TREAL
    call fftw_execute_dft_r2c(fft_array(fft%slot)%planf, in, out)
#else
    call fftw_execute_dft(fft_array(fft%slot)%planf, in, out)
#endif
    call fft_operation_count(fft)

    POP_SUB(X(fft_forward_1d))
  end subroutine X(fft_forward_1d)

  ! ---------------------------------------------------------
  subroutine X(fft_backward_3d)(fft, in, out, norm)
    type(fft_t),     intent(in)    :: fft
    CMPLX,           intent(inout) :: in(:,:,:)
    R_TYPE,          intent(out)   :: out(:,:,:)
    FLOAT, optional, intent(out)   :: norm 
    
    integer :: ii, jj, kk, slot
    FLOAT :: scaling_factor
    type(profile_t), save :: prof_bw
    logical :: scale
#ifdef HAVE_CLFFT
    CMPLX, allocatable :: cout(:, :, :)
    type(accel_mem_t) :: rsbuffer, fsbuffer
#endif

    PUSH_SUB(X(fft_backward_3d))
    
    call profiling_in(prof_bw,TOSTRING(X(FFT_BACKWARD)))

    scale = .true.

    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
#ifdef R_TREAL
      call fftw_execute_dft_c2r(fft_array(slot)%planb, in, out)
#else
      call fftw_execute_dft(fft_array(slot)%planb, in, out)
#endif
    case (FFTLIB_NFFT)
      scale = .false. ! the result is already scaled
      call X(nfft_backward)(fft_array(slot)%nfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%nfft%norm

    case (FFTLIB_PNFFT)
      scale = .false. ! the result is already scaled
      call X(pnfft_backward)(fft_array(slot)%pnfft, in(:,:,:), out(:,:,:))
      if(present(norm)) norm = fft_array(slot)%pnfft%norm

    case (FFTLIB_PFFT)
      if (all(fft_array(slot)%fs_n /= 0)) then
        ASSERT(fft_array(slot)%fs_data(1,1,1) == in(1,1,1))
      end if  
      if (all(fft_array(slot)%rs_n /= 0)) then
        ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == out(1,1,1))
      end if
#ifdef HAVE_PFFT
      call pfft_execute(fft_array(slot)%planb)
#endif
    case(FFTLIB_ACCEL)
#ifdef HAVE_CLFFT

      call accel_create_buffer(rsbuffer, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%rs_n(1:3)))
      call accel_create_buffer(fsbuffer, ACCEL_MEM_READ_WRITE, TYPE_CMPLX, product(fft_array(slot)%fs_n(1:3)))

      call accel_write_buffer(fsbuffer, product(fft_array(slot)%fs_n(1:3)), in)

      call accel_finish()

      call clfftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, accel%command_queue, &
        fsbuffer%mem, rsbuffer%mem, cl_status)
      if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftEnqueueTransform')

      call accel_finish()

      SAFE_ALLOCATE(cout(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)))

      call accel_read_buffer(rsbuffer, product(fft_array(slot)%rs_n(1:3)), cout)

      call accel_finish()

      out(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3)) = &
        cout(1:fft_array(slot)%rs_n(1), 1:fft_array(slot)%rs_n(2), 1:fft_array(slot)%rs_n(3))

      call accel_release_buffer(rsbuffer)
      call accel_release_buffer(fsbuffer)
      
      scale = .false. ! scaling is done by the library
#endif
    case default
      call messages_write('Invalid FFT library.')
      call messages_fatal()
    end select

    if(scale) then
      ! multiply by 1/(N1*N2*N3)
      ! separate divisions to avoid integer overflow
      scaling_factor = M_ONE/TOFLOAT(fft_array(slot)%rs_n_global(1))
      scaling_factor = scaling_factor/TOFLOAT(fft_array(slot)%rs_n_global(2))
      scaling_factor = scaling_factor/TOFLOAT(fft_array(slot)%rs_n_global(3))
      !$omp parallel do private(jj,ii)
      do kk = 1, fft_array(slot)%rs_n(3)
        do jj = 1, fft_array(slot)%rs_n(2)
          do ii = 1, fft_array(slot)%rs_n(1)
            out(ii, jj, kk) = out(ii, jj, kk)*scaling_factor
          end do
        end do
      end do
      !$omp end parallel do
    end if

    call fft_operation_count(fft)

    call profiling_out(prof_bw)

    POP_SUB(X(fft_backward_3d))
  end subroutine X(fft_backward_3d)

  ! ---------------------------------------------------------

  subroutine X(fft_backward_accel)(fft, in, out)
    type(fft_t),       intent(in)    :: fft
    type(accel_mem_t), intent(in)    :: in
    type(accel_mem_t), intent(inout) :: out

    integer :: slot
    type(profile_t), save :: prof_bw
#ifdef HAVE_CLFFT
    integer                    :: bsize
    integer(8)                 :: tmp_buf_size
    type(accel_mem_t)         :: tmp_buf
#endif

    PUSH_SUB(X(fft_backward_accel))
    
    call profiling_in(prof_bw,TOSTRING(X(FFT_BACKWARD_ACCEL)))

    slot = fft%slot
    ASSERT(fft_array(slot)%library == FFTLIB_ACCEL)

#ifdef HAVE_CUDA

#ifdef R_TREAL
    call cuda_fft_execute_z2d(fft_array(slot)%cuda_plan_bw, in%mem, out%mem)
#else
    call cuda_fft_execute_z2z_backward(fft_array(slot)%cuda_plan_bw, in%mem, out%mem)
#endif
    call fft_operation_count(fft)
    call accel_finish()

#endif
    
#ifdef HAVE_CLFFT
    call clfftGetTmpBufSize(fft_array(slot)%cl_plan_bw, tmp_buf_size, cl_status)
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftGetTmpBufSize')

    if(tmp_buf_size > 0) then
      call accel_create_buffer(tmp_buf, ACCEL_MEM_READ_WRITE, TYPE_BYTE, int(tmp_buf_size, 4))
    end if

    if(tmp_buf_size > 0) then
      call clfftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, accel%command_queue, &
        in%mem, out%mem, tmp_buf%mem, cl_status)
    else
      call clfftEnqueueTransform(fft_array(slot)%cl_plan_bw, CLFFT_FORWARD, accel%command_queue, &
        in%mem, out%mem, cl_status)
    end if
    
    if(cl_status /= CLFFT_SUCCESS) call clfft_print_error(cl_status, 'clfftEnqueueTransform')
    
    call fft_operation_count(fft)

    call accel_finish()

    if(tmp_buf_size > 0) call accel_release_buffer(tmp_buf)

#endif

    call profiling_out(prof_bw)

    POP_SUB(X(fft_backward_accel))
  end subroutine X(fft_backward_accel)

  ! ---------------------------------------------------------
  subroutine X(fft_backward_1d)(fft, in, out)
    type(fft_t), intent(in)    :: fft
    CMPLX,       intent(inout) :: in(:)
    R_TYPE,      intent(out)   :: out(:)
 
    PUSH_SUB(X(fft_backward_1d))
    
#ifdef R_TREAL
    call fftw_execute_dft_c2r(fft_array(fft%slot)%planb, in, out)
#else
    call fftw_execute_dft(fft_array(fft%slot)%planb, in, out)
#endif

    call fft_operation_count(fft)

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft_array(fft%slot)%rs_n_global(1), R_TOTYPE(M_ONE) / fft_array(fft%slot)%rs_n_global(1), out)

    POP_SUB(X(fft_backward_1d))
  end subroutine X(fft_backward_1d)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
