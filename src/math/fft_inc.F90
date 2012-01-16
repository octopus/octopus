!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011 J. Alberdi, P. Garcia Risueño, M. Oliveira
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
!! $Id: fftw3.F90 8581 2011-11-17 15:47:34Z micael $


! ---------------------------------------------------------
subroutine X(fft_forward)(fft, in, out)
    type(fft_t),     intent(in)  :: fft
    R_TYPE,          intent(in)  :: in(:,:,:)
    CMPLX,           intent(out) :: out(:,:,:)

    integer :: ii, jj, kk, slot, n1, n2, n3
    type(profile_t), save :: prof_fw

    PUSH_SUB(X(fft_forward))

    call profiling_in(prof_fw, "FFT_FW")

    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
      ii = min(1, fft_array(slot)%rs_n(1))
      jj = min(1, fft_array(slot)%rs_n(2))
      kk = min(1, fft_array(slot)%rs_n(3))
      call fftw_execute_dft(fft_array(slot)%planf, in(ii,jj,kk), out(ii,jj,kk))
    case (FFTLIB_PFFT)
      ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == in(1,1,1))
      ASSERT(fft_array(slot)%fs_data(1,1,1) == out(1,1,1))
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
      call pfft_execute(fft_array(slot)%pfft_planf)
#endif
    end select

    call profiling_out(prof_fw)

    POP_SUB(X(fft_forward))
  end subroutine X(fft_forward)

  ! ---------------------------------------------------------
  subroutine X(fft_forward1)(fft, in, out)
    type(fft_t), intent(in)  :: fft
    R_TYPE,      intent(in)  :: in(:)
    CMPLX,       intent(out) :: out(:)

    PUSH_SUB(X(fft_forward1))

    call fftw_execute_dft(fft_array(fft%slot)%planf, in(1), out(1))

    POP_SUB(X(fft_forward1))
  end subroutine X(fft_forward1)

  ! ---------------------------------------------------------
  subroutine X(fft_backward)(fft, in, out)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in(:,:,:)
    R_TYPE,      intent(out) :: out(:,:,:)

    integer :: ii, jj, kk, scaling_factor, slot
    type(profile_t), save :: prof_bw

    PUSH_SUB(X(fft_backward))
    
    call profiling_in(prof_bw,"FFT_BW")

    slot = fft%slot
    select case (fft_array(slot)%library)
    case (FFTLIB_FFTW)
      call fftw_execute_dft(fft_array(slot)%planb, in(1,1,1), out(1,1,1))
    case (FFTLIB_PFFT)
      ASSERT(fft_array(slot)%fs_data(1,1,1) == in(1,1,1))
      ASSERT(fft_array(slot)%X(rs_data)(1,1,1) == out(1,1,1))
#ifdef HAVE_PFFT
      call pfft_execute(fft_array(slot)%pfft_planb)
#endif
    end select

    ! multiply by 1/(N1*N2*N2)
    scaling_factor = fft_array(slot)%rs_n_global(1)*fft_array(slot)%rs_n_global(2)*fft_array(slot)%rs_n_global(3)
    !$omp parallel do
    do kk = 1, fft_array(slot)%rs_n(3)
      do jj = 1, fft_array(slot)%rs_n(2)
        do ii = 1, fft_array(slot)%rs_n(1)
          out(ii, jj, kk) = out(ii, jj, kk)/scaling_factor
        end do
      end do
    end do
    !$omp end parallel do

    call profiling_out(prof_bw)

    POP_SUB(X(fft_backward))
  end subroutine X(fft_backward)

  ! ---------------------------------------------------------
  subroutine X(fft_backward1)(fft, in, out)
    type(fft_t), intent(in)  :: fft
    CMPLX,       intent(in)  :: in(:)
    R_TYPE,      intent(out) :: out(:)
 
    PUSH_SUB(X(fft_backward1))
    
    call fftw_execute_dft(fft_array(fft%slot)%planb, in(1), out(1))

    ! multiply by 1/(N1*N2*N2)
    call lalg_scal(fft_array(fft%slot)%rs_n_global(1), R_TOTYPE(M_ONE) / fft_array(fft%slot)%rs_n_global(1), out)

    POP_SUB(X(fft_backward1))
  end subroutine X(fft_backward1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
