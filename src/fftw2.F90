!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module fft
  use global
  use linalg
  
  implicit none

  ! global constants
  integer, parameter :: &
       fft_real    = 0, &
       fft_complex = 1

  ! fftw constants WARNING -> they should be private
  integer, parameter :: &
       fftw_forward         = -1, &
       fftw_backward        =  1, &
       fftw_real_to_complex = -1, &
       fftw_complex_to_real =  1, &
       fftw_estimate        =  0, &
       fftw_measure         =  1, &
       fftw_in_place        =  8, &
       fftw_out_of_place    =  0, &
       fftw_use_wisdom      = 16, &
       fftw_threadsafe      =128

  type fft_type
    integer :: slot                ! in which slot do we have this fft

    integer :: n(3)                ! size of the fft
    integer :: is_real             ! is the fft real or complex
    integer(POINTER_SIZE) :: planf ! the plan for forward transforms
    integer(POINTER_SIZE) :: planb ! the plan for backward transforms
  end type fft_type

  integer, private :: fft_refs(FFT_MAX)
  type(fft_type), private :: fft_array(FFT_MAX)
  logical, private :: fft_optimize

contains
  ! initialize the table
  subroutine fft_all_init()
    integer :: i

    call oct_parse_logical("FFTOptimize", .true., fft_optimize)
    
    do i = 1, FFT_MAX
      fft_refs(i) = NULL
    end do
  end subroutine fft_all_init

  ! delete all plans
  subroutine fft_all_end()
    integer :: i

    do i = 1, FFT_MAX
      if(fft_refs(i).ne.NULL) then
        call fftw_f77_destroy_plan(fft_array(i)%planf)
        call fftw_f77_destroy_plan(fft_array(i)%planb)
        fft_refs(i) = NULL
      end if
    end do
  end subroutine fft_all_end

  subroutine fft_init(n, is_real, fft)
    integer, intent(inout) :: n(3)
    integer, intent(in) :: is_real
    type(fft_type), intent(out) :: fft

    integer :: i, j

    ! optimize dimensions in non-periodic directions
    do i = conf%periodic_dim+1, conf%dim 
      if(n(i).ne.1 .and. fft_optimize) &
           call oct_fft_optimize(n(i), 7, 1) ! always ask for an odd number
    end do    

    ! find out if fft has already been allocated
    j = 0
    do i = FFT_MAX, 1, -1
      if(fft_refs(i).ne.NULL) then
        if(n(1) == fft_array(i)%n(1).and.n(2) == fft_array(i)%n(2).and.n(3) == fft_array(i)%n(3).and. &
             is_real == fft_array(i)%is_real) then
          fft = fft_array(i)             ! return a copy
          fft_refs(i) = fft_refs(i) + 1  ! increment the ref count
          return
        end if
      else
        j = i
      end if
    end do

    if(j == 0) then
      message(1) = "Not enough slots for FFTs"
      message(2) = "Please increase FFT_MAX in fft.F90 and recompile"
      call write_fatal(2)
    end if

    ! j now contains an empty slot
    fft_refs(j)          = 1
    fft_array(j)%slot    = j
    fft_array(j)%n       = n
    fft_array(j)%is_real = is_real
    if(is_real == fft_real) then
      call rfftwnd_f77_create_plan(fft_array(j)%planf, conf%dim, n, &
           fftw_real_to_complex + fftw_forward,  fftw_measure + fftw_threadsafe)
      call rfftwnd_f77_create_plan(fft_array(j)%planb, conf%dim, n, &
           fftw_complex_to_real + fftw_backward, fftw_measure + fftw_threadsafe)
    else
      call  fftwnd_f77_create_plan(fft_array(j)%planf, conf%dim, n, &
           fftw_forward,  fftw_measure + fftw_threadsafe)
      call  fftwnd_f77_create_plan(fft_array(j)%planb, conf%dim, n, &
           fftw_backward, fftw_measure + fftw_threadsafe)
    end if
    fft = fft_array(j)

    write(message(1), '(a,i4,a,i4,a,i4,a,i2)') "Info: FFT allocated with size (", &
         n(1), ",", n(2), ",", n(3), ") in slot ", j
    call write_info(1)

  end subroutine fft_init

  subroutine fft_copy(fft_i, fft_o)
    type(fft_type), intent( in) :: fft_i
    type(fft_type), intent(out) :: fft_o
    
    ASSERT(fft_i%slot>=1.and.fft_i%slot<=FFT_MAX)
    ASSERT(fft_refs(fft_i%slot) > 0)

    fft_o = fft_i
    fft_refs(fft_i%slot) = fft_refs(fft_i%slot) + 1
  end subroutine fft_copy

  subroutine fft_end(fft)
    type(fft_type), intent(inout) :: fft

    integer :: i

    i = fft%slot
    if(fft_refs(i)==NULL) then
      message(1) = "Trying to deallocate fft that has not been allocated"
      call write_warning(1)
    else
      if(fft_refs(i) > 1) then
        fft_refs(i) = fft_refs(i) - 1
      else
        fft_refs(i) = NULL
        call fftw_f77_destroy_plan(fft_array(i)%planf)
        call fftw_f77_destroy_plan(fft_array(i)%planb)

        write(message(1), '(a,i4,a,i4,a,i4,a,i2)') "Info: FFT deallocated from slot ", i
        call write_info(1)
      end if
    end if

  end subroutine fft_end

  subroutine fft_getdim_real(fft, d)
    type(fft_type), intent(in) :: fft
    integer, intent(out) :: d(3)

    d = fft%n
  end subroutine fft_getdim_real

  subroutine fft_getdim_complex(fft, d)
    type(fft_type), intent(in) :: fft
    integer, intent(out) :: d(3)

    d = fft%n
    if(fft%is_real == fft_real)  d(1) = d(1)/2 + 1
  end subroutine fft_getdim_complex

  ! these routines simply call fftw
  ! first the real to complex versions
  subroutine dfft_forward(fft, r, c)
    type(fft_type), intent(in) :: fft
    real(r8), intent(in)     :: r(:,:,:)
    complex(r8), intent(out) :: c(:,:,:)

    call rfftwnd_f77_one_real_to_complex(fft%planf, r, c)
  end subroutine dfft_forward

  subroutine dfft_backward(fft, c, r)
    type(fft_type), intent(in) :: fft
    complex(r8), intent(in) :: c(:,:,:)
    real(r8), intent(out)   :: r(:,:,:)

    integer :: n

    call rfftwnd_f77_one_complex_to_real(fft%planb, c, r)

    ! multiply by 1/(N1*N2*N2)
    n = fft%n(1)*fft%n(2)*fft%n(3)
    call dscal(n, M_ONE/real(n, r8), r(1, 1, 1), 1)

  end subroutine dfft_backward

  ! first the complex versions
  subroutine zfft_forward(fft, r, c)
    type(fft_type), intent(in) :: fft
    complex(r8), intent(in)  :: r(:,:,:)
    complex(r8), intent(out) :: c(:,:,:)

    call fftwnd_f77_one(fft%planf, r, c)
  end subroutine zfft_forward

  subroutine zfft_backward(fft, c, r)
    type(fft_type), intent(in) :: fft
    complex(r8), intent(in)  :: c(:,:,:)
    complex(r8), intent(out) :: r(:,:,:)

    integer :: n

    call fftwnd_f77_one(fft%planb, c, r)

    ! multiply by 1/(N1*N2*N2)
    n = fft%n(1)*fft%n(2)*fft%n(3)
    call zscal(n, M_z1/real(n, r8), r(1, 1, 1), 1)

  end subroutine zfft_backward

  ! convert between array index and G vector
  function pad_feq(i, n, mode)
    integer, intent(in) :: i,n
    logical, intent(in) :: mode
    integer :: pad_feq 
    
    if(mode) then      ! index to frequency number
      if( i <= n/2 + 1 ) then
        pad_feq = i - 1
      else
        pad_feq = i - n -1
      endif
    else
      if( i >= 0 ) then
        pad_feq = i + 1
      else
        pad_feq = i + n + 1
      endif
    endif
    
    return
  end function pad_feq

end module fft
