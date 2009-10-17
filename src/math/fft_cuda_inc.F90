!! Copyright (C) 2008  X. Andrade
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
!! $Id: fft_cuda_inc.F90 3946 2008-03-25 23:09:38Z xavier $

module fft_t_m
  use c_pointer_m

  implicit none  

  private

  public :: fft_t

  type fft_t
    private
    type(c_ptr) :: opaque
  end type fft_t

end module

module fft_m
  use global_m
  use messages_m
  use datasets_m
  use loct_math_m
  use parser_m
  use lalg_basic_m

  use fft_t_m

  implicit none

  private
  public ::           &
       fft_t,         &
       fft_all_init,  &
       fft_all_end,   &
       fft_init,      &
       fft_end,       &
       fft_copy,      &
       pad_feq,       &
       dfft_forward,  &
       zfft_forward,  &
       dfft_backward, &
       zfft_backward 
  
  ! global constants
  integer, public, parameter ::                &
       fft_real    = 0,                        &
       fft_complex = 1
  
  interface
    
    subroutine fft_all_init
    end subroutine fft_all_init
    
    subroutine fft_all_end
    end subroutine fft_all_end
    
    subroutine fft_init(n, is_real, fft, optimize)
      use fft_t_m
      integer,           intent(inout) :: n(MAX_DIM)
      integer,           intent(in)    :: is_real
      type(fft_t),       intent(out)   :: fft
      logical, optional, intent(in)    :: optimize
    end subroutine fft_init
    
    subroutine fft_copy(fft_i, fft_o)
      use fft_t_m
      type(fft_t), intent(in)  :: fft_i
      type(fft_t), intent(out) :: fft_o
    end subroutine fft_copy
    
    subroutine fft_end(fft)
      use fft_t_m
      type(fft_t), intent(inout) :: fft
    end subroutine fft_end

    subroutine fft_getdim_real(fft, d)
      use fft_t_m
      type(fft_t), intent(in) :: fft
      integer,    intent(out) :: d(MAX_DIM)
    end subroutine fft_getdim_real

    subroutine fft_getdim_complex(fft, d)
      use fft_t_m
      type(fft_t), intent(in)  :: fft
      integer,     intent(out) :: d(MAX_DIM)
    end subroutine fft_getdim_complex

    subroutine dfft_forward(fft, r, c)
      use fft_t_m
      type(fft_t), intent(in)  :: fft
      FLOAT,       intent(in)  :: r(:,:,:)
      CMPLX,       intent(out) :: c(:,:,:)
    end subroutine dfft_forward
    
    subroutine dfft_backward(fft, c, r)
      use fft_t_m
      type(fft_t), intent(in) :: fft
      CMPLX, intent(in)  :: c(:, :, :)
      FLOAT, intent(out) :: r(:, :, :)
    end subroutine dfft_backward
    
    subroutine zfft_forward(fft, in, out)
      use fft_t_m
      type(fft_t), intent(in) :: fft
      CMPLX, intent(in)  :: in(:,:,:)
      CMPLX, intent(out) :: out(:,:,:)
    end subroutine zfft_forward
    
    subroutine zfft_backward(fft, in, out)
      use fft_t_m
      type(fft_t), intent(in) :: fft
      CMPLX, intent(in)  ::  in(:, :, :)
      CMPLX, intent(out) :: out(:, :, :)
    end subroutine zfft_backward

    integer function pad_feq(i, n, mode)
      integer, intent(in) :: i,n
      logical, intent(in) :: mode
    end function pad_feq
  end interface

  end module fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
