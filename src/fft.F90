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

#include "config_F90.h"

#ifdef HAVE_FFTW
  module fft
  use global

  implicit none


! fftw constants
  integer, parameter :: fftw_forward = -1
  integer, parameter :: fftw_backward = 1
  integer, parameter :: fftw_real_to_complex = -1
  integer, parameter :: fftw_complex_to_real = 1
  integer, parameter :: fftw_estimate = 0
  integer, parameter :: fftw_measure = 1
  integer, parameter :: fftw_in_place = 8
  integer, parameter :: fftw_out_of_place = 0
  integer, parameter :: fftw_use_wisdom = 16
  integer, parameter :: fftw_threadsafe = 128


  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes real FFT transforms using the fftw package         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! input:                                                                      !
! n1,n2,n3		grid size                                             !
! re			real array                                            !
! co			complex array (hermitian)                             !
! mode			determines transform mode:                            !
!			mode < 0 : forward (real to complex)                  !
!			mode > 0 : backward (complex to real, with the        !
!					(1/(n1*n2*n3)) factor                 !
!                                                                             !
! output:                                                                     !
! re,co                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rfft3d(re, co, n1, n2, n3, plan, mode)

    integer, intent(in)               :: n1,n2,n3,mode
    integer(POINTER_SIZE), intent(in) :: plan
    real(r8), intent(inout)           :: re(n1,n2,n3)
    complex(r8), intent(inout)        :: co(n1/2+1,n2,n3)
    
    if(mode<0) then
      call rfftwnd_f77_one_real_to_complex(plan,re,co)
    else
      call rfftwnd_f77_one_complex_to_real(plan,co,re)
      call dscal(n1*n2*n3, 1.0_r8/real(n1*n2*n3, r8), re, 1)
    endif
    
    return
  end subroutine rfft3d

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

#endif
