!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module cube_function
  use global
  use mesh
#if defined(HAVE_FFT)
  use fft
#endif

  implicit none

  type dcf
    integer :: n(3)   ! the linear dimensions of the cube

    real(r8),    pointer :: RS(:,:,:)
    complex(r8), pointer :: FS(:,:,:)

#if defined(HAVE_FFT)
    integer :: nx ! = n(1)/2 + 1, first dimension of the FS array
    type(fft_type), pointer :: fft
#endif
  end type dcf

  type zcf
    integer :: n(3)   ! the linear dimensions of the cube

    complex(r8), pointer :: RS(:,:,:)
    complex(r8), pointer :: FS(:,:,:)

#if defined(HAVE_FFT)
    integer :: nx ! = n(1),  first dimension of the FS array
    type(fft_type), pointer :: fft
#endif
  end type zcf

contains

#include "undef.F90"
#include "real.F90"
#include "cf_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "cf_inc.F90"

end module cube_function
  
