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

#ifdef HAVE_FFT
  ! this routine computes
  ! cf_o = cf_o + exp(-k vec) cf_i
  subroutine cf_phase_factor(m, vec, cf_i, cf_o)
    type(mesh_type), intent(in) :: m
    real(r8), intent(in)        :: vec(3)
    type(dcf), intent(in)       :: cf_i
    type(dcf), intent(inout)    :: cf_o
  
    complex(r8) :: k(3)
    integer     :: n(3), ix, iy, iz, ixx, iyy, izz
  
    ASSERT(all(cf_i%n == cf_o%n))
    ASSERT(associated(cf_i%FS).and.associated(cf_o%FS))

    k = M_z0
    k(1:conf%dim) = M_zI * ((2.0_r8*M_Pi)/(cf_i%n(1:conf%dim)*m%h(1:conf%dim)))

    n  = cf_i%n
    do iz = 1, n(3)
      izz = pad_feq(iz, n(3), .true.)
      do iy = 1, n(2)
        iyy = pad_feq(iy, n(2), .true.)
        do ix = 1, cf_i%nx
          ixx = pad_feq(ix, n(1), .true.)
          
          cf_o%FS(ix, iy, iz) = cf_o%FS(ix, iy, iz) + &
               exp( -(k(1)*vec(1)*ixx + k(2)*vec(2)*iyy + k(3)*vec(3)*izz) ) * cf_i%FS(ix, iy, iz)
        end do
      end do
    end do
  end subroutine cf_phase_factor
#endif


#include "undef.F90"
#include "real.F90"
#include "cf_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "cf_inc.F90"

end module cube_function
  
