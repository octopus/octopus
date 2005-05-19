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
!!
!! $Id$

#include "global.h"

module cube_function
  use global
  use mesh
#if defined(HAVE_FFT)
  use fft
#endif

  implicit none

  private
  public :: dcf, zcf, &
            cf_phase_factor, &
            dcf_new, zcf_new, &
            dcf_new_from, zcf_new_from, &
            dcf_alloc_rs, zcf_alloc_rs, &
            dcf_free_rs, zcf_free_rs, &
            dcf_alloc_fs, zcf_alloc_fs, &
            dcf_free_fs, zcf_free_fs, &
            dcf_free, zcf_free, &
            dcf_fft_init, zcf_fft_init, &
            dcf_RS2FS, zcf_RS2FS, &
            dcf_FS2RS, zcf_FS2RS, &
            dcf_FS_lapl, zcf_FS_lapl, &
            dcf_FS_grad, zcf_FS_grad, &
            cf_surface_average

  type dcf
    integer :: n(3)   ! the linear dimensions of the cube

    FLOAT, pointer :: RS(:,:,:)
    CMPLX, pointer :: FS(:,:,:)

#if defined(HAVE_FFT)
    integer :: nx ! = n(1)/2 + 1, first dimension of the FS array
    type(fft_type), pointer :: fft
#endif
  end type dcf

  type zcf
    integer :: n(3)   ! the linear dimensions of the cube

    CMPLX, pointer :: RS(:,:,:)
    CMPLX, pointer :: FS(:,:,:)

#if defined(HAVE_FFT)
    integer :: nx ! = n(1),  first dimension of the FS array
    type(fft_type), pointer :: fft
#endif
  end type zcf

contains

  ! This function calculates the surface average of any function.
  ! WARNING: Some more careful testing should be done on this.
  FLOAT function cf_surface_average(cf) result(x)
    type(dcf), intent(in)       :: cf

    integer ix, iy, iz, npoints
    x = M_ZERO

    do iy = 2, cf%n(2) - 1
       do iz = 2, cf%n(3) - 1
          x = x + (cf%RS(1, iy, iz) + cf%RS(cf%n(1), iy, iz))
       enddo
    enddo

    do ix = 2, cf%n(1) - 1
       do iz = 2, cf%n(3) - 1
          x = x + (cf%RS(ix, 1, iz) + cf%RS(ix, cf%n(2), iz))
       enddo
    enddo

    do ix = 2, cf%n(1) - 1
       do iy = 2, cf%n(2) - 1
          x = x + (cf%RS(ix, iy, 1) + cf%RS(ix, iy, cf%n(3)))
       enddo
    enddo

    do iz = 2, cf%n(3) - 1
       x = x + cf%RS(1, 1, iz) + cf%RS(cf%n(1), 1, iz) + &
               cf%RS(1, cf%n(2), iz) + cf%RS(cf%n(1), cf%n(2), 1) 
    enddo

    do iy = 2, cf%n(2) - 1
       x = x + cf%RS(1, iy, 1) + cf%RS(cf%n(1), iy, 1) + &
               cf%RS(1, iy, cf%n(3)) + cf%RS(cf%n(1), iy, cf%n(3)) 
    enddo

    do ix = 2, cf%n(1) - 1
       x = x + cf%RS(ix, 1, 1) + cf%RS(ix, cf%n(2), 1) + &
               cf%RS(ix, 1, cf%n(3)) + cf%RS(ix, cf%n(2), cf%n(3)) 
    enddo

    x = x + cf%RS(1, 1, 1)             + cf%RS(cf%n(1), 1, 1) + &
            cf%RS(1, cf%n(2), 1)       + cf%RS(cf%n(1), cf%n(2), 1) + &
            cf%RS(1, 1, cf%n(3))       + cf%RS(cf%n(1), 1, cf%n(3)) + &
            cf%RS(1, cf%n(2), cf%n(3)) + cf%RS(cf%n(1), cf%n(2), cf%n(3))

    npoints = 6*(cf%n(1)-2)**2 + 12*(cf%n(1)-2) + 8
    x = x/npoints

  end function cf_surface_average

#ifdef HAVE_FFT
  ! this routine computes
  ! cf_o = cf_o + exp(-k vec) cf_i
  subroutine cf_phase_factor(m, vec, cf_i, cf_o)
    type(mesh_type), intent(in) :: m
    FLOAT, intent(in)        :: vec(3)
    type(dcf), intent(in)       :: cf_i
    type(dcf), intent(inout)    :: cf_o
  
    CMPLX :: k(3)
    integer     :: n(3), ix, iy, iz, ixx, iyy, izz
  
    ASSERT(all(cf_i%n == cf_o%n))
    ASSERT(associated(cf_i%FS).and.associated(cf_o%FS))

    k = M_z0
    k(1:conf%dim) = M_zI * ((M_TWO*M_Pi)/(cf_i%n(1:conf%dim)*m%h(1:conf%dim)))

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
  
