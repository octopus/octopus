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

module functions
  use global
  use lib_oct_parser
  use lib_basic_alg
  use mesh_function
  use cube_function
  use derivatives
#if defined(HAVE_FFT)
  use fft
#endif

  implicit none

  integer, parameter ::     &
       REAL_SPACE = 0,      &
       RECIPROCAL_SPACE = 1

  integer, private :: derivatives_space

  type(der_discr_type) :: f_der

#if defined(HAVE_FFT)
  type(dcf), private :: dcf_der, dcf_aux  ! these auxiliary variables are used to calculate
  type(zcf), private :: zcf_der, zcf_aux  ! derivatives in fourier space
#endif

contains

  subroutine functions_init(m)
    type(mesh_type), intent(inout) :: m
    integer :: norder, j

    call push_sub('functions_init')

#ifdef HAVE_FFT
  call loct_parse_int('DerivativesSpace', REAL_SPACE, derivatives_space)
  if(derivatives_space < 0 .or. derivatives_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", derivatives_space, &
         "' is not a valid DerivativesSpace"
    message(2) = '(0 <= DerivativesSpace <=1)'
    call write_fatal(2)
  end if
#else
  derivatives_space = REAL_SPACE
#endif

  if(derivatives_space == REAL_SPACE) then
    call loct_parse_int('OrderDerivatives', 4, norder)
    m%der_order   = norder

    call derivatives_init(f_der)
    message(1) = 'Info: Derivatives calculated in real-space'
#if defined(HAVE_FFT)
  else
    call dcf_new(m%l, dcf_der)
    call dcf_fft_init(dcf_der)
    call dcf_new_from(dcf_aux, dcf_der)

    call zcf_new(m%l, zcf_der)
    call zcf_fft_init(zcf_der)
    call zcf_new_from(zcf_aux, zcf_der)
    
    message(1) = 'Info: Derivatives calculated in reciprocal-space'
#endif
  end if

  call write_info(1)

  call pop_sub()
end subroutine functions_init

subroutine functions_end(m)
  type(mesh_type), intent(inout) :: m
  integer :: j
  if(derivatives_space == REAL_SPACE) then
    call mf_end_der(m)
#if defined(HAVE_FFT)
  else
    call dcf_free(dcf_der)
    call zcf_free(zcf_der)
    call dcf_free(dcf_aux)
    call zcf_free(zcf_aux)
#endif
  end if
end subroutine functions_end

#include "undef.F90"
#include "real.F90"
#include "f_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "f_inc.F90"

end module functions
