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
       FOURIER_SPACE = 1

  type f_der_type
    type(mesh_type), pointer :: m            ! a pointer to mesh

    integer                  :: space        ! derivatives calculated in real or fourier space

    ! derivatives in real space
    integer                  :: n_ghost(3)   ! ghost points to add in each dimension
    type(der_discr_type)     :: der_discr    ! discretization of the derivatives

    ! derivatives in fourier space
#if defined(HAVE_FFT)
    type(dcf) :: dcf_der, dcf_aux            ! auxiliary variables
    type(zcf) :: zcf_der, zcf_aux            ! derivatives in fourier space
#endif
  end type f_der_type

contains

  ! ---------------------------------------------------------
  subroutine f_der_init(m, f_der)
    type(mesh_type),     pointer     :: m
    type(f_der_type), intent(out) :: f_der

    integer :: norder, j

    call push_sub('f_der_init')

    f_der%m => m ! keep a working pointer to the underlying mesh

#ifdef HAVE_FFT
    call loct_parse_int('DerivativesSpace', REAL_SPACE, f_der%space)
    if((f_der%space.ne.REAL_SPACE).and.(f_der%space.ne.FOURIER_SPACE)) then
      write(message(1), '(a,i5,a)') "Input: '", f_der%space, &
         "' is not a valid DerivativesSpace"
      message(2) = '(DerivativesSpace = real_space | fourier_space)'
      call write_fatal(2)
    end if
#else
    f_der%space = REAL_SPACE
#endif
    
    if(f_der%space == REAL_SPACE) then
      call derivatives_init(m, f_der%der_discr, f_der%n_ghost)
      message(1) = 'Info: Derivatives calculated in real-space'
#if defined(HAVE_FFT)
    else
      if(f_der%m%use_curvlinear) then
        message(1) = "When using curvilinear coordinates you must use"
        message(2) = "DerivativesSpace = real_space"
        call write_fatal(2)
      end if

      message(1) = 'Info: Derivatives calculated in reciprocal-space'
#endif
    end if

    call write_info(1)
    
    call pop_sub()
  end subroutine f_der_init


  ! ---------------------------------------------------------
  subroutine f_der_build(f_der)
    type(f_der_type), intent(inout) :: f_der
    
    call push_sub('f_der_build')

    if(f_der%space == REAL_SPACE) then
      call derivatives_build(f_der%der_discr)
#if defined(HAVE_FFT)
    else
      call dcf_new(f_der%m%l, f_der%dcf_der)
      call dcf_fft_init(f_der%dcf_der)
      call dcf_new_from(f_der%dcf_aux, f_der%dcf_der)
      
      call zcf_new(f_der%m%l, f_der%zcf_der)
      call zcf_fft_init(f_der%zcf_der)
      call zcf_new_from(f_der%zcf_aux, f_der%zcf_der)
#endif
    end if

    call pop_sub()
  end subroutine f_der_build


  ! ---------------------------------------------------------
  subroutine f_der_end(f_der)
    type(f_der_type), intent(inout) :: f_der

    call push_sub('f_der_end')

    ASSERT(associated(f_der%m))
    ASSERT(f_der%space==REAL_SPACE.or.f_der%space==FOURIER_SPACE)

    if(f_der%space == REAL_SPACE) then
      call derivatives_end(f_der%der_discr)
#if defined(HAVE_FFT)
    else
      call dcf_free(f_der%dcf_der)
      call zcf_free(f_der%zcf_der)
      call dcf_free(f_der%dcf_aux)
      call zcf_free(f_der%zcf_aux)
#endif
    end if

    nullify(f_der%m)

    call pop_sub()
  end subroutine f_der_end

#include "undef.F90"
#include "real.F90"
#include "f_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "f_inc.F90"

end module functions
