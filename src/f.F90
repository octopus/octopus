
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

module functions_m
  use global_m
  use messages_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use mesh_m
  use mesh_function_m
  use cube_function_m
  use derivatives_m
  use varinfo_m
#if defined(HAVE_FFT)
  use fft_m
#endif

  implicit none

  private
  public  ::                    &
    f_der_t,                    &
    f_der_init,                 &
    f_der_build,                &
    f_der_end,                  &
    dmf2cf, zmf2cf,             &
    dcf2mf, zcf2mf,             &
    dcf_FS2mf, zcf_FS2mf,       &
    df_laplacian, zf_laplacian, &
    df_gradient, zf_gradient,   &
    df_divergence,              &
    zf_divergence,              &
    df_curl, zf_curl,           &
    df_angular_momentum,        &
    zf_angular_momentum,        &
    df_l2, zf_l2

  integer, public, parameter :: &
    REAL_SPACE = 0,             &
    FOURIER_SPACE = 1

  type f_der_t
    type(mesh_t), pointer :: m            ! a pointer to mesh

    integer               :: space        ! derivatives calculated in real or fourier space

    ! derivatives in real space
    integer               :: n_ghost(3)   ! ghost points to add in each dimension
    type(der_discr_t)     :: der_discr    ! discretization of the derivatives

    ! derivatives in fourier space
#if defined(HAVE_FFT)
    type(dcf_t) :: dcf_der, dcf_aux            ! auxiliary variables
    type(zcf_t) :: zcf_der, zcf_aux            ! derivatives in fourier space
#endif
  end type f_der_t

contains

  ! ---------------------------------------------------------
  subroutine f_der_init(f_der, sb, use_curvlinear)
    type(f_der_t),     intent(out) :: f_der
    type(simul_box_t), intent(in)  :: sb
    logical,           intent(in)  :: use_curvlinear

    call push_sub('f.f_der_init')

#ifdef HAVE_FFT
    !%Variable DerivativesSpace
    !%Type integer
    !%Default real_space
    !%Section Mesh::Derivatives
    !%Description
    !% Defines in which space the gradients and the Laplacian are calculated.
    !%Option real_space 0
    !% Derivatives are calculated in real-space using finite differences.
    !%Option fourier_space 1
    !% Derivatives are calculated in
    !% reciprocal space. Obviously this case implies cyclic boundary conditions,
    !% so be careful.
    !%End
    call loct_parse_int(check_inp('DerivativesSpace'), REAL_SPACE, f_der%space)
    if(.not.varinfo_valid_option('DerivativesSpace', f_der%space)) call input_error('DerivativesSpace')
#else
    f_der%space = REAL_SPACE
#endif
    call messages_print_var_option(stdout, "DerivativesSpace", f_der%space)

    if(f_der%space == REAL_SPACE) then
      call derivatives_init(sb, f_der%der_discr, f_der%n_ghost, use_curvlinear)
#if defined(HAVE_FFT)
    else
      if(use_curvlinear) then
        message(1) = "When using curvilinear coordinates you must use"
        message(2) = "DerivativesSpace = real_space"
        call write_fatal(2)
      end if
#endif
    end if

    call pop_sub()
  end subroutine f_der_init


  ! ---------------------------------------------------------
  subroutine f_der_build(sb, m, f_der)
    type(simul_box_t),    intent(in)    :: sb
    type(mesh_t), target, intent(in)    :: m
    type(f_der_t),        intent(inout) :: f_der

    call push_sub('f.f_der_build')

    f_der%m => m ! keep a working pointer to the underlying mesh

    if(f_der%space == REAL_SPACE) then
      call derivatives_build(m, f_der%der_discr)
#if defined(HAVE_FFT)
    else
      call dcf_new(m%l, f_der%dcf_der)
      call dcf_fft_init(f_der%dcf_der, sb)
      call dcf_new_from(f_der%dcf_aux, f_der%dcf_der)

      call zcf_new(m%l, f_der%zcf_der)
      call zcf_fft_init(f_der%zcf_der, sb)
      call zcf_new_from(f_der%zcf_aux, f_der%zcf_der)
#endif
    end if

    call pop_sub()
  end subroutine f_der_build


  ! ---------------------------------------------------------
  subroutine f_der_end(f_der)
    type(f_der_t), intent(inout) :: f_der

    call push_sub('f.f_der_end')

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

end module functions_m
