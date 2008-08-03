!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
  use cube_function_m
  use datasets_m
  use derivatives_m
  use global_m
  use lalg_basic_m
  use loct_math_m
  use loct_parser_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use simul_box_m
  use varinfo_m
  use fft_m

  implicit none

  private
  public  ::                    &
    f_der_t,                    &
    f_der_init,                 &
    f_der_build,                &
    f_der_end,                  &
    dmf2cf, zmf2cf,             &
    dcf2mf, zcf2mf,             &
    df_gradient, zf_gradient,   &
    df_divergence,              &
    zf_divergence,              &
    df_multipoles,              &
    zf_multipoles,              &
    df_curl, zf_curl,           &
    df_angular_momentum,        &
    zf_angular_momentum,        &
    df_l2, zf_l2

  public  ::                    &
    dcf_FS2mf, zcf_FS2mf

  type f_der_t
    type(mesh_t), pointer :: m            ! a pointer to mesh

    ! derivatives in real space
    integer               :: n_ghost(3)   ! ghost points to add in each dimension
    type(der_discr_t)     :: der_discr    ! discretization of the derivatives
  end type f_der_t

contains

  ! ---------------------------------------------------------
  subroutine f_der_init(f_der, sb, use_curvlinear)
    type(f_der_t),     intent(out) :: f_der
    type(simul_box_t), intent(in)  :: sb
    logical,           intent(in)  :: use_curvlinear

    call derivatives_init(sb, f_der%der_discr, f_der%n_ghost, use_curvlinear)

  end subroutine f_der_init


  ! ---------------------------------------------------------
  subroutine f_der_build(sb, m, f_der)
    type(simul_box_t),    intent(in)    :: sb
    type(mesh_t), target, intent(in)    :: m
    type(f_der_t),        intent(inout) :: f_der

    call push_sub('f.f_der_build')

    f_der%m => m ! keep a working pointer to the underlying mesh
    call derivatives_build(m, f_der%der_discr)

    call pop_sub()
  end subroutine f_der_build


  ! ---------------------------------------------------------
  subroutine f_der_end(f_der)
    type(f_der_t), intent(inout) :: f_der

    call push_sub('f.f_der_end')

    ASSERT(associated(f_der%m))

     call derivatives_end(f_der%der_discr)

    call pop_sub()
  end subroutine f_der_end

#include "undef.F90"
#include "real.F90"
#include "functions_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "functions_inc.F90"

end module functions_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
