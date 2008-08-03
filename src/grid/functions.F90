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
    dmf2cf, zmf2cf,             &
    dcf2mf, zcf2mf,             &
    df_multipoles,              &
    zf_multipoles,              &
    df_angular_momentum,        &
    zf_angular_momentum,        &
    df_l2, zf_l2,               &
    dcf_FS2mf, zcf_FS2mf

contains

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
