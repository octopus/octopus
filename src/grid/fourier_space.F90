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

module fourier_space_m
  use cube_function_m
  use global_m
  use mesh_m
  use messages_m
  use fft_m
  use simul_box_m

  implicit none
  private
  public ::                     &
    dcf_alloc_rs,               &
    zcf_alloc_rs,               &
    dcf_alloc_fs,               & 
    zcf_alloc_fs,               &
    dcf_free_fs,                &
    zcf_free_fs,                &
    dcf_fft_init,               &
    zcf_fft_init,               &
    dcf_RS2FS,                  &
    zcf_RS2FS,                  &
    dcf_FS2RS,                  &
    zcf_FS2RS

contains

#include "undef.F90"
#include "real.F90"
#include "fourier_space_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "fourier_space_inc.F90"

end module fourier_space_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
