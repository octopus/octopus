!! Copyright (C) 2008 X. Andrade
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

module gridhier_m
  use curvilinear_m
  use datasets_m
  use derivatives_m
  use geometry_m
  use global_m
  use parser_m
  use mesh_m
  use messages_m
  use multigrid_m
  use profiling_m

  implicit none

  private
  public ::                         &
    dgridhier_t,                    &
    zgridhier_t,                    &
    dmg_pointer_t,                  &
    zmg_pointer_t,                  &
    gridhier_init,                  &
    gridhier_end

  type dmg_pointer_t
    FLOAT, pointer :: p(:)
  end type dmg_pointer_t

  type zmg_pointer_t
    CMPLX, pointer :: p(:)
  end type zmg_pointer_t

  type dgridhier_t
    type(dmg_pointer_t), pointer :: level(:)
  end type dgridhier_t

  type zgridhier_t
    type(zmg_pointer_t), pointer :: level(:)
  end type zgridhier_t

  interface gridhier_init
    module procedure dgridhier_init, zgridhier_init
  end interface
  
  interface gridhier_end
    module procedure dgridhier_end, zgridhier_end
  end interface

contains

#include "undef.F90"
#include "real.F90"
#include "gridhier_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "gridhier_inc.F90"

end module gridhier_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
