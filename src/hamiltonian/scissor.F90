!! Copyright (C) 2009 X. Andrade
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
!! $Id: scissor.F90 4866 2009-01-14 23:21:40Z xavier $

#include "global.h"

module scissor_m
  use global_m
  use lalg_basic_m
  use mesh_m
  use mesh_function_m
  use states_m

  implicit none 

  private

  public ::                  &
       scissor_t,            &
       scissor_nullify,      &
       scissor_init,         &
       dscissor_apply,       &
       zscissor_apply,       &
       scissor_end

  type scissor_t
    logical                 :: apply
    FLOAT                   :: gap
    FLOAT, pointer          :: dpsi(:, :, :, :)
    CMPLX, pointer          :: zpsi(:, :, :, :)
    type(states_t), pointer :: st            ! only for the dimension components
  end type scissor_t
  
  interface scissor_init
    module procedure dscissor_init, zscissor_init
  end interface

contains
  subroutine scissor_nullify(this)
    type(scissor_t), intent(out) :: this
    
    this%apply = .false.
  end subroutine scissor_nullify

  subroutine scissor_end(this)
    type(scissor_t), intent(inout) :: this
    
    this%apply = .false.
    nullify(this%dpsi, this%zpsi)
  end subroutine scissor_end

#include "undef.F90"
#include "real.F90"
#include "scissor_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "scissor_inc.F90"

end module scissor_m
