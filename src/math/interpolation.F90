!! Copyright (C) 2014 X. Andrade
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Id$
  
#include <global.h>

module interpolation_m
  use global_m
  use messages_m
  use profiling_m
  use qshep_m

  implicit none

  private
  
  public ::                                     &
    interpolation_t,                            &
    interpolation_init,                         &
    interpolation_end,                          &
    interpolation_point,                        &
    interpolation_point_t,                      &
    interpolation_point_end,                    &
    interpolation_point_evaluate


  type interpolation_t
    private

    integer          :: ndim
    integer          :: npoints
    FLOAT, pointer   :: points(:, :)
  end type interpolation_t

  type interpolation_point_t
    private

    type(interpolation_t), pointer :: interpolation
    FLOAT, pointer                 :: position(:)
  end type interpolation_point_t

  interface interpolation_point_evaluate
    module procedure dinterpolation_point_evaluate
    module procedure zinterpolation_point_evaluate
  end interface interpolation_point_evaluate

contains

  subroutine interpolation_init(this, ndim, npoints, points)
    type(interpolation_t), intent(out)   :: this
    integer,              intent(in)     :: ndim
    integer,              intent(in)     :: npoints
    FLOAT,                intent(in)     :: points(:, :)

    PUSH_SUB(interpolation_init)

    this%ndim = ndim
    this%npoints = npoints
    
    SAFE_ALLOCATE(this%points(1:this%npoints, 1:this%ndim))
    
    this%points(1:this%npoints, 1:this%ndim) = points(1:this%npoints, 1:this%ndim)
    
    POP_SUB(interpolation_init)

  end subroutine interpolation_init

  ! -----------------------------------------------

  subroutine interpolation_end(this)
    type(interpolation_t), intent(inout) :: this

    PUSH_SUB(interpolation_end)

    SAFE_DEALLOCATE_P(this%points)

    POP_SUB(interpolation_end)

  end subroutine interpolation_end

  ! -----------------------------------------------

  subroutine interpolation_point(this, position, interp_point)
    type(interpolation_t), target, intent(in)    :: this
    FLOAT,                         intent(in)    :: position(:)
    type(interpolation_point_t),   intent(out)   :: interp_point

    PUSH_SUB(interpolation_point)

    interp_point%interpolation => this

    ASSERT(ubound(position, dim = 1) >= this%ndim)

    SAFE_ALLOCATE(interp_point%position(1:this%ndim))

    interp_point%position(1:this%ndim) = position(1:this%ndim)

    POP_SUB(interpolation_point)

  end subroutine interpolation_point

  ! -----------------------------------------------

  subroutine interpolation_point_end(this)
    type(interpolation_point_t),   intent(out)   :: this
    
    PUSH_SUB(interpolation_point_end)

    SAFE_DEALLOCATE_P(this%position)

    POP_SUB(interpolation_point_end)

  end subroutine interpolation_point_end
  
#include "undef.F90"
#include "real.F90"
#include "interpolation_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "interpolation_inc.F90"

end module interpolation_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
