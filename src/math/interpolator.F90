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

module interpolator_m
  use global_m
  use messages_m
  use profiling_m
  use qshep_m

  implicit none

  private
  
  public ::                                    &
    interpolator_t,                            &
    interpolator_init,                         &
    interpolator_end,                          &
    interpolator_interpolate

  type interpolator_t
    integer          :: ndim
    integer          :: npoints
    FLOAT, pointer   :: points(:, :)
  end type interpolator_t

  interface interpolator_interpolate
    module procedure dinterpolator_interpolate
    module procedure zinterpolator_interpolate
  end interface interpolator_interpolate

contains

  subroutine interpolator_init(this, ndim, npoints, points)
    type(interpolator_t), intent(out)   :: this
    integer,              intent(in)    :: ndim
    integer,              intent(in)    :: npoints
    FLOAT,                intent(in)    :: points(:, :)

    this%ndim = ndim
    this%npoints = npoints
    
    SAFE_ALLOCATE(this%points(1:this%npoints, 1:this%ndim))
    
    this%points(1:this%npoints, 1:this%ndim) = points(1:this%npoints, 1:this%ndim)
    
  end subroutine interpolator_init

  ! -----------------------------------------------

  subroutine interpolator_end(this)
    type(interpolator_t), intent(inout) :: this

    SAFE_DEALLOCATE_P(this%points)

  end subroutine interpolator_end

#include "undef.F90"
#include "real.F90"
#include "interpolator_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "interpolator_inc.F90"

end module interpolator_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
