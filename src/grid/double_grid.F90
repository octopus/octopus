!! Copyright (C) 2002-2006 X. Andrade
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
!! $Id: double_grid.F90 2711 2007-02-13 17:36:18Z xavier $

#include "global.h"

module double_grid_m

  use curvlinear_m
  use datasets_m
  use geometry_m
  use global_m
  use math_m
  use mesh_m
  use messages_m
  use mesh_function_m
  use loct_m
  use loct_parser_m
  use par_vec_m
  use profiling_m
  use simul_box_m
  use species_m
  use splines_m
  use submesh_m

#ifdef USE_OMP
  use omp_lib
#endif

  implicit none
  
  private
  
  public ::              &
       double_grid_t,    &
       double_grid_init, &
       double_grid_end,  &
       double_grid_apply_local,     &
       double_grid_apply_non_local, &
       double_grid_get_rmax,        &
       double_grid_get_hmax,        &
       double_grid_enlarge

  type double_grid_t
     private
     integer :: order
     integer :: npoints
     integer :: spacing_divisor
     integer :: interpolation_min
     integer :: interpolation_max
     integer :: nn
     logical :: use_double_grid
     FLOAT, pointer :: co(:)
  end type double_grid_t

  type(profile_t), save :: double_grid_local_prof, double_grid_nonlocal_prof

contains
  
  subroutine double_grid_init(this, sb)
    type(double_grid_t), intent(out) :: this
    type(simul_box_t),   intent(in)  :: sb

    call push_sub('double_grid.double_grid_init')

    this%spacing_divisor = 3
    this%nn = (this%spacing_divisor - 1)/2
   
    !%Variable DoubleGrid
    !%Type logical
    !%Default no
    !%Section Mesh
    !%Description
    !% Enables or disables the use of a double-grid technique to
    !% increase the precision of the application of the
    !% pseudopotentials.
    !%End
    if (sb%dim == 3) then 
      call loct_parse_logical(datasets_check('DoubleGrid'), .false., this%use_double_grid)
    else
      this%use_double_grid = .false.
    end if

    if(this%use_double_grid) call messages_devel_version('Double grid')

    !%Variable DoubleGridOrder
    !%Type integer 
    !%Default 9
    !%Section Mesh
    !%Description
    !% Order of the interpolation used for the double grid. Must be
    !% an odd number. Low-order interpolation schemes are not
    !% recommended. The default is to use 9th-order interpolation.
    !%End
    call loct_parse_int(datasets_check('DoubleGridOrder'), 9, this%order)
    
    ASSERT(mod(this%order,2) == 1)
    
    this%interpolation_min = -this%order/2
    this%interpolation_max =  this%order/2 + 1
    this%npoints = this%interpolation_max - this%interpolation_min + 1
    
    ALLOCATE(this%co(this%interpolation_min:this%interpolation_max), this%npoints)

    call calc_coefficients

    call pop_sub()

    contains

      subroutine calc_coefficients
        FLOAT, allocatable :: points(:)
        integer :: ii

        call push_sub('double_grid.double_grid_init.calc_coefficients')

        ALLOCATE(points(this%interpolation_min:this%interpolation_max), this%npoints)
        
        do ii = this%interpolation_min,  this%interpolation_max
          points(ii) = M_THREE*ii - M_ONE
        end do

        call interpolation_coefficients(this%npoints, points, M_ZERO, this%co)

        deallocate(points)
        call pop_sub()
        
      end subroutine calc_coefficients

  end subroutine double_grid_init

  subroutine double_grid_end(this)
    type(double_grid_t), intent(inout) :: this
 
    call push_sub('double_grid.double_grid_end')
    deallocate(this%co)
    call pop_sub()

  end subroutine double_grid_end
  
  FLOAT function double_grid_get_hmax(this, mesh) result(hmax)
    type(double_grid_t), intent(in) :: this
    type(mesh_t),        intent(in) :: mesh

    call push_sub('double_grid.double_grid_get_hmax')

    hmax = maxval(mesh%h(1:MAX_DIM))
      
    if(this%use_double_grid)  hmax = hmax / this%spacing_divisor    

    call pop_sub()
  end function double_grid_get_hmax

  FLOAT function double_grid_get_rmax(this, s, m) result(rmax)
    type(double_grid_t),     intent(in) :: this
    type(species_t),         intent(in) :: s
    type(mesh_t),            intent(in) :: m
    
    call push_sub('double_grid.double_grid_get_rmax')

    rmax = spline_cutoff_radius(s%ps%vl, s%ps%projectors_sphere_threshold)
    if(this%use_double_grid) then 
      rmax = rmax + this%interpolation_max * maxval(m%h(1:3))
    end if

    call pop_sub()
  end function double_grid_get_rmax

  integer function double_grid_enlarge(this)
    type(double_grid_t),     intent(in) :: this
    
    call push_sub('double_grid.double_grid_enlarge')

    double_grid_enlarge = 0 
    if(this%use_double_grid) double_grid_enlarge = this%interpolation_max * this%nn

    call pop_sub()
  end function double_grid_enlarge

#define profiler double_grid_local_prof
#define profiler_label "DOUBLE_GRID_LOCAL"
#define double_grid_apply double_grid_apply_local
#define calc_pot(vv) vv = spline_eval(s%ps%vl, r)

#include "double_grid_apply.F90"

#undef calc_pot
#undef double_grid_apply
#undef profiler
#undef profiler_label

#define profiler double_grid_nonlocal_prof
#define profiler_label "DOUBLE_GRID_NL"
#define double_grid_apply double_grid_apply_non_local
#define calc_pot(vv) call species_real_nl_projector(s, x, l, lm, ic, vv, tmp)

#include "double_grid_apply.F90"

end module double_grid_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
