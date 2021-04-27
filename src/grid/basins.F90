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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module basins_oct_m
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use par_vec_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: &
     basins_t,       &
     basins_init,    &
     basins_end,     &
     basins_analyze, &
     basins_write

  type basins_t
    private
    integer, allocatable, public :: map(:)

    integer              :: number
    integer, allocatable :: position(:)
    FLOAT,   allocatable :: val(:)
    FLOAT,   allocatable :: volume(:)

    FLOAT,   allocatable :: population(:)
  end type basins_t

contains

  !----------------------------------------------------------------
  subroutine basins_init(this, mesh)
    type(basins_t), intent(out) :: this
    type(mesh_t),   intent(in)  :: mesh
    
    PUSH_SUB(basins_init)

    if(mesh%parallel_in_domains) &
      call messages_experimental("Bader basins parallel in domains")
    
    SAFE_ALLOCATE(this%map(1:mesh%np))
    this%map(1:mesh%np) = -1

    POP_SUB(basins_init)
  end subroutine basins_init


  !----------------------------------------------------------------
  subroutine basins_end(this)
    type(basins_t), intent(inout) :: this

    PUSH_SUB(basins_end)

    ASSERT(allocated(this%map))
    SAFE_DEALLOCATE_A(this%map)

    SAFE_DEALLOCATE_A(this%position)
    SAFE_DEALLOCATE_A(this%val)
    SAFE_DEALLOCATE_A(this%volume)
    SAFE_DEALLOCATE_A(this%population)
    
    POP_SUB(basins_end)
  end subroutine basins_end


  !----------------------------------------------------------------
  subroutine basins_analyze(this, mesh, f, rho, threshold)
    type(basins_t), intent(inout) :: this
    type(mesh_t),   intent(in)    :: mesh
    FLOAT,          intent(in)    :: f(:)
    FLOAT,          intent(in)    :: rho(:, :)
    FLOAT,          intent(in)    :: threshold

    integer :: jj, xmax, ymax, zmax
    integer :: cur_color, dum

    PUSH_SUB(basins_analyze)

    ASSERT(allocated(this%map))

    cur_color = 0

    xmax = 1; ymax = 1; zmax = 1
    if(mesh%sb%dim < 3) zmax = 0
    if(mesh%sb%dim < 2) ymax = 0

    do jj = 1, mesh%np
      if(this%map(jj)  ==  -1) then
        dum = steep_fill(jj)
      end if
    end do

    call analyze()

    POP_SUB(basins_analyze)

  contains

    !----------------------------------------------------------------
    recursive integer function steep_fill(ii) result(color)
      integer, intent(in) :: ii

      integer :: ii_max

      ! No PUSH/POP_SUB in a recursive routine, as the recursion limit of 50 can be reached.
      
      if(this%map(ii) >= 0) then
        color = this%map(ii)
        return
      end if
      this%map(ii) = -2

      ii_max = get_max(ii, M_ZERO)
      if(ii_max == -1) then ! this is required to get ring attractors
        ii_max = get_max(ii, threshold)
      end if

      if(ii_max /= -1) then
        color = steep_fill(ii_max)
      else
        color = cur_color
        cur_color = cur_color + 1
      end if
      this%map(ii) = color

    end function steep_fill


    !----------------------------------------------------------------
    integer function get_max(ii, threshold)
      integer, intent(in) :: ii
      FLOAT,   intent(in) :: threshold

      FLOAT   :: f_max
      integer :: xx, yy, zz, index
      integer :: point(MAX_DIM), point2(MAX_DIM)

      PUSH_SUB(basins_analyze.get_max)

      point = 0
      call mesh_local_index_to_coords(mesh, ii, point)

      f_max   = f(ii)
      get_max = -1
      do xx = -xmax, xmax
        do yy = -ymax, ymax
          do zz = -zmax, zmax
            if(xx==0.and.yy==0.and.zz==0) cycle
              
            point2(:) = point(:)
            point2(1) = point2(1) + xx
            point2(2) = point2(2) + yy
            point2(3) = point2(3) + zz

            index = mesh_local_index_from_coords(mesh, point2)

            if(index <= 0 .or. index > mesh%np) cycle
            if(this%map(index) == -2) cycle

            if(f_max <= f(index) + threshold) then
              f_max   = f(index)
              get_max = index
            end if

          end do
        end do
      end do     

      POP_SUB(basins_analyze.get_max)
    end function get_max


    !----------------------------------------------------------------
    subroutine analyze()
      integer :: ii, jj, ii_max
      FLOAT :: f_max

      PUSH_SUB(basins_analyze.analyze)

      this%number = maxval(this%map) + 1
      if(this%number <= 0) then
        message(1) = "Internal error analysing basins of attraction"
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(this%position  (1:this%number))
      SAFE_ALLOCATE(this%val       (1:this%number))
      SAFE_ALLOCATE(this%volume    (1:this%number))
      SAFE_ALLOCATE(this%population(1:this%number))

      this%position(:) = -1
      this%val(:)    = M_ZERO
      this%volume(:)   = M_ZERO
      this%population(:)   = M_ZERO
      do ii = 1, this%number
        ii_max = -1
        f_max  = -huge(f_max)

        do jj = 1, mesh%np
          if(this%map(jj) /= ii-1) cycle
          if(f_max <= f(jj)) then
            ii_max = jj
            f_max  = f(jj)
          end if

          if(mesh%use_curvilinear) then
            this%volume(ii) = this%volume(ii) + mesh%vol_pp(jj)
            this%population(ii) = this%population(ii) + mesh%vol_pp(jj)*sum(rho(jj, :))
          else
            this%volume(ii) = this%volume(ii) + mesh%volume_element
            this%population(ii) = this%population(ii) + mesh%volume_element*sum(rho(jj, :))
          end if
        end do
        
        this%position(ii) = ii_max
        this%val(ii)      = f_max
      end do

      POP_SUB(basins_analyze.analyze)
    end subroutine analyze

  end subroutine basins_analyze


  !----------------------------------------------------------------
  subroutine basins_write(this, mesh, iunit)
    type(basins_t), intent(in) :: this
    type(mesh_t),   intent(in) :: mesh
    integer,        intent(in) :: iunit

    integer :: ii
    type(unit_t) :: unit_vol
    FLOAT :: xx(1:mesh%sb%dim)

    PUSH_SUB(basins_write)

    unit_vol = units_out%length**mesh%sb%dim

    write(iunit, '(a,i5)') 'Number of basins = ', this%number
    write(iunit, '(1x)')

    do ii = 1, this%number
      write(iunit, '(a,i5)') '# ', ii
      xx = units_from_atomic(units_out%length, mesh_x_global(mesh, this%position(ii)))
      write(iunit, '(a,3(f12.6,a), a)') '  position = (', &
         xx(1), ',', xx(2), ',', xx(3), ') ', units_abbrev(units_out%length)
      write(iunit, '(a,f12.6)') '  value = ', this%val(ii)
      write(iunit, '(a,f12.6,a,a)') '  volume = ', units_from_atomic(unit_vol, this%volume(ii)), ' ', units_abbrev(unit_vol)
      write(iunit, '(a,f12.6)') '  population = ', this%population(ii)
    end do

    POP_SUB(basins_write)
  end subroutine basins_write

end module basins_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
