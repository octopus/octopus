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
!! $Id: mesh.F90 2781 2007-03-23 10:58:32Z lorenzen $

#include "global.h"

module basins_m
  use global_m
  use mesh_m
  use index_m
  use messages_m
  use par_vec_m
  use profiling_m

  implicit none

  private
  public :: &
     basins_t,       &
     basins_init,    &
     basins_end,     &
     basins_analyze, &
     basins_write

  type basins_t
    integer, pointer :: map(:)

    integer          :: number
    integer, pointer :: position(:)
    FLOAT,   pointer :: value(:)
    FLOAT,   pointer :: volume(:)

    FLOAT,   pointer :: population(:)
  end type basins_t

contains

  !----------------------------------------------------------------
  subroutine basins_init(this, mesh)
    type(basins_t), intent(out) :: this
    type(mesh_t),   intent(in)  :: mesh
    
    PUSH_SUB(basins_init)

    SAFE_ALLOCATE(this%map(1:mesh%np))
    this%map(1:mesh%np) = -1

    POP_SUB(basins_init)
  end subroutine basins_init


  !----------------------------------------------------------------
  subroutine basins_end(this)
    type(basins_t), intent(inout) :: this

    PUSH_SUB(basins_end)

    ASSERT(associated(this%map))
    SAFE_DEALLOCATE_P(this%map)

    if(associated(this%position)) then
      SAFE_DEALLOCATE_P(this%position)
      SAFE_DEALLOCATE_P(this%value)
      SAFE_DEALLOCATE_P(this%volume)
      SAFE_DEALLOCATE_P(this%population)
    end if
    
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

    ASSERT(associated(this%map))

    message(1) = "Info: Calculating basins of attraction"
    call messages_info(1)

    cur_color = 0

    xmax = 1; ymax = 1; zmax = 1
    if(mesh%sb%dim < 3) zmax = 0
    if(mesh%sb%dim < 2) ymax = 0

    do jj = 1, mesh%np
      if(this%map(jj) .eq. -1) then
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

      PUSH_SUB(basins_analyze.steep_fill)

      if(this%map(ii) >= 0) then
        color = this%map(ii)
        POP_SUB(basins_analyze.steep_fill)
        return
      end if
      this%map(ii) = -2

      ii_max = get_max(ii, M_ZERO)
      if(ii_max == -1) then ! this is required to get ring attractors
        ii_max = get_max(ii, threshold)
      end if

      if(ii_max.ne.-1) then
        color = steep_fill(ii_max)
      else
        color = cur_color
        cur_color = cur_color + 1
      end if
      this%map(ii) = color

      POP_SUB(basins_analyze.steep_fill)
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
      if(mesh%parallel_in_domains) then
        ! When running in parallel, get global number of point i.
        call index_to_coords(mesh%idx, mesh%sb%dim, mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ii - 1), point)
      else
        call index_to_coords(mesh%idx, mesh%sb%dim, ii, point)
      end if

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

            index = index_from_coords(mesh%idx, mesh%sb%dim, point2)

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
      SAFE_ALLOCATE(this%value     (1:this%number))
      SAFE_ALLOCATE(this%volume    (1:this%number))
      SAFE_ALLOCATE(this%population(1:this%number))

      this%position(:) = -1
      this%value(:)    = M_ZERO
      this%volume(:)   = M_ZERO
      this%population(:)   = M_ZERO
      do ii = 1, this%number
        ii_max = -1
        f_max  = -huge(f_max)

        do jj = 1, mesh%np
          if(this%map(jj) .ne. ii-1) cycle
          if(f_max <= f(jj)) then
            ii_max = jj
            f_max  = f(jj)
          end if
          this%volume(ii) = this%volume(ii) + mesh%vol_pp(jj)
          this%population(ii) = this%population(ii) + mesh%vol_pp(jj)*sum(rho(jj, :))
        end do
        
        this%position(ii) = ii_max
        this%value(ii)    = f_max
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

    PUSH_SUB(basins_write)

    write(iunit, '(a,i5)') 'Number of basins = ', this%number
    write(iunit, '(1x)')

    do ii = 1, this%number
      write(iunit, '(a,i5)') '# ', ii
      write(iunit, '(a,3(f12.6,a))') '  position = (', &
         mesh%x(this%position(ii), 1), ',', &
         mesh%x(this%position(ii), 2), ',', &
         mesh%x(this%position(ii), 3), ')'
      write(iunit, '(a,f12.6)') '  value = ', this%value(ii)
      write(iunit, '(a,f12.6)') '  volume = ', this%volume(ii)
      write(iunit, '(a,f12.6)') '  population = ', this%population(ii)
    end do

    POP_SUB(basins_write)
  end subroutine basins_write

end module basins_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
