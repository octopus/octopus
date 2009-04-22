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
!! $Id: gridhier_inc.F90 4205 2008-05-29 07:54:35Z xavier $

subroutine X(gridhier_init)(this, mgrid, add_points_for_boundaries)
  type(X(gridhier_t)),  intent(out) :: this
  type(multigrid_t),    intent(in)  :: mgrid
  logical,              intent(in)  :: add_points_for_boundaries
  
  integer :: cl, l
  call push_sub('poisson_multigrid.gridhier_init')
  
  cl = mgrid%n_levels
  
  SAFE_ALLOCATE(this%level(0:cl))
  
  do l = 0, cl
    if(add_points_for_boundaries) then
      SAFE_ALLOCATE(this%level(l)%p(1:mgrid%level(l)%mesh%np_part))
      this%level(l)%p(1:mgrid%level(l)%mesh%np_part) = M_ZERO
    else
      SAFE_ALLOCATE(this%level(l)%p(1:mgrid%level(l)%mesh%np))
      this%level(l)%p(1:mgrid%level(l)%mesh%np) = M_ZERO
    end if
  end do
  
  call pop_sub()
end subroutine X(gridhier_init)

! ---------------------------------------------------------
subroutine X(gridhier_end)(this, mgrid)
  type(X(gridhier_t)),  intent(inout) :: this
  type(multigrid_t),    intent(in)    :: mgrid

  integer :: cl, l
  call push_sub('poisson_multigrid.gridhier_end')
  
  cl = mgrid%n_levels
  
  do l = 0, cl
    SAFE_DEALLOCATE_P(this%level(l)%p)
  end do
  SAFE_DEALLOCATE_P(this%level)
  
  call pop_sub()
end subroutine X(gridhier_end)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
