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

subroutine X(gridhier_init)(this, base_der, np_part_size)
  type(X(gridhier_t)),          intent(out) :: this
  type(derivatives_t),  target, intent(in)  :: base_der
  logical,                      intent(in)  :: np_part_size
  
  integer :: cl, ll, np
  type(derivatives_t), pointer :: der

  call push_sub('poisson_multigrid.Xgridhier_init')
  
  cl = multigrid_number_of_levels(base_der)

  SAFE_ALLOCATE(this%level(0:cl))
  
  der => base_der
  do ll = 0, cl
    if(np_part_size) then
      np = der%mesh%np_part
    else
      np = der%mesh%np
    end if

    SAFE_ALLOCATE(this%level(ll)%p(1:np))
    this%level(ll)%p(1:np) = M_ZERO
    
    der => der%coarser
  end do
  
  call pop_sub('poisson_multigrid.Xgridhier_init')
end subroutine X(gridhier_init)

! ---------------------------------------------------------
subroutine X(gridhier_end)(this)
  type(X(gridhier_t)),  intent(inout) :: this

  integer :: ll
  call push_sub('poisson_multigrid.Xgridhier_end')
  
  do ll = 0, ubound(this%level, dim = 1)
    SAFE_DEALLOCATE_P(this%level(ll)%p)
  end do

  SAFE_DEALLOCATE_P(this%level)
  
  call pop_sub('poisson_multigrid.Xgridhier_end')
end subroutine X(gridhier_end)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
