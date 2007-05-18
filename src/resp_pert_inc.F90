!! Copyright (C) 2007 the octopus team
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
!! $Id: response.F90 2548 2006-11-06 21:42:27Z xavier $

subroutine X(resp_pert_apply) (this, gr, f_in, f_out)
  type(resp_pert_t), intent(in)    :: this
  type(grid_t),      intent(inout) :: gr
  R_TYPE,            intent(inout) :: f_in(:)
  R_TYPE,            intent(out)   :: f_out(:)

  R_TYPE, allocatable :: lf(:,:)

  ASSERT(this%dir.ne.-1)

  select case(this%resp_type)
  case(RESP_PERTURBATION_ELECTRIC)
    f_out(1:NP) = f_in(1:NP) * gr%m%x(1:NP, this%dir)

  case(RESP_PERTURBATION_MAGNETIC)
    ! I believe this should be NP and not NP_PART
    ALLOCATE(lf(gr%m%np_part, gr%sb%dim), gr%m%np_part*gr%sb%dim)
    
    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%f_der, f_in, lf)
    f_out(1:NP) = lf(1:NP, this%dir)/M_TWO
    
    deallocate(lf)
  end select
  
end subroutine X(resp_pert_apply)
