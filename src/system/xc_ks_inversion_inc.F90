!! Copyright (C) 2010 H. Appel
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
!! $Id: $

! ---------------------------------------------------------
subroutine X(xc_ks_inversion_calc)(ks_inversion, gr, hm, st, ex, ec, vxc)
  use xc_functl_m

  type(xc_ks_inversion_t),  intent(in)    :: ks_inversion
  type(grid_t),             intent(inout) :: gr
  type(hamiltonian_t),      intent(inout) :: hm
  type(states_t),           intent(inout) :: st
  FLOAT,                    intent(inout) :: ex, ec
  FLOAT, optional,          intent(inout) :: vxc(:,:) ! vxc(gr%mesh%np, st%d%nspin)

  if(ks_inversion%level == XC_KS_INVERSION_NONE) return

  call push_sub('xc_ks_inversion_inc.Xxc_ks_inversion_calc')
  
  ! compute ks inversion
  select case (ks_inversion%level)
    ! adiabatic ks inversion
  case(XC_KS_INVERSION_ADIABATIC)
    ! TODO: compute density for input states
    !       call invertvxc_iter with auxilary hamiltonian and states
    !       update vxc
  end select
  
  call pop_sub()

end subroutine X(xc_ks_inversion_calc)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
