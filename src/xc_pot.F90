!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

subroutine X(xc_pot) (xcs, m, f_der, st, vxc, ex, ec, ip, qtot)
  type(xc_type), intent(inout) :: xcs
  type(mesh_type), intent(IN) :: m
  type(f_der_type), intent(inout) :: f_der
  type(states_type), intent(inout) :: st
  FLOAT, intent(out)    :: vxc(m%np, st%d%nspin), ex, ec
  FLOAT, intent(in) :: ip, qtot

  integer :: i

  ! for fxc != vxc...
  ! fxc is always LDA!!!!
!!$  FLOAT, allocatable, save :: save_vxc(:,:)
!!$  logical, save :: first_time = .true.

  call push_sub('xc_pot')
  
  ! "ip" variable is only meaningfull if LB94 potential is to be used. Should
  ! bear the opposite of the value of the last eigenvalue ( = ionization potential)
  ! If "self-consistent" LB94 is not to be used, should be 1/32
  ! ip = M_ONE/CNST(32.0)
  
  ex  = M_ZERO
  ec  = M_ZERO

  do i = 0, N_XC_FAMILIES - 1
    ! The XC_FAMILY_OEP is handled separatly by h_calc_vhxc
    if(btest(xcs%family, i)) then
      select case(ibset(0, i))
      case(XC_FAMILY_LDA)
        call xc_get_lda(xcs, m, st, vxc, ex, ec)
      case(XC_FAMILY_GGA)
        call xc_get_gga(xcs, m, f_der, st, vxc, ex, ec, ip, qtot)
!!$      case(XC_FAMILY_MGGA)
!!$        call X(xc_mgga) (xcs%functl, xcs, m, nst, st%d%nspin, psi, occ, eigenval, &
!!$             rho, vx, ex)
      end select
    end if
  end do

  ! Warning: For vxc != vxc
!!$  if(first_time) then
!!$    first_time = .false.
!!$    allocate(save_vxc(m%np, st%d%nspin))
!!$    save_vxc = vx + vc
!!$
!!$    ! now, we get the LDA xc
!!$    xcs%x_family = XC_FAMILY_LDA
!!$    xcs%x_func   = X_FUNC_LDA_NREL
!!$    call xc_lda(xcs%x_func, m, st, vx, ex)
!!$
!!$    xcs%c_family = XC_FAMILY_LDA
!!$    xcs%c_func   = C_FUNC_LDA_PZ
!!$    call xc_lda(xcs%c_func, m, st, vc, ec)
!!$
!!$    save_vxc = save_vxc - vx - vc
!!$  end if
!!$  vx = vx + save_vxc

  call pop_sub()
end subroutine X(xc_pot)
