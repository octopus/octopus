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

subroutine R_FUNC(xc_pot) (xcs, m, st, hartr, vxc, ex, ec, ip, qtot)
  type(xc_type), intent(inout) :: xcs
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  type(hartree_type), intent(inout) :: hartr
  real(r8), intent(out)    :: vxc(m%np, st%nspin), ex, ec
  real(r8), intent(in) :: ip, qtot

  real(r8), allocatable :: vaux(:, :) 

  ! for fxc != vxc...
  ! fxc is always LDA!!!!
!!$  real(r8), allocatable, save :: save_vxc(:,:)
!!$  logical, save :: first_time = .true.

  call push_sub('xc_pot')

  ! "ip" variable is only meaningfull if LB94 potential is to be used. Should
  ! bear the opposite of the value of the last eigenvalue ( = ionization potential)
  ! If "self-consistent" LB94 is not to be used, should be 1/32
  ! ip = M_ONE/32.0_r8
  
  allocate(vaux(m%np, st%nspin))
  vxc = M_ZERO; vaux = M_ZERO
  ex  = M_ZERO; ec = M_ZERO

  select case(xcs%x_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%x_func, xcs%nlcc, m, st, vxc, ex)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%x_func, xcs%nlcc, m, st, vxc, ex, &
                ip, qtot, xcs%lb94_modified, xcs%lb94_beta, xcs%lb94_threshold)
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%x_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vx, ex)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%x_func, xcs%nlcc, m, st, hartr, vxc, ex)
  end select

  select case(xcs%c_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%c_func, xcs%nlcc, m, st, vaux, ec)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%c_func, xcs%nlcc, m, st, vaux, ec, &
                ip, qtot, xcs%lb94_modified, xcs%lb94_beta, xcs%lb94_threshold)
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%c_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vc, ec)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%c_func, xcs%nlcc, m, st, hartr, vaux, ec)
  end select

  vxc = vxc + vaux

  ! Warning: For vxc != vxc
!!$  if(first_time) then
!!$    first_time = .false.
!!$    allocate(save_vxc(m%np, st%nspin))
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

  deallocate(vaux)
  call pop_sub()
end subroutine R_FUNC(xc_pot)
