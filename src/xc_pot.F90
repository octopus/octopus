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

subroutine R_FUNC(xc_pot) (xcs, m, st, hartr, vxc, ex, ec)
  type(xc_type), intent(inout) :: xcs
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  type(hartree_type), intent(inout) :: hartr
  real(r8), intent(out)    :: vxc(m%np, st%nspin), ex, ec

  real(r8), allocatable :: vaux(:, :) 

  ! for fxc != vxc...
  ! fxc is always LDA!!!!
!!$  real(r8), allocatable, save :: save_vxc(:,:)
!!$  logical, save :: first_time = .true.

  sub_name = 'xc_pot'; call push_sub()

  allocate(vaux(m%np, st%nspin))
  vxc = 0.0_r8; vaux = 0.0_r8
  ex  = 0.0_r8; ec = 0.0_r8

  select case(xcs%x_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%x_func, m, st, vxc, ex)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%x_func, m, st, vxc, ex)
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%x_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vx, ex)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%x_func, m, st, hartr, vxc, ex)
  end select

  select case(xcs%c_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%c_func, m, st, vaux, ec)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%c_func, m, st, vaux, ec)
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%c_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vc, ec)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%c_func, m, st, hartr, vaux, ec)
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
  call pop_sub(); return
end subroutine R_FUNC(xc_pot)

subroutine R_FUNC(xc_lda) (func, m, st, pot, energy)
  integer, intent(in) :: func
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: energy, pot(m%np, st%nspin)

  real(r8) :: d(st%nspin), p(st%nspin), d1, d2, e
  integer  :: i

  sub_name = 'xc_lda'; call push_sub()

  energy = 0._r8
  do i = 1, m%np
    d(:) = st%rho(i, :)

    ! WARNING: is this OK?
    if(st%nlcc) then    
      d(1) = d(1) + st%rho_core(i)/st%spin_channels
    end if

    select case(func)
    case(X_FUNC_LDA_NREL)
      call xc_x_lda(.false., st%nspin, d, p, e) 
    case(X_FUNC_LDA_REL)
      call xc_x_lda(.true.,  st%nspin, d, p, e) 
    case(C_FUNC_LDA_PZ)
      call xc_c_pz(st%nspin, d, p, e)
    case(C_FUNC_LDA_PW92)
      call xc_c_pw92(st%nspin, d, p, e)
    end select

    ! WARNING: this is probably not OK.
    energy = energy + d(1) * e * m%vol_pp
    pot(i, :) = p(:)
  end do

  call pop_sub(); return
end subroutine R_FUNC(xc_lda)
