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

subroutine R_FUNC(xc_pot) (xcs, m, st, hartr, vx, vc, ex, ec, vx_off, vc_off)
  type(xc_type), intent(inout) :: xcs
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  type(hartree_type), intent(inout) :: hartr
  real(r8), intent(out) :: vx(m%np, st%nspin), vc(m%np, st%nspin), ex, ec
  R_TYPE, pointer :: vx_off(:), vc_off(:)

  ! for fxc != vxc...
  ! fxc is always LDA!!!!
!!$  real(r8), allocatable, save :: save_vxc(:,:)
!!$  logical, save :: first_time = .true.

  sub_name = 'xc_pot'; call push_sub()

  vx = 0.0_r8; vc = 0.0_r8
  ex = 0.0_r8; ec = 0.0_r8
  if(st%ispin == 3) then
    vx_off = R_TOTYPE(0._r8)
    vc_off = R_TOTYPE(0._r8)
  end if

  select case(xcs%x_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%x_func, m, st, vx, ex, vx_off)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%x_func, m, st, vx, ex)
#ifdef HAVE_LAPACK
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%x_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vx, ex)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%x_func, m, st, hartr, vx, ex)
#endif
  end select

  select case(xcs%c_family)
  case(XC_FAMILY_ZER)
  case(XC_FAMILY_LDA)
    call R_FUNC(xc_lda) (xcs%c_func, m, st, vc, ec, vc_off)
  case(XC_FAMILY_GGA)
    call xc_gga(xcs%c_func, m, st, vc, ec)
#ifdef HAVE_LAPACK
!  case(XC_FAMILY_MGGA)
!    call R_FUNC(xc_mgga) (xcs%c_func, xcs, m, nst, st%nspin, psi, occ, eigenval, &
!        rho, vc, ec)
  case(XC_FAMILY_KLI)
    call R_FUNC(xc_kli) (xcs%c_func, m, st, hartr, vc, ec)
#endif
  end select

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

  call pop_sub()
  return
end subroutine R_FUNC(xc_pot)

subroutine R_FUNC(xc_lda) (func, m, st, pot, energy, pot_off)
  integer, intent(in) :: func
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: pot(m%np, st%nspin), energy
  R_TYPE, pointer, optional :: pot_off(:)

  real(r8) :: d(st%nspin), p(st%nspin), d1, d2, e
  R_TYPE   :: ri, a(2,2)
  integer  :: i

  energy = 0._r8
  do i = 1, m%np
    if(st%ispin == 3) then
      ri = st%R_FUNC(rho_off)(i)

      if(R_ABS(ri) < 1e-8_r8) then
        d(1:2) = st%rho(i, 1:2)
        a(1,1) = R_TOTYPE(1._r8); a(1,2) = R_TOTYPE(0._r8)
        a(2,1) = R_TOTYPE(0._r8); a(2,2) = R_TOTYPE(1._r8)
      else
        d1 = st%rho(i, 1) + st%rho(i, 2)
        d2 = st%rho(i, 1) - st%rho(i, 2)
        d2 = sqrt(d2*d2 + 4._r8*R_ABS(ri)**2)
        
        ! eigenvalues of density matrix
        d(1) = (d1 + d2) / 2._r8
        !d(2) = (d1 - d2) / 2._r8
        
        ! normalization
        d1 = sqrt( (st%rho(i, 1) - d(1))**2 + R_ABS(ri)**2 )
        !d2 = sqrt( (st%rho(i, 2) - d(2))**2 + R_ABS(ri)**2 )
        
        ! eigenfunctions
        a(1,1) = (st%rho(i, 1) - d(1)) / d1
        a(2,1) = -ri / d1
        a(2,2) =  a(1,1)
        a(1,2) = -a(2,1)
      end if
      
    else
      d(:) = st%rho(i, :)
    end if

    if(st%nlcc) then    
      d(:) = d(:) + st%rho_core(i)/st%nspin
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

    energy = energy + sum(d) * e * m%vol_pp

    if(st%ispin == 3) then ! rotate bak potential
      pot(i, 1) = R_ABS(a(1,1))**2 * p(1) + R_ABS(a(1,2))**2 * p(2)
      pot(i, 2) = R_ABS(a(2,1))**2 * p(1) + R_ABS(a(2,2))**2 * p(2)
      pot_off(i) = a(1,1)*R_CONJ(a(2,1))*p(1) + a(1,2)*R_CONJ(a(2,2))*p(2)
    else
      pot(i, :) = p(:)
    end if
  end do

  return
end subroutine R_FUNC(xc_lda)
