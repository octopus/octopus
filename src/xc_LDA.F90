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

subroutine R_FUNC(xc_lda) (func, nlcc, m, st, pot, energy)
  integer, intent(in) :: func
  logical, intent(in) :: nlcc
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: energy, pot(m%np, st%nspin)

  real(r8) :: d(st%spin_channels), dtot, dpol, vpol, p(st%nspin), pd(st%spin_channels), e, &
              dummy_x, dummy_vx(st%spin_channels)
  real(r8), parameter :: tiny = 1.0e-12_r8
  integer  :: i, is

  sub_name = 'xc_lda'; call push_sub()

  energy = M_ZERO
  do i = 1, m%np
    if(st%ispin==SPINORS) then
     dtot = st%rho(i, 1) + st%rho(i, 2)
     dpol = sqrt( (st%rho(i, 1)-st%rho(i, 2))**2 + M_FOUR*(st%rho(i, 3)**2+st%rho(i, 4)**2) )
     d(1) = M_HALF * ( dtot + dpol )
     d(2) = M_HALF * ( dtot - dpol )
    else
     do is = 1, st%spin_channels
        d(is) = st%rho(i, is)
     enddo
    endif
    if(nlcc) then    
      d(1:st%spin_channels) = d(1:st%spin_channels) + st%rho_core(i)/st%spin_channels
    end if

    select case(func)
    case(X_FUNC_LDA_NREL)
      call exchng(0, st%spin_channels, d, e, pd)
    case(X_FUNC_LDA_REL)
      call exchng(1, st%spin_channels, d, e, pd) 
    case(C_FUNC_LDA_PZ)
      call pzxc(1, st%spin_channels, d, dummy_x, e, dummy_vx, pd)
    case(C_FUNC_LDA_PW92)
      call pw92c(st%spin_channels, d, e, pd)
    end select

    if (st%ispin==SPINORS) then
       vpol  = (pd(1)-pd(2)) * (st%rho(i, 1)-st%rho(i, 2)) / (dpol+tiny)
       p(1) = M_HALF * ( pd(1) + pd(2) + vpol )
       p(2) = M_HALF * ( pd(1) + pd(2) - vpol )
       p(3) = (pd(1)-pd(2)) * st%rho(i, 3) / (dpol+tiny)
       p(4) = (pd(1)-pd(2)) * st%rho(i, 4) / (dpol+tiny)
    else
       do is = 1, st%nspin
          p(is) = pd(is)
       enddo
    endif

    energy = energy + sum(d(1:st%spin_channels)) * e * m%vol_pp
    pot(i, :) = p(:)
  end do

  call pop_sub(); return
end subroutine R_FUNC(xc_lda)
