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

subroutine xc_get_lda(xcs, m, st, vxc, ex, ec)
  type(xc_type),     intent(in)  :: xcs
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(in)  :: st
  FLOAT,          intent(out) :: ex, ec, vxc(m%np, st%d%nspin)
  
  FLOAT, allocatable :: d(:), p(:), pd(:), pd1(:)
  FLOAT :: dtot, dpol, vpol, e
  integer  :: i, ixc, is, spin_channels

  call push_sub('xc_lda')

  spin_channels = st%d%spin_channels

  allocate(d(spin_channels), p(st%d%nspin), pd(spin_channels), pd1(spin_channels))

  do i = 1, m%np
    if(st%d%ispin==SPINORS) then
      dtot = st%rho(i, 1) + st%rho(i, 2)
      dpol = sqrt( (st%rho(i, 1)-st%rho(i, 2))**2 + M_FOUR*(st%rho(i, 3)**2+st%rho(i, 4)**2) )
      d(1) = max(M_HALF * ( dtot + dpol ), M_ZERO)
      d(2) = max(M_HALF * ( dtot - dpol ), M_ZERO)
    else
      do is = 1, spin_channels
        d(is) = max(st%rho(i, is), M_ZERO)
      enddo
    endif
    if(xcs%nlcc) then
      d(1:spin_channels) = d(1:spin_channels) + st%rho_core(i)/spin_channels
    end if

    pd = M_ZERO
    functl_loop: do ixc = 1, XC_LDA_N
      if(.not.btest(xcs%lda_functl, ixc)) cycle

      ! call xc library
      call xc_lda(xcs%lda_conf(ixc), d(1), e, pd1(1))

      if(ixc == XC_LDA_X) then
        ex = ex + sum(d(1:spin_channels)) * e * m%vol_pp(i)
      else
        ec = ec + sum(d(1:spin_channels)) * e * m%vol_pp(i)
      end if

      pd = pd + pd1
    end do functl_loop

    if (st%d%ispin==SPINORS) then
      vpol  = (pd(1)-pd(2)) * (st%rho(i, 1)-st%rho(i, 2)) / (dpol+tiny)
      p(1) = M_HALF * ( pd(1) + pd(2) + vpol )
      p(2) = M_HALF * ( pd(1) + pd(2) - vpol )
      p(3) = (pd(1)-pd(2)) * st%rho(i, 3) / (dpol+tiny)
      p(4) = (pd(1)-pd(2)) * st%rho(i, 4) / (dpol+tiny)
    else
      do is = 1, st%d%nspin
        p(is) = pd(is)
      enddo
    endif
    
    vxc(i, :) = vxc(i, :) + p(:)
  end do

  deallocate(d, p, pd, pd1)

  call pop_sub()
end subroutine xc_get_lda
