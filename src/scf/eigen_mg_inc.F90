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
!! $Id: eigen_mg_inc.F90 4195 2008-05-25 18:15:35Z xavier $

! ---------------------------------------------------------
subroutine X(eigen_solver_mg) (gr, st, h, pre, tol, niter, converged, diff, verbose)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  type(preconditioner_t), intent(in) :: pre
  FLOAT,               intent(in)    :: tol
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged
  FLOAT,     optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,   optional, intent(in)    :: verbose

  integer  :: ik, p, iter, maxter, conv, conv_, idim, ns, ip, is, inb
  R_TYPE   :: aa, rta, rtb, rtc, alpha, beta, dh, sigma, vv
  FLOAT    :: bb, res
  R_TYPE, allocatable :: hpsi(:, :), hdiag(:, :)

  call push_sub('eigen_cg.eigen_solver_mg')

  ALLOCATE(hpsi(1:NP, 1:st%d%dim), NP*st%d%dim)
  ALLOCATE(hdiag(1:NP, 1:st%d%dim), NP*st%d%dim)

  do ik = 1, st%d%nik
    
    call X(hpsi_diag) (h, gr, hdiag, ik)

    do p = conv + 1, st%nst

      do iter = 1, niter*100
        
        call X(hpsi)(h, gr, st%X(psi)(:,:, p, ik), hpsi, p, ik)
        
        bb = X(states_nrm2)(gr%m, st%d%dim, st%X(psi)(:, :, p, ik))
        call lalg_scal(NP, M_ONE/bb, st%X(psi)(:, 1, p, ik))
        
        aa = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, p, ik), hpsi)

        st%eigenval(p, ik) = aa

        res = X(states_residue)(gr%m, st%d%dim, hpsi, st%eigenval(p, ik), st%X(psi)(:, :, p, ik))
 
        if(res < tol) exit

!        print*, st%eigenval(p, ik), res
       
        bb = M_ONE
        
        do ip = 1, NP
          
          vv = sqrt(gr%m%vol_pp(ip))

          ! apply the hamiltonian in the point

          alpha = M_ZERO
          do is = 1, gr%LAP%n
            inb = ip + gr%LAP%ri(is, gr%LAP%rimap(ip))
            alpha = alpha + gr%LAP%w_re(is, 1)*st%X(psi)(inb, 1, p, ik)
          end do

          alpha = -M_HALF*alpha + (h%vhxc(ip, 1) + h%ep%vpsl(ip))*st%X(psi)(ip, 1, p, ik)

          alpha = alpha*vv

          beta = st%X(psi)(ip, 1, p, ik)*vv
          dh = hdiag(ip, 1)
          
          rta = alpha - dh*beta
          rtb = bb*dh - aa
          rtc = aa*beta - bb*alpha
          
          sigma = -M_TWO*rtc/(rtb + sqrt(rtb**2 - M_FOUR*rta*rtc))
          
          st%X(psi)(ip, 1, p, ik) = st%X(psi)(ip, 1, p, ik) - sigma/vv
          
          aa = aa - M_TWO*sigma*alpha + sigma**2*dh
          bb = bb - M_TWO*sigma*beta  + sigma**2
        end do
        
        st%eigenval(p, ik) = aa/bb
        
      end do
      
      niter = iter

      if(present(diff)) then
        diff(p, ik) = res
      end if

    end do
  end do

  deallocate(hpsi)

  call pop_sub()
end subroutine X(eigen_solver_mg)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
