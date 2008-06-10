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
  type(grid_t),           intent(inout) :: gr
  type(states_t),         intent(inout) :: st
  type(hamiltonian_t),    intent(inout) :: h
  type(preconditioner_t), intent(in) :: pre
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(inout) :: converged
  FLOAT,                  intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,   optional,    intent(in)    :: verbose

  integer  :: ik, ist, ist2, iter, maxter, conv, conv_, idim, ns, ip, is, inb
  R_TYPE   :: rta, rtb, rtc, alpha, dh, vv, s1, s2, s3
  FLOAT    :: ww, bb
  R_TYPE, allocatable :: hpsi(:, :), hdiag(:, :), cc(:, :), sigma(:), beta(:), aa(:)

  call push_sub('eigen_cg.eigen_solver_mg')

  ALLOCATE(hpsi(1:NP, 1:st%d%dim), NP*st%d%dim)
  ALLOCATE(hdiag(1:NP, 1:st%d%dim), NP*st%d%dim)
  ALLOCATE(cc(1:st%nst, 1:st%nst), st%nst**2)
  ALLOCATE(aa(1:st%nst), st%nst)
  ALLOCATE(sigma(1:st%nst), st%nst)
  ALLOCATE(beta(1:st%nst), st%nst)

  cc = M_Z0

  ww = 5.0

  do ik = 1, st%d%nik
    
    call X(hpsi_diag) (h, gr, hdiag, ik)
    
    do iter = 1, niter*100

      do ist = 1, st%nst

        call X(hpsi)(h, gr, st%X(psi)(:,:, ist, ik), hpsi, ist, ik)
        
        bb = X(states_nrm2)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik))
        call lalg_scal(NP, M_ONE/bb, st%X(psi)(:, 1, ist, ik))
        
        cc(ist, 1:ist - 1) = cc(ist, 1:ist - 1)/bb

        aa(ist) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), hpsi)

        st%eigenval(ist, ik) = aa(ist)

        diff(ist, ik) = X(states_residue)(gr%m, st%d%dim, hpsi, st%eigenval(ist, ik), st%X(psi)(:, :, ist, ik))
 
        if(mod(iter, 100) ==  1) print*, iter, ist, st%eigenval(ist, ik), diff(ist, ik)

      end do
      
      do ist = 1, st%nst

        cc(ist, ist) = M_ONE
        do ist2 = 1, ist - 1
          cc(ist, ist2) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, ist, ik), st%X(psi)(:, :, ist2, ik))
        end do

      end do
      
      do ip = 1, NP

        vv = sqrt(gr%m%vol_pp(ip))
        dh = hdiag(ip, 1)
          
        do ist = 1, st%nst

          beta(ist) = st%X(psi)(ip, 1, ist, ik)*vv
          
          ! apply the hamiltonian in the point

          alpha = M_ZERO
          do is = 1, gr%LAP%n
            inb = ip + gr%LAP%ri(is, gr%LAP%rimap(ip))
            alpha = alpha + gr%LAP%w_re(is, 1)*st%X(psi)(inb, 1, ist, ik)
          end do

          alpha = -M_HALF*alpha + (h%vhxc(ip, 1) + h%ep%vpsl(ip))*st%X(psi)(ip, 1, ist, ik)

          alpha = alpha*vv

          s1 = M_ZERO
          s2 = M_ZERO
          s3 = M_ZERO
          do ist2 = 1, ist - 1 
            s1 = s1 + ww*beta(ist2)**2/cc(ist2, ist2)
            s2 = s2 + ww*beta(ist2)*cc(ist, ist2)/cc(ist2, ist2)
            s3 = s3 + ww*cc(ist, ist2)**2/cc(ist2, ist2)
          end do

          rta = alpha - dh*beta(ist) + s2 - beta(ist)*s1
          rtb = cc(ist, ist)*dh - aa(ist) + cc(ist, ist)*s1 - s3
          rtc = aa(ist)*beta(ist) - cc(ist, ist)*alpha + beta(ist)*s3 - s2
          
          sigma(ist) = -M_TWO*rtc/(rtb + sqrt(rtb**2 - M_FOUR*rta*rtc))
          
          st%X(psi)(ip, 1, ist, ik) = st%X(psi)(ip, 1, ist, ik) - sigma(ist)/vv
          
          aa(ist) = aa(ist) - M_TWO*sigma(ist)*alpha + sigma(ist)**2*dh
          cc(ist, ist) = cc(ist, ist) - M_TWO*sigma(ist)*beta(ist)  + sigma(ist)**2

          do ist2 = 1, ist - 1
            cc(ist, ist2) = cc(ist, ist2) - sigma(ist)*beta(ist2) - sigma(ist2)*beta(ist) + sigma(ist)*sigma(ist2)
          end do
          
        end do
        
      end do
      
    end do

    niter = iter

  end do

  call pop_sub()
end subroutine X(eigen_solver_mg)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
