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
!!
!! $Id$

! CONJUGATE-GRADIENTS METHOD.
subroutine eigen_solver_cg2(m, f_der, st, h, tol, niter, converged, errorflag, diff, reorder, verbose)
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(IN)    :: h
  FLOAT,                  intent(in)    :: tol
  integer,                intent(inout) :: niter
  integer,                intent(out)   :: errorflag
  integer,                intent(inout) :: converged
  FLOAT,        optional, intent(out)   :: diff(1:st%nst,1:st%d%nik)
  logical,      optional, intent(in)    :: reorder
  logical,      optional, intent(in)    :: verbose

  R_TYPE, allocatable :: h_psi(:,:), g(:,:), g0(:,:), &
       cg(:,:), ppsi(:,:)

  R_TYPE :: es(2), a0, b0, gg, gg0, gg1, gamma, theta, norma
  FLOAT :: cg0, e0, res
  integer  :: ik, np, moved, p, j, iter, i, maxter
  logical  :: reord = .true., verbose_

  call push_sub('eigen_solver_cg2')

  verbose_ = .false.; if(present(verbose)) verbose_ = verbose
  if(verbose_) then
    message(1) = " "
    message(2) = stars
    message(3) = "Diagonalization with the conjugate gradients algorithm."
    write(message(4),'(a,e8.2)') '  Tolerance: ',tol
    write(message(5),'(a,i6)')   '  Maximum number of iterations per eigenstate:', niter
    message(6) = ""
    call write_info(6)
  endif

  np = m%np

  if(present(reorder)) reord = reorder
  maxter = niter
  niter = 0
  moved = 0

  allocate(h_psi(np, st%d%dim), g(np, st%d%dim), g0(np, st%d%dim), &
       cg(np, st%d%dim), ppsi(np, st%d%dim))

  ! Start of main loop, which runs over all the eigenvectors searched
  ik_loop: do ik = 1, st%d%nik
    eigenfunction_loop : do p = converged + 1, st%nst

      if(verbose_) then
        write(message(2),'(a,i4,a)') 'Eigenstate # ',p,':'
      endif

      ! Orthogonalize starting eigenfunctions to those already calculated...
      call X(states_gram_schmidt)(p, m, st%d%dim, &
           st%X(psi)(1:np, 1:st%d%dim, 1:p, ik), start=p)

      ! Calculate starting gradient: |hpsi> = H|psi>
      call X(Hpsi)(h, m, f_der, st%X(psi)(:,:, p, ik) , h_psi, ik)

      ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
      st%eigenval(p, ik) = R_REAL(X(states_dotp) (m, st%d%dim, st%X(psi)(:,:, p, ik), h_psi))

      ! Starts iteration for this band
      iter_loop: do iter = 1, maxter

        ! if approximate inverse preconditioner....
        !call pre(hpsi%val   , g%val   ) 
        !call pre(psi(m)%val , ppsi%val)
        g = h_psi
        ppsi = st%X(psi)(:,:, p, ik)
        
        es(1) = X(states_dotp) (m, st%d%dim, st%X(psi)(:,:, p, ik), g)
        es(2) = X(states_dotp) (m, st%d%dim, st%X(psi)(:,:, p, ik), ppsi)
        es(1) = es(1)/es(2)
        g = g - es(1)*ppsi
        
        ! Orthogonalize to lowest eigenvalues (already calculated)
        do j = 1, p - 1
          a0 = X(states_dotp) (m, st%d%dim, st%X(psi)(:,:, j, ik), g)
          g(:,:) = g(:,:) - a0 * st%X(psi)(:,:, j, ik)
        end do
        
        if(iter .ne. 1) gg1 = X(states_dotp) (m, st%d%dim, g, g0)
        
        ! Approximate inverse preconditioner...
        !call ipre(g, g0)
        g0 = g
        
        gg = X(states_dotp) (m, st%d%dim, g, g0)
        
        ! Starting or following iterations...
        if(iter .eq. 1) then
          gg0 = gg
          cg(:,:) = g(:,:)
        else
          !gamma = gg/gg0        ! (Fletcher-Reeves)
          gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
          gg0 = gg
          cg(:,:) = gamma*cg(:,:)
          cg(:,:) = cg(:,:) + g(:,:)

          norma = gamma*cg0*sin(theta)
          cg(:,:) = cg(:,:) - norma * st%X(psi)(:,:, p, ik)
        end if
        
        ! cg contains now the conjugate gradient
        cg0 = X(states_nrm2) (m, st%d%dim, cg(:,:))
        call X(Hpsi) (h, m, f_der, cg, ppsi, ik)
        
        ! Line minimization.
        a0 = X(states_dotp) (m, st%d%dim, st%X(psi)(:,:, p, ik), ppsi)
        a0 = M_TWO * a0 / cg0
        b0 = X(states_dotp) (m, st%d%dim, cg(:,:), ppsi)
        b0 = b0/cg0**2
        e0 = st%eigenval(p, ik)
        theta = atan(R_REAL(a0/(e0 - b0)))/M_TWO
        es(1) = ((e0-b0)*cos(M_TWO*theta) + a0*sin(M_TWO*theta) + e0 + b0) / M_TWO
        es(2) =(-(e0-b0)*cos(M_TWO*theta) - a0*sin(M_TWO*theta) + e0 + b0) / M_TWO
        
        ! Choose the minimum solutions. 
        if (R_REAL(es(2)) < R_REAL(es(1))) then
          theta = theta + M_PI/M_TWO
        end if
        st%eigenval(p, ik) = min(R_REAL(es(1)), R_REAL(es(2)))
        
        ! Upgrade psi...
        a0 = cos(theta)
        b0 = sin(theta)/cg0
        ! This does the sum: st%X(psi)(:, :, p, ik) = a0*st%X(psi)(:, :, p, ik) + b0*cg(:, :)
        ! It can crash in Intel compiler version 8 otherwise.
        do j = 1, st%d%dim
           call lalg_scal(np, a0, st%X(psi)(:, j, p, ik))
           call lalg_axpy(np, b0, cg(:, j), st%X(psi)(:, j, p, ik))
        enddo
        
        ! Calculate H|psi>
        h_psi = a0*h_psi + b0*ppsi

        res = X(states_residue)(m, st%d%dim, h_psi, st%eigenval(p, ik), &
             st%X(psi)(:, :, p, ik))
        ! Test convergence.
        if(res < tol) then
          converged = converged + 1
          exit iter_loop
        end if
        
      end do iter_loop

      if(verbose_) then
        if(res<tol) then
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2)),"     converged. Iterations:", iter, '   [Res = ',res,']'
        else
          write(message(1),'(a,a,i5,a,e8.2,a)') trim(message(2))," not converged. Iterations:", iter, '   [Res = ',res,']'
        endif
        call write_info(1)
      endif

      niter = niter + iter + 1

      if(present(diff)) then
        diff(p, ik) = res
      end if

      ! Reordering... (this should be improved)
      if(p>1 .and. reord) then
        if(st%eigenval(p, ik) - st%eigenval(p-1, ik) .lt. -M_TWO*tol) then
          do i = p - 2, 1, -1
            if (st%eigenval(p, ik) - st%eigenval(i, ik) .gt. M_TWO*tol) exit
          end do
          i = i + 1
          moved = moved + 1
          e0 = st%eigenval(p, ik)
          cg = st%X(psi)(:,:, p, ik)
          do j = p, i + 1 , -1
            st%eigenval(j, ik) = st%eigenval(j-1, ik)
            st%X(psi)(:,:, j, ik) = st%X(psi)(:,:, j-1, ik)
          end do
          st%eigenval(i, ik) = e0
          st%X(psi)(:,:, i, ik) = cg
        end if
      end if
      ! End of reordering
    end do eigenfunction_loop
  end do ik_loop

  ! Deallocation of variables
  deallocate(h_psi, g, g0, cg, ppsi)

  errorflag = 0

  if(verbose_) then
    message(1) = stars
    message(2) = " "
    call write_info(2)
  endif
 
  call pop_sub()
end subroutine eigen_solver_cg2


