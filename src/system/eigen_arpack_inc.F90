!!!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
	
	
subroutine X(eigen_solver_arpack)(gr, st, hm, tol_, niter, ncv, converged, ik, diff)
  type(grid_t),        intent(in)    :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  FLOAT,               intent(in)    :: tol_
  integer,             intent(inout) :: niter
  integer,             intent(in)    :: ncv
  integer,             intent(inout) :: converged
  integer,             intent(in)    :: ik
  FLOAT,     optional, intent(out)   :: diff(1:st%nst)
	
  logical, allocatable :: select(:)
  R_TYPE, allocatable  :: ax(:),  resid(:), v(:, :),   &
                          workd(:), workev(:), workl(:), zd(:), &
                          psi(:,:)
                     
  integer :: ldv, nev, iparam(11), ipntr(14), ido, n, lworkl, info, ierr, &
             i, j, ishfts, maxitr, mode1, ist
  FLOAT :: tol, sigmar, sigmai
  FLOAT, allocatable :: rwork(:), d(:, :) 
  CMPLX :: sigma 	
	!!!!WARNING: No support for spinors, yet. No support for complex wavefunctions.
  PUSH_SUB('eigen_arpack.eigen_solver_arpack')
	

	
#if defined(HAVE_MPI)
  if(gr%mesh%parallel_in_domains) then
    message(1) = 'Error: Arpack-Solver not parallelized for domain decomposition.'
    call messages_fatal(1)
    !  FIXME: Need to adjust mesh%x and mesh%vol_pp occurences in the code below
    !         appropriately for domain decomposition. Also parallelization
    !         of the vectors has to be taken care of.
  end if
#endif
	
  ldv = gr%mesh%np
  n = gr%mesh%np
  nev = st%nst
  lworkl  = 3*ncv**2+6*ncv
  SAFE_ALLOCATE(ax(ldv))
  SAFE_ALLOCATE(d(ncv+1, 3))
  SAFE_ALLOCATE(resid(ldv))
  SAFE_ALLOCATE(v(ldv, ncv))
  SAFE_ALLOCATE(workd(3*ldv))
  SAFE_ALLOCATE(workev(3*ncv))
  SAFE_ALLOCATE(workl(lworkl))
  SAFE_ALLOCATE(select(ncv))
  
  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
  
#if defined(R_TCOMPLEX)
  SAFE_ALLOCATE(rwork(ncv))
  SAFE_ALLOCATE(zd(ncv+1))
#endif
	
  select = .true.
  tol    = tol_
  ido    = 0
  info = 1
  
  do i = 1, gr%mesh%np
!      resid(i) = sum(st%X(psi)(i, 1, 1:st%nst, ik))*sqrt(gr%mesh%vol_pp(1))
      resid(i) = R_TOTYPE(M_ONE)
  end do

	
  ishfts = 1
  maxitr = niter
  mode1 = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
	
  do
#if defined(R_TCOMPLEX)      
    call znaupd  ( ido, 'I', n, 'SR', nev, tol, resid, ncv, &
               v, ldv, iparam, ipntr, workd, workl, lworkl, &
               rwork,info )
#else 
    call dnaupd  ( ido, 'I', n, 'SR', nev, tol, resid, ncv, &
               v, ldv, iparam, ipntr, workd, workl, lworkl, & 
               info )
#endif      
      
    if( abs(ido).ne.1) exit
    call av (n, workd(ipntr(1)), workd(ipntr(2)))
  end do
  ! If info is larger than zero, it may not be an error (i.e., not all eigenvectors
  ! were converged)
  if(info .lt. 0) then
    write(message(1),'(a,i5)') 'Error with ARPACK _naupd, info = ', info
    write(message(2),'(a)')    'Check the documentation of _naupd.'
    call messages_fatal(2)
  end if

#if defined(R_TCOMPLEX) 
  call zneupd  (.true., 'A', select, zd, v, ldv, sigma, &
        workev, 'I', n, 'SR', nev, tol, resid, ncv, & 
        v, ldv, iparam, ipntr, workd, workl, lworkl, &
        rwork, ierr)
        d(:,1)=real(zd(:))
        d(:,2)=aimag(zd(:))
        d(:,3)=M_ZERO
        
#else	
  call dneupd ( .true., 'A', select, d, d(1,2), v, ldv, &
       sigmar, sigmai, workev, 'I', n, 'SR', nev, tol, &
       resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
       lworkl, ierr )
#endif
       
  if(ierr .ne. 0) then
    write(message(1),'(a,i5)') 'Error with ARPACK _neupd, info = ', info
    write(message(2),'(a)')    'Check the documentation of _neupd.'
    call messages_fatal(2)
  end if

  ! This sets the number of converged eigenvectors.
  converged =  iparam(5)

  call dmout(6, converged, 3, d, ncv+1, -6, 'Ritz values (Real, Imag) and residual residuals')

  ! This sets niter to the number of matrix-vector operations.
  niter = iparam(9)
  do j = 1, converged
    do i = 1, gr%mesh%np
!       st%X(psi)(i, 1, j, ik) = v(i, j)/sqrt(gr%mesh%vol_pp(1))
      psi(i,1) = v(i, j)/sqrt(gr%mesh%volume_element) 
    end do
    call states_set_state(st, gr%mesh, j, ik, psi)
    
    print *,"st", j, "norm", sqrt(X(mf_dotp)(gr%mesh, st%d%dim, psi, psi, dotu = .true.))!,&
!      sqrt(sum(psi(:, 1)*psi(:, 1))*gr%mesh%volume_element), sqrt(X(mf_integrate)(gr%mesh,psi(:,1)**2))
        
    st%eigenval(j, ik) = d(j, 1)
    if(associated(st%zeigenval%Im))then 
      st%zeigenval%Im(j, ik) = d(j, 2)
    end if

    if(abs(workl(ipntr(11)+j-1))< M_EPSILON) then
      diff(j) = M_ZERO
    else
      diff(j) = workl(ipntr(11)+j-1)
    end if
  end do

  !Fill unconverged states
  do j = converged + 1, st%nst
    do i = 1, gr%mesh%np
!       st%X(psi)(i, 1, j, ik) = R_TOTYPE(M_ONE)
      psi(i,1) = R_TOTYPE(M_ONE) 
    end do
    call states_set_state(st, gr%mesh, j, ik, psi)

    st%eigenval(j, ik) = M_HUGE
    if(associated(st%zeigenval%Im))then 
      st%zeigenval%Im(j, ik) = M_HUGE
    end if
    diff(j) = M_HUGE
  end do


  SAFE_DEALLOCATE_A(ax)
  SAFE_DEALLOCATE_A(d)
  SAFE_DEALLOCATE_A(resid)
  SAFE_DEALLOCATE_A(v)
  SAFE_DEALLOCATE_A(workd)
  SAFE_DEALLOCATE_A(workev)
  SAFE_DEALLOCATE_A(workl)
  SAFE_DEALLOCATE_A(select)
  
  SAFE_DEALLOCATE_A(psi)
  
#if defined(R_TCOMPLEX)
  SAFE_DEALLOCATE_A(rwork)
  SAFE_DEALLOCATE_A(zd)  
#endif



   POP_SUB('eigen_arpack.eigen_solver_arpack')
contains

  ! ---------------------------------------------------------
  subroutine av (n, v, w)
    integer :: n
    R_TYPE :: v(n), w(n)
    integer :: i, NP, NP_PART
    R_TYPE, allocatable :: psi(:, :), hpsi(:, :)
    
    NP = gr%mesh%np
    NP_PART = gr%mesh%np_part

    SAFE_ALLOCATE(psi(NP_PART, hm%d%dim))
    SAFE_ALLOCATE(hpsi(NP_PART, hm%d%dim))

    do i = 1, NP
!       print *,i, NP, NP_PART, ist
      psi(i, 1) = v(i)/sqrt(gr%mesh%vol_pp(1))
    end do
    do i = NP+1, NP_PART
      psi(i, 1) = M_ZERO
    end do
! psi(1:NP,1) = v(:)
! psi(NP+1:NP_PART,1) = M_ZERO
    
    
    call X(hamiltonian_apply) (hm, gr%der, psi, hpsi, 1, ik)
    
! w(:) = hpsi(1:NP,1)
        
    do i = 1, NP
      w(i) = hpsi(i, 1)*sqrt(gr%mesh%vol_pp(1))
    end do


    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

  end subroutine av

end subroutine X(eigen_solver_arpack)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End: