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
	
	
subroutine X(eigen_solver_arpack)(arpack, gr, st, hm, tol_, niter, converged, ik, diff)
  type(eigen_arpack_t),intent(in)    :: arpack
  type(grid_t),        intent(in)    :: gr
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  FLOAT,               intent(in)    :: tol_
  integer,             intent(inout) :: niter
  integer,             intent(inout) :: converged
  integer,             intent(in)    :: ik
  FLOAT,     optional, intent(out)   :: diff(:) !< (1:st%nst)
	
  logical, allocatable :: select(:)
  R_TYPE, allocatable  :: resid(:), v(:, :), &
                          workd(:), workev(:), workl(:), zd(:), &
                          psi(:,:), hpsi(:,:)
                     
  integer :: ldv, nev, iparam(11), ipntr(14), ido, n, lworkl, info, ierr, &
             i, j, ishfts, maxitr, mode1, ist, idim, ncv
  FLOAT :: tol, sigmar, sigmai, resid_sum, tmp
  FLOAT, allocatable :: rwork(:), d(:, :)
  CMPLX :: sigma, eps_temp
  integer :: mpi_comm
  character(len=2) :: which
  	
	!!!!WARNING: No support for spinors, yet. 
  PUSH_SUB(eigen_arpack.eigen_solver_arpack)

  !Enable debug info
  if(in_debug_mode) call arpack_debug(conf%debug_level)
  
  mpi_comm = mpi_world%comm
  if (gr%mesh%parallel_in_domains) mpi_comm = gr%mesh%mpi_grp%comm
  
  ncv = arpack%arnoldi_vectors
  n = gr%mesh%np
  ldv = gr%mesh%np
  nev = st%nst
  lworkl  = 3*ncv**2+6*ncv

  SAFE_ALLOCATE(d(ncv+1, 3))
  SAFE_ALLOCATE(resid(ldv))       !residual vector 
  SAFE_ALLOCATE(v(ldv, ncv))      !Arnoldi basis vectors / Eigenstates
  SAFE_ALLOCATE(workd(3*ldv))
  SAFE_ALLOCATE(workev(3*ncv))
  SAFE_ALLOCATE(workl(lworkl))
  SAFE_ALLOCATE(select(ncv))
  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
  
#if defined(R_TCOMPLEX)
  SAFE_ALLOCATE(rwork(ncv))
  SAFE_ALLOCATE(zd(ncv+1))
#endif
  which = arpack%sort	
  select(:) = .true.
  tol  = tol_
  ido  = 0
  info = arpack%init_resid ! 0. random resid vector 
                           ! 1. calculate resid vector 
                           ! 2. resid vector constant = 1 
  
  if(info == 1) then !Calculate the residual vector
    print*, 'eigen_solver_arpack info1, allocate for calculating residual.'
    SAFE_ALLOCATE(hpsi(1:gr%mesh%np_part, 1:st%d%dim))
  
    resid(:) = R_TOTYPE(M_ZERO)
    do ist = 1, st%nst
      call states_get_state(st, gr%mesh, ist, ik, psi)      
      do idim = 1, st%d%dim
       call X(hamiltonian_apply) (hm, gr%der, psi, hpsi, idim, ik)
       ! XXX this will hardly work because tmp is not necessarily written to...
       ! In fact as of lately, it is never written to as the mentioned sorting trick is not in use
       !if (st%eigenval(ist, ik) > CNST(1e3)) tmp = st%eigenval(ist, ik) -  CNST(1e3) ! compensate the ugly sorting trick
       resid(1:ldv) = resid(1:ldv) + hpsi(1:ldv, idim) - tmp * psi(1:ldv, idim)
       if(associated(st%zeigenval%Im)) resid(1:ldv) = resid(1:ldv) - M_zI * st%zeigenval%Im(ist, ik) * psi(1:ldv, idim)
      end do
    end do
    !resid(:) = resid(:) * sqrt(gr%mesh%volume_element)
    SAFE_DEALLOCATE_A(hpsi)

    resid_sum = abs(sum(resid(:)**2))
    print *,"residual", resid_sum
    if(resid_sum < M_EPSILON .or. resid_sum > M_HUGE) then
      resid(:) = R_TOTYPE(M_ONE)
    end if
    
  else
    resid(:) = R_TOTYPE(M_ONE)
  end if
  
!   do i = 1, ldv
! !      resid(i) = sum(st%X(psi)(i, 1, 1:st%nst, ik))*sqrt(gr%mesh%vol_pp(1))
!       resid(i) = R_TOTYPE(M_ONE)
!   end do
!   ishfts = 1
!   maxitr = niter
!   mode1 = 1
  iparam(1) = 1
  iparam(3) = niter
  iparam(7) = 1

  do
#if defined(R_TCOMPLEX)
    if(arpack%use_parpack) then
#if defined(HAVE_PARPACK)
      call pznaupd  ( mpi_comm, &
            ido, 'I', n, which, nev, tol, resid, ncv, &
            v, ldv, iparam, ipntr, workd, workl, lworkl, &
            rwork, info)

#endif
    else
      call znaupd  ( & 
            ido, 'I', n, which, nev, tol, resid, ncv, &
            v, ldv, iparam, ipntr, workd, workl, lworkl, &
            rwork, info)
    end if

#else 
    if(arpack%use_parpack) then
#if defined(HAVE_PARPACK)
      call pdnaupd  ( mpi_comm, &
            ido, 'I', n, which, nev, tol, resid, ncv, &
            v, ldv, iparam, ipntr, workd, workl, lworkl, &
            info )  
#endif
    else 
      call dnaupd  ( & 
            ido, 'I', n, which, nev, tol, resid, ncv, &
            v, ldv, iparam, ipntr, workd, workl, lworkl, &
            info)
    end if
#endif      
      
    if( abs(ido).ne.1) exit
    
    !!!call av (arpack, ldv, workd(ipntr(1)), workd(ipntr(2))) ! calculate H * psi
    call av (arpack, n, workd(ipntr(1)), workd(ipntr(2))) ! calculate H * psi
    
  end do
  
  !Error Check
  call arpack_check_error('naupd', info)
  
 

#if defined(R_TCOMPLEX) 
  if(arpack%use_parpack) then
#if defined(HAVE_PARPACK) 
    call pzneupd  (mpi_comm,&
          .true., 'A', select, zd, v, ldv, sigma, &
          workev, 'I', n, which, nev, tol, resid, ncv, & 
          v, ldv, iparam, ipntr, workd, workl, lworkl, &
          rwork, info)
          d(:,1)=real(zd(:))
          d(:,2)=aimag(zd(:))
          d(:,3)=M_ZERO
#endif
  else
    call zneupd  (&
          .true., 'A', select, zd, v, ldv, sigma, &
          workev, 'I', n,  which, nev, tol, resid, ncv, & 
          v, ldv, iparam, ipntr, workd, workl, lworkl, &
          rwork, info)
          d(:,1)=real(zd(:))
          d(:,2)=aimag(zd(:))
          d(:,3)=M_ZERO
  end if    
        
#else	
  if(arpack%use_parpack) then
#if defined(HAVE_PARPACK)
    call pdneupd (mpi_comm,&
         .true., 'A', select, d, d(1,2), v, ldv, &
         sigmar, sigmai, workev, 'I', n, which, nev, tol, &
         resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
         lworkl, info )
 
#endif
  else
    call dneupd (&
         .true., 'A', select, d, d(1,2), v, ldv, &
         sigmar, sigmai, workev, 'I', n, which, nev, tol, &
         resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
         lworkl, info )
  end if

#endif

  !Error Check    
  call arpack_check_error('neupd', info) 

  ! This sets niter to the number of matrix-vector operations.
  niter = iparam(9)
  
  ! The number of converged eigenvectors.
  converged =  iparam(5)
    
  do j = 1, converged
!     do i = 1, n
!       psi(i,1) = v(i, j)!/sqrt(gr%mesh%volume_element) 
!     end do
!     do i = n + 1, gr%mesh%np_part
!       psi(i,1) = R_TOTYPE(M_ZERO) 
!     end do
    psi(1:n, 1) = v(1:n, j)
    psi(n+1:gr%mesh%np_part,1) = R_TOTYPE(M_ZERO) 
    
    call states_set_state(st, gr%mesh, j, ik, psi)

    eps_temp = (d(j, 1) + M_zI * d(j, 2)) / arpack%rotation

    st%eigenval(j, ik) = real(eps_temp)
    if(associated(st%zeigenval%Im)) then
      st%zeigenval%Im(j, ik) = aimag(eps_temp)
    end if

    if(abs(workl(ipntr(11)+j-1))< M_EPSILON) then
      diff(j) = M_ZERO
    else
      diff(j) = workl(ipntr(11)+j-1)
    end if
  end do

  !Fill unconverged states with (nice) garbage  
  ! or maybe we should go with whatever we have
  !do j = converged + 1, st%nst
  !  do i = 1, gr%mesh%np
  !    psi(i,1) = R_TOTYPE(M_ONE) 
  !  end do
  !  call states_set_state(st, gr%mesh, j, ik, psi)

  !  st%eigenval(j, ik) = M_HUGE
  !  if(associated(st%zeigenval%Im))then 
  !    st%zeigenval%Im(j, ik) = M_HUGE
  !  end if
  !  diff(j) = M_HUGE
  !end do



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

   POP_SUB(eigen_arpack.eigen_solver_arpack)
contains

  ! ---------------------------------------------------------
  subroutine av (arpack, n, v, w)
    type(eigen_arpack_t), intent(in) :: arpack
    integer,              intent(in) :: n
    R_TYPE,               intent(in) :: v(n)
    R_TYPE,               intent(out):: w(n)
    
    integer :: i, NP, NP_PART
    R_TYPE, allocatable :: psi(:, :), hpsi(:, :)
    
    PUSH_SUB(X(eigen_solver_arpack).av)

    NP = gr%mesh%np
    NP_PART = gr%mesh%np_part

    ASSERT(n == NP .or. n == NP_PART)

    SAFE_ALLOCATE(psi(NP_PART, hm%d%dim))
    SAFE_ALLOCATE(hpsi(NP_PART, hm%d%dim))

!     do i = 1, NP
!       psi(i, 1) = v(i)!/sqrt(gr%mesh%volume_element)
!     end do
!     do i = NP+1, NP_PART
!       psi(i, 1) = M_ZERO
!     end do
    
    psi(1:n,1) = v(1:n)
    psi(n+1:NP_PART, 1) = M_ZERO
    
    call X(hamiltonian_apply) (hm, gr%der, psi, hpsi, 1, ik)
    
    w(1:n) = arpack%rotation * hpsi(1:n,1) ! XXX only works if complex
    
 !    do i = 1, NP
 !       w(i) = hpsi(i, 1)!*sqrt(gr%mesh%volume_element)
 !     end do


    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)

    POP_SUB(X(eigen_solver_arpack).av)
  end subroutine av

end subroutine X(eigen_solver_arpack)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
