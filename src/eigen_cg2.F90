! CONJUGATE-GRADIENTS METHOD.
!
! Method-dependent parameters (optional):
!
! logical, intent(in), optional :: reorder
!       This flag decides wether the eigenvectors/values are ordered before
!       output. By default, they are.
subroutine eigen_solver_cg2(sys, h, st, tol, niter, converged, errorflag, diff, reorder)
  type(system_type), intent(IN)      :: sys
  type(hamiltonian_type), intent(IN) :: h
  type(states_type), intent(inout)   :: st
  real(r8), intent(in)               :: tol
  integer, intent(inout)             :: niter
  integer, intent(out)               :: errorflag, converged
  real(r8), intent(out), optional    :: diff(:,:)
  logical, intent(in), optional      :: reorder

  R_TYPE, allocatable :: h_psi(:,:), g(:,:), g0(:,:), &
       cg(:,:), ppsi(:,:), tmp_wf(:,:)

  R_TYPE :: es(2), a0, b0, gg, gg0, gg1, gamma, theta, norma
  real(r8) :: cg0, e0
  integer  :: ik, np, moved, p, j, iter, i, maxter
  logical  :: reord = .true.

  sub_name = 'eigen_solver_cg2'; call push_sub()

  np = sys%m%np

  if(present(reorder)) reord = reorder
  maxter = niter
  niter = 0
  converged = 0
  moved = 0

  allocate(h_psi(np, st%dim), g(np, st%dim), g0(np, st%dim), &
       cg(0:np, st%dim), ppsi(np, st%dim))
  cg(0,:) = R_TOTYPE(0._r8)

  ! Start of main loop, which runs over all the eigenvectors searched
  ik_loop: do ik = 1, st%nik
    eigenfunction_loop : do p = 1, st%nst
      if(present(diff)) then
        allocate(tmp_wf(np, st%dim))
        tmp_wf = st%R_FUNC(psi)(1:,:, p, ik)
      end if

      ! Orthogonalize starting eigenfunctions to those already calculated...
      call R_FUNC(states_gram_schmidt)(p, sys%m, st%dim, &
           st%R_FUNC(psi)(:,:, 1:p, ik), start=p)
      
      ! Calculate starting gradient: |hpsi> = H|psi>
      call R_FUNC(Hpsi)(h, sys, ik, st%R_FUNC(psi)(:,:, p, ik) , h_psi)

      ! Calculates starting eigenvalue: e(p) = <psi(p)|H|psi>
      st%eigenval(p, ik) = R_REAL(R_FUNC(states_dotp) (sys%m, st%dim, 
           st%R_FUNC(psi)(1:,:, p, ik), h_psi))

      ! Starts iteration for this band
      iter_loop: do iter = 1, maxter

        ! if approximate inverse preconditioner....
        !call pre(hpsi%val   , g%val   ) 
        !call pre(psi(m)%val , ppsi%val)
        g = h_psi
        ppsi = st%R_FUNC(psi)(1:,:, p, ik)
        
        es(1) = R_FUNC(states_dotp) (sys%m, st%dim, st%R_FUNC(psi)(1:,:, p, ik), g)
        es(2) = R_FUNC(states_dotp) (sys%m, st%dim, st%R_FUNC(psi)(1:,:, p, ik), ppsi)
        es(1) = es(1)/es(2)
        g = g - es(1)*ppsi
        
        ! Orthogonalize to lowest eigenvalues (already calculated)
        do j = 1, p - 1
          a0 = R_FUNC(states_dotp) (sys%m, st%dim, st%R_FUNC(psi)(1:,:, j, ik), g)
          g(:,:) = g(:,:) - a0 * st%R_FUNC(psi)(1:,:, j, ik)
        end do
        
        if(iter .ne. 1) gg1 = R_FUNC(states_dotp) (sys%m, st%dim, g, g0)
        
        ! Approximate inverse preconditioner...
        !call ipre(g, g0)
        g0 = g
        
        gg = R_FUNC(states_dotp) (sys%m, st%dim, g, g0)
        
        ! Starting or following iterations...
        if(iter .eq. 1) then
          gg0 = gg
          cg(1:, :) = g(:,:)
        else
          !gamma = gg/gg0        ! (Fletcher-Reeves)
          gamma = (gg - gg1)/gg0   ! (Polack-Ribiere)
          gg0 = gg
          cg = gamma*cg
          cg(1:,:) = cg(1:,:) + g(:,:)
          norma = gamma*cg0*sin(theta)
          cg(1:,:) = cg(1:,:) - norma * st%R_FUNC(psi)(1:,:, p, ik)
        end if
        
        ! cg contains now the conjugate gradient
        cg0 = R_FUNC(states_nrm2) (sys%m, st%dim, cg(1:,:))
        call R_FUNC(Hpsi) (h, sys, ik, cg, ppsi)
        
        ! Minimization.
        a0 = R_FUNC(states_dotp) (sys%m, st%dim, st%R_FUNC(psi)(1:,:, p, ik), ppsi)
        a0 = 2.0_r8 * a0 / cg0
        b0 = R_FUNC(states_dotp) (sys%m, st%dim, cg(1:,:), ppsi)
        b0 = b0/cg0**2
        e0 = st%eigenval(p, ik)
        theta = atan(a0/(e0 - b0))/2.0_r8
        es(1) = ((e0-b0)*cos(2.0_r8*theta) + a0*sin(2.0_r8*theta) + e0 + b0) / 2.0_r8
        es(2) =(-(e0-b0)*cos(2.0_r8*theta) - a0*sin(2.0_r8*theta) + e0 + b0) / 2.0_r8
        
        ! Choose the minimum solutions. 
        if (R_REAL(es(2)) < R_REAL(es(1))) then
          theta = theta + M_PI/2.0_r8
        end if
        st%eigenval(p, ik) = min(R_REAL(es(1)), R_REAL(es(2)))
        
        ! Upgrade psi...
        a0 = cos(theta)
        b0 = sin(theta)/cg0
        st%R_FUNC(psi)(:,:, p, ik) = a0*st%R_FUNC(psi)(:,:, p, ik) + b0*cg(:,:)
        
        ! Test convergence.
        if(abs(st%eigenval(p, ik) - e0) < tol) then
          converged = converged + 1
          exit iter_loop
        end if
        
        ! Calculate H|psi>
        h_psi = a0*h_psi + b0*ppsi
        
      end do iter_loop

      niter = niter + iter + 1

      if(present(diff)) then
        tmp_wf = st%R_FUNC(psi)(1:np,:, p, ik) - tmp_wf
        diff(p, ik) = R_FUNC(mesh_nrm2) (sys%m, tmp_wf)
        deallocate(tmp_wf)
      end if

      ! Reordering... (this should be improved)
      if(p>1 .and. reord) then
        if(st%eigenval(p, ik) - st%eigenval(p-1, ik) .lt. -2.0_r8*tol) then
          do i = p - 2, 1, -1
            if (st%eigenval(p, ik) - st%eigenval(i, ik) .gt. 2.0_r8*tol) exit
          end do
          i = i + 1
          moved = moved + 1
          e0 = st%eigenval(p, ik)
          cg = st%R_FUNC(psi)(:,:, p, ik)
          do j = p, i + 1 , -1
            st%eigenval(j, ik) = st%eigenval(j-1, ik)
            st%R_FUNC(psi)(:,:, j, ik) = st%R_FUNC(psi)(:,:, j-1, ik)
          end do
          st%eigenval(i, ik) = e0
          st%R_FUNC(psi)(:,:, i, ik) = cg
        end if
      end if
      ! End of reordering
    end do eigenfunction_loop
  end do ik_loop

  ! Deallocation of variables
  deallocate(h_psi, g, g0, cg, ppsi)
  
  call pop_sub()
  return
end subroutine eigen_solver_cg2
