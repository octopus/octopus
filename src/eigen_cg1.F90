! This subroutine performs NCG cycles of conjugated-gradients minimization    !
subroutine eigen_solver_cg1(ncg, sys, h, st, diff)
  integer, intent(in) :: ncg
  type(system_type), intent(IN) :: sys
  type(hamiltonian_type), intent(IN) :: h
  type(states_type), intent(inout) :: st
  real(r8), intent(out), optional :: diff(st%nst, st%nik)

  integer :: ik, iter, p, q, i, np, idim
  real(r8) :: r
  R_TYPE :: xkHxk, xkxk, Rk, pkpk, pkHxk, pkHpki,     &
       uk, alpha, Ak, Bk, Ck, pkHpk, s, xkpk, gkgk
  R_TYPE, allocatable :: xk(:,:), hxk(:,:), gk(:,:), pk(:,:), tmp_wf(:,:)

  sub_name = 'eigen_solver_cg1'; call push_sub()

  allocate(xk(0:sys%m%np, st%dim), hxk(sys%m%np, st%dim), &
       gk(sys%m%np, st%dim), pk(0:sys%m%np, st%dim))

  np = sys%m%np ! shortcuts

  xkHxk=0.0_r8; xkxk=0.0_r8; Rk=0.0_r8
  gkgk=0.0_r8;  xkpk=0.0_r8; pkpk=0.0_r8;  pkHxk=0.0_r8
  pkHpk=0.0_r8; uk=0.0_r8;   alpha=0.0_r8; Ak=0.0_r8
  Bk=0.0_r8;    Ck=0.0_r8
  xk=0.0_r8;    hxk=0.0_r8;  gk=0.0_r8;    pk=0.0_r8

! The cycles are performed
  do p = 1, st%nst
    do ik = 1, st%nik
      if(present(diff)) then
        allocate(tmp_wf(np, st%dim))
        tmp_wf(1:np,:) = st%R_FUNC(psi)(1:np,:, p, ik)
      end if

      ! orthogonalize state p to all previous ones
      call R_FUNC(states_gram_schmidt)(p, sys%m, st%dim, &
           st%R_FUNC(psi)(:,:, 1:p, ik), start = p)

      xk(1:np,:) = st%R_FUNC(psi)(1:np,:, p, ik)
      call R_FUNC(Hpsi) (h, sys, ik, xk, hxk)

      xkHxk      = R_REAL(R_FUNC(states_dotp) (sys%m, st%dim, xk(1:np,:), hxk))
      xkxk       = 1._r8
      Rk         = xkHxk
      do iter = 1, ncg
        gk = 2._r8*(hxk - Rk*xk(1:np,:))/xkxk
        do q = 1, p - 1
          s = R_FUNC(states_dotp) (sys%m, st%dim, st%R_FUNC(psi)(1:np,:, q, ik), gk)
          do idim = 1, st%dim
            call R_FUNC(axpy) (np, -s, st%R_FUNC(psi)(1:np, idim, q, ik), 1, gk(1:np, idim), 1)
          end do
        end do
        r = R_FUNC(states_nrm2) (sys%m, st%dim, gk)**2
        
        select case (iter)
        case(1)
          pk(1:np,:) = -gk
        case default
          uk = r/gkgk
          pk(1:np,:) = -gk + uk*pk(1:np,:)
        end select
        
        gkgk = r
        xkpk = R_FUNC(states_dotp) (sys%m, st%dim, xk(1:np,:), pk(1:np,:))
        pkpk = R_FUNC(states_nrm2) (sys%m, st%dim, pk(1:np,:))**2
        pkHxk= R_FUNC(states_dotp) (sys%m, st%dim, pk(1:np,:), hxk)
        call R_FUNC(Hpsi) (h, sys, ik, pk, gk)
        pkHpk = R_FUNC(states_dotp) (sys%m, st%dim, pk(1:np,:), gk)

        Ak = pkHpk*xkpk - pkHxk*pkpk
        Bk = pkHpk*xkxk - xkHxk*pkpk
        Ck = pkHxk*xkxk - xkHxk*xkpk
        alpha = (-Bk + sqrt(Bk*Bk - 4._r8*Ak*Ck))/(2._r8*Ak)

        xk = xk + alpha*pk
        hxk = hxk + alpha*gk
        xkxk = R_FUNC(states_nrm2) (sys%m, st%dim, xk(1:np,:))**2
        xkHxk= R_FUNC(states_dotp) (sys%m, st%dim, xk(1:np,:), hxk)
        Rk = xkHxk/xkxk
        
        s = sqrt(xkxk)
        xk(1:np,:) = xk(1:np,:)/s
        st%R_FUNC(psi)(1:np,:, p, ik) = xk(1:np,:)
      end do

      if(present(diff)) then
        tmp_wf = st%R_FUNC(psi)(1:np,:, p, ik) - tmp_wf
        diff(p, ik) = R_FUNC(mesh_nrm2) (sys%m, tmp_wf)
        deallocate(tmp_wf)
      end if
    enddo
  enddo
     
  deallocate(xk, hxk, gk, pk)

  call pop_sub()
  return
end subroutine eigen_solver_cg1
