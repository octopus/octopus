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

!!! This file handles the evaluation of the OEP potential, in the KLI or full OEP
!!! as described in S. Kuemmel and J. Perdew, PRL 90, 043004 (2003)

!!! This file has to be outside the module xc, for it requires the Hpsi

subroutine X(h_xc_oep)(xcs, m, h, sys, st, vxc, ex, ec)
  type(xc_type),     intent(in)    :: xcs
  type(mesh_type),   intent(in)    :: m
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in)    :: sys
  type(states_type), intent(inout) :: st
  real(r8),          intent(inout) :: vxc(m%np, st%nspin), ex, ec
  
  type(xc_oep_type) :: oep
  real(r8) :: e
  integer :: is, ist, ixc, ifunc

  call push_sub('h_xc_oep')

  ! this routine is only prepared for finite systems, and ispin = 1, 2
  if(st%ispin > 2 .or. st%nik>st%ispin) then
    message(1) = "OEP only works for finite systems and collinear spin!"
    call write_fatal(1)
  end if

  ! initialize oep structure
  allocate(oep%eigen_type(st%nst), oep%eigen_index(st%nst))
  allocate(oep%vxc(m%np), oep%lxc(m%np, st%nst), oep%uxc_bar(st%nst))

  ! obtain the spin factors
  call xc_oep_SpinFactor(oep, st%nspin)

  spin: do is = 1, min(st%nspin, 2)
    oep%lxc = M_ZERO
    e       = M_ZERO

    ! get lxc
    functl_loop: do ixc = 0, N_X_FUNCTL+N_C_FUNCTL-1
      if(.not.btest(xcs%functl, ixc)) cycle
      ifunc = ibset(0, ixc)
      if(.not.( &
           ifunc == X_FUNC_OEP_X    .or. &
           ifunc == X_FUNC_OEP_SIC  .or. &
           ifunc == C_FUNC_OEP_SIC)) cycle
      
      select case(ifunc)
      case(X_FUNC_OEP_X)
        call X(oep_x) (m, st, is, oep, e)
        
      case(X_FUNC_OEP_SIC)
        call X(oep_x_sic) (xcs, m, st, is, oep, e)
      case(C_FUNC_OEP_SIC)
        call X(oep_c_sic) (xcs, m, st, is, oep, e)
      end select
      
      if(ixc < N_X_FUNCTL) then
        ex = ex + e
      else
        ec = ec + e
      end if
    end do functl_loop
  
    ! get the HOMO state
    call xc_oep_AnalizeEigen(oep, st, is)

    ! calculate uxc_bar for the occupied states
    do ist = 1, st%nst
      oep%uxc_bar(ist) = sum(R_REAL(st%X(psi)(:, 1, ist, is) * oep%lxc(:, ist)))*m%vol_pp
    end do

    ! solve the KLI equation
    oep%vxc = M_ZERO
    call X(xc_KLI_solve) (m, st, is, oep, xcs%oep_level)
    if(xcs%oep_level == 2) call X(h_xc_oep_solve)(m, h, sys, st, is, vxc(:,is), oep)

    vxc(:, is) = vxc(:, is) + oep%vxc(:)
  end do spin

  deallocate(oep%eigen_type, oep%eigen_index)
  deallocate(oep%vxc, oep%lxc, oep%uxc_bar)

  call pop_sub()
end subroutine X(h_xc_OEP)

subroutine X(h_xc_oep_solve) (m, h, sys, st, is, vxc, oep)
  type(mesh_type),   intent(in)    :: m
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in)    :: sys
  type(states_type), intent(in)    :: st
  integer,           intent(in)    :: is
  real(r8),          intent(inout) :: vxc(m%np)
  type(xc_oep_type), intent(inout) :: oep

  integer :: iter, i, ist
  real(r8) :: vxc_bar
  real(r8), allocatable :: s(:), vxc_old(:)
  R_TYPE, allocatable :: b(:), psi(:,:)

  allocate(b(m%np), s(m%np), psi(m%np, 1), vxc_old(m%np))

  vxc_old = vxc(:)
  do iter = 1, 30
    ! fix xc potential (needed for Hpsi)
    vxc(:) = vxc_old(:) + oep%vxc(:)

    ! iteration ver all states
    s = M_ZERO
    do ist = 1, st%nst
      ! evaluate right-hand side
      vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:))*m%vol_pp
      b(:) = (oep%vxc(:)*R_CONJ(st%X(psi)(:, 1, ist, is))  - oep%lxc(:, ist))  &
           - (vxc_bar - oep%uxc_bar(ist))*R_CONJ(st%X(psi)(:, 1, ist, is))
      
      ! initialize psi to something
      !call X(states_random)(m, psi(:,1))
      psi(:,1) = b(:)
      
      ! and we now solve the equation [h-eps_i] psi_i = b_i
      call get_psi()
      
      ! calculate this funny function s
      s(:) = s(:) + M_TWO*R_REAL(psi(:,1)*st%X(psi)(:, 1, ist, is))
    end do
    print *, iter, "s = ", dmf_nrm2(m, s)
    oep%vxc(:) = oep%vxc(:) + 2._r8*s(:)

    do ist = 1, st%nst
      if(oep%eigen_type(ist) == 2) then
        vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:))*m%vol_pp
        oep%vxc(:) = oep%vxc(:) - (vxc_bar - oep%uxc_bar(ist))
      end if
    end do

  end do

  vxc(:) = vxc_old(:)
  deallocate(b, s, psi, vxc_old)

contains
  ! This subroutine uses conjugated gradients to solve the linear equation
  ! (H - eps_i) psi_i = b_i, with psi_i orthogonal to psi_{KS i}
  ! as in H.R. Schwarz, "Numerische Mathematik", pag. 603
  ! Thanks to Heiko Appel for the initial code
  subroutine get_psi()
    integer :: iter
    R_TYPE, allocatable :: res(:,:), p(:,:), x(:,:), z(:,:), tmp(:,:)
    R_TYPE :: r, ek, qk, sprod, spold, norm
    
    allocate(res(m%np, 1), p(m%np, 1), x(m%np, 1), z(m%np, 1), tmp(m%np, 1))
    
    ! Orthogonalize starting psi to phi
    r = X(states_dotp) (m, st%dim, psi, st%X(psi)(:, :, ist, is))
    psi = psi - r*st%X(psi)(:, :, ist, is)

    ! Calculate starting gradient: |hpsi> = (H-eps_i)|psi> + b
    call X(Hpsi)(h, m, psi, res, sys, 1)
    res(:,1) = res(:,1) - st%eigenval(ist, 1)*psi(:,1) + b(:)

    ! orthogonalize direction to phi
    r = X(states_dotp) (m, st%dim, res, st%X(psi)(:, :, ist, is))
    res = res - r*st%X(psi)(:, :, ist, is)

    p = -res
    spold = X(states_nrm2)(m, st%dim, res)**2

    iter_loop: do iter = 1, 50
      ! verify
      call X(Hpsi)(h, m, psi, tmp, sys, 1)
      tmp(:,1) = tmp(:,1) - st%eigenval(ist, 1)*psi(:,1) + b(:)
      if(X(states_nrm2)(m, st%dim, tmp).lt.1e-6_8) exit iter_loop

      if(iter > 1) then
        r = X(states_nrm2)(m, st%dim, res)**2
        ek = r/spold
        spold = r

        p  = -res + ek*p
      end if

      call X(Hpsi)(h, m, p, z, sys, 1)
      z(:,1) = z(:,1) - st%eigenval(ist, 1)*p(:,1)
      r = X(states_dotp) (m, st%dim, z, st%X(psi)(:, :, ist, is))
      z = z - r*st%X(psi)(:, :, ist, is)

      qk    = spold/X(states_dotp)(m, st%dim, p, z)
      psi   = psi + qk*p
      res   = res + qk*z

    end do iter_loop

    deallocate(res, p, x, z, tmp)
  end subroutine get_psi
end subroutine X(h_xc_oep_solve)
