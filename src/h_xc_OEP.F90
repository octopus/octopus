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
!!! This is why it needs the xc_functl module. I prefer to put it here since
!!! the rest of the hamiltonian module does not know about the gory details
!!! of how xc is defined and calculated.

subroutine X(h_xc_oep)(xcs, m, f_der, h, st, vxc, ex, ec)
  use xc_functl
  type(xc_type),          intent(IN)    :: xcs
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(hamiltonian_type), intent(IN)    :: h
  type(states_type),      intent(inout) :: st
  FLOAT,                  intent(inout) :: vxc(m%np, st%d%nspin)
  FLOAT,                  intent(inout) :: ex, ec
  
  type(xc_oep_type) :: oep
  FLOAT :: e
  integer :: is, ist, ixc, ifunc, ierr

  if(xcs%oep_level.eq.XC_OEP_NONE) return

  call push_sub('h_xc_oep')

  ! this routine is only prepared for finite systems, and ispin = 1, 2
  if(st%d%ispin > SPIN_POLARIZED .or. st%d%nik>st%d%ispin) then
    message(1) = "OEP only works for finite systems and collinear spin!"
    call write_fatal(1)
  end if


  ! initialize oep structure
  allocate(oep%eigen_type(st%nst))
  allocate(oep%eigen_index(st%nst))
  allocate(oep%lxc(m%np, st%st_start:st%st_end))
  allocate(oep%uxc_bar(st%nst))
  allocate(oep%vxc(m%np))
  

  ! obtain the spin factors
  call xc_oep_SpinFactor(oep, st%d%nspin)

  ! this part handles the (pure) orbital functionals
  spin: do is = 1, min(st%d%nspin, 2)
    oep%lxc = M_ZERO

    ! get lxc
    functl_loop: do ixc = 1, 2
      if(xcs%functl(ixc)%family.ne.XC_FAMILY_OEP) cycle

      e = M_ZERO
      select case(xcs%functl(ixc)%id)
      case(XC_OEP_X)
        call X(oep_x) (m, f_der, st, is, oep, e)
        ex = ex + e
      end select
    end do functl_loop
  
    ! SIC a la PZ is handled here
    if(xcs%sic_correction.ne.0) then
      call X(oep_sic) (xcs, m, f_der, st, is, oep, vxc, ex, ec)
    end if
    
    ! get the HOMO state
    call xc_oep_AnalizeEigen(oep, st, is)

    ! calculate uxc_bar for the occupied states
    do ist = st%st_start, st%st_end
      oep%uxc_bar(ist) = sum(R_REAL(st%X(psi)(:, 1, ist, is) * oep%lxc(:, ist))*m%vol_pp(:))
    end do
#if defined(HAVE_MPI)
    if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      call mpi_barrier(mpi_comm_world, ierr)
      do ist = 1, st%nst
         call mpi_bcast(oep%uxc_bar(ist), 1, MPI_FLOAT, st%node(ist), MPI_COMM_WORLD, ierr)
      enddo
    endif
#endif
    
    ! solve the KLI equation
    oep%vxc = M_ZERO
    call X(xc_KLI_solve) (m, st, is, oep, xcs%oep_level)

    ! if asked, solve the full OEP equation
    if(xcs%oep_level == XC_OEP_FULL) then
      call X(h_xc_oep_solve)(m, f_der, h, st, is, vxc(:,is), oep)
    end if

    vxc(:, is) = vxc(:, is) + oep%vxc(:)
  end do spin

  deallocate(oep%eigen_type, oep%eigen_index, oep%vxc, oep%lxc, oep%uxc_bar)
  call pop_sub()
end subroutine X(h_xc_OEP)


subroutine X(h_xc_oep_solve) (m, f_der, h, st, is, vxc, oep)
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(hamiltonian_type), intent(IN)    :: h
  type(states_type),      intent(IN)    :: st
  integer,                intent(in)    :: is
  FLOAT,                  intent(inout) :: vxc(m%np)
  type(xc_oep_type),      intent(inout) :: oep

  integer :: iter, ist
  FLOAT :: vxc_bar
  FLOAT, allocatable :: s(:), vxc_old(:)
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
      vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:) * m%vol_pp(:))
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
    oep%vxc(:) = oep%vxc(:) + M_TWO*s(:)

    do ist = 1, st%nst
      if(oep%eigen_type(ist) == 2) then
        vxc_bar = sum(R_ABS(st%X(psi)(:, 1, ist, is))**2 * oep%vxc(:) * m%vol_pp(:))
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
    R_TYPE :: r, ek, qk, spold
    
    allocate(res(m%np, 1), p(m%np, 1), x(m%np, 1), z(m%np, 1), tmp(m%np, 1))
    
    ! Orthogonalize starting psi to phi
    r = X(states_dotp) (m, st%d%dim, psi, st%X(psi)(:, :, ist, is))
    psi = psi - r*st%X(psi)(:, :, ist, is)

    ! Calculate starting gradient: |hpsi> = (H-eps_i)|psi> + b
    call X(Hpsi)(h, m, f_der, psi, res, 1)
    res(:,1) = res(:,1) - st%eigenval(ist, 1)*psi(:,1) + b(:)

    ! orthogonalize direction to phi
    r = X(states_dotp) (m, st%d%dim, res, st%X(psi)(:, :, ist, is))
    res = res - r*st%X(psi)(:, :, ist, is)

    p = -res
    spold = X(states_nrm2)(m, st%d%dim, res)**2

    iter_loop: do iter = 1, 50
      ! verify
      call X(Hpsi)(h, m, f_der, psi, tmp, 1)
      tmp(:,1) = tmp(:,1) - st%eigenval(ist, 1)*psi(:,1) + b(:)
      if(X(states_nrm2)(m, st%d%dim, tmp).lt.1e-6_8) exit iter_loop

      if(iter > 1) then
        r = X(states_nrm2)(m, st%d%dim, res)**2
        ek = r/spold
        spold = r

        p  = -res + ek*p
      end if

      call X(Hpsi)(h, m, f_der, p, z, 1)
      z(:,1) = z(:,1) - st%eigenval(ist, 1)*p(:,1)
      r = X(states_dotp) (m, st%d%dim, z, st%X(psi)(:, :, ist, is))
      z = z - r*st%X(psi)(:, :, ist, is)

      qk    = spold/X(states_dotp)(m, st%d%dim, p, z)
      psi   = psi + qk*p
      res   = res + qk*z

    end do iter_loop

    deallocate(res, p, x, z, tmp)
  end subroutine get_psi
end subroutine X(h_xc_oep_solve)
