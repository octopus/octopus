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

! calculates the eigenvalues of the real orbitals
subroutine X(hamiltonian_eigenval)(h, m, st)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  type(states_type),      intent(inout) :: st

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('hamiltonian_eigenval')
  allocate(Hpsi(m%np, st%d%dim))

  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      call X(hpsi) (h, m, st%X(psi)(:, :, ist, ik), hpsi, ik)
      e = X(states_dotp)(m, st%d%dim, st%X(psi)(1:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
end subroutine X(hamiltonian_eigenval)

subroutine X(Hpsi) (h, m, psi, hpsi, ik, t)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  integer,                intent(in)  :: ik
  R_TYPE,                 intent(IN)  :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  FLOAT, optional,        intent(in)  :: t

  call push_sub('Hpsi')

  call X(kinetic) (h, m, psi, hpsi, ik)
  call X(vlpsi)   (h, m, psi, hpsi, ik)
  if(h%ep%nvnl > 0) call X(vnlpsi)  (h, m, psi, hpsi, ik)

  ! Relativistic corrections...
  select case(h%reltype)
  case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
  case(SPIN_ORBIT)
    call zso (h, m, psi, hpsi, h%d%dim, ik)
#endif
  case default
    message(1) = 'Error: Internal.'
    call write_fatal(1)
  end select

  if(present(t)) then
    if (h%d%cdft) then
      message(1) = "TDCDFT not yet implemented"
      call write_fatal(1)
    end if
    call X(vlasers)  (h, m, psi, hpsi, t)
    call X(vborders) (h, m, psi, hpsi)
  elseif (h%d%cdft) then
    call X(current_extra_terms) (h, m, psi, hpsi, ik)
  end if

  call pop_sub()
end subroutine X(Hpsi)

#if defined(R_TCOMPLEX)
subroutine X(magnus) (h, m, psi, hpsi, ik, vmagnus)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  integer,                intent(in)  :: ik
  R_TYPE,                 intent(IN)  :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  FLOAT,                  intent(IN)  :: vmagnus(m%np, h%d%nspin, 2)

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim
  ! We will assume, for the moment, no spinors.

  call push_sub('magnus')

  allocate(auxpsi(m%np, h%d%dim), aux2psi(m%np, h%d%dim))

  call X(kinetic) (h, m, psi, hpsi, ik)

  auxpsi = hpsi
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, m, psi, auxpsi, ik)
  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(:, 1) = hpsi(:, 1) -  M_zI*vmagnus(:, 1, 1)*auxpsi(:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      hpsi(:, 1) = hpsi(:, 1) - M_zI*vmagnus(:, 1, 1)*auxpsi(:, 1)
    else
      hpsi(:, 1) = hpsi(:, 1) - M_zI*vmagnus(:, 2, 1)*auxpsi(:, 1)
    end if
  end select

  select case(h%d%ispin)
  case(UNPOLARIZED)
    auxpsi(:, 1) = vmagnus(:, 1, 1)*psi(:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      auxpsi(:, 1) = vmagnus(:, 1, 1)*psi(:, 1)
    else
      auxpsi(:, 1) = vmagnus(:, 2, 1) *psi(:, 1)
    end if
  end select
  call X(kinetic) (h, m, auxpsi, aux2psi, ik)
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, m, auxpsi, aux2psi, ik)
  hpsi(:, 1) = hpsi(:, 1) + M_zI*aux2psi(:, 1)

  do idim = 1, h%d%dim
      hpsi(:, idim) = hpsi(:, idim) + h%Vpsl(:)*psi(:,idim)
  end do

  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(:, 1) = hpsi(:, 1) + vmagnus(:, 1, 2)*psi(:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then
      hpsi(:, 1) = hpsi(:, 1) + vmagnus(:, 1, 2)*psi(1:, 1)
    else
      hpsi(:, 1) = hpsi(:, 1) + vmagnus(:, 2, 2)*psi(1:, 1)
    end if
  end select
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, m, psi, Hpsi, ik)
  call X(vborders) (h, m, psi, hpsi)

  deallocate(auxpsi, aux2psi)
  call pop_sub()
end subroutine X(magnus)
#endif

subroutine X(kinetic) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  R_TYPE,                 intent(IN)  :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  integer :: ik

  integer :: idim
#if defined(COMPLEX_WFNS)
  integer :: i
  R_TYPE, allocatable :: grad(:,:)
  FLOAT :: k2
#endif  

  call push_sub('kinetic')

  if(conf%periodic_dim>0) then
#if defined(COMPLEX_WFNS)
    allocate(grad(3, m%np))
    k2 = sum(h%d%kpoints(:, ik)**2)
    do idim = 1, h%d%dim
      call X(f_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
      call X(f_gradient)  (m, psi(:, idim), grad(:, :))
      do i = 1, m%np
        Hpsi(i, idim) = -M_HALF*(Hpsi(i, idim) &
             + M_TWO*M_zI*sum(h%d%kpoints(:, ik)*grad(:, i)) &
             - k2*psi(i, idim))
      end do
    end do
    deallocate(grad)
#else
    message(1) = "Real wavefunction for ground state not yet implemented for polymers:"
    message(2) = "Reconfigure with --enable-complex, and remake"
    call write_fatal(2)
#endif
  
  else
    do idim = 1, h%d%dim
      call X(f_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
    end do
    call lalg_scal(m%np, h%d%dim, R_TOTYPE(-M_HALF), Hpsi)
  end if


  call pop_sub()
end subroutine X(kinetic)

subroutine X(current_extra_terms) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  R_TYPE,                 intent(IN)  :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  integer,                intent(in)  :: ik

  integer :: idim, k, ispin
  FLOAT :: b(conf%dim), a(conf%dim), r(3)
  FLOAT, allocatable :: div(:)
  R_TYPE, allocatable :: grad(:,:)

  call push_sub('current_extra_terms')

  b = M_ZERO
  if (associated(h%ep%b)) then
    b = b + h%ep%b/P_C
  end if

  do idim = 1, conf%dim
    select case(h%d%ispin)
    case(UNPOLARIZED)
      hpsi(:, 1) = hpsi(:, 1) + h%ahxc(idim, :, 1)*h%ahxc(idim, :, 1)*psi(1:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(:, 1) = hpsi(:, 1) + h%ahxc(idim, :, 1)*h%ahxc(idim, :, 1)*psi(1:, 1)
      else
        hpsi(:, 1) = hpsi(:, 1) + h%ahxc(idim, :, 2)*h%ahxc(idim, :, 2)*psi(1:, 1)
      end if
    case (SPINORS)
      ! not implemented yet
    end select
  end do

  allocate(div(m%np))
  select case (h%d%ispin)
  case(UNPOLARIZED)
    call df_divergence(m, h%ahxc(:, :, 1), div)
    hpsi(:, 1) = hpsi(:, 1) - M_zI * div*psi(1:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then ! we have a spin down
      call df_divergence(m, h%ahxc(:, :, 1), div)
    else
      call df_divergence(m, h%ahxc(:, :, 2), div)
    end if
    hpsi(:, 1) = hpsi(:, 1) - M_zI * div*psi(1:, 1)
  case(SPINORS)
    ! not implemented yet
  end select
  deallocate(div)

  allocate(grad(3, m%np))
  select case (h%d%ispin)
  case(UNPOLARIZED)
    call X(f_gradient)(m, psi(:, 1), grad)
    do k = 1, m%np
      hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(:, k, 1), grad(:, k)) 
    end do
  case(SPIN_POLARIZED)
    call X(f_gradient)(m, psi(:, 1), grad)
    do k = 1, m%np
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(:, k, 1), grad(:, k))
      else
        hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(:, k, 2), grad(:, k))
      end if
    end do
  case(SPINORS)
    ! not implemented yet
  end select

  if (.not.all(b == M_ZERO)) then
    do k = 1, m%np
      call mesh_xyz(m, k, r)
      a = -M_HALF*(/r(2)*b(3) - r(3)*b(2), r(3)*b(1) - r(1)*b(3), r(1)*b(2) - r(2)*b(1)/)

      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(a, a)*psi(k, :)

      select case(h%d%ispin)
      case(UNPOLARIZED)
        hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(a, h%ahxc(:, k, 1))*psi(k, 1)
      case(SPIN_POLARIZED)
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(a, h%ahxc(:, k, 1))*psi(k, 1)
        else
          hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(a, h%ahxc(:, k, 2))*psi(k, 1)
        end if
      case (SPINORS)
        ! not implemented yet
      end select

      select case(h%d%ispin)
      case(UNPOLARIZED)
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a, grad(:, k))
      case(SPIN_POLARIZED)
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a, h%ahxc(:, k, 1))*psi(k, 1)
        else
          hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(a, h%ahxc(:, k, 2))*psi(k, 1)
        end if
      case (SPINORS)
        ! not implemented yet
      end select
 
    end do
  end if
  deallocate(grad)

  call pop_sub()
end subroutine X(current_extra_terms)

subroutine X(vnlpsi) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  integer,                intent(in)    :: ik

  integer :: idim, ikbc, jkbc, ivnl
  R_TYPE :: uVpsi
  R_TYPE, allocatable :: lpsi(:), lHpsi(:)
  type(nonlocal_op), pointer :: nlop
  R_TYPE, allocatable :: tmp(:)

  call push_sub('vnlpsi')

  do ivnl = 1, h%ep%nvnl
    nlop => h%ep%vnl(ivnl)
    allocate(tmp(nlop%n))
    
    do_dim: do idim = 1, h%d%dim
      allocate(lpsi(nlop%n), lhpsi(nlop%n))
      if (conf%periodic_dim==0) then
        lpsi(:) = psi(nlop%jxyz(:), idim)
      else
        lpsi(:) = nlop%phases(:,ik)*psi(nlop%jxyz(:), idim)
      end if
      lHpsi(:) = M_z0
    
      do ikbc = 1, nlop%c
        do jkbc = 1, nlop%c
          tmp   = R_TOTYPE(nlop%uv(:, ikbc))
          uvpsi = lalg_dot(nlop%n, tmp, lpsi)*m%vol_pp*nlop%uvu(ikbc, jkbc)
          if (conf%periodic_dim==0) then
            tmp = R_TOTYPE(nlop%uv(:, jkbc))
            call lalg_axpy(nlop%n, uvpsi, tmp, lHpsi)
          else
#           ifdef R_TCOMPLEX
              tmp = R_CONJ(nlop%phases(:, ik)*nlop%uv(:, jkbc))
              call lalg_axpy(nlop%n, uvpsi, tmp, lHpsi)
#           else
              message(1) = "Real wavefunction for ground state not yet implemented for polymers:"
              message(2) = "Reconfigure with --enable-complex, and remake"
              call write_fatal(2)
#           endif
          end if
        end do
      end do
      hpsi(nlop%jxyz(:), idim) = hpsi(nlop%jxyz(:), idim) + lhpsi(:)
      deallocate(lpsi, lhpsi)
    end do do_dim
    
    deallocate(tmp)
  end do

  call pop_sub()
end subroutine X(vnlpsi)

subroutine X(vlpsi) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)

  integer :: idim

  call push_sub('vlpsi')

  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%vpsl(:))*psi(:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then ! we have a spin down
      hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%vpsl(:))*psi(:, 1)
    else
      hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 2) + h%vpsl(:))*psi(:, 1)
    end if
  case(SPINORS)
    hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%vpsl(:))*psi(:, 1) + &
                 (h%vhxc(:, 3) + M_zI*h%vhxc(:, 4))*psi(:, 2)
    hpsi(:, 2) = hpsi(:, 2) + (h%vhxc(:, 2) + h%vpsl(:))*psi(:, 2) + &
                 (h%vhxc(:, 3) - M_zI*h%vhxc(:, 4))*psi(:, 1)
  end select

  call pop_sub()
end subroutine X(vlpsi)

subroutine X(vlasers) (h, m, psi, hpsi, t)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)

  FLOAT, intent(in) :: t

  integer :: k, idim
  FLOAT :: x(3), f(3)
  R_TYPE, allocatable :: grad(:,:)

  call push_sub('vlasers')

  if(h%ep%no_lasers > 0) then
    select case(h%gauge)
    case(1) ! length gauge
      call epot_laser_field(h%ep, t, f)

      do k = 1, m%np
        call mesh_xyz(m, k, x)
        hpsi(k,:) = hpsi(k,:) + sum(x*f) * psi(k,:)
      end do

    case(2) ! velocity gauge
      call epot_laser_vector_field(h%ep, t, f)
      allocate(grad(3, m%np))
      do idim = 1, h%d%dim
        call X(f_gradient)(m, psi(:, idim), grad)
        do k = 1, m%np
          hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
               sum(f**2)/M_TWO * psi(k, idim)
        end do
      end do
      deallocate(grad)
    end select
  end if

  call pop_sub()
end subroutine X(vlasers)

subroutine X(vborders) (h, m, psi, hpsi)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)

  integer :: idim

  call push_sub('vborders')

  if(h%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, h%d%dim
       hpsi(:, idim) = hpsi(:, idim) + M_zI*h%ab_pot(:)*psi(1:, idim)
    end do
  end if

  call pop_sub()
end subroutine X(vborders)

subroutine X(h_calc_vhxc)(h, m, st, calc_eigenval)
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type),        intent(IN)    :: m
  type(states_type),      intent(inout) :: st
  logical,      optional, intent(in)    :: calc_eigenval

  integer :: is, idim
  FLOAT, allocatable :: rho_aux(:), vhartree(:), j_aux(:,:), ahartree(:,:)
  CMPLX, allocatable :: ztmp1(:), ztmp2(:)

  ! next 2 lines are for an RPA calculation
  !FLOAT, save, pointer ::  RPA_Vhxc(:, :)
  !logical, save :: RPA_first = .true.

  call push_sub('h_calc_vhxc')

  if(.not. h%ip_app) then ! No Hartree or xc if independent electrons.
    h%epot = M_ZERO
    h%vhxc = M_ZERO
    if (h%d%cdft) h%ahxc = M_ZERO

    ! now we add the hartree contribution from the density
    allocate(vhartree(m%np))
    vhartree = M_ZERO
    if(h%d%spin_channels == 1) then
      call poisson_solve(m, vhartree, st%rho(:, 1))
    else
      allocate(rho_aux(m%np))                    ! need an auxiliary array to
      rho_aux(:) = st%rho(:, 1)                  ! calculate the total density
      do is = 2, h%d%spin_channels
        rho_aux(:) = rho_aux(:) + st%rho(:, is)
      end do
      call poisson_solve(m, vhartree, rho_aux) ! solve the poisson equation
      deallocate(rho_aux)
    end if

    select case (h%d%ispin)
    case(UNPOLARIZED)
      h%vhxc(:, 1) = h%vhxc(:, 1) + vhartree
    case (SPIN_POLARIZED)
      h%vhxc = h%vhxc + reshape(vhartree, (/m%np, 2/), vhartree)
    case (SPINORS)
      h%vhxc(:, 1:2) = h%vhxc(:, 1:2) + reshape(vhartree, (/m%np, 2/), vhartree)
    end select

    ! We first add 1/2 int n vH, to then subtract int n (vxc + vH)
    ! this yields the correct formula epot = - int n (vxc + vH/2)
    do is = 1, h%d%spin_channels
      h%epot = h%epot + M_HALF*dmf_dotp(m, st%rho(:, is), vhartree)
    end do
    deallocate(vhartree)

    ! now we add the hartree contribution from the current
    if (h%d%cdft) then
      allocate(ahartree(conf%dim, m%np))
      ahartree = M_ZERO
      if(h%d%spin_channels == 1) then
        do idim = 1, conf%dim
          call poisson_solve(m, ahartree(idim, :), st%j(idim, :, 1))
        end do
      else
        allocate(j_aux(conf%dim, m%np))         ! need an auxiliary array to
        j_aux(:,:) = st%j(:, :, 1)              ! calculate the total current
        do is = 2, h%d%spin_channels
          j_aux(:,:) = j_aux(:,:) + st%j(:, :, is)
        end do
        do idim = 1, conf%dim
          call poisson_solve(m, ahartree(idim, :), j_aux(idim, :)) ! solve the poisson equation
        end do
        deallocate(j_aux)
      end if

      select case (h%d%ispin)
      case(UNPOLARIZED)
        h%ahxc(:, :, 1) = h%ahxc(:, :, 1) + ahartree(:, :)
      case (SPIN_POLARIZED)
        h%ahxc = h%ahxc + reshape(ahartree, (/conf%dim, m%np, 2/), ahartree)
      case (SPINORS)
        h%ahxc(:, :, 1:2) = h%ahxc(:, :, 1:2) + reshape(ahartree, (/conf%dim, m%np, 2/), ahartree)
      end select

      ! We first add 1/2 int j.aH, to then subtract int j.(axc + aH)
      ! this yields the correct formula epot = - int j.(axc + aH/2)
      ! WARNING 1: the axc we store is in fact axc/c. So, to get the energy rigth, we have to multiply by c. 
      ! WARNING 2: axc is not yet implemented, so we will just add -1/2 int j.aH
      do is = 1, h%d%spin_channels
        do idim = 1, conf%dim
          h%epot = h%epot - M_HALF*P_c*dmf_dotp(m, st%j(idim, :, is), ahartree(idim, :))
        end do
      end do

      deallocate(ahartree)
    end if

    ! next 3 lines are for an RPA calculation
    !if(RPA_first) then
    !  allocate(RPA_Vhxc(m%np, h%d%nspin))
    !  RPA_Vhxc = h%vhxc

    ! now we calculate the xc terms
    call X(xc_pot)(h%xc, m, st, h%vhxc, h%ex, h%ec, &
         -minval(st%eigenval(st%nst, :)), st%qtot)

    ! The OEP family has to handle specially
    if(iand(h%xc%family, XC_FAMILY_OEP).ne.0) then
      call X(h_xc_oep)(h%xc, m, h, st, h%vhxc, h%ex, h%ec)
    end if

    ! next 5 lines are for an RPA calculation
    !  RPA_Vhxc = h%vhxc - RPA_Vhxc  ! RPA_Vhxc now includes the xc potential
    !  RPA_first = .false.
    !else
    !  h%vhxc = h%vhxc + RPA_Vhxc
    !end if

    select case(h%d%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1)) &
           - dmf_dotp(m, st%rho(:, 2), h%vhxc(:, 2))
    case(SPINORS)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1)) &
           - dmf_dotp(m, st%rho(:, 2), h%vhxc(:, 2))

      allocate(ztmp1(m%np), ztmp2(m%np))
      ztmp1 = st%rho(:, 3) + M_zI*st%rho(:, 4)
      ztmp2 = h%vhxc(:, 3) - M_zI*h%vhxc(:, 4)
      ! WARNING missing real() ???
      h%epot = h%epot - M_TWO*zmf_dotp(m, ztmp1, ztmp2)
      deallocate(ztmp1, ztmp2)
    end select

  end if

  ! this, I think, belongs here
  if(present(calc_eigenval)) call X(hamiltonian_eigenval) (h, m, st)

  call pop_sub()
end subroutine X(h_calc_vhxc)
