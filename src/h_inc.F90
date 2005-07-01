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

! calculates the eigenvalues of the real orbitals
subroutine X(hamiltonian_eigenval)(h, m, f_der, st)
  type(hamiltonian_type), intent(in)    :: h
  type(mesh_type),        intent(in)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('hamiltonian_eigenval')
  allocate(Hpsi(m%np, st%d%dim))

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      call X(hpsi) (h, m, f_der, st%X(psi)(:, :, ist, ik), hpsi, ik)
      e = X(states_dotp)(m, st%d%dim, st%X(psi)(1:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
end subroutine X(hamiltonian_eigenval)


! ---------------------------------------------------------
subroutine X(Hpsi) (h, m, f_der, psi, hpsi, ik, t)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out)   :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  FLOAT, optional,        intent(in)    :: t

  call push_sub('Hpsi')

  call X(kinetic) (h, m, f_der, psi, hpsi, ik)
  call X(vlpsi)   (h, m, psi, hpsi, ik)
  if(h%ep%nvnl > 0) call X(vnlpsi)  (h, m, psi, hpsi, ik)

  ! Relativistic corrections...
  select case(h%reltype)
  case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
  case(SPIN_ORBIT)
    call zso (h, m, psi, hpsi)
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
    call X(vlasers)  (h, m, f_der, psi, hpsi, t)
    call X(vborders) (h, psi, hpsi)
  elseif (h%d%cdft) then
    call X(current_extra_terms) (h, m, f_der, psi, hpsi, ik)
  end if

  call pop_sub()
end subroutine X(Hpsi)

! ---------------------------------------------------------
subroutine X(magnus) (h, m, f_der, psi, hpsi, ik, vmagnus)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  type(f_der_type),       intent(inout) :: f_der
  integer,                intent(in)  :: ik
  R_TYPE,                 intent(IN)  :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  FLOAT,                  intent(IN)  :: vmagnus(m%np, h%d%nspin, 2)

  R_TYPE, allocatable :: auxpsi(:, :), aux2psi(:, :)
  integer :: idim
  ! We will assume, for the moment, no spinors.

  call push_sub('magnus')

  allocate(auxpsi(m%np, h%d%dim), aux2psi(m%np, h%d%dim))

  call X(kinetic) (h, m, f_der, psi, hpsi, ik)

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
  call X(kinetic) (h, m, f_der, auxpsi, aux2psi, ik)
  if (h%ep%nvnl > 0) call X(vnlpsi)  (h, m, auxpsi, aux2psi, ik)
  hpsi(:, 1) = hpsi(:, 1) + M_zI*aux2psi(:, 1)

  do idim = 1, h%d%dim
      hpsi(:, idim) = hpsi(:, idim) + h%ep%Vpsl(:)*psi(:,idim)
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
  call X(vborders) (h, psi, hpsi)

  deallocate(auxpsi, aux2psi)
  call pop_sub()
end subroutine X(magnus)


! ---------------------------------------------------------
subroutine X(kinetic) (h, m, f_der, psi, hpsi, ik)
  type(hamiltonian_type), intent(IN)  :: h
  type(mesh_type),        intent(IN)  :: m
  type(f_der_type),       intent(inout) :: f_der
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
    allocate(grad(m%np, conf%dim))
    k2 = sum(h%d%kpoints(:, ik)**2)
    do idim = 1, h%d%dim
      call X(f_laplacian) (f_der, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
      call X(f_gradient)  (f_der, psi(:, idim), grad(:, :))
      do i = 1, m%np
        Hpsi(i, idim) = -M_HALF*(Hpsi(i, idim) &
             + M_TWO*M_zI*sum(h%d%kpoints(:, ik)*grad(i, :)) &
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
      call X(f_laplacian) (f_der, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
    end do
    call lalg_scal(m%np, h%d%dim, R_TOTYPE(-M_HALF), Hpsi)
  end if

  if (h%em_app) then
    call lalg_scal(m%np, h%d%dim, R_TOTYPE(M_ONE/h%m_ratio), Hpsi)
  end if


  call pop_sub()
end subroutine X(kinetic)

! ---------------------------------------------------------
subroutine X(current_extra_terms) (h, m, f_der, psi, hpsi, ik)
  type(hamiltonian_type), intent(in)    :: h
  type(mesh_type),        intent(in)    :: m
  type(f_der_type),       intent(inout) :: f_der
  R_TYPE,                 intent(in)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(out)   :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  integer,                intent(in)    :: ik

  integer :: idim, k
  FLOAT, allocatable :: div(:)
  R_TYPE, allocatable :: grad(:,:)

  call push_sub('current_extra_terms')

  do idim = 1, conf%dim
    select case(h%d%ispin)
    case(UNPOLARIZED)
      hpsi(:, 1) = hpsi(:, 1) + h%ahxc(:, idim, 1)*h%ahxc(:, idim, 1)*psi(:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(:, 1) = hpsi(:, 1) + h%ahxc(:, idim, 1)*h%ahxc(:, idim, 1)*psi(:, 1)
      else
        hpsi(:, 1) = hpsi(:, 1) + h%ahxc(:, idim, 2)*h%ahxc(:, idim, 2)*psi(:, 1)
      end if
    case (SPINORS)
      ! not implemented yet
    end select
  end do

  allocate(div(m%np))
  select case (h%d%ispin)
  case(UNPOLARIZED)
    call df_divergence(f_der, h%ahxc(:, :, 1), div)
    hpsi(:, 1) = hpsi(:, 1) - M_zI * div*psi(1:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then ! we have a spin down
      call df_divergence(f_der, h%ahxc(:, :, 1), div)
    else
      call df_divergence(f_der, h%ahxc(:, :, 2), div)
    end if
    hpsi(:, 1) = hpsi(:, 1) - M_zI * div*psi(1:, 1)
  case(SPINORS)
    ! not implemented yet
  end select
  deallocate(div)

  allocate(grad(m%np, conf%dim))
  select case (h%d%ispin)
  case(UNPOLARIZED)
    call X(f_gradient)(f_der, psi(:, 1), grad)
    do k = 1, m%np
      hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(k,:, 1), grad(k, :))
    end do
  case(SPIN_POLARIZED)
    call X(f_gradient)(f_der, psi(:, 1), grad)
    do k = 1, m%np
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(k, :, 1), grad(k, :))
      else
        hpsi(k, 1) = hpsi(k, 1) - M_zI * dot_product(h%ahxc(k, :, 2), grad(k, :))
      end if
    end do
  case(SPINORS)
    ! not implemented yet
  end select

  if (associated(h%ep%a)) then
    do k = 1, m%np
      hpsi(k, :) = hpsi(k, :) + M_HALF*dot_product(h%ep%a(k, :), h%ep%a(k, :))*psi(k, :)

      select case(h%d%ispin)
      case(UNPOLARIZED)
        hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(h%ep%a(k, :), h%ahxc(k, :, 1))*psi(k, 1)
      case(SPIN_POLARIZED)
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(h%ep%a(k, :), h%ahxc(k, :, 1))*psi(k, 1)
        else
          hpsi(k, 1) = hpsi(k, 1) + M_TWO*dot_product(h%ep%a(k, :), h%ahxc(k, :, 2))*psi(k, 1)
        end if
      case (SPINORS)
        ! not implemented yet
      end select

      select case(h%d%ispin)
      case(UNPOLARIZED)
        hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%ep%a(k, :), grad(k, :))
      case(SPIN_POLARIZED)
        if(modulo(ik+1, 2) == 0) then ! we have a spin down
          hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%ep%a(k, :), h%ahxc(k, :, 1))*psi(k, 1)
        else
          hpsi(k, 1) = hpsi(k, 1) - M_zI*dot_product(h%ep%a(k, :), h%ahxc(k, :, 2))*psi(k, 1)
        end if
      case (SPINORS)
        ! not implemented yet
      end select
 
    end do
  end if
  deallocate(grad)

  call pop_sub()
end subroutine X(current_extra_terms)


! ---------------------------------------------------------
subroutine X(vnlpsi) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(in)    :: h
  type(mesh_type),        intent(in)    :: m
  R_TYPE,                 intent(in)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)
  integer,                intent(in)    :: ik

  integer :: idim, ikbc, jkbc, ivnl
  R_TYPE  :: uVpsi
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
          tmp   = R_TOTYPE(nlop%uv(:, ikbc)*m%vol_pp(nlop%jxyz(:)))
          uvpsi = lalg_dot(nlop%n, tmp, lpsi)*nlop%uvu(ikbc, jkbc)
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


! ---------------------------------------------------------
subroutine X(vlpsi) (h, m, psi, hpsi, ik)
  type(hamiltonian_type), intent(in)    :: h
  type(mesh_type),        intent(in)    :: m
  integer,                intent(in)    :: ik
  R_TYPE,                 intent(in)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)

  integer :: idim
  R_TYPE, allocatable :: lhpsi(:,:)

  call push_sub('vlpsi')

  select case(h%d%ispin)
  case(UNPOLARIZED)
    hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%ep%vpsl(:))*psi(:, 1)
  case(SPIN_POLARIZED)
    if(modulo(ik+1, 2) == 0) then ! we have a spin down
      hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%ep%vpsl(:))*psi(:, 1)
    else
      hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 2) + h%ep%vpsl(:))*psi(:, 1)
    end if
  case(SPINORS)
    hpsi(:, 1) = hpsi(:, 1) + (h%vhxc(:, 1) + h%ep%vpsl(:))*psi(:, 1) + &
                 (h%vhxc(:, 3) + M_zI*h%vhxc(:, 4))*psi(:, 2)
    hpsi(:, 2) = hpsi(:, 2) + (h%vhxc(:, 2) + h%ep%vpsl(:))*psi(:, 2) + &
                 (h%vhxc(:, 3) - M_zI*h%vhxc(:, 4))*psi(:, 1)
  end select

  if (associated(h%ep%e)) then
    do idim = 1, h%d%dim
      hpsi(:, idim) = hpsi(:, idim) + h%ep%v*psi(:, idim)
    end do
  end if

  if (associated(h%ep%b)) then
    allocate(lhpsi(m%np, h%d%dim))
    select case (h%d%ispin)
    case (UNPOLARIZED)
    case (SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        lhpsi(:, 1) = - M_HALF/P_C*sqrt(dot_product(h%ep%b, h%ep%b))*psi(:, 1)

        hpsi(:, 1) = hpsi(:, 1) 
      else
        lhpsi(:, 1) = + M_HALF/P_C*sqrt(dot_product(h%ep%b, h%ep%b))*psi(:, 1)

        hpsi(:, 1) = hpsi(:, 1) 
      end if
    case (SPINORS)
      lhpsi(:, 1) = M_HALF/P_C*( h%ep%b(3)*psi(:, 1) + (h%ep%b(1) - M_zI*h%ep%b(2))*psi(:, 2))
      lhpsi(:, 2) = M_HALF/P_C*(-h%ep%b(3)*psi(:, 2) + (h%ep%b(1) + M_zI*h%ep%b(2))*psi(:, 1))
    end select
    if (h%em_app) then
      call lalg_scal(m%np, h%d%dim, R_TOTYPE(h%g_ratio/h%m_ratio), lhpsi)
    end if
    hpsi = hpsi + lhpsi
    deallocate(lhpsi)
  end if

  call pop_sub()
end subroutine X(vlpsi)


! ---------------------------------------------------------
subroutine X(vlasers) (h, m, f_der, psi, hpsi, t)
  type(hamiltonian_type), intent(IN)    :: h
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  R_TYPE,                 intent(IN)    :: psi(:,:)  !  psi(m%np, h%d%dim)
  R_TYPE,                 intent(inout) :: Hpsi(:,:) !  Hpsi(m%np, h%d%dim)

  FLOAT, intent(in) :: t

  integer :: k, idim
  FLOAT :: x(conf%dim), v, a(conf%dim)
  R_TYPE, allocatable :: grad(:,:)

  call push_sub('vlasers')

  if(h%ep%no_lasers > 0) then
    select case(h%gauge)
    case(1) ! length gauge

      do k = 1, m%np
        call epot_laser_scalar_pot(h%ep, m%x(k,:), t, v)

        hpsi(k,:) = hpsi(k,:) + v * psi(k,:)
      end do

    case(2) ! velocity gauge

      call epot_laser_vector_pot(h%ep, t, a)
      allocate(grad(m%np, conf%dim))
      do idim = 1, h%d%dim
        call X(f_gradient)(f_der, psi(:, idim), grad)
        do k = 1, m%np
          hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(a(:)*grad(k,:)) + &
               sum(a**2)/M_TWO * psi(k, idim)
        end do
      end do
      deallocate(grad)
    end select
  end if

  call pop_sub()
end subroutine X(vlasers)


! ---------------------------------------------------------
subroutine X(vborders) (h, psi, hpsi)
  type(hamiltonian_type), intent(IN)    :: h
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
