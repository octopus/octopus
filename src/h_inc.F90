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
subroutine X(hamiltonian_eigenval)(h, st, sys)
  type(hamiltonian_type), intent(IN) :: h
  type(states_type), intent(inout) :: st
  type(system_type), intent(in) :: sys

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('hamiltonian_eigenval')
  allocate(Hpsi(sys%m%np, st%dim))

  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      call X(Hpsi) (h, sys%m, st, sys, ik, st%X(psi)(:, :, ist, ik), Hpsi)
      e = X(states_dotp)(sys%m, st%dim, st%X(psi)(1:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
end subroutine X(hamiltonian_eigenval)

subroutine X(Hpsi) (h, m, st, sys, ik, psi, Hpsi, t)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(in) :: sys ! this is necessary due to the nl part of the PP
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)
  real(r8), intent(in), optional :: t

  call push_sub('Hpsi')

  call X(kinetic) (h, ik, m, st, psi, Hpsi)
  call X(vlpsi)   (h, m, st, ik, psi, Hpsi)
  if(sys%nlpp) call X(vnlpsi)  (ik, m, st, sys, psi, Hpsi)

  ! Relativistic corrections...
  select case(h%reltype)
  case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
  case(SPIN_ORBIT)
    call zso (h, m, sys%natoms, sys%atom, st%dim, ik, psi, Hpsi)
#endif
  case default
    message(1) = 'Error: Internal.'
    call write_fatal(1)
  end select

  if(present(t)) then
    call X(vlasers)  (h, m, st, psi, Hpsi, t)
    call X(vborders) (h, m, st, psi, Hpsi)
  endif

  call pop_sub()
end subroutine X(Hpsi)

subroutine X(kinetic) (h, ik, m, st, psi, Hpsi)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)
  integer :: ik

  R_TYPE, allocatable :: grad(:,:)
  real(r8) :: k2
  integer :: idim, i, np, dim

  call push_sub('kinetic')
  np = m%np
  dim = st%dim

  if(conf%periodic_dim>0) then
#if defined(COMPLEX_WFNS)
    allocate(grad(3, m%np))
    k2 = sum(st%kpoints(:, ik)**2)
    do idim = 1, dim
      call X(f_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
      call X(f_gradient)  (m, psi(:, idim), grad(:, :))
      do i = 1, m%np
        Hpsi(i, idim) = -M_HALF*(Hpsi(i, idim) &
             + M_TWO*M_zI*sum(st%kpoints(:, ik)*grad(:, i)) &
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
    do idim = 1, dim
      call X(f_laplacian) (m, psi(:, idim), Hpsi(:, idim), cutoff_ = M_TWO*h%cutoff)
    end do
    call X(scal)(m%np*dim, R_TOTYPE(-M_HALF), Hpsi(1, 1), 1)
  end if


  call pop_sub()
end subroutine X(kinetic)

subroutine X(vnlpsi) (ik, m, st, sys, psi, Hpsi)
  integer, intent(in) :: ik
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: i, is, idim, ia, ikbc, jkbc, k, l, lm, add_lm
  real(r8) :: x(conf%dim)
  R_TYPE :: phase
  R_TYPE :: uVpsi
  R_TYPE, allocatable :: lpsi(:), lHpsi(:)
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
#ifdef R_TCOMPLEX
  complex(r8), allocatable :: conj(:)
#endif

  call push_sub('vnlpsi')

  ! Ionic pseudopotential
  do_atm: do ia = 1, sys%natoms
    atm => sys%atom(ia)
    spec => atm%spec
    ! do we have a pseudopotential, or a local pot?
    if(spec%local) cycle do_atm

    ! calculate the phase factors exp(ik*x) 
    ! they are necessary to get the correct periodicity 
    ! for the nonlocal part of the pseudopotential
    if (conf%periodic_dim/=0) then
      allocate(atm%phases(atm%mps,st%nik))
      do i=1,atm%mps
        call mesh_xyz(m, atm%jxyz(i), x)
        do k=1,st%nik
          atm%phases(i,k)=exp(M_zI*sum(st%kpoints(:,k)*x(:)))
        end do
      end do
    end if

#   ifdef R_TCOMPLEX
      allocate(conj(atm%mps))
#   endif

    do_dim: do idim = 1, st%dim
      allocate(lpsi(atm%mps), lHpsi(atm%mps))
      if (conf%periodic_dim==0) then
        lpsi(:) = psi(atm%jxyz(:), idim)
      else
        lpsi(:) = atm%phases(:,ik)*psi(atm%jxyz(:), idim)
      end if
      lHpsi(:) = M_z0
      add_lm = 1
      do_l: do l = 0, spec%ps%l_max
        if (l == spec%ps%L_loc) then
          add_lm = add_lm + (2*l + 1)
          cycle do_l
        end if
        
        do_m: do lm = -l, l
          do ikbc = 1, spec%ps%kbc
            do jkbc = 1, spec%ps%kbc
              uvpsi = R_DOT(atm%mps, atm%X(uv)(1, add_lm, ikbc), 1, lpsi(1), 1) * &
                   m%vol_pp*atm%X(uvu)(add_lm, ikbc, jkbc) 
              if (conf%periodic_dim==0) then
                call X(axpy) (atm%mps, uvpsi, atm%X(uv)(1, add_lm, jkbc), 1, lHpsi(1), 1)
              else
#               ifdef R_TCOMPLEX
                  conj = R_CONJ(atm%phases(:,ik))*atm%X(uv)(:, add_lm, jkbc)
                  call X(axpy) (atm%mps, uvpsi, conj(1), 1, lHpsi(1), 1)
#               else
                  message(1) = "Real wavefunction for ground state not yet implemented for polymers:"
                  message(2) = "Reconfigure with --enable-complex, and remake"
                  call write_fatal(2)
#               endif
              end if
            end do
          end do
          
          add_lm = add_lm + 1
        end do do_m
      end do do_l
      Hpsi(atm%jxyz(:), idim) = Hpsi(atm%jxyz(:), idim) + lHpsi(:)
      deallocate(lpsi, lHpsi)
    end do do_dim

#   ifdef R_TCOMPLEX
      deallocate(conj)
#   endif
  end do do_atm

  call pop_sub()
end subroutine X(vnlpsi)

subroutine X(vlpsi) (h, m, st, ik, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in) :: ik
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: is, idim, np, dim

  call push_sub('vlpsi')

  np = m%np
  dim = st%dim

    do idim = 1, dim
      Hpsi(:, idim) = Hpsi(:, idim) + h%Vpsl*psi(1:,idim)
    end do

    select case(st%ispin)
    case(UNPOLARIZED)
      hpsi(:, 1) = hpsi(:, 1) + h%vhxc(:, 1)*psi(1:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(:, 1) = hpsi(:, 1) + h%vhxc(:, 1)*psi(1:, 1)
      else
        hpsi(:, 1) = hpsi(:, 1) + h%vhxc(:, 2)*psi(1:, 1)
      end if
    case(SPINORS)
      hpsi(:, 1) = hpsi(:, 1) + h%vhxc(:, 1)*psi(1:, 1) + &
                   (h%vhxc(:, 3) - M_zI*h%vhxc(:, 4))*psi(1:, 2)
      hpsi(:, 2) = hpsi(:, 2) + h%vhxc(:, 2)*psi(1:, 2) + &
                   (h%vhxc(:, 3) + M_zI*h%vhxc(:, 4))*psi(1:, 1)
    end select

  call pop_sub()
end subroutine X(vlpsi)

subroutine X(vlasers) (h, m, st, psi, Hpsi, t)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)
  real(r8), intent(in) :: t

  integer :: is, i, k, idim, np, dim
  real(r8) :: x(3), f(3)
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
      do idim = 1, st%dim
        call X(f_gradient)(m, psi(:, idim), grad)
        do k = 1, m%np
          hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
               sum(f**2)/2._r8 * psi(k, idim)
        end do
      end do
      deallocate(grad)
    end select
  end if
  
  call pop_sub()
end subroutine X(vlasers)

subroutine X(vborders) (h, m, st, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: idim

  call push_sub('vborders')

  if(h%ab .eq. IMAGINARY_ABSORBING) then
    do idim = 1, st%dim
       hpsi(:, idim) = hpsi(:, idim) + M_zI*h%ab_pot(:)*psi(1:, idim)
    end do
  end if

  call pop_sub()
end subroutine X(vborders)

subroutine X(h_calc_vhxc)(h, m, st, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type),        intent(in)    :: m
  type(states_type),      intent(inout) :: st
  type(system_type),      intent(in), optional :: sys

  integer :: is
  real(r8), allocatable :: rho_aux(:), vhartree(:)

  call push_sub('h_calc_vhxc')

  if(.not. h%ip_app) then ! No Hartree or xc if independent electrons.
    h%epot = M_ZERO

    ! first we calculate the xc terms
    call X(xc_pot)(h%xc, m, st, h%vhxc, h%ex, h%ec, &
         -minval(st%eigenval(st%nst, :)), st%qtot)

    select case(h%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1)) &
           - dmf_dotp(m, st%rho(:, 2), h%vhxc(:, 2))
    case(SPINORS)
      h%epot = h%epot - dmf_dotp(m, st%rho(:, 1), h%vhxc(:, 1)) &
           - dmf_dotp(m, st%rho(:, 2), h%vhxc(:, 2))
      h%epot = h%epot - M_TWO*zmf_dotp(m, st%rho(:, 3) + M_zI*st%rho(:, 4), &
           h%vhxc(:, 3) - M_zI* h%vhxc(:, 4))
    end select

    ! now we add the hartree contribution
    allocate(vhartree(m%np))
    vhartree = M_ZERO
    if(st%spin_channels == 1) then
      call poisson_solve(m, vhartree, st%rho(:, 1))
    else
      allocate(rho_aux(m%np))                    ! need an auxiliary array to
      rho_aux(:) = st%rho(:, 1)                  ! calculate the total density
      do is = 2, st%spin_channels
        rho_aux(:) = rho_aux(:) + st%rho(:, is)
      end do
      call poisson_solve(m, vhartree, rho_aux) ! solve the poisson equation
      deallocate(rho_aux)
    end if

    select case (h%ispin)
    case(UNPOLARIZED)
       h%vhxc(:, 1) = h%vhxc(:, 1) + vhartree
    case (SPIN_POLARIZED)
       h%vhxc = h%vhxc + reshape(vhartree, (/m%np, 2/), vhartree)
    case (SPINORS)
       h%vhxc(:, 1:2) = h%vhxc(:, 1:2) + reshape(vhartree, (/m%np, 2/), vhartree)
    end select

    do is = 1, st%spin_channels
      h%epot = h%epot - M_HALF*dmf_dotp(m, st%rho(:, is), vhartree)
    end do

    deallocate(vhartree)
  end if

  ! this, I think, belongs here
  if(present(sys)) call X(hamiltonian_eigenval) (h, st, sys)

  call pop_sub()
end subroutine X(h_calc_vhxc)
