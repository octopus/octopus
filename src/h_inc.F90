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
subroutine R_FUNC(hamiltonian_eigenval)(h, st, sys)
  type(hamiltonian_type), intent(IN) :: h
  type(states_type), intent(inout) :: st
  type(system_type), intent(in) :: sys

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  sub_name = 'hamiltonian_eigenval'; call push_sub()
  allocate(Hpsi(sys%m%np, st%dim))

  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      call R_FUNC(Hpsi) (h, sys%m, st, sys, ik, st%R_FUNC(psi)(:, :, ist, ik), Hpsi)
      e = R_FUNC(states_dotp)(sys%m, st%dim, st%R_FUNC(psi)(1:, :, ist, ik), Hpsi)
      st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
  return
end subroutine R_FUNC(hamiltonian_eigenval)

subroutine R_FUNC(Hpsi) (h, m, st, sys, ik, psi, Hpsi, t)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(in) :: sys ! this is necessary due to the nl part of the PP
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(0:m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)
  real(r8), intent(in), optional :: t

  sub_name = 'Hpsi'; call push_sub()

  call R_FUNC(kinetic) (m, st, psi, Hpsi)
  call R_FUNC(vlpsi)   (h, m, st, ik, psi, Hpsi)
#if defined(THREE_D)
  call R_FUNC(vnlpsi)  (m, st, sys, psi, Hpsi)
#endif
  if(present(t)) then
    call R_FUNC(vlasers)  (h, m, st, psi, Hpsi, t)
    call R_FUNC(vborders) (h, m, st, psi, Hpsi)
  endif

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

  call pop_sub(); return
end subroutine R_FUNC(Hpsi)

subroutine R_FUNC(kinetic) (m, st, psi, Hpsi)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  R_TYPE, intent(IN) :: psi(0:m%np, st%dim)
  R_TYPE, intent(out) :: Hpsi(m%np, st%dim)

  integer :: idim, i, np, dim
  R_TYPE :: d

  sub_name = 'kinetic'; call push_sub()
  np = m%np
  dim = st%dim

# if defined(POLYMERS)
    R_TYPE, allocatable :: grad(:,:)
    real(r8) :: k2

    allocate(grad(3, m%np))
    k2 = sum(st%kpoints(:, ik)**2)/2._r8
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (m, psi(:, idim), lapl=Hpsi(:, idim), grad=grad(:,:))

      do i = 1, m%np
        Hpsi(i, idim) = Hpsi(i, idim) &
             + 2._r8*M_zI*sum(st%kpoints(:, ik)*grad(:, i)) &
             - k2*psi(i, idim)
      end do
    end do
    deallocate(grad)
# else
    d = R_TOTYPE(-1.0_r8/2.0_r8)
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (m, psi(:, idim), lapl=Hpsi(:, idim), alpha=d )
    end do
#endif

  call pop_sub(); return
end subroutine R_FUNC(kinetic)

#if defined(THREE_D)
subroutine R_FUNC(vnlpsi) (m, st, sys, psi, Hpsi)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(0:m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: is, idim, ia, ikbc, jkbc, l, lm, add_lm
  R_TYPE :: uVpsi
  R_TYPE, allocatable :: lpsi(:), lHpsi(:)
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  R_TYPE, external :: R_DOT

  sub_name = 'vnlpsi'; call push_sub()

  ! Ionic pseudopotential
  do_atm: do ia = 1, sys%natoms
    atm => sys%atom(ia)
    spec => atm%spec
    ! do we have a pseudopotential, or a local pot?
    if(spec%local) cycle do_atm

    do_dim: do idim = 1, st%dim
      allocate(lpsi(atm%mps), lHpsi(atm%mps))
      lpsi(:) = psi(atm%jxyz(:), idim)
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
               uvpsi = R_DOT(atm%mps, atm%R_FUNC(uv)(:, add_lm, ikbc), 1, lpsi(:), 1) * &
                       m%vol_pp*atm%R_FUNC(uvu)(add_lm, ikbc, jkbc) 
               call R_FUNC(axpy) (atm%mps, uvpsi, atm%R_FUNC(uv)(:, add_lm, jkbc), 1, lHpsi(:), 1)
            end do
          end do

          add_lm = add_lm + 1
        end do do_m
      end do do_l
      Hpsi(atm%jxyz(:), idim) = Hpsi(atm%jxyz(:), idim) + lHpsi(:)
      deallocate(lpsi, lHpsi)
    end do do_dim
  end do do_atm

  call pop_sub(); return
end subroutine R_FUNC(vnlpsi)
#endif

subroutine R_FUNC(vlpsi) (h, m, st, ik, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in) :: ik
  R_TYPE, intent(in) :: psi(0:m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: is, idim, np, dim

  sub_name = 'vlpsi'; call push_sub()

  np = m%np
  dim = st%dim

    do idim = 1, dim
      Hpsi(:, idim) = Hpsi(:, idim) + (h%Vpsl + h%Vhartree)*psi(1:,idim)
    end do

    select case(st%ispin)
    case(UNPOLARIZED)
      hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1)
    case(SPIN_POLARIZED)
      if(modulo(ik+1, 2) == 0) then ! we have a spin down
        hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1)
      else
        hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 2)*psi(1:, 1)
      end if
    case(SPINORS)
      hpsi(:, 1) = hpsi(:, 1) + h%vxc(:, 1)*psi(1:, 1) + &
                   (h%vxc(:, 3) - M_zI*h%vxc(:, 4))*psi(1:, 2)
      hpsi(:, 2) = hpsi(:, 2) + h%vxc(:, 2)*psi(1:, 2) + &
                   (h%vxc(:, 3) + M_zI*h%vxc(:, 4))*psi(1:, 1)
    end select

  call pop_sub(); return
end subroutine R_FUNC(vlpsi)

subroutine R_FUNC(vlasers) (h, m, st, psi, Hpsi, t)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(0:m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)
  real(r8), intent(in) :: t

  integer :: is, i, k, idim, np, dim
  real(r8) :: x(3), f(3)
  R_TYPE, allocatable :: grad(:,:)

  sub_name = 'vlasers'; call push_sub()

    if(h%no_lasers > 0) then
      select case(h%gauge)
      case(1) ! length gauge
        call laser_field(h%no_lasers, h%lasers, t, f)

        do k = 1, m%np
          call mesh_xyz(m, k, x)
          hpsi(k,:) = hpsi(k,:) + sum(x*f) * psi(k,:)
        end do

      case(2) ! velocity gauge
        call laser_vector_field(h%no_lasers, h%lasers, t, f)
        allocate(grad(3, m%np))
        do idim = 1, st%dim
           call R_FUNC(mesh_derivatives)(m, psi(:, idim), grad=grad)
           do k = 1, m%np
             hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
                  sum(f**2)/2._r8 * psi(k, idim)
           end do
        end do
        deallocate(grad)
      end select
    end if

  call pop_sub(); return
end subroutine R_FUNC(vlasers)

subroutine R_FUNC(vborders) (h, m, st, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  R_TYPE, intent(in) :: psi(0:m%np, st%dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, st%dim)

  integer :: idim

  sub_name = 'vborders'; call push_sub()

  if(h%ab .eq. 1) then
    do idim = 1, st%dim
       hpsi(:, idim) = hpsi(:, idim) + M_zI*h%ab_pot(:)*psi(1:, idim)
    end do
  end if

  call pop_sub(); return
end subroutine R_FUNC(vborders)

subroutine R_FUNC(hamiltonian_setup)(h, m, st, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type), intent(in) :: m
  type(states_type), intent(inout) :: st
  type(system_type), intent(in), optional :: sys

  integer :: is

  if(h%ip_app) return

  sub_name = 'hamiltonian_setup'; call push_sub()

  call hartree_solve(h%hart, m, h%vhartree, st%rho(:, 1))
  h%epot = - M_HALF*dmesh_dotp(m, st%rho(:, 1), h%vhartree)

  call R_FUNC(xc_pot)(h%xc, m, st, h%hart, h%vxc, h%ex, h%ec)
  select case(h%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot + dmesh_dotp(m, st%rho(:, 1), h%vxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot + dmesh_dotp(m, M_HALF*(st%rho(:, 1)+st%rho(:, 2)), &
                                          h%vxc(:, 1)) +                            &
                        dmesh_dotp(m, M_HALF*(st%rho(:, 1)-st%rho(:, 2)), &
                                          h%vxc(:, 2))
    case(SPINORS)
      h%epot = h%epot + dmesh_dotp(m, M_HALF*(st%rho(:, 1)+st%rho(:, 4)), &
                                          h%vxc(:, 1))
      h%epot = h%epot + dmesh_dotp(m, M_HALF*(st%rho(:, 1)-st%rho(:, 4)), &
                                          h%vxc(:, 2))
      h%epot = h%epot + M_TWO*zmesh_dotp(m, M_HALF*(st%rho(:, 2) + M_zI*st%rho(:, 3)), &
                                                    h%vxc(:, 3) - M_zI*h%vxc(:, 4))

  end select

  ! this, I think, belongs here
  if(present(sys)) call R_FUNC(hamiltonian_eigenval) (h, st, sys)

  call pop_sub(); return
end subroutine R_FUNC(hamiltonian_setup)
