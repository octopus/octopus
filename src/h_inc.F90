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
subroutine R_FUNC(hamiltonian_eigenval)(h, sys, st_start, st_end)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(inout) :: sys
  integer, intent(in) :: st_start, st_end

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  sub_name = 'hamiltonian_eigenval'; call push_sub()
  allocate(Hpsi(sys%m%np, sys%st%dim))

  do ik = 1, sys%st%nik
    do ist = st_start, st_end
      call R_FUNC(Hpsi) (h, sys, ik, sys%st%R_FUNC(psi)(:, :, ist, ik), Hpsi)
      e = R_FUNC(states_dotp)(sys%m, sys%st%dim, sys%st%R_FUNC(psi)(1:, :, ist, ik), Hpsi)
      sys%st%eigenval(ist, ik) = R_REAL(e)
    end do
  end do

  deallocate(Hpsi)
  call pop_sub()
  return
end subroutine R_FUNC(hamiltonian_eigenval)

subroutine R_FUNC(Hpsi) (h, sys, ik, psi, Hpsi, t)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(out) :: Hpsi(sys%m%np, sys%st%dim)
  real(r8), intent(in), optional :: t

  sub_name = 'Hpsi'; call push_sub()

  call R_FUNC(kinetic) (sys, psi, Hpsi)
  call R_FUNC(vlpsi)   (h, sys, ik, psi, Hpsi)
  call R_FUNC(vnlpsi)  (sys, psi, Hpsi)
  if(present(t)) then
    call R_FUNC(vlasers)  (h, sys, psi, Hpsi, t)
    call R_FUNC(vborders) (h, sys, psi, Hpsi)
  endif

  ! Relativistic corrections...
  select case(h%reltype)
  case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
  case(SPIN_ORBIT)
    call zso      (h, sys, ik, psi, Hpsi)
#endif
  case default
    message(1) = 'Error: Internal.'
    call write_fatal(1)
  end select

  call pop_sub(); return
end subroutine R_FUNC(Hpsi)

subroutine R_FUNC(kinetic) (sys, psi, Hpsi)
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(out) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: idim, i, np, dim

  sub_name = 'kinetic'; call push_sub()
  np = sys%m%np
  dim = sys%st%dim

# if defined(POLYMERS)
    R_TYPE, allocatable :: grad(:,:)
    real(r8) :: k2

    allocate(grad(3,sys%m%np))
    k2 = sum(sys%st%kpoints(:, ik)**2)/2._r8
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (sys%m, psi(:, idim), lapl=Hpsi(:, idim), grad=grad(:,:))

      do i = 1, sys%m%np
        Hpsi(i, idim) = Hpsi(i, idim) &
             + 2._r8*M_zI*sum(sys%st%kpoints(:, ik)*grad(:, i)) &
             - k2*psi(i, idim)
      end do
    end do
    deallocate(grad)
# else
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (sys%m, psi(:, idim), lapl=Hpsi(:, idim), &
                                     alpha = R_TOTYPE(-1.0_r8/2.0_r8) )
    end do
#endif

  call pop_sub(); return
end subroutine R_FUNC(kinetic)

subroutine R_FUNC(vnlpsi) (sys, psi, Hpsi)
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)

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

    do_dim: do idim = 1, sys%st%dim
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
                       sys%m%vol_pp*atm%R_FUNC(uvu)(add_lm, ikbc, jkbc) 
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

subroutine R_FUNC(vlpsi) (h, sys, ik, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in) :: sys
  integer, intent(in) :: ik
  R_TYPE, intent(in) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: is, idim, np, dim

  sub_name = 'vlpsi'; call push_sub()

  np = sys%m%np
  dim = sys%st%dim

    do idim = 1, dim
      Hpsi(:, idim) = Hpsi(:, idim) + (h%Vpsl + h%Vhartree)*psi(1:,idim)
    end do

    select case(sys%st%ispin)
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

subroutine R_FUNC(vlasers) (h, sys, psi, Hpsi, t)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in) :: sys
  R_TYPE, intent(in) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)
  real(r8), intent(in) :: t

  integer :: is, i, k, idim, np, dim
  real(r8) :: x(3), f(3)
  R_TYPE, allocatable :: grad(:,:)

  sub_name = 'vlasers'; call push_sub()

    if(h%no_lasers > 0) then
      select case(h%gauge)
      case(1) ! length gauge
        call laser_field(h%no_lasers, h%lasers, t, f)

        do k = 1, sys%m%np
          call mesh_xyz(sys%m, k, x)
          hpsi(k,:) = hpsi(k,:) + sum(x*f) * psi(k,:)
        end do

      case(2) ! velocity gauge
        call laser_vector_field(h%no_lasers, h%lasers, t, f)
        allocate(grad(3, sys%m%np))
        do idim = 1, sys%st%dim
           call R_FUNC(mesh_derivatives)(sys%m, psi(:, idim), grad=grad)
           do k = 1, sys%m%np
             hpsi(k, idim) = hpsi(k, idim) - M_zI * sum(f(:)*grad(:, k)) + &
                  sum(f**2)/2._r8 * psi(k, idim)
           end do
        end do
        deallocate(grad)
      end select
    end if

  call pop_sub(); return
end subroutine R_FUNC(vlasers)

subroutine R_FUNC(vborders) (h, sys, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in) :: sys
  R_TYPE, intent(in) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: idim

  sub_name = 'vborders'; call push_sub()

  if(h%ab .eq. 1) then
    do idim = 1, sys%st%dim
       hpsi(:, idim) = hpsi(:, idim) + M_zI*h%ab_pot(:)*psi(1:, idim)
    end do
  end if

  call pop_sub(); return
end subroutine R_FUNC(vborders)

subroutine R_FUNC(hamiltonian_setup)(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: is
  sub_name = 'hamiltonian_setup'; call push_sub()

  if(h%ip_app) return

  call hartree_solve(h%hart, sys%m, h%vhartree, sys%st%rho(:, 1))
  h%epot = - HALF*dmesh_dotp(sys%m, sys%st%rho(:, 1), h%vhartree)

  call R_FUNC(xc_pot)(h%xc, sys%m, sys%st, h%hart, h%vxc, h%ex, h%ec)
  select case(h%ispin)
    case(UNPOLARIZED)
      h%epot = h%epot + dmesh_dotp(sys%m, sys%st%rho(:, 1), h%vxc(:, 1))
    case(SPIN_POLARIZED)
      h%epot = h%epot + dmesh_dotp(sys%m, HALF*(sys%st%rho(:, 1)+sys%st%rho(:, 2)), &
                                          h%vxc(:, 1)) +                            &
                        dmesh_dotp(sys%m, HALF*(sys%st%rho(:, 1)-sys%st%rho(:, 2)), &
                                          h%vxc(:, 2))
    case(SPINORS)
      h%epot = h%epot + dmesh_dotp(sys%m, HALF*(sys%st%rho(:, 1)+sys%st%rho(:, 4)), &
                                          h%vxc(:, 1))
      h%epot = h%epot + dmesh_dotp(sys%m, HALF*(sys%st%rho(:, 1)-sys%st%rho(:, 4)), &
                                          h%vxc(:, 2))
      h%epot = h%epot + TWO*zmesh_dotp(sys%m, HALF*(sys%st%rho(:, 2) + M_zI*sys%st%rho(:, 3)), &
                                                    h%vxc(:, 3) - M_zI*h%vxc(:, 4))

  end select

  call pop_sub(); return
end subroutine R_FUNC(hamiltonian_setup)
