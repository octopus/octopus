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

subroutine R_FUNC(Hpsi) (h, sys, ik, psi, Hpsi)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(out) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: is, idim, np, dim, a, lm
  R_TYPE :: uVpsi
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  
  ! for spin-orbit coupling
  real(r8), allocatable :: Vtot(:), dVtot(:,:)

  !sub_name = 'Hpsi'; call push_sub()

  np = sys%m%np
  dim = sys%st%dim

  call kinetic_energy(sys%m)

  ! Local potential
  do idim = 1, dim
    Hpsi(:, idim) = Hpsi(:, idim) + (h%Vpsl + h%Vhartree)*psi(1:,idim)
  end do

  select case(sys%st%ispin)
  case(1) ! dim = 1
    Hpsi(:, 1) = Hpsi(:, 1) + h%Vxc(:, 1)*psi(1:, 1)
  case(2) ! dim = 1
    if(modulo(ik, 2) == 0) then ! we have a spin down
      Hpsi(:, 1) = Hpsi(:, 1) + h%Vxc(:, 1)*psi(1:, 1)
    else ! spin down
      Hpsi(:, 1) = Hpsi(:, 1) + h%Vxc(:, 2)*psi(1:, 1)
    end if
  case(3) ! dim = 2
    Hpsi(:, 1) = Hpsi(:, 1) + &
         h%Vxc(:, 1)*psi(1:, 1) + h%R_FUNC(Vxc_off)(:)*psi(1:, 2)
    Hpsi(:, 2) = Hpsi(:, 2) + &
         h%Vxc(:, 2)*psi(1:, 2) + R_CONJ(h%R_FUNC(Vxc_off)(:))*psi(1:, 1)
  end select

  ! Non-local part
  call R_FUNC(vnlpsi) (sys, psi, Hpsi)

  ! spin-orbit coupling
#if defined(COMPLEX_WFNS)
  if(h%soc) then
    allocate(Vtot(sys%m%np), dVtot(3, sys%m%np))
    deallocate(Vtot, dVtot)
  end if
#endif

  !call pop_sub()

contains

  subroutine kinetic_energy(m)
    type(mesh_type), intent(IN) :: m

    integer :: idim, i

# if defined(POLYMERS)
    R_TYPE, allocatable :: grad(:,:)
    real(r8) :: k2

    allocate(grad(3,m%np))
    k2 = sum(sys%st%kpoints(:, ik)**2)/2._r8
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (m, psi(:, idim), lapl=Hpsi(:, idim), grad=grad(:,:))

      do i = 1, m%np
        Hpsi(i, idim) = Hpsi(i, idim) &
             + 2._r8*M_zI*sum(sys%st%kpoints(:, ik)*grad(:, i)) &
             - k2*psi(i, idim)
      end do
    end do
    deallocate(grad)
# else
    do idim = 1, dim
      call R_FUNC(mesh_derivatives) (m, psi(:, idim), lapl=Hpsi(:, idim))
    end do
#endif

    Hpsi = - Hpsi/2.0_r8
  end subroutine kinetic_energy

end subroutine R_FUNC(Hpsi)

subroutine R_FUNC(vnlpsi) (sys, psi, Hpsi)
  type(system_type), intent(IN) :: sys
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: is, idim, np, dim, a, lm, i, j, l, m
  R_TYPE :: uVpsi
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  
  np = sys%m%np
  dim = sys%st%dim

  ! Ionic pseudopotential
  do a = 1, sys%natoms
    atm => sys%atom(a)
    spec => atm%spec

    ! do we have a pseudopotential, or a local pot?
    if(.not.spec%local) then
      do idim = 1, dim
        lm = 1
        do l = 0, spec%ps%l_max
           do m = -l, l
              do i = 1, spec%ps%kbc
              do j = 1, spec%ps%kbc
                 uVpsi = sum(atm%uV(:, lm, i)*psi(atm%Jxyz(:), idim))*sys%m%vol_pp * atm%uVu(lm, i, j)
                 Hpsi(atm%Jxyz(:), idim) = Hpsi(atm%Jxyz(:), idim) + uVpsi*atm%uV(:, lm, j)
              end do
              end do
              lm = lm + 1
           end do
        end do
      end do
    end if
  enddo

  !call pop_sub()
end subroutine R_FUNC(vnlpsi)

subroutine R_FUNC(hamiltonian_setup)(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  real(r8), allocatable :: v_aux(:,:)
  R_TYPE, pointer :: v_aux2(:)
  integer :: i

  h%epot = 0._r8 ! The energy coming from the potentials

  if(.not.h%ip_app) then
    call hartree_solve(h%hart, sys%m, h%Vhartree, sys%st%rho)
    do i = 1, sys%st%nspin
      h%epot = h%epot - 0.5_r8*dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vhartree)
    end do
    
    allocate(v_aux(h%np, sys%st%nspin))
    if(h%ispin == 3) then
      allocate(v_aux2(h%np))
    else
      nullify(v_aux2)
    end if

    call R_FUNC(xc_pot)(h%xc, sys%m, sys%st, h%hart, &
         h%Vxc, v_aux, h%ex, h%ec, h%R_FUNC(Vxc_off),  v_aux2)

    h%Vxc = h%Vxc + v_aux
    if(h%ispin == 3) then
      h%R_FUNC(Vxc_off) = h%R_FUNC(Vxc_off) + v_aux2
      deallocate(v_aux2); nullify(v_aux2)
      h%epot = h%epot - sys%m%vol_pp*2._r8* &
           sum(real(R_CONJ(sys%st%R_FUNC(rho_off)(:))*h%R_FUNC(Vxc_off)(:), r8))
    end if

    do i = 1, sys%st%nspin
      h%epot = h%epot - dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vxc(:, i))
    end do
    deallocate(v_aux)
    
  end if

end subroutine R_FUNC(hamiltonian_setup)
