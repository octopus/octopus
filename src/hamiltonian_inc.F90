! calculates the eigenvalues of the real orbitals
subroutine R_FUNC(hamiltonian_eigenval)(h, sys, st_start, st_end)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(inout) :: sys
  integer, intent(in) :: st_start, st_end

  R_TYPE, allocatable :: Hpsi(:,:)
  integer :: ik, ist

  sub_name = 'hamiltonian_eigenval'; call push_sub()
  allocate(Hpsi(sys%m%np, sys%st%dim))

  do ik = 1, sys%st%nik
    do ist = st_start, st_end
      call R_FUNC(Hpsi) (h, sys, ik, sys%st%R_FUNC(psi)(:, :, ist, ik), Hpsi)
      sys%st%eigenval(ist, ik) = R_REAL(R_FUNC(states_dotp)(sys%m, sys%st%dim, sys%st%R_FUNC(psi)(1:, :, ist, ik), Hpsi))
    end do
  end do

  call pop_sub()
  return
end subroutine R_FUNC(hamiltonian_eigenval)

subroutine R_FUNC(Hpsi) (h, sys, ik, psi, Hpsi)
  type(hamiltonian_type), intent(IN) :: h
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: ik
  R_TYPE, intent(IN) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(out) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: idim, np, dim, a, lm
  R_TYPE :: uVpsi
  type(atom_type), pointer :: atm
  type(specie_type), pointer :: spec
  
  !sub_name = 'Hpsi'; call push_sub()

  np = sys%m%np
  dim = sys%st%dim

  ! Kinetic part -\nabla^2/2
  do idim = 1, dim
    call R_FUNC(mesh_derivatives) (sys%m, psi(:, idim), lapl=Hpsi(:, idim))
  end do
  Hpsi = - Hpsi/2.0_r8

  ! Local potential
  ! TODO: one has to implement the k vectors here
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
  case(4) ! dim = 2  ---- WARNING: not tested!!!
    Hpsi(:, 1) = Hpsi(:, 1) + &
         h%Vxc(:, 1)*psi(1:, 1) + h%Vxc(:, 3)*psi(1:, 2)
    Hpsi(:, 2) = Hpsi(:, 2) + &
         h%Vxc(:, 2)*psi(1:, 2) + h%Vxc(:, 4)*psi(1:, 1)
  end select

  ! Ionic pseudopotential
  do a = 1, sys%natoms
    atm => sys%atom(a)
    spec => atm%spec

    ! do we have a pseudopotential, or a local pot?
    if(.not.spec%local) then
      do idim = 1, dim
        do lm = 1, (spec%ps%L_max + 1)**2
          uVpsi = sum(atm%uV(:, lm)*psi(atm%Jxyz(:), idim))*sys%m%vol_pp &
               * atm%uVu(lm)
          Hpsi(atm%Jxyz(:), idim) = Hpsi(atm%Jxyz(:), idim) + uVpsi*atm%uV(:, lm)
        end do
      end do
    end if
  enddo

end subroutine R_FUNC(Hpsi)

subroutine R_FUNC(hamiltonian_setup)(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  real(r8), allocatable :: v_aux(:,:)
  integer :: i

  h%epot = 0._r8 ! The energy coming from the potentials

  if(.not.h%ip_app) then
    call hartree_solve(h%hart, sys%m, h%Vhartree, sys%st%rho)
    do i = 1, sys%st%nspin
      h%epot = h%epot - 0.5_r8*dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vhartree)
    end do
    
    allocate(v_aux(h%np, sys%st%nspin))
    call R_FUNC(xc_pot)(h%xc, sys%m, sys%st, h%hart, h%rho_core, &
         h%Vxc, v_aux, h%ex, h%ec)
    h%Vxc = h%Vxc + v_aux
    do i = 1, sys%st%nspin
      h%epot = h%epot - dmesh_dotp(sys%m, sys%st%rho(:, i), h%Vxc(:, i))
    end do
    deallocate(v_aux)
    
  end if

end subroutine R_FUNC(hamiltonian_setup)
