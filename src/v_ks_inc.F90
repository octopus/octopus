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

subroutine X(h_calc_vhxc)(ks, h, m, f_der, st, calc_eigenval)
  type(v_ks_type),        intent(inout) :: ks
  type(hamiltonian_type), intent(inout) :: h
  type(mesh_type),        intent(IN)    :: m
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  logical,      optional, intent(in)    :: calc_eigenval

  integer :: is, idim
  FLOAT, allocatable :: rho_aux(:), vhartree(:), j_aux(:,:), ahartree(:,:), rho(:, :)
  CMPLX, allocatable :: ztmp1(:), ztmp2(:)

  ! next 2 lines are for an RPA calculation
  !FLOAT, save, pointer ::  RPA_Vhxc(:, :)
  !logical, save :: RPA_first = .true.

  ! No Hartree or xc if independent electrons.
  if(ks%ip_app) then
    if(present(calc_eigenval)) call X(hamiltonian_eigenval) (h, m, f_der, st)
    return
  end if

  call push_sub('h_calc_vhxc')

  h%epot = M_ZERO
  h%vhxc = M_ZERO
  if (h%d%cdft) h%ahxc = M_ZERO
  
  ! now we add the hartree contribution from the density
  allocate(vhartree(m%np))
  vhartree = M_ZERO
  if(h%d%spin_channels == 1) then
    call dpoisson_solve(m, f_der, vhartree, st%rho(:, 1))
  else
    allocate(rho_aux(m%np))                    ! need an auxiliary array to
    rho_aux(:) = st%rho(:, 1)                  ! calculate the total density
    do is = 2, h%d%spin_channels
      rho_aux(:) = rho_aux(:) + st%rho(:, is)
    end do
    call dpoisson_solve(m, f_der, vhartree, rho_aux) ! solve the poisson equation
    deallocate(rho_aux)
  end if
  if (h%em_app) then
    call lalg_scal(m%np, M_ONE/h%e_ratio, vhartree)
  end if
  
  h%vhxc(:, 1) = h%vhxc(:, 1) + vhartree(:)
  if(h%d%ispin > UNPOLARIZED) h%vhxc(:, 2) = h%vhxc(:, 2) + vhartree
  
  ! We first add 1/2 int n vH, to then subtract int n (vxc + vH)
  ! this yields the correct formula epot = - int n (vxc + vH/2)
  do is = 1, h%d%spin_channels
    h%epot = h%epot + M_HALF*dmf_dotp(m, st%rho(:, is), vhartree)
  end do
  deallocate(vhartree)
  
  ! now we add the hartree contribution from the current
  if (h%d%cdft) then
    allocate(ahartree(m%np, conf%dim))
    ahartree = M_ZERO
    if(h%d%spin_channels == 1) then
      do idim = 1, conf%dim
        call dpoisson_solve(m, f_der, ahartree(:, idim), st%j(:, idim, 1))
      end do
    else
      allocate(j_aux(conf%dim, m%np))         ! need an auxiliary array to
      j_aux(:,:) = st%j(:,:, 1)               ! calculate the total current
      do is = 2, h%d%spin_channels
        j_aux(:,:) = j_aux(:,:) + st%j(:, :, is)
      end do
      do idim = 1, conf%dim
        call dpoisson_solve(m, f_der, ahartree(:,idim), j_aux(:, idim)) ! solve the poisson equation
      end do
      deallocate(j_aux)
    end if

    h%ahxc(:,:, 1) = h%ahxc(:,:, 1) + ahartree(:,:)
    if(h%d%ispin > UNPOLARIZED) h%ahxc(:,:, 2) = h%ahxc(:,:, 2) + ahartree(:,:)
    
    ! We first add 1/2 int j.aH, to then subtract int j.(axc + aH)
    ! this yields the correct formula epot = - int j.(axc + aH/2)
    ! WARNING 1: the axc we store is in fact axc/c. So, to get the energy rigth, we have to multiply by c. 
    do is = 1, h%d%spin_channels
      do idim = 1, conf%dim
        h%epot = h%epot + M_HALF*P_c*dmf_dotp(m, st%j(idim, :, is), ahartree(:, idim))
      end do
    end do
    
    deallocate(ahartree)
  end if

  ! next 3 lines are for an RPA calculation
  !if(RPA_first) then
  !  allocate(RPA_Vhxc(m%np, h%d%nspin))
  !  RPA_Vhxc = h%vhxc
  
  ! now we calculate the xc terms
  h%ex = M_ZERO
  h%ec = M_ZERO
  h%exc_j = M_ZERO
  allocate(rho(m%np, st%d%nspin))
  if(associated(st%rho_core)) then
     do is = 1, st%d%spin_channels
        rho(:, is) = st%rho(:, is) + st%rho_core(:)/st%d%spin_channels
     enddo
  else
     rho = st%rho
  endif
  if (h%d%cdft) then
    call xc_get_vxc_and_axc(ks%xc, m, f_der, rho, st%j, st%d%ispin, h%vhxc, h%ahxc, &
       h%ex, h%ec, h%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
  else
    call xc_get_vxc(ks%xc, m, f_der, rho, st%d%ispin, h%vhxc, h%ex, h%ec, &
       -minval(st%eigenval(st%nst, :)), st%qtot)
  end if
  deallocate(rho)

  ! The OEP family has to handle specially
  call X(xc_oep_calc)(ks%oep, ks%xc, m, f_der, h, st, h%vhxc, h%ex, h%ec)
  
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
  
  ! this, I think, belongs here
  if(present(calc_eigenval)) call X(hamiltonian_eigenval) (h, m, f_der, st)

  call pop_sub()
end subroutine X(h_calc_vhxc)
