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

subroutine X(oep_x) (m, f_der, st, is, oep, ex)
  type(mesh_type),   intent(IN)    :: m
  type(f_der_type),  intent(in)    :: f_der
  type(states_type), intent(IN)    :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,             intent(out)   :: ex

  integer :: i, j
  FLOAT :: r
  R_TYPE, allocatable :: F_ij(:), lx(:) 

  allocate(F_ij(m%np), lx(m%np))

  do i = 1, st%nst
    if(st%occ(i, is) .le. small) cycle ! we only need the occupied states
    
    lx = R_TOTYPE(M_ZERO)
    do j = i, st%nst
      if(st%occ(j, is) .le. small) cycle
      
      call get_F_ij(i, j, is, F_ij)

      ! lx will the be used to get the energy
      r = M_ONE
      if(i.ne.j) r = M_TWO
      lx(:) = lx(:) - r*oep%socc*st%occ(j, is)*F_ij(:)*st%dpsi(:, 1, j, is)

      oep%X(lxc)(:, i) = oep%X(lxc)(:, i) - &
         oep%socc*st%occ(j, is)*st%dpsi(:, 1, j, is)*F_ij(:)
      if(i.ne.j) then
        oep%X(lxc)(:, j) = oep%X(lxc)(:, j) - &
           oep%socc*st%occ(i, is)*st%dpsi(:, 1, i, is)*R_CONJ(F_ij(:))
      end if
    end do

    r = sum(R_REAL(st%X(psi)(:, 1, i, is)*lx(:))*m%vol_pp(:))
    ex = ex + M_HALF*oep%socc*oep%sfact*st%occ(i, is)*r
  end do

  deallocate(lx, F_ij)

contains
#if defined(R_TREAL)
  ! calculates for real wave-functions
  !    F_ij = \int d^3 r2 \frac{\phi_i,is(r2) \phi_j,is(r2)}{|r-r2|}
  subroutine get_F_ij(i, j, is, F_ij)
    integer, intent(in)  :: i, j, is
    FLOAT,   intent(out) :: F_ij(:)

    FLOAT, allocatable :: rho_ij(:)

    allocate(rho_ij(1:m%np))
      
    rho_ij(:) = st%dpsi(:, 1, i, is)*st%dpsi(:, 1, j, is)

    F_ij = M_ZERO
    call dpoisson_solve(m, f_der, F_ij, rho_ij)

    deallocate(rho_ij)

  end subroutine get_F_ij

#else
  ! WARNING: should call zpoisson_solve
  ! calculates for complex wave-functions
  !    F_ij = \int d^3 r2 \frac{\phi_i,is^*(r2) \phi_j,is(r2)}{|r-r2|}
  subroutine get_F_ij(i, j, is, F_ij)
    integer, intent(in)  :: i, j, is
    CMPLX,   intent(out) :: F_ij(:)

    FLOAT, allocatable :: rho_ij(:), pot(:)

    allocate(rho_ij(m%np), pot(m%np))
    
    ! first we solve for the real part
    rho_ij(:) = real(st%zpsi(:, 1, i, is))*real(st%zpsi(:, 1, j, is)) + &
       aimag(st%zpsi(:, 1, i, is))*aimag(st%zpsi(:, 1, j, is))

    pot(:) = M_ZERO
    call dpoisson_solve(m, f_der, pot, rho_ij)

    ! this is the real part of F_ij
    F_ij(:) = pot(:)

    ! and now the imaginary part
    rho_ij(:) = real(st%zpsi(:, 1, i, is))*aimag(st%zpsi(:, 1, j, is)) - &
       aimag(st%zpsi(:, 1, i, is))*real(st%zpsi(:, 1, j, is))

    pot(:) = M_ZERO
    call dpoisson_solve(m, f_der, pot, rho_ij)
    deallocate(rho_ij)

    F_ij(:) = F_ij(:) + M_zI*pot(:)

    deallocate(rho_ij, pot)
  
  end subroutine get_F_ij
#endif


end subroutine X(oep_x)
