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

subroutine X(oep_x) (m, st, is, oep, ex)
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(in)    :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,          intent(out)   :: ex

  integer :: i, j
  FLOAT :: r
  R_TYPE, allocatable :: lx(:)
  FLOAT, allocatable :: rho_ij(:)
#if defined(R_TREAL)
  FLOAT, allocatable :: pot(:)             ! For real
#else
  FLOAT, allocatable :: pot_r(:), pot_i(:) ! For complex
#endif  

  allocate(lx(m%np))

  do i = 1, st%nst
    lx = R_TOTYPE(M_ZERO)
    if(st%occ(i, is) .gt. small) then ! we only need the occupied states
      do j = 1, st%nst
        if(st%occ(j, is) .gt. small) then
          allocate(rho_ij(1:m%np))

#ifdef R_TREAL
          allocate(pot(m%np))
          pot = M_ZERO

          rho_ij(:) = st%dpsi(:, 1, i, is)*st%dpsi(:, 1, j, is)
          call poisson_solve(m, pot, rho_ij)
          deallocate(rho_ij)
          lx(:) = lx(:) - oep%socc*st%occ(j, is)*pot(:)*st%dpsi(:, 1, j, is)
          deallocate(pot)
#else
          allocate(pot_r(m%np), pot_i(m%np))
          pot_r = M_ZERO
          pot_i = M_ZERO

          rho_ij(:) = real(st%zpsi(:, 1, i, is))*real(st%zpsi(:, 1, j, is)) + &
               aimag(st%zpsi(:, 1, i, is))*aimag(st%zpsi(:, 1, j, is))
          call poisson_solve(m, pot_r, rho_ij)
          ! and now the imaginary part
          rho_ij(:) = real(st%zpsi(:, 1, i, is))*aimag(st%zpsi(:, 1, j, is)) - &
               aimag(st%zpsi(:, 1, i, is))*real(st%zpsi(:, 1, j, is))
          call poisson_solve(m, pot_i, rho_ij)
          deallocate(rho_ij)
          lx(:) = lx(:) - st%occ(j, is)*oep%socc* &
               (pot_r(:) + M_zI*pot_i(:))*conjg(st%zpsi(:, 1, j, is))
          deallocate(pot_r, pot_i)
#endif
            
        end if
      end do

      oep%lxc(:, i) = oep%lxc(:, i) + lx(:)

      r = sum(R_REAL(st%X(psi)(:, 1, i, is) * lx(:)))*m%vol_pp
      ex = ex + M_HALF*oep%socc*oep%sfact*st%occ(i, is)*r
    end if
  end do
  
  deallocate(lx)
end subroutine X(oep_x)
