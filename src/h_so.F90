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

subroutine zso (h, m, psi, hpsi, dim, ik)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  integer, intent(in) :: dim, ik
  R_TYPE, intent(in) :: psi(:, :)
  R_TYPE, intent(inout) :: Hpsi(:, :)

  type(nonlocal_op), pointer :: nlop
  integer :: is, ia, i, j,  mps,add_lm, ikbc,jkbc, idim, l, lm, ivnl
  CMPLX, allocatable :: tpsi(:, :), tHpsi(:, :)
  type(specie_type), pointer :: spec
  CMPLX :: uvpsi

  call push_sub('zso')
  do ivnl = 1, h%ep%nvnl
     nlop => h%ep%vnl(ivnl)
     mps = nlop%n
     allocate(tpsi(mps, 2), thpsi(mps, 2))
     tpsi(:, 1) = psi(nlop%jxyz(:), 1)
     tpsi(:, 2) = psi(nlop%jxyz(:), 2)
     thpsi = M_z0
     do ikbc = 1, nlop%c
        do jkbc = 1, nlop%c
          ! WARNING: Not every efficient. Has to be changed, and checked.
          stop 'Does not work due to vol_pp'
          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 1), &
          !   tpsi(:, 1))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 2))
          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 1), &
          !   tpsi(:, 2))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 1))
 
          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 2), &
          !   tpsi(:, 1))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, M_zI*uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 2))
          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 2), &
          !   tpsi(:, 2))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, -M_zI*uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 1))

          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 3), &
          !   tpsi(:, 1))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 1))
          !uvpsi = lalg_dot(mps, nlop%so_luv(:, ikbc, 3), &
          !   tpsi(:, 2))*m%vol_pp*nlop%so_uvu(ikbc, jkbc)
          call lalg_axpy(mps, -uvpsi/2, nlop%so_uv(:, jkbc), tHpsi(:, 2))
        end do
     end do
     hpsi(nlop%jxyz(:), 1) = hpsi(nlop%jxyz(:), 1) + thpsi(:, 1)
     hpsi(nlop%jxyz(:), 2) = hpsi(nlop%jxyz(:), 2) + thpsi(:, 2)
     deallocate(tpsi, thpsi)
   end do

  call pop_sub()
end subroutine zso
