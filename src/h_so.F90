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

subroutine zso (h, m, psi, hpsi, natoms, atom, dim, ik)
  type(hamiltonian_type), intent(in) :: h
  type(mesh_type), intent(in) :: m
  integer, intent(in) :: dim, natoms, ik
  type(atom_type), intent(in) :: atom(natoms)
  R_TYPE, intent(in) :: psi(m%np, dim)
  R_TYPE, intent(inout) :: Hpsi(m%np, dim)

  integer :: is, ia, i, j,  mps,add_lm, ikbc,jkbc, idim, l, lm 
  CMPLX, allocatable :: tpsi(:, :), tHpsi(:, :)
  type(specie_type), pointer :: spec
  CMPLX :: uvpsi

  call push_sub('zso')

  atm: do ia = 1, natoms
     spec => atom(ia)%spec
     allocate(tpsi(atom(ia)%mps, 2), tHpsi(atom(ia)%mps, 2))
     tpsi(:, 1) = psi(atom(ia)%jxyz(:), 1)
     tpsi(:, 2) = psi(atom(ia)%jxyz(:), 2)
     tHpsi = M_z0
        add_lm = 1
        do l = 0, spec%ps%l_max
           do lm = -l, l
              do ikbc = 1, spec%ps%kbc
                 do jkbc = 1, spec%ps%kbc
                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 1), 1, &
                               tpsi(1, 1), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 2), 1)
                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 1), 1, &
                               tpsi(1, 2), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 1), 1)
 
                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 2), 1, &
                               tpsi(1, 1), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, M_zI*uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 2), 1)
                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 2), 1, &
                               tpsi(1, 2), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, -M_zI*uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 1), 1)

                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 3), 1, &
                               tpsi(1, 1), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 1), 1)
                     uvpsi = la_dot(atom(ia)%mps, atom(ia)%so_luv(1, add_lm, ikbc, 3), 1, &
                               tpsi(1, 2), 1)*m%vol_pp*atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call la_axpy(atom(ia)%mps, -uvpsi/2, atom(ia)%so_uv(1, add_lm, jkbc), 1, &
                            tHpsi(1, 2), 1)
                  enddo
              enddo
              add_lm = add_lm + 1
           enddo
        enddo
        Hpsi(atom(ia)%jxyz(:), 1) = Hpsi(atom(ia)%jxyz(:), 1) + tHpsi(:, 1)
        Hpsi(atom(ia)%jxyz(:), 2) = Hpsi(atom(ia)%jxyz(:), 2) + tHpsi(:, 2)
     deallocate(tpsi, tHpsi)
  end do atm

  call pop_sub()
end subroutine zso
