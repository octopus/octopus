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

subroutine zso (h, sys, ik, psi, Hpsi)
  type(hamiltonian_type), intent(in) :: h
  type(system_type), intent(in) :: sys
  integer, intent(in) :: ik
  R_TYPE, intent(in) :: psi(0:sys%m%np, sys%st%dim)
  R_TYPE, intent(inout) :: Hpsi(sys%m%np, sys%st%dim)

  integer :: is, ia, i, j,  mps,add_lm, ikbc,jkbc, idim, l, lm 
  complex(r8), allocatable :: tpsi(:, :), tHpsi(:, :)
  type(specie_type), pointer :: spec
  complex(r8) :: uvpsi
  R_TYPE, external :: R_DOT

  sub_name = 'zso'; call push_sub()

  atm: do ia = 1, sys%natoms
     spec => sys%atom(ia)%spec
     allocate(tpsi(sys%atom(ia)%mps, 2), tHpsi(sys%atom(ia)%mps, 2))
     tpsi(:, 1) = psi(sys%atom(ia)%jxyz(:), 1)
     tpsi(:, 2) = psi(sys%atom(ia)%jxyz(:), 2)
     tHpsi = M_z0
        add_lm = 1
        do l = 0, spec%ps%l_max
           do lm = -l, l
              do ikbc = 1, spec%ps%kbc
                 do jkbc = 1, spec%ps%kbc
                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 1), 1, &
                               tpsi(:, 1), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,       uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 2), 1)
                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 1), 1, &
                               tpsi(:, 2), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,       uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 1), 1)
 
                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 2), 1, &
                               tpsi(:, 1), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,   M_zI*uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 2), 1)
                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 2), 1, &
                               tpsi(:, 2), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,  -M_zI*uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 1), 1)

                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 3), 1, &
                               tpsi(:, 1), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,       uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 1), 1)
                     uvpsi = R_DOT(sys%atom(ia)%mps, sys%atom(ia)%so_luv(:, add_lm, ikbc, 3), 1, &
                               tpsi(:, 2), 1)*sys%m%vol_pp*sys%atom(ia)%so_uvu(add_lm, ikbc, jkbc)
                     call zaxpy(sys%atom(ia)%mps,      -uvpsi/2, sys%atom(ia)%so_uv(:, add_lm, jkbc), 1, &
                            tHpsi(:, 2), 1)
                  enddo
              enddo
              add_lm = add_lm + 1
           enddo
        enddo
        Hpsi(sys%atom(ia)%jxyz(:), 1) = Hpsi(sys%atom(ia)%jxyz(:), 1) + tHpsi(:, 1)
        Hpsi(sys%atom(ia)%jxyz(:), 2) = Hpsi(sys%atom(ia)%jxyz(:), 2) + tHpsi(:, 2)
     deallocate(tpsi, tHpsi)
  end do atm

  call pop_sub(); return
end subroutine zso
