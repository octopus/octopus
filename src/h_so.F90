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

  integer :: is, ia, i, j
  real(r8) :: x, y, z, r, xa(3)
  complex(r8), allocatable :: gradpsi(:, :, :), lpsi(:, :, :), tpsi(:, :)

  sub_name = 'zso'; call push_sub()

  allocate(gradpsi(3,   sys%m%np, 2),  &
           lpsi   (3,   sys%m%np, 2),  &
           tpsi   (   0:sys%m%np, 2))


  call zmesh_derivatives(sys%m, psi(:, 1), grad = gradpsi(:, :, 1), alpha = -M_zI)
  call zmesh_derivatives(sys%m, psi(:, 2), grad = gradpsi(:, :, 2), alpha = -M_zI)
  do_atm: do ia = 1, sys%natoms
    do i = 1, sys%m%np
       call mesh_r(sys%m, i, r, sys%atom(ia)%x, xa)
       lpsi(1, i, :) = xa(2)*gradpsi(3, i, :) - xa(3)*gradpsi(2, i, :)
       lpsi(2, i, :) = xa(3)*gradpsi(1, i, :) - xa(1)*gradpsi(3, i, :)
       lpsi(3, i, :) = xa(1)*gradpsi(2, i, :) - xa(2)*gradpsi(1, i, :)
    end do
    tpsi(0, 1:2) = M_Z0
    call zcopy(sys%m%np,        lpsi(1, 1:sys%m%np, 2), 1, tpsi(1:, 1), 1)
    call zcopy(sys%m%np,        lpsi(1, 1:sys%m%np, 1), 1, tpsi(1:, 2), 1)
    call zaxpy(sys%m%np, -M_zI, lpsi(2, 1:sys%m%np, 2), 1, tpsi(1:, 1), 1)
    call zaxpy(sys%m%np,  M_zI, lpsi(2, 1:sys%m%np, 1), 1, tpsi(1:, 2), 1)
    call zaxpy(sys%m%np,  M_z1, lpsi(3, 1:sys%m%np, 1), 1, tpsi(1:, 1), 1)
    call zaxpy(sys%m%np, -M_z1, lpsi(3, 1:sys%m%np, 2), 1, tpsi(1:, 2), 1)
    call zscal(2*sys%m%np+2, R_TOTYPE(1.0/2.0), tpsi, 1)
    call znlso(sys%atom(ia), tpsi(:, 1), Hpsi(:, 1)) 
    call znlso(sys%atom(ia), tpsi(:, 2), Hpsi(:, 2)) 
  end do do_atm

  deallocate(gradpsi, lpsi, tpsi)

  call pop_sub(); return
  contains

  subroutine znlso (atm, psi, Hpsi)
    type(atom_type), intent(IN) :: atm
    complex(r8), intent(IN) :: psi(0:sys%m%np)
    complex(r8), intent(inout) :: Hpsi(sys%m%np)

    integer :: ikbc, jkbc, l, lm, add_lm
    complex(r8) :: uVpsi
    type(specie_type), pointer :: spec

    sub_name = 'zsopsi'; call push_sub()

    spec => atm%spec
    ! do we have a pseudopotential, or a local pot?
    if(spec%local) return

    add_lm = 1
    l_loop: do l = 0, spec%ps%l_max
       if (l == spec%ps%L_loc) then
          add_lm = add_lm + (2*l + 1)
          cycle l_loop
       end if
       m_loop: do lm = -l, l
          do ikbc = 1, spec%ps%kbc
             do jkbc = 1, spec%ps%kbc
                uVpsi = sum(atm%uV(:, add_lm, ikbc)*psi(atm%Jxyz(:)))*sys%m%vol_pp * &
                     atm%uvuso(add_lm, ikbc, jkbc)
                Hpsi(atm%Jxyz(:)) = Hpsi(atm%Jxyz(:)) + uVpsi*atm%uV(:, add_lm, jkbc)
             end do
          end do
          add_lm = add_lm + 1
       end do m_loop
    end do l_loop

    call pop_sub(); return
  end subroutine znlso

end subroutine zso
