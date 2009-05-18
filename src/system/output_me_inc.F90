!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any %later version.
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
!!
!! $Id: unocc.F90 5320 2009-04-23 15:29:57Z marques $


! ---------------------------------------------------------
! Prints out the mult	ipole matrix elements between KS states.
! It prints the states to the file opened in iunit.
! It prints the (l,m) multipole moment, for
! the Kohn-Sham states in the irreducible subspace ik.
! ---------------------------------------------------------
subroutine X(output_me_ks_multipoles)(fname, st, gr, l, m, ik)
  character(len=*), intent(in) :: fname
  type(states_t),   intent(in) :: st
  type(grid_t),     intent(in) :: gr
  integer,          intent(in) :: l, m, ik
  
  integer :: ii, jj, i, iunit
  FLOAT, allocatable :: multipole(:, :)
  FLOAT :: r, x(MAX_DIM), ylm
  R_TYPE :: multip_element
  
  call push_sub('output_me_inc.output_me_ks_multipoles')
  
  iunit = io_open(file = fname, action = 'write')

  write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
  write(iunit, fmt = '(a,i2,a,i2)') '# l =', l, '; m =', m
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  if(l>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_out%length%abbrev)//']^', l
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_out%length%abbrev)//']'
  end if
  
  SAFE_ALLOCATE(multipole(1:gr%mesh%np_part, 1:st%d%dim))
  multipole = M_ZERO
  do ii = 1, st%d%dim
    do i = 1, gr%mesh%np
      call mesh_r(gr%mesh, i, r, x = x)
      call loct_ylm(1, x(1), x(2), x(3), l, m, ylm)
      multipole(i, ii) = r**l * ylm
    end do
  end do
  
  do ii = 1, st%nst
    do jj = 1, st%nst
      multip_element = X(mf_dotp) (gr%mesh, st%d%dim, &
           R_CONJ(st%X(psi)(:, :, ii, 1)), &
           st%X(psi)(:, :, jj, 1) * multipole(:, :))
      multip_element = multip_element / units_out%length%factor**l

      write(iunit, fmt='(f20.12)', advance = 'no') R_REAL(multip_element)
#if defined(R_TCOMPLEX)
      write(iunit, fmt='(a1,f20.12)', advance = 'no') ',', R_AIMAG(multip_element)
#endif
      write(iunit, fmt='(a)', advance = 'no') '   '
    end do
    write(iunit, '(a)') ''
  end do
  
  SAFE_DEALLOCATE_A(multipole)
  call io_close(iunit)

  call pop_sub()
end subroutine X(output_me_ks_multipoles)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
