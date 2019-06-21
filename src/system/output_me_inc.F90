!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! Now we have three "ks_multipoles" routines, for the 1, 2, and 3D
! cases. Eventually they can probably be merged back into one, once
! we write the 2D case in more general form, using "cylindrical multipoles",
! analogous to the spherical ones.

! ---------------------------------------------------------
!> Prints out the multipole matrix elements between KS states.
!!
!! It prints the states to the file opened in iunit.
!! It prints the (ll,mm) multipole moment, for
!! the Kohn-Sham states in the irreducible subspace ik.
!!
! ---------------------------------------------------------
subroutine X(output_me_ks_multipoles)(fname, st, gr, ll, mm, ik)
  character(len=*), intent(in) :: fname
  type(states_t),   intent(in) :: st
  type(grid_t),     intent(in) :: gr
  integer,          intent(in) :: ll, mm, ik
  
  integer :: ist, jst, ip, iunit
  FLOAT, allocatable :: multipole(:)
  FLOAT :: rr, xx(MAX_DIM), ylm
  R_TYPE :: multip_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  iunit = io_open(file = fname, action = 'write')

  write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
  write(iunit, fmt = '(a,i2,a,i2)') '# l =', ll, '; m =', mm
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  if(ll>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_abbrev(units_out%length))//']^', ll
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  end if
  
  SAFE_ALLOCATE(multipole(1:gr%mesh%np))

  multipole = M_ZERO
  do ip = 1, gr%mesh%np
    call mesh_r(gr%mesh, ip, rr, coords = xx)
    call loct_ylm(1, xx(1), xx(2), xx(3), ll, mm, ylm)
    multipole(ip) = rr**ll * ylm
  end do
  
  ASSERT(.not. st%parallel_in_states)
  ! how to do this properly? states_matrix
  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:gr%mesh%np, 1) = psij(1:gr%mesh%np, 1)*multipole(1:gr%mesh%np)

      multip_element = X(mf_dotp)(gr%mesh, st%d%dim, psii, psij)

      multip_element = units_from_atomic(units_out%length**ll, multip_element)

      write(iunit, fmt='(f20.12)', advance = 'no') R_REAL(multip_element)
#if defined(R_TCOMPLEX)
      write(iunit, fmt='(a1,f20.12)', advance = 'no') ',', R_AIMAG(multip_element)
#endif
      write(iunit, fmt='(a)', advance = 'no') '   '
    end do
    write(iunit, '(a)') ''
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)
  SAFE_DEALLOCATE_A(multipole)
  call io_close(iunit)

  POP_SUB(X(output_me_ks_multipoles))
end subroutine X(output_me_ks_multipoles)



! ---------------------------------------------------------
!> Prints out the dipole matrix elements (X or Y) between single
!! orbitals, in the 1d case.
!!
!! It prints the states to the file opened in iunit.
!! It prints the moment, for single orbital states
!! irreducible subspace ik. It only prints the first order moments
!! X or Y. Eventually it should print the circular multipoles of
!! arbitrary order, similar to the 3D case.
!!
!! The argument dir should be 1 (X) or 2 (Y).
! ---------------------------------------------------------
subroutine X(output_me_ks_multipoles2d)(fname, st, gr, dir, ik)
  character(len=*), intent(in) :: fname
  type(states_t),   intent(in) :: st
  type(grid_t),     intent(in) :: gr
  integer,          intent(in) :: dir, ik
  
  integer :: ist, jst, ip, iunit
  FLOAT, allocatable :: dipole(:)
  FLOAT :: rr, xx(MAX_DIM)
  R_TYPE :: dipole_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles2d))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  iunit = io_open(file = fname, action = 'write')

  select case(dir)
  case(1)
    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | X | Phi_j>' 
  case(2)
    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | Y | Phi_j>' 
  end select
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  
  SAFE_ALLOCATE(dipole(1:gr%mesh%np))

  dipole = M_ZERO
  do ip = 1, gr%mesh%np
    call mesh_r(gr%mesh, ip, rr, coords = xx)
    dipole(ip) = xx(dir)
  end do
  
  ASSERT(.not. st%parallel_in_states)
  ! how to do this properly? states_matrix
  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:gr%mesh%np, 1) = psij(1:gr%mesh%np, 1) * dipole(1:gr%mesh%np)

      dipole_element = X(mf_dotp)(gr%mesh, st%d%dim, psii, psij)

      dipole_element = units_from_atomic(units_out%length, dipole_element)

      write(iunit, fmt='(f20.12)', advance = 'no') R_REAL(dipole_element)
#if defined(R_TCOMPLEX)
      write(iunit, fmt='(a1,f20.12)', advance = 'no') ',', R_AIMAG(dipole_element)
#endif
      write(iunit, fmt='(a)', advance = 'no') '   '
    end do
    write(iunit, '(a)') ''
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)
  SAFE_DEALLOCATE_A(dipole)
  call io_close(iunit)

  POP_SUB(X(output_me_ks_multipoles2d))
end subroutine X(output_me_ks_multipoles2d)


! ---------------------------------------------------------
!> Prints out the multipole matrix elements (X**l) between single
!! orbitals, in the 1d case.
!!
!! It prints the states to the file opened in iunit.
!! It prints the moment of ll-th order, for single orbital states
!! irreducible subspace ik.
! ---------------------------------------------------------
subroutine X(output_me_ks_multipoles1d)(fname, st, gr, ll, ik)
  character(len=*), intent(in) :: fname
  type(states_t),   intent(in) :: st
  type(grid_t),     intent(in) :: gr
  integer,          intent(in) :: ll, ik

  integer :: iunit, ip, ist, jst
  FLOAT, allocatable :: dipole(:)
  FLOAT :: rr, xx(MAX_DIM)
  R_TYPE :: dipole_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles1d))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  iunit = io_open(file = fname, action = 'write')

  write(iunit, fmt = '(a)') '# Dipole matrix elements file: <Phi_i | X**l | Phi_j>' 
  write(iunit, fmt = '(a,i2)') '# l =', ll
  write(iunit, fmt = '(a,i4)') '# ik =', ik
  if(ll>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_abbrev(units_out%length))//']^', ll
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  end if
  
  SAFE_ALLOCATE(dipole(1:gr%mesh%np))


  dipole = M_ZERO
  do ip = 1, gr%mesh%np
    call mesh_r(gr%mesh, ip, rr, coords = xx)
    dipole(ip) = xx(1)**ll
  end do
  
  ASSERT(.not. st%parallel_in_states)
  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:gr%mesh%np, 1) = psij(1:gr%mesh%np, 1) * dipole (1:gr%mesh%np)

      dipole_element = X(mf_dotp)(gr%mesh, st%d%dim, psii, psij)

      dipole_element = units_from_atomic(units_out%length**ll, dipole_element)

      write(iunit, fmt='(f20.12)', advance = 'no') R_REAL(dipole_element)
#if defined(R_TCOMPLEX)
      write(iunit, fmt='(a1,f20.12)', advance = 'no') ',', R_AIMAG(dipole_element)
#endif
      write(iunit, fmt='(a)', advance = 'no') '   '
    end do
    write(iunit, '(a)') ''
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)
  SAFE_DEALLOCATE_A(dipole)
  call io_close(iunit)

  POP_SUB(X(output_me_ks_multipoles1d))
end subroutine X(output_me_ks_multipoles1d)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
