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
subroutine X(output_me_ks_multipoles)(fname, namespace, st, gr, ll, mm, ik)
  character(len=*),    intent(in) :: fname
  type(namespace_t),   intent(in) :: namespace
  type(states_elec_t), intent(in) :: st
  type(grid_t),        intent(in) :: gr
  integer,             intent(in) :: ll, mm, ik
  
  integer :: ist, jst, ip, iunit
  FLOAT, allocatable :: multipole(:)
  FLOAT :: rr, xx(1:gr%sb%dim), ylm
  R_TYPE :: multip_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  iunit = io_open(fname, namespace, action = 'write')

  write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
  write(iunit, fmt = '(a,i2,a,i2)') '# l =', ll, '; m =', mm
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  if(ll>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_abbrev(units_out%length))//']^', ll
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  end if
  
  SAFE_ALLOCATE(multipole(1:gr%mesh%np))

  do ip = 1, gr%mesh%np
    call mesh_r(gr%mesh, ip, rr, coords = xx)
    call loct_ylm(1, xx(1), xx(2), xx(3), ll, mm, ylm)
    multipole(ip) = rr**ll * ylm
  end do
  
  ASSERT(.not. st%parallel_in_states)
  ! how to do this properly? states_elec_matrix
  do ist = 1, st%nst

    call states_elec_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_elec_get_state(st, gr%mesh, jst, 1, psij)

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
subroutine X(output_me_ks_multipoles2d)(fname, namespace, st, mesh, dir, ik)
  character(len=*),    intent(in) :: fname
  type(namespace_t),   intent(in) :: namespace
  type(states_elec_t), intent(in) :: st
  type(mesh_t),        intent(in) :: mesh
  integer,             intent(in) :: dir, ik
    
  integer :: ist, jst, ip, iunit
  FLOAT, allocatable :: dipole(:)
  R_TYPE :: dipole_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles2d))

  SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:mesh%np, 1:st%d%dim))
  
  iunit = io_open(fname, namespace, action = 'write')

  select case(dir)
  case(1)
    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | X | Phi_j>' 
  case(2)
    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | Y | Phi_j>' 
  end select
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  
  SAFE_ALLOCATE(dipole(1:mesh%np))

  do ip = 1, mesh%np
    dipole(ip) = mesh%x(ip, dir)
  end do
  
  ASSERT(.not. st%parallel_in_states)
  ! how to do this properly? states_elec_matrix
  do ist = 1, st%nst

    call states_elec_get_state(st, mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_elec_get_state(st, mesh, jst, 1, psij)

      psij(1:mesh%np, 1) = psij(1:mesh%np, 1) * dipole(1:mesh%np)

      dipole_element = X(mf_dotp)(mesh, st%d%dim, psii, psij)

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
subroutine X(output_me_ks_multipoles1d)(fname, namespace, st, mesh, ll, ik)
  character(len=*),    intent(in) :: fname
  type(namespace_t),   intent(in) :: namespace
  type(states_elec_t), intent(in) :: st
  type(mesh_t),        intent(in) :: mesh
  integer,             intent(in) :: ll, ik

  integer :: iunit, ip, ist, jst
  FLOAT, allocatable :: dipole(:)
  R_TYPE :: dipole_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles1d))

  SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:mesh%np, 1:st%d%dim))
  
  iunit = io_open(fname, namespace, action = 'write')

  write(iunit, fmt = '(a)') '# Dipole matrix elements file: <Phi_i | X**l | Phi_j>' 
  write(iunit, fmt = '(a,i2)') '# l =', ll
  write(iunit, fmt = '(a,i4)') '# ik =', ik
  if(ll>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_abbrev(units_out%length))//']^', ll
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  end if
  
  SAFE_ALLOCATE(dipole(1:mesh%np))

  do ip = 1, mesh%np
    dipole(ip) = mesh%x(ip, 1)**ll
  end do
  
  ASSERT(.not. st%parallel_in_states)
  do ist = 1, st%nst

    call states_elec_get_state(st, mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_elec_get_state(st, mesh, jst, 1, psij)

      psij(1:mesh%np, 1) = psij(1:mesh%np, 1) * dipole (1:mesh%np)

      dipole_element = X(mf_dotp)(mesh, st%d%dim, psii, psij)

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

! ---------------------------------------------------------
!> Prints out the dipole matrix elements between KS states.
!!
!! It prints the states to the file opened in iunit.
!!
! ---------------------------------------------------------
subroutine X(output_me_dipole)(this, fname, namespace, space, st, gr, hm, ions, ik)
  type(output_me_t),   intent(in) :: this
  character(len=*),    intent(in) :: fname
  type(namespace_t),   intent(in) :: namespace
  type(space_t),       intent(in) :: space
  type(states_elec_t), intent(in) :: st
  type(grid_t),        intent(in) :: gr
  type(hamiltonian_elec_t), intent(in) :: hm
  type(ions_t),        intent(in) :: ions
  integer,             intent(in) :: ik
  
  integer :: ist, jst, ip, iunit, idir, idim, ispin
  R_TYPE :: dip_element
  R_TYPE, allocatable :: psii(:, :), psij(:, :), gpsii(:,:,:)

  PUSH_SUB(X(output_me_dipole))

  ASSERT(.not. st%parallel_in_states)
  if(st%d%ispin == SPINORS) then
    call messages_not_implemented("Dipole matrix elements with spinors", namespace=namespace)
  end if

  ispin = st%d%get_spin_index(ik)

  SAFE_ALLOCATE(psii(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(gpsii(1:gr%mesh%np, 1:space%dim, 1:st%d%dim))
  
  do idir = 1, space%dim

    iunit = io_open(trim(fname)//index2axis(idir), namespace, action = 'write')

    write(iunit, '(a)') '# Dipole matrix elements file: |<Phi_i | r | Phi_j>|' 
    write(iunit, '(a,i4)')      '# ik =', ik
    write(iunit, '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  
    do ist = this%st_start, this%st_end

      call states_elec_get_state(st, gr%mesh, ist, ik, psii)

      if (.not. space%is_periodic()) then

        do idim = 1, st%d%dim  
          do ip = 1, gr%mesh%np
            gpsii(ip, idir, idim) = psii(ip, idim)*gr%mesh%x(ip, idir)
          end do
        end do

      else

        do idim = 1, st%d%dim
           call boundaries_set(gr%der%boundaries, psii(:, idim))
        end do

        !We need the phase here as the routines for the nonlocal contributions assume that the wavefunctions have a phase.
#ifdef R_TCOMPLEX
        if (allocated(hm%hm_base%phase)) then
          call states_elec_set_phase(st%d, psii, hm%hm_base%phase(1:gr%mesh%np_part, ik), gr%mesh%np_part, .false.)
        end if
#endif

        do idim = 1, st%d%dim
          call X(derivatives_grad)(gr%der, psii(:, idim), gpsii(:, :, idim), set_bc = .false.)
        end do

        !A nonlocal contribution from the MGGA potential must be included
        !This must be done first, as this is like a position-dependent mass 
        if (family_is_mgga_with_exc(hm%xc)) then
          do idim = 1, st%d%dim
            !$omp parallel do
            do ip = 1, gr%mesh%np
              gpsii(ip, idir, idim) = (M_ONE + M_TWO*hm%vtau(ip, ispin))*gpsii(ip, idir, idim)
            end do
            !$omp end parallel do
          end do
        end if

        !A nonlocal contribution from the pseudopotential must be included
        call X(projector_commute_r_allatoms_alldir)(hm%ep%proj, ions, gr%mesh, st%d%dim, &
                    gr%der%boundaries, ik, psii, gpsii) 
        
        !A nonlocal contribution from the scissor must be included
        if(hm%scissor%apply) then
          call scissor_commute_r(hm%scissor, gr%mesh, ik, psii, gpsii)
        end if

        if(hm%lda_u_level /= DFT_U_NONE) then
          call X(lda_u_commute_r)(hm%lda_u, gr%mesh, st%d, namespace, ik, psii, gpsii, allocated(hm%hm_base%phase))
        end if
      end if

      do jst = this%st_start, this%st_end

        call states_elec_get_state(st, gr%mesh, jst, ik, psij)
#ifdef R_TCOMPLEX
        if (allocated(hm%hm_base%phase)) then
          call states_elec_set_phase(st%d, psij, hm%hm_base%phase(1:gr%mesh%np, ik), gr%mesh%np, .false.)
        end if
#endif

        dip_element = X(mf_dotp)(gr%mesh, st%d%dim, gpsii(:, idir, :), psij)
        if (space%is_periodic()) then
          if(abs(st%eigenval(ist, ik) - st%eigenval(jst, ik)) > CNST(1e-5)) then  
            dip_element = -dip_element/((st%eigenval(ist, ik) - st%eigenval(jst, ik)))
          else
            dip_element = R_TOTYPE(M_ZERO)
          end if
        end if
        dip_element = units_from_atomic(units_out%length, dip_element)

        write(iunit, fmt='(f20.12)', advance = 'no') abs(dip_element)
        write(iunit, fmt='(a)', advance = 'no') '   '
      end do
      write(iunit, '(a)') ''
    end do
    call io_close(iunit)

  end do

  SAFE_DEALLOCATE_A(gpsii)
  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)

  POP_SUB(X(output_me_dipole))
end subroutine X(output_me_dipole)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
