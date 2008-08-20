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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$


! ---------------------------------------------------------
! This routine fills state psi with an atomic orbital -- provided
! by the pseudopotential structure in geo.
!
! Which state will be placed is determined by index "n". Each atom
! provides niwfs pseudo-orbitals (this number is given in geo%atom(ia)%spec%niwfs
! for atom number ia). This number is actually multiplied by two in case
! of spin-unrestricted or spinors calculations.
!
! The pseudo-orbitals are placed in order in the following way (Natoms
! is the total number of atoms).
!
! n = 1 => first orbital of atom 1,
! n = 2 => first orbital of atom 2.
! n = 3 => first orbital of atom 3.
! ....
! n = Natoms => first orbital of atom Natoms
! n = Natoms + 1 = > second orbital of atom 1
! ....
!
! If at some point in this loop an atom pseudo cannot provide the corresponding
! orbital (because the niws orbitals have been exhausted), it moves on to the following
! atom.
!
! In the spinors case, it changes a bit:
!
! n = 1 => first spin-up orbital of atom 1, assigned to the spin-up component of the spinor.
! n = 2 => first spin-down orbital of atom 1, assigned to the spin-down component of the spinor.
! n = 3 => first spin-up orbital of atom 2, assigned to the spin-up component of the spinor.
! ....
! ---------------------------------------------------------
Subroutine X(lcao_initial_wf) (n, m, geo, sb, psi, ispin, ik, kpoints, err)
  integer,                  intent(in)    :: n
  type(mesh_t),             intent(in)    :: m
  type(simul_box_t),        intent(in)    :: sb
  type(geometry_t), target, intent(in)    :: geo
  R_TYPE,                   intent(inout) :: psi(:, :)
  integer,                  intent(in)    :: ispin
  integer,                  intent(in)    :: ik
  FLOAT,                    intent(in)    :: kpoints(:)
  integer,                  intent(out)   :: err

  type(species_t), pointer :: s
  type(periodic_copy_t)   :: pc
  integer :: norbs, ia, i, icell, j, idim, k, wf_dim
  FLOAT :: x(MAX_DIM), pos(MAX_DIM)
  
  
  call push_sub('lcao_inc.Xlcao_initial_wf')

  err = 0

  norbs = 0
  do ia = 1, geo%natoms
    norbs = norbs + geo%atom(ia)%spec%niwfs
  end do

  wf_dim = 1
  if (ispin == SPINORS) then
    wf_dim = 2
    norbs = norbs*2
  end if

  if ((n > norbs) .or. (n < 1)) then
    err = 1
    call pop_sub()
    return
  end if

  psi(1:m%np, 1:wf_dim) = R_TOTYPE(M_ZERO)
  
  select case(ispin)

    case(UNPOLARIZED, SPIN_POLARIZED)
      ! The index "i" goes over all the orbitals supplied by the pseudopotentials, and the
      ! orbitals in placed in psi whenever it matches "n". The index "j" runs over the orbitals
      ! of each atom; whenever j is larger than the number of orbitals that can actually be supplied
      ! by the atom, the atom is skipped.
      i = 1
      j = 0
      do
        j = j + 1
        do ia = 1, geo%natoms
          s => geo%atom(ia)%spec
          if(j > s%niwfs) cycle
          if(n == i) then
            call periodic_copy_init(pc, sb, geo%atom(ia)%x, range = maxval(sb%lsize(1:sb%periodic_dim)))
            do icell = 1, periodic_copy_num(pc)
              pos = periodic_copy_position(pc, sb, icell)
              do k = 1, m%np
                x(1:calc_dim) = m%x(k, 1:calc_dim) - pos(1:calc_dim)
                psi(k, 1) =  psi(k, 1) + exp(M_zI*sum(kpoints(1:calc_dim)*x(1:calc_dim)))* &
                  R_TOTYPE(species_get_iwf(s, j, calc_dim, states_spin_channel(ispin, ik, 1), x(1:calc_dim)))
              end do
            end do
            call periodic_copy_end(pc)
            call X(states_normalize_orbital)(m, wf_dim, psi)
            call pop_sub()
            return
          end if
          i = i + 1
        end do
      end do

    case(SPINORS)

      i = 1
      j = 0
      do
        j = j + 1
        do ia = 1, geo%natoms
          s => geo%atom(ia)%spec
          if(j > s%niwfs) cycle
          do idim = 1, 2
            if(n == i) then
              call periodic_copy_init(pc, sb, geo%atom(ia)%x, range = maxval(sb%lsize(1:sb%periodic_dim)))
              do icell = 1, periodic_copy_num(pc)
                pos = periodic_copy_position(pc, sb, icell)
                do k = 1, m%np
                  x(1:calc_dim) = m%x(k, 1:calc_dim) - geo%atom(ia)%x(1:calc_dim)
                  psi(k, idim) =  &
                       R_TOTYPE(species_get_iwf(s, j, calc_dim, idim, x(1:calc_dim)))
                end do
              end do
              call periodic_copy_end(pc)
              call X(states_normalize_orbital)(m, wf_dim, psi)
              call pop_sub()
              return
            end if
            i = i + 1
          end do
        end do
      end do

  end select

  call pop_sub()
end subroutine X(lcao_initial_wf)

! ---------------------------------------------------------
subroutine X(lcao_init) (this, gr, geo, st, norbs)
  type(lcao_t),         intent(inout) :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(states_t),       intent(in)    :: st
  integer,              intent(in)    :: norbs

  integer :: ik, n, ierr
  integer :: kstart, kend

  call push_sub('lcao_inc.Xlcao_init')

  kstart = this%st%d%kpt%start
  kend = this%st%d%kpt%end

  do ik = kstart, kend
    do n = 1, this%st%nst
      call X(lcao_initial_wf)(n, gr%m, geo, gr%sb, this%st%X(psi)(:, :, n, ik), this%st%d%ispin, ik, st%d%kpoints(:, ik), ierr)
      if(ierr .ne. 0) then
        write(message(1),'(a)') 'Internal error in lcao_wf'
        call write_fatal(1)
      end if
    end do
  end do

  ! Allocation of variables
  ALLOCATE(this%X(s)(norbs, norbs, kstart:kend), norbs**2*this%st%d%kpt%nlocal)

  do ik = kstart, kend
    call states_blockt_mul(gr%m, this%st, this%st%st_start, this%st%st_end, this%st%st_start, this%st%st_end, &
         this%st%X(psi)(:, :, :, ik), this%st%X(psi)(:, :, :, ik), this%X(s)(:, :, ik), symm = .true.)
  end do

  call pop_sub()
end subroutine X(lcao_init)

! ---------------------------------------------------------
subroutine X(lcao_wf) (this, st, gr, h, start)
  type(lcao_t),        intent(in)    :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(hamiltonian_t), intent(in)    :: h
  integer,             intent(in)    :: start
  integer :: dim, nst, ik, n1, n2, idim, norbs, lcao_start
  integer :: sts, ste, ii
  R_TYPE, allocatable :: hpsi(:, :, :)
  FLOAT, allocatable :: ev(:)
  type(batch_t) :: psib, hpsib
  R_TYPE, allocatable :: hamilt(:, :, :)
  integer :: kstart, kend

  call push_sub('lcao_inc.Xlcao_wf')
  
  norbs = this%st%nst
  dim = st%d%dim
  nst = st%nst
  kstart = this%st%d%kpt%start
  kend = this%st%d%kpt%end

  ALLOCATE(hpsi(NP, dim, st%d%block_size), NP*dim*st%d%block_size)
  ALLOCATE(hamilt(norbs, norbs, kstart:kend), norbs**2*this%st%d%kpt%nlocal)

  do ik = kstart, kend
    do sts = 1, norbs, st%d%block_size
      ste = min(norbs, sts + st%d%block_size - 1)

      call batch_init(psib, sts, ste, this%st%X(psi)(:, :, sts:, ik))
      call batch_init(hpsib, sts, ste, hpsi)
      call X(hpsi_batch)(h, gr, psib, hpsib, ik)
      call batch_end(psib)
      call batch_end(hpsib)

      ii = 1
      do n1 = sts, ste
        do n2 = n1, norbs
          hamilt(n1, n2, ik) = X(mf_dotp)(gr%m, dim, hpsi(:, :, ii), this%st%X(psi)(:, :, n2, ik))
          hamilt(n2, n1, ik) = R_CONJ(hamilt(n1, n2, ik))
        end do
        ii = ii + 1
      end do

    end do
  end do

  deallocate(hpsi)

  ALLOCATE(ev(norbs), norbs)

  lcao_start = start
  if(st%parallel_in_states .and. st%st_start > start) lcao_start = st%st_start

  do ik =  kstart, kend
    call lalg_geneigensolve(norbs, hamilt(1:norbs, 1:norbs, ik), this%X(s) (1:norbs, 1:norbs, ik), ev)

    st%eigenval(start:nst, ik) = ev(start:nst)
  
    ! Change of base
    do n1 = lcao_start, st%st_end
      do idim = 1, dim
        st%X(psi)(1:NP, idim, n1, ik) = R_TOTYPE(M_ZERO)
        do n2 = 1, norbs
          call lalg_axpy(NP, hamilt(n2, n1, ik), this%st%X(psi)(:, idim, n2, ik), st%X(psi)(:, idim, n1, ik))
        end do
      end do
    end do

  end do

  deallocate(ev, hamilt)

  call pop_sub()
end subroutine X(lcao_wf)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
