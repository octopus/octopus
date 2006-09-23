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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$


! ---------------------------------------------------------
subroutine X(lcao_initial_wf) (n, m, geo, psi, ispin, ik, err)
  integer,                  intent(in)  :: n
  type(mesh_t),             intent(in)  :: m
  type(geometry_t), target, intent(in)  :: geo
  R_TYPE,                   intent(out) :: psi(:, :)
  integer,                  intent(in)  :: ispin
  integer,                  intent(in)  :: ik
  integer,                  intent(out) :: err

  integer :: norbs, ia, i, j, idim, k, wf_dim
  type(specie_t), pointer :: s
  FLOAT :: x(MAX_DIM), r

  call push_sub('lcao_inc.Xlcao_initial_wf')

  err = 0
  psi = R_TOTYPE(M_ZERO)

  norbs = 0
  do ia = 1, geo%natoms
    norbs = norbs + geo%atom(ia)%spec%niwfs
  end do

  wf_dim = 1
  if (ispin == SPINORS) then
    wf_dim = 2
    norbs = norbs * 2
  end if

  if ((n > norbs) .or. (n < 1)) then
    err = 1
    call pop_sub()
    return
  end if

  idim = 1
  i = 1; j = 0
  do
    j = j + 1
    do ia = 1, geo%natoms
      s => geo%atom(ia)%spec
      do idim = 1, wf_dim
        if(j > s%niwfs) cycle
        if(n == i) then
          do k = 1, m%np
            x(1:calc_dim) = m%x(k, 1:calc_dim) - geo%atom(ia)%x(1:calc_dim)
            psi(k, idim) =  R_TOTYPE(specie_get_iwf(s, j, calc_dim, states_spin_channel(ispin, ik, idim), x(1:calc_dim)))
          end do
          r = X(states_nrm2)(m, wf_dim, psi)
          psi = psi/r
          call pop_sub()
          return
        end if
        i = i + 1
      end do
    end do
  end do

  call pop_sub()
end subroutine X(lcao_initial_wf)


! ---------------------------------------------------------
subroutine X(lcao_init) (lcao_data, gr, geo, h, norbs)
  type(lcao_t), target, intent(inout) :: lcao_data
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  integer,              intent(in)    :: norbs

  type(states_t), pointer :: st
  integer :: ik, n1, n2, n, ierr
  R_TYPE, allocatable :: hpsi(:,:)

  call push_sub('lcao_inc.Xlcao_init')

  st => lcao_data%st

  do ik = 1, st%d%nik
    do n = 1, st%nst
      call X(lcao_initial_wf) (n, gr%m, geo, st%X(psi)(:, :, n, ik), st%d%ispin, ik, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Internal error in lcao_wf.'
        call write_fatal(1)
      end if
    end do
  end do

  ! Allocation of variables
  ALLOCATE(lcao_data%X(hamilt) (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
  ALLOCATE(lcao_data%X(s)      (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
  ALLOCATE(lcao_data%X(k)      (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
  ALLOCATE(lcao_data%X(v)      (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)

  ! Overlap and kinetic+so matrices.
  ALLOCATE(hpsi(NP_PART, st%d%dim), NP_PART*st%d%dim)
  do ik = 1, st%d%nik
    do n1 = 1, st%nst
      call X(kinetic) (h, gr, st%X(psi)(:, :, n1, ik), hpsi(:, :), ik)
      ! Relativistic corrections...
      select case(h%reltype)
      case(NOREL)
#ifdef R_TCOMPLEX
      case(SPIN_ORBIT)
        call zso (h, gr, st%zpsi(:, :, n1, ik), hpsi(:, :), ik)
#endif
      end select

      do n2 = n1, st%nst
        lcao_data%X(k)(n1, n2, ik) = X(states_dotp)(gr%m, st%d%dim, hpsi, st%X(psi)(:, : ,n2, ik))
        lcao_data%X(k)(n2, n1, ik) = lcao_data%X(k)(n1, n2, ik)
        lcao_data%X(s)(n1, n2, ik) = X(states_dotp)(gr%m, st%d%dim, st%X(psi)(:, :, n1, ik), &
             st%X(psi)(:, : ,n2, ik))
        lcao_data%X(s)(n2, n1, ik) = lcao_data%X(s)(n1, n2, ik)
      end do

    end do
  end do
  deallocate(hpsi)

  call pop_sub()
end subroutine X(lcao_init)


! ---------------------------------------------------------
subroutine X(lcao_wf) (lcao_data, st, m, h, start)
  type(lcao_t),        intent(inout) :: lcao_data
  type(states_t),      intent(inout) :: st
  type(mesh_t),        intent(in)    :: m
  type(hamiltonian_t), intent(in)    :: h
  integer,             intent(in)    :: start

  integer :: np, dim, nst, ik, n1, n2, idim, norbs
  R_TYPE, allocatable :: hpsi(:,:)
  FLOAT, allocatable :: ev(:)

  call push_sub('lcao_inc.Xlcao_wf')

  norbs = lcao_data%st%nst
  np = m%np
  dim = st%d%dim
  nst = st%nst

  ! Hamiltonian and overlap matrices.
  ALLOCATE(hpsi(m%np_part, dim), m%np_part*dim)
  do ik = 1, st%d%nik
    do n1 = 1, lcao_data%st%nst
      hpsi = M_ZERO
      call X(vlpsi) (h, m, lcao_data%st%X(psi)(:,:, n1, ik), hpsi(:,:), ik)
      if (h%ep%nvnl > 0) call X(vnlpsi) (h, m, lcao_data%st%X(psi)(:,:, n1, ik), hpsi(:,:), ik)
      do n2 = n1, lcao_data%st%nst
        lcao_data%X(v) (n1, n2, ik) = X(states_dotp)(m, dim, hpsi, lcao_data%st%X(psi)(1:, : ,n2, ik))
        lcao_data%X(hamilt) (n1, n2, ik) = lcao_data%X(k) (n1, n2, ik) + lcao_data%X(v) (n1 , n2, ik)
        lcao_data%X(hamilt) (n2, n1, ik) = R_CONJ(lcao_data%X(hamilt) (n1, n2, ik))
      end do
    end do
  end do

  do ik = 1, st%d%nik
    ALLOCATE(ev(norbs), norbs)
    call lalg_geneigensolve(norbs, lcao_data%X(hamilt) (1:norbs, 1:norbs, ik), &
         lcao_data%X(s) (1:norbs, 1:norbs, ik), ev)

    do n1 = start, nst
      st%eigenval(n1, ik) = ev(n1)
      st%X(psi)(:, :, n1, ik) = R_TOTYPE(M_ZERO)
    end do
    deallocate(ev)

    ! Change of base
    do n1 = start, nst
      do idim = 1, dim
        do n2 = 1, norbs
          call lalg_axpy(np, lcao_data%X(hamilt) (n2, n1, ik), lcao_data%st%X(psi)(:, idim, n2, ik), &
               st%X(psi)(:, idim, n1, ik))
        end do
      end do
    end do
    
  end do

  deallocate(hpsi)
  call pop_sub()
end subroutine X(lcao_wf)
