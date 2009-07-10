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
! calculates the eigenvalues of the real orbitals
subroutine X(calculate_eigenvalues)(hm, gr, st, time)
  type(hamiltonian_t), intent(inout) :: hm
  type(grid_t) ,       intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  FLOAT,  optional,    intent(in)    :: time

  R_TYPE, allocatable :: hpsi(:, :, :), psi(:, :, :)
  integer :: ik, ist, minst, maxst, idim
  type(batch_t) :: psib, hpsib
  type(profile_t), save :: prof

  call push_sub('energy_inc.Xcalculate_eigenvalues')
  call profiling_in(prof, "EIGENVALUE_CALC")

  if(in_debug_mode) then
    write(message(1), '(a)') 'Debug: Calculating eigenvalues.'
    call write_info(1)
  end if

  if(gr%sb%open_boundaries.and.calc_mode_is(CM_GS)) then
    ! For open boundaries we know the eigenvalues.
    st%eigenval = st%ob_eigenval
  else
    ! FIXME: for TD open boundaries this is wrong. But the GS case like above
    ! is also wrong.
    ! The correct way to calculate the eigenvalue here is:
    !      / Psi_L | H_LL H_LC 0    | Psi_L \
    ! e = <  Psi_C | H_CL H_CC H_CR | Psi_C  >
    !      \ Psi_R | 0    H_RC H_RR | Psi_R /
    ! But I am not sure how to calculate this right now.

    do ik = st%d%kpt%start, st%d%kpt%end

      SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size))
      hpsi(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size) = R_TOTYPE(M_ZERO)
      if(st%np_size) then
        SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim, 1:st%d%block_size))
        psi(1:gr%mesh%np, 1:st%d%dim, 1:st%d%block_size) = R_TOTYPE(M_ZERO)
      end if

      do minst = st%st_start, st%st_end, st%d%block_size
        maxst = min(st%st_end, minst + st%d%block_size - 1)

        if(st%np_size) then
          call batch_init(psib, st%d%dim, minst, maxst, psi)
          call batch_set(psib, gr%mesh%np,  st%X(psi)(:, :, minst:, ik))
        else
          call batch_init(psib, st%d%dim, minst, maxst, st%X(psi)(:, :, minst:, ik))
        end if

        call batch_init(hpsib, st%d%dim, minst, maxst, hpsi)

        if(present(time)) then
          call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik, time)
        else
          call X(hamiltonian_apply_batch)(hm, gr, psib, hpsib, ik)
        end if

        call batch_end(psib)
        call batch_end(hpsib)
        
        do ist = minst, maxst 
          st%eigenval(ist, ik) = &
            R_REAL(X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), hpsi(:, :, ist - minst + 1)))
        end do

      end do

      if(st%np_size) then
        SAFE_DEALLOCATE_A(psi)
      end if
      SAFE_DEALLOCATE_A(hpsi)

    end do

  end if

  call profiling_out(prof)
  call pop_sub()
end subroutine X(calculate_eigenvalues)

! ---------------------------------------------------------
FLOAT function X(electronic_kinetic_energy)(hm, gr, st) result(t0)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist, idim
  R_TYPE, allocatable :: tpsi(:, :), psi(:, :)
  FLOAT, allocatable  :: t(:, :)

  call push_sub('energy_inc.Xelectronic_kinetic_energy')

  if(st%np_size) then
    SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:st%d%dim))
  end if
  SAFE_ALLOCATE(tpsi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(t(st%st_start:st%st_end, 1:st%d%nik))

  t = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      tpsi = R_TOTYPE(M_ZERO)
      if(st%np_size) then
        do idim = 1, st%d%dim
          call lalg_copy(gr%mesh%np, st%X(psi)(:, idim, ist, ik), psi(:, idim))
        end do
        call X(hamiltonian_apply)(hm, gr, psi, tpsi, ist, ik, kinetic_only = .true.)
      else
        call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, ist, ik), tpsi, ist, ik, kinetic_only = .true.)
      end if
      t(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), tpsi)
    end do
  end do
  
  t0 = states_eigenvalues_sum(st, t)
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(tpsi)
  SAFE_DEALLOCATE_A(t)
  call pop_sub()
end function X(electronic_kinetic_energy)

! ---------------------------------------------------------
FLOAT function X(electronic_external_energy)(hm, gr, st) result(v)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: vpsi(:, :)
  FLOAT, allocatable :: t(:, :)

  call push_sub('energy_inc.Xelectronic_external_energy')

  SAFE_ALLOCATE(vpsi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(t(st%st_start:st%st_end, 1:st%d%nik))
  t = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      vpsi = R_TOTYPE(M_ZERO)
      call X(vexternal) (hm, gr, st%X(psi)(:, :, ist, ik), vpsi, ik)
      t(ist, ik) = X(mf_dotp) (gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), vpsi)
    end do
  end do

  v = states_eigenvalues_sum(st, t)

  SAFE_DEALLOCATE_A(vpsi)
  SAFE_DEALLOCATE_A(t)
  call pop_sub()
end function X(electronic_external_energy)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
