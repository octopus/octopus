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
subroutine X(calculate_eigenvalues)(hm, gr, st, t)
  type(hamiltonian_t), intent(inout) :: hm
  type(grid_t) ,       intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  FLOAT,  optional,    intent(in)    :: t

  R_TYPE, allocatable :: Hpsi(:,:)
  R_TYPE :: e
  integer :: ik, ist

  call push_sub('energy_inc.Xcalculate_eigenvalues')
  SAFE_ALLOCATE(Hpsi(1:gr%mesh%np, 1:st%d%dim))

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
      do ist = st%st_start, st%st_end
        if(present(t)) then
          call X(hamiltonian_apply) (hm, gr, st%X(psi)(:, :, ist, ik), hpsi, ist, ik, t)
        else
          call X(hamiltonian_apply) (hm, gr, st%X(psi)(:, :, ist, ik), hpsi, ist, ik)
        end if
        e = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), Hpsi)
        st%eigenval(ist, ik) = R_REAL(e)
      end do
    end do

#ifdef HAVE_MPI
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      message(1) = "Fixme: eigenvalues are not allgathered."
      call write_warning(1)
    end if
#endif
  end if

  SAFE_DEALLOCATE_A(Hpsi)
  call pop_sub()
end subroutine X(calculate_eigenvalues)

! ---------------------------------------------------------
FLOAT function X(electronic_kinetic_energy)(hm, gr, st) result(t0)
  type(hamiltonian_t), intent(in)    :: hm
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st

  integer :: ik, ist
  R_TYPE, allocatable :: tpsi(:, :)
  FLOAT, allocatable  :: t(:, :)

  call push_sub('energy_inc.Xelectronic_kinetic_energy')

  SAFE_ALLOCATE(tpsi(1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(t(st%st_start:st%st_end, 1:st%d%nik))
  t = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      tpsi = R_TOTYPE(M_ZERO)
      call X(hamiltonian_apply)(hm, gr, st%X(psi)(:, :, ist, ik), tpsi, ist, ik, kinetic_only = .true.)
      t(ist, ik) = X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, ist, ik), tpsi)
    end do
  end do
  
  t0 = states_eigenvalues_sum(st, t)

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
