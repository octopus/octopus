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
!! $Id: subspace_inc.F90 6258 2009-12-26 20:07:24Z xavier $

! ---------------------------------------------------------
! This routine diagonalises the Hamiltonian in the subspace defined by the states.
subroutine X(subspace_diag)(der, st, hm, ik, eigenval, psi, diff)
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(:, :, :)
  FLOAT, optional,     intent(out)   :: diff(:)

  R_TYPE, allocatable :: h_subspace(:, :), f(:, :, :)
  integer             :: ist, ist2, jst, size, idim
  FLOAT               :: nrm2
  type(profile_t),     save    :: diagon_prof
  type(batch_t) :: psib, hpsib, whole_psib

  PUSH_SUB(X(subspace_diag))
  call profiling_in(diagon_prof, "SUBSPACE_DIAG")

#ifdef HAVE_MPI
  if(st%parallel_in_states) then
    if(present(diff)) then
      call X(subspace_diag_par_states)(der, st, hm, ik, eigenval, psi, diff)
    else
      call X(subspace_diag_par_states)(der, st, hm, ik, eigenval, psi)
    end if
  else
#endif

    SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(f(1:der%mesh%np, 1:st%d%dim, 1:st%d%block_size))

    ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
    do ist = st%st_start, st%st_end, st%d%block_size
      size = min(st%d%block_size, st%st_end - ist + 1)

      call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))

      call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)

      call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)

      call batch_end(psib)

      call batch_init(whole_psib, hm%d%dim, ist, st%nst, psi(:, :, ist:st%nst))
      call X(mesh_batch_dotp_matrix)(der%mesh, hpsib, whole_psib, h_subspace)
      call batch_end(whole_psib)
      call batch_end(hpsib)
      
    end do

    do ist = st%st_start, st%st_end
      do jst = ist, st%st_end
        h_subspace(jst, ist) = R_CONJ(h_subspace(ist, jst))
      end do
    end do

    ! Diagonalize the Hamiltonian in the subspace.
    call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

    ! Calculate the new eigenfunctions as a linear combination of the
    ! old ones.
    call batch_init(whole_psib, hm%d%dim, 1, st%nst, psi(:, :, :))
    call X(mesh_batch_rotate)(der%mesh, whole_psib, h_subspace)
    call batch_end(whole_psib)

    ! Renormalize.
    do ist = st%st_start, st%st_end
      nrm2 = X(mf_nrm2)(der%mesh, st%d%dim, psi(:, :, ist))
      do idim = 1, st%d%dim
        call lalg_scal(der%mesh%np, M_ONE/nrm2, psi(:, idim, ist))
      end do
    end do

    ! Recalculate the residues if requested by the diff argument.
    if(present(diff)) then 
      
      do ist = st%st_start, st%st_end, st%d%block_size
        size = min(st%d%block_size, st%st_end - ist + 1)

        call batch_init(psib, hm%d%dim, ist, ist + size - 1, psi(:, :, ist:))
        call batch_init(hpsib, hm%d%dim, ist, ist + size - 1, f)
        
        call X(hamiltonian_apply_batch)(hm, der, psib, hpsib, ik)

        call batch_end(psib)
        call batch_end(hpsib)
        
        do ist2 = ist, ist + size - 1
          diff(ist2) = X(states_residue)(der%mesh, st%d%dim, f(:, :, ist2 - ist + 1), eigenval(ist2), psi(:, :, ist2))
        end do
      end do

    end if

    SAFE_DEALLOCATE_A(f)
    SAFE_DEALLOCATE_A(h_subspace)

#ifdef HAVE_MPI
  end if
#endif

  call profiling_out(diagon_prof)
  POP_SUB(X(subspace_diag))

end subroutine X(subspace_diag)

#ifdef HAVE_MPI
! --------------------------------------------------------- 
! This routine diagonalises the Hamiltonian in the subspace defined by
! the states; this version is aware of parallelization in states but
! consumes more memory.
!
subroutine X(subspace_diag_par_states)(der, st, hm, ik, eigenval, psi, diff)
  type(derivatives_t), intent(in)    :: der
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  FLOAT,               intent(out)   :: eigenval(:)
  R_TYPE,              intent(inout) :: psi(1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end)
  FLOAT, optional,     intent(out)   :: diff(1:st%nst)

  R_TYPE, allocatable :: h_subspace(:,:), ff(:,:,:)
  integer             :: ist
  FLOAT               :: nrm2
#if defined(HAVE_MPI)
  integer             :: tmp
  FLOAT               :: ldiff(st%lnst)
#endif

  PUSH_SUB(X(subspace_diag_par_states))

  SAFE_ALLOCATE(h_subspace(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(ff(1:der%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

  ! Calculate the matrix representation of the Hamiltonian in the subspace <psi|H|psi>.
  do ist = st%st_start, st%st_end
    call X(hamiltonian_apply)(hm, der, psi(:, :, ist), ff(:, :, ist), ist, ik)
  end do
  call states_blockt_mul(der%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       psi(:, :, :), ff, h_subspace, symm=.true.)

  ! Diagonalize the Hamiltonian in the subspace.
  call lalg_eigensolve(st%nst, h_subspace, eigenval(:))

  ! The new states are given by the eigenvectors of the matrix.
  ff(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end) = psi(1:der%mesh%np, 1:st%d%dim, st%st_start:st%st_end)

  call states_block_matr_mul(der%mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
       ff, h_subspace, psi(:, :, :))

  ! Renormalize.
  do ist = st%st_start, st%st_end
    nrm2 = X(mf_nrm2)(der%mesh, st%d%dim, psi(:, :, ist))
    psi(1:der%mesh%np, 1:st%d%dim, ist) = psi(1:der%mesh%np, 1:st%d%dim, ist)/nrm2
  end do

  ! Recalculate the residues if requested by the diff argument.
  if(present(diff)) then 
    do ist = st%st_start, st%st_end
      call X(hamiltonian_apply)(hm, der, psi(:, :, ist) , ff(:, :, st%st_start), ist, ik)
      diff(ist) = X(states_residue)(der%mesh, st%d%dim, ff(:, :, st%st_start), eigenval(ist), &
           psi(:, :, ist))
    end do

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      ldiff = diff(st%st_start:st%st_end)
      call lmpi_gen_allgatherv(st%lnst, ldiff, tmp, diff(:), st%mpi_grp)
    end if
#endif
  end if
  
  SAFE_DEALLOCATE_A(ff)
  SAFE_DEALLOCATE_A(h_subspace)
  
  POP_SUB(X(subspace_diag_par_states))
  
end subroutine X(subspace_diag_par_states) 
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
