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
!! $Id: target_excited_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_excited(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    integer :: ierr, ip
    PUSH_SUB(target_init_excited)

    message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
    call messages_info(1)

    call states_look (trim(restart_dir)//GS_DIR, gr%mesh%mpi_grp, ip, ip, tg%st%nst, ierr)
    tg%st%st_start = 1
    tg%st%st_end   = tg%st%nst

    SAFE_DEALLOCATE_P(tg%st%occ)
    SAFE_DEALLOCATE_P(tg%st%eigenval)
    SAFE_DEALLOCATE_P(tg%st%node)

    SAFE_ALLOCATE(     tg%st%occ(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(tg%st%eigenval(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(    tg%st%node(1:tg%st%nst))
    if(tg%st%d%ispin == SPINORS) then
      SAFE_DEALLOCATE_P(tg%st%spin)
      SAFE_ALLOCATE(tg%st%spin(1:3, 1:tg%st%nst, 1:tg%st%d%nik))
    end if
    call states_allocate_wfns(tg%st, gr%mesh, TYPE_CMPLX)
    tg%st%node(:)  = 0

    call restart_read(trim(restart_dir)//GS_DIR, tg%st, gr, ierr, exact = .true.)
    call excited_states_init(tg%est, tg%st, "oct-excited-state-target") 

    POP_SUB(target_init_excited)
  end subroutine target_init_excited


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_excited()
    PUSH_SUB(target_end_excited)

    POP_SUB(target_end_excited)
  end subroutine target_end_excited


  ! ----------------------------------------------------------------------
  subroutine target_output_excited(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    PUSH_SUB(target_output_excited)
    
    call loct_mkdir(trim(dir))
    call output_states(tg%est%st, gr, geo, trim(dir)//'/st', outp)
    call excited_states_output(tg%est, trim(dir))

    POP_SUB(target_output_excited)
  end subroutine target_output_excited
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_excited(tg, gr, psi) result(j1)
    type(target_t),   intent(inout) :: tg
    type(grid_t),     intent(inout) :: gr
    type(states_t),   intent(inout) :: psi

    PUSH_SUB(target_j1_excited)

    j1 = abs(zstates_mpdotp(gr%mesh, tg%est, psi))**2

    POP_SUB(target_j1_excited)
  end function target_j1_excited


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_excited(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out

    CMPLX :: zdet
    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    integer :: ik, ist, jst, ia, ib, n_pairs, nst, kpoints, jj
    PUSH_SUB(target_chi_excited)

    n_pairs = tg%est%n_pairs
    kpoints = psi_in%d%nik
    nst = psi_in%nst

    SAFE_ALLOCATE(cI(1:n_pairs))
    SAFE_ALLOCATE(dI(1:n_pairs))
    SAFE_ALLOCATE(mat(1:tg%est%st%nst, 1:nst, 1:psi_in%d%nik))
    SAFE_ALLOCATE(mm(1:nst, 1:nst, 1:kpoints, 1:n_pairs))
    SAFE_ALLOCATE(mk(1:gr%mesh%np_part, 1:psi_in%d%dim))
    SAFE_ALLOCATE(lambda(1:n_pairs, 1:n_pairs))

    call zstates_matrix(gr%mesh, tg%est%st, psi_in, mat)

    do ia = 1, n_pairs
      cI(ia) = tg%est%weight(ia)
      call zstates_matrix_swap(mat, tg%est%pair(ia))
      mm(1:nst, 1:nst, 1:kpoints, ia) = mat(1:nst, 1:kpoints, 1:kpoints)
      dI(ia) = zstates_mpdotp(gr%mesh, tg%est%st, psi_in, mat)
      if(abs(dI(ia)) > CNST(1.0e-12)) then
        do ik = 1, kpoints
          zdet = lalg_inverter(nst, mm(1:nst, 1:nst, ik, ia))
        end do
      end if
      call zstates_matrix_swap(mat, tg%est%pair(ia))
    end do

    do ia = 1, n_pairs
      do ib = 1, n_pairs
        lambda(ia, ib) = conjg(cI(ib)) * cI(ia) * conjg(dI(ia)) * dI(ib)
      end do
    end do

    select case(psi_in%d%ispin)
    case(UNPOLARIZED)
      write(message(1), '(a)') 'Internal error in target.target_chi: unpolarized.'
      call messages_fatal(1)

    case(SPIN_POLARIZED)
      ASSERT(chi_out%d%nik  ==  2)
      do ik = 1, kpoints
        do ist = chi_out%st_start, chi_out%st_end
          chi_out%zpsi(:, :, ist, ik) = M_z0
          do ia = 1, n_pairs
            if(ik /= tg%est%pair(ia)%sigma) cycle
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle
            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle
              mk = M_z0
              do jst = 1, nst
                if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
                mk(:, :) = mk(:, :) + conjg(mm(ist, jst, ik, ib)) * tg%est%st%zpsi(:, :, jj, ik)
              end do
              call lalg_axpy(gr%mesh%np_part, psi_in%d%dim, M_z1, lambda(ib, ia) * mk(:, :), chi_out%zpsi(:, :, ist, ik))
            end do
          end do
        end do
      end do
        
    case(SPINORS)
      ASSERT(chi_out%d%nik  ==  1)

      do ist = chi_out%st_start, chi_out%st_end
        chi_out%zpsi(:, :, ist, 1) = M_z0

        do ia = 1, n_pairs
          if(abs(dI(ia)) < CNST(1.0e-12)) cycle

          do ib = 1, n_pairs
            if(abs(dI(ib)) < CNST(1.0e-12)) cycle

            mk = M_z0
            do jst = 1, nst
              if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
              mk(:, :) = mk(:, :) + conjg(mm(ist, jst, 1, ib)) * tg%est%st%zpsi(:, :, jj, 1)
            end do

            call lalg_axpy(gr%mesh%np_part, 2, M_z1, lambda(ib, ia) * mk(:, :), chi_out%zpsi(:, :, ist, 1))
          end do
        end do
      end do

    end select

    SAFE_DEALLOCATE_A(cI)
    SAFE_DEALLOCATE_A(dI)
    SAFE_DEALLOCATE_A(mat)
    SAFE_DEALLOCATE_A(mm)
    SAFE_DEALLOCATE_A(mk)
    SAFE_DEALLOCATE_A(lambda)
    POP_SUB(target_chi_excited)
  end subroutine target_chi_excited


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
