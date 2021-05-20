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


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_excited(mesh, namespace, space, tg, td, restart, kpoints)
    type(mesh_t),      intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(target_t),    intent(inout) :: tg
    type(td_t),        intent(in)    :: td
    type(restart_t),   intent(in)    :: restart
    type(kpoints_t),   intent(in)    :: kpoints

    integer :: ierr, nik, dim

    PUSH_SUB(target_init_excited)

    message(1) =  'Info: TargetOperator is a linear combination of Slater determinants.'
    call messages_info(1)

    tg%move_ions = ion_dynamics_ions_move(td%ions_dyn)
    tg%dt = td%dt

    call states_elec_look(restart, nik, dim, tg%st%nst, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to read states information."
      call messages_fatal(1)
    end if
    tg%st%st_start = 1
    tg%st%st_end   = tg%st%nst

    SAFE_DEALLOCATE_A(tg%st%occ)
    SAFE_DEALLOCATE_A(tg%st%eigenval)
    SAFE_DEALLOCATE_A(tg%st%node)

    SAFE_ALLOCATE(     tg%st%occ(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(tg%st%eigenval(1:tg%st%nst, 1:tg%st%d%nik))
    SAFE_ALLOCATE(    tg%st%node(1:tg%st%nst))
    if(tg%st%d%ispin == SPINORS) then
      SAFE_DEALLOCATE_A(tg%st%spin)
      SAFE_ALLOCATE(tg%st%spin(1:3, 1:tg%st%nst, 1:tg%st%d%nik))
    end if
    call states_elec_allocate_wfns(tg%st, mesh, TYPE_CMPLX)
    tg%st%node(:)  = 0

    call states_elec_load(restart, namespace, space, tg%st, mesh, kpoints, ierr)
    if (ierr /= 0) then
      message(1) = "Unable to read wavefunctions."
      call messages_fatal(1)
    end if

    call excited_states_init(tg%est, tg%st, "oct-excited-state-target", namespace) 

    POP_SUB(target_init_excited)
  end subroutine target_init_excited


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_excited()
    PUSH_SUB(target_end_excited)

    POP_SUB(target_end_excited)
  end subroutine target_end_excited


  ! ----------------------------------------------------------------------
  subroutine target_output_excited(tg, namespace, space, gr, dir, ions, hm, outp)
    type(target_t),      intent(in)  :: tg
    type(namespace_t),   intent(in)  :: namespace
    type(space_t),       intent(in)    :: space
    type(grid_t),        intent(in)  :: gr
    character(len=*),    intent(in)  :: dir
    type(ions_t),        intent(in)  :: ions
    type(hamiltonian_elec_t), intent(in)  :: hm
    type(output_t),      intent(in)  :: outp

    PUSH_SUB(target_output_excited)
    
    call io_mkdir(trim(dir), namespace)
    call output_states(outp, namespace, space, trim(dir)//'/st', tg%est%st, gr, ions, hm, -1)
    call excited_states_output(tg%est, trim(dir), namespace)

    POP_SUB(target_output_excited)
  end subroutine target_output_excited
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_excited(tg, namespace, gr, psi) result(j1)
    type(target_t),      intent(in) :: tg
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(in) :: gr
    type(states_elec_t), intent(in) :: psi

    PUSH_SUB(target_j1_excited)

    j1 = abs(zstates_elec_mpdotp(namespace, gr%mesh, tg%est, psi))**2

    POP_SUB(target_j1_excited)
  end function target_j1_excited


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_excited(tg, namespace, gr, psi_in, chi_out)
    type(target_t),         intent(in)    :: tg
    type(namespace_t),      intent(in)    :: namespace
    type(grid_t),           intent(in)    :: gr
    type(states_elec_t),    intent(in)    :: psi_in
    type(states_elec_t),    intent(inout) :: chi_out

    CMPLX, allocatable :: cI(:), dI(:), mat(:, :, :), mm(:, :, :, :), mk(:, :), lambda(:, :)
    CMPLX, allocatable :: zpsi(:, :), zchi(:, :)
    integer :: ik, ist, jst, ia, ib, n_pairs, nst, kpoints, jj, idim, ip
    PUSH_SUB(target_chi_excited)

    n_pairs = tg%est%n_pairs
    kpoints = psi_in%d%nik
    nst = psi_in%nst

    
    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:psi_in%d%dim))
    SAFE_ALLOCATE(zchi(1:gr%mesh%np, 1:psi_in%d%dim))
    SAFE_ALLOCATE(cI(1:n_pairs))
    SAFE_ALLOCATE(dI(1:n_pairs))
    SAFE_ALLOCATE(mat(1:tg%est%st%nst, 1:nst, 1:psi_in%d%nik))
    SAFE_ALLOCATE(mm(1:nst, 1:nst, 1:kpoints, 1:n_pairs))
    SAFE_ALLOCATE(mk(1:gr%mesh%np_part, 1:psi_in%d%dim))
    SAFE_ALLOCATE(lambda(1:n_pairs, 1:n_pairs))

    call zstates_elec_matrix(tg%est%st, psi_in, gr%mesh, mat)

    do ia = 1, n_pairs
      cI(ia) = tg%est%weight(ia)
      call zstates_elec_matrix_swap(mat, tg%est%pair(ia))
      mm(1:nst, 1:nst, 1:kpoints, ia) = mat(1:nst, 1:kpoints, 1:kpoints)
      dI(ia) = zstates_elec_mpdotp(namespace, gr%mesh, tg%est%st, psi_in, mat)
      if(abs(dI(ia)) > CNST(1.0e-12)) then
        do ik = 1, kpoints
          call lalg_inverter(nst, mm(1:nst, 1:nst, ik, ia))
        end do
      end if
      call zstates_elec_matrix_swap(mat, tg%est%pair(ia))
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
          
          zchi(1:gr%mesh%np, 1:psi_in%d%dim) = M_z0

          do ia = 1, n_pairs
            if(ik /= tg%est%pair(ia)%kk) cycle
            if(abs(dI(ia)) < CNST(1.0e-12)) cycle
            do ib = 1, n_pairs
              if(abs(dI(ib)) < CNST(1.0e-12)) cycle
              mk = M_z0

              do jst = 1, nst
                if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
                call states_elec_get_state(tg%est%st, gr%mesh, jj, ik, zpsi)

                do idim = 1, psi_in%d%dim
                  do ip = 1, gr%mesh%np
                    mk(ip, idim) = mk(ip, idim) + conjg(mm(ist, jst, ik, ib))*zpsi(ip, idim)
                  end do
                end do
              end do

              call lalg_axpy(gr%mesh%np_part, psi_in%d%dim, M_z1, lambda(ib, ia)*mk(:, :), zchi)

            end do
          end do

          call states_elec_set_state(chi_out, gr%mesh, ist, ik, zchi)
          
        end do
      end do
        
    case(SPINORS)
      ASSERT(chi_out%d%nik  ==  1)

      do ist = chi_out%st_start, chi_out%st_end
        
        zchi(1:gr%mesh%np, 1:psi_in%d%dim) = M_z0
        
        do ia = 1, n_pairs
          if(abs(dI(ia)) < CNST(1.0e-12)) cycle

          do ib = 1, n_pairs
            if(abs(dI(ib)) < CNST(1.0e-12)) cycle

            mk = M_z0
            do jst = 1, nst
              if(jst  ==  tg%est%pair(ib)%i) jj = tg%est%pair(ia)%a
              call states_elec_get_state(tg%est%st, gr%mesh, jj, ik, zpsi)
              
              do idim = 1, psi_in%d%dim
                do ip = 1, gr%mesh%np
                  mk(ip, idim) = mk(ip, idim) + conjg(mm(ist, jst, 1, ib))*zpsi(ip, idim)
                end do
              end do
            end do

            call lalg_axpy(gr%mesh%np_part, 2, M_z1, lambda(ib, ia)*mk(:, :), zchi)
          end do
        end do

        call states_elec_set_state(chi_out, gr%mesh, ist, ik, zchi)
        
      end do

    end select

    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zchi)
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
