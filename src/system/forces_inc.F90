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

subroutine X(forces_gather)(geo, force)
  type(geometry_t), intent(in)    :: geo
  R_TYPE,           intent(inout) :: force(:, :)
  
  R_TYPE,  allocatable :: force_local(:, :)
  integer, allocatable :: recv_count(:), recv_displ(:)

  PUSH_SUB(X(forces_gather))

  call profiling_in(prof_comm, "FORCES_COMM")
    
  ! each node has a piece of the force array, they have to be
  ! collected by all nodes
  
  ASSERT(ubound(force, dim = 1) == geo%space%dim)

  SAFE_ALLOCATE(recv_count(1:geo%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(recv_displ(1:geo%atoms_dist%mpi_grp%size))
  SAFE_ALLOCATE(force_local(1:geo%space%dim, 1:max(1, geo%atoms_dist%nlocal)))
  
  recv_count(1:geo%atoms_dist%mpi_grp%size) = geo%space%dim*geo%atoms_dist%num(0:geo%atoms_dist%mpi_grp%size - 1)
  recv_displ(1:geo%atoms_dist%mpi_grp%size) = geo%space%dim*(geo%atoms_dist%range(1, 0:geo%atoms_dist%mpi_grp%size - 1) - 1)
  
  if(geo%atoms_dist%nlocal > 0) then
    force_local(1:geo%space%dim, 1:geo%atoms_dist%nlocal) = force(1:geo%space%dim, geo%atoms_dist%start:geo%atoms_dist%end)
  endif

#ifdef HAVE_MPI  
  call MPI_Allgatherv(&
    force_local(1, 1), geo%space%dim*geo%atoms_dist%nlocal, R_MPITYPE, &
    force(1, 1), recv_count(1), recv_displ(1), R_MPITYPE, &
    geo%atoms_dist%mpi_grp%comm, mpi_err)
#endif
  
  SAFE_DEALLOCATE_A(recv_count)
  SAFE_DEALLOCATE_A(recv_displ)
  SAFE_DEALLOCATE_A(force_local)

  call profiling_out(prof_comm)
  POP_SUB(X(forces_gather))
end subroutine X(forces_gather)

!---------------------------------------------------------------------------
subroutine X(forces_from_local_potential)(gr, geo, ep, st, time, gdensity, force)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(in)    :: time
  R_TYPE,                         intent(in)    :: gdensity(:, :)
  R_TYPE,                         intent(inout) :: force(:, :)

  FLOAT,  allocatable :: vloc(:)
  R_TYPE, pointer     :: zvloc(:)
  integer             :: ip, idir, iatom
 
  PUSH_SUB(X(forces_from_local_potential))

  SAFE_ALLOCATE(vloc(1:gr%mesh%np))
  SAFE_ALLOCATE(zvloc(1:gr%mesh%np))
  
  do iatom = geo%atoms_dist%start, geo%atoms_dist%end

    if(.not.simul_box_in_box(gr%mesh%sb, geo, geo%atom(iatom)%x) .and. ep%ignore_external_ions) then
      force(1:gr%mesh%sb%dim, iatom) = M_ZERO
      cycle
    end if
    
    vloc(1:gr%mesh%np) = M_ZERO
    
    call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, vloc, time)

    forall(ip = 1:gr%mesh%np) zvloc(ip) = vloc(ip)

    do idir = 1, gr%mesh%sb%dim
      force(idir, iatom) = -X(mf_dotp)(gr%mesh, zvloc, gdensity(:, idir))
    end do

  end do

  if(geo%atoms_dist%parallel) call X(forces_gather)(geo, force)
  !if(geo%atoms_dist%parallel .and. geo%atoms_dist%nlocal > 0) call X(forces_gather)(geo, force)

  SAFE_DEALLOCATE_A(vloc)
  SAFE_DEALLOCATE_P(zvloc)

  POP_SUB(X(forces_from_local_potential))
end subroutine X(forces_from_local_potential)

!---------------------------------------------------------------------------
! Ref: Kikuji Hirose, Tomoya Ono, Yoshitaka Fujimoto, and Shigeru Tsukamoto,
! First-principles calculations in real-space formalism: Electronic configurations
! and transport properties of nanostructures, Imperial College Press (2005)
! Section 1.6, page 12
subroutine X(forces_from_potential)(gr, geo, ep, st, time)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(in)    :: time
 
  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  FLOAT,  allocatable :: grad_rho(:, :), force(:, :)
  CMPLX :: phase
#ifdef HAVE_MPI
  FLOAT, allocatable  :: force_local(:, :), grad_rho_local(:, :)
#endif

  PUSH_SUB(X(forces_from_potential))

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho = M_ZERO
  SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call states_get_state(st, gr%mesh, ist, iq, psi)

      do idim = 1, st%d%dim
        call X(derivatives_set_bc)(gr%der, psi(:, idim))

        ikpoint = states_dim_get_kpoint_index(st%d, iq)
        if(simul_box_is_periodic(gr%sb) .and. .not. kpoints_point_is_gamma(gr%sb%kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, ikpoint)

          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:gr%sb%dim)*gr%mesh%x(ip, 1:gr%sb%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
          end do
        endif

        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)

        ff = st%d%kweights(iq) * st%occ(ist, iq) * M_TWO
        do idir = 1, gr%mesh%sb%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff*R_REAL(R_CONJ(psi(ip, idim))*grad_psi(ip, idir, idim))
          end do
        end do

      end do

      call profiling_count_operations(np*st%d%dim*gr%mesh%sb%dim*(2 + R_MUL))

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, gr%mesh%sb%dim

          force(idir, iatom) = force(idir, iatom) - M_TWO * st%d%kweights(iq) * st%occ(ist, iq) * &
            R_REAL(X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, psi, grad_psi(:, idir, :)))

        end do
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(grad_psi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, "FORCES_COMM")
    call comm_allreduce(st%st_kpt_mpi_grp%comm, force, dim = (/gr%mesh%sb%dim, geo%natoms/))
    call comm_allreduce(st%st_kpt_mpi_grp%comm, grad_rho, dim = (/np, gr%mesh%sb%dim/))
    call profiling_out(prof_comm)
  end if
#endif

  do iatom = 1, geo%natoms
    geo%atom(iatom)%f(1:gr%mesh%sb%dim) = geo%atom(iatom)%f(1:gr%mesh%sb%dim) + force(1:gr%mesh%sb%dim, iatom)
  end do

  call dforces_from_local_potential(gr, geo, ep, st, time, grad_rho, force)

  do iatom = 1, geo%natoms
    do idir = 1, gr%mesh%sb%dim
      geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + force(idir, iatom)
    end do
  end do

  SAFE_DEALLOCATE_A(force)
  POP_SUB(X(forces_from_potential))
end subroutine X(forces_from_potential)

! --------------------------------------------------------------------------------
subroutine X(forces_born_charges)(gr, geo, ep, st, time, lr, lr2, lr_dir, born_charges)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(inout) :: ep
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(in)    :: time
  type(lr_t),                     intent(inout) :: lr
  type(lr_t),                     intent(inout) :: lr2
  integer,                        intent(in)    :: lr_dir
  type(born_charges_t),           intent(out)   :: born_charges

  ! lr, lr2 should be the wfns from electric perturbation in the lr_dir direction
  ! lr is for +omega, lr2 is for -omega.
  ! for each atom, Z*(i,j) = dF(j)/dE(i)

  integer :: iatom, ist, iq, idim, idir, np, np_part, ip, ikpoint
  FLOAT :: ff, kpoint(1:MAX_DIM)
  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: dl_psi(:, :)
  R_TYPE, allocatable :: dl_psi2(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi2(:, :, :)
  CMPLX,  allocatable :: grad_rho(:, :), force(:, :)
  CMPLX :: phase
#ifdef HAVE_MPI
  CMPLX, allocatable  :: force_local(:, :), grad_rho_local(:, :)
#endif

  PUSH_SUB(X(forces_born_charges))

  ! need all to calculate Born charges
  ASSERT(lr_dir > 0 .and. lr_dir <= gr%mesh%sb%dim)

  np      = gr%mesh%np
  np_part = gr%mesh%np_part

  SAFE_ALLOCATE(grad_dl_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_dl_psi2(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_psi(1:np, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho = M_ZERO
  SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi(1:np_part, 1:st%d%dim))
  SAFE_ALLOCATE(dl_psi2(1:np_part, 1:st%d%dim))

  !THE NON-LOCAL PART (parallel in states and k-points)
  do iq = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        call states_get_state(st, gr%mesh, idim, ist, iq, psi(:, idim))
        call X(derivatives_set_bc)(gr%der, psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr%X(dl_psi)(:, idim, ist, iq), dl_psi(:, idim))
        call X(derivatives_set_bc)(gr%der, dl_psi(:, idim))
        call lalg_copy(gr%mesh%np_part, lr2%X(dl_psi)(:, idim, ist, iq), dl_psi2(:, idim))
        call X(derivatives_set_bc)(gr%der, dl_psi2(:, idim))

        ikpoint = states_dim_get_kpoint_index(st%d, iq)
        if(simul_box_is_periodic(gr%sb) .and. .not. kpoints_point_is_gamma(gr%sb%kpoints, ikpoint)) then

          kpoint = M_ZERO
          kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, ikpoint)

          do ip = 1, np_part
            phase = exp(-M_zI*sum(kpoint(1:gr%sb%dim)*gr%mesh%x(ip, 1:gr%sb%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
            dl_psi(ip, idim) = phase*dl_psi(ip, idim)
            dl_psi2(ip, idim) = phase*dl_psi2(ip, idim)
          end do
        endif
       
        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi(:, idim), grad_dl_psi(:, :, idim), set_bc = .false.)
        call X(derivatives_grad)(gr%der, dl_psi2(:, idim), grad_dl_psi2(:, :, idim), set_bc = .false.)

        !accumulate to calculate the gradient of the density
        ff = st%d%kweights(iq) * st%occ(ist, iq)
        do idir = 1, gr%mesh%sb%dim
          do ip = 1, np
            grad_rho(ip, idir) = grad_rho(ip, idir) + ff * &
              (R_CONJ(grad_psi(ip, idir, idim)) * dl_psi(ip, idim) + R_CONJ(psi(ip, idim)) * grad_dl_psi(ip, idir, idim) &
              + R_CONJ(dl_psi2(ip, idim)) * grad_psi(ip, idir, idim) + R_CONJ(grad_dl_psi2(ip, idir, idim)) * psi(ip, idim))
          end do
        end do
      end do

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, gr%mesh%sb%dim

          force(idir, iatom) = force(idir, iatom) - st%d%kweights(iq) * st%occ(ist, iq) * &
            (X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, grad_psi(:, idir, :), dl_psi) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, psi, grad_dl_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, dl_psi2, grad_psi(:, idir, :)) &
            + X(projector_matrix_element)(ep%proj(iatom), st%d%dim, iq, grad_dl_psi2(:, idir, :), psi))
          
        end do
      end do

    end do
  end do
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(dl_psi)
  SAFE_DEALLOCATE_A(dl_psi2)
  SAFE_DEALLOCATE_A(grad_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi)
  SAFE_DEALLOCATE_A(grad_dl_psi2)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call profiling_in(prof_comm, "FORCES_COMM")
    call comm_allreduce(st%st_kpt_mpi_grp%comm, force, dim = (/gr%mesh%sb%dim, geo%natoms/))
    call comm_allreduce(st%st_kpt_mpi_grp%comm, grad_rho, dim = (/np, gr%mesh%sb%dim/))
    call profiling_out(prof_comm)
  end if
#endif
  
  do iatom = 1, geo%natoms
    born_charges%charge(lr_dir, 1:gr%mesh%sb%dim, iatom) = force(1:gr%mesh%sb%dim, iatom)
  end do

  call zforces_from_local_potential(gr, geo, ep, st, time, grad_rho, force)

  SAFE_DEALLOCATE_A(grad_rho)

  do iatom = 1, geo%natoms
    born_charges%charge(lr_dir, lr_dir, iatom) = born_charges%charge(lr_dir, lr_dir, iatom) + species_zval(geo%atom(iatom)%spec)
    do idir = 1, gr%mesh%sb%dim
      born_charges%charge(lr_dir, idir, iatom) = born_charges%charge(lr_dir, idir, iatom) + force(idir, iatom)
    enddo
  enddo

  SAFE_DEALLOCATE_A(force)
  POP_SUB(X(forces_born_charges))
  
end subroutine X(forces_born_charges)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
