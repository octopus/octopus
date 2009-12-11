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

subroutine X(forces_from_potential)(gr, geo, ep, st, time, lr, lr2, lr_dir, Born_charges)
  type(grid_t),                   intent(inout) :: gr
  type(geometry_t),               intent(inout) :: geo
  type(epot_t),                   intent(in)    :: ep
  type(states_t),                 intent(inout) :: st
  FLOAT,                          intent(in)    :: time
  type(lr_t),           optional, intent(inout) :: lr
  type(lr_t),           optional, intent(inout) :: lr2
  integer,              optional, intent(in)    :: lr_dir
  type(Born_charges_t), optional, intent(out)   :: Born_charges
  ! provide these optional arguments to calculate Born effective charges rather than forces
  ! lr, lr2 should be the wfns from electric perturbation in the lr_dir direction
  ! lr is for +omega, lr2 is for -omega.
  ! for each atom, Z*(i,j) = dF(j)/dE(i)

  integer :: iatom, ist, ik, idim, idir, np, np_part, ip
  FLOAT :: ff

  R_TYPE, allocatable :: psi(:, :)
  R_TYPE, allocatable :: dl_psi(:, :)
  R_TYPE, allocatable :: dl_psi2(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi2(:, :, :)
  FLOAT,  allocatable :: vloc(:)
  CMPLX,  allocatable :: grad_rho(:, :), force(:, :), zvloc(:)
  CMPLX :: phase
#ifdef HAVE_MPI
  integer, allocatable :: recv_count(:), recv_displ(:)
  CMPLX, allocatable  :: force_local(:, :), grad_rho_local(:, :)
  type(profile_t), save :: prof
#endif

  call push_sub('forces_inc.Xforces_from_potential')

  ASSERT(present(lr) .eqv. present(lr_dir))
  ASSERT(present(lr) .eqv. present(lr2))
  ASSERT(present(lr) .eqv. present(Born_charges))
  ! need all to calculate Born charges
  if(present(lr_dir)) then
    ASSERT(lr_dir > 0 .and. lr_dir <= gr%mesh%sb%dim)
  end if

  np = gr%mesh%np
  np_part = gr%mesh%np_part

  if(present(lr)) then
    SAFE_ALLOCATE( grad_dl_psi(1:np_part, 1:gr%mesh%sb%dim, 1:st%d%dim))
    SAFE_ALLOCATE(grad_dl_psi2(1:np_part, 1:gr%mesh%sb%dim, 1:st%d%dim))
  endif
  SAFE_ALLOCATE(grad_psi(1:np_part, 1:gr%mesh%sb%dim, 1:st%d%dim))
  SAFE_ALLOCATE(grad_rho(1:np, 1:gr%mesh%sb%dim))
  grad_rho(1:np, 1:gr%mesh%sb%dim) = M_ZERO
  SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
  force = M_ZERO

  ! even if there is no fine mesh, we need to make another copy
  SAFE_ALLOCATE(psi(1:np_part, 1:st%d%dim))
  if(present(lr)) then
    SAFE_ALLOCATE(dl_psi(1:np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dl_psi2(1:np_part, 1:st%d%dim))
  endif

  !THE NON-LOCAL PART (parallel in states and k-points)
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        call lalg_copy(gr%mesh%np_part, st%X(psi)(:, idim, ist, ik), psi(:, idim))
        call X(derivatives_set_bc)(gr%der, psi(:, idim))
        
        if (present(lr)) then
          call lalg_copy(gr%mesh%np_part, lr%X(dl_psi)(:, idim, ist, ik), dl_psi(:, idim))
          call X(derivatives_set_bc)(gr%der, dl_psi(:, idim))
          call lalg_copy(gr%mesh%np_part, lr2%X(dl_psi)(:, idim, ist, ik), dl_psi2(:, idim))
          call X(derivatives_set_bc)(gr%der, dl_psi2(:, idim))
        endif

        if(simul_box_is_periodic(gr%sb) .and. .not. kpoint_is_gamma(st%d, ik)) then
          do ip = 1, np_part
            phase = exp(-M_zI*sum(st%d%kpoints(1:gr%sb%dim, ik)*gr%mesh%x(ip, 1:gr%sb%dim)))
            psi(ip, idim) = phase*psi(ip, idim)
            if(present(lr)) then
              dl_psi(ip, idim) = phase*dl_psi(ip, idim)
              dl_psi2(ip, idim) = phase*dl_psi2(ip, idim)
            endif
          end do
        endif
       
        ! calculate the gradients of the wavefunctions
        ! and set boundary conditions in preparation for applying projectors
        call X(derivatives_grad)(gr%der, psi(:, idim), grad_psi(:, :, idim), set_bc = .false.)
        do idir = 1, gr%mesh%sb%dim
          call X(derivatives_set_bc)(gr%der, grad_psi(:, idir, idim))
        enddo

        if (present(lr)) then
          call X(derivatives_grad)(gr%der, dl_psi(:, idim), grad_dl_psi(:, :, idim), set_bc = .false.)
          call X(derivatives_grad)(gr%der, dl_psi2(:, idim), grad_dl_psi2(:, :, idim), set_bc = .false.)

          do idir = 1, gr%mesh%sb%dim
            call X(derivatives_set_bc)(gr%der, grad_dl_psi(:, idir, idim))
            call X(derivatives_set_bc)(gr%der, grad_dl_psi2(:, idir, idim))
          enddo
        endif

        !accumulate to calculate the gradient of the density
        if (present(lr)) then
          ff = st%d%kweights(ik) * st%occ(ist, ik)
          do idir = 1, gr%mesh%sb%dim
            do ip = 1, np
              grad_rho(ip, idir) = grad_rho(ip, idir) + ff * &
                   (R_CONJ(grad_psi(ip, idir, idim)) * dl_psi(ip, idim) + R_CONJ(psi(ip, idim)) * grad_dl_psi(ip, idir, idim) &
                   + R_CONJ(dl_psi2(ip, idim)) * grad_psi(ip, idir, idim) + R_CONJ(grad_dl_psi2(ip, idir, idim)) * psi(ip, idim))
            end do
          end do
        else
          ff = st%d%kweights(ik) * st%occ(ist, ik) * M_TWO
          do idir = 1, gr%mesh%sb%dim
            do ip = 1, np
              grad_rho(ip, idir) = grad_rho(ip, idir) + &
                   ff * R_CONJ(psi(ip, idim)) * grad_psi(ip, idir, idim)
            end do
          end do
        endif
      end do

      call profiling_count_operations(np*st%d%dim*gr%mesh%sb%dim*(2 + R_MUL))
      ! probably this is not accurate anymore

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(projector_is_null(ep%proj(iatom))) cycle
        do idir = 1, gr%mesh%sb%dim

          if(present(lr)) then
            force(idir, iatom) = force(idir, iatom) - st%d%kweights(ik) * st%occ(ist, ik) * &
               (X(psia_project_psib)(ep%proj(iatom), st%d%dim, grad_psi(:, idir, :), dl_psi, ik) &
              + X(psia_project_psib)(ep%proj(iatom), st%d%dim, psi, grad_dl_psi(:, idir, :), ik) &
              + X(psia_project_psib)(ep%proj(iatom), st%d%dim, dl_psi2, grad_psi(:, idir, :), ik) &
              + X(psia_project_psib)(ep%proj(iatom), st%d%dim, grad_dl_psi2(:, idir, :), psi, ik))
          else
            force(idir, iatom) = force(idir, iatom) - M_TWO * st%d%kweights(ik) * st%occ(ist, ik) * &
              X(psia_project_psib)(ep%proj(iatom), st%d%dim, psi, grad_psi(:, idir, :), ik)
          endif

        end do
      end do

    end do
  end do
  
  SAFE_DEALLOCATE_A(psi)
  if(present(lr)) then
    SAFE_DEALLOCATE_A(dl_psi)
    SAFE_DEALLOCATE_A(dl_psi2)
  endif

  SAFE_DEALLOCATE_A(grad_psi)
  if(present(lr)) then
    SAFE_DEALLOCATE_A(grad_dl_psi)
    SAFE_DEALLOCATE_A(grad_dl_psi2)
  endif

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then

    call profiling_in(prof, "FORCES_MPI")

    !reduce the force
    SAFE_ALLOCATE(force_local(1:gr%mesh%sb%dim, 1:geo%natoms))
    force_local = force
    call MPI_Allreduce(force_local, force, gr%mesh%sb%dim*geo%natoms, MPI_CMPLX, MPI_SUM, st%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(force_local)

    !reduce the gradient of the density
    SAFE_ALLOCATE(grad_rho_local(1:np, 1:gr%mesh%sb%dim))
    call lalg_copy(np, gr%mesh%sb%dim, grad_rho, grad_rho_local)
    call MPI_Allreduce(grad_rho_local(1, 1), grad_rho(1, 1), np*gr%mesh%sb%dim, MPI_CMPLX, &
      MPI_SUM, st%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(grad_rho_local)

    call profiling_out(prof)

  end if

  if(st%d%kpt%parallel) then
    
    call profiling_in(prof, "FORCES_MPI")
    
    !reduce the force
    SAFE_ALLOCATE(force_local(1:gr%mesh%sb%dim, 1:geo%natoms))
    force_local = force
    call MPI_Allreduce(force_local(1, 1), force(1, 1), gr%mesh%sb%dim*geo%natoms, MPI_CMPLX, &
      MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(force_local)
    
    !reduce the gradient of the density
    SAFE_ALLOCATE(grad_rho_local(1:np, 1:gr%mesh%sb%dim))
    call lalg_copy(np, gr%mesh%sb%dim, grad_rho, grad_rho_local)
    call MPI_Allreduce(grad_rho_local(1, 1), grad_rho(1, 1), np*gr%mesh%sb%dim, MPI_CMPLX, &
      MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(grad_rho_local)

    call profiling_out(prof)

  end if

#endif
  
  do iatom = 1, geo%natoms
    if(present(lr)) then
      Born_charges%charge(lr_dir, 1:gr%mesh%sb%dim, iatom) = force(1:gr%mesh%sb%dim, iatom)
    else
      geo%atom(iatom)%f(1:gr%mesh%sb%dim) = geo%atom(iatom)%f(1:gr%mesh%sb%dim) + real(force(1:gr%mesh%sb%dim, iatom))
    endif
  end do

  ! THE LOCAL PART (parallel in atoms)

  SAFE_ALLOCATE(vloc(1:np))
  SAFE_ALLOCATE(zvloc(1:np))
  
  do iatom = geo%atoms%start, geo%atoms%end
    
    vloc(1:np) = M_ZERO
    
    call epot_local_potential(ep, gr%der, gr%dgrid, psolver, geo, iatom, vloc, time)

    forall(ip = 1:np) zvloc(ip) = vloc(ip)

    do idir = 1, gr%mesh%sb%dim
      force(idir, iatom) = -zmf_dotp(gr%mesh, zvloc, grad_rho(:, idir))
    end do

  end do

  SAFE_DEALLOCATE_A(vloc)
  SAFE_DEALLOCATE_A(zvloc)
  SAFE_DEALLOCATE_A(grad_rho)

#ifdef HAVE_MPI
  if(geo%atoms%parallel) then

    call profiling_in(prof, "FORCES_MPI")

    ! each node has a piece of the force array, they have to be
    ! collected by all nodes, MPI_Allgatherv does precisely this (if
    ! we get the arguments right).

    SAFE_ALLOCATE(recv_count(1:geo%atoms%mpi_grp%size))
    SAFE_ALLOCATE(recv_displ(1:geo%atoms%mpi_grp%size))
    SAFE_ALLOCATE(force_local(1:gr%mesh%sb%dim, 1:geo%atoms%nlocal))

    recv_count(1:geo%atoms%mpi_grp%size) = gr%mesh%sb%dim*geo%atoms%num(0:geo%atoms%mpi_grp%size - 1)
    recv_displ(1:geo%atoms%mpi_grp%size) = gr%mesh%sb%dim*(geo%atoms%range(1, 0:geo%atoms%mpi_grp%size - 1) - 1)

    force_local(1:gr%mesh%sb%dim, 1:geo%atoms%nlocal) = force(1:gr%mesh%sb%dim, geo%atoms%start:geo%atoms%end)

    call MPI_Allgatherv(&
         force_local(1, 1), gr%mesh%sb%dim*geo%atoms%nlocal, MPI_CMPLX, &
         force(1, 1), recv_count(1), recv_displ(1), MPI_CMPLX, &
         geo%atoms%mpi_grp%comm, mpi_err)

    SAFE_DEALLOCATE_A(recv_count)
    SAFE_DEALLOCATE_A(recv_displ)
    SAFE_DEALLOCATE_A(force_local)

    call profiling_out(prof)

  end if
#endif

  if(present(lr)) then
    do iatom = 1, geo%natoms
      Born_charges%charge(lr_dir, lr_dir, iatom) = Born_charges%charge(lr_dir, lr_dir, iatom) + &
        species_zval(geo%atom(iatom)%spec)
      do idir = 1, gr%mesh%sb%dim
        Born_charges%charge(lr_dir, idir, iatom) = Born_charges%charge(lr_dir, idir, iatom) + force(idir, iatom)
      enddo
    enddo
  else
    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + real(force(idir, iatom))
      end do
    end do
  endif

  SAFE_DEALLOCATE_A(force)
  call pop_sub()
  
end subroutine X(forces_from_potential)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
