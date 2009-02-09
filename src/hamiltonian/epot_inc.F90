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

subroutine X(calc_forces_from_potential)(gr, geo, ep, st, time, lr, lr_dir)
  type(grid_t),     intent(inout)  :: gr
  type(geometry_t), intent(inout)  :: geo
  type(epot_t),     intent(in)     :: ep
  type(states_t),   intent(inout)  :: st
  FLOAT,            intent(in)     :: time
  type(lr_t), optional, intent(inout) :: lr
  integer,    optional, intent(in) :: lr_dir    
  ! provide these optional arguments to calculate Born effective charges rather than forces
  ! lr should be the wfns from electric perturbation in the lr_dir direction
  ! for each atom, Z*(i,j) = dF(j)/dE(i)

  integer :: iatom, ist, ik, idim, idir, np, ip, factor

  R_TYPE, pointer     :: psi(:, :)
  R_TYPE, pointer     :: dl_psi(:, :)
  R_TYPE, allocatable :: grad_psi(:, :, :)
  R_TYPE, allocatable :: grad_dl_psi(:, :, :)
  FLOAT,  allocatable :: vloc(:)
  CMPLX,  allocatable :: grad_rho(:, :), force(:, :)
#ifdef HAVE_MPI
  integer, allocatable :: recv_count(:), recv_displ(:)
  CMPLX, allocatable  :: force_local(:, :), grad_rho_local(:, :)
  type(profile_t), save :: prof
#endif

  call push_sub('epot_inc.Xcalc_forces_from_potential')

  ASSERT(present(lr) .eqv. present(lr_dir))
  ! need both to calculate Born charges

  np = gr%fine%m%np

  if(present(lr)) then
    ALLOCATE(grad_dl_psi(np, 1:NDIM, st%d%dim), np*NDIM*st%d%dim)
  endif
  ALLOCATE(grad_psi(np, 1:NDIM, st%d%dim), np*NDIM*st%d%dim)
  ALLOCATE(grad_rho(np, MAX_DIM), np*MAX_DIM)
  grad_rho(1:np, 1:NDIM) = M_ZERO
  ALLOCATE(force(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
  force = M_ZERO

  if(gr%have_fine_mesh) then
     ALLOCATE(psi(gr%fine%m%np_part, st%d%dim), gr%fine%m%np_part*st%d%dim)
     if(present(lr)) then
       ALLOCATE(dl_psi(gr%fine%m%np_part, st%d%dim), gr%fine%m%np_part*st%d%dim)
     endif
  endif

  !THE NON-LOCAL PART (parallel in states and k-points)
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      
      if(gr%have_fine_mesh) then
        do idim = 1, st%d%dim
          call X(multigrid_coarse2fine)(gr%fine, st%X(psi)(:, idim, ist, ik), psi(:, idim))

          if (present(lr)) then
            call X(multigrid_coarse2fine)(gr%fine, lr%X(dl_psi)(:, idim, ist, ik), dl_psi(:, idim))
          endif
        end do
      else
        psi => st%X(psi)(:, :, ist, ik)

        if (present(lr)) then
          dl_psi => lr%X(dl_psi)(:, :, ist, ik)
        endif
      end if

      ! calculate the gradients of the wave-functions
      do idim = 1, st%d%dim

        call X(derivatives_grad)(gr%fine%der, psi(:, idim), grad_psi(:, :, idim))

        if (present(lr)) then
          call X(derivatives_grad)(gr%fine%der, dl_psi(:, idim), grad_dl_psi(:, :, idim))
        endif

        !accumulate to calculate the gradient of the density
        ! is -sigma needed here too? is Z* linear in density? no, yes
        if (present(lr)) then
          forall (idir = 1:NDIM, ip = 1:np) &
            grad_rho(ip, idir) = grad_rho(ip, idir) + st%d%kweights(ik) * st%occ(ist, ik) * M_TWO * &
            (R_CONJ(dl_psi(ip, idim)) * grad_psi(ip, idir, idim) + R_CONJ(psi(ip, idim)) * grad_dl_psi(ip, idir, idim))
        else
          forall (idir = 1:NDIM, ip = 1:np) grad_rho(ip, idir) = grad_rho(ip, idir) + &
               st%d%kweights(ik)*st%occ(ist, ik)*M_TWO*R_CONJ(psi(ip, idim))*grad_psi(ip, idir, idim)
        endif
      end do

      call profiling_count_operations(np*st%d%dim*NDIM*(2 + R_MUL))
      ! probably this is not accurate anymore

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(ep%proj(iatom)%type == 0) cycle
        do idir = 1, NDIM

          if(present(lr)) then
            force(idir, iatom) = force(idir, iatom) - M_TWO * st%occ(ist, ik) * &
              (X(psia_project_psib)(ep%proj_fine(iatom), st%d%dim, dl_psi, grad_psi(:, idir, :), ik) &
              + X(psia_project_psib)(ep%proj_fine(iatom), st%d%dim, psi, grad_dl_psi(:, idir, :), ik))
          else
            force(idir, iatom) = force(idir, iatom) - M_TWO * st%occ(ist, ik) * &
              X(psia_project_psib)(ep%proj_fine(iatom), st%d%dim, psi, grad_psi(:, idir, :), ik)
          endif

        end do
      end do

    end do
  end do

  if(gr%have_fine_mesh) then
    deallocate(psi)
    if(present(lr)) then
      deallocate(dl_psi)
    endif
  endif

  deallocate(grad_psi)
  if(present(lr)) then
    deallocate(grad_dl_psi)
  endif

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then

    call profiling_in(prof, "FORCES_MPI")

    !reduce the force
    ALLOCATE(force_local(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
    force_local = force
    call MPI_Allreduce(force_local, force, MAX_DIM*geo%natoms, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    deallocate(force_local)

    !reduce the gradient of the density
    ALLOCATE(grad_rho_local(np, MAX_DIM), np*MAX_DIM)
    call lalg_copy(np, MAX_DIM, grad_rho, grad_rho_local)
    call MPI_Allreduce(grad_rho_local, grad_rho, np*MAX_DIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    deallocate(grad_rho_local)

    call profiling_out(prof)

  end if

  if(st%d%kpt%parallel) then
    
    call profiling_in(prof, "FORCES_MPI")
    
    !reduce the force
    ALLOCATE(force_local(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
    force_local = force
    call MPI_Allreduce(force_local, force, MAX_DIM*geo%natoms, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
    deallocate(force_local)
    
    !reduce the gradient of the density
    ALLOCATE(grad_rho_local(np, MAX_DIM), np*MAX_DIM)
    call lalg_copy(np, MAX_DIM, grad_rho, grad_rho_local)
    call MPI_Allreduce(grad_rho_local, grad_rho, np*MAX_DIM, MPI_FLOAT, MPI_SUM, st%d%kpt%mpi_grp%comm, mpi_err)
    deallocate(grad_rho_local)

    call profiling_out(prof)

  end if

#endif
  
  do iatom = 1, geo%natoms
    if(present(lr)) then
      geo%atom(iatom)%Born_charge(lr_dir, 1:MAX_DIM) = force(1:MAX_DIM, iatom)
      ! the force has other contributions from Hartree and Vxc, but the Born charge does not 
    else
      geo%atom(iatom)%f(1:MAX_DIM) = geo%atom(iatom)%f(1:MAX_DIM) + real(force(1:MAX_DIM, iatom))
    endif
  end do

  write(*,*) 'done with non-local part of force'

  ! THE LOCAL PART (parallel in atoms)

  ALLOCATE(vloc(1:np), np)
  
  do iatom = geo%atoms%start, geo%atoms%end
    
    vloc(1:np) = M_ZERO
    
    call epot_local_potential(ep, gr, gr%fine%m, geo, iatom, vloc, time)
    
    do idir = 1, NDIM
      force(idir, iatom) = -dmf_dotp(gr%fine%m, real(grad_rho(:, idir)), vloc) &
                           - M_zI * dmf_dotp(gr%fine%m, aimag(grad_rho(:, idir)), vloc)
    end do

  end do

  deallocate(vloc)

#ifdef HAVE_MPI
  if(geo%atoms%parallel) then

    call profiling_in(prof, "FORCES_MPI")

    ! each node has a piece of the force array, they have to be
    ! collected by all nodes, MPI_Allgatherv does precisely this (if
    ! we get the arguments right).

    ALLOCATE(recv_count(geo%atoms%mpi_grp%size), geo%atoms%mpi_grp%size)
    ALLOCATE(recv_displ(geo%atoms%mpi_grp%size), geo%atoms%mpi_grp%size)
    ALLOCATE(force_local(1:MAX_DIM, 1:geo%atoms%nlocal), MAX_DIM*geo%atoms%nlocal)

    recv_count(1:geo%atoms%mpi_grp%size) = MAX_DIM*geo%atoms%num(0:geo%atoms%mpi_grp%size - 1)
    recv_displ(1:geo%atoms%mpi_grp%size) = MAX_DIM*(geo%atoms%range(1, 0:geo%atoms%mpi_grp%size - 1) - 1)

    force_local(1:MAX_DIM, 1:geo%atoms%nlocal) = force(1:MAX_DIM, geo%atoms%start:geo%atoms%end)

    call MPI_Allgatherv(&
         force_local, MAX_DIM*geo%atoms%nlocal, MPI_FLOAT, &
         force, recv_count(1), recv_displ, MPI_FLOAT, &
         geo%atoms%mpi_grp%comm, mpi_err)

    deallocate(recv_count, recv_displ, force_local)

    call profiling_out(prof)

  end if
#endif

  if(present(lr)) then
    do iatom = 1, geo%natoms
      geo%atom(iatom)%Born_charge(lr_dir, lr_dir) = geo%atom(iatom)%Born_charge(lr_dir, lr_dir) + geo%atom(iatom)%spec%Z_val
      do idir = 1, MAX_DIM
        geo%atom(iatom)%Born_charge(lr_dir, idir) = geo%atom(iatom)%Born_charge(lr_dir, idir) + force(idir, iatom)
!        write(*,'(a,i4,a,f4.1,a,i4,a,i4,a,)') 'iatom = ', iatom, ', Z = ', geo%atom(iatom)%spec%Z_val, ', lr_dir = ', lr_dir, ', idir = ', idir
!        write(*,'(a,f10.6,a,f10.6)') 'Re Born = ', real(geo%atom(iatom)%Born_charge(lr_dir, idir)), ', Im Born = ', aimag(geo%atom(iatom)%Born_charge(lr_dir, idir))
      enddo
    enddo
  else
    forall (iatom = 1:geo%natoms, idir = 1:MAX_DIM) geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + real(force(idir, iatom))
  endif

  deallocate(force)
  call pop_sub()
  
end subroutine X(calc_forces_from_potential)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
