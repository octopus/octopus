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

subroutine X(calc_forces_from_potential)(gr, geo, ep, st, time)
  type(grid_t),     intent(inout)  :: gr
  type(geometry_t), intent(inout)  :: geo
  type(epot_t),     intent(in)     :: ep
  type(states_t),   intent(inout)  :: st
  FLOAT,            intent(in)     :: time

  integer :: iatom, ist, ik, idim, idir, np

  R_TYPE :: psi_proj_gpsi
  R_TYPE, allocatable :: gpsi(:, :, :)
  R_TYPE, pointer     :: psi(:, :)
  FLOAT,  allocatable :: grho(:, :), vloc(:)
  FLOAT,  allocatable :: force(:, :)
#ifdef HAVE_MPI
  integer, allocatable :: recv_count(:), recv_displ(:)
  FLOAT, allocatable  :: force_local(:, :), grho_local(:, :)
  type(profile_t), save :: prof
#endif

  np = gr%fine%m%np

  ALLOCATE(gpsi(np, 1:NDIM, st%d%dim), np*NDIM*st%d%dim)
  ALLOCATE(grho(np, MAX_DIM), np*MAX_DIM)

  ALLOCATE(force(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)

  if(gr%have_fine_mesh) ALLOCATE(psi(gr%fine%m%np_part, st%d%dim), gr%fine%m%np_part*st%d%dim)

  force = M_ZERO

  !$omp parallel workshare
  grho(1:np, 1:NDIM) = M_ZERO
  !$omp end parallel workshare    

  !THE NON-LOCAL PART (parallel in states)
  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      
      if(gr%have_fine_mesh) then
        do idim = 1, st%d%dim
          call X(multigrid_coarse2fine)(gr%fine, st%X(psi)(:, idim, ist, ik), psi(:, idim))
        end do
      else
        psi => st%X(psi)(:, :, ist, ik)
      end if

      ! calculate the gradient of the wave-function
      do idim = 1, st%d%dim

        call X(derivatives_grad)(gr%fine%f_der%der_discr, psi(:, idim), gpsi(:, :, idim))

        !accumulate to calculate the gradient of the density
        do idir = 1, NDIM
          grho(1:np, idir) = grho(1:np, idir) + &
               st%d%kweights(ik)*st%occ(ist, ik)*M_TWO*R_REAL(psi(1:np, idim)*R_CONJ(gpsi(1:np, idir, idim)))
        end do
      end do

      call profiling_count_operations(np*st%d%dim*NDIM*(2 + R_MUL))

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(ep%proj(iatom)%type == 0) cycle
        do idir = 1, NDIM
          
          psi_proj_gpsi = X(psia_project_psib)(ep%proj_fine(iatom), st%d%dim, psi, gpsi(:, idir, :), ik)
          
          force(idir, iatom) = force(idir, iatom) - M_TWO*st%occ(ist, ik)*R_REAL(psi_proj_gpsi)

        end do
      end do

    end do
  end do

  if(gr%have_fine_mesh) deallocate(psi)
  deallocate(gpsi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then

    call profiling_in(prof, "FORCES_MPI")

    !reduce the force
    ALLOCATE(force_local(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
    force_local = force
    call MPI_Allreduce(force_local, force, MAX_DIM*geo%natoms, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    deallocate(force_local)

    !reduce the gradient of the density
    ALLOCATE(grho_local(np, MAX_DIM), np*MAX_DIM)
    call lalg_copy(np, MAX_DIM, grho, grho_local)
    call MPI_Allreduce(grho_local, grho, np*MAX_DIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    deallocate(grho_local)

    call profiling_out(prof)

  end if
#endif
  
  do iatom = 1, geo%natoms
    geo%atom(iatom)%f(1:MAX_DIM) = geo%atom(iatom)%f(1:MAX_DIM) + force(1:MAX_DIM, iatom)
  end do

  ! THE LOCAL PART (parallel in atoms)

  ALLOCATE(vloc(1:np), np)
  
  do iatom = geo%atoms%start, geo%atoms%end
    
    vloc(1:np) = M_ZERO
    
    call epot_local_potential(ep, gr, gr%fine%m, geo, geo%atom(iatom), vloc, time)
    
    do idir = 1, NDIM
      force(idir, iatom) = -dmf_dotp(gr%fine%m, grho(:, idir), vloc)
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

  do iatom = 1, geo%natoms
    geo%atom(iatom)%f(1:MAX_DIM) = geo%atom(iatom)%f(1:MAX_DIM) + force(1:MAX_DIM, iatom)
  end do
  
  deallocate(force)
  
end subroutine X(calc_forces_from_potential)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
