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

  integer :: iatom, ist, ik, idim, idir

  R_TYPE :: psi_proj_gpsi
  R_TYPE, allocatable :: gpsi(:, :, :)
  FLOAT,  allocatable :: grho(:, :), vloc(:)
  FLOAT,  allocatable :: force(:, :)
#ifdef HAVE_MPI
  FLOAT, allocatable  :: force_local(:, :)
#endif

  ALLOCATE(gpsi(gr%m%np, 1:NDIM, st%d%dim), gr%m%np*NDIM*st%d%dim)
  ALLOCATE(grho(NP, MAX_DIM), NP*MAX_DIM)

  ALLOCATE(force(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)

  force = M_ZERO

  !$omp parallel workshare
  grho(1:NP, 1:NDIM) = M_ZERO
  !$omp end parallel workshare    

  !the non-local part
  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end

      ! calculate the gradient of the wave-function
      do idim = 1, st%d%dim
        call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)(:, idim, ist, ik), gpsi(:, :, idim))

        !accumulate to calculate the gradient of the density
        do idir = 1, NDIM
          grho(1:NP, idir) = grho(1:NP, idir) + st%d%kweights(ik)*st%occ(ist, ik) * M_TWO * &
               R_REAL(st%X(psi)(1:NP, idim, ist, ik)*R_CONJ(gpsi(1:NP, idir, idim)))
        end do
      end do

      call profiling_count_operations(NP*st%d%dim*NDIM*(2 + R_MUL))

      ! iterate over the projectors
      do iatom = 1, geo%natoms
        if(ep%p(iatom)%type == 0) cycle
        do idir = 1, NDIM
          
          psi_proj_gpsi = &
            X(psia_project_psib)(gr%m, ep%p(iatom), st%d%dim, st%X(psi)(:, :, ist, ik), gpsi(:, idir, :), ep%reltype, ik)
          
          force(idir, iatom) = force(idir, iatom) - M_TWO * st%occ(ist, ik) * R_REAL(psi_proj_gpsi)

        end do
      end do

    end do
  end do

  deallocate(gpsi)

  !now add the local part
  ALLOCATE(vloc(1:NP), NP)
  
  do iatom = 1, geo%natoms
    
    !$omp parallel workshare
    vloc(1:NP) = M_ZERO
    !$omp end parallel workshare
    
    call epot_local_potential(ep, gr, geo, geo%atom(iatom), vloc, time)
    
    do idir = 1, NDIM
      force(idir, iatom) = force(idir, iatom) - dmf_dotp(gr%m, grho(1:NP, idir), vloc(1:NP))
    end do

  end do
  
  deallocate(vloc)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ALLOCATE(force_local(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)
    force_local = force
    call MPI_Allreduce(force_local, force, MAX_DIM*geo%natoms, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    deallocate(force_local)
  end if
#endif
  
  do iatom = 1, geo%natoms
    geo%atom(iatom)%f(1:MAX_DIM) = geo%atom(iatom)%f(1:MAX_DIM) + force(1:MAX_DIM, iatom)
  end do
  
  deallocate(force)
  
end subroutine X(calc_forces_from_potential)

!This function calculates |cpsi> = [x,V_nl] |psi>

subroutine X(conmut_vnl_r)(gr, geo, ep, dim, idir, iatom, psi, cpsi, ik)
  type(grid_t),          intent(in)     :: gr
  type(geometry_t),      intent(in)     :: geo
  type(epot_t), target,  intent(in)     :: ep
  integer,               intent(in)     :: dim
  integer,               intent(in)     :: idir
  integer,               intent(in)     :: iatom
  R_TYPE,                intent(in)     :: psi(:)
  R_TYPE,                intent(out)    :: cpsi(:,:)
  integer,               intent(in)     :: ik

  integer ::  n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), pxlpsi(:,:), xplpsi(:,:)
  integer, pointer :: jxyz(:)
  FLOAT,   pointer :: smx(:, :)

  call push_sub('epot_inc.Xconmut_vnl_r')


  cpsi(1:gr%m%np, 1:dim) = M_ZERO

  n_s = ep%p(iatom)%sphere%ns
  jxyz => ep%p(iatom)%sphere%jxyz
  smx => ep%p(iatom)%sphere%x

  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(xplpsi(n_s, dim), n_s*dim)
  ALLOCATE(pxlpsi(n_s, dim), n_s*dim)

  do idim = 1, dim
    lpsi(1:n_s, idim) = psi(jxyz(1:n_s))
  end do

  ! x V_nl |psi>
  call X(project_sphere)(gr%m, ep%p(iatom), dim, lpsi, xplpsi, ep%reltype)
  do idim = 1, dim
    xplpsi(1:n_s, idim) = smx(1:n_s, idir) * xplpsi(1:n_s, idim)
  end do

  ! V_nl x |psi>
  do idim = 1, dim
    lpsi(1:n_s, idim) = smx(1:n_s, idir) * lpsi(1:n_s, idim)
  end do
  call X(project_sphere)(gr%m, ep%p(iatom), dim, lpsi, pxlpsi, ep%reltype)

  ! |cpsi> = x V_nl |psi> - V_nl x |psi> 
  do idim = 1, dim
    cpsi(jxyz(1:n_s), idim) = cpsi(jxyz(1:n_s), idim) + xplpsi(1:n_s, idim) - pxlpsi(1:n_s, idim)
  end do

  deallocate(lpsi, xplpsi, pxlpsi)

  call pop_sub()

end subroutine X(conmut_vnl_r)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
