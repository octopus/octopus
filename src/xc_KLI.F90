!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

subroutine X(xc_KLI_solve) (m, st, is, oep, oep_level)
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(in)    :: st
  type(xc_oep_type), intent(inout) :: oep
  integer,           intent(in)    :: oep_level, is

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE)
#endif
  integer :: i, j, n, ierr, request
  FLOAT, allocatable :: d(:)
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:)
  FLOAT, allocatable :: Ma(:,:), x(:,:), y(:,:)
  R_TYPE, allocatable :: phi1(:, :), phi2(:, :)
  R_TYPE :: occ

  ! some intermediate quantities
  ! vxc contains the Slater part!
  allocate(rho_sigma(m%np))
#if defined(HAVE_MPI)
  allocate(d(m%np))
  do i = 1, m%np
    d(i) = max(sum(oep%socc*st%occ(:, is)*R_ABS(st%X(psi)(i, 1, :, is))**2), CNST(1e-20))
  enddo
  if(st%st_end - st%st_start + 1 < st%nst) then
    call mpi_allreduce(d(1), rho_sigma(1), m%np, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    rho_sigma = d
  endif
  do i = 1, m%np
     d(i)   = oep%socc * sum(st%occ(:, is)*R_REAL(oep%lxc(i, :)*st%X(psi)(i, 1, :, is)))/rho_sigma(i)
  enddo
  if(st%st_end - st%st_start + 1 < st%nst) then
    call mpi_allreduce(d(1), oep%vxc(1), m%np, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
  else
    oep%vxc = d
  endif
  deallocate(d)
#else
  do i = 1, m%np
    rho_sigma(i) = max(sum(oep%socc*st%occ(:, is)*R_ABS(st%X(psi)(i, 1, :, is))**2), CNST(1e-20))
    oep%vxc(i)   = oep%socc* &
       sum(st%occ(:, is)*R_REAL(oep%lxc(i, :)*st%X(psi)(i, 1, :, is)))/rho_sigma(i)
  enddo
#endif

  if(oep_level == XC_OEP_SLATER) then
    deallocate(rho_sigma)
    return
  end if

  n = oep%eigen_n

  allocate(v_bar_S(st%nst))
  do i = st%st_start, st%st_end
    if(st%occ(i, is) .gt. small) then
      v_bar_S(i) = sum(R_ABS(st%X(psi)(:, 1, i, is))**2 * oep%vxc(:) * m%vol_pp(:))
    end if
  end do
#if defined(HAVE_MPI)
  ! Broadcast the vector v_bar_S to all processors
  if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      do i = 1, st%nst
         call mpi_bcast(v_bar_S(i), 1, MPI_FLOAT, st%node(i), MPI_COMM_WORLD, ierr)
      enddo
  endif
#endif

  allocate(phi1(m%np, st%d%dim), phi2(m%np, st%d%dim))
  phi1 = M_ZERO; phi2 = M_ZERO


  ! If there is more than one state, so solve linear equation.
  linear_equation: if(n > 0) then 


    allocate(x(n, 1))
    x = M_ZERO
    allocate(Ma(n, n), y(n, 1))
    do i = 1, n
#if defined(HAVE_MPI)
     if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
       if(st%node(oep%eigen_index(i)).ne.0) then
         call mpi_isend(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
                        0, i, MPI_COMM_WORLD, request, ierr)
       endif
       if(mpiv%node.eq.0) then
         if(st%node(oep%eigen_index(i)).ne.0) then
            call mpi_irecv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
                           i, MPI_COMM_WORLD, request, ierr)
            call mpi_wait(request, status, ierr)
         else
            phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
         endif
       endif
       call mpi_barrier(MPI_COMM_WORLD, ierr)
     else
       phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
     endif
#else
       phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
#endif
      do j = i, n
#if defined(HAVE_MPI)
       if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
        if(st%node(oep%eigen_index(j)).ne.0) then
           call mpi_isend(st%X(psi)(1, 1, oep%eigen_index(j), is), st%d%dim*m%np, R_MPITYPE, &
                          0, j, MPI_COMM_WORLD, request, ierr)
        endif
        if(mpiv%node.eq.0) then
           if(st%node(oep%eigen_index(j)).ne.0) then
              call mpi_irecv(phi2(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(j)), &
                             j, MPI_COMM_WORLD, request, ierr)
              call mpi_wait(request, status, ierr)
           else
              phi2(:, :) = st%X(psi)(:, :, oep%eigen_index(j), is)
           endif
        endif
        call mpi_barrier(MPI_COMM_WORLD, ierr)
       else
        phi2(:, :) = st%X(psi)(:, :, oep%eigen_index(j), is)
       endif
#else
        phi2(:, :) = st%X(psi)(:, :, oep%eigen_index(j), is)
#endif
        Ma(i, j) = - sum(R_ABS(phi1(:, 1))**2 * R_ABS(phi2(:, 1))**2 / rho_sigma(:) * m%vol_pp(:))
        Ma(j,i) = Ma(i,j)
      end do
      Ma(i,i) = 1 + Ma(i,i)
      
      y(i, 1) = v_bar_S(oep%eigen_index(i)) - oep%uxc_bar(oep%eigen_index(i))
    end do

#if defined(HAVE_MPI)     
    ! Only node zero will hold the good x
    if(mpiv%node == 0) then
#endif
    call lalg_linsyssolve(n, 1, Ma, y, x)
#if defined(HAVE_MPI)
    endif
    call mpi_barrier(mpi_comm_world, ierr)
#endif
    deallocate(Ma, y)

    do i = 1, n
       ! add contribution of low lying states. Once again, we have to send every state
       ! that is not in node zero to node 0
#if defined(HAVE_MPI)
       if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
        if(st%node(oep%eigen_index(i)).ne.0) then
          call mpi_isend(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
                         0, i, MPI_COMM_WORLD, request, ierr)
        endif
        if(mpiv%node.eq.0) then
           if(st%node(oep%eigen_index(i)).ne.0) then
              call mpi_irecv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
                             i, MPI_COMM_WORLD, request, ierr)
              call mpi_wait(request, status, ierr)
           else
              phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
           endif
        endif
        call mpi_barrier(MPI_COMM_WORLD, ierr)

        if(st%node(oep%eigen_index(i)).ne.0) then
           call mpi_isend(st%occ(oep%eigen_index(i), is), 1, MPI_FLOAT, &
                          0, i, MPI_COMM_WORLD, request, ierr)
        endif
        if(mpiv%node.eq.0) then
           if(st%node(oep%eigen_index(i)).ne.0) then
              call mpi_irecv(occ, 1, MPI_FLOAT, st%node(oep%eigen_index(i)), &
                             i, MPI_COMM_WORLD, status, ierr)
              call mpi_wait(request, status, ierr)
           else
              occ = st%occ(oep%eigen_index(i), is)
           endif
        endif
        call mpi_barrier(MPI_COMM_WORLD, ierr)

       else
        phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
        occ = st%occ(oep%eigen_index(i), is)
       endif
#else
        phi1(:, :) = st%X(psi)(:, :, oep%eigen_index(i), is)
        occ = st%occ(oep%eigen_index(i), is)
#endif
      oep%vxc(:) = oep%vxc(:) + oep%socc * occ * x(i,1) * R_ABS(phi1(:, 1))**2 / rho_sigma(:)
    end do
    deallocate(x)


  end if linear_equation
  ! The previous stuff is only needed if n>0.


#if defined(HAVE_MPI)
  ! Since only node 0 holds the true vxc, broadcast it to all processors.
  call mpi_barrier(mpi_comm_world, ierr)
  if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      call mpi_bcast(oep%vxc(1), m%np, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
  endif
  call mpi_barrier(mpi_comm_world, ierr)
#endif

    
  deallocate(v_bar_S, rho_sigma)
end subroutine X(xc_KLI_solve)
