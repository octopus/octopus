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
!!
!! $Id$

! ---------------------------------------------------------
subroutine X(xc_KLI_solve) (m, st, is, oep)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  integer,        intent(in)    :: is
  type(xc_oep_t), intent(inout) :: oep

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE), ierr
#endif
  integer :: i, j, n
  FLOAT, allocatable :: d(:)
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:)
  FLOAT, allocatable :: Ma(:,:), x(:,:), y(:,:)
  R_TYPE, allocatable :: phi1(:,:), phi2(:,:)
  R_TYPE :: occ

  call push_sub('xc_KLI.xc_KLI_solve')

  ! some intermediate quantities
  ! vxc contains the Slater part!
  ALLOCATE(rho_sigma(m%np), m%np)

  do i = 1, m%np
    rho_sigma(i) = max(sum(oep%socc*st%occ(st%st_start:st%st_end, is) * &
      R_ABS(st%X(psi)(i, 1, st%st_start:st%st_end, is))**2), CNST(1e-20))
  end do
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ALLOCATE(d(m%np), m%np)
    call MPI_Allreduce(rho_sigma(:), d(:), m%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, ierr)
    rho_sigma(1:m%np) = d(1:m%np)
  end if
#endif

  do i = 1, m%np
    oep%vxc(i)   = oep%socc* sum(st%occ(st%st_start:st%st_end, is) * &
      R_REAL(oep%X(lxc)(i, st%st_start:st%st_end)*st%X(psi)(i, 1, st%st_start:st%st_end, is)))/rho_sigma(i)
  end do
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(oep%vxc(:),   d(:), m%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, ierr)
    oep%vxc(1:m%np)   = d(1:m%np)
    deallocate(d)
  end if
#endif

  if(oep%level == XC_OEP_SLATER) then
    deallocate(rho_sigma)
    return
  end if

  n = oep%eigen_n

  ALLOCATE(v_bar_S(st%nst), st%nst)
  do i = st%st_start, st%st_end
    if(st%occ(i, is) .gt. small) then
      v_bar_S(i) = dmf_dotp(m, (R_ABS(st%X(psi)(1:m%np, 1, i, is)))**2 , oep%vxc)
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Broadcast the vector v_bar_S to all processors
    if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
      do i = 1, st%nst
        call MPI_Bcast(v_bar_S(i), 1, MPI_FLOAT, st%node(i), st%mpi_grp%comm, ierr)
      end do
    end if
  end if
#endif

  ! If there is more than one state, so solve linear equation.
  linear_equation: if(n > 0) then
    ALLOCATE(phi1(m%np, st%d%dim), m%np*st%d%dim)
    ALLOCATE(phi2(m%np, st%d%dim), m%np*st%d%dim)
    phi1 = M_ZERO; phi2 = M_ZERO

    ALLOCATE(d(m%np), m%np)
    ALLOCATE(x(n, 1), n*1)
    x = M_ZERO

    ALLOCATE(Ma(n, n), n*n)
    ALLOCATE(y(n, 1), n*1)

    do i = 1, n

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        if(.not.mpi_grp_is_root(st%mpi_grp)) then
          if(st%node(oep%eigen_index(i)) == st%mpi_grp%rank) then
            call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
              0, i, st%mpi_grp%comm, ierr)
          end if
        else ! master node
          if(st%node(oep%eigen_index(i)).ne.0) then
            call MPI_Recv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
              i, st%mpi_grp%comm, status, ierr)
          else
            phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
          end if
        end if
        call MPI_Barrier(st%mpi_grp%comm, ierr)
      else
        phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
      endif
#else
      phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
#endif

      if(mpi_grp_is_root(st%mpi_grp)) then
        d(1:m%np) = (R_REAL(phi1(1:m%np, 1))**2 + &
          R_AIMAG(phi1(1:m%np, 1))**2) / rho_sigma(1:m%np)
      end if

      do j = i, n
#if defined(HAVE_MPI)
        if(st%parallel_in_states) then
          if(.not.mpi_grp_is_root(st%mpi_grp)) then
            if(st%node(oep%eigen_index(j)) == st%mpi_grp%rank) then
              call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(j), is), st%d%dim*m%np, R_MPITYPE, &
                0, j, st%mpi_grp%comm, ierr)
            end if
          else
            if(st%node(oep%eigen_index(j)).ne.0) then
              call MPI_Recv(phi2(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(j)), &
                j, st%mpi_grp%comm, status, ierr)
            else
              phi2(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(j), is)
            end if
          end if
          call MPI_Barrier(st%mpi_grp%comm, ierr)
        else
          phi2(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(j), is)
        end if
#else
        phi2(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(j), is)
#endif

        if(mpi_grp_is_root(st%mpi_grp)) then
          Ma(i, j) = - dmf_dotp(m, d, (R_REAL(phi2(1:m%np, 1))**2 + R_AIMAG(phi2(1:m%np, 1))**2) )
          Ma(j,i) = Ma(i,j)
        end if
      end do

      if(mpi_grp_is_root(st%mpi_grp)) then
        Ma(i,i) = M_ONE + Ma(i,i)
        y(i, 1) = v_bar_S(oep%eigen_index(i)) - oep%uxc_bar(oep%eigen_index(i))
      end if
    end do

    if(mpi_grp_is_root(st%mpi_grp)) then
      call lalg_linsyssolve(n, 1, Ma, y, x)
    end if
    deallocate(Ma, y)

#if defined(HAVE_MPI)
    if(st%parallel_in_states)  call MPI_Barrier(st%mpi_grp%comm, ierr)
#endif

    do i = 1, n
      ! add contribution of low lying states. Once again, we have to send every state
      ! that is not in node zero to node 0
#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        if(.not.mpi_grp_is_root(st%mpi_grp)) then
          if(st%node(oep%eigen_index(i)) == st%mpi_grp%rank) then
            call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
              0, i, st%mpi_grp%comm, ierr)
          end if
        else ! master
          if(st%node(oep%eigen_index(i)).ne.0) then
            call MPI_Recv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
              i, st%mpi_grp%comm, status, ierr)
          else
            phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
          end if
        end if
        call MPI_Barrier(st%mpi_grp%comm, ierr)
      else
        phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
      end if
#else
      phi1(1:m%np, 1:st%d%dim) = st%X(psi)(1:m%np, 1:st%d%dim, oep%eigen_index(i), is)
#endif
      occ = st%occ(oep%eigen_index(i), is)

      if(mpi_grp_is_root(st%mpi_grp)) then
        oep%vxc(1:m%np) = oep%vxc(1:m%np) + &
          oep%socc * occ * x(i,1) * R_ABS(phi1(1:m%np, 1))**2 / rho_sigma(1:m%np)
      end if
    end do

    deallocate(d, x, phi1, phi2)
  end if linear_equation
  ! The previous stuff is only needed if n>0.

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Since only node 0 holds the true vxc, broadcast it to all processors.
    call MPI_Bcast(oep%vxc(:), m%np, MPI_FLOAT, 0, st%mpi_grp%comm, ierr)
  end if
#endif

  deallocate(v_bar_S, rho_sigma)
  call pop_sub()
end subroutine X(xc_KLI_solve)
