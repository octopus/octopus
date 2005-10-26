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

subroutine X(xc_KLI_solve) (m, st, is, oep)
  type(mesh_type),   intent(in)    :: m
  type(states_type), intent(in)    :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE), ierr
#endif
  integer :: i, j, n
  FLOAT, allocatable :: d(:)
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:)
  FLOAT, allocatable :: Ma(:,:), x(:,:), y(:,:)
  R_TYPE, allocatable :: phi1(:,:), phi2(:,:)
  R_TYPE :: occ

  ! some intermediate quantities
  ! vxc contains the Slater part!
  allocate(rho_sigma(m%np))

  do i = 1, m%np
    rho_sigma(i) = max(sum(oep%socc*st%occ(st%st_start:st%st_end, is) * &
       R_ABS(st%X(psi)(i, 1, st%st_start:st%st_end, is))**2), CNST(1e-20))
  end do
#if defined(HAVE_MPI)
  allocate(d(m%np))

  call MPI_Allreduce(rho_sigma(:), d(:), m%np, MPI_FLOAT, MPI_SUM, st%comm, ierr)
  rho_sigma(1:m%np) = d(1:m%np)
#endif

  do i = 1, m%np
    oep%vxc(i)   = oep%socc* sum(st%occ(st%st_start:st%st_end, is) * &
       R_REAL(oep%X(lxc)(i, st%st_start:st%st_end)*st%X(psi)(i, 1, st%st_start:st%st_end, is)))/rho_sigma(i)
  end do
#if defined(HAVE_MPI)
  call MPI_Allreduce(oep%vxc(:),   d(:), m%np, MPI_FLOAT, MPI_SUM, st%comm, ierr)
  oep%vxc(1:m%np)   = d(1:m%np)

  deallocate(d)
#endif

  if(oep%level == XC_OEP_SLATER) then
    deallocate(rho_sigma)
    return
  end if

  n = oep%eigen_n

  allocate(v_bar_S(st%nst))
  do i = st%st_start, st%st_end
    if(st%occ(i, is) .gt. small) then
      v_bar_S(i) = sum(R_ABS(st%X(psi)(1:m%np, 1, i, is))**2 * oep%vxc(1:m%np) * m%vol_pp(1:m%np))
    end if
  end do

#if defined(HAVE_MPI)
  ! Broadcast the vector v_bar_S to all processors
  if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
    do i = 1, st%nst
      call MPI_Bcast(v_bar_S(i), 1, MPI_FLOAT, st%node(i), st%comm, ierr)
    end do
  end if
#endif

  ! If there is more than one state, so solve linear equation.
  linear_equation: if(n > 0) then
    allocate(phi1(m%np, st%d%dim), phi2(m%np, st%d%dim))
    phi1 = M_ZERO; phi2 = M_ZERO

    allocate(d(m%np), x(n, 1))
    x = M_ZERO

    allocate(Ma(n, n), y(n, 1))

    do i = 1, n

#if defined(HAVE_MPI)
      if(mpiv%node.ne.0) then
        if(st%node(oep%eigen_index(i)) == mpiv%node) then
          call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
             0, i, st%comm, status, ierr)
        end if
      else ! master node
        if(st%node(oep%eigen_index(i)).ne.0) then
          call MPI_Recv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
             i, st%comm, status, ierr)
        else
          phi1(:,:) = st%X(psi)(:,:, oep%eigen_index(i), is)
        end if
      end if
      call MPI_Barrier(st%comm, ierr)
#else
      phi1(:,:) = st%X(psi)(:,:, oep%eigen_index(i), is)
#endif

      if(mpiv%node == 0) then
        d(1:m%np) = (R_REAL(phi1(1:m%np, 1))**2 + &
           R_AIMAG(phi1(1:m%np, 1))**2) / rho_sigma(1:m%np) * m%vol_pp(1:m%np)
      end if

      do j = i, n
#if defined(HAVE_MPI)
        if(mpiv%node.ne.0) then
          if(st%node(oep%eigen_index(j)) == mpiv%node) then
            call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(j), is), st%d%dim*m%np, R_MPITYPE, &
               0, j, st%comm, status, ierr)
          end if
        else
          if(st%node(oep%eigen_index(j)).ne.0) then
            call MPI_Recv(phi2(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(j)), &
               j, st%comm, status, ierr)
          else
            phi2(:,:) = st%X(psi)(:,:, oep%eigen_index(j), is)
          end if
        end if
        call MPI_Barrier(st%comm, ierr)
#else
        phi2(:,:) = st%X(psi)(:,:, oep%eigen_index(j), is)
#endif

        if(mpiv%node == 0) then
          Ma(i, j) = - sum(d(1:m%np) * (R_REAL(phi2(1:m%np, 1))**2 + R_AIMAG(phi2(1:m%np, 1))**2) )
          Ma(j,i) = Ma(i,j)
        end if
      end do

      if(mpiv%node == 0) then
        Ma(i,i) = M_ONE + Ma(i,i)
        y(i, 1) = v_bar_S(oep%eigen_index(i)) - oep%uxc_bar(oep%eigen_index(i))
      end if
    end do

    if(mpiv%node == 0) then
      call lalg_linsyssolve(n, 1, Ma, y, x)
    end if
    deallocate(Ma, y)

#if defined(HAVE_MPI)
    call MPI_Barrier(st%comm, ierr)
#endif

    do i = 1, n
      ! add contribution of low lying states. Once again, we have to send every state
      ! that is not in node zero to node 0
#if defined(HAVE_MPI)
      if(mpiv%node.ne.0) then
        if(st%node(oep%eigen_index(i)) == mpiv%node) then
          call MPI_Send(st%X(psi)(1, 1, oep%eigen_index(i), is), st%d%dim*m%np, R_MPITYPE, &
             0, i, st%comm, status, ierr)
        end if
      else ! master
        if(st%node(oep%eigen_index(i)).ne.0) then
          call MPI_Recv(phi1(1, 1), st%d%dim*m%np, R_MPITYPE, st%node(oep%eigen_index(i)), &
             i, st%comm, status, ierr)
        else
          phi1(:,:) = st%X(psi)(:,:, oep%eigen_index(i), is)
        end if
      end if
      call MPI_Barrier(st%comm, ierr)
#else
      phi1(:,:) = st%X(psi)(:,:, oep%eigen_index(i), is)
#endif
      occ = st%occ(oep%eigen_index(i), is)

      if(mpiv%node == 0) then
        oep%vxc(1:m%np) = oep%vxc(1:m%np) + &
           oep%socc * occ * x(i,1) * R_ABS(phi1(1:m%np, 1))**2 / rho_sigma(1:m%np)
      end if
    end do

    deallocate(d, x, phi1, phi2)
  end if linear_equation
  ! The previous stuff is only needed if n>0.

#if defined(HAVE_MPI)
  ! Since only node 0 holds the true vxc, broadcast it to all processors.
  call MPI_Bcast(oep%vxc(:), m%np, MPI_FLOAT, 0, st%comm, ierr)
#endif

  deallocate(v_bar_S, rho_sigma)
end subroutine X(xc_KLI_solve)
