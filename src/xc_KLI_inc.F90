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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine X(xc_KLI_solve) (m, st, is, oep)
  type(mesh_t),   intent(in)    :: m
  type(states_t), intent(in)    :: st
  integer,        intent(in)    :: is
  type(xc_oep_t), intent(inout) :: oep

  integer :: i, j, n, kssi, kssj, proc
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:), sqphi(:, :, :), d(:)
  FLOAT, allocatable :: Ma(:,:), x(:,:), y(:,:)
  R_TYPE :: occ

  call profiling_in(C_PROFILING_XC_KLI)
  call push_sub('xc_KLI.xc_KLI_solve')

  ! some intermediate quantities
  ! vxc contains the Slater part!
  ALLOCATE(rho_sigma(m%np), m%np)
  ALLOCATE(sqphi(m%np, st%d%dim, st%nst), m%np*st%d%dim*st%nst)

  do i = st%st_start, st%st_end
    sqphi(1:m%np, 1:st%d%dim, i) = R_REAL(st%X(psi)(1:m%np, 1:st%d%dim, i, is))**2 + &
                                   R_AIMAG(st%X(psi)(1:m%np, 1:st%d%dim, i, is))**2
  end do

  do i = 1, m%np
    rho_sigma(i) = max(sum(oep%socc*st%occ(st%st_start:st%st_end, is) * &
      sqphi(i, 1, st%st_start:st%st_end)), CNST(1e-20))
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ALLOCATE(d(m%np), m%np)
    call MPI_Allreduce(rho_sigma(1), d(1), m%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    rho_sigma(1:m%np) = d(1:m%np)
  end if
#endif

  do i = 1, m%np
    oep%vxc(i) = M_ZERO
    do j = st%st_start, st%st_end
      oep%vxc(i) = oep%vxc(i) + oep%socc * st%occ(j, is) * oep%X(lxc)(i, j) * st%X(psi)(i, 1, j, is)
    end do
    oep%vxc(i) = oep%vxc(i) / rho_sigma(i)
  end do
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(oep%vxc(1),   d(1), m%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    oep%vxc(1:m%np)   = d(1:m%np)
    deallocate(d)
  end if
#endif

  if(oep%level == XC_OEP_SLATER) then
    deallocate(rho_sigma, sqphi)
    call pop_sub()
    call profiling_out(C_PROFILING_XC_KLI)
    return
  end if

  n = oep%eigen_n

  ALLOCATE(v_bar_S(st%nst), st%nst)
  do i = st%st_start, st%st_end
    if(st%occ(i, is) .gt. small) then
      v_bar_S(i) = dmf_dotp(m, sqphi(:, 1, i) , oep%vxc)
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Broadcast the vector v_bar_S  and sqphi to all processors
    do i = 1, st%nst
      call MPI_Bcast(v_bar_S(i), 1, MPI_FLOAT, st%node(i), st%mpi_grp%comm, mpi_err)
    end do
    do i = 1, n
      kssi = oep%eigen_index(i)
      call MPI_Bcast(sqphi(1, 1, kssi), m%np, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
    end do
  end if
#endif

  ! If there is more than one state, so solve linear equation.
  linear_equation: if(n > 0) then
    ALLOCATE(d(m%np), m%np)
    ALLOCATE(x(n, 1), n*1)
    ALLOCATE(Ma(n, n), n*n)
    ALLOCATE(y(n, 1), n*1)
    x = M_ZERO; y = M_ZERO; Ma = M_ZERO; d = M_ZERO
    proc = st%mpi_grp%rank

    i_loop: do i = 1, n
      kssi = oep%eigen_index(i)
      if(proc .eq. st%node(kssi)) then
        d(1:m%np) = sqphi(1:m%np, 1, kssi) / rho_sigma(1:m%np)
        j_loop: do j = i, n
          kssj = oep%eigen_index(j)
          Ma(i, j) = - dmf_dotp(m, d, sqphi(:, 1, kssj) )
        end do j_loop
        Ma(i,i) = M_ONE + Ma(i,i)
        y(i, 1) = v_bar_S(kssi) - oep%uxc_bar(kssi)
      end if
    end do i_loop

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      do i = 1, n
        kssi = oep%eigen_index(i)
        call MPI_Bcast(y(i, 1), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        do j = 1, n
           call MPI_Bcast(Ma(i, j), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        end do
     end do
    end if
#endif

    do i = 1, n
      do j = i, n
        Ma(j, i) = Ma(i, j)
      end do
    end do

    call lalg_linsyssolve(n, 1, Ma, y, x)

    do i = 1, n
      kssi = oep%eigen_index(i)
      occ = st%occ(kssi, is)
      oep%vxc(1:m%np) = oep%vxc(1:m%np) + &
        oep%socc * occ * x(i,1) * sqphi(1:m%np, 1, kssi ) / rho_sigma(1:m%np)
    end do

    deallocate(d, x, Ma, y)

  end if linear_equation
  ! The previous stuff is only needed if n>0.

  deallocate(v_bar_S, rho_sigma, sqphi)
  call pop_sub()
  call profiling_out(C_PROFILING_XC_KLI)
end subroutine X(xc_KLI_solve)

