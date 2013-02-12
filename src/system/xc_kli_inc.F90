!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
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
subroutine X(xc_KLI_solve) (mesh, st, is, oep)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  integer,        intent(in)    :: is
  type(xc_oep_t), intent(inout) :: oep

  integer :: ist, ip, jst, eigen_n, kssi, kssj, proc
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:), sqphi(:, :, :), dd(:)
  FLOAT, allocatable :: Ma(:,:), xx(:,:), yy(:,:)
  R_TYPE :: occ

  call profiling_in(C_PROFILING_XC_KLI)
  PUSH_SUB(X(xc_KLI_solve))
  ! some intermediate quantities
  ! vxc contains the Slater part!
  SAFE_ALLOCATE(rho_sigma(1:mesh%np))
  SAFE_ALLOCATE(sqphi(1:mesh%np, 1:st%d%dim, 1:st%nst))

  do ist = st%st_start, st%st_end
    sqphi(1:mesh%np, 1:st%d%dim, ist) = R_REAL (st%X(psi)(1:mesh%np, 1:st%d%dim, ist, is))**2 + &
                                        R_AIMAG(st%X(psi)(1:mesh%np, 1:st%d%dim, ist, is))**2
  end do

  do ip = 1, mesh%np
    rho_sigma(ip) = max(sum(oep%socc * st%occ(st%st_start:st%st_end, is) * &
      sqphi(ip, 1, st%st_start:st%st_end)), CNST(1e-20))
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    call MPI_Allreduce(rho_sigma(1), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    rho_sigma(1:mesh%np) = dd(1:mesh%np)
  end if
#endif

  do ip = 1, mesh%np
    oep%vxc(ip,1) = M_ZERO
    do ist = st%st_start, st%st_end
      oep%vxc(ip,1) = oep%vxc(ip,1) + oep%socc * st%occ(ist, is) * oep%X(lxc)(ip, ist, is) * st%X(psi)(ip, 1, ist, is)
    end do
    oep%vxc(ip,1) = oep%vxc(ip,1) / rho_sigma(ip)
  end do
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(oep%vxc(1,1), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    oep%vxc(1:mesh%np,1) = dd(1:mesh%np)
    SAFE_DEALLOCATE_A(dd)
  end if
#endif
  if(oep%level == XC_OEP_SLATER) then
    SAFE_DEALLOCATE_A(rho_sigma)
    SAFE_DEALLOCATE_A(sqphi)
    call profiling_out(C_PROFILING_XC_KLI)
    POP_SUB(X(xc_KLI_solve))
    return
  end if
  eigen_n = oep%eigen_n

  SAFE_ALLOCATE(v_bar_S(1:st%nst))
  do ist = st%st_start, st%st_end
    if(st%occ(ist, is) .gt. small) then
      v_bar_S(ist) = dmf_dotp(mesh, sqphi(:, 1, ist) , oep%vxc(:,1))
    end if
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Broadcast the vector v_bar_S  and sqphi to all processors
    do ist = 1, st%nst
      call MPI_Bcast(v_bar_S(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
    end do
    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      call MPI_Bcast(sqphi(1, 1, kssi), mesh%np, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
    end do
  end if
#endif
  ! If there is more than one state, then solve linear equation.
  linear_equation: if(eigen_n > 0) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    SAFE_ALLOCATE(xx(1:eigen_n, 1:1))
    SAFE_ALLOCATE(Ma(1:eigen_n, 1:eigen_n))
    SAFE_ALLOCATE(yy(1:eigen_n, 1:1))
    xx = M_ZERO
    yy = M_ZERO
    Ma = M_ZERO
    dd = M_ZERO
    proc = st%mpi_grp%rank

    i_loop: do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      if(proc .eq. st%node(kssi)) then
        dd(1:mesh%np) = sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
        j_loop: do jst = ist, eigen_n
          kssj = oep%eigen_index(jst)
          Ma(ist, jst) = - dmf_dotp(mesh, dd, sqphi(:, 1, kssj) )
        end do j_loop
        Ma(ist, ist) = M_ONE + Ma(ist, ist)
        yy(ist, 1) = v_bar_S(kssi) - oep%uxc_bar(kssi,is)
      end if
    end do i_loop

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      do ist = 1, eigen_n
        kssi = oep%eigen_index(ist)
        call MPI_Bcast(yy(ist, 1), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        do jst = 1, eigen_n
           call MPI_Bcast(Ma(ist, jst), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        end do
     end do
    end if
#endif

    do ist = 1, eigen_n
      do jst = ist, eigen_n
        Ma(jst, ist) = Ma(ist, jst)
      end do
    end do

    call lalg_linsyssolve(eigen_n, 1, Ma, yy, xx)

    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      occ = st%occ(kssi, is)
      oep%vxc(1:mesh%np,1) = oep%vxc(1:mesh%np,1) + &
        oep%socc * occ * xx(ist, 1) * sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
    end do

    SAFE_DEALLOCATE_A(dd)
    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(Ma)
    SAFE_DEALLOCATE_A(yy)

  end if linear_equation
  ! The previous stuff is only needed if eigen_n>0.

  SAFE_DEALLOCATE_A(v_bar_S)
  SAFE_DEALLOCATE_A(rho_sigma)
  SAFE_DEALLOCATE_A(sqphi)
  POP_SUB(X(xc_KLI_solve))
  call profiling_out(C_PROFILING_XC_KLI)
end subroutine X(xc_KLI_solve)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

