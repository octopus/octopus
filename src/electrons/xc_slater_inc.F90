!! Copyright (C) 2020 N. Tancogne-Dejean
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!!
!------------------------------------------------------------
subroutine X(slater_calc) (namespace, mesh, space, exxop, st, kpoints, ex, vxc)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(space_t),               intent(in)    :: space
  type(exchange_operator_t),   intent(in)    :: exxop
  type(states_elec_t),         intent(in)    :: st
  type(kpoints_t),             intent(in)    :: kpoints
  FLOAT,                       intent(inout) :: ex
  FLOAT, optional,             intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  integer :: ist, ip, isp, ik
  FLOAT :: rr, nn, mm, alpha, betar, betai, alpha2, beta2
  R_TYPE, allocatable:: bij(:,:)
  R_TYPE, allocatable :: psi(:,:), xpsi(:,:)
  CMPLX :: tmp
  FLOAT, allocatable :: tmp_vxc(:,:)
  FLOAT :: nup, ndn, sqmod_updn, weight, eig
  FLOAT :: global_b(4), local_b(4), local_v(4), global_v(4)
  type(states_elec_t) :: xst


  PUSH_SUB(X(slater_calc))

  !We first apply the exchange operator to all the states
  call xst%nullify()
  eig = M_ZERO
  call X(exchange_operator_compute_potentials)(exxop, namespace, space, mesh, st, xst, kpoints)

  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(xpsi(1:mesh%np, 1:st%d%dim))
  if(st%d%ispin == SPINORS) then
    SAFE_ALLOCATE(bij(1:mesh%np, 1:3))
    bij(:,:) = M_ZERO
  end if

  if(present(vxc) .and. st%d%ispin /= SPINORS) then
    SAFE_ALLOCATE(tmp_vxc(1:mesh%np, 1:st%d%spin_channels))
    tmp_vxc(:,:) = M_ZERO
  end if

  !We now construct the right-hand side of the equation defining the Slater potential
  do ik = st%d%kpt%start, st%d%kpt%end
    isp = st%d%get_spin_index(ik)
    do ist = st%st_start, st%st_end
      weight = st%occ(ist, ik) * st%d%kweights(ik)  
      if(abs(weight) < M_EPSILON) cycle

      call states_elec_get_state(xst, mesh, ist, ik, xpsi)
      call states_elec_get_state(st, mesh, ist, ik, psi)
      !Here we accumulate the result for the potential
      if(present(vxc)) then
        if( st%d%ispin /= SPINORS) then
          do ip = 1, mesh%np
            !If there is no density at this point, we simply ignore it
            tmp_vxc(ip, isp) = tmp_vxc(ip, isp) + weight* R_REAL(R_CONJ(psi(ip,1))*xpsi(ip,1)) / (st%rho(ip, isp) + M_EPSILON)
          end do
        else
          do ip = 1, mesh%np
            bij(ip, 1) = bij(ip, 1) + weight * R_REAL(xpsi(ip, 1)*R_CONJ(psi(ip, 1)))
            bij(ip, 2) = bij(ip, 2) + weight * R_REAL(xpsi(ip, 2)*R_CONJ(psi(ip, 2)))
            bij(ip, 3) = bij(ip, 3) + weight * xpsi(ip, 1) * R_CONJ(psi(ip, 2))
          end do
          !The last component is simply the complex conjuguate of bij(3), so we do not compute it.
        end if
      end if

      ! get the contribution to the exchange energy
      eig = eig + M_HALF * weight * R_REAL(X(mf_dotp)(mesh, st%d%dim, xpsi, psi))
    end do
  end do

  call states_elec_end(xst)

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    ! sum all contributions to the exchange energy
    call comm_allreduce(st%st_kpt_mpi_grp, eig)
    if(present(vxc)) then
      if(st%d%ispin == SPINORS) then
        call comm_allreduce(st%st_kpt_mpi_grp, bij)
      else
        call comm_allreduce(st%st_kpt_mpi_grp, tmp_vxc)
      end if
    end if
  end if
#endif
  ex = ex + eig

  if(present(vxc) .and. st%d%ispin /= SPINORS) then
    call lalg_axpy(mesh%np, st%d%spin_channels, M_ONE, tmp_vxc, vxc)
  end if
  SAFE_DEALLOCATE_A(tmp_vxc)

  !We now construct the full potential from the density and the quantity b_ij
  !Note that the final potential is real and stored as v_upup, v_downdown, Re(v_updown), and Im(updown)
  !See Eq. 22 in SI of PRB 98, 035140 (2018)
  if(present(vxc) .and. st%d%ispin == SPINORS) then
    do ip = 1, mesh%np

      nn = st%rho(ip, 1) + st%rho(ip, 2) + M_EPSILON
      sqmod_updn = st%rho(ip, 3)**2 + st%rho(ip, 4)**2

      ! 1/(2nD), where n is the charge density and D = n_uu*n_dd - n_ud*n_du
      ! where D is the determinant of the spin-density matrix
      rr = (st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn)
      if(abs(rr) < CNST(1.0e-13) ) then !The matrix is singular
        !For determining the potential, we go to the local frame given by the local magnetization
        !In this frame, the spin density matrix has a single non-zero element on the diagonal (because the matrix is singular)
        !This allows us to determine the potential in this frame, and to rotate it back in the original frame
        !We start by finding the non-zero component in the local rotated frame
        mm = sqrt((st%rho(ip, 1) - st%rho(ip, 2))**2 + M_FOUR*sqmod_updn)
        if((st%rho(ip, 1) - st%rho(ip, 2)) > 0) then
          nup = max(M_HALF*(nn + mm), M_ZERO)
          ndn = max(M_HALF*(nn - mm), M_ZERO)
        else
          nup = max(M_HALF*(nn - mm), M_ZERO)
          ndn = max(M_HALF*(nn + mm), M_ZERO)
        end if

        !Coefficients of the rotation matrix
        call get_rotation_matrix(st%rho(ip, :), alpha, betar, betai)
        alpha2 = alpha**2
        beta2 = betar**2 + betai**2

        !We rotate to the local frame
        !We add here the factor of two that is originally in the Sylvester equation
        !In the non singular case, this is cancelled by a factor one half
        global_b(1) = M_TWO*TOFLOAT(bij(ip, 1))
        global_b(2) = M_TWO*TOFLOAT(bij(ip, 2))
        global_b(3) = M_TWO*TOFLOAT(bij(ip, 3))
        global_b(4) = M_TWO*R_AIMAG(bij(ip, 3))
        call rotate_to_local(global_b, alpha, betar, betai, alpha2, beta2, local_b)
       
        local_v(1:4) = M_ZERO
        if(nup > ndn) then !We are up
          local_v(1) = local_b(1) / (M_TWO * nup)
          local_v(3) = local_b(3) / nup
          local_v(4) = local_b(4) / nup
        else !We are down
          local_v(2) = local_b(2) / (M_TWO * ndn)
          local_v(3) = local_b(3) / ndn
          local_v(4) = local_b(4) / ndn
        end if

        !We rotate to the original frame and we accumulate the result
        call rotate_to_global(local_v, alpha, betar, betai, alpha2, beta2, global_v)
        vxc(ip, 1:4) = vxc(ip, 1:4) + global_v(1:4)   

      else

        rr = M_ONE/(nn * rr)

        vxc(ip, 1) = vxc(ip, 1) + rr * ( &
           (nn * st%rho(ip, 2) - sqmod_updn) * R_REAL(bij(ip, 1)) + sqmod_updn * R_REAL(bij(ip, 2)) &
          - M_TWO * st%rho(ip,2) * ( st%rho(ip,3) * R_REAL(bij(ip,3)) + st%rho(ip,4) * R_AIMAG(bij(ip,3))))

        vxc(ip, 2) = vxc(ip, 2) + rr * ( &
           (nn * st%rho(ip, 1) - sqmod_updn) * R_REAL(bij(ip, 2)) + sqmod_updn * R_REAL(bij(ip, 1)) &
          - M_TWO * st%rho(ip,1) * ( st%rho(ip,3) * R_REAL(bij(ip,3)) + st%rho(ip,4) * R_AIMAG(bij(ip,3)))) 

        tmp = -TOCMPLX(st%rho(ip, 3), st%rho(ip,4)) * (st%rho(ip, 2) * bij(ip, 1) + st%rho(ip, 1) * bij(ip,2)) &
              + (M_TWO *st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn) * bij(ip, 3) &
              + (TOCMPLX(st%rho(ip, 3),st%rho(ip,4)))**2 * R_CONJ(bij(ip,3))

        vxc(ip, 3) = vxc(ip, 3) + rr * real(tmp)
        vxc(ip, 4) = vxc(ip, 4) + rr * aimag(tmp)
      end if
    end do 
  end if
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(xpsi)
  SAFE_DEALLOCATE_A(bij)

  POP_SUB(X(slater_calc))
end subroutine X(slater_calc)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
