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
! ---------------------------------------------------------
subroutine X(xc_slater_calc)(namespace, psolver, mesh, st, ex, vxc)
  type(namespace_t),        intent(in)    :: namespace
  type(poisson_t),          intent(in)    :: psolver
  type(mesh_t),             intent(in)    :: mesh
  type(states_elec_t),      intent(inout) :: st
  FLOAT,                    intent(inout) :: ex
  FLOAT, optional,          intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  FLOAT :: eig
  integer :: isp
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(XC_SLATER)))
  PUSH_SUB(X(xc_slater_calc))

  eig = M_ZERO
  if(st%d%ispin == SPIN_POLARIZED) then
    !At the moment we treat only spins and not k-points
    do isp = 1, st%d%nspin
      call X(slater) (namespace, mesh, psolver, st, isp, eig, vxc)
    end do
  else
    call X(slater) (namespace, mesh, psolver, st, 1, eig, vxc)
  end if
  ex = ex + eig

  POP_SUB(X(xc_slater_calc))
  call profiling_out(prof)
end subroutine X(xc_slater_calc)

!------------------------------------------------------------
!> This routine is adapted from the X(oep_x) routine
!------------------------------------------------------------
subroutine X(slater) (namespace, mesh, psolver, st, isp, ex, vxc)
  type(namespace_t),           intent(in)    :: namespace
  type(mesh_t),                intent(in)    :: mesh
  type(poisson_t),             intent(in)    :: psolver
  type(states_elec_t), target, intent(in)    :: st
  integer,                     intent(in)    :: isp
  FLOAT,                       intent(inout) :: ex
  FLOAT, optional,             intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  integer :: ii, jst, ist, i_max, node_to, node_fr, ist_s, ist_r, idm, ip
  integer, allocatable :: recv_stack(:), send_stack(:)
  FLOAT :: rr, socc, nn, mm, alpha, betar, betai, alpha2, beta2
  R_TYPE, allocatable:: bij(:,:)
  R_TYPE, allocatable :: rho_ij(:), pot_ij(:), psi(:,:), wf_ist(:,:)
  CMPLX :: tmp
  FLOAT, allocatable :: tmp_vxc(:)
  FLOAT :: global_b(4), local_b(4), local_v(4), global_v(4)
  FLOAT :: nup, ndn, sqmod_updn


#if defined(HAVE_MPI)
  integer :: send_req, status(MPI_STATUS_SIZE)
#endif

  PUSH_SUB(X(slater))

  socc = M_ONE / st%smear%el_per_state

  SAFE_ALLOCATE(pot_ij(1:mesh%np))
  SAFE_ALLOCATE(rho_ij(1:mesh%np))
  SAFE_ALLOCATE(recv_stack(1:st%nst+1))
  SAFE_ALLOCATE(send_stack(1:st%nst+1))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(wf_ist(1:mesh%np, 1:st%d%dim))
  if(st%d%ispin == SPINORS) then
    SAFE_ALLOCATE(bij(1:mesh%np, 1:3))
    bij(:,:) = M_ZERO
  end if

  if(present(vxc) .and. st%d%ispin /= SPINORS) then
    SAFE_ALLOCATE(tmp_vxc(1:mesh%np))
    tmp_vxc(:) = M_ZERO
  end if

  ! This is the maximum number of blocks for each processor
  i_max = int((st%mpi_grp%size + 2)/2) - 1

  do ii = 0, i_max
    ! node where to send the wavefunctions
    node_to = mod(st%mpi_grp%rank + ii, st%mpi_grp%size)
    if(ii == i_max .and. mod(st%mpi_grp%size, 2) == 0 .and. node_to < st%mpi_grp%size/2) then
      node_to = -1
    end if

    ! node from which we receive the wavefunctions
    node_fr = st%mpi_grp%rank - ii
    if(node_fr < 0) node_fr = st%mpi_grp%size + node_fr
    node_fr = mod(node_fr, st%mpi_grp%size)
    if(ii == i_max .and. mod(st%mpi_grp%size, 2) == 0 .and. st%mpi_grp%rank < st%mpi_grp%size/2) then
      node_fr = -1
    end if

    ! check which wavefunctions we have to send/recv, and put them in a stack
    recv_stack(:) = -1
    ist_r = 1
    send_stack(:) = -1
    ist_s = 1

    do jst = 1, st%nst
      if(st%node(jst) == node_fr) then
        recv_stack(ist_r) = jst
        ist_r = ist_r + 1
      end if
      if(node_to /= -1 .and. st%node(jst) == st%mpi_grp%rank) then
        send_stack(ist_s) = jst
        ist_s = ist_s + 1
      end if
    end do

    ! now we do a loop over all states
    ist_r = 0
    if(node_fr < 0) ist_r = st%nst + 1
    ist_s = 0
    if(node_to < 0) ist_s = st%nst + 1
    do
      ! increment send counter
      if(ist_s <= st%nst) ist_s = ist_s  + 1

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ! send wavefunction
        send_req = 0
        if((send_stack(ist_s) > 0).and.(node_to /= st%mpi_grp%rank)) then
          call states_elec_get_state(st, mesh, send_stack(ist_s), isp, psi)
          call MPI_Isend(psi, mesh%np*st%d%dim, R_MPITYPE, &
            node_to, send_stack(ist_s), st%mpi_grp%comm, send_req, mpi_err)
        end if
      end if
#endif

      ! increment receive counter
      if(ist_r <= st%nst) ist_r = ist_r  + 1

      ! receive wavefunction
      if(recv_stack(ist_r) > 0) then
        if(node_fr == st%mpi_grp%rank) then
          call states_elec_get_state(st, mesh, send_stack(ist_r), isp, wf_ist)
#if defined(HAVE_MPI)
        else
          if(st%parallel_in_states) then
            call MPI_Recv(wf_ist, mesh%np*st%d%dim, R_MPITYPE, &
              node_fr, recv_stack(ist_r), st%mpi_grp%comm, status, mpi_err)
          end if
#endif
        end if
      end if

#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        if(send_req /= 0) call MPI_Wait(send_req, status, mpi_err)
        send_req = 0
      end if
#endif

      if(recv_stack(ist_r) > 0) then
        ! this is where we calculate the elements of the matrix
        ist = recv_stack(ist_r)
        do jst = st%st_start, st%st_end

          if(abs(st%occ(ist, isp)) < M_EPSILON .or. abs(st%occ(jst, isp)) < M_EPSILON) cycle

          if((st%node(ist) == st%mpi_grp%rank).and.(jst < ist)) cycle
          if((st%occ(ist, isp) <= M_EPSILON).or.(st%occ(jst, isp) <= M_EPSILON)) cycle

          call states_elec_get_state(st, mesh, jst, isp, psi)
          !The co-density is summed over the spin indices before solving the Poisson equation
          rho_ij(1:mesh%np) = R_CONJ(wf_ist(1:mesh%np, 1))*psi(1:mesh%np, 1)
          do idm = 2, st%d%dim
            rho_ij(1:mesh%np) = rho_ij(1:mesh%np) &
                           + R_CONJ(wf_ist(1:mesh%np, idm))*psi(1:mesh%np, idm)
          end do
          pot_ij(1:mesh%np) = R_TOTYPE(M_ZERO)
          call X(poisson_solve)(psolver, pot_ij, rho_ij, all_nodes=.false.)

          rr = socc * st%occ(ist, isp) * st%occ(jst, isp)
          if(ist /= jst .and. st%d%ispin /= SPINORS) rr = rr * M_TWO

          !Here we accumulate the result for the potential
          if(present(vxc)) then
            if( st%d%ispin /= SPINORS) then
              do ip = 1, mesh%np
                !If there is no density at this point, we simply ignore it
                tmp_vxc(ip) = tmp_vxc(ip) - rr * R_REAL(R_CONJ(rho_ij(ip)) * pot_ij(ip)) / (st%rho(ip, isp) + M_EPSILON)
              end do
            else
              !This is given by Eq. 19 in SI of PRB 98, 035140 (2018)
              if(ist /= jst) then
                bij(1:mesh%np, 1) = bij(1:mesh%np, 1) - M_TWO * rr &
                      * R_REAL(wf_ist(1:mesh%np, 1)*R_CONJ(psi(1:mesh%np, 1)) * pot_ij(1:mesh%np))
                bij(1:mesh%np, 2) = bij(1:mesh%np, 2) - M_TWO * rr &
                      * R_REAL(wf_ist(1:mesh%np, 2)*R_CONJ(psi(1:mesh%np, 2)) * pot_ij(1:mesh%np))
                !As we only compute the terms ist >= jst, we get a symmetric form
                bij(1:mesh%np, 3) = bij(1:mesh%np, 3) - rr &
                      * ( wf_ist(1:mesh%np, 1)* R_CONJ(psi(1:mesh%np, 2)) * pot_ij(1:mesh%np) &
                        + psi(1:mesh%np, 1)   * R_CONJ(wf_ist(1:mesh%np, 2) * pot_ij(1:mesh%np)))
              else
                bij(1:mesh%np, 1) = bij(1:mesh%np, 1) - rr &
                      * R_REAL(wf_ist(1:mesh%np, 1)*R_CONJ(psi(1:mesh%np, 1)) * pot_ij(1:mesh%np))
                bij(1:mesh%np, 2) = bij(1:mesh%np, 2) - rr &
                      * R_REAL(wf_ist(1:mesh%np, 2)*R_CONJ(psi(1:mesh%np, 2)) * pot_ij(1:mesh%np))
                bij(1:mesh%np, 3) = bij(1:mesh%np, 3) - rr &
                      * wf_ist(1:mesh%np, 1) * R_CONJ(psi(1:mesh%np, 2)) * pot_ij(1:mesh%np)
              end if
              !The last component is simply the complex conjuguate of bij(3), so we do not compute it.
            end if
          end if
          ! get the contribution (ist, jst) to the exchange energy
          ex = ex - M_HALF * rr * R_REAL(X(mf_dotp)(mesh, rho_ij, pot_ij))
        end do
      end if

      ! all done?
      if((send_stack(ist_s) < 0) .and. (recv_stack(ist_r) < 0)) exit
    end do
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! sum all contributions to the exchange energy
    call MPI_Allreduce(ex, rr, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    ex = rr
    if(present(vxc)) then
      if(st%d%ispin == SPINORS) then
        call comm_allreduce(st%mpi_grp, bij)
      else
        call comm_allreduce(st%mpi_grp, tmp_vxc)
      end if
    end if
  end if
#endif
  if(present(vxc) .and. st%d%ispin /= SPINORS) then
    vxc(:, isp) = vxc(:, isp) + tmp_vxc
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
        !In this frame, the spin density matrix has a single non-zero element on the diagonal (because wthe matrix is singular)
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
  

  SAFE_DEALLOCATE_A(recv_stack)
  SAFE_DEALLOCATE_A(send_stack)
  SAFE_DEALLOCATE_A(pot_ij)
  SAFE_DEALLOCATE_A(rho_ij)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(wf_ist)
  SAFE_DEALLOCATE_A(bij)

  POP_SUB(X(slater))
end subroutine X(slater)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
