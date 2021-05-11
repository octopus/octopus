!! Copyright (C) 2021 N. Tancogne-Dejean
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
subroutine X(x_fbe_calc)(psolver, mesh, der, st, ex, vxc)
  type(poisson_t),          intent(in)    :: psolver
  type(mesh_t),             intent(in)    :: mesh
  type(derivatives_t),      intent(in)    :: der
  type(states_elec_t),      intent(inout) :: st
  FLOAT,                    intent(inout) :: ex
  FLOAT, optional,          intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  FLOAT :: eig
  integer :: isp
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(X_FBE)))
  PUSH_SUB(X(x_fbe_calc))

  ASSERT(st%d%ispin /= SPINORS)

  eig = M_ZERO
  !At the moment we treat only spins and not k-points
  do isp = 1, st%d%nspin
    call X(fbe)(mesh, der, psolver, st, isp, eig, vxc)
  end do
  ex = ex + eig

  POP_SUB(X(x_fbe_calc))
  call profiling_out(prof)
end subroutine X(x_fbe_calc)

!------------------------------------------------------------
!> This routine is adapted from the X(slater) routine
!------------------------------------------------------------
subroutine X(fbe) (mesh, der, psolver, st, isp, ex, vxc)
  type(mesh_t),                intent(in)    :: mesh
  type(derivatives_t),         intent(in)    :: der
  type(poisson_t),             intent(in)    :: psolver
  type(states_elec_t), target, intent(inout) :: st
  integer,                     intent(in)    :: isp
  FLOAT,                       intent(inout) :: ex
  FLOAT, optional,             intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)

  integer :: ii, jst, ist, i_max, node_to, node_fr, ist_s, ist_r, idm, ip
  integer, allocatable :: recv_stack(:), send_stack(:)
  FLOAT :: rr, socc, tmp
  R_TYPE, allocatable :: rho_ij(:), pot_ij(:), psi(:,:), wf_ist(:,:), grad_rho_ij(:,:)
  FLOAT, allocatable :: tmp_vxc(:,:), tmp_pot(:,:), div(:), grad_rho(:,:)


#if defined(HAVE_MPI)
  integer :: send_req, status(MPI_STATUS_SIZE)
#endif

  PUSH_SUB(X(fbe))

  socc = M_ONE / st%smear%el_per_state

  SAFE_ALLOCATE(pot_ij(1:mesh%np))
  SAFE_ALLOCATE(grad_rho_ij(1:mesh%np, 1:mesh%sb%dim))
  SAFE_ALLOCATE(rho_ij(1:mesh%np_part))
  SAFE_ALLOCATE(recv_stack(1:st%nst+1))
  SAFE_ALLOCATE(send_stack(1:st%nst+1))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(wf_ist(1:mesh%np, 1:st%d%dim))

  if(present(vxc)) then
    SAFE_ALLOCATE(tmp_pot(1:mesh%np, 1:mesh%sb%dim))
    SAFE_ALLOCATE(tmp_vxc(1:mesh%np_part, 1:mesh%sb%dim))
    tmp_vxc(:,:) = M_ZERO
    SAFE_ALLOCATE(div(1:mesh%np_part))
    div(:) = M_ZERO
    !We need the gradient of the density to evaluate the potential
    SAFE_ALLOCATE(grad_rho(1:mesh%np, 1:mesh%sb%dim))
    call dderivatives_grad(der, st%rho(:,1), grad_rho)
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
          rho_ij(1:mesh%np) = R_CONJ(wf_ist(1:mesh%np, 1))*psi(1:mesh%np, 1)

          pot_ij(1:mesh%np) = R_TOTYPE(M_ZERO)
          call X(poisson_solve)(psolver, pot_ij, rho_ij, all_nodes=.false.)

          rr = socc * st%occ(ist, isp) * st%occ(jst, isp)
          if(ist /= jst) rr = rr * M_TWO

          !Here we accumulate the result for the potential
          if(present(vxc)) then
            !There are two contributitions, the gradient of the Slater potential, and the 
            !Poisson equation with the gradient of the exchange hole

            !We first compute the Slater contribution
            do ip = 1, mesh%np
              !If there is no density at this point, we simply ignore it
              div(ip) = div(ip) - rr * R_REAL(R_CONJ(rho_ij(ip)) * pot_ij(ip)) / (st%rho(ip, isp) + M_EPSILON)
            end do

            call X(derivatives_grad)(der, rho_ij(:), grad_rho_ij(:, :))

            do idm = 1, mesh%sb%dim
              do ip = 1, mesh%np
                !For numerical reasons, we need to compute the gradient of the \rho_ij/n this way
                !else, we have severe numerical troubles at the border of the simulation box
                !especially in the case of one electron only (\rho_ij = n)
                tmp = R_REAL(R_CONJ(grad_rho_ij(ip,idm)) * pot_ij(ip)) * st%rho(ip, isp) &
                     -R_REAL(R_CONJ(rho_ij(ip)) * pot_ij(ip)) * grad_rho(ip, idm) 

                tmp_vxc(ip, idm) = tmp_vxc(ip, idm)  - rr * tmp / ( st%rho(ip, isp)**2 + M_EPSILON)
              end do
            end do

          end if
           

          ! not an energy functional
          ex = M_ZERO 
        end do
      end if

      ! all done?
      if((send_stack(ist_s) < 0) .and. (recv_stack(ist_r) < 0)) exit
    end do
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    if(present(vxc)) then
      call comm_allreduce(st%mpi_grp, tmp_vxc)
      call comm_allreduce(st%mpi_grp, div)
    end if
  end if
#endif
  if(present(vxc)) then
    !First we add the part from the Slater potential
    vxc(1:mesh%np, isp) = div(1:mesh%np)

    call dderivatives_div(der, tmp_vxc(:, :), div(:)) 
    tmp_vxc(:, 1) = M_ZERO
    call dpoisson_solve(psolver, tmp_vxc(:, 1), div(:)) 
    vxc(1:mesh%np, isp) = vxc(1:mesh%np, isp) + tmp_vxc(1:mesh%np, 1) / M_PI / M_FOUR
    SAFE_DEALLOCATE_A(div)
  end if
  SAFE_DEALLOCATE_A(tmp_vxc)

  SAFE_DEALLOCATE_A(recv_stack)
  SAFE_DEALLOCATE_A(send_stack)
  SAFE_DEALLOCATE_A(pot_ij)
  SAFE_DEALLOCATE_A(rho_ij)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(wf_ist)

  POP_SUB(X(fbe))
end subroutine X(fbe)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
