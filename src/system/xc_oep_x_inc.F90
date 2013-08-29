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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$


!------------------------------------------------------------
!> The parallelization of this routine is done in the following way:
!!
!! We have to calculate the sum
!! \f$ lxc = \sum_{j>i} l_{ij} \f$
!! where the states i and j are divided in blocks and scattered among
!! the processors. Each processor will calculate a sub-block of the
!! matrix \f$ l_{ij} \f$. Examples of the partitioning (for 3 and 4 blocks/processors)
!! \verbatim
!! (1  2  1)    (1  2  3  1)
!! (-  2  3)    (-  2  3  4)
!! (-  -  3)    (-  -  3  4)
!!  -  -  -  -  (-  -  -  4)
!! \endverbatim
!!  where the numbers indicate the processor that will do the work
!------------------------------------------------------------

subroutine X(oep_x) (der, st, is, jdm, lxc, ex, exx_coef)
  type(derivatives_t), intent(inout) :: der
  type(states_t), target, intent(in) :: st
  integer,        intent(in)    :: is
  integer,        intent(in)    :: jdm
  R_TYPE,         intent(inout) :: lxc(:, st%st_start:, :) !< (1:der%mesh%np, :st%st_end, nspin)
  FLOAT,          intent(inout) :: ex
  FLOAT,          intent(in)    :: exx_coef !< amount of EXX (for hybrids)

  integer :: ii, jst, ist, i_max, node_to, node_fr, ist_s, ist_r, isp, idm
  integer, allocatable :: recv_stack(:), send_stack(:)
  FLOAT :: rr, socc
  R_TYPE, pointer     :: send_buffer(:)
  R_TYPE, allocatable :: rho_ij(:), F_ij(:), psi(:), wf_ist(:)

#if defined(HAVE_MPI)
  integer :: send_req, status(MPI_STATUS_SIZE)
#endif

  call profiling_in(C_PROFILING_XC_EXX)
  PUSH_SUB(X(oep_x))

  if(der%mesh%sb%kpoints%reduced%npoints > 1) call messages_not_implemented("exchange operator with k-points")

  socc = M_ONE / st%smear%el_per_state
  !
  ! distinguish between 'is' being the spin_channel index (collinear)
  ! and being the spinor (noncollinear)
  if (st%d%ispin==SPINORS) then
    isp = 1
    idm = is
  else
    isp = is
    idm = 1
  end if 
  ! Note: we assume that st%occ is known in all nodes

  SAFE_ALLOCATE(F_ij(1:der%mesh%np))
  SAFE_ALLOCATE(rho_ij(1:der%mesh%np))
  SAFE_ALLOCATE(send_buffer(1:der%mesh%np))
  SAFE_ALLOCATE(recv_stack(1:st%nst+1))
  SAFE_ALLOCATE(send_stack(1:st%nst+1))
  SAFE_ALLOCATE(psi(der%mesh%np))
  SAFE_ALLOCATE(wf_ist(der%mesh%np))

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
          call states_get_state(st, der%mesh, jdm, send_stack(ist_s), isp, psi)
          call MPI_Isend(psi, der%mesh%np, R_MPITYPE, &
            node_to, send_stack(ist_s), st%mpi_grp%comm, send_req, mpi_err)
        end if
      end if
#endif

      ! increment receive counter
      if(ist_r <= st%nst) ist_r = ist_r  + 1

      ! receive wavefunction
      if(recv_stack(ist_r) > 0) then
        if(node_fr == st%mpi_grp%rank) then
          call states_get_state(st, der%mesh, jdm, send_stack(ist_r), isp, wf_ist)
          if(st%cmplxscl%space) wf_ist = R_CONJ(wf_ist)
#if defined(HAVE_MPI)
        else
          if(st%parallel_in_states) then
            call MPI_Recv(wf_ist, der%mesh%np, R_MPITYPE, &
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
        send_buffer(1:der%mesh%np) = R_TOTYPE(M_ZERO)
        do jst = st%st_start, st%st_end

          if((st%node(ist) == st%mpi_grp%rank).and.(jst < ist).and..not.(st%d%ispin==SPINORS)) cycle
          if((st%occ(ist, isp) <= M_EPSILON).or.(st%occ(jst, isp) <= M_EPSILON)) cycle

          call states_get_state(st, der%mesh, jdm, jst, isp, psi)
          rho_ij(1:der%mesh%np) = R_CONJ(wf_ist(1:der%mesh%np))*psi(1:der%mesh%np)
          F_ij(1:der%mesh%np) = R_TOTYPE(M_ZERO)
          call X(poisson_solve)(psolver, F_ij, rho_ij, all_nodes=.false.)

          ! this quantity has to be added to lxc(1:der%mesh%np, ist)
          call states_get_state(st, der%mesh, idm, jst, isp, psi)
          send_buffer(1:der%mesh%np) = send_buffer(1:der%mesh%np) + &
            socc*st%occ(jst, isp)*F_ij(1:der%mesh%np)*R_CONJ(psi(1:der%mesh%np))

          ! if off-diagonal, then there is another contribution
          ! note that the wf jst is always in this node
          if((ist /= jst).and..not.(st%d%ispin==SPINORS)) then
            lxc(1:der%mesh%np, jst, is) = lxc(1:der%mesh%np, jst, is) - &
              exx_coef * socc * st%occ(ist, isp) * R_CONJ(F_ij(1:der%mesh%np)*wf_ist(1:der%mesh%np))
          end if
          ! get the contribution (ist, jst) to the exchange energy
          rr = M_ONE
          if(ist /= jst .and. .not.(st%d%ispin==SPINORS)) rr = M_TWO

          ex = ex - exx_coef* M_HALF * rr * &
              st%occ(ist, isp) * socc*st%occ(jst, isp) * &
              R_REAL(X(mf_dotp)(der%mesh, psi(1:der%mesh%np), wf_ist(:)*F_ij(:)))
        end do

        if(st%node(ist) == st%mpi_grp%rank) then
          ! either add the contribution ist
          lxc(1:der%mesh%np, ist, is) = lxc(1:der%mesh%np, ist, is) - exx_coef * send_buffer(1:der%mesh%np)

#if defined(HAVE_MPI)
        else
          if(st%parallel_in_states) then
            ! or send it to the node that has wf ist
            call MPI_Isend(send_buffer(1), der%mesh%np, R_MPITYPE, &
             node_fr, ist, st%mpi_grp%comm, send_req, mpi_err)
          end if
#endif
        end if

      end if

#if defined(HAVE_MPI)
      ! now we have to receive the contribution to lxc from the node to
      ! which we sent the wavefunction ist
      if(st%parallel_in_states) then
        if((node_to >= 0) .and. (send_stack(ist_s) > 0) .and. (node_to /= st%mpi_grp%rank)) then
          call MPI_Recv(psi(:), der%mesh%np, R_MPITYPE, &
            node_to, send_stack(ist_s), st%mpi_grp%comm, status, mpi_err)

          lxc(1:der%mesh%np, send_stack(ist_s), is) = lxc(1:der%mesh%np, send_stack(ist_s), is) - &
            exx_coef * psi(1:der%mesh%np)
        end if

        if(send_req /= 0) call MPI_Wait(send_req, status, mpi_err)
      end if
#endif

      ! all done?
      if((send_stack(ist_s) < 0) .and. (recv_stack(ist_r) < 0)) exit
    end do
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! sum all contributions to the exchange energy
    call MPI_Allreduce(ex, rr, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    ex = rr
  end if
#endif

  SAFE_DEALLOCATE_A(recv_stack)
  SAFE_DEALLOCATE_A(send_stack)
  SAFE_DEALLOCATE_A(F_ij)
  SAFE_DEALLOCATE_A(rho_ij)
  SAFE_DEALLOCATE_P(send_buffer)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(wf_ist)

  call profiling_out(C_PROFILING_XC_EXX)
  POP_SUB(X(oep_x))
end subroutine X(oep_x)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
