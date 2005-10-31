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


!------------------------------------------------------------
! The paralellization of this routine is done in the following way:
! We have to calculate the sum
!    lxc = sum_{j>i) l_ij
! where the states i and j are divided in blocks and scattered among
! the processors. Each processor will calculate a sub-block of the
! matrix l_ij. Examples of the partitioning (for 3 and 4 blocks/processors)
!
! (1  2  1)    (1  2  3  1)
! (   2  3)    (   2  3  4)
! (      3)    (      3  4)
!              (         4)
!
!  where the numbers indicate the processor that will do the work
!------------------------------------------------------------

subroutine X(oep_x) (gr, st, is, oep, ex)
  type(grid_type),   intent(inout) :: gr
  type(states_type), intent(in)    :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,             intent(out)   :: ex

  integer :: i, j, ist, jst, i_max, node_to, node_fr, ist_s, ist_r
  integer, allocatable :: recv_stack(:), send_stack(:)
  FLOAT :: r
  R_TYPE, pointer     :: wf_ist(:), send_buffer(:)
  R_TYPE, allocatable :: rho_ij(:), F_ij(:)

#if defined(HAVE_MPI)
  R_TYPE,  pointer     :: recv_buffer(:)
  integer :: ierr, send_req, status(MPI_STATUS_SIZE)
#endif

  ! Note: we assume that st%occ is known in all nodes

  call push_sub('xc_OEP_x.oep_x')

  allocate(F_ij(NP), rho_ij(NP), send_buffer(NP))
  allocate(recv_stack(st%nst+1), send_stack(st%nst+1))

#if defined(HAVE_MPI)
  allocate(recv_buffer(NP))
#endif

  ! This is the maximum number of blocks for each processor
  i_max = int((st%numprocs + 2)/2) - 1

  do i = 0, i_max
    ! node where to send the wavefunctions
    node_to = mod(st%rank + i, st%numprocs)
    if(i==i_max .and. mod(st%numprocs,2)==0 .and. node_to < st%numprocs/2) then
      node_to = -1
    end if

    ! node from which we receive the wave-functions
    node_fr = st%rank - i
    if(node_fr < 0) node_fr = st%numprocs + node_fr
    node_fr = mod(node_fr, st%numprocs)
    if(i==i_max .and. mod(st%numprocs, 2)==0 .and. st%rank < st%numprocs/2) then
      node_fr = -1
    end if

    ! check which wave-functions we have to send/recv, and put them in a stack
    recv_stack(:) = -1; ist_r = 1
    send_stack(:) = -1; ist_s = 1
    do j = 1, st%nst
      if(st%node(j) == node_fr) then
        recv_stack(ist_r) = j
        ist_r = ist_r + 1
      end if
      if(node_to.ne.-1 .and. st%node(j) == st%rank) then
        send_stack(ist_s) = j
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
      ! send wave-function
      send_req = 0
      if((send_stack(ist_s) > 0).and.(node_to.ne.st%rank)) then
        call MPI_Isend(st%X(psi)(:, 1, send_stack(ist_s), is), NP, R_MPITYPE, &
          node_to, send_stack(ist_s), st%comm, send_req, ierr)
      end if
#endif

      ! increment receive counter
      if(ist_r <= st%nst) ist_r = ist_r  + 1

      ! receive wave-function
      if(recv_stack(ist_r) > 0) then
        if(node_fr == st%rank) then
          wf_ist => st%X(psi)(1:NP, 1, recv_stack(ist_r), is)
#if defined(HAVE_MPI)
        else
          call MPI_Recv(recv_buffer, NP, R_MPITYPE, &
            node_fr, recv_stack(ist_r), st%comm, status, ierr)
          wf_ist => recv_buffer(1:NP)
#endif
        end if
      end if

#if defined(HAVE_MPI)
      if(send_req.ne.0) call MPI_Wait(send_req, status, ierr)
      send_req = 0
#endif

      if(recv_stack(ist_r) > 0) then
        ! this is where we calculate the elements of the matrix
        ist = recv_stack(ist_r)
        send_buffer(1:NP) = R_TOTYPE(M_ZERO)
        do jst = st%st_start, st%st_end
          if((st%node(ist) == st%rank).and.(jst < ist)) cycle
          if((st%occ(ist, is).le.small).or.(st%occ(jst, is).le.small)) cycle

          rho_ij(1:NP) = R_CONJ(wf_ist(1:NP))*st%X(psi)(1:NP, 1, jst, is)
          F_ij(1:NP) = R_TOTYPE(M_ZERO)
          call X(poisson_solve)(gr, F_ij, rho_ij)

          ! this quantity has to be added to oep%X(lxc)(1:NP, ist)
          send_buffer(1:NP) = send_buffer(1:NP) + &
            oep%socc*st%occ(jst, is)*F_ij(1:NP)*R_CONJ(st%X(psi)(1:NP, 1, jst, is))

          ! if off-diagonal, then there is another contribution
          ! not that the wf jst is always in this node
          if(ist.ne.jst) then
            oep%X(lxc)(1:NP, jst) = oep%X(lxc)(1:NP, jst) - &
              oep%socc * st%occ(ist, is) * R_CONJ(F_ij(1:NP)*wf_ist(1:NP))
          end if

          ! get the contribution (ist, jst) to the exchange energy
          r = M_ONE
          if(ist.ne.jst) r = M_TWO

          ex = ex - M_HALF * r * oep%sfact * oep%socc*st%occ(ist, is) * oep%socc*st%occ(jst, is) * &
            sum(R_REAL(wf_ist(1:NP) * F_ij(1:NP) * R_CONJ(st%X(psi)(1:NP, 1, jst, is))) * gr%m%vol_pp(1:NP))
        end do

        if(st%node(ist) == st%rank) then
          ! either add the contribution ist
          oep%X(lxc)(1:NP, ist) = oep%X(lxc)(1:NP, ist) - send_buffer(1:NP)

#if defined(HAVE_MPI)
        else
          ! or send it to the node that has wf ist
          call MPI_Isend(send_buffer(:), NP, R_MPITYPE, &
            node_fr, ist, st%comm, send_req, ierr)
#endif
        end if

      end if

#if defined(HAVE_MPI)
      ! now we have to receive the contribution to lxc from the node to
      ! which we sent the wave-function ist
      if((node_to >= 0).and.(send_stack(ist_s) > 0).and.(node_to.ne.st%rank)) then
        call MPI_Recv(recv_buffer(:), NP, R_MPITYPE, &
          node_to, send_stack(ist_s), st%comm, status, ierr)

        oep%X(lxc)(1:NP, send_stack(ist_s)) = oep%X(lxc)(1:NP, send_stack(ist_s)) - &
          recv_buffer(1:NP)
      end if

      if(send_req.ne.0) call MPI_Wait(send_req, status, ierr)
#endif

      ! all done?
      if((send_stack(ist_s) < 0).and.(recv_stack(ist_r) < 0)) exit
    end do
  end do

#if defined(HAVE_MPI)
  ! sum all contributions to the exchange energy
  call MPI_Allreduce(ex, r, 1, MPI_FLOAT, MPI_SUM, st%comm, ierr)
  ex = r

  deallocate(recv_buffer)
#endif

  deallocate(recv_stack, send_stack)
  deallocate(F_ij, rho_ij, send_buffer)

  call pop_sub()

end subroutine X(oep_x)
