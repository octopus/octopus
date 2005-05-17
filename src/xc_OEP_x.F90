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

subroutine X(oep_x) (m, f_der, st, is, oep, ex)
  type(mesh_type),   intent(IN)    :: m
  type(f_der_type),  intent(in)    :: f_der
  type(states_type), intent(IN)    :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,             intent(out)   :: ex

  integer :: i, j, ist, jst, i_max, node_to, node_fr, ist_s, ist_r
  integer, allocatable :: recv_stack(:), send_stack(:)
  FLOAT :: r
  R_TYPE, pointer     :: wf_ist(:), send_buffer(:)
  R_TYPE, allocatable :: rho_ij(:), F_ij(:)
#if defined(HAVE_MPI)
  R_TYPE, pointer :: recv_buffer(:)
  integer :: ierr, req, status
#endif

  ! Note: we assume that st%occ is known in all nodes

  call push_sub('oep_x')

  allocate(F_ij(m%np), rho_ij(m%np), send_buffer(m%np))

#if defined(HAVE_MPI)
  allocate(recv_buffer(m%np))
#endif

  i_max = int((mpiv%numprocs + 2)/2) - 1
  allocate(recv_stack(st%nst+1), send_stack(st%nst+1))
  do i = 0, i_max
    ! node where to send the wavefunctions
    node_to = mod(mpiv%node + i, mpiv%numprocs)
    if(i==i_max .and. mod(mpiv%numprocs,2)==0 .and. node_to < mpiv%numprocs/2) then
      node_to = -1
    end if

    ! node from which we receive the wave-functions
    node_fr = mpiv%node - i
    if(node_fr < 0) node_fr = mpiv%numprocs + node_fr
    node_fr = mod(node_fr, mpiv%numprocs)
    if(i==i_max .and. mod(mpiv%numprocs, 2)==0 .and. mpiv%node < mpiv%numprocs/2) then
      node_fr = -1
    end if

    ! build receive and send stacks
    recv_stack(:) = -1; ist_r = 1
    send_stack(:) = -1; ist_s = 1
    do j = 1, st%nst
      if(st%node(j) == node_fr) then
        recv_stack(ist_r) = j
        ist_r = ist_r + 1
      end if
      if(st%node(j) == mpiv%node) then
        send_stack(ist_s) = j
        ist_s = ist_s + 1
      end if
    end do

    ! now we do a loop over all states
    ist_r = 0
    if(node_fr <0) ist_r = st%nst + 1
    ist_s = 0
    if(node_to < 0) ist_s = st%nst + 1
    do
      ! send wave-function
      if(ist_s <= st%nst) ist_s = ist_s  + 1

#if defined(HAVE_MPI)          
      if((send_stack(ist_s) > 0).and.(node_to.ne.mpiv%node)) then
        call MPI_ISend(st%X(psi)(:, 1, send_stack(ist_s), is), m%np, R_MPITYPE, &
            node_to, send_stack(ist_s), MPI_COMM_WORLD, req, ierr)
      end if
#endif

      ! receive wave-function
      if(ist_r <= st%nst) ist_r = ist_r  + 1

      if(recv_stack(ist_r) > 0) then
        if(node_fr == mpiv%node) then
          wf_ist => st%X(psi)(:, 1, recv_stack(ist_r), is)
#if defined(HAVE_MPI)
        else
          call MPI_Recv(recv_buffer(:), m%np, R_MPITYPE, &
              node_fr, recv_stack(ist_r), MPI_COMM_WORLD, status, ierr)
          wf_ist => recv_buffer(:)
#endif
        end if

        ! calculate stuff
        ist = recv_stack(ist_r)
        send_buffer(:) = R_TOTYPE(M_ZERO)
        do jst = st%st_start, st%st_end
          if((st%node(ist) == mpiv%node).and.(jst < ist)) cycle
          if((st%occ(ist, is).le.small).or.(st%occ(jst, is).le.small)) cycle

          rho_ij(:) = R_CONJ(wf_ist(:))*st%X(psi)(:, 1, jst, is)
          F_ij = R_TOTYPE(M_ZERO)
          call X(poisson_solve)(m, f_der, F_ij, rho_ij)

          send_buffer(:) = send_buffer(:) + &
              oep%socc*st%occ(jst, is)*F_ij(:)*R_CONJ(st%X(psi)(:, 1, jst, is))

          if(ist.ne.jst) then
            oep%X(lxc)(:, jst) = oep%X(lxc)(:, jst) - &
                oep%socc * st%occ(ist, is) * R_CONJ(F_ij(:)*wf_ist(:))
          end if

          r = M_ONE
          if(ist.ne.jst) r = M_TWO
          ex = ex - M_HALF * r * oep%sfact * oep%socc*st%occ(ist, is) * oep%socc*st%occ(jst, is) * &
              sum(R_REAL(wf_ist(:) * F_ij(:) * R_CONJ(st%X(psi)(:, 1, jst, is))) * m%vol_pp(:))
        end do

        if(st%node(ist) == mpiv%node) then
          oep%X(lxc)(:, ist) = oep%X(lxc)(:, ist) - send_buffer(:)
#if defined(HAVE_MPI)
        else
          call MPI_ISend(send_buffer(:), m%np, R_MPITYPE, &
              node_fr, 999+ist, MPI_COMM_WORLD, req, ierr)
#endif
        end if
        
      end if

#if defined(HAVE_MPI)
      if((send_stack(ist_s) > 0).and.(node_to.ne.mpiv%node)) then
        call MPI_Recv(send_buffer(:), m%np, R_MPITYPE, &
            node_to, 999+send_stack(ist_s), MPI_COMM_WORLD, status, ierr)

        oep%X(lxc)(:, send_stack(ist_s)) = oep%X(lxc)(:, send_stack(ist_s)) - &
            send_buffer(:)
      end if
#endif

      if((send_stack(ist_s) < 0).and.(recv_stack(ist_r) < 0)) exit
    end do
  end do

#if defined(HAVE_MPI)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(ex, r, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
  ex = r

  deallocate(recv_buffer)
#endif

  deallocate(recv_stack, send_stack)
  deallocate(F_ij, rho_ij, send_buffer)

  call pop_sub()

end subroutine X(oep_x)
