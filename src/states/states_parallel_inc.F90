!! Copyright (C) 2015 X. Andrade
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

subroutine X(states_parallel_gather)(st, dims, psi)
  type(states_t), intent(in)    :: st
  integer,        intent(in)    :: dims(2)
  R_TYPE,         intent(inout) :: psi(:, :, :)

  !no PUSH_SUB, called too often
  
  call profiling_in(prof_gather, 'STATES_GATHER')

  if(st%parallel_in_states) then
    !this should really be an allgather, we use the simpler allreduce
    !for the moment to get it working
    psi(1:st%st_start - 1, 1:dims(1), 1:dims(2)) = CNST(0.0)
    psi(st%st_end + 1:st%nst, 1:dims(1), 1:dims(2)) = CNST(0.0)
    call comm_allreduce(st%mpi_grp%comm, psi)
  end if

  call profiling_out(prof_gather)
end subroutine X(states_parallel_gather)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
