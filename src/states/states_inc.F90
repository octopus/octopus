!! Copyright (C) 2011 X. Andrade
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
!! $Id: states_inc.F90 7770 2011-04-26 00:51:17Z dstrubbe $

subroutine X(states_get_state2)(st, mesh, ist, iqn, psi, left)
  type(states_t),    intent(in)    :: st
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: ist       !< current state
  integer,           intent(in)    :: iqn       !< current k-point
  R_TYPE,            intent(out)   :: psi(:, :)
  logical, optional, intent (in)   :: left
  
  integer :: idim

  PUSH_SUB(X(states_get_state2))

  do idim =  1, st%d%dim
    call X(states_get_state1)(st, mesh, idim, ist, iqn, psi(:, idim), left)
  end do

  POP_SUB(X(states_get_state2))
end subroutine X(states_get_state2)

! ------------------------------------------------------------

subroutine X(states_get_state1)(st, mesh, idim, ist, iqn, psi, left)
  type(states_t),    intent(in)    :: st
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: idim   !< current dimension
  integer,           intent(in)    :: ist    !< current k-point
  integer,           intent(in)    :: iqn
  R_TYPE,            intent(out)   :: psi(:)
  logical, optional, intent(in)    :: left 

  PUSH_SUB(X(states_get_state1))

  if (optional_default(left, .false.)) then
    ASSERT(st%have_left_states)
    call batch_get_state(st%psibL(st%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)
  else
    call batch_get_state(st%psib(st%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)
  end if


  POP_SUB(X(states_get_state1))
end subroutine X(states_get_state1)

! ------------------------------------------------------------

subroutine X(states_set_state2)(st, mesh, ist, iqn, psi, left)
  type(states_t),    intent(inout) :: st
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: ist       !< current dimension
  integer,           intent(in)    :: iqn       !< current k-point
  R_TYPE,            intent(in)    :: psi(:, :)
  logical, optional, intent(in)    :: left
  
  integer :: idim

  PUSH_SUB(X(states_set_state2))

  do idim =  1, st%d%dim
    call X(states_set_state1)(st, mesh, idim, ist, iqn, psi(:, idim), left)
  end do

  POP_SUB(X(states_set_state2))
end subroutine X(states_set_state2)

! ------------------------------------------------------------

subroutine X(states_set_state1)(st, mesh, idim, ist, iqn, psi, left)
  type(states_t),    intent(inout) :: st
  type(mesh_t),      intent(in)    :: mesh
  integer,           intent(in)    :: idim   !< current dimension
  integer,           intent(in)    :: ist    !< current state
  integer,           intent(in)    :: iqn    !< current k-point
  R_TYPE,            intent(in)    :: psi(:)
  logical, optional, intent(in)    :: left

  PUSH_SUB(X(states_set_state1))
  
  if (optional_default(left, .false.)) then
    ASSERT(st%have_left_states)
    call batch_set_state(st%psibL(st%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)
  else
    call batch_set_state(st%psib(st%iblock(ist, iqn), iqn), (/ist, idim/), mesh%np, psi)
  end if
  
  POP_SUB(X(states_set_state1))
end subroutine X(states_set_state1)

! ------------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
