!! Copyright (C) 2008 X. Andrade
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

!--------------------------------------------------------------
subroutine X(batch_init_contiguous)(this, dim, st_start, st_end, psi)
  class(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: dim
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, intent(in)    :: psi(:, :, st_start:)

  integer :: ist

  PUSH_SUB(X(batch_init_contiguous))

  ASSERT(st_end >= st_start)

  call batch_init_empty(this, dim, st_end - st_start + 1)

  this%X(psicont) => psi(:, :, st_start:)

  ASSERT(ubound(psi, dim = 3) >= st_end)

  do ist = st_start, st_end
    call this%add_state(ist, psi(:, :, ist))
  end do

  this%type_of = R_TYPE_VAL

  POP_SUB(X(batch_init_contiguous))
end subroutine X(batch_init_contiguous)

!--------------------------------------------------------------
subroutine X(batch_add_state)(this, ist, psi)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: ist
  R_TYPE, target, intent(in)    :: psi(:, :)

  integer :: idim, ii

  PUSH_SUB(X(batch_add_state))

  ASSERT(this%current <= this%nst)

  this%states(this%current)%ist    =  ist
  this%states(this%current)%X(psi) => psi

  ! now we also populate the linear array
  do idim = 1, this%dim
    ii = this%dim*(this%current - 1) + idim
    this%states_linear(ii)%X(psi) => psi(:, idim)
    this%ist_idim_index(ii, 1) = ist
    this%ist_idim_index(ii, 2) = idim
  end do

  this%max_size = max(this%max_size, ubound(this%states(this%current)%X(psi), dim = 1))
  
  this%current = this%current + 1

  POP_SUB(X(batch_add_state))
end subroutine X(batch_add_state)

!--------------------------------------------------------------

subroutine X(batch_add_state_linear)(this, psi)
  class(batch_t),  intent(inout) :: this
  R_TYPE,  target, intent(in)    :: psi(:)

  PUSH_SUB(X(batch_add_state_linear))

  ASSERT(this%current <= this%nst_linear)
  this%states_linear(this%current)%X(psi) => psi
  this%ist_idim_index(this%current, 1) = this%current

  this%max_size = max(this%max_size, ubound(this%states_linear(this%current)%X(psi), dim = 1))
  
  this%current = this%current + 1

  POP_SUB(X(batch_add_state_linear))
end subroutine X(batch_add_state_linear)


!--------------------------------------------------------------
subroutine X(batch_allocate)(this, st_start, st_end, np, mirror, special)
  class(batch_t),    intent(inout) :: this
  integer,           intent(in)    :: st_start
  integer,           intent(in)    :: st_end
  integer,           intent(in)    :: np
  logical, optional, intent(in)    :: mirror     !< If .true., this batch will keep a copy when packed. Default: .false.
  logical, optional, intent(in)    :: special    !< If .true., the allocation will be handled in C (to use pinned memory for GPUs)

  integer :: ist, nst

  PUSH_SUB(X(batch_allocate))

  nst = st_end - st_start + 1
  if(optional_default(special, .false.)) then
    call c_f_pointer(X(allocate_hardware_aware)(np*this%dim*nst), this%X(psicont), [np,this%dim,nst])
    this%special_memory = .true.
  else
    SAFE_ALLOCATE(this%X(psicont)(1:np, 1:this%dim, 1:nst))
  end if

  this%is_allocated = .true.
  this%mirror = optional_default(mirror, .false.)  
  
  do ist = st_start, st_end
    call this%add_state(ist, this%X(psicont)(:, :, ist - st_start + 1))
  end do

  POP_SUB(X(batch_allocate))
end subroutine X(batch_allocate)

!--------------------------------------------------------------
subroutine X(batch_allocate_temporary)(this)
  class(batch_t),  intent(inout) :: this

  integer :: ist, idim

  PUSH_SUB(X(batch_allocate_temporary))

  ASSERT(.not. associated(this%X(psicont)))
  
  if(this%special_memory) then
    call c_f_pointer(X(allocate_hardware_aware)(this%max_size*this%dim*this%nst), this%X(psicont), [this%max_size,this%dim,this%nst])
  else
    SAFE_ALLOCATE(this%X(psicont)(1:this%max_size, 1:this%dim, 1:this%nst))
  end if
  
  do ist = 1, this%nst
    this%states(ist)%X(psi) => this%X(psicont)(:, :, ist)
    do idim = 1, this%dim
      this%states_linear((ist - 1)*this%dim + idim)%X(psi) => this%X(psicont)(:, idim, ist)
    end do
  end do

  POP_SUB(X(batch_allocate_temporary))
end subroutine X(batch_allocate_temporary)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
