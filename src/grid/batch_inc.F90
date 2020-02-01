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
subroutine X(batch_init_with_memory_3)(this, dim, st_start, st_end, psi)
  class(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: dim
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, contiguous, intent(in)    :: psi(:, :, st_start:)

  PUSH_SUB(X(batch_init_with_memory_3))

  ASSERT(st_end >= st_start)

  call batch_init_empty(this, dim, st_end - st_start + 1, ubound(psi, dim=1))

  this%type_of = R_TYPE_VAL
  this%X(ff) => psi(:, :, st_start:)
  this%X(ff_linear)(1:this%np, 1:this%nst_linear) => this%X(ff)(:, :, :)

  ASSERT(ubound(psi, dim = 3) >= st_end)

  call X(batch_build_indices)(this, st_start, st_end)

  POP_SUB(X(batch_init_with_memory_3))
end subroutine X(batch_init_with_memory_3)

subroutine X(batch_init_with_memory_2)(this, dim, st_start, st_end, psi)
  class(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: dim
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, contiguous, intent(in)    :: psi(:, :)

  R_TYPE, pointer :: psip(:, :, :)

  PUSH_SUB(X(batch_init_with_memory_2))

  ASSERT(st_end == st_start .or. dim == 1)

  psip(1:ubound(psi, dim=1), 1:dim, st_start:st_end) => psi(:, :)

  call X(batch_init_with_memory_3)(this, dim, st_start, st_end, psip)

  POP_SUB(X(batch_init_with_memory_2))
end subroutine X(batch_init_with_memory_2)

subroutine X(batch_init_with_memory_1)(this, psi)
  class(batch_t),             intent(out)   :: this
  R_TYPE, target, contiguous, intent(in)    :: psi(:)

  R_TYPE, pointer :: psip(:, :, :)
  PUSH_SUB(X(batch_init_with_memory_1))

  psip(1:ubound(psi, dim=1), 1:1, 1:1) => psi(:)
  call X(batch_init_with_memory_3)(this, 1, 1, 1, psip)

  POP_SUB(X(batch_init_with_memory_1))
end subroutine X(batch_init_with_memory_1)

!--------------------------------------------------------------
subroutine X(batch_build_indices)(this, st_start, st_end)
  class(batch_t), intent(inout) :: this
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end

  integer :: idim, ii, ist

  PUSH_SUB(X(batch_build_indices))

  do ist = st_start, st_end
    ! now we also populate the linear array
    do idim = 1, this%dim
      ii = this%dim*(ist - st_start) + idim
      this%ist_idim_index(ii, 1) = ist
      this%ist_idim_index(ii, 2) = idim
    end do
    this%ist(ist - st_start + 1) = ist
  end do

  POP_SUB(X(batch_build_indices))
end subroutine X(batch_build_indices)


!--------------------------------------------------------------
subroutine X(batch_allocate)(this)
  class(batch_t),    intent(inout) :: this

  integer :: ist

  PUSH_SUB(X(batch_allocate))

  if(this%special_memory) then
    call c_f_pointer(X(allocate_hardware_aware)(this%np*this%dim*this%nst), this%X(ff), &
      [this%np,this%dim,this%nst])
  else
    SAFE_ALLOCATE(this%X(ff)(1:this%np, 1:this%dim, 1:this%nst))
  end if
  this%type_of = R_TYPE_VAL
  this%X(ff_linear)(1:this%np, 1:this%nst_linear) => this%X(ff)(:, :, :)

  this%is_allocated = .true.
  this%own_memory = .true.
  
  POP_SUB(X(batch_allocate))
end subroutine X(batch_allocate)

subroutine X(batch_init)(this, dim, st_start, st_end, np, mirror, special)
  class(batch_t),    intent(inout) :: this
  integer,           intent(in)    :: dim
  integer,           intent(in)    :: st_start
  integer,           intent(in)    :: st_end
  integer,           intent(in)    :: np
  logical, optional, intent(in)    :: mirror     !< If .true., this batch will keep a copy when packed. Default: .false.
  logical, optional, intent(in)    :: special    !< If .true., the allocation will be handled in C (to use pinned memory for GPUs)

  PUSH_SUB(X(batch_init))

  call batch_init_empty(this, dim, st_end - st_start + 1, np)
  this%mirror = optional_default(mirror, .false.)
  this%special_memory = optional_default(special, .false.)

  call this%X(allocate)()
  call X(batch_build_indices)(this, st_start, st_end)

  POP_SUB(X(batch_init))
end subroutine X(batch_init)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
