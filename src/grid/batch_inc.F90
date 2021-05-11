!! Copyright (C) 2008 X. Andrade, 2020 S. Ohlmann
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
  this%X(ff_linear)(1:this%np, 1:this%nst_linear) => this%X(ff)

  ASSERT(ubound(psi, dim=3) >= st_end)
  ASSERT(ubound(psi, dim=2) == dim)

  call batch_build_indices(this, st_start, st_end)

  POP_SUB(X(batch_init_with_memory_3))
end subroutine X(batch_init_with_memory_3)

subroutine X(batch_init_with_memory_2)(this, dim, st_start, st_end, psi)
  class(batch_t),  intent(out)   :: this
  integer,        intent(in)    :: dim
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  R_TYPE, target, contiguous, intent(in)    :: psi(:, :)

  R_TYPE, pointer, contiguous :: psip(:, :, :)

  PUSH_SUB(X(batch_init_with_memory_2))

  ASSERT(st_end == st_start .or. dim == 1)

  psip(1:ubound(psi, dim=1), 1:dim, st_start:st_end) => psi(:, :)

  call X(batch_init_with_memory_3)(this, dim, st_start, st_end, psip)

  POP_SUB(X(batch_init_with_memory_2))
end subroutine X(batch_init_with_memory_2)

subroutine X(batch_init_with_memory_1)(this, psi)
  class(batch_t),             intent(out)   :: this
  R_TYPE, target, contiguous, intent(in)    :: psi(:)

  R_TYPE, pointer, contiguous :: psip(:, :, :)
  PUSH_SUB(X(batch_init_with_memory_1))

  psip(1:ubound(psi, dim=1), 1:1, 1:1) => psi(:)
  call X(batch_init_with_memory_3)(this, 1, 1, 1, psip)

  POP_SUB(X(batch_init_with_memory_1))
end subroutine X(batch_init_with_memory_1)


!--------------------------------------------------------------
subroutine X(batch_allocate_unpacked_host)(this)
  class(batch_t),    intent(inout) :: this

  PUSH_SUB(X(batch_allocate_unpacked_host))

  if(this%special_memory) then
    call c_f_pointer(X(allocate_hardware_aware)(this%np*this%dim*this%nst), this%X(ff), &
      [this%np,this%dim,this%nst])
  else
    SAFE_ALLOCATE(this%X(ff)(1:this%np, 1:this%dim, 1:this%nst))
  end if
  this%X(ff_linear)(1:this%np, 1:this%nst_linear) => this%X(ff)

  this%is_allocated = .true.

  POP_SUB(X(batch_allocate_unpacked_host))
end subroutine X(batch_allocate_unpacked_host)

!--------------------------------------------------------------
subroutine X(batch_allocate_packed_host)(this)
  class(batch_t),    intent(inout) :: this

  PUSH_SUB(X(batch_allocate_packed_host))

  if(this%special_memory) then
    call c_f_pointer(X(allocate_hardware_aware)(this%pack_size(1)*this%pack_size(2)), this%X(ff_pack), this%pack_size)
  else
    SAFE_ALLOCATE(this%X(ff_pack)(1:this%pack_size(1), 1:this%pack_size(2)))
  end if

  POP_SUB(X(batch_allocate_packed_host))
end subroutine X(batch_allocate_packed_host)

!--------------------------------------------------------------
subroutine X(batch_init)(this, dim, st_start, st_end, np, special, packed)
  class(batch_t),    intent(inout) :: this
  integer,           intent(in)    :: dim
  integer,           intent(in)    :: st_start
  integer,           intent(in)    :: st_end
  integer,           intent(in)    :: np
  logical, optional, intent(in)    :: special    !< If .true., the allocation will be handled in C (to use pinned memory for GPUs)
  logical, optional, intent(in)    :: packed    !< If .true., the allocation will be handled in C (to use pinned memory for GPUs)

  PUSH_SUB(X(batch_init))

  call batch_init_empty(this, dim, st_end - st_start + 1, np)
  this%special_memory = optional_default(special, .false.)
  this%type_of = R_TYPE_VAL
  call batch_build_indices(this, st_start, st_end)

  if(optional_default(packed, .false.)) then
    call this%X(allocate_packed_host)()
    this%status_of = BATCH_PACKED
    this%status_host = BATCH_PACKED
    this%host_buffer_count = this%host_buffer_count + 1
  else
    call this%X(allocate_unpacked_host)()
  end if

  this%own_memory = .true.

  POP_SUB(X(batch_init))
end subroutine X(batch_init)

!--------------------------------------------------------------
subroutine X(batch_pack_copy)(this)
  class(batch_t),    intent(inout) :: this

  integer :: ist, ip, sp, ep, bsize
  type(profile_t), save :: prof_copy
  ! no push_sub, called too frequently

  call profiling_in(prof_copy, TOSTRING(X(BATCH_PACK_COPY)))

  bsize = hardware%X(block_size)

  !$omp parallel do private(ep, ist, ip)
  do sp = 1, this%pack_size(2), bsize
    ep = min(sp + bsize - 1, this%pack_size(2))
    do ist = 1, this%nst_linear
      do ip = sp, ep
        this%X(ff_pack)(ist, ip) = this%X(ff_linear)(ip, ist)
      end do
    end do
  end do

  call profiling_count_transfers(this%nst_linear*this%pack_size(2), this%type())
  call profiling_out(prof_copy)
end subroutine X(batch_pack_copy)

!--------------------------------------------------------------
subroutine X(batch_unpack_copy)(this)
  class(batch_t),    intent(inout) :: this

  integer :: ist, ip
  type(profile_t), save :: prof_copy
  ! no push_sub, called too frequently

  call profiling_in(prof_copy, TOSTRING(X(BATCH_UNPACK_COPY)))

  !$omp parallel do private(ist)
  do ip = 1, this%pack_size(2)
    do ist = 1, this%nst_linear
      this%X(ff_linear)(ip, ist) = this%X(ff_pack)(ist, ip)
    end do
  end do
  call profiling_count_transfers(this%nst_linear*this%pack_size(2), this%type())
  call profiling_out(prof_copy)
end subroutine X(batch_unpack_copy)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
