!! Copyright (C) 2009 X. Andrade
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

! ---------------------------------------------------

subroutine X(subarray_gather)(this, array, subarray)
  type(subarray_t),    intent(in)  :: this
  R_TYPE,              intent(in)  :: array(:)
  R_TYPE,              intent(out) :: subarray(:)

  type(profile_t), save :: prof
  integer :: iblock, ii

  call profiling_in(prof, "SUBARRAY_GATHER")

  do iblock = 1, this%nblocks
    forall(ii = 1:this%blength(iblock)) subarray(this%dest(iblock) + ii) = array(this%offsets(iblock) + ii - 1)
  end do

  call profiling_count_transfers(this%npoints, array(1))

  call profiling_out(prof)
end subroutine X(subarray_gather)

! ---------------------------------------------------

#if defined(R_TREAL) || defined(R_TCOMPLEX)
subroutine X(subarray_gather_batch)(this, arrayb, subarrayb)
  type(subarray_t),    intent(in)    :: this
  type(batch_t),       intent(in)    :: arrayb
  type(batch_t),       intent(inout) :: subarrayb

  type(profile_t), save :: prof
  integer :: iblock, ii, ist, bsize
  R_TYPE  :: aa
  type(accel_mem_t) :: blength_buff
  type(accel_mem_t) :: offsets_buff
  type(accel_mem_t) :: dest_buff

  PUSH_SUB(X(subarray_gather_batch))

  call profiling_in(prof, "SUBARRAY_GATHER_BATCH")


  ASSERT(arrayb%status() == subarrayb%status())
  call arrayb%check_compatibility_with(subarrayb)
    
  select case(arrayb%status())
  case(BATCH_DEVICE_PACKED)

    call accel_create_buffer(blength_buff, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nblocks)
    call accel_create_buffer(offsets_buff, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nblocks)
    call accel_create_buffer(dest_buff, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nblocks)

    call accel_write_buffer(blength_buff, this%nblocks, this%blength)
    call accel_write_buffer(offsets_buff, this%nblocks, this%offsets)
    call accel_write_buffer(dest_buff, this%nblocks, this%dest)
    
    call accel_set_kernel_arg(kernel_subarray_gather, 0, blength_buff)
    call accel_set_kernel_arg(kernel_subarray_gather, 1, offsets_buff)
    call accel_set_kernel_arg(kernel_subarray_gather, 2, dest_buff)
    call accel_set_kernel_arg(kernel_subarray_gather, 3, arrayb%pack%buffer)
    call accel_set_kernel_arg(kernel_subarray_gather, 4, log2(arrayb%pack%size_real(1)))
    call accel_set_kernel_arg(kernel_subarray_gather, 5, subarrayb%pack%buffer)
    call accel_set_kernel_arg(kernel_subarray_gather, 6, log2(subarrayb%pack%size_real(1)))

    bsize = accel_kernel_workgroup_size(kernel_subarray_gather)/subarrayb%pack%size_real(1)

    call accel_kernel_run(kernel_subarray_gather, &
      (/subarrayb%pack%size_real(1), bsize, this%nblocks/), (/subarrayb%pack%size_real(1), bsize, 1/))

    call accel_finish()
    
    call accel_release_buffer(blength_buff)
    call accel_release_buffer(offsets_buff)
    call accel_release_buffer(dest_buff)
    
  case(BATCH_PACKED)
    do iblock = 1, this%nblocks
      forall(ii = 1:this%blength(iblock))
        forall(ist = 1:arrayb%pack%size(1))
          subarrayb%X(ff_pack)(ist, this%dest(iblock) + ii) = arrayb%X(ff_pack)(ist, this%offsets(iblock) + ii - 1)
        end forall
      end forall
    end do
    
  case(BATCH_NOT_PACKED)
    !$omp parallel do private(iblock, ii)
    do ist = 1, arrayb%nst_linear
      do iblock = 1, this%nblocks
        forall(ii = 1:this%blength(iblock))
          subarrayb%X(ff_linear)(this%dest(iblock) + ii, ist) = &
            arrayb%X(ff_linear)(this%offsets(iblock) + ii - 1, ist)
        end forall
      end do
    end do
    
  end select

  ! Avoid warning: 'aa' is used uninitialized; it is just to define which type
  aa = R_TOTYPE(M_ZERO)
  call profiling_count_transfers(arrayb%nst_linear*this%npoints, aa)

  call profiling_out(prof)
  POP_SUB(X(subarray_gather_batch))
end subroutine X(subarray_gather_batch)
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
