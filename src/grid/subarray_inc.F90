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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: subarray_inc.F90 3030 2007-06-25 16:45:05Z marques $

! ---------------------------------------------------

subroutine X(subarray_gather)(this, array, subarray)
    type(subarray_t),    intent(in)  :: this
    R_TYPE,              intent(in)  :: array(:)
    R_TYPE,              intent(out) :: subarray(:)

    type(profile_t), save :: prof
    integer :: iblock, ii, isa

    call profiling_in(prof, "SUBARRAY_GATHER")

    isa = 0
    do iblock = 1, this%nblocks
      forall(ii = 1:this%blength(iblock)) subarray(isa + ii) = array(this%offsets(iblock) + ii - 1)
      isa = isa + this%blength(iblock)
    end do

    call profiling_count_transfers(this%npoints, array(1))

    call profiling_out(prof)
end subroutine X(subarray_gather)

#if defined(R_TREAL) || defined(R_TCOMPLEX)
subroutine X(subarray_gather_batch)(this, arrayb, subarrayb)
    type(subarray_t),    intent(in)    :: this
    type(batch_t),       intent(in)    :: arrayb
    type(batch_t),       intent(inout) :: subarrayb

    type(profile_t), save :: prof
    integer :: iblock, ii, isa, ist
    R_TYPE  :: aa
#ifdef HAVE_OPENCL
    type(opencl_mem_t) :: blength_buff
    type(opencl_mem_t) :: offsets_buff
#endif

    call profiling_in(prof, "SUBARRAY_GATHER_BATCH")


    ASSERT(batch_status(arrayb) == batch_status(subarrayb))
    
    select case(batch_status(arrayb))
#ifdef HAVE_OPENCL
    case(BATCH_CL_PACKED)
      call opencl_create_buffer(blength_buff, CL_MEM_READ_ONLY, TYPE_INTEGER, this%nblocks)
      call opencl_create_buffer(offsets_buff, CL_MEM_READ_ONLY, TYPE_INTEGER, this%nblocks)
      call opencl_write_buffer(blength_buff, this%nblocks, this%blength)
      call opencl_write_buffer(offsets_buff, this%nblocks, this%offsets)
      
      call opencl_set_kernel_arg(kernel_subarray_gather, 0, this%nblocks)
      call opencl_set_kernel_arg(kernel_subarray_gather, 1, blength_buff)
      call opencl_set_kernel_arg(kernel_subarray_gather, 2, offsets_buff)
      call opencl_set_kernel_arg(kernel_subarray_gather, 3, arrayb%pack%buffer)
      call opencl_set_kernel_arg(kernel_subarray_gather, 4, log2(arrayb%pack%size_real(1)))
      call opencl_set_kernel_arg(kernel_subarray_gather, 5, subarrayb%pack%buffer)
      call opencl_set_kernel_arg(kernel_subarray_gather, 6, log2(subarrayb%pack%size_real(1)))
      
      call opencl_kernel_run(kernel_subarray_gather, (/subarrayb%pack%size_real(1)/), (/subarrayb%pack%size_real(1)/))
      
      call opencl_release_buffer(blength_buff)
      call opencl_release_buffer(offsets_buff)
#endif
    case(BATCH_PACKED)
      isa = 0
      do iblock = 1, this%nblocks
        forall(ii = 1:this%blength(iblock))
          forall(ist = 1:arrayb%pack%size(1))
            subarrayb%pack%X(psi)(ist, isa + ii) = arrayb%pack%X(psi)(ist, this%offsets(iblock) + ii - 1)
          end forall
        end forall
        isa = isa + this%blength(iblock)
      end do
      
    case(BATCH_NOT_PACKED)
      !$omp parallel do private(isa, iblock, ii)
      do ist = 1, arrayb%nst_linear
        isa = 0
        do iblock = 1, this%nblocks
          forall(ii = 1:this%blength(iblock))
            subarrayb%states_linear(ist)%X(psi)(isa + ii) = &
              arrayb%states_linear(ist)%X(psi)(this%offsets(iblock) + ii - 1)
          end forall
          isa = isa + this%blength(iblock)
        end do
      end do

    end select

    call profiling_count_transfers(arrayb%nst_linear*this%npoints, aa)

    call profiling_out(prof)
end subroutine X(subarray_gather_batch)
#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
