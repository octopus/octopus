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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: batch.F90 4298 2008-06-18 15:03:00Z dstrubbe $

#include "global.h"

module batch_m
  use blas_m
  use c_pointer_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use octcl_kernel_m
  use datasets_m
  use global_m
  use hardware_m
  use lalg_adv_m
  use lalg_basic_m
  use parser_m
  use math_m
  use messages_m
  use opencl_m
  use profiling_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::                         &
    batch_state_t,                  &
    batch_state_l_t,                &
    batch_pack_t,                   &
    batch_t,                        &
    batch_init,                     &
    batch_copy,                     &
    batch_end,                      &
    batch_add_state,                &
    dbatch_new,                     &
    zbatch_new,                     &
    batch_delete,                   &
    batch_set,                      &
    batch_set_zero,                 &
    batch_is_packed,                &
    batch_pack_was_modified,        &
    batch_is_ok,                    &
    batch_axpy,                     &
    batch_copy_data,                &
    batch_pack,                     &
    batch_unpack,                   &
    batch_sync,                     &
    batch_status,                   &
    batch_type,                     &
    batch_inv_index,                &
    batch_linear_index,             &
    batch_set_state,                &
    batch_get_state,                &
    batch_get_points,               &
    batch_set_points,               &
    batch_pack_size

  !--------------------------------------------------------------
  type batch_state_t
    FLOAT, pointer :: dpsi(:, :)
    CMPLX, pointer :: zpsi(:, :)
    integer        :: ist
  end type batch_state_t

  type batch_state_l_t
    FLOAT, pointer :: dpsi(:)
    CMPLX, pointer :: zpsi(:)
  end type batch_state_l_t
  
  type batch_pack_t
    integer                        :: size(1:2)
    integer                        :: size_real(1:2)
    FLOAT, pointer                 :: dpsi(:, :)
    CMPLX, pointer                 :: zpsi(:, :)
    type(opencl_mem_t)             :: buffer
  end type batch_pack_t
  
  type batch_t
    type(batch_state_t), pointer   :: states(:)
    integer                        :: nst
    integer                        :: current
    integer                        :: dim

    integer                        :: ndims
    integer,             pointer   :: index(:, :)

    logical                        :: is_allocated

    !> We also need a linear array with the states in order to calculate derivatives, etc.
    integer                        :: nst_linear
    type(batch_state_l_t), pointer :: states_linear(:)

    !> If the memory is contiguous, we can perform some operations faster.
    FLOAT,               pointer   :: dpsicont(:, :, :)
    CMPLX,               pointer   :: zpsicont(:, :, :)
    integer                        :: status
    integer                        :: in_buffer_count !< whether there is a copy in the opencl buffer
    logical                        :: dirty     !< if this is true, the buffer has different data
    type(batch_pack_t)             :: pack
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure  batch_init_empty
    module procedure  batch_init_empty_linear
    module procedure dbatch_init_contiguous
    module procedure zbatch_init_contiguous
  end interface

  interface batch_add_state
    module procedure dbatch_add_state
    module procedure zbatch_add_state
    module procedure dbatch_add_state_linear
    module procedure zbatch_add_state_linear
  end interface

  interface batch_set
    module procedure dbatch_set
    module procedure zbatch_set
  end interface

  interface batch_axpy
    module procedure dbatch_axpy_const
    module procedure zbatch_axpy_const
    module procedure dbatch_axpy_vec
    module procedure zbatch_axpy_vec
  end interface

  interface batch_set_state
    module procedure dbatch_set_state1
    module procedure zbatch_set_state1
    module procedure dbatch_set_state2
    module procedure zbatch_set_state2
  end interface

  interface batch_get_state
    module procedure dbatch_get_state1
    module procedure zbatch_get_state1
    module procedure dbatch_get_state2
    module procedure zbatch_get_state2
  end interface

  interface batch_get_points
    module procedure dbatch_get_points
    module procedure zbatch_get_points
    module procedure batch_get_points_cl
  end interface batch_get_points

  interface batch_set_points
    module procedure dbatch_set_points
    module procedure zbatch_set_points
    module procedure batch_set_points_cl
  end interface batch_set_points

  type(profile_t), save :: axpy_prof, get_points_prof, set_points_prof

  integer, public, parameter :: &
    BATCH_NOT_PACKED = 0,       &
    BATCH_PACKED     = 1,       &
    BATCH_CL_PACKED  = 2

  integer, parameter :: CL_PACK_MAX_BUFFER_SIZE = 4 !< this value controls the size (in number of wave-functions)
                                                    !! of the buffer used to copy states to the opencl device.

contains

  !--------------------------------------------------------------
  subroutine batch_end(this, copy)
    type(batch_t),           intent(inout) :: this
    logical,       optional, intent(in)    :: copy

    PUSH_SUB(batch_end)

    if(batch_is_packed(this)) then
      call batch_unpack(this, copy)
    end if

    if(this%is_allocated) then
      call batch_delete(this)
    end if

    nullify(this%dpsicont)
    nullify(this%zpsicont)

    SAFE_DEALLOCATE_P(this%states)
    SAFE_DEALLOCATE_P(this%states_linear)
    
    SAFE_DEALLOCATE_P(this%index)

    POP_SUB(batch_end)
  end subroutine batch_end

  !--------------------------------------------------------------
  subroutine batch_delete(this)
    type(batch_t),  intent(inout) :: this
    
    integer :: ii
    
    PUSH_SUB(batch_delete)

    this%is_allocated = .false.

    do ii = 1, this%nst
      nullify(this%states(ii)%dpsi)
      nullify(this%states(ii)%zpsi)
    end do
    
    do ii = 1, this%nst_linear
      nullify(this%states_linear(ii)%dpsi)
      nullify(this%states_linear(ii)%zpsi)
    end do
    
    this%current = 1
    
    SAFE_DEALLOCATE_P(this%dpsicont)
    SAFE_DEALLOCATE_P(this%zpsicont)
    
    POP_SUB(batch_delete)
  end subroutine batch_delete
  
  !--------------------------------------------------------------

  subroutine batch_init_empty (this, dim, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nst
    
    integer :: ist

    PUSH_SUB(batch_init_empty)

    this%is_allocated = .false.
    this%nst = nst
    this%dim = dim
    this%current = 1
    nullify(this%dpsicont, this%zpsicont)
    
    SAFE_ALLOCATE(this%states(1:nst))
    do ist = 1, nst
      nullify(this%states(ist)%dpsi)
      nullify(this%states(ist)%zpsi)
    end do
    
    this%nst_linear = nst*dim
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
    end do
    
    this%in_buffer_count = 0
    this%status = BATCH_NOT_PACKED

    nullify(this%pack%dpsi)
    nullify(this%pack%zpsi)

    this%ndims = 2
    SAFE_ALLOCATE(this%index(1:this%nst_linear, 1:this%ndims))

    POP_SUB(batch_init_empty)
  end subroutine batch_init_empty


  !--------------------------------------------------------------
  !> When we are interested in batches of 1D functions
  subroutine batch_init_empty_linear(this, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: nst
    
    integer :: ist

    PUSH_SUB(batch_init_empty_linear)

    this%is_allocated = .false.
    this%nst = 0
    this%dim = 0
    this%current = 1
    nullify(this%dpsicont, this%zpsicont)
    nullify(this%states)

    this%nst_linear = nst
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
    end do

    this%in_buffer_count = 0
    this%status = BATCH_NOT_PACKED

    nullify(this%pack%dpsi)
    nullify(this%pack%zpsi)

    this%ndims = 1
    SAFE_ALLOCATE(this%index(1:this%nst_linear, 1:this%ndims))

    POP_SUB(batch_init_empty_linear)
  end subroutine batch_init_empty_linear


  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    type(batch_t), intent(in)   :: this

    integer :: ist
    
    PUSH_SUB(batch_is_ok)
    
    ok = (this%nst_linear >= 1) .and. associated(this%states_linear)
    ok = ubound(this%states_linear, dim = 1) == this%nst_linear
    if(ok) then
      do ist = 1, this%nst_linear
        ok = ok.and. &
          (associated(this%states_linear(ist)%dpsi) .or. associated(this%states_linear(ist)%zpsi))
      end do
    end if

    POP_SUB(batch_is_ok)
  end function batch_is_ok

  !--------------------------------------------------------------

  subroutine batch_copy(bin, bout, reference)
    type(batch_t),           intent(in)    :: bin
    type(batch_t),           intent(out)   :: bout
    logical,       optional, intent(in)    :: reference !< If .true. (the default) the copy points to the same memory.
                                                        !! If .false. new memory is allocated.

    integer :: ii, np

    PUSH_SUB(batch_copy)

    call batch_init_empty(bout, bin%dim, bin%nst)

    if(optional_default(reference, .true.)) then
      
      do ii = 1, bout%nst
        if(associated(bin%states(ii)%dpsi)) bout%states(ii)%dpsi => bin%states(ii)%dpsi
        if(associated(bin%states(ii)%zpsi)) bout%states(ii)%zpsi => bin%states(ii)%zpsi
      end do
      
      do ii = 1, bout%nst_linear
        if(associated(bin%states_linear(ii)%dpsi)) bout%states_linear(ii)%dpsi => bin%states_linear(ii)%dpsi
        if(associated(bin%states_linear(ii)%zpsi)) bout%states_linear(ii)%zpsi => bin%states_linear(ii)%zpsi
      end do

      bout%current             = bin%current
      bout%status              = bin%status
      bout%in_buffer_count     = bin%in_buffer_count
      bout%dirty               = bin%dirty
      bout%pack%size(1:2)      = bin%pack%size(1:2)
      bout%pack%size_real(1:2) = bin%pack%size_real(1:2)

      if(opencl_is_enabled()) then
        bout%pack%buffer = bin%pack%buffer
      end if

    else      

      if(batch_type(bin) == TYPE_FLOAT) then
        np = 0
        do ii = 1, bin%nst_linear
          np = max(np, ubound(bin%states_linear(ii)%dpsi, dim = 1))
        end do

        call dbatch_new(bout, 1, bin%nst, np)
      else
        np = 0
        do ii = 1, bin%nst_linear
          np = max(np, ubound(bin%states_linear(ii)%zpsi, dim = 1))
        end do

        call zbatch_new(bout, 1, bin%nst, np)
      end if

!      if(batch_is_packed(bin)) call batch_pack(bout, copy = .false.)

    end if

    do ii = 1, bout%nst
      bout%states(ii)%ist = bin%states(ii)%ist
    end do

    bout%index(1:bin%nst, 1:bin%ndims) = bin%index(1:bin%nst, 1:bin%ndims)

    POP_SUB(batch_copy)
  end subroutine batch_copy

  ! ----------------------------------------------------
  !> THREADSAFE
  type(type_t) pure function batch_type(this) result(btype)
    type(batch_t),      intent(in)    :: this

    if(associated(this%states_linear(1)%dpsi)) btype = TYPE_FLOAT
    if(associated(this%states_linear(1)%zpsi)) btype = TYPE_CMPLX
  end function batch_type

  ! ----------------------------------------------------
  !> THREADSAFE
  integer pure function batch_status(this) result(bstatus)
    type(batch_t),      intent(in)    :: this

    bstatus = this%status
  end function batch_status
  
  ! ----------------------------------------------------

  logical pure function batch_is_packed(this) result(in_buffer)
    type(batch_t),      intent(in)    :: this

    in_buffer = this%in_buffer_count > 0
  end function batch_is_packed

  ! ----------------------------------------------------

  integer pure function batch_max_size(this) result(size)
    type(batch_t),      intent(in)    :: this

    integer :: ist

    size = 0
    do ist = 1, this%nst_linear
      if(associated(this%states_linear(ist)%dpsi)) then
        size = max(size, ubound(this%states_linear(ist)%dpsi, dim = 1))
      else
        size = max(size, ubound(this%states_linear(ist)%zpsi, dim = 1))
      end if
    end do

  end function batch_max_size

  ! ----------------------------------------------------

  subroutine batch_pack_was_modified(this)
    type(batch_t),      intent(inout) :: this

    this%dirty = .true.
  end subroutine batch_pack_was_modified

  ! ----------------------------------------------------

  integer function batch_pack_size(this) result(size)
    type(batch_t),      intent(inout) :: this

    size = batch_max_size(this)
    if(opencl_is_enabled()) size = opencl_padded_size(size)
    size = size*pad_pow2(this%nst_linear)*types_get_size(batch_type(this))

  end function batch_pack_size

  ! ----------------------------------------------------

  subroutine batch_pack(this, copy)
    type(batch_t),      intent(inout) :: this
    logical, optional,  intent(in)    :: copy

    logical :: copy_
    type(profile_t), save :: prof, prof_copy

    PUSH_SUB(batch_pack)

    call profiling_in(prof, "BATCH_PACK")
    ASSERT(batch_is_ok(this))

    copy_ = .true.
    if(present(copy)) copy_ = copy

    if(.not. batch_is_packed(this)) then
      this%pack%size(1) = pad_pow2(this%nst_linear)
      this%pack%size(2) = batch_max_size(this)

      if(opencl_is_enabled()) this%pack%size(2) = opencl_padded_size(this%pack%size(2))

      this%pack%size_real = this%pack%size
      if(batch_type(this) == TYPE_CMPLX) this%pack%size_real(1) = 2*this%pack%size_real(1)

      if(opencl_is_enabled()) then
#ifdef HAVE_OPENCL
        this%status = BATCH_CL_PACKED

        call opencl_create_buffer(this%pack%buffer, CL_MEM_READ_WRITE, batch_type(this), product(this%pack%size))
#endif
      else
        this%status = BATCH_PACKED
        if(batch_type(this) == TYPE_FLOAT) then
          SAFE_ALLOCATE(this%pack%dpsi(1:this%pack%size(1), 1:this%pack%size(2)))
        else
          SAFE_ALLOCATE(this%pack%zpsi(1:this%pack%size(1), 1:this%pack%size(2)))
        end if
      end if
      
      if(copy_) then
        call profiling_in(prof_copy, "BATCH_PACK_COPY")

        this%dirty = .false.
        if(opencl_is_enabled()) then

          call batch_write_to_opencl_buffer(this)

        else
          call pack_copy()
        end if

        call profiling_out(prof_copy)
      else
        this%dirty = .true.
      end if
    end if

    INCR(this%in_buffer_count, 1)

    call profiling_out(prof)

    POP_SUB(batch_pack)

  contains

    subroutine pack_copy()
      integer :: ist, ip, sp, ep, bsize
      
 
      if(batch_type(this) == TYPE_FLOAT) then

        bsize = hardware%dblock_size
      
        !$omp parallel do private(ep, ist, ip)
        do sp = 1, this%pack%size(2), bsize
          ep = min(sp + bsize - 1, this%pack%size(2))
          forall(ist = 1:this%nst_linear)
            forall(ip = sp:ep)
              this%pack%dpsi(ist, ip) = this%states_linear(ist)%dpsi(ip)
            end forall
          end forall
        end do

        call profiling_count_transfers(this%nst_linear*this%pack%size(2), M_ONE)
      else

        bsize = hardware%zblock_size

        !$omp parallel do private(ep, ist, ip)
        do sp = 1, this%pack%size(2), bsize
          ep = min(sp + bsize - 1, this%pack%size(2))
          forall(ist = 1:this%nst_linear)
            forall(ip = sp:ep)
              this%pack%zpsi(ist, ip) = this%states_linear(ist)%zpsi(ip)
            end forall
          end forall
        end do

        call profiling_count_transfers(this%nst_linear*this%pack%size(2), M_ZI)
     end if

    end subroutine pack_copy

  end subroutine batch_pack

  ! ----------------------------------------------------

  subroutine batch_unpack(this, copy)
    type(batch_t),      intent(inout) :: this
    logical, optional,  intent(in)    :: copy

    logical :: copy_
    type(profile_t), save :: prof

    PUSH_SUB(batch_unpack)

    call profiling_in(prof, "BATCH_UNPACK")

    if(batch_is_packed(this)) then

      if(this%in_buffer_count == 1) then
        copy_ = .true.
        if(present(copy)) copy_ = copy
        
        if(copy_) call batch_sync(this)
        
        ! now deallocate
        this%status = BATCH_NOT_PACKED
        
        if(opencl_is_enabled()) then
#ifdef HAVE_OPENCL
          call opencl_release_buffer(this%pack%buffer)
#endif
        else
          SAFE_DEALLOCATE_P(this%pack%dpsi)
          SAFE_DEALLOCATE_P(this%pack%zpsi)
        end if
      end if
      
      INCR(this%in_buffer_count, -1)
    end if

    call profiling_out(prof)

    POP_SUB(batch_unpack)
    
  end subroutine batch_unpack

  ! ----------------------------------------------------

  subroutine batch_sync(this)
    type(batch_t),      intent(inout) :: this
    
    type(profile_t), save :: prof

    PUSH_SUB(batch_sync)

    if(batch_is_packed(this) .and. this%dirty) then
      call profiling_in(prof, "BATCH_UNPACK_COPY")
      
      if(opencl_is_enabled()) then
        call batch_read_from_opencl_buffer(this)
      else
        call unpack_copy()
      end if
      
      call profiling_out(prof)
      this%dirty = .false.
    end if
    
    POP_SUB(batch_sync)
    
  contains

    subroutine unpack_copy()
      integer :: ist, ip

      if(batch_type(this) == TYPE_FLOAT) then

         !$omp parallel do private(ist)
         do ip = 1, this%pack%size(2)
           forall(ist = 1:this%nst_linear)
             this%states_linear(ist)%dpsi(ip) = this%pack%dpsi(ist, ip) 
           end forall
         end do

        call profiling_count_transfers(this%nst_linear*this%pack%size(2), M_ONE)
      else
         
        !$omp parallel do private(ist)
        do ip = 1, this%pack%size(2)
          forall(ist = 1:this%nst_linear)
            this%states_linear(ist)%zpsi(ip) = this%pack%zpsi(ist, ip) 
          end forall
        end do

        call profiling_count_transfers(this%nst_linear*this%pack%size(2), M_ZI)
      end if
    end subroutine unpack_copy

  end subroutine batch_sync

  ! ----------------------------------------------------

  subroutine batch_write_to_opencl_buffer(this)
    type(batch_t),      intent(inout)  :: this

#ifdef HAVE_OPENCL
    integer :: ist, ist2, unroll
    type(opencl_mem_t) :: tmp
    type(profile_t), save :: prof_pack
    type(cl_kernel) :: kernel

    PUSH_SUB(batch_write_to_opencl_buffer)

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call opencl_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call opencl_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if

    else
      ! we copy to a temporary array and then we re-arrange data

      if(batch_type(this) == TYPE_FLOAT) then
        kernel = dpack
      else
        kernel = zpack
      end if
      
      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), unroll*this%pack%size(2))
      
      do ist = 1, this%nst_linear, unroll
        
        ! copy a number 'unroll' of states to the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)

          if(batch_type(this) == TYPE_FLOAT) then
            call opencl_write_buffer(tmp, ubound(this%states_linear(ist2)%dpsi, dim = 1), this%states_linear(ist2)%dpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          else
            call opencl_write_buffer(tmp, ubound(this%states_linear(ist2)%zpsi, dim = 1), this%states_linear(ist2)%zpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          end if
        end do

        ! now call an opencl kernel to rearrange the data
        call opencl_set_kernel_arg(kernel, 0, this%pack%size(1))
        call opencl_set_kernel_arg(kernel, 1, ist - 1)
        call opencl_set_kernel_arg(kernel, 2, tmp)
        call opencl_set_kernel_arg(kernel, 3, this%pack%buffer)

        call profiling_in(prof_pack, "CL_PACK")
        call opencl_kernel_run(kernel, (/this%pack%size(2), unroll/), (/opencl_max_workgroup_size()/unroll, unroll/))

        if(batch_type(this) == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack%size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack%size(2), M_ZI)
        end if

        call opencl_finish()
        call profiling_out(prof_pack)

      end do

      call opencl_release_buffer(tmp)

    end if

    POP_SUB(batch_write_to_opencl_buffer)
#endif
  end subroutine batch_write_to_opencl_buffer

  ! ------------------------------------------------------------------

  subroutine batch_read_from_opencl_buffer(this)
    type(batch_t),      intent(inout) :: this

#ifdef HAVE_OPENCL
    integer :: ist, ist2, unroll
    type(opencl_mem_t) :: tmp
    type(cl_kernel) :: kernel
    type(profile_t), save :: prof_unpack

    PUSH_SUB(batch_read_from_opencl_buffer)

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call opencl_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call opencl_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if
    else

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      ! we use a kernel to move to a temporary array and then we read
      call opencl_create_buffer(tmp, CL_MEM_READ_ONLY, batch_type(this), unroll*this%pack%size(2))

      if(batch_type(this) == TYPE_FLOAT) then
        kernel = dunpack
      else
        kernel = zunpack
      end if

      do ist = 1, this%nst_linear, unroll
        call opencl_set_kernel_arg(kernel, 0, this%pack%size(1))
        call opencl_set_kernel_arg(kernel, 1, ist - 1)
        call opencl_set_kernel_arg(kernel, 2, this%pack%buffer)
        call opencl_set_kernel_arg(kernel, 3, tmp)

        call profiling_in(prof_unpack, "CL_UNPACK")
        call opencl_kernel_run(kernel, (/unroll, this%pack%size(2)/), (/unroll, opencl_max_workgroup_size()/unroll/))

        if(batch_type(this) == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack%size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack%size(2), M_ZI)
        end if

        call opencl_finish()
        call profiling_out(prof_unpack)

        ! copy a number 'unroll' of states from the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)
          
          if(batch_type(this) == TYPE_FLOAT) then
            call opencl_read_buffer(tmp, ubound(this%states_linear(ist2)%dpsi, dim = 1), this%states_linear(ist2)%dpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          else
            call opencl_read_buffer(tmp, ubound(this%states_linear(ist2)%zpsi, dim = 1), this%states_linear(ist2)%zpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          end if
        end do

      end do

      call opencl_release_buffer(tmp)
    end if

    POP_SUB(batch_read_from_opencl_buffer)
#endif
  end subroutine batch_read_from_opencl_buffer


!--------------------------------------------------------------

  subroutine batch_set_zero(this)
    type(batch_t),     intent(inout) :: this

    integer :: ist_linear

    PUSH_SUB(batch_set_zero)

    call batch_pack_was_modified(this)

    if(batch_is_packed(this) .and. opencl_is_enabled()) then
#ifdef HAVE_OPENCL
      call opencl_set_buffer_to_zero(this%pack%buffer, batch_type(this), product(this%pack%size))
#endif
    else if(batch_is_packed(this) .and. batch_type(this) == TYPE_FLOAT) then
      this%pack%dpsi = M_ZERO
    else if(batch_is_packed(this) .and. batch_type(this) == TYPE_CMPLX) then
      this%pack%zpsi = M_ZERO
    else
      do ist_linear = 1, this%nst_linear
        if(associated(this%states_linear(ist_linear)%dpsi)) then
          this%states_linear(ist_linear)%dpsi = M_ZERO
        else
          this%states_linear(ist_linear)%zpsi = M_ZERO
        end if
      end do

    end if

    POP_SUB(batch_set_zero)
  end subroutine batch_set_zero

! --------------------------------------------------------------

subroutine batch_copy_data(np, xx, yy)
  integer,           intent(in)    :: np
  type(batch_t),     intent(in)    :: xx
  type(batch_t),     intent(inout) :: yy

  integer :: ist
  type(profile_t), save :: prof
  integer :: localsize

  PUSH_SUB(batch_copy_data)
  call profiling_in(prof, "BATCH_COPY_DATA")

  ASSERT(batch_type(yy) == batch_type(xx))
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))

  select case(batch_status(xx))
  case(BATCH_CL_PACKED)

#ifdef HAVE_OPENCL
    call opencl_set_kernel_arg(kernel_copy, 0, xx%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 1, log2(xx%pack%size_real(1)))
    call opencl_set_kernel_arg(kernel_copy, 2, yy%pack%buffer)
    call opencl_set_kernel_arg(kernel_copy, 3, log2(yy%pack%size_real(1)))
    
    localsize = opencl_max_workgroup_size()/yy%pack%size_real(1)
    call opencl_kernel_run(kernel_copy, (/yy%pack%size_real(1), pad(np, localsize)/), (/yy%pack%size_real(1), localsize/))
    
    call opencl_finish()
#endif

  case(BATCH_PACKED)
    if(batch_type(yy) == TYPE_FLOAT) then
      call blas_copy(np*xx%pack%size(1), xx%pack%dpsi(1, 1), 1, yy%pack%dpsi(1, 1), 1)
    else
      call blas_copy(np*xx%pack%size(1), xx%pack%zpsi(1, 1), 1, yy%pack%zpsi(1, 1), 1)
    end if

  case(BATCH_NOT_PACKED)
    !$omp parallel do private(ist)
    do ist = 1, yy%nst_linear
      if(batch_type(yy) == TYPE_CMPLX) then
        call blas_copy(np, xx%states_linear(ist)%zpsi(1), 1, yy%states_linear(ist)%zpsi(1), 1)
      else
        call blas_copy(np, xx%states_linear(ist)%dpsi(1), 1, yy%states_linear(ist)%dpsi(1), 1)
      end if
    end do

  end select

  call batch_pack_was_modified(yy)

  call profiling_out(prof)
  POP_SUB(batch_copy_data)
end subroutine batch_copy_data

! ------------------------------------------------------
integer function batch_inv_index(this, cind) result(index)
  type(batch_t),     intent(in)    :: this
  integer,           intent(in)    :: cind(:)

  do index = 1, this%nst_linear
    if(all(cind(1:this%ndims) == this%index(index, 1:this%ndims))) exit
  end do

  ASSERT(index <= this%nst_linear)

end function batch_inv_index

! ------------------------------------------------------

integer pure function batch_linear_index(this, cind) result(index)
  type(batch_t),     intent(in)    :: this
  integer,           intent(in)    :: cind(:)
  
  if(ubound(cind, dim = 1) == 1) then
    index = cind(1)
  else
    index = (cind(1) - 1)*this%dim + cind(2)
  end if

end function batch_linear_index


! --------------------------------------------------------------

subroutine batch_get_points_cl(this, sp, ep, psi, ldpsi)
  type(batch_t),       intent(in)    :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(opencl_mem_t),  intent(inout) :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif

  PUSH_SUB(batch_get_points_cl)
  call profiling_in(get_points_prof, "GET_POINTS")

  select case(batch_status(this))
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_CL_PACKED)

    tsize = types_get_size(batch_type(this))/types_get_size(TYPE_FLOAT)
    offset = this%index(1, 1) - 1

#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(kernel, 'points.cl', 'get_points')

    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, sp)
    call opencl_set_kernel_arg(kernel_ref, 1, ep)
    call opencl_set_kernel_arg(kernel_ref, 2, offset*tsize)
    call opencl_set_kernel_arg(kernel_ref, 3, this%nst_linear*tsize)
    call opencl_set_kernel_arg(kernel_ref, 4, this%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 5, this%pack%size_real(1))
    call opencl_set_kernel_arg(kernel_ref, 6, psi)
    call opencl_set_kernel_arg(kernel_ref, 7, ldpsi*tsize)

    call opencl_kernel_run(kernel_ref, (/this%pack%size_real(1), ep -&
      & sp + 1/), (/this%pack%size_real(1), 1/))
#endif
  end select

  call profiling_out(get_points_prof)

  POP_SUB(batch_get_points_cl)
end subroutine batch_get_points_cl

! --------------------------------------------------------------

subroutine batch_set_points_cl(this, sp, ep, psi, ldpsi)
  type(batch_t),       intent(inout) :: this
  integer,             intent(in)    :: sp  
  integer,             intent(in)    :: ep
  type(opencl_mem_t),  intent(in)    :: psi
  integer,             intent(in)    :: ldpsi

  integer :: tsize, offset
#ifdef HAVE_OPENCL
  type(octcl_kernel_t), save :: kernel
  type(cl_kernel)         :: kernel_ref
#endif

  PUSH_SUB(batch_set_points_cl)
  call profiling_in(set_points_prof, "SET_POINTS")

  call batch_pack_was_modified(this)

  select case(batch_status(this))
  case(BATCH_NOT_PACKED, BATCH_PACKED)
    call messages_not_implemented('batch_get_points_cl for non-CL batches')

  case(BATCH_CL_PACKED)

    tsize = types_get_size(batch_type(this))&
      &/types_get_size(TYPE_FLOAT)
    offset = this%index(1, 1) - 1
#ifdef HAVE_OPENCL
    call octcl_kernel_start_call(kernel, 'points.cl', 'set_points')
    
    kernel_ref = octcl_kernel_get_ref(kernel)

    call opencl_set_kernel_arg(kernel_ref, 0, sp)
    call opencl_set_kernel_arg(kernel_ref, 1, ep)
    call opencl_set_kernel_arg(kernel_ref, 2, offset*tsize)
    call opencl_set_kernel_arg(kernel_ref, 3, this%nst_linear*tsize)
    call opencl_set_kernel_arg(kernel_ref, 4, psi)
    call opencl_set_kernel_arg(kernel_ref, 5, ldpsi*tsize)
    call opencl_set_kernel_arg(kernel_ref, 6, this%pack%buffer)
    call opencl_set_kernel_arg(kernel_ref, 7, this%pack%size_real(1))

    call opencl_kernel_run(kernel_ref, (/this%pack%size_real(1), ep -&
      & sp + 1/), (/this%pack%size_real(1), 1/))
#endif
  end select

  call profiling_out(set_points_prof)

  POP_SUB(batch_set_points_cl)
end subroutine batch_set_points_cl

#include "real.F90"
#include "batch_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_inc.F90"
#include "undef.F90"


end module batch_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
