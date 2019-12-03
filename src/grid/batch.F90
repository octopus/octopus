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

#include "global.h"

module batch_oct_m
  use accel_oct_m
  use allocate_hardware_aware_oct_m
  use blas_oct_m
  use global_oct_m
  use hardware_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                         &
    batch_state_t,                  &
    batch_state_l_t,                &
    batch_pack_t,                   &
    batch_t,                        &
    batch_init,                     &
    batch_copy_data,                &
    batch_end,                      &
    batch_add_state,                &
    dbatch_allocate,                &
    zbatch_allocate,                &
    batch_deallocate,               &
    batch_is_packed,                &
    batch_is_ok,                    &
    batch_pack,                     &
    batch_unpack,                   &
    batch_status,                   &
    batch_type,                     &
    batch_inv_index,                &
    batch_ist_idim_to_linear,       &
    batch_linear_to_ist,            &
    batch_linear_to_idim,           &
    batch_pack_size,                &
    batch_remote_access_start,      &
    batch_remote_access_stop
  
  !--------------------------------------------------------------
  type batch_state_t
    ! Components are public by default
    FLOAT,      pointer :: dpsi(:, :)
    CMPLX,      pointer :: zpsi(:, :)
    integer        :: ist
  end type batch_state_t

  type batch_state_l_t
    ! Components are public by default
    FLOAT,      pointer :: dpsi(:)
    CMPLX,      pointer :: zpsi(:)
  end type batch_state_l_t
  
  type batch_pack_t
    ! Components are public by default
    integer                        :: size(1:2)
    integer                        :: size_real(1:2)
    FLOAT,      allocatable        :: dpsi(:, :)
    CMPLX,      allocatable        :: zpsi(:, :)
    type(accel_mem_t)             :: buffer
  end type batch_pack_t
  
  type batch_t
    private
    type(batch_state_t),   pointer, public :: states(:)
    integer,                        public :: nst
    integer                                :: current
    integer,                        public :: dim
    integer                                :: max_size

    integer                                :: ndims
    integer,               pointer         :: ist_idim_index(:, :)

    logical                                :: is_allocated
    logical                                :: mirror !< keep a copy of the batch data in unpacked form

    !> We also need a linear array with the states in order to calculate derivatives, etc.
    integer,                        public :: nst_linear
    type(batch_state_l_t), pointer, public :: states_linear(:)

    !> If the memory is contiguous, we can perform some operations faster.
    FLOAT,                 pointer, public :: dpsicont(:, :, :)
    CMPLX,                 pointer, public :: zpsicont(:, :, :)
    integer                                :: status
    integer                                :: in_buffer_count !< whether there is a copy in the opencl buffer
    type(batch_pack_t),             public :: pack
    type(type_t)                           :: type !< only available if the batched is packed
    logical :: special_memory
  contains
    procedure :: batches_are_compatible
    procedure :: copy => batch_copy
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure  batch_init_empty
    module procedure  batch_init_empty_linear
    module procedure dbatch_init_contiguous
    module procedure zbatch_init_contiguous
  end interface batch_init

  interface batch_add_state
    module procedure dbatch_add_state
    module procedure zbatch_add_state
    module procedure dbatch_add_state_linear
    module procedure zbatch_add_state_linear
  end interface batch_add_state

  integer, public, parameter :: &
    BATCH_NOT_PACKED     = 0,   &
    BATCH_PACKED         = 1,   &
    BATCH_DEVICE_PACKED  = 2

  integer, parameter :: CL_PACK_MAX_BUFFER_SIZE = 4 !< this value controls the size (in number of wave-functions)
                                                    !! of the buffer used to copy states to the opencl device.

contains

  !--------------------------------------------------------------
  subroutine batch_end(this, copy)
    class(batch_t),          intent(inout) :: this
    logical,       optional, intent(in)    :: copy

    PUSH_SUB(batch_end)

    if(this%is_allocated .and. batch_is_packed(this)) then
      !deallocate directly to avoid unnecessary copies
      this%status = BATCH_NOT_PACKED
      this%in_buffer_count = 1
      
      if(accel_is_enabled()) then
        call accel_release_buffer(this%pack%buffer)
      else
        SAFE_DEALLOCATE_A(this%pack%dpsi)
        SAFE_DEALLOCATE_A(this%pack%zpsi)
      end if
    else if(batch_is_packed(this)) then
      call batch_unpack(this, copy, force = .true.)
    end if

    if(this%is_allocated) then
      call batch_deallocate(this)
    end if

    nullify(this%dpsicont)
    nullify(this%zpsicont)

    SAFE_DEALLOCATE_P(this%states)
    SAFE_DEALLOCATE_P(this%states_linear)
    
    SAFE_DEALLOCATE_P(this%ist_idim_index)

    POP_SUB(batch_end)
  end subroutine batch_end

  !--------------------------------------------------------------
  subroutine batch_deallocate(this)
    class(batch_t),  intent(inout) :: this
    
    integer :: ii
    
    PUSH_SUB(batch_deallocate)

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
    
    if(this%special_memory) then
      if(associated(this%dpsicont)) then
        call deallocate_hardware_aware(c_loc(this%dpsicont(1,1,1)))
        nullify(this%dpsicont)
      end if
      if(associated(this%zpsicont)) then
        call deallocate_hardware_aware(c_loc(this%zpsicont(1,1,1)))
        nullify(this%zpsicont)
      end if
    else
      SAFE_DEALLOCATE_P(this%dpsicont)
      SAFE_DEALLOCATE_P(this%zpsicont)
    end if
    
    POP_SUB(batch_deallocate)
  end subroutine batch_deallocate

  !--------------------------------------------------------------

  subroutine batch_deallocate_temporary(this)
    type(batch_t),  intent(inout) :: this
    
    integer :: ii
    
    PUSH_SUB(batch_deallocate_temporary)

    do ii = 1, this%nst
      nullify(this%states(ii)%dpsi)
      nullify(this%states(ii)%zpsi)
    end do
    
    do ii = 1, this%nst_linear
      nullify(this%states_linear(ii)%dpsi)
      nullify(this%states_linear(ii)%zpsi)
    end do
    
    if(this%special_memory) then
      if(associated(this%dpsicont)) then
        call deallocate_hardware_aware(c_loc(this%dpsicont(1,1,1)))
        nullify(this%dpsicont)
      end if
      if(associated(this%zpsicont)) then
        call deallocate_hardware_aware(c_loc(this%zpsicont(1,1,1)))
        nullify(this%zpsicont)
      end if
    else
      SAFE_DEALLOCATE_P(this%dpsicont)
      SAFE_DEALLOCATE_P(this%zpsicont)

      nullify(this%dpsicont)
      nullify(this%zpsicont)
    end if
        
    POP_SUB(batch_deallocate_temporary)
  end subroutine batch_deallocate_temporary
  
  !--------------------------------------------------------------
  subroutine batch_allocate_temporary(this)
    type(batch_t),  intent(inout) :: this

    PUSH_SUB(batch_allocate_temporary)
    
    if(batch_type(this) == TYPE_FLOAT) then
      call dbatch_allocate_temporary(this)
    else if(batch_type(this) == TYPE_CMPLX) then
      call zbatch_allocate_temporary(this)
    end if

    POP_SUB(batch_allocate_temporary)
  end subroutine batch_allocate_temporary

  !--------------------------------------------------------------
  subroutine batch_init_empty (this, dim, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim
    integer,       intent(in)    :: nst
    
    integer :: ist

    PUSH_SUB(batch_init_empty)

    this%is_allocated = .false.
    this%mirror = .false.
    this%special_memory = .false.
    this%nst = nst
    this%dim = dim
    this%current = 1
    this%type = TYPE_NONE
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

    this%max_size = 0
    this%in_buffer_count = 0
    this%status = BATCH_NOT_PACKED

    this%ndims = 2
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))

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
    this%mirror = .false.
    this%special_memory = .false.
    this%nst = 0
    this%dim = 0
    this%current = 1
    this%type = TYPE_NONE
    nullify(this%dpsicont, this%zpsicont)
    nullify(this%states)

    this%nst_linear = nst
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
    end do

    this%max_size = 0
    this%in_buffer_count = 0
    this%status = BATCH_NOT_PACKED

    this%ndims = 1
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))

    POP_SUB(batch_init_empty_linear)
  end subroutine batch_init_empty_linear


  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    class(batch_t), intent(in)   :: this

    integer :: ist
    logical :: all_assoc(1:2)
    
    ! no push_sub, called too frequently
    
    ok = (this%nst_linear >= 1) .and. associated(this%states_linear)
    ok = ok .and. ubound(this%states_linear, dim = 1) == this%nst_linear
    if(ok .and. .not. batch_is_packed(this)) then
      ! ensure that either all real are associated, or all cplx are associated
      all_assoc = .true.
      do ist = 1, this%nst_linear
        all_assoc(1) = all_assoc(1) .and. associated(this%states_linear(ist)%dpsi)
        all_assoc(2) = all_assoc(2) .and. associated(this%states_linear(ist)%zpsi)
      end do

      ok = ok .and. (count(all_assoc, dim = 1) == 1)
   end if

  end function batch_is_ok

  !--------------------------------------------------------------

  subroutine batch_copy(bin, bout, pack, copy_data)
    class(batch_t),          intent(in)    :: bin
    class(batch_t),          intent(out)   :: bout
    logical,       optional, intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(bin)
    logical,       optional, intent(in)    :: copy_data  !< If .true. the new batch will be packed. Default: .false.

    integer :: ii, np
    logical :: pack_, copy_data_

    PUSH_SUB(batch_copy)

    call batch_init_empty(bout, bin%dim, bin%nst)

    copy_data_ = optional_default(copy_data, .false.)

    bout%type = bin%type

    if(batch_type(bin) == TYPE_FLOAT) then

      np = 0
      do ii = 1, bin%nst_linear
        np = max(np, ubound(bin%states_linear(ii)%dpsi, dim = 1))
      end do

      call dbatch_allocate(bout, 1, bin%nst, np)

    else if(batch_type(bin) == TYPE_CMPLX) then
      np = 0
      do ii = 1, bin%nst_linear
        np = max(np, ubound(bin%states_linear(ii)%zpsi, dim = 1))
      end do

      call zbatch_allocate(bout, 1, bin%nst, np)

    end if

    pack_ = batch_is_packed(bin)
    if(present(pack)) pack_ = pack
    if(pack_) call batch_pack(bout, copy = .false.)

    do ii = 1, bout%nst
      bout%states(ii)%ist = bin%states(ii)%ist
    end do

    bout%ist_idim_index(1:bin%nst_linear, 1:bin%ndims) = bin%ist_idim_index(1:bin%nst_linear, 1:bin%ndims)

    if(copy_data_) call batch_copy_data(np, bin, bout)
    
    POP_SUB(batch_copy)
  end subroutine batch_copy

  ! ----------------------------------------------------
  !> THREADSAFE
  type(type_t) pure function batch_type(this) result(btype)
    class(batch_t),      intent(in)    :: this

    if(.not. batch_is_packed(this)) then
      if(associated(this%states_linear(1)%dpsi)) btype = TYPE_FLOAT
      if(associated(this%states_linear(1)%zpsi)) btype = TYPE_CMPLX
    else
      btype = this%type
    end if
     
  end function batch_type

  ! ----------------------------------------------------
  !> THREADSAFE
  integer pure function batch_status(this) result(bstatus)
    class(batch_t),      intent(in)    :: this

    bstatus = this%status
  end function batch_status
  
  ! ----------------------------------------------------

  logical pure function batch_is_packed(this) result(in_buffer)
    class(batch_t),      intent(in)    :: this

    in_buffer = this%in_buffer_count > 0
  end function batch_is_packed

  ! ----------------------------------------------------

  integer pure function batch_max_size(this) result(size)
    class(batch_t),      intent(in)    :: this

    size = this%max_size
  end function batch_max_size

  ! ----------------------------------------------------

  integer function batch_pack_size(this) result(size)
    class(batch_t),      intent(inout) :: this

    size = batch_max_size(this)
    if(accel_is_enabled()) size = accel_padded_size(size)
    size = size*pad_pow2(this%nst_linear)*types_get_size(batch_type(this))

  end function batch_pack_size

  ! ----------------------------------------------------

  subroutine batch_pack(this, copy)
    class(batch_t),      intent(inout) :: this
    logical,   optional, intent(in)    :: copy

    logical :: copy_
    type(profile_t), save :: prof, prof_copy

    ! no push_sub, called too frequently

    call profiling_in(prof, "BATCH_PACK")
    ASSERT(batch_is_ok(this))

    copy_ = .true.
    if(present(copy)) copy_ = copy

    if(.not. batch_is_packed(this)) then
      this%type = batch_type(this)
      this%pack%size(1) = pad_pow2(this%nst_linear)
      this%pack%size(2) = batch_max_size(this)

      if(accel_is_enabled()) this%pack%size(2) = accel_padded_size(this%pack%size(2))

      this%pack%size_real = this%pack%size
      if(type_is_complex(batch_type(this))) this%pack%size_real(1) = 2*this%pack%size_real(1)

      if(accel_is_enabled()) then
        this%status = BATCH_DEVICE_PACKED
        call accel_create_buffer(this%pack%buffer, ACCEL_MEM_READ_WRITE, batch_type(this), product(this%pack%size))
      else
        this%status = BATCH_PACKED
        if(batch_type(this) == TYPE_FLOAT) then
          SAFE_ALLOCATE(this%pack%dpsi(1:this%pack%size(1), 1:this%pack%size(2)))
        else if(batch_type(this) == TYPE_CMPLX) then
          SAFE_ALLOCATE(this%pack%zpsi(1:this%pack%size(1), 1:this%pack%size(2)))
        end if
      end if
      
      if(copy_) then
        call profiling_in(prof_copy, "BATCH_PACK_COPY")
        if(accel_is_enabled()) then
          call batch_write_to_opencl_buffer(this)
        else
          call pack_copy()
        end if

        call profiling_out(prof_copy)
      end if

      if(this%is_allocated .and. .not. this%mirror) call batch_deallocate_temporary(this)

    end if

    INCR(this%in_buffer_count, 1)

    call profiling_out(prof)

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

      else if(batch_type(this) == TYPE_CMPLX) then

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

      end if

      call profiling_count_transfers(this%nst_linear*this%pack%size(2), batch_type(this))
      
    end subroutine pack_copy

  end subroutine batch_pack

  ! ----------------------------------------------------

  subroutine batch_unpack(this, copy, force)
    class(batch_t),     intent(inout) :: this
    logical, optional,  intent(in)    :: copy
    logical, optional,  intent(in)    :: force  !< if force = .true., unpack independently of the counter

    logical :: copy_
    type(profile_t), save :: prof

    PUSH_SUB(batch_unpack)

    call profiling_in(prof, "BATCH_UNPACK")

    if(batch_is_packed(this)) then

      if(this%in_buffer_count == 1 .or. optional_default(force, .false.)) then

        if(this%is_allocated .and. .not. this%mirror) call batch_allocate_temporary(this)
        
        copy_ = .true.
        if(present(copy)) copy_ = copy
        if(this%is_allocated .and. .not. this%mirror) copy_ = .true.
        
        if(copy_) call batch_sync(this)
        
        ! now deallocate
        this%status = BATCH_NOT_PACKED
        this%in_buffer_count = 1

        if(accel_is_enabled()) then
          call accel_release_buffer(this%pack%buffer)
        else
          SAFE_DEALLOCATE_A(this%pack%dpsi)
          SAFE_DEALLOCATE_A(this%pack%zpsi)
        end if
      end if
      
      INCR(this%in_buffer_count, -1)
    end if

    call profiling_out(prof)

    POP_SUB(batch_unpack)

  end subroutine batch_unpack

  ! ----------------------------------------------------

  subroutine batch_sync(this)
    class(batch_t),      intent(inout) :: this
    
    type(profile_t), save :: prof

    PUSH_SUB(batch_sync)

    if(batch_is_packed(this)) then
      call profiling_in(prof, "BATCH_UNPACK_COPY")
      
      if(accel_is_enabled()) then
        call batch_read_from_opencl_buffer(this)
      else
        call unpack_copy()
      end if
      
      call profiling_out(prof)
    end if
    
    POP_SUB(batch_sync)
    
  contains

    subroutine unpack_copy()
      integer :: ist, ip

      if(batch_type(this) == TYPE_FLOAT) then

        do ist = 1, this%nst_linear
          ASSERT(associated(this%states_linear(ist)%dpsi))
        end do
        
        !$omp parallel do private(ist)
        do ip = 1, this%pack%size(2)
          forall(ist = 1:this%nst_linear)
            this%states_linear(ist)%dpsi(ip) = this%pack%dpsi(ist, ip) 
          end forall
        end do
        
      else if(batch_type(this) == TYPE_CMPLX) then

        do ist = 1, this%nst_linear
          ASSERT(associated(this%states_linear(ist)%zpsi))
        end do
        
        !$omp parallel do private(ist)
        do ip = 1, this%pack%size(2)
          forall(ist = 1:this%nst_linear)
            this%states_linear(ist)%zpsi(ip) = this%pack%zpsi(ist, ip) 
          end forall
        end do
        
      end if
      
      call profiling_count_transfers(this%nst_linear*this%pack%size(2), batch_type(this))
      
    end subroutine unpack_copy

  end subroutine batch_sync

  ! ----------------------------------------------------

  subroutine batch_write_to_opencl_buffer(this)
    class(batch_t),      intent(inout)  :: this

    integer :: ist, ist2, unroll
    type(accel_mem_t) :: tmp
    type(profile_t), save :: prof_pack
    type(accel_kernel_t), pointer :: kernel

    PUSH_SUB(batch_write_to_opencl_buffer)

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call accel_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else if(batch_type(this) == TYPE_CMPLX) then
        call accel_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      else
        ASSERT(.false.)
      end if

    else
      ! we copy to a temporary array and then we re-arrange data

      if(batch_type(this) == TYPE_FLOAT) then
        kernel => dpack
      else
        kernel => zpack
      end if
      
      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      call accel_create_buffer(tmp, ACCEL_MEM_READ_ONLY, batch_type(this), unroll*this%pack%size(2))
      
      do ist = 1, this%nst_linear, unroll
        
        ! copy a number 'unroll' of states to the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)

          if(batch_type(this) == TYPE_FLOAT) then
            call accel_write_buffer(tmp, ubound(this%states_linear(ist2)%dpsi, dim = 1), this%states_linear(ist2)%dpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          else
            call accel_write_buffer(tmp, ubound(this%states_linear(ist2)%zpsi, dim = 1), this%states_linear(ist2)%zpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          end if
        end do

        ! now call an opencl kernel to rearrange the data
        call accel_set_kernel_arg(kernel, 0, this%pack%size(1))
        call accel_set_kernel_arg(kernel, 1, this%pack%size(2))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, tmp)
        call accel_set_kernel_arg(kernel, 4, this%pack%buffer)

        call profiling_in(prof_pack, "CL_PACK")
        call accel_kernel_run(kernel, (/this%pack%size(2), unroll/), (/accel_max_workgroup_size()/unroll, unroll/))

        if(batch_type(this) == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack%size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack%size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_pack)

      end do

      call accel_release_buffer(tmp)

    end if

    POP_SUB(batch_write_to_opencl_buffer)
  end subroutine batch_write_to_opencl_buffer

  ! ------------------------------------------------------------------

  subroutine batch_read_from_opencl_buffer(this)
    class(batch_t),      intent(inout) :: this

    integer :: ist, ist2, unroll
    type(accel_mem_t) :: tmp
    type(accel_kernel_t), pointer :: kernel
    type(profile_t), save :: prof_unpack

    PUSH_SUB(batch_read_from_opencl_buffer)

    ASSERT(batch_is_ok(this))

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(batch_type(this) == TYPE_FLOAT) then
        call accel_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call accel_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if
    else

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      ! we use a kernel to move to a temporary array and then we read
      call accel_create_buffer(tmp, ACCEL_MEM_WRITE_ONLY, batch_type(this), unroll*this%pack%size(2))

      if(batch_type(this) == TYPE_FLOAT) then
        kernel => dunpack
      else
        kernel => zunpack
      end if

      do ist = 1, this%nst_linear, unroll
        call accel_set_kernel_arg(kernel, 0, this%pack%size(1))
        call accel_set_kernel_arg(kernel, 1, this%pack%size(2))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, this%pack%buffer)
        call accel_set_kernel_arg(kernel, 4, tmp)

        call profiling_in(prof_unpack, "CL_UNPACK")
        call accel_kernel_run(kernel, (/unroll, this%pack%size(2)/), (/unroll, accel_max_workgroup_size()/unroll/))

        if(batch_type(this) == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack%size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack%size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_unpack)

        ! copy a number 'unroll' of states from the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)
          
          if(batch_type(this) == TYPE_FLOAT) then
            call accel_read_buffer(tmp, ubound(this%states_linear(ist2)%dpsi, dim = 1), this%states_linear(ist2)%dpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          else
            call accel_read_buffer(tmp, ubound(this%states_linear(ist2)%zpsi, dim = 1), this%states_linear(ist2)%zpsi, &
              offset = (ist2 - ist)*this%pack%size(2))
          end if
        end do

      end do

      call accel_release_buffer(tmp)
    end if

    POP_SUB(batch_read_from_opencl_buffer)
  end subroutine batch_read_from_opencl_buffer

! ------------------------------------------------------
integer function batch_inv_index(this, cind) result(index)
  class(batch_t),     intent(in)    :: this
  integer,            intent(in)    :: cind(:)

  do index = 1, this%nst_linear
    if(all(cind(1:this%ndims) == this%ist_idim_index(index, 1:this%ndims))) exit
  end do

  ASSERT(index <= this%nst_linear)

end function batch_inv_index

! ------------------------------------------------------

integer pure function batch_ist_idim_to_linear(this, cind) result(index)
  class(batch_t),     intent(in)    :: this
  integer,            intent(in)    :: cind(:)
  
  if(ubound(cind, dim = 1) == 1) then
    index = cind(1)
  else
    index = (cind(1) - 1)*this%dim + cind(2)
  end if

end function batch_ist_idim_to_linear

! ------------------------------------------------------

integer pure function batch_linear_to_ist(this, linear_index) result(ist)
  class(batch_t),     intent(in)    :: this
  integer,            intent(in)    :: linear_index
  
  ist = this%ist_idim_index(linear_index, 1)

end function batch_linear_to_ist

! ------------------------------------------------------

integer pure function batch_linear_to_idim(this, linear_index) result(idim)
  class(batch_t),     intent(in)    :: this
  integer,            intent(in)    :: linear_index
  
  idim = this%ist_idim_index(linear_index, 2)
  
end function batch_linear_to_idim

! ------------------------------------------------------

subroutine batch_remote_access_start(this, mpi_grp, rma_win)
  class(batch_t),  intent(inout) :: this
  type(mpi_grp_t), intent(in)    :: mpi_grp
  integer,         intent(out)   :: rma_win

  PUSH_SUB(batch_remote_access_start)

  ASSERT(.not. accel_is_enabled())
  
  if(mpi_grp%size > 1) then
    call batch_pack(this)
    
    if(batch_type(this) == TYPE_CMPLX) then
#ifdef HAVE_MPI2
      call MPI_Win_create(this%pack%zpsi(1, 1), int(product(this%pack%size)*types_get_size(batch_type(this)), MPI_ADDRESS_KIND), &
        types_get_size(batch_type(this)), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
    end if
    
    if(batch_type(this) == TYPE_FLOAT) then
#ifdef HAVE_MPI2
      call MPI_Win_create(this%pack%dpsi(1, 1), int(product(this%pack%size)*types_get_size(batch_type(this)), MPI_ADDRESS_KIND), &
        types_get_size(batch_type(this)), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
    end if
    
  else
    rma_win = -1
  end if
  
  POP_SUB(batch_remote_access_start)
end subroutine batch_remote_access_start

! ------------------------------------------------------

subroutine batch_remote_access_stop(this, rma_win)
  class(batch_t),     intent(inout) :: this
  integer,            intent(inout) :: rma_win
  
  PUSH_SUB(batch_remote_access_stop)

  if(rma_win /= -1) then
#ifdef HAVE_MPI2
    call MPI_Win_free(rma_win, mpi_err)
#endif
    call batch_unpack(this)
  end if
  
  POP_SUB(batch_remote_access_stop)
end subroutine batch_remote_access_stop

! --------------------------------------------------------------

subroutine batch_copy_data(np, xx, yy)
  integer,           intent(in)    :: np
  class(batch_t),    intent(in)    :: xx
  class(batch_t),    intent(inout) :: yy

  integer :: ist, dim2, dim3
  type(profile_t), save :: prof
  integer :: localsize

  PUSH_SUB(batch_copy_data)
  call profiling_in(prof, "BATCH_COPY_DATA")

  call xx%batches_are_compatible(yy)

  select case(batch_status(xx))
  case(BATCH_DEVICE_PACKED)
    call accel_set_kernel_arg(kernel_copy, 0, np)
    call accel_set_kernel_arg(kernel_copy, 1, xx%pack%buffer)
    call accel_set_kernel_arg(kernel_copy, 2, log2(xx%pack%size_real(1)))
    call accel_set_kernel_arg(kernel_copy, 3, yy%pack%buffer)
    call accel_set_kernel_arg(kernel_copy, 4, log2(yy%pack%size_real(1)))
    
    localsize = accel_kernel_workgroup_size(kernel_copy)/yy%pack%size_real(1)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
    
    call accel_kernel_run(kernel_copy, (/yy%pack%size_real(1), dim2, dim3/), (/yy%pack%size_real(1), localsize, 1/))
    
    call accel_finish()

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

  call profiling_out(prof)
  POP_SUB(batch_copy_data)
end subroutine batch_copy_data

! --------------------------------------------------------------

subroutine batches_are_compatible(xx, yy)
  class(batch_t),    intent(in) :: xx
  class(batch_t),    intent(in) :: yy

  PUSH_SUB(batches_are_compatible)

  ASSERT(batch_type(yy) == batch_type(xx))
  ASSERT(xx%nst_linear == yy%nst_linear)
  ASSERT(batch_status(xx) == batch_status(yy))
  ASSERT(xx%dim == yy%dim)

  POP_SUB(batches_are_compatible)

end subroutine batches_are_compatible

#include "real.F90"
#include "batch_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "batch_inc.F90"
#include "undef.F90"

end module batch_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
