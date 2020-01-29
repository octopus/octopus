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
  use iso_c_binding
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
    batch_init
  
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
    integer                                :: status_of
    integer                                :: in_buffer_count !< whether there is a copy in the opencl buffer
    type(batch_pack_t),             public :: pack
    type(type_t)                           :: type_of !< only available if the batched is packed
    logical :: special_memory
  contains
    procedure :: dbatch_add_state
    procedure :: zbatch_add_state
    procedure :: dbatch_add_state_linear
    procedure :: zbatch_add_state_linear
    generic   :: add_state => dbatch_add_state, zbatch_add_state, dbatch_add_state_linear, zbatch_add_state_linear
    procedure ::  dallocate => dbatch_allocate
    procedure ::  zallocate => zbatch_allocate
    procedure :: check_compatibility_with => batch_check_compatibility_with
    procedure :: clone_to => batch_clone_to
    procedure :: clone_to_array => batch_clone_to_array
    procedure :: copy_to => batch_copy_to
    procedure :: copy_data_to => batch_copy_data_to
    procedure :: deallocate => batch_deallocate
    procedure :: do_pack => batch_do_pack
    procedure :: do_unpack => batch_do_unpack
    procedure :: end => batch_end
    procedure :: inv_index => batch_inv_index
    procedure :: is_ok => batch_is_ok
    procedure :: is_packed => batch_is_packed
    procedure :: ist_idim_to_linear => batch_ist_idim_to_linear
    procedure :: linear_to_idim => batch_linear_to_idim
    procedure :: linear_to_ist => batch_linear_to_ist
    procedure :: pack_size => batch_pack_size
    procedure :: remote_access_start => batch_remote_access_start
    procedure :: remote_access_stop => batch_remote_access_stop
    procedure :: status => batch_status
    procedure :: type => batch_type
    procedure :: type_as_int => batch_type_as_integer
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure  batch_init_empty
    module procedure dbatch_init_contiguous
    module procedure zbatch_init_contiguous
    module procedure dbatch_init_contiguous_2d
    module procedure zbatch_init_contiguous_2d
    module procedure dbatch_init_single
    module procedure zbatch_init_single
  end interface batch_init

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

    if(this%is_allocated .and. this%is_packed()) then
      !deallocate directly to avoid unnecessary copies
      this%status_of = BATCH_NOT_PACKED
      this%in_buffer_count = 1
      
      if(accel_is_enabled()) then
        call accel_release_buffer(this%pack%buffer)
      else
        SAFE_DEALLOCATE_A(this%pack%dpsi)
        SAFE_DEALLOCATE_A(this%pack%zpsi)
      end if
    else if(this%is_packed()) then
      call this%do_unpack(copy, force = .true.)
    end if

    if(this%is_allocated) then
      call this%deallocate()
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
    
    if(this%type() == TYPE_FLOAT) then
      call dbatch_allocate_temporary(this)
    else if(this%type() == TYPE_CMPLX) then
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
    this%type_of = TYPE_NONE
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
    this%status_of = BATCH_NOT_PACKED

    this%ndims = 2
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))

    POP_SUB(batch_init_empty)
  end subroutine batch_init_empty

  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    class(batch_t), intent(in)   :: this

    integer :: ist
    logical :: all_assoc(1:2)
    
    ! no push_sub, called too frequently
    
    ok = (this%nst_linear >= 1) .and. associated(this%states_linear)
    ok = ok .and. ubound(this%states_linear, dim = 1) == this%nst_linear
    if(ok .and. .not. this%is_packed()) then
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

  subroutine batch_clone_to(this, dest, pack, copy_data)
    class(batch_t),              intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest
    logical,        optional,    intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,        optional,    intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.

    PUSH_SUB(batch_clone_to)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE(batch_t, dest)
    else
      message(1) = "Internal error: destination batch in batch_clone_to has been previously allocated."
      call messages_fatal(1)
    end if

    call this%copy_to(dest, pack, copy_data)

    POP_SUB(batch_clone_to)
  end subroutine batch_clone_to

  !--------------------------------------------------------------

  subroutine batch_clone_to_array(this, dest, n_batches, pack, copy_data)
    class(batch_t),              intent(in)    :: this
    class(batch_t), allocatable, intent(out)   :: dest(:)
    integer,                     intent(in)    :: n_batches
    logical,        optional,    intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,        optional,    intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.

    integer :: ib

    PUSH_SUB(batch_clone_to_array)

    if (.not. allocated(dest)) then
      SAFE_ALLOCATE_TYPE_ARRAY(batch_t, dest, (1:n_batches))
    else
      message(1) = "Internal error: destination batch in batch_clone_to_array has been previously allocated."
      call messages_fatal(1)
    end if

    do ib = 1, n_batches
      call this%copy_to(dest(ib), pack, copy_data)
    end do

    POP_SUB(batch_clone_to_array)
  end subroutine batch_clone_to_array
  
  !--------------------------------------------------------------

  subroutine batch_copy_to(this, dest, pack, copy_data)
    class(batch_t),          intent(in)    :: this
    class(batch_t),          intent(out)   :: dest
    logical,       optional, intent(in)    :: pack       !< If .false. the new batch will not be packed. Default: batch_is_packed(this)
    logical,       optional, intent(in)    :: copy_data  !< If .true. the batch data will be copied to the destination batch. Default: .false.

    integer :: ii, np

    PUSH_SUB(batch_copy_to)

    call batch_init_empty(dest, this%dim, this%nst)

    dest%type_of = this%type_of

    if(this%type() == TYPE_FLOAT) then

      np = 0
      do ii = 1, this%nst_linear
        np = max(np, ubound(this%states_linear(ii)%dpsi, dim = 1))
      end do

      call dest%dallocate(1, this%nst, np)

    else if(this%type() == TYPE_CMPLX) then
      np = 0
      do ii = 1, this%nst_linear
        np = max(np, ubound(this%states_linear(ii)%zpsi, dim = 1))
      end do

      call dest%zallocate(1, this%nst, np)

    else
      message(1) = "Internal error: unknown batch type in batch_copy_to."
      call messages_fatal(1)
    end if

    if(optional_default(pack, this%is_packed())) call dest%do_pack(copy = .false.)

    do ii = 1, dest%nst
      dest%states(ii)%ist = this%states(ii)%ist
    end do

    dest%ist_idim_index(1:this%nst_linear, 1:this%ndims) = this%ist_idim_index(1:this%nst_linear, 1:this%ndims)

    if(optional_default(copy_data, .false.)) call this%copy_data_to(np, dest)
    
    POP_SUB(batch_copy_to)
  end subroutine batch_copy_to

  ! ----------------------------------------------------
  !> THREADSAFE
  type(type_t) pure function batch_type(this) result(btype)
    class(batch_t),      intent(in)    :: this

    if(.not. this%is_packed()) then
      if(associated(this%states_linear(1)%dpsi)) btype = TYPE_FLOAT
      if(associated(this%states_linear(1)%zpsi)) btype = TYPE_CMPLX
    else
      btype = this%type_of
    end if
     
  end function batch_type

  ! ----------------------------------------------------
  !> For debuging purpose only
  integer pure function batch_type_as_integer(this) result(itype)
    class(batch_t),      intent(in)    :: this

    type(type_t) :: btype

    itype = 0
    btype = this%type()
    if( btype == TYPE_FLOAT ) itype = 1
    if( btype == TYPE_CMPLX ) itype = 2

  end function batch_type_as_integer

  ! ----------------------------------------------------
  !> THREADSAFE
  integer pure function batch_status(this) result(bstatus)
    class(batch_t),      intent(in)    :: this

    bstatus = this%status_of
  end function batch_status
  
  ! ----------------------------------------------------

  logical pure function batch_is_packed(this) result(in_buffer)
    class(batch_t),      intent(in)    :: this

    in_buffer = this%in_buffer_count > 0
  end function batch_is_packed

  ! ----------------------------------------------------

  integer function batch_pack_size(this) result(size)
    class(batch_t),      intent(inout) :: this

    size = this%max_size
    if(accel_is_enabled()) size = accel_padded_size(size)
    size = size*pad_pow2(this%nst_linear)*types_get_size(this%type())

  end function batch_pack_size

  ! ----------------------------------------------------

  subroutine batch_do_pack(this, copy)
    class(batch_t),      intent(inout) :: this
    logical,   optional, intent(in)    :: copy

    logical :: copy_
    type(profile_t), save :: prof, prof_copy

    ! no push_sub, called too frequently

    call profiling_in(prof, "BATCH_DO_PACK")
    ASSERT(this%is_ok())

    copy_ = .true.
    if(present(copy)) copy_ = copy

    if(.not. this%is_packed()) then
      this%type_of = this%type()
      this%pack%size(1) = pad_pow2(this%nst_linear)
      this%pack%size(2) = this%max_size

      if(accel_is_enabled()) this%pack%size(2) = accel_padded_size(this%pack%size(2))

      this%pack%size_real = this%pack%size
      if(type_is_complex(this%type())) this%pack%size_real(1) = 2*this%pack%size_real(1)

      if(accel_is_enabled()) then
        this%status_of = BATCH_DEVICE_PACKED
        call accel_create_buffer(this%pack%buffer, ACCEL_MEM_READ_WRITE, this%type(), product(this%pack%size))
      else
        this%status_of = BATCH_PACKED
        if(this%type() == TYPE_FLOAT) then
          SAFE_ALLOCATE(this%pack%dpsi(1:this%pack%size(1), 1:this%pack%size(2)))
        else if(this%type() == TYPE_CMPLX) then
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
      
 
      if(this%type() == TYPE_FLOAT) then

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

      else if(this%type() == TYPE_CMPLX) then

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

      else
        message(1) = "Internal error: unknown batch type in batch_do_pack."
        call messages_fatal(1)
      end if

      call profiling_count_transfers(this%nst_linear*this%pack%size(2), this%type())
      
    end subroutine pack_copy

  end subroutine batch_do_pack

  ! ----------------------------------------------------

  subroutine batch_do_unpack(this, copy, force)
    class(batch_t),     intent(inout) :: this
    logical, optional,  intent(in)    :: copy
    logical, optional,  intent(in)    :: force  !< if force = .true., unpack independently of the counter

    logical :: copy_
    type(profile_t), save :: prof

    PUSH_SUB(batch_do_unpack)

    call profiling_in(prof, "BATCH_DO_UNPACK")

    if(this%is_packed()) then

      if(this%in_buffer_count == 1 .or. optional_default(force, .false.)) then

        if(this%is_allocated .and. .not. this%mirror) call batch_allocate_temporary(this)
        
        copy_ = .true.
        if(present(copy)) copy_ = copy
        if(this%is_allocated .and. .not. this%mirror) copy_ = .true.
        
        if(copy_) call batch_sync(this)
        
        ! now deallocate
        this%status_of = BATCH_NOT_PACKED
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

    POP_SUB(batch_do_unpack)

  end subroutine batch_do_unpack

  ! ----------------------------------------------------

  subroutine batch_sync(this)
    class(batch_t),      intent(inout) :: this
    
    type(profile_t), save :: prof

    PUSH_SUB(batch_sync)

    if(this%is_packed()) then
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

      if(this%type() == TYPE_FLOAT) then

        do ist = 1, this%nst_linear
          ASSERT(associated(this%states_linear(ist)%dpsi))
        end do
        
        !$omp parallel do private(ist)
        do ip = 1, this%pack%size(2)
          forall(ist = 1:this%nst_linear)
            this%states_linear(ist)%dpsi(ip) = this%pack%dpsi(ist, ip) 
          end forall
        end do
        
      else if(this%type() == TYPE_CMPLX) then

        do ist = 1, this%nst_linear
          ASSERT(associated(this%states_linear(ist)%zpsi))
        end do
        
        !$omp parallel do private(ist)
        do ip = 1, this%pack%size(2)
          forall(ist = 1:this%nst_linear)
            this%states_linear(ist)%zpsi(ip) = this%pack%zpsi(ist, ip) 
          end forall
        end do
        
      else
        message(1) = "Internal error: unknown batch type in batch_sync."
        call messages_fatal(1)
      end if
      
      call profiling_count_transfers(this%nst_linear*this%pack%size(2), this%type())
      
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

    ASSERT(this%is_ok())

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(this%type() == TYPE_FLOAT) then
        call accel_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else if(this%type() == TYPE_CMPLX) then
        call accel_write_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      else
        ASSERT(.false.)
      end if

    else
      ! we copy to a temporary array and then we re-arrange data

      if(this%type() == TYPE_FLOAT) then
        kernel => dpack
      else
        kernel => zpack
      end if
      
      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      call accel_create_buffer(tmp, ACCEL_MEM_READ_ONLY, this%type(), unroll*this%pack%size(2))
      
      do ist = 1, this%nst_linear, unroll
        
        ! copy a number 'unroll' of states to the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)

          if(this%type() == TYPE_FLOAT) then
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

        if(this%type() == TYPE_FLOAT) then
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

    ASSERT(this%is_ok())

    if(this%nst_linear == 1) then
      ! we can copy directly
      if(this%type() == TYPE_FLOAT) then
        call accel_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%dpsi, dim = 1), this%states_linear(1)%dpsi)
      else
        call accel_read_buffer(this%pack%buffer, ubound(this%states_linear(1)%zpsi, dim = 1), this%states_linear(1)%zpsi)
      end if
    else

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack%size(1))

      ! we use a kernel to move to a temporary array and then we read
      call accel_create_buffer(tmp, ACCEL_MEM_WRITE_ONLY, this%type(), unroll*this%pack%size(2))

      if(this%type() == TYPE_FLOAT) then
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

        if(this%type() == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack%size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack%size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_unpack)

        ! copy a number 'unroll' of states from the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)
          
          if(this%type() == TYPE_FLOAT) then
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
    call this%do_pack()
    
    if(this%type() == TYPE_CMPLX) then
#ifdef HAVE_MPI2
      call MPI_Win_create(this%pack%zpsi(1, 1), int(product(this%pack%size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
        types_get_size(this%type()), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
    else if (this%type() == TYPE_FLOAT) then
#ifdef HAVE_MPI2
      call MPI_Win_create(this%pack%dpsi(1, 1), int(product(this%pack%size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
        types_get_size(this%type()), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
    else
      message(1) = "Internal error: unknown batch type in batch_remote_access_start."
      call messages_fatal(1)
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
    call this%do_unpack()
  end if
  
  POP_SUB(batch_remote_access_stop)
end subroutine batch_remote_access_stop

! --------------------------------------------------------------

subroutine batch_copy_data_to(this, np, dest)
  class(batch_t),    intent(in)    :: this
  integer,           intent(in)    :: np
  class(batch_t),    intent(inout) :: dest

  integer :: ist, dim2, dim3
  type(profile_t), save :: prof
  integer :: localsize

  PUSH_SUB(batch_copy_data_to)
  call profiling_in(prof, "BATCH_COPY_DATA_TO")

  call this%check_compatibility_with(dest)

  select case(this%status())
  case(BATCH_DEVICE_PACKED)
    call accel_set_kernel_arg(kernel_copy, 0, np)
    call accel_set_kernel_arg(kernel_copy, 1, this%pack%buffer)
    call accel_set_kernel_arg(kernel_copy, 2, log2(this%pack%size_real(1)))
    call accel_set_kernel_arg(kernel_copy, 3, dest%pack%buffer)
    call accel_set_kernel_arg(kernel_copy, 4, log2(dest%pack%size_real(1)))
    
    localsize = accel_kernel_workgroup_size(kernel_copy)/dest%pack%size_real(1)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
    
    call accel_kernel_run(kernel_copy, (/dest%pack%size_real(1), dim2, dim3/), (/dest%pack%size_real(1), localsize, 1/))
    
    call accel_finish()

  case(BATCH_PACKED)
    if(dest%type() == TYPE_FLOAT) then
      call blas_copy(np*this%pack%size(1), this%pack%dpsi(1, 1), 1, dest%pack%dpsi(1, 1), 1)
    else
      call blas_copy(np*this%pack%size(1), this%pack%zpsi(1, 1), 1, dest%pack%zpsi(1, 1), 1)
    end if

  case(BATCH_NOT_PACKED)
    !$omp parallel do private(ist)
    do ist = 1, dest%nst_linear
      if(dest%type() == TYPE_CMPLX) then
        call blas_copy(np, this%states_linear(ist)%zpsi(1), 1, dest%states_linear(ist)%zpsi(1), 1)
      else
        call blas_copy(np, this%states_linear(ist)%dpsi(1), 1, dest%states_linear(ist)%dpsi(1), 1)
      end if
    end do

  end select

  call profiling_out(prof)
  POP_SUB(batch_copy_data_to)
end subroutine batch_copy_data_to

! --------------------------------------------------------------

subroutine batch_check_compatibility_with(this, target, only_check_dim)
  class(batch_t),    intent(in) :: this
  class(batch_t),    intent(in) :: target
  logical, optional, intent(in) :: only_check_dim

  PUSH_SUB(batch_check_compatibility_with)

  ASSERT(this%type() == target%type())
  if(.not. optional_default(only_check_dim, .false.)) then
    ASSERT(this%nst_linear == target%nst_linear)
  end if
  ASSERT(this%status() == target%status())
  ASSERT(this%dim == target%dim)

  POP_SUB(batch_check_compatibility_with)

end subroutine batch_check_compatibility_with

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
