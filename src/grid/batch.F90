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
    batch_t,                        &
    batch_init,                     &
    dbatch_init,                    &
    zbatch_init
  
  type batch_t
    private
    integer,                        public :: nst
    integer                                :: current
    integer,                        public :: dim
    integer                                :: max_size

    integer                                :: ndims
    integer,               pointer         :: ist_idim_index(:, :)
    integer,           allocatable, public :: ist(:)

    logical                                :: is_allocated
    logical                                :: mirror !< keep a copy of the batch data in unpacked form

    !> We also need a linear array with the states in order to calculate derivatives, etc.
    integer,                        public :: nst_linear

    integer                                :: status_of
    integer                                :: in_buffer_count !< whether there is a copy in the opencl buffer
    logical :: special_memory


    !> unpacked variables; linear variables are pointers with different shapes
    FLOAT, pointer, contiguous, public :: dff(:, :, :)
    CMPLX, pointer, contiguous, public :: zff(:, :, :)
    FLOAT, pointer, contiguous, public :: dff_linear(:, :)
    CMPLX, pointer, contiguous, public :: zff_linear(:, :)
    !> packed variables; only rank-2 arrays due to padding to powers of 2
    FLOAT, pointer, contiguous, public :: dff_pack(:, :)
    CMPLX, pointer, contiguous, public :: zff_pack(:, :)

    integer                   , public :: pack_size(1:2)
    integer                   , public :: pack_size_real(1:2)

    type(accel_mem_t)         , public    :: ff_device


    integer, public :: layout !< either BATCH_NOT_PACKED or BATCH_PACKED
    integer, public :: location !< either BATCH_HOST or BATCH_DEVICE
    type(type_t) :: type_of !< either TYPE_FLOAT or TYPE_COMPLEX

  contains
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
    procedure :: pack_total_size => batch_pack_total_size
    procedure :: remote_access_start => batch_remote_access_start
    procedure :: remote_access_stop => batch_remote_access_stop
    procedure :: status => batch_status
    procedure :: type => batch_type
    procedure :: type_as_int => batch_type_as_integer
  end type batch_t

  !--------------------------------------------------------------
  interface batch_init
    module procedure dbatch_init_with_memory_3
    module procedure zbatch_init_with_memory_3
    module procedure dbatch_init_with_memory_2
    module procedure zbatch_init_with_memory_2
    module procedure dbatch_init_with_memory_1
    module procedure zbatch_init_with_memory_1
  end interface batch_init

  integer, public, parameter :: &
    BATCH_NOT_PACKED     = 0,   &
    BATCH_PACKED         = 1,   &
    BATCH_DEVICE_PACKED  = 2,   &
    BATCH_HOST           = 3,   &
    BATCH_DEVICE         = 4

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
        call accel_release_buffer(this%ff_device)
      else
        if(associated(this%dff_pack)) then
          call deallocate_hardware_aware(c_loc(this%dff_pack(1,1)))
        end if
        if(associated(this%zff_pack)) then
          call deallocate_hardware_aware(c_loc(this%zff_pack(1,1)))
        end if
        nullify(this%dff_pack)
        nullify(this%zff_pack)
      end if
    else if(this%is_packed()) then
      call this%do_unpack(copy, force = .true.)
    end if

    if(this%is_allocated) then
      call this%deallocate()
    end if

    SAFE_DEALLOCATE_P(this%ist_idim_index)
    SAFE_DEALLOCATE_A(this%ist)

    POP_SUB(batch_end)
  end subroutine batch_end

  !--------------------------------------------------------------
  subroutine batch_deallocate(this)
    class(batch_t),  intent(inout) :: this
    
    PUSH_SUB(batch_deallocate)

    this%is_allocated = .false.

    this%current = 1
    
    if(this%special_memory) then
      if(associated(this%dff)) then
        call deallocate_hardware_aware(c_loc(this%dff(1,1,1)))
      end if
      if(associated(this%zff)) then
        call deallocate_hardware_aware(c_loc(this%zff(1,1,1)))
      end if
    else
      SAFE_DEALLOCATE_P(this%dff)
      SAFE_DEALLOCATE_P(this%zff)
    end if
    nullify(this%dff)
    nullify(this%dff_linear)
    nullify(this%zff)
    nullify(this%zff_linear)
    
    POP_SUB(batch_deallocate)
  end subroutine batch_deallocate

  !--------------------------------------------------------------

  subroutine batch_deallocate_temporary(this)
    type(batch_t),  intent(inout) :: this
    
    PUSH_SUB(batch_deallocate_temporary)

    if(this%special_memory) then
      if(associated(this%dff)) then
        call deallocate_hardware_aware(c_loc(this%dff(1,1,1)))
        nullify(this%dff)
        nullify(this%dff_linear)
      end if
      if(associated(this%zff)) then
        call deallocate_hardware_aware(c_loc(this%zff(1,1,1)))
        nullify(this%zff)
        nullify(this%zff_linear)
      end if
    else
      SAFE_DEALLOCATE_P(this%dff)
      SAFE_DEALLOCATE_P(this%zff)
      nullify(this%dff_linear)
      nullify(this%zff_linear)
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
    
    PUSH_SUB(batch_init_empty)

    this%is_allocated = .false.
    this%mirror = .false.
    this%special_memory = .false.
    this%nst = nst
    this%dim = dim
    this%current = 1
    this%type_of = TYPE_NONE
    
    this%nst_linear = nst*dim

    this%max_size = 0
    this%in_buffer_count = 0
    this%status_of = BATCH_NOT_PACKED

    this%ndims = 2
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))
    SAFE_ALLOCATE(this%ist(1:this%nst))

    nullify(this%dff, this%zff, this%dff_linear, this%zff_linear)
    nullify(this%dff_pack, this%zff_pack)
    this%layout = BATCH_NOT_PACKED
    this%location = BATCH_HOST

    POP_SUB(batch_init_empty)
  end subroutine batch_init_empty

  !--------------------------------------------------------------
  logical function batch_is_ok(this) result(ok)
    class(batch_t), intent(in)   :: this

    ! no push_sub, called too frequently
    
    ok = this%nst_linear >= 1
    ok = ok .and. (associated(this%dff) .or.  associated(this%zff) .or. &
                   associated(this%dff_pack) .or.  associated(this%zff_pack))
    ok = ok .and. (this%type() == TYPE_CMPLX .or. this%type() == TYPE_FLOAT)
    if(ok .and. .not. this%is_packed()) then
      if(this%type() == TYPE_FLOAT) then
        ok = ok .and. associated(this%dff_linear)
        ok = ok .and. ubound(this%dff_linear, dim=2) == this%nst_linear
      else if(this%type() == TYPE_CMPLX) then
        ok = ok .and. associated(this%zff_linear)
        ok = ok .and. ubound(this%zff_linear, dim=2) == this%nst_linear
      else
        ok = .false.
      end if
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

    integer :: np

    PUSH_SUB(batch_copy_to)

    call batch_init_empty(dest, this%dim, this%nst)

    dest%type_of = this%type_of

    if(this%type() == TYPE_FLOAT) then

      np = ubound(this%dff_linear, dim=1)
      call dest%dallocate(1, this%nst, np)

    else if(this%type() == TYPE_CMPLX) then

      np = ubound(this%zff_linear, dim=1)
      call dest%zallocate(1, this%nst, np)

    else
      message(1) = "Internal error: unknown batch type in batch_copy_to."
      call messages_fatal(1)
    end if

    if(optional_default(pack, this%is_packed())) call dest%do_pack(copy = .false.)

    dest%ist_idim_index(1:this%nst_linear, 1:this%ndims) = this%ist_idim_index(1:this%nst_linear, 1:this%ndims)
    dest%ist(1:this%nst) = this%ist(1:this%nst)

    if(optional_default(copy_data, .false.)) call this%copy_data_to(np, dest)
    
    POP_SUB(batch_copy_to)
  end subroutine batch_copy_to

  ! ----------------------------------------------------
  !> THREADSAFE
  type(type_t) pure function batch_type(this) result(btype)
    class(batch_t),      intent(in)    :: this

    btype = this%type_of
     
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

  integer function batch_pack_total_size(this) result(size)
    class(batch_t),      intent(inout) :: this

    size = this%max_size
    if(accel_is_enabled()) size = accel_padded_size(size)
    size = size*pad_pow2(this%nst_linear)*types_get_size(this%type())

  end function batch_pack_total_size

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
      this%pack_size(1) = pad_pow2(this%nst_linear)
      this%pack_size(2) = this%max_size

      if(accel_is_enabled()) this%pack_size(2) = accel_padded_size(this%pack_size(2))

      this%pack_size_real = this%pack_size
      if(type_is_complex(this%type())) this%pack_size_real(1) = 2*this%pack_size_real(1)

      if(accel_is_enabled()) then
        this%status_of = BATCH_DEVICE_PACKED
        call accel_create_buffer(this%ff_device, ACCEL_MEM_READ_WRITE, this%type(), product(this%pack_size))
      else
        this%status_of = BATCH_PACKED
        this%layout = BATCH_PACKED
        ! always use hardware aware memory here
        if(this%type() == TYPE_FLOAT) then
          call c_f_pointer(dallocate_hardware_aware(this%pack_size(1)*this%pack_size(2)), this%dff_pack, this%pack_size)
        else if(this%type() == TYPE_CMPLX) then
          call c_f_pointer(zallocate_hardware_aware(this%pack_size(1)*this%pack_size(2)), this%zff_pack, this%pack_size)
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
        do sp = 1, this%pack_size(2), bsize
          ep = min(sp + bsize - 1, this%pack_size(2))
          forall(ist = 1:this%nst_linear)
            forall(ip = sp:ep)
              this%dff_pack(ist, ip) = this%dff_linear(ip, ist)
            end forall
          end forall
        end do

      else if(this%type() == TYPE_CMPLX) then

        bsize = hardware%zblock_size

        !$omp parallel do private(ep, ist, ip)
        do sp = 1, this%pack_size(2), bsize
          ep = min(sp + bsize - 1, this%pack_size(2))
          forall(ist = 1:this%nst_linear)
            forall(ip = sp:ep)
              this%zff_pack(ist, ip) = this%zff_linear(ip, ist)
            end forall
          end forall
        end do

      else
        message(1) = "Internal error: unknown batch type in batch_do_pack."
        call messages_fatal(1)
      end if

      call profiling_count_transfers(this%nst_linear*this%pack_size(2), this%type())
      
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
          call accel_release_buffer(this%ff_device)
        else
          if(associated(this%dff_pack)) then
            call deallocate_hardware_aware(c_loc(this%dff_pack(1,1)))
          end if
          if(associated(this%zff_pack)) then
            call deallocate_hardware_aware(c_loc(this%zff_pack(1,1)))
          end if
          nullify(this%dff_pack)
          nullify(this%zff_pack)
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

        !$omp parallel do private(ist)
        do ip = 1, this%pack_size(2)
          forall(ist = 1:this%nst_linear)
            this%dff_linear(ip, ist) = this%dff_pack(ist, ip)
          end forall
        end do
        
      else if(this%type() == TYPE_CMPLX) then

        !$omp parallel do private(ist)
        do ip = 1, this%pack_size(2)
          forall(ist = 1:this%nst_linear)
            this%zff_linear(ip, ist) = this%zff_pack(ist, ip)
          end forall
        end do
        
      else
        message(1) = "Internal error: unknown batch type in batch_sync."
        call messages_fatal(1)
      end if
      
      call profiling_count_transfers(this%nst_linear*this%pack_size(2), this%type())
      
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
        call accel_write_buffer(this%ff_device, ubound(this%dff_linear, dim=1), this%dff_linear(:, 1))
      else if(this%type() == TYPE_CMPLX) then
        call accel_write_buffer(this%ff_device, ubound(this%zff_linear, dim=1), this%zff_linear(:, 1))
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
      
      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack_size(1))

      call accel_create_buffer(tmp, ACCEL_MEM_READ_ONLY, this%type(), unroll*this%pack_size(2))
      
      do ist = 1, this%nst_linear, unroll
        
        ! copy a number 'unroll' of states to the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)

          if(this%type() == TYPE_FLOAT) then
            call accel_write_buffer(tmp, ubound(this%dff_linear, dim=1), this%dff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          else
            call accel_write_buffer(tmp, ubound(this%zff_linear, dim=1), this%zff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          end if
        end do

        ! now call an opencl kernel to rearrange the data
        call accel_set_kernel_arg(kernel, 0, this%pack_size(1))
        call accel_set_kernel_arg(kernel, 1, this%pack_size(2))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, tmp)
        call accel_set_kernel_arg(kernel, 4, this%ff_device)

        call profiling_in(prof_pack, "CL_PACK")
        call accel_kernel_run(kernel, (/this%pack_size(2), unroll/), (/accel_max_workgroup_size()/unroll, unroll/))

        if(this%type() == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack_size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack_size(2), M_ZI)
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
        call accel_read_buffer(this%ff_device, ubound(this%dff_linear, dim=1), this%dff_linear(:, 1))
      else
        call accel_read_buffer(this%ff_device, ubound(this%zff_linear, dim=1), this%zff_linear(:, 1))
      end if
    else

      unroll = min(CL_PACK_MAX_BUFFER_SIZE, this%pack_size(1))

      ! we use a kernel to move to a temporary array and then we read
      call accel_create_buffer(tmp, ACCEL_MEM_WRITE_ONLY, this%type(), unroll*this%pack_size(2))

      if(this%type() == TYPE_FLOAT) then
        kernel => dunpack
      else
        kernel => zunpack
      end if

      do ist = 1, this%nst_linear, unroll
        call accel_set_kernel_arg(kernel, 0, this%pack_size(1))
        call accel_set_kernel_arg(kernel, 1, this%pack_size(2))
        call accel_set_kernel_arg(kernel, 2, ist - 1)
        call accel_set_kernel_arg(kernel, 3, this%ff_device)
        call accel_set_kernel_arg(kernel, 4, tmp)

        call profiling_in(prof_unpack, "CL_UNPACK")
        call accel_kernel_run(kernel, (/unroll, this%pack_size(2)/), (/unroll, accel_max_workgroup_size()/unroll/))

        if(this%type() == TYPE_FLOAT) then
          call profiling_count_transfers(unroll*this%pack_size(2), M_ONE)
        else
          call profiling_count_transfers(unroll*this%pack_size(2), M_ZI)
        end if

        call accel_finish()
        call profiling_out(prof_unpack)

        ! copy a number 'unroll' of states from the buffer
        do ist2 = ist, min(ist + unroll - 1, this%nst_linear)
          
          if(this%type() == TYPE_FLOAT) then
            call accel_read_buffer(tmp, ubound(this%dff_linear, dim=1), this%dff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
          else
            call accel_read_buffer(tmp, ubound(this%zff_linear, dim=1), this%zff_linear(:, ist2), &
              offset = (ist2 - ist)*this%pack_size(2))
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
      call MPI_Win_create(this%zff_pack(1, 1), int(product(this%pack_size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
        types_get_size(this%type()), MPI_INFO_NULL, mpi_grp%comm, rma_win, mpi_err)
#endif
    else if (this%type() == TYPE_FLOAT) then
#ifdef HAVE_MPI2
      call MPI_Win_create(this%dff_pack(1, 1), int(product(this%pack_size)*types_get_size(this%type()), MPI_ADDRESS_KIND), &
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
    call accel_set_kernel_arg(kernel_copy, 1, this%ff_device)
    call accel_set_kernel_arg(kernel_copy, 2, log2(this%pack_size_real(1)))
    call accel_set_kernel_arg(kernel_copy, 3, dest%ff_device)
    call accel_set_kernel_arg(kernel_copy, 4, log2(dest%pack_size_real(1)))
    
    localsize = accel_kernel_workgroup_size(kernel_copy)/dest%pack_size_real(1)

    dim3 = np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(np, localsize))
    
    call accel_kernel_run(kernel_copy, (/dest%pack_size_real(1), dim2, dim3/), (/dest%pack_size_real(1), localsize, 1/))
    
    call accel_finish()

  case(BATCH_PACKED)
    if(dest%type() == TYPE_FLOAT) then
      call blas_copy(np*this%pack_size(1), this%dff_pack(1, 1), 1, dest%dff_pack(1, 1), 1)
    else
      call blas_copy(np*this%pack_size(1), this%zff_pack(1, 1), 1, dest%zff_pack(1, 1), 1)
    end if

  case(BATCH_NOT_PACKED)
    !$omp parallel do private(ist)
    do ist = 1, dest%nst_linear
      if(dest%type() == TYPE_CMPLX) then
        call blas_copy(np, this%zff_linear(1, ist), 1, dest%zff_linear(1, ist), 1)
      else
        call blas_copy(np, this%dff_linear(1, ist), 1, dest%dff_linear(1, ist), 1)
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
