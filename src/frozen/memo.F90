#include "global.h"

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

#define DICT_TEMPLATE_NAME data

module memo_oct_m

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX

  implicit none

  private

  public ::        &
    MEMO_NAME_LEN

  public ::           &
    MEMO_OK,          &
    MEMO_KEY_ERROR,   &
    MEMO_EMPTY_ERROR

  public :: &
    memo_t

  public ::    &
    memo_in,   &
    memo_init, &
    memo_set,  &
    memo_get,  &
    memo_del,  &
    memo_copy, &
    memo_end

#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER

  integer, parameter :: MEMO_NAME_LEN = 63

  integer, parameter :: MEMO_OK          = DATA_DICT_OK
  integer, parameter :: MEMO_KEY_ERROR   = DATA_DICT_KEY_ERROR
  integer, parameter :: MEMO_EMPTY_ERROR = DATA_DICT_EMPTY_ERROR

  integer, parameter :: DATA_STAT_NULL = 0
  integer, parameter :: DATA_STAT_INIT = 1

  type :: data_t
    private
    integer                                  :: stat = DATA_STAT_NULL
    real(kind=wp), dimension(:), allocatable :: data
  end type data_t

  type :: memo_t
    private
    type(data_dict_t) :: dict
  end type memo_t

  interface data_new
    module procedure data_new_type
    module procedure data_new_copy
  end interface data_new

  interface data_init
    module procedure data_init_type
    module procedure data_init_copy
  end interface data_init

  interface memo_set
    module procedure memo_set_r0
    module procedure memo_set_r1
  end interface memo_set

  interface memo_get
    module procedure memo_get_r0
    module procedure memo_get_r1
  end interface memo_get

  interface memo_del
    module procedure memo_del_nv
    module procedure memo_del_r0
    module procedure memo_del_r1
  end interface memo_del

contains

#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY

  ! ---------------------------------------------------------
  function data_new_type(data) result(this)
    real(kind=wp), dimension(:), intent(in) :: data
    
    type(data_t), pointer :: this

    PUSH_SUB(data_new_type)

    nullify(this)
    SAFE_ALLOCATE(this)
    call data_init(this, data)

    POP_SUB(data_new_type)
  end function data_new_type

  ! ---------------------------------------------------------
  function data_new_copy(that) result(this)
    type(data_t), intent(in) :: that
    
    type(data_t), pointer :: this

    PUSH_SUB(data_new_copy)

    ASSERT(data_alloc(that))
    nullify(this)
    SAFE_ALLOCATE(this)
    call data_init(this, that)

    POP_SUB(data_new_copy)
  end function data_new_copy

  ! ---------------------------------------------------------
  subroutine data_del(this)
    type(data_t), pointer, intent(inout) :: this

    PUSH_SUB(data_del)

    if(associated(this))then
      call data_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(data_del)
  end subroutine data_del

  ! ---------------------------------------------------------
  function data_alloc(this) result(that)
    type(data_t), intent(in) :: this
    
    logical :: that

    PUSH_SUB(data_alloc)

    that = .false.
    select case(this%stat)
    case(DATA_STAT_NULL)
      ASSERT(.not.allocated(this%data))
    case(DATA_STAT_INIT)
      ASSERT(allocated(this%data))
      that = .true.
    case default
      ASSERT(.false.)
    end select

    POP_SUB(data_alloc)
  end function data_alloc

  ! ---------------------------------------------------------
  subroutine data_init_type(this, data)
    type(data_t),                intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: data

    PUSH_SUB(data_init_type)

    ASSERT(size(data)>0)
    this%stat = DATA_STAT_INIT
    SAFE_ALLOCATE(this%data(size(data)))
    this%data(:) = data(:)

    POP_SUB(data_init_type)
  end subroutine data_init_type

  ! ---------------------------------------------------------
  subroutine data_init_copy(this, that)
    type(data_t), intent(out) :: this
    type(data_t), intent(in)  :: that

    PUSH_SUB(data_init_copy)

    ASSERT(data_alloc(that))
    call data_copy(this, that)

    POP_SUB(data_init_copy)
  end subroutine data_init_copy

  ! ---------------------------------------------------------
  subroutine data_get(this, data)
    type(data_t),                 target, intent(in)  :: this
    real(kind=wp), dimension(:), pointer, intent(out) :: data

    PUSH_SUB(data_get)

    ASSERT(data_alloc(this))
    data => this%data

    POP_SUB(data_get)
  end subroutine data_get

  ! ---------------------------------------------------------
  subroutine data_copy(this, that)
    type(data_t), intent(inout) :: this
    type(data_t), intent(in)    :: that

    PUSH_SUB(data_copy)

    call data_end(this)
    select case(that%stat)
    case(DATA_STAT_NULL)
      ASSERT(.not.data_alloc(that))
    case(DATA_STAT_INIT)
      ASSERT(data_alloc(that))
      call data_init(this, that%data)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(data_copy)
  end subroutine data_copy

  ! ---------------------------------------------------------
  subroutine data_end(this)
    type(data_t), intent(inout) :: this

    PUSH_SUB(data_end)

    select case(this%stat)
    case(DATA_STAT_NULL)
      ASSERT(.not.data_alloc(this))
    case(DATA_STAT_INIT)
      ASSERT(data_alloc(this))
      SAFE_DEALLOCATE_A(this%data)
    case default
      ASSERT(.false.)
    end select
    this%stat = DATA_STAT_NULL
    
    POP_SUB(data_end)
  end subroutine data_end

  ! ---------------------------------------------------------
  subroutine memo_init(this, size)
    type(memo_t),      intent(out) :: this
    integer, optional, intent(in)  :: size

    PUSH_SUB(memo_init)

    call data_dict_init(this%dict, size)

    POP_SUB(memo_init)
  end subroutine memo_init

  ! ---------------------------------------------------------
  function memo_in(this, key) result(has)
    type(memo_t),     intent(in) :: this
    character(len=*), intent(in) :: key

    logical :: has

    PUSH_SUB(memo_in)

    has = data_dict_has_key(this%dict, key)
    
    POP_SUB(memo_in)
  end function memo_in

  ! ---------------------------------------------------------
  subroutine memo_set_r0(this, key, val)
    type(memo_t),     intent(inout) :: this
    character(len=*), intent(in)    :: key
    real(kind=wp),    intent(in)    :: val

    real(kind=wp), dimension(1) :: ival
    type(data_t),       pointer :: data
    
    PUSH_SUB(memo_set_r0)

    ival(1) = val
    data => data_new(ival)
    call data_dict_set(this%dict, key, data)
    nullify(data)

    POP_SUB(memo_set_r0)
  end subroutine memo_set_r0

  ! ---------------------------------------------------------
  subroutine memo_set_r1(this, key, val)
    type(memo_t),                intent(inout) :: this
    character(len=*),            intent(in)    :: key
    real(kind=wp), dimension(:), intent(in)    :: val

    type(data_t), pointer :: data
    
    PUSH_SUB(memo_set_r1)

    data => data_new(val)
    call data_dict_set(this%dict, key, data)
    nullify(data)

    POP_SUB(memo_set_r1)
  end subroutine memo_set_r1

  ! ---------------------------------------------------------
  subroutine memo_get_r0(this, key, val, ierr)
    type(memo_t),      intent(in)  :: this
    character(len=*),  intent(in)  :: key
    real(kind=wp),     intent(out) :: val
    integer, optional, intent(out) :: ierr

    type(data_t),                pointer :: data
    real(kind=wp), dimension(:), pointer :: ival
    integer                              :: jerr

    PUSH_SUB(memo_get_r0)

    nullify(data, ival)
    call data_dict_get(this%dict, key, data, jerr)
    if(jerr==MEMO_OK)then
      ASSERT(associated(data))
      call data_get(data, ival)
      ASSERT(associated(ival))
      ASSERT(size(ival)==1)
      val = ival(1)
      nullify(ival)
    end if
    nullify(data)
    if(present(ierr)) ierr = jerr
    
    POP_SUB(memo_get_r0)
  end subroutine memo_get_r0

  ! ---------------------------------------------------------
  subroutine memo_get_r1(this, key, val, ierr)
    type(memo_t),                intent(in)  :: this
    character(len=*),            intent(in)  :: key
    real(kind=wp), dimension(:), intent(out) :: val
    integer,           optional, intent(out) :: ierr

    type(data_t),                pointer :: data
    real(kind=wp), dimension(:), pointer :: ival
    integer                              :: jerr

    PUSH_SUB(memo_get_r1)

    nullify(data, ival)
    call data_dict_get(this%dict, key, data, jerr)
    if(jerr==MEMO_OK)then
      ASSERT(associated(data))
      call data_get(data, ival)
      ASSERT(associated(ival))
      ASSERT(size(ival)==size(val))
      val = ival
      nullify(ival)
    end if
    nullify(data)
    if(present(ierr)) ierr = jerr
    
    POP_SUB(memo_get_r1)
  end subroutine memo_get_r1

  ! ---------------------------------------------------------
  subroutine memo_del_nv(this, key, ierr)
    type(memo_t),      intent(inout) :: this
    character(len=*),  intent(in)    :: key
    integer, optional, intent(out)   :: ierr

    type(data_t), pointer :: data
    integer               :: jerr

    PUSH_SUB(memo_del_nv)

    nullify(data)
    call data_dict_del(this%dict, key, data, jerr)
    if(jerr==MEMO_OK)then
      ASSERT(associated(data))
      call data_del(data)
    end if
    nullify(data)
    if(present(ierr)) ierr = jerr

    POP_SUB(memo_del_nv)
  end subroutine memo_del_nv

  ! ---------------------------------------------------------
  subroutine memo_del_r0(this, key, val, ierr)
    type(memo_t),      intent(inout) :: this
    character(len=*),  intent(in)    :: key
    real(kind=wp),     intent(out)   :: val
    integer, optional, intent(out)   :: ierr

    type(data_t),                pointer :: data
    real(kind=wp), dimension(:), pointer :: ival
    integer                              :: jerr

    PUSH_SUB(memo_del)

    nullify(data, ival)
    call data_dict_del(this%dict, key, data, jerr)
    if(jerr==MEMO_OK)then
      ASSERT(associated(data))
      call data_get(data, ival)
      ASSERT(associated(ival))
      ASSERT(size(ival)==1)
      val = ival(1)
      nullify(ival)
      call data_del(data)
    end if
    nullify(data)
    if(present(ierr)) ierr = jerr

    POP_SUB(memo_del_r0)
  end subroutine memo_del_r0

  ! ---------------------------------------------------------
  subroutine memo_del_r1(this, key, val, ierr)
    type(memo_t),                intent(inout) :: this
    character(len=*),            intent(in)    :: key
    real(kind=wp), dimension(:), intent(out)   :: val
    integer,           optional, intent(out)   :: ierr

    type(data_t),                pointer :: data
    real(kind=wp), dimension(:), pointer :: ival
    integer                              :: jerr

    PUSH_SUB(memo_del_r1)

    nullify(data, ival)
    call data_dict_del(this%dict, key, data, jerr)
    if(jerr==MEMO_OK)then
      ASSERT(associated(data))
      call data_get(data, ival)
      ASSERT(associated(ival))
      ASSERT(size(ival)==size(val))
      val = ival
      nullify(ival)
      call data_del(data)
    end if
    nullify(data)
    if(present(ierr)) ierr = jerr

    POP_SUB(memo_del_r1)
  end subroutine memo_del_r1

  ! ---------------------------------------------------------
  subroutine memo_copy(this, that)
    type(memo_t), intent(inout) :: this
    type(memo_t), intent(in)    :: that

    type(data_dict_iterator_t)   :: iter
    character(len=MEMO_NAME_LEN) :: name
    type(data_t),        pointer :: data
    integer                      :: ierr

    PUSH_SUB(memo_copy)

    call memo_end(this)
    call data_dict_init(iter, that%dict)
    do
      nullify(data)
      call data_dict_next(iter, name, data, ierr)
      if(ierr/=MEMO_OK)exit
      ASSERT(associated(data))
      call data_dict_set(this%dict, trim(adjustl(name)), data_new(data))
    end do
    call data_dict_end(iter)
    nullify(data)
    
    POP_SUB(memo_copy)
  end subroutine memo_copy

  ! ---------------------------------------------------------
  subroutine memo_end(this)
    type(memo_t), intent(inout) :: this

    type(data_t), pointer :: data
    integer               :: ierr

    PUSH_SUB(memo_end)

    do
      nullify(data)
      call data_dict_pop(this%dict, data, ierr)
      if(ierr/=MEMO_OK)exit
      ASSERT(associated(data))
      call data_del(data)
    end do
    nullify(data)
    ASSERT(data_dict_len(this%dict)==0)
    call data_dict_end(this%dict)
    
    POP_SUB(memo_end)
  end subroutine memo_end

end module memo_oct_m

#undef DICT_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

