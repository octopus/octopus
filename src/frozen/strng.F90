#include "global.h"

module strng_oct_m

  use global_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m


  implicit none

  private

  public ::            &
    STRNG_OK,          &
    STRNG_SIZE_ERROR,  &
    STRNG_INDEX_ERROR

  public ::  &
    strng_t

  public ::           &
    strng_iterator_t

  public ::       &
    operator(==), &
    operator(/=)

  public ::        &
    strng_new,     &
    strng_del,     &
    strng_init,    &
    strng_len,     &
    strng_hash,    &
    strng_set,     &
    strng_get,     &
    strng_pop,     &
    strng_append,  &
    strng_tolower, &
    strng_toupper, &
    strng_copy,    &
    strng_end

  real(kind=wp), parameter :: STRNG_GROWTH_FACTOR = 1.1_wp
  integer,       parameter :: STRNG_INIT_LEN      = 63

  integer, parameter :: STRNG_OK          = 0
  integer, parameter :: STRNG_SIZE_ERROR  =-1
  integer, parameter :: STRNG_INDEX_ERROR =-2

  character(len=*), parameter :: ALPH_LOWER="abcdefghijklmnopqrstuvwxyz"
  character(len=*), parameter :: ALPH_UPPER="ABCDEFGHIJKLMNOPQRSTUVWXYZ"

  type :: strng_t
    private
    integer                          :: len  = 0
    integer                          :: size = 0
    character, dimension(:), pointer :: value =>null()
  end type strng_t

  type :: strng_iterator_t
    private
    integer                          :: pos   = 0
    integer                          :: len   = 0
    character, dimension(:), pointer :: value =>null()
  end type strng_iterator_t

  interface operator(==)
    module procedure strng_equal
  end interface operator(==)

  interface operator(/=)
    module procedure strng_not_equal
  end interface operator(/=)

  interface strng_new
    module procedure strng_new_null
    module procedure strng_new_char
    module procedure strng_new_strng
  end interface strng_new

  interface strng_init
    module procedure strng_init_null
    module procedure strng_init_char
    module procedure strng_init_strng
  end interface strng_init

  interface strng_get
    module procedure strng_get_char
    module procedure strng_get_strng
  end interface strng_get

  interface strng_append
    module procedure strng_append_char
    module procedure strng_append_strng
  end interface strng_append

contains

  ! ---------------------------------------------------------
  elemental function strng_isdef(this) result(is)
    type(strng_t), intent(in) :: this

    logical :: is

    is = ((this%size>0).and.associated(this%value))

  end function strng_isdef

  ! ---------------------------------------------------------
  elemental function strng_len(this) result(len)
    type(strng_t), intent(in) :: this

    integer :: len

    len = this%len

  end function strng_len

  ! ---------------------------------------------------------
  elemental function strng_equal(this, that) result(eqv)
    type(strng_t), intent(in) :: this
    type(strng_t), intent(in) :: that

    logical :: eqv

    integer :: i

    eqv = .false.
    if(strng_isdef(this).and.strng_isdef(that))then
      if(this%len==that%len)then
        eqv = .true.
        do i = 1, this%len
          if(this%value(i)/=that%value(i))then
            eqv = .false.
            exit
          end if
        end do
      end if
    end if

  end function strng_equal

  ! ---------------------------------------------------------
  elemental function strng_not_equal(this, that) result(eqv)
    type(strng_t), intent(in) :: this
    type(strng_t), intent(in) :: that

    logical :: eqv

    eqv=(.not.strng_equal(this, that))

  end function strng_not_equal

  ! ---------------------------------------------------------
  !Daniel J. Bernstein Hash Function
  ! ---------------------------------------------------------
  elemental function strng_hash(this, size) result(hash)
    type(strng_t),    intent(in) :: this
    integer, optional, intent(in) :: size

    integer :: hash

    integer :: indx

    hash = 5381
    do indx = 1, this%len
      hash = ieor(33*hash, iachar(this%value(indx)))
    end do
    if(present(size)) hash = modulo(hash, size) + 1

  end function strng_hash

  ! ---------------------------------------------------------
  elemental subroutine strng_nullify(this)
    type(strng_t), intent(inout) :: this

    nullify(this%value)
    this%size = 0
    this%len = 0

  end subroutine strng_nullify

  ! ---------------------------------------------------------
  subroutine strng_allocate(this, size)
    type(strng_t),     intent(out) :: this
    integer, optional, intent(in)  :: size

    real(kind=wp) :: tmp
    integer       :: pwr

    PUSH_SUB(strng_allocate)

    call strng_nullify(this)
    this%size = STRNG_INIT_LEN
    if(present(size))then
      if(STRNG_INIT_LEN<size)then
        tmp = log(real(size,kind=wp)) - log(real(STRNG_INIT_LEN,kind=wp))
        tmp = tmp / log(STRNG_GROWTH_FACTOR)
        pwr = max(ceiling(tmp), 1)
        this%size = ceiling((STRNG_GROWTH_FACTOR**pwr)*real(STRNG_INIT_LEN,kind=wp))
      end if
    end if
    SAFE_ALLOCATE(this%value(this%size))

    POP_SUB(strng_allocate)
  end subroutine strng_allocate

  ! ---------------------------------------------------------
  subroutine strng_reallocate(this, extra)
    type(strng_t), intent(inout) :: this
    integer,       intent(in)    :: extra

    type(strng_t) :: buff
    integer       :: size

    PUSH_SUB(strng_reallocate)

    if(strng_isdef(this))then
      size = int(STRNG_GROWTH_FACTOR * real(this%len+extra,kind=wp))
      if(this%size<size)then
        call strng_shallow_copy(buff, this)
        call strng_allocate(this, size)
        ASSERT(buff%size<this%size)
        call strng_append_strng(this, buff)
        call strng_nullify(buff)
      end if
    else
      call strng_allocate(this)
    end if

    POP_SUB(strng_reallocate)
  end subroutine strng_reallocate

  ! ---------------------------------------------------------
  subroutine strng_inew(this)
    type(strng_t), pointer :: this

    PUSH_SUB(strng_inew)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(strng_inew)
  end subroutine strng_inew

  ! ---------------------------------------------------------
  subroutine strng_idel(this)
    type(strng_t), pointer :: this

    PUSH_SUB(strng_idel)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(strng_idel)
  end subroutine strng_idel

  ! ---------------------------------------------------------
  subroutine strng_new_null(this)
    type(strng_t), pointer :: this

    PUSH_SUB(strng_new_null)

    call strng_inew(this)
    call strng_init(this)

    POP_SUB(strng_new_null)
  end subroutine strng_new_null

  ! ---------------------------------------------------------
  subroutine strng_new_char(this, value)
    type(strng_t),   pointer     :: this
    character(len=*), intent(in) :: value

    PUSH_SUB(strng_new_char)

    call strng_inew(this)
    call strng_init(this, trim(adjustl(value)))

    POP_SUB(strng_new_char)
  end subroutine strng_new_char

  ! ---------------------------------------------------------
  subroutine strng_new_strng(this, that)
    type(strng_t), pointer     :: this
    type(strng_t),  intent(in) :: that

    PUSH_SUB(strng_new_strng)

    call strng_inew(this)
    call strng_init(this, that)

    POP_SUB(strng_new_strng)
  end subroutine strng_new_strng

  ! ---------------------------------------------------------
  subroutine strng_del(this)
    type(strng_t), pointer :: this

    PUSH_SUB(strng_del)

    call strng_end(this)
    call strng_idel(this)

    POP_SUB(strng_del)
  end subroutine strng_del

  ! ---------------------------------------------------------
  subroutine strng_init_null(this)
    type(strng_t), intent(out) :: this

    PUSH_SUB(strng_init_null)

    call strng_allocate(this)

    POP_SUB(strng_init_null)
  end subroutine strng_init_null

  ! ---------------------------------------------------------
  subroutine strng_init_char(this, value)
    type(strng_t),    intent(out) :: this
    character(len=*), intent(in)  :: value

    PUSH_SUB(strng_init_char)

    call strng_allocate(this, len_trim(adjustl(value)))
    call strng_append_char(this, trim(adjustl(value)))

    POP_SUB(strng_init_char)
  end subroutine strng_init_char

  ! ---------------------------------------------------------
  subroutine strng_init_strng(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that

    PUSH_SUB(strng_init_strng)

    call strng_copy(this, that)

    POP_SUB(strng_init_strng)
  end subroutine strng_init_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_set(this, i, char, ierr)
    type(strng_t),     intent(inout) :: this
    integer,           intent(in)    :: i
    character,         intent(in)    :: char
    integer, optional, intent(out)   :: ierr

    if(present(ierr)) ierr = STRNG_INDEX_ERROR
    if(strng_isdef(this))then
      if((0<i).and.(i<=this%len))then
        this%value(i) = char
        if(present(ierr)) ierr = STRNG_OK
      end if
    end if

  end subroutine strng_set

  ! ---------------------------------------------------------
  elemental subroutine strng_get_char(this, i, value, ierr)
    type(strng_t),     intent(in)  :: this
    integer,           intent(in)  :: i
    character,         intent(out) :: value
    integer, optional, intent(out) :: ierr

    value = ""
    if(present(ierr)) ierr = STRNG_INDEX_ERROR
    if(strng_isdef(this))then
      if((0<i).and.(i<=this%len))then
        value = this%value(i)
        if(present(ierr)) ierr = STRNG_OK
      end if
    end if

  end subroutine strng_get_char

  ! ---------------------------------------------------------
  !elemental
  pure subroutine strng_get_strng(this, value, ierr)
    type(strng_t),     intent(in)  :: this
    character(len=*),  intent(out) :: value
    integer, optional, intent(out) :: ierr

    integer :: idx

    value = ""
    if(present(ierr)) ierr = STRNG_SIZE_ERROR
    if(strng_isdef(this))then
      if(this%len<=len(value))then
        do idx = 1, this%len
          value(idx:idx) = this%value(idx)
        end do
        if(present(ierr)) ierr = STRNG_OK
      end if
    end if

  end subroutine strng_get_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_pop(this, char)
    type(strng_t), intent(inout) :: this
    character,     intent(out)   :: char

    char = ""
    if(this%len>0)then
      char = this%value(this%len)
      this%value(this%len) = ""
      this%len = this%len - 1
    end if

  end subroutine strng_pop

  ! ---------------------------------------------------------
  subroutine strng_append_char(this, char)
    type(strng_t),    intent(inout) :: this
    character(len=*), intent(in)    :: char

    integer :: idx

    PUSH_SUB(strng_append_char)

    call strng_reallocate(this, len(char))
    do idx = 1, len(char)
      this%value(this%len+idx) = char(idx:idx)
    end do
    this%len = this%len + len(char)

    POP_SUB(strng_append_char)
  end subroutine strng_append_char

  ! ---------------------------------------------------------
  subroutine strng_append_strng(this, strng)
    type(strng_t), intent(inout) :: this
    type(strng_t), intent(in)    :: strng

    PUSH_SUB(strng_append_strng)

    call strng_reallocate(this, strng%len)
    this%value(this%len+1:this%len+strng%len) = strng%value(1:strng%len)
    this%len = this%len + strng%len

    POP_SUB(strng_append_strng)
  end subroutine strng_append_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_tolower(this)
    type(strng_t),   intent(inout) :: this

    integer :: idx, pos

    do idx = 1, this%len
      pos = scan(ALPH_UPPER, this%value(idx))
      if(pos>0) this%value(idx) = ALPH_LOWER(pos:pos)
    end do

  end subroutine strng_tolower

  ! ---------------------------------------------------------
  elemental subroutine strng_toupper(this)
    type(strng_t),   intent(inout) :: this

    integer :: idx, pos

    do idx = 1, this%len
      pos = scan(ALPH_LOWER, this%value(idx))
      if(pos>0) this%value(idx) = ALPH_UPPER(pos:pos)
    end do

  end subroutine strng_toupper

  ! ---------------------------------------------------------
  subroutine strng_shallow_copy(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that

    PUSH_SUB(strng_shallow_copy)

    this%len = that%len
    this%size = that%size
    this%value => that%value

    POP_SUB(strng_shallow_copy)
  end subroutine strng_shallow_copy

  ! ---------------------------------------------------------
  subroutine strng_copy(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that

    PUSH_SUB(strng_copy)

    call strng_allocate(this, that%size)
    call strng_append_strng(this, that)

    POP_SUB(strng_copy)
  end subroutine strng_copy

  ! ---------------------------------------------------------
  subroutine strng_end(this)
    type(strng_t), intent(inout) :: this

    PUSH_SUB(strng_end)

    SAFE_DEALLOCATE_P(this%value)
    nullify(this%value)
    this%size = 0
    this%len = 0

    POP_SUB(strng_end)
  end subroutine strng_end

end module strng_oct_m

!! Local Variables:
!! mode: f90
!! End:
