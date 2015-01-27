#include "global.h"

module strng_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  implicit none

  private

  public ::       &
    operator(==), &
    operator(/=)

  public ::         &
    strng_init,    &
    strng_len,     &
    strng_hash,    &
    strng_get,     &
    strng_set,     &
    strng_append,  &
    strng_tolower, &
    strng_toupper, &
    strng_copy,    &
    strng_end

  real(kind=wp), parameter :: STRNG_GROWTH_FACTOR = 1.1_wp
  integer,       parameter :: STRNG_INIT_LEN      = 63

  integer, public, parameter :: STRNG_OK          = 0
  integer, public, parameter :: STRNG_SIZE_ERROR  =-1
  integer, public, parameter :: STRNG_INDEX_ERROR =-2

  character(len=*), parameter :: ALPH_LOWER="abcdefghijklmnopqrstuvwxyz"
  character(len=*), parameter :: ALPH_UPPER="ABCDEFGHIJKLMNOPQRSTUVWXYZ"

  type, public :: strng_t
    private
    integer                          :: len  = 0
    integer                          :: size = 0
    character, dimension(:), pointer :: value =>null()
  end type strng_t

  type, public :: strng_iterator_t
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

  interface strng_init
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
    !
    logical :: is
    !
    is=((this%size>0).and.associated(this%value))
    return
  end function strng_isdef

  ! ---------------------------------------------------------
  elemental function strng_len(this) result(len)
    type(strng_t), intent(in) :: this
    !
    integer :: len
    !
    len=this%len
    return
  end function strng_len

  ! ---------------------------------------------------------
  elemental function strng_equal(this, that) result(eqv)
    type(strng_t), intent(in) :: this
    type(strng_t), intent(in) :: that
    !
    logical :: eqv
    !
    integer :: i
    !
    eqv=.false.
    if(strng_isdef(this).and.strng_isdef(that))then
      if(this%len==that%len)then
        eqv=.true.
        do i = 1, this%len
          if(this%value(i)/=that%value(i))then
            eqv=.false.
            exit
          end if
        end do
      end if
    end if
    return
  end function strng_equal

  ! ---------------------------------------------------------
  elemental function strng_not_equal(this, that) result(eqv)
    type(strng_t), intent(in) :: this
    type(strng_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(.not.strng_equal(this, that))
    return
  end function strng_not_equal

  ! ---------------------------------------------------------
  !Daniel J. Bernstein Hash Function
  ! ---------------------------------------------------------
  elemental function strng_hash(this, size) result(hash)
    type(strng_t),    intent(in) :: this
    integer, optional, intent(in) :: size
    !
    integer :: hash
    !
    integer :: i
    !
    hash=5381
    do i = 1, this%len
      hash = ieor(33*hash, iachar(this%value(i)))
    end do
    if(present(size))hash=modulo(hash, size)+1
    return
  end function strng_hash

  ! ---------------------------------------------------------
  elemental subroutine strng_nullify(this)
    type(strng_t), intent(inout) :: this
    !
    nullify(this%value)
    this%size=0
    this%len=0
    return
  end subroutine strng_nullify

  ! ---------------------------------------------------------
  subroutine strng_allocate(this, size)
    type(strng_t),    intent(out) :: this
    integer, optional, intent(in)  :: size
    !
    integer :: n
    !
    PUSH_SUB(strng_init_size)
    call strng_nullify(this)
    this%size=STRNG_INIT_LEN
    if(present(size))then
      if(STRNG_INIT_LEN<size)then
        n=max(ceiling((log(real(size,kind=wp))-log(real(STRNG_INIT_LEN,kind=wp)))/log(STRNG_GROWTH_FACTOR)),1)
        this%size=ceiling((STRNG_GROWTH_FACTOR**n)*real(STRNG_INIT_LEN,kind=wp))
      end if
    end if
    SAFE_ALLOCATE(this%value(this%size))
    POP_SUB(strng_init_size)
    return
  end subroutine strng_allocate

  ! ---------------------------------------------------------
  subroutine strng_reallocate(this, extra)
    type(strng_t), intent(inout) :: this
    integer,        intent(in)    :: extra
    !
    type(strng_t) :: buff
    integer        :: size
    !
    PUSH_SUB(strng_reallocate)
    if(strng_isdef(this))then
      size=int(STRNG_GROWTH_FACTOR*real(this%len+extra,kind=wp))
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
    return
  end subroutine strng_reallocate

  ! ---------------------------------------------------------
  subroutine strng_init_char(this, value)
    type(strng_t),             intent(out) :: this
    character(len=*), optional, intent(in)  :: value
    !
    PUSH_SUB(strng_init_char)
    if(present(value))then
      call strng_allocate(this, len_trim(adjustl(value)))
      call strng_append_char(this,trim(adjustl(value)))
    else
      call strng_allocate(this)
    end if
    POP_SUB(strng_init_char)
    return
  end subroutine strng_init_char

  ! ---------------------------------------------------------
  subroutine strng_init_strng(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that
    !
    PUSH_SUB(strng_init_strng)
    call strng_copy(this, that)
    POP_SUB(strng_init_strng)
    return
  end subroutine strng_init_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_pop(this, char)
    type(strng_t), intent(inout) :: this
    character,      intent(out)   :: char
    !
    char=""
    if(this%len>0)then
      char=this%value(this%len)
      this%value(this%len)=""
      this%len=this%len-1
    end if
    return
  end subroutine strng_pop

  ! ---------------------------------------------------------
  elemental subroutine strng_get_char(this, i, value, ierr)
    type(strng_t),    intent(in)  :: this
    integer,           intent(in)  :: i
    character,         intent(out) :: value
    integer, optional, intent(out) :: ierr
    !
    value=""
    if(present(ierr))ierr=STRNG_INDEX_ERROR
    if(strng_isdef(this))then
      if((0<i).and.(i<=this%len))then
        value=this%value(i)
        if(present(ierr))ierr=STRNG_OK
      end if
    end if
    return
  end subroutine strng_get_char

  ! ---------------------------------------------------------
  elemental subroutine strng_get_strng(this, value, ierr)
    type(strng_t),    intent(in)  :: this
    character(len=*),  intent(out) :: value
    integer, optional, intent(out) :: ierr
    !
    integer :: i
    !
    value=""
    if(present(ierr))ierr=STRNG_SIZE_ERROR
    if(strng_isdef(this))then
      if(this%len<=len(value))then
        do i = 1, this%len
          value(i:i)=this%value(i)
        end do
        if(present(ierr))ierr=STRNG_OK
      end if
    end if
    return
  end subroutine strng_get_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_set(this, i, char, ierr)
    type(strng_t),    intent(inout) :: this
    integer,           intent(in)    :: i
    character,         intent(in)    :: char
    integer, optional, intent(out)   :: ierr
    !
    if(present(ierr))ierr=STRNG_INDEX_ERROR
    if(strng_isdef(this))then
      if((0<i).and.(i<=this%len))then
        this%value(i)=char
        if(present(ierr))ierr=STRNG_OK
      end if
    end if
    return
  end subroutine strng_set

  ! ---------------------------------------------------------
  subroutine strng_append_char(this, char)
    type(strng_t),   intent(inout) :: this
    character(len=*), intent(in)    :: char
    !
    integer :: i
    !
    PUSH_SUB(strng_append_char)
    call strng_reallocate(this, len(char))
    do i = 1, len(char)
      this%value(this%len+i)=char(i:i)
    end do
    this%len=this%len+len(char)
    POP_SUB(strng_append_char)
    return
  end subroutine strng_append_char

  ! ---------------------------------------------------------
  subroutine strng_append_strng(this, strng)
    type(strng_t), intent(inout) :: this
    type(strng_t), intent(in)    :: strng
    !
    PUSH_SUB(strng_append_strng)
    call strng_reallocate(this, strng%len)
    this%value(this%len+1:this%len+strng%len)=strng%value(1:strng%len)
    this%len=this%len+strng%len
    POP_SUB(strng_append_strng)
    return
  end subroutine strng_append_strng

  ! ---------------------------------------------------------
  elemental subroutine strng_tolower(this)
    type(strng_t),   intent(inout) :: this
    !
    integer :: i, p
    !
    do i = 1, this%len
      p=scan(ALPH_UPPER, this%value(i))
      if(p>0)this%value(i)=ALPH_LOWER(p:p)
    end do
    return
  end subroutine strng_tolower

  ! ---------------------------------------------------------
  elemental subroutine strng_toupper(this)
    type(strng_t),   intent(inout) :: this
    !
    integer :: i, p
    !
    do i = 1, this%len
      p=scan(ALPH_LOWER, this%value(i))
      if(p>0)this%value(i)=ALPH_UPPER(p:p)
    end do
    return
  end subroutine strng_toupper

  ! ---------------------------------------------------------
  subroutine strng_shallow_copy(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that
    !
    PUSH_SUB(strng_shallow_copy)
    this%len=that%len
    this%size=that%size
    this%value=>that%value
    POP_SUB(strng_shallow_copy)
    return
  end subroutine strng_shallow_copy

  ! ---------------------------------------------------------
  subroutine strng_copy(this, that)
    type(strng_t), intent(out) :: this
    type(strng_t), intent(in)  :: that
    !
    PUSH_SUB(strng_copy)
    call strng_allocate(this, that%size)
    call strng_append_strng(this, that)
    POP_SUB(strng_copy)
    return
  end subroutine strng_copy

  ! ---------------------------------------------------------
  subroutine strng_end(this)
    type(strng_t), intent(inout) :: this
    !
    PUSH_SUB(strng_end)
    SAFE_DEALLOCATE_P(this%value)
    nullify(this%value)
    this%size=0
    this%len=0
    POP_SUB(strng_end)
    return
  end subroutine strng_end

end module strng_m

!! Local Variables:
!! mode: f90
!! End:
