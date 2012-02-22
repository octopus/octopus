#include "global.h"

module json_m

  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public :: json_isdef
  public :: json_init
  public :: json_end

  public :: scan
  public :: json_len
  public :: json_string
  public :: json_write
  public :: json_get
  public :: json_set
  public :: json_append
  public :: json_next
  public :: operator(==)

  FLOAT, parameter :: kind_parm=CNST(1.0)

  integer, public, parameter :: wp=kind(kind_parm)

  real(kind=wp), parameter :: JSON_STRING_GROWTH_FACTOR = 1.1_wp
  integer,       parameter :: JSON_STRING_INIT_LEN      = 63

  real(kind=wp), parameter :: JSON_TABLE_GROWTH_FACTOR = 1.5_wp
  integer,       parameter :: JSON_TABLE_INIT_LEN=127

  integer, public, parameter :: JSON_OK          = 0
  integer, public, parameter :: JSON_UNDEF_ERROR =-1
  integer, public, parameter :: JSON_TYPE_ERROR  =-2
  integer, public, parameter :: JSON_SIZE_ERROR  =-3
  integer, public, parameter :: JSON_INDEX_ERROR =-4
  integer, public, parameter :: JSON_KEY_ERROR   =-5

  integer, parameter :: JSON_UNDEF_TYPE   =-1
  integer, parameter :: JSON_NULL_TYPE    = 0
  integer, parameter :: JSON_LOGICAL_TYPE = 1
  integer, parameter :: JSON_INTEGER_TYPE = 2
  integer, parameter :: JSON_REAL_TYPE    = 3
  integer, parameter :: JSON_STRING_TYPE  = 4
  integer, parameter :: JSON_ARRAY_TYPE   = 5
  integer, parameter :: JSON_OBJECT_TYPE  = 6

  character, parameter :: backslash=achar(92)
  character, parameter :: space=achar(32)
  character, parameter :: bspace=achar(8)
  character, parameter :: tab=achar(9)
  character, parameter :: newline=achar(10)
  character, parameter :: vt=achar(11)
  character, parameter :: formfeed=achar(12)
  character, parameter :: creturn=achar(13)

  type, public :: json_null_t
    private
    integer :: type=JSON_UNDEF_TYPE
  end type json_null_t

  type, public :: json_logical_t
    private
    integer :: type=JSON_UNDEF_TYPE
    logical :: value=.false.
  end type json_logical_t

  type, public :: json_integer_t
    private
    integer :: type=JSON_UNDEF_TYPE
    integer :: value=0
  end type json_integer_t

  type, public :: json_real_t
    private
    integer       :: type=JSON_UNDEF_TYPE
    real(kind=wp) :: value=0.0_wp
  end type json_real_t

  type, public :: json_string_iterator_t
    private
    integer                          :: pos=0
    integer                          :: len=0
    character, dimension(:), pointer :: value
  end type json_string_iterator_t

  type, public :: json_string_t
    private
    integer                          :: type=JSON_UNDEF_TYPE
    integer                          :: len=0
    integer                          :: size=0
    character, dimension(:), pointer :: value
  end type json_string_t

  type, public :: json_value_t
    private
    integer                       :: type=JSON_UNDEF_TYPE
    type(json_null_t),    pointer :: null
    type(json_logical_t), pointer :: logical
    type(json_integer_t), pointer :: integer
    type(json_real_t),    pointer :: real
    type(json_string_t),  pointer :: string
    type(json_array_t),   pointer :: array
    type(json_object_t),  pointer :: object
  end type json_value_t

  type :: json_value_node_t
    private
    type(json_value_t),      pointer :: value
    type(json_value_node_t), pointer :: next
  end type json_value_node_t

  type, public :: json_array_iterator_t
    private
    type(json_value_node_t), pointer :: node
  end type json_array_iterator_t

  type, public :: json_array_t
    private
    integer                          :: type=JSON_UNDEF_TYPE
    integer                          :: size=0
    type(json_value_node_t), pointer :: head
    type(json_value_node_t), pointer :: tail
  end type json_array_t

  type, public :: json_member_t
    private
    integer                      :: type=JSON_UNDEF_TYPE
    type(json_string_t), pointer :: ident
    type(json_value_t),  pointer :: value
  end type json_member_t

  type :: json_member_node_t
    private
    type(json_member_t),      pointer :: member
    type(json_member_node_t), pointer :: next
  end type json_member_node_t

  type :: json_table_node_t
    private
    type(json_member_node_t), pointer :: head
  end type json_table_node_t

  type, public :: json_object_iterator_t
    private
    integer                                        :: pos=0
    integer                                        :: size=0
    type(json_member_node_t),              pointer :: node
    type(json_table_node_t), dimension(:), pointer :: table
  end type json_object_iterator_t

  type, public :: json_object_t
    private
    integer                                        :: type=JSON_UNDEF_TYPE
    integer                                        :: size=0
    integer                                        :: used=0
    type(json_table_node_t), dimension(:), pointer :: table
  end type json_object_t

  type, public :: json_t
    private
    integer                      :: type=JSON_UNDEF_TYPE
    type(json_array_t),  pointer :: array
    type(json_object_t), pointer :: object
  end type json_t

  interface operator(==)
    module procedure json_null_equal
    module procedure json_logical_equal
    module procedure json_integer_equal
    module procedure json_real_equal
    module procedure json_string_equal
    module procedure json_array_equal
    module procedure json_member_equal
    module procedure json_object_equal
    module procedure json_value_equal
    module procedure json_json_equal
  end interface operator(==)

  interface scan
    module procedure json_string_scan_char_string
    module procedure json_string_scan_string_char
    module procedure json_string_scan_string_string
  end interface scan

  interface json_isdef
    module procedure json_null_isdef
    module procedure json_logical_isdef
    module procedure json_integer_isdef
    module procedure json_real_isdef
    module procedure json_string_iterator_isdef
    module procedure json_string_isdef
    module procedure json_array_iterator_isdef
    module procedure json_array_isdef
    module procedure json_member_isdef
    module procedure json_object_iterator_isdef
    module procedure json_object_isdef
    module procedure json_value_isdef
    module procedure json_json_isdef
  end interface json_isdef

  interface json_init
    module procedure json_null_init
    module procedure json_logical_init
    module procedure json_integer_init_integer
    module procedure json_integer_init_string
    module procedure json_integer_init_json_string
    module procedure json_real_init_real
    module procedure json_real_init_string
    module procedure json_real_init_json_string
    module procedure json_string_iterator_init
    module procedure json_string_init
    module procedure json_array_iterator_init
    module procedure json_array_init
    module procedure json_array_init_logical
    module procedure json_array_init_integer
    module procedure json_array_init_real
    module procedure json_array_init_string
    module procedure json_member_init
    module procedure json_object_iterator_init
    module procedure json_object_init
    module procedure json_value_init_null
    module procedure json_value_init_logical
    module procedure json_value_init_integer
    module procedure json_value_init_real
    module procedure json_value_init_string
    module procedure json_value_init_array
    module procedure json_value_init_object
    module procedure json_json_array_init
    module procedure json_json_object_init
  end interface json_init

  interface json_end
    module procedure json_null_end
    module procedure json_logical_end
    module procedure json_integer_end
    module procedure json_real_end
    module procedure json_string_iterator_end
    module procedure json_string_end
    module procedure json_array_iterator_end
    module procedure json_array_end
    module procedure json_member_end
    module procedure json_object_iterator_end
    module procedure json_object_end
    module procedure json_value_end
    module procedure json_json_end
  end interface json_end

  interface json_len
    module procedure json_string_len
    module procedure json_array_len
    module procedure json_object_len
  end interface json_len

  interface json_string
    module procedure json_null_string
    module procedure json_logical_string
    module procedure json_integer_string
    module procedure json_real_string
    module procedure json_string_string
    module procedure json_array_string
    module procedure json_member_string
    module procedure json_object_string
    module procedure json_value_string
    module procedure json_json_string
  end interface json_string

  interface json_write
    module procedure json_null_write
    module procedure json_logical_write
    module procedure json_integer_write
    module procedure json_real_write
    module procedure json_string_write
    module procedure json_array_write
    module procedure json_member_write
    module procedure json_object_write
    module procedure json_value_write
    module procedure json_json_write
  end interface json_write

  interface json_set
    module procedure json_string_set
    module procedure json_array_set_value
    module procedure json_array_set_null
    module procedure json_array_set_logical
    module procedure json_array_set_integer
    module procedure json_array_set_real
    module procedure json_array_set_string
    module procedure json_array_set_array
    module procedure json_array_set_array_logical
    module procedure json_array_set_array_integer
    module procedure json_array_set_array_real
    module procedure json_array_set_array_string
    module procedure json_array_set_object
    module procedure json_object_set_member
    module procedure json_object_set_value
    module procedure json_object_set_null
    module procedure json_object_set_logical
    module procedure json_object_set_integer
    module procedure json_object_set_real
    module procedure json_object_set_string
    module procedure json_object_set_array
    module procedure json_object_set_array_logical
    module procedure json_object_set_array_integer
    module procedure json_object_set_array_real
    module procedure json_object_set_array_string
    module procedure json_object_set_object
  end interface json_set

  interface json_get
    module procedure json_null_get
    module procedure json_logical_get
    module procedure json_integer_get
    module procedure json_real_get
    module procedure json_string_get_string
    module procedure json_string_get_char
    module procedure json_array_get_value
    module procedure json_array_get_null
    module procedure json_array_get_logical
    module procedure json_array_get_integer
    module procedure json_array_get_real
    module procedure json_array_get_string
    module procedure json_array_get_array
    module procedure json_array_get_array_logical
    module procedure json_array_get_array_integer
    module procedure json_array_get_array_real
    module procedure json_array_get_array_string
    module procedure json_array_get_object
    module procedure json_object_get_value
    module procedure json_object_get_null
    module procedure json_object_get_logical
    module procedure json_object_get_integer
    module procedure json_object_get_real
    module procedure json_object_get_string
    module procedure json_object_get_array
    module procedure json_object_get_array_logical
    module procedure json_object_get_array_integer
    module procedure json_object_get_array_real
    module procedure json_object_get_array_string
    module procedure json_object_get_object
    module procedure json_json_get_array
    module procedure json_json_get_object
  end interface json_get

  interface json_append
    module procedure json_string_append
    module procedure json_array_append_value
    module procedure json_array_append_null
    module procedure json_array_append_logical
    module procedure json_array_append_integer
    module procedure json_array_append_real
    module procedure json_array_append_string
    module procedure json_array_append_array
    module procedure json_array_append_array_logical
    module procedure json_array_append_array_integer
    module procedure json_array_append_array_real
    module procedure json_array_append_array_string
    module procedure json_array_append_object
  end interface json_append

  interface json_extend
    module procedure json_string_extend_char
    module procedure json_string_extend_string
  end interface json_extend

  interface json_pop
    module procedure json_string_pop
    module procedure json_array_pop
    module procedure json_object_pop
  end interface json_pop

  interface json_next
    module procedure json_string_iterator_next
    module procedure json_array_iterator_next
    module procedure json_object_iterator_next
  end interface json_next

contains

  subroutine json_write_string(string, unit)
    character(len=*),  intent(in) :: string
    integer, optional, intent(in) :: unit
    !
    if(present(unit))then
      write(unit=unit, fmt="(a)", advance="no") string
    else
      write(unit=*, fmt="(a)", advance="no") string
    end if
    return
  end subroutine json_write_string

  subroutine json_write_line(unit)
    integer, optional, intent(in) :: unit
    !
    if(present(unit))then
      write(unit=unit, fmt=*)
    else
      write(unit=*, fmt=*)
    end if
    return
  end subroutine json_write_line

  elemental function json_null_isdef(this) result(is)
    type(json_null_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_NULL_TYPE)
    return
  end function json_null_isdef

  elemental subroutine json_null_init(this)
    type(json_null_t), intent(out) :: this
    !
    this%type=JSON_NULL_TYPE
    return
  end subroutine json_null_init

  elemental subroutine json_null_end(this)
    type(json_null_t), intent(inout) :: this
    !
    this%type=JSON_UNDEF_TYPE
    return
  end subroutine json_null_end

  elemental subroutine json_null_get(this, ierr)
    type(json_null_t), intent(in)  :: this
    integer,           intent(out) :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_null_isdef(this))ierr=JSON_OK
    return
  end subroutine json_null_get

  elemental function json_null_equal(this_1, this_2) result(eqv)
    type(json_null_t), intent(in) :: this_1
    type(json_null_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_null_isdef(this_1).and.json_null_isdef(this_2))eqv=.true.
    return
  end function json_null_equal

  subroutine json_null_string(this, string)
    type(json_null_t),   intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    if(json_null_isdef(this))call json_string_extend_char(string, "null")
    return
  end subroutine json_null_string

  subroutine json_null_write(this, unit)
    type(json_null_t), intent(in) :: this
    integer, optional, intent(in) :: unit
    !
    if(json_null_isdef(this)) call json_write_string("null", unit)
    return
  end subroutine json_null_write

  elemental function json_logical_isdef(this) result(is)
    type(json_logical_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_LOGICAL_TYPE)
    return
  end function json_logical_isdef

  elemental subroutine json_logical_init(this, value)
    type(json_logical_t), intent(out) :: this
    logical,              intent(in)  :: value
    !
    this%value=value
    this%type=JSON_LOGICAL_TYPE
    return
  end subroutine json_logical_init

  elemental subroutine json_logical_end(this)
    type(json_logical_t), intent(inout) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%value=.false.
    return
  end subroutine json_logical_end

  elemental function json_logical_equal(this_1, this_2) result(eqv)
    type(json_logical_t), intent(in) :: this_1
    type(json_logical_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_logical_isdef(this_1).and.json_logical_isdef(this_2))&
      eqv=(this_1%value.eqv.this_2%value)
    return
  end function json_logical_equal

  elemental subroutine json_logical_get(this, value, ierr)
    type(json_logical_t), intent(in)  :: this
    logical,              intent(out) :: value
    integer,              intent(out) :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_logical_isdef(this))then
      value=this%value
      ierr=JSON_OK
    end if
    return
  end subroutine json_logical_get

  subroutine json_logical_string(this, string)
    type(json_logical_t), intent(in)    :: this
    type(json_string_t),  intent(inout) :: string
    !
    if(this%value)then
      call json_string_extend_char(string, "true")
    else
      call json_string_extend_char(string, "false")
    end if
    return
  end subroutine json_logical_string

  subroutine json_logical_write(this, unit)
    type(json_logical_t), intent(in) :: this
    integer,    optional, intent(in) :: unit
    !
    if(json_logical_isdef(this))then
      if(this%value)then
        call json_write_string("true", unit)
      else
        call json_write_string("false", unit)
      end if
    end if
    return
  end subroutine json_logical_write

  elemental function json_integer_isdef(this) result(is)
    type(json_integer_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_INTEGER_TYPE)
    return
  end function json_integer_isdef

  elemental subroutine json_integer_init_integer(this, value)
    type(json_integer_t), intent(out) :: this
    integer,              intent(in)  :: value
    !
    this%value=value
    this%type=JSON_INTEGER_TYPE
    return
  end subroutine json_integer_init_integer

  elemental subroutine json_integer_init_string(this, value)
    type(json_integer_t), intent(out) :: this
    character(len=*),     intent(in)  :: value
    !
    read(unit=value, fmt=*) this%value
    this%type=JSON_INTEGER_TYPE
    return
  end subroutine json_integer_init_string

  pure subroutine json_integer_init_json_string(this, value)
    type(json_integer_t), intent(out) :: this
    type(json_string_t),  intent(in)  :: value
    !
    character(len=value%len) :: buff
    integer                  :: ierr
    !
    call json_string_get_string(value, buff, ierr)
    if(ierr==JSON_OK)call json_integer_init_string(this, trim(adjustl(buff)))
    return
  end subroutine json_integer_init_json_string

  elemental subroutine json_integer_end(this)
    type(json_integer_t), intent(inout) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%value=0
    return
  end subroutine json_integer_end

  elemental function json_integer_len(this) result(len)
    type(json_integer_t), intent(in) :: this
    !
    integer :: len
    !
    if(this%value>0)then
      len=floor(log10(real(abs(this%value),kind=wp)))+1
    elseif(this%value<0)then
      len=floor(log10(real(abs(this%value),kind=wp)))+2
    else
      len=1
    end if
    return
  end function json_integer_len

  elemental function json_integer_equal(this_1, this_2) result(eqv)
    type(json_integer_t), intent(in) :: this_1
    type(json_integer_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_integer_isdef(this_1).and.json_integer_isdef(this_2))&
      eqv=(this_1%value==this_2%value)
    return
  end function json_integer_equal

  elemental subroutine json_integer_get(this, value, ierr)
    type(json_integer_t), intent(in)  :: this
    integer,              intent(out) :: value
    integer,              intent(out) :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_integer_isdef(this))then
      value=this%value
      ierr=JSON_OK
    end if
    return
  end subroutine json_integer_get

  subroutine json_integer_string(this, string)
    type(json_integer_t), intent(in)    :: this
    type(json_string_t),  intent(inout) :: string
    !
    character(len=json_integer_len(this)) :: buff
    character(len=6)                      :: frmt
    !
    if(json_integer_isdef(this))then
      write(unit=frmt, fmt="(a2,i3.3,a1)") "(i", len(buff), ")"
      write(unit=buff, fmt=frmt) this%value
      call json_string_extend_char(string, trim(adjustl(buff)))
    end if
    return
  end subroutine json_integer_string

  subroutine json_integer_write(this, unit)
    type(json_integer_t), intent(in) :: this
    integer,    optional, intent(in) :: unit
    !
    character(len=json_integer_len(this)) :: buff
    character(len=6)                      :: frmt
    !
    if(json_integer_isdef(this))then
      write(unit=frmt, fmt="(a2,i3.3,a1)") "(i", len(buff), ")"
      write(unit=buff, fmt=frmt) this%value
      call json_write_string(trim(adjustl(buff)), unit)
    end if
    return
  end subroutine json_integer_write

  elemental function json_real_isdef(this) result(is)
    type(json_real_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_REAL_TYPE)
    return
  end function json_real_isdef

  elemental subroutine json_real_init_real(this, value)
    type(json_real_t), intent(out) :: this
    real(kind=wp),     intent(in)  :: value
    !
    this%value=value
    this%type=JSON_REAL_TYPE
    return
  end subroutine json_real_init_real

  elemental subroutine json_real_init_string(this, value)
    type(json_real_t), intent(out) :: this
    character(len=*),  intent(in)  :: value
    !
    read(unit=value, fmt=*) this%value
    this%type=JSON_REAL_TYPE
    return
  end subroutine json_real_init_string

  pure subroutine json_real_init_json_string(this, value)
    type(json_real_t),   intent(out) :: this
    type(json_string_t), intent(in)  :: value
    !
    character(len=value%len) :: buff
    integer                  :: ierr
    !
    call json_string_get_string(value, buff, ierr)
    if(ierr==JSON_OK)call json_real_init_string(this, trim(adjustl(buff)))
    return
  end subroutine json_real_init_json_string

  elemental subroutine json_real_end(this)
    type(json_real_t), intent(inout) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%value=0.0_wp
    return
  end subroutine json_real_end

  elemental function json_real_equal(this_1, this_2) result(eqv)
    type(json_real_t), intent(in) :: this_1
    type(json_real_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_real_isdef(this_1).and.json_real_isdef(this_2))&
      eqv=abs(this_1%value-this_2%value)<2.0_wp*spacing(max(abs(this_1%value),abs(this_2%value)))
    return
  end function json_real_equal

  elemental subroutine json_real_get(this, value, ierr)
    type(json_real_t), intent(in)  :: this
    real(kind=wp),     intent(out) :: value
    integer,           intent(out) :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_real_isdef(this))then
      value=this%value
      ierr=JSON_OK
    end if
    return
  end subroutine json_real_get

  subroutine json_real_string(this, string)
    type(json_real_t),   intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    character(len=precision(this%value)+9) :: buff
    character(len=11)                      :: fmt
    integer                                :: p
    !
    if(json_real_isdef(this))then
      p=precision(this%value)
      write(unit=fmt, fmt="(a3,i3.3,a1,i3.3,a1)") "(es", p+9, ".", p+1, ")"
      write(unit=buff, fmt=fmt) this%value
      call json_string_extend_char(string, trim(adjustl(buff)))
    end if
    return
  end subroutine json_real_string

  subroutine json_real_write(this, unit)
    type(json_real_t), intent(in) :: this
    integer, optional, intent(in) :: unit
    !
    character(len=precision(this%value)+9) :: buff
    character(len=11)                      :: fmt
    integer                                :: p
    !
    if(json_real_isdef(this))then
      p=precision(this%value)
      write(unit=fmt, fmt="(a3,i3.3,a1,i3.3,a1)") "(es", p+9, ".", p+1, ")"
      write(unit=buff, fmt=fmt) this%value
      call json_write_string(trim(adjustl(buff)), unit)
    end if
    return
  end subroutine json_real_write

  elemental subroutine json_string_iterator_nullify(this)
    type(json_string_iterator_t), intent(out) :: this
    !
    this%pos=0
    this%len=0
    this%value=>null()
    return
  end subroutine json_string_iterator_nullify

  elemental function json_string_iterator_isdef(this) result(is)
    type(json_string_iterator_t), intent(in) :: this
    !
    logical :: is
    !
    is=associated(this%value)
    return
  end function json_string_iterator_isdef

  subroutine json_string_iterator_init(this, string)
    type(json_string_iterator_t), intent(out) :: this
    type(json_string_t),          intent(in)  :: string
    !
    call json_string_iterator_nullify(this)
    this%value=>null()
    if(json_string_isdef(string))then
      this%pos=1
      this%len=string%len
      this%value=>string%value
    end if
    return
  end subroutine json_string_iterator_init

  elemental subroutine json_string_iterator_end(this)
    type(json_string_iterator_t), intent(inout) :: this
    !
    this%pos=0
    this%len=0
    this%value=>null()
    return
  end subroutine json_string_iterator_end

  elemental subroutine json_string_iterator_next(this, char)
    type(json_string_iterator_t), intent(inout) :: this
    character,                    intent(out)   :: char
    !
    if(this%pos<=this%len)then
      char=this%value(this%pos)
      this%pos=this%pos+1
    end if
    return
  end subroutine json_string_iterator_next

  elemental subroutine json_string_nullify(this)
    type(json_string_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%len=0
    this%size=0
    this%value=>null()
    return
  end subroutine json_string_nullify

  elemental function json_string_isdef(this) result(is)
    type(json_string_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_STRING_TYPE)
    return
  end function json_string_isdef

  subroutine json_string_init(this, value)
    type(json_string_t),        intent(out) :: this
    character(len=*), optional, intent(in)  :: value
    !
    call json_string_nullify(this)
    this%len=0
    this%size=JSON_STRING_INIT_LEN
    SAFE_ALLOCATE(this%value(this%size))
    this%type=JSON_STRING_TYPE
    if(present(value))&
      call json_string_extend_char(this, value)
    return
  end subroutine json_string_init

  subroutine json_string_end(this)
    type(json_string_t), intent(inout) :: this
    !
    this%type=JSON_UNDEF_TYPE
    SAFE_DEALLOCATE_P(this%value)
    this%value=>null()
    this%size=0
    this%len=0
    return
  end subroutine json_string_end

  subroutine json_string_reallocate(this, extra)
    type(json_string_t), intent(inout) :: this
    integer,             intent(in)    :: extra
    !
    character, dimension(:), pointer :: buff
    real(kind=wp)                    :: need
    integer                          :: i, n
    !
    if(json_string_isdef(this))then
      need=real(this%len+extra, kind=wp)
      if(this%size<int(JSON_STRING_GROWTH_FACTOR*need))then
        n=max(ceiling((log(need)-log(real(this%size,kind=wp)))/log(JSON_STRING_GROWTH_FACTOR)),1)
        this%size=ceiling((JSON_STRING_GROWTH_FACTOR**n)*real(this%size,kind=wp))
        SAFE_ALLOCATE(buff(this%size))
        forall(i=1:this%len)buff(i)=this%value(i)
        SAFE_DEALLOCATE_P(this%value)
        this%value=>buff
      end if
    end if
    return
  end subroutine json_string_reallocate

  elemental function json_string_len(this) result(len)
    type(json_string_t), intent(in) :: this
    !
    integer :: len
    !
    len=this%len
    return
  end function json_string_len

  elemental function json_string_equal(this_1, this_2) result(eqv)
    type(json_string_t), intent(in) :: this_1
    type(json_string_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    integer :: i
    !
    eqv=.false.
    if(json_string_isdef(this_1).and.json_string_isdef(this_2))then
      if(this_1%len==this_2%len)then
        eqv=.true.
        do i = 1, this_1%len
          if(this_1%value(i)/=this_2%value(i))then
            eqv=.false.
            exit
          end if
        end do
      end if
    end if
    return
  end function json_string_equal

  pure function json_string_scan_char_string(string, set, back) result(i)
    character(len=*),    intent(in) :: string
    type(json_string_t), intent(in) :: set
    logical,   optional, intent(in) :: back
    !
    integer :: i
    !
    character(len=set%len) :: buff
    integer                :: ierr
    !
    i=0
    if(json_string_isdef(set))then
      call json_string_get_string(set, buff, ierr)
      if(ierr==JSON_OK) i=scan(string, buff, back)
    end if
    return
  end function json_string_scan_char_string

  pure function json_string_scan_string_char(string, set, back) result(i)
    type(json_string_t), intent(in) :: string
    character(len=*),    intent(in) :: set
    logical,   optional, intent(in) :: back
    !
    integer :: i
    !
    character(len=string%len) :: buff
    integer                   :: ierr
    !
    i=0
    if(json_string_isdef(string))then
      call json_string_get_string(string, buff, ierr)
      if(ierr==JSON_OK) i=scan(buff, set, back)
    end if
    return
  end function json_string_scan_string_char

  pure function json_string_scan_string_string(string, set, back) result(i)
    type(json_string_t), intent(in) :: string
    type(json_string_t), intent(in) :: set
    logical,   optional, intent(in) :: back
    !
    integer :: i
    !
    character(len=string%len) :: strb
    character(len=set%len)    :: setb
    integer                   :: ierr
    !
    i=0
    if(json_string_isdef(set).and.json_string_isdef(string))then
      call json_string_get_string(string, strb, ierr)
      if(ierr==JSON_OK)then
        call json_string_get_string(set, setb, ierr)
        if(ierr==JSON_OK) i=scan(strb, setb, back)
      end if
    end if
    return
  end function json_string_scan_string_string

  elemental subroutine json_string_escape(ichar, ochar, escape)
    character, intent(in)  :: ichar
    character, intent(out) :: ochar
    logical,   intent(out) :: escape
    !
    escape=.true.
    select case(ichar)
    case('"') 
      ochar='"'
    case(backslash)
      ochar=backslash
    case("/")
      ochar="/"
    case(bspace)
      ochar="b"
    case(formfeed)
      ochar="f"
    case(newline)
      ochar="n"
    case(creturn)
      ochar="r"
    case(tab)
      ochar="t"
    !case("u") !no unicode...
    case default
      escape=.false.
      ochar=ichar
    end select
    return
  end subroutine json_string_escape

  subroutine json_string_string(this, string)
    type(json_string_t), intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    character :: char
    integer   :: i
    logical   :: escape
    !
    if(json_string_isdef(this))then
      call json_string_append(string, '"')
      do i = 1, this%len
        call json_string_escape(this%value(i), char, escape)
        if(escape)call json_string_append(string, backslash)
        call json_string_append(string, char)
      end do
      call json_string_extend_string(string, this)
      call json_string_append(string, '"')
    end if
    return
  end subroutine json_string_string

  subroutine json_string_write(this, unit)
    type(json_string_t), intent(in) :: this
    integer,   optional, intent(in) :: unit
    !
    character :: char     
    integer   :: i
    logical   :: escape
    !
    if(json_string_isdef(this))then
      call json_write_string('"', unit)
      do i = 1, this%len
        call json_string_escape(this%value(i), char, escape)
        if(escape)call json_write_string(backslash, unit)
        call json_write_string(char, unit)
      end do
      call json_write_string('"', unit)
    end if
    return
  end subroutine json_string_write

  elemental subroutine json_string_set(this, i, char, ierr)
    type(json_string_t), intent(inout) :: this
    integer,             intent(in)    :: i
    character,           intent(in)    :: char
    integer,             intent(out)   :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_string_isdef(this))then
      ierr=JSON_INDEX_ERROR
      if((0<i).and.(i<=this%len))then
        this%value(i)=char
        ierr=JSON_OK
      end if
    end if
    return
  end subroutine json_string_set

  elemental subroutine json_string_get_string(this, string, ierr)
    type(json_string_t), intent(in)  :: this
    character(len=*),    intent(out) :: string
    integer,             intent(out) :: ierr
    !
    integer :: i
    !
    string=""
    ierr=JSON_UNDEF_ERROR
    if(json_string_isdef(this))then
      ierr=JSON_SIZE_ERROR
      if(this%len<=len(string))then
        do i = 1, this%len
          string(i:i)=this%value(i)
        end do
        ierr=JSON_OK
      end if
    end if
    return
  end subroutine json_string_get_string

  elemental subroutine json_string_get_char(this, i, char, ierr)
    type(json_string_t), intent(in)  :: this
    integer,             intent(in)  :: i
    character,           intent(out) :: char
    integer,             intent(out) :: ierr
    !
    ierr=JSON_UNDEF_ERROR
    if(json_string_isdef(this))then
      ierr=JSON_INDEX_ERROR
      if((0<i).and.(i<=this%len))then
        char=this%value(i)
        ierr=JSON_OK
      end if
    end if
    return
  end subroutine json_string_get_char

  elemental subroutine json_string_pop(this, char)
    type(json_string_t), intent(inout) :: this
    character,           intent(out)   :: char
    !
    char=this%value(this%len)
    this%value(this%len)=""
    this%len=this%len-1
    return
  end subroutine json_string_pop

  subroutine json_string_append(this, char)
    type(json_string_t), intent(inout) :: this
    character,           intent(in)    :: char
    !
    call json_string_reallocate(this, 1)
    this%len=this%len+1
    this%value(this%len)=char
    return
  end subroutine json_string_append

  subroutine json_string_extend_char(this, char)
    type(json_string_t), intent(inout) :: this
    character(len=*),    intent(in)    :: char
    !
    integer :: i
    !
    call json_string_reallocate(this, len(char))
    do i = 1, len(char)
      this%value(this%len+i)=char(i:i)
    end do
    this%len=this%len+len(char)
    return
  end subroutine json_string_extend_char

  subroutine json_string_extend_string(this, string)
    type(json_string_t), intent(inout) :: this
    type(json_string_t), intent(in)    :: string
    !
    call json_string_reallocate(this, string%len)
    this%value(this%len+1:this%len+string%len)=string%value(1:string%len)
    this%len=this%len+string%len
    return
  end subroutine json_string_extend_string

  elemental subroutine json_array_iterator_nullify(this)
    type(json_array_iterator_t), intent(out) :: this
    !
    this%node=>null()
    return
  end subroutine json_array_iterator_nullify

  elemental function json_array_iterator_isdef(this) result(is)
    type(json_array_iterator_t), intent(in) :: this
    !
    logical :: is
    !
    is=associated(this%node)
    return
  end function json_array_iterator_isdef

  subroutine json_array_iterator_init(this, array)
    type(json_array_iterator_t), intent(out) :: this
    type(json_array_t),          intent(in)  :: array
    !
    call json_array_iterator_nullify(this)
    this%node=>null()
    if(json_array_isdef(array).and.(array%size>0)) this%node=>array%head
    return
  end subroutine json_array_iterator_init

  elemental subroutine json_array_iterator_end(this)
    type(json_array_iterator_t), intent(inout) :: this
    !
    this%node=>null()
    return
  end subroutine json_array_iterator_end

  subroutine json_array_iterator_next(this, value)
    type(json_array_iterator_t), intent(inout) :: this
    type(json_value_t),          pointer       :: value
    !
    value=>null()
    if(associated(this%node))then
      value=>this%node%value
      this%node=>this%node%next
    end if
    return
  end subroutine json_array_iterator_next

  elemental subroutine json_array_nullify(this)
    type(json_array_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%size=0
    this%head=>null()
    this%tail=>null()
    return
  end subroutine json_array_nullify

  elemental function json_array_isdef(this) result(is)
    type(json_array_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_ARRAY_TYPE)
    return
  end function json_array_isdef

  elemental subroutine json_array_init(this)
    type(json_array_t), intent(out) :: this
    !
    call json_array_nullify(this)
    this%size=0
    this%head=>null()
    this%type=JSON_ARRAY_TYPE
    return
  end subroutine json_array_init

  subroutine json_array_init_logical(this, vals)
    type(json_array_t),    intent(out) :: this
    logical, dimension(:), intent(in)  :: vals
    !
    integer :: i
    !
    call json_array_init(this)
    do i = 1, size(vals)
      call json_array_append_logical(this, vals(i))
    end do
    return
  end subroutine json_array_init_logical

  subroutine json_array_init_integer(this, vals)
    type(json_array_t),    intent(out) :: this
    integer, dimension(:), intent(in)  :: vals
    !
    integer :: i
    !
    call json_array_init(this)
    do i = 1, size(vals)
      call json_array_append_integer(this, vals(i))
    end do
    return
  end subroutine json_array_init_integer

  subroutine json_array_init_real(this, vals)
    type(json_array_t),          intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: vals
    !
    integer :: i
    !
    call json_array_init(this)
    do i = 1, size(vals)
      call json_array_append_real(this, vals(i))
    end do
    return
  end subroutine json_array_init_real

  subroutine json_array_init_string(this, vals)
    type(json_array_t),             intent(out) :: this
    character(len=*), dimension(:), intent(in)  :: vals
    !
    integer :: i
    !
    call json_array_init(this)
    do i = 1, size(vals)
      call json_array_append_string(this, vals(i))
    end do
    return
  end subroutine json_array_init_string

  subroutine json_array_end(this)
    type(json_array_t), intent(inout) :: this
    !
    type(json_value_node_t), pointer :: node
    !
    this%type=JSON_UNDEF_TYPE
    do
      node=>this%head
      if(.not.associated(node))exit
      this%size=this%size-1
      call json_value_end(node%value)
      SAFE_DEALLOCATE_P(node%value)
      node%value=>null()
      this%head=>node%next
      SAFE_DEALLOCATE_P(node)
    end do
    this%size=0
    this%head=>null()
    return
  end subroutine json_array_end

  elemental function json_array_len(this) result(size)
    type(json_array_t), intent(in) :: this
    !
    integer :: size
    !
    size=this%size
    return
  end function json_array_len

  function json_array_equal(this_1, this_2) result(eqv)
    type(json_array_t), intent(in) :: this_1
    type(json_array_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    type(json_array_iterator_t) :: iter_1, iter_2
    type(json_value_t), pointer :: value_1, value_2
    !
    eqv=.false.
    if(json_array_isdef(this_1).and.json_array_isdef(this_2))then
      if(this_1%size==this_2%size)then
        eqv=.true.
        call json_array_iterator_init(iter_1, this_1)
        call json_array_iterator_init(iter_2, this_2)
        call json_array_iterator_next(iter_1, value_1)
        call json_array_iterator_next(iter_2, value_2)
        do while(associated(value_1).and.associated(value_2))
          eqv=eqv.and.json_value_equal(value_1, value_2)
          if(.not.eqv)exit
          call json_array_iterator_next(iter_1, value_1)
          call json_array_iterator_next(iter_2, value_2)
        end do
      end if
      call json_array_iterator_end(iter_1)
      call json_array_iterator_end(iter_2)
    end if
    return
  end function json_array_equal

  subroutine json_array_string(this, string)
    type(json_array_t),  intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    type(json_array_iterator_t) :: iter
    type(json_value_t), pointer :: value
    !
    if(json_array_isdef(this))then
      call json_string_append(string, "[")
      call json_array_iterator_init(iter, this)
      call json_array_iterator_next(iter, value)
      if(associated(value))then
        do
          call json_value_string(value, string)
          call json_array_iterator_next(iter, value)
          if(.not.associated(value))exit
          call json_string_append(string, ",")
        end do
      end if
      call json_array_iterator_end(iter)
      call json_string_append(string, "]")
    end if
    return
  end subroutine json_array_string

  subroutine json_array_write(this, unit, level, separator, count)
    type(json_array_t), intent(in) :: this
    integer,   optional, intent(in) :: unit
    integer,   optional, intent(in) :: level
    character, optional, intent(in) :: separator
    integer,   optional, intent(in) :: count
    !
    type(json_array_iterator_t) :: iter
    type(json_value_t), pointer :: value
    character                   :: sep
    integer                     :: lvl, cnt
    !
    if(json_array_isdef(this))then
      lvl=0
      cnt=1
      sep=" "
      if(present(level))lvl=level
      if(present(count))cnt=count
      if(present(separator))sep=separator
      call json_write_string("[", unit)
      call json_array_iterator_init(iter, this)
      call json_array_iterator_next(iter, value)
      if(associated(value))then
        do
          call json_write_line(unit)
          call json_write_string(repeat(sep, cnt*(lvl+1)), unit)
          call json_value_write(value, unit, lvl+1, separator, count)
          call json_array_iterator_next(iter, value)
          if(.not.associated(value))exit
          call json_write_string(",", unit)
        end do
        call json_write_line(unit)
        call json_write_string(repeat(sep, cnt*lvl), unit)
      end if
      call json_write_string("]", unit)
      call json_array_iterator_end(iter)
    end if
    return
  end subroutine json_array_write

  subroutine json_array_append_value(this, value)
    type(json_array_t),         intent(inout) :: this
    type(json_value_t), target, intent(in)    :: value
    !
    type(json_value_node_t), pointer :: node
    !
    if(json_value_isdef(value).and.json_array_isdef(this))then
      SAFE_ALLOCATE(node)
      node%value=>value
      node%next=>null()
      if(associated(this%head))then
        this%tail%next=>node
      else
        this%head=>node
      end if
      this%tail=>node
      this%size=this%size+1
    end if
    return
  end subroutine json_array_append_value

  subroutine json_array_append_null(this)
    type(json_array_t), intent(inout) :: this
    !
    type(json_value_t), pointer :: json_value
    type(json_null_t),  pointer :: type_value
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_null_init(type_value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_null(json_value, type_value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_null

  subroutine json_array_append_logical(this, value)
    type(json_array_t), intent(inout) :: this
    logical,            intent(in)    :: value
    !
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: type_value
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_logical_init(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_logical(json_value, type_value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_logical

  subroutine json_array_append_integer(this, value)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: value
    !
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: type_value
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_integer_init_integer(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_integer(json_value, type_value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_integer

  subroutine json_array_append_real(this, value)
    type(json_array_t), intent(inout) :: this
    real(kind=wp),      intent(in)    :: value
    !
    type(json_value_t), pointer :: json_value
    type(json_real_t),  pointer :: type_value
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_real_init_real(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_real(json_value, type_value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_real

  subroutine json_array_append_string(this, value)
    type(json_array_t), intent(inout) :: this
    character(len=*),   intent(in)    :: value
    !
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: type_value
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_string_init(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_string(json_value, type_value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_string

  subroutine json_array_append_array(this, value)
    type(json_array_t), intent(inout) :: this
    type(json_array_t), intent(in)    :: value
    !
    type(json_value_t), pointer :: json_value
    !
    if(json_array_isdef(value).and.json_array_isdef(this))then
      SAFE_ALLOCATE(json_value)
      call json_value_init_array(json_value, value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_array

  subroutine json_array_append_array_logical(this, vals)
    type(json_array_t),    intent(inout) :: this
    logical, dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_logical(json_array, vals)
      call json_array_append_array(this, json_array)
    end if
    return
  end subroutine json_array_append_array_logical

  subroutine json_array_append_array_integer(this, vals)
    type(json_array_t),    intent(inout) :: this
    integer, dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_integer(json_array, vals)
      call json_array_append_array(this, json_array)
    end if
    return
  end subroutine json_array_append_array_integer

  subroutine json_array_append_array_real(this, vals)
    type(json_array_t),          intent(inout) :: this
    real(kind=wp), dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_real(json_array, vals)
      call json_array_append_array(this, json_array)
    end if
    return
  end subroutine json_array_append_array_real

  subroutine json_array_append_array_string(this, vals)
    type(json_array_t),             intent(inout) :: this
    character(len=*), dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_string(json_array, vals)
      call json_array_append_array(this, json_array)
    end if
    return
  end subroutine json_array_append_array_string

  subroutine json_array_append_object(this, value)
    type(json_array_t),  intent(inout) :: this
    type(json_object_t), intent(in)    :: value
    !
    type(json_value_t), pointer :: json_value
    !
    if(json_object_isdef(value).and.json_array_isdef(this))then
      SAFE_ALLOCATE(json_value)
      call json_value_init_object(json_value, value)
      call json_array_append_value(this, json_value)
    end if
    return
  end subroutine json_array_append_object

  subroutine json_array_set_value(this, i, value, ierr)
    type(json_array_t),         intent(inout) :: this
    integer,                    intent(in)    :: i
    type(json_value_t), target, intent(in)    :: value
    integer,                    intent(out)   :: ierr
    !
    type(json_value_node_t), pointer :: node
    integer                          :: idx
    !
    ierr=JSON_UNDEF_ERROR
    if(json_value_isdef(value).and.json_array_isdef(this))then
      ierr=JSON_INDEX_ERROR
      if((0<i).and.(i<=this%size))then
        idx=1
        node=>this%head
        do while(associated(node))
          if(idx==i)exit
          idx=idx+1
          node=>node%next
        end do
        if(associated(node))then
          call json_value_end(node%value)
          SAFE_DEALLOCATE_P(node%value)
          node%value=>value
          ierr=JSON_OK
        end if
      end if
    end if
    return
  end subroutine json_array_set_value

  subroutine json_array_set_null(this, i, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    integer,            intent(out)   :: ierr
    !
    type(json_value_t), pointer :: json_value
    type(json_null_t),  pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_null_init(type_value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_null(json_value, type_value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_null

  subroutine json_array_set_logical(this, i, value, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    logical,            intent(in)    :: value
    integer,            intent(out)   :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_logical_init(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_logical(json_value, type_value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_logical

  subroutine json_array_set_integer(this, i, value, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    integer,            intent(in)    :: value
    integer,            intent(out)   :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_integer_init_integer(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_integer(json_value, type_value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_integer

  subroutine json_array_set_real(this, i, value, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    real(kind=wp),      intent(in)    :: value
    integer,            intent(out)   :: ierr
    !
    type(json_value_t), pointer :: json_value
    type(json_real_t),  pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_real_init_real(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_real(json_value, type_value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_real

  subroutine json_array_set_string(this, i, value, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    character(len=*),   intent(in)    :: value
    integer,            intent(out)   :: ierr
    !
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(type_value)
      call json_string_init(type_value, value)
      SAFE_ALLOCATE(json_value)
      call json_value_init_string(json_value, type_value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_string

  subroutine json_array_set_array(this, i, value, ierr)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: i
    type(json_array_t), intent(in)    :: value
    integer,            intent(out)   :: ierr
    !
    type(json_value_t), pointer :: json_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(value).and.json_array_isdef(this))then
      SAFE_ALLOCATE(json_value)
      call json_value_init_array(json_value, value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_array

  subroutine json_array_set_array_logical(this, i, vals, ierr)
    type(json_array_t),    intent(inout) :: this
    integer,               intent(in)    :: i
    logical, dimension(:), intent(in)    :: vals
    integer,               intent(out)   :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_logical(json_array, vals)
      call json_array_set_array(this, i, json_array, ierr)
    end if
    return
  end subroutine json_array_set_array_logical

  subroutine json_array_set_array_integer(this, i, vals, ierr)
    type(json_array_t),    intent(inout) :: this
    integer,               intent(in)    :: i
    integer, dimension(:), intent(in)    :: vals
    integer,               intent(out)   :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_integer(json_array, vals)
      call json_array_set_array(this, i, json_array, ierr)
    end if
    return
  end subroutine json_array_set_array_integer

  subroutine json_array_set_array_real(this, i, vals, ierr)
    type(json_array_t),          intent(inout) :: this
    integer,                     intent(in)    :: i
    real(kind=wp), dimension(:), intent(in)    :: vals
    integer,                     intent(out)   :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_real(json_array, vals)
      call json_array_set_array(this, i, json_array, ierr)
    end if
    return
  end subroutine json_array_set_array_real

  subroutine json_array_set_array_string(this, i, vals, ierr)
    type(json_array_t),             intent(inout) :: this
    integer,                        intent(in)    :: i
    character(len=*), dimension(:), intent(in)    :: vals
    integer,                        intent(out)   :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      SAFE_ALLOCATE(json_array)
      call json_array_init_string(json_array, vals)
      call json_array_set_array(this, i, json_array, ierr)
    end if
    return
  end subroutine json_array_set_array_string

  subroutine json_array_set_object(this, i, value, ierr)
    type(json_array_t),  intent(inout) :: this
    integer,             intent(in)    :: i
    type(json_object_t), intent(in)    :: value
    integer,   optional, intent(out)   :: ierr
    !
    type(json_value_t), pointer :: json_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(value).and.json_array_isdef(this))then
      SAFE_ALLOCATE(json_value)
      call json_value_init_object(json_value, value)
      call json_array_set_value(this, i, json_value, ierr)
    end if
    return
  end subroutine json_array_set_object

  subroutine json_array_get_self_logical(this, vals, ierr)
    type(json_array_t),    intent(in)  :: this
    logical, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: type_value
    type(json_array_iterator_t)   :: iter
    integer                       :: i
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      ierr=JSON_SIZE_ERROR
      if(this%size==size(vals))then
        call json_array_iterator_init(iter, this)
        if(json_array_iterator_isdef(iter))then
          ierr=JSON_OK
          do i = 1, this%size
            json_value=>null()
            call json_array_iterator_next(iter, json_value)
            if(.not.associated(json_value))then
              ierr=JSON_SIZE_ERROR
              exit
            end if
            call json_value_get_logical(json_value, type_value, ierr)
            if(ierr/=JSON_OK)exit
            call json_logical_get(type_value, vals(i), ierr)
            if(ierr/=JSON_OK)exit
          end do
        end if
        call json_array_iterator_end(iter)
      end if
    end if
    return
  end subroutine json_array_get_self_logical

  subroutine json_array_get_self_integer(this, vals, ierr)
    type(json_array_t),    intent(in)  :: this
    integer, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: type_value
    type(json_array_iterator_t)   :: iter
    integer                       :: i
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      ierr=JSON_SIZE_ERROR
      if(this%size==size(vals))then
        call json_array_iterator_init(iter, this)
        if(json_array_iterator_isdef(iter))then
          ierr=JSON_OK
          do i = 1, this%size
            json_value=>null()
            call json_array_iterator_next(iter, json_value)
            if(.not.associated(json_value))then
              ierr=JSON_SIZE_ERROR
              exit
            end if
            call json_value_get_integer(json_value, type_value, ierr)
            if(ierr/=JSON_OK)exit
            call json_integer_get(type_value, vals(i), ierr)
            if(ierr/=JSON_OK)exit
          end do
        end if
        call json_array_iterator_end(iter)
      end if
    end if
    return
  end subroutine json_array_get_self_integer

  subroutine json_array_get_self_real(this, vals, ierr)
    type(json_array_t),          intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: vals
    integer,                     intent(out) :: ierr
    !
    type(json_value_t), pointer :: json_value
    type(json_real_t),  pointer :: type_value
    type(json_array_iterator_t) :: iter
    integer                     :: i
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      ierr=JSON_SIZE_ERROR
      if(this%size==size(vals))then
        call json_array_iterator_init(iter, this)
        if(json_array_iterator_isdef(iter))then
          ierr=JSON_OK
          do i = 1, this%size
            json_value=>null()
            call json_array_iterator_next(iter, json_value)
            if(.not.associated(json_value))then
              ierr=JSON_SIZE_ERROR
              exit
            end if
            call json_value_get_real(json_value, type_value, ierr)
            if(ierr/=JSON_OK)exit
            call json_real_get(type_value, vals(i), ierr)
            if(ierr/=JSON_OK)exit
          end do
        end if
        call json_array_iterator_end(iter)
      end if
    end if
    return
  end subroutine json_array_get_self_real

  subroutine json_array_get_self_string(this, vals, ierr)
    type(json_array_t),             intent(in)  :: this
    character(len=*), dimension(:), intent(out) :: vals
    integer,                        intent(out) :: ierr
    !
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: type_value
    type(json_array_iterator_t)  :: iter
    integer                      :: i
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      ierr=JSON_SIZE_ERROR
      if(this%size==size(vals))then
        call json_array_iterator_init(iter, this)
        if(json_array_iterator_isdef(iter))then
          ierr=JSON_OK
          do i = 1, this%size
            json_value=>null()
            call json_array_iterator_next(iter, json_value)
            if(.not.associated(json_value))then
              ierr=JSON_SIZE_ERROR
              exit
            end if
            call json_value_get_string(json_value, type_value, ierr)
            if(ierr/=JSON_OK)exit
            call json_string_get_string(type_value, vals(i), ierr)
            if(ierr/=JSON_OK)exit
          end do
        end if
        call json_array_iterator_end(iter)
      end if
    end if
    return
  end subroutine json_array_get_self_string

  subroutine json_array_get_value(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    type(json_value_t), pointer     :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_node_t), pointer :: node
    integer                          :: idx
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      ierr=JSON_INDEX_ERROR
      if((0<i).and.(i<=this%size))then
        idx=1
        node=>this%head
        do while(associated(node))
          if(idx==i)exit
          idx=idx+1
          node=>node%next
        end do
        value=>node%value
      end if
      if(associated(value))ierr=JSON_OK
    end if
    return
  end subroutine json_array_get_value

  subroutine json_array_get_null(this, i, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    integer,            intent(out) :: ierr
    !
    type(json_value_t), pointer :: json_value
    type(json_null_t),  pointer :: type_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        type_value=>null()
        call json_value_get_null(json_value, type_value, ierr)
        if(ierr==JSON_OK) call json_null_get(type_value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_null

  subroutine json_array_get_logical(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    logical,            intent(out) :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_logical(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_logical_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_logical

  subroutine json_array_get_integer(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    integer,            intent(out) :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_integer(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_integer_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_integer

  subroutine json_array_get_real(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    real(kind=wp),      intent(out) :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_t), pointer :: json_value
    type(json_real_t),  pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_real(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_real_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_real

  subroutine json_array_get_string(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    character(len=*),   intent(out) :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_string(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_string_get_string(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_string

  subroutine json_array_get_array(this, i, value, ierr)
    type(json_array_t), intent(in)  :: this
    integer,            intent(in)  :: i
    type(json_array_t), pointer     :: value
    integer,            intent(out) :: ierr
    !
    type(json_value_t),  pointer :: json_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        value=>null()
        call json_value_get_array(json_value, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_array

  subroutine json_array_get_array_logical(this, i, vals, ierr)
    type(json_array_t),    intent(in)  :: this
    integer,               intent(in)  :: i
    logical, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_array=>null()
      call json_array_get_array(this, i, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_logical(json_array, vals, ierr)
    end if
    return
  end subroutine json_array_get_array_logical

  subroutine json_array_get_array_integer(this, i, vals, ierr)
    type(json_array_t),    intent(in)  :: this
    integer,               intent(in)  :: i
    integer, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_array=>null()
      call json_array_get_array(this, i, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_integer(json_array, vals, ierr)
    end if
    return
  end subroutine json_array_get_array_integer

  subroutine json_array_get_array_real(this, i, vals, ierr)
    type(json_array_t),          intent(in)  :: this
    integer,                     intent(in)  :: i
    real(kind=wp), dimension(:), intent(out) :: vals
    integer,                     intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_array=>null()
      call json_array_get_array(this, i, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_real(json_array, vals, ierr)
    end if
    return
  end subroutine json_array_get_array_real

  subroutine json_array_get_array_string(this, i, vals, ierr)
    type(json_array_t),             intent(in)  :: this
    integer,                        intent(in)  :: i
    character(len=*), dimension(:), intent(out) :: vals
    integer,                        intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_array=>null()
      call json_array_get_array(this, i, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_string(json_array, vals, ierr)
    end if
    return
  end subroutine json_array_get_array_string

  subroutine json_array_get_object(this, i, value, ierr)
    type(json_array_t),  intent(in)  :: this
    integer,             intent(in)  :: i
    type(json_object_t), pointer     :: value
    integer,             intent(out) :: ierr
    !
    type(json_value_t), pointer :: json_value
    !
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this))then
      json_value=>null()
      call json_array_get_value(this, i, json_value, ierr)
      if(ierr==JSON_OK)then
        value=>null()
        call json_value_get_object(json_value, value, ierr)
      end if
    end if
    return
  end subroutine json_array_get_object

  subroutine json_array_pop(this, value)
    type(json_array_t),     intent(inout) :: this
    type(json_value_t),     pointer       :: value
    !
    type(json_value_node_t), pointer :: node
    !
    if(associated(this%head))then
      node=>this%head
      this%head=>node%next
      this%size=this%size-1
      value=>node%value
      if(.not.associated(this%head))this%tail=>null()
      SAFE_DEALLOCATE_P(node)
    else
      value=>null()
    end if
    return
  end subroutine json_array_pop

  elemental subroutine json_member_nullify(this)
    type(json_member_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%ident=>null()
    this%value=>null()
    return
  end subroutine json_member_nullify

  elemental function json_member_isdef(this) result(is)
    type(json_member_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type/=JSON_UNDEF_TYPE)
    return
  end function json_member_isdef

  subroutine json_member_init(this, ident, value)
    type(json_member_t),         intent(out) :: this
    type(json_string_t), target, intent(in)  :: ident
    type(json_value_t),  target, intent(in)  :: value
    !
    call json_member_nullify(this)
    this%ident=>ident
    this%value=>value
    this%type=json_value_type(this%value)
    return
  end subroutine json_member_init

  subroutine json_member_end(this)
    type(json_member_t), intent(inout) :: this
    !
    if(associated(this%ident))then
      call json_string_end(this%ident)
      SAFE_DEALLOCATE_P(this%ident)
    end if
    this%ident=>null()
    if(associated(this%value))then
      call json_value_end(this%value)
      SAFE_DEALLOCATE_P(this%value)
    end if
    this%value=>null()
    this%type=JSON_UNDEF_TYPE
    return
  end subroutine json_member_end

  function json_member_equal(this_1, this_2) result(eqv)
    type(json_member_t), intent(in) :: this_1
    type(json_member_t), intent(in) :: this_2
    !
    logical :: eqv
    integer :: i
    !
    eqv=.false.
    if(json_member_isdef(this_1).and.json_member_isdef(this_2))then
      eqv=json_string_equal(this_1%ident, this_2%ident)
      do i = 1, this_1%ident%len
        write(unit=*,fmt="(a1)",advance="no") this_1%ident%value(i)
      end do
      write(unit=*,fmt=*)
      do i = 1, this_2%ident%len
        write(unit=*,fmt="(a1)",advance="no") this_2%ident%value(i)
      end do
      write(unit=*,fmt=*)
      eqv=eqv.and.json_value_equal(this_1%value, this_2%value)
    end if
    return
  end function json_member_equal

  subroutine json_member_string(this, string)
    type(json_member_t), intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    if(json_member_isdef(this))then
      call json_string_string(this%ident, string)
      call json_string_append(string, ":")
      call json_value_string(this%value, string)
    end if
    return
  end subroutine json_member_string

  subroutine json_member_write(this, unit, level, separator, count)
    type(json_member_t), intent(in) :: this
    integer,   optional, intent(in) :: unit
    integer,   optional, intent(in) :: level
    character, optional, intent(in) :: separator
    integer,   optional, intent(in) :: count
    !
    if(json_member_isdef(this))then
      call json_string_write(this%ident, unit)
      call json_write_string(": ", unit)
      call json_value_write(this%value, unit, level, separator, count)
    end if
    return
  end subroutine json_member_write

  elemental subroutine json_object_iterator_nullify(this)
    type(json_object_iterator_t), intent(out) :: this
    !
    this%pos=0
    this%size=0
    this%node=>null()
    this%table=>null()
    return
  end subroutine json_object_iterator_nullify

  elemental function json_object_iterator_isdef(this) result(is)
    type(json_object_iterator_t), intent(in) :: this
    !
    logical :: is
    !
    is=associated(this%table)
    return
  end function json_object_iterator_isdef

  subroutine json_object_iterator_init(this, object)
    type(json_object_iterator_t), intent(out) :: this
    type(json_object_t),          intent(in)  :: object
    !
    integer :: i
    !
    call json_object_iterator_nullify(this)
    if(object%used>0)then
      this%size=object%size
      this%table=>object%table
      do i = 1, object%size
        if(associated(object%table(i)%head))then
          this%pos=i
          this%node=>object%table(i)%head
          exit
        end if
      end do
    end if
    return
  end subroutine json_object_iterator_init

  elemental subroutine json_object_iterator_end(this)
    type(json_object_iterator_t), intent(inout) :: this
    !
    this%pos=0
    this%node=>null()
    this%table=>null()
    return
  end subroutine json_object_iterator_end

  subroutine json_object_iterator_next(this, member)
    type(json_object_iterator_t), intent(inout) :: this
    type(json_member_t),          pointer       :: member
    !
    integer :: i
    !
    member=>null()
    if(associated(this%node))then
      member=>this%node%member
      if(associated(this%node%next))then
        this%node=>this%node%next
      else
        this%node=>null()
        do i = this%pos+1, this%size
          if(associated(this%table(i)%head))then
            this%pos=i
            this%node=>this%table(i)%head
            exit
          end if
        end do
      end if
    end if
    return
  end subroutine json_object_iterator_next

  elemental subroutine json_object_nullify(this)
    type(json_object_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%size=0
    this%used=0
    this%table=>null()
    return
  end subroutine json_object_nullify

  elemental function json_object_isdef(this) result(is)
    type(json_object_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type==JSON_OBJECT_TYPE)
    return
  end function json_object_isdef

  subroutine json_object_init(this)
    type(json_object_t), intent(out) :: this
    !
    integer :: i
    !
    call json_object_nullify(this)
    this%used=0
    this%size=JSON_TABLE_INIT_LEN
    SAFE_ALLOCATE(this%table(this%size))
    forall(i=1:this%size)this%table(i)%head=>null()
    this%type=JSON_OBJECT_TYPE
    return
  end subroutine json_object_init

  subroutine json_object_end(this)
    type(json_object_t), intent(inout) :: this
    !
    type(json_member_node_t), pointer :: node
    integer                           :: i
    !
    this%type=JSON_UNDEF_TYPE
    do i = 1, this%size
      do
        node=>this%table(i)%head
        if(.not.associated(node))exit
        call json_member_end(node%member)
        SAFE_DEALLOCATE_P(node%member)
        node%member=>null()
        this%table(i)%head=>node%next
        this%used=this%used-1
        SAFE_DEALLOCATE_P(node)
        node=>null()
      end do
    end do
    this%used=0
    this%size=0
    SAFE_DEALLOCATE_P(this%table)
    this%table=>null()
    return
  end subroutine json_object_end

  subroutine json_object_reallocate(this)
    type(json_object_t), intent(inout) :: this
    !
    type(json_object_t)               :: buff
    type(json_member_node_t), pointer :: node
    real(kind=wp)                     :: need
    integer                           :: i, n
    !
    need=real(this%used+1, kind=wp)
    if(this%size<int(JSON_TABLE_GROWTH_FACTOR*need))then
      buff%used=0
      n=max(ceiling((log(need)-log(real(this%size,kind=wp)))/log(JSON_TABLE_GROWTH_FACTOR)),1)
      buff%size=ceiling((JSON_TABLE_GROWTH_FACTOR**n)*real(this%size,kind=wp))
      SAFE_ALLOCATE(buff%table(buff%size))
      forall(i=1:buff%size)buff%table(i)%head=>null()
      do i = 1, this%size
        do
          node=>this%table(i)%head
          if(.not.associated(node))exit
          call json_object_set_member(buff, node%member)
          this%table(i)%head=>node%next
          this%used=this%used-1
          SAFE_DEALLOCATE_P(node)
          node=>null()
        end do
      end do
      this%size=buff%size
      SAFE_DEALLOCATE_P(this%table)
      this%table=>buff%table
    end if
    return
  end subroutine json_object_reallocate

  elemental function json_object_len(this) result(len)
    type(json_object_t), intent(in) :: this
    !
    integer :: len
    !
    len=this%used
    return
  end function json_object_len

  function json_object_equal(this_1, this_2) result(eqv)
    type(json_object_t), intent(in) :: this_1
    type(json_object_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    type(json_object_iterator_t) :: iter
    type(json_member_t), pointer :: member
    type(json_value_t),  pointer :: value
    integer                      :: ierr
    !
    eqv=.false.
    if(json_object_isdef(this_1).and.json_object_isdef(this_2))then
      if(this_1%used==this_2%used)then
        eqv=.true.
        call json_object_iterator_init(iter, this_1)
        call json_object_iterator_next(iter, member)
        do while(associated(member))
          value=>null()
          call json_object_get_value(this_2, member%ident, value, ierr)
          if(ierr/=JSON_OK)then
            eqv=.false.
            exit
          end if
          if(.not.json_value_equal(member%value, value))then
            eqv=.false.
            exit
          end if
          call json_object_iterator_next(iter, member)
        end do
      end if
      call json_object_iterator_end(iter)
    end if
    return
  end function json_object_equal

  !Daniel J. Bernstein Hash Function
  elemental function json_object_hash(string, size) result(hash)
    type(json_string_t), intent(in) :: string
    integer,             intent(in) :: size
    !
    integer :: hash
    !
    integer :: i
    !
    hash=5381
    do i = 1, string%len
      hash = ieor(33*hash, iachar(string%value(i)))
    end do
    hash=modulo(hash, size)+1
    return
  end function json_object_hash

  subroutine json_object_string(this, string)
    type(json_object_t), intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    type(json_member_node_t), pointer :: node
    integer                           :: i
    !
    if(json_object_isdef(this))then
      call json_string_append(string, "{")
      do i = 1, this%size
        node=>this%table(i)%head
        do while(associated(node))
          call json_member_string(node%member, string)
          if(.not.associated(node%next))exit
          node=>node%next
          call json_string_append(string, ",")
        end do
      end do
      call json_string_append(string, "}")
    end if
    return
  end subroutine json_object_string

  subroutine json_object_write(this, unit, level, separator, count)
    type(json_object_t), intent(in) :: this
    integer,   optional, intent(in) :: unit
    integer,   optional, intent(in) :: level
    character, optional, intent(in) :: separator
    integer,   optional, intent(in) :: count
    !
    type(json_object_iterator_t) :: iter
    type(json_member_t), pointer :: member
    character                    :: sep
    integer                      :: lvl, cnt
    !
    if(json_object_isdef(this))then
      lvl=0
      cnt=1
      sep=" "
      if(present(level))lvl=level
      if(present(count))cnt=count
      if(present(separator))sep=separator
      call json_write_string("{", unit)
      call json_object_iterator_init(iter, this)
      call json_object_iterator_next(iter, member)
      if(associated(member))then
        do
          call json_write_line(unit)
          call json_write_string(repeat(sep, cnt*(lvl+1)), unit)
          call json_member_write(member, unit, lvl+1, separator, count)
          call json_object_iterator_next(iter, member)
          if(.not.associated(member))exit
          call json_write_string(",", unit)
        end do
        call json_write_line(unit)
        call json_write_string(repeat(sep, cnt*lvl), unit)
      end if
      call json_object_iterator_end(iter)
      call json_write_string("}", unit)
    end if
    return
  end subroutine json_object_write

  subroutine json_object_pop(this, member)
    type(json_object_t), intent(inout) :: this
    type(json_member_t), pointer       :: member
    !
    type(json_member_node_t), pointer :: node
    integer                           :: i
    !
    member=>null()
    do i = 1, this%size
      if(associated(this%table(i)%head))then
        node=>this%table(i)%head
        this%table(i)%head=>node%next
        this%used=this%used-1
        member=>node%member
        SAFE_DEALLOCATE_P(node)
        node=>null()
        exit
      end if
    end do
    return
  end subroutine json_object_pop

  subroutine json_object_set_member(this, member)
    type(json_object_t),         intent(inout) :: this
    type(json_member_t), target, intent(in)    :: member
    !
    type(json_member_node_t), pointer :: node
    integer                           :: n
    !
    if(json_member_isdef(member))then
      call json_object_reallocate(this)
      n=json_object_hash(member%ident, this%size)
      if(associated(this%table(n)%head))then
        node=>this%table(n)%head
        do
          if(json_string_equal(node%member%ident, member%ident))then
            call json_member_end(node%member)
            SAFE_DEALLOCATE_P(node%member)
            node%member=>member
            exit
          end if
          if(associated(node%next))then
            node=>node%next
          else
            SAFE_ALLOCATE(node%next)
            node=>node%next
            node%member=>member
            node%next=>null()
            this%used=this%used+1
            exit
          end if
        end do
      else
        SAFE_ALLOCATE(node)
        node%member=>member
        node%next=>null()
        this%table(n)%head=>node
        this%used=this%used+1
      end if
    end if
    return
  end subroutine json_object_set_member

  subroutine json_object_set_value(this, ident, value)
    type(json_object_t), intent(inout) :: this
    type(json_string_t), intent(in)    :: ident
    type(json_value_t),  intent(in)    :: value
    !
    type(json_member_t), pointer :: member
    !
    if(json_value_isdef(value))then
      SAFE_ALLOCATE(member)
      call json_member_init(member, ident, value)
      call json_object_set_member(this, member)
    end if
    return
  end subroutine json_object_set_value

  subroutine json_object_set_null(this, ident)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: json_value
    type(json_null_t),   pointer :: type_value
    !
    SAFE_ALLOCATE(string)
    call json_string_init(string, ident)
    SAFE_ALLOCATE(type_value)
    call json_null_init(type_value)
    SAFE_ALLOCATE(json_value)
    call json_value_init_null(json_value, type_value)
    call json_object_set_value(this, string, json_value)
    return
  end subroutine json_object_set_null

  subroutine json_object_set_logical(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    logical,             intent(in)    :: value
    !
    type(json_string_t),  pointer :: string
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: type_value
    !
    SAFE_ALLOCATE(string)
    call json_string_init(string, ident)
    SAFE_ALLOCATE(type_value)
    call json_logical_init(type_value, value)
    SAFE_ALLOCATE(json_value)
    call json_value_init_logical(json_value, type_value)
    call json_object_set_value(this, string, json_value)
    return
  end subroutine json_object_set_logical

  subroutine json_object_set_integer(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    integer,             intent(in)    :: value
    !
    type(json_string_t),  pointer :: string
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: type_value
    !
    SAFE_ALLOCATE(string)
    call json_string_init(string, ident)
    SAFE_ALLOCATE(type_value)
    call json_integer_init_integer(type_value, value)
    SAFE_ALLOCATE(json_value)
    call json_value_init_integer(json_value, type_value)
    call json_object_set_value(this, string, json_value)
    return
  end subroutine json_object_set_integer

  subroutine json_object_set_real(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    real(kind=wp),       intent(in)    :: value
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: json_value
    type(json_real_t),   pointer :: type_value
    !
    SAFE_ALLOCATE(string)
    call json_string_init(string, ident)
    SAFE_ALLOCATE(type_value)
    call json_real_init_real(type_value, value)
    SAFE_ALLOCATE(json_value)
    call json_value_init_real(json_value, type_value)
    call json_object_set_value(this, string, json_value)
    return
  end subroutine json_object_set_real

  subroutine json_object_set_string(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    character(len=*),    intent(in)    :: value
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: type_value
    !
    SAFE_ALLOCATE(string)
    call json_string_init(string, ident)
    SAFE_ALLOCATE(type_value)
    call json_string_init(type_value, value)
    SAFE_ALLOCATE(json_value)
    call json_value_init_string(json_value, type_value)
    call json_object_set_value(this, string, json_value)
    return
  end subroutine json_object_set_string

  subroutine json_object_set_array(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    type(json_array_t),  intent(in)    :: value
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: json_value
    !
    if(json_array_isdef(value))then
      SAFE_ALLOCATE(string)
      call json_string_init(string, ident)
      SAFE_ALLOCATE(json_value)
      call json_value_init_array(json_value, value)
      call json_object_set_value(this, string, json_value)
    end if
    return
  end subroutine json_object_set_array

  subroutine json_object_set_array_logical(this, ident, vals)
    type(json_object_t),   intent(inout) :: this
    character(len=*),      intent(in)    :: ident
    logical, dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    SAFE_ALLOCATE(json_array)
    call json_array_init_logical(json_array, vals)
    call json_object_set_array(this, ident, json_array)
    return
  end subroutine json_object_set_array_logical

  subroutine json_object_set_array_integer(this, ident, vals)
    type(json_object_t),   intent(inout) :: this
    character(len=*),      intent(in)    :: ident
    integer, dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    SAFE_ALLOCATE(json_array)
    call json_array_init_integer(json_array, vals)
    call json_object_set_array(this, ident, json_array)
    return
  end subroutine json_object_set_array_integer

  subroutine json_object_set_array_real(this, ident, vals)
    type(json_object_t),         intent(inout) :: this
    character(len=*),            intent(in)    :: ident
    real(kind=wp), dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    SAFE_ALLOCATE(json_array)
    call json_array_init_real(json_array, vals)
    call json_object_set_array(this, ident, json_array)
    return
  end subroutine json_object_set_array_real

  subroutine json_object_set_array_string(this, ident, vals)
    type(json_object_t),            intent(inout) :: this
    character(len=*),               intent(in)    :: ident
    character(len=*), dimension(:), intent(in)    :: vals
    !
    type(json_array_t), pointer :: json_array
    !
    SAFE_ALLOCATE(json_array)
    call json_array_init_string(json_array, vals)
    call json_object_set_array(this, ident, json_array)
    return
  end subroutine json_object_set_array_string

  subroutine json_object_set_object(this, ident, value)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: ident
    type(json_object_t), intent(in)    :: value
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: json_value
    !
    if(json_object_isdef(value))then
      SAFE_ALLOCATE(string)
      call json_string_init(string, ident)
      SAFE_ALLOCATE(json_value)
      call json_value_init_object(json_value, value)
      call json_object_set_value(this, string, json_value)
    end if
    return
  end subroutine json_object_set_object

  subroutine json_object_get_value(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    type(json_string_t), intent(in)  :: ident
    type(json_value_t),  pointer     :: value
    integer,             intent(out) :: ierr
    !
    type(json_member_node_t), pointer :: node
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this).and.json_string_isdef(ident))then
      ierr=JSON_KEY_ERROR
      node=>this%table(json_object_hash(ident, this%size))%head
      do while(associated(node))
        if(json_string_equal(node%member%ident, ident))then
          if(json_value_isdef(node%member%value))value=>node%member%value
          exit
        end if
        node=>node%next
      end do
      if(associated(value))ierr=JSON_OK
    end if
    return
  end subroutine json_object_get_value

  subroutine json_object_get_null(this, ident, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    integer,             intent(out) :: ierr
    !
    type(json_string_t)         :: json_ident
    type(json_value_t), pointer :: json_value
    type(json_null_t),  pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_null(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_null_get(json_tpval, ierr)
      end if
    end if
    return
  end subroutine json_object_get_null

  subroutine json_object_get_logical(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    logical,             intent(out) :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)           :: json_ident
    type(json_value_t),   pointer :: json_value
    type(json_logical_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_logical(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_logical_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_object_get_logical

  subroutine json_object_get_integer(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    integer,             intent(out) :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)           :: json_ident
    type(json_value_t),   pointer :: json_value
    type(json_integer_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_integer(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_integer_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_object_get_integer

  subroutine json_object_get_real(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    real(kind=wp),       intent(out) :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)         :: json_ident
    type(json_value_t), pointer :: json_value
    type(json_real_t),  pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_real(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_real_get(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_object_get_real

  subroutine json_object_get_string(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    character(len=*),    intent(out) :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)          :: json_ident
    type(json_value_t),  pointer :: json_value
    type(json_string_t), pointer :: json_tpval
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK)then
        json_tpval=>null()
        call json_value_get_string(json_value, json_tpval, ierr)
        if(ierr==JSON_OK) call json_string_get_string(json_tpval, value, ierr)
      end if
    end if
    return
  end subroutine json_object_get_string

  subroutine json_object_get_array(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    type(json_array_t),  pointer     :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)         :: json_ident
    type(json_value_t), pointer :: json_value
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK) call json_value_get_array(json_value, value, ierr)
    end if
    return
  end subroutine json_object_get_array

  subroutine json_object_get_array_logical(this, ident, vals, ierr)
    type(json_object_t),   intent(in)  :: this
    character(len=*),      intent(in)  :: ident
    logical, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      json_array=>null()
      call json_object_get_array(this, ident, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_logical(json_array, vals, ierr)
    end if
    return
  end subroutine json_object_get_array_logical

  subroutine json_object_get_array_integer(this, ident, vals, ierr)
    type(json_object_t),   intent(in)  :: this
    character(len=*),      intent(in)  :: ident
    integer, dimension(:), intent(out) :: vals
    integer,               intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      json_array=>null()
      call json_object_get_array(this, ident, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_integer(json_array, vals, ierr)
    end if
    return
  end subroutine json_object_get_array_integer

  subroutine json_object_get_array_real(this, ident, vals, ierr)
    type(json_object_t),         intent(in)  :: this
    character(len=*),            intent(in)  :: ident
    real(kind=wp), dimension(:), intent(out) :: vals
    integer,                     intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      json_array=>null()
      call json_object_get_array(this, ident, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_real(json_array, vals, ierr)
    end if
    return
  end subroutine json_object_get_array_real

  subroutine json_object_get_array_string(this, ident, vals, ierr)
    type(json_object_t),            intent(in)  :: this
    character(len=*),               intent(in)  :: ident
    character(len=*), dimension(:), intent(out) :: vals
    integer,                        intent(out) :: ierr
    !
    type(json_array_t), pointer :: json_array
    !
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      json_array=>null()
      call json_object_get_array(this, ident, json_array, ierr)
      if((ierr==JSON_OK).and.json_array_isdef(json_array))&
        call json_array_get_self_string(json_array, vals, ierr)
    end if
    return
  end subroutine json_object_get_array_string

  subroutine json_object_get_object(this, ident, value, ierr)
    type(json_object_t), intent(in)  :: this
    character(len=*),    intent(in)  :: ident
    type(json_object_t), pointer     :: value
    integer,             intent(out) :: ierr
    !
    type(json_string_t)         :: json_ident
    type(json_value_t), pointer :: json_value
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this))then
      call json_string_init(json_ident, ident)
      json_value=>null()
      call json_object_get_value(this, json_ident, json_value, ierr)
      call json_string_end(json_ident)
      if(ierr==JSON_OK) call json_value_get_object(json_value, value, ierr)
    end if
    return
  end subroutine json_object_get_object

  elemental subroutine json_value_nullify(this)
    type(json_value_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%null=>null()
    this%logical=>null()
    this%integer=>null()
    this%real=>null()
    this%string=>null()
    this%array=>null()
    this%object=>null()
    return
  end subroutine json_value_nullify

  elemental function json_value_isdef(this) result(is)
    type(json_value_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type/=JSON_UNDEF_TYPE)
    return
  end function json_value_isdef

  elemental function json_value_type(this) result(id)
    type(json_value_t), intent(in) :: this
    !
    integer :: id
    !
    id=this%type
    return
  end function json_value_type

  subroutine json_value_init_null(this, value)
    type(json_value_t),        intent(out) :: this
    type(json_null_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%null=>value
    this%type=JSON_NULL_TYPE
    return
  end subroutine json_value_init_null

  subroutine json_value_init_logical(this, value)
    type(json_value_t),           intent(out) :: this
    type(json_logical_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%logical=>value
    this%type=JSON_LOGICAL_TYPE
    return
  end subroutine json_value_init_logical

  subroutine json_value_init_integer(this, value)
    type(json_value_t),           intent(out) :: this
    type(json_integer_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%integer=>value
    this%type=JSON_INTEGER_TYPE
    return
  end subroutine json_value_init_integer

  subroutine json_value_init_real(this, value)
    type(json_value_t),        intent(out) :: this
    type(json_real_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%real=>value
    this%type=JSON_REAL_TYPE
    return
  end subroutine json_value_init_real

  subroutine json_value_init_string(this, value)
    type(json_value_t),          intent(out) :: this
    type(json_string_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%string=>value
    this%type=JSON_STRING_TYPE
    return
  end subroutine json_value_init_string

  subroutine json_value_init_array(this, value)
    type(json_value_t),         intent(out) :: this
    type(json_array_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%array=>value
    this%type=JSON_ARRAY_TYPE
    return
  end subroutine json_value_init_array

  subroutine json_value_init_object(this, value)
    type(json_value_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: value
    !
    call json_value_nullify(this)
    this%object=>value
    this%type=JSON_OBJECT_TYPE
    return
  end subroutine json_value_init_object

  subroutine json_value_end(this)
    type(json_value_t), intent(inout) :: this
    !
    select case(this%type)
    case(JSON_NULL_TYPE)
      call json_null_end(this%null)
      SAFE_DEALLOCATE_P(this%null)
      this%null=>null()
    case(JSON_LOGICAL_TYPE)
      call json_logical_end(this%logical)
      SAFE_DEALLOCATE_P(this%logical)
      this%logical=>null()
    case(JSON_INTEGER_TYPE)
      call json_integer_end(this%integer)
      SAFE_DEALLOCATE_P(this%integer)
      this%integer=>null()
    case(JSON_REAL_TYPE)
      call json_real_end(this%real)
      SAFE_DEALLOCATE_P(this%real)
      this%real=>null()
    case(JSON_STRING_TYPE)
      call json_string_end(this%string)
      SAFE_DEALLOCATE_P(this%string)
      this%string=>null()
    case(JSON_ARRAY_TYPE)
      call json_array_end(this%array)
      SAFE_DEALLOCATE_P(this%array)
      this%array=>null()
    case(JSON_OBJECT_TYPE)
      call json_object_end(this%object)
      SAFE_DEALLOCATE_P(this%object)
      this%object=>null()
    end select
    this%type=JSON_UNDEF_TYPE
    return
  end subroutine json_value_end

  function json_value_equal(this_1, this_2) result(eqv)
    type(json_value_t), intent(in) :: this_1
    type(json_value_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_value_isdef(this_1).and.json_value_isdef(this_2))then
      if(this_1%type==this_2%type)then
        select case(this_1%type)
        case(JSON_NULL_TYPE)
          eqv=json_null_equal(this_1%null, this_2%null)
        case(JSON_LOGICAL_TYPE)
          eqv=json_logical_equal(this_1%logical, this_2%logical)
        case(JSON_INTEGER_TYPE)
          eqv=json_integer_equal(this_1%integer, this_2%integer)
        case(JSON_REAL_TYPE)
          eqv=json_real_equal(this_1%real, this_2%real)
        case(JSON_STRING_TYPE)
          eqv=json_string_equal(this_1%string, this_2%string)
        case(JSON_ARRAY_TYPE)
          eqv=json_array_equal(this_1%array, this_2%array)
        case(JSON_OBJECT_TYPE)
          eqv=json_object_equal(this_1%object, this_2%object)
        end select
      end if
    end if
    return
  end function json_value_equal

  subroutine json_value_string(this, string)
    type(json_value_t),  intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    select case(this%type)
    case(JSON_NULL_TYPE)
      call json_null_string(this%null, string)
    case(JSON_LOGICAL_TYPE)
      call json_logical_string(this%logical, string)
    case(JSON_INTEGER_TYPE)
      call json_integer_string(this%integer, string)
    case(JSON_REAL_TYPE)
      call json_real_string(this%real, string)
    case(JSON_STRING_TYPE)
      call json_string_string(this%string, string)
    case(JSON_ARRAY_TYPE)
      call json_array_string(this%array, string)
    case(JSON_OBJECT_TYPE)
      call json_object_string(this%object, string)
    end select
    return
  end subroutine json_value_string

  subroutine json_value_get_null(this, value, ierr)
    type(json_value_t), intent(in)  :: this
    type(json_null_t),  pointer     :: value
    integer,            intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_null_isdef(this%null))then
      value=>this%null
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_null

  subroutine json_value_get_logical(this, value, ierr)
    type(json_value_t),   intent(in)  :: this
    type(json_logical_t), pointer     :: value
    integer,              intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_logical_isdef(this%logical))then
      value=>this%logical
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_logical

  subroutine json_value_get_integer(this, value, ierr)
    type(json_value_t),   intent(in)  :: this
    type(json_integer_t), pointer     :: value
    integer,              intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_integer_isdef(this%integer))then
      value=>this%integer
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_integer

  subroutine json_value_get_real(this, value, ierr)
    type(json_value_t), intent(in)  :: this
    type(json_real_t),  pointer     :: value
    integer,            intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_real_isdef(this%real))then
      value=>this%real
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_real

  subroutine json_value_get_string(this, value, ierr)
    type(json_value_t),  intent(in)  :: this
    type(json_string_t), pointer     :: value
    integer,             intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_string_isdef(this%string))then
      value=>this%string
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_string

  subroutine json_value_get_array(this, value, ierr)
    type(json_value_t), intent(in)  :: this
    type(json_array_t), pointer     :: value
    integer,            intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this%array))then
      value=>this%array
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_array

  subroutine json_value_get_object(this, value, ierr)
    type(json_value_t),  intent(in)  :: this
    type(json_object_t), pointer     :: value
    integer,             intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this%object))then
      value=>this%object
      ierr=JSON_OK
    end if
    return
  end subroutine json_value_get_object

  subroutine json_value_write(this, unit, level, separator, count)
    type(json_value_t),  intent(in) :: this
    integer,   optional, intent(in) :: unit
    integer,   optional, intent(in) :: level
    character, optional, intent(in) :: separator
    integer,   optional, intent(in) :: count
    !
    select case(this%type)
    case(JSON_NULL_TYPE)
      call json_null_write(this%null, unit)
    case(JSON_LOGICAL_TYPE)
      call json_logical_write(this%logical, unit)
    case(JSON_INTEGER_TYPE)
      call json_integer_write(this%integer, unit)
    case(JSON_REAL_TYPE)
      call json_real_write(this%real, unit)
    case(JSON_STRING_TYPE)
      call json_string_write(this%string, unit)
    case(JSON_ARRAY_TYPE)
      call json_array_write(this%array, unit, level, separator, count)
    case(JSON_OBJECT_TYPE)
      call json_object_write(this%object, unit, level, separator, count)
    end select
    return
  end subroutine json_value_write

  elemental subroutine json_json_nullify(this)
    type(json_t), intent(out) :: this
    !
    this%type=JSON_UNDEF_TYPE
    this%array=>null()
    this%object=>null()
    return
  end subroutine json_json_nullify

  elemental function json_json_isdef(this) result(is)
    type(json_t), intent(in) :: this
    !
    logical :: is
    !
    is=(this%type/=JSON_UNDEF_TYPE)
    return
  end function json_json_isdef

  subroutine json_json_array_init(this, value)
    type(json_t),               intent(out) :: this
    type(json_array_t), target, intent(in)  :: value
    !
    call json_json_nullify(this)
    this%array=>value
    this%type=JSON_ARRAY_TYPE
    return
  end subroutine json_json_array_init

  subroutine json_json_object_init(this, value)
    type(json_t),                intent(out) :: this
    type(json_object_t), target, intent(in)  :: value
    !
    call json_json_nullify(this)
    this%object=>value
    this%type=JSON_OBJECT_TYPE
    return
  end subroutine json_json_object_init

  subroutine json_json_end(this)
    type(json_t), intent(inout) :: this
    !
    select case(this%type)
    case(JSON_ARRAY_TYPE)
      call json_array_end(this%array)
      SAFE_DEALLOCATE_P(this%array)
      this%array=>null()
    case(JSON_OBJECT_TYPE)
      call json_object_end(this%object)
      SAFE_DEALLOCATE_P(this%object)
      this%object=>null()
    end select
    this%type=JSON_UNDEF_TYPE
    return
  end subroutine json_json_end

  function json_json_equal(this_1, this_2) result(eqv)
    type(json_t), intent(in) :: this_1
    type(json_t), intent(in) :: this_2
    !
    logical :: eqv
    !
    eqv=.false.
    if(json_json_isdef(this_1).and.json_json_isdef(this_2))then
      if(this_1%type==this_2%type)then
        select case(this_1%type)
        case(JSON_ARRAY_TYPE)
          eqv=json_array_equal(this_1%array, this_2%array)
        case(JSON_OBJECT_TYPE)
          eqv=json_object_equal(this_1%object, this_2%object)
        end select
      end if
    end if
    return
  end function json_json_equal

  subroutine json_json_string(this, string)
    type(json_t),        intent(in)    :: this
    type(json_string_t), intent(inout) :: string
    !
    select case(this%type)
    case(JSON_ARRAY_TYPE)
      call json_array_string(this%array, string)
    case(JSON_OBJECT_TYPE)
      call json_object_string(this%object, string)
    end select
    return
  end subroutine json_json_string

  subroutine json_json_write(this, unit, separator, count)
    type(json_t),        intent(in) :: this
    integer,   optional, intent(in) :: unit
    character, optional, intent(in) :: separator
    integer,   optional, intent(in) :: count
    !
    !application/json
    select case(this%type)
    case(JSON_ARRAY_TYPE)
      call json_array_write(this%array, unit, separator=separator, count=count)
    case(JSON_OBJECT_TYPE)
      call json_object_write(this%object, unit, separator=separator, count=count)
    end select
    call json_write_line(unit)
    return
  end subroutine json_json_write

  subroutine json_json_get_array(this, value, ierr)
    type(json_t),       intent(in)  :: this
    type(json_array_t), pointer     :: value
    integer,            intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_array_isdef(this%array))then
      value=>this%array
      ierr=JSON_OK
    end if
    return
  end subroutine json_json_get_array

  subroutine json_json_get_object(this, value, ierr)
    type(json_t),        intent(in)  :: this
    type(json_object_t), pointer     :: value
    integer,             intent(out) :: ierr
    !
    value=>null()
    ierr=JSON_UNDEF_ERROR
    if(json_object_isdef(this%object))then
      value=>this%object
      ierr=JSON_OK
    end if
    return
  end subroutine json_json_get_object

end module json_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
