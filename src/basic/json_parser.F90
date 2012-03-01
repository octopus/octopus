#include "global.h"

module json_parser_m

  use global_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::            &
    json_parser_init,  &
    json_parser_parse, &
    json_parser_error, &
    json_parser_end

  integer, parameter :: LINE_LEN=1024

  character, parameter :: backslash="\" !"
  character, parameter :: space=" "
  character, parameter :: bspace=achar(8)
  character, parameter :: tab=achar(9)
  character, parameter :: newline=achar(10)
  character, parameter :: vt=achar(11)
  character, parameter :: formfeed=achar(12)
  character, parameter :: creturn=achar(13)

  character(len=*), parameter :: whitespace=space//tab//newline//vt//formfeed//creturn
  character(len=*), parameter :: ualph="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  character(len=*), parameter :: lalph="abcdefghijklmnopqrstuvwxyz"

  character(len=*), parameter :: numbers="-+.0123456789"

  type, public :: json_parser_t
    private
    integer             :: unit=-1
    integer             :: line=0
    integer             :: bpos=0
    integer             :: ierr=0
    type(json_string_t) :: buff
  end type json_parser_t

contains

  subroutine json_parser_init(parser, unit, iostat)
    type(json_parser_t), intent(inout) :: parser
    integer,             intent(in)    :: unit
    integer,             intent(out)   :: iostat
    !
    parser%unit=unit
    parser%line=1
    parser%bpos=1
    parser%ierr=0
    call json_parser_readline(parser, iostat)
    return
  end subroutine json_parser_init

  subroutine json_parser_end(parser)
    type(json_parser_t), intent(inout) :: parser
    !
    parser%unit=-1
    parser%line=0
    parser%bpos=0
    parser%ierr=0
    call json_end(parser%buff)
    return
  end subroutine json_parser_end

  elemental function json_parser_is_space(char) result(iss)
    character, intent(in) :: char
    !
    logical :: iss
    !
    iss=(scan(whitespace, char)/=0)
    return
  end function json_parser_is_space

  elemental function json_parser_to_lower(uchr) result(lchr)
    character, intent(in) :: uchr
    !
    character :: lchr
    !
    integer :: i
    !
    i=index(ualph, uchr)
    if(i>0)then
      lchr=lalph(i:i)
    else
      lchr=uchr
    end if
    return
  end function json_parser_to_lower

  elemental function json_parser_to_upper(lchr) result(uchr)
    character, intent(in) :: lchr
    !
    character :: uchr
    !
    integer :: i
    !
    i=index(lalph, lchr)
    if(i>0)then
      uchr=ualph(i:i)
    else
      uchr=lchr
    end if
    return
  end function json_parser_to_upper

  subroutine json_parser_readline(parser, iostat)
    type(json_parser_t), intent(inout) :: parser
    integer,             intent(out)   :: iostat
    !
    character(len=LINE_LEN) :: buff
    !
    read(unit=parser%unit,fmt="(a)",iostat=iostat) buff
    parser%ierr=iostat
    if(iostat==0)then
      call json_end(parser%buff)
      call json_init(parser%buff, trim(buff))
      parser%line=parser%line+1
      parser%bpos=1
    end if
    return
  end subroutine json_parser_readline

  subroutine json_parser_get_char(parser, char, iostat)
    type(json_parser_t), intent(inout) :: parser
    character,           intent(out)   :: char
    integer,             intent(out)   :: iostat
    !
    integer :: ierr
    !
    if(parser%bpos>json_len(parser%buff)) call json_parser_readline(parser, iostat)
    if(iostat==0)then
      call json_get(parser%buff, parser%bpos, char, ierr)
      parser%bpos=parser%bpos+1
    end if
    return
  end subroutine json_parser_get_char

  subroutine json_parser_peek_char(parser, char, iostat)
    type(json_parser_t), intent(inout) :: parser
    character,           intent(out)   :: char
    integer,             intent(out)   :: iostat
    !
    call json_parser_get_char(parser, char, iostat)
    if(iostat==0)parser%bpos=parser%bpos-1
    return
  end subroutine json_parser_peek_char

  subroutine json_parser_get_nonblank_char(parser, char, iostat)
    type(json_parser_t), intent(inout) :: parser
    character,           intent(out)   :: char
    integer,             intent(out)   :: iostat
    !
    do
      call json_parser_get_char(parser, char, iostat)
      if(iostat/=0)exit
      if(.not.json_parser_is_space(char))exit
    end do
    return
  end subroutine json_parser_get_nonblank_char

  subroutine json_parser_peek_nonblank_char(parser, char, iostat)
    type(json_parser_t), intent(inout) :: parser
    character,           intent(out)   :: char
    integer,             intent(out)   :: iostat
    !
    call json_parser_get_nonblank_char(parser, char, iostat)
    if(iostat==0)parser%bpos=parser%bpos-1
    return
  end subroutine json_parser_peek_nonblank_char

  subroutine json_parser_parse_literal(parser, string, ok, iostat)
    type(json_parser_t), intent(inout) :: parser
    character(len=*),    intent(in)    :: string
    logical,             intent(out)   :: ok
    integer,             intent(out)   :: iostat
    !
    character :: char
    integer   :: i
    !
    ok=.false.
    call json_parser_get_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(json_parser_to_lower(char)==string(1:1)))then
      ok=.true.
      do i = 2, len(string)
        call json_parser_get_char(parser, char, iostat)
        if((iostat/=0).or.(json_parser_to_lower(char)/=string(i:i)))then
          ok=.false.
          exit
        end if
      end do
    end if
    return
  end subroutine json_parser_parse_literal

  subroutine json_parser_parse_null(parser, value, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_null_t),   intent(out)   :: value
    integer,             intent(out)   :: iostat
    !
    character(len=*), parameter :: data="null"
    !
    logical :: ok
    !
    call json_end(value)
    call json_parser_parse_literal(parser, data, ok, iostat)
    if((iostat==0).and.(ok)) call json_init(value)
    return
  end subroutine json_parser_parse_null

  subroutine json_parser_parse_false(parser, value, iostat)
    type(json_parser_t),  intent(inout) :: parser
    type(json_logical_t), intent(inout) :: value
    integer,              intent(out)   :: iostat
    !
    character(len=*), parameter :: data="false"
    !
    logical :: ok
    !
    call json_end(value)
    call json_parser_parse_literal(parser, data, ok, iostat)
    if((iostat==0).and.(ok)) call json_init(value, .false.)
    return
  end subroutine json_parser_parse_false

  subroutine json_parser_parse_true(parser, value, iostat)
    type(json_parser_t),  intent(inout) :: parser
    type(json_logical_t), intent(inout) :: value
    integer,              intent(out)   :: iostat
    !
    character(len=*), parameter :: data="true"
    !
    logical :: ok
    !
    call json_end(value)
    call json_parser_parse_literal(parser, data, ok, iostat)
    if((iostat==0).and.(ok)) call json_init(value, .true.)
    return
  end subroutine json_parser_parse_true

  subroutine json_parser_parse_number(parser, string, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_string_t), intent(out)   :: string
    integer,             intent(out)   :: iostat
    !
    character :: char
    !
    call json_end(string)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(scan(numbers,char)/=0))then
      call json_init(string)
      do
        call json_parser_get_char(parser, char, iostat)
        if((iostat/=0).or.(scan(numbers//"eE",char)==0))exit
        call json_append(string, char)
        call json_parser_peek_char(parser, char, iostat)
        if((iostat/=0).or.(scan(numbers//"eE",char)==0))exit
      end do
      if(iostat/=0)call json_end(string)
    end if
    return
  end subroutine json_parser_parse_number

  subroutine json_parser_handle_special(parser, char, iostat)
    type(json_parser_t), intent(inout) :: parser
    character,           intent(out)   :: char
    integer,             intent(out)   :: iostat
    !
    character :: ichr
    !
    call json_parser_peek_char(parser, ichr, iostat)
    if(iostat==0)then
      select case(ichr)
      case('"') 
        char='"'
        call json_parser_get_char(parser, ichr, iostat)
      case(backslash) 
        char=backslash
        call json_parser_get_char(parser, ichr, iostat)
      case("/")
        char="/"
        call json_parser_get_char(parser, ichr, iostat)
      case("b")
        char=bspace
        call json_parser_get_char(parser, ichr, iostat)
      case("f")
        char=formfeed
        call json_parser_get_char(parser, ichr, iostat)
      case("n")
        char=newline
        call json_parser_get_char(parser, ichr, iostat)
      case("r")
        char=creturn
        call json_parser_get_char(parser, ichr, iostat)
      case("t")
        char=tab
        call json_parser_get_char(parser, ichr, iostat)
      !case("u") no unicode...
      case default
        char=backslash
      end select
    end if
    return
  end subroutine json_parser_handle_special

   subroutine json_parser_parse_string(parser, string, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_string_t), intent(out)   :: string
    integer,             intent(out)   :: iostat
    !
    character :: char
    !
    call json_end(string)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(char=='"'))then
      call json_parser_get_char(parser, char, iostat)
      call json_init(string)
      do
        call json_parser_get_char(parser, char, iostat)
        if(iostat/=0)exit
        select case(char)
        case(backslash)
          call json_parser_handle_special(parser, char, iostat)
        case('"')
          exit
        end select
        if(iostat/=0)exit
        call json_append(string, char)
      end do
      if((iostat/=0).or.(char/='"'))call json_end(string)
    end if
    return
  end subroutine json_parser_parse_string

  recursive subroutine json_parser_parse_value(parser, value, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_value_t),  intent(out)   :: value
    integer,             intent(out)   :: iostat
    !
    type(json_null_t),    pointer :: nvalue
    type(json_logical_t), pointer :: lvalue
    type(json_integer_t), pointer :: ivalue
    type(json_real_t) ,   pointer :: rvalue
    type(json_string_t),  pointer :: string
    type(json_array_t),   pointer :: array
    type(json_object_t),  pointer :: object
    character                     :: char
    !
    call json_end(value)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if(iostat==0)then
      select case(json_parser_to_lower(char))
      case("n")
        SAFE_ALLOCATE(nvalue)
        call json_parser_parse_null(parser, nvalue, iostat)
        if((iostat==0).and.(json_isdef(nvalue)))then
          call json_init(value, nvalue)
        else
          call json_end(nvalue)
          SAFE_DEALLOCATE_P(nvalue)
        end if
        nullify(nvalue)
      case("f")
        SAFE_ALLOCATE(lvalue)
        call json_parser_parse_false(parser, lvalue, iostat)
        if((iostat==0).and.(json_isdef(lvalue)))then
          call json_init(value, lvalue)
        else
          call json_end(lvalue)
          SAFE_DEALLOCATE_P(lvalue)
        end if
        nullify(lvalue)
      case("t")
        SAFE_ALLOCATE(lvalue)
        call json_parser_parse_true(parser, lvalue, iostat)
        if((iostat==0).and.(json_isdef(lvalue)))then
          call json_init(value, lvalue)
        else
          call json_end(lvalue)
          SAFE_DEALLOCATE_P(lvalue)
        end if
        nullify(lvalue)
      case("-","+",".","0","1","2","3","4","5","6","7","8","9")
        SAFE_ALLOCATE(string)
        call json_parser_parse_number(parser, string, iostat)
        if((iostat==0).and.(json_isdef(string)))then
          if(scan(string,".eE")==0)then
            SAFE_ALLOCATE(ivalue)
            call json_init(ivalue, string)
            if(json_isdef(ivalue))then
              call json_init(value, ivalue)
            else
              SAFE_DEALLOCATE_P(ivalue)
            end if
            nullify(ivalue)
          else
            SAFE_ALLOCATE(rvalue)
            call json_init(rvalue, string)
            if(json_isdef(rvalue))then
              call json_init(value, rvalue)
            else
              SAFE_DEALLOCATE_P(rvalue)
            end if
            nullify(rvalue)
          end if
        end if
        call json_end(string)
        SAFE_DEALLOCATE_P(string)
        nullify(string)
      case('"')
        SAFE_ALLOCATE(string)
        call json_parser_parse_string(parser, string, iostat)
        if((iostat==0).and.(json_isdef(string)))then
          call json_init(value, string)
        else
          call json_end(string)
          SAFE_DEALLOCATE_P(string)
        end if
        nullify(string)
      case("[")
        SAFE_ALLOCATE(array)
        call json_parser_parse_array(parser, array, iostat)
        if((iostat==0).and.(json_isdef(array)))then
          call json_init(value, array)
        else
          call json_end(array)
          SAFE_DEALLOCATE_P(array)
        end if
        nullify(array)
      case("{")
        SAFE_ALLOCATE(object)
        call json_parser_parse_object(parser, object, iostat)
        if((iostat==0).and.(json_isdef(object)))then
          call json_init(value, object)
        else
          call json_end(object)
          SAFE_DEALLOCATE_P(object)
        end if
        nullify(object)
      end select
    end if
    return
  end subroutine json_parser_parse_value

  recursive subroutine json_parser_parse_array(parser, array, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_array_t),  intent(out)   :: array
    integer,             intent(out)   :: iostat
    !
    type(json_value_t), pointer :: value
    character                   :: char
    !
    call json_end(array)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(char=="["))then
      call json_parser_get_nonblank_char(parser, char, iostat)
      call json_parser_peek_nonblank_char(parser, char, iostat)
      if(iostat==0)then
        call json_init(array)
        if(char=="]")then
          call json_parser_get_nonblank_char(parser, char, iostat)
        else
          do
            SAFE_ALLOCATE(value)
            call json_parser_parse_value(parser, value, iostat)
            if(json_isdef(value))then
              call json_append(array, value)
              call json_parser_get_nonblank_char(parser, char, iostat)
              if((iostat/=0).or.(char/=","))exit
            else
              call json_end(value)
              SAFE_DEALLOCATE_P(value)
              call json_end(array)
              exit
            end if
            nullify(value)
          end do
        end if
        if((.not.json_isdef(array)).or.(iostat/=0).or.(char/="]"))call json_end(array)
      end if
    end if
    return
  end subroutine json_parser_parse_array

  recursive subroutine json_parser_parse_member(parser, member, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_member_t), intent(out)   :: member
    integer,             intent(out)   :: iostat
    !
    type(json_string_t), pointer :: string
    type(json_value_t),  pointer :: value
    character                    :: char
    !
    call json_end(member)
    SAFE_ALLOCATE(string)
    call json_parser_parse_string(parser, string, iostat)
    if((iostat==0).and.(json_isdef(string)))then
      call json_parser_get_nonblank_char(parser, char, iostat)
      if((iostat==0).and.(char==':'))then
        SAFE_ALLOCATE(value)
        call json_parser_parse_value(parser, value, iostat)
        if((iostat==0).and.(json_isdef(value)))then
          call json_init(member, string, value)
        else
          call json_end(value)
          SAFE_DEALLOCATE_P(value)
        end if
        nullify(value)
      end if
    end if
    if((iostat/=0).or.(.not.json_isdef(member)))then
      call json_end(string)
      SAFE_DEALLOCATE_P(string)
    end if
    nullify(string)
    return
  end subroutine json_parser_parse_member

  recursive subroutine json_parser_parse_object(parser, object, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_object_t), intent(out)   :: object
    integer,             intent(out)   :: iostat
    !
    type(json_member_t), pointer :: member
    character                    :: char
    !
    call json_end(object)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(char=="{"))then
      call json_parser_get_nonblank_char(parser, char, iostat)
      call json_parser_peek_nonblank_char(parser, char, iostat)
      if(iostat==0)then
        call json_init(object)
        if(char=="}")then
          call json_parser_get_nonblank_char(parser, char, iostat)
        else
          do
            SAFE_ALLOCATE(member)
            call json_parser_parse_member(parser, member, iostat)
            if((iostat==0).and.(json_isdef(member)))then
              call json_set(object, member)
              call json_parser_get_nonblank_char(parser, char, iostat)
              if((iostat/=0).or.(char/=","))exit
            else
              call json_end(member)
              SAFE_DEALLOCATE_P(member)
              call json_end(object)
              exit
            end if
            nullify(member)
          end do
        end if
        if((.not.json_isdef(object)).or.(iostat/=0).or.(char/="}")) call json_end(object)
      end if
    end if
    return
  end subroutine json_parser_parse_object

  subroutine json_parser_parse(parser, json, iostat)
    type(json_parser_t), intent(inout) :: parser
    type(json_t),        intent(out)   :: json
    integer,             intent(out)   :: iostat
    !
    type(json_array_t),  pointer :: array
    type(json_object_t), pointer :: object
    character                    :: char
    !
    call json_end(json)
    call json_parser_peek_nonblank_char(parser, char, iostat)
    if((iostat==0).and.(char=="["))then
      SAFE_ALLOCATE(array)
      call json_parser_parse_array(parser, array, iostat)
      if((iostat==0).and.(json_isdef(array)))then
        call json_init(json, array)
      else
        SAFE_DEALLOCATE_P(array)
      end if
      nullify(array)
    else if((iostat==0).and.(char=="{"))then
      SAFE_ALLOCATE(object)
      call json_parser_parse_object(parser, object, iostat)
      if((iostat==0).and.(json_isdef(object)))then
        call json_init(json, object)
      else
        SAFE_DEALLOCATE_P(object)
      end if
      nullify(object)
    end if
    return
  end subroutine json_parser_parse

  subroutine json_parser_error(parser, unit)
    type(json_parser_t), intent(in) :: parser
    integer,   optional, intent(in) :: unit
    !
    character(len=LINE_LEN) :: buff
    integer                 :: ierr
    !
    if(parser%ierr/=0)then
      if(present(unit))then
        write(unit=unit, fmt="(a,i3)") "I/O Error: ", parser%ierr
      else
        write(unit=*, fmt="(a,i3)") "I/O Error: ", parser%ierr
      end if
    else
      if(present(unit))then
        write(unit=unit, fmt="(a,i3)") "Parsing Error at line: ", parser%line
      else
        write(unit=*, fmt="(a,i3)") "Parsing Error at line: ", parser%line
      end if
    end if
    call json_get(parser%buff, buff, ierr)
    if(present(unit))then
      write(unit=unit, fmt="(a)") trim(buff)
      write(unit=unit, fmt="(a)") repeat("-",parser%bpos-1)//"^"
    else
      write(unit=*, fmt="(a)") trim(buff)
      write(unit=*, fmt="(a)") repeat("-",parser%bpos-1)//"^"
    end if
    return
  end subroutine json_parser_error

end module json_parser_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
