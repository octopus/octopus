!! Copyright (C) 2003-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module parser_oct_m
  use global_oct_m
  use loct_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use unit_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::              &
    block_t,             &
    parser_init,         &
    parser_end,          &
    parse_init,          &
    parse_putsym,        &
    parse_end,           &
    parse_is_defined,    &
    parse_variable,      &
    parse_block,         &
    parse_block_end,     &
    parse_block_n,       &
    parse_block_cols,    &
    parse_block_integer, &
    parse_block_float,   &
    parse_block_cmplx,   &
    parse_block_string,  &
    parse_block_logical, &
    parse_expression,    &
    parse_array

  type :: block_t
    private
    integer, pointer :: p
  end type block_t

  ! The following characters should not be allowed in variable names.
  character(len=27), parameter, public :: parser_varname_excluded_characters = '|!''"#$%&/\()=?{}+-*^.,;:<> '

  interface parse_init
    integer function oct_parse_init(file_out, mpiv_node)
      implicit none
      character(len=*), intent(in) :: file_out
      integer,          intent(in) :: mpiv_node
    end function oct_parse_init
  end interface parse_init

  interface parse_putsym
    subroutine oct_parse_putsym_int(sym, i)
      implicit none
      character(len=*), intent(in) :: sym
      integer,          intent(in) :: i
    end subroutine oct_parse_putsym_int

    subroutine oct_parse_putsym_double(sym, d)
      implicit none
      character(len=*), intent(in) :: sym
      real(8),          intent(in) :: d
    end subroutine oct_parse_putsym_double
  end interface parse_putsym

  interface parse_input_file
    integer function oct_parse_input(file_in, set_used)
      implicit none
      character(len=*), intent(in) :: file_in
      integer,          intent(in) :: set_used
    end function oct_parse_input
  end interface parse_input_file

  interface parse_environment
    subroutine oct_parse_environment(prefix)
      implicit none
      character(len=*), intent(in) :: prefix
    end subroutine oct_parse_environment
  end interface parse_environment

  interface parse_end
    subroutine oct_parse_end()
      implicit none
    end subroutine oct_parse_end
  end interface parse_end

  interface sym_output_table
    subroutine oct_sym_output_table(only_unused, mpiv_node)
      implicit none
      integer, intent(in) :: only_unused
      integer, intent(in) :: mpiv_node
    end subroutine oct_sym_output_table
  end interface sym_output_table

  interface parse_isdef
    integer function oct_parse_isdef(name)
      implicit none
      character(len=*), intent(in) :: name
    end function oct_parse_isdef
  end interface parse_isdef

  interface
    subroutine oct_parse_int(name, def, res)
      implicit none
      character(len=*), intent(in)  :: name
      integer(8),       intent(in)  :: def
      integer(8),       intent(out) :: res
    end subroutine oct_parse_int

    subroutine oct_parse_double(name, def, res)
      implicit none
      character(len=*), intent(in)  :: name
      real(8),          intent(in)  :: def
      real(8),          intent(out) :: res
    end subroutine oct_parse_double

    subroutine oct_parse_complex(name, def, res)
      implicit none
      character(len=*), intent(in)  :: name
      complex(8),       intent(in)  :: def
      complex(8),       intent(out) :: res
    end subroutine oct_parse_complex

    subroutine oct_parse_string(name, def, res)
      implicit none
      character(len=*), intent(in)  :: name, def
      character(len=*), intent(out) :: res
    end subroutine oct_parse_string

    integer function oct_parse_block(name, blk)
      import block_t
      implicit none
      character(len=*), intent(in)  :: name
      type(block_t),    intent(out) :: blk
    end function oct_parse_block
  end interface

  interface parse_variable
    module procedure parse_integer
    module procedure parse_integer8
    module procedure parse_integer48
    module procedure parse_integer84
    module procedure parse_logical
    module procedure parse_string
    module procedure parse_cmplx
    module procedure oct_parse_double_unit
  end interface parse_variable

  interface parse_block_end
    subroutine oct_parse_block_end(blk)
      import block_t
      implicit none
      type(block_t), intent(inout) :: blk
    end subroutine oct_parse_block_end
  end interface parse_block_end

  interface parse_block_n
    integer function oct_parse_block_n(blk)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
    end function oct_parse_block_n
  end interface parse_block_n

  interface parse_block_cols
    integer function oct_parse_block_cols(blk, line)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
      integer,       intent(in) :: line
    end function oct_parse_block_cols
  end interface parse_block_cols

  interface parse_block_integer
    subroutine oct_parse_block_int(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in)  :: blk
      integer,       intent(in)  :: l, c
      integer,       intent(out) :: res
    end subroutine oct_parse_block_int

    subroutine oct_parse_block_int8(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in)    :: blk
      integer, intent(in)          :: l, c
      integer(8), intent(out)      :: res
    end subroutine oct_parse_block_int8
  end interface parse_block_integer

  interface parse_block_float
    subroutine oct_parse_block_double(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in)  :: blk
      integer,       intent(in)  :: l, c
      real(8),       intent(out) :: res
    end subroutine oct_parse_block_double

    module procedure oct_parse_block_double_unit
  end interface parse_block_float

  interface parse_block_cmplx
    subroutine oct_parse_block_complex(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in)  :: blk
      integer,       intent(in)  :: l, c
      complex(8),    intent(out) :: res
    end subroutine oct_parse_block_complex
  end interface parse_block_cmplx

  interface parse_block_string
    subroutine oct_parse_block_string(blk, l, c, res)
      import block_t
      implicit none
      type(block_t),    intent(in)  :: blk
      integer,          intent(in)  :: l, c
      character(len=*), intent(out) :: res
    end subroutine oct_parse_block_string
  end interface parse_block_string

  ! ---------------------------------------------------------
  !> The public subroutine parse_expression accepts two
  !! possible interfaces, one which assumes that the variables
  !! in the expression are "x(:)", "r" and "t", and another
  !! one which permits to set one variable to whichever string.
  !! Examples of usage:
  !!
  !! call parse_expression(f_re, f_im, ndim, x(:), r, t, &
  !!   "0.5*0.01*r^2")
  !!
  !! call parse_expression(f_re, f_im, "t", t, "cos(0.01*t)")
  ! ---------------------------------------------------------

  interface
    subroutine oct_parse_expression(re, im, ndim, x, r, t, pot)
      implicit none
      real(8),          intent(in)  :: x, r, t
      integer,          intent(in)  :: ndim
      real(8),          intent(out) :: re, im
      character(len=*), intent(in)  :: pot
    end subroutine oct_parse_expression
  end interface

  interface parse_expression
    subroutine oct_parse_expression1(re, im, c, x, string)
      implicit none
      real(8),          intent(out) :: re, im
      character(len=*), intent(in)  :: c
      real(8),          intent(in)  :: x
      character(len=*), intent(in)  :: string
    end subroutine oct_parse_expression1

    module procedure oct_parse_expression_vec
  end interface

contains

  ! ---------------------------------------------------------
  subroutine parser_init()

    integer :: ierr
    logical :: file_exists

    ! check files are present
    inquire(file=trim(conf%share)//'/variables', exist=file_exists)
    if(.not. file_exists) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Error initializing parser'
      write(stderr,'(a)') 'Cannot open variables file: '//trim(conf%share)//'/variables'
      call parse_fatal()
    end if

    inquire(file='inp', exist=file_exists)
    if(.not. file_exists) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Error initializing parser'
      write(stderr,'(a)') 'Cannot open input file!'
      write(stderr,'(a)') 'Please provide an input file with name inp in the current workdir'
      call parse_fatal()
    end if

    ! initialize the parser
    if(mpi_grp_is_root(mpi_world)) call loct_mkdir('exec')
    ierr = parse_init('exec/parser.log', mpi_world%rank)
    if(ierr /= 0) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Error initializing parser: cannot write to exec/parser.log.'
      write(stderr,'(a)') 'Do you have write permissions in this directory?'
      call parse_fatal()
    end if

    ! read in option definitions
    ierr = parse_input_file(trim(conf%share)//'/variables', set_used = 1)
    if(ierr /= 0) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Error initializing parser'
      write(stderr,'(a)') 'Cannot open variables file: '//trim(conf%share)//'/variables'
      call parse_fatal()
    end if

    ! setup standard input
    ierr = parse_input_file('inp', set_used = 0)
    if(ierr /= 0) then
      write(stderr,'(a)') '*** Fatal Error (description follows)'
      write(stderr,'(a)') 'Error initializing parser'
      write(stderr,'(a)') 'Cannot open input file!'
      write(stderr,'(a)') 'Please provide an input file with name inp in the current workdir'
      call parse_fatal()
    end if

    ! parse OCT_ prefixed variables from environment
    call parse_environment("OCT_")

  end subroutine parser_init


  ! ---------------------------------------------------------
  subroutine parser_end()

    call sym_output_table(only_unused = 1, mpiv_node = mpi_world%rank)
    call parse_end()

  end subroutine parser_end

  ! ---------------------------------------------------------

  logical function parse_is_defined(namespace, name) result(isdef)
    type(namespace_t), intent(in) :: namespace
    character(len=*),  intent(in) :: name

    isdef = parse_isdef(parse_get_full_name(namespace, name)) /= 0

  end function parse_is_defined

  ! ---------------------------------------------------------

  subroutine parse_integer(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: def
    integer,           intent(out)   :: res

    integer(8) :: res8

    call parse_check_varinfo(name)
    call oct_parse_int(parse_get_full_name(namespace, name), int(def, 8), res8)

    res = int(res8)

  end subroutine parse_integer

  ! ---------------------------------------------------------

  subroutine parse_integer8(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    integer(8),        intent(in)    :: def
    integer(8),        intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_int(parse_get_full_name(namespace, name), def, res)

  end subroutine parse_integer8

  ! ---------------------------------------------------------

  subroutine parse_integer48(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: def
    integer(8),        intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_int(parse_get_full_name(namespace, name), int(def, 8), res)

  end subroutine parse_integer48

  ! ---------------------------------------------------------

  subroutine parse_integer84(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    integer(8),        intent(in)    :: def
    integer,           intent(out)   :: res

    integer(8) :: res8

    call parse_check_varinfo(name)
    call oct_parse_int(parse_get_full_name(namespace, name), def, res8)

    res = int(res8)

  end subroutine parse_integer84

  ! ---------------------------------------------------------

  subroutine parse_string(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    character(len=*),  intent(in)    :: def
    character(len=*),  intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_string(parse_get_full_name(namespace, name), def, res)

  end subroutine parse_string

  ! ---------------------------------------------------------
  !> logical is a FORTRAN type, so we emulate the routine with integers
  subroutine parse_logical(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    logical,           intent(in)    :: def
    logical,           intent(out)   :: res

    integer(8) :: idef, ires

    call parse_check_varinfo(name)

    idef = 0
    if(def) idef = 1

    call oct_parse_int(parse_get_full_name(namespace, name), idef, ires)
    res = (ires /= 0)

  end subroutine parse_logical

  ! ---------------------------------------------------------

  subroutine parse_cmplx(namespace, name, def, res)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    complex(8),        intent(in)    :: def
    complex(8),        intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_complex(parse_get_full_name(namespace, name), def, res)

  end subroutine parse_cmplx

  ! ---------------------------------------------------------

  integer function parse_block(namespace, name, blk, check_varinfo_)
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: name
    type(block_t),     intent(out)   :: blk
    logical, optional, intent(in)    :: check_varinfo_

    logical check_varinfo

    check_varinfo = .true.
    if(present(check_varinfo_)) check_varinfo = check_varinfo_

    if(check_varinfo) then
      call parse_check_varinfo(name)
    end if
    parse_block = oct_parse_block(parse_get_full_name(namespace, name), blk)

  end function parse_block

  ! ---------------------------------------------------------

  subroutine parse_block_logical(blk, l, c, res)
    type(block_t), intent(in)  :: blk
    integer,       intent(in)  :: l, c
    logical,       intent(out) :: res

    integer :: ires

    call oct_parse_block_int(blk, l, c, ires)
    res = (ires /= 0)

  end subroutine parse_block_logical

  ! ---------------------------------------------------------

  subroutine oct_parse_double_unit(namespace, name, def, res, unit)
    type(namespace_t),      intent(in)  :: namespace
    character(len=*),       intent(in)  :: name
    real(8),                intent(in)  :: def
    real(8),                intent(out) :: res
    type(unit_t), optional, intent(in)  :: unit

    call parse_check_varinfo(name)

    if(present(unit)) then
      call oct_parse_double(parse_get_full_name(namespace, name), units_from_atomic(unit, def), res)
      res = units_to_atomic(unit, res)
    else
      call oct_parse_double(parse_get_full_name(namespace, name), def, res)
    end if

  end subroutine oct_parse_double_unit

  ! ---------------------------------------------------------

  subroutine oct_parse_block_double_unit(blk, l, c, res, unit)
    type(block_t), intent(in)  :: blk
    integer,       intent(in)  :: l, c
    real(8),       intent(out) :: res
    type(unit_t),  intent(in)  :: unit

    call oct_parse_block_double(blk, l, c, res)
    res = units_to_atomic(unit, res)

  end subroutine oct_parse_block_double_unit

  ! ---------------------------------------------------------
  subroutine oct_parse_expression_vec(re, im, ndim, x, r, t, pot)
    real(8),          intent(out) :: re, im
    integer,          intent(in)  :: ndim
    real(8),          intent(in)  :: x(:), r, t
    character(len=*), intent(in)  :: pot

    real(8) :: xc(1:MAX_DIM)

    xc = M_ZERO
    xc(1:ndim) = x(1:ndim)
    call oct_parse_expression(re, im, ndim, xc(1), r, t, pot)
  end subroutine oct_parse_expression_vec


  ! ----------------------------------------------------------------------
  !> A very primitive way to "preprocess" a string that contains reference
  !! to the elements of a two-dimensional array, substituting them with
  !! the values of the array x. This way the string can be processed by
  !! the parser later.
  subroutine parse_array(inp_string, x, arraychar)
    character(len=*), intent(inout) :: inp_string
    FLOAT,            intent(in)    :: x(:, :)
    character(len=1), intent(in)    :: arraychar

    integer              :: i, m, n_atom, coord, string_length
    character (LEN=100)  :: v_string

    string_length = len(inp_string)
    do i = 1, string_length - 1
       if(inp_string(i:i+1) == arraychar//"[") then
          m = 0
          if(inp_string(i+3:i+3) == ",") m = 1
          if(inp_string(i+4:i+4) == ",") m = 2
          if(m == 0) then
             write(stderr, '(a)') "*** Fatal Error (description follows)"
             write(stderr, '(a)') "Attempting to parse a string with array elements larger than 99"
             call parse_fatal()
          end if
          read(inp_string(i+2:i+1+m),*) n_atom
          read(inp_string(i+3+m:i+3+m),*) coord
          write(v_string,*) x(n_atom, coord)
          inp_string = inp_string(:i-1) // "(" // trim(v_string) // ")" // inp_string(i+5+m:)
       end if
    end do

  end subroutine parse_array

  ! ----------------------------------------------------------------------

  subroutine parse_check_varinfo(varname)
    character(len=*), intent(in) :: varname

    if(.not. varinfo_exists(varname)) then
      write(stderr,'(a)') "*** Fatal Internal Error (description follows)"
      write(stderr,'(a)') 'Attempting to parse undocumented variable '//trim(varname)//'.'
      call parse_fatal()
    end if

  end subroutine parse_check_varinfo

  ! ----------------------------------------------------------------------
  !> Given a namespace and a variable name, this function will iterate over all
  !! namespace ancestors contained in the namespace, until it finds one for
  !! which the variable is defined. If it finds such namespace it returns the
  !! variable name prefixed with the namespace. If it does not find any suitable
  !! namespace it returns the variable name without any prefix.
  !!
  !! To make it clear what we mean by all namespace ancestors contained in a
  !! given namesspace, lets suppose that we have the following namespace:
  !!
  !!   "A.B.C"
  !!
  !! Clearly "A" and "B" are ancestors of "C", but also the full path to "B",
  !! that is "A.B", is an ancestor. For practical purposes we will also consider
  !! that "C" is an ancestor of itself. So "C", "B.C" and "A.B.C" are also
  !! ancestors.
  !!
  !! The order in which the function iterates over the possible namespace
  !! ancestors is crucial, as it effectively determines namespace precedence in
  !! the input file. The order is such that it goes from the rigth-most ancestor
  !! to the left-most, and form the more complete path to the least complete. So
  !! for the above example, the order will be the following:
  !!
  !!   "A.B.C", "B.C", "C", "A.B", "B", "A"
  !!
  !! as "C" is the right-most ancestor while "A" is the left-most and "A.B.C" is
  !! the most complete path while "C" is the least complete.
  function parse_get_full_name(namespace, varname) result(name)
    type(namespace_t), target, intent(in)  :: namespace
    character(len=*),          intent(in)  :: varname
    character(len=:),          allocatable :: name

    logical :: found
    integer :: is
    type(namespace_t), pointer :: ancestor
    character(len=MAX_NAMESPACE_LEN) :: ancestor_name

    found = .false.

    ! Loop over all ancestors, starting from the right-most
    ancestor => namespace
    do while (associated(ancestor) .and. .not. found)

      ! Loop over all paths to this ancestor, starting from the most complete path
      ancestor_name = ancestor%get()
      is = -1
      do while (len_trim(ancestor_name) > 0 .and. is /= 0 .and. .not. found)
        ! Check if the current path is found in the input file
        name = trim(ancestor_name) // "." // trim(varname)
        found = parse_isdef(trim(name)) /= 0

        ! Remove the left-most namespace ("is" will be zero if there is only one namespace left)
        is = index(ancestor_name, ".")
        ancestor_name = ancestor_name(is+1:)
      end do
      ancestor => ancestor%parent
    end do

    ! If no suitable namespace found, just return the variable name
    if (.not. found) name = varname

  end function parse_get_full_name

  ! ----------------------------------------------------------------------
  subroutine parse_fatal()

#ifdef HAVE_MPI
    if(mpi_world%comm /= -1) call MPI_Abort(mpi_world%comm, 999, mpi_err)
#endif
    stop

  end subroutine parse_fatal

end module parser_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
