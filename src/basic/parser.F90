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
  use unit_oct_m
  use varinfo_oct_m
  
  implicit none

  private
  public ::              &
    parser_t,            &
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

  type :: parser_t
    private
    integer :: dummy
  end type parser_t
  
  type :: block_t
    private
    integer, pointer :: p
  end type block_t

  interface parse_init
    integer function oct_parse_init(file_out, mpiv_node)
      implicit none
      character(len=*), intent(in)  :: file_out
      integer, intent(in) :: mpiv_node
    end function oct_parse_init
  end interface parse_init

  interface parse_putsym
    subroutine oct_parse_putsym_int(sym, i)
      implicit none
      character(len=*), intent(in)  :: sym
      integer, intent(in) :: i
    end subroutine oct_parse_putsym_int
    subroutine oct_parse_putsym_double(sym, d)
      implicit none
      character(len=*), intent(in)  :: sym
      real(8), intent(in) :: d
    end subroutine oct_parse_putsym_double
    module procedure oct_parse_putsym_double4
  end interface parse_putsym

  interface parse_input_file
    integer function oct_parse_input(file_in, set_used)
      implicit none
      character(len=*), intent(in)  :: file_in
      integer,          intent(in)  :: set_used
    end function oct_parse_input
  end interface parse_input_file

  interface parse_environment
    subroutine oct_parse_environment(prefix)
      implicit none
      character(len=*), intent(in)  :: prefix
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
      character(len=*), intent(in) :: name
      integer(8), intent(in)       :: def
      integer(8), intent(out)      :: res
    end subroutine oct_parse_int

    subroutine oct_parse_double(name, def, res)
      implicit none
      character(len=*), intent(in)  :: name
      real(8),          intent(in)  :: def
      real(8),          intent(out) :: res
    end subroutine oct_parse_double

    subroutine oct_parse_complex(name, def, res)
      implicit none
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex
    
    subroutine oct_parse_string(name, def, res)
      implicit none
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string

    integer function oct_parse_block(name, blk)
      import block_t
      implicit none
      character(len=*), intent(in) :: name
      type(block_t), intent(out) :: blk
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
    module procedure oct_parse_double4_unit
    module procedure oct_parse_double8_unit
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
      integer, intent(in) :: line
    end function oct_parse_block_cols
  end interface parse_block_cols

  interface parse_block_integer
    subroutine oct_parse_block_int(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int
  end interface parse_block_integer

  interface parse_block_float
    subroutine oct_parse_block_double(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double
    module procedure oct_parse_block_double4
    module procedure oct_parse_block_double4_unit
    module procedure oct_parse_block_double8_unit
  end interface parse_block_float

  interface parse_block_cmplx
    subroutine oct_parse_block_complex(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex
    module procedure oct_parse_block_complex4
  end interface parse_block_cmplx

  interface parse_block_string
    subroutine oct_parse_block_string(blk, l, c, res)
      import block_t
      implicit none
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
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
    module procedure oct_parse_expression_vec4
    module procedure oct_parse_expression14
  end interface

contains

  ! ---------------------------------------------------------
  subroutine parser_init(self)
    type(parser_t), intent(out) :: self
    
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
  subroutine parser_end(self)
    type(parser_t), intent(inout) :: self
    
    call sym_output_table(only_unused = 1, mpiv_node = mpi_world%rank)
    call parse_end()

  end subroutine parser_end

  ! ---------------------------------------------------------  

  logical function parse_is_defined(self, name) result(isdef)
    type(parser_t), intent(in)   :: self
    character(len=*), intent(in) :: name

    isdef = parse_isdef(name) /= 0
    
  end function parse_is_defined

  ! ---------------------------------------------------------  
  
  subroutine parse_integer(name, def, res)
    character(len=*), intent(in)    :: name
    integer,          intent(in)    :: def
    integer,          intent(out)   :: res

    integer(8) :: res8
    
    call parse_check_varinfo(name)
    call oct_parse_int(name, int(def, 8), res8)

    res = int(res8)
    
  end subroutine parse_integer

  ! ---------------------------------------------------------

  subroutine parse_integer8(name, def, res)
    character(len=*), intent(in)    :: name
    integer(8),       intent(in)    :: def
    integer(8),       intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_int(name, def, res)
    
  end subroutine parse_integer8

  ! ---------------------------------------------------------  
  
  subroutine parse_integer48(name, def, res)
    character(len=*), intent(in)    :: name
    integer,          intent(in)    :: def
    integer(8),       intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_int(name, int(def, 8), res)
    
  end subroutine parse_integer48

  ! ---------------------------------------------------------  
  
  subroutine parse_integer84(name, def, res)
    character(len=*), intent(in)    :: name
    integer(8),       intent(in)    :: def
    integer,          intent(out)   :: res

    integer(8) :: res8
    
    call parse_check_varinfo(name)
    call oct_parse_int(name, def, res8)

    res = int(res8)
    
  end subroutine parse_integer84

  ! ---------------------------------------------------------
  
  subroutine parse_string(name, def, res)
    character(len=*), intent(in)    :: name
    character(len=*), intent(in)    :: def
    character(len=*), intent(out)   :: res
    
    call parse_check_varinfo(name)
    call oct_parse_string(name, def, res)
    
  end subroutine parse_string
  
  ! ---------------------------------------------------------
  !> logical is a FORTRAN type, so we emulate the routine with integers
  subroutine parse_logical(name, def, res)
    character(len=*), intent(in)    :: name
    logical,          intent(in)    :: def
    logical,          intent(out)   :: res

    integer(8) :: idef, ires

    call parse_check_varinfo(name)
    
    idef = 0
    if(def) idef = 1

    call oct_parse_int(name, idef, ires)
    res = (ires /= 0)

  end subroutine parse_logical

  ! ---------------------------------------------------------
  
  subroutine parse_cmplx(name, def, res)
    character(len=*), intent(in)    :: name
    complex(8),       intent(in)    :: def
    complex(8),       intent(out)   :: res

    call parse_check_varinfo(name)
    call oct_parse_complex(name, def, res)
    
  end subroutine parse_cmplx

  ! ---------------------------------------------------------
  
  integer function parse_block(self, name, blk, check_varinfo_)
    type(parser_t),    intent(in)    :: self
    character(len=*),  intent(in)    :: name
    type(block_t),     intent(out)   :: blk
    logical, optional, intent(in)    :: check_varinfo_

    logical check_varinfo

    check_varinfo = .true.
    if(present(check_varinfo_)) check_varinfo = check_varinfo_

    if(check_varinfo) call parse_check_varinfo(name)
    parse_block = oct_parse_block(name, blk)

  end function parse_block

  ! ---------------------------------------------------------

  subroutine parse_block_logical(blk, l, c, res)
    type(block_t), intent(in) :: blk
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res

    integer :: ires

    call oct_parse_block_int(blk, l, c, ires)
    res = (ires /= 0)

  end subroutine parse_block_logical

  !> The code may want to compile in single-precision mode.
  !! As I did not want to change the parser library, these
  !! driver functions just convert their arguments.

  ! ---------------------------------------------------------
  subroutine oct_parse_putsym_double4(sym, d4)
    character(len=*), intent(in) :: sym
    real(4), intent(in) :: d4

    call oct_parse_putsym_double(sym, real(d4, 8))
  end subroutine oct_parse_putsym_double4


  ! ---------------------------------------------------------
  subroutine oct_parse_double4(name, def4, res4)
    character(len=*), intent(in) :: name
    real(4), intent(in)          :: def4
    real(4), intent(out)         :: res4

    real(8) :: res8
    call oct_parse_double(name, real(def4, 8), res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_double4

  ! ---------------------------------------------------------

  subroutine oct_parse_double4_unit(name, def4, res4, unit)
    character(len=*),       intent(in)  :: name
    real(4),                intent(in)  :: def4
    real(4),                intent(out) :: res4
    type(unit_t), optional, intent(in)  :: unit

    real(8) :: res8

    call parse_check_varinfo(name)

    if(present(unit)) then
      call oct_parse_double(name, units_from_atomic(unit, real(def4, 8)), res8)
      res4 = real(units_to_atomic(unit, res8), kind=4)
    else
      call oct_parse_double(name, real(def4, 8), res8)
      res4 = real(res8, kind=4)
    end if
    
  end subroutine oct_parse_double4_unit

  ! ---------------------------------------------------------

  subroutine oct_parse_double8_unit(name, def, res, unit)
    character(len=*), intent(in)  :: name
    real(8),          intent(in)  :: def
    real(8),          intent(out) :: res
    type(unit_t), optional, intent(in)  :: unit

    call parse_check_varinfo(name)
    
    if(present(unit)) then
      call oct_parse_double(name, units_from_atomic(unit, def), res)
      res = units_to_atomic(unit, res)
    else
      call oct_parse_double(name, def, res)
    end if
    
  end subroutine oct_parse_double8_unit

  ! ---------------------------------------------------------
  subroutine oct_parse_block_double4(blk, l, c, res4)
    type(block_t), intent(in) :: blk
    integer, intent(in)          :: l, c
    real(4), intent(out)         :: res4

    real(8) :: res8
    call oct_parse_block_double(blk, l, c, res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_block_double4

  ! ---------------------------------------------------------

  subroutine oct_parse_block_double4_unit(blk, l, c, res4, unit)
    type(block_t), intent(in)  :: blk
    integer,       intent(in)  :: l, c
    real(4),       intent(out) :: res4
    type(unit_t),  intent(in)  :: unit

    real(8) :: res8
    call oct_parse_block_double(blk, l, c, res8)
    res4 = real(units_to_atomic(unit, res8), kind=4)
  end subroutine oct_parse_block_double4_unit

  ! ---------------------------------------------------------

  subroutine oct_parse_block_double8_unit(blk, l, c, res, unit)
    type(block_t), intent(in)  :: blk
    integer,       intent(in)  :: l, c
    real(8),       intent(out) :: res
    type(unit_t),  intent(in)  :: unit

    call oct_parse_block_double(blk, l, c, res)
    res = units_to_atomic(unit, res)

  end subroutine oct_parse_block_double8_unit

  ! ---------------------------------------------------------
  subroutine oct_parse_block_complex4(blk, l, c, res4)
    type(block_t), intent(in) :: blk
    integer, intent(in)          :: l, c
    complex(4), intent(out)      :: res4

    complex(8) :: res8
    call oct_parse_block_complex(blk, l, c, res8)
    res4 = cmplx(res8, kind=4)
  end subroutine oct_parse_block_complex4

  ! ---------------------------------------------------------
  subroutine oct_parse_expression_vec(re, im, ndim, x, r, t, pot)
    real(8), intent(out) :: re, im
    integer, intent(in)  :: ndim
    real(8), intent(in)  :: x(:), r, t
    character(len=*), intent(in) :: pot

    real(8) :: xc(1:MAX_DIM)

    xc = M_ZERO
    xc(1:ndim) = x(1:ndim)
    call oct_parse_expression(re, im, ndim, xc(1), r, t, pot)
  end subroutine oct_parse_expression_vec

  ! ---------------------------------------------------------
  subroutine oct_parse_expression_vec4(re, im, ndim, x, r, t, pot)
    real(4), intent(out) :: re, im
    integer, intent(in)  :: ndim
    real(4), intent(in)  :: x(:), r, t
    character(len=*), intent(in) :: pot

    real(8) :: xc(1:MAX_DIM)
    real(8) :: re8, im8

    xc = M_ZERO
    xc(1:ndim) = real(x(1:ndim), 8)
    call oct_parse_expression(re8, im8, ndim, xc(1), real(r, 8), real(t, 8), pot)
    re = real(re8, 4)
    im = real(im8, 4)
  end subroutine oct_parse_expression_vec4

  ! ---------------------------------------------------------
  subroutine oct_parse_expression14(re, im, c, x, string)
    real(4), intent(out) :: re, im
    character(len=*), intent(in) :: c
    real(4), intent(in) :: x
    character(len=*), intent(in) :: string
    real(8) :: re8, im8
    call oct_parse_expression1(re8, im8, c, real(x, 8), string)
    re = real(re8, 4)
    im = real(im8, 4)
  end subroutine oct_parse_expression14


  ! ----------------------------------------------------------------------
  !> A very primitive way to "preprocess" a string that contains reference
  !! to the elements of a two-dimensional array, substituting them with
  !! the values of the array x. This way the string can be processed by
  !! the parser later.
  subroutine parse_array(inp_string, x, arraychar)
    character(len=*), intent(inout)  :: inp_string
    FLOAT, intent(in) :: x(:, :)
    character(len=1), intent(in) :: arraychar
    integer              :: i,m,n_atom,coord,string_length
    character (LEN=100)  :: v_string

    string_length = len(inp_string)
    do i=1, string_length - 1
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
