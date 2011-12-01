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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

! This module is only supposed to be used within this file.
module block_t_m
  implicit none
  
  private

  type, public :: block_t
    private
    integer, pointer :: p
  end type block_t

end module block_t_m


module parser_m
  use block_t_m
  use global_m
  use loct_m
  use mpi_m
  use unit_m

  implicit none

  ! Define which routines can be seen from the outside.
  private
  public ::              &
    block_t,             &   ! This is defined in block_t_m above
    parser_init,         &
    parser_end,          &
    parse_init,          &
    parse_putsym,        &
    parse_end,           &
    parse_isdef,         &
    parse_integer,       &
    parse_float,         &
    parse_cmplx,         &
    parse_string,        &
    parse_logical,       &
    parse_block,         &
    parse_block_end,     &
    parse_block_n,       &
    parse_block_cols,    &
    parse_block_integer, &
    parse_block_float,   &
    parse_block_cmplx,   &
    parse_block_string,  &
    parse_block_logical, &
    parse_expression

  interface parse_init
    integer function oct_parse_init(file_out, mpiv_node)
      character(len=*), intent(in)  :: file_out
      integer, intent(in) :: mpiv_node
    end function oct_parse_init
  end interface

  interface parse_putsym
    subroutine oct_parse_putsym_int(sym, i)
      character(len=*), intent(in)  :: sym
      integer, intent(in) :: i
    end subroutine oct_parse_putsym_int
    subroutine oct_parse_putsym_double(sym, d)
      character(len=*), intent(in)  :: sym
      real(8), intent(in) :: d
    end subroutine oct_parse_putsym_double
    module procedure oct_parse_putsym_double4
  end interface

  interface parse_input_file
    integer function oct_parse_input(file_in)
      character(len=*), intent(in)  :: file_in
    end function oct_parse_input
  end interface

  interface parse_end
    subroutine oct_parse_end()
    end subroutine oct_parse_end
  end interface

  interface parse_isdef
    integer function oct_parse_isdef(name)
      character(len=*), intent(in) :: name
    end function oct_parse_isdef
  end interface

  interface parse_integer
    subroutine oct_parse_int(name, def, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_int
  end interface

  interface parse_float
    subroutine oct_parse_double(name, def, res)
      character(len=*), intent(in)  :: name
      real(8),          intent(in)  :: def
      real(8),          intent(out) :: res
    end subroutine oct_parse_double
    module procedure oct_parse_double4_unit
    module procedure oct_parse_double8_unit
  end interface

  interface parse_cmplx
    subroutine oct_parse_complex(name, def, res)
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex
    module procedure oct_parse_complex4
  end interface

  interface parse_string
    subroutine oct_parse_string(name, def, res)
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string
  end interface

  interface parse_block
    integer function oct_parse_block(name, blk)
      use block_t_m
      character(len=*), intent(in) :: name
      type(block_t), intent(out) :: blk
    end function oct_parse_block
  end interface

  interface parse_block_end
    subroutine oct_parse_block_end(blk)
      use block_t_m
      type(block_t), intent(in) :: blk
    end subroutine oct_parse_block_end
  end interface

  interface parse_block_n
    integer function oct_parse_block_n(blk)
      use block_t_m
      type(block_t), intent(in) :: blk
    end function oct_parse_block_n
  end interface

  interface parse_block_cols
    integer function oct_parse_block_cols(blk, line)
      use block_t_m
      type(block_t), intent(in) :: blk
      integer, intent(in) :: line
    end function oct_parse_block_cols
  end interface

  interface parse_block_integer
    subroutine oct_parse_block_int(blk, l, c, res)
      use block_t_m
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int
  end interface

  interface parse_block_float
    subroutine oct_parse_block_double(blk, l, c, res)
      use block_t_m
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double
    module procedure oct_parse_block_double4
    module procedure oct_parse_block_double4_unit
    module procedure oct_parse_block_double8_unit
  end interface

  interface parse_block_cmplx
    subroutine oct_parse_block_complex(blk, l, c, res)
      use block_t_m
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex
    module procedure oct_parse_block_complex4
  end interface

  interface parse_block_string
    subroutine oct_parse_block_string(blk, l, c, res)
      use block_t_m
      type(block_t), intent(in) :: blk
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
    end subroutine oct_parse_block_string
  end interface

  ! ---------------------------------------------------------
  ! The public subroutine parse_expression accepts two
  ! possible interfaces, one which assumes that the variables
  ! in the expression are "x(:)", "r" and "t", and another
  ! one which permits to set one variable to whichever string.
  ! Examples of usage:
  !
  ! call parse_expression(f_re, f_im, ndim, x(:), r, t, &
  !   "0.5*0.01*r^2")
  !
  ! call parse_expression(f_re, f_im, "t", t, "cos(0.01*t)")
  ! ---------------------------------------------------------

  interface
    subroutine oct_parse_expression(re, im, ndim, x, r, t, pot)
      real(8),          intent(in)  :: x, r, t
      integer,          intent(in)  :: ndim
      real(8),          intent(out) :: re, im
      character(len=*), intent(in)  :: pot
    end subroutine oct_parse_expression
  end interface

  interface parse_expression
    subroutine oct_parse_expression1(re, im, c, x, string)
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
  subroutine parser_init
    integer :: ierr
    
    ! initialize the parser
    if(mpi_grp_is_root(mpi_world)) call loct_mkdir('exec')
    ierr = parse_init('exec/parser.log', mpi_world%rank)
    if(ierr .ne. 0) then
      write(0,'(a)') '*** Fatal Error (description follows)'
      write(0,'(a)') 'Error initializing liboct'
      write(0,'(a)') 'Do you have write permissions in this directory?'
#ifdef HAVE_MPI
      call MPI_Finalize(mpi_err)
#endif
      stop
    end if

    ! read in default variables
    ierr = parse_input_file(trim(conf%share)//'/variables')
    if(ierr .ne. 0) then
      write(0,'(a)') '*** Fatal Error (description follows)'
      write(0,'(a)') 'Cannot open variables file: '//trim(conf%share)//'/variables'
#ifdef HAVE_MPI
      call MPI_Finalize(mpi_err)
#endif
      stop
    end if

    ! setup standard input
    ierr = parse_input_file('inp')
    if(ierr .ne. 0) then 
      write(0,'(a)') '*** Fatal Error (description follows)' 
      write(0,'(a)') 'Error initializing liboct' 
      write(0,'(a)') 'Cannot open input file!' 
      write(0,'(a)') 'Please provide an input file with name inp in the current workdir'
#ifdef HAVE_MPI
      call MPI_Finalize(mpi_err)
#endif
      stop
    end if


  end subroutine parser_init


  ! ---------------------------------------------------------
  subroutine parser_end

    call parse_end()

  end subroutine parser_end


  ! ---------------------------------------------------------
  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine parse_logical(name, def, res)
    character(len=*), intent(in) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res

    integer :: idef, ires

    idef = 0
    if(def) idef = 1

    call oct_parse_int(name, idef, ires)
    res = (ires .ne. 0)

  end subroutine parse_logical


  ! ---------------------------------------------------------
  subroutine parse_block_logical(blk, l, c, res)
    type(block_t), intent(in) :: blk
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res

    integer :: ires

    call oct_parse_block_int(blk, l, c, ires)
    res = (ires .ne. 0)

  end subroutine parse_block_logical

  ! The code may want to compile in single-precision mode.
  ! As I did not want to change the parser library, these
  ! driver functions just convert their arguments.

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
    type(unit_t),     intent(in)  :: unit

    call oct_parse_double(name, units_from_atomic(unit, def), res)

    res = units_to_atomic(unit, res)
    
  end subroutine oct_parse_double8_unit

  ! ---------------------------------------------------------
  subroutine oct_parse_complex4(name, def4, res4)
    character(len=*), intent(in) :: name
    complex(4), intent(in)       :: def4
    complex(4), intent(out)      :: res4

    complex(8) :: res8
    call oct_parse_complex(name, cmplx(def4, kind=8), res8)
    res4 = real(res8, kind=4)
  end subroutine oct_parse_complex4


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

end module parser_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
