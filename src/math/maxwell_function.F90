!! Copyright (C) 2019 R. Jestaedt, H. Appel, F. Bonafe, M. Oliveira, N. Tancogne-Dejean
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

!>--------------------------------------------------------------
!! This module defines "Maxwell functions", to be used by
!! the Maxwell module especially for Maxwell incident waves
!!--------------------------------------------------------------
module maxwell_function_oct_m
  use iso_c_binding
  use fft_oct_m
  use global_oct_m
  use io_oct_m
  use loct_math_oct_m
  use math_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use splines_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use namespace_oct_m

  implicit none

  public :: mxf_t,                       &
            mxf_init,                    &
            mxf_init_const_wave,         &
            mxf_init_const_phase,        &
            mxf_init_gaussian_wave,      &
            mxf_init_cosinoidal_wave,    &
            mxf_init_fromexpr,           &
            mxf,                         &
            mxf_read,                    &
            mxf_is_empty

  integer, public, parameter ::  &
    MXF_EMPTY            =  10001,  &
    MXF_CONST_WAVE       =  10002,  &
    MXF_CONST_PHASE      =  10004,  &
    MXF_GAUSSIAN_WAVE    =  10005,  &
    MXF_COSINOIDAL_WAVE  =  10006,  &
    MXF_LOGISTIC_WAVE    =  10007,  &
    MXF_TRAPEZOIDAL_WAVE =  10008,  &
    MXF_FROM_FILE        =  10009,  &
    MXF_NUMERICAL        =  10010,  &
    MXF_FROM_EXPR        =  10011,  &
    MXF_FOURIER_SERIES   =  10012,  &
    MXF_ZERO_FOURIER     =  10013

  type mxf_t
    integer :: mode              = MXF_EMPTY
    FLOAT   :: k_vector(MAX_DIM) = M_ZERO
    FLOAT   :: r0(MAX_DIM)       = M_ZERO  !< vector at the maximum of the pulse
    FLOAT   :: width             = M_ZERO  !< the width of the pulse
    FLOAT   :: a0                = M_ZERO
    FLOAT   :: dx                = M_ZERO !< the space-discretization value.
    FLOAT   :: init_x            = M_ZERO
    FLOAT   :: final_x           = M_ZERO
    FLOAT   :: gr                = M_ZERO
    integer :: niter             = 0
    integer :: nfreqs            = 0

    type(spline_t)         :: amplitude
    character(len=200)     :: expression
    FLOAT, pointer :: val(:)    => NULL()
    FLOAT, pointer :: valww(:)  => NULL()
    type(fft_t) :: fft_handler
  end type mxf_t

  interface mxf
    module procedure mxft
  end interface mxf

contains

  !------------------------------------------------------------
  !> This function initializes "f" from the MXFunctions block.
  subroutine mxf_read(f, namespace, function_name, ierr)
    type(mxf_t),       intent(inout) :: f
    type(namespace_t), intent(in)    :: namespace
    character(len=*),  intent(in)    :: function_name
    integer,           intent(out)   :: ierr  !< Error code, 0 on success.

    type(block_t) :: blk
    integer :: nrows, ncols, i, function_type, idim
    character(len=100) :: row_name, function_expression
    FLOAT :: a0, r0(MAX_DIM), gr, width, k_vector(MAX_DIM)

    PUSH_SUB(mxf_read)

    !%Variable MaxwellFunctions
    !%Type block
    !%Section Time-Dependent
    !%Description
    !% This block specifies the shape of a "spatial-dependent function", such as the
    !% envelope needed when using the <tt>MaxwellFunctions</tt> block. Each line in the block
    !% specifies one function. The first element of each line will be a string
    !% that defines the name of the function. The second element specifies which type
    !% of function we are using; in the following we provide an example for each of the
    !% possible types:
    !%
    !%Option mxf_const_wave 10002
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_const_wave | kx | ky | kz | x0 | y0 | z0
    !% <br>%</tt>
    !%
    !% The function is constant plane wave <math> f(x,y,z) = a0 * cos( kx*(x-x0) + ky*(y-y0) + kz*(z-z0) ) </math>
    !%
    !%Option mxf_const_phase 10004
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_const_phase | kx | ky | kz | x0 | y0 | z0  
    !% <br>%</tt>
    !%
    !% The function is a constant phase of <math> f(x,y,z) = a0 * (kx*x0 + ky*y0 + kz*z0)
    !%
    !%Option mxf_gaussian_wave 10005
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_gaussian_wave | kx | ky | kz | x0 | y0 | z0 | width
    !% <br>%</tt>
    !% 
    !% The function is a Gaussian, <math> f(x,y,z) = a0 * \exp( -( kx*(x-x0) + ky*(y-y0) + kz*(z-z0) )^2 / (2 width^2) ) </math>
    !%
    !%Option mxf_cosinoidal_wave 10006
    !% 
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_cosinoidal_wave | kx | ky | kz | x0 | y0 | z0 | width
    !% <br>%</tt>
    !%
    !% <math> f(t) =  \cos( \frac{\pi}{2} \frac{kx*(x-x0)+ky*(y-y0)+kz*(z-z0)-2\width}{width} + pi )  </math>
    !%
    !% If <math> | k_x*x+k_y*y+k_z*z - x_0 | > \xi_0 </math>, then <math> f(x,y,z) = 0 </math>.
    !%
    !%Option mxf_logistic_wave 10007
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_logistic_wave | kx | ky | kz | x0 | y0 | z0 | gr | width 
    !% <br>%</tt>
    !% 
    !% The function is a logistic function, <math> f(x,y,z) = a0 * 1/(1+exp(gr*(kx*(x-x0)+ky*(y-y0)+kz*(kz*(z-z0))+width/2))) * 1/(1+exp(-gr*(kx*(x-x0)+ky*(y-y0)+kz*(kz*(z-z0))-width/2)))  </math>
    !%
    !%Option mxf_trapezoidal_wave 10008
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_trapezoidal_wave | kx | ky | kz | x0 | y0 | z0 | gr | width 
    !% <br>%</tt>
    !% 
    !% The function is a logistic function, <math> f(x,y,z) = a0 * ( ( 1-gr*(k*(r-r0)-width/2)*Theta(k*(r-r0)-width/2) ) * Theta(-(k*(r-r0)+width/2+1/gr))
    !%                                                             + (-1+gr*(k*(r-r0)+width/2)*Theta(k*(r-r0)+width/2) ) * Theta(-(k*(r-r0)-width/2+1/gr)) </math>
    !%
    !%Option mxf_from_expr 10011
    !%
    !% <tt>%MaxwellFunctions
    !% <br>&nbsp;&nbsp; "function-name" | mxf_from_expr | "expression"
    !% <br>%</tt>
    !%
    !% The temporal shape of the field is given as an expression (e.g., <tt>cos(2.0*x-3*y+4*z)</tt>. The 
    !% letter <i>x</i>, <i>y</i>, <i>z</i> means spatial coordinates, obviously. 
    !% The expression is used to construct the function <i>f</i>
    !% that defines the field.
    !%End
    ierr = -3
    if(parse_block(namespace, 'MaxwellFunctions', blk) /= 0) then
      ierr = -1
      POP_SUB(mxf_read)
      return
    end if

    nrows = parse_block_n(blk)
    row_loop: do i = 1, nrows
      call parse_block_string(blk, i-1, 0, row_name)
      if(trim(row_name) == trim(function_name)) then

        ncols = parse_block_cols(blk, i-1)
        call parse_block_integer(blk, i-1, 1, function_type)

        a0 = M_ONE; r0 = M_ZERO; width = M_ZERO; k_vector = M_ZERO
        select case(function_type)
        case (MXF_CONST_WAVE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call mxf_init_const_wave(f, a0, k_vector, r0)
        case (MXF_CONST_PHASE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call mxf_init_const_phase(f, a0, k_vector, r0)
        case (MXF_GAUSSIAN_WAVE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call parse_block_float(blk, i-1, 8, width, units_inp%length)
          call mxf_init_gaussian_wave(f, a0, k_vector, r0, width)
        case (MXF_COSINOIDAL_WAVE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call parse_block_float(blk, i-1, 8, width, units_inp%length)
          call mxf_init_cosinoidal_wave(f, a0, k_vector, r0, width)
        case (MXF_LOGISTIC_WAVE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call parse_block_float(blk, i-1, 8, gr, units_inp%length)
          call parse_block_float(blk, i-1, 9, width, units_inp%length)
          call mxf_init_logistic_wave(f, a0, k_vector, r0, gr, width)
        case (MXF_TRAPEZOIDAL_WAVE)
          do idim = 1, 3
            call parse_block_float(blk, i-1, 1+idim, k_vector(idim), unit_one/units_inp%length)
          end do
          do idim = 1, 3
            call parse_block_float(blk, i-1, 4+idim, r0(idim), units_inp%length)
          end do
          call parse_block_float(blk, i-1, 8, gr, units_inp%length)
          call parse_block_float(blk, i-1, 9, width, units_inp%length)
          call mxf_init_trapezoidal_wave(f, a0, k_vector, r0, gr, width)
        case (MXF_FROM_EXPR)
          call parse_block_string(blk, i-1, 2, function_expression)
          call mxf_init_fromexpr(f, trim(function_expression))
        case default
          ierr = -2
          call parse_block_end(blk)
          POP_SUB(mxf_read)
          return
        end select

        ierr = 0
        exit row_loop
      end if
    end do row_loop

    call parse_block_end(blk)

    POP_SUB(mxf_read)
  end subroutine mxf_read
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init(f)
    type(mxf_t), intent(inout) :: f

    PUSH_SUB(mxf_init)

    f%mode = MXF_EMPTY
    f%niter = 0
    f%dx = M_ZERO
    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init)
  end subroutine mxf_init
  !------------------------------------------------------------


  !------------------------------------------------------------
  logical function mxf_is_empty(f)
    type(mxf_t), intent(in) :: f

    PUSH_SUB(mxf_is_empty)
    mxf_is_empty = (f%mode == MXF_EMPTY)

    POP_SUB(mxf_is_empty)
  end function mxf_is_empty
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_const_wave(f, a0, k_vector, r0)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM)

    PUSH_SUB(mxf_init_const_wave)

    f%mode = MXF_CONST_WAVE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_const_wave)
  end subroutine mxf_init_const_wave
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_const_phase(f, a0, k_vector, r0)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM)

    PUSH_SUB(mxf_init_const_phase)

    f%mode = MXF_CONST_PHASE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_const_phase)
  end subroutine mxf_init_const_phase
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_gaussian_wave(f, a0, k_vector, r0, width)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM), width

    PUSH_SUB(mxf_init_gaussian_wave)

    f%mode = MXF_GAUSSIAN_WAVE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    f%width = width
    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_gaussian_wave)
  end subroutine mxf_init_gaussian_wave
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_cosinoidal_wave(f, a0, k_vector, r0, width)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM), width

    PUSH_SUB(mxf_init_cosinoidal_wave)

    f%mode = MXF_COSINOIDAL_WAVE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    f%width = width
    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_cosinoidal_wave)
  end subroutine mxf_init_cosinoidal_wave
  !------------------------------------------------------------
    

  !------------------------------------------------------------
  subroutine mxf_init_logistic_wave(f, a0, k_vector, r0, gr, width)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM), gr, width

    PUSH_SUB(mxf_init_logistic_wave)

    f%mode = MXF_LOGISTIC_WAVE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    f%gr = gr
    f%width = width

    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_logistic_wave)
  end subroutine
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_trapezoidal_wave(f, a0, k_vector, r0, gr, width)
    type(mxf_t), intent(inout) :: f
    FLOAT,       intent(in)    :: a0, k_vector(MAX_DIM), r0(MAX_DIM), gr, width

    PUSH_SUB(mxf_init_trapezoidal_wave)

    f%mode = MXF_TRAPEZOIDAL_WAVE
    f%a0 = a0
    f%k_vector = k_vector
    f%r0 = r0
    f%gr = gr
    f%width = width

    nullify(f%val)
    nullify(f%valww)

    POP_SUB(mxf_init_trapezoidal_wave)
  end subroutine
  !------------------------------------------------------------


  !------------------------------------------------------------
  subroutine mxf_init_fromexpr(f, expression)
    type(mxf_t),      intent(inout) :: f
    character(len=*), intent(in)    :: expression

    PUSH_SUB(mxf_init_fromexpr)

    f%mode       = MXF_FROM_EXPR
    f%expression = trim(expression)
    nullify(f%val)
    nullify(f%valww)
    
    POP_SUB(mxf_init_fromexpr)
  end subroutine mxf_init_fromexpr
  !------------------------------------------------------------


  !------------------------------------------------------------
  FLOAT function mxft(f, x) result(y)
    type(mxf_t), intent(in) :: f
    FLOAT,       intent(in) :: x(:)

    FLOAT :: r, xx, limit_1, limit_2, limit_3, limit_4
    integer :: xdim

    ! no push_sub because it is called too frequently

    xdim = size(x)

    select case(f%mode)
    case (MXF_CONST_WAVE)

      y = f%a0 * cos( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) )

    case (MXF_CONST_PHASE)

      y = f%a0 * sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim)))

    case (MXF_GAUSSIAN_WAVE)

       r = exp( - ( ( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) &
            / sqrt(sum(f%k_vector(1:xdim)**2)) )**2 / (M_TWO*f%width**2) ) )
      y = f%a0 * r * cos( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) )

    case (MXF_COSINOIDAL_WAVE)

      r = M_ZERO
      if(abs( sum( f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))/sqrt(sum(f%k_vector(1:xdim)**2)) ) ) <= f%width) then
         r = - cos((M_Pi/M_TWO) * ((sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) &
              / sqrt(sum(f%k_vector(1:xdim)**2)) - M_TWO*f%width) / f%width))
        r = r * cos( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) )
      end if
      y = f%a0 * r

    case (MXF_LOGISTIC_WAVE)

       r = M_ONE/(M_ONE + exp(f%gr*(sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) &
            / sqrt(sum(f%k_vector(1:xdim)**2)) - f%width/M_TWO))) &
            + M_ONE/(M_ONE + exp(-f%gr*(sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) &
            / sqrt(sum(f%k_vector(1:xdim)**2)) + f%width/M_TWO))) - M_ONE
      y = f%a0 * r * cos( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) )

    case (MXF_TRAPEZOIDAL_WAVE)

      xx = sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim)))/sqrt(sum(f%k_vector(1:xdim)**2))
      limit_1 = - f%width/M_TWO - M_ONE/f%gr
      limit_2 = - f%width/M_TWO
      limit_3 =   f%width/M_TWO
      limit_4 =   f%width/M_TWO + M_ONE/f%gr
      if ( ( xx > limit_1 ) .and. ( xx <= limit_2 ) ) then
        r = M_ONE + f%gr * (sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim)))/sqrt(sum(f%k_vector(1:xdim)**2)) + f%width/M_TWO)
      else if ( ( xx > limit_2 ) .and. ( xx <= limit_3 ) ) then
        r = M_ONE
      else if ( ( xx > limit_3 ) .and. ( xx <= limit_4 ) ) then
        r = M_ONE - f%gr * (sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim)))/sqrt(sum(f%k_vector(1:xdim)**2)) - f%width/M_TWO)
      else
        r = M_ZERO
      end if
      y = f%a0 * r * cos( sum(f%k_vector(1:xdim)*(x(:) - f%r0(1:xdim))) )

    case default

      y = M_ZERO

    end select

  end function mxft
  !------------------------------------------------------------


end module maxwell_function_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
