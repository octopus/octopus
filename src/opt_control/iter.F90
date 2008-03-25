!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: opt_control.F90 2875 2007-04-30 16:54:15Z acastro $

#include "global.h"

module opt_control_iter_m
  use datasets_m
  use varinfo_m
  use global_m
  use loct_parser_m
  use opt_control_constants_m
  use opt_control_parameters_m
  use grid_m
  use states_m
  use h_sys_output_m
  use messages_m

  implicit none

  private
  public :: oct_iterator_t,           &
            oct_iterator_init,        &
            oct_iterator_end,         &
            iteration_manager,        &
            iteration_manager_direct, &
            iterator_write

  type oct_iterator_t
    FLOAT              :: eps
    integer            :: ctr_iter_max
    integer            :: ctr_iter
    FLOAT, pointer     :: convergence(:,:)
    FLOAT              :: bestJ1, bestJ1_fluence, bestJ1_J
    integer            :: bestJ1_ctr_iter
    logical            :: maximize
    type(oct_control_parameters_t) :: best_par
  end type oct_iterator_t

contains

  ! ---------------------------------------------------------
  subroutine oct_iterator_init(iterator, par, maximize)
    type(oct_iterator_t), intent(inout)        :: iterator
    type(oct_control_parameters_t), intent(in) :: par
    logical, intent(in)                        :: maximize

    call push_sub('iter.oct_iter_init')

    !%Variable OCTEps
    !%Type float
    !%Section Optimal Control
    !%Default 1.0e-6
    !%Description
    !% Define the convergence threshold. It computes the difference between the "input"
    !% field in the iterative procedure, and the "output" field. If this difference is
    !% less then OCTEps the iteration is stopped. This difference is defined as:
    !% 
    !% <math>
    !% D[\epsilon^{i},\epsilon^{o}] = \int_0^T dt \vert \epsilon^{i}(t)-\epsilon^{o}(t)\vert^2\,.
    !% </math>
    !%
    !% (If there are several control fields, this difference is defined as the sum over
    !% all the individual differences).
    !%
    !% Whenever this condition is satisfied, it means that we have reached a solution point
    !% point of the QOCT equations, i.e. a critical point of the QOCT functional (not
    !% necessarily a maximum, and not necessarily the global maximum). 
    !%End
    call loct_parse_float(check_inp('OCTEps'), CNST(1.0e-6), iterator%eps)
    if(iterator%eps < M_ZERO) iterator%eps = tiny(CNST(1.0))

    !%Variable OCTMaxIter
    !%Type integer
    !%Section Optimal Control
    !%Default 10
    !%Description
    !% OCTMaxIter defines the maximum number of iterations.
    !% Typical values range from 10-100.
    !%End
    call loct_parse_int(check_inp('OCTMaxIter'), 10, iterator%ctr_iter_max)

    if( iterator%ctr_iter_max < 0 .and. iterator%eps < M_ZERO ) then
      message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
      call write_fatal(1)
    end if
    if(iterator%ctr_iter_max < 0) iterator%ctr_iter_max = huge(iterator%ctr_iter_max)

    iterator%ctr_iter = 0

    iterator%maximize = maximize

    ALLOCATE(iterator%convergence(5, 0:iterator%ctr_iter_max), (iterator%ctr_iter_max+1)*5)
    iterator%convergence = M_ZERO

    if(maximize) then
      iterator%bestJ1        = -CNST(1.0e20)
    else
      iterator%bestj1        = CNST(1.0e20)
    end if
    iterator%bestJ1_fluence  = M_ZERO
    iterator%bestJ1_J        = M_ZERO
    iterator%bestJ1_ctr_iter = 0

    call parameters_copy(iterator%best_par, par)

    call pop_sub()
  end subroutine oct_iterator_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_iterator_end(iterator)
    type(oct_iterator_t), intent(inout) :: iterator

    call push_sub('iter.oct_iterator_end')

    deallocate(iterator%convergence)
    nullify(iterator%convergence)
    call parameters_end(iterator%best_par)

    call pop_sub()
  end subroutine oct_iterator_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical function iteration_manager(j1, par, par_prev, iterator) result(stoploop)
    FLOAT, intent(in) :: j1
    type(oct_control_parameters_t), intent(in)  :: par
    type(oct_control_parameters_t), intent(in)  :: par_prev
    type(oct_iterator_t), intent(inout) :: iterator

    FLOAT :: fluence, jfunctional, j2
    logical :: bestj1

    call push_sub('iter.iteration_manager')
    
    stoploop = .false.

    fluence = parameters_fluence(par)
    j2 = parameters_j2(par)
    jfunctional = j1 + j2

    iterator%convergence(1, iterator%ctr_iter) = jfunctional
    iterator%convergence(2, iterator%ctr_iter) = j1
    iterator%convergence(3, iterator%ctr_iter) = j2
    ! WARNING: this does not consider the possibility of different 
    ! alphas for different control parameters.
    iterator%convergence(4, iterator%ctr_iter) = par%alpha(1)
    iterator%convergence(5, iterator%ctr_iter) = parameters_diff(par, par_prev)
    
    if(iterator%ctr_iter .eq. iterator%ctr_iter_max) then
      message(1) = "Info: Maximum number of iterations reached"
      call write_info(1)
      stoploop = .true.
    end if

    if( (iterator%eps > M_ZERO) .and. &
        (iterator%convergence(5, iterator%ctr_iter) < iterator%eps) .and. &
        (iterator%ctr_iter > 0 ) ) then
      message(1) = "Info: Convergence threshold reached"
      call write_info(1)
      stoploop = .true.
    end if

    write(message(1), '(a,i5)') 'Optimal control iteration #', iterator%ctr_iter
    call messages_print_stress(stdout, trim(message(1)))

    write(message(1), '(6x,a,f12.5)')  " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')  " => J        = ", jfunctional
    write(message(3), '(6x,a,f12.5)')  " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')  " => Fluence  = ", fluence
    write(message(5), '(6x,a,f12.5)')  " => Penalty  = ", par%alpha(1)
    write(message(6), '(6x,a,es12.2)') " => D[e,e']  = ", iterator%convergence(5, iterator%ctr_iter)
    if(iterator%ctr_iter .ne. 0) then
      call write_info(6)
    else
      call write_info(5)
    end if
    call messages_print_stress(stdout)


    if(iterator%maximize) then
      bestj1 = (j1 > iterator%bestj1)
    else
      bestj1 = (j1 < iterator%bestj1)
    end if

    ! store field with best J1
    if(bestj1) then
      iterator%bestJ1          = j1
      iterator%bestJ1_J        = jfunctional
      iterator%bestJ1_fluence  = fluence       
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      call parameters_end(iterator%best_par)
      call parameters_copy(iterator%best_par, par)
    end if

    iterator%ctr_iter = iterator%ctr_iter + 1
    
    call pop_sub()
  end function iteration_manager
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iteration_manager_direct(j1, par, iterator, dx)
    FLOAT, intent(in) :: j1
    type(oct_control_parameters_t), intent(in)  :: par
    type(oct_iterator_t), intent(inout) :: iterator
    FLOAT, optional, intent(in) :: dx


    FLOAT :: j, j2, fluence
    call push_sub('iter.iteration_manager_direct')

    fluence = parameters_fluence(par)
    j2 = - par%alpha(1) * (fluence - par%targetfluence)
    j  = j1 + j2

    iterator%convergence(1, iterator%ctr_iter) = j
    iterator%convergence(2, iterator%ctr_iter) = j1
    iterator%convergence(3, iterator%ctr_iter) = j2
    ! WARNING: this does not consider the possibility of different 
    ! alphas for different control parameters.
    iterator%convergence(4, iterator%ctr_iter) = par%alpha(1)
    if(present(dx)) then
      iterator%convergence(5, iterator%ctr_iter) = dx
    else
      iterator%convergence(5, iterator%ctr_iter) = M_ZERO
    end if

    write(message(1), '(a,i5)') 'Optimal control iteration #', iterator%ctr_iter
    call messages_print_stress(stdout, trim(message(1)))

    write(message(1), '(6x,a,f12.5)')    " => J1       = ", j1
    write(message(2), '(6x,a,f12.5)')    " => J        = ", j
    write(message(3), '(6x,a,f12.5)')    " => J2       = ", j2
    write(message(4), '(6x,a,f12.5)')    " => Fluence  = ", fluence
    call write_info(4)
    if(present(dx)) then
      write(message(1), '(6x,a,f12.5)')  " => Delta    = ", dx
      call write_info(1)
    end if
    call messages_print_stress(stdout)

    ! store field with best J1
    if(j1 > iterator%bestJ1) then
      iterator%bestJ1          = j1
      iterator%bestJ1_J        = j
      iterator%bestJ1_fluence  = fluence       
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      call parameters_end(iterator%best_par)
      call parameters_copy(iterator%best_par, par)
    end if

    call pop_sub()
  end subroutine iteration_manager_direct
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iterator_write(iterator, par)
    type(oct_iterator_t),           intent(in) :: iterator
    type(oct_control_parameters_t), intent(in) :: par

    character(len=80)  :: filename
    call push_sub('iter.iterator_write')

    write(filename,'(a,i3.3)') 'opt-control/laser.', iterator%ctr_iter
    call parameters_write(filename, par)

    call pop_sub()
  end subroutine iterator_write

end module opt_control_iter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
