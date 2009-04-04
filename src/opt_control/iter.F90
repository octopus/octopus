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
  use io_m
  use global_m
  use loct_parser_m
  use opt_control_parameters_m
  use grid_m
  use states_m
  use h_sys_output_m
  use messages_m
  use profiling_m

  implicit none

  private
  public :: oct_iterator_t,           &
            oct_iterator_init,        &
            oct_iterator_end,         &
            iteration_manager,        &
            iteration_manager_direct, &
            iterator_write,           &
            oct_iterator_bestpar,     &
            oct_iterator_current,     &
            oct_iterator_maxiter,     &
            oct_iterator_tolerance


  type oct_iterator_t
    private
    FLOAT              :: eps
    integer            :: ctr_iter_max
    integer            :: ctr_iter
    FLOAT              :: bestJ1
    FLOAT              :: bestJ1_fluence
    FLOAT              :: bestJ1_J
    integer            :: bestJ1_ctr_iter
    type(oct_control_parameters_t), pointer :: best_par
    integer            :: convergence_iunit
  end type oct_iterator_t

contains


  ! ---------------------------------------------------------
  subroutine oct_iterator_init(iterator, par)
    type(oct_iterator_t), intent(inout)        :: iterator
    type(oct_control_parameters_t), intent(in) :: par

    call push_sub('iter.oct_iter_init')

    !%Variable OCTEps
    !%Type float
    !%Section Calculation Modes::Optimal Control
    !%Default 1.0e-6
    !%Description
    !% Define the convergence threshold. It computes the difference between the "input"
    !% field in the iterative procedure, and the "output" field. If this difference is
    !% less than OCTEps the iteration is stopped. This difference is defined as:
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
    call loct_parse_float(datasets_check('OCTEps'), CNST(1.0e-6), iterator%eps)
    if(iterator%eps < M_ZERO) iterator%eps = tiny(CNST(1.0))

    !%Variable OCTMaxIter
    !%Type integer
    !%Section Calculation Modes::Optimal Control
    !%Default 10
    !%Description
    !% OCTMaxIter defines the maximum number of iterations.
    !% Typical values range from 10-100.
    !%End
    call loct_parse_int(datasets_check('OCTMaxIter'), 10, iterator%ctr_iter_max)

    if( iterator%ctr_iter_max < 0 .and. iterator%eps < M_ZERO ) then
      message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
      call write_fatal(1)
    end if
    if(iterator%ctr_iter_max < 0) iterator%ctr_iter_max = huge(iterator%ctr_iter_max)

    iterator%ctr_iter = 0
    iterator%bestJ1        = -CNST(1.0e20)
    iterator%bestJ1_fluence  = M_ZERO
    iterator%bestJ1_J        = M_ZERO
    iterator%bestJ1_ctr_iter = 0

    ALLOCATE(iterator%best_par, 1)
    call parameters_copy(iterator%best_par, par)

    iterator%convergence_iunit = io_open('opt-control/convergence', action='write')
    write(iterator%convergence_iunit, '(5(a))') '# iteration', '  J[Psi,chi,epsilon]', &
                                                '            J_1[Psi]', &
                                                '        J_2[epsilon]', &
                                                '               Delta'

    write(iterator%convergence_iunit, '(a)') &
     '###########################################################################################'


    call pop_sub()
  end subroutine oct_iterator_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_iterator_end(iterator)
    type(oct_iterator_t), intent(inout) :: iterator

    call push_sub('iter.oct_iterator_end')

    call parameters_write('opt-control/laser.bestJ1', iterator%best_par)

    call parameters_end(iterator%best_par)
    deallocate(iterator%best_par)
    nullify(iterator%best_par)
    call io_close(iterator%convergence_iunit)

    call pop_sub()
  end subroutine oct_iterator_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  logical function iteration_manager(j1, par, par_prev, iterator) result(stoploop)
    FLOAT, intent(in) :: j1
    type(oct_control_parameters_t), intent(in)  :: par
    type(oct_control_parameters_t), intent(in)  :: par_prev
    type(oct_iterator_t), intent(inout) :: iterator

    FLOAT :: fluence, jfunctional, j2, delta
    logical :: bestj1

    call push_sub('iter.iteration_manager')
    
    stoploop = .false.

    fluence = parameters_fluence(par)
    j2 = parameters_j2(par)
    jfunctional = j1 + j2
    delta = parameters_diff(par, par_prev)
    
    if(iterator%ctr_iter .eq. iterator%ctr_iter_max) then
      message(1) = "Info: Maximum number of iterations reached"
      call write_info(1)
      stoploop = .true.
    end if

    if( (iterator%eps > M_ZERO) .and. &
        (delta < iterator%eps) .and. &
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
    write(message(5), '(6x,a,f12.5)')  " => Penalty  = ", parameters_alpha(par, 1)
    write(message(6), '(6x,a,es12.2)') " => D[e,e']  = ", delta
    if(iterator%ctr_iter .ne. 0) then
      call write_info(6)
    else
      call write_info(5)
    end if
    call messages_print_stress(stdout)

    bestj1 = (j1 > iterator%bestj1)
    ! store field with best J1
    if(bestj1) then
      iterator%bestJ1          = j1
      iterator%bestJ1_J        = jfunctional
      iterator%bestJ1_fluence  = fluence       
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      call parameters_end(iterator%best_par)
      call parameters_copy(iterator%best_par, par)
    end if

    write(iterator%convergence_iunit, '(i11,4f20.8)')                &
      iterator%ctr_iter, jfunctional, j1, j2, parameters_diff(par, par_prev)

    iterator%ctr_iter = iterator%ctr_iter + 1

    call pop_sub()
  end function iteration_manager
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine iteration_manager_direct(j, par, iterator, dx)
    FLOAT, intent(in) :: j
    type(oct_control_parameters_t), intent(in)  :: par
    type(oct_iterator_t), intent(inout) :: iterator
    FLOAT, optional, intent(in) :: dx

    FLOAT :: j1, j2, fluence, delta
    call push_sub('iter.iteration_manager_direct')

    fluence = parameters_fluence(par)
    j2 = parameters_j2(par)
    j1 = j - j2

    if(present(dx)) then
      delta = dx
    else
      delta = M_ZERO
    end if

    if(iterator%ctr_iter .eq. 0) then
      write(message(1), '(a)') 'Initial-guess field'
      call messages_print_stress(stdout, trim(message(1)))
    else
      write(message(1), '(a,i5)') 'Function evaluation #', iterator%ctr_iter
      call messages_print_stress(stdout, trim(message(1)))
    end if

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

    write(iterator%convergence_iunit, '(i11,4f20.8)')                &
      iterator%ctr_iter, j, j1, j2, delta

    iterator%ctr_iter = iterator%ctr_iter + 1

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
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine oct_iterator_bestpar(par, iterator)
    type(oct_iterator_t), intent(inout)     :: iterator
    type(oct_control_parameters_t), pointer :: par
    par => iterator%best_par
  end subroutine oct_iterator_bestpar
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function oct_iterator_current(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_current = iterator%ctr_iter
  end function oct_iterator_current
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  integer function oct_iterator_maxiter(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_maxiter = iterator%ctr_iter_max
  end function oct_iterator_maxiter
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  FLOAT function oct_iterator_tolerance(iterator)
    type(oct_iterator_t), intent(in)     :: iterator
    oct_iterator_tolerance = iterator%eps
  end function oct_iterator_tolerance
  ! ---------------------------------------------------------


end module opt_control_iter_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
