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


  ! ---------------------------------------------------------
  subroutine oct_iterator_init(iterator, oct)
    type(oct_iterator_t), intent(inout) :: iterator
    type(oct_t), intent(in)             :: oct

    call push_sub('opt_control_iter.oct_iter_init')

    !%Variable OCTEps
    !%Type float
    !%Section Optimal Control
    !%Default 0.001
    !%Description
    !% Define the convergence threshold.
    !% For the monotonically convergent scheme: If the increase of the 
    !% target functional is less then OCTEps the iteration is stopped.
    !% Example
    !% OCTEps = 0.00001
    !%End
    call loct_parse_float(check_inp('OCTEps'), CNST(1.0e-3), iterator%eps)

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

    iterator%old_functional = -CNST(1e10)
    iterator%ctr_iter       = 0

    ALLOCATE(iterator%convergence(4,0:iterator%ctr_iter_max),(iterator%ctr_iter_max+1)*4)

    iterator%bestJ1          = M_ZERO
    iterator%bestJ1_fluence  = M_ZERO
    iterator%bestJ1_J        = M_ZERO
    iterator%bestJ1_ctr_iter = 0

    call pop_sub()
  end subroutine oct_iterator_init


  ! ---------------------------------------------------------
  subroutine oct_iterator_end(iterator)
    type(oct_iterator_t), intent(inout) :: iterator

    call push_sub('opt_control_iter.oct_iter_end')

    deallocate(iterator%convergence)
    nullify(iterator%convergence)

    call pop_sub()
  end subroutine oct_iterator_end


  ! ---------------------------------------------------------
  logical function iteration_manager(oct, gr, par, &
                                     td, psi, target, iterator) result(stoploop)
    type(oct_t), intent(in)             :: oct
    type(grid_t), intent(in)            :: gr
    type(oct_control_parameters_t), intent(in)  :: par
    type(td_t), intent(in)              :: td
    type(states_t), intent(in)          :: psi
    type(target_t), intent(in)          :: target
    type(oct_iterator_t), intent(inout) :: iterator

    FLOAT :: fluence, overlap, jfunctional, j2
    character(len=80)  :: filename

    call push_sub('opt_control.iteration_manager')
    
    stoploop = .false.

    overlap = j1(oct, gr%m, td%max_iter, psi, target)
    fluence = laser_fluence(par)
    j2 = j2_functional(par)
    jfunctional = overlap - j2

    iterator%convergence(1,iterator%ctr_iter) = jfunctional
    iterator%convergence(2,iterator%ctr_iter) = overlap
    iterator%convergence(3,iterator%ctr_iter) = j2
    iterator%convergence(4,iterator%ctr_iter) = tdf(par%td_penalty(1), 1)
    
    if((iterator%ctr_iter .eq. iterator%ctr_iter_max) .or. &
       (iterator%eps>M_ZERO.and.abs(jfunctional-iterator%old_functional) < iterator%eps)) then

      if((iterator%ctr_iter .eq. iterator%ctr_iter_max)) then
        message(1) = "Info: Maximum number of iterations reached"
        call write_info(1)
      endif

      if(iterator%eps > M_ZERO .and. abs(jfunctional-iterator%old_functional) < iterator%eps ) then
        message(1) = "Info: Convergence threshold reached"
        call write_info(1)
      endif
      
      stoploop = .TRUE.
    end if

    write(message(1), '(a,i3)') 'Optimal control iteration #', iterator%ctr_iter
    call messages_print_stress(stdout, trim(message(1)))

    if(oct%mode_fixed_fluence) then
      write(message(1), '(6x,a,f10.5)') " => J1       = ", overlap
      write(message(2), '(6x,a,f10.5)') " => J        = ", jfunctional
      write(message(3), '(6x,a,f10.5)') " => Fluence  = ", fluence
      write(message(4), '(6x,a,f10.5)') " => penalty  = ", &
            real(tdf(par%td_penalty(1), 1), REAL_PRECISION)
      call write_info(4)
    else
      write(message(1), '(6x,a,f10.5)') " => J1       = ", overlap
      write(message(2), '(6x,a,f10.5)') " => J        = ", jfunctional
      write(message(3), '(6x,a,f10.5)') " => Fluence  = ", fluence
      call write_info(3)
    end if
    call messages_print_stress(stdout)

    ! store field with best J1
    if(overlap > iterator%bestJ1) then
      iterator%bestJ1          = overlap
      iterator%bestJ1_J        = jfunctional
      iterator%bestJ1_fluence  = fluence       
      iterator%bestJ1_ctr_iter = iterator%ctr_iter
      ! dump to disc
      write(filename,'(a)') 'opt-control/laser.bestJ1'
      call parameters_write(filename, par)
    end if

    iterator%ctr_iter = iterator%ctr_iter + 1
    iterator%old_functional = jfunctional
    
    call pop_sub()
  end function iteration_manager


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
