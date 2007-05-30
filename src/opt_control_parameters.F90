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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, no_parameters, td)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer, intent(in) :: no_parameters   
    type(td_t),   intent(in) :: td
    integer :: j

    call push_sub('opt_control_parameters.parameters_init')

    cp%no_parameters = no_parameters
    cp%dt = td%dt
    cp%ntiter = 2*td%max_iter
    ALLOCATE(cp%f(cp%no_parameters), cp%no_parameters)
    do j = 1, cp%no_parameters
      call tdf_init_numerical(cp%f(j), cp%ntiter, M_HALF*cp%dt)
    end do

    ALLOCATE(cp%laser_pol(MAX_DIM, cp%no_parameters), MAX_DIM*cp%no_parameters)
    cp%laser_pol = M_z0

    call pop_sub()
  end subroutine parameters_init


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(epot_t), intent(in) :: ep
    integer :: j

    call push_sub('opt_control_parameters.parameters_set')

    do j = 1, cp%no_parameters
      cp%f(j) = ep%lasers(j)%f
      cp%laser_pol(:, j) = ep%lasers(j)%pol(:)
    end do

    call pop_sub()
  end subroutine parameters_set


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(epot_t), intent(inout) :: ep

    integer :: j

    call push_sub('opt_control_paramters.parameters_to_h')

    do j = 1, cp%no_parameters
      ep%lasers(j)%f = cp%f(j)
      ep%lasers(j)%pol(:) = cp%laser_pol(:, j)
    end do

    call pop_sub()
  end subroutine parameters_to_h


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    integer :: j

    call push_sub('opt_control_parameters.parameters_end')

    do j = 1, cp%no_parameters
      call tdf_end(cp%f(j))
    end do
    deallocate(cp%f); nullify(cp%f)
    deallocate(cp%laser_pol)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, cp)
    character(len=*), intent(in) :: filename
    type(oct_control_parameters_t), intent(in) :: cp

    integer :: i, iunit
    FLOAT :: t

    call push_sub('opt_control_parameters.parameters_write')

    iunit = io_open(filename, action='write')
    do i = 1, cp%ntiter + 1
      t = (i-1)*cp%dt*M_HALF
      write(iunit, '(3es20.8e3)') t, tdf(cp%f(1), t)
    end do
    call io_close(iunit)

    call pop_sub()
  end subroutine parameters_write


  ! ---------------------------------------------------------
  subroutine update_field(oct, penalty, iter, cp, gr, td, psi, chi)
    type(oct_t), intent(in)    :: oct
    type(oct_penalty_t), intent(in) :: penalty
    integer, intent(in)        :: iter
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in)   :: gr
    type(td_t), intent(in)     :: td
    type(states_t), intent(in) :: psi
    type(states_t), intent(in) :: chi
    
    CMPLX :: d1
    CMPLX :: d2(NDIM)
    integer :: ik, p, dim, i, pol
    FLOAT :: value

    CMPLX, allocatable :: rpsi(:, :)
    
    call push_sub('opt_control.update_field')
    
    ALLOCATE(rpsi(gr%m%np_part, psi%d%dim), gr%m%np_part*psi%d%dim)

    ! TODO This should be a product between Slater determinants.
    d2 = M_z0 
    do ik = 1, psi%d%nik
      do p  = psi%st_start, psi%st_end
        do pol = 1, NDIM
          do dim = 1, psi%d%dim
            rpsi(:, dim) = psi%zpsi(:, dim, p, ik)*cp%laser_pol(pol, 1)*gr%m%x(:, pol)
          end do
          d2(pol) = zstates_dotp(gr%m, psi%d%dim, chi%zpsi(:, :, p, ik), rpsi)
        end do
      end do
    end do
    deallocate(rpsi)

    d1 = M_z1
    if(oct%algorithm_type .eq. oct_algorithm_zbr98) d1 = zstates_mpdotp(gr%m, psi, chi)

    value = aimag(d1*d2(1))/tdf(penalty%td_penalty(1), 2*iter+1)
    call tdf_set_numerical(cp%f(1), 2*iter+1, value)
    i = int(sign(M_ONE, td%dt))
    if(iter==0.or.iter==td%max_iter) then
      call tdf_set_numerical(cp%f(1), 2*iter+i+1, value)
      call tdf_set_numerical(cp%f(1), 2*iter+2*i+1, value)
    else
      value = M_HALF*(M_THREE*tdf(cp%f(1), 2*iter+1) - tdf(cp%f(1), 2*iter-2*i+1))
      call tdf_set_numerical(cp%f(1), 2*iter+i+1, value)
      value = M_HALF*( M_FOUR*tdf(cp%f(1), 2*iter+1) - M_TWO*tdf(cp%f(1), 2*iter-2*i+1))
      call tdf_set_numerical(cp%f(1), 2*iter+2*i+1, value)
    end if

    call pop_sub()
  end subroutine update_field

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
