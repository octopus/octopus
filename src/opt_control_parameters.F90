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

#include "global.h"


  ! ---------------------------------------------------------
  subroutine parameters_init(cp, gr, td)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    call push_sub('opt_control_parameters.parameters_init')

    ALLOCATE(cp%laser(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter))
    cp%laser(:, :) = M_ZERO

  end subroutine parameters_init


  ! ---------------------------------------------------------
  subroutine parameters_set(cp, gr, td, ep)
    type(oct_control_parameters_t), intent(inout) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    type(epot_t), intent(in) :: ep

    integer :: jj, kk
    FLOAT   :: t
    FLOAT, allocatable :: field(:)

    call push_sub('opt_control_parameters.parameters_set')

    ALLOCATE(field(1:NDIM), NDIM)
    cp%laser = M_ZERO;

      do jj = 0, 2*td%max_iter
        t = td%dt*(jj-1)/M_TWO
        !i = int(abs(M_TWO*t/l(1)%dt) + M_HALF)
        do kk=1, ep%no_lasers
          call laser_field(gr%sb, ep%no_lasers, ep%lasers, t, field)
          cp%laser(:,jj) = cp%laser(:,jj) + field(:)
        end do
      enddo

    deallocate(field)
    call pop_sub()
  end subroutine parameters_set


  ! ---------------------------------------------------------
  subroutine parameters_to_h(cp, gr, td, ep)
    type(oct_control_parameters_t), intent(in) :: cp
    type(grid_t), intent(in) :: gr
    type(td_t),   intent(in) :: td
    type(epot_t), intent(inout) :: ep

    call push_sub('opt_control_paramters.parameters_to_h')

    ! TODO Probably is this better if this was not a pointer?
    ep%lasers(1)%numerical => cp%laser

    call pop_sub()
  end subroutine parameters_to_h


  ! ---------------------------------------------------------
  subroutine parameters_end(cp)
    type(oct_control_parameters_t), intent(inout) :: cp
    call push_sub('opt_control_parameters.parameters_end')

    deallocate(cp%laser)
    nullify(cp%laser)

    call pop_sub()
  end subroutine parameters_end


  ! ---------------------------------------------------------
  subroutine parameters_write(filename, steps, n_dim, las, dt)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: steps
    integer,          intent(in) :: n_dim
    FLOAT,            intent(in) :: las(1:n_dim,0:2*steps)
    FLOAT,            intent(in) :: dt
    integer :: i, iunit
    FLOAT   :: tgrid(0:2*steps)

    call push_sub('opt_control_parameters.parameters_write')
    
    call t_lookup(2*steps+1,dt/CNST(2.0),tgrid)
    iunit = io_open(filename, action='write')
    do i = 0, 2*steps, 2
       write(iunit, '(4es30.16e4)') tgrid(i), las(:, i)
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

    CMPLX, allocatable :: rpsi(:, :)
    
    call push_sub('opt_control.update_field')
    
    ! case SWITCH
    ! check dipole moment function - > velocity gauge

    ALLOCATE(rpsi(gr%m%np_part, psi%d%dim), gr%m%np_part*psi%d%dim)

    ! TODO This should be a product between Slater determinants if we want
    ! to make it work for many-particle systems.      
    d2 = M_z0 
    do ik = 1, psi%d%nik
      do p  = psi%st_start, psi%st_end
        do pol = 1, NDIM
          do dim = 1, psi%d%dim
            rpsi(:, dim) = psi%zpsi(:, dim, p, ik)*oct%laser_pol(pol)*gr%m%x(:, pol)
          end do
          d2(pol) = zstates_dotp(gr%m, psi%d%dim, chi%zpsi(:, :, p, ik), rpsi)
        end do
      end do
    end do
    deallocate(rpsi)

    d1 = M_z1
    if(oct%algorithm_type .eq. oct_algorithm_zbr98) d1 = zstates_mpdotp(gr%m, psi, chi)

    ! Q: How to distinguish between the cases ?
    !if((NDIM-dof).eq.1) then
    !l(1:NDIM, 2*iter) = M_ONE/tdpenalty(1:NDIM,2*iter)*real(m_z1/(m_z2*m_zI) * &
    !(laser_pol(1:NDIM)*(d1*d2(1:NDIM)+conjg(d1)*conjg(d2(1:NDIM)))))
    !endif
    !if((NDIM-dof).eq.0) then
    cp%laser(1:NDIM, 2*iter) = aimag(d1*d2(1:NDIM))/penalty%tdpenalty(1:NDIM,2*iter)
    !endif
    ! extrapolate to t+-dt/2
    i = int(sign(M_ONE, td%dt))
    if(iter==0.or.iter==td%max_iter) then
      cp%laser(1:NDIM, 2*iter+  i) = cp%laser(1:NDIM, 2*iter)
      cp%laser(1:NDIM, 2*iter+2*i) = cp%laser(1:NDIM, 2*iter)
    else
      cp%laser(1:NDIM, 2*iter+  i) = M_HALF*(M_THREE*cp%laser(1:NDIM, 2*iter) -       cp%laser(1:NDIM, 2*iter-2*i))
      cp%laser(1:NDIM, 2*iter+2*i) = M_HALF*( M_FOUR*cp%laser(1:NDIM, 2*iter) - M_TWO*cp%laser(1:NDIM, 2*iter-2*i))
    end if

    call pop_sub()
  end subroutine update_field
