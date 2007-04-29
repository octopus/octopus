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
!! $Id: opt_control.F90 2868 2007-04-26 17:11:47Z acastro $

#include "global.h"

  ! ---------------------------------------------------------
  ! Gets the fluence of the laser field, defined as:
  ! laser_fluence = \sum_{pol} \integrate_0^T
  !                 laserin(t, pol)^2 dt
  ! ---------------------------------------------------------
  FLOAT function laser_fluence(laserin, dt)
    FLOAT, intent(in) :: laserin(:,:)
    FLOAT, intent(in) :: dt
    call push_sub('opt_control_aux.calc_fluence')

    laser_fluence = SUM(laserin**2) * abs(dt)/M_TWO      

    call pop_sub()
  end function laser_fluence


  ! ---------------------------------------------------------
  subroutine write_field(filename, steps, n_dim, las, dt)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: steps
    integer,          intent(in) :: n_dim
    FLOAT,            intent(in) :: las(1:n_dim,0:2*steps)
    FLOAT,            intent(in) :: dt
    integer :: i, iunit
    FLOAT   :: tgrid(0:2*steps)

    call push_sub('opt_control.write_field')
    
    call t_lookup(2*steps+1,dt/CNST(2.0),tgrid)
    iunit = io_open(filename, action='write')
    do i = 0, 2*steps, 2
       write(iunit, '(4es30.16e4)') tgrid(i), las(:, i)
    end do
    call io_close(iunit)

    call pop_sub()
  end subroutine write_field


  ! ---------------------------------------------------------
  subroutine write_fieldw(filename, ndim, steps, las, dt)
    ! in w=(2pi f) space
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: ndim
    integer,          intent(in) :: steps
    FLOAT,            intent(in) :: las(1:ndim,0:2*steps)
    FLOAT,            intent(in) :: dt
  
    integer :: i, iunit
    FLOAT   :: wgrid(0:2*steps)

    call push_sub('opt_control.write_fieldw')

    call w_lookup(2*steps+1, dt, wgrid)

    iunit = io_open(filename, action='write')
    do i = 0, 2*steps
       write(iunit, '(4es30.16e4)') wgrid(i), las(:, i)
    end do
    call io_close(iunit)
    
    call pop_sub()
  end subroutine write_fieldw



  ! ---------------------------------------------------------
  ! calculate chi = \hat{O} psi
  ! do loop <target_st|Psi> for all States
  subroutine target_calc(oct, method, gr, targetst, psi_in, chi_out)
    type(oct_t),       intent(in)  :: oct
    integer,           intent(in)  :: method
    type(grid_t),      intent(in)  :: gr
    type(states_t),    intent(in)  :: targetst, psi_in
    type(states_t),    intent(inout) :: chi_out
    
    CMPLX   :: olap
    integer :: ik, p, dim
  
    call push_sub('opt_control.target_calc')

    if(oct%targetmode==oct_targetmode_static) then
      if(oct%totype.eq.oct_tg_local) then ! only zr98 and wg05
        do ik = 1, psi_in%d%nik
          do p  = psi_in%st_start, psi_in%st_end
            do dim = 1, psi_in%d%dim
              ! multiply orbtials with local operator
              ! FIXME: for multiple particles 1,1,1 -> dim,p,ik
              chi_out%zpsi(:,dim,p,ik) = targetst%zpsi(:, 1, 1 , 1)*psi_in%zpsi(:, dim, p, ik)
            end do
          end do
        end do
      else ! totype nonlocal (all other totypes)
        do ik = 1, psi_in%d%nik
          do p  = psi_in%st_start, psi_in%st_end
            olap = zstates_dotp(gr%m, targetst%d%dim, targetst%zpsi(:, :, p, ik), &
                                                      psi_in%zpsi(:, :, p, ik))
            if(method == oct_algorithm_zr98) &
              chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
            if(method == oct_algorithm_zbr98) then
              chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:, :, p, ik)
            end if
            if(method == oct_algorithm_wg05) &
              chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
          end do
        end do
      end if
    else
      ! time-dependent target
      message(1) = 'Info: Time-dependent Target selected'
      call write_info(1)
      chi_out%zpsi = M_z0
    end if

!      do ik = 1, psi%d%nik
!        do p  = psi%st_start, psi%st_end
!          olap = M_z0     
!          olap = zstates_dotp(gr%m, psi_in%d%dim,  targetst%zpsi(:,:, p, ik), psi_in%zpsi(:,:, p, ik))
!          if(method == 'ZR98') &
!            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
!          if(method == 'ZBR98') &
!            chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:,:,p,ik)
!          if(method == 'WG05') &
!            chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:,:,p,ik)
!        end do
!      end do
      
    call pop_sub()
  end subroutine target_calc

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
