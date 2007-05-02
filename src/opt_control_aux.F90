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
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function j2_functional(oct, laser, dt) result(j2)
    type(oct_t), intent(in) :: oct
    FLOAT, intent(in)       :: laser(:,:)
    FLOAT, intent(in)       :: dt
    j2 = SUM(oct%tdpenalty * laser**2) * abs(dt)/M_TWO
  end function j2_functional


  ! ---------------------------------------------------------
  FLOAT function overlap_function(oct, m, td_fitness, max_iter, psi, target_st)
    type(oct_t), intent(in)    :: oct
    type(mesh_t), intent(in)   :: m
    FLOAT, intent(in)          :: td_fitness(:)
    integer, intent(in)        :: max_iter
    type(states_t), intent(in) :: psi, target_st

    call push_sub('opt_control.overlap_function')

    if(oct%targetmode==oct_targetmode_td) then
       ! 1/T * int(<Psi| O | Psi>)
       overlap_function = SUM(td_fitness) / real(max_iter, REAL_PRECISION) 
    else
      overlap_function = abs(zstates_mpdotp(m, psi, target_st))
    end if

    call pop_sub()
  end function overlap_function


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
  subroutine target_calc(oct, gr, targetst, psi_in, chi_out)
    type(oct_t),       intent(in)  :: oct
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
            if(oct%algorithm_type == oct_algorithm_zr98) &
              chi_out%zpsi(:,:,p,ik) = olap*targetst%zpsi(:, :, p, ik)
            if(oct%algorithm_type == oct_algorithm_zbr98) then
              chi_out%zpsi(:,:,p,ik) = targetst%zpsi(:, :, p, ik)
            end if
            if(oct%algorithm_type == oct_algorithm_wg05) &
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


  ! ---------------------------------------------------------
  subroutine calc_tdfitness(tg, gr, psi, tdtarget, merit)
    type(td_target_t), intent(in) :: tg(:)
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    CMPLX,             intent(in) :: tdtarget(:)
    FLOAT,             intent(out):: merit
    integer             :: jj, ik, p, dim

    call push_sub('opt_control.calc_tdfitness')

    merit = M_ZERO
    do jj=1, size(tg)
      if(tg(jj)%type.eq.oct_tgtype_local) then
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                zmf_integrate(gr%m, tdtarget(:)* &
                abs(psi%zpsi(:,dim,ik,p))**2)
            end do
          end do
        end do
      else
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                abs(zmf_integrate(gr%m, tdtarget(:)* &
                conjg(psi%zpsi(:,dim,ik,p))))**2
            end do
          end do
        end do
      end if
    end do

    call pop_sub()
  end subroutine calc_tdfitness

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
