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


  ! ---------------------------------------------------------
  ! Calculates the J1 functional, i.e.:
  ! <Psi(T)|\hat{O}|Psi(T) in the time-independent
  ! case, or else \int_0^T dt <Psi(t)|\hat{O}(t)|Psi(t) in 
  ! the time-dependent case.
  ! ---------------------------------------------------------
  FLOAT function j1(oct, m, td_fitness, max_iter, psi, target_st)
    type(oct_t), intent(in)    :: oct
    type(mesh_t), intent(in)   :: m
    FLOAT, intent(in)          :: td_fitness(:)
    integer, intent(in)        :: max_iter
    type(states_t), intent(in) :: psi, target_st

    call push_sub('opt_control.overlap_function')

    if(oct%targetmode==oct_targetmode_td) then
       ! 1/T * int(<Psi| O | Psi>)
       j1 = sum(td_fitness) / real(max_iter, REAL_PRECISION) 
    else
      j1 = abs(zstates_mpdotp(m, psi, target_st))**2
    end if

    call pop_sub()
  end function j1


  ! ---------------------------------------------------------
  ! calculate chi = \hat{O} psi
  ! do loop <target_st|Psi> for all States
  subroutine target_calc(oct, gr, target, psi_in, chi_out)
    type(oct_t),       intent(in)  :: oct
    type(grid_t),      intent(in)  :: gr
    type(states_t),    intent(in)  :: psi_in
    type(target_t),    intent(in)  :: target

    type(states_t),    intent(inout) :: chi_out
    
    CMPLX   :: olap
    integer :: ik, p, dim
  
    call push_sub('opt_control.target_calc')

    if(oct%targetmode==oct_targetmode_static) then
      if(target%totype.eq.oct_tg_local) then ! only zr98 and wg05
        do ik = 1, psi_in%d%nik
          do p  = psi_in%st_start, psi_in%st_end
            do dim = 1, psi_in%d%dim
              chi_out%zpsi(:,dim,p,ik) = target%st%zpsi(:, 1, 1 , 1)*psi_in%zpsi(:, dim, p, ik)
            end do
          end do
        end do
      else ! totype nonlocal 
        olap = zstates_mpdotp(gr%m, target%st, psi_in)
        do ik = 1, psi_in%d%nik
          do p  = psi_in%st_start, psi_in%st_end
            if(oct%algorithm_type == oct_algorithm_zr98) &
              chi_out%zpsi(:,:,p,ik) = olap*target%st%zpsi(:, :, p, ik)
            if(oct%algorithm_type == oct_algorithm_zbr98) then
              chi_out%zpsi(:,:,p,ik) = target%st%zpsi(:, :, p, ik)
            end if
            if(oct%algorithm_type == oct_algorithm_wg05) &
              chi_out%zpsi(:,:,p,ik) = olap*target%st%zpsi(:, :, p, ik)
          end do
        end do
      end if
    else
      ! time-dependent target
      message(1) = 'Info: Time-dependent Target selected'
      call write_info(1)
      chi_out%zpsi = M_z0
    end if
      
    call pop_sub()
  end subroutine target_calc


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
