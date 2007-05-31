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
  ! Gets the fluence of the laser field, defined as:
  ! laser_fluence = \sum_{pol} \integrate_0^T
  !                 laserin(t, pol)^2 dt
  ! ---------------------------------------------------------
  FLOAT function laser_fluence(par)
    type(oct_control_parameters_t), intent(in) :: par
    integer :: i, j
    FLOAT :: t
    call push_sub('opt_control_aux.calc_fluence')

    ! WARNING: This is probably very inefficient; there should be functions in
    ! the tdf module taken care of integrating functions.
    laser_fluence = M_ZERO
    do j = 1, par%no_parameters
      do i = 1, par%ntiter+1
        t = (i-1) * par%dt
        laser_fluence = laser_fluence + abs(tdf(par%f(j), i))**2 
      end do
    end do
    laser_fluence = laser_fluence * par%dt

    call pop_sub()
  end function laser_fluence


  ! ---------------------------------------------------------
  ! Gets the J2 functional (which is the fluence, but weighted
  ! by a penalty function.
  ! ---------------------------------------------------------
  FLOAT function j2_functional(oct, penalty, par) result(j2)
    type(oct_t), intent(in)         :: oct
    type(oct_penalty_t), intent(in) :: penalty
    integer :: i, j
    FLOAT :: t
    type(oct_control_parameters_t), intent(in) :: par
    j2 = M_ZERO
    do j = 1, par%no_parameters
      do i = 1, par%ntiter + 1
        t = (i-1) * par%dt
        j2 = j2 + tdf(penalty%td_penalty(j), i)*abs(tdf(par%f(j), i))**2 
      end do
    end do
    j2 = j2 * par%dt
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
  subroutine write_fieldw(filename, ndim, steps, las, dt)
    ! in w=(2pi f) space
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: ndim
    integer,          intent(in) :: steps
    FLOAT,            intent(in) :: las(1:ndim,0:steps)
    FLOAT,            intent(in) :: dt
  
    integer :: i, iunit
    FLOAT   :: wgrid(0:steps)

    call push_sub('opt_control.write_fieldw')

    call w_lookup(steps+1, dt, wgrid)

    iunit = io_open(filename, action='write')
    do i = 0, steps
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
  subroutine calc_tdfitness(tdt, gr, psi, merit)
    type(td_target_set_t), intent(inout) :: tdt
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    FLOAT,             intent(out):: merit
    integer             :: jj, ik, p, dim

    call push_sub('opt_control.calc_tdfitness')

    merit = M_ZERO
    do jj = 1, tdt%no_tdtargets
      if(tdt%tdtg(jj)%type.eq.oct_tgtype_local) then
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                zmf_integrate(gr%m, tdt%tdtarget(:)* &
                abs(psi%zpsi(:,dim,ik,p))**2)
            end do
          end do
        end do
      else
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              merit = merit + &
                abs(zmf_integrate(gr%m, tdt%tdtarget(:)* &
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
