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
!! $Id: td_trans_src.F90 3030 2008-01-03 10:30:05Z nitsche $

! Calculation of the source term for the modified Crank-Nicholson propagator.

#include "global.h"

module td_trans_src_m
  use datasets_m
  use global_m
  use grid_m
  use lalg_adv_m
  use lalg_basic_m
  use loct_parser_m
  use loct_m
  use math_m
  use messages_m
  use nl_operator_m
  use states_m
  use system_m
  use td_trans_intf_m
  use varinfo_m
  use derivatives_m
  use td_trans_mem_m
  use io_m

  implicit none

  private
  public ::             &
    source_init,        &
    source_end,         &
    calc_source_wf,     &
    calc_source_wf_sp

contains

  ! ---------------------------------------------------------
  ! Read extended states for time t=0
  subroutine read_psi0_ext(dir, st, psi0, energy, np, gr)
    character(len=*), intent(in) :: dir ! directory with the extended states
    type(states_t),   intent(in) :: st  ! states
    CMPLX,     intent(out)   :: psi0(:, :, :, :, :, :) ! np, (lead=1;center=2), ndim, nst, nik, NLEADS
    FLOAT,     intent(out)   :: energy ! H(lead)*psi(lead)+H(lead,C)*psi(C) == energy*psi(lead)
    integer,   intent(in)    :: np      ! length of the extended block
    type(grid_t),  intent(in):: gr

    CMPLX, allocatable :: tmp(:,:,:,:) ! NP+2*np, ndim, nst, nik
    integer :: ierr

    call push_sub('td_trans_src.read_psi0_ext')
    ALLOCATE(tmp(NP+2*np, st%d%dim, st%st_start:st%st_end,  st%d%nik), (NP+2*np)*st%d%dim*st%lnst*st%d%nik*NLEADS )

    ! try to read from file
    psi0 = M_z0
    call read_binary(NP+2*np, tmp, 3, ierr, trim(dir))
    psi0(1:np, 1, 1, 1, 1, LEFT)  = tmp(1:np,1,1,1) ! np, (lead=1;center=2), ndim, nst, nik
    psi0(1:np, 2, 1, 1, 1, LEFT)  = tmp(np+1:2*np,1,1,1) ! np, (lead=1;center=2), ndim, nst, nik
    psi0(1:np, 1, 1, 1, 1, RIGHT) = tmp(NP+np+1:NP+2*np,1,1,1) ! np, (lead=1;center=2), ndim, nst, nik
    psi0(1:np, 2, 1, 1, 1, RIGHT) = tmp(NP+1:NP+np,1,1,1) ! np, (lead=1;center=2), ndim, nst, nik 

    ! FIXME
    energy = 3.5*(M_PI/(gr%sb%lsize(2)+gr%sb%h(2)))**2/M_EIGHT

    deallocate(tmp)
    call pop_sub()
  end subroutine read_psi0_ext

  ! ---------------------------------------------------------
  ! Allocate memory for extended groundstate and calculate eigenstates of the lead
  subroutine source_init(st, src_factor, st_psi0, dt, energy, max_iter, inp, gr)
    type(states_t),  intent(in)  :: st
    CMPLX, pointer      :: src_factor(:)       ! max_iter
    CMPLX, pointer      :: st_psi0(:, :, :, :, : ,:)  ! np, (lead=1;center=2), ndim, nst, nik, nleads, nleads 
    FLOAT, intent(in)   :: dt      ! timestep
    FLOAT, intent(out)  :: energy
    integer, intent(in) :: max_iter, inp
    type(grid_t),  intent(in)  :: gr

    CMPLX     :: factor
    integer   :: im

    call push_sub('td_trans_src.source_init')

    ALLOCATE(src_factor(0:max_iter), max_iter+1)
    ALLOCATE(st_psi0(inp, 2, 1, st%st_start:st%st_end, st%d%nik, NLEADS), inp*2*1*st%lnst*st%d%nik*NLEADS )
    
    ! read the extended blocks of the wave funtion
    ! TODO: for all states
    st_psi0 = M_z0
    call read_psi0_ext('ext_eigenstate.obf', st, st_psi0, energy, inp, gr)
    ! and do some precalculations
    src_factor(:) = M_z0
    factor = (M_z1 - M_zI*M_HALF*dt*energy) / (M_z1 + M_zI*M_HALF*dt*energy)
    src_factor(0) = M_z1 / (M_z1 + M_zI*M_HALF*dt*energy)
    do im=1, max_iter
      src_factor(im) = src_factor(im-1) * factor
    end do
    call pop_sub()
  end subroutine source_init

  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term
  ! f(m)*offdiag*psi_lead + i dt/2 sum[k=0,m](f(m-k)(mem(k)+mem(k-1))psi_c(0))
  subroutine calc_source_wf(max_iter, m, il, offdiag, factor, mem, dt, np, psi0, src_wf)
    integer, intent(in)   :: max_iter ! the maximum timestep
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: il   ! which lead
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: factor(0:np) ! factor(m) = (1 - i dt/2 en)^m / (1 + i dt/2 en)^(m+1)
    CMPLX,   intent(in)   :: mem(1:np, 1:np, 0:max_iter) ! the memory coefficients
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: psi0(:, :)! np, (lead=1; center=2)
    CMPLX,   intent(out)  :: src_wf(:) ! the resulting wave function
    
    CMPLX,   allocatable   :: temp(:)
    integer                :: k

    call push_sub('td_trans_src.calc_source_wf')

    ALLOCATE(temp(np), np)
! FIXME: use recursive relation
    ! factor(m)*V^T*psi_lead -> src_wf
    src_wf(:) = factor(m)*psi0(:,1)
    if (il.eq.LEFT) then
      call ztrmv('U','N','N',np,offdiag,np,src_wf,1)
    else
      call ztrmv('L','N','N',np,offdiag,np,src_wf,1)
    end if
    temp(:) = M_zI*M_HALF*dt*psi0(:,2)
    do k=0, m
      call zsymv('U', np, factor(m-k), mem(:,:,k), np, temp, 1, M_z1, src_wf, 1)
      if (k.gt.0) then
        call zsymv('U', np, factor(m-k), mem(:,:,k-1), np, temp, 1, M_z1, src_wf, 1)
      end if
    end do

    deallocate(temp)
    call pop_sub()
  end subroutine calc_source_wf

  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term (for sparse mem-coefficients)
  ! f(m)*offdiag*psi_lead + i dt/2 sum[k=0,m](f(m-k)(mem(k)+mem(k-1))psi_c(0))
  subroutine calc_source_wf_sp(max_iter, m, il, offdiag, factor, sp_mem, dt, np, order, dim, psi0,&
                               mem_s, mapping, src_wf)
    integer, intent(in)   :: max_iter ! the maximum timestep
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: order
    integer, intent(in)   :: dim
    integer, intent(in)   :: il   ! which lead
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: factor(0:np) ! factor(m) = (1 - i dt/2 en)^m / (1 + i dt/2 en)^(m+1)
    CMPLX,   intent(in)   :: sp_mem(1:np*order, 0:max_iter)
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: psi0(:, :)! np, (lead=1; center=2)
    CMPLX,   intent(in)   :: mem_s(:, :, :)
    integer, intent(in)   :: mapping(:)   ! the mapping
    CMPLX,   intent(out)  :: src_wf(:) ! the resulting wave function
    
    CMPLX,   allocatable   :: tmem(:, :)

    call push_sub('td_trans_src.calc_source_wf_sp')

    ALLOCATE(tmem(np, np), np**2)
!FIXME: when implementing the recursive relation fix this also
   ! call make_full_matrix(np, order, dim, sp_mem, mem_s, tmem, mapping)
    call calc_source_wf(max_iter, m, il, offdiag, factor, tmem, dt, np, psi0, src_wf)
    deallocate(tmem)

    deallocate(tmem)
    call pop_sub()
  end subroutine calc_source_wf_sp


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine source_end(src_factor, st_psi0)
    CMPLX, pointer :: src_factor(:)   ! max_iter
    CMPLX, pointer :: st_psi0(:, :, :, :, :, :)

    call push_sub('td_trans_src.source_end')

    if(associated(src_factor)) then
      deallocate(src_factor)
      nullify(src_factor)
    end if

    if(associated(st_psi0)) then
      deallocate(st_psi0)
      nullify(st_psi0)
    end if

    call pop_sub()
  end subroutine source_end
end module td_trans_src_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
