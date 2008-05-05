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
  subroutine read_psi0_ext(st, psi0, np, gr)
    type(states_t),   intent(in) :: st  ! states
    CMPLX,     intent(out)   :: psi0(:, :, :, :, :, :) ! np, (lead=1;center=2), ndim, nst, nik, NLEADS
    integer,   intent(in)    :: np      ! length of the extended block
    type(grid_t),  intent(in):: gr

    CMPLX, allocatable :: tmp(:,:,:,:) ! NP+2*np, ndim, nst, nik
    integer :: ierr, s, ist, ik
    character(len=100) :: filename

    call push_sub('td_trans_src.read_psi0_ext')
    s = (NP+2*np)*st%d%dim
    ALLOCATE(tmp(NP+2*np, st%d%dim, st%st_start:st%st_end,  st%d%nik), s )

    ! try to read from file
    psi0 = M_z0
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        write(filename, '(a,i6.6,a,i6.6,a)') 'ext_eigenstate-', ist, '-', ik, '.obf'
        call read_binary(s, tmp(1:NP+2*np, 1:st%d%dim, ist, ik), 3, ierr, filename)
      end do
    end do
    
!    call read_binary(s, tmp(1:NP+2*np, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik), 3, ierr, trim(dir))
    psi0(1:np, 1, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, LEFT)  = &
        tmp(1:np,1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik)
    psi0(1:np, 2, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, LEFT)  = &
        tmp(np+1:2*np, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik)
    psi0(1:np, 1, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, RIGHT) = &
        tmp(NP+np+1:NP+2*np, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik)
    psi0(1:np, 2, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, RIGHT) = &
        tmp(NP+1:NP+np, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik)
    deallocate(tmp)
    call pop_sub()
  end subroutine read_psi0_ext

  ! ---------------------------------------------------------
  ! Allocate memory for extended groundstate and calculate eigenstates of the lead
  subroutine source_init(st, src_prev, st_psi0, dt, inp, gr)
    type(states_t),  intent(in)  :: st
    CMPLX, pointer      :: src_prev(:, :, :, :, :)! intf%np, ndim, nst, nik, NLEADS
    CMPLX, pointer      :: st_psi0(:, :, :, :, : ,:)  ! np, (lead=1;center=2), ndim, nst, nik, nleads, nleads 
    FLOAT, intent(in)   :: dt      ! timestep
    integer, intent(in) :: inp
    type(grid_t),  intent(in)  :: gr

    call push_sub('td_trans_src.source_init')

    ALLOCATE(src_prev(inp, 1, st%st_start:st%st_end, st%d%nik, NLEADS), inp*NLEADS)
    ALLOCATE(st_psi0(inp, 2, 1, st%st_start:st%st_end, st%d%nik, NLEADS), inp*2*1*st%lnst*st%d%nik*NLEADS )
    
    ! read the extended blocks of the wave funtion
    ! TODO: for all states
    st_psi0 = M_z0
    src_prev = M_z0
    call read_psi0_ext(st, st_psi0, inp, gr)
    ! and do some precalculations

    call pop_sub()
  end subroutine source_init

  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term (recursive variant)
  ! S(m) = f*u(m)*u(m-1)*S(m-1) + dt**2/2*lambda(m,0)/u(m)*f0*(Q(m)+Q(m-1))*psi_c(0)
  subroutine calc_source_wf(maxiter, m, np, il, offdiag, mem, dt, psi0, u, f0, factor, lambda, src)
    integer, intent(in)   :: maxiter
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: il   ! which lead
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: mem(1:np, 1:np) ! the effective memory coefficient
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: psi0(:, :)! np, (lead=1; center=2)
    CMPLX,   intent(in)   :: u(0:maxiter)
    CMPLX,   intent(in)   :: f0
    CMPLX,   intent(in)   :: factor
    CMPLX,   intent(in)   :: lambda
    CMPLX,   intent(inout):: src(:) ! the old wave function in, the new out
    
    CMPLX                  :: tmp, alpha
    integer                :: k

    call push_sub('td_trans_src.calc_source_wf')

    if (m.eq.0) then
      src(:) = psi0(:,1)
      if (il.eq.LEFT) then
        call ztrmv('U','N','N',np,offdiag,np,src,1)
      else
        call ztrmv('L','N','N',np,offdiag,np,src,1)
      end if
      tmp = -M_zI*dt*u(0)*f0
      call zsymv('U', np, tmp*M_zI*M_HALF*dt, mem, np, psi0(:,2), 1, tmp,  src, 1)
    else
      tmp = factor*u(m)*u(m-1)
      alpha = M_HALF*dt**2*f0*lambda/u(m)
      call zsymv('U', np, alpha, mem, np, psi0(:,2), 1, tmp, src, 1)
    end if

    call pop_sub()
  end subroutine calc_source_wf

  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term (for sparse mem-coefficients)
  subroutine calc_source_wf_sp(maxiter, m, np, il, offdiag, sp_mem, dt, order, dim, psi0,&
                               mem_s, mapping, u, f0, factor, lambda, src)
    integer, intent(in)   :: maxiter ! the maximum timestep
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: order
    integer, intent(in)   :: dim
    integer, intent(in)   :: il   ! which lead
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: sp_mem(:)
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: psi0(:, :)! np, (lead=1; center=2)
    CMPLX,   intent(in)   :: mem_s(:, :, :)
    integer, intent(in)   :: mapping(:)   ! the mapping
    CMPLX,   intent(in)   :: u(0:maxiter)
    CMPLX,   intent(in)   :: f0
    CMPLX,   intent(in)   :: factor
    CMPLX,   intent(in)   :: lambda
    CMPLX,   intent(inout):: src(:) ! the old wave function in, the new out
    
    CMPLX,   allocatable   :: tmem(:, :)

    call push_sub('td_trans_src.calc_source_wf_sp')

    ALLOCATE(tmem(np, np), np**2)
    call make_full_matrix(np, order, dim, sp_mem, mem_s, tmem, mapping)
    call calc_source_wf(maxiter, m, np, il, offdiag, tmem, dt, psi0, u, f0, factor, lambda, src)
    deallocate(tmem)

    call pop_sub()
  end subroutine calc_source_wf_sp


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine source_end(src_prev, st_psi0)
    CMPLX, pointer :: src_prev(:, :, :, :, :)
    CMPLX, pointer :: st_psi0(:, :, :, :, :, :)

    call push_sub('td_trans_src.source_end')

    if(associated(src_prev)) then
      deallocate(src_prev)
      nullify(src_prev)
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
