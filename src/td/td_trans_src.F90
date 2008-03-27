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
  use io_m

  implicit none

  private
  public ::             &
    source_init,        &
    source_end,         &
    calc_source_wf,     &
    calc_source_wf_new, &
    deriv_coeffs

contains

  ! ---------------------------------------------------------
  ! one sided derivative (1d)
  subroutine deriv_coeffs(order, dx, coeffs, direction)
    integer,                intent(in)    :: order
    FLOAT,                  intent(in)    :: dx   ! h(1)
    FLOAT,                  intent(out)   :: coeffs(:)
    integer,                intent(in)    :: direction

    integer :: k, j, morder
    FLOAT, allocatable :: cc(:,:,:)

    call push_sub('td_trans_src.deriv_coeffs')

    ASSERT(order >= 1)
    ASSERT(direction.eq.1.or.direction.eq.-1)

    morder = 2*order
    ALLOCATE(cc(0:morder, 0:morder, 0:1), (morder+1)*(morder+1)*(1+1))
    call weights(1, morder, cc, direction)

    if (direction.eq.1) then
      coeffs(1) = cc(0, morder, 1) / dx**2
      k = 1
      do j = 1, morder
        k = k + 1
        coeffs(k) = cc( j,   morder, 1) / dx**2
      end do
    else
      coeffs(morder+1) = cc(0, morder, 1) / dx**2
      k = morder+1
      do j = 1, morder
        k = k - 1
        coeffs(k) = cc( j,   morder, 1) / dx**2
      end do
    end if

    deallocate(cc)

    call pop_sub()
  end subroutine deriv_coeffs

  ! ---------------------------------------------------------
  ! Read extended states for time t=0
  subroutine read_psi0_ext(dir, st, psi0, np, gr)
    character(len=*), intent(in) :: dir ! directory with the extended states
    type(states_t),   intent(in) :: st  ! states
    CMPLX,     intent(out)   :: psi0(:, :, :, :, :, :) ! np, (lead=1;center=2), ndim, nst, nik, NLEADS
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

    call pop_sub()
  end subroutine read_psi0_ext

  ! ---------------------------------------------------------
  ! Allocate memory for extended groundstate and calculate eigenstates of the lead
  subroutine source_init(st, src_factor, st_sincos, st_psi0, st_phase, dt, energy, qx, max_iter, inp, nst, nik, gr, order, np)
    type(states_t),  intent(in)  :: st
    CMPLX, pointer      :: src_factor(:)       ! max_iter
    CMPLX, pointer      :: st_sincos(:, :, :)  ! np, (sin=1;cos=2;ext_sin=3;ext_cos=4), nleads
    CMPLX, pointer      :: st_psi0(:, :, :, :, : ,:)  ! np, (lead=1;center=2), ndim, nst, nik, nleads, nleads 
    FLOAT, pointer      :: st_phase(:, :, :)   ! nleads, phase at the interface
    FLOAT, intent(in)   :: dt, energy, qx      ! timestep, energy
    integer, intent(in) :: max_iter, inp, nst, nik, order, np
    type(grid_t),  intent(in)  :: gr

    CMPLX     :: qlr, factor
    integer   :: im, id, length
    FLOAT     :: en, spacing, x, sqrt2
    integer, pointer :: lxyz(:,:)
    FLOAT :: lsize(3)
    FLOAT, allocatable :: coeffs(:,:)
    CMPLX, allocatable :: ext_psi0(:,:)
    !integer     :: ierr

    call push_sub('td_trans_src.source_init')

    sqrt2 = sqrt(M_TWO)
    length = 2*order+1
    ALLOCATE(src_factor(0:max_iter), max_iter+1)
    ALLOCATE(st_sincos(inp, 4, NLEADS), inp*4*NLEADS)
    ALLOCATE(st_psi0(inp, 2, 1, st%st_start:st%st_end, st%d%nik, NLEADS), inp*2*1*st%lnst*st%d%nik*NLEADS )
    ALLOCATE(ext_psi0(np+2, 1) , np+2)
    ALLOCATE(st_phase(nst, nik, NLEADS), NLEADS*nst*nik)
    ALLOCATE(coeffs(length,NLEADS),length*NLEADS)
    
    lxyz => gr%m%lxyz
    lsize = gr%sb%lsize
    spacing = gr%sb%h(1)
    en = energy
    src_factor(:) = M_z0

    factor = (M_z1 - M_zI*M_HALF*dt*en) / (M_z1 + M_zI*M_HALF*dt*en)
    src_factor(0) = M_z1 / (M_z1 + M_zI*M_HALF*dt*en)
    do im=1, max_iter
      src_factor(im) = src_factor(im-1) * factor
    end do
    ! read the extended blocks of the wave funtion
    ! left block
    st_psi0 = M_z0
    !call read_binary(NP+2, ext_psi0, 3, ierr, '/home/nitsche/octopus_n/oct_extended_initial_state.obf')
    !call read_psi0_ext('/home/nitsche/octopus_n/oct_extended_initial_state.obf', st, st_psi0, inp, gr)
    !st_psi0(1, 1, 1, 1, 1 ,1) = ext_psi0(1,1) ! np, ndim, nst, nik, nleads, (lead=1;center=2), nleads 
    !st_psi0(1, 2, 1, 1, 1 ,1) = ext_psi0(2,1) ! np, ndim, nst, nik, nleads, (lead=1;center=2), nleads 
    !st_psi0(1, 1, 1, 1, 1 ,2) = ext_psi0(np+2,1) ! np, ndim, nst, nik, nleads, (lead=1;center=2), nleads 
    !st_psi0(1, 2, 1, 1, 1 ,2) = ext_psi0(np+1,1) ! np, ndim, nst, nik, nleads, (lead=1;center=2), nleads 

    ! now generate the extended eigenstate with separated x-wavefunction
    ! yet only for a box psi=psi(x)*psi(y)*psi(z)
    ! e_tot = e_x + e_y + e_z, we need e_x for q (for box-sized leads)
    ! FIXME: (probably) run trough all possible linear combinations
    ! yet only the groundstate of the separated function is multiplied

    ! subtract the transversal energy from the total energy to get the longitudinal part
    ! which is the energy used for the transport
    do id=2, gr%sb%dim ! FIXME: when the last gridpoint is not the border, recalculate transversal energy
      en = en - M_HALF*(M_PI/(M_TWO*(lsize(id))))**2
      !en = en - M_HALF*(M_PI/(M_TWO*(lsize(id)+spacing)))**2
      !en = en - (M_ONE-cos(M_PI/(M_TWO*(lsize(id)+spacing))))
    end do
    if(en.le.M_ZERO) then
      write(message(1), '(a,f14.6,a)') "The input energy : '", energy, &
                  "' is smaller than the groundstate energy for the given system."
      call write_fatal(1)
    end if
    do im=1, inp
      ! LEFT interface: shift to the right interface for easier calculation
      x = (lxyz(NP-inp+im,1)-1)*spacing - lsize(1)
      st_sincos(im, 1, LEFT) = sqrt2*sin(qx*x)
      st_sincos(im, 2, LEFT) = sqrt2*cos(qx*x)
      x = lxyz(im,1)*spacing + lsize(1)
      st_sincos(im, 3, LEFT) = sqrt2*sin(qx*x)
      st_sincos(im, 4, LEFT) = sqrt2*cos(qx*x)
      ! RIGHT interface: shift to the left interface for easier calculation
      x = (lxyz(im,1)+1)*spacing + lsize(1)
      st_sincos(im, 1, RIGHT) = sqrt2*sin(qx*x)
      st_sincos(im, 2, RIGHT) = sqrt2*cos(qx*x)
      ! in front of the interface
      x = lxyz(NP-inp+im,1)*spacing - lsize(1)
      st_sincos(im, 3, RIGHT) = sqrt2*sin(qx*x)
      st_sincos(im, 4, RIGHT) = sqrt2*cos(qx*x)
      ! now multiply the other seperated wavefunction(s) (f(y)*g(z) or f(y,z))
      do id=2, gr%sb%dim ! for box-sized leads
        qlr = M_PI/(M_TWO*lsize(id))
        st_sincos(im, 1:2, LEFT ) = st_sincos(im, 1:2, LEFT )*cos(qlr*lxyz(inp-im+1,id)*spacing)
        st_sincos(im, 3:4, LEFT ) = st_sincos(im, 3:4, LEFT )*cos(qlr*lxyz(im,id)*spacing)
        st_sincos(im, 1:2, RIGHT) = st_sincos(im, 1:2, RIGHT)*cos(qlr*lxyz(NP-inp+im,id)*spacing)
        st_sincos(im, 3:4, RIGHT) = st_sincos(im, 3:4, RIGHT)*cos(qlr*lxyz(NP+1-im,id)*spacing)
      end do
    end do
    ! TODO 
    ! 1. calculate ALL needed eigenstates for leads (yet only the groundstate)
    ! 2. diagonalize eq. (17) to calculate extended eigenstates
    ! for now just set the initial eigenstate (of the extended system)

    ! Matching algorithm
    ! 1. calculate phaseshift of the lead eigenfunction (eq. (41) in paper)
    ! 1.1 calculate the derivatives at the boundaries
    ! 2. match at interface in order to rescale the central wave function

    deallocate(coeffs, ext_psi0)
    call pop_sub()
  end subroutine source_init


  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term
  ! f(m)*offdiag*st_ext - i dt/2 sum[k=0,m](f(m-k)(mem(k)+mem(k-1))R(cos(ph_s)Z1-sin(ph_s)Z2))
  subroutine calc_source_wf(max_iter, m, il, ph_s, offdiag, factor, mem, dt, np, st_sincos, src_wf)
    integer, intent(in)   :: max_iter ! the maximum timestep
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: il   ! which lead
    FLOAT,   intent(in)   :: ph_s ! phase shift
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: factor(0:np) ! factor(m) = (1 - i dt/2 en)^m / (1 + i dt/2 en)^(m+1)
    CMPLX,   intent(in)   :: mem(1:np, 1:np, 0:max_iter) ! the memory coefficients
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: st_sincos(:, :) ! the sin and cos wave function
    CMPLX,   intent(out)  :: src_wf(:) ! the resulting wave function
    
    CMPLX,   allocatable   :: temp(:)
    integer                :: k

    call push_sub('td_trans_src.calc_source_wf')

    ALLOCATE(temp(np), np)

    ! factor(m)*V^T*st_ext -> src_wf
    src_wf(:) = factor(m)*(cos(ph_s)*st_sincos(:, 1) + sin(ph_s)*st_sincos(:, 2))
    if (il.eq.LEFT) then
      call ztrmv('U','N','N',np,offdiag,np,src_wf,1)
    else
      call ztrmv('L','N','N',np,offdiag,np,src_wf,1)
    end if
!write(*,*) trim(lead_name(il))//' lead wf', cos(ph_s)*st_sincos(:,3) + sin(ph_s)*st_sincos(:,4)
    temp(:) = -M_zI*M_HALF*dt*(cos(ph_s)*st_sincos(:,3) + sin(ph_s)*st_sincos(:,4))
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
  ! calculates the source-wavefunction for the source-term
  ! f(m)*offdiag*psi_lead + i dt/2 sum[k=0,m](f(m-k)(mem(k)+mem(k-1))psi_c(0))
  subroutine calc_source_wf_new(max_iter, m, il, offdiag, factor, mem, dt, np, psi0, src_wf)
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

    call push_sub('td_trans_src.calc_source_wf_new')

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
  end subroutine calc_source_wf_new

  ! ---------------------------------------------------------
  ! calculates the source-wavefunction for the source-term (for sparse mem-coefficients)
  ! f(m)*offdiag*psi_lead + i dt/2 sum[k=0,m](f(m-k)(mem(k)+mem(k-1))psi_c(0))
  subroutine calc_source_wf_sp(max_iter, m, il, offdiag, factor, sp_mem, dt, np, psi0,&
                               mem_s, order, dim, mapping, src_wf)
    integer, intent(in)   :: max_iter ! the maximum timestep
    integer, intent(in)   :: m    ! the m-th timestep
    integer, intent(in)   :: np   ! intf%np, the size of the wave functions
    integer, intent(in)   :: il   ! which lead
    CMPLX,   intent(in)   :: offdiag(:, :) ! the matrix V^T
    CMPLX,   intent(in)   :: factor(0:np) ! factor(m) = (1 - i dt/2 en)^m / (1 + i dt/2 en)^(m+1)
    CMPLX,   intent(in)   :: sp_mem(:)
!    CMPLX,   intent(in)   :: mem(1:np, 1:np, 0:max_iter) ! the memory coefficients
    FLOAT,   intent(in)   :: dt ! time step
    CMPLX,   intent(in)   :: psi0(:, :)! np, (lead=1; center=2)
    CMPLX,   intent(in)   :: mem_s(:, :, :)
    integer, intent(in)   :: order
    integer, intent(in)   :: dim
    integer, intent(in)   :: mapping(:)   ! the mapping
    CMPLX,   intent(out)  :: src_wf(:) ! the resulting wave function
    
    CMPLX,   allocatable   :: tmem(:, :)

    call push_sub('td_trans_src.calc_source_wf_sp')

    ALLOCATE(tmem(np, np), np**2)
    !call make_full_matrix(np, order, dim, sp_mem, mem_s, tmem, mapping)
    call calc_source_wf_new(max_iter, m, il, offdiag, factor, tmem, dt, np, psi0, src_wf)
    deallocate(tmem)

    deallocate(tmem)
    call pop_sub()
  end subroutine calc_source_wf_sp


  ! ---------------------------------------------------------
  ! Free arrays.
  subroutine source_end(src_factor, st_sincos, st_phase, st_psi0)
    CMPLX, pointer :: src_factor(:)   ! max_iter
    CMPLX, pointer :: st_sincos(:, :, :)
    FLOAT, pointer :: st_phase(:, :, :)
    CMPLX, pointer :: st_psi0(:, :, :, :, :, :)

    call push_sub('td_trans_src.source_end')

    if(associated(src_factor)) then
      deallocate(src_factor)
      nullify(src_factor)
    end if

    if(associated(st_sincos)) then
      deallocate(st_sincos)
      nullify(st_sincos)
    end if

    if(associated(st_phase)) then
      deallocate(st_phase)
      nullify(st_phase)
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
