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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module potential_interpolation_m
  use grid_m
  use global_m
  use lalg_basic_m
  use loct_pointer_m
  use math_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use restart_m

  implicit none

  private
  public ::                                       &
    potential_interpolation_t,                    &
    potential_interpolation_nullify,              &
    potential_interpolation_copy,                 &
    potential_interpolation_init,                 &
    potential_interpolation_end,                  &
    potential_interpolation_run_zero_iter,        &
    potential_interpolation_new,                  &
    potential_interpolation_set,                  &
    potential_interpolation_get,                  &
    potential_interpolation_diff,                 &
    potential_interpolation_interpolate,          &
    potential_interpolation_dump,                 &
    potential_interpolation_load

  type potential_interpolation_t
    FLOAT, pointer      :: v_old(:, :, :) => null()
    FLOAT, pointer      :: Imv_old(:, :, :) => null()
  end type potential_interpolation_t

contains
  
  ! ---------------------------------------------------------
  subroutine potential_interpolation_nullify(this)
    type(potential_interpolation_t), intent(out) :: this
    
    PUSH_SUB(potential_interpolation_nullify)
    
    this%v_old => null()
    this%imv_old => null()

    POP_SUB(potential_interpolation_nullify)
  end subroutine potential_interpolation_nullify
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_copy(vkso, vksi)
    type(potential_interpolation_t), intent(inout) :: vkso
    type(potential_interpolation_t), intent(in)    :: vksi
    
    PUSH_SUB(potential_interpolation_copy)

    call loct_pointer_copy(vkso%v_old, vksi%v_old)
    call loct_pointer_copy(vkso%Imv_old, vksi%Imv_old)

    POP_SUB(potential_interpolation_copy)
  end subroutine potential_interpolation_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_init(potential_interpolation, cmplxscl, np, nspin)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    logical, intent(in) :: cmplxscl
    integer, intent(in) :: np, nspin
    
    PUSH_SUB(potential_interpolation_init)

    SAFE_ALLOCATE(potential_interpolation%v_old(1:np, 1:nspin, 0:3))
    potential_interpolation%v_old(:, :, :) = M_ZERO
    
    if(cmplxscl) then
      SAFE_ALLOCATE(potential_interpolation%Imv_old(1:np, 1:nspin, 0:3))
      potential_interpolation%Imv_old(:, :, :) = M_ZERO
    end if

    POP_SUB(potential_interpolation_init)
  end subroutine potential_interpolation_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_end(potential_interpolation)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    
    PUSH_SUB(potential_interpolation_end)

    ASSERT(associated(potential_interpolation%v_old)) 
    SAFE_DEALLOCATE_P(potential_interpolation%v_old)
    SAFE_DEALLOCATE_P(potential_interpolation%Imv_old) 

    POP_SUB(potential_interpolation_end)
  end subroutine potential_interpolation_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_run_zero_iter(potential_interpolation, np, nspin, vhxc, imvhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: vhxc(:, :)
    FLOAT, optional,   intent(in)    :: imvhxc(:, :)

    integer :: i, ispin, ip
    
    PUSH_SUB(potential_interpolation_run_zero_iter)

    forall(i = 1:3, ispin = 1:nspin, ip = 1:np)
      potential_interpolation%v_old(ip, ispin, i) = vhxc(ip, ispin)
    end forall
    
    if(present(imvhxc)) then
      forall(i = 1:3, ispin = 1:nspin, ip = 1:np)
        potential_interpolation%imv_old(ip, ispin, i) = imvhxc(ip, ispin)
      end forall
    end if

    POP_SUB(potential_interpolation_run_zero_iter)
  end subroutine potential_interpolation_run_zero_iter
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_new(potential_interpolation, np, nspin, time, dt, vhxc, imvhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: time, dt
    FLOAT,             intent(in)    :: vhxc(:, :)
    FLOAT, optional,   intent(in)    :: imvhxc(:, :)

    PUSH_SUB(potential_interpolation_new)

    call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, 2), potential_interpolation%v_old(:, :, 3))
    call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, 1), potential_interpolation%v_old(:, :, 2))
    call lalg_copy(np, nspin, vhxc(:, :),     potential_interpolation%v_old(:, :, 1))
    call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
      potential_interpolation%v_old(:, :, 1:3), time, potential_interpolation%v_old(:, :, 0))

    if(present(imvhxc)) then
      call lalg_copy(np, nspin, potential_interpolation%Imv_old(:, :, 2), potential_interpolation%Imv_old(:, :, 3))
      call lalg_copy(np, nspin, potential_interpolation%Imv_old(:, :, 1), potential_interpolation%Imv_old(:, :, 2))
      call lalg_copy(np, nspin, imvhxc(:, :),     potential_interpolation%Imv_old(:, :, 1))
      call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
        potential_interpolation%imv_old(:, :, 1:3), time, potential_interpolation%imv_old(:, :, 0))
    end if

    POP_SUB(potential_interpolation_new)
  end subroutine potential_interpolation_new
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_set(potential_interpolation, np, nspin, i, vhxc, imvhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(inout) :: imvhxc(:, :)

    PUSH_SUB(potential_interpolation_set)

    call lalg_copy(np, nspin, vhxc, potential_interpolation%v_old(:, :, i))
    
    if(present(imvhxc)) then
      call lalg_copy(np, nspin, imvhxc, potential_interpolation%imv_old(:, :, i))
    end if


    POP_SUB(potential_interpolation_set)
  end subroutine potential_interpolation_set
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_get(potential_interpolation, np, nspin, i, vhxc, imvhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(inout) :: imvhxc(:, :)

    PUSH_SUB(potential_interpolation_set)

    call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, i), vhxc)
    
    if(present(imvhxc)) then
      call lalg_copy(np, nspin, potential_interpolation%imv_old(:, :, i), imvhxc)
    end if

    POP_SUB(potential_interpolation_set)
  end subroutine potential_interpolation_get
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine potential_interpolation_diff(vksold, gr, nspin, vhxc, i, diff)
    type(potential_interpolation_t), intent(inout) :: vksold
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: nspin
    FLOAT,             intent(inout) :: vhxc(:, :)
    integer,           intent(in)    :: i
    FLOAT,             intent(out)   :: diff

    integer :: is
    FLOAT :: d
    FLOAT, allocatable :: dtmp(:)
    PUSH_SUB(potential_interpolation_diff)

    SAFE_ALLOCATE(dtmp(1:gr%mesh%np))

    diff = M_ZERO
    do is = 1, nspin
      dtmp(1:gr%mesh%np) = vhxc(1:gr%mesh%np, is) - vksold%v_old(1:gr%mesh%np, is, i)
      d = dmf_nrm2(gr%mesh, dtmp)
      if(d > diff) diff = d
    end do

    SAFE_DEALLOCATE_A(dtmp)
    POP_SUB(potential_interpolation_diff)
  end subroutine potential_interpolation_diff
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine potential_interpolation_interpolate(potential_interpolation, order, time, dt, t, vhxc, imvhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: order
    FLOAT,             intent(in)    :: time, dt, t
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(inout) :: imvhxc(:, :)


    PUSH_SUB(potential_interpolation_interpolate)

    select case(order)
    case(3)
      call interpolate( (/time, time - dt, time - M_TWO*dt/), potential_interpolation%v_old(:, :, 0:2), t, vhxc(:, :))
      if(present(imvhxc)) then
        call interpolate( (/time, time - dt, time - M_TWO*dt/), potential_interpolation%Imv_old(:, :, 0:2), t, imvhxc(:, :))
      end if
    case(2)
      call interpolate( (/time, time-dt/), potential_interpolation%v_old(:, :, 0:1), t, vhxc(:, :))
      if(present(imvhxc)) then
        call interpolate( (/time, time-dt/), potential_interpolation%imv_old(:, :, 0:1), t, imvhxc(:, :))
      end if
    case default
       ASSERT(.false.)
    end select

    POP_SUB(potential_interpolation_interpolate)
  end subroutine potential_interpolation_interpolate
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_dump(potential_interpolation, restart, gr, cmplxscl, nspin, err2)
    type(potential_interpolation_t), intent(in)    :: potential_interpolation
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    !type(mpi_grp_t) :: mpi_grp
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(potential_interpolation_dump)

    if (cmplxscl) then
      SAFE_ALLOCATE(zv_old(1:gr%mesh%np))
    end if

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a6,i2.2,i3.3)') 'vprev_', ii, is
        if (cmplxscl) then
          zv_old = potential_interpolation%v_old(1:gr%mesh%np, is, ii) &
            + M_zI*potential_interpolation%Imv_old(1:gr%mesh%np, is, ii)
          call zrestart_write_mesh_function(restart, filename, gr%mesh, zv_old, err, use_mpi_grp = .true.)
        else
          call drestart_write_mesh_function(restart, filename, gr%mesh, &
            potential_interpolation%v_old(1:gr%mesh%np, is, ii), err, use_mpi_grp = .true.)
        end if
        ! the unit is energy actually, but this only for restart, and can be kept in atomic units
        ! for simplicity
        if (err /= 0) err2 = err2 + 1
      end do
    end do

    if (cmplxscl) then
      SAFE_DEALLOCATE_A(zv_old)
    end if
    POP_SUB(potential_interpolation_dump)
  end subroutine potential_interpolation_dump
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_load(potential_interpolation, restart, gr, cmplxscl, nspin, err2)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(potential_interpolation_load)

    if (cmplxscl) then
      SAFE_ALLOCATE(zv_old(1:gr%mesh%np))
    end if

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a,i2.2,i3.3)') 'vprev_', ii, is
        if(cmplxscl) then
          call zrestart_read_mesh_function(restart, trim(filename), gr%mesh, zv_old(1:gr%mesh%np), err)
          potential_interpolation%v_old(1:gr%mesh%np, is, ii)   =  real(zv_old(1:gr%mesh%np))
          potential_interpolation%imv_old(1:gr%mesh%np, is, ii) = aimag(zv_old(1:gr%mesh%np))
        else
          call drestart_read_mesh_function(restart, trim(filename), gr%mesh, &
            potential_interpolation%v_old(1:gr%mesh%np, is, ii), err)
        end if
        if (err /= 0) then
          err2 = err2 + 1
          message(1) = "Unable to read VKS restart file '" // trim(filename) // "'"
          call messages_warning(1)
        endif
      end do
    end do

    if (cmplxscl) then
      SAFE_DEALLOCATE_A(zv_old)
    end if
    POP_SUB(potential_interpolation_load)
  end subroutine potential_interpolation_load
  ! ---------------------------------------------------------

end module potential_interpolation_m
  
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
