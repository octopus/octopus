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

#include "global.h"

module potential_interpolation_oct_m
  use grid_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use loct_pointer_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use profiling_oct_m
  use restart_oct_m

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
  end type potential_interpolation_t

  integer :: interpolation_steps

contains
  
  ! ---------------------------------------------------------
  subroutine potential_interpolation_nullify(this)
    type(potential_interpolation_t), intent(out) :: this
    
    PUSH_SUB(potential_interpolation_nullify)
    
    this%v_old => null()

    POP_SUB(potential_interpolation_nullify)
  end subroutine potential_interpolation_nullify
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_copy(vkso, vksi)
    type(potential_interpolation_t), intent(inout) :: vkso
    type(potential_interpolation_t), intent(in)    :: vksi
    
    PUSH_SUB(potential_interpolation_copy)

    call loct_pointer_copy(vkso%v_old, vksi%v_old)

    POP_SUB(potential_interpolation_copy)
  end subroutine potential_interpolation_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_init(potential_interpolation, np, nspin, order)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer, intent(in) :: np, nspin
    integer, optional, intent(in) :: order
    
    PUSH_SUB(potential_interpolation_init)

    interpolation_steps = optional_default(order, 3)

    SAFE_ALLOCATE(potential_interpolation%v_old(1:np, 1:nspin, 0:interpolation_steps))
    potential_interpolation%v_old(:, :, :) = M_ZERO
    
    POP_SUB(potential_interpolation_init)
  end subroutine potential_interpolation_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_end(potential_interpolation)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    
    PUSH_SUB(potential_interpolation_end)

    ASSERT(associated(potential_interpolation%v_old)) 
    SAFE_DEALLOCATE_P(potential_interpolation%v_old)

    POP_SUB(potential_interpolation_end)
  end subroutine potential_interpolation_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_run_zero_iter(potential_interpolation, np, nspin, vhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: vhxc(:, :)

    integer :: i, ispin, ip
    
    PUSH_SUB(potential_interpolation_run_zero_iter)

    forall(i = 1:interpolation_steps, ispin = 1:nspin, ip = 1:np)
      potential_interpolation%v_old(ip, ispin, i) = vhxc(ip, ispin)
    end forall
    
    POP_SUB(potential_interpolation_run_zero_iter)
  end subroutine potential_interpolation_run_zero_iter
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_new(potential_interpolation, np, nspin, time, dt, vhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: time, dt
    FLOAT,             intent(in)    :: vhxc(:, :)

    FLOAT, allocatable :: times(:)
    integer :: j

    PUSH_SUB(potential_interpolation_new)

    SAFE_ALLOCATE(times(interpolation_steps))
    do j = 1, interpolation_steps
      times(j) = time - j*dt
    end do

    do j = interpolation_steps, 2, -1
      call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, j-1), potential_interpolation%v_old(:, :, j))
    end do
    call lalg_copy(np, nspin, vhxc(:, :),     potential_interpolation%v_old(:, :, 1))
    call interpolate( times, potential_interpolation%v_old(:, :, 1:interpolation_steps), &
         time, potential_interpolation%v_old(:, :, 0))

    SAFE_DEALLOCATE_A(times)
    POP_SUB(potential_interpolation_new)
  end subroutine potential_interpolation_new
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_set(potential_interpolation, np, nspin, i, vhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)

    PUSH_SUB(potential_interpolation_set)

    call lalg_copy(np, nspin, vhxc, potential_interpolation%v_old(:, :, i))

    POP_SUB(potential_interpolation_set)
  end subroutine potential_interpolation_set
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_get(potential_interpolation, np, nspin, i, vhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)

    PUSH_SUB(potential_interpolation_get)

    call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, i), vhxc)
    
    POP_SUB(potential_interpolation_get)
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
  subroutine potential_interpolation_interpolate(potential_interpolation, order, time, dt, t, vhxc)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: order
    FLOAT,             intent(in)    :: time, dt, t
    FLOAT,             intent(inout) :: vhxc(:, :)

    integer :: j
    FLOAT, allocatable :: times(:)

    PUSH_SUB(potential_interpolation_interpolate)

    SAFE_ALLOCATE(times(interpolation_steps))
    do j = 1, interpolation_steps
      times(j) = time - (j-1)*dt
    end do

    call interpolate( times(1:order), potential_interpolation%v_old(:, :, 0:order-1), t, vhxc(:, :))

    SAFE_DEALLOCATE_A(times)
    POP_SUB(potential_interpolation_interpolate)
  end subroutine potential_interpolation_interpolate
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_dump(potential_interpolation, restart, gr, nspin, err2)
    type(potential_interpolation_t), intent(in)    :: potential_interpolation
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    !type(mpi_grp_t) :: mpi_grp
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(potential_interpolation_dump)

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a6,i2.2,i3.3)') 'vprev_', ii, is
        call drestart_write_mesh_function(restart, filename, gr%mesh, &
          potential_interpolation%v_old(1:gr%mesh%np, is, ii), err)
        ! the unit is energy actually, but this only for restart, and can be kept in atomic units
        ! for simplicity
        if (err /= 0) err2 = err2 + 1
      end do
    end do

    POP_SUB(potential_interpolation_dump)
  end subroutine potential_interpolation_dump
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_load(potential_interpolation, restart, gr, nspin, err2)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(potential_interpolation_load)

    err2 = 0
    do ii = 1, interpolation_steps - 1
      do is = 1, nspin
        write(filename,'(a,i2.2,i3.3)') 'vprev_', ii, is
        call drestart_read_mesh_function(restart, trim(filename), gr%mesh, &
          potential_interpolation%v_old(1:gr%mesh%np, is, ii), err)
        if (err /= 0) then
          err2 = err2 + 1
          message(1) = "Unable to read VKS restart file '" // trim(filename) // "'"
          call messages_warning(1)
        end if
      end do
    end do

    POP_SUB(potential_interpolation_load)
  end subroutine potential_interpolation_load
  ! ---------------------------------------------------------

end module potential_interpolation_oct_m
  
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
