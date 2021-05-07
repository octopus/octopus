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
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m

  implicit none

  private
  public ::                                       &
    potential_interpolation_t,                    &
    potential_interpolation_copy,                 &
    potential_interpolation_init,                 &
    potential_interpolation_end,                  &
    potential_interpolation_run_zero_iter,        &
    potential_interpolation_new,                  &
    potential_interpolation_set,                  &
    potential_interpolation_get,                  &
    potential_interpolation_interpolate,          &
    potential_interpolation_dump,                 &
    potential_interpolation_load

  type potential_interpolation_t
    private
    FLOAT, allocatable, public :: v_old(:, :, :)
    FLOAT, allocatable, public :: vtau_old(:, :, :)
    logical                    :: mgga_with_exc
  end type potential_interpolation_t

  integer :: interpolation_steps

contains

  ! ---------------------------------------------------------
  subroutine potential_interpolation_copy(vkso, vksi)
    type(potential_interpolation_t), intent(inout) :: vkso
    type(potential_interpolation_t), intent(in)    :: vksi
    
    PUSH_SUB(potential_interpolation_copy)

    SAFE_ALLOCATE_SOURCE_A(vkso%v_old, vksi%v_old)
    SAFE_ALLOCATE_SOURCE_A(vkso%vtau_old, vksi%vtau_old)

    POP_SUB(potential_interpolation_copy)
  end subroutine potential_interpolation_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_init(potential_interpolation, np, nspin, mgga_with_exc, order)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer, intent(in) :: np, nspin
    logical, intent(in) :: mgga_with_exc
    integer, optional, intent(in) :: order
    
    PUSH_SUB(potential_interpolation_init)

    interpolation_steps = optional_default(order, 3)

    SAFE_ALLOCATE(potential_interpolation%v_old(1:np, 1:nspin, 0:interpolation_steps))
    potential_interpolation%v_old(:, :, :) = M_ZERO
    
    potential_interpolation%mgga_with_exc = mgga_with_exc

    if(potential_interpolation%mgga_with_exc) then

      SAFE_ALLOCATE(potential_interpolation%vtau_old(1:np, 1:nspin, 0:interpolation_steps))
      potential_interpolation%vtau_old(:, :, :) = M_ZERO

   end if

    POP_SUB(potential_interpolation_init)
  end subroutine potential_interpolation_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_end(potential_interpolation)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    
    PUSH_SUB(potential_interpolation_end)

    ASSERT(allocated(potential_interpolation%v_old)) 
    SAFE_DEALLOCATE_A(potential_interpolation%v_old)
    SAFE_DEALLOCATE_A(potential_interpolation%vtau_old)

    POP_SUB(potential_interpolation_end)
  end subroutine potential_interpolation_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_run_zero_iter(potential_interpolation, np, nspin, vhxc, vtau)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: vhxc(:, :)
    FLOAT, optional,   intent(in)    :: vtau(:, :)

    integer :: i, ispin, ip

    PUSH_SUB(potential_interpolation_run_zero_iter)

    do i = 1, interpolation_steps
      do ispin = 1, nspin
        do ip = 1, np
          potential_interpolation%v_old(ip, ispin, i) = vhxc(ip, ispin)
        end do
      end do
    end do

   if(present(vtau)) then
      do i = 1, interpolation_steps
        do ispin = 1, nspin
          do ip = 1, np
            potential_interpolation%vtau_old(ip, ispin, i) = vtau(ip, ispin)
          end do
        end do
      end do
    end if

    POP_SUB(potential_interpolation_run_zero_iter)
  end subroutine potential_interpolation_run_zero_iter
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_new(potential_interpolation, np, nspin, time, dt, vhxc, vtau)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: time, dt
    FLOAT,             intent(in)    :: vhxc(:, :)
    FLOAT, optional,   intent(in)    :: vtau(:, :)

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

    if(present(vtau)) then
      do j = interpolation_steps, 2, -1
        call lalg_copy(np, nspin, potential_interpolation%vtau_old(:, :, j-1), potential_interpolation%vtau_old(:, :, j))
      end do
      call lalg_copy(np, nspin, vtau(:, :),     potential_interpolation%vtau_old(:, :, 1))
      call interpolate( times, potential_interpolation%vtau_old(:, :, 1:interpolation_steps), &
           time, potential_interpolation%vtau_old(:, :, 0))
    end if

    SAFE_DEALLOCATE_A(times)
    POP_SUB(potential_interpolation_new)
  end subroutine potential_interpolation_new
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_set(potential_interpolation, np, nspin, i, vhxc, vtau)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(in)    :: vtau(:, :)

    PUSH_SUB(potential_interpolation_set)

    call lalg_copy(np, nspin, vhxc, potential_interpolation%v_old(:, :, i))
    
    if(present(vtau)) then
      call lalg_copy(np, nspin, vtau, potential_interpolation%vtau_old(:, :, i))
    end if

    POP_SUB(potential_interpolation_set)
  end subroutine potential_interpolation_set
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_get(potential_interpolation, np, nspin, i, vhxc, vtau)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: np, nspin
    integer,           intent(in)    :: i
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(inout) :: vtau(:, :)

    PUSH_SUB(potential_interpolation_get)

    call lalg_copy(np, nspin, potential_interpolation%v_old(:, :, i), vhxc)

    if(present(vtau)) then
      call lalg_copy(np, nspin, potential_interpolation%vtau_old(:, :, i), vtau)
    end if

    POP_SUB(potential_interpolation_get)
  end subroutine potential_interpolation_get
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine potential_interpolation_interpolate(potential_interpolation, order, time, dt, t, vhxc, vtau)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    integer,           intent(in)    :: order
    FLOAT,             intent(in)    :: time, dt, t
    FLOAT,             intent(inout) :: vhxc(:, :)
    FLOAT, optional,   intent(inout) :: vtau(:, :)

    integer :: j
    FLOAT, allocatable :: times(:)

    PUSH_SUB(potential_interpolation_interpolate)

    SAFE_ALLOCATE(times(interpolation_steps))
    do j = 1, interpolation_steps
      times(j) = time - (j-1)*dt
    end do

    call interpolate( times(1:order), potential_interpolation%v_old(:, :, 0:order-1), t, vhxc(:, :))

    if(present(vtau)) then
      call interpolate( times(1:order), potential_interpolation%vtau_old(:, :, 0:order-1), t, vtau(:, :))
    end if

    SAFE_DEALLOCATE_A(times)
    POP_SUB(potential_interpolation_interpolate)
  end subroutine potential_interpolation_interpolate
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_dump(potential_interpolation, space, restart, mesh, nspin, err2)
    type(potential_interpolation_t), intent(in)    :: potential_interpolation
    type(restart_t),   intent(in)    :: restart
    type(space_t),     intent(in)    :: space
    type(mesh_t),      intent(in)    :: mesh
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    !type(mpi_grp_t) :: mpi_grp

    PUSH_SUB(potential_interpolation_dump)

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a6,i2.2,i3.3)') 'vprev_', ii, is
        call drestart_write_mesh_function(restart, space, filename, mesh, &
          potential_interpolation%v_old(1:mesh%np, is, ii), err)
        ! the unit is energy actually, but this only for restart, and can be kept in atomic units
        ! for simplicity
        if (err /= 0) err2 = err2 + 1
      end do
    end do


    if(potential_interpolation%mgga_with_exc) then
      err2 = 0
      do ii = 1, 2
        do is = 1, nspin
          write(filename,'(a6,i2.2,i3.3)') 'vtauprev_', ii, is
          call drestart_write_mesh_function(restart, space, filename, mesh, &
            potential_interpolation%vtau_old(1:mesh%np, is, ii), err)
          ! the unit is energy actually, but this only for restart, and can be kept in atomic units
          ! for simplicity
          if (err /= 0) err2 = err2 + 1
        end do
      end do
    end if

    POP_SUB(potential_interpolation_dump)
  end subroutine potential_interpolation_dump
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine potential_interpolation_load(potential_interpolation, namespace, space, restart, mesh, nspin, err2)
    type(potential_interpolation_t), intent(inout) :: potential_interpolation
    type(namespace_t), intent(in)    :: namespace
    type(space_t),     intent(in)    :: space
    type(restart_t),   intent(in)    :: restart
    type(mesh_t),      intent(in)    :: mesh
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename

    PUSH_SUB(potential_interpolation_load)

    err2 = 0
    do ii = 1, interpolation_steps - 1
      do is = 1, nspin
        write(filename,'(a,i2.2,i3.3)') 'vprev_', ii, is
        call drestart_read_mesh_function(restart, space, filename, mesh, &
          potential_interpolation%v_old(1:mesh%np, is, ii), err)
        if (err /= 0) then
          err2 = err2 + 1
          message(1) = "Unable to read VKS restart file '" // trim(filename) // "'"
          call messages_warning(1, namespace=namespace)
        end if
      end do
    end do

    if(potential_interpolation%mgga_with_exc) then
      err2 = 0
      do ii = 1, interpolation_steps - 1
        do is = 1, nspin
          write(filename,'(a,i2.2,i3.3)') 'vtauprev_', ii, is
          call drestart_read_mesh_function(restart, space, filename, mesh, &
            potential_interpolation%vtau_old(1:mesh%np, is, ii), err)
          if (err /= 0) then
            err2 = err2 + 1
            message(1) = "Unable to read VKS restart file '" // trim(filename) // "'"
            call messages_warning(1, namespace=namespace)
          end if
        end do
      end do
    end if

    POP_SUB(potential_interpolation_load)
  end subroutine potential_interpolation_load
  ! ---------------------------------------------------------

end module potential_interpolation_oct_m
  
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
