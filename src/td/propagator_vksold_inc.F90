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
!! $Id: propagator_vksold_inc.F90 $

  ! ---------------------------------------------------------
  subroutine vksinterp_nullify(this)
    type(vksinterp_t), intent(out) :: this
    PUSH_SUB(vksinterp_nullify)
    this%v_old => null()
    this%imv_old => null()
    POP_SUB(vksinterp_nullify)
  end subroutine vksinterp_nullify
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_copy(vkso, vksi)
    type(vksinterp_t), intent(inout) :: vkso
    type(vksinterp_t), intent(in)    :: vksi
    PUSH_SUB(vksinterp_copy)

    call loct_pointer_copy(vkso%v_old, vksi%v_old)
    call loct_pointer_copy(vkso%Imv_old, vksi%Imv_old)

    POP_SUB(vksinterp_copy)
  end subroutine vksinterp_copy
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_init(vksinterp, cmplxscl, np, nspin)
    type(vksinterp_t), intent(inout) :: vksinterp
    logical, intent(in) :: cmplxscl
    integer, intent(in) :: np, nspin
    PUSH_SUB(vksinterp_init)

    SAFE_ALLOCATE(vksinterp%v_old(1:np, 1:nspin, 0:3))
    vksinterp%v_old(:, :, :) = M_ZERO
    if(cmplxscl) then
      SAFE_ALLOCATE(vksinterp%Imv_old(1:np, 1:nspin, 0:3))
      vksinterp%Imv_old(:, :, :) = M_ZERO
    end if

    POP_SUB(vksinterp_init)
  end subroutine vksinterp_init
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_end(vksinterp)
    type(vksinterp_t), intent(inout) :: vksinterp
    PUSH_SUB(vksinterp_end)

    ASSERT(associated(vksinterp%v_old)) 
    SAFE_DEALLOCATE_P(vksinterp%v_old)
    SAFE_DEALLOCATE_P(vksinterp%Imv_old) 

    POP_SUB(vksinterp_end)
  end subroutine vksinterp_end
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_run_zero_iter(vksinterp, cmplxscl, np, nspin, vhxc, imvhxc)
    type(vksinterp_t), intent(inout) :: vksinterp
    logical, intent(in) :: cmplxscl
    integer, intent(in) :: np, nspin
    FLOAT,   intent(in) :: vhxc(:, :), imvhxc(:, :)

    integer :: i, ispin, ip
    PUSH_SUB(vksinterp_run_zero_iter)

    forall(i = 1:3, ispin = 1:nspin, ip = 1:np)
      vksinterp%v_old(ip, ispin, i) = vhxc(ip, ispin)
    end forall
    if(cmplxscl) then
      forall(i = 1:3, ispin = 1:nspin, ip = 1:np)
        vksinterp%imv_old(ip, ispin, i) = imvhxc(ip, ispin)
      end forall
    end if

    POP_SUB(vksinterp_run_zero_iter)
  end subroutine vksinterp_run_zero_iter
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_new(vksinterp, cmplxscl, np, nspin, vhxc, imvhxc, time, dt)
    type(vksinterp_t), intent(inout) :: vksinterp
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(in)    :: vhxc(:, :), imvhxc(:, :)
    FLOAT,             intent(in)    :: time, dt

    PUSH_SUB(vksinterp_new)

    call lalg_copy(np, nspin, vksinterp%v_old(:, :, 2), vksinterp%v_old(:, :, 3))
    call lalg_copy(np, nspin, vksinterp%v_old(:, :, 1), vksinterp%v_old(:, :, 2))
    call lalg_copy(np, nspin, vhxc(:, :),     vksinterp%v_old(:, :, 1))
    call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
      vksinterp%v_old(:, :, 1:3), time, vksinterp%v_old(:, :, 0))
    if(cmplxscl) then
      call lalg_copy(np, nspin, vksinterp%Imv_old(:, :, 2), vksinterp%Imv_old(:, :, 3))
        call lalg_copy(np, nspin, vksinterp%Imv_old(:, :, 1), vksinterp%Imv_old(:, :, 2))
        call lalg_copy(np, nspin, imvhxc(:, :),     vksinterp%Imv_old(:, :, 1))
        call interpolate( (/time - dt, time - M_TWO*dt, time - M_THREE*dt/), &
          vksinterp%imv_old(:, :, 1:3), time, vksinterp%imv_old(:, :, 0))
      end if

    POP_SUB(vksinterp_new)
  end subroutine vksinterp_new
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_set(vksinterp, cmplxscl, np, nspin, vhxc, imvhxc, i)
    type(vksinterp_t), intent(inout) :: vksinterp
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(inout) :: vhxc(:, :), imvhxc(:, :)
    integer,           intent(in)    :: i

    PUSH_SUB(vksinterp_set)

    call lalg_copy(np, nspin, vhxc, vksinterp%v_old(:, :, i))
    if(cmplxscl) then
      call lalg_copy(np, nspin, imvhxc, vksinterp%imv_old(:, :, i))
    end if


    POP_SUB(vksinterp_set)
  end subroutine vksinterp_set
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_get(vksinterp, cmplxscl, np, nspin, vhxc, imvhxc, i)
    type(vksinterp_t), intent(inout) :: vksinterp
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(inout) :: vhxc(:, :), imvhxc(:, :)
    integer,           intent(in)    :: i

    PUSH_SUB(vksinterp_set)

    call lalg_copy(np, nspin, vksinterp%v_old(:, :, i), vhxc)
    if(cmplxscl) then
      call lalg_copy(np, nspin, vksinterp%imv_old(:, :, i), imvhxc)
    end if

    POP_SUB(vksinterp_set)
  end subroutine vksinterp_get
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine vksinterp_diff(vksold, gr, nspin, vhxc, i, diff)
    type(vksinterp_t), intent(inout) :: vksold
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: nspin
    FLOAT,             intent(inout) :: vhxc(:, :)
    integer,           intent(in)    :: i
    FLOAT,             intent(out)   :: diff

    integer :: is
    FLOAT :: d
    FLOAT, allocatable :: dtmp(:)
    PUSH_SUB(vksinterp_diff)

    SAFE_ALLOCATE(dtmp(1:gr%mesh%np))

    diff = M_ZERO
    do is = 1, nspin
      dtmp(1:gr%mesh%np) = vhxc(1:gr%mesh%np, is) - vksold%v_old(1:gr%mesh%np, is, i)
      d = dmf_nrm2(gr%mesh, dtmp)
      if(d > diff) diff = d
    end do

    SAFE_DEALLOCATE_A(dtmp)
    POP_SUB(vksinterp_diff)
  end subroutine vksinterp_diff
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine vksinterp_interpolate(vksinterp, order, cmplxscl, np, nspin, vhxc, imvhxc, time, dt, t)
    type(vksinterp_t), intent(inout) :: vksinterp
    integer,           intent(in)    :: order
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: np, nspin
    FLOAT,             intent(inout) :: vhxc(:, :), imvhxc(:, :)
    FLOAT,             intent(in)    :: time, dt, t

    PUSH_SUB(vksinterp_interpolate)

    select case(order)
    case(3)
      call interpolate( (/time, time - dt, time - M_TWO*dt/), vksinterp%v_old(:, :, 0:2), t, vhxc(:, :))
      if(cmplxscl) &
        call interpolate( (/time, time - dt, time - M_TWO*dt/), vksinterp%Imv_old(:, :, 0:2), t, imvhxc(:, :))
    case(2)
      call interpolate( (/time, time-dt/), vksinterp%v_old(:, :, 0:1), t, vhxc(:, :))
      if(cmplxscl) &
        call interpolate( (/time, time-dt/), vksinterp%imv_old(:, :, 0:1), t, imvhxc(:, :))
    case default
       ASSERT(.false.)
    end select

    POP_SUB(vksinterp_interpolate)
  end subroutine vksinterp_interpolate
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_dump(vksinterp, restart, gr, cmplxscl, nspin, err2)
    type(vksinterp_t), intent(in)    :: vksinterp
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(vksinterp_dump)

    if (cmplxscl) then
      SAFE_ALLOCATE(zv_old(gr%mesh%np))
    end if

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a6,i2.2,i3.3)') 'vprev_', ii, is
        if (cmplxscl) then
          zv_old = vksinterp%v_old(1:gr%mesh%np, is, ii) + M_zI * vksinterp%Imv_old(1:gr%mesh%np, is, ii)
          call zrestart_write_mesh_function(restart, filename, gr%mesh, zv_old, err)
        else
          call drestart_write_mesh_function(restart, filename, gr%mesh, vksinterp%v_old(1:gr%mesh%np, is, ii), err)
        end if
        ! the unit is energy actually, but this only for restart, and can be kept in atomic units
        ! for simplicity
        if (err /= 0) err2 = err2 + 1
      end do
    end do

    if (cmplxscl) then
      SAFE_DEALLOCATE_A(zv_old)
    end if
    POP_SUB(vksinterp_dump)
  end subroutine vksinterp_dump
  ! ---------------------------------------------------------

  ! ---------------------------------------------------------
  subroutine vksinterp_load(vksinterp, restart, gr, cmplxscl, nspin, err2)
    type(vksinterp_t), intent(inout) :: vksinterp
    type(restart_t),   intent(in)    :: restart
    type(grid_t),      intent(in)    :: gr
    logical,           intent(in)    :: cmplxscl
    integer,           intent(in)    :: nspin
    integer,           intent(out)   :: err2

    integer :: ii, is, err
    character(len=256) :: filename
    CMPLX, allocatable :: zv_old(:)
    PUSH_SUB(vksinterp_load)

    if (cmplxscl) then
      SAFE_ALLOCATE(zv_old(1:gr%mesh%np))
    end if

    err2 = 0
    do ii = 1, 2
      do is = 1, nspin
        write(filename,'(a,i2.2,i3.3)') 'vprev_', ii, is
        if(cmplxscl) then
          call zrestart_read_mesh_function(restart, trim(filename), gr%mesh, zv_old(1:gr%mesh%np), err)
          vksinterp%v_old(1:gr%mesh%np, is, ii)   =  real(zv_old(1:gr%mesh%np))
          vksinterp%imv_old(1:gr%mesh%np, is, ii) = aimag(zv_old(1:gr%mesh%np))
        else
          call drestart_read_mesh_function(restart, trim(filename), gr%mesh, vksinterp%v_old(1:gr%mesh%np, is, ii), err)
        end if
        if (err /= 0) err2 = err2 + 1
      end do
    end do

    if (cmplxscl) then
      SAFE_DEALLOCATE_A(zv_old)
    end if
    POP_SUB(vksinterp_load)
  end subroutine vksinterp_load
  ! ---------------------------------------------------------

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
