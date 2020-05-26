!! Copyright (C) 2008 X. Andrade, & ARW
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

module thermal_gradient_oct_m
  use global_oct_m
  use grid_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_dim_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                               &
    thermal_gradient_t,                        &
    thermal_gradient_nullify,                  &
    thermal_gradient_init,                     &
    thermal_gradient_init_vec_pot,             &
    thermal_gradient_is_applied,               &
    thermal_gradient_set_vec_pot,              &
    thermal_gradient_set_vec_pot_vel,          &
    thermal_gradient_get_vec_pot,              &
    thermal_gradient_get_vec_pot_vel,          &
    thermal_gradient_get_vec_pot_acc,          &
    thermal_gradient_propagate,                &
    thermal_gradient_propagate_vel,            &
    thermal_gradient_get_energy,               &
    thermal_gradient_dump,                     &
    thermal_gradient_load,                     &
    thermal_gradient_end,                      &
    thermal_gradient_get_force

  type thermal_gradient_t
    FLOAT   :: vecpot(1:MAX_DIM)   
    FLOAT   :: vecpot_vel(1:MAX_DIM)
    FLOAT   :: vecpot_acc(1:MAX_DIM)    
    FLOAT   :: vecpot_kick(1:MAX_DIM)
    FLOAT   :: force(1:MAX_DIM)
    FLOAT   :: wp2
    integer :: ndim
    logical :: with_thermal_gradient
 !   integer :: dynamics
    FLOAT   :: kicktime 
  end type thermal_gradient_t

contains

  subroutine thermal_gradient_nullify(this)
    type(thermal_gradient_t), intent(out) :: this

    PUSH_SUB(thermal_gradient_nullify)
    this%with_thermal_gradient = .false.

    POP_SUB(thermal_gradient_nullify)
  end subroutine thermal_gradient_nullify

  ! ---------------------------------------------------------
  subroutine thermal_gradient_init(this, sb)
    type(thermal_gradient_t),     intent(out)   :: this
    type(simul_box_t),       intent(in)    :: sb

    integer :: ii, iop
    type(block_t) :: blk

    PUSH_SUB(thermal_gradient_init)

    this%with_thermal_gradient = .false.
    this%vecpot = M_ZERO
    this%vecpot_vel = M_ZERO
    this%vecpot_acc = M_ZERO
    this%vecpot_kick = M_ZERO
    this%force = M_ZERO
    this%ndim = sb%dim

    !%Variable ThermalGradient
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The thermal gradient is used to include a uniform (but time-dependent)
    !% external thermal gradient in a time-dependent run for
    !% a periodic system. An optional second row specifies the initial
    !% value for the gradient of the temperature over the temperature (which is set
    !% to zero by default). By default this field is not included.
    !% If <tt>KPointsUseSymmetries = yes</tt>, then <tt>SymmetryBreakDir</tt>
    !% must be set in the same direction.
    !%End
    ! Read the thermal gradient
 
    if(parse_block('ThermalGradient', blk) == 0) then

      this%with_thermal_gradient = .true.

      do ii = 1, this%ndim
        call parse_block_float(blk, 0, ii - 1, this%vecpot_kick(ii))
      end do

      call parse_block_end(blk)
      if(.not. simul_box_is_periodic(sb)) then
        message(1) = "ThermalGradient is intended for periodic systems."
        call messages_warning(1)
      end if

      if(sb%kpoints%use_symmetries) then
        do iop = 1, symmetries_number(sb%symm)
          if(iop == symmetries_identity_index(sb%symm)) cycle
          if(.not. symm_op_invariant_cart(sb%symm%ops(iop), this%vecpot, CNST(1e-5))) then
            message(1) = "ThermalGradient breaks (at least) one of the symmetries used to reduce the k-points."
            !message(2) = "Set SymmetryBreakDir equal to GaugeVectorField."
            call messages_fatal(2)
          end if
        end do
      end if

    end if


    this%kicktime = M_ZERO
    !call parse_variable('GaugeFieldDelay', M_ZERO, this%kicktime)

    if(abs(this%kicktime) <= M_EPSILON) then
       this%vecpot(1:this%ndim) = this%vecpot_kick(1:this%ndim)
    endif

    POP_SUB(thermal_gradient_init)
  end subroutine thermal_gradient_init


  ! ---------------------------------------------------------
  subroutine thermal_gradient_end(this)
    type(thermal_gradient_t),     intent(inout) :: this

    PUSH_SUB(thermal_gradient_end)

    this%with_thermal_gradient = .false.

    POP_SUB(thermal_gradient_end)
  end subroutine thermal_gradient_end


  ! ---------------------------------------------------------
  logical pure function thermal_gradient_is_applied(this) result(is_applied)
    type(thermal_gradient_t),  intent(in) :: this

    is_applied = this%with_thermal_gradient
  end function thermal_gradient_is_applied


  ! ---------------------------------------------------------
  subroutine thermal_gradient_set_vec_pot(this, vec_pot)
    type(thermal_gradient_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(:) !< (this%ndim)

    PUSH_SUB(thermal_gradient_set_vec_pot)
    this%vecpot(1:this%ndim) = vec_pot(1:this%ndim)

    POP_SUB(thermal_gradient_set_vec_pot)
  end subroutine thermal_gradient_set_vec_pot


  ! ---------------------------------------------------------
  subroutine thermal_gradient_set_vec_pot_vel(this, vec_pot_vel)
    type(thermal_gradient_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(:) !< (this%ndim)

    PUSH_SUB(thermal_gradient_set_vec_pot_vel)
    this%vecpot_vel(1:this%ndim) = vec_pot_vel(1:this%ndim)

    POP_SUB(thermal_gradient_set_vec_pot_vel)
  end subroutine thermal_gradient_set_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine thermal_gradient_get_vec_pot(this, vec_pot)
    type(thermal_gradient_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot(:) !< (this%ndim)

    PUSH_SUB(thermal_gradient_get_vec_pot)
    vec_pot(1:this%ndim) = this%vecpot(1:this%ndim)

    POP_SUB(thermal_gradient_get_vec_pot)
  end subroutine thermal_gradient_get_vec_pot


  ! ---------------------------------------------------------
  subroutine thermal_gradient_get_vec_pot_vel(this, vec_pot_vel)
    type(thermal_gradient_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_vel(:) !< (this%ndim)

    PUSH_SUB(thermal_gradient_get_vec_pot_vel)
    vec_pot_vel(1:this%ndim) = this%vecpot_vel(1:this%ndim)

    POP_SUB(thermal_gradient_get_vec_pot_vel)
  end subroutine thermal_gradient_get_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine thermal_gradient_get_vec_pot_acc(this, vec_pot_acc)
    type(thermal_gradient_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_acc(:) !< (this%ndim)

    PUSH_SUB(thermal_gradient_get_vec_pot_acc)
    vec_pot_acc(1:this%ndim) = this%vecpot_acc(1:this%ndim)

    POP_SUB(thermal_gradient_get_vec_pot_acc)
  end subroutine thermal_gradient_get_vec_pot_acc

  ! ---------------------------------------------------------
  subroutine thermal_gradient_propagate(this, dt, time)
    type(thermal_gradient_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: dt
    FLOAT,                intent(in)    :: time
    
    logical, save :: warning_shown = .false.
    integer :: idim

    PUSH_SUB(thermal_gradient_propagate)
  
    this%vecpot_acc(1:this%ndim) = this%force(1:this%ndim)

    ! apply kick, in case kicktime=0 the kick has already been applied
    if(this%kicktime > M_ZERO .and. time-dt <= this%kicktime .and. time >= this%kicktime )  then
      this%vecpot(1:this%ndim) = this%vecpot(1:this%ndim) +  this%vecpot_kick(1:this%ndim)
      call messages_write('     ----------------  Applying gauge kick  ----------------')
      call messages_info()
    endif
     
    this%vecpot(1:this%ndim) = this%vecpot(1:this%ndim) + dt * this%vecpot_vel(1:this%ndim) + &
      M_HALF * dt**2 * this%force(1:this%ndim)

    !In the case of a kick, the induced field could not be higher than the initial kick
    do idim = 1, this%ndim
      if(.not. warning_shown .and. this%vecpot_kick(idim) /= M_ZERO .and.  &
         abs(this%vecpot(idim))> abs(this%vecpot_kick(idim))*1.01 .and. .not. this%kicktime > M_ZERO ) then
       
        warning_shown = .true.

        write(message(1),'(a)') 'It seems that the thermal-field might be diverging. You should probably check'
        write(message(2),'(a)') 'the simulation parameters, in particular the number of k-points.'
        call messages_warning(2)
      end if
    end do
    POP_SUB(thermal_gradient_propagate)
  end subroutine thermal_gradient_propagate

  ! ---------------------------------------------------------
  subroutine thermal_gradient_propagate_vel(this, dt)
    type(thermal_gradient_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: dt

    PUSH_SUB(thermal_gradient_propagate_vel)

    this%vecpot_vel(1:this%ndim) = this%vecpot_vel(1:this%ndim) + &
      M_HALF * dt * (this%vecpot_acc(1:this%ndim) + this%force(1:this%ndim))

    POP_SUB(thermal_gradient_propagate_vel)
  end subroutine thermal_gradient_propagate_vel


  ! ---------------------------------------------------------
  subroutine thermal_gradient_init_vec_pot(this, sb, st)
    type(thermal_gradient_t),  intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(states_t),       intent(in)    :: st
    
    PUSH_SUB(thermal_gradient_init_vec_pot)

    this%wp2 = M_FOUR*M_PI*st%qtot/sb%rcell_volume

    write (message(1), '(a,f12.6,a)') "Info: Electron-gas plasmon frequency", &
         units_from_atomic(units_out%energy, sqrt(this%wp2)), " ["//trim(units_abbrev(units_out%energy))//"]"
    call messages_info(1)

    POP_SUB(thermal_gradient_init_vec_pot)
  end subroutine thermal_gradient_init_vec_pot

  ! ---------------------------------------------------------
  FLOAT function thermal_gradient_get_energy(this, sb) result(energy)
    type(thermal_gradient_t),  intent(in)    :: this
    type(simul_box_t),    intent(in)    :: sb

    PUSH_SUB(thermal_gradient_get_energy)
    energy = sb%rcell_volume / (CNST(8.0) * M_PI * P_c**2) * sum(this%vecpot_vel(1:this%ndim)**2)

    POP_SUB(thermal_gradient_get_energy)
  end function thermal_gradient_get_energy


  ! ---------------------------------------------------------
  subroutine thermal_gradient_dump(restart, tfield, ierr)
    type(restart_t),      intent(in)  :: restart
    type(thermal_gradient_t),  intent(in)  :: tfield
    integer,              intent(out) :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)
    
    PUSH_SUB(thermal_gradient_dump)

    ierr = 0
    
    if (restart_skip(restart)) then
      POP_SUB(thermal_gradient_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing thermal gradient restart."
      call messages_info(1)
    end if

    SAFE_ALLOCATE(vecpot(1:tfield%ndim, 1:2))
    vecpot = M_ZERO
    call thermal_gradient_get_vec_pot(tfield, vecpot(:, 1))
    call thermal_gradient_get_vec_pot_vel(tfield, vecpot(:, 2))

    call drestart_write_binary(restart, "thermal_gradient", 2*tfield%ndim, vecpot, err)
    SAFE_DEALLOCATE_A(vecpot)
    if (err /= 0) ierr = ierr + 1

    if (debug%info) then
      message(1) = "Debug: Writing thermal gradient restart done."
      call messages_info(1)
    end if

    POP_SUB(thermal_gradient_dump)
  end subroutine thermal_gradient_dump


  ! ---------------------------------------------------------
  subroutine thermal_gradient_load(restart, tfield, ierr)
    type(restart_t),      intent(in)    :: restart
    type(thermal_gradient_t),  intent(inout) :: tfield
    integer,              intent(out)   :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)
    
    PUSH_SUB(thermal_gradient_load)

    ierr = 0
    
    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(thermal_gradient_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading thermal gradient restart."
      call messages_info(1)
    end if

    SAFE_ALLOCATE(vecpot(1:tfield%ndim, 1:2))
    call drestart_read_binary(restart, "thermal_gradient", 2*tfield%ndim, vecpot, err)
    if (err /= 0) ierr = ierr + 1

    call thermal_gradient_set_vec_pot(tfield, vecpot(:,1))
    call thermal_gradient_set_vec_pot_vel(tfield, vecpot(:,2))
    SAFE_DEALLOCATE_A(vecpot)
    
    if (debug%info) then
      message(1) = "Debug: Reading thermal gradient restart done."
      call messages_info(1)
    end if

    POP_SUB(thermal_gradient_load)
  end subroutine thermal_gradient_load

  ! ---------------------------------------------------------

  subroutine thermal_gradient_get_force(this, gr, st)
    type(thermal_gradient_t),  intent(inout)    :: this
    type(grid_t),         intent(in)    :: gr
    type(states_t),       intent(in)    :: st

    integer :: idir,ispin,istot

    PUSH_SUB(thermal_gradient_get_force)

    this%force(1:gr%sb%dim) = M_ZERO
    

!    select case(this%dynamics)
!    case(OPTION__GAUGEFIELDDYNAMICS__NONE)
!      this%force(1:gr%sb%dim) = M_ZERO 

    !case(OPTION__GAUGEFIELDDYNAMICS__POLARIZATION)
    !  istot = 1
    ! !  if (st%d%nspin > 1) istot = 2
    !   do idir = 1, gr%sb%periodic_dim
    !     this%force(idir) = M_ZERO
    !     do ispin = 1, istot                      
    !       this%force(idir) = this%force(idir) - &
    !                            CNST(4.0)*M_PI*P_c/gr%sb%rcell_volume*dmf_integrate(gr%mesh, st%current(:, idir, ispin))
    !     end do
    !   end do

!    case default
!      ASSERT(.false.)
!    end select

    POP_SUB(thermal_gradient_get_force)
  end subroutine thermal_gradient_get_force


end module thermal_gradient_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
