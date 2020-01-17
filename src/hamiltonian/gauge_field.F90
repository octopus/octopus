!! Copyright (C) 2008 X. Andrade
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

module gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use symmetries_oct_m
  use symm_op_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                               &
    gauge_field_t,                        &
    gauge_field_nullify,                  &
    gauge_field_init,                     &
    gauge_field_init_vec_pot,             &
    gauge_field_is_applied,               &
    gauge_field_set_vec_pot,              &
    gauge_field_set_vec_pot_vel,          &
    gauge_field_get_vec_pot,              &
    gauge_field_get_vec_pot_vel,          &
    gauge_field_get_vec_pot_acc,          &
    gauge_field_propagate,                &
    gauge_field_propagate_vel,            &
    gauge_field_get_energy,               &
    gauge_field_dump,                     &
    gauge_field_load,                     &
    gauge_field_end,                      &
    gauge_field_get_force

  type gauge_field_t
    private
    FLOAT   :: vecpot(1:MAX_DIM)
    FLOAT   :: vecpot_vel(1:MAX_DIM)
    FLOAT   :: vecpot_acc(1:MAX_DIM)
    FLOAT   :: vecpot_kick(1:MAX_DIM)
    FLOAT   :: force(1:MAX_DIM)
    FLOAT   :: wp2
    integer :: ndim
    logical, public :: with_gauge_field
    integer :: dynamics
    FLOAT   :: kicktime
  end type gauge_field_t

contains

  subroutine gauge_field_nullify(this)
    type(gauge_field_t), intent(out) :: this

    PUSH_SUB(gauge_field_nullify)
    this%with_gauge_field = .false.

    POP_SUB(gauge_field_nullify)
  end subroutine gauge_field_nullify

  ! ---------------------------------------------------------
  subroutine gauge_field_init(this, namespace, sb)
    type(gauge_field_t),     intent(out)   :: this
    type(namespace_t),       intent(in)    :: namespace
    type(simul_box_t),       intent(in)    :: sb

    integer :: ii, iop
    type(block_t) :: blk

    PUSH_SUB(gauge_field_init)

    this%with_gauge_field = .false.
    this%vecpot = M_ZERO
    this%vecpot_vel = M_ZERO
    this%vecpot_acc = M_ZERO
    this%vecpot_kick = M_ZERO
    this%force = M_ZERO
    this%ndim = sb%dim

    !%Variable GaugeFieldDynamics
    !%Type integer
    !%Default polarization
    !%Section Hamiltonian
    !%Description
    !% This variable select the dynamics of the gauge field used to
    !% apply a finite electric field to periodic systems in
    !% time-dependent runs.
    !%Option none 0
    !% The gauge field does not have dynamics. The induced polarization field is zero.
    !%Option polarization 1
    !% The gauge field follows the dynamic described in
    !% Bertsch et al, Phys. Rev. B 62 7998 (2000).
    !%End

    call parse_variable(namespace, 'GaugeFieldDynamics', OPTION__GAUGEFIELDDYNAMICS__POLARIZATION, this%dynamics)

    !%Variable GaugeFieldPropagate
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% Propagate the gauge field with initial condition set by GaugeVectorField or zero if not specified
    !%End

    call parse_variable(namespace, 'GaugeFieldPropagate', .false., this%with_gauge_field)

    !%Variable GaugeVectorField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The gauge vector field is used to include a uniform (but time-dependent)
    !% external electric field in a time-dependent run for
    !% a periodic system. An optional second row specifies the initial
    !% value for the time derivative of the gauge field (which is set
    !% to zero by default). By default this field is not included.
    !% If <tt>KPointsUseSymmetries = yes</tt>, then <tt>SymmetryBreakDir</tt>
    !% must be set in the same direction.
    !% This is used with utility <tt>oct-dielectric_function</tt>
    !% according to GF Bertsch, J-I Iwata, A Rubio, and K Yabana,
    !% <i>Phys. Rev. B</i> <b>62</b>, 7998-8002 (2000).
    !%End
    ! Read the initial gauge vector field

    if(parse_block(namespace, 'GaugeVectorField', blk) == 0) then

      this%with_gauge_field = .true.

      do ii = 1, this%ndim
        call parse_block_float(blk, 0, ii - 1, this%vecpot_kick(ii))
      end do

      call parse_block_end(blk)
      if(.not. simul_box_is_periodic(sb)) then
        message(1) = "GaugeVectorField is intended for periodic systems."
        call messages_warning(1, namespace=namespace)
      end if

      if(sb%kpoints%use_symmetries) then
        do iop = 1, symmetries_number(sb%symm)
          if(iop == symmetries_identity_index(sb%symm)) cycle
          if(.not. symm_op_invariant_cart(sb%symm%ops(iop), this%vecpot_kick, CNST(1e-5))) then
            message(1) = "The GaugeVectorField breaks (at least) one of the symmetries used to reduce the k-points."
            message(2) = "Set SymmetryBreakDir equal to GaugeVectorField."
            call messages_fatal(2, namespace=namespace)
          end if
        end do
      end if

    end if

    !%Variable GaugeFieldDelay
    !%Type float
    !%Default 0.
    !%Section Hamiltonian
    !%Description
    !% The application of the gauge field acts as a probe of the system. For dynamical
    !% systems one can apply this probe with a delay relative to the start of the simulation.
    !%End

    call parse_variable(namespace, 'GaugeFieldDelay', M_ZERO, this%kicktime)

    if(abs(this%kicktime) <= M_EPSILON) then
      this%vecpot(1:this%ndim) = this%vecpot_kick(1:this%ndim)
    endif

    POP_SUB(gauge_field_init)
  end subroutine gauge_field_init


  ! ---------------------------------------------------------
  subroutine gauge_field_end(this)
    type(gauge_field_t),     intent(inout) :: this

    PUSH_SUB(gauge_field_end)

    this%with_gauge_field = .false.

    POP_SUB(gauge_field_end)
  end subroutine gauge_field_end


  ! ---------------------------------------------------------
  logical pure function gauge_field_is_applied(this) result(is_applied)
    type(gauge_field_t),  intent(in) :: this

    is_applied = this%with_gauge_field
  end function gauge_field_is_applied


  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(:) !< (this%ndim)

    PUSH_SUB(gauge_field_set_vec_pot)
    this%vecpot(1:this%ndim) = vec_pot(1:this%ndim)

    POP_SUB(gauge_field_set_vec_pot)
  end subroutine gauge_field_set_vec_pot


  ! ---------------------------------------------------------
  subroutine gauge_field_set_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(:) !< (this%ndim)

    PUSH_SUB(gauge_field_set_vec_pot_vel)
    this%vecpot_vel(1:this%ndim) = vec_pot_vel(1:this%ndim)

    POP_SUB(gauge_field_set_vec_pot_vel)
  end subroutine gauge_field_set_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot(:) !< (this%ndim)

    PUSH_SUB(gauge_field_get_vec_pot)
    vec_pot(1:this%ndim) = this%vecpot(1:this%ndim)

    POP_SUB(gauge_field_get_vec_pot)
  end subroutine gauge_field_get_vec_pot


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_vel(:) !< (this%ndim)

    PUSH_SUB(gauge_field_get_vec_pot_vel)
    vec_pot_vel(1:this%ndim) = this%vecpot_vel(1:this%ndim)

    POP_SUB(gauge_field_get_vec_pot_vel)
  end subroutine gauge_field_get_vec_pot_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_get_vec_pot_acc(this, vec_pot_acc)
    type(gauge_field_t),  intent(in)  :: this
    FLOAT,                intent(out) :: vec_pot_acc(:) !< (this%ndim)

    PUSH_SUB(gauge_field_get_vec_pot_acc)
    vec_pot_acc(1:this%ndim) = this%vecpot_acc(1:this%ndim)

    POP_SUB(gauge_field_get_vec_pot_acc)
  end subroutine gauge_field_get_vec_pot_acc

  ! ---------------------------------------------------------
  subroutine gauge_field_propagate(this, dt, time, namespace)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: dt
    FLOAT,                intent(in)    :: time
    type(namespace_t),    intent(in)    :: namespace

    logical, save :: warning_shown = .false.
    integer :: idim

    PUSH_SUB(gauge_field_propagate)

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

        write(message(1),'(a)') 'It seems that the gauge-field might be diverging. You should probably check'
        write(message(2),'(a)') 'the simulation parameters, in particular the number of k-points.'
        call messages_warning(2, namespace=namespace)
      end if
    end do
    POP_SUB(gauge_field_propagate)
  end subroutine gauge_field_propagate

  ! ---------------------------------------------------------
  subroutine gauge_field_propagate_vel(this, dt)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: dt

    PUSH_SUB(gauge_field_propagate_vel)

    this%vecpot_vel(1:this%ndim) = this%vecpot_vel(1:this%ndim) + &
      M_HALF * dt * (this%vecpot_acc(1:this%ndim) + this%force(1:this%ndim))

    POP_SUB(gauge_field_propagate_vel)
  end subroutine gauge_field_propagate_vel


  ! ---------------------------------------------------------
  subroutine gauge_field_init_vec_pot(this, sb, st)
    type(gauge_field_t),  intent(inout) :: this
    type(simul_box_t),    intent(in)    :: sb
    type(states_elec_t),  intent(in)    :: st

    PUSH_SUB(gauge_field_init_vec_pot)

    this%wp2 = M_FOUR*M_PI*st%qtot/sb%rcell_volume

    write (message(1), '(a,f12.6,a)') "Info: Electron-gas plasmon frequency", &
      units_from_atomic(units_out%energy, sqrt(this%wp2)), " ["//trim(units_abbrev(units_out%energy))//"]"
    call messages_info(1)

    POP_SUB(gauge_field_init_vec_pot)
  end subroutine gauge_field_init_vec_pot

  ! ---------------------------------------------------------
  FLOAT function gauge_field_get_energy(this, sb) result(energy)
    type(gauge_field_t),  intent(in)    :: this
    type(simul_box_t),    intent(in)    :: sb

    PUSH_SUB(gauge_field_get_energy)
    energy = sb%rcell_volume / (CNST(8.0) * M_PI * P_c**2) * sum(this%vecpot_vel(1:this%ndim)**2)

    POP_SUB(gauge_field_get_energy)
  end function gauge_field_get_energy


  ! ---------------------------------------------------------
  subroutine gauge_field_dump(restart, gfield, ierr)
    type(restart_t),      intent(in)  :: restart
    type(gauge_field_t),  intent(in)  :: gfield
    integer,              intent(out) :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)

    PUSH_SUB(gauge_field_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(gauge_field_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing gauge field restart."
      call messages_info(1)
    end if

    SAFE_ALLOCATE(vecpot(1:gfield%ndim, 1:2))
    vecpot = M_ZERO
    call gauge_field_get_vec_pot(gfield, vecpot(:, 1))
    call gauge_field_get_vec_pot_vel(gfield, vecpot(:, 2))

    call drestart_write_binary(restart, "gauge_field", 2*gfield%ndim, vecpot, err)
    SAFE_DEALLOCATE_A(vecpot)
    if (err /= 0) ierr = ierr + 1

    if (debug%info) then
      message(1) = "Debug: Writing gauge field restart done."
      call messages_info(1)
    end if

    POP_SUB(gauge_field_dump)
  end subroutine gauge_field_dump


  ! ---------------------------------------------------------
  subroutine gauge_field_load(restart, gfield, ierr)
    type(restart_t),      intent(in)    :: restart
    type(gauge_field_t),  intent(inout) :: gfield
    integer,              intent(out)   :: ierr

    integer :: err
    FLOAT, allocatable :: vecpot(:,:)

    PUSH_SUB(gauge_field_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(gauge_field_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading gauge field restart."
      call messages_info(1)
    end if

    SAFE_ALLOCATE(vecpot(1:gfield%ndim, 1:2))
    call drestart_read_binary(restart, "gauge_field", 2*gfield%ndim, vecpot, err)
    if (err /= 0) ierr = ierr + 1

    call gauge_field_set_vec_pot(gfield, vecpot(:,1))
    call gauge_field_set_vec_pot_vel(gfield, vecpot(:,2))
    SAFE_DEALLOCATE_A(vecpot)

    if (debug%info) then
      message(1) = "Debug: Reading gauge field restart done."
      call messages_info(1)
    end if

    POP_SUB(gauge_field_load)
  end subroutine gauge_field_load

  ! ---------------------------------------------------------

  subroutine gauge_field_get_force(this, gr, st)
    type(gauge_field_t),  intent(inout) :: this
    type(grid_t),         intent(in)    :: gr
    type(states_elec_t),  intent(in)    :: st

    integer :: idir,ispin,istot

    PUSH_SUB(gauge_field_get_force)


    select case(this%dynamics)
    case(OPTION__GAUGEFIELDDYNAMICS__NONE)
      this%force(1:gr%sb%dim) = M_ZERO

    case(OPTION__GAUGEFIELDDYNAMICS__POLARIZATION)
      istot = 1
      if (st%d%nspin > 1) istot = 2
      do idir = 1, gr%sb%periodic_dim
        this%force(idir) = M_ZERO
        do ispin = 1, istot
          this%force(idir) = this%force(idir) - &
            CNST(4.0)*M_PI*P_c/gr%sb%rcell_volume*dmf_integrate(gr%mesh, st%current(:, idir, ispin))
        end do
      end do

    case default
      ASSERT(.false.)
    end select

    POP_SUB(gauge_field_get_force)
  end subroutine gauge_field_get_force


end module gauge_field_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
