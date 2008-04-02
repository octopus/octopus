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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: gauge_field.F90 3988 2008-03-31 15:06:50Z fnog $

#include "global.h"

module gauge_field_m
  use datasets_m
  use derivatives_m
  use functions_m
  use global_m
  use grid_m
  use io_m
  use lalg_basic_m
  use loct_parser_m
  use splines_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use simul_box_m
  use units_m
  use logrid_m
  use ps_m
  use specie_m
  use specie_pot_m
  use solids_m
  use geometry_m
  use states_m
  use submesh_m
  use lasers_m
  use profiling_m
  use mpi_m
  use mpi_debug_m
  use varinfo_m
  use projector_m

  implicit none

  private

  public ::                               &
    gauge_field_t,                        &
    gauge_field_init,                     &
    gauge_field_is_applied,               &
    gauge_field_set_vec_pot,     &
    gauge_field_set_vec_pot_vel, &
    gauge_field_get_vec_pot,     &
    gauge_field_get_vec_pot_vel, &
    gauge_field_get_vec_pot_acc, &
    gauge_field_propagate,                &
    gauge_field_propagate_vel,            &
    gauge_field_get_force,                &
    gauge_field_end

  type gauge_field_t
    private
    FLOAT   :: vecpot(1:MAX_DIM)   
    FLOAT   :: vecpot_vel(1:MAX_DIM)
    FLOAT   :: vecpot_acc(1:MAX_DIM)      
    logical :: with_gauge_field    
  end type gauge_field_t

contains

  ! ---------------------------------------------------------
  subroutine gauge_field_init(this, sb)
    type(gauge_field_t),     intent(out)   :: this
    type(simul_box_t),       intent(in)    :: sb

    integer :: ii
    type(block_t) :: blk

    call push_sub('epot.epot_init')

    this%with_gauge_field = .false.
    this%vecpot = M_ZERO
    this%vecpot_vel = M_ZERO
    this%vecpot_acc = M_ZERO
    
    !%Variable GaugeVectorField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The gauge vector field is used to include a uniform (but time dependent)
    !% external electric field in a time dependent run for a periodic system
    !% By default this field is kept null.
    !%End
    ! Read the initial gauge vector field

    if(simul_box_is_periodic(sb)) then
      if(loct_parse_block(check_inp('GaugeVectorField'), blk) == 0) then

        this%with_gauge_field = .true.

    	do ii = 1, sb%dim
          call loct_parse_block_float(blk, 0, ii - 1, this%vecpot(ii))
	end do

	call loct_parse_block_end(blk)

      end if
    end if

    call pop_sub()
  end subroutine gauge_field_init

  ! ---------------------------------------------------------

  subroutine gauge_field_end(this)
    type(gauge_field_t),     intent(inout) :: this

    this%with_gauge_field = .false.
  end subroutine gauge_field_end
  
  logical pure function gauge_field_is_applied(this) result(is_applied)
    type(gauge_field_t),  intent(in) :: this

    is_applied = this%with_gauge_field
  end function gauge_field_is_applied

  subroutine gauge_field_set_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(1:MAX_DIM)

    this%vecpot = vec_pot 
  end subroutine gauge_field_set_vec_pot

  subroutine gauge_field_set_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(1:MAX_DIM)

    this%vecpot_vel = vec_pot_vel
  end subroutine gauge_field_set_vec_pot_vel

  function gauge_field_get_vec_pot(this) result(vec_pot)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot(1:MAX_DIM)

    vec_pot = this%vecpot
  end function gauge_field_get_vec_pot

  function gauge_field_get_vec_pot_vel(this) result(vec_pot_vel)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_vel(1:MAX_DIM)

    vec_pot_vel = this%vecpot_vel
  end function gauge_field_get_vec_pot_vel

  function gauge_field_get_vec_pot_acc(this) result(vec_pot_acc)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_acc(1:MAX_DIM)

    vec_pot_acc = this%vecpot_acc
  end function gauge_field_get_vec_pot_acc

  subroutine gauge_field_propagate(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: force(1:MAX_DIM)
    FLOAT,                intent(in)    :: dt

    this%vecpot_acc = force

    this%vecpot = this%vecpot + dt*this%vecpot_vel + M_HALF*dt**2*force
  end subroutine gauge_field_propagate

  subroutine gauge_field_propagate_vel(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: force(1:MAX_DIM)
    FLOAT,                intent(in)    :: dt

    this%vecpot_vel = this%vecpot_vel + M_HALF*dt*(this%vecpot_acc + force)
  end subroutine gauge_field_propagate_vel

  function gauge_field_get_force(this, gr, st) result(force)
    type(gauge_field_t),  intent(in)    :: this
    type(grid_t),         intent(inout) :: gr
    type(states_t),       intent(inout) :: st
    
    integer :: ik, ist, idir, idim
    CMPLX, allocatable :: gpsi(:, :, :)
    FLOAT :: force(1:MAX_DIM), n_el
#ifdef HAVE_MPI
    FLOAT :: force_tmp(1:MAX_DIM)
#endif

    ALLOCATE(gpsi(gr%m%np, 1:NDIM, st%d%dim), gr%m%np*NDIM*st%d%dim)

    n_el = M_ZERO
    do ik = 1, st%d%spin_channels
      n_el = n_el + dmf_integrate(gr%m, st%rho(1:NP, ik))
    end do

    force(1:MAX_DIM) = -this%vecpot(1:MAX_DIM)*n_el
    
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:, idim, ist, ik), gpsi(:, :, idim))
          
          do idir = 1, gr%sb%dim
            force(idir) = force(idir) + &
                 st%d%kweights(ik)*st%occ(ist, ik)/M_ZI*zstates_dotp(gr%m, st%d%dim, st%zpsi(:, :, ist, ik), gpsi(:, idir, :))
          end do
          
        end do
      end do
    end do
    
#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(force, force_tmp, 3, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      force_tmp = force
    end if
#endif

    deallocate(gpsi)

    force = force*M_FOUR*M_PI*P_c/gr%sb%rcell_volume
    
  end function gauge_field_get_force

end module gauge_field_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
