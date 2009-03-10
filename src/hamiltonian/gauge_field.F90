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
  use species_m
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
    gauge_force_t,                        &   
    gauge_field_t,                        &
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
    gauge_field_get_force,                &
    gauge_field_get_energy,               &
    gauge_field_apply,                    &
    gauge_field_end

  type gauge_force_t
    private
    FLOAT   :: vecpot(1:MAX_DIM)   
  end type gauge_force_t

  type gauge_field_t
    private
    FLOAT   :: vecpot(1:MAX_DIM)   
    FLOAT   :: vecpot_vel(1:MAX_DIM)
    FLOAT   :: vecpot_acc(1:MAX_DIM)    
    FLOAT   :: wp2
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
    !% The gauge vector field is used to include a uniform (but time
    !% dependent) external electric field in a time dependent run for
    !% a periodic system. An optional second row specifies the initial
    !% value for the time derivative of the gauge_field (that is set
    !% to zero by default).  By default this field is not included.
    !%End
    
    ! Read the initial gauge vector field
    
    if(loct_parse_block(datasets_check('GaugeVectorField'), blk) == 0) then
      
      this%with_gauge_field = .true.
      
      do ii = 1, sb%dim
        call loct_parse_block_float(blk, 0, ii - 1, this%vecpot(ii))
      end do
      
      call loct_parse_block_end(blk)
      
    end if

    call pop_sub()
  end subroutine gauge_field_init

  ! ---------------------------------------------------------

  subroutine gauge_field_end(this)
    type(gauge_field_t),     intent(inout) :: this

    this%with_gauge_field = .false.
  end subroutine gauge_field_end

  ! ---------------------------------------------------------
  
  logical pure function gauge_field_is_applied(this) result(is_applied)
    type(gauge_field_t),  intent(in) :: this

    is_applied = this%with_gauge_field
  end function gauge_field_is_applied

  ! ---------------------------------------------------------

  subroutine gauge_field_set_vec_pot(this, vec_pot)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot(1:MAX_DIM)

    this%vecpot = vec_pot 
  end subroutine gauge_field_set_vec_pot

  ! ---------------------------------------------------------

  subroutine gauge_field_set_vec_pot_vel(this, vec_pot_vel)
    type(gauge_field_t),  intent(inout) :: this
    FLOAT,                intent(in)    :: vec_pot_vel(1:MAX_DIM)

    this%vecpot_vel = vec_pot_vel
  end subroutine gauge_field_set_vec_pot_vel

  ! ---------------------------------------------------------

  function gauge_field_get_vec_pot(this) result(vec_pot)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot(1:MAX_DIM)

    vec_pot = this%vecpot
  end function gauge_field_get_vec_pot

  ! ---------------------------------------------------------

  function gauge_field_get_vec_pot_vel(this) result(vec_pot_vel)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_vel(1:MAX_DIM)

    vec_pot_vel = this%vecpot_vel
  end function gauge_field_get_vec_pot_vel

  ! ---------------------------------------------------------

  function gauge_field_get_vec_pot_acc(this) result(vec_pot_acc)
    type(gauge_field_t),  intent(in) :: this
    FLOAT :: vec_pot_acc(1:MAX_DIM)

    vec_pot_acc = this%vecpot_acc
  end function gauge_field_get_vec_pot_acc

  ! ---------------------------------------------------------

  subroutine gauge_field_propagate(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    type(gauge_force_t),  intent(in)    :: force
    FLOAT,                intent(in)    :: dt

    this%vecpot_acc(1:MAX_DIM) = force%vecpot(1:MAX_DIM)

    this%vecpot = this%vecpot + dt*this%vecpot_vel + M_HALF*dt**2*force%vecpot
  end subroutine gauge_field_propagate

  ! ---------------------------------------------------------

  subroutine gauge_field_propagate_vel(this, force, dt)
    type(gauge_field_t),  intent(inout) :: this
    type(gauge_force_t),  intent(in)    :: force
    FLOAT,                intent(in)    :: dt

    this%vecpot_vel = this%vecpot_vel + M_HALF*dt*(this%vecpot_acc + force%vecpot)
  end subroutine gauge_field_propagate_vel

  ! ---------------------------------------------------------

  subroutine gauge_field_init_vec_pot(this, m, sb, st, dt)
    type(gauge_field_t),  intent(inout) :: this
    type(mesh_t),         intent(in)    :: m
    type(simul_box_t),    intent(in)    :: sb
    type(states_t),       intent(in)    :: st
    FLOAT,                intent(in)    :: dt
    
    integer :: ik
    FLOAT :: n_el

    n_el = M_ZERO
    do ik = 1, st%d%spin_channels
      n_el = n_el + dmf_integrate(m, st%rho(1:m%np, ik))
    end do
    
    this%wp2 = M_FOUR*M_PI*n_el/sb%rcell_volume

    write (message(1), '(a,f12.6,a)') "Info: Electron gas plasmon frequency", &
         sqrt(this%wp2)*units_inp%energy%factor, " ["//trim(units_inp%energy%abbrev)//"]"
    call write_info(1)

  end subroutine gauge_field_init_vec_pot

  ! ---------------------------------------------------------

  subroutine gauge_field_get_force(this, gr, geo, pj, phases, st, force)
    type(gauge_field_t),  intent(inout) :: this
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(projector_t),    intent(in)    :: pj(:)
    CMPLX,                intent(in)    :: phases(:, :)
    type(states_t),       intent(inout) :: st
    type(gauge_force_t),  intent(out)   :: force

    integer :: ik, ist, idir, idim, iatom
    CMPLX, allocatable :: gpsi(:, :, :), cpsi(:, :), epsi(:, :)
#ifdef HAVE_MPI
    FLOAT :: force_tmp(1:MAX_DIM)
#endif

    ALLOCATE(epsi(gr%mesh%np_part, st%d%dim), gr%mesh%np_part*st%d%dim)
    ALLOCATE(gpsi(gr%mesh%np, 1:NDIM, st%d%dim), gr%mesh%np*NDIM*st%d%dim)
    ALLOCATE(cpsi(gr%mesh%np, st%d%dim), gr%mesh%np*st%d%dim)

    force%vecpot(1:MAX_DIM) = M_ZERO
    
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end

        do idim = 1, st%d%dim

          call zset_bc(gr%der, st%zpsi(:, idim, ist, ik))

          if(simul_box_is_periodic(gr%sb)) then
            epsi(1:NP_PART, idim) = phases(1:NP_PART, ik)*st%zpsi(1:NP_PART, idim, ist, ik)
          end if

          call zderivatives_grad(gr%der, epsi(:, idim), gpsi(:, :, idim), set_bc = .false.)

        end do

        do idir = 1, gr%sb%dim
          do iatom = 1, geo%natoms
            if(species_is_ps(geo%atom(iatom)%spec)) then
              call zprojector_commute_r(pj(iatom), gr, st%d%dim, idir, ik, epsi(:, 1), cpsi(:, :))
              gpsi(1:NP, idir, 1:st%d%dim) = gpsi(1:NP, idir, 1:st%d%dim) + cpsi(1:NP, 1:st%d%dim)
            end if
          end do
        end do
        
        do idir = 1, gr%sb%dim
          force%vecpot(idir) = force%vecpot(idir) + M_FOUR*M_PI*P_c/gr%sb%rcell_volume*st%d%kweights(ik)*st%occ(ist, ik)*&
               aimag(zmf_dotp(gr%mesh, st%d%dim, epsi, gpsi(:, idir, :)))
        end do
        
      end do
    end do
    
#ifdef HAVE_MPI
    if(st%parallel_in_states) then
      call MPI_Allreduce(force%vecpot, force_tmp, MAX_DIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
      force%vecpot = force_tmp
    end if
#endif

    force%vecpot(1:MAX_DIM) = force%vecpot(1:MAX_DIM) - this%wp2*this%vecpot(1:MAX_DIM)

    deallocate(gpsi, cpsi)

  end subroutine gauge_field_get_force

  ! ---------------------------------------------------------

  FLOAT function gauge_field_get_energy(this, sb) result(energy)
    type(gauge_field_t),  intent(in)    :: this
    type(simul_box_t),    intent(in)    :: sb

    energy = sb%rcell_volume/(M_EIGHT*M_PI*P_c**2)*sum(this%vecpot_vel(1:MAX_DIM)**2)

  end function gauge_field_get_energy

  ! ---------------------------------------------------------
  subroutine gauge_field_apply(this, gr, dim, psi, grad, hpsi)
    type(gauge_field_t), intent(in)    :: this
    type(grid_t),        intent(inout) :: gr
    integer,             intent(in)    :: dim
    CMPLX,               intent(in)    :: psi(:,:)  ! psi(NP_PART, hm%d%dim)
    CMPLX,               intent(in)    :: grad(:, :, :)
    CMPLX,               intent(inout) :: hpsi(:,:)
    
    integer :: ip, idim, a2
    FLOAT :: vecpot(1:MAX_DIM)
    
    call push_sub('hamiltonian_inc.Xvgauge')
    
    ASSERT(gauge_field_is_applied(this))
    
    vecpot = gauge_field_get_vec_pot(this)/P_c
    a2 = sum(vecpot(1:MAX_DIM)**2)
    
    forall(idim = 1:dim, ip = 1:NP)
      hpsi(ip, idim) = hpsi(ip, idim) + M_HALF*a2*psi(ip, idim) + M_zI*dot_product(vecpot(1:MAX_DIM), grad(ip, 1:MAX_DIM, idim))
    end forall
    
    call pop_sub()
  end subroutine gauge_field_apply
  
end module gauge_field_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
