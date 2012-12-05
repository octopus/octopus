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
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module magnetic_m
  use derivatives_m
  use geometry_m
  use global_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use poisson_m
  use profiling_m
  use species_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                    &
    magnetic_density,          &
    magnetic_moment,           &
    write_magnetic_moments,    &
    magnetic_local_moments,    &
    calc_physical_current,     &
    magnetic_induced

contains

  ! ---------------------------------------------------------
  subroutine magnetic_density(mesh, st, rho, md)
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st
    FLOAT,          intent(in)  :: rho(:,:) !< (np, st%d%nspin)
    FLOAT,          intent(out) :: md(:,:)  !< (np, 3)

    PUSH_SUB(magnetic_density)

    select case (st%d%ispin)
    case (UNPOLARIZED)
      md = M_ZERO

    case (SPIN_POLARIZED)
      md = M_ZERO
      md(1:mesh%np, 3) = rho(1:mesh%np, 1) - rho(1:mesh%np, 2)

    case (SPINORS)
      md(1:mesh%np, 1) =  M_TWO*rho(1:mesh%np, 3)
      md(1:mesh%np, 2) = -M_TWO*rho(1:mesh%np, 4)
      md(1:mesh%np, 3) = rho(1:mesh%np, 1) - rho(1:mesh%np, 2)
    end select

    POP_SUB(magnetic_density)
  end subroutine magnetic_density


  ! ---------------------------------------------------------
  subroutine magnetic_moment(mesh, st, rho, mm)
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st
    FLOAT,          intent(in)  :: rho(:,:)
    FLOAT,          intent(out) :: mm(3)

    FLOAT, allocatable :: md(:,:)

    PUSH_SUB(states_magnetic_moment)

    SAFE_ALLOCATE(md(1:mesh%np, 1:MAX_DIM))
    call magnetic_density(mesh, st, rho, md)

    mm(1) = dmf_integrate(mesh, md(:, 1))
    mm(2) = dmf_integrate(mesh, md(:, 2))
    mm(3) = dmf_integrate(mesh, md(:, 3))

    SAFE_DEALLOCATE_A(md)

    POP_SUB(states_magnetic_moment)
  end subroutine magnetic_moment


  ! ---------------------------------------------------------
  subroutine write_magnetic_moments(iunit, mesh, st, geo, lmm_r)
    integer,          intent(in) :: iunit
    type(mesh_t),     intent(in) :: mesh
    type(states_t),   intent(in) :: st
    type(geometry_t), intent(in) :: geo
    FLOAT,            intent(in) :: lmm_r
    
    integer :: ia
    FLOAT :: mm(max(mesh%sb%dim, 3))
    FLOAT, allocatable :: lmm(:,:)
    
    PUSH_SUB(write_magnetic_moments)
    
    call magnetic_moment(mesh, st, st%rho, mm)
    SAFE_ALLOCATE(lmm(1:max(mesh%sb%dim, 3), 1:geo%natoms))
    call magnetic_local_moments(mesh, st, geo, st%rho, lmm_r, lmm)
    
    if(mpi_grp_is_root(mpi_world)) then
      
      write(iunit, '(a)') 'Total Magnetic Moment:'
      if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
        write(iunit, '(a,f10.6)') ' mz = ', mm(3)
      else if(st%d%ispin == SPINORS) then ! non-collinear
        write(iunit, '(1x,3(a,f10.6,3x))') 'mx = ', mm(1),'my = ', mm(2),'mz = ', mm(3)
      end if
      
      write(iunit, '(a,a,a,f7.3,a)') 'Local Magnetic Moments (sphere radius [', &
        trim(units_abbrev(units_out%length)),'] = ', units_from_atomic(units_out%length, lmm_r), '):'
      if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
        write(iunit,'(a,6x,14x,a)') ' Ion','mz'
        do ia = 1, geo%natoms
          write(iunit,'(i4,a10,f15.6)') ia, trim(species_label(geo%atom(ia)%spec)), lmm(3, ia)
        end do
      else if(st%d%ispin == SPINORS) then ! non-collinear
        write(iunit,'(a,8x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
        do ia = 1, geo%natoms
          write(iunit,'(i4,a10,9f15.6)') ia, trim(species_label(geo%atom(ia)%spec)), lmm(1:mesh%sb%dim, ia)
        end do
      end if
      
    end if
    
    SAFE_DEALLOCATE_A(lmm)
    
    POP_SUB(write_magnetic_moments)
  end subroutine write_magnetic_moments

  ! ---------------------------------------------------------
  subroutine magnetic_local_moments(mesh, st, geo, rho, rr, lmm)
    type(mesh_t),     intent(in)  :: mesh
    type(states_t),   intent(in)  :: st
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(in)  :: rho(:,:)
    FLOAT,            intent(in)  :: rr
    FLOAT,            intent(out) :: lmm(max(mesh%sb%dim, 3), geo%natoms)

    integer :: ia, ip, idir
    FLOAT :: ri
    FLOAT, allocatable :: md(:, :), aux(:, :)

    PUSH_SUB(magnetic_local_moments)

    SAFE_ALLOCATE(md (1:mesh%np, 1:max(mesh%sb%dim, 3)))
    SAFE_ALLOCATE(aux(1:mesh%np, 1:max(mesh%sb%dim, 3)))

    call magnetic_density(mesh, st, rho, md)
    lmm = M_ZERO
    do ia = 1, geo%natoms
      aux = M_ZERO
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, ri, origin = geo%atom(ia)%x)
        if (ri > rr) cycle
        aux(ip, 1:max(mesh%sb%dim, 3)) = md(ip, 1:max(mesh%sb%dim, 3))
      end do

      do idir = 1, max(mesh%sb%dim, 3)
        lmm(idir, ia) = dmf_integrate(mesh, aux(1:mesh%np, idir))
      enddo

    end do

    SAFE_DEALLOCATE_A(md)
    SAFE_DEALLOCATE_A(aux)

    POP_SUB(magnetic_local_moments)
  end subroutine magnetic_local_moments


  ! ---------------------------------------------------------
  subroutine calc_physical_current(der, st, jj)
    type(derivatives_t),  intent(in)    :: der
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)   :: jj(:,:,:)

    PUSH_SUB(calc_physical_current)

    ! Paramagnetic contribution to the physical current
    call states_calc_quantities(der, st, paramagnetic_current = jj)

    ! \todo
    ! Diamagnetic contribution to the physical current

    POP_SUB(calc_physical_current)
  end subroutine calc_physical_current


  ! ---------------------------------------------------------
  !> This subroutine receives as input a current, and produces
  !! as an output the vector potential that it induces.
  !! \warning There is probably a problem for 2D. For 1D none of this makes sense?
  subroutine magnetic_induced(der, st, a_ind, b_ind)
    type(derivatives_t),  intent(in)    :: der
    type(states_t),       intent(inout) :: st
    FLOAT,                intent(out)   :: a_ind(:, :) ! a_ind(gr%mesh%np_part, gr%mesh%sb%dim)
    FLOAT,                intent(out)   :: b_ind(:, :)
    ! if gr%mesh%sb%dim=3, b_ind(gr%mesh%np_part, gr%mesh%sb%dim)
    ! if gr%mesh%sb%dim=2, b_ind(gr%mesh%np_part, 1)

    integer :: idir
    FLOAT, allocatable :: jj(:, :, :)

    PUSH_SUB(magnetic_induced)

    ! If the states are real, we should never have reached here, but
    ! just in case we return zero.
    if(states_are_real(st)) then
      a_ind = M_ZERO
      b_ind = M_ZERO
      POP_SUB(magnetic_induced)
      return
    end if

    SAFE_ALLOCATE(jj(1:der%mesh%np_part, 1:der%mesh%sb%dim, 1:st%d%nspin))
    call states_calc_quantities(der, st, paramagnetic_current = jj)

    a_ind = M_ZERO
    do idir = 1, der%mesh%sb%dim
      call dpoisson_solve(psolver, a_ind(:, idir), jj(:, idir, 1))
    end do
    ! This minus sign is introduced here because the current that has been used
    ! before is the "number-current density", and not the "charge-current density",
    ! and therefore there is a minus sign missing (electrons are negative charges...)
    a_ind = - a_ind / P_C

    call dderivatives_curl(der, a_ind, b_ind)

    SAFE_DEALLOCATE_A(jj)
    POP_SUB(magnetic_induced)
  end subroutine magnetic_induced


end module magnetic_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
