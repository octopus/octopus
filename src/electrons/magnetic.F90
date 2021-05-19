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

module magnetic_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                    &
    magnetic_density,          &
    magnetic_moment,           &
    write_magnetic_moments,    &
    magnetic_local_moments,    &
    calc_physical_current,     &
    magnetic_induced,          &
    magnetic_total_magnetization, &
    write_total_xc_torque     ,&
    calc_xc_torque

contains

  ! ---------------------------------------------------------
  subroutine magnetic_density(mesh, st, rho, md)
    type(mesh_t),        intent(in)  :: mesh
    type(states_elec_t), intent(in)  :: st
    FLOAT,               intent(in)  :: rho(:,:) !< (np, st%d%nspin)
    FLOAT,               intent(out) :: md(:,:)  !< (np, 3)

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
    type(mesh_t),        intent(in)  :: mesh
    type(states_elec_t), intent(in)  :: st
    FLOAT,               intent(in)  :: rho(:,:)
    FLOAT,               intent(out) :: mm(3)

    FLOAT, allocatable :: md(:,:)

    PUSH_SUB(states_elec_magnetic_moment)

    SAFE_ALLOCATE(md(1:mesh%np, 1:3))
    call magnetic_density(mesh, st, rho, md)

    mm(1) = dmf_integrate(mesh, md(:, 1), reduce = .false.)
    mm(2) = dmf_integrate(mesh, md(:, 2), reduce = .false.)
    mm(3) = dmf_integrate(mesh, md(:, 3), reduce = .false.)

    if(mesh%parallel_in_domains) then
      call mesh%allreduce(mm)
    end if

    SAFE_DEALLOCATE_A(md)

    POP_SUB(states_elec_magnetic_moment)
  end subroutine magnetic_moment


  ! ---------------------------------------------------------
  subroutine write_magnetic_moments(iunit, mesh, st, ions, boundaries, lmm_r)
    integer,             intent(in) :: iunit
    type(mesh_t),        intent(in) :: mesh
    type(states_elec_t), intent(in) :: st
    type(ions_t),        intent(in) :: ions
    type(boundaries_t),  intent(in) :: boundaries
    FLOAT,               intent(in) :: lmm_r
    
    integer :: ia
    FLOAT :: mm(max(mesh%sb%dim, 3))
    FLOAT, allocatable :: lmm(:,:)
    
    PUSH_SUB(write_magnetic_moments)
    
    call magnetic_moment(mesh, st, st%rho, mm)
    SAFE_ALLOCATE(lmm(1:max(mesh%sb%dim, 3), 1:ions%natoms))
    call magnetic_local_moments(mesh, st, ions, boundaries, st%rho, lmm_r, lmm)
    
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
        do ia = 1, ions%natoms
          write(iunit,'(i4,a10,f15.6)') ia, trim(species_label(ions%atom(ia)%species)), lmm(3, ia)
        end do
      else if(st%d%ispin == SPINORS) then ! non-collinear
        write(iunit,'(a,8x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
        do ia = 1, ions%natoms
          write(iunit,'(i4,a10,9f15.6)') ia, trim(species_label(ions%atom(ia)%species)), lmm(:, ia)
        end do
      end if
      
    end if
    
    SAFE_DEALLOCATE_A(lmm)
    
    POP_SUB(write_magnetic_moments)
  end subroutine write_magnetic_moments

  ! ---------------------------------------------------------
  subroutine magnetic_local_moments(mesh, st, ions, boundaries, rho, rr, lmm)
    type(mesh_t),         intent(in)  :: mesh
    type(states_elec_t),  intent(in)  :: st
    type(ions_t),         intent(in)  :: ions
    type(boundaries_t),   intent(in) :: boundaries
    FLOAT,                intent(in)  :: rho(:,:)
    FLOAT,                intent(in)  :: rr
    FLOAT,                intent(out) :: lmm(max(mesh%sb%dim, 3), ions%natoms)

    integer :: ia, idir, is
    FLOAT, allocatable :: md(:, :)
    type(submesh_t) :: sphere
    FLOAT :: cosqr, sinqr
    CMPLX, allocatable :: phase_spiral(:)

    PUSH_SUB(magnetic_local_moments)

    SAFE_ALLOCATE(md (1:mesh%np, 1:max(mesh%sb%dim, 3)))
    
    call magnetic_density(mesh, st, rho, md)
    lmm = M_ZERO
    do ia = 1, ions%natoms
      call submesh_init(sphere, ions%space, mesh, ions%latt, ions%pos(:, ia), rr)

      if(boundaries%spiral) then 
        SAFE_ALLOCATE(phase_spiral(1:sphere%np))
        do is = 1, sphere%np
          phase_spiral(is) = exp(+M_zI*sum((sphere%x(is,:) - mesh%x(sphere%map(is),:)) &
                                       *boundaries%spiral_q(1:mesh%sb%dim)))
        end do

        if(mesh%sb%dim>= 3) then
          lmm(1,ia) = M_ZERO
          lmm(2,ia) = M_ZERO
 
          do is = 1, sphere%np
            !There is a factor of 1/2 in phase_spiral
            cosqr = TOFLOAT(phase_spiral(is))
            sinqr = aimag(phase_spiral(is))
            lmm(1,ia) = lmm(1,ia)+md(sphere%map(is),1)*cosqr - md(sphere%map(is),2)*sinqr
            lmm(2,ia) = lmm(2,ia)+md(sphere%map(is),1)*sinqr + md(sphere%map(is),2)*cosqr
          end do
          lmm(1,ia) = lmm(1,ia)*mesh%volume_element
          lmm(2,ia) = lmm(2,ia)*mesh%volume_element
          lmm(3,ia) = dsm_integrate_frommesh(mesh, sphere, md(1:mesh%np,3), reduce = .false.)
        else
          ASSERT(.not.boundaries%spiral)
        end if

        SAFE_DEALLOCATE_A(phase_spiral)
      else
        do idir = 1, max(mesh%sb%dim, 3)
          lmm(idir, ia) = dsm_integrate_frommesh(mesh, sphere, md(1:mesh%np,idir), reduce = .false.)
        end do
      end if
      
      call submesh_end(sphere) 
    end do

    if(mesh%parallel_in_domains) then
      call mesh%allreduce(lmm)
    end if 
    
    SAFE_DEALLOCATE_A(md)

    POP_SUB(magnetic_local_moments)
  end subroutine magnetic_local_moments

  ! ---------------------------------------------------------
  subroutine magnetic_total_magnetization(mesh, st, qq, trans_mag)
    type(mesh_t),        intent(in)  :: mesh
    type(states_elec_t), intent(in)  :: st
    FLOAT,               intent(in)  :: qq(:)
    CMPLX,               intent(out) :: trans_mag(6)

    integer :: ip
    CMPLX, allocatable :: tmp(:,:)
    FLOAT, allocatable :: md(:, :)
    CMPLX :: expqr
    type(profile_t), save :: prof

    PUSH_SUB(magnetic_total_magnetization)

    call profiling_in(prof, "TOTAL_MAGNETIZATION")

    SAFE_ALLOCATE(tmp(1:mesh%np, 1:6))
    SAFE_ALLOCATE(md (1:mesh%np, 1:max(mesh%sb%dim, 3)))

    call magnetic_density(mesh, st, st%rho, md)
    do ip = 1, mesh%np
      expqr = exp(-M_zI*sum(mesh%x(ip, 1:mesh%sb%dim)*qq(1:mesh%sb%dim)))
      tmp(ip,1) = expqr*md(ip,1)
      tmp(ip,2) = expqr*md(ip,2)
      tmp(ip,3) = expqr*md(ip,3)
      tmp(ip,4) = conjg(expqr)*md(ip,1)
      tmp(ip,5) = conjg(expqr)*md(ip,2)
      tmp(ip,6) = conjg(expqr)*md(ip,3)
    end do

    !TODO: combine the reductions here
    do ip = 1, 6
      trans_mag(ip) = zmf_integrate(mesh, tmp(:,ip))
    end do
     

    SAFE_DEALLOCATE_A(md)
    SAFE_DEALLOCATE_A(tmp)

    call profiling_out(prof)

    POP_SUB(magnetic_total_magnetization)
  end subroutine magnetic_total_magnetization


  ! ---------------------------------------------------------
  !TODO: We should remove this routine and use st%current. NTD
  subroutine calc_physical_current(der, st, kpoints, jj)
    type(derivatives_t),  intent(in)    :: der
    type(states_elec_t),  intent(inout) :: st
    type(kpoints_t),      intent(in)    :: kpoints
    FLOAT,                intent(out)   :: jj(:,:,:)

    PUSH_SUB(calc_physical_current)

    ! Paramagnetic contribution to the physical current
    call states_elec_calc_quantities(der, st, kpoints, .false., paramagnetic_current = jj)

    ! \todo
    ! Diamagnetic contribution to the physical current

    POP_SUB(calc_physical_current)
  end subroutine calc_physical_current


  ! ---------------------------------------------------------
  !> This subroutine receives as input a current, and produces
  !! as an output the vector potential that it induces.
  !! \warning There is probably a problem for 2D. For 1D none of this makes sense?
  subroutine magnetic_induced(der, st, psolver, kpoints, a_ind, b_ind)
    type(derivatives_t),  intent(in)    :: der
    type(states_elec_t),  intent(inout) :: st
    type(poisson_t),      intent(in)    :: psolver
    type(kpoints_t),      intent(in)    :: kpoints
    FLOAT,                intent(out)   :: a_ind(:, :) !< a_ind(der%mesh%np_part, der%dim)
    FLOAT,                intent(out)   :: b_ind(:, :)
    !< if der%dim=3, b_ind(der%mesh%np_part, der%dim)
    !< if der%dim=2, b_ind(der%mesh%np_part, 1)

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

    SAFE_ALLOCATE(jj(1:der%mesh%np_part, 1:der%dim, 1:st%d%nspin))
    call states_elec_calc_quantities(der, st, kpoints, .false., paramagnetic_current = jj)

    !We sum the current for up and down, valid for collinear and noncollinear spins
    if(st%d%nspin > 1) then
      do idir = 1, der%dim
        jj(:, idir, 1) = jj(:, idir, 1) + jj(:, idir, 2)
      end do 
    end if

    a_ind = M_ZERO
    do idir = 1, der%dim
      call dpoisson_solve(psolver, a_ind(:, idir), jj(:, idir, 1))
    end do
    ! This minus sign is introduced here because the current that has been used
    ! before is the "number-current density", and not the "charge-current density",
    ! and therefore there is a minus sign missing (electrons are negative charges...)
    a_ind(1:der%mesh%np, 1:der%dim) = - a_ind(1:der%mesh%np, 1:der%dim) / P_C

    call dderivatives_curl(der, a_ind, b_ind)

    SAFE_DEALLOCATE_A(jj)
    POP_SUB(magnetic_induced)
  end subroutine magnetic_induced

  subroutine write_total_xc_torque(iunit, mesh, vxc, st)
    integer,                  intent(in) :: iunit
    type(mesh_t),             intent(in) :: mesh
    FLOAT,                    intent(in) :: vxc(:,:)
    type(states_elec_t),      intent(in) :: st

    FLOAT, allocatable :: torque(:,:)
    FLOAT :: tt(3)

    PUSH_SUB(write_total_xc_torque)

    SAFE_ALLOCATE(torque(1:mesh%np, 1:3))

    call calc_xc_torque(mesh, vxc, st, torque)

    tt(1) = dmf_integrate(mesh, torque(:,1))
    tt(2) = dmf_integrate(mesh, torque(:,2))
    tt(3) = dmf_integrate(mesh, torque(:,3))
   
    if(mpi_grp_is_root(mpi_world)) then
      write(iunit, '(a)') 'Total xc torque:'
      write(iunit, '(1x,3(a,es10.3,3x))') 'Tx = ', tt(1),'Ty = ', tt(2),'Tz = ', tt(3)
    end if

    SAFE_DEALLOCATE_A(torque)

    POP_SUB(write_total_xc_torque)
  end subroutine write_total_xc_torque

  ! ---------------------------------------------------------
  subroutine calc_xc_torque(mesh, vxc, st, torque)
    type(mesh_t),             intent(in) :: mesh
    FLOAT,                    intent(in) :: vxc(:,:)
    type(states_elec_t),      intent(in) :: st
    FLOAT,                 intent(inout) :: torque(:,:)

    FLOAT :: mag(3), Bxc(3)
    integer :: ip

    PUSH_SUB(calc_xc_torque)

    ASSERT(st%d%ispin == SPINORS)

    do ip = 1, mesh%np
      mag(1) =  M_TWO * st%rho(ip, 3)
      mag(2) = -M_TWO * st%rho(ip, 4)
      mag(3) = st%rho(ip, 1) - st%rho(ip, 2)
      Bxc(1) = -M_TWO * vxc(ip, 3)
      Bxc(2) =  M_TWO * vxc(ip, 4)
      Bxc(3) = -(vxc(ip, 1) - vxc(ip, 2))
      torque(ip, 1) = mag(2) * Bxc(3) - mag(3) * Bxc(2)
      torque(ip, 2) = mag(3) * Bxc(1) - mag(1) * Bxc(3)
      torque(ip, 3) = mag(1) * Bxc(2) - mag(2) * Bxc(1)
    end do

    POP_SUB(calc_xc_torque)
  end subroutine calc_xc_torque



end module magnetic_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
