!!! Copyright (C) 2008-2010 David Strubbe
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
!! $Id: kdotp_calc.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module kdotp_calc_m
  use global_m
  use hamiltonian_m
  use linear_response_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use pert_m
  use profiling_m
  use states_m
  use states_calc_m
  use sternheimer_m
  use system_m
  use utils_m

  implicit none

  private
  public ::                        &
    dcalc_eff_mass_inv,            &
    zcalc_eff_mass_inv,            &
    zcalc_band_velocity,           &
    zcalc_dipole_periodic,         &
    dkdotp_add_diagonal,           &
    zkdotp_add_diagonal,           &
    kdotp_wfs_tag

contains

! ---------------------------------------------------------
  character(len=100) function kdotp_wfs_tag(dir, dir2) result(str)
    integer,           intent(in) :: dir 
    integer, optional, intent(in) :: dir2

    PUSH_SUB(kdotp_wfs_tag)

    write(str, '(2a)') "wfs_", index2axis(dir)
    if(present(dir2)) write(str, '(3a)') trim(str), "_", index2axis(dir2)

    POP_SUB(kdotp_wfs_tag)

  end function kdotp_wfs_tag

! ---------------------------------------------------------
! v = (dE_nk/dk)/hbar = -Im < u_nk | -i grad | u_nk >
! This is identically zero for real wavefunctions.
subroutine zcalc_band_velocity(sys, hm, pert, velocity)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: hm
  type(pert_t),        intent(inout) :: pert
  FLOAT,               intent(out)   :: velocity(:,:,:)

  integer :: ik, ist, idir
  CMPLX, allocatable :: psi(:, :), pertpsi(:,:)
#ifdef HAVE_MPI
  FLOAT, allocatable :: vel_temp(:,:,:)
#endif

  PUSH_SUB(zkdotp_calc_band_velocity)

  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  SAFE_ALLOCATE(pertpsi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  velocity(:, :, :) = M_ZERO

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    do ist = sys%st%st_start, sys%st%st_end

      call states_get_state(sys%st, sys%gr%mesh, ist, ik, psi)

      do idir = 1, sys%gr%sb%periodic_dim
        call pert_setup_dir(pert, idir)
        call zpert_apply(pert, sys%gr, sys%geo, hm, ik, psi, pertpsi)
        velocity(ik, ist, idir) = -aimag(zmf_dotp(sys%gr%mesh, sys%st%d%dim, psi, pertpsi))
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(pertpsi)

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states .or. sys%st%d%kpt%parallel) then
    SAFE_ALLOCATE(vel_temp(1:sys%st%d%nik, 1:sys%st%nst, 1:sys%gr%sb%periodic_dim))

    call MPI_Allreduce(velocity, vel_temp, sys%st%d%nik*sys%st%nst*sys%gr%sb%periodic_dim, &
      MPI_FLOAT, MPI_SUM, sys%st%st_kpt_mpi_grp%comm, mpi_err)

    velocity(:,:,:) = vel_temp(:,:,:)
    SAFE_DEALLOCATE_A(vel_temp)
  endif
#endif

  POP_SUB(zkdotp_calc_band_velocity)
end subroutine zcalc_band_velocity

! ---------------------------------------------------------
! This routine cannot be used with d/dk wavefunctions calculated
! by perturbation theory.
subroutine zcalc_dipole_periodic(sys, lr, dipole)
  type(system_t),         intent(inout) :: sys
  type(lr_t),             intent(in)    :: lr(:,:)
  FLOAT,                  intent(out)   :: dipole(:)

  integer idir, ist, ik, idim
  type(mesh_t), pointer :: mesh
  CMPLX :: term, moment
  mesh => sys%gr%mesh

  ! mu_i = sum(m occ, k) <u_mk(0)|(-id/dk_i|u_mk(0)>)
  !      = Im sum(m occ, k) <u_mk(0)|(d/dk_i|u_mk(0)>)

  PUSH_SUB(zcalc_dipole_periodic)

  do idir = 1, sys%gr%sb%periodic_dim
    moment = M_ZERO

    do ik = 1, sys%st%d%nik
      term = M_ZERO

      do ist = 1, sys%st%nst
        do idim = 1, sys%st%d%dim
          term = term + zmf_dotp(mesh, sys%st%zpsi(1:mesh%np, idim, ist, ik), &
            lr(idir, 1)%zdl_psi(1:mesh%np, idim, ist, ik))
        enddo
      enddo

      moment = moment + term * sys%st%d%kweights(ik) * sys%st%smear%el_per_state
    enddo

    dipole(idir) = -moment
  enddo

  POP_SUB(zcalc_dipole_periodic)

end subroutine zcalc_dipole_periodic

#include "undef.F90"
#include "real.F90"
#include "kdotp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "kdotp_calc_inc.F90"

end module kdotp_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
