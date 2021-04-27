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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module kdotp_calc_oct_m
  use comm_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use linear_response_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use pert_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use utils_oct_m

  implicit none

  private
  public ::                        &
    dcalc_eff_mass_inv,            &
    zcalc_eff_mass_inv,            &
    zcalc_band_velocity,           &
    dkdotp_add_occ,                &
    zkdotp_add_occ,                &
    dkdotp_add_diagonal,           &
    zkdotp_add_diagonal,           &
    kdotp_wfs_tag

contains

! ---------------------------------------------------------
  character(len=100) function kdotp_wfs_tag(dir, dir2) result(str)
    integer,           intent(in) :: dir 
    integer, optional, intent(in) :: dir2

    PUSH_SUB(kdotp_wfs_tag)

    str = "wfs_" // index2axis(dir)
    if(present(dir2)) str = trim(str) // "_" // index2axis(dir2)

    POP_SUB(kdotp_wfs_tag)

  end function kdotp_wfs_tag

! ---------------------------------------------------------
!> v = (dE_nk/dk)/hbar = -Im < u_nk | -i grad | u_nk >
!! This is identically zero for real wavefunctions.
subroutine zcalc_band_velocity(namespace, space, gr, st, hm, ions, pert, velocity)
  type(namespace_t),        intent(in)    :: namespace
  type(space_t),            intent(in)    :: space
  type(grid_t),             intent(in)    :: gr
  type(states_elec_t),      intent(in)    :: st
  type(hamiltonian_elec_t), intent(inout) :: hm
  type(ions_t),             intent(in)    :: ions
  type(pert_t),             intent(inout) :: pert
  FLOAT,                    intent(out)   :: velocity(:,:,:)

  integer :: ik, ist, idir
  CMPLX, allocatable :: psi(:, :), pertpsi(:,:)
  type(profile_t), save :: prof

  PUSH_SUB(zkdotp_calc_band_velocity)

  call profiling_in(prof, "CALC_BAND_VELOCITY")

  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(pertpsi(1:gr%mesh%np, 1:st%d%dim))

  velocity(:, :, :) = M_ZERO

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call states_elec_get_state(st, gr%mesh, ist, ik, psi)

      do idir = 1, space%periodic_dim
        call pert_setup_dir(pert, idir)
        call zpert_apply(pert, namespace, gr, ions, hm, ik, psi, pertpsi)
        velocity(idir, ist, ik) = -aimag(zmf_dotp(gr%mesh, st%d%dim, psi, pertpsi))
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(pertpsi)

  call comm_allreduce(st%st_kpt_mpi_grp, velocity)

  call profiling_out(prof)

  POP_SUB(zkdotp_calc_band_velocity)
end subroutine zcalc_band_velocity

#include "undef.F90"
#include "real.F90"
#include "kdotp_calc_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "kdotp_calc_inc.F90"

end module kdotp_calc_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
