!! Copyright (C) 2017 Johannes Flick
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

module photon_mode_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use lalg_adv_oct_m
  use linear_response_oct_m
  use linear_solver_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use scf_tol_oct_m
  use varinfo_oct_m
  use xc_oct_m
  use XC_F90(lib_m)
  use xc_functl_oct_m

  implicit none

  private
  public ::                      &
    photon_mode_t,               &
    photon_mode_init,            &
    photon_mode_end

  type photon_mode_t
    ! All components are public by default
    integer               :: nmodes
    FLOAT, allocatable    :: omega_array(:), lambda_array(:)  ! frequencies and interaction strength
    FLOAT, allocatable    :: pol_array(:,:)                   ! polarization of the photon field
    FLOAT, allocatable    :: pol_dipole_array(:,:)            ! polarization*dipole operator
    FLOAT                 :: ex                               ! photon exchange energy
    FLOAT                 :: pt_number            ! number of photons in mode
    FLOAT, allocatable    :: correlator(:,:)      ! correlation function <n(r)(ad+a)>
    type(lr_t)            :: lr                   !< to solve the equation H psi = b
  end type photon_mode_t

contains

  ! ---------------------------------------------------------
  subroutine photon_mode_init(this, namespace, gr)
    type(photon_mode_t),  intent(out)   :: this
    type(namespace_t),    intent(in)    :: namespace 
    type(grid_t),         intent(inout) :: gr

    type(block_t)         :: blk
    integer               :: ii, ncols

    PUSH_SUB(photon_mode_init)

    !%Variable PhotonModes
    !%Type block
    !%Section Hamiltonian::XC
    !%Description
    !% Syntax:
    !%PhotonModes
    !%omega1 | lambda1| PolX1 | PolY1 | PolZ1
    !%...
    !%
    !%End
    !% frequency of the cavity mode
    !% coupling strength
    !% polarization of the cavity mode in (x,y,z) direction

    this%nmodes = 0
    if(parse_block(namespace, 'PhotonModes', blk) == 0) then

       this%nmodes = parse_block_n(blk)
       SAFE_ALLOCATE(this%omega_array(1:this%nmodes))
       SAFE_ALLOCATE(this%lambda_array(1:this%nmodes))
       SAFE_ALLOCATE(this%pol_array(1:this%nmodes,3))
       SAFE_ALLOCATE(this%pol_dipole_array(1:gr%mesh%np,1:this%nmodes))
       do ii = 1, this%nmodes
          ncols = parse_block_cols(blk, ii-1)
          call parse_block_float(blk, ii-1, 0, this%omega_array(ii))   !row, column
          call parse_block_float(blk, ii-1, 1, this%lambda_array(ii))  !row, column
          call parse_block_float(blk, ii-1, 2, this%pol_array(ii,1))   !row, column
          call parse_block_float(blk, ii-1, 3, this%pol_array(ii,2))   !row, column
          call parse_block_float(blk, ii-1, 4, this%pol_array(ii,3))   !row, column
          this%pol_dipole_array(1:gr%mesh%np,ii) = this%pol_array(ii,1)*gr%mesh%x(1:gr%mesh%np, 1) &
                                        + this%pol_array(ii,2)*gr%mesh%x(1:gr%mesh%np, 2) &
                                        + this%pol_array(ii,3)*gr%mesh%x(1:gr%mesh%np, 3)
       end do
      call parse_block_end(blk)
    else
      call messages_write('You need to specify the photon modes!')
      call messages_fatal()
    end if

    this%ex = M_ZERO
    this%pt_number = M_ZERO

    SAFE_ALLOCATE(this%correlator(1:gr%mesh%np,1))
    this%correlator = M_ZERO

    POP_SUB(photon_mode_init)
  end subroutine photon_mode_init

  ! ---------------------------------------------------------

  subroutine photon_mode_end(this)
    type(photon_mode_t), intent(inout) :: this

    PUSH_SUB(photon_mode_end)

    SAFE_DEALLOCATE_A(this%correlator)

    SAFE_DEALLOCATE_A(this%omega_array)
    SAFE_DEALLOCATE_A(this%lambda_array)

    SAFE_DEALLOCATE_A(this%pol_array)
    SAFE_DEALLOCATE_A(this%pol_dipole_array)

    POP_SUB(photon_mode_end)
  end subroutine photon_mode_end


end module photon_mode_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
