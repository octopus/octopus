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
  use global_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                      &
    photon_mode_t,               &
    photon_mode_init,            &
    photon_mode_end

  type photon_mode_t
    ! All components are public by default
    integer               :: nmodes             !< Number of photon modes
    integer               :: dim                !< Dimensionality of the electronic system
    FLOAT, allocatable    :: omega(:)           !< Mode frequencies
    FLOAT, allocatable    :: lambda(:)          !< Interaction strength
    FLOAT, allocatable    :: pol(:,:)           !< Polarization of the photon field
    FLOAT, allocatable    :: pol_dipole(:,:)    !< Polarization*dipole operator
    FLOAT                 :: ex                 !< Photon exchange energy
    FLOAT, allocatable    :: number(:)          !< Number of photons in mode
    FLOAT, allocatable    :: correlator(:,:)    !< Correlation function <n(r)(ad+a)>
    FLOAT                 :: n_electrons        !< Number of electrons
  end type photon_mode_t

contains

  ! ---------------------------------------------------------
  subroutine photon_mode_init(this, namespace, mesh, dim, n_electrons)
    type(photon_mode_t),  intent(out) :: this
    type(namespace_t),    intent(in)  :: namespace
    type(mesh_t),         intent(in)  :: mesh
    integer,              intent(in)  :: dim
    FLOAT,                intent(in)  :: n_electrons

    type(block_t)         :: blk
    integer               :: ii, idir, ncols

    PUSH_SUB(photon_mode_init)

    this%dim = dim
    this%n_electrons = n_electrons

    !%Variable PhotonModes
    !%Type block
    !%Section Hamiltonian::XC
    !%Description
    !% Each line of the block should specify one photon mode. The syntax is the following:
    !%
    !%PhotonModes
    !%omega1 | lambda1| PolX1 | PolY1 | PolZ1
    !%...
    !%
    !%End
    !%
    !% The first column is the mode frequency, in units of energy.
    !% The second column is the coupling strength, in units of energy.
    !% The remaining columns specify the polarization direction of the mode.
    !% If the polarization vector should be normalized to one. If that is not the case
    !% the code will normalize it.
    this%nmodes = 0
    if(parse_block(namespace, 'PhotonModes', blk) == 0) then

       this%nmodes = parse_block_n(blk)
       SAFE_ALLOCATE(this%omega(1:this%nmodes))
       SAFE_ALLOCATE(this%lambda(1:this%nmodes))
       SAFE_ALLOCATE(this%pol(1:this%nmodes, this%dim))
       SAFE_ALLOCATE(this%pol_dipole(1:mesh%np, 1:this%nmodes))

       do ii = 1, this%nmodes
          ncols = parse_block_cols(blk, ii-1)

          ! Sanity check
          if (ncols /= this%dim + 2) then
            call messages_input_error(namespace, 'PhotonModes', 'Incorrect number of columns')
          end if

          ! Read line
          call parse_block_float(blk, ii-1, 0, this%omega(ii), units_inp%energy)  ! frequency
          call parse_block_float(blk, ii-1, 1, this%lambda(ii), units_inp%energy) ! coupling strength
          do idir = 1, this%dim
            call parse_block_float(blk, 0, idir + 1, this%pol(ii, idir)) ! polarization vector components
          end do

          ! Normalize polarization vector
          this%pol = this%pol/sqrt(dot_product(this%pol(ii,:), this%pol(ii,:)))

          ! Calculate polarization dipole
          this%pol_dipole(:, ii) = M_ZERO
          do idir = 1, this%dim
            this%pol_dipole(1:mesh%np, ii) = this%pol_dipole(1:mesh%np, ii) + this%pol(ii, idir)*mesh%x(1:mesh%np, idir)
          end do

        end do
      call parse_block_end(blk)
    else
      call messages_write('You need to specify the photon modes!')
      call messages_fatal(namespace=namespace)
    end if

    this%ex = M_ZERO
    SAFE_ALLOCATE(this%number(1:this%nmodes))
    this%number = M_ZERO

    SAFE_ALLOCATE(this%correlator(1:mesh%np, 1:this%nmodes))
    this%correlator = M_ZERO

    POP_SUB(photon_mode_init)
  end subroutine photon_mode_init

  ! ---------------------------------------------------------
  subroutine photon_mode_end(this)
    type(photon_mode_t), intent(inout) :: this

    PUSH_SUB(photon_mode_end)

    SAFE_DEALLOCATE_A(this%correlator)

    SAFE_DEALLOCATE_A(this%omega)
    SAFE_DEALLOCATE_A(this%lambda)
    SAFE_DEALLOCATE_A(this%number)

    SAFE_DEALLOCATE_A(this%pol)
    SAFE_DEALLOCATE_A(this%pol_dipole)

    POP_SUB(photon_mode_end)
  end subroutine photon_mode_end


end module photon_mode_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
