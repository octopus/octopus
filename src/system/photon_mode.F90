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
  use io_oct_m
  use lalg_adv_oct_m
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
    FLOAT, allocatable    :: pt_number(:)                     ! number of photons in mode
    FLOAT, allocatable    :: correlator(:,:)                  ! correlation function <n(r)(ad+a)>
  end type photon_mode_t

contains

  ! ---------------------------------------------------------
  subroutine photon_mode_init(this, namespace, gr)
    type(photon_mode_t),  intent(out)   :: this
    type(namespace_t),    intent(in)    :: namespace 
    type(grid_t),         intent(inout) :: gr

    type(block_t)         :: blk
    integer               :: ii, ncols, iunit, ierr
    character(MAX_PATH_LEN) :: filename
    logical               :: file_exists

    PUSH_SUB(photon_mode_init)

    this%nmodes = 0

    !%Variable PhotonmodesFilename
    !%Type string
    !%Default "photonmodes"
    !%Section Linear Response::Casida
    !%Description
    !% Filename for photon modes in text format
    !%  - first line contains 2 integers: number of photon modes and number of
    !%    columns
    !%  - each further line contains the given number of floats for one photon
    !%    mode
    !%End
    call parse_variable(namespace, 'PhotonmodesFilename', 'photonmodes', filename)
    inquire(file=trim(filename), exist=file_exists)
    if(file_exists) then
      if(mpi_grp_is_root(mpi_world)) then
        message(1) = 'Opening '//trim(filename)
        call messages_info(1)
        ! open file on root
        iunit = io_open(trim(filename), namespace, action='read', form='formatted')

        ! get dimensions from first line
        read(iunit, *) this%nmodes, ncols

        write(message(1), '(3a,i7,a,i3,a)') 'Reading file ', trim(filename), ' with ', &
          this%nmodes, ' photon modes and ', ncols, ' columns.'
        call messages_info(1)

        SAFE_ALLOCATE(this%omega_array(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda_array(1:this%nmodes))
        SAFE_ALLOCATE(this%pol_array(1:this%nmodes,3))

        ! now read in all modes
        do ii = 1, this%nmodes
          if(ncols == 5) then
            read(iunit, *) this%omega_array(ii), this%lambda_array(ii), &
              this%pol_array(ii,1), this%pol_array(ii,2), this%pol_array(ii,3)
          else
            ! error if not 5 columns
            message(1) = 'Error: unexpected number of columns in '//filename
            call messages_fatal(1)
          end if
        end do
        call io_close(iunit)
      end if
#ifdef HAVE_MPI
      ! broadcast first array dimensions, then allocate and broadcast arrays
      call MPI_Bcast(this%nmodes, 1, MPI_INTEGER, 0, mpi_world%comm, ierr)
      call MPI_Bcast(ncols, 1, MPI_INTEGER, 0, mpi_world%comm, ierr)
      if(.not. mpi_grp_is_root(mpi_world)) then
        SAFE_ALLOCATE(this%omega_array(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda_array(1:this%nmodes))
        SAFE_ALLOCATE(this%pol_array(1:this%nmodes,3))
      end if
      call MPI_Bcast(this%omega_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
      call MPI_Bcast(this%lambda_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
      call MPI_Bcast(this%pol_array(1,1), this%nmodes*3, MPI_FLOAT, 0, mpi_world%comm, ierr)
#endif
    end if


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
      write(message(1), '(a,i5,a)') 'Info: Happy to have ', this%nmodes, ' photon modes with us.'
    else
      call messages_write('You need to specify the photon modes!')
      call messages_fatal()
    end if

    this%ex = M_ZERO
    SAFE_ALLOCATE(this%pt_number(1:this%nmodes))
    this%pt_number = M_ZERO

    SAFE_ALLOCATE(this%correlator(1:gr%mesh%np,1:this%nmodes))
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
    SAFE_DEALLOCATE_A(this%pt_number)

    SAFE_DEALLOCATE_A(this%pol_array)
    SAFE_DEALLOCATE_A(this%pol_dipole_array)

    POP_SUB(photon_mode_end)
  end subroutine photon_mode_end


end module photon_mode_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
