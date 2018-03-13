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
  use hamiltonian_oct_m
  use io_oct_m
  use lalg_adv_oct_m
  use linear_response_oct_m
  use linear_solver_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use states_oct_m
  use states_dim_oct_m
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
    FLOAT                 :: omega, lambda        ! omega and lambda for photons
    FLOAT                 :: pol_x, pol_y, pol_z  ! polarization of the photon field
    FLOAT                 :: ex                   ! photon exchange energy
    FLOAT                 :: pt_number            ! number of photons in mode
    FLOAT, pointer        :: pol_dipole(:,:)      ! polarization*dipole operator
    FLOAT, pointer        :: correlator(:,:)      ! correlation function <n(r)(ad+a)>
    type(lr_t)            :: lr         !< to solve the equation H psi = b
    integer               :: nmodes
    FLOAT, pointer        :: omega_array(:), lambda_array(:)  ! the arrays are only temporary, later merge, OEP has to be adapted
    FLOAT, pointer        :: pol_array(:,:)             ! polarization of the photon field
    FLOAT, pointer        :: q0_array(:)
    FLOAT, pointer        :: p0_array(:)
    logical               :: has_arrays                 ! will become obsolete after single mode removal
    logical               :: has_q0_p0
  end type photon_mode_t

contains

  ! ---------------------------------------------------------
  subroutine photon_mode_init(this, gr)
    type(photon_mode_t),  intent(out)   :: this
    type(grid_t),         intent(inout) :: gr

    type(block_t)         :: blk
    integer               :: ii, ncols, iunit, ierr
    character(MAX_PATH_LEN) :: filename
    logical               :: file_exists

    PUSH_SUB(photon_mode_init)

    this%nmodes = 0
    this%has_arrays = .false.
    this%has_q0_p0 = .false.

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
    call parse_variable('PhotonmodesFilename', 'photonmodes', filename)
    inquire(file=trim(filename), exist=file_exists)
    if(file_exists) then
      this%has_arrays = .true.
      if(mpi_grp_is_root(mpi_world)) then
        message(1) = 'Opening '//trim(filename)
        call messages_info(1)
        ! open file on root
        iunit = io_open(trim(filename), action='read', form='formatted')

        ! get dimensions from first line
        read(iunit, *) this%nmodes, ncols

        write(message(1), '(3a,i7,a,i3,a)') 'Reading file ', trim(filename), ' with ', &
          this%nmodes, ' photon modes and ', ncols, ' columns.'
        call messages_info(1)

        SAFE_ALLOCATE(this%omega_array(1:this%nmodes))
        SAFE_ALLOCATE(this%lambda_array(1:this%nmodes))
        SAFE_ALLOCATE(this%pol_array(1:this%nmodes,3))
        if (ncols > 5) then
          this%has_q0_p0 = .true.
          SAFE_ALLOCATE(this%q0_array(1:this%nmodes))
          SAFE_ALLOCATE(this%p0_array(1:this%nmodes))
        end if

        ! now read in all modes
        do ii = 1, this%nmodes
          if(ncols == 5) then
            read(iunit, *) this%omega_array(ii), this%lambda_array(ii), &
              this%pol_array(ii,1), this%pol_array(ii,2), this%pol_array(ii,3)
          else if(ncols == 7) then
            read(iunit, *) this%omega_array(ii), this%lambda_array(ii), &
              this%pol_array(ii,1), this%pol_array(ii,2), this%pol_array(ii,3), &
              this%q0_array(ii), this%p0_array(ii)
          else
            ! error if not 5 or 7 columns
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
        if(ncols > 5) then
          this%has_q0_p0 = .true.
          SAFE_ALLOCATE(this%q0_array(1:this%nmodes))
          SAFE_ALLOCATE(this%p0_array(1:this%nmodes))
        end if
      end if
      call MPI_Bcast(this%omega_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
      call MPI_Bcast(this%lambda_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
      call MPI_Bcast(this%pol_array(1,1), this%nmodes*3, MPI_FLOAT, 0, mpi_world%comm, ierr)
      if(this%has_q0_p0) then
        call MPI_Bcast(this%q0_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
        call MPI_Bcast(this%p0_array(1), this%nmodes, MPI_FLOAT, 0, mpi_world%comm, ierr)
      end if
#endif
    end if


    !%Variable PhotonModes
    !%Type block
    !%Section Hamiltonian::XC
    !%Description
    !% Syntax:
    !%PhotonModes
    !%omega1 | lambda 1| PolX1 | PolY1 | PolZ1 | q0 | p0
    !%...
    !%
    !%End
    !% frequency, coupling strength, pol in (x,y,z), q(t0), p(t0)
    ! todo extend to more modes and put defaults
    if(.not. file_exists .and. parse_block('PhotonModes', blk) == 0) then
      this%has_arrays = .true.
      this%nmodes = parse_block_n(blk)
      SAFE_ALLOCATE(this%omega_array(1:this%nmodes))
      SAFE_ALLOCATE(this%lambda_array(1:this%nmodes))
      SAFE_ALLOCATE(this%pol_array(1:this%nmodes,3))
      do ii = 1, this%nmodes
         ncols = parse_block_cols(blk, ii-1)
         call parse_block_float(blk, ii-1, 0, this%omega_array(ii))   !row, column
         call parse_block_float(blk, ii-1, 1, this%lambda_array(ii))  !row, column
         call parse_block_float(blk, ii-1, 2, this%pol_array(ii,1))   !row, column
         call parse_block_float(blk, ii-1, 3, this%pol_array(ii,2))   !row, column
         call parse_block_float(blk, ii-1, 4, this%pol_array(ii,3))   !row, column
         if (ncols > 5) then
           this%has_q0_p0 = .true.
           SAFE_ALLOCATE(this%q0_array(1:this%nmodes))
           SAFE_ALLOCATE(this%p0_array(1:this%nmodes))
           call parse_block_float(blk, ii-1, 5, this%q0_array(ii))   !row, column
           call parse_block_float(blk, ii-1, 6, this%p0_array(ii))   !row, column
         end if
      end do
      call parse_block_end(blk)
      write(message(1), '(a,i5,a)') 'Info: Happy to have ', this%nmodes, ' photon modes with us.'
    else
      !The following syntax is depreciated, but still in for downwards compatibility, Todo remove, OEP only works with this syntax
      !%Variable OEPOmega
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% Frequency of the single photon mode
      !%End
      call parse_variable('OEPOmega', M_ONE, this%omega)

      !%Variable OEPLambda
      !%Type float
      !%Default 0.01
      !%Section Hamiltonian::XC
      !%Description
      !% Electron-photon coupling strenght of the single photon mode
      !%End
      call parse_variable('OEPLambda', CNST(1e-2), this%lambda)

      !%Variable OEPPolX
      !%Type float
      !%Default 1.0
      !%Section Hamiltonian::XC
      !%Description
      !% polarization component of photon field in x-direction
      !%End
      call parse_variable('OEPPolX', M_ONE, this%pol_x)

      !%Variable OEPPolY
      !%Type float
      !%Default 0.0
      !%Section Hamiltonian::XC
      !%Description
      !% polarization component of photon field in y-direction
      !%End
      call parse_variable('OEPPolY', M_ZERO, this%pol_y)

      !%Variable OEPPolZ
      !%Type float
      !%Default 0.0
      !%Section Hamiltonian::XC
      !%Description
      !% polarization component of photon field in z-direction
      !%End
      call parse_variable('OEPPolZ', M_ZERO, this%pol_z)
    end if

    this%ex = M_ZERO
    this%pt_number = M_ZERO

    SAFE_ALLOCATE(this%correlator(1:gr%mesh%np,1))
    this%correlator = M_ZERO

    SAFE_ALLOCATE(this%pol_dipole(1:gr%mesh%np,1))
    this%pol_dipole(1:gr%mesh%np,1) = this%pol_x*gr%mesh%x(1:gr%mesh%np, 1) &
                                      + this%pol_y*gr%mesh%x(1:gr%mesh%np, 2) &
                                      + this%pol_z*gr%mesh%x(1:gr%mesh%np, 3)

    POP_SUB(photon_mode_init)
  end subroutine photon_mode_init

  ! ---------------------------------------------------------

  subroutine photon_mode_end(this)
    type(photon_mode_t), intent(inout) :: this

    PUSH_SUB(photon_mode_end)

    SAFE_DEALLOCATE_P(this%pol_dipole)
    SAFE_DEALLOCATE_P(this%correlator)

    if(this%has_arrays) then
      SAFE_DEALLOCATE_P(this%omega_array)
      SAFE_DEALLOCATE_P(this%lambda_array)
      SAFE_DEALLOCATE_P(this%pol_array)
      if (this%has_q0_p0) then 
        SAFE_DEALLOCATE_P(this%q0_array)
        SAFE_DEALLOCATE_P(this%p0_array)
      end if
    end if

    POP_SUB(photon_mode_end)
  end subroutine photon_mode_end


end module photon_mode_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
