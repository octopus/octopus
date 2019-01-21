!! Copyright (C) 2015-2018 Xavier Andrade
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

module ps_xml_oct_m
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use profiling_oct_m
  use pseudo_oct_m
  use ps_in_grid_oct_m

  implicit none

  private
  public ::     &
    ps_xml_t,       &
    ps_xml_init,    &
    ps_xml_end

  type ps_xml_t
    logical            :: initialized
    logical            :: kleinman_bylander
    logical            :: nlcc
    logical            :: has_density
    integer            :: atomic_number
    FLOAT              :: mass
    FLOAT              :: valence_charge
    integer            :: lmax
    integer            :: llocal
    integer            :: nchannels
    integer            :: grid_size
    integer            :: nwavefunctions
    FLOAT, allocatable :: grid(:)
    FLOAT, allocatable :: weights(:)
    FLOAT, allocatable :: potential(:, :)
    FLOAT, allocatable :: wavefunction(:, :)
    FLOAT, allocatable :: projector(:, :, :)
    FLOAT, allocatable :: dij(:, :, :)
    FLOAT, allocatable :: nlcc_density(:)
    FLOAT, allocatable :: density(:)
    integer, allocatable :: wf_n(:)
    integer, allocatable :: wf_l(:)
    FLOAT, allocatable :: wf_occ(:)
    type(pseudo_t)     :: pseudo
  end type ps_xml_t

contains

  ! ---------------------------------------------------------
  subroutine ps_xml_init(this, filename, fmt, ierr)
    type(ps_xml_t),   intent(inout) :: this
    character(len=*), intent(in)    :: filename
    integer,          intent(in)    :: fmt
    integer,          intent(out)   :: ierr

    integer :: ll, ii, ic, jc
    type(pseudo_t) :: pseudo

    PUSH_SUB(ps_xml_init)

    this%initialized = .false.
    
    call pseudo_init(pseudo, filename, fmt, ierr)

    if(ierr == PSEUDO_STATUS_FILE_NOT_FOUND) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not found")
      call messages_fatal()
    end if

    if(ierr == PSEUDO_STATUS_UNKNOWN_FORMAT) then
      call messages_write("Cannot determine the format for pseudopotential file '" // trim(filename) // "'")
      call messages_fatal()
    end if

    if(ierr == PSEUDO_STATUS_UNSUPPORTED_TYPE_ULTRASOFT) then
      call messages_write("Ultrasoft pseudopotential file '" // trim(filename) // "' not supported")
      call messages_fatal()
    end if

    if(ierr == PSEUDO_STATUS_UNSUPPORTED_TYPE_PAW) then
      call messages_write("PAW pseudopotential file '" // trim(filename) // "' not supported")
      call messages_fatal()
    end if
    
    if(ierr == PSEUDO_STATUS_UNSUPPORTED_TYPE) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not supported")
      call messages_fatal()
    end if
    
    if(ierr == PSEUDO_STATUS_FORMAT_NOT_SUPPORTED) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not supported")
      call messages_fatal()
    end if

    this%initialized = .true.
    this%valence_charge = pseudo_valence_charge(pseudo)
    this%mass = pseudo_mass(pseudo)
    this%lmax = pseudo_lmax(pseudo)
    this%llocal = pseudo_llocal(pseudo) 
    this%nchannels = pseudo_nchannels(pseudo)
   
    this%kleinman_bylander = (pseudo_type(pseudo) == PSEUDO_TYPE_KLEINMAN_BYLANDER)

    this%grid_size = pseudo_mesh_size(pseudo)

    SAFE_ALLOCATE(this%grid(1:this%grid_size))
    SAFE_ALLOCATE(this%weights(1:this%grid_size))
    
    call pseudo_grid(pseudo, this%grid(1))
    call pseudo_grid_weights(pseudo, this%weights(1))
    
    if(.not. this%kleinman_bylander) then

      SAFE_ALLOCATE(this%potential(1:this%grid_size, 0:this%lmax))
      SAFE_ALLOCATE(this%wavefunction(1:this%grid_size, 0:this%lmax))
      
      do ll = 0, this%lmax
        if(.not. pseudo_has_radial_function(pseudo, ll)) then
          call messages_write("The pseudopotential file '"//trim(filename)//"' does not contain")
          call messages_new_line()
          call messages_write("the wave functions. Octopus cannot use it.")
          call messages_fatal()
        end if

        call pseudo_radial_function(pseudo, ll, this%wavefunction(1, ll))
        call pseudo_radial_potential(pseudo, ll, this%potential(1, ll))
      end do

    else

      SAFE_ALLOCATE(this%potential(1:this%grid_size, this%llocal:this%llocal))
      
      call pseudo_local_potential(pseudo, this%potential(1, this%llocal))

      SAFE_ALLOCATE(this%projector(1:this%grid_size, 0:this%lmax, 1:this%nchannels))

      this%projector = CNST(0.0)
      
      do ll = 0, this%lmax
        if(this%llocal == ll) cycle
        do ic = 1, this%nchannels
          call pseudo_projector(pseudo, ll, ic, this%projector(1, ll, ic))
        end do
      end do

      SAFE_ALLOCATE(this%dij(0:this%lmax, 1:this%nchannels, 1:this%nchannels))

      this%dij = CNST(0.0)
      do ll = 0, this%lmax
        do ic = 1, this%nchannels
          do jc = 1, this%nchannels
            this%dij(ll, ic, jc) = pseudo_dij(pseudo, ll, ic, jc)
          end do
        end do
      end do

      this%nwavefunctions = pseudo_nwavefunctions(pseudo)
      
      SAFE_ALLOCATE(this%wavefunction(1:this%grid_size, 1:this%nwavefunctions))
      SAFE_ALLOCATE(this%wf_n(1:this%nwavefunctions))
      SAFE_ALLOCATE(this%wf_l(1:this%nwavefunctions))
      SAFE_ALLOCATE(this%wf_occ(1:this%nwavefunctions))
      
      do ii = 1, this%nwavefunctions
        call pseudo_wavefunction(pseudo, ii, this%wf_n(ii), this%wf_l(ii), this%wf_occ(ii), this%wavefunction(1, ii))
      end do
      
    end if

    this%has_density = pseudo_has_density(pseudo)

    if(this%has_density) then
      SAFE_ALLOCATE(this%density(1:this%grid_size))
      call pseudo_density(pseudo, this%density(1))
    end if
    
    this%nlcc = pseudo_has_nlcc(pseudo)
    if(this%nlcc) then
      SAFE_ALLOCATE(this%nlcc_density(1:this%grid_size))
      call pseudo_nlcc_density(pseudo, this%nlcc_density(1))
    end if

    if(.not. this%kleinman_bylander) call ps_xml_check_normalization(this)

    this%pseudo = pseudo
    
    POP_SUB(ps_xml_init)
  end subroutine ps_xml_init

  ! ---------------------------------------------------------
  !> checks normalization of the pseudo wavefunctions
  subroutine ps_xml_check_normalization(this)
    type(ps_xml_t), intent(in) :: this
    
    integer :: ll, ip
    FLOAT   :: nrm, rr

    PUSH_SUB(ps_xml_check_normalization)

    !  checking normalization of the wavefunctions
    do ll = 0, this%lmax
      nrm = M_ZERO
      do ip = 1, this%grid_size
        rr = this%grid(ip)
        nrm = nrm + this%wavefunction(ip, ll)**2*this%weights(ip)*rr**2
      end do
      nrm = sqrt(nrm)

      nrm = abs(nrm - M_ONE)
      if (nrm > CNST(1.0e-5)) then
        write(message(1), '(a,i2,a)') "Eigenstate for l = ", ll, ' is not normalized'
        write(message(2), '(a, f12.6,a)') '(abs(1 - norm) = ', nrm, ')'
        call messages_warning(2)
      end if

    end do
    
    POP_SUB(ps_xml_check_normalization)
  end subroutine ps_xml_check_normalization
  
  ! ---------------------------------------------------------
  subroutine ps_xml_end(this)
    type(ps_xml_t), intent(inout) :: this

    PUSH_SUB(ps_xml_end)

    call pseudo_end(this%pseudo)

    SAFE_DEALLOCATE_A(this%grid)
    SAFE_DEALLOCATE_A(this%weights)
    SAFE_DEALLOCATE_A(this%potential)
    SAFE_DEALLOCATE_A(this%wavefunction)
    SAFE_DEALLOCATE_A(this%projector)
    SAFE_DEALLOCATE_A(this%dij)
    SAFE_DEALLOCATE_A(this%nlcc_density)
    SAFE_DEALLOCATE_A(this%density)
    SAFE_DEALLOCATE_A(this%wf_n)
    SAFE_DEALLOCATE_A(this%wf_l)
    SAFE_DEALLOCATE_A(this%wf_occ)

    POP_SUB(ps_xml_end)
  end subroutine ps_xml_end

end module ps_xml_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
