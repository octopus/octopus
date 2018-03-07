!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module ps_qso_oct_m
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use profiling_oct_m
  use pseudo_oct_m
  use ps_in_grid_oct_m
  use xml_oct_m

  implicit none

  private
  public ::     &
    ps_qso_t,       &
    ps_qso_init,    &
    ps_qso_end

  type ps_qso_t
    logical            :: oncv
    logical            :: nlcc
    integer            :: atomic_number
    FLOAT              :: mass
    FLOAT              :: valence_charge
    integer            :: lmax
    integer            :: llocal
    FLOAT              :: mesh_spacing
    integer            :: nchannels
    integer            :: grid_size
    FLOAT, allocatable :: potential(:, :)
    FLOAT, allocatable :: wavefunction(:, :)
    FLOAT, allocatable :: projector(:, :, :)
    FLOAT, allocatable :: dij(:, :, :)
    FLOAT, allocatable :: nlcc_density(:)
  end type ps_qso_t

contains

  ! ---------------------------------------------------------
  subroutine ps_qso_init(this, filename, ierr)
    type(ps_qso_t),   intent(inout) :: this
    character(len=*), intent(in)    :: filename
    integer,          intent(out)   :: ierr

    integer :: ll, ii, ic, jc
    type(pseudo_t) :: pseudo

    PUSH_SUB(ps_qso_init)

    call pseudo_init(pseudo, filename, ierr)

    if(ierr /= 0) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not found")
      call messages_fatal()
    end if

    this%valence_charge = pseudo_valence_charge(pseudo)
    this%mesh_spacing = pseudo_mesh_spacing(pseudo)
    this%mass = pseudo_mass(pseudo)
    this%lmax = pseudo_lmax(pseudo)
    this%llocal = pseudo_llocal(pseudo) 
    this%nchannels = pseudo_nchannels(pseudo)
   
    this%oncv = (pseudo_type(pseudo) == PSEUDO_TYPE_KLEINMAN_BYLANDER)

    this%grid_size = pseudo_mesh_size(pseudo)
    
    if(.not. this%oncv) then

      SAFE_ALLOCATE(this%potential(1:this%grid_size, 0:this%lmax))
      SAFE_ALLOCATE(this%wavefunction(1:this%grid_size, 0:this%lmax))
      
      do ll = 0, this%lmax
        call pseudo_radial_potential(pseudo, ll, this%potential(1, ll))
        call pseudo_radial_function(pseudo, ll, this%wavefunction(1, ll))
      end do

    else

      SAFE_ALLOCATE(this%potential(1:this%grid_size, -1:-1))
      
      call pseudo_local_potential(pseudo, this%potential(1, -1))

      SAFE_ALLOCATE(this%projector(1:this%grid_size, 0:3, 1:2))
      
      do ll = 0, this%lmax
        do ic = 1, this%nchannels
          call pseudo_projector(pseudo, ll, ic, this%projector(1, ll, ic))
        end do
      end do

      SAFE_ALLOCATE(this%dij(0:this%lmax, 1:this%nchannels, 1:this%nchannels))

      this%dij = CNST(0.0)
      do ll = 0, this%lmax
        do ic = 1, this%nchannels
          do jc = 1, this%nchannels
            call pseudo_dij(pseudo, ll, ic, jc, this%dij(ll, ic, jc))
          end do
        end do
      end do

    end if

    this%nlcc = pseudo_has_nlcc(pseudo)
    if(this%nlcc) then
      SAFE_ALLOCATE(this%nlcc_density(1:this%grid_size))
      call pseudo_nlcc_density(pseudo, this%nlcc_density(1))
    end if

    call pseudo_end(pseudo)
    
    if(.not. this%oncv) call ps_qso_check_normalization(this)

    POP_SUB(ps_qso_init)
  end subroutine ps_qso_init

  ! ---------------------------------------------------------
  !> checks normalization of the pseudo wavefunctions
  subroutine ps_qso_check_normalization(this)
    type(ps_qso_t), intent(in) :: this
    
    integer :: ll, ip
    FLOAT   :: nrm, rr

    PUSH_SUB(ps_qso_check_normalization)

    !  checking normalization of the wavefunctions
    do ll = 0, this%lmax
      nrm = 0.0
      do ip = 1, this%grid_size
        rr = (ip - 1)*this%mesh_spacing
        nrm = nrm + this%wavefunction(ip, ll)**2*this%mesh_spacing*rr**2
      end do
      nrm = sqrt(nrm)

      nrm = abs(nrm - M_ONE)
      if (nrm > CNST(1.0e-5)) then
        write(message(1), '(a,i2,a)') "Eigenstate for l = ", ll, ' is not normalized'
        write(message(2), '(a, f12.6,a)') '(abs(1 - norm) = ', nrm, ')'
        call messages_warning(2)
      end if

    end do
      
    POP_SUB(ps_qso_check_normalization)
  end subroutine ps_qso_check_normalization
  
  ! ---------------------------------------------------------
  subroutine ps_qso_end(this)
    type(ps_qso_t), intent(inout) :: this

    PUSH_SUB(ps_qso_end)

    SAFE_DEALLOCATE_A(this%potential)
    SAFE_DEALLOCATE_A(this%wavefunction)
    SAFE_DEALLOCATE_A(this%projector)
    SAFE_DEALLOCATE_A(this%dij)
    SAFE_DEALLOCATE_A(this%nlcc_density)

    POP_SUB(ps_qso_end)
  end subroutine ps_qso_end

end module ps_qso_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
