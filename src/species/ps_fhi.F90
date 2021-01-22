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

module ps_fhi_oct_m
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use ps_cpi_file_oct_m
  use ps_cpi_oct_m
  use ps_fhi_file_oct_m
  use ps_in_grid_oct_m

  implicit none

  private
  public ::     &
    ps_fhi_t,       &
    ps_fhi_init,    &
    ps_fhi_end,     &
    ps_fhi_process

  !> remember that the FHI format is basically the CPI format with a header
  type ps_fhi_t
    ! Components are public by default
    type(ps_fhi_file_t), allocatable, private :: fhi_file !< This just includes the extra header
    type(ps_cpi_file_t), allocatable, private :: cpi_file !< This includes the real pseudopotential
    type(ps_in_grid_t),  allocatable          :: ps_grid  !< the pseudopotential in the grid
    type(valconf_t),     allocatable, private :: conf
  end type ps_fhi_t

contains

  ! ---------------------------------------------------------
  subroutine ps_fhi_init(ps_fhi, filename, namespace)
    type(ps_fhi_t),    intent(inout) :: ps_fhi
    character(len=*),  intent(in)    :: filename
    type(namespace_t), intent(in)    :: namespace

    integer :: iunit
    logical :: found

    PUSH_SUB(ps_fhi_init)

    SAFE_ALLOCATE(ps_fhi%fhi_file)
    SAFE_ALLOCATE(ps_fhi%cpi_file)
    SAFE_ALLOCATE(ps_fhi%ps_grid)
    SAFE_ALLOCATE(ps_fhi%conf)

    inquire(file = filename, exist = found)

    if(.not.found) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not found")
      call messages_fatal(namespace=namespace)
    end if

    iunit = io_open(filename, action='read', form='formatted', status='old')
    call ps_fhi_file_read(iunit, ps_fhi%fhi_file, namespace)
    call ps_cpi_file_read(iunit, ps_fhi%cpi_file)
    call io_close(iunit)

    call ps_cpi_file_to_grid(ps_fhi%cpi_file, ps_fhi%ps_grid)

    POP_SUB(ps_fhi_init)
  end subroutine ps_fhi_init

  
  ! ---------------------------------------------------------
  subroutine ps_fhi_end(ps_fhi)
    type(ps_fhi_t), intent(inout) :: ps_fhi

    PUSH_SUB(ps_fhi_end)

    SAFE_DEALLOCATE_A(ps_fhi%fhi_file)

    call ps_cpi_file_end(ps_fhi%cpi_file)
    SAFE_DEALLOCATE_A(ps_fhi%cpi_file)
    SAFE_DEALLOCATE_A(ps_fhi%conf)

    call ps_in_grid_end(ps_fhi%ps_grid)
    SAFE_DEALLOCATE_A(ps_fhi%ps_grid)

    POP_SUB(ps_fhi_end)
  end subroutine ps_fhi_end


  ! ---------------------------------------------------------
  subroutine ps_fhi_process(ps_fhi, lmax, lloc, namespace)
    type(ps_fhi_t),    intent(inout) :: ps_fhi
    integer,           intent(in)    :: lmax, lloc
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(ps_fhi_process)

    if(lmax /= ps_fhi%fhi_file%lmax) then
      message(1) = "Inconsistency in pseudopotential :"
      write(message(2),'(a,i2,a,i2)') "  Input file says lmax = ", lmax, &
        " but ps file says lmax = ", ps_fhi%fhi_file%lmax
      call messages_warning(2, namespace=namespace)
    end if
    if(lloc /= ps_fhi%fhi_file%lloc) then
      message(1) = "Inconsistency in pseudopotential :"
      write(message(2),'(a,i2,a,i2)') "  Input file says lloc = ", lloc, &
        " but ps file says lloc = ", ps_fhi%fhi_file%lloc
      call messages_warning(2, namespace=namespace)
    end if

    ! check norm of rphi
    call ps_in_grid_check_rphi(ps_fhi%ps_grid, namespace)

    ! Fix the local potential. Final argument is the core radius
    call ps_in_grid_vlocal(ps_fhi%ps_grid, lloc, M_THREE, namespace)

    ! Calculate kb cosines and norms
    call ps_in_grid_kb_cosines(ps_fhi%ps_grid, lloc)

    ! Define the KB-projector cut-off radii
    call ps_in_grid_cutoff_radii(ps_fhi%ps_grid, lloc)

    ! Calculate KB-projectors
    call ps_in_grid_kb_projectors(ps_fhi%ps_grid)

    POP_SUB(ps_fhi_process)
  end subroutine ps_fhi_process

end module ps_fhi_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
