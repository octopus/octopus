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

module ps_cpi_oct_m
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use logrid_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use ps_cpi_file_oct_m
  use ps_in_grid_oct_m

  implicit none

  private
  public ::              &
    ps_cpi_t,            &
    ps_cpi_init,         &
    ps_cpi_end,          &
    ps_cpi_file_to_grid, &
    ps_cpi_process

  type ps_cpi_t
    ! Components are public by default
    type(ps_cpi_file_t), allocatable, private :: cpi_file
    type(ps_in_grid_t),  allocatable          :: ps_grid

    type(valconf_t),     allocatable, private :: conf    !< what to do with this?
  end type ps_cpi_t

contains

  ! ---------------------------------------------------------
  subroutine ps_cpi_init(ps_cpi, filename, namespace)
    type(ps_cpi_t),    intent(inout) :: ps_cpi
    character(len=*),  intent(in)    :: filename
    type(namespace_t), intent(in)    :: namespace

    integer :: iunit
    logical :: found

    PUSH_SUB(ps_cpi_init)

    ! allocate data
    SAFE_ALLOCATE(ps_cpi%cpi_file)
    SAFE_ALLOCATE(ps_cpi%ps_grid)
    SAFE_ALLOCATE(ps_cpi%conf)

    inquire(file = filename, exist = found)
    if(.not.found) then
      call messages_write("Pseudopotential file '" // trim(filename) // "' not found")
      call messages_fatal(namespace=namespace)
    end if
    
    iunit = io_open(filename, action='read', form='formatted', status='old')
    call ps_cpi_file_read(iunit, ps_cpi%cpi_file)
    call io_close(iunit)

    ! Fills the valence configuration data.
    !call valconf_guess(ps_cpi%conf, )

    call ps_cpi_file_to_grid(ps_cpi%cpi_file, ps_cpi%ps_grid)

    POP_SUB(ps_cpi_init)
  end subroutine ps_cpi_init

  
  ! ---------------------------------------------------------
  subroutine ps_cpi_end(ps_cpi)
    type(ps_cpi_t), intent(inout) :: ps_cpi

    call ps_in_grid_end (ps_cpi%ps_grid)
    call ps_cpi_file_end(ps_cpi%cpi_file)

    SAFE_DEALLOCATE_A(ps_cpi%cpi_file)
    SAFE_DEALLOCATE_A(ps_cpi%ps_grid)
    SAFE_DEALLOCATE_A(ps_cpi%conf)

  end subroutine ps_cpi_end


  !----------------------------------------------------------------
  subroutine ps_cpi_file_to_grid(cpi_file, ps_grid)
    type(ps_cpi_file_t), intent(in)  :: cpi_file
    type(ps_in_grid_t),  intent(out) :: ps_grid

    ! Initializes the pseudo in the logarithmic grid.
    call ps_in_grid_init(ps_grid,                      &
      LOGRID_CPI, cpi_file%a, cpi_file%rofi(2), cpi_file%nr,  &
      cpi_file%no_l_channels, 0)
    
    ps_grid%zval        = cpi_file%zval
    ps_grid%vps(:,:)    = cpi_file%vps(:,:)
    ps_grid%rphi(:,:,1) = cpi_file%rphi(:,:)
    ps_grid%rphi(:,:,2) = cpi_file%rphi(:,:)
    ps_grid%rphi(:,:,3) = cpi_file%rphi(:,:)

    ps_grid%core_corrections = cpi_file%core_corrections
    if(ps_grid%core_corrections) then
      ps_grid%chcore(:) = cpi_file%chcore(:)
    end if

  end subroutine ps_cpi_file_to_grid


  ! ---------------------------------------------------------
  subroutine ps_cpi_process(ps_cpi, lloc, namespace)
    type(ps_cpi_t),    intent(inout) :: ps_cpi
    integer,           intent(in)    :: lloc
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(ps_cpi_process)

    ! check norm of rphi
    call ps_in_grid_check_rphi(ps_cpi%ps_grid, namespace)

    ! Fix the local potential. Final argument is the core radius
    call ps_in_grid_vlocal(ps_cpi%ps_grid, lloc, M_THREE, namespace)

    ! Calculate kb cosines and norms
    call ps_in_grid_kb_cosines(ps_cpi%ps_grid, lloc)

    ! Define the KB-projector cut-off radii
    call ps_in_grid_cutoff_radii(ps_cpi%ps_grid, lloc)

    ! Calculate KB-projectors
    call ps_in_grid_kb_projectors(ps_cpi%ps_grid)

    POP_SUB(ps_cpi_process)
  end subroutine ps_cpi_process

end module ps_cpi_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
