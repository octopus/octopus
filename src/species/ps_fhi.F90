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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: tm.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module ps_fhi_m
  use atomic_m
  use global_m
  use io_m
  use messages_m
  use profiling_m
  use ps_cpi_file_m
  use ps_cpi_m
  use ps_fhi_file_m
  use ps_in_grid_m

  implicit none

  private
  public ::     &
    ps_fhi_t,       &
    ps_fhi_init,    &
    ps_fhi_end,     &
    ps_fhi_process

  ! remember that the FHI format is basically the CPI format with a header
  type ps_fhi_t
    type(ps_fhi_file_t), pointer :: fhi_file ! This just includes the extra header
    type(ps_cpi_file_t), pointer :: cpi_file ! This includes the real pseudopotential
    type(ps_in_grid_t),  pointer :: ps_grid  ! the pseudopotential in the grid
    type(valconf_t),     pointer :: conf     ! 
  end type ps_fhi_t

contains

  ! ---------------------------------------------------------
  subroutine ps_fhi_init(ps_fhi, filename)
    type(ps_fhi_t),   intent(inout) :: ps_fhi
    character(len=*), intent(in)    :: filename

    character(len=256) :: filename2
    integer :: iunit
    logical :: found

    PUSH_SUB(ps_fhi_init)

    SAFE_ALLOCATE(ps_fhi%fhi_file)
    SAFE_ALLOCATE(ps_fhi%cpi_file)
    SAFE_ALLOCATE(ps_fhi%ps_grid)
    SAFE_ALLOCATE(ps_fhi%conf)

    ! Find out where the hell the file is.
    filename2 = trim(filename) // '.fhi'
    inquire(file=filename2, exist=found)
    if(.not.found) then
      filename2 = trim(conf%share) // "/PP/FHI/" // trim(filename) // ".fhi"
      inquire(file=filename2, exist=found)
      if(.not.found) then
        message(1) = "Pseudopotential file '" // trim(filename) // ".fhi' not found"
        call messages_fatal(1)
      end if
    end if

    message(1) = "Reading pseudopotential from file:"
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call messages_info(2)

    iunit = io_open(filename2, action='read', form='formatted', status='old', is_tmp=.true.)
    call ps_fhi_file_read(iunit, ps_fhi%fhi_file)
    call ps_cpi_file_read(iunit, ps_fhi%cpi_file)
    call io_close(iunit)

    call ps_cpi_file_to_grid(ps_fhi%cpi_file, ps_fhi%ps_grid)

    POP_SUB(ps_fhi_init)
  end subroutine ps_fhi_init

  
  ! ---------------------------------------------------------
  subroutine ps_fhi_end(ps_fhi)
    type(ps_fhi_t), intent(inout) :: ps_fhi

    SAFE_DEALLOCATE_P(ps_fhi%fhi_file)
    SAFE_DEALLOCATE_P(ps_fhi%cpi_file)
    SAFE_DEALLOCATE_P(ps_fhi%ps_grid)
    SAFE_DEALLOCATE_P(ps_fhi%conf)
  end subroutine ps_fhi_end


  ! ---------------------------------------------------------
  subroutine ps_fhi_process(ps_fhi, lmax, lloc)
    type(ps_fhi_t), intent(inout) :: ps_fhi
    integer,        intent(in)    :: lmax, lloc

    PUSH_SUB(ps_fhi_process)

    if(lmax.ne.ps_fhi%fhi_file%lmax) then
      message(1) = "Inconsistency in pseudopotential :"
      write(message(2),'(a,i2,a,i2)') "  Input file says lmax = ", lmax, &
        " but ps file says lmax = ", ps_fhi%fhi_file%lmax
      call messages_warning(2)
    end if
    if(lloc.ne.ps_fhi%fhi_file%lloc) then
      message(1) = "Inconsistency in pseudopotential :"
      write(message(2),'(a,i2,a,i2)') "  Input file says lloc = ", lloc, &
        " but ps file says lloc = ", ps_fhi%fhi_file%lloc
      call messages_warning(2)
    end if

    ! check norm of rphi
    call ps_in_grid_check_rphi(ps_fhi%ps_grid)

    ! Fix the local potential. Final argument is the core radius
    call ps_in_grid_vlocal(ps_fhi%ps_grid, lloc, M_THREE)

    ! Calculate kb cosines and norms
    call ps_in_grid_kb_cosines(ps_fhi%ps_grid, lloc)

    ! Define the KB-projector cut-off radii
    call ps_in_grid_cutoff_radii(ps_fhi%ps_grid, lloc)

    ! Calculate KB-projectors
    call ps_in_grid_kb_projectors(ps_fhi%ps_grid)

    POP_SUB(ps_fhi_process)
  end subroutine ps_fhi_process

end module ps_fhi_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
