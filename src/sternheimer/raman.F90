!! Copyright (C) 2004 Xavier Andrade
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
!! $Id: em_resp.F90 2686 2007-02-03 22:10:51Z xavier $

#include "global.h"

module raman_m
  use datasets_m
  use em_resp_calc_m
  use geometry_m
  use global_m
  use grid_m
  use h_sys_output_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use linear_response_m
  use parser_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use phonons_lr_m
  use restart_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m

  implicit none

  private
  public :: &
       raman_run
  
contains

  ! ---------------------------------------------------------
  subroutine raman_run(sys, hm, fromscratch)
    type(system_t),         intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: hm
    logical,                intent(in)    :: fromscratch

    type(lr_t)          :: psi_elec(MAX_DIM)
    type(lr_t)          :: psi_vib(1)

    integer :: idir, ierr, sigma
    character(len=100) :: dirname, str_tmp

    integer :: nnm, inm

    !CONSTRUCT
    PUSH_SUB(raman_run)

    call messages_devel_version("Raman response")

    call restart_look_and_read(sys%st, sys%gr, sys%geo)

    ! read electric perturbation
    sigma = 1
    do idir = 1, sys%gr%sb%dim
      call lr_init(psi_elec(idir))
      call lr_allocate(psi_elec(idir), sys%st, sys%gr%mesh)

      str_tmp =  em_wfs_tag(idir, 1)
      write(dirname,'(2a)') EM_RESP_DIR, trim(wfs_tag_sigma(str_tmp, sigma))
      call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, ierr, lr = psi_elec(idir))
      
      if(ierr.ne.0) then
        message(1) = "Could not load response wavefunctions from '"//trim(tmpdir)//dirname
        call write_fatal(1)
      end if

    end do

    call lr_init(psi_vib(1))
    call lr_allocate(psi_vib(1), sys%st, sys%gr%mesh)

    !
    ! Calculate the number of vibrational modes (here we should
    ! consider the possibility of linear molecules and not only two
    ! atoms).
    !
    if(sys%geo%natoms > 2) then 
      nnm = (sys%geo%natoms - 2)*sys%gr%sb%dim
    else
      nnm = (sys%geo%natoms - 2)*sys%gr%sb%dim + 1
    end if

    ! iterate over normal modes 
    do inm = 1, nnm
      !read vibronic perturbation
      call restart_read(trim(restart_dir)//VIB_MODES_DIR//trim(wfs_tag_sigma(phn_nm_wfs_tag(inm), 1)), &
           sys%st, sys%gr, sys%geo, ierr, lr = psi_vib(1))
      if(ierr .ne. 0) then
        message(1) = "Could not load vibrational response wavefunctions."
        call write_fatal(1)
      end if
    end do

    call lr_dealloc(psi_elec(1))
    call lr_dealloc(psi_vib(1))

    call states_deallocate_wfns(sys%st)

    POP_SUB(raman_run)
  end subroutine raman_run

end module raman_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
