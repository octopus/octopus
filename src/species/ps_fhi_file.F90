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

module ps_fhi_file_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::                 &
    ps_fhi_file_t,          &
    ps_fhi_file_read

  ! First, the contents of the file.
  type ps_fhi_file_t
    ! This is the general ABINIT header
    character(len=256) :: title
    FLOAT              :: znucl  ! charge of the nucleus
    FLOAT              :: zion   ! valence charge
    integer            :: pspdat ! date of creation of PP (DDMMYY)
    integer            :: pspcod ! code for the pseudopotential (6 for .fhi)
    integer            :: pspxc  ! exchange-correlation used to generate the psp
    integer            :: lmax   ! Maximum l to use
    integer            :: lloc   ! local part of the pseudo to use
    integer            :: mmax   ! Maximum number of points in real space grid
    FLOAT              :: r2well ! ??

    ! this is specific for FHI
    FLOAT              :: rchrg  ! The core charge becomes zero beyond rchrg
    FLOAT              :: fchrg  !
    FLOAT              :: qchrg  !
  end type ps_fhi_file_t
  
contains

  ! ---------------------------------------------------------
  subroutine ps_fhi_file_read(unit, psf)
    integer,             intent(in)    :: unit
    type(ps_fhi_file_t), intent(inout) :: psf

    character(len=3) :: line

    PUSH_SUB(read_file_data)

    read(unit, *) psf%title
    read(unit, *) psf%znucl, psf%zion, psf%pspdat
    read(unit, *) psf%pspcod, psf%pspxc, psf%lmax, psf%lloc, psf%mmax, psf%r2well

    if(psf%pspcod.ne.6) then
      message(1) = "Inconsistency in pseudopotential file:"
      write(message(2),'(a,i2)') "  expecting pspcod = 6, but found ", psf%pspcod
      call messages_fatal(2)
    end if
    
    read(unit, '(a3)') line
    if(line.ne.'4--') then  ! have non-local core corrections
      backspace(unit)
      read(unit, *) psf%rchrg, psf%fchrg, psf%qchrg
    else
      psf%rchrg = M_ZERO
      psf%fchrg = M_ZERO
      psf%qchrg = M_ZERO
    end if

    ! skip next 3 lines (#5, #6, #7)
    read(unit, '(a3)') line
    read(unit, '(a3)') line
    read(unit, '(a3)') line
    
    POP_SUB(read_file_data)
  end subroutine ps_fhi_file_read

end module ps_fhi_file_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
