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

module ps_fhi_file_oct_m
  use global_oct_m
  use messages_oct_m
  use namespace_oct_m

  implicit none

  private

  public ::                 &
    ps_fhi_file_t,          &
    ps_fhi_file_read

  ! First, the contents of the file.
  type ps_fhi_file_t
    ! Components are public by default

    ! This is the general ABINIT header
    character(len=256), private :: title
    FLOAT,              private :: znucl  ! charge of the nucleus
    FLOAT,              private :: zion   ! valence charge
    integer,            private :: pspdat ! date of creation of PP (DDMMYY)
    integer,            private :: pspcod ! code for the pseudopotential (6 for .fhi)
    integer,            private :: pspxc  ! exchange-correlation used to generate the psp
    integer                     :: lmax   ! Maximum l to use
    integer                     :: lloc   ! local part of the pseudo to use
    integer,            private :: mmax   ! Maximum number of points in real space grid
    FLOAT,              private :: r2well ! ??

    ! this is specific for FHI
    FLOAT,              private :: rchrg  ! The core charge becomes zero beyond rchrg
    FLOAT,              private :: fchrg  !
    FLOAT,              private :: qchrg  !
  end type ps_fhi_file_t

contains

  ! ---------------------------------------------------------
  subroutine ps_fhi_file_read(unit, psf, namespace)
    integer,             intent(in)    :: unit
    type(ps_fhi_file_t), intent(inout) :: psf
    type(namespace_t),   intent(in)    :: namespace

    character(len=3) :: line

    PUSH_SUB(read_file_data)

    read(unit, *) psf%title
    read(unit, *) psf%znucl, psf%zion, psf%pspdat
    read(unit, *) psf%pspcod, psf%pspxc, psf%lmax, psf%lloc, psf%mmax, psf%r2well

    if(psf%pspcod /= 6) then
      message(1) = "Inconsistency in pseudopotential file:"
      write(message(2),'(a,i2)') "  expecting pspcod = 6, but found ", psf%pspcod
      call messages_fatal(2, namespace=namespace)
    end if

    read(unit, '(a3)') line
    if(line /= '4--') then  ! have non-local core corrections
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

end module ps_fhi_file_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
