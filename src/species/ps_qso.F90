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
!! $Id: ps_qso.F90 12895 2015-02-07 20:32:18Z dstrubbe $

#include "global.h"

module ps_qso_m
  use atomic_m
  use global_m
  use io_m
  use messages_m
  use profiling_m
  use ps_in_grid_m
  use xml_file_m

  implicit none

  private
  public ::     &
    ps_qso_t,       &
    ps_qso_init,    &
    ps_qso_end

  type ps_qso_t
    integer            :: atomic_number
    FLOAT              :: mass
    FLOAT              :: valence_charge
    integer            :: lmax
    integer            :: llocal
    FLOAT              :: mesh_spacing
    integer            :: grid_size
    FLOAT, allocatable :: potential(:, :)
    FLOAT, allocatable :: projector(:, :)
  end type ps_qso_t

contains

  ! ---------------------------------------------------------
  subroutine ps_qso_init(this, filename)
    type(ps_qso_t),   intent(inout) :: this
    character(len=*), intent(in)    :: filename

    character(len=256) :: filename2
    integer :: iunit, l, ll, size, ierr, ii
    logical :: found
    logical, allocatable :: found_l(:)    
    type(xml_file_t) :: qso_file
    type(xml_tag_t)  :: tag
    
    PUSH_SUB(ps_qso_init)

    ierr = xml_file_init(qso_file, trim(filename)// '.xml')

    if(ierr /= 0) ierr = xml_file_init(qso_file, trim(filename)// '.XML')

    if(ierr /= 0) ierr = xml_file_init(qso_file, trim(conf%share) // "/PP/qso/" // trim(filename) // ".xml")

    if(ierr /= 0) then
      call messages_write("Pseudopotential file '" // trim(filename) // ".xml' not found")
      call messages_fatal()
    end if    

    ierr = xml_file_get_tag_value(qso_file, 'lmax', this%lmax)
    ierr = xml_file_get_tag_value(qso_file, 'llocal', this%llocal)
    ierr = xml_file_get_tag_value(qso_file, 'mass', this%mass)
    ierr = xml_file_get_tag_value(qso_file, 'valence_charge', this%valence_charge)
    ierr = xml_file_get_tag_value(qso_file, 'mesh_spacing', this%mesh_spacing)

    do ii = 0, this%lmax
      call xml_file_tag(qso_file, 'projector', ii, tag)
      ierr = xml_tag_get_attribute_value(tag, 'l', ll)
      ierr = xml_tag_get_attribute_value(tag, 'size', size)

      if(ii == 0) then
        this%grid_size = size
        SAFE_ALLOCATE(this%potential(1:size, 0:this%lmax))
        SAFE_ALLOCATE(this%projector(1:size, 0:this%lmax))
      else
        ASSERT(size == this%grid_size)
      end if
      
      call xml_tag_get_tag_value_array(tag, 'radial_potential', size, this%potential(:, ll))
      call xml_tag_get_tag_value_array(tag, 'radial_function', size, this%projector(:, ll))

      call xml_tag_end(tag)

    end do

    call messages_write('QSO pseudopotential for '//trim(filename)//':', new_line = .true.)
    call messages_write('  maximum angular momentum component = ')
    call messages_write(this%lmax, new_line = .true.)
    call messages_write('  local angular momentum component   = ')
    call messages_write(this%llocal, new_line = .true.)
    call messages_info()
    
    call xml_file_end(qso_file)

    POP_SUB(ps_qso_init)
  end subroutine ps_qso_init

  ! ---------------------------------------------------------
  subroutine ps_qso_end(this)
    type(ps_qso_t), intent(inout) :: this

    PUSH_SUB(ps_qso_end)

    SAFE_DEALLOCATE_A(this%potential)
    SAFE_DEALLOCATE_A(this%projector)

    POP_SUB(ps_qso_end)
  end subroutine ps_qso_end

end module ps_qso_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
