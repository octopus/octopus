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
!! $Id$

#include "global.h"

module vibrations_m
  use geometry_m
  use global_m
  use io_m
  use lalg_adv_m
  use messages_m
  use profiling_m
  use simul_box_m
  use species_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public :: &
       vibrations_t, &
       vibrations_init, &
       vibrations_end,  &
       vibrations_normalize_dyn_matrix, &
       vibrations_diag_dyn_matrix, &
       vibrations_get_index, &
       vibrations_get_atom,  &
       vibrations_get_dir,   &
       vibrations_output
  
  type vibrations_t
    integer :: num_modes
    integer :: ndim
    integer :: natoms
    FLOAT, pointer :: dyn_matrix(:,:), normal_mode(:,:), freq(:)

    FLOAT :: disp
    FLOAT :: total_mass
  end type vibrations_t

contains

  ! ---------------------------------------------------------
  subroutine vibrations_init(this, geo, sb)
    type(vibrations_t),     intent(out) :: this
    type(geometry_t),    intent(inout) :: geo
    type(simul_box_t),   intent(inout) :: sb

    integer :: iatom

    PUSH_SUB(vibrations_init)

    this%ndim = sb%dim
    this%natoms = geo%natoms
    this%num_modes = geo%natoms*sb%dim
    SAFE_ALLOCATE(this%dyn_matrix(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%normal_mode(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%freq(1:this%num_modes))

    this%total_mass = M_ZERO
    do iatom = 1, geo%natoms
      this%total_mass = this%total_mass + species_weight(geo%atom(iatom)%spec)
    end do

    POP_SUB(vibrations_init)
  end subroutine vibrations_init


  ! ---------------------------------------------------------
  subroutine vibrations_end(this)
    type(vibrations_t),     intent(inout) :: this

    PUSH_SUB(vibrations_end)

    SAFE_DEALLOCATE_P(this%dyn_matrix)
    SAFE_DEALLOCATE_P(this%freq)
    SAFE_DEALLOCATE_P(this%normal_mode)

    POP_SUB(vibrations_end)
  end subroutine vibrations_end

  ! ---------------------------------------------------------

  subroutine vibrations_normalize_dyn_matrix(this, geo)
    type(vibrations_t),      intent(inout) :: this
    type(geometry_t),        intent(inout) :: geo

    FLOAT :: factor
    integer :: iatom, idir, jatom, jdir, imat, jmat

    PUSH_SUB(vibrations_normalize_dyn_matrix)

    do iatom = 1, this%natoms
      do idir = 1, this%ndim
        
        imat = vibrations_get_index(this, iatom, idir)

        do jatom = 1, this%natoms
          do jdir = 1, this%ndim
            
            jmat = vibrations_get_index(this, jatom, jdir)

            factor = this%total_mass/sqrt(species_weight(geo%atom(iatom)%spec)) / &
              sqrt(species_weight(geo%atom(jatom)%spec))
            
            this%dyn_matrix(imat, jmat) = this%dyn_matrix(imat, jmat) * factor

          end do
        end do

      end do
    end do

    POP_SUB(vibrations_normalize_dyn_matrix)
  end subroutine vibrations_normalize_dyn_matrix


  ! ---------------------------------------------------------
  subroutine vibrations_diag_dyn_matrix(this)
    type(vibrations_t),      intent(inout) :: this
    
    PUSH_SUB(vibrations_diag_dyn_matrix)

    this%normal_mode = M_ZERO
    
    this%normal_mode = this%dyn_matrix
    call lalg_eigensolve(this%num_modes, this%normal_mode, this%freq)

    this%freq(1:this%num_modes) = this%freq(1:this%num_modes) / this%total_mass

    POP_SUB(vibrations_diag_dyn_matrix)
  end subroutine vibrations_diag_dyn_matrix


  ! ---------------------------------------------------------
  integer pure function vibrations_get_index(this, iatom, idim)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: iatom
    integer,            intent(in) :: idim

    vibrations_get_index = (iatom - 1)*this%ndim + idim
  end function vibrations_get_index


  ! ---------------------------------------------------------
  integer pure function vibrations_get_atom(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_atom = 1 + (index - 1)/ this%ndim 
  end function vibrations_get_atom


  ! ---------------------------------------------------------
  integer pure function vibrations_get_dir(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_dir =  1 + mod(index - 1, this%ndim)
  end function vibrations_get_dir


  ! ---------------------------------------------------------
  subroutine vibrations_output(this, suffix)
    type(vibrations_t),   intent(in) :: this
    character (len=*),    intent(in) :: suffix
    
    integer :: iunit, i, j, iatom, jatom, idir, jdir, imat, jmat

    PUSH_SUB(vibrations_output)

    ! create directory for output
    call io_mkdir(VIB_MODES_DIR)

    ! output dynamic matrix
    iunit = io_open(VIB_MODES_DIR//'dynamical_matrix'//trim(suffix), action='write')

    do iatom = 1, this%natoms
      do idir = 1, this%ndim
        imat = vibrations_get_index(this, iatom, idir)

        do jatom = 1, this%natoms
          do jdir = 1, this%ndim
            jmat = vibrations_get_index(this, jatom, jdir)
            write(iunit, '(i6, i3, i6, i3, f16.4)') iatom, idir, jatom, jdir, &
                 units_from_atomic(unit_invcm, this%dyn_matrix(imat, jmat))
          end do
        end do

      end do
    end do

    call io_close(iunit)

    ! output frequencies and eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_frequencies'//trim(suffix), action='write')
    do i = 1, this%num_modes
      write(iunit, '(i6,f14.5)') i, units_from_atomic(unit_invcm, sqrt(abs(this%freq(i))))
    end do
    call io_close(iunit)

    ! output eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_modes'//trim(suffix), action='write')
    do i = 1, this%num_modes
      write(iunit, '(i6)', advance='no') i
      do j = 1, this%num_modes
        write(iunit, '(es14.5)', advance='no') this%normal_mode(j, i)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)

    POP_SUB(vibrations_output)
  end subroutine vibrations_output

end module vibrations_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
