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

module phonons_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use lib_adv_alg_m
  use mesh_m
  use messages_m
  use output_m
  use system_m

  implicit none

  private
  public :: &
       phonons_t, &
       phonons_init, &
       phonons_end,  &
       phonons_normalize_dm, &
       phonons_diagonalize_dm, &
       phn_idx, &
       phonons_output
  
  type phonons_t
    integer :: dim
    integer :: ndim
    integer :: natoms
    FLOAT, pointer :: dm(:,:), vec(:,:), freq(:)

    FLOAT :: disp
    FLOAT :: total_mass
  end type phonons_t

contains

  ! ---------------------------------------------------------
  subroutine phonons_init(ph, sys)
    type(phonons_t),     intent(out) :: ph
    type(system_t),      intent(inout) :: sys

    integer :: iatom

    ph%ndim = sys%gr%sb%dim
    ph%natoms = sys%geo%natoms
    ph%dim = sys%geo%natoms*sys%gr%sb%dim
    ALLOCATE(ph%dm(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%vec(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%freq(ph%dim), ph%dim)

    ph%total_mass = M_ZERO
    do iatom = 1, sys%geo%natoms
      ph%total_mass = ph%total_mass + sys%geo%atom(iatom)%spec%weight
    end do

  end subroutine phonons_init


  ! ---------------------------------------------------------
  subroutine phonons_end(ph)
    type(phonons_t),     intent(inout) :: ph

    deallocate(ph%dm)
    deallocate(ph%freq)
    deallocate(ph%vec)

  end subroutine phonons_end

  ! ---------------------------------------------------------

  subroutine phonons_normalize_dm(ph, geo)
    type(phonons_t),      intent(inout) :: ph
    type(geometry_t),     intent(inout) :: geo

    FLOAT :: factor
    integer :: iatom, idir, jatom, jdir

    do iatom = 1, ph%natoms
      do idir = 1, ph%ndim
        
        do jatom = 1, ph%natoms
          do jdir = 1, ph%ndim

            factor = ph%total_mass/sqrt(geo%atom(iatom)%spec%weight)/sqrt(geo%atom(jatom)%spec%weight)
            
            ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir)) = &
                 ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir)) * factor

          end do
        end do

      end do
    end do

  end subroutine phonons_normalize_dm

  subroutine phonons_diagonalize_dm(ph)
    type(phonons_t),      intent(inout) :: ph
    
    ph%vec = M_ZERO
    
    ! diagonalize DM
    call lalg_eigensolve(ph%dim, ph%dm, ph%vec, ph%freq)

    ph%freq(1:ph%dim) = ph%freq(1:ph%dim) / ph%total_mass

  end subroutine phonons_diagonalize_dm

  integer function phn_idx(ph, iatom, idim)
    type(phonons_t), intent(in) :: ph
    integer,         intent(in) :: iatom, idim
    phn_idx = (iatom-1)*ph%ndim + idim
  end function phn_idx

  subroutine phonons_output(ph, suffix)
    type(phonons_t),   intent(in) :: ph
    character (len=*), intent(in) :: suffix
    
    integer :: iunit, i, j, iatom, jatom, idir, jdir


    ! create directory for output
    call io_mkdir('phonons')

    ! output dynamic matrix
    call io_assign(iunit)
    iunit = io_open('phonons/DM'//trim(suffix), action='write')

    do iatom = 1, ph%natoms
      do idir = 1, ph%ndim

        do jatom = 1, ph%natoms
          do jdir = 1, ph%ndim
            
            write(iunit, '(f14.3)', advance='no') &
                 ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir)) * 219474.63 ! output cm^-1
          end do
        end do
        write(iunit, '(1x)')
      end do
    end do

    call io_close(iunit)

    ! output phonon frequencies and eigenvectors
    iunit = io_open('phonons/freq'//trim(suffix), action='write')
    do i = 1, ph%dim
      write(iunit, *) i, sqrt(abs(ph%freq(i))) * 219474.63 ! output cm^-1
    end do
    call io_close(iunit)

    ! output phonon eigenvectors
    iunit = io_open('phonons/vec'//trim(suffix), action='write')
    do i = 1, ph%dim
      write(iunit, '(i6)', advance='no') i
      do j = 1, ph%dim
        write(iunit, '(es14.5)', advance='no') ph%vec(j, i)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
    
  end subroutine phonons_output

end module phonons_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
