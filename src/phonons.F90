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
  use datasets_m
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lib_adv_alg_m
  use lib_oct_parser_m
  use mesh_m
  use messages_m
  use output_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use units_m
  use v_ks_m

  implicit none

  private
  public :: &
       phonons_t, &
       phonons_init, &
       phonons_end,  &
       phonons_diagonalize_dm, &
       phonons_index, &
       phonons_output
  
  type phonons_t
    integer :: dim
    integer :: ndim
    integer :: natoms
    FLOAT, pointer :: dm(:,:), vec(:,:), freq(:)

    FLOAT :: disp
  end type phonons_t

contains

  ! ---------------------------------------------------------
  subroutine phonons_init(ph, sys)
    type(phonons_t),     intent(out) :: ph
    type(system_t),      intent(inout) :: sys

    ph%ndim = sys%gr%sb%dim
    ph%natoms = sys%geo%natoms
    ph%dim = sys%geo%natoms*sys%gr%sb%dim
    ALLOCATE(ph%dm(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%vec(ph%dim, ph%dim), ph%dim*ph%dim)
    ALLOCATE(ph%freq(ph%dim), ph%dim)

  end subroutine phonons_init


  ! ---------------------------------------------------------
  subroutine phonons_end(ph)
    type(phonons_t),     intent(inout) :: ph

    deallocate(ph%dm)
    deallocate(ph%freq)

  end subroutine phonons_end

  ! ---------------------------------------------------------
  subroutine phonons_diagonalize_dm(ph)
    type(phonons_t),      intent(inout) :: ph
    
    ph%vec = M_ZERO
    
    ! diagonalize DM
    call lalg_eigensolve(ph%dim, ph%dm, ph%vec, ph%freq)
    
  end subroutine phonons_diagonalize_dm

  integer function phonons_index(ph, iatom, idim)
    type(phonons_t), intent(in) :: ph
    integer,         intent(in) :: iatom, idim
    phonons_index = (iatom-1)*ph%ndim + idim
  end function phonons_index

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
            
            write(iunit, '(es14.5)', advance='no') &
                 ph%dm(phonons_index(ph, iatom, idir), phonons_index(ph, jatom, jdir))
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
