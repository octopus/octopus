!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$

#include "global.h"

module lda_u_io_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use lda_u_oct_m
  use messages_oct_m
  use mpi_oct_m
  use species_oct_m
  use states_oct_m
  use types_oct_m  
  use unit_oct_m
  use unit_system_oct_m
 
  implicit none

  private

  public ::                             &
       lda_u_write_occupation_matrices, &
       lda_u_write_effectiveU,          &
       lda_u_write_U

  character(len=1), parameter :: &
    l_notation(0:3) = (/ 's', 'p', 'd', 'f' /)

contains

 !> Prints the occupation matrices at the end of the scf calculation.
 subroutine lda_u_write_occupation_matrices(dir, this, st)
   type(lda_u_t),     intent(in)    :: this
   character(len=*),  intent(in)    :: dir
   type(states_t),    intent(in)    :: st

   integer :: iunit, ios, ispin, im, imp
   FLOAT :: hubbardl
 
   PUSH_SUB(lda_u_write_occupation_matrices)

   if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
   iunit = io_open(trim(dir) // "/occ_matrices", action='write')
   write(iunit,'(a)') ' Occupation matrices '

   do ios = 1, this%norbsets
     do ispin = 1,st%d%nspin 
        write(iunit,'(a, i4, a, i4)') ' Orbital set ', ios, ' spin ', ispin
        do im = 1, this%orbsets(ios)%norbs
          write(iunit,'(1x)',advance='no') 

          if (states_are_real(st)) then
            do imp = 1, this%orbsets(ios)%norbs-1
              write(iunit,'(f14.8)',advance='no') this%dn(im,imp,ispin,ios)  
            end do
            write(iunit,'(f14.8)') this%dn(im,this%orbsets(ios)%norbs,ispin,ios)
          else
            do imp = 1, this%orbsets(ios)%norbs-1
              write(iunit,'(f14.8,f14.8)',advance='no') this%zn(im,imp,ispin,ios)
            end do
            write(iunit,'(f14.8,f14.8)') this%zn(im,this%orbsets(ios)%norbs,ispin,ios) 
          end if
        end do
     end do !ispin
   end do !iatom
   call io_close(iunit)

   end if

   POP_SUB(lda_u_write_occupation_matrices)
 end subroutine lda_u_write_occupation_matrices

 !--------------------------------------------------------- 
 subroutine lda_u_write_effectiveU(dir, this, st)
   type(lda_u_t),     intent(in)    :: this
   character(len=*),  intent(in)    :: dir
   type(states_t),    intent(in)    :: st

   integer :: iunit

   PUSH_SUB(lda_u_write_effectiveU)

   if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
     iunit = io_open(trim(dir) // "/effectiveU", action='write')
     call lda_u_write_U(this, iunit)
     call io_close(iunit)
   end if

   POP_SUB(lda_u_write_effectiveU)
 end subroutine lda_u_write_effectiveU


 !--------------------------------------------------------- 
 subroutine lda_u_write_U(this, iunit)
   type(lda_u_t),     intent(in) :: this
   integer,           intent(in) :: iunit

   integer :: ios

   PUSH_SUB(lda_u_write_U)
  
   if(mpi_grp_is_root(mpi_world)) then

     write(iunit, '(a,a,a,f7.3,a)') 'Effective Hubbard U [', &
       trim(units_abbrev(units_out%energy)),']:'
     write(iunit,'(a,6x,14x,a)') ' Orbital',  'U'
     do ios = 1, this%norbsets
      if(this%orbsets(ios)%nn /= 0 ) then
        write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                      this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), this%orbsets(ios)%Ueff  
      else
        write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                              l_notation(this%orbsets(ios)%ll), this%orbsets(ios)%Ueff
      end if
     end do
   end if
 
   POP_SUB(lda_u_write_U)
 end subroutine lda_u_write_U

end module lda_u_io_oct_m
