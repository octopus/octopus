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
  use atomic_orbital_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use io_function_oct_m
  use lda_u_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use restart_oct_m
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
       lda_u_write_U,                   &
       lda_u_write_magnetization,       &
       lda_u_load,                      &
       lda_u_dump

contains

 !> Prints the occupation matrices at the end of the scf calculation.
 subroutine lda_u_write_occupation_matrices(dir, this, st)
   type(lda_u_t),     intent(in)    :: this
   character(len=*),  intent(in)    :: dir
   type(states_t),    intent(in)    :: st

   integer :: iunit, ios, ispin, im, imp
 
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

   if(this%level == DFT_U_ACBN0) then
     iunit = io_open(trim(dir) // "/renorm_occ_matrices", action='write')
     write(iunit,'(a)') ' Renormalized occupation matrices '

     do ios = 1, this%norbsets
       do ispin = 1,st%d%nspin
          write(iunit,'(a, i4, a, i4)') ' Orbital set ', ios, ' spin ', ispin
          do im = 1, this%orbsets(ios)%norbs
            write(iunit,'(1x)',advance='no')

            if (states_are_real(st)) then
              do imp = 1, this%orbsets(ios)%norbs-1
                write(iunit,'(f14.8)',advance='no') this%dn_alt(im,imp,ispin,ios)
              end do
              write(iunit,'(f14.8)') this%dn_alt(im,this%orbsets(ios)%norbs,ispin,ios)
            else
              do imp = 1, this%orbsets(ios)%norbs-1
                write(iunit,'(f14.8,f14.8)',advance='no') this%zn_alt(im,imp,ispin,ios)
              end do
              write(iunit,'(f14.8,f14.8)') this%zn_alt(im,this%orbsets(ios)%norbs,ispin,ios)
            end if
          end do
       end do !ispin
     end do !iatom
     call io_close(iunit)
   end if

   end if

   POP_SUB(lda_u_write_occupation_matrices)
 end subroutine lda_u_write_occupation_matrices

 !--------------------------------------------------------- 
 subroutine lda_u_write_effectiveU(dir, this)
   type(lda_u_t),     intent(in)    :: this
   character(len=*),  intent(in)    :: dir

   integer :: iunit, ios

   PUSH_SUB(lda_u_write_effectiveU)

   if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
     iunit = io_open(trim(dir) // "/effectiveU", action='write')
     call lda_u_write_U(this, iunit)

     write(iunit, '(a,a,a,f7.3,a)') 'Hubbard U [', &
       trim(units_abbrev(units_out%energy)),']:'
     write(iunit,'(a,6x,14x,a)') ' Orbital',  'U'
     do ios = 1, this%norbsets
       if(this%orbsets(ios)%ndim == 1) then
         if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
         else
           write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
         end if
      else
        if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                        int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
         else
           write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), &
                         int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Ubar)
         end if
      end if
     end do


     write(iunit, '(a,a,a,f7.3,a)') 'Hund J [', &
       trim(units_abbrev(units_out%energy)),']:'
     write(iunit,'(a,6x,14x,a)') ' Orbital',  'J'
     do ios = 1, this%norbsets
       if(this%orbsets(ios)%ndim == 1) then
         if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
         else
           write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
         end if
      else
        if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                        int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
         else
           write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), &
                         int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Jbar)
         end if
      end if
     end do

     call io_close(iunit)
   end if

   POP_SUB(lda_u_write_effectiveU)
 end subroutine lda_u_write_effectiveU


 !--------------------------------------------------------- 
 subroutine lda_u_write_magnetization(dir, this, geo, mesh, st)
   type(lda_u_t),     intent(in)    :: this
   character(len=*),  intent(in)    :: dir
   type(geometry_t),  intent(in)    :: geo
   type(mesh_t),      intent(in)    :: mesh
   type(states_t),    intent(in)    :: st

   integer :: iunit, ia, ios, im
   FLOAT, allocatable :: mm(:,:)

   if( .not. mpi_grp_is_root(mpi_world)) return

   PUSH_SUB(lda_u_write_magnetization)

   call io_mkdir(dir)
    iunit = io_open(trim(dir)//"/magnetization.xsf", action='write', position='asis')

    if(this%nspins > 1) then
      SAFE_ALLOCATE(mm(1:geo%natoms, 1:mesh%sb%dim))
      mm(1:geo%natoms, 1:mesh%sb%dim) = M_ZERO
      !We compute the magnetization vector for each orbital set
      do ios = 1, this%norbsets
        ia = this%orbsets(ios)%iatom
        do im = 1, this%orbsets(ios)%norbs
          if (states_are_real(st)) then
            mm(ia, 3) = mm(ia, 3) + this%dn(im,im,1,ios) - this%dn(im,im,2,ios) 
          else
            mm(ia, 3) = mm(ia, 3) + real(this%zn(im,im,1,ios) - this%zn(im,im,2,ios))
            !Spinors
            if(this%nspins /= this%spin_channels) then
              mm(ia, 1) = mm(ia, 1) + 2*real(this%zn(im,im,3,ios))
              mm(ia, 2) = mm(ia, 2) - 2*aimag(this%zn(im,im,3,ios))
            end if
          end if  
        end do !im
      end do ! ios
      call write_xsf_geometry(iunit, geo, mesh, forces = mm)
      SAFE_DEALLOCATE_A(mm)
    end if

    call io_close(iunit)

   POP_SUB(lda_u_write_magnetization)
 end subroutine lda_u_write_magnetization

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
       if(this%orbsets(ios)%ndim == 1) then 
         if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)  
         else
           write(iunit,'(i4,a10, 3x, a1, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
         end if
      else
        if(this%orbsets(ios)%nn /= 0 ) then
           write(iunit,'(i4,a10, 2x, i1, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                        this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), &
                        int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
         else
           write(iunit,'(i4,a10, 3x, a1, i1, a2, f15.6)') ios, trim(species_label(this%orbsets(ios)%spec)), &
                                l_notation(this%orbsets(ios)%ll), &
                         int(M_TWO*(this%orbsets(ios)%jj)), '/2', units_from_atomic(units_out%energy, this%orbsets(ios)%Ueff)
         end if 
      end if
     end do
   end if
 
   POP_SUB(lda_u_write_U)
 end subroutine lda_u_write_U

  ! ---------------------------------------------------------
  subroutine lda_u_dump(restart, this, st, ierr, iter)
    type(restart_t),      intent(in)  :: restart
    type(lda_u_t),        intent(in)  :: this
    type(states_t),       intent(in)  :: st
    integer,              intent(out) :: ierr
    integer, optional,    intent(in)  :: iter

    integer :: err, occsize
    FLOAT, allocatable :: Ueff(:), docc(:)
    CMPLX, allocatable :: zocc(:)

    PUSH_SUB(lda_u_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(lda_u_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing LDA+U restart."
      call messages_info(1)
    end if

    occsize = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
    if(this%level == DFT_U_ACBN0) occsize = occsize*2
 
    if (states_are_real(st)) then
      SAFE_ALLOCATE(docc(1:occsize))
      docc = M_ZERO
      call dlda_u_get_occupations(this, docc)
      call drestart_write_binary(restart, "lda_u_occ", occsize, docc, err)      
      if (err /= 0) ierr = ierr + 1
      SAFE_DEALLOCATE_A(docc)
    else
      SAFE_ALLOCATE(zocc(1:occsize))
      zocc = M_ZERO
      call zlda_u_get_occupations(this, zocc)
      call zrestart_write_binary(restart, "lda_u_occ", occsize, zocc, err)
      if (err /= 0) ierr = ierr + 1
      SAFE_DEALLOCATE_A(zocc)
    end if


    if(this%level == DFT_U_ACBN0) then
      SAFE_ALLOCATE(Ueff(1:this%norbsets))
      Ueff = M_ZERO
      call lda_u_get_effectiveU(this, Ueff(:))
      call drestart_write_binary(restart, "lda_u_Ueff", this%norbsets, Ueff, err)
      SAFE_DEALLOCATE_A(Ueff)
      if (err /= 0) ierr = ierr + 1
    end if

    if (debug%info) then
      message(1) = "Debug: Writing LDA+U restart done."
      call messages_info(1)
    end if

    POP_SUB(lda_u_dump)
  end subroutine lda_u_dump


 ! ---------------------------------------------------------
  subroutine lda_u_load(restart, this, st, ierr)
    type(restart_t),      intent(in)    :: restart
    type(lda_u_t),        intent(inout) :: this
    type(states_t),       intent(in)    :: st
    integer,              intent(out)   :: ierr

    integer :: err, occsize
    FLOAT, allocatable :: Ueff(:), docc(:)
    CMPLX, allocatable :: zocc(:)

    PUSH_SUB(lda_u_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(lda_u_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading LDA+U restart."
      call messages_info(1)
    end if

    !We have to read the effective U first, as we call lda_u_uptade_potential latter
    if(this%level == DFT_U_ACBN0) then
      SAFE_ALLOCATE(Ueff(1:this%norbsets))
      call drestart_read_binary(restart, "lda_u_Ueff", this%norbsets, Ueff, err)
      if (err /= 0) ierr = ierr + 1
      call lda_u_set_effectiveU(this, Ueff)
      SAFE_DEALLOCATE_A(Ueff)
    end if


    occsize = this%maxnorbs*this%maxnorbs*this%nspins*this%norbsets
    if(this%level == DFT_U_ACBN0) occsize = occsize*2

    if (states_are_real(st)) then
      SAFE_ALLOCATE(docc(1:occsize))
      call drestart_read_binary(restart, "lda_u_occ", occsize, docc, err) 
      if (err /= 0) ierr = ierr + 1
      call dlda_u_set_occupations(this, docc)
      call dlda_u_update_potential(this, st)
      SAFE_DEALLOCATE_A(docc)
    else
      SAFE_ALLOCATE(zocc(1:occsize))
      call zrestart_read_binary(restart, "lda_u_occ", occsize, zocc, err)
      if (err /= 0) ierr = ierr + 1
      call zlda_u_set_occupations(this, zocc)
      call zlda_u_update_potential(this, st)
      SAFE_DEALLOCATE_A(zocc)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading LDA+U restart done."
      call messages_info(1)
    end if

    POP_SUB(lda_u_load)
  end subroutine lda_u_load

end module lda_u_io_oct_m
