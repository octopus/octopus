!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module restart
use global
use io
use states
use mesh
use output

implicit none

private
public :: restart_write, restart_load

integer, parameter :: RESTART_PLAIN  = 1, &
                      RESTART_NETCDF = 2

contains

subroutine restart_write(dir, st, m, ierr, iter, nvprev, vprev)
  character(len=*),  intent(in) :: dir
  type(states_type), intent(in) :: st
  type(mesh_type), intent(in)   :: m
  integer, intent(out)          :: ierr
  integer, intent(in), optional :: iter
  integer, intent(in), optional :: nvprev
  FLOAT, intent(in), optional   :: vprev(:, :, :) ! (m%np, st%d%nspin, nvprev)

  logical :: re
  integer :: iunit, err, ik, ist, idim, i, is, how
  character(len=6) :: filename, ext

  call push_sub('states_write')

  ierr = 0
  call loct_mkdir(dir)

  MPINODE0 then
    call io_assign(iunit)
    open(unit = iunit, file=trim(dir)//'/data', iostat = err, &
       form='formatted', action='write', status='replace', recl = 200)
    if(err .ne. 0) then
       ierr = 1; return
    endif
  endif

  re = .false.
  if(associated(st%dpsi)) then
    re = .true.
  elseif(.not.associated(st%zpsi)) then
    ierr = 2; return
  endif

  ! Fix the restart format...
  how = output_fill_how("Plain"); ext = " "
#if defined(HAVE_NETCDF)
  if(st%restart_format == RESTART_NETCDF) then
    how = output_fill_how("NETCDF"); ext =".ncdf"
  endif
#endif

  i = 1
  MPINODE0 write(iunit,'(a)') &
'#     #kpoint            #st            #dim    filename           eigenvalue[a.u.]'
  MPINODE0 write(iunit,'(a)') '%Wavefunctions'
  do ik = 1, st%d%nik
     do ist = 1, st%nst
        do idim = 1, st%d%dim
           write(filename,'(i6.6)') i
           MPINODE0 write(unit = iunit, fmt = *) ik, ' | ', ist, ' | ', idim, ' | "', &
                                                 trim(filename)//trim(ext), '" | ', st%eigenval(ist, ik)
           if(st%st_start <= ist .and. st%st_end >= ist) then
              if(re) then
                 call doutput_function(how, trim(dir), filename, m, st%dpsi(:, idim, ist, ik), M_ONE)
              else
                 call zoutput_function(how, trim(dir), filename, m, st%zpsi(:, idim, ist, ik), M_ONE)
              endif
           endif
           i = i + 1
        enddo
     enddo
  enddo
  MPINODE0 write(iunit,'(a)') '%'
  if(present(iter)) then
     MPINODE0 write(iunit,'(a,i5)') 'Iter = ', iter
  endif

  if(present(vprev) .and. mpiv%node == 0) then
    do i = 1, nvprev
       do is = 1, st%d%nspin
          write(filename,'(a1,i2.2,i3.3)') 'p',i, is
          call doutput_function(how, trim(dir), filename, m, vprev(1:m%np, is, i), M_ONE)
       enddo
    enddo
  endif

#if defined(HAVE_MPI)
  call mpi_barrier(MPI_COMM_WORLD, i) ! Since some processors did more than others...
#endif
  MPINODE0 call io_close(iunit)
  call pop_sub()
end subroutine restart_write

subroutine restart_load(dir, st, m, ierr, iter, nvprev, vprev)
  character(len=*),  intent(in)    :: dir
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(in)      :: m
  integer, intent(out)             :: ierr
  integer, intent(out), optional   :: iter
  integer, intent(in), optional    :: nvprev
  FLOAT, intent(inout), optional   :: vprev(:, :, :) ! (m%np, st%d%nspin, nvprev)

  logical :: re
  integer :: iunit, err, ik, ist, idim, i, is
  character(len=50) :: blockname
  character(len=12) :: filename
  character(len=1) :: char
  logical, allocatable :: filled(:, :, :)

  call push_sub('restart_load')

  ierr = 0
  call io_assign(iunit)
  open(unit = iunit, file=trim(dir)//'/data', iostat = err, &
       form='formatted', action='read', status='old', recl = 200)
  if(err .ne. 0) then
    ierr = 1; return
  endif

  re = .false.
  if(associated(st%dpsi)) then
    re = .true.
  elseif(.not.associated(st%zpsi)) then
    return
  endif

  allocate(filled(st%d%dim, st%st_start:st%st_end, st%d%nik)); filled = .false.

  read(iunit, *); read(iunit, *) ! Skip two lines...
  do
     read(unit = iunit, fmt = '(a)', iostat = i) char
     if(char=='%') exit
     backspace(unit = iunit)
     read(unit = iunit, iostat = i, fmt = *) ik, char, ist, char, idim, char, filename, char, st%eigenval(ist, ik)
     if(index_is_wrong()) cycle
     if(st%st_start <= ist .and. st%st_end >= ist) then
        if(re) then
            call dinput_function(trim(dir)//'/'//trim(filename), m, st%dpsi(:, idim, ist, ik), err)
        else
            call zinput_function(trim(dir)//'/'//trim(filename), m, st%zpsi(:, idim, ist, ik), err)
        endif
        if(err<=0) filled(idim, ist, ik) = .true.
     endif
  enddo

  if(present(iter)) then
    read(unit = iunit, fmt = *) filename, filename, iter
  endif
  if(present(vprev)) then
     do i = 1, nvprev
        do is = 1, st%d%nspin
           write(filename,'(a1,i2.2,i3.3)') 'p',i, is
           call dinput_function(trim(dir)//'/'//trim(filename), m, vprev(1:m%np, is, i), err)
           ! If problems, try with the netcdf files. However, it will always try to read first
           ! the plain files, which may be a problem if both are present, and the "good" one is the netcdf.
           ! This has to be fixed.
           if(err>0) call dinput_function(trim(dir)//'/'//trim(filename)//'.ncdf', m, vprev(1:m%np, is, i), err)
           if(err>0) then 
             write(message(1),'(a)') 'Problem reading potential in restart file for extrapolation'
             write(message(2),'(a,i3)') 'Error code = ',err
             call write_warning(2)
           endif
        enddo
     enddo
  endif

  if(any(.not.filled)) then
     call fill()
     ierr = -1
  endif

  deallocate(filled)
  call io_close(iunit)
  call pop_sub()

  contains

  subroutine fill() ! Put random function in orbitals that could not be read.
     do ik = 1, st%d%nik
        do ist = 1, st%nst
           do idim = 1, st%d%dim
              if(filled(idim, ist, ik)) cycle
              write(message(1),'(a,3i4)') 'Randomizing wavefunction: #dim, #ist #ik = ', idim, ist, ik
              call write_warning(1)
           enddo
        enddo
     enddo
  end subroutine fill
  logical function index_is_wrong() ! .true. if the index (idim, ist, ik) is not present in st structure...
    if(idim > st%d%dim .or. idim < 1 .or.   &
       ist  > st%nst   .or. ist  < 1 .or.   &
       ik   > st%d%nik .or. ik   < 1) then
      index_is_wrong = .true.
    else
      index_is_wrong = .false.
    endif
  end function index_is_wrong
end subroutine restart_load


end module restart
