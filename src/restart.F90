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
!!
!! $Id$

#include "global.h"

module restart
  use lib_oct_parser
  use global
  use messages
  use syslabels
  use io
  use states
  use mesh
  use grid
  use output
#ifdef HAVE_MPI
  use mpi_mod
#endif
#if defined(HAVE_MPI) && !defined(MPI_H)
  use mpi
#endif

  implicit none


#if defined(HAVE_MPI) && defined(MPI_H)
# include "mpif.h"
#endif


  private
  public ::                            &
       restart_init, clean_stop,       &
       drestart_write, zrestart_write, &
       drestart_read, zrestart_read,   &
       restart_format, restart_look

  integer, parameter ::    &
       RESTART_PLAIN  = 1, &
       RESTART_NETCDF = 2
  integer :: restart_format


contains

  ! returns true if a file named stop exists
  function clean_stop()
    logical clean_stop, file_exists

    clean_stop = .false.
    inquire(file='stop', exist=file_exists)
    if(file_exists) then
       message(1) = 'Clean STOP'
       message(2) = "(don't forget to remove the file 'stop' ;)"
       call write_warning(2)
       clean_stop = .true.
    end if

    return
  end function clean_stop


  subroutine restart_init
    integer :: i

    ! read restart format information
    call loct_parse_int(check_inp('RestartFileFormat'), RESTART_PLAIN, i)
    if (i<RESTART_PLAIN .or. i>RESTART_NETCDF) then
       write(message(1),'(a,i4,a)') "Input: '", i,"' is not a valid RestartFileFormat"
       message(2) = '(RestartFileFormat = plain | netcdf)'
       call write_fatal(2)
    end if

    ! Fix the restart format...
    restart_format = output_fill_how("Plain")
#if defined(HAVE_NETCDF)
    if(i == RESTART_NETCDF) then
       restart_format = output_fill_how("NETCDF")
    end if
#endif

  end subroutine restart_init

subroutine restart_look (dir, m, kpoints, dim, nst, ierr)
  character(len=*),  intent(in)    :: dir
  type(mesh_type),   intent(in)    :: m
  integer, intent(out) :: kpoints, dim, nst, ierr

  character(len=256)   :: line
  character(len=12)    :: filename
  character(len=1)     :: char
  integer :: iunit, iunit2, err, i, ist, idim, ik
  FLOAT :: occ, eigenval

  ierr = 0
  iunit  = iopar_open(m, trim(dir)//'/wfns', action='read', status='old', die=.false.)
  if(iunit < 0) then
      ierr = -1
      return
  end if
  iunit2 = iopar_open(m, trim(dir)//'/occs', action='read', status='old', die=.false.)
  if(iunit2 < 0) then
     call iopar_close(m, iunit)
     ierr = -1
     return
  end if

  ! Skip two lines.
  call iopar_read(m, iunit, line, err); call iopar_read(m, iunit, line, err)
  call iopar_read(m, iunit2, line, err); call iopar_read(m, iunit2, line, err)

  kpoints = 1
  dim = 1
  nst = 1
  do
    call iopar_read(m, iunit, line, i)
    read(line, '(a)') char
    if(i.ne.0.or.char=='%') exit
    read(line, *) ik, char, ist, char, idim, char, filename
    if(ik > kpoints) kpoints = ik
    if(idim == 2)    dim     = 2
    if(ist>nst)      nst     = ist
    call iopar_read(m, iunit2, line, err)
    read(line, *) occ, char, eigenval
  end do

  call iopar_close(m, iunit)
  call iopar_close(m, iunit2)
end subroutine restart_look



#include "undef.F90"
#include "real.F90"
#include "restart_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "restart_inc.F90"

end module restart
