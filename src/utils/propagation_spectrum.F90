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

program propagation_spectrum
  use command_line_m
  use datasets_m
  use global_m
  use io_m
  use loct_m
  use messages_m
  use parser_m
  use profiling_m
  use spectrum_m
  use unit_m
  use unit_system_m

  implicit none

  integer :: in_file(3), out_file(3), eq_axes, nspin, lmax, time_steps, ierr
  logical :: calculate_tensor
  type(spec_t) :: spectrum
  type(unit_system_t) :: file_units

  ! Initialize stuff
  call global_init()
  call getopt_init(ierr)
  if(ierr.eq.0) call getopt_propagation_spectrum
  call parser_init()
  call messages_init()

  call datasets_init(1)
  call io_init()

  call unit_system_init()

  call spectrum_init(spectrum)

  select case (spectrum%spectype)
    case (SPECTRUM_ABSORPTION)
      call read_files('multipoles')
      call calculate('cross_section')
    case (SPECTRUM_ENERGYLOSS)
      call calculate_ftchd('ftchd', 'dynamic_structure_factor')
    case default
      write(message(1), '(a)') 'No PropagationSpectrumType defined,'
      write(message(2), '(a)') 'cannot calculate the spectrum.'
      call messages_warning(2)
  end select

  call io_end()
  call datasets_end()
  call messages_end()
  call parser_end()
  call global_end()

  contains

    !/*----------------------------------------------------------------------------
    ! Two spectra can be calculated with these routines: photoabsorption spectrum
    ! or the dynamic struture factor. This is controlled by the input variable
    ! PropagationSpectrumType. In the photoabsorption case the calculation
    ! goes as follows: start the search for the file(s) that must be processed.
    ! In a future version of octopus, this could be controlled by a
    ! command-line option. The routine sets the eq_axes and calculate_tensor
    ! variables, as well as the in_file unit numbers.
    !
    ! The options are:
    ! (i)  A file called "multipoles" is found. In this case, no other file is
    !      considered.
    !      (i.1) If this file signals three equivalent axes, the full tensor is
    !            calculated, and placed in "cross_section_tensor".
    !      (i.2) If this file signals fewer than three equivalent axes, the full
    !            tensor cannot be calculated, and instead a "cross_section_vector"
    !            will be generated.
    ! (ii) A file called "multipoles.1" is found. In this case, the program will
    !      always try to generate the full tensor (calculate_tensor = .true.).
    !      The file "cross_section_tensor" should be generated at the end.
    !      Other files are searched for, depending on the equivalent axes that
    !      are written in the "multipoles.1" file.
    !      (ii.1) Three equivalent axes. No other file is searched for.
    !      (ii.2) Two equivalent axes. File "multipoles.2" is searched for; the
    !             program ends if it is not found.
    !      (ii.3) No equivalent axes. Files "multipoles.2" and "multipoles.3" are
    !             searched for; the program ends if they are not found.
    ! In the case of dynamic structure factor calculation the program looks for
    ! ftchds.cos and ftchds.sin. In the future versions it should also look for
    ! the files ftchds.lXX_mYY, which would be used to obtain the directionally
    ! averaged dynamic structure factor.
    !---------------------------------------------------------------------------*/!
    subroutine read_files(fname)
      character(len=*), intent(in) :: fname

      type(kick_t) :: kick
      FLOAT :: dt

      PUSH_SUB(read_files)
      
      in_file(1) = io_open(trim(fname), action='read', status='old', die=.false.)
      if(in_file(1) < 0) in_file(1) = io_open('td.general/'//trim(fname), action='read', status='old', die=.false.)
      if(in_file(1) >= 0) then
        write(message(1),'(3a)') 'File "', trim(fname), '" found. This will be the only file to be processed.'
        write(message(2),'(a)')  '(If more than one file is to be used, the files should be called'
        write(message(3),'(5a)') '"', trim(fname), '.1", "', trim(fname), '.2", etc.'
        write(message(4),'(a)')
        call messages_info(4)
        
        ! OK, so we only have one file. Now we have to see what it has inside.
        call spectrum_mult_info(in_file(1), nspin, kick, time_steps, dt, file_units, lmax=lmax)
        eq_axes = kick%pol_equiv_axes
        if(eq_axes == 3) then
          calculate_tensor = .true.
          write(message(1),'(3a)') 'The file "', trim(fname), '" tells me that the system has three equivalent axes.'
          write(message(2),'(a)')  'I will calculate the full tensor, written in file "XXXX_tensor".'
          call messages_info(2)
        else if(eq_axes == 2) then
          write(message(1),'(3a)') 'The file "', trim(fname), '" tells me that the system has two equivalent axes.'
          write(message(2),'(a)')  'However, I am only using this file; cannot calculate the full tensor.'
          write(message(3),'(a)')  'A file "XXXX_vector" will be generated instead.'
          call messages_warning(3)
          calculate_tensor = .false.
        else
          write(message(1),'(3a)') 'The file "', trim(fname), '" tells me that the system has no usable symmetry. '
          write(message(2),'(a)')  'However, I am only using this file; cannot calculate the full tensor.'
          write(message(3),'(a)')  'A file "XXXX_vector" will be generated instead.'
          call messages_warning(3)
          calculate_tensor = .false.
        end if
        
      else  ! We will try to load more fname.1 files...
        
        ! In this case, we will always want the full tensor
        calculate_tensor = .true.
        
        in_file(1) = io_open(trim(fname)//'.1', action='read', status='old', die=.false.)
        if(in_file(1) < 0) in_file(1) = io_open('td.general/'//trim(fname)//'.1', action='read', status='old', die=.false.)
        if(in_file(1) < 0) then ! Could not find proper files. Die and complain.
          write(message(1),'(5a)') 'No "', trim(fname), '" or "', trim(fname), '.1" file found. At least one of those'
          write(message(2),'(a)')  'should be visible.'
          call messages_fatal(2)
        end if
        
        call spectrum_mult_info(in_file(1), nspin, kick, time_steps, dt, file_units, lmax=lmax)
        eq_axes = kick%pol_equiv_axes
        
        if(eq_axes == 3) then
          write(message(1),'(3a)') 'The file "', trim(fname), '.1" tells me that the system has three equivalent axes.'
          write(message(2),'(a)') 'I will calculate the full tensor, written in file "cross_section_tensor".'
          call messages_info(2)
          
        else if(eq_axes == 2) then
          in_file(2) = io_open(trim(fname)//'.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) in_file(2) = io_open('td.general/'//trim(fname)//'.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) then
            write(message(1),'(3a)') 'The file "', trim(fname), '.1" tells me that the system has two equivalent axes,'
            write(message(2),'(3a)') 'but I cannot find a "', trim(fname), '.2".'
            call messages_fatal(2)
          end if
          write(message(1),'(5a)') 'Found two files, "', trim(fname), '.1" and "', trim(fname), '.2".'
          write(message(2),'(a)')  'Two polarization axes are equivalent. I will generate the full tensor.'
          call messages_info(2)
          
        else ! No equivalent axes
          in_file(2) = io_open(trim(fname)//'.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) in_file(2) = io_open('td.general/'//trim(fname)//'.2', action='read', status='old', die=.false.)
          if(in_file(2) < 0) then
            write(message(1),'(3a)') 'The file "', trim(fname), '.1" tells me that the system has three inequivalent axes,'
            write(message(2),'(3a)') 'but I cannot find a "', trim(fname), '.2".'
            call messages_fatal(2)
          end if
          in_file(3) = io_open(trim(fname)//'.3', action='read', status='old', die=.false.)
          if(in_file(3) < 0) in_file(3) = io_open('td.general/'//trim(fname)//'.3', action='read', status='old', die=.false.)
          if(in_file(3) < 0) then
            write(message(1),'(3a)') 'The file "', trim(fname), '.1" tells me that the system has three inequivalent axes,'
            write(message(2),'(3a)') 'but I cannot find a "', trim(fname), '.3".'
            call messages_fatal(2)
          end if
          write(message(1),'(7a)') 'Found three files, "', trim(fname), '.1", "', trim(fname), '.2" and "', trim(fname), '.3".'
          write(message(2),'(a)')  'No symmetry information will be used.'
          call messages_info(2)
        end if
        
      end if

      POP_SUB(read_files)
    end subroutine read_files


    !----------------------------------------------------------------------------   
    subroutine calculate(fname)
      character(len=*), intent(in) :: fname

      integer :: ii, jj
      character(len=150), allocatable :: filename(:)

      PUSH_SUB(calculate)

      if(.not.calculate_tensor) then

        out_file(1) = io_open(trim(fname)//'_vector', action='write')
        call spectrum_cross_section(in_file(1), out_file(1), spectrum)
        call io_close(in_file(1))
        call io_close(out_file(1))
        
      else
        select case(eq_axes)
        case(0, 1); jj = 3
        case(2);    jj = 2
        case(3);    jj = 1
        end select

        SAFE_ALLOCATE(filename(1:jj))
        do ii = 1, jj
          write(filename(ii),'(2a,i1)') trim(fname), '_vector.',ii
          out_file(ii) = io_open(trim(filename(ii)), action='write')
          call spectrum_cross_section(in_file(ii), out_file(ii), spectrum)
          call io_close(in_file(ii))
          call io_close(out_file(ii))
          in_file(ii)  = io_open(trim(filename(ii)), action='read', status='old')
        end do

        out_file(1) = io_open(trim(fname)//'_tensor', action='write')
        call spectrum_cross_section_tensor(spectrum, out_file(1), in_file(1:jj))
        do ii = 1, jj
          call io_close(in_file(ii))
        end do
        call io_close(out_file(1))
        
      end if

      POP_SUB(calculate)
    end subroutine calculate


    !----------------------------------------------------------------------------   
    subroutine calculate_ftchd(fname_in, fname_out)
      character(len=*), intent(in) :: fname_in, fname_out

      PUSH_SUB(calculate_ftchd)

      ! read files
      in_file(1) = io_open(trim(fname_in) // '.sin', action='read', status='old', die=.false.)
      in_file(2) = io_open(trim(fname_in) // '.cos', action='read', status='old', die=.false.)

      out_file(1) = io_open(trim(fname_out), action='write')
      call spectrum_dyn_structure_factor(in_file(1), in_file(2), out_file(1), spectrum)
      call io_close(in_file(1))
      call io_close(out_file(1))

      POP_SUB(calculate_ftchd)

    end subroutine calculate_ftchd

end program propagation_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
