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


! ---------------------------------------------------------
integer function X(restart_write_function)(dir, filename, m, f) result(ierr)
  character(len=*), intent(in) :: dir
  character(len=*), intent(in) :: filename
  type(mesh_type),  intent(in) :: m
  R_TYPE,           intent(in) :: f(:)

  call push_sub('restart_write_function')

  ierr = X(output_function) (restart_format, trim(dir), trim(filename), &
     m, f(:), M_ONE)

  call pop_sub()
end function X(restart_write_function)


! ---------------------------------------------------------
integer function X(restart_read_function)(dir, filename, m, f) result(ierr)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(mesh_type),  intent(in)  :: m
  R_TYPE,           intent(out) :: f(:)

  call push_sub('restart_read_function')

  ierr = X(input_function) (trim(dir)//'/'//trim(filename), m, f(:))

  ! If problems, try with the netcdf files
  if(ierr>0) ierr = X(input_function) (trim(dir)//'/'//trim(filename)//'.ncdf', m, f(:))

  call pop_sub()
end function X(restart_read_function)


! ---------------------------------------------------------
integer function X(restart_write) (dir, st, m, iter) result(ierr)
  character(len=*),  intent(in) :: dir
  type(states_type), intent(in) :: st
  type(mesh_type),   intent(in) :: m
  integer, intent(in), optional :: iter

  integer :: iunit, err, ik, ist, idim, i, is
  character(len=6) :: filename

  call push_sub('restart_write')

  call loct_mkdir(dir)

  ierr = 0
  if(mpiv%node==0) then
    call io_assign(iunit)
    open(unit = iunit, file=trim(dir)//'/data', iostat = err, &
       form='formatted', action='write', status='replace', recl = 200)

    if(err .ne. 0) then
      ierr = -1

      call io_free(iunit)
      call pop_sub()
      return
    end if
  end if

  ASSERT(associated(st%X(psi)))

  i = 1
  if(mpiv%node==0) then 
    write(iunit,'(2a)') '#     #kpoint            #st            #dim    ', &
       'filename           occupations           eigenvalue[a.u.]'
    write(iunit,'(a)') '%Wavefunctions'
  end if

  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim
        write(filename,'(i6.6)') i
        
        if(mpiv%node==0) write(unit = iunit, fmt=*) ik, ' | ', ist, ' | ', idim, ' | "', &
           trim(filename), '" | ', st%occ(ist,ik), ' | ', st%eigenval(ist, ik)
        
        if(st%st_start <= ist .and. st%st_end >= ist) then
          if(X(restart_write_function)(dir, filename, m, st%X(psi) (:, idim, ist, ik)) == 0) then
            ierr = ierr + 1
          end if
        end if

        i = i + 1
      end do
    end do
  end do
  
  if(mpiv%node==0)  then 
    write(iunit,'(a)') '%'
    if(present(iter)) write(iunit,'(a,i5)') 'Iter = ', iter
  end if

#if defined(HAVE_MPI)
  call mpi_barrier(MPI_COMM_WORLD, i) ! Since some processors did more than others...
#endif

  if(mpiv%node==0) call io_close(iunit)
  call pop_sub()
end function X(restart_write)


! returns
! <0 => Fatal error
! =0 => read all wave-functions
! >0 => could only read x wavefunctions
integer function X(restart_read) (dir, st, m, iter) result(ierr)
  character(len=*),  intent(in)    :: dir
  type(states_type), intent(inout) :: st
  type(mesh_type),   intent(in)    :: m
  integer, optional, intent(out)   :: iter

  integer :: iunit, err, ik, ist, idim, i, is
  character(len=50) :: blockname
  character(len=12) :: filename
  character(len=1) :: char
  logical, allocatable :: filled(:, :, :)

  call push_sub('restart_read')

  ierr = 0
  call io_assign(iunit)
  open(unit = iunit, file=trim(dir)//'/data', iostat = err, &
       form='formatted', action='read', status='old', recl = 200)

  if(err .ne. 0) then
    ierr = -1

    call io_free(iunit)
    call pop_sub()
    return
  endif

  ASSERT(associated(st%X(psi)))

  allocate(filled(st%d%dim, st%st_start:st%st_end, st%d%nik)); filled = .false.

  read(iunit, *); read(iunit, *) ! Skip two lines...
  do
    read(unit = iunit, fmt = '(a)', iostat = i) char
    if(i.ne.0.or.char=='%') exit
    backspace(unit = iunit)
    read(unit = iunit, iostat = i, fmt = *) ik, char, ist, char, idim, char, filename, char, &
       st%occ(ist,ik), st%eigenval(ist, ik)

    if(index_is_wrong()) cycle
    if(st%st_start <= ist .and. st%st_end >= ist) then
      if(X(restart_read_function) (dir, filename, m, st%X(psi) (:, idim, ist, ik))<=0) then
        filled(idim, ist, ik) = .true.
        ierr = ierr + 1
      end if
    end if
  end do

  if(present(iter)) then
    read(unit = iunit, fmt = *) filename, filename, iter
  endif

!  if(present(vprev)) then
!     do i = 1, nvprev
!        do is = 1, st%d%nspin
!           write(filename,'(a1,i2.2,i3.3)') 'p',i, is
!           call dinput_function(trim(dir)//'/'//trim(filename), m, vprev(1:m%np, is, i), err)
!           ! If problems, try with the netcdf files. However, it will always try to read first
!           ! the plain files, which may be a problem if both are present, and the "good" one is the netcdf.
!           ! This has to be fixed.
!           if(err>0) call dinput_function(trim(dir)//'/'//trim(filename)//'.ncdf', m, vprev(1:m%np, is, i), err)
!           if(err>0) then 
!             write(message(1),'(a)') 'Problem reading potential in restart file for extrapolation'
!             write(message(2),'(a,i3)') 'Error code = ',err
!             call write_warning(2)
!           endif
!        enddo
!     enddo
!  endif

  if(any(.not.filled)) call fill()
  if(ierr == (st%st_end-st%st_start)*st%d%nik*st%d%dim) ierr = 0 ! Alles OK

  deallocate(filled)
  call io_close(iunit)
  call pop_sub()

contains

  subroutine fill() ! Put random function in orbitals that could not be read.
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          if(filled(idim, ist, ik)) cycle
          write(message(1),'(a,3i4)') 'Randomizing wavefunction: #dim, #ist, #ik = ', idim, ist, ik
          call write_warning(1)

          call states_generate_random(st, m, ist, ist)
        end do
      end do
    end do
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

end function X(restart_read)
