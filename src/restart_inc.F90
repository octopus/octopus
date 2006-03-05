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


! ---------------------------------------------------------
subroutine X(restart_write_function)(dir, filename, gr, f, ierr, size)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(grid_t),  intent(in)  :: gr
  integer,          intent(out) :: ierr
  integer,          intent(in)  :: size
  R_TYPE,           intent(in)  :: f(size)

  call push_sub('restart_inc.Xrestart_write_function')

  call X(output_function) (restart_format, trim(dir), trim(filename), &
    gr%m, gr%sb, f(:), M_ONE, ierr, is_tmp=.true.)

  call pop_sub()
end subroutine X(restart_write_function)


! ---------------------------------------------------------
subroutine X(restart_read_function)(dir, filename, m, f, ierr)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(mesh_t),  intent(in)  :: m
  R_TYPE,           intent(out) :: f(1:m%np)
  integer,          intent(out) :: ierr

  call push_sub('restart_inc.Xrestart_read_function')

  ! try first to load plain binary files
  call X(input_function) (trim(dir)//'/'//trim(filename), m, f(:), ierr, is_tmp=.true.)

  ! if we do not succeed try NetCDF
  if(ierr>0) call X(input_function) (trim(dir)//'/'//trim(filename)//'.ncdf', m, f(:), ierr, is_tmp=.true.)

  call pop_sub()
end subroutine X(restart_read_function)


! ---------------------------------------------------------
subroutine X(restart_write) (dir, st, gr, ierr, iter)
  character(len=*),  intent(in)  :: dir
  type(states_t), intent(in)  :: st
  type(grid_t),   intent(in)  :: gr
  integer,           intent(out) :: ierr
  integer, optional, intent(in)  :: iter

  integer :: iunit, iunit2, iunit_mesh, err, ik, ist, idim, i
  character(len=40) :: filename, mformat

  call push_sub('restart_inc.Xrestart_write')

  ASSERT(associated(st%X(psi)))

  mformat = '(f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)'
  ierr = 0
  if(mpi_grp_is_root(mpi_world)) then
    call io_mkdir(dir, is_tmp=.true.)

    iunit = io_open(trim(dir)//'/wfns', action='write', is_tmp=.true.)
    write(iunit,'(a)') '#     #kpoint            #st            #dim    filename'
    write(iunit,'(a)') '%Wavefunctions'

    iunit2 = io_open(trim(dir)//'/occs', action='write', is_tmp=.true.)
    write(iunit2,'(a)') '# occupations           eigenvalue[a.u.]        K-Points'
    write(iunit2,'(a)') '%Occupations_Eigenvalues_K-Points'

    iunit_mesh = io_open(trim(dir)//'/mesh', action='write', is_tmp=.true.)
    write(iunit_mesh,'(a)') '# This file contains the necessary information to generate the'
    write(iunit_mesh,'(a)') '# mesh with which the functions in this directory were calculated,'
    write(iunit_mesh,'(a)') '# except for the geometry of the system.'
    call curvlinear_dump(gr%cv, iunit_mesh)
    call simul_box_dump(gr%sb, iunit_mesh)
    call io_close(iunit_mesh)
  end if

  i = 1
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim
        write(filename,'(i10.10)') i

        if(mpi_grp_is_root(mpi_world)) then
          write(unit=iunit,  fmt=*) ik, ' | ', ist, ' | ', idim, ' | "', trim(filename), '"'
          write(unit=iunit2, fmt=mformat) st%occ(ist,ik), ' | ', st%eigenval(ist, ik), ' | ', &
            st%d%kpoints(1,ik), ' | ', st%d%kpoints(2,ik), ' | ', st%d%kpoints(3,ik)
        end if

        if(st%st_start <= ist .and. st%st_end >= ist) then
          call X(restart_write_function)(dir, filename, gr, st%X(psi) (:, idim, ist, ik), err, size(st%X(psi),1))
          if(err == 0) ierr = ierr + 1
        end if
#if defined(HAVE_MPI)
        call TS(MPI_Barrier)(MPI_COMM_WORLD, err) ! now we all wait
#endif
        i = i + 1
      end do
    end do
  end do

  if(ierr == st%d%nik*(st%st_end - st%st_start + 1)*st%d%dim) ierr = 0 ! Alles OK

  if(mpi_grp_is_root(mpi_world)) then
    write(iunit,'(a)') '%'
    if(present(iter)) write(iunit,'(a,i7)') 'Iter = ', iter
    write(iunit2, '(a)') '%'

    call io_close(iunit)
    call io_close(iunit2)
  end if

#if defined(HAVE_MPI)
  call TS(MPI_Barrier)(MPI_COMM_WORLD, err) ! Since some processors did more than others...
#endif

  call pop_sub()
end subroutine X(restart_write)


! ---------------------------------------------------------
! returns
! <0 => Fatal error
! =0 => read all wave-functions
! >0 => could only read x wavefunctions
subroutine X(restart_read) (dir, st, gr, ierr, iter)
  character(len=*),  intent(in)  :: dir
  type(states_t), intent(inout)  :: st
  type(grid_t),      intent(in)  :: gr
  integer,           intent(out) :: ierr
  integer, optional, intent(out) :: iter

  integer              :: iunit, iunit2, iunit_mesh, err, ik, ist, idim, i
  character(len=12)    :: filename
  character(len=1)     :: char
  logical, allocatable :: filled(:, :, :)
  character(len=256)   :: line
  character(len=50)    :: str

  R_TYPE, allocatable  :: phi(:)
  type(mesh_t)         :: old_mesh
  type(curvlinear_t)   :: old_cv
  type(simul_box_t)    :: old_sb
  logical              :: mesh_change, full_interpolation

  call push_sub('restart_inc.Xrestart_read')

  ! sanity check
  ASSERT(associated(st%X(psi)))

  ierr = 0

  write(str, '(a,i5)') 'Loading restart information'
  call messages_print_stress(stdout, trim(str))

  ! open files to read
  call open_files()
  if(ierr.ne.0) then
    write(message(1),'(a)') 'Could not load any previous restart information.'
    call write_info(1)
    call messages_print_stress(stdout)
    call pop_sub()
    return
  end if

  ! Reads out the previous mesh info.

  call read_previous_mesh()

  ! now we really start
  ALLOCATE(filled(st%d%dim, st%st_start:st%st_end, st%d%nik), st%d%dim*(st%st_end-st%st_start+1)*st%d%nik)
  filled = .false.

  ! Skip two lines.
  call iopar_read(gr%m, iunit,  line, err); call iopar_read(gr%m, iunit,  line, err)
  call iopar_read(gr%m, iunit2, line, err); call iopar_read(gr%m, iunit2, line, err)

  do
    call iopar_read(gr%m, iunit, line, i)
    read(line, '(a)') char
    if(i.ne.0.or.char=='%') exit

    call iopar_backspace(gr%m, iunit)

    call iopar_read(gr%m, iunit, line, err)
    read(line, *) ik, char, ist, char, idim, char, filename
    if(index_is_wrong()) then
      call iopar_read(gr%m, iunit2, line, err)
      cycle
    end if

    call iopar_read(gr%m, iunit2, line, err)
    read(line, *) st%occ(ist, ik), char, st%eigenval(ist, ik)

    if(ist >= st%st_start .and. ist <= st%st_end) then
      if(.not.mesh_change) then
        call X(restart_read_function) (dir, filename, gr%m, st%X(psi) (:, idim, ist, ik), err)
        if(err <= 0) then
          filled(idim, ist, ik) = .true.
          ierr = ierr + 1
        end if
      else
        call X(restart_read_function) (dir, filename, old_mesh, phi, err)
        call X(mf_interpolate)(old_mesh, gr%m, full_interpolation, phi, st%X(psi) (:, idim, ist, ik))
        if(err <= 0) then
          filled(idim, ist, ik) = .true.
          ierr = ierr + 1
        end if
      endif
    end if
  end do

  if(present(iter)) then
    call iopar_read(gr%m, iunit, line, err)
    read(line, *) filename, filename, iter
  end if

  if(any(.not.filled)) call fill()
  if(ierr == 0) then
    ierr = -1 ! no files read
#if !defined(HAVE_MPI)
    write(message(1),'(a)') 'No files could be read. No restart information can be used.'
    call write_info(1)
#endif
  else
    ! Everything o. k.
    if(ierr == (st%st_end - st%st_start + 1)*st%d%nik*st%d%dim) then
      ierr = 0
#if !defined(HAVE_MPI)
      write(message(1),'(a)') 'All the needed files were succesfully read.'
      call write_info(1)
#endif
    else
#if !defined(HAVE_MPI)
      write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' files out of ', &
        (st%st_end - st%st_start + 1)*st%d%nik*st%d%dim, 'could be read.'
      call write_info(1)
#endif
    endif
  end if

  deallocate(filled)
  call iopar_close(gr%m, iunit)
  call iopar_close(gr%m, iunit2)

  if(mesh_change) call interpolation_end()

  call messages_print_stress(stdout)
  call pop_sub()

contains

  subroutine interpolation_init
    call mesh_init_stage_1(old_sb, old_mesh, gr%geo, old_cv, gr%f_der%n_ghost)
    call mesh_init_stage_2(old_sb, old_mesh, gr%geo, old_cv)
    call mesh_init_stage_3(old_mesh, gr%geo, old_cv)
    ALLOCATE(phi(old_mesh%np_global), old_mesh%np_global)
  end subroutine interpolation_init

  subroutine interpolation_end
    call mesh_end(old_mesh)
    deallocate(phi)
  end subroutine interpolation_end

  subroutine read_previous_mesh

    mesh_change = .false.
    full_interpolation = .true.

#if defined(HAVE_MPI)
    return ! For the moment, avoid the complications of parallel stuff.
#endif
    if(present(iter)) return ! No intepolation, in case we are in the td part.

    iunit_mesh  = iopar_open(gr%m, trim(dir)//'/mesh', action='read', status='old', die=.false.)
    if(iunit_mesh < 0) return

    read(iunit_mesh, *); read(iunit_mesh, *); read(iunit_mesh, *)
    call curvlinear_init_from_file(old_cv, iunit_mesh)
    call simul_box_init_from_file(old_sb, gr%geo, iunit_mesh)
    call iopar_close(gr%m, iunit_mesh)

    if( .not. (old_cv.eq.gr%cv) ) then
      mesh_change = .true.
      return
    end if
    if( .not. (old_sb.eq.gr%sb) ) then
      mesh_change = .true.
      ! First, check wether the spacings are the same.
      if(old_sb%h .app. gr%sb%h) full_interpolation = .false.
    end if

    if(mesh_change) then
      if(.not.full_interpolation) then
        write(message(1),'(a)') 'The functions stored in "tmp/restart_gs" were obtained with' 
        write(message(2),'(a)') 'a different simulation box. The possible missing regions will be'
        write(message(3),'(a)') 'padded with zeros.'
        call write_info(3)
      else
        write(message(1),'(a)') 'The functions stored in "tmp/restart_gs" were obtained with' 
        write(message(2),'(a)') 'a different mesh. The values of the functions for the current'
        write(message(3),'(a)') 'calculations will be interpolated/extrapolated.'
        call write_info(3)
      end if
      call interpolation_init()
    endif
    
  end subroutine read_previous_mesh

  ! ---------------------------------------------------------
  subroutine open_files
    iunit  = iopar_open(gr%m, trim(dir)//'/wfns', action='read', status='old', die=.false.)
    if(iunit < 0) then
      ierr = -1
      return
    end if

    iunit2 = iopar_open(gr%m, trim(dir)//'/occs', action='read', status='old', die=.false.)
    if(iunit2 < 0) then
      call iopar_close(gr%m, iunit)
      ierr = -1
    end if
  end subroutine open_files

  ! ---------------------------------------------------------
  subroutine fill() ! Put random function in orbitals that could not be read.
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          if(filled(idim, ist, ik)) cycle
          write(message(1),'(a,3i4)') 'Randomizing wavefunction: #dim, #ist, #ik = ', idim, ist, ik
          call write_warning(1)

          call states_generate_random(st, gr%m, ist, ist)
          st%occ(ist, ik) = M_ZERO
        end do
      end do
    end do
  end subroutine fill

  ! ---------------------------------------------------------
  logical function index_is_wrong() ! .true. if the index (idim, ist, ik) is not present in st structure...
    if(idim > st%d%dim .or. idim < 1 .or.   &
      ist   > st%nst   .or. ist  < 1 .or.   &
      ik    > st%d%nik .or. ik   < 1) then
      index_is_wrong = .true.
    else
      index_is_wrong = .false.
    end if
  end function index_is_wrong

end subroutine X(restart_read)
