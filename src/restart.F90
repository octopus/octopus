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

module restart_m
  use lib_oct_m
  use lib_oct_parser_m
  use global_m
  use messages_m
  use datasets_m
  use io_m
  use states_m
  use curvlinear_m
  use simul_box_m
  use mesh_m
  use mesh_function_m
  use grid_m
  use output_m
  use mpi_m
  use mpi_debug_m
  use varinfo_m
  use math_m

  implicit none

  private
  public ::             &
    restart_init,       &
    clean_stop,         &
    restart_write,      &
    restart_read,       &
    restart_format,     &
    restart_look

  integer, parameter :: &
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
      call write_warning(1)
      clean_stop = .true.
#if defined(HAVE_MPI)
      call MPI_Barrier(mpi_world%comm, mpi_err)
#endif
      if(mpi_grp_is_root(mpi_world)) call loct_rm('stop')
    end if

    return
  end function clean_stop


  ! ---------------------------------------------------------
  ! read restart format information
  subroutine restart_init
    integer :: i

    call push_sub('restart.restart_init')

    !%Variable RestartFileFormat
    !%Type integer
    !%Default restart_plain
    !%Section Generalities::IO
    !%Description
    !% Determines in which format the restart file should be written
    !%Option restart_plain 1
    !% Binary (platform dependent) format
    !%Option restart_netcdf 2
    !% NetCDF (platform independent) format. This requires the NETCDF library.
    !%End
    call loct_parse_int(check_inp('RestartFileFormat'), RESTART_PLAIN, i)
    if(.not.varinfo_valid_option('RestartFileFormat', i)) call input_error('RestartFileFormat')
    call messages_print_var_option(stdout, "RestartFileFormat", i)

    ! Fix the restart format...
    restart_format = output_fill_how("Plain")
#if defined(HAVE_NETCDF)
    if(i == RESTART_NETCDF) then
      restart_format = output_fill_how("NETCDF")
    end if
#endif

    call pop_sub()
  end subroutine restart_init


  ! ---------------------------------------------------------
  subroutine restart_look (dir, m, kpoints, dim, nst, ierr)
    character(len=*), intent(in) :: dir
    type(mesh_t),     intent(in) :: m
    integer,         intent(out) :: kpoints, dim, nst, ierr

    character(len=256) :: line
    character(len=12)  :: filename
    character(len=1)   :: char
    integer :: iunit, iunit2, err, i, ist, idim, ik
    FLOAT :: occ, eigenval

    call push_sub('restart.restart_look')

    ierr = 0
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp = .true., grp = m%mpi_grp)
    if(iunit < 0) then
      ierr = -1
      return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp = .true., grp = m%mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = m%mpi_grp)
      ierr = -1
      return
    end if

    ! Skip two lines.
    call iopar_read(m%mpi_grp, iunit, line, err); call iopar_read(m%mpi_grp, iunit, line, err)
    call iopar_read(m%mpi_grp, iunit2, line, err); call iopar_read(m%mpi_grp, iunit2, line, err)

    kpoints = 1
    dim = 1
    nst = 1
    do
      call iopar_read(m%mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit
      read(line, *) ik, char, ist, char, idim, char, filename
      if(ik > kpoints) kpoints = ik
      if(idim == 2)    dim     = 2
      if(ist>nst)      nst     = ist
      call iopar_read(m%mpi_grp, iunit2, line, err)
      read(line, *) occ, char, eigenval
    end do

    call io_close(iunit, grp = m%mpi_grp)
    call io_close(iunit2, grp = m%mpi_grp)
    call pop_sub()
  end subroutine restart_look

  ! ---------------------------------------------------------
  subroutine restart_write(dir, st, gr, ierr, iter)
    character(len=*),  intent(in)  :: dir
    type(states_t),    intent(in)  :: st
    type(grid_t),      intent(in)  :: gr
    integer,           intent(out) :: ierr
    integer, optional, intent(in)  :: iter

    integer :: iunit, iunit2, iunit_mesh, err, ik, ist, idim, i
    character(len=40) :: filename, mformat

    call push_sub('restart.restart_write')

    ASSERT((associated(st%dpsi) .and. st%d%wfs_type == M_REAL) .or. (associated(st%zpsi) .and. st%d%wfs_type == M_CMPLX))

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
            if (st%d%wfs_type == M_REAL) then
              call drestart_write_function(dir, filename, gr, st%dpsi(:, idim, ist, ik), err, size(st%dpsi,1))
            else
              call zrestart_write_function(dir, filename, gr, st%zpsi(:, idim, ist, ik), err, size(st%zpsi,1))
            end if
            if(err == 0) ierr = ierr + 1
          end if
#if defined(HAVE_MPI)
          call MPI_Barrier(MPI_COMM_WORLD, mpi_err) ! now we all wait
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
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err) ! Since some processors did more than others...
#endif

    call pop_sub()
  end subroutine restart_write


  ! ---------------------------------------------------------
  ! returns
  ! <0 => Fatal error
  ! =0 => read all wave-functions
  ! >0 => could only read x wavefunctions
  subroutine restart_read(dir, st, gr, ierr, iter)
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

    FLOAT, allocatable   :: dphi(:)
    CMPLX, allocatable   :: zphi(:)
    type(mesh_t)         :: old_mesh
    type(curvlinear_t)   :: old_cv
    type(simul_box_t)    :: old_sb
    logical              :: mesh_change, full_interpolation

    call push_sub('restart.restart_read')

    ! sanity check
    ASSERT((associated(st%dpsi) .and. st%d%wfs_type == M_REAL) .or. (associated(st%zpsi) .and. st%d%wfs_type == M_CMPLX))


    ierr = 0

    write(str, '(a,i5)') 'Loading restart information'
    call messages_print_stress(stdout, trim(str))

    ! open files to read
    iunit  = io_open(trim(dir)//'/wfns', action='read', status='old', die=.false., is_tmp = .true., grp = gr%m%mpi_grp)
    if(iunit < 0) then
      ierr = -1
      return
    end if
    iunit2 = io_open(trim(dir)//'/occs', action='read', status='old', die=.false., is_tmp = .true., grp = gr%m%mpi_grp)
    if(iunit2 < 0) then
      call io_close(iunit, grp = gr%m%mpi_grp)
      ierr = -1
    end if

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
    call iopar_read(gr%m%mpi_grp, iunit,  line, err); call iopar_read(gr%m%mpi_grp, iunit,  line, err)
    call iopar_read(gr%m%mpi_grp, iunit2, line, err); call iopar_read(gr%m%mpi_grp, iunit2, line, err)

    do
      call iopar_read(gr%m%mpi_grp, iunit, line, i)
      read(line, '(a)') char
      if(i.ne.0.or.char=='%') exit

      call iopar_backspace(gr%m%mpi_grp, iunit)

      call iopar_read(gr%m%mpi_grp, iunit, line, err)
      read(line, *) ik, char, ist, char, idim, char, filename
      if(index_is_wrong()) then
        call iopar_read(gr%m%mpi_grp, iunit2, line, err)
        cycle
      end if

      call iopar_read(gr%m%mpi_grp, iunit2, line, err)
      read(line, *) st%occ(ist, ik), char, st%eigenval(ist, ik)

      if(ist >= st%st_start .and. ist <= st%st_end) then
        if(.not.mesh_change) then
          if (st%d%wfs_type == M_REAL) then
            call drestart_read_function(dir, filename, gr%m, st%dpsi(:, idim, ist, ik), err)
          else
            call zrestart_read_function(dir, filename, gr%m, st%zpsi(:, idim, ist, ik), err)
          end if
          if(err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          end if
        else
          if (st%d%wfs_type == M_REAL) then
            call drestart_read_function(dir, filename, old_mesh, dphi, err)
            call dmf_interpolate(old_mesh, gr%m, full_interpolation, dphi, st%dpsi(:, idim, ist, ik))
          else
            call zrestart_read_function(dir, filename, old_mesh, zphi, err)
            call zmf_interpolate(old_mesh, gr%m, full_interpolation, zphi, st%zpsi(:, idim, ist, ik))
          end if
          if(err <= 0) then
            filled(idim, ist, ik) = .true.
            ierr = ierr + 1
          end if
        endif
      end if
    end do

    if(present(iter)) then
      call iopar_read(gr%m%mpi_grp, iunit, line, err)
      read(line, *) filename, filename, iter
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      call MPI_Allreduce(ierr, err, 1, MPI_INTEGER, MPI_SUM, st%mpi_grp%comm, mpi_err)
      ierr = err
    end if
#endif

    call fill()
    if(ierr == 0) then
      ierr = -1 ! no files read
      write(message(1),'(a)') 'No files could be read. No restart information can be used.'
      call write_info(1)
    else
      ! Everything o. k.
      if(ierr == st%nst*st%d%nik*st%d%dim) then
        ierr = 0
        write(message(1),'(a)') 'All the needed files were succesfully read.'
        call write_info(1)
      else
        write(message(1),'(a,i4,a,i4,a)') 'Only ', ierr,' files out of ', &
             st%nst*st%d%nik*st%d%dim, ' could be read.'
        call write_info(1)
      endif
    end if

    deallocate(filled)
    call io_close(iunit, grp = gr%m%mpi_grp)
    call io_close(iunit2, grp = gr%m%mpi_grp)

    if(mesh_change) call interpolation_end()

    call messages_print_stress(stdout)
    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine interpolation_init
      call mesh_init_stage_1(old_sb, old_mesh, gr%geo, old_cv, gr%f_der%n_ghost)
      call mesh_init_stage_2(old_sb, old_mesh, gr%geo, old_cv)
      call mesh_init_stage_3(old_mesh, gr%geo, old_cv)
      if (st%d%wfs_type == M_REAL) then
        ALLOCATE(dphi(old_mesh%np_global), old_mesh%np_global)
      else
        ALLOCATE(zphi(old_mesh%np_global), old_mesh%np_global)
      end if
    end subroutine interpolation_init


    ! ---------------------------------------------------------
    subroutine interpolation_end
      call mesh_end(old_mesh)
      if (st%d%wfs_type == M_REAL) then
        deallocate(dphi)
      else
        deallocate(zphi)
      end if
    end subroutine interpolation_end


    ! ---------------------------------------------------------
    subroutine read_previous_mesh

      mesh_change = .false.
      full_interpolation = .true.

#if defined(HAVE_MPI)
      return ! For the moment, avoid the complications of parallel stuff.
#endif
      if(present(iter)) return ! No intepolation, in case we are in the td part.

      iunit_mesh  = io_open(trim(dir)//'/mesh', action='read', status='old', die=.false., grp = gr%m%mpi_grp)
      if(iunit_mesh < 0) return

      read(iunit_mesh, *); read(iunit_mesh, *); read(iunit_mesh, *)
      call curvlinear_init_from_file(old_cv, iunit_mesh)
      call simul_box_init_from_file(old_sb, iunit_mesh)
      call io_close(iunit_mesh, grp = gr%m%mpi_grp)

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

  end subroutine restart_read

#include "undef.F90"
#include "real.F90"
#include "restart_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "restart_inc.F90"

end module restart_m
