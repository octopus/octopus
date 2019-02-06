!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

! ---------------------------------------------------------
!
!> Reads a mesh function from file filename, and puts it into ff. If
!! the map argument is passed, the subroutine will reorder the values
!! in the file according to it, missing values will be filled with
!! zeros. (For the moment this is only implemented for the obf format.)
!!
!! On output, ierr signals how everything went:
!! ierr > 0 => Error. The function ff was not read: \n
!!              1 : illegal filename (must have ".obf" or ".ncdf" extension). \n
!!              2 : file could not be successfully opened. \n
!!              3 : file opened, but error reading. \n
!!              4 : The number of points/mesh dimensions do not coincide. \n
!!              5 : Format error (one or several warnings are written) \n
!! ierr = 0 => Success. \n
!! ierr < 0 => Success, but some kind of type conversion was necessary. The value
!!             of ierr is then: \n
!!             -1 : function in file is real, sp. \n
!!             -2 : function in file is complex, sp. \n
!!             -3 : function in file is real, dp. \n
!!             -4 : function in file is complex, dp. \n
! ---------------------------------------------------------

subroutine X(io_function_input)(filename, mesh, ff, ierr, map)
  character(len=*),  intent(in)    :: filename
  type(mesh_t),      intent(in)    :: mesh
  R_TYPE,            intent(inout) :: ff(:)
  integer,           intent(out)   :: ierr
  integer, optional, intent(in)    :: map(:)
  !

#if defined(HAVE_MPI)
  R_TYPE, allocatable :: ff_global(:)
#endif
  !
  PUSH_SUB(X(io_function_input))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    ! Only root reads. Therefore, only root needs a buffer
    ! ff_global for the whole function.
    SAFE_ALLOCATE(ff_global(1:1))
    if(mpi_grp_is_root(mesh%mpi_grp)) then
      SAFE_DEALLOCATE_A(ff_global)
      SAFE_ALLOCATE(ff_global(1:mesh%np_global))
      call X(io_function_input_global)(filename, mesh, ff_global, ierr, map)
    end if
    if(debug%info) call messages_debug_newlines(2)

    ! Only root knows if the file was successfully read.
    ! Now, it tells everybody else.
    call mpi_debug_in(mesh%vp%comm, C_MPI_BCAST)
    call MPI_Bcast(ierr, 1, MPI_INTEGER, mesh%vp%root, mesh%vp%comm, mpi_err)
    call mpi_debug_out(mesh%vp%comm, C_MPI_BCAST)

    ! Only scatter, when successfully read the file(s).
    if(ierr <= 0) then
      call vec_scatter(mesh%vp, mesh%vp%root, ff_global, ff)
    end if

    SAFE_DEALLOCATE_A(ff_global)
#else
    ! internal error
    ASSERT(.false.) 
#endif
  else
    call X(io_function_input_global)(filename, mesh, ff, ierr, map)
  end if

  POP_SUB(X(io_function_input))

end subroutine X(io_function_input)


! ---------------------------------------------------------
subroutine X(io_function_input_global)(filename, mesh, ff, ierr, map)
  character(len=*),  intent(in)    :: filename
  type(mesh_t),      intent(in)    :: mesh
  R_TYPE,            intent(inout) :: ff(:)
  integer,           intent(out)   :: ierr
  integer, optional, intent(in)    :: map(:)

  integer :: ip, np, ii, jj, kk, file_size
  integer(8) :: dims(3)
  FLOAT, allocatable :: x_in(:, :)
  FLOAT, allocatable :: x_out(:, :)
  type(cube_t) :: cube
  type(cube_function_t) :: cf
  
  R_TYPE, pointer :: read_ff(:)

  call profiling_in(read_prof, "DISK_READ")
  PUSH_SUB(X(io_function_input_global))

  ierr = 0

  select case(trim(io_get_extension(filename)))
  case("obf")

    if(present(map)) then

      call io_binary_get_info(io_workpath(filename), np, file_size, ierr)

      if (ierr == 0) then
        SAFE_ALLOCATE(read_ff(1:np))

        call io_binary_read(io_workpath(filename), np, read_ff, ierr)
        call profiling_count_transfers(np, read_ff(1))
        
        if (ierr == 0) then
          ff(1:mesh%np_global) = M_ZERO
          do ip = 1, min(np, ubound(map, dim = 1))
            if(map(ip) > 0) ff(map(ip)) = read_ff(ip)
          end do
        end if

        SAFE_DEALLOCATE_P(read_ff)
      end if

    else
      call io_binary_read(io_workpath(filename), mesh%np_global, ff, ierr)
      call profiling_count_transfers(mesh%np_global, ff(1))
    end if

  case default
    ierr = 1
  end select

  POP_SUB(X(io_function_input_global))
  call profiling_out(read_prof)
end subroutine X(io_function_input_global)

! ---------------------------------------------------------

subroutine X(io_function_output_vector)(how, dir, fname, mesh, ff, vector_dim, unit, ierr, &
  geo, grp, root, is_global, vector_dim_labels)
  integer(8),                 intent(in)  :: how
  character(len=*),           intent(in)  :: dir
  character(len=*),           intent(in)  :: fname
  type(mesh_t),               intent(in)  :: mesh
  R_TYPE,           target,   intent(in)  :: ff(:, :)
  integer,                    intent(in)  :: vector_dim
  type(unit_t),               intent(in)  :: unit
  integer,                    intent(out) :: ierr
  type(geometry_t), optional, intent(in)  :: geo
  type(mpi_grp_t),  optional, intent(in)  :: grp !< the group that shares the same data, must contain the domains group
  integer,          optional, intent(in)  :: root !< which process is going to write the data
  logical,          optional, intent(in)  :: is_global !< Input data is mesh%np_global? And, thus, it has not be gathered
  character(len=*), optional, intent(in)  :: vector_dim_labels(:)

  integer :: ivd
  integer(8) :: how_seq
  character(len=MAX_PATH_LEN) :: full_fname
  R_TYPE, pointer :: ff_global(:, :)
  logical :: i_am_root, is_global_
  integer :: root_, comm

  PUSH_SUB(X(io_function_output_vector))

  if(present(vector_dim_labels)) then
    ASSERT(ubound(vector_dim_labels, dim = 1) >= vector_dim)
  end if

  ASSERT(vector_dim < 10)

  ierr = 0
  is_global_ = optional_default(is_global, .false.)

  if (is_global_) then
    ASSERT(ubound(ff, dim = 1) == mesh%np_global .or. ubound(ff, dim = 1) == mesh%np_part_global)
  else
    ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)
  end if

  i_am_root = .true.
#ifdef HAVE_MPI
  comm = MPI_COMM_NULL
#endif
  root_ = optional_default(root, 0)

  if(mesh%parallel_in_domains) then
    comm = mesh%vp%comm

    i_am_root = (mesh%vp%rank == root_)

    if (.not. is_global_) then
      if(bitand(how, OPTION__OUTPUTFORMAT__BOUNDARY_POINTS) /= 0) then
        call messages_not_implemented("OutputFormat = boundary_points with domain parallelization")
        SAFE_ALLOCATE(ff_global(1:mesh%np_part_global, 1:vector_dim))
        ! FIXME: needs version of vec_gather that includes boundary points. See ticket #127
      else
        SAFE_ALLOCATE(ff_global(1:mesh%np_global, 1:vector_dim))
      end if

      !note: here we are gathering data that we won`t write if grp is
      !present, but to avoid it we will have to find out all if the
      !processes are members of the domain line where the root of grp
      !lives

      do ivd = 1, vector_dim
#ifdef HAVE_MPI        
        call vec_gather(mesh%vp, root_, ff_global(:, ivd), ff(:, ivd))
#endif
      end do
      
    else
      ff_global => ff
    end if
  else
    ff_global => ff
  end if

  if(present(grp)) then
    i_am_root = i_am_root .and. (grp%rank == root_)
    comm = grp%comm
  end if

  if(i_am_root) then

    how_seq = how
    
    if(bitand(how, OPTION__OUTPUTFORMAT__VTK) /= 0) call out_vtk()
    how_seq = bitand(how_seq, not(OPTION__OUTPUTFORMAT__VTK)) ! remove from the list of formats

    if(how_seq /= 0) then !perhaps there is nothing left to do
      do ivd = 1, vector_dim
        if(present(vector_dim_labels)) then
          full_fname = trim(fname)//'-'//vector_dim_labels(ivd)
        else
          write(full_fname, '(2a,i1)') trim(fname), '-', ivd
        end if
        
        call X(io_function_output_global)(how_seq, dir, full_fname, mesh, ff_global(:, ivd), unit, ierr, geo)
      end do
    end if

  end if

#ifdef HAVE_MPI
  if(comm /= MPI_COMM_NULL .and. comm /= 0 .and. .not. is_global_) then
    ! I have to broadcast the error code
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, comm, mpi_err)
    ! Add a barrier to ensure that the process are synchronized
    call MPI_Barrier(comm, mpi_err)
  end if
#endif
  
  if(mesh%parallel_in_domains .and. .not. is_global_) then
    SAFE_DEALLOCATE_P(ff_global)
  else
    nullify(ff_global)
  end if

  POP_SUB(X(io_function_output_vector))

contains

  subroutine out_vtk()
    type(cube_t) :: cube
    type(cube_function_t), allocatable :: cf(:)
    character(len=MAX_PATH_LEN) :: filename
    FLOAT :: dk(3)  
    integer :: ii

    PUSH_SUB(X(io_function_output_vector).out_vtk)

    call cube_init(cube, mesh%idx%ll, mesh%sb)

    SAFE_ALLOCATE(cf(1:vector_dim))

    do ivd = 1, vector_dim
      call cube_function_null(cf(ivd))
      call X(cube_function_alloc_RS)(cube, cf(ivd))
      call X(mesh_to_cube)(mesh, ff_global(:, ivd), cube, cf(ivd))
    end do

    filename = io_workpath(trim(dir)//'/'//trim(fname)//".vtk")

    forall (ii = 1:3) dk(ii)= units_from_atomic(units_out%length, mesh%spacing(ii))

    call X(vtk_out_cf_vector)(filename, fname, ierr, cf, vector_dim, cube, dk, unit)

    do ivd = 1, vector_dim
      call X(cube_function_free_RS)(cube, cf(ivd))
    end do

    call cube_end(cube)

    SAFE_DEALLOCATE_A(cf)

    POP_SUB(X(io_function_output_vector).out_vtk)    
  end subroutine out_vtk

end subroutine X(io_function_output_vector)

! ---------------------------------------------------------
subroutine X(io_function_output) (how, dir, fname, mesh, ff, unit, ierr, geo, grp, root, is_global)
  integer(8),                 intent(in)  :: how
  character(len=*),           intent(in)  :: dir
  character(len=*),           intent(in)  :: fname
  type(mesh_t),               intent(in)  :: mesh
  R_TYPE,           target,   intent(in)  :: ff(:)
  type(unit_t),               intent(in)  :: unit
  integer,                    intent(out) :: ierr
  type(geometry_t), optional, intent(in)  :: geo
  type(mpi_grp_t),  optional, intent(in)  :: grp !< the group that shares the same data, must contain the domains group
  integer,          optional, intent(in)  :: root !< which process is going to write the data
  logical,          optional, intent(in)  :: is_global !< Input data is mesh%np_global? And, thus, it has not be gathered

  logical :: is_global_
#if defined(HAVE_MPI)
  logical :: i_am_root
  integer :: root_, comm
  R_TYPE, pointer :: ff_global(:)
#endif

  PUSH_SUB(X(io_function_output))
  ierr = 0
  is_global_ = optional_default(is_global, .false.)
  if (is_global_) then
    ASSERT(ubound(ff, dim = 1) == mesh%np_global .or. ubound(ff, dim = 1) == mesh%np_part_global)
  else
    ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)
  end if

#if defined(HAVE_MPI)

  i_am_root = .true.
  comm = MPI_COMM_NULL
  root_ = optional_default(root, 0)
  if(mesh%parallel_in_domains) then
    comm = mesh%vp%comm
    i_am_root = (mesh%vp%rank == root_)
    if (.not. is_global_) then
      if(bitand(how, OPTION__OUTPUTFORMAT__BOUNDARY_POINTS) /= 0) then
        call messages_not_implemented("OutputFormat = boundary_points with domain parallelization")
        SAFE_ALLOCATE(ff_global(1:mesh%np_part_global))
        ! FIXME: needs version of vec_gather that includes boundary points. See ticket #127
      else
        SAFE_ALLOCATE(ff_global(1:mesh%np_global))
      end if

      !note: here we are gathering data that we won`t write if grp is
      !present, but to avoid it we will have to find out all if the
      !processes are members of the domain line where the root of grp
      !lives
      call vec_gather(mesh%vp, root_, ff_global, ff)
    else
      ff_global => ff
    end if
  else
    ff_global => ff
  end if

  if(present(grp)) then
    i_am_root = i_am_root .and. (grp%rank == root_)
    comm = grp%comm
  end if

  if(i_am_root) then
    call X(io_function_output_global)(how, dir, fname, mesh, ff_global, unit, ierr, geo = geo)
  end if
  if(comm /= MPI_COMM_NULL .and. comm /= 0 .and. .not. is_global_) then
    ! I have to broadcast the error code
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, comm, mpi_err)
    ! Add a barrier to ensure that the process are synchronized
    call MPI_Barrier(comm, mpi_err)
  end if

  if(mesh%parallel_in_domains .and. .not. is_global_) then
    SAFE_DEALLOCATE_P(ff_global)
  else
    nullify(ff_global)
  end if

#else

  ! serial mode
  ASSERT(.not. mesh%parallel_in_domains)
  call X(io_function_output_global)(how, dir, fname, mesh, ff, unit, ierr, geo = geo)

#endif

  POP_SUB(X(io_function_output))
end subroutine X(io_function_output)


! ---------------------------------------------------------
subroutine X(io_function_output_global) (how, dir, fname, mesh, ff, unit, ierr, geo)
  integer(8),                 intent(in)  :: how
  character(len=*),           intent(in)  :: dir, fname
  type(mesh_t),               intent(in)  :: mesh
  R_TYPE,                     intent(in)  :: ff(:)  !< (mesh%np_global or mesh%np_part_global)
  type(unit_t),               intent(in)  :: unit
  integer,                    intent(out) :: ierr
  type(geometry_t), optional, intent(in)  :: geo

  character(len=512) :: filename
  character(len=20)  :: mformat, mformat2, mfmtheader
  integer            :: iunit, ip, idir, jj, np_max
  FLOAT              :: x0

  call profiling_in(write_prof, "DISK_WRITE")
  PUSH_SUB(X(io_function_output_global))

  call io_mkdir(dir)

! Define the format
  mformat    = '(99es23.14E3)'
  mformat2   = '(i12,99es34.24E3)'
  mfmtheader = '(a,a10,5a23)'

  ASSERT(how > 0)
  ASSERT(ubound(ff, dim = 1) >= mesh%np_global)

  np_max = mesh%np_global
  ! should we output boundary points?
  if(bitand(how, OPTION__OUTPUTFORMAT__BOUNDARY_POINTS) /= 0) then
    if(ubound(ff, dim = 1) >= mesh%np_part_global) then
      np_max = mesh%np_part_global
    else
      write(message(1),'(2a)') trim(fname), ': not outputting boundary points; they are not available'
      call messages_warning(1)
      ! FIXME: in this case, one could allocate an array of the larger size and apply boundary conditions
    endif
  endif

  if(bitand(how, OPTION__OUTPUTFORMAT__BINARY)     /= 0) call out_binary()
  if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_X)     /= 0) call out_axis (1, 2, 3) ! x ; y=0,z=0
  if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Y)     /= 0) call out_axis (2, 1, 3) ! y ; x=0,z=0
  if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Z)     /= 0) call out_axis (3, 1, 2) ! z ; x=0,y=0
  if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_X)    /= 0) call out_plane(1, 2, 3) ! x=0; y; z;
  if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Y)    /= 0) call out_plane(2, 1, 3) ! y=0; x; z;
  if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Z)    /= 0) call out_plane(3, 1, 2) ! z=0; x; y;
  if(bitand(how, OPTION__OUTPUTFORMAT__MESH_INDEX) /= 0) call out_mesh_index()
  if(bitand(how, OPTION__OUTPUTFORMAT__CUBE)       /= 0) call out_cube()
  if(bitand(how, OPTION__OUTPUTFORMAT__VTK) /= 0) call out_vtk()

  POP_SUB(X(io_function_output_global))
  call profiling_out(write_prof)

contains

  ! ---------------------------------------------------------
  subroutine out_binary()
    character(len=512) :: workdir

    PUSH_SUB(X(io_function_output_global).out_binary)

    workdir = io_workpath(dir)
    call io_binary_write(trim(workdir)//'/'//trim(fname)//'.obf', np_max, ff, ierr)

    call profiling_count_transfers(np_max, ff(1))
    POP_SUB(X(io_function_output_global).out_binary)
  end subroutine out_binary


  ! ---------------------------------------------------------
  subroutine out_axis(d1, d2, d3)
    integer, intent(in) :: d1, d2, d3
    
    integer :: ixvect(MAX_DIM)
    FLOAT   :: xx(1:MAX_DIM)
    R_TYPE  :: fu

    PUSH_SUB(X(io_function_output_global).out_axis)

    filename = trim(dir)//'/'//trim(fname)//"."//index2axis(d2)//"=0,"//index2axis(d3)//"=0"
    iunit = io_open(filename, action='write')

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d1), 'Re', 'Im'
    do ip = 1, np_max
      call index_to_coords(mesh%idx, ip, ixvect)

      if(ixvect(d2)==0.and.ixvect(d3)==0) then
        xx = units_from_atomic(units_out%length, mesh_x_global(mesh, ip))
        fu = units_from_atomic(unit, ff(ip))
        write(iunit, mformat, iostat=ierr) xx(d1), fu
      end if
    end do

    call io_close(iunit)
    POP_SUB(X(io_function_output_global).out_axis)
  end subroutine out_axis


  ! ---------------------------------------------------------
  subroutine out_plane(d1, d2, d3)
    integer, intent(in) :: d1, d2, d3

    integer :: ix, iy, iz, jdim
    integer :: ixvect(MAX_DIM)
    integer :: ixvect_test(MAX_DIM)
    FLOAT   :: xx(1:MAX_DIM)
    R_TYPE  :: fu

    PUSH_SUB(X(io_function_output_global).out_plane)

    filename = trim(dir)//'/'//trim(fname)//"."//index2axis(d1)//"=0"
    iunit = io_open(filename, action='write')

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d2), index2axis(d3), 'Re', 'Im'

! here we find the indices for coordinate 0 along all directions apart from d2
! and d3, to get a plane. Do the same as for ix, but with all the other
! dimensions.

    ixvect=1
    do jdim=1,mesh%sb%dim
      if (jdim==d2 .or. jdim==d3) cycle

      do ix = mesh%idx%nr(1, jdim), mesh%idx%nr(2, jdim)
! NOTE: MJV: how could this return anything but ix=0? Answ: if there is a shift in origin
        ixvect_test = 1
        ixvect_test(jdim) = ix
        ip = index_from_coords(mesh%idx, ixvect_test)
        if(ip /= 0) then 
          call index_to_coords(mesh%idx, ip, ixvect_test)
          if(ixvect_test(jdim) == 0) exit
        end if
      end do
      ixvect(jdim) = ix
    end do ! loop over dimensions

    ! have found ix such that coordinate d1 is 0 for this value of ix
    ! ixvect is prepared for all dimensions apart from d2 and d3
    
    do iy = mesh%idx%nr(1, d2), mesh%idx%nr(2, d2)
      write(iunit, mformat, iostat=ierr)
      do iz = mesh%idx%nr(1, d3), mesh%idx%nr(2, d3)

        ixvect(d2) = iy
        ixvect(d3) = iz
        ip = index_from_coords(mesh%idx, ixvect)

        if(ip <= np_max .and. ip > 0) then
          xx = units_from_atomic(units_out%length, mesh_x_global(mesh, ip))
          fu = units_from_atomic(unit, ff(ip))
          write(iunit, mformat, iostat=ierr)  &
            xx(d2), xx(d3), fu
        end if
      end do
    end do

    call io_close(iunit)

    POP_SUB(X(io_function_output_global).out_plane)
  end subroutine out_plane

  ! ---------------------------------------------------------
  subroutine out_mesh_index()
    FLOAT :: xx(1:MAX_DIM)
    R_TYPE :: fu

    integer :: idir

    PUSH_SUB(X(io_function_output_global).out_mesh_index)

    iunit = io_open(trim(dir)//'/'//trim(fname)//".mesh_index", action='write')

    write(iunit, mfmtheader, iostat=ierr) '#', 'Index', 'x', 'y', 'z', 'Re', 'Im'
    xx = mesh_x_global(mesh, 1)
    x0 = xx(1)
    if(ierr == 0) write(iunit, mformat, iostat=ierr)

    do ip = 1, np_max
      xx = mesh_x_global(mesh, ip)
       if (ierr == 0 .and. x0 /= xx(1)) then
          write(iunit, mformat, iostat=ierr)      ! write extra lines for gnuplot grid mode
          x0 = xx(1)
       end if
       fu = units_from_atomic(unit, ff(ip))
       if(ierr==0) write(iunit, mformat2, iostat=ierr) ip, &
         (units_from_atomic(units_out%length, xx(idir)), idir = 1, 3), fu
    end do

    if(ierr == 0) write(iunit, mformat, iostat=ierr)
    call io_close(iunit)

    POP_SUB(X(io_function_output_global).out_mesh_index)
  end subroutine out_mesh_index

  ! ---------------------------------------------------------
  !> see http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/
  !! Writes only real part
  subroutine out_cube()

    integer :: ix, iy, iz, idir, idir2, iatom
    integer :: int_unit(3)
    FLOAT   :: offset(MAX_DIM)
    type(cube_t) :: cube
    type(cube_function_t) :: cf
    character(len=8) :: fmt

    PUSH_SUB(X(io_function_output_global).out_cube)

    ASSERT(present(geo))

    ! put values in a nice cube
    call cube_init(cube, mesh%idx%ll, mesh%sb)
    call cube_function_null(cf)
    call X(cube_function_alloc_RS) (cube, cf)
    call X(mesh_to_cube) (mesh, ff, cube, cf)

    ! the offset is different in periodic directions
    offset = M_ZERO
    offset(1:3) = units_from_atomic(units_out%length, -matmul(mesh%sb%rlattice_primitive(1:3,1:3), mesh%sb%lsize(1:3)))

    do idir = mesh%sb%periodic_dim+1, 3
      offset(idir) = units_from_atomic(units_out%length, -(cube%rs_n_global(idir) - 1)/2*mesh%spacing(idir))
    end do

    iunit = io_open(trim(dir)//'/'//trim(fname)//".cube", action='write')

    write(iunit, '(2a)') 'Generated by octopus ', trim(conf%version)
    write(iunit, '(4a)') 'git: ', trim(conf%git_commit), " build: ",  trim(conf%build_time)
    write(iunit, '(i5,3f12.6)') geo%natoms, offset(1:3)

    ! According to http://gaussian.com/cubegen/
    ! If N1<0 the input cube coordinates are assumed to be in Bohr, otherwise, they are interpreted as Angstroms. 
    int_unit(1:3) = 1
    if (units_out%length%abbrev == "b") then
      int_unit(1) = -1
    end if

    do idir = 1, 3
      write(iunit, '(i5,3f12.6)') int_unit(idir)*cube%rs_n_global(idir), (units_from_atomic(units_out%length, &
        mesh%spacing(idir)*mesh%sb%rlattice_primitive(idir2, idir)), idir2 = 1, 3)
    end do
    do iatom = 1, geo%natoms
      write(iunit, '(i5,4f12.6)') int(species_z(geo%atom(iatom)%species)),  M_ZERO, &
        (units_from_atomic(units_out%length, geo%atom(iatom)%x(idir)), idir = 1, 3)
    end do

    do ix = 1, cube%rs_n_global(1)
      do iy = 1, cube%rs_n_global(2)
        do iz = 1, cube%rs_n_global(3), 6

          if(iz + 6 - 1 <= cube%rs_n_global(3)) then
            write(iunit,'(6e14.6)') units_from_atomic(unit, R_REAL(cf%X(RS)(ix, iy, iz:iz + 6 - 1)))
          else
            write(fmt, '(a,i1,a)') '(', cube%rs_n_global(3) - iz + 1, 'e14.6)'
            write(iunit, trim(fmt)) units_from_atomic(unit, R_REAL(cf%X(RS)(ix, iy, iz:cube%rs_n_global(3))))
          end if
        
        end do
      end do
    end do

    call io_close(iunit)

    call X(cube_function_free_RS)(cube, cf)
    call cube_end(cube)

    POP_SUB(X(io_function_output_global).out_cube)
  end subroutine out_cube

  ! ---------------------------------------------------------

  subroutine out_vtk()
    type(cube_t) :: cube
    type(cube_function_t) :: cf
    FLOAT :: dk(3), pnt(3)
    integer :: i, i1, i2, i3 
    FLOAT, ALLOCATABLE :: points(:,:,:,:)
    
    PUSH_SUB(X(io_function_output_global).out_vtk)

    forall (i = 1:3) dk(i)= units_from_atomic(units_out%length, mesh%spacing(i))
    
    call cube_init(cube, mesh%idx%ll, mesh%sb, spacing = dk )
    call cube_function_null(cf)
    call X(cube_function_alloc_RS) (cube, cf)
    call X(mesh_to_cube) (mesh, ff, cube, cf)

    filename = io_workpath(trim(dir)//'/'//trim(fname)//".vtk")
   

    if(mesh%sb%nonorthogonal) then
      ! non-orthogonal grid
      SAFE_ALLOCATE(points(cube%rs_n_global(1),cube%rs_n_global(2),cube%rs_n_global(3),3))
      
      do i1 =1 , cube%rs_n_global(1)
        do i2 =1 , cube%rs_n_global(2)
          do i3 =1 , cube%rs_n_global(3)
            pnt(1:3) =(/cube%Lrs(i1, 1),cube%Lrs(i2, 2),cube%Lrs(i3, 3)/)
            points(i1,i2,i3, 1:3) = matmul(mesh%sb%rlattice_primitive(1:3,1:3), pnt(1:3))
          end do
        end do
      end do
      
      call X(vtk_out_cf_structured)(filename, fname, ierr, cf, cube, unit, points)     
      SAFE_DEALLOCATE_A(points) 
    else  
      !Ordinary grid
      call X(vtk_out_cf)(filename, fname, ierr, cf, cube, dk(:), unit)
    end if  

    call X(cube_function_free_RS)(cube, cf)
    call cube_end(cube)

    POP_SUB(X(io_function_output_global).out_vtk)
  end subroutine out_vtk

end subroutine X(io_function_output_global)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
