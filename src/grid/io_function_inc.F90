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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

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
!!              5 : Format or NetCDF error (one or several warnings are written) \n
!! ierr = 0 => Success. \n
!! ierr < 0 => Success, but some kind of type conversion was necessary. The value
!!             of ierr is then: \n
!!             -1 : function in file is real, sp. \n
!!             -2 : function in file is complex, sp. \n
!!             -3 : function in file is real, dp. \n
!!             -4 : function in file is complex, dp. \n
! ---------------------------------------------------------

subroutine X(io_function_input)(filename, mesh, ff, ierr, is_tmp, map)
  character(len=*),  intent(in)    :: filename
  type(mesh_t),      intent(in)    :: mesh
  R_TYPE,            intent(inout) :: ff(:)
  integer,           intent(out)   :: ierr
  logical, optional, intent(in)    :: is_tmp
  integer, optional, intent(in)    :: map(:)
  !
  logical :: is_tmp_

#if defined(HAVE_MPI)
  R_TYPE, allocatable :: ff_global(:)
#endif
  !
  PUSH_SUB(X(io_function_input))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  is_tmp_ = .false.
  if(present(is_tmp)) is_tmp_ = is_tmp

  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    ! Only root reads. Therefore, only root needs a buffer
    ! ff_global for the whole function.
    SAFE_ALLOCATE(ff_global(1:1))
    if(mpi_grp_is_root(mesh%mpi_grp)) then
      SAFE_DEALLOCATE_A(ff_global)
      SAFE_ALLOCATE(ff_global(1:mesh%np_global))
      call X(io_function_input_global)(filename, mesh, ff_global, ierr, is_tmp_, map)
    end if
    if(in_debug_mode) call messages_debug_newlines(2)

    ! Only root knows if the file was successfully read.
    ! Now, it tells everybody else.
    call mpi_debug_in(mesh%vp%comm, C_MPI_BCAST)
    call MPI_Bcast(ierr, 1, MPI_INTEGER, mesh%vp%root, mesh%vp%comm, mpi_err)
    call mpi_debug_out(mesh%vp%comm, C_MPI_BCAST)

    ! Only scatter, when successfully read the file(s).
    if(ierr.le.0) then
      call X(vec_scatter)(mesh%vp, mesh%vp%root, ff_global, ff)
    end if

    SAFE_DEALLOCATE_A(ff_global)
#else
    ! internal error
    ASSERT(.false.) 
#endif
  else
    call X(io_function_input_global)(filename, mesh, ff, ierr, is_tmp_, map)
  end if

  POP_SUB(X(io_function_input))

end subroutine X(io_function_input)


! ---------------------------------------------------------
subroutine X(io_function_input_global)(filename, mesh, ff, ierr, is_tmp, map)
  character(len=*),  intent(in)    :: filename
  type(mesh_t),      intent(in)    :: mesh
  R_TYPE,            intent(inout) :: ff(:)
  integer,           intent(out)   :: ierr
  logical,           intent(in)    :: is_tmp
  integer, optional, intent(in)    :: map(:)

  integer :: ip, np, ii, jj, kk
  integer(8) :: dims(3)
  FLOAT, allocatable :: x_in(:, :)
  FLOAT, allocatable :: x_out(:, :)
  type(cube_t) :: cube
  type(cube_function_t) :: cf
  
#if defined(HAVE_NETCDF)
  character(len=512) :: file
  integer :: ncid, status
  integer :: function_kind
#if defined(R_TCOMPLEX)
  type(cube_function_t) :: re, im
#endif
#endif
  R_TYPE, pointer :: read_ff(:)

  call profiling_in(read_prof, "DISK_READ")
  PUSH_SUB(X(io_function_input_global))

  ierr = 0

#if defined(HAVE_NETCDF)
  function_kind = X(output_kind)*kind(ff(1)) ! +4 for real, single; +8 for real, double;
  ! -4 for complex, single, -8 for real, double
#endif

  select case(trim(io_get_extension(filename)))
#if defined(HAVE_NETCDF)
  case("ncdf")
    ASSERT(.not. present(map))
    file = io_workpath(filename, is_tmp=is_tmp)
    status = nf90_open(trim(file), NF90_WRITE, ncid)
    if(status.ne.NF90_NOERR) then
      ierr = 2
    else
      call cube_init(cube, mesh%idx%ll, mesh%sb)
      call cube_function_null(cf)
      call X(cube_function_alloc_RS)(cube, cf)
#if defined(R_TCOMPLEX)
      call cube_function_null(re)
      call cube_function_null(im)
      call dcube_function_alloc_RS(cube, re)
      call dcube_function_alloc_RS(cube, im)
      call read_netcdf()
      cf%zRS = re%dRS + M_zI*im%dRS
      call X(cube_to_mesh) (cube, cf, mesh, ff)
      call dcube_function_free_RS(cube, re)
      call dcube_function_free_RS(cube, im)
#else
      call read_netcdf()
      call X(cube_to_mesh)(cube, cf, mesh, ff)
#endif
      call X(cube_function_free_RS)(cube, cf)
      call cube_end(cube)
    end if
#endif
  case("obf")

    if(present(map)) then

      call io_binary_get_info(filename, np, ierr)

      SAFE_ALLOCATE(read_ff(1:np))

      call io_binary_read(filename, np, read_ff, ierr)
      call profiling_count_transfers(np, read_ff(1))

      ff(1:mesh%np_global) = M_ZERO
      do ip = 1, min(np, ubound(map, dim = 1))
        if(map(ip) > 0) ff(map(ip)) = read_ff(ip)
      end do

      SAFE_DEALLOCATE_P(read_ff)

    else
      call io_binary_read(filename, mesh%np_global, ff, ierr)
      call profiling_count_transfers(mesh%np_global, ff(1))
    end if

  case("csv")
    if (mesh%sb%box_shape .ne. PARALLELEPIPED) then
      message(1) = "Box shape must be parallelepiped when a .csv file is used."
      call messages_fatal(1)
    end if 

    call cube_init(cube, mesh%idx%ll, mesh%sb)
    call cube_function_null(cf)
    call X(cube_function_alloc_RS)(cube, cf)
    
    call io_csv_get_info(filename, dims, ierr)
    
    if (ierr .ne. 0) then
      message(1) = "Could not read file "//trim(filename)//""
      call messages_fatal(1)
    end if
    
    SAFE_ALLOCATE(read_ff(1:dims(1)*dims(2)*dims(3)))
    call io_csv_read(filename, dims(1)*dims(2)*dims(3), read_ff, ierr)

    if (ierr .ne. 0) then
      message(1) = "Could not read file "//trim(filename)//""
      call messages_fatal(1)
    end if

    if(mesh%sb%dim == 1) then
      SAFE_ALLOCATE(x_in(1:dims(1), 1:1))
      SAFE_ALLOCATE(x_out(1:cube%rs_n_global(1), 1:1))
      
      do ii = 1, dims(1)
        x_in(ii,:) = (/ real(ii-1)*(real(cube%rs_n_global(1)-1)/(dims(1)-1)) /)
      end do

      
      do ii = 1, cube%rs_n_global(1)
        x_out(ii,1) = real(ii-1)
      end do
        
      call X(mf_interpolate_points)(2, int(dims(1),4), x_in(:,:),&
           &read_ff, cube%rs_n_global(1), x_out(:,:), ff)
           
      SAFE_DEALLOCATE_A(x_in)
      SAFE_DEALLOCATE_A(x_out)
      
    else if(mesh%sb%dim == 2) then
      SAFE_ALLOCATE(x_in(1:(dims(1)*dims(2)), 1:2))
      SAFE_ALLOCATE(x_out(1:(cube%rs_n_global(1)*cube%rs_n_global(2)), 1:2))
      
      do ii = 1, dims(2)
        do jj = 1, dims(1)
          x_in((ii-1)*dims(1) + jj,:) = & 
              &(/ real(jj-1)*(real(cube%rs_n_global(1)-1)/(dims(1)-1)),&
              &   real(ii-1)*(real(cube%rs_n_global(2)-1)/(dims(2)-1)) /)
        end do
      end do
      
      do ii = 1, cube%rs_n_global(2)
        do jj = 1, cube%rs_n_global(1)
          x_out((ii-1)*cube%rs_n_global(1) + jj,:) = (/ real(jj-1), real(ii-1) /)
        end do
      end do
      call X(mf_interpolate_points)(2, int(dims(1)*dims(2),4), x_in(:,:),&
           &read_ff, cube%rs_n_global(1)*cube%rs_n_global(2), x_out(:,:), ff)
           
      SAFE_DEALLOCATE_A(x_in)
      SAFE_DEALLOCATE_A(x_out)
    else if(mesh%sb%dim == 3) then
      SAFE_ALLOCATE(x_in(1:(dims(1)*dims(2)*dims(3)), 1:3))
      SAFE_ALLOCATE(x_out(1:(cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)), 1:3))
      
      do ii = 1, dims(3)
        do jj = 1, dims(2)
          do kk = 1, dims(1)
            x_in((ii-1)*dims(1)*dims(2) + (jj-1)*dims(1) + kk,:) = &
              &(/ real(kk-1)*(real(cube%rs_n_global(1)-1)/(dims(1)-1))  ,&
              &   real(jj-1)*(real(cube%rs_n_global(2)-1)/(dims(2)-1))  ,&
              &   real(ii-1)*(real(cube%rs_n_global(3)-1)/(dims(3)-1))  /)
          end do
        end do
      end do
    
      do ii = 1, cube%rs_n_global(3)
        do jj = 1, cube%rs_n_global(2)
          do kk = 1, cube%rs_n_global(1)
            x_out((ii-1)*cube%rs_n_global(1)*cube%rs_n_global(2) + (jj-1)*cube%rs_n_global(1) + kk,:) = &
              &(/ real(kk-1), real(jj-1), real(ii-1) /)
            end do
          end do
        end do
      call X(mf_interpolate_points)(3, int(dims(1)*dims(2)*dims(3),4), x_in(:,:),&
           &read_ff, cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3), x_out(:,:), ff)
      SAFE_DEALLOCATE_A(x_in)
      SAFE_DEALLOCATE_A(x_out)
    end if
    
    cf%X(RS) = reshape(ff, (/ cube%rs_n_global(1), cube%rs_n_global(2), cube%rs_n_global(3) /))
    
    call X(cube_to_mesh) (cube, cf, mesh, ff)
    call X(cube_function_free_RS)(cube, cf)
    call cube_end(cube)
    
    SAFE_DEALLOCATE_P(read_ff)
  case default
    ierr = 1
  end select

  POP_SUB(X(io_function_input_global))
  call profiling_out(read_prof)

#if defined(HAVE_NETCDF)

contains

  ! ---------------------------------------------------------
  subroutine read_netcdf()
    integer :: data_id, data_im_id, &
      dim_data_id(MAX_DIM), ndim(MAX_DIM), xtype, file_kind
    FLOAT, allocatable :: xx(:, :, :)

    PUSH_SUB(X(io_function_input_global).read_netcdf)

    !Inquire about dimensions
    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_1", dim_data_id(1))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_2", dim_data_id(2))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_3", dim_data_id(3))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(1), len = ndim(3))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(2), len = ndim(2))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(3), len = ndim(1))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if((ndim(1) .ne. cube%rs_n_global(1)) .or. &
      (ndim(2) .ne. cube%rs_n_global(2)) .or. &
      (ndim(3) .ne. cube%rs_n_global(3))) then
      ierr = 12
      POP_SUB(X(io_function_input_global).read_netcdf)
      return
    end if

    if(status == NF90_NOERR) then
      status = nf90_inq_varid (ncid, "rdata", data_id)
      call ncdf_error('nf90_inq_varid', status, file, ierr)
    end if
    status = nf90_inq_varid(ncid, "idata", data_im_id)
    if(status == 0) then
      file_kind = -1
    else
      file_kind = 1
    end if
    status = 0

    if(status == NF90_NOERR) then
      status = nf90_inquire_variable (ncid, data_id, xtype = xtype)
      call ncdf_error('nf90_inquire_variable', status, file, ierr)
    end if

    if(xtype == NF90_FLOAT) then
      file_kind = file_kind*4
    else
      file_kind = file_kind*8
    end if
    if(file_kind .ne. function_kind) then
      select case(file_kind)
      case(4)
        ierr = -1
      case(-4)
        ierr = -2
      case(8)
        ierr = -3
      case(-8)
        ierr = -4
      end select
    end if

    SAFE_ALLOCATE(xx(1:cube%rs_n_global(3), 1:cube%rs_n_global(2), 1:cube%rs_n_global(1)))
#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      select case(mesh%sb%dim)
      case(1)
        status = nf90_get_var (ncid, data_id, xx(1, 1, :))
      case(2)
        status = nf90_get_var (ncid, data_id, xx(1, :, :))
      case(3)
        status = nf90_get_var (ncid, data_id, xx)
      end select
      call transpose3(xx, re%dRS)
      call ncdf_error('nf90_get_var', status, file, ierr)
    end if
    if(file_kind<0) then
      if(status == NF90_NOERR) then
        select case(mesh%sb%dim)
        case(1)
          status = nf90_get_var (ncid, data_im_id, xx(1, 1, :))
        case(2)
          status = nf90_get_var (ncid, data_im_id, xx(1, :, :))
        case(3)
          status = nf90_get_var (ncid, data_im_id, xx)
        end select
        call transpose3(xx, im%dRS)
        call ncdf_error('nf90_get_var', status, file, ierr)
      end if
    else
      im%dRS = M_ZERO
    end if
#else
    if(status == NF90_NOERR) then
      select case(mesh%sb%dim)
      case(1)
        status = nf90_get_var (ncid, data_id, xx(1, 1, :))
      case(2)
        status = nf90_get_var (ncid, data_id, xx(1, :, :))
      case(3)
        status = nf90_get_var (ncid, data_id, xx)
      end select
      call transpose3(xx, cf%dRS)
      call ncdf_error('nf90_get_var', status, file, ierr)
    end if
#endif
    SAFE_DEALLOCATE_A(xx)

    status = nf90_close(ncid)
    POP_SUB(X(io_function_input_global).read_netcdf)
  end subroutine read_netcdf

#endif

end subroutine X(io_function_input_global)


! ---------------------------------------------------------
subroutine X(io_function_output) (how, dir, fname, mesh, ff, unit, ierr, is_tmp, geo, grp)
  integer,                    intent(in)  :: how
  character(len=*),           intent(in)  :: dir
  character(len=*),           intent(in)  :: fname
  type(mesh_t),               intent(in)  :: mesh
  R_TYPE,           target,   intent(in)  :: ff(:)
  type(unit_t),               intent(in)  :: unit
  integer,                    intent(out) :: ierr
  logical,          optional, intent(in)  :: is_tmp
  type(geometry_t), optional, intent(in)  :: geo
  type(mpi_grp_t),  optional, intent(in)  :: grp !< the group that shares the same data, must contain the domains group
  !
  logical :: is_tmp_, i_am_root
  integer :: comm
#if defined(HAVE_MPI)
  R_TYPE, pointer :: ff_global(:)
#endif
  !
  PUSH_SUB(X(io_function_output))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  is_tmp_ = .false.
  if(present(is_tmp)) is_tmp_ = is_tmp

#if defined(HAVE_MPI)

  i_am_root = .true.
  comm = MPI_COMM_NULL

  if(mesh%parallel_in_domains) then
    SAFE_ALLOCATE(ff_global(1:mesh%np_global))

    !note: here we are gathering data that we won`t write if grp is
    !present, but to avoid it we will have to find out all if the
    !processes are members of the domain line where the root of grp
    !lives
    call X(vec_gather)(mesh%vp, 0, ff_global, ff)

    i_am_root = (mesh%vp%rank == 0)
    comm = mesh%vp%comm
  else
    ff_global => ff
  end if

  if(present(grp)) then
    i_am_root = i_am_root .and. (grp%rank == 0)
    comm = grp%comm
  end if

  if(i_am_root) then
    call X(io_function_output_global)(how, dir, fname, mesh, ff_global, unit, ierr, is_tmp = is_tmp_, geo = geo)
  end if

  if(comm /= MPI_COMM_NULL) then
    ! I have to broadcast the error code
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, comm, mpi_err)
    ! Add a barrier to ensure that the process are synchronized
    call MPI_Barrier(comm, mpi_err)
  end if

  if(mesh%parallel_in_domains) then
    SAFE_DEALLOCATE_P(ff_global)
  else
    nullify(ff_global)
  end if

#else

  ! serial mode
  ASSERT(.not. mesh%parallel_in_domains)
  call X(io_function_output_global)(how, dir, fname, mesh, ff, unit, ierr, is_tmp = is_tmp_, geo = geo)

#endif

  POP_SUB(X(io_function_output))
end subroutine X(io_function_output)


! ---------------------------------------------------------
subroutine X(io_function_output_global) (how, dir, fname, mesh, ff, unit, ierr, is_tmp, geo)
  integer,                    intent(in)  :: how
  character(len=*),           intent(in)  :: dir, fname
  type(mesh_t),               intent(in)  :: mesh
  R_TYPE,                     intent(in)  :: ff(:)  ! ff(mesh%np_global)
  type(unit_t),               intent(in)  :: unit
  integer,                    intent(out) :: ierr
  logical,                    intent(in)  :: is_tmp
  type(geometry_t), optional, intent(in)  :: geo

  character(len=512) :: filename
  character(len=20)  :: mformat, mformat2, mfmtheader
  integer            :: iunit, ip, idir, jj, np_max
  FLOAT              :: x0

  call profiling_in(write_prof, "DISK_WRITE")
  PUSH_SUB(X(io_function_output_global))

  call io_mkdir(dir, is_tmp=is_tmp)

! Define the format; check if code is single precision or double precision
!  FIXME: this may need to be expanded for MAX_DIM > 3
#if defined(SINGLE_PRECISION)
  mformat    = '(4es15.6E3)'
  mformat2   = '(i6,5es15.6E3)'
  mfmtheader = '(a,a7,5a15)'
#else
  mformat    = '(4es23.14E3)'
  mformat2   = '(i6,5es34.24E3)'
  mfmtheader = '(a,a10,5a23)'
#endif

  if(how == 0) then
    message(1) = "Internal error: cannot call io_function with outp%how = 0."
    call messages_fatal(1)
  endif

  np_max = mesh%np_global
  ! should we output boundary points?
  if(iand(how, C_OUTPUT_HOW_BOUNDARY_POINTS)  .ne.0) np_max = mesh%np_part_global

  if(iand(how, C_OUTPUT_HOW_BINARY)    .ne.0) call out_binary()
  if(iand(how, C_OUTPUT_HOW_AXIS_X)    .ne.0) call out_axis (1, 2, 3) ! x ; y=0,z=0
  if(iand(how, C_OUTPUT_HOW_AXIS_Y)    .ne.0) call out_axis (2, 1, 3) ! y ; x=0,z=0
  if(iand(how, C_OUTPUT_HOW_AXIS_Z)    .ne.0) call out_axis (3, 1, 2) ! z ; x=0,y=0
  if(iand(how, C_OUTPUT_HOW_PLANE_X)   .ne.0) call out_plane(1, 2, 3) ! x=0; y; z;
  if(iand(how, C_OUTPUT_HOW_PLANE_Y)   .ne.0) call out_plane(2, 1, 3) ! y=0; x; z;
  if(iand(how, C_OUTPUT_HOW_PLANE_Z)   .ne.0) call out_plane(3, 1, 2) ! z=0; x; y;
  if(iand(how, C_OUTPUT_HOW_MESH_INDEX).ne.0) call out_mesh_index()
  if(iand(how, C_OUTPUT_HOW_DX)        .ne.0) call out_dx()
  if(iand(how, C_OUTPUT_HOW_XCRYSDEN)  .ne.0) then
    call out_xcrysden(.true.)
#ifdef R_TCOMPLEX
    call out_xcrysden(.false.)
#endif
  endif
  if(iand(how, C_OUTPUT_HOW_CUBE)      .ne.0) call out_cube()

  if(iand(how, C_OUTPUT_HOW_MATLAB).ne.0) then
#if defined(R_TCOMPLEX)
    do jj = 1, 3 ! re, im, abs
#else
    do jj = 1, 1 ! only real part
#endif
      if(iand(how, C_OUTPUT_HOW_PLANE_X).ne.0) call out_matlab(how, 1, 2, 3, jj) ! x=0; y; z; 
      if(iand(how, C_OUTPUT_HOW_PLANE_Y).ne.0) call out_matlab(how, 2, 1, 3, jj) ! y=0; x; z;
      if(iand(how, C_OUTPUT_HOW_PLANE_Z).ne.0) call out_matlab(how, 3, 1, 2, jj) ! z=0; x; y;
    end do
    if(iand(how, C_OUTPUT_HOW_MESHGRID).ne.0) then
      do jj = 4, 5 ! meshgrid
        if(iand(how, C_OUTPUT_HOW_PLANE_X).ne.0) call out_matlab(how, 1, 2, 3, jj) ! x=0; y; z; 
        if(iand(how, C_OUTPUT_HOW_PLANE_Y).ne.0) call out_matlab(how, 2, 1, 3, jj) ! y=0; x; z;
        if(iand(how, C_OUTPUT_HOW_PLANE_Z).ne.0) call out_matlab(how, 3, 1, 2, jj) ! z=0; x; y;
      end do
    end if
  end if

#if defined(HAVE_NETCDF)
  if(iand(how, C_OUTPUT_HOW_NETCDF)    .ne.0) call out_netcdf()
#endif

  POP_SUB(X(io_function_output_global))
  call profiling_out(write_prof)

contains

  ! ---------------------------------------------------------
  subroutine out_binary()
    character(len=512) :: workdir

    PUSH_SUB(X(io_function_output_global).out_binary)

    workdir = io_workpath(dir, is_tmp=is_tmp)
    call io_binary_write(trim(workdir)//'/'//trim(fname)//'.obf', mesh%np_global, ff, ierr)

    call profiling_count_transfers(mesh%np_global, ff(1))
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
    iunit = io_open(filename, action='write', is_tmp=is_tmp)

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d1), 'Re', 'Im'
    do ip = 1, np_max
      call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixvect)

      if(ixvect(d2)==0.and.ixvect(d3)==0) then
        xx = units_from_atomic(units_out%length, mesh_x_global(mesh, ip))
        fu = units_from_atomic(unit, ff(ip))
        write(iunit, mformat, iostat=ierr) xx(d1), R_REAL(fu), R_AIMAG(fu)
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
    iunit = io_open(filename, action='write', is_tmp=is_tmp)

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d2), index2axis(d3), 'Re', 'Im'
    write(iunit, mformat)

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
        ip = index_from_coords(mesh%idx, mesh%sb%dim, ixvect_test)
        if(ip /= 0) then 
          call index_to_coords(mesh%idx, mesh%sb%dim, ip, ixvect_test)
          if(ixvect_test(jdim) == 0) exit
        end if
      end do
      ixvect(jdim) = ix
    end do ! loop over dimensions

    ! have found ix such that coordinate d1 is 0 for this value of ix
    ! ixvect is prepared for all dimensions apart from d2 and d3
    
    do iy = mesh%idx%nr(1, d2), mesh%idx%nr(2, d2)
      write(iunit, *)
      do iz = mesh%idx%nr(1, d3), mesh%idx%nr(2, d3)

        ixvect(d2) = iy
        ixvect(d3) = iz
        ip = index_from_coords(mesh%idx, mesh%sb%dim, ixvect)

        if(ip <= mesh%np_global .and. ip > 0) then
          xx = units_from_atomic(units_out%length, mesh_x_global(mesh, ip))
          fu = units_from_atomic(unit, ff(ip))
          write(iunit, mformat, iostat=ierr)  &
            xx(d2), xx(d3), R_REAL(fu), R_AIMAG(fu)
        end if
      end do
    end do

    write(iunit, mformat, iostat=ierr)
    call io_close(iunit)

    POP_SUB(X(io_function_output_global).out_plane)
  end subroutine out_plane


  ! ---------------------------------------------------------
  subroutine out_matlab(how, d1, d2, d3, out_what)
    integer, intent(in) :: how, d1, d2, d3, out_what

    integer :: ix, iy, record_length
    integer :: min_d2, min_d3, max_d2, max_d3
    FLOAT, allocatable :: out_vec(:)
    R_TYPE  :: fu

    PUSH_SUB(X(io_function_output_global).out_matlab)
    
    min_d2 = mesh%idx%nr(1, d2) + mesh%idx%enlarge(d2)
    max_d2 = mesh%idx%nr(2, d2) - mesh%idx%enlarge(d2)
    min_d3 = mesh%idx%nr(1, d3) + mesh%idx%enlarge(d3)
    max_d3 = mesh%idx%nr(2, d3) - mesh%idx%enlarge(d3)    
    
    if(iand(how, C_OUTPUT_HOW_BOUNDARY_POINTS).ne.0) then
      min_d2 = mesh%idx%nr(1, d2)
      max_d2 = mesh%idx%nr(2, d2)
      min_d3 = mesh%idx%nr(1, d3)
      max_d3 = mesh%idx%nr(2, d3)
    end if

    select case(out_what)      
    case(1:3)
      filename = trim(dir)//'/'//trim(fname)//"."//index2axis(d1)//"=0.matlab"
#if defined(R_TCOMPLEX)
      filename = trim(filename)//"."//trim(index2label(out_what))
#endif
    case(4)
      filename = trim(dir)//'/meshgrid.'//index2axis(d1)//"=0."//trim(index2axis(d3)) ! meshgrid d3
    case(5)
      filename = trim(dir)//'/meshgrid.'//index2axis(d1)//"=0."//trim(index2axis(d2)) ! meshgrid d2
    end select

    record_length = (max_d3 - min_d3 + 1)*23  ! 23 because of F23.14 below
    iunit = io_open(filename, action='write', is_tmp=is_tmp, recl=record_length)


    SAFE_ALLOCATE(out_vec(min_d3:max_d3))
    
    do ix = min_d2, max_d2

      out_vec = M_ZERO

      do iy = min_d3, max_d3
        
        select case(d1)
        case(1)
          ip = mesh%idx%lxyz_inv( 0, ix, iy)    ! plane_x
        case(2)
          ip = mesh%idx%lxyz_inv(ix,  0, iy)    ! plane_y
        case(3)
          ip = mesh%idx%lxyz_inv(ix, iy,  0)    ! plane_z
        end select

        select case(out_what)
        case(4)
          out_vec(iy) = mesh%x(ip, d2)      ! meshgrid d2 (this is swapped wrt. 
        case(5)
          out_vec(iy) = mesh%x(ip, d3)      ! meshgrid d3  to the filenames)
        end select

        if (ip < 1 .or. ip > np_max) cycle
        
        fu = units_from_atomic(unit, ff(ip))

        select case(out_what)
        case(1)
          out_vec(iy) = R_REAL(fu)  ! real part
        case(2)
          out_vec(iy) = R_AIMAG(fu) ! imaginary part
        case(3)
          out_vec(iy) = R_ABS(fu)   ! absolute value
        end select
        
      end do

      ! now we write to the disk
      write(iunit,'(32767f23.14)') (out_vec(iy), iy = min_d3, max_d3)

    end do

    SAFE_DEALLOCATE_A(out_vec)
    call io_close(iunit)

    POP_SUB(X(io_function_output_global).out_matlab)
  end subroutine out_matlab


  ! ---------------------------------------------------------
  subroutine out_mesh_index()
    FLOAT :: xx(1:MAX_DIM)
    R_TYPE :: fu

    integer :: idir

    PUSH_SUB(X(io_function_output_global).out_mesh_index)

    iunit = io_open(trim(dir)//'/'//trim(fname)//".mesh_index", action='write', is_tmp=is_tmp)

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
         (units_from_atomic(units_out%length, xx(idir)), idir = 1, 3), R_REAL(fu), R_AIMAG(fu)
    end do

    if(ierr == 0) write(iunit, mformat, iostat=ierr)
    call io_close(iunit)

    POP_SUB(X(io_function_output_global).out_mesh_index)
  end subroutine out_mesh_index


  ! ---------------------------------------------------------
  !> Writes real and imaginary parts
  subroutine out_dx()
    integer :: ix, iy, iz, idir
    FLOAT   :: offset(MAX_DIM)
    character(len=40) :: nitems
    type(cube_t) :: cube
    type(cube_function_t) :: cf

    PUSH_SUB(X(io_function_output_global).out_dx)

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

    ! just for nice formatting of the output
    write(nitems,*) cube%rs_n_global(1)*cube%rs_n_global(2)*cube%rs_n_global(3)
    nitems=trim(adjustl(nitems))

    iunit = io_open(trim(dir)//'/'//trim(fname)//".dx", action='write', is_tmp=is_tmp)

    write(iunit, '(a,3i7)') 'object 1 class gridpositions counts', cube%rs_n_global(1:3)
    write(iunit, '(a,3f12.6)') ' origin', offset(1:3)
    write(iunit, '(a,3f12.6)') ' delta ', (units_from_atomic(units_out%length, &
                                           mesh%spacing(1)*mesh%sb%rlattice_primitive(idir, 1)), idir = 1, 3)
    write(iunit, '(a,3f12.6)') ' delta ', (units_from_atomic(units_out%length, &
                                           mesh%spacing(2)*mesh%sb%rlattice_primitive(idir, 2)), idir = 1, 3)
    write(iunit, '(a,3f12.6)') ' delta ', (units_from_atomic(units_out%length, &
                                           mesh%spacing(3)*mesh%sb%rlattice_primitive(idir, 3)), idir = 1, 3)
    write(iunit, '(a,3i7)') 'object 2 class gridconnections counts', cube%rs_n_global(1:3)
#if defined(R_TREAL)
    write(iunit, '(a,a,a)') 'object 3 class array type float rank 0 items ', nitems, ' data follows'
#else
    write(iunit, '(a,a,a)') 'object 3 class array type float category complex rank 0 items ', nitems, ' data follows'
#endif
    do ix = 1, cube%rs_n_global(1)
      do iy = 1, cube%rs_n_global(2)
        do iz = 1, cube%rs_n_global(3)
          write(iunit,'(2f25.15)') units_from_atomic(unit, cf%X(RS)(ix, iy, iz))
        end do
      end do
    end do
    write(iunit, '(a)') 'object "regular positions regular connections" class field'
    write(iunit, '(a)') ' component "positions" value 1'
    write(iunit, '(a)') ' component "connections" value 2'
    write(iunit, '(a)') ' component "data" value 3'
    write(iunit, '(a)') 'end'

    call io_close(iunit)
    call cube_end(cube)
    call X(cube_function_free_RS)(cube, cf)

    POP_SUB(X(io_function_output_global).out_dx)
  end subroutine out_dx


  ! ---------------------------------------------------------
  !> see http://local.wasp.uwa.edu.au/~pbourke/dataformats/cube/
  !! Writes only real part
  subroutine out_cube()
    integer :: ix, iy, iz, idir, idir2, iatom
    FLOAT   :: offset(MAX_DIM)
    type(cube_t) :: cube
    type(cube_function_t) :: cf

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

    iunit = io_open(trim(dir)//'/'//trim(fname)//".cube", action='write', is_tmp=is_tmp)

    write(iunit, '(2a)') 'Generated by octopus ', trim(conf%version)
    write(iunit, '(4a)') 'svn: ', trim(conf%latest_svn), " build: ",  trim(conf%build_time)
    write(iunit, '(i5,3f12.6)') geo%natoms, offset(1:3)
    do idir = 1, 3
      write(iunit, '(i5,3f12.6)') cube%rs_n_global(idir), (units_from_atomic(units_out%length, &
        mesh%spacing(idir)*mesh%sb%rlattice_primitive(idir2, idir)), idir2 = 1, 3)
    end do
    do iatom = 1, geo%natoms
      write(iunit, '(i5,4f12.6)') int(species_z(geo%atom(iatom)%spec)),  M_ZERO, &
        (units_from_atomic(units_out%length, geo%atom(iatom)%x(idir)), idir = 1, 3)
    end do

    do ix = 1, cube%rs_n_global(1)
      do iy = 1, cube%rs_n_global(2)
        do iz = 1, cube%rs_n_global(3)
          write(iunit,'(e14.6)', advance='no') units_from_atomic(unit, R_REAL(cf%X(RS)(ix, iy, iz)))
          if(mod(iz-1, 6) == 5)  write(iunit, '(1x)')
        end do
        write(iunit, '(1x)')
      end do
    end do

    call io_close(iunit)
    call cube_end(cube)
    call X(cube_function_free_RS)(cube, cf)

    POP_SUB(X(io_function_output_global).out_cube)
  end subroutine out_cube


  ! ---------------------------------------------------------
  !> For format specification see:
  !! http://www.xcrysden.org/doc/XSF.html#__toc__11
  !! XCrySDen can only read 3D output, though it could be
  !! extended to plot a function on a 2D plane.
  !! Writes real part unless write_real = false and called in complex version
  subroutine out_xcrysden(write_real)
    logical, intent(in) :: write_real

    integer :: ix, iy, iz, idir2, ix2, iy2, iz2, my_n(3)
    FLOAT :: lattice_vectors(3,3)
    FLOAT :: offset(3)
    type(cube_t) :: cube
    type(cube_function_t) :: cf
    character(len=80) :: fname_ext

    PUSH_SUB(X(io_function_output_global).out_xcrysden)

#ifdef R_TCOMPLEX
    if(write_real) then
      fname_ext = trim(fname) // '.real'
    else
      fname_ext = trim(fname) // '.imag'
    endif
#else
    fname_ext = trim(fname)
#endif

    if(mesh%sb%dim .ne. 3) then
      write(message(1), '(a)') 'Cannot output function '//trim(fname_ext)//' in XCrySDen format except in 3D.'
      call messages_warning(1)
      return
    endif

    ! put values in a nice cube
    call cube_init(cube, mesh%idx%ll, mesh%sb)
    call cube_function_null(cf)
    call X(cube_function_alloc_RS) (cube, cf)
    call X(mesh_to_cube) (mesh, ff, cube, cf)

    ! Note that XCrySDen uses "general" not "periodic" grids
    ! mesh%idx%ll is "general" in aperiodic directions,
    ! but "periodic" in periodic directions.
    ! Making this assignment, the output grid is entirely "general"
    my_n(1:mesh%sb%periodic_dim) = mesh%idx%ll(1:mesh%sb%periodic_dim) + 1
    my_n(mesh%sb%periodic_dim + 1:3) = mesh%idx%ll(mesh%sb%periodic_dim + 1:3)

    ! This differs from mesh%sb%rlattice if it is not an integer multiple of the spacing
    do idir = 1, 3
      do idir2 = 1, 3
        lattice_vectors(idir, idir2) = mesh%spacing(idir) * (my_n(idir) - 1) * mesh%sb%rlattice_primitive(idir2, idir)
      enddo
    enddo
    
    iunit = io_open(trim(dir)//'/'//trim(fname_ext)//".xsf", action='write', is_tmp=is_tmp)

    ASSERT(present(geo))
    call write_xsf_geometry(iunit, geo, mesh)

    write(iunit, '(a)') 'BEGIN_BLOCK_DATAGRID3D'
    write(iunit, '(4a)') 'units: coords = ', trim(units_abbrev(units_out%length)), &
                            ', function = ', trim(units_abbrev(unit))
    write(iunit, '(a)') 'DATAGRID_3D_function'
    write(iunit, '(3i7)') my_n(1:3)
    write(iunit, '(a)') '0.0 0.0 0.0'

    do idir = 1, 3
      write(iunit, '(3f12.6)') (units_from_atomic(units_out%length, &
        lattice_vectors(idir2, idir)), idir2 = 1, 3)
    enddo

    do iz = 1, my_n(3)
      do iy = 1, my_n(2)
        do ix = 1, my_n(1)
          ! this is about "general" grids also
          if (ix == mesh%idx%ll(1) + 1) then
            ix2 = 1
          else
            ix2 = ix
          endif

          if (iy == mesh%idx%ll(2) + 1) then
            iy2 = 1
          else
            iy2 = iy
          endif

          if (iz == mesh%idx%ll(3) + 1) then
            iz2 = 1
          else
            iz2 = iz
          endif

#ifdef R_TCOMPLEX
          if(.not. write_real) then
            write(iunit,'(2f25.15)') aimag(units_from_atomic(unit, cf%X(RS)(ix2, iy2, iz2)))
          else
            write(iunit,'(2f25.15)') real(units_from_atomic(unit, cf%X(RS)(ix2, iy2, iz2)))
          endif
#else
          write(iunit,'(2f25.15)') units_from_atomic(unit, cf%X(RS)(ix2, iy2, iz2))
#endif
        end do
      end do
    end do

    write(iunit, '(a)') 'END_DATAGRID3D'
    write(iunit, '(a)') 'END_BLOCK_DATAGRID3D'

    call io_close(iunit)
    call cube_end(cube)
    call X(cube_function_free_RS)(cube, cf)

    POP_SUB(X(io_function_output_global).out_xcrysden)
  end subroutine out_xcrysden


#if defined(HAVE_NETCDF)
  ! ---------------------------------------------------------
  subroutine out_netcdf()

    type(cube_t) :: cube
    type(cube_function_t) :: cf



    PUSH_SUB(X(io_function_output_global).out_netcdf)


    ! put values in a nice cube
    call cube_init(cube, mesh%idx%ll, mesh%sb)
    call cube_function_null(cf)
    call X(cube_function_alloc_RS) (cube, cf)
    call X(mesh_to_cube) (mesh, ff, cube, cf)

    filename = io_workpath(trim(dir)//'/'//trim(fname)//".ncdf", is_tmp=is_tmp)

     
    call X(out_cf_netcdf)(filename, ierr, cf, cube, mesh%sb%dim,& 
          units_from_atomic(units_out%length, mesh%spacing(1:mesh%sb%dim)), transpose = .true.)

    call cube_end(cube)
    call X(cube_function_free_RS)(cube, cf)

    POP_SUB(X(io_function_output_global).out_netcdf)
  end subroutine out_netcdf



#endif /*defined(HAVE_NETCDF)*/

end subroutine X(io_function_output_global)

#if defined(HAVE_NETCDF)
  ! --------------------------------------------------------- 
  !>  Writes a cube_function in netcdf format
  ! ---------------------------------------------------------
  subroutine X(out_cf_netcdf)(filename, ierr, cf, cube, sb_dim, spacing, transpose)
    character(len=*),      intent(in) :: filename        !< the file name
    integer,               intent(out):: ierr            !< error message   
    type(cube_function_t), intent(in) :: cf              !< the cube_function to be written 
    type(cube_t),          intent(in) :: cube            !< the underlying cube mesh
    integer,               intent(in) :: sb_dim          !< the simulation box dimensions aka sb%dim
    FLOAT,                 intent(in) :: spacing(:)      !< the mesh spacing already converted to units_out
    logical,               intent(in) :: transpose       !< whether we want the function cf(x,y,z) to be saved as cf(z,y,x)


    integer :: ncid, status, data_id, pos_id, dim_min
    integer :: dim_data_id(3), dim_pos_id(2)

    REAL_SINGLE :: pos(2, 3)
    FLOAT, allocatable :: xx(:, :, :)
    
#if defined(R_TCOMPLEX)
    integer :: data_im_id
#endif

    PUSH_SUB(X(out_cf_netcdf))

    ierr = 0


    status = nf90_create(trim(filename), NF90_CLOBBER, ncid)
    if(status.ne.NF90_NOERR) then
      ierr = 2
      POP_SUB(X(out_cf_netcdf))
      return
    end if

    ! dimensions
    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_1", cube%rs_n_global(3), dim_data_id(1))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_2", cube%rs_n_global(2), dim_data_id(2))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_3", cube%rs_n_global(1), dim_data_id(3))
      call ncdf_error('nf90_der_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "pos_1", 2, dim_pos_id(1))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "pos_2", 3, dim_pos_id(2))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    dim_min = 3 - sb_dim + 1

#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_FLOAT, dim_data_id(dim_min:3), data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if
#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_DOUBLE, dim_data_id(dim_min:3), data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if
#endif
#if defined(R_TCOMPLEX)
#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_FLOAT, dim_data_id(dim_min:3), data_im_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if

#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_DOUBLE, dim_data_id(dim_min:3), data_im_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if
#endif
#endif
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "pos", NF90_FLOAT,  dim_pos_id,  pos_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if

    ! attributes
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_id, "field", "rdata, scalar")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    end if
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_id, "positions", "pos, regular")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    end if
#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_im_id, "field", "idata, scalar")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    end if
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_im_id, "positions", "pos, regular")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    end if
#endif

    ! end definitions
    status = nf90_enddef (ncid)

    ! data
    pos(:,:) = M_ZERO
    pos(1, 1:sb_dim) = &
      real( - (cube%rs_n_global(1:sb_dim) - 1)/2*spacing(1:sb_dim), 4)
    pos(2, 1:sb_dim) = real(spacing(1:sb_dim), 4)

    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, pos_id, pos(:,:))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if

    SAFE_ALLOCATE(xx(1:cube%rs_n_global(3), 1:cube%rs_n_global(2), 1:cube%rs_n_global(1)))
#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      if (transpose) then
        call transpose3(real(cf%X(RS), REAL_PRECISION), xx)
      else 
        xx = real(cf%X(RS))
      end if    
      call write_variable(ncid, data_id, status, sb_dim, xx)
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
    if(status == NF90_NOERR) then
      if ( transpose ) then 
        call transpose3(aimag(cf%X(RS)), xx)
      else
        xx = aimag(cf%X(RS))
      end if
      call write_variable(ncid, data_im_id, status, sb_dim, xx)
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
#else
    if(status == NF90_NOERR) then
      if ( transpose ) then 
        call transpose3(cf%X(RS), xx)
      else             
        xx=cf%X(RS)
      end if 
      call write_variable(ncid, data_id, status, sb_dim, xx)
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
#endif
    SAFE_DEALLOCATE_A(xx)

    ! close
    status = nf90_close(ncid)

    POP_SUB(X(out_cf_netcdf))
    
    contains

    ! ---------------------------------------------------------
    subroutine write_variable(ncid, data_id, status, sb_dim, xx)
      integer, intent(in)  :: ncid, data_id
      integer, intent(out) :: status
      integer, intent(in)  :: sb_dim
      FLOAT,   intent(in)  :: xx(:,:,:)

      PUSH_SUB(X(out_cf_netcdf).write_variable)

      select case(sb_dim)
      case(1)
        status = nf90_put_var (ncid, data_id, xx(1,1,:))
      case(2)
        status = nf90_put_var (ncid, data_id, xx(1,:,:))
      case(3)
        status = nf90_put_var (ncid, data_id, xx)
      end select
    
      POP_SUB(X(out_cf_netcdf).write_variable)
    end subroutine write_variable
    

  end subroutine X(out_cf_netcdf)



#endif /*defined(HAVE_NETCDF)*/



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
