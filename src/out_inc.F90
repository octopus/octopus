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
! Reads a function from file filename, and puts it into f. The input file
! may be a "plain" file (no extension), or a netcdf file ".ncdf" extension.
! On output, ierr signals how everything went:
! ierr > 0 => Error. The function f was not read:
!              1 : illegal filename (must have no extension, or ".ncdf" extension.
!              2 : file could not be succesfully opened.
!              3 : file opened, but error reading.
!              4 : The number of points/mesh dimensions do not coincide.
!              5 : NetCDF error (one or several warnings are emitted)
! ierr = 0 => Success.
! ierr < 0 => Success, but some kind of type conversion was necessary. The
!             of ierr is then:
!             -1 : function in file is real, sp.
!             -2 : function in file is complex, sp.
!             -3 : function in file is real, dp.
!             -4 : function in file is complex, dp.
! ---------------------------------------------------------

subroutine X(input_function)(filename, m, f, ierr, is_tmp)
  character(len=*),     intent(in)  :: filename
  type(mesh_t),      intent(in)  :: m
  R_TYPE,               intent(out) :: f(:)
  integer,              intent(out) :: ierr
  logical, optional,    intent(in)  :: is_tmp

  logical :: is_tmp_ = .false.

#if defined(HAVE_MPI)
  integer             :: mpi_err
  R_TYPE, allocatable :: f_global(:)
#endif

  call push_sub('out_inc.Xinput_function')

  if(present(is_tmp)) is_tmp_ = is_tmp

  if(m%parallel_in_domains) then
#if defined(HAVE_MPI)
    ! Only root reads. Therefore, only root needs a buffer
    ! f_global for the whole function.
    ALLOCATE(f_global(1), 1)
    if(mpi_grp_is_root(m%mpi_grp)) then
      deallocate(f_global); ALLOCATE(f_global(m%np_global), m%np_global)
      call X(input_function_global)(filename, m, f_global, ierr, is_tmp_)
    end if
    if(in_debug_mode) call write_debug_newlines(2)

    ! Only root knows if the file was succesfully read.
    ! Now, it tells everybody else.
    call MPI_Debug_IN(m%vp%comm, C_MPI_BCAST)
    call MPI_Bcast(ierr, 1, MPI_INTEGER, m%vp%root, m%vp%comm, mpi_err)
    call MPI_Debug_OUT(m%vp%comm, C_MPI_BCAST)

    ! Only scatter, when successfully read the file(s).
    if(ierr.le.0) then
      call X(vec_scatter)(m%vp, f_global, f)
    end if

    deallocate(f_global)
#else
    ASSERT(.false.) ! internal error
#endif
  else
    call X(input_function_global)(filename, m, f, ierr, is_tmp_)
  end if

  call pop_sub()

end subroutine X(input_function)


! ---------------------------------------------------------
subroutine X(input_function_global)(filename, m, f, ierr, is_tmp)
  character(len=*),     intent(in)  :: filename
  type(mesh_t),      intent(in)  :: m
  R_TYPE,               intent(out) :: f(:)
  integer,              intent(out) :: ierr
  logical,              intent(in)  :: is_tmp

  integer :: iunit, i, function_kind, file_kind

#if defined(HAVE_NETCDF)
  type(X(cf)) :: c
#if defined(R_TCOMPLEX)
  type(dcf) :: re, im
#endif
#endif

  call push_sub('out_inc.Xinput_function_global')

  ierr = 0
  function_kind = X(output_kind)*kind(f(1)) ! +4 for real, single; +8 for real, double;
  ! -4 for complex, single, -8 for real, double

  select case(trim(io_get_extension(filename)))
  case("")
     call plain()
#if defined(HAVE_NETCDF)
  case("ncdf")
#if defined(R_TCOMPLEX)
     call X(cf_new)(m%l, c); call dcf_new(m%l, re); call dcf_new(m%l, im)
     call X(cf_alloc_RS)(c); call dcf_alloc_RS(re); call dcf_alloc_RS(im)
     call dx_cdf()
     c%RS = re%RS + M_zI*im%RS
     call X(cf2mf) (m, c, f)
     call X(cf_free)(c); call dcf_free(re); call dcf_free(im)
#else
     call X(cf_new)(m%l, c)
     call X(cf_alloc_RS)(c)
     call dx_cdf()
     call X(cf2mf) (m, c, f)
     call X(cf_free)(c)
#endif
#endif

  case default
     ierr = 1
  end select

  call pop_sub()

contains


  ! ---------------------------------------------------------
  subroutine plain()
    integer                 :: file_np
    real(4),    allocatable :: rs(:)
    real(8),    allocatable :: rd(:)
    complex(4), allocatable :: cs(:)
    complex(8), allocatable :: cd(:)

    iunit = io_open(filename, action='read', status='old', form='unformatted', die=.false., is_tmp=is_tmp)

    if(iunit< 0) then
      ierr = 2
      return
    end if

    read(unit = iunit, iostat = i) file_kind, file_np
    if(i.ne.0) then
      ierr = 3
    else if (file_np .ne. m%np_global) then
      ierr = 4
    end if
    
    if(ierr==0) then
      if(file_kind == function_kind) then
        read(unit = iunit) f(1:m%np_global)
        
      else ! Adequate conversions....
        select case(file_kind)
        case(doutput_kind*4) ! Real, single precision
          ALLOCATE(rs(m%np_global), m%np_global)
          read(unit = iunit) rs(1:m%np_global)
          f = rs
          deallocate(rs)
          ierr = -1

        case(zoutput_kind*4) ! Complex, single precision
          ALLOCATE(cs(m%np_global), m%np_global)
          read(unit = iunit) cs(1:m%np_global)
          f = cs
          deallocate(cs)
          ierr = -2

        case(doutput_kind*8) ! Real, double precision
          ALLOCATE(rd(m%np_global), m%np_global)
          read(unit = iunit) rd(1:m%np_global)
          f = rd
          deallocate(rd)
          ierr = -3

        case(zoutput_kind*8) ! Complex, double precision
          ALLOCATE(cd(m%np_global), m%np_global)
          read(unit = iunit) cd(1:m%np_global)
          f = cd
          deallocate(cd)
          ierr = -4
          
        end select
      end if
    end if

    call io_close(iunit)
  end subroutine plain


  ! ---------------------------------------------------------
#if defined(HAVE_NETCDF)

  ! ---------------------------------------------------------
  subroutine dx_cdf()
    integer :: ncid, ndims, nvars, natts, status, data_id, data_im_id, pos_id, &
         dim_data_id(3), dim_pos_id(2), ndim(3), xtype
    real(r4)           :: pos(2, 3)
    logical            :: function_is_complex = .false.
    character(len=512) :: file

    file = io_workpath(filename, is_tmp=is_tmp)

    status = nf90_open(trim(file), NF90_WRITE, ncid)
    if(status.ne.NF90_NOERR) then
       ierr = 2
       return
    end if

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
       status = nf90_inquire_dimension (ncid, dim_data_id(1), len = ndim(1))
       call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if(status == NF90_NOERR) then
       status = nf90_inquire_dimension (ncid, dim_data_id(2), len = ndim(2))
       call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if(status == NF90_NOERR) then
       status = nf90_inquire_dimension (ncid, dim_data_id(3), len = ndim(3))
       call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    end if
    if((ndim(1) .ne. c%n(1)) .or. &
         (ndim(2) .ne. c%n(2)) .or. &
         (ndim(3) .ne. c%n(3))) then
       ierr = 12; return
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
       case(4);  ierr = -1
       case(-4); ierr = -2
       case(8);  ierr = -3
       case(-8); ierr = -4
       end select
    end if

#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
       status = nf90_get_var (ncid, data_id, re%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
       call ncdf_error('nf90_get_var', status, file, ierr)
    end if
    if(file_kind<0) then
       if(status == NF90_NOERR) then
          status = nf90_get_var (ncid, data_im_id, im%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
          call ncdf_error('nf90_get_var', status, file, ierr)
       end if
    end if
#else
    if(status == NF90_NOERR) then
       status = nf90_get_var (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
       call ncdf_error('nf90_get_var', status, file, ierr)
    end if
#endif

    status = nf90_close(ncid)
  end subroutine dx_cdf

#endif

end subroutine X(input_function_global)


! ---------------------------------------------------------
subroutine X(output_function) (how, dir, fname, mesh, sb, f, u, ierr, is_tmp)
  integer,              intent(in)  :: how
  character(len=*),     intent(in)  :: dir, fname
  type(mesh_t),      intent(in)  :: mesh
  type(simul_box_t), intent(in)  :: sb
  R_TYPE,               intent(in)  :: f(:)  ! f(mesh%np)
  FLOAT,                intent(in)  :: u
  integer,              intent(out) :: ierr
  logical, optional,    intent(in)  :: is_tmp

  logical :: is_tmp_ = .false.

#if defined(HAVE_MPI)
  R_TYPE, allocatable :: f_global(:)
  integer             :: mpi_err
#endif

  call push_sub('out_inc.Xoutput_function')

  if(present(is_tmp)) is_tmp_ = is_tmp

  if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
    ALLOCATE(f_global(mesh%np_global), mesh%np_global)

    call X(vec_gather)(mesh%vp, f_global, f)

    if(mesh%vp%rank.eq.mesh%vp%root) then
      call X(output_function_global)(how, dir, fname, mesh, sb, f_global, u, ierr, is_tmp_)
    end if

    ! I have to broadcast the error code
    call MPI_Bcast(ierr, 1, MPI_INTEGER, mesh%vp%root, mesh%vp%comm, mpi_err)

    if(in_debug_mode) call write_debug_newlines(2)

    deallocate(f_global)
#else
    ASSERT(.false.)
#endif
  else
    call X(output_function_global)(how, dir, fname, mesh, sb, f, u, ierr, is_tmp_)
  end if

  call pop_sub()
end subroutine X(output_function)


! ---------------------------------------------------------
subroutine X(output_function_global) (how, dir, fname, m, sb, f, u, ierr, is_tmp)
  integer,              intent(in)  :: how
  character(len=*),     intent(in)  :: dir, fname
  type(mesh_t),         intent(in)  :: m
  type(simul_box_t),    intent(in)  :: sb
  R_TYPE,               intent(in)  :: f(:)  ! f(m%np_global)
  FLOAT,                intent(in)  :: u
  integer,              intent(out) :: ierr
  logical,              intent(in)  :: is_tmp


  character(len=256) :: filename
  character(len=20)  :: mformat, mformat2, mfmtheader
  integer            :: iunit, i, j, np_max
  FLOAT              :: x0
  logical            :: gnuplot_mode = .false.

  call push_sub('out_inc.Xoutput_function_global')

  call io_mkdir(dir)

! Define the format; check if code is single precision or double precision
#if defined(SINGLE_PRECISION)
    mformat    = '(4es15.6)'
    mformat2   = '(i6,5es15.6)'
    mfmtheader = '(a,a7,5a15)'
#else
    mformat    = '(4es23.14)'
    mformat2   = '(i6,5es34.24)'
    mfmtheader = '(a,a10,5a23)'
#endif

  if(iand(how, output_gnuplot)   .ne.0) then
    gnuplot_mode = .true.
    mformat    = '(4f23.14)'
    mformat2   = '(i6,5f34.24)'
  end if

  np_max = m%np_global
  ! should we output boundary points?
  if(iand(how, boundary_points)  .ne.0) np_max = m%np_part_global

  if(iand(how, output_plain)     .ne.0) call plain()
  if(iand(how, output_axis_x)    .ne.0) call out_axis (1, 2, 3) ! x ; y=0,z=0
  if(iand(how, output_axis_y)    .ne.0) call out_axis (2, 1, 3) ! y ; x=0,z=0
  if(iand(how, output_axis_z)    .ne.0) call out_axis (3, 1, 2) ! z ; x=0,y=0
  if(iand(how, output_plane_x)   .ne.0) call out_plane(1, 2, 3) ! x=0; y; z;
  if(iand(how, output_plane_y)   .ne.0) call out_plane(2, 1, 3) ! y=0; x; z;
  if(iand(how, output_plane_z)   .ne.0) call out_plane(3, 1, 2) ! z=0; x; y;
  if(iand(how, output_mesh_index).ne.0) call out_mesh_index()
  if(iand(how, output_dx)        .ne.0) call dx()

  if(iand(how, output_matlab).ne.0) then
    do j = 1, 3 ! re, im, abs
      if(iand(how, output_plane_x).ne.0) call out_matlab(how, 1, 2, 3, j) ! x=0; y; z; 
      if(iand(how, output_plane_y).ne.0) call out_matlab(how, 2, 1, 3, j) ! y=0; x; z;
      if(iand(how, output_plane_z).ne.0) call out_matlab(how, 3, 1, 2, j) ! z=0; x; y;
    end do
    if(iand(how, output_meshgrid).ne.0) then
      do j = 4, 5 ! meshgrid
        if(iand(how, output_plane_x).ne.0) call out_matlab(how, 1, 2, 3, j) ! x=0; y; z; 
        if(iand(how, output_plane_y).ne.0) call out_matlab(how, 2, 1, 3, j) ! y=0; x; z;
        if(iand(how, output_plane_z).ne.0) call out_matlab(how, 3, 1, 2, j) ! z=0; x; y;
      end do
    end if
  end if

#if defined(HAVE_NETCDF)
  if(iand(how, output_dx_cdf)    .ne.0) call dx_cdf()
#endif

  call pop_sub()

contains


  ! ---------------------------------------------------------
  subroutine plain()

    iunit = io_open(trim(dir)//'/'//trim(fname), action='write', form='unformatted', is_tmp=is_tmp)

    write(unit=iunit, iostat=ierr) X(output_kind)*kind(f(1)), m%np_global
    write(unit=iunit, iostat=ierr) f(1:m%np_global)
    call io_close(iunit)

  end subroutine plain


  ! ---------------------------------------------------------
  subroutine out_axis(d1, d2, d3)
    integer, intent(in) :: d1, d2, d3

    filename = trim(dir)//'/'//trim(fname)//"."//index2axis(d2)//"=0,"//index2axis(d3)//"=0"
    iunit = io_open(filename, action='write', is_tmp=is_tmp)

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d1), 'Re', 'Im'
    do i = 1, np_max
      if(m%Lxyz(i, d2)==0.and.m%Lxyz(i, d3)==0) then     
        write(iunit, mformat, iostat=ierr) m%x_global(i, d1), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do

    call io_close(iunit)

  end subroutine out_axis


  ! ---------------------------------------------------------
  subroutine out_plane(d1, d2, d3)
    integer, intent(in) :: d1, d2, d3

    filename = trim(dir)//'/'//trim(fname)//"."//index2axis(d1)//"=0"
    iunit = io_open(filename, action='write', is_tmp=is_tmp)

    write(iunit, mfmtheader, iostat=ierr) '#', index2axis(d2), index2axis(d3), 'Re', 'Im'
    x0 = m%x_global(1, d2)
    if(gnuplot_mode) write(iunit, mformat)

    do i = 1, np_max
      if(gnuplot_mode.and.x0 /= m%x_global(i, d2)) then
        write(iunit, mformat, iostat=ierr)  ! write extra lines for gnuplot grid mode
        x0 = m%x_global(i, d2)
      end if
      if(m%Lxyz(i, d1) == 0) then  ! are we in the plane?
          write(iunit, mformat, iostat=ierr)  &
            m%x_global(i, d2), m%x_global(i, d3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do

    if(gnuplot_mode) write(iunit, mformat, iostat=ierr)
    call io_close(iunit)

  end subroutine out_plane


  ! ---------------------------------------------------------
  subroutine out_matlab(how, d1, d2, d3, out_what)
    integer, intent(in) :: how, d1, d2, d3, out_what

    integer :: ix, iy, record_length
    integer :: min_d2, min_d3, max_d2, max_d3
    FLOAT, allocatable :: out_vec(:)
    
    min_d2 = m%nr(1, d2) + m%enlarge(d2)
    max_d2 = m%nr(2, d2) - m%enlarge(d2)
    min_d3 = m%nr(1, d3) + m%enlarge(d3)
    max_d3 = m%nr(2, d3) - m%enlarge(d3)    
    
    if(iand(how, boundary_points).ne.0) then
      min_d2 = m%nr(1, d2); max_d2 = m%nr(2, d2); 
      min_d3 = m%nr(1, d3); max_d3 = m%nr(2, d3); 
    end if

    select case(out_what)      
    case(1:3);  filename = &
      trim(dir)//'/'//trim(fname)//"."//index2axis(d1)//"=0.matlab."//trim(index2label(out_what))
    case(4);    filename = &
      trim(dir)//'/meshgrid.'//index2axis(d1)//"=0."//trim(index2axis(d3))      ! meshgrid d3
    case(5);    filename = &
      trim(dir)//'/meshgrid.'//index2axis(d1)//"=0."//trim(index2axis(d2))      ! meshgrid d2
    end select

    record_length = (max_d3 - min_d3 + 1)*23  ! 23 because of F23.14 below
    iunit = io_open(filename, action='write', is_tmp=is_tmp, recl=record_length)


    ALLOCATE(out_vec(min_d3:max_d3), max_d3-min_d3 + 1)
    
    do ix = min_d2, max_d2

      out_vec = M_ZERO

      do iy = min_d3, max_d3
        
        select case(d1)
        case(1); i = m%Lxyz_inv( 0, ix, iy)    ! plane_x
        case(2); i = m%Lxyz_inv(ix,  0, iy)    ! plane_y
        case(3); i = m%Lxyz_inv(ix, iy,  0)    ! plane_z
        end select

        select case(out_what)
        case(4); out_vec(iy) = m%x(i, d2)      ! meshgrid d2 (this is swapped wrt. 
        case(5); out_vec(iy) = m%x(i, d3)      ! meshgrid d3  to the filenames)
        end select

        if (i < 1 .or. i > np_max) cycle
        
        select case(out_what)
        case(1); out_vec(iy) = R_REAL(f(i))/u  ! real part
        case(2); out_vec(iy) = R_AIMAG(f(i))/u ! imaginary part
        case(3); out_vec(iy) = R_ABS(f(i))/u   ! absolute value
        end select
        
      end do

      ! now we write to the disk
      write(iunit,'(32767f23.14)') (out_vec(iy), iy = min_d3, max_d3)

    end do

    deallocate(out_vec)
    call io_close(iunit)

  end subroutine out_matlab


  ! ---------------------------------------------------------
  subroutine out_mesh_index()

    iunit = io_open(trim(dir)//'/'//trim(fname)//".mesh_index", action='write', is_tmp=is_tmp)

    write(iunit, mfmtheader, iostat=ierr) '#', 'Index', 'x', 'y', 'z', 'Re', 'Im'
    x0 = m%x_global(1,1)
    if(ierr == 0.and.gnuplot_mode) write(iunit, mformat, iostat=ierr)

    do i= 1, np_max
       if (ierr == 0.and.gnuplot_mode.and.x0 /= m%x_global(i, 1)) then
          write(iunit, mformat, iostat=ierr)      ! write extra lines for gnuplot grid mode
          x0 = m%x_global(i, 1)
       end if
       if(ierr==0) write(iunit, mformat2, iostat=ierr) i, m%x_global(i,1),  &
            m%x_global(i,2), m%x_global(i,3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
    end do

    if(ierr == 0.and.gnuplot_mode) write(iunit, mformat, iostat=ierr)
    call io_close(iunit)
  end subroutine out_mesh_index


  ! ---------------------------------------------------------
  subroutine dx()
    integer :: ix, iy, iz
    FLOAT   :: offset(3)
    character(len=40) :: nitems
    type(X(cf)) :: c

    ! put values in a nice cube
    call X(cf_new) (m%l, c)
    call X(cf_alloc_RS) (c)
    call X(mf2cf) (m, f, c)

    ! the offset is different in periodic directions
    do i = 1, sb%periodic_dim
      offset(i)=-(c%n(i))/2 * m%h(i) / units_out%length%factor
    end do
    do i = sb%periodic_dim+1, 3
      offset(i)=-(c%n(i) - 1)/2 * m%h(i) / units_out%length%factor
    end do

! just for nice formatting of the output
    write(nitems,*)c%n(1)*c%n(2)*c%n(3)
    nitems=trim(adjustl(nitems))

    iunit = io_open(trim(dir)//'/'//trim(fname)//".dx", action='write', is_tmp=is_tmp)

    write(iunit, '(a,3i7)') 'object 1 class gridpositions counts', c%n(:)
    write(iunit, '(a,3f12.6)') ' origin', offset(:)
    write(iunit, '(a,f12.6,a)') ' delta ',m%h(1) / units_out%length%factor, '    0.000000    0.000000'
    write(iunit, '(a,f12.6,a)') ' delta     0.000000',m%h(2) / units_out%length%factor, '    0.000000'
    write(iunit, '(a,f12.6)') ' delta     0.000000    0.000000',m%h(3) / units_out%length%factor
    write(iunit, '(a,3i7)') 'object 2 class gridconnections counts',c%n(:)
#if defined(R_TREAL)
    write(iunit, '(a,a,a)') 'object 3 class array type float rank 0 items ',nitems,' data follows'
#else
    write(iunit, '(a,a,a)') 'object 3 class array type float category complex rank 0 items ',&
                             nitems,' data follows'
#endif
    do ix = 1, c%n(1)
      do iy = 1, c%n(2)
        do iz = 1, c%n(3)
          write(iunit,'(2e20.10)') c%RS(ix, iy, iz)
        end do
      end do
    end do
    write(iunit, '(a)') 'object "regular positions regular connections" class field'
    write(iunit, '(a)') ' component "positions" value 1'
    write(iunit, '(a)') ' component "connections" value 2'
    write(iunit, '(a)') ' component "data" value 3'
    write(iunit, '(a)') 'end'

    call io_close(iunit)
    call X(cf_free) (c)

  end subroutine dx


#if defined(HAVE_NETCDF)
  ! ---------------------------------------------------------
  subroutine dx_cdf()
    integer :: ncid, status, data_id, data_im_id, pos_id, dim_data_id(3), dim_pos_id(2)
    real(r4) :: pos(2, 3)
    type(X(cf)) :: c

    ierr = 0

    ! put values in a nice cube
    call X(cf_new) (m%l, c)
    call X(cf_alloc_RS) (c)
    call X(mf2cf) (m, f, c)

    filename = io_workpath(trim(dir)//'/'//trim(fname)//".ncdf", is_tmp=is_tmp);

    status = nf90_create(trim(filename), NF90_CLOBBER, ncid)
    if(status.ne.NF90_NOERR) then
      ierr = 2
      return
    end if

    ! dimensions
    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_1", c%n(1), dim_data_id(1))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_2", c%n(2), dim_data_id(2))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    end if

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_3", c%n(3), dim_data_id(3))
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

#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_FLOAT, dim_data_id, data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if
#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_DOUBLE, dim_data_id, data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if
#endif
#if defined(R_TCOMPLEX)
#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_FLOAT, dim_data_id, data_im_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    end if

#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_DOUBLE, dim_data_id, data_im_id)
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
    pos(1,:) = real(-(c%n(:) - 1)/2 * m%h(:) / units_out%length%factor, 4)
    pos(2,:) = real(m%h(:) / units_out%length%factor, 4)

    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, pos_id, pos(:,:))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if

#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_id, real(c%RS, PRECISION), map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_im_id, aimag(c%RS), map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
#else
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    end if
#endif

    ! close
    status = nf90_close(ncid)
    call X(cf_free) (c)

  end subroutine dx_cdf

#endif

end subroutine X(output_function_global)
