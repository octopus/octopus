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
subroutine X(input_function)(filename, m, f, ierr)
  character(len=*), intent(in)  :: filename
  type(mesh_type),  intent(in)  :: m
  R_TYPE,           intent(out) :: f(:)
  integer,          intent(out) :: ierr

  integer :: iunit, i, function_kind, file_kind

  ierr = 0
  function_kind = X(output_kind)*kind(f(1)) ! +4 for real, single; +8 for real, double;
                                            ! -4 for complex, single, -8 for real, double

  select case(trim(get_extension(filename)))
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

contains


  ! ---------------------------------------------------------
  subroutine plain()
    integer                 :: file_np
    real(4),    allocatable :: rs(:)
    real(8),    allocatable :: rd(:)
    complex(4), allocatable :: cs(:)
    complex(8), allocatable :: cd(:)

    call io_assign(iunit)
    iunit = io_open(filename, action='read', status='old', form='unformatted', die=.false.)
    
    if(iunit< 0) then
      ierr = 2
      return
    end if

    read(unit = iunit, iostat = i) file_kind, file_np
    if(i.ne.0) then
      ierr = 3
    else if (file_np .ne. m%np) then
      ierr = 4
    end if

    if(ierr==0) then
      if(file_kind == function_kind) then
        read(unit = iunit) f(1:m%np)
        
      else ! Adequate conversions....
        select case(file_kind)
        case(doutput_kind*4) ! Real, single precision
          allocate(rs(m%np))
          read(unit = iunit) rs(1:m%np)
          f = rs
          deallocate(rs)
          ierr = -1
          
        case(zoutput_kind*4) ! Complex, single precision
          allocate(cs(m%np))
          read(unit = iunit) cs(1:m%np)
          f = cs
          deallocate(cs)
          ierr = -2
          
        case(doutput_kind*8) ! Real, double precision
          allocate(rd(m%np))
          read(unit = iunit) rd(1:m%np)
          f = rd
          deallocate(rd)
          ierr = -3
          
        case(zoutput_kind*8) ! Complex, double precision
          allocate(cd(m%np))
          read(unit = iunit) cd(1:m%np)
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

    file = io_workpath(filename)

    status = nf90_open(trim(file), NF90_WRITE, ncid)
    if(status.ne.NF90_NOERR) then
      ierr = 2
      return
    end if

    !Inquire about dimensions
    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_1", dim_data_id(1))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_2", dim_data_id(2))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_inq_dimid (ncid, "dim_3", dim_data_id(3))
      call ncdf_error('nf90_inq_dimid', status, file, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(1), len = ndim(1))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    endif
    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(2), len = ndim(2))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    endif
    if(status == NF90_NOERR) then
      status = nf90_inquire_dimension (ncid, dim_data_id(3), len = ndim(3))
      call ncdf_error('nf90_inquire_dimension', status, file, ierr)
    endif
    if((ndim(1) .ne. c%n(1)) .or. &
       (ndim(2) .ne. c%n(2)) .or. &
       (ndim(3) .ne. c%n(3))) then
      ierr = 12; return
    endif

    if(status == NF90_NOERR) then
      status = nf90_inq_varid (ncid, "rdata", data_id)
      call ncdf_error('nf90_inq_varid', status, file, ierr)
    endif
    status = nf90_inq_varid(ncid, "idata", data_im_id)
    if(status == 0) then
       file_kind = -1
    else
       file_kind = 1
    endif
    status = 0
       
    if(status == NF90_NOERR) then
      status = nf90_inquire_variable (ncid, data_id, xtype = xtype)
      call ncdf_error('nf90_inquire_variable', status, file, ierr)
    endif

    if(xtype == NF90_FLOAT) then
      file_kind = file_kind*4
    else
      file_kind = file_kind*8
    endif
    if(file_kind .ne. function_kind) then
      select case(file_kind)
        case(4);  ierr = -1
        case(-4); ierr = -2
        case(8);  ierr = -3
        case(-8); ierr = -4
      end select
    endif

#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      status = nf90_get_var (ncid, data_id, re%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_get_var', status, file, ierr)
    endif
    if(file_kind<0) then
      if(status == NF90_NOERR) then
        status = nf90_get_var (ncid, data_im_id, im%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
        call ncdf_error('nf90_get_var', status, file, ierr)
      endif
    endif
#else
    if(status == NF90_NOERR) then
      status = nf90_get_var (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_get_var', status, file, ierr)
    endif
#endif

    status = nf90_close(ncid)
  end subroutine dx_cdf

#endif

end subroutine X(input_function)


! ---------------------------------------------------------
subroutine X(output_function) (how, dir, fname, m, f, u, ierr)
  integer,          intent(in)  :: how
  character(len=*), intent(in)  :: dir, fname
  type(mesh_type),  intent(in)  :: m
  R_TYPE,           intent(in)  :: f(:)  ! f(m%np)
  FLOAT,            intent(in)  :: u
  integer,          intent(out) :: ierr
  
  integer :: i
  FLOAT   :: x0
  character(len=20)  :: mformat, mformat2, mfmtheader
  logical            :: gnuplot_mode = .false.

  call io_mkdir(dir)

! Define the format; check if code is single precision or double precision
#if defined(SINGLE_PRECISION)
    mformat    = '(4es15.6)'
    mformat2   = '(i4,5es15.6)'
    mfmtheader = '(a,a7,5a15)'
#else
    mformat    = '(4es23.14)'
    mformat2   = '(i4,5es23.14)'
    mfmtheader = '(a,a10,5a23)'
#endif

  if(iand(how, output_gnuplot)   .ne.0) gnuplot_mode = .true.
  if(iand(how, output_plain)     .ne.0) call plain()
  if(iand(how, output_axis_x)    .ne.0) call axis_x()
  if(iand(how, output_axis_y)    .ne.0) call axis_y()
  if(iand(how, output_axis_z)    .ne.0) call axis_z()
  if(iand(how, output_plane_x)   .ne.0) call plane_x()
  if(iand(how, output_plane_y)   .ne.0) call plane_y()
  if(iand(how, output_plane_z)   .ne.0) call plane_z()
  if(iand(how, output_mesh_index).ne.0) call out_mesh_index()
  if(iand(how, output_dx)        .ne.0) call dx()
#if defined(HAVE_NETCDF)
  if(iand(how, output_dx_cdf)    .ne.0) call dx_cdf()
#endif

contains
  ! ---------------------------------------------------------
  subroutine plain()
    integer :: iunit

    iunit = io_open(trim(dir)//'/'//trim(fname), action='write', form='unformatted')

    write(unit=iunit) X(output_kind)*kind(f(1)), m%np
    write(unit=iunit) f(1:m%np)
    call io_close(iunit)

  end subroutine plain


  ! ---------------------------------------------------------
  subroutine axis_x()
    integer :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".y=0,z=0", action='write')

    write(iunit, mfmtheader) '#', 'x', 'Re', 'Im'
    do i = 1, m%np
      if(m%Lxyz(i, 2)==0.and.m%Lxyz(i, 3)==0) then
        write(iunit, mformat) m%x(i,1), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    call io_close(iunit)

  end subroutine axis_x


  ! ---------------------------------------------------------
  subroutine axis_y()
    integer :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".x=0,z=0", action='write')

    write(iunit, mfmtheader) '#', 'y', 'Re', 'Im'
    do i = 1, m%np
      if(m%Lxyz(i, 1)==0.and.m%Lxyz(i, 3)==0) then
        write(iunit, mformat) m%x(i,2), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    call io_close(iunit)

  end subroutine axis_y


  ! ---------------------------------------------------------
  subroutine axis_z()
    integer :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".x=0,y=0", action='write')

    write(iunit, mfmtheader) '#', 'z', 'Re', 'Im'
    do i = 1, m%np
      if(m%Lxyz(i, 1)==0.and.m%Lxyz(i, 2)==0) then
        write(iunit, mformat) m%x(i,3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    call io_close(iunit)

  end subroutine axis_z


  ! ---------------------------------------------------------
  subroutine plane_x()
    integer  :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".x=0", action='write')

    write(iunit, mfmtheader) '#', 'x', 'y', 'Re', 'Im'
    x0 = m%x(1,2)
    if(gnuplot_mode) write(iunit, mformat)
    
    do i = 1, m%np
      if(gnuplot_mode.and.x0 /= m%x(i, 2)) then
        write(iunit, mformat)      ! write extra lines for gnuplot grid mode
        x0 = m%x(i, 2)
      end if
      if(m%Lxyz(i,1)==0) then
        write(iunit, mformat) m%x(i,2), m%x(i,3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    
    if(gnuplot_mode) write(iunit, mformat)
    call io_close(iunit)
    
  end subroutine plane_x


  ! ---------------------------------------------------------
   subroutine plane_y()
    integer  :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".y=0", action='write')

    write(iunit, MFMTHEADER) '#', 'x', 'z', 'Re', 'Im'
    x0 = m%x(1,1)
    if(gnuplot_mode) write(iunit, mformat)
    
    do i = 1, m%np
      if(gnuplot_mode.and.x0 /= m%x(i, 1)) then
        write(iunit, mformat)      ! write extra lines for gnuplot grid mode
        x0 = m%x(i, 1)
      end if
      if(m%Lxyz(i,2)==0) then
        write(iunit, mformat) m%x(i,1), m%x(i,3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    
    if(gnuplot_mode) write(iunit, mformat)
    call io_close(iunit)
    
  end subroutine plane_y


  ! ---------------------------------------------------------
  subroutine plane_z()
    integer  :: iunit, i

    iunit = io_open(trim(dir)//'/'//trim(fname)//".z=0", action='write')

    write(iunit, MFMTHEADER) '#', 'x', 'y', 'Re', 'Im'
    x0 = m%x(1,1)
    if(gnuplot_mode) write(iunit, mformat)
    
    do i = 1, m%np
      if(gnuplot_mode.and.x0 /= m%x(i, 1)) then
        write(iunit, mformat)      ! write extra lines for gnuplot grid mode
        x0 = m%x(i, 1)
      end if
       
      if(m%Lxyz(i,3) == 0) then
        write(iunit, mformat) m%x(i,1), m%x(i,2), R_REAL(f(i))/u, R_AIMAG(f(i))/u
      end if
    end do
    
    if(gnuplot_mode) write(iunit, mformat)
    call io_close(iunit)
    
  end subroutine plane_z


  ! ---------------------------------------------------------
  subroutine out_mesh_index()
    integer  :: iunit, i
    
    iunit = io_open(trim(dir)//'/'//trim(fname)//".mesh_index", action='write')

    write(iunit, mfmtheader, iostat=ierr) '#', 'Index', 'x', 'y', 'z', 'Re', 'Im'
    x0 = m%x(1,1)
    if(ierr == 0.and.gnuplot_mode) write(iunit, mformat, iostat=ierr)
       
    do i=1,m%np
       if (ierr == 0.and.gnuplot_mode.and.x0 /= m%x(i, 1)) then
          write(iunit, mformat, iostat=ierr)      ! write extra lines for gnuplot grid mode
          x0 = m%x(i, 1)
       endif
       if(ierr==0) write(iunit, mformat2, iostat=ierr) i, m%x(i,1),  &
            m%x(i,2), m%x(i,3), R_REAL(f(i))/u, R_AIMAG(f(i))/u
    enddo

    if(ierr == 0.and.gnuplot_mode) write(iunit, mformat, iostat=ierr)
    call io_close(iunit)
  end subroutine out_mesh_index


  ! ---------------------------------------------------------
  subroutine dx()
    integer :: iunit, ix, iy, iz
    FLOAT   :: offset(3)
    character(len=40) :: nitems
    type(X(cf)) :: c

    ! put values in a nice cube
    call X(cf_new) (m%l, c)
    call X(cf_alloc_RS) (c)
    call X(mf2cf) (m, f, c)
    
    ! the offset is different in periodic directions
    do i = 1,conf%periodic_dim
      offset(i)=-(c%n(i))/2 * m%h(i) / units_out%length%factor
    end do
    do i = conf%periodic_dim+1, 3
      offset(i)=-(c%n(i) - 1)/2 * m%h(i) / units_out%length%factor
    end do

! just for nice formatting of the output
    write(nitems,*)c%n(1)*c%n(2)*c%n(3)
    nitems=trim(adjustl(nitems))

    iunit = io_open(trim(dir)//'/'//trim(fname)//".dx", action='write')
    
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
    character(len=512) :: filename
    integer :: ncid, status, data_id, data_im_id, pos_id, dim_data_id(3), dim_pos_id(2)
    real(r4) :: pos(2, 3)
    type(X(cf)) :: c

    ierr = 0

    ! put values in a nice cube
    call X(cf_new) (m%l, c)
    call X(cf_alloc_RS) (c)
    call X(mf2cf) (m, f, c)

    filename = io_workpath(trim(dir)//'/'//trim(fname)//".ncdf");

    status = nf90_create(trim(filename), NF90_CLOBBER, ncid)
    if(status.ne.NF90_NOERR) then
      ierr = 2
      return
    end if

    ! dimensions
    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_1", c%n(1), dim_data_id(1))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_2", c%n(2), dim_data_id(2))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "dim_3", c%n(3), dim_data_id(3))
      call ncdf_error('nf90_der_dim', status, filename, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "pos_1", 2, dim_pos_id(1))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    endif

    if(status == NF90_NOERR) then
      status = nf90_def_dim (ncid, "pos_2", 3, dim_pos_id(2))
      call ncdf_error('nf90_def_dim', status, filename, ierr)
    endif

#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_FLOAT, dim_data_id, data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    endif
#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "rdata", NF90_DOUBLE, dim_data_id, data_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    endif
#endif
#if defined(R_TCOMPLEX)
#if defined(SINGLE_PRECISION)
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_FLOAT, dim_data_id, data_im_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    endif

#else
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "idata", NF90_DOUBLE, dim_data_id, data_im_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    endif
#endif
#endif
    if(status == NF90_NOERR) then
      status = nf90_def_var (ncid, "pos", NF90_FLOAT,  dim_pos_id,  pos_id)
      call ncdf_error('nf90_def_var', status, filename, ierr)
    endif

    ! attributes
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_id, "field", "rdata, scalar")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    endif
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_id, "positions", "pos, regular")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    endif
#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_im_id, "field", "idata, scalar")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    endif
    if(status == NF90_NOERR) then
      status = nf90_put_att (ncid, data_im_id, "positions", "pos, regular")
      call ncdf_error('nf90_put_att', status, filename, ierr)
    endif
#endif

    ! end definitions
    status = nf90_enddef (ncid)

    ! data
    pos(1,:) = real(-(c%n(:) - 1)/2 * m%h(:) / units_out%length%factor, 4)
    pos(2,:) = real(m%h(:) / units_out%length%factor, 4)

    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, pos_id, pos(:,:))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    endif

#if defined(R_TCOMPLEX)
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_id, real(c%RS, PRECISION), map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    endif
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_im_id, aimag(c%RS), map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    endif
#else
    if(status == NF90_NOERR) then
      status = nf90_put_var (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /))
      call ncdf_error('nf90_put_var', status, filename, ierr)
    endif
#endif

    ! close
    status = nf90_close(ncid)
    call X(cf_free) (c)

  end subroutine dx_cdf

#endif

end subroutine X(output_function)
