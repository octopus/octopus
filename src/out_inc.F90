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

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reads a function from file filename, and puts it into f. The input file
! may be a "plain" file (no extension), or a netcdf file ".ncdf" extension.
! On output, ierr signals how everything went:
! ierr > 0 => Error. The function f was not read:
!             1 : illegal filename (must have no extension, or ".ncdf" extension.
!             2 : file could not be succesfully opened.
!             3 : file opened, but error reading.
!             4 : The number of points/mesh dimensions do not coincide. 
!             5 : NetCDF error (one or several warnings are emitted)
! ierr = 0 => Success.
! ierr < 0 => Success, but some kind of type conversion was necessary. The
!             of ierr is then:
!             -1 : function in file is real, sp.
!             -2 : function in file is complex, sp.
!             -3 : function in file is real, dp.
!             -4 : function in file is complex, dp.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/!
subroutine X(input_function)(filename, m, f, ierr)
  character(len=*), intent(in)  :: filename
  type(mesh_type),  intent(in)  :: m
  R_TYPE,           intent(out) :: f(:)
  integer, intent(out)          :: ierr

  integer :: iunit, i, function_kind, file_kind
  type(X(cf)) :: c
#if defined(R_TCOMPLEX)
  type(dcf) :: re, im
#endif

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

  subroutine plain
    integer :: file_np
    real(4),    allocatable :: rs(:)
    real(8),    allocatable :: rd(:)
    complex(4), allocatable :: cs(:)
    complex(8), allocatable :: cd(:)
    call io_assign(iunit)
    open(unit = iunit, file = trim(filename), &
         status = 'old', action = 'read', form = 'unformatted', iostat = i)
      if(i.ne.0) then; ierr = 2; return; endif
    read(unit = iunit, iostat = i) file_kind, file_np
      if(i.ne.0) then; ierr = 3; return; endif
      if(file_np .ne. m%np) then; ierr = 4; return; endif
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
        case default
      end select
    endif
    call io_close(iunit)
  end subroutine plain

#if defined(HAVE_NETCDF)
! I define this macro in order to call ncdf_error every time.
#define NCDFCALL(x, y) status = x y; call ncdf_error(#x, status, filename, ierr)
  subroutine dx_cdf()
    integer :: ncid, ndims, nvars, natts, status, data_id, data_im_id, pos_id, &
               dim_data_id(3), dim_pos_id(2), ndim(3), xtype
    real(r4) :: pos(2, 3)
    logical :: function_is_complex = .false.

    NCDFCALL(nf90_open, (trim(filename), NF90_WRITE, ncid))

    !Inquire about dimensions
    NCDFCALL(nf90_inq_dimid, (ncid, "dim_1", dim_data_id(1)))
    NCDFCALL(nf90_inq_dimid, (ncid, "dim_2", dim_data_id(2)))
    NCDFCALL(nf90_inq_dimid, (ncid, "dim_3", dim_data_id(3)))
    NCDFCALL(nf90_inquire_dimension, (ncid, dim_data_id(1), len = ndim(1)))
    NCDFCALL(nf90_inquire_dimension, (ncid, dim_data_id(2), len = ndim(2)))
    NCDFCALL(nf90_inquire_dimension, (ncid, dim_data_id(3), len = ndim(3)))
    if((ndim(1) .ne. c%n(1)) .or. &
       (ndim(2) .ne. c%n(2)) .or. &
       (ndim(3) .ne. c%n(3))) then
      ierr = 12; return
    endif

    NCDFCALL(nf90_inq_varid, (ncid, "rdata", data_id))
    status = nf90_inq_varid(ncid, "idata", data_im_id)
    if(status == 0) then
       file_kind = -1
    else
       file_kind = 1
    endif
       
    NCDFCALL(nf90_inquire_variable, (ncid, data_id, xtype = xtype))

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
    NCDFCALL(nf90_get_var, (ncid, data_id, re%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
    if(file_kind<0) then
       NCDFCALL(nf90_get_var, (ncid, data_im_id, im%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
    endif
#else
    NCDFCALL(nf90_get_var, (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
#endif

    NCDFCALL(nf90_close, (ncid))
  end subroutine dx_cdf

#undef NCDFCALL
#endif

end subroutine X(input_function)

subroutine X(output_function) (how, dir, fname, m, f, u)
  integer,          intent(in) :: how
  character(len=*), intent(in) :: dir, fname
  type(mesh_type),  intent(IN) :: m
  R_TYPE,           intent(IN) :: f(:)  ! f(m%np)
  FLOAT,            intent(in) :: u
  
  integer :: iunit, i
  type(X(cf)) :: c
  character(len=20) :: mformat

  ! do not bother with errors
  call loct_mkdir(trim(dir))

  if(iand(how, output_plain) .ne.0) call plain()
  if(how==output_plain) return ! Return if no other output is needed.

  call X(cf_new) (m%l, c)
  call X(cf_alloc_RS) (c)
  call X(mf2cf) (m, f, c)

! Define the format; check if code is single precision or double precision
#if defined(SINGLE_PRECISION)
    mformat = '(4es15.6)'
#else
    mformat = '(4es23.14)'
#endif

  if(iand(how, output_axis_x) .ne.0) call axis_x()
  if(iand(how, output_axis_y) .ne.0) call axis_y()
  if(iand(how, output_axis_z) .ne.0) call axis_z()
  if(iand(how, output_plane_x).ne.0) call plane_x()
  if(iand(how, output_plane_y).ne.0) call plane_y()
  if(iand(how, output_plane_z).ne.0) call plane_z()
  if(iand(how, output_dx)     .ne.0) call dx()
#if defined(HAVE_NETCDF)
  if(iand(how, output_dx_cdf) .ne.0) call dx_cdf()
#endif

  call X(cf_free) (c)

contains
#define MFMTHEADER '(a,a3,a12,a12,a12)'

  subroutine plain()
    call io_assign(iunit)
    open(unit = iunit, file = trim(dir) // "/" // trim(fname), status = 'unknown', form = 'unformatted')
    write(unit = iunit) X(output_kind)*kind(f(1)), m%np
    write(unit = iunit) f(1:m%np)
    call io_close(iunit)
  end subroutine plain

  subroutine axis_x()
    integer  :: ix
    FLOAT :: x

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".y=0,z=0", status='unknown')
    do ix = 1, c%n(1)
       x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
       write(iunit, mformat) x, R_REAL(c%RS(ix, c%n(2)/2 + 1, c%n(3)/2 + 1))/u, &
            R_AIMAG(c%RS(ix, c%n(2)/2 + 1, c%n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_x

  subroutine axis_y()
    integer  :: iy
    FLOAT :: y

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".x=0,z=0", status='unknown')
    do iy = 1, c%n(2)
       y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
       write(iunit, mformat) y, R_REAL(c%RS(c%n(1)/2 + 1, iy, c%n(3)/2 + 1))/u, &
            R_AIMAG(c%RS(c%n(1)/2 + 1, iy, c%n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_y

  subroutine axis_z()
    integer  :: iz
    FLOAT :: z

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".x=0,y=0", status='unknown')
    do iz = 1, c%n(3)
       z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
       write(iunit, mformat) z, R_REAL(c%RS(c%n(1)/2 + 1, c%n(2)/2 + 1, iz))/u, &
            R_AIMAG(c%RS(c%n(1)/2 + 1, c%n(2)/2 + 1, iz))/u
    enddo
    call io_close(iunit)
  end subroutine axis_z

  subroutine plane_x()
    integer  :: iy, iz
    FLOAT :: y, z

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".x=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'y', 'z', 'Re', 'Im'
    do iy = 1, c%n(2)
      y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
      do iz = 1, c%n(3)
        z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, mformat) y, z, R_REAL(c%RS(c%n(1)/2 + 1, iy, iz))/u, &
             R_AIMAG(c%RS(c%n(1)/2 + 1, iy, iz))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_x

  subroutine plane_y()
    integer  :: ix, iz
    FLOAT :: x, z

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".y=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'x', 'z', 'Re', 'Im'
    do ix = 1, c%n(1)
      x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iz = 1, c%n(3)
        z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, mformat) x, z, R_REAL(c%RS(ix, c%n(2)/2 + 1, iz))/u, &
             R_AIMAG(c%RS(ix, c%n(2)/2 + 1, iz))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_y

  subroutine plane_z()
    integer  :: ix, iy
    FLOAT :: x, y

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".z=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'x', 'y', 'Re', 'Im'
    do ix = 1, c%n(1)
      x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iy = 1, c%n(2)
        y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
        write(iunit, mformat) x, y, R_REAL(c%RS(ix, iy, c%n(3)/2 + 1))/u, &
             R_AIMAG(c%RS(ix, iy, c%n(3)/2 + 1))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_z

  subroutine dx()
    integer :: ix, iy, iz
    FLOAT :: offset(3)
    character(LEN=40) :: nitems
    
! the offset is different in periodic directions
    do i=1,conf%periodic_dim
      offset(i)=-(c%n(i))/2 * m%h(i) / units_out%length%factor
    end do
    do i=conf%periodic_dim+1,3
      offset(i)=-(c%n(i) - 1)/2 * m%h(i) / units_out%length%factor
    end do

! just for nice formatting of the output
    write(nitems,*)c%n(1)*c%n(2)*c%n(3)
    nitems=TRIM(ADJUSTL(nitems))

    call io_assign(iunit)
    open(iunit, file=trim(dir) // "/" // trim(fname) // ".dx", status='unknown')
    
    write(iunit, '(a,3i7)') 'object 1 class gridpositions counts',c%n(:)
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
  end subroutine dx

#if defined(HAVE_NETCDF)
! I define this macro in order to call ncdf_error every time.
#define NCDFCALL(x, y) status = x y; call ncdf_error(#x, status, filename, ierr)
  subroutine dx_cdf()
    character(len=200) :: filename
    integer :: ncid, status, data_id, data_im_id, pos_id, dim_data_id(3), dim_pos_id(2), ierr
    real(r4) :: pos(2, 3)

    filename = trim(dir) // "/" // trim(fname) // ".ncdf"
    NCDFCALL(nf90_create, (trim(filename), NF90_CLOBBER, ncid))

    ! dimensions
    NCDFCALL(nf90_def_dim, (ncid, "dim_1", c%n(1), dim_data_id(1)))
    NCDFCALL(nf90_def_dim, (ncid, "dim_2", c%n(2), dim_data_id(2)))
    NCDFCALL(nf90_def_dim, (ncid, "dim_3", c%n(3), dim_data_id(3)))

    NCDFCALL(nf90_def_dim, (ncid, "pos_1", 2, dim_pos_id(1)))
    NCDFCALL(nf90_def_dim, (ncid, "pos_2", 3, dim_pos_id(2)))

#if defined(SINGLE_PRECISION)
    NCDFCALL(nf90_def_var, (ncid, "rdata", NF90_FLOAT, dim_data_id, data_id))
#else
    NCDFCALL(nf90_def_var, (ncid, "rdata", NF90_DOUBLE, dim_data_id, data_id))
#endif
#if defined(R_TCOMPLEX)
#if defined(SINGLE_PRECISION)
    NCDFCALL(nf90_def_var, (ncid, "idata", NF90_FLOAT, dim_data_id, data_im_id))
#else
    NCDFCALL(nf90_def_var, (ncid, "idata", NF90_DOUBLE, dim_data_id, data_im_id))
#endif
#endif
    NCDFCALL(nf90_def_var, (ncid, "pos", NF90_FLOAT,  dim_pos_id,  pos_id))

    ! attributes
    NCDFCALL(nf90_put_att, (ncid, data_id, "field", "rdata, scalar"))
    NCDFCALL(nf90_put_att, (ncid, data_id, "positions", "pos, regular"))
#if defined(R_TCOMPLEX)
    NCDFCALL(nf90_put_att, (ncid, data_im_id, "field", "idata, scalar"))
    NCDFCALL(nf90_put_att, (ncid, data_im_id, "positions", "pos, regular"))
#endif

    ! end definitions
    NCDFCALL(nf90_enddef, (ncid))

    ! data
    pos(1,:) = real(-(c%n(:) - 1)/2 * m%h(:) / units_out%length%factor, 4)
    pos(2,:) = real(m%h(:) / units_out%length%factor, 4)

    NCDFCALL(nf90_put_var, (ncid, pos_id, pos(:,:)))

#if defined(R_TCOMPLEX)
    NCDFCALL(nf90_put_var, (ncid, data_id, real(c%RS, PRECISION), map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
    NCDFCALL(nf90_put_var, (ncid, data_im_id, aimag(c%RS), map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
#else
    NCDFCALL(nf90_put_var, (ncid, data_id, c%RS, map = (/c%n(3)*c%n(2), c%n(2), 1 /)))
#endif

    ! close
    NCDFCALL(nf90_close, (ncid))
  end subroutine dx_cdf

#undef NCDFCALL
#endif

end subroutine X(output_function)
