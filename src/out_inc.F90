subroutine R_FUNC(output_function) (outp, dir, fname, m, f, u)
  type(output_type), intent(IN) :: outp
  character(len=*), intent(IN) :: dir, fname
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)
  real(r8), intent(in) :: u
  
  integer :: iunit, i
  type(X(cf)) :: c

  ! do not bother with errors
  call oct_mkdir(trim(dir))

  call X(cf_new) (m%l, c)
  call X(cf_alloc_RS) (c)
  call R_FUNC(mf2cf) (m, f, c)

  if(iand(outp%how, output_axis_x) .ne.0) call axis_x()
  if(iand(outp%how, output_axis_y) .ne.0) call axis_y()
  if(iand(outp%how, output_axis_z) .ne.0) call axis_z()
  if(iand(outp%how, output_plane_x).ne.0) call plane_x()
  if(iand(outp%how, output_plane_y).ne.0) call plane_y()
  if(iand(outp%how, output_plane_z).ne.0) call plane_z()
  if(iand(outp%how, output_dx)     .ne.0) call dx()
#if defined(HAVE_NETCDF) && defined(R_TREAL)
  if(iand(outp%how, output_dx_cdf) .ne.0) call dx_cdf()
#endif

  call X(cf_free) (c)

contains
#define MFORMAT    '(4e20.8)'
#define MFMTHEADER '(a,a3,a12,a12,a12)'

  subroutine axis_x()
    integer  :: ix
    real(r8) :: x

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".y=0,z=0", status='unknown')
    do ix = 1, c%n(1)
       x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
       write(iunit, MFORMAT) x, R_REAL(c%RS(ix, c%n(2)/2 + 1, c%n(3)/2 + 1))/u, &
            R_AIMAG(c%RS(ix, c%n(2)/2 + 1, c%n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_x

  subroutine axis_y()
    integer  :: iy
    real(r8) :: y

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0,z=0", status='unknown')
    do iy = 1, c%n(2)
       y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
       write(iunit, MFORMAT) y, R_REAL(c%RS(c%n(1)/2 + 1, iy, c%n(3)/2 + 1))/u, &
            R_AIMAG(c%RS(c%n(1)/2 + 1, iy, c%n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_y

  subroutine axis_z()
    integer  :: iz
    real(r8) :: z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0,y=0", status='unknown')
    do iz = 1, c%n(3)
       z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
       write(iunit, MFORMAT) z, R_REAL(c%RS(c%n(1)/2 + 1, c%n(2)/2 + 1, iz))/u, &
            R_AIMAG(c%RS(c%n(1)/2 + 1, c%n(2)/2 + 1, iz))/u
    enddo
    call io_close(iunit)
  end subroutine axis_z

  subroutine plane_x()
    integer  :: iy, iz
    real(r8) :: y, z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'y', 'z', 'Re', 'Im'
    do iy = 1, c%n(2)
      y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
      do iz = 1, c%n(3)
        z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, MFORMAT) y, z, R_REAL(c%RS(c%n(1)/2 + 1, iy, iz))/u, &
             R_AIMAG(c%RS(c%n(1)/2 + 1, iy, iz))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_x

  subroutine plane_y()
    integer  :: ix, iz
    real(r8) :: x, z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".y=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'x', 'z', 'Re', 'Im'
    do ix = 1, c%n(1)
      x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iz = 1, c%n(3)
        z = (iz - c%n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, MFORMAT) x, z, R_REAL(c%RS(ix, c%n(2)/2 + 1, iz))/u, &
             R_AIMAG(c%RS(ix, c%n(2)/2 + 1, iz))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_y

  subroutine plane_z()
    integer  :: ix, iy
    real(r8) :: x, y

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".z=0", status='unknown')
    write(iunit, MFMTHEADER) '#', 'x', 'y', 'Re', 'Im'
    do ix = 1, c%n(1)
      x = (ix - c%n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iy = 1, c%n(2)
        y = (iy - c%n(2)/2 - 1)*m%h(2)/units_out%length%factor
        write(iunit, MFORMAT) x, y, R_REAL(c%RS(ix, iy, c%n(3)/2 + 1))/u, &
             R_AIMAG(c%RS(ix, iy, c%n(3)/2 + 1))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_z

  subroutine dx()
    integer :: ix, iy, iz
    real(r8) :: offset(3)
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
    open(iunit, file=trim(dir)+"/"+trim(fname)+".dx", status='unknown')
    
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

#if defined(HAVE_NETCDF) && defined(R_TREAL)
  subroutine dx_cdf()
    integer :: ncid, status, data_id, pos_id, dim_data_id(3), dim_pos_id(2)
    real(r4) :: pos(2, 3)

    status = nf90_create(trim(dir)+"/"+trim(fname)+".ncdf", NF90_CLOBBER, ncid)
    if (status /= NF90_NOERR) then
      call ncdf_error("nf90_create", status); return
    end if

    ! dimensions
    status = nf90_def_dim(ncid, "dim_1", c%n(1), dim_data_id(1))
    status = nf90_def_dim(ncid, "dim_2", c%n(2), dim_data_id(2))
    status = nf90_def_dim(ncid, "dim_3", c%n(3), dim_data_id(3))

    status = nf90_def_dim(ncid, "pos_1", 2, dim_pos_id(1))
    status = nf90_def_dim(ncid, "pos_2", 3, dim_pos_id(2))

    ! variables
    status = nf90_def_var(ncid, fname, NF90_DOUBLE, dim_data_id, data_id)
    status = nf90_def_var(ncid, "pos", NF90_FLOAT,  dim_pos_id,  pos_id)

    ! attributes
    status = nf90_put_att(ncid, data_id, "field", trim(fname)+", scalar")
    status = nf90_put_att(ncid, data_id, "positions", "pos, regular")

    ! end definitions
    status = nf90_enddef(ncid)

    ! data
    pos(1,:) = real(-(c%n(:) - 1)/2 * m%h(:) / units_out%length%factor, r4)
    pos(2,:) = real(m%h(:) / units_out%length%factor, r4)

    status = nf90_put_var(ncid, pos_id, pos(:,:))
    ! we have to transpose the matrix: stupid Fortran!
    status = nf90_put_var(ncid, data_id, c%RS, map=(/c%n(1)*c%n(2), c%n(1), 1/))

    ! close
    status = nf90_close(ncid)
    
  end subroutine dx_cdf

  subroutine ncdf_error(func, status)
    character(len=*), intent(in) :: func
    integer, intent(in) :: status

    message(1) = "NETCDF error in function'"+trim(func)+"'"
    write(message(2), '(6x,a,i4)')'error code = ', status
    call write_warning(2)

  end subroutine ncdf_error
#endif

end subroutine R_FUNC(output_function)

