subroutine R_FUNC(output_function) (outp, dir, fname, m, f, u)
  type(output_type), intent(IN) :: outp
  character(len=*), intent(IN) :: dir, fname
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(IN) :: f(m%np)
  real(r8), intent(in) :: u
  
  integer :: iunit, i
  R_TYPE, allocatable :: c(:,:,:)

  ! do not bother with errors
  call oct_mkdir(C_string(trim(dir)))

  allocate(c(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
  c = R_TOTYPE(0._r8)
  call R_FUNC(mesh_to_cube) (m, f, c)

  if(iand(outp%how, output_axis_x).ne.0) call axis_x()
  if(conf%dim > 1.and.iand(outp%how, output_axis_y).ne.0) call axis_y()
  if(conf%dim > 2.and.iand(outp%how, output_axis_z).ne.0) call axis_z()
  if(conf%dim > 2.and.iand(outp%how, output_plane_x).ne.0) call plane_x()
  if(conf%dim > 2.and.iand(outp%how, output_plane_y).ne.0) call plane_y()
  if(conf%dim > 1.and.iand(outp%how, output_plane_z).ne.0) call plane_z()
  if(iand(outp%how, output_dx).ne.0) call dx()

  deallocate(c)

contains

  subroutine axis_x()
    integer  :: ix
    real(r8) :: x

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".y=0,z=0", status='unknown')
    do ix = 1, m%fft_n(1)
       x = (ix - m%fft_n(1)/2 - 1)*m%h(1)/units_out%length%factor
       write(iunit, *) x, R_REAL(c(ix, m%fft_n(2)/2 + 1, m%fft_n(3)/2 + 1))/u, &
            R_AIMAG(c(ix, m%fft_n(2)/2 + 1, m%fft_n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_x

  subroutine axis_y()
    integer  :: iy
    real(r8) :: y

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0,z=0", status='unknown')
    do iy = 1, m%fft_n(2)
       y = (iy - m%fft_n(2)/2 - 1)*m%h(2)/units_out%length%factor
       write(iunit, *) y, R_REAL(c(m%fft_n(1)/2 + 1, iy, m%fft_n(3)/2 + 1))/u, &
            R_AIMAG(c(m%fft_n(1)/2 + 1, iy, m%fft_n(3)/2 + 1))/u
    enddo
    call io_close(iunit)
  end subroutine axis_y

  subroutine axis_z()
    integer  :: iz
    real(r8) :: z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0,y=0", status='unknown')
    do iz = 1, m%fft_n(3)
       z = (iz - m%fft_n(3)/2 - 1)*m%h(3)/units_out%length%factor
       write(iunit, *) z, R_REAL(c(m%fft_n(1)/2 + 1, m%fft_n(2)/2 + 1, iz))/u, &
            R_AIMAG(c(m%fft_n(1)/2 + 1, m%fft_n(2)/2 + 1, iz))/u
    enddo
    call io_close(iunit)
  end subroutine axis_z

  subroutine plane_x()
    integer  :: iy, iz
    real(r8) :: y, z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0", status='unknown')
    do iy = 1, m%fft_n(2)
      y = (iy - m%fft_n(2)/2 - 1)*m%h(2)/units_out%length%factor
      do iz = 1, m%fft_n(3)
        z = (iz - m%fft_n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, *) y, z, R_REAL(c(m%fft_n(1)/2 + 1, iy, iz))/u, &
             R_AIMAG(c(m%fft_n(1)/2 + 1, iy, iz))/u
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
    do ix = 1, m%fft_n(1)
      x = (ix - m%fft_n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iz = 1, m%fft_n(3)
        z = (iz - m%fft_n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, *) x, z, R_REAL(c(ix, m%fft_n(2)/2 + 1, iz))/u, &
             R_AIMAG(c(ix, m%fft_n(2)/2 + 1, iz))/u
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
    do ix = 1, m%fft_n(1)
      x = (ix - m%fft_n(1)/2 - 1)*m%h(1)/units_out%length%factor
      do iy = 1, m%fft_n(2)
        y = (iy - m%fft_n(2)/2 - 1)*m%h(2)/units_out%length%factor
        write(iunit, *) x, y, R_REAL(c(ix, iy, m%fft_n(3)/2 + 1))/u, &
             R_AIMAG(c(ix, iy, m%fft_n(3)/2 + 1))/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_z

  subroutine dx()
    integer :: ix, iy, iz

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".dx", status='unknown')
    
    write(iunit, '(a)') 'header = marker "Start data\n"'
    write(iunit,'(a,i4,a,i4,a,i4)') 'grid = ', m%fft_n(1), ' x ', m%fft_n(2), ' x ', m%fft_n(3)
    write(iunit,'(a)') 'format = ascii'
    write(iunit,'(a)') 'interleaving = record'
    write(iunit,'(a)') 'majority = row'
    write(iunit,'(a)') 'field = field0'
    write(iunit,'(a)') 'type = float'
    write(iunit,'(a)') 'dependency = positions'
    
    write(iunit,'(a,6(a,f12.6))') 'positions = regular, regular, regular', &
         ',', -(m%fft_n(1) - 1)/2 * m%h(1) / units_out%length%factor, &
         ',', m%h(1) / units_out%length%factor, &
         ',', -(m%fft_n(2) - 1)/2 * m%h(2) / units_out%length%factor, &
         ',', m%h(2) / units_out%length%factor, &
         ',', -(m%fft_n(3) - 1)/2 * m%h(3) / units_out%length%factor, &
         ',', m%h(3) / units_out%length%factor

    write(iunit, '(1x)')
    write(iunit, '(a)') 'end'
    write(iunit, '(a)') 'Start data'
    do ix = 1, m%fft_n(1)
      do iy = 1, m%fft_n(2)
        do iz = 1, m%fft_n(3)
          write(iunit, '(e17.10)') c(ix, iy, iz)
        end do
      end do
    end do
    
    call io_close(iunit)
  end subroutine dx

end subroutine R_FUNC(output_function)

