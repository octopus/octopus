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

#if defined(ONE_D)
  call io_assign(iunit)
  open(iunit, file=trim(dir)+"/"+trim(fname)+".dat", status='unknown')
  do i = 1, m%np
    write(iunit, *) m%lx(i)/units_out%length%factor, f(i)/u
  end do
  call io_close(i)

#elif defined(THREE_D)

  allocate(c(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
  c = R_TOTYPE(0._r8)
  call R_FUNC(mesh_to_cube) (m, f, c)

  if(iand(outp%how, output_plane_x).ne.0) call plane_x()
  if(iand(outp%how, output_plane_y).ne.0) call plane_y()
  if(iand(outp%how, output_plane_z).ne.0) call plane_z()
  if(iand(outp%how, output_dx     ).ne.0) call dx()

  deallocate(c)
#endif

contains
  subroutine plane_x()
    integer  :: iy, iz
    real(r8) :: y, z

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".x=0", status='unknown')
    do iy = 1, m%fft_n(2)
      y = (iy - m%fft_n(2)/2 - 1)*m%h(2)/units_out%length%factor
      do iz = 1, m%fft_n(3)
        z = (iz - m%fft_n(3)/2 - 1)*m%h(3)/units_out%length%factor
        write(iunit, *) y, z, c(m%fft_n(1)/2 + 1, iy, iz)/u
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
        write(iunit, *) x, z, c(ix, m%fft_n(2)/2 + 1, iz)/u
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
        write(iunit, *) x, y, c(ix, iy, m%fft_n(3)/2 + 1)/u
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)
  end subroutine plane_z

  subroutine dx()
    character(len=12) :: posxc, posyc, poszc, h1c, h2c, h3c
    character(len=200) :: pwd, line, datafile
    integer :: ix, iy, iz
    real(r8) :: pos

    ! first write the block data
    call io_assign(iunit)
    datafile = trim(dir) + "/" + trim(fname) + ".dxdata"
    open(iunit, file=trim(datafile), status='unknown')
    do ix = 1, m%fft_n(1)
      do iy = 1, m%fft_n(2)
        do iz = 1, m%fft_n(3)
          write(iunit, '(e17.10)') c(ix, iy, iz)
        end do
      end do
    end do
    call io_close(iunit)
    
    ! get working directory
    call clear_str(pwd); call oct_getcwd(pwd)
    open(iunit, file=trim(dir)+"/"+trim(fname)+".general", status='unknown')
    
    write(line,'(a)')              'file = '+trim(pwd(1:len_trim(pwd)-1))+'/'+trim(datafile)
    write(iunit,'(a)') trim(line)
    write(line,'(a,i4,a,i4,a,i4)') 'grid = ', m%fft_n(1), ' x ', m%fft_n(2), ' x ', m%fft_n(3);
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'format = ascii'
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'interleaving = record'
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'majority = row'
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'field = field0'
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'type = float'
    write(iunit,'(a)') trim(line)
    write(line,'(a)')              'dependency = positions'
    write(iunit,'(a)') trim(line)
    
    pos = -(m%fft_n(1) - 1)/2 * m%h(1); write(posxc,'(f12.6)') pos / units_out%length%factor
    pos = -(m%fft_n(1) - 1)/2 * m%h(2); write(posyc,'(f12.6)') pos / units_out%length%factor
    pos = -(m%fft_n(1) - 1)/2 * m%h(3); write(poszc,'(f12.6)') pos / units_out%length%factor
    write(h1c,'(f12.6)') m%h(1) / units_out%length%factor
    write(h2c,'(f12.6)') m%h(2) / units_out%length%factor
    write(h3c,'(f12.6)') m%h(3) / units_out%length%factor
    write(line,'(a)') &
         'positions = regular, regular, regular, '+ &
         posxc + ', ' + &
         h1c   + ', ' + &
         posyc + ', ' + &
         h2c   + ', ' + &
         poszc + ', ' + &
         h3c
    write(iunit,'(a)') trim(line)
    
    call io_close(iunit)
  end subroutine dx

end subroutine R_FUNC(output_function)

