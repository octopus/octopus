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

subroutine PES_mask_init(v, m, st)
  type(PES_mask_type), intent(out)   :: v
  type(mesh_type),     intent(inout) :: m
  type(states_type),   intent(IN)    :: st

  message(1) = 'Info: Calculating PES using mask technique'
  call write_info(1)
      
  ! alloc ffts in case they are not allocated yet
  call fft_init(m%l, fft_complex, v%fft)

  ! setup arrays to be used
  allocate(&
       v%k(m%l(1), m%l(2), m%l(3), st%dim, st%st_start:st%st_end, st%nik), &
       v%r(m%l(1), m%l(2), m%l(3),         st%st_start:st%st_end, st%nik))

  v%k = M_z0
  v%r = M_ZERO
   
end subroutine PES_mask_init

subroutine PES_mask_end(v)
  type(PES_mask_type), intent(inout) :: v

  if(associated(v%k)) then
    call fft_end(v%fft)
    deallocate(v%k, v%r)
    nullify   (v%k, v%r)
  end if

end subroutine PES_mask_end

subroutine PES_mask_doit(v, m, st, dt, mask)
  type(PES_mask_type), intent(inout) :: v
  type(mesh_type),     intent(IN)    :: m
  type(states_type),   intent(IN)    :: st
  FLOAT,               intent(in)    :: dt
  FLOAT,               pointer       :: mask(:)

  integer :: j, idim, ist, ik, ix, iy, iz, ix3(3), ixx(3)
  CMPLX, allocatable :: wf1(:,:,:), wf2(:,:,:)
  FLOAT :: temp(3), vec

  ! propagate wave-function in momentum space
  temp(:) = M_TWO*M_PI/(m%l(:)*m%h(:))
  do ix = 1, m%l(1)
    ixx(1) = pad_feq(ix, m%l(1), .true.)
    do iy = 1, m%l(2)
      ixx(2) = pad_feq(iy, m%l(2), .true.)
      do iz = 1, m%l(3)
        ixx(3) = pad_feq(iz, m%l(3), .true.)

        vec = sum((temp(:) * ixx(:))**2) / M_TWO
        v%k(ix, iy, iz,:,:,:) = v%k(ix, iy, iz,:,:,:)*exp(-M_zI*dt*vec)
      enddo
    enddo
  enddo
      
  ! we now add the contribution from this timestep
  allocate( &
       wf1(m%l(1), m%l(2), m%l(3)), &
       wf2(m%l(1), m%l(2), m%l(3)))
  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      do idim = 1, st%dim

        wf1 = M_z0

        do j = 1, m%np
          ix3(:) = m%Lxyz(:,j) + m%l(:)/2 + 1
          wf1(ix3(1), ix3(2), ix3(3)) = mask(j)*st%zpsi(j, idim, ist, ik)
        end do

        ! and add to our density (sum for idim, also)
        v%r(:,:,:, ist, ik) = v%r(:,:,:, ist, ik) + abs(wf1)**2

        ! now we FT
        call zfft_forward(v%fft, wf1, wf2)

        ! and add to our spectrum
        v%k(:,:,:, idim, ist, ik) = v%k(:,:,:, idim, ist, ik) + wf2

      end do
    end do
  end do

  deallocate(wf1, wf2)

  return
end subroutine PES_mask_doit

subroutine PES_mask_output(v, m, st, file)
  type(PES_mask_type), intent(IN) :: v
  type(mesh_type),     intent(IN) :: m
  type(states_type),   intent(IN) :: st
  character(len=*),    intent(in) :: file

  FLOAT, allocatable :: spis(:,:,:), arpis(:,:,:)
  FLOAT :: vec, temp(3)
  integer, allocatable :: npoints(:), ar_npoints(:)
  integer :: p, ik, ii, ix, iy, iz, ixx(3), iunit
  character(len=100) :: fn
  
  integer,  parameter :: n = 600, ar_n = 90
  FLOAT, parameter :: step = CNST(0.005)

  allocate( &
       spis (n,    st%st_start:st%st_end, st%nik), &
       arpis(ar_n, st%st_start:st%st_end, st%nik))
  allocate(npoints(n), ar_npoints(ar_n))
  spis = M_ZERO; arpis = M_ZERO
  npoints = 0;  ar_npoints = 0
  
  temp(:) = M_TWO*M_PI/(m%l(:)*m%h(:))
  do ix = 1, m%l(1)
    ixx(1) = pad_feq(ix, m%l(1), .true.)
    do iy = 1, m%l(2)
      ixx(2) = pad_feq(iy, m%l(2), .true.)
      do iz = 1, m%l(3)
        ixx(3) = pad_feq(iz, m%l(3), .true.)
        
        if(ixx(1).ne.0 .or. ixx(2).ne.0 .or. ixx(3).ne.0) then
          ! power spectrum
          vec = sum((temp(:) * ixx(:))**2) / M_TWO
          ii = nint(vec / step) + 1
          if(ii <= n) then
            do ik = 1, st%nik
              do p = st%st_start, st%st_end
                spis(ii, p, ik) = spis(ii, p, ik) + st%occ(p, ik) * &
                     sum(abs(v%k(ix, iy, iz, :, p, ik))**2)
                npoints(ii) = npoints(ii) + 1
              end do
            end do
          end if
          
          ! angle resolved (assumes the pol is in the x direction)
          if(ixx(3)==0 .and. (ixx(1).ne.0 .or. ixx(2).ne.0)) then
            vec = atan2(real(ixx(2), PRECISION), real(ixx(1), PRECISION))
            ii  = nint(abs(vec)*(ar_n-1)/M_PI) + 1
            if(ii <= ar_n) then ! should always be true
              do ik = 1, st%nik
                do p = st%st_start, st%st_end
                  arpis(ii, p, ik) = arpis(ii, p, ik) + &
                       st%occ(p, ik)*sum(abs(v%k(ix, iy, iz, :, p, ik))**2)
                  ar_npoints(ii) = ar_npoints(ii) + 1
                end do
              end do
            end if
          end if
        end if
        
      end do
    end do
  end do

  ! first output power spectra
  do ik = 1, st%nik
    do p = st%st_start, st%st_end
      call io_assign(iunit)
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_power.', ik, '.',p
      open(iunit, status='unknown', file=trim(fn))
      do ix = 1, n
        if(npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*step/units_out%energy%factor, spis(ix, p, ik), npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do
  call io_assign(iunit)
  write(fn, '(a,a)') trim(file), '_power.sum'
  open(iunit, status='unknown', file=trim(fn))
  do ix = 1, n
    if(npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*step/units_out%energy%factor, sum(spis(ix, :, :)), npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now output ar spectra
  do ik = 1, st%nik
    do p = st%st_start, st%st_end
      call io_assign(iunit)
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar.', ik, '.', p
      open(iunit, status='unknown', file=trim(fn))
      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, PRECISION), &
               arpis(ix, p, ik), ar_npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do
  call io_assign(iunit);
  write(fn, '(a,a)') trim(file), '_ar.sum'
  open(iunit, status='unknown', file=trim(fn))
  do ix = 1, ar_n
    if(ar_npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, PRECISION), &
           sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now we do the ar spectrum in real space
  arpis = M_ZERO
  ar_npoints = 0
  do ix = 1, m%l(1)
    ixx(1) = ix - m%l(1)/2 - 1
    do iy = 1, m%l(2)
      ixx(2) = iy - m%l(2)/2 - 1
      do iz = 1, m%l(3)
        ixx(3) = iz - m%l(3)/2 - 1
        
        ! angle resolved
        if(ixx(3)==0.and.(ixx(1).ne.0 .or. ixx(2).ne.0)) then
          vec = atan2(real(ixx(2), PRECISION), real(ixx(1), PRECISION))
          ii  = nint(abs(vec)*(ar_n-1)/M_PI) + 1
          if(ii <= ar_n) then ! should always be true
            do ik = 1, st%nik
              do p = st%st_start, st%st_end
                arpis(ii, p, ik) = arpis(ii, p, ik) + st%occ(p, ik)*v%r(ix, iy, iz, p, ik)
                ar_npoints(ii) = ar_npoints(ii) + 1
              end do
            end do
          end if
        end if
        
      end do
    end do
  end do

  ! now output real ar spectra
  do ik = 1, st%nik
    do p = st%st_start, st%st_end
      call io_assign(iunit)
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar_r.', ik, '.', p
      open(iunit, status='unknown', file=trim(fn))
      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, PRECISION), &
               arpis(ix, p, ik), ar_npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do
  call io_assign(iunit)
  write(fn, '(a,a)') trim(file), '_ar_r.sum'
  open(iunit, status='unknown', file=trim(fn))
  do ix = 1, ar_n
    if(ar_npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, PRECISION), &
           sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)
  
  deallocate(spis, arpis)
  deallocate(npoints, ar_npoints)

end subroutine PES_mask_output
