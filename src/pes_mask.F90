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
  type(PES_mask_type), intent(out) :: v
  type(mesh_type), intent(inout) :: m
  type(states_type), intent(IN) :: st

  message(1) = 'Info: Calculating PES using mask technique'
  call write_info(1)
      
  ! alloc ffts in case they are not allocated yet
  call mesh_alloc_ffts(m, 1)

  ! setup arrays to be used
  allocate(&
       v%k(m%fft_n(1), m%fft_n(2), m%fft_n(3), st%dim, st%st_start:st%st_end, st%nik), &
       v%r(m%fft_n(1), m%fft_n(2), m%fft_n(3),         st%st_start:st%st_end, st%nik))

  v%k = M_z0
  v%r = 0._r8
   
end subroutine PES_mask_init

subroutine PES_mask_end(v)
  type(PES_mask_type), intent(inout) :: v

  if(associated(v%k)) then
    deallocate(v%k, v%r)
    nullify   (v%k, v%r)
  end if

end subroutine PES_mask_end

subroutine PES_mask_doit(v, m, st, dt, mask)
  type(PES_mask_type), intent(inout) :: v
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(IN) :: dt
  real(r8), pointer :: mask(:)

  integer :: j, idim, ist, ik, ix, iy, iz, ix3(3), ixx(3)
  complex(r8), allocatable :: wf1(:,:,:), wf2(:,:,:)
  real(r8) :: temp(3), vec

  ! propagate wave-function in momentum space
  temp(:) = 2.0_r8*M_PI/(m%fft_n(:)*m%h(:))
  do ix = 1, m%fft_n(1)
    ixx(1) = pad_feq(ix, m%fft_n(1), .true.)
    do iy = 1, m%fft_n(2)
      ixx(2) = pad_feq(iy, m%fft_n(2), .true.)
      do iz = 1, m%fft_n(3)
        ixx(3) = pad_feq(iz, m%fft_n(3), .true.)

        vec = sum((temp(:) * ixx(:))**2) / 2._r8
        v%k(ix, iy, iz,:,:,:) = v%k(ix, iy, iz,:,:,:)*exp(-M_zI*dt*vec)
      enddo
    enddo
  enddo
      
  ! we now add the contribution from this timestep
  allocate( &
       wf1(m%fft_n(1), m%fft_n(2), m%fft_n(3)), &
       wf2(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      do idim = 1, st%dim

        wf1 = M_z0

        do j = 1, m%np
          ix3(:) = m%Lxyz(:,j) + m%fft_n(:)/2 + 1
          wf1(ix3(1), ix3(2), ix3(3)) = mask(j)*st%zpsi(j, idim, ist, ik)
        end do

        ! and add to our density (sum for idim, also)
        v%r(:,:,:, ist, ik) = v%r(:,:,:, ist, ik) + abs(wf1)**2

        ! now we FT
        call fftwnd_f77_one(m%zplanf, wf1, wf2)

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
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  character(len=*), intent(in) :: file

  real(r8), allocatable :: spis(:,:,:), arpis(:,:,:)
  real(r8) :: vec, temp(3)
  integer, allocatable :: npoints(:), ar_npoints(:)
  integer :: p, ik, ii, ix, iy, iz, ixx(3), iunit
  character(len=100) :: fn
  
  integer,  parameter :: n = 600, ar_n = 90
  real(r8), parameter :: step = 0.005_r8

  allocate( &
       spis (n,    st%st_start:st%st_end, st%nik), &
       arpis(ar_n, st%st_start:st%st_end, st%nik))
  allocate(npoints(n), ar_npoints(ar_n))
  spis = 0._r8; arpis = 0._r8
  npoints = 0;  ar_npoints = 0
  
  temp(:) = 2.0_r8*M_PI/(m%fft_n(:)*m%h(:))
  do ix = 1, m%fft_n(1)
    ixx(1) = pad_feq(ix, m%fft_n(1), .true.)
    do iy = 1, m%fft_n(2)
      ixx(2) = pad_feq(iy, m%fft_n(2), .true.)
      do iz = 1, m%fft_n(3)
        ixx(3) = pad_feq(iz, m%fft_n(3), .true.)
        
        if(ixx(1).ne.0 .or. ixx(2).ne.0 .or. ixx(3).ne.0) then
          ! power spectrum
          vec = sum((temp(:) * ixx(:))**2) / 2._r8
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
            vec = atan2(real(ixx(2), r8), real(ixx(1), r8))
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
          write(iunit, *)  (ix-1)*180._r8/real(ar_n-1, r8), &
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
      write(iunit, *)  (ix-1)*180._r8/real(ar_n-1, r8), &
           sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now we do the ar spectrum in real space
  arpis = 0._r8
  ar_npoints = 0
  do ix = 1, m%fft_n(1)
    ixx(1) = ix - m%fft_n(1)/2 - 1
    do iy = 1, m%fft_n(2)
      ixx(2) = iy - m%fft_n(2)/2 - 1
      do iz = 1, m%fft_n(3)
        ixx(3) = iz - m%fft_n(3)/2 - 1
        
        ! angle resolved
        if(ixx(3)==0.and.(ixx(1).ne.0 .or. ixx(2).ne.0)) then
          vec = atan2(real(ixx(2), r8), real(ixx(1), r8))
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
          write(iunit, *)  (ix-1)*180._r8/real(ar_n-1, r8), &
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
      write(iunit, *)  (ix-1)*180._r8/real(ar_n-1, r8), &
           sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)
  
  deallocate(spis, arpis)
  deallocate(npoints, ar_npoints)

end subroutine PES_mask_output
