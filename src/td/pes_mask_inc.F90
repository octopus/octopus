!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine PES_mask_init(v, m, sb, st)
  type(PES_mask_t),  intent(out)   :: v
  type(mesh_t),      intent(inout) :: m
  type(simul_box_t), intent(in)    :: sb
  type(states_t),    intent(in)    :: st

  call push_sub('pes_mask_inc.PES_mask_init')

  message(1) = 'Info: Calculating PES using mask technique.'
  call write_info(1)

  ! alloc ffts in case they are not allocated yet
  call fft_init(m%idx%ll, fft_complex, v%fft, optimize = .not.simul_box_is_periodic(sb))

  ! setup arrays to be used
  SAFE_ALLOCATE(v%k(1:m%idx%ll(1), 1:m%idx%ll(2), 1:m%idx%ll(3), 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
  SAFE_ALLOCATE(v%r(1:m%idx%ll(1), 1:m%idx%ll(2), 1:m%idx%ll(3),             st%st_start:st%st_end, 1:st%d%nik))

  v%k = M_z0
  v%r = M_ZERO

  call pop_sub()
end subroutine PES_mask_init


! ---------------------------------------------------------
subroutine PES_mask_end(v)
  type(PES_mask_t), intent(inout) :: v

  call push_sub('pes_mask_inc.PES_mask_end')

  if(associated(v%k)) then
    call fft_end(v%fft)
    SAFE_DEALLOCATE_P(v%k)
    SAFE_DEALLOCATE_P(v%r)
  end if

  call pop_sub()
end subroutine PES_mask_end


! ---------------------------------------------------------
subroutine PES_mask_doit(v, m, st, dt, mask)
  type(PES_mask_t), intent(inout) :: v
  type(mesh_t),     intent(in)    :: m
  type(states_t),   intent(in)    :: st
  FLOAT,            intent(in)    :: dt
  FLOAT,            pointer       :: mask(:)

  integer :: j, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  CMPLX, allocatable :: wf1(:,:,:), wf2(:,:,:)
  FLOAT :: temp(MAX_DIM), vec

  call push_sub('pes_mask_inc.PES_mask_doit')

  ! propagate wavefunction in momentum space
  temp(:) = M_TWO*M_PI/(m%idx%ll(:)*m%h(:))
  do ix = 1, m%idx%ll(1)
    ixx(1) = pad_feq(ix, m%idx%ll(1), .true.)
    do iy = 1, m%idx%ll(2)
      ixx(2) = pad_feq(iy, m%idx%ll(2), .true.)
      do iz = 1, m%idx%ll(3)
        ixx(3) = pad_feq(iz, m%idx%ll(3), .true.)

        vec = sum((temp(:) * ixx(:))**2) / M_TWO
        v%k(ix, iy, iz,:,:,:) = v%k(ix, iy, iz,:,:,:)*exp(-M_zI*dt*vec)
      end do
    end do
  end do

  ! we now add the contribution from this timestep
  SAFE_ALLOCATE(wf1(1:m%idx%ll(1), 1:m%idx%ll(2), 1:m%idx%ll(3)))
  SAFE_ALLOCATE(wf2(1:m%idx%ll(1), 1:m%idx%ll(2), 1:m%idx%ll(3)))
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        wf1 = M_z0

        do j = 1, m%np
          ix3(:) = m%idx%Lxyz(j,:) + m%idx%ll(:)/2 + 1
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

  SAFE_DEALLOCATE_A(wf1)
  SAFE_DEALLOCATE_A(wf2)
  call pop_sub()
end subroutine PES_mask_doit


! ---------------------------------------------------------
subroutine PES_mask_output(v, m, st, file)
  type(PES_mask_t), intent(in) :: v
  type(mesh_t),     intent(in) :: m
  type(states_t),   intent(in) :: st
  character(len=*), intent(in) :: file

  FLOAT, allocatable :: spis(:,:,:), arpis(:,:,:)
  FLOAT :: vec, temp(MAX_DIM)
  integer, allocatable :: npoints(:), ar_npoints(:)
  integer :: p, ik, ii, ix, iy, iz, ixx(MAX_DIM), iunit
  character(len=100) :: fn

  integer,  parameter :: n = 600, ar_n = 90
  FLOAT, parameter :: step = CNST(0.005)

  call push_sub('pes_mask_inc.PES_mask_output')

  SAFE_ALLOCATE( spis(1:n,    st%st_start:st%st_end, 1:st%d%nik))
  SAFE_ALLOCATE(arpis(1:ar_n, st%st_start:st%st_end, 1:st%d%nik))
  SAFE_ALLOCATE(npoints(1:n))
  SAFE_ALLOCATE(ar_npoints(1:ar_n))
  spis = M_ZERO; arpis = M_ZERO
  npoints = 0;  ar_npoints = 0

  temp(:) = M_TWO*M_PI/(m%idx%ll(:)*m%h(:))
  do ix = 1, m%idx%ll(1)
    ixx(1) = pad_feq(ix, m%idx%ll(1), .true.)
    do iy = 1, m%idx%ll(2)
      ixx(2) = pad_feq(iy, m%idx%ll(2), .true.)
      do iz = 1, m%idx%ll(3)
        ixx(3) = pad_feq(iz, m%idx%ll(3), .true.)

        if(ixx(1).ne.0 .or. ixx(2).ne.0 .or. ixx(3).ne.0) then
          ! power spectrum
          vec = sum((temp(:) * ixx(:))**2) / M_TWO
          ii = nint(vec / step) + 1
          if(ii <= n) then
            do ik = st%d%kpt%start, st%d%kpt%end
              do p = st%st_start, st%st_end
                spis(ii, p, ik) = spis(ii, p, ik) + st%occ(p, ik) * &
                  sum(abs(v%k(ix, iy, iz, :, p, ik))**2)
                npoints(ii) = npoints(ii) + 1
              end do
            end do
          end if

          ! angle-resolved (assumes the pol is in the x-direction)
          if(ixx(3)==0 .and. (ixx(1).ne.0 .or. ixx(2).ne.0)) then
            vec = atan2(real(ixx(2), REAL_PRECISION), real(ixx(1), REAL_PRECISION))
            ii  = nint(abs(vec)*(ar_n-1)/M_PI) + 1
            if(ii <= ar_n) then ! should always be true
              do ik = st%d%kpt%start, st%d%kpt%end
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
  do ik = st%d%kpt%start, st%d%kpt%end
    do p = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_power.', ik, '.',p
      iunit = io_open(fn, action='write')

      do ix = 1, n
        if(npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*step/units_out%energy%factor, spis(ix, p, ik), npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_power.sum'
  iunit = io_open(fn, action='write')

  do ix = 1, n
    if(npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*step/units_out%energy%factor, sum(spis(ix, :, :)), npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now output ar spectra
  do ik = st%d%kpt%start, st%d%kpt%end
    do p = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar.', ik, '.', p
      iunit = io_open(fn, action='write')

      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, REAL_PRECISION), &
            arpis(ix, p, ik), ar_npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_ar.sum'
  iunit = io_open(fn, action='write')

  do ix = 1, ar_n
    if(ar_npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, REAL_PRECISION), &
        sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now we do the ar spectrum in real space
  arpis = M_ZERO
  ar_npoints = 0
  do ix = 1, m%idx%ll(1)
    ixx(1) = ix - m%idx%ll(1)/2 - 1
    do iy = 1, m%idx%ll(2)
      ixx(2) = iy - m%idx%ll(2)/2 - 1
      do iz = 1, m%idx%ll(3)
        ixx(3) = iz - m%idx%ll(3)/2 - 1

        ! angle-resolved
        if(ixx(3)==0.and.(ixx(1).ne.0 .or. ixx(2).ne.0)) then
          vec = atan2(real(ixx(2), REAL_PRECISION), real(ixx(1), REAL_PRECISION))
          ii  = nint(abs(vec)*(ar_n-1)/M_PI) + 1
          if(ii <= ar_n) then ! should always be true
            do ik = st%d%kpt%start, st%d%kpt%end
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
  do ik = st%d%kpt%start, st%d%kpt%end
    do p = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar_r.', ik, '.', p
      iunit = io_open(fn, action='write')

      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, REAL_PRECISION), &
            arpis(ix, p, ik), ar_npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_ar_r.sum'
  iunit = io_open(fn, action='write')
  do ix = 1, ar_n
    if(ar_npoints(ix) > 0) then
      write(iunit, *)  (ix-1)*CNST(180.0)/real(ar_n-1, REAL_PRECISION), &
        sum(arpis(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  SAFE_DEALLOCATE_A(spis)
  SAFE_DEALLOCATE_A(arpis)
  SAFE_DEALLOCATE_A(npoints)
  SAFE_DEALLOCATE_A(ar_npoints)

  call pop_sub()
end subroutine PES_mask_output

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
