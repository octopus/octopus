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
subroutine PES_mask_init(mask, mesh, sb, st)
  type(PES_mask_t),  intent(out)   :: mask
  type(mesh_t),      intent(inout) :: mesh
  type(simul_box_t), intent(in)    :: sb
  type(states_t),    intent(in)    :: st

  integer :: ll(MAX_DIM)

  call push_sub('pes_mask_inc.PES_mask_init')

  message(1) = 'Info: Calculating PES using mask technique.'
  call write_info(1)

  ! allocate FFTs in case they are not allocated yet
  call fft_init(mesh%idx%ll, fft_complex, mask%fft, optimize = .not.simul_box_is_periodic(sb))

  ll(1:MAX_DIM) = mesh%idx%ll(1:MAX_DIM)

  ! setup arrays to be used
  SAFE_ALLOCATE(mask%k(1:ll(1),1:ll(2),1:ll(3),1:st%d%dim,st%st_start:st%st_end,1:st%d%nik))
  SAFE_ALLOCATE(mask%r(1:ll(1),1:ll(2),1:ll(3), st%st_start:st%st_end, 1:st%d%nik))

  mask%k = M_z0
  mask%r = M_ZERO

  call pop_sub('pes_mask_inc.PES_mask_init')
end subroutine PES_mask_init


! ---------------------------------------------------------
subroutine PES_mask_end(mask)
  type(PES_mask_t), intent(inout) :: mask

  call push_sub('pes_mask_inc.PES_mask_end')

  if(associated(mask%k)) then
    call fft_end(mask%fft)
    SAFE_DEALLOCATE_P(mask%k)
    SAFE_DEALLOCATE_P(mask%r)
  end if

  call pop_sub('pes_mask_inc.PES_mask_end')
end subroutine PES_mask_end


! ---------------------------------------------------------
subroutine PES_mask_calc(mask, mesh, st, dt, mask_fn)
  type(PES_mask_t), intent(inout) :: mask
  type(mesh_t),     intent(in)    :: mesh
  type(states_t),   intent(in)    :: st
  FLOAT,            intent(in)    :: dt
  FLOAT,            pointer       :: mask_fn(:)

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  CMPLX, allocatable :: wf1(:,:,:), wf2(:,:,:)
  FLOAT :: temp(MAX_DIM), vec

  call push_sub('pes_mask_inc.PES_mask_calc')

  ! propagate wavefunction in momentum space
  temp(:) = M_TWO * M_PI / (mesh%idx%ll(:) * mesh%spacing(:))
  do ix = 1, mesh%idx%ll(1)
    ixx(1) = pad_feq(ix, mesh%idx%ll(1), .true.)
    do iy = 1, mesh%idx%ll(2)
      ixx(2) = pad_feq(iy, mesh%idx%ll(2), .true.)
      do iz = 1, mesh%idx%ll(3)
        ixx(3) = pad_feq(iz, mesh%idx%ll(3), .true.)

        vec = sum((temp(:) * ixx(:))**2) / M_TWO
        mask%k(ix, iy, iz,:,:,:) = mask%k(ix, iy, iz,:,:,:) * exp(-M_zI * dt * vec)
      end do
    end do
  end do

  ! we now add the contribution from this timestep
  SAFE_ALLOCATE(wf1(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))
  SAFE_ALLOCATE(wf2(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        wf1 = M_z0

        do ip = 1, mesh%np
          ix3(:) = mesh%idx%Lxyz(ip, :) + mesh%idx%ll(:)/2 + 1
          wf1(ix3(1), ix3(2), ix3(3)) = mask_fn(ip) * st%zpsi(ip, idim, ist, ik)
        end do

        ! and add to our density (sum for idim, also)
        mask%r(:,:,:, ist, ik) = mask%r(:,:,:, ist, ik) + abs(wf1)**2

        ! now we FT
        call zfft_forward(mask%fft, wf1, wf2)

        ! and add to our spectrum
        mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf2

      end do
    end do
  end do

  SAFE_DEALLOCATE_A(wf1)
  SAFE_DEALLOCATE_A(wf2)
  call pop_sub('pes_mask_inc.PES_mask_calc')
end subroutine PES_mask_calc


! ---------------------------------------------------------
subroutine PES_mask_output(mask, mesh, st, file)
  type(PES_mask_t), intent(in) :: mask
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st
  character(len=*), intent(in) :: file

  FLOAT, allocatable :: spis(:,:,:), arpes(:,:,:)
  FLOAT :: vec, temp(MAX_DIM)
  integer, allocatable :: npoints(:), ar_npoints(:)
  integer :: ist, ik, ii, ix, iy, iz, ixx(MAX_DIM), iunit
  character(len=100) :: fn

  integer,  parameter :: nn = 600, ar_n = 90
  FLOAT, parameter :: step = CNST(0.005)

  call push_sub('pes_mask_inc.PES_mask_output')

  SAFE_ALLOCATE( spis(1:nn,   st%st_start:st%st_end, 1:st%d%nik))
  SAFE_ALLOCATE(arpes(1:ar_n, st%st_start:st%st_end, 1:st%d%nik))
  SAFE_ALLOCATE(npoints(1:nn))
  SAFE_ALLOCATE(ar_npoints(1:ar_n))
  spis = M_ZERO
  arpes = M_ZERO
  npoints = 0
  ar_npoints = 0

  temp(:) = M_TWO * M_PI / (mesh%idx%ll(:) * mesh%spacing(:))
  do ix = 1, mesh%idx%ll(1)
    ixx(1) = pad_feq(ix, mesh%idx%ll(1), .true.)
    do iy = 1, mesh%idx%ll(2)
      ixx(2) = pad_feq(iy, mesh%idx%ll(2), .true.)
      do iz = 1, mesh%idx%ll(3)
        ixx(3) = pad_feq(iz, mesh%idx%ll(3), .true.)

        if(ixx(1).ne.0 .or. ixx(2).ne.0 .or. ixx(3).ne.0) then
          ! power spectrum
          vec = sum((temp(:) * ixx(:))**2) / M_TWO
          ii = nint(vec / step) + 1
          if(ii <= nn) then
            do ik = st%d%kpt%start, st%d%kpt%end
              do ist = st%st_start, st%st_end
                spis(ii, ist, ik) = spis(ii, ist, ik) + st%occ(ist, ik) * &
                  sum(abs(mask%k(ix, iy, iz, :, ist, ik))**2)
                npoints(ii) = npoints(ii) + 1
              end do
            end do
          end if

          ! angle-resolved (assumes the pol is in the x-direction)
          if(ixx(3)==0 .and. (ixx(1).ne.0 .or. ixx(2).ne.0)) then
            vec = atan2(real(ixx(2), REAL_PRECISION), real(ixx(1), REAL_PRECISION))
            ii  = nint(abs(vec) * (ar_n - 1)/M_PI) + 1
            if(ii <= ar_n) then ! should always be true
              do ik = st%d%kpt%start, st%d%kpt%end
                do ist = st%st_start, st%st_end
                  arpes(ii, ist, ik) = arpes(ii, ist, ik) + &
                    st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik))**2)
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
    do ist = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_power.', ik, '.', ist
      iunit = io_open(fn, action='write')

      do ix = 1, nn
        if(npoints(ix) > 0) then
          write(iunit, *)  units_from_atomic(units_out%energy, (ix - 1) * step), spis(ix, ist, ik), npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_power.sum'
  iunit = io_open(fn, action='write')

  do ix = 1, nn
    if(npoints(ix) > 0) then
      write(iunit, *)  units_from_atomic(units_out%energy, (ix - 1) * step), sum(spis(ix, :, :)), npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now output ar spectra
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar.', ik, '.', ist
      iunit = io_open(fn, action='write')

      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0) / real(ar_n-1, REAL_PRECISION), &
            arpes(ix, ist, ik), ar_npoints(ix)
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
        sum(arpes(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now we do the ar spectrum in real space
  arpes = M_ZERO
  ar_npoints = 0
  do ix = 1, mesh%idx%ll(1)
    ixx(1) = ix - mesh%idx%ll(1)/2 - 1
    do iy = 1, mesh%idx%ll(2)
      ixx(2) = iy - mesh%idx%ll(2)/2 - 1
      do iz = 1, mesh%idx%ll(3)
        ixx(3) = iz - mesh%idx%ll(3)/2 - 1

        ! angle-resolved
        if(ixx(3) == 0 .and. (ixx(1) .ne. 0 .or. ixx(2) .ne. 0)) then
          vec = atan2(real(ixx(2), REAL_PRECISION), real(ixx(1), REAL_PRECISION))
          ii  = nint(abs(vec) * (ar_n-1) / M_PI) + 1
          if(ii <= ar_n) then ! should always be true
            do ik = st%d%kpt%start, st%d%kpt%end
              do ist = st%st_start, st%st_end
                arpes(ii, ist, ik) = arpes(ii, ist, ik) + st%occ(ist, ik) * mask%r(ix, iy, iz, ist, ik)
                ar_npoints(ii) = ar_npoints(ii) + 1
              end do
            end do
          end if
        end if

      end do
    end do
  end do

  ! now output real angle-resolved spectra
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar_r.', ik, '.', ist
      iunit = io_open(fn, action='write')

      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *) (ix - 1) * CNST(180.0) / real(ar_n - 1, REAL_PRECISION), &
            arpes(ix, ist, ik), ar_npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_ar_r.sum'
  iunit = io_open(fn, action='write')
  do ix = 1, ar_n
    if(ar_npoints(ix) > 0) then
      write(iunit, *) (ix - 1) * CNST(180.0) / real(ar_n - 1, REAL_PRECISION), &
        sum(arpes(ix, :, :)), ar_npoints(ix)
    end if
  end do
  call io_close(iunit)

  SAFE_DEALLOCATE_A(spis)
  SAFE_DEALLOCATE_A(arpes)
  SAFE_DEALLOCATE_A(npoints)
  SAFE_DEALLOCATE_A(ar_npoints)

  call pop_sub('pes_mask_inc.PES_mask_output')
end subroutine PES_mask_output

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
