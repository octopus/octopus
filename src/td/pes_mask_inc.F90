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
subroutine PES_mask_init(mask, mesh, sb, st, hm, max_iter,dt)
  type(PES_mask_t),    intent(out)   :: mask
  type(mesh_t),        intent(inout) :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(states_t),      intent(in)    :: st
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: max_iter
  FLOAT,               intent(in)    :: dt

  integer :: ll(MAX_DIM), il, it
  FLOAT :: field(MAX_DIM)

  PUSH_SUB(PES_mask_init)

  message(1) = 'Info: Calculating PES using mask technique.'
  call messages_info(1)

  ! allocate FFTs in case they are not allocated yet
  call fft_init(mesh%idx%ll,sb%dim,fft_complex,mask%fft, optimize = .not.simul_box_is_periodic(sb))

  ll(1:MAX_DIM) = mesh%idx%ll(1:MAX_DIM)

  ! setup arrays to be used
  SAFE_ALLOCATE(mask%k(1:ll(1),1:ll(2),1:ll(3),1:st%d%dim,1:st%nst,1:st%d%nik))
  SAFE_ALLOCATE(mask%r(1:ll(1),1:ll(2),1:ll(3), 1:st%nst, 1:st%d%nik))
  SAFE_ALLOCATE(mask%vec_pot(1:max_iter,1:MAX_DIM))

  mask%k = M_z0
  mask%r = M_ZERO
  mask%vec_pot=M_ZERO
  field=M_ZERO

  if (st%parallel_in_states) then
    message(1)= "Info: PES_mask parallelization on states experimental"
    call messages_info(1)

    if(mesh%parallel_in_domains) then
      write(message(1),'(a)') "PES_mask: simultaneous parallelization on mesh and states not supported"
      write(message(2),'(a)') "Modify ParallelizationStrategy and rerun." 
      call messages_fatal(2) 
      end if    
  endif



  ! Precalculate the potetial vector for all the simulation time
  do it = 2, max_iter
    do il = 1, hm%ep%no_lasers
      select case(laser_kind(hm%ep%lasers(il)))
      case(E_FIELD_ELECTRIC)
        call   laser_field(hm%ep%lasers(il), field(1:mesh%sb%dim), it*dt)
      case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
        write(message(1),'(a)') 'PES with mask method does not work with magnetic fields yet'
        call messages_fatal(2)
     end select
      mask%vec_pot(it,:)=mask%vec_pot(it,:)+field(:) !Sum up all the fields
   end do
    mask%vec_pot(it,:)= mask%vec_pot(it,:)+mask%vec_pot(it-1,:)
  end do

  mask%vec_pot(:,:)=-mask%vec_pot(:,:)*dt


  !%Variable PESMaskSpectEnergyMax 
  !%Type float
  !%Section Time-Dependent::PES
  !%Description
  !% The maximum energy for the PES spectrum (default 30 a.u.).
  !%End
  call parse_float(datasets_check('PESMaskSpectEnergyMax'),&
       units_to_atomic(units_inp%energy,CNST(30.0)),mask%energyMax)

  !%Variable PESMaskSpectEnergyStep 
  !%Type float
  !%Section Time-Dependent::PES
  !%Description
  !% The PES spectrum energy step (default 0.05 a.u.).
  !%End
  call parse_float(datasets_check('PESMaskSpectEnergyStep'),&
       units_to_atomic(units_inp%energy,CNST(0.05)),mask%energyStep)



  POP_SUB(PES_mask_init)
end subroutine PES_mask_init


! ---------------------------------------------------------
subroutine PES_mask_end(mask)
  type(PES_mask_t), intent(inout) :: mask

  PUSH_SUB(PES_mask_end)

  if(associated(mask%k)) then
    call fft_end(mask%fft)
    SAFE_DEALLOCATE_P(mask%k)
    SAFE_DEALLOCATE_P(mask%r)
    SAFE_DEALLOCATE_P(mask%vec_pot)
  end if

  POP_SUB(PES_mask_end)
end subroutine PES_mask_end


! ---------------------------------------------------------
subroutine PES_mask_calc(mask, mesh, st, dt, mask_fn,hm,geo,iter)
  type(PES_mask_t), intent(inout) :: mask
  type(mesh_t),     intent(in)    :: mesh
  type(states_t),   intent(in)    :: st
  integer,          intent(in)    :: iter
  FLOAT,            intent(in)    :: dt
  FLOAT,            pointer       :: mask_fn(:) !namely hm%ab_pot
  type(hamiltonian_t),   intent(in)    :: hm
  type(geometry_t), intent(in)    :: geo

  integer :: ip, idim, ist, ik, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
  CMPLX, allocatable :: wf1(:,:,:), wf2(:,:,:), psi(:)
  FLOAT :: temp(MAX_DIM), vec

  integer :: ip_local,size

#if defined(HAVE_MPI)
  integer :: status(MPI_STATUS_SIZE)
  integer :: iproc,dataSize
#endif

  PUSH_SUB(PES_mask_calc)


  ! propagate wavefunction in momentum space in presence of a laser field 
  temp(:) = M_TWO * M_PI / (mesh%idx%ll(:) * mesh%spacing(:))
  do ix = 1, mesh%idx%ll(1)
    ixx(1) = pad_feq(ix, mesh%idx%ll(1), .true.)
    do iy = 1, mesh%idx%ll(2)
      ixx(2) = pad_feq(iy, mesh%idx%ll(2), .true.)
      do iz = 1, mesh%idx%ll(3)
        ixx(3) = pad_feq(iz, mesh%idx%ll(3), .true.)

        vec = sum((temp(1:mesh%sb%dim) * ixx(1:mesh%sb%dim)&
             -mask%vec_pot(iter,1:mesh%sb%dim))**2) / M_TWO

        mask%k(ix, iy, iz,:,:,:) = mask%k(ix, iy, iz,:,:,:) * exp(-M_zI * dt * vec)
      end do
    end do
  end do


  ! we now add the contribution from this timestep
  SAFE_ALLOCATE(psi(1:mesh%np))
  SAFE_ALLOCATE(wf1(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))
  SAFE_ALLOCATE(wf2(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))

  size = (mesh%idx%ll(1))*(mesh%idx%ll(2))*(mesh%idx%ll(3)) 

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        wf1 = M_z0
        wf2 = M_z0

        call states_get_state(st, mesh, idim, ist, ik, psi)

        do ip_local = 1, mesh%np

          ! Convert from local to global mesh index
          ip = index_from_coords(mesh%idx,mesh%sb%dim, int(mesh%x(ip_local, :)/mesh%spacing(:)))

          ix3(:) = mesh%idx%lxyz(ip, :) + mesh%idx%ll(:)/2 + 1 

          wf1(ix3(1), ix3(2), ix3(3)) = mask_fn(ip_local)*psi(ip_local)
          
        end do

#if defined(HAVE_MPI)

        if(mesh%parallel_in_domains)then

           !send all the wavefunctions to root node
           if(mesh%mpi_grp%rank .gt. 0 ) then
              call MPI_Send(wf1,size,& 
                   MPI_CMPLX,0,666, mesh%mpi_grp%comm, mpi_err)
           else
              do iproc= 1, mesh%mpi_grp%size-1
                 call MPI_Recv(wf2,size,&
                      MPI_CMPLX,iproc,666, mesh%mpi_grp%comm, status, mpi_err)
                ! add contribute for other nodes
                 wf1 = wf1 +wf2 
              end do
 
           end if
        end if

#endif

        ! and add to our density (sum for idim, also)
        mask%r(:,:,:, ist, ik) = mask%r(:,:,:, ist, ik) + abs(wf1)**2

        ! now we FT
        call zfft_forward(mask%fft, wf1, wf2)

        ! and add to our spectrum
        mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf2

      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi)

#ifdef HAVE_MPI
 ! wait for all processors to finish!
  if(st%mpi_grp%size .gt. 1 ) then
     call MPI_Barrier(st%mpi_grp%comm, mpi_err)
   end if
#endif

  SAFE_DEALLOCATE_A(wf1)
  SAFE_DEALLOCATE_A(wf2)

  POP_SUB(PES_mask_calc)
end subroutine PES_mask_calc


!---------------------------------------------------------------------------
! Collect the states from all the nodes when the code run parallel on states
! --------------------------------------------------------------------------
subroutine PES_mask_collect(mask, st,mesh)
  type(PES_mask_t), intent(inout) :: mask
  type(states_t),   intent(in) :: st
  type(mesh_t),     intent(in)    :: mesh

#ifdef HAVE_MPI
  FLOAT, allocatable :: wfr(:,:,:)
  CMPLX, allocatable :: wf(:,:,:)
  integer :: iproc, size
  integer :: idim, ist, ik
  integer :: status(MPI_STATUS_SIZE)
#endif

 PUSH_SUB(PES_mask_collect)

#ifdef HAVE_MPI

  SAFE_ALLOCATE(wf(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))
  SAFE_ALLOCATE(wfr(1:mesh%idx%ll(1), 1:mesh%idx%ll(2), 1:mesh%idx%ll(3)))

  size = product(mesh%idx%ll(1:3))

  do ik = 1, st%d%nik
     do ist = 1, st%nst

	if(in_debug_mode) then

          if(st%mpi_grp%rank .gt. 0 ) then
            wfr = mask%r(:, :, :, ist, ik) 
            call MPI_Send(wfr(1, 1, 1), size, MPI_FLOAT, 0, 2, st%mpi_grp%comm, mpi_err)
          else
            do iproc = 1, st%mpi_grp%size-1 
              call MPI_Recv(wfr(1, 1, 1), size, MPI_FLOAT, iproc, 2, st%mpi_grp%comm, status, mpi_err)
              mask%r(:,:,:, ist, ik) = mask%r(:,:,:, ist, ik) + wfr
            end do
          end if
          
        end if

        do idim = 1, st%d%dim

           !send all the wavefunctions to root node
           if(st%mpi_grp%rank .gt. 0 ) then
              wf= mask%k(:, :, :, idim, ist, ik) 
              call MPI_Send(wf(1, 1, 1), size, MPI_CMPLX, 0, 1, st%mpi_grp%comm, mpi_err)
           else
              !root node collects all the data
              do iproc= 1, st%mpi_grp%size-1 
                 call MPI_Recv(wf(1, 1, 1), size, MPI_CMPLX,iproc, 1, st%mpi_grp%comm, status, mpi_err)
                 mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) + wf
              end do
           end if

        end do
     end do
  end do

  !wait for all thred to finish 
  if(st%mpi_grp%size .gt. 1 ) then
     call MPI_Barrier(st%mpi_grp%comm, mpi_err)
  end if


  SAFE_DEALLOCATE_A(wf)
  SAFE_DEALLOCATE_A(wfr)

#endif

  POP_SUB(PES_mask_collect)
end subroutine PES_mask_collect




! ---------------------------------------------------------
subroutine PES_mask_output(mask, mesh, st, file)
  type(PES_mask_t), intent(in) :: mask
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st
  character(len=*), intent(in) :: file

  FLOAT, allocatable :: spis(:,:,:), arpes(:,:,:)
  FLOAT :: vec, temp(MAX_DIM),rdens
  FLOAT, allocatable :: npoints(:), ar_npoints(:)
  integer :: ist, ik, ii, ix, iy, iz, ixx(MAX_DIM), iunit
  character(len=100) :: fn

  integer,  parameter ::  ar_n = 90
  integer :: nn 
  FLOAT  :: step

  step=mask%energyStep
  nn  = nint(mask%energyMax/step)

  PUSH_SUB(PES_mask_output)

  SAFE_ALLOCATE( spis(1:nn,   1:st%nst, 1:st%d%nik))
  SAFE_ALLOCATE(arpes(1:ar_n, 1:st%nst, 1:st%d%nik))
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
          ! the power spectrum
          vec = sum((temp(:) * ixx(:))**2) / M_TWO
          ii = nint(vec / step) + 1
          if(ii <= nn) then
            do ik = 1,st%d%nik
              do ist = 1, st%nst
                spis(ii, ist, ik) = spis(ii, ist, ik) + st%occ(ist, ik) * &
                  sum(abs(mask%k(ix, iy, iz, :, ist, ik))**2)
                npoints(ii) = npoints(ii) + st%occ(ist, ik) ! count only occupied states
              end do
            end do
          end if

          ! angle-resolved (assumes the pol is in the x-direction)
          if(ixx(3)==0 .and. (ixx(1).ne.0 .or. ixx(2).ne.0)) then
            vec = atan2(real(ixx(2), REAL_PRECISION), real(ixx(1), REAL_PRECISION))
            ii  = nint(abs(vec) * (ar_n - 1)/M_PI) + 1
            if(ii <= ar_n) then ! should always be true
              do ik = 1, st%d%kpt%nglobal
                do ist = 1, st%nst
                  arpes(ii, ist, ik) = arpes(ii, ist, ik) + &
                    st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik))**2)
                  ar_npoints(ii) = ar_npoints(ii) +  st%occ(ist, ik) ! count only occupied states
                end do
              end do
            end if
          end if
        end if

      end do
    end do
  end do

  ! first output power spectra
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      write(fn , '(a,a,i1.1,a,i2.2)') trim(file), '_power.', ik, '.', ist
      iunit = io_open(fn, action='write')

      do ix = 1, nn
        if(npoints(ix) > 0) then
          write(iunit, *)  units_from_atomic(units_out%energy, (ix - 1) * step), spis(ix, ist, ik)/npoints(ix), npoints(ix)
        end if
      end do
      call io_close(iunit)
    end do
  end do

  write(fn, '(a,a)') trim(file), '_power.sum'
  iunit = io_open(fn, action='write')

  do ix = 1, nn
    if(npoints(ix) > 0) then
      write(iunit, *)  units_from_atomic(units_out%energy, (ix - 1) * step), sum(spis(ix, :, :))/npoints(ix), npoints(ix)
    end if
  end do
  call io_close(iunit)

  ! now output ar spectra
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      write(fn, '(a,a,i1.1,a,i2.2)') trim(file), '_ar.', ik, '.', ist
      iunit = io_open(fn, action='write')

      do ix = 1, ar_n
        if(ar_npoints(ix) > 0) then
          write(iunit, *)  (ix-1)*CNST(180.0) / real(ar_n-1, REAL_PRECISION), &
            arpes(ix, ist, ik)/ar_npoints(ix), ar_npoints(ix)
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
        sum(arpes(ix, :, :))/ar_npoints(ix), ar_npoints(ix)
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
            do ik = 1, st%d%nik
              do ist = 1, st%nst
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
  do ik = 1, st%d%nik
    do ist = 1, st%nst
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

if(in_debug_mode) then

  !Now the total PES density cut along the x_axis
  !since this function is mainly for debug purposes I don`t see a reason to expand this part any further 
  write(fn, '(a,a)') trim(file), '_den.sum'
  iunit = io_open(fn, action='write')
  do ix = 1, mesh%idx%ll(1)
     iy= mesh%idx%ll(2)/2+1
     iz= mesh%idx%ll(3)/2+1

     rdens=M_ZERO
     do ik = 1, st%d%nik
        do ist = 1, st%nst
           rdens = rdens+ mask%r(ix,iy,iz,ist,ik)*st%occ(ist, ik)               
        end do
     end do 

     write(iunit,*)&
        units_from_atomic(units_inp%length,(ix-mesh%idx%ll(1)/2-1)*mesh%spacing(1)),rdens,&
        sum( mask%r(ix,iy,iz,:,:))
           
  end do
end if

  call io_close(iunit)

  SAFE_DEALLOCATE_A(spis)
  SAFE_DEALLOCATE_A(arpes)
  SAFE_DEALLOCATE_A(npoints)
  SAFE_DEALLOCATE_A(ar_npoints)

  POP_SUB(PES_mask_output)
end subroutine PES_mask_output



! ---------------------------------------------------------
subroutine PES_mask_restart_write(mask, mesh, st)
  type(PES_mask_t), intent(in) :: mask
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st

  character(len=80) :: filename, dir ,path
  integer :: itot, ik, ist, idim , np, ierr
  integer :: ll(MAX_DIM)

  PUSH_SUB(PES_mask_restart_write)

  ll(1:MAX_DIM) = mesh%idx%ll(1:MAX_DIM)
  np =ll(1)*ll(2)*ll(3); 

  dir=trim(tmpdir)//'td/'

  itot = 1

 !assumes that only the main thread is allowed to dump the restart info
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim

        write(filename,'(i10.10)') itot

        path=trim(dir)//'pes_'//trim(filename)//'.obf'
        
        call io_binary_write(path,np, mask%k(:,:,:, idim, ist, ik), ierr)


        itot = itot + 1
      end do
    end do
  end do



  POP_SUB(PES_mask_restart_write)
end subroutine PES_mask_restart_write

! ---------------------------------------------------------
subroutine PES_mask_restart_read(mask, mesh, st)
  type(PES_mask_t), intent(inout) :: mask
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st

  character(len=80) :: filename, dir ,path

  integer :: itot, ik, ist, idim , np, ierr
  integer :: ll(MAX_DIM)

  PUSH_SUB(PES_mask_restart_read)

  ll(1:MAX_DIM) = mesh%idx%ll(1:MAX_DIM)
  np =ll(1)*ll(2)*ll(3); 


  dir=trim(tmpdir)//'td/'

  itot = 1
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim
         write(filename,'(i10.10)') itot

        path=trim(dir)//'pes_'//trim(filename)//'.obf'
        
        call io_binary_read(path,np, mask%k(:,:,:, idim, ist, ik), ierr)

        itot = itot + 1
      end do
    end do
  end do

  POP_SUB(PES_mask_restart_read)
end subroutine PES_mask_restart_read



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
