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


! ---------------------------------------------------------
! Write the photoelectron wavefunctions in real space
! ---------------------------------------------------------
  subroutine PES_mask_output_states(st, gr, geo, dir, outp, mask)

    type(states_t),         intent(in) :: st
    type(grid_t),           intent(in) :: gr
    type(geometry_t),       intent(in)    :: geo
    character(len=*),       intent(in)    :: dir
    type(output_t),   intent(in)    :: outp
    type(PES_mask_t),       intent(in)    :: mask

    integer :: ik, ist, idim, idir, is, ierr, ip
    character(len=80) :: fname
    type(unit_t) :: fn_unit

    integer :: ip_local, ix, iy, iz, ix3(MAX_DIM), ixx(MAX_DIM)
    CMPLX, allocatable :: PsiAB(:,:,:,:),wf1(:,:,:)
    FLOAT,allocatable :: RhoAB(:,:) 
    FLOAT :: temp(MAX_DIM), vec
    FLOAT :: dd
    integer :: il,ll(MAX_DIM)
    type(mesh_t):: mesh   

    type(batch_t)        :: psib
    type(density_calc_t) :: dens_calc

    

    PUSH_SUB(PES_mask_output_states)

    mesh= gr%mesh
!    ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
    SAFE_ALLOCATE(wf1(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(PsiAB(1:mesh%np_part,1:st%d%dim,1:st%nst,1:st%d%nik))
    SAFE_ALLOCATE(RhoAB(1:mesh%np_part,1:st%d%nspin))
                 
  RhoAB= M_ZERO
  
  !Calculate the pes density \Psi_A + \Psi_B on the simulation Box  
  call density_calc_init(dens_calc, st, gr, RhoAB)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        wf1 = M_z0

        ! back to real space
        call PES_mask_K_to_X(mask,mesh,mask%k(:,:,:, idim, ist, ik),wf1)        


!        do ip_local = 1, mesh%np_global
!
!          ! Convert from local to global mesh index
!          ip= index_from_coords(mesh%idx,mesh%sb%dim,nint(mesh%x(ip_local,:)/mask%spacing(:)))
!
!          ix3(:) =  mesh%idx%Lxyz(ip, :) + mask%ll(:)/2 + 1 
!          
!          if (mask%mode .ne. MODE_EXACT) then 
!            PsiAB(ip_local, idim, ist, ik) = st%zpsi(ip_local, idim, ist, ik) + wf1(ix3(1), ix3(2), ix3(3)) 
!          else 
!            PsiAB(ip_local, idim, ist, ik) =  wf1(ix3(1), ix3(2), ix3(3)) 
!          end if 
! 
!      end do

        call PES_mask_square_to_mesh(mask,mesh,PsiAB(:, idim, ist, ik),wf1,MaskHow = 3,Const = M_ONE)

        if (mask%mode .ne. MODE_EXACT) then 
          PsiAB(:, idim, ist, ik) = PsiAB(:, idim, ist, ik) + st%zpsi(:, idim, ist, ik) 
        end if
 
      end do
    end do
     
    call batch_init(psib, st%d%dim, st%st_start, st%st_end, PsiAB(:, :, st%st_start:, ik))
    call density_calc_accumulate(dens_calc, ik, psib) 
    call batch_end(psib)

  end do
   
  call density_calc_end(dens_calc)


  ! THE OUTPUT 
  if(iand(outp%what, C_OUTPUT_PES_DENSITY) .ne. 0) then
    fn_unit = units_out%length**(-gr%mesh%sb%dim)
    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(a)') 'pes_den'
      else
        write(fname, '(a,i1)') 'pes_den-sp', is
      endif
      call dio_function_output(outp%how, dir, fname, gr%fine%mesh, &
        RhoAB(:, is), fn_unit, ierr, is_tmp = .false., geo = geo, grp = st%mpi_grp)
    end do
  end if

  if(iand(outp%what, C_OUTPUT_PES_WFS).ne.0) then
    fn_unit = sqrt(units_out%length**(-gr%mesh%sb%dim))
    do ist = st%st_start, st%st_end
!        if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i4.4,a,i1)') 'pes_wf-k', ik, '-st', ist, '-sp', idim
              else
                write(fname, '(a,i3.3,a,i4.4)')      'pes_wf-k', ik, '-st', ist
              endif
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i4.4,a,i1)')        'pes_wf-st', ist, '-sp', idim
              else
                write(fname, '(a,i4.4)')             'pes_wf-st', ist
              endif
            endif
              
            call zio_function_output(outp%how, dir, fname, gr%mesh, &
                   PsiAB(1:, idim, ist, ik), fn_unit, ierr, is_tmp = .false., geo = geo)

          end do
        end do
 !       end if
    end do
  end if


  SAFE_DEALLOCATE_A(wf1)
  SAFE_DEALLOCATE_A(PsiAB)
  SAFE_DEALLOCATE_A(RhoAB)

  POP_SUB(PES_mask_output_states)
  end subroutine PES_mask_output_states

! ---------------------------------------------------------
!
! Calculates the momentum resolved photoelectron probability
!            P(k) = \sum_i |\Psi_{B,i}(k)|^2 
!
! ---------------------------------------------------------
subroutine PES_mask_create_full_map(mask, st, PESK, wfAk)
  type(PES_mask_t),     intent(in)  :: mask
  type(states_t),       intent(in)  :: st
  FLOAT,                intent(out) :: PESK(:,:,:)
  CMPLX, optional,      intent(in) :: wfAk(:,:,:,:,:,:)

  integer :: ist, ik, ii, kx, ky, kz,idim
  FLOAT   :: scale


  PUSH_SUB(PES_mask_create_full_map)

  PESK = M_ZERO


  do kx = 1, mask%ll(1)
    do ky = 1, mask%ll(2)
      do kz = 1, mask%ll(3)

        do ik = 1,st%d%nik
          do ist = 1, st%nst
            if(present(wfAk))then
              PESK(kx,ky,kz) = PESK(kx,ky,kz) + st%occ(ist, ik) *&
                       sum(abs(mask%k(kx, ky, kz, :, ist, ik) + wfAk(kx,ky,kz,:, ist, ik)  )**2)
            else
              PESK(kx,ky,kz) = PESK(kx,ky,kz) + st%occ(ist, ik) * sum(abs(mask%k(kx, ky, kz, :, ist, ik)  )**2)
            end if
          end do
        end do

      end do
    end do
  end do

  ! This is needed in order to normalize the Fourier integral 
  scale = M_ONE
  do idim=1, mask%mesh%sb%dim
    scale = scale *( mask%spacing(idim)/sqrt(M_TWO*M_PI))**2
  end do
  PESK = PESK *scale 
   


  POP_SUB(PES_mask_create_full_map)
  end subroutine PES_mask_create_full_map

! --------------------------------------------------------
!
!  Qshep interpolation helper function initialization.
!  Generates the linearized version of PESK (cube_f) and the associated
!  qshep interpolator opbject (interp)
!
! ---------------------------------------------------------
subroutine PES_mask_interpolator_init(PESK, Lk, dim, cube_f, interp)
  FLOAT,                intent(in)   :: PESK(:,:,:)
  FLOAT,                intent(in)   :: Lk(:)
  integer,              intent(in)   :: dim
  FLOAT, pointer,       intent(inout):: cube_f(:)
  type(qshep_t),        intent(out)  :: interp
  
  integer :: np, ii, ll(MAX_DIM), ix, iy, iz
  FLOAT   :: KK(MAX_DIM)
  FLOAT, allocatable ::  kx(:),ky(:),kz(:)

  PUSH_SUB(PES_mask_interpolator_init)

  ll = 1
  do ii = 1, dim
    ll(ii) = size(PESK,ii) 
  end do 
  

  np = ll(1)*ll(2)*ll(3)  

  !check dim
  if (dim .lt. 2 .or. dim .gt. 3) then
    write (message(1),*) "This  interpolator works only for 2 <= dim <= 3" 
    call messages_fatal(1)
  end if

  
  SAFE_ALLOCATE(cube_f(np))    

  SAFE_ALLOCATE(kx(np))    
  SAFE_ALLOCATE(ky(np))    
  SAFE_ALLOCATE(kz(np))    

  cube_f = M_ZERO
  kx = M_ZERO
  ky = M_ZERO
  kz = M_ZERO


  ii=1
  do ix = 1, ll(1)
    KK(1) = Lk(ix)
    do iy = 1, ll(2)
      KK(2) = Lk(iy)
      do iz = 1, ll(3)
        KK(3) = Lk(iz)

        cube_f(ii) =  PESK(ix,iy,iz)

        kx(ii) = KK(1)
        ky(ii) = KK(2)
        kz(ii) = KK(3)

        ii = ii +1
      end do 
    end do
  end do 
  


  select case(dim)
  case (2)
    call init_qshep(interp, np, cube_f, kx, ky) 
  case (3)
    call init_qshep(interp, np, cube_f, kx, ky, kz) 
  end select
 
  SAFE_DEALLOCATE_A(kx)    
  SAFE_DEALLOCATE_A(ky)    
  SAFE_DEALLOCATE_A(kz)    
  

  POP_SUB(PES_mask_interpolator_init)
  end subroutine PES_mask_interpolator_init

! ---------------------------------------------------------
!  Destroy the interpolation objects
! ---------------------------------------------------------
subroutine PES_mask_interpolator_end(cube_f, interp)
  FLOAT, pointer,       intent(inout):: cube_f(:)
  type(qshep_t)                      :: interp


  PUSH_SUB(PES_mask_interpolator_end)
  
  call kill_qshep(interp)

  SAFE_DEALLOCATE_P(cube_f)

  POP_SUB(PES_mask_interpolator_end)
  end subroutine PES_mask_interpolator_end




! ---------------------------------------------------------
subroutine PES_mask_dump_full_map(mask, st, outp, file, dir)
  type(PES_mask_t),     intent(in) :: mask
  type(states_t),       intent(in) :: st
  character(len=*),     intent(in) :: file
  type(output_t), intent(in) :: outp
  integer,              intent(in) :: dir

  integer :: ist, ik, ii, ix, iy, iz, iunit,idim
  FLOAT ::  KK(MAX_DIM),temp, scale



  PUSH_SUB(PES_mask_dump_full_map)

  iunit = io_open(file, action='write')
  
  scale = M_ONE
  do idim=1, mask%mesh%sb%dim
    scale = scale *( mask%spacing(idim)/sqrt(M_TWO*M_PI))**2
  end do


!! INCOMPLETE !!
  do ix = 1, mask%ll(1)
    KK(1) = mask%Lk(ix)
    do iy = 1, mask%ll(2)
      KK(2) = mask%Lk(iy)
!      do iz = 1, mask%ll(3)
!        KK(3) = mask%Lk(iz)

       iz = 1
       temp = M_ZERO

       do ik = 1,st%d%nik
          do ist = 1, st%nst
             temp = temp + st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik) )**2)
           end do
       end do
      
       temp = temp * scale
         
       write(iunit, '(es18.11,2x,es18.11,2x,es18.12,i10,2x,i10,2x,es18.12,2x,i10)')& 
                    KK(1), KK(2), temp , ix , iy,  sum(KK(1:mask%mesh%sb%dim)**2) / M_TWO&
                     ,int((sum(KK(1:mask%mesh%sb%dim)**2) / M_TWO) / mask%energyStep) + 1

    end do
    write(iunit, *)  

  end do

  call io_close(iunit)


  POP_SUB(PES_mask_dump_full_map)
  end subroutine PES_mask_dump_full_map

! ---------------------------------------------------------
subroutine PES_mask_dump_full_mapM(PESK, file, Lk, dim, dir)
  FLOAT,                intent(in) :: PESK(:,:,:)
  FLOAT,                intent(in) :: Lk(:)
  character(len=*),     intent(in) :: file
  integer,              intent(in) :: dim
  integer,              intent(in) :: dir

  integer :: ist, ik, ii, ix, iy, iz, iunit,idim
  FLOAT ::  KK(MAX_DIM),temp, scale
  integer :: ll(MAX_DIM)


  PUSH_SUB(PES_mask_dump_full_mapM)

  iunit = io_open(file, action='write')

  ll = 1
  do ii = 1, dim
    ll(ii) = size(PESK,ii) 
  end do 
  


!! INCOMPLETE !!
  do ix = 1, ll(1)
    KK(1) = Lk(ix)
    do iy = 1, ll(2)
      KK(2) = Lk(iy)

      temp = M_ZERO
      do iz = 1, ll(3)
        KK(3) = Lk(iz)

        temp = temp + PESK(ix,iy,iz)
      end do  

      temp = temp * abs(lk(2)-lk(1))
        
      write(iunit, '(es18.11,2x,es18.11,2x,es18.12,i10,2x,i10,2x,es18.12,2x,i10)')& 
                    KK(1), KK(2), temp , ix , iy,  sum(KK(1:dim)**2) / M_TWO&
                     ,int((sum(KK(1:dim)**2) / M_TWO)) + 1

    end do
    write(iunit, *)  

  end do

  call io_close(iunit)


  POP_SUB(PES_mask_dump_full_mapM)
  end subroutine PES_mask_dump_full_mapM

! ---------------------------------------------------------
subroutine PES_mask_dump_power_totalM(PESK, file, Lk, dim, Emax, Estep, interpolate)
  FLOAT,                intent(in) :: PESK(:,:,:)
  FLOAT,                intent(in) :: Lk(:)
  character(len=*),     intent(in) :: file
  integer,              intent(in) :: dim
  FLOAT,                intent(in) :: Emax
  FLOAT,                intent(in) :: Estep
  logical,              intent(in) :: interpolate

  integer :: ist, ik, ii, ix, iy, iz, iunit,idim
  FLOAT ::  KK(MAX_DIM),vec

  integer :: nn
  FLOAT  :: step,DE
  FLOAT, allocatable :: npoints(:), pes(:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, pointer :: cube_f(:)
  type(qshep_t) :: interp

  FLOAT :: Dtheta, Dphi, theta, phi,EE
  integer :: np, Ntheta, Nphi, ith, iph
  integer :: ll(MAX_DIM)


  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_dump")

  PUSH_SUB(PES_mask_dump_power_totalM)

  ll = 1
  do ii = 1, dim
    ll(ii) = size(PESK,ii) 
  end do 

  step= Estep
  nn  = int(Emax/step)


  Ntheta = 36
  Dtheta = M_TWO*M_PI/Ntheta

  Nphi = 18
  Dphi = M_PI/Nphi

  SAFE_ALLOCATE(pes(1:nn))
  pes = M_ZERO

  SAFE_ALLOCATE(npoints(1:nn))
  npoints = M_ZERO

  !in 1D we do not interpolate 
  if ( (.not. interpolate) .or.  (dim .eq. 1) ) then 

    do ix = 1,ll(1)
      KK(1) = Lk(ix)
      do iy = 1, ll(2)
        KK(2) = Lk(iy)
        do iz = 1, ll(3)
          KK(3) = Lk(iz)

          if(KK(1).ne.0 .or. KK(2).ne.0 .or. KK(3).ne.0) then
            ! the power spectrum
            vec = sum(KK(1:dim)**2) / M_TWO
            ii = int(vec / step) + 1

            if(ii <= nn) then

              pes(ii) = pes(ii)+PESK(ix,iy,iz)
              npoints(ii) = npoints(ii) + M_ONE


            end if
          end if

        end do 
      end do
    end do 
   
    do ii = 2, nn
       EE = (ii-1)*step
       !Multiply for the correct Jacobian factor
       pes(ii) = pes(ii)*sqrt(M_TWO*EE)**(dim - 2) 
       pes(ii) = pes(ii) / npoints(ii)
    end do  
 

  ! Interpolate the output
  else

    call PES_mask_interpolator_init(PESK, Lk, dim, cube_f, interp)

    select case(dim)
    case(2)

      do ii = 1, nn
        EE = (ii-1)*step
        do ith = 0, Ntheta
          theta = ith*Dtheta
          KK(1) = sqrt(2*EE)*cos(theta) 
          KK(2) = sqrt(2*EE)*sin(theta)
          pes(ii) = pes(ii) + qshep_interpolate(interp, cube_f, KK(1:2))
        end do
      end do

      pes = pes * Dtheta

    case(3)

      do ii = 1, nn
        EE = (ii-1)*step
        do ith = 0, Ntheta
          theta = ith * Dtheta
          do iph = 0, Nphi
            phi = iph * Dphi
            KK(1) = sqrt(M_TWO*EE)*sin(phi)*cos(theta) 
            KK(2) = sqrt(M_TWO*EE)*sin(phi)*sin(theta)
            KK(3) = sqrt(M_TWO*EE)*cos(phi)
            pes(ii) = pes(ii) + qshep_interpolate(interp, cube_f, KK(1:3))
          end do 
        end do
        pes(ii) = pes(ii) * sqrt(M_TWO*EE)    
      end do

      pes = pes * Dtheta * Dphi

    end select

    call PES_mask_interpolator_end(cube_f, interp)

  end if 


  if (interpolate) then 
    call  PES_mask_write_power_total(file, step, pes)
  else 
    call  PES_mask_write_power_total(file, step, pes, npoints)
  end if



  SAFE_DEALLOCATE_A(pes)
  SAFE_DEALLOCATE_A(npoints)

  POP_SUB(PES_mask_dump_power_totalM)

  call profiling_out(prof)

  end subroutine PES_mask_dump_power_totalM


! ---------------------------------------------------------
subroutine PES_mask_dump_power_total(mask, st, file, wfAk)
  type(PES_mask_t),     intent(in) :: mask
  type(states_t),       intent(in) :: st
  character(len=*),     intent(in) :: file
  CMPLX, optional,      intent(in) :: wfAk(:,:,:,:,:,:)

  integer :: ist, ik, ii, ix, iy, iz, iunit,idim
  FLOAT ::  KK(MAX_DIM),vec

  integer :: nn
  FLOAT  :: step,DE
  FLOAT, allocatable :: npoints(:), pes(:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:), kx(:),ky(:),kz(:)
  FLOAT :: Dtheta, Dphi, theta, phi,EE
  type(qshep_t) :: interp
  integer :: np, Ntheta, Nphi, ith, iph


  type(profile_t), save :: prof
  call profiling_in(prof, "PESMASK_dump")

  PUSH_SUB(PES_mask_dump_power_total)


  step=mask%energyStep
  nn  = int(mask%energyMax/step)

  DE= (mask%Lk(2)-mask%Lk(1))**2/M_TWO
  
  np = mask%ll(1)*mask%ll(2)*mask%ll(3)  

  Ntheta = 36
  Dtheta = M_TWO*M_PI/Ntheta

  Nphi = 18
  Dphi = M_PI/Nphi

  SAFE_ALLOCATE(pes(1:nn))
  pes = M_ZERO

  SAFE_ALLOCATE(npoints(1:nn))
  npoints = M_ZERO

  !in 1D we do not interpolate 
  if ((.not. mask%interpolate_out) .or. (mask%mesh%sb%dim .eq. 1) ) then 

    do ix = 1, mask%ll(1)
      KK(1) = mask%Lk(ix)
      do iy = 1, mask%ll(2)
        KK(2) = mask%Lk(iy)
        do iz = 1, mask%ll(3)
          KK(3) = mask%Lk(iz)

          if(KK(1).ne.0 .or. KK(2).ne.0 .or. KK(3).ne.0) then
            ! the power spectrum
            vec = sum(KK(1:mask%mesh%sb%dim)**2) / M_TWO
            ii = int(vec / step) + 1

            if(ii <= nn) then
              do ik = 1,st%d%nik
                do ist = 1, st%nst
                  if(present(wfAk))then
                    pes(ii) = pes(ii) + st%occ(ist, ik) *& 
                             sum(abs(mask%k(ix, iy, iz, :, ist, ik) + wfAk(ix,iy,iz,:, ist, ik)  )**2)
                  else 
                    pes(ii) = pes(ii) + st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik)  )**2)
                  end if 
                end do
              end do
              npoints(ii) = npoints(ii) + M_ONE
            end if
            !Multiply for the correct Jacobian factor
!            pes(ii) = pes(ii)*sqrt(2*vec)**(mask%mesh%sb%dim - 2) 
           end if

         end do 
       end do
     end do 
  


  ! Interpolate the output
  else

    select case(mask%mesh%sb%dim)
    case(2)
      
      SAFE_ALLOCATE(cube_f(np))    
      SAFE_ALLOCATE(kx(np))    
      SAFE_ALLOCATE(ky(np))    

      cube_f = M_ZERO
      kx = M_ZERO
      ky = M_ZERO
      npoints = M_ONE

      ii=1
      do ix = 1, mask%ll(1)
        KK(1) = mask%Lk(ix)
        do iy = 1, mask%ll(2)
          KK(2) = mask%Lk(iy)
          do iz = 1, mask%ll(3)
            KK(3) = mask%Lk(iz)

            do ik = 1,st%d%nik
              do ist = 1, st%nst
                if (present(wfAk)) then
                  cube_f(ii) = cube_f(ii) + st%occ(ist, ik) *&
                               sum(abs(mask%k(ix, iy, iz, :, ist, ik) + wfAk(ix,iy,iz,:, ist, ik)  )**2)
                else
                  cube_f(ii) = cube_f(ii) + st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik)  )**2)
                end if
              end do
            end do

            kx(ii) = KK(1)
            ky(ii) = KK(2)

            ii = ii +1
          end do 
        end do
      end do 

      call init_qshep(interp, np, cube_f, kx, ky)

      do ii = 1, nn
        EE = ii*step
        do ith = 0, Ntheta
          theta = ith*Dtheta
          KK(1) = sqrt(2*EE)*cos(theta) 
          KK(2) = sqrt(2*EE)*sin(theta)
          pes(ii) = pes(ii) + qshep_interpolate(interp, cube_f, KK(1:2))
        end do
!        pes(ii) = pes(ii) * sqrt(2*EE)      
        pes(ii) = pes(ii)      
      end do

      pes = pes * Dtheta


      call kill_qshep(interp)
   
      SAFE_DEALLOCATE_A(cube_f)    
      SAFE_DEALLOCATE_A(kx)    
      SAFE_DEALLOCATE_A(ky)    

    case(3)

      SAFE_ALLOCATE(cube_f(np))    
      SAFE_ALLOCATE(kx(np))    
      SAFE_ALLOCATE(ky(np))    
      SAFE_ALLOCATE(kz(np))    

      cube_f = M_ZERO
      kx = M_ZERO
      ky = M_ZERO
      kz = M_ZERO
      npoints = M_ONE

      ii=1
      do ix = 1, mask%ll(1)
        KK(1) = mask%Lk(ix)
        do iy = 1, mask%ll(2)
          KK(2) = mask%Lk(iy)
          do iz = 1, mask%ll(3)
            KK(3) = mask%Lk(iz)

            do ik = 1,st%d%nik
              do ist = 1, st%nst
                if (present(wfAk)) then
                  cube_f(ii) = cube_f(ii) + st%occ(ist, ik) *&
                               sum(abs(mask%k(ix, iy, iz, :, ist, ik) + wfAk(ix,iy,iz,:, ist, ik)  )**2)
                else
                  cube_f(ii) = cube_f(ii) + st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik)  )**2)
                end if
              end do
            end do

            kx(ii) = KK(1)
            ky(ii) = KK(2)
            kz(ii) = KK(3)

            ii = ii +1
          end do 
        end do
      end do 

      call init_qshep(interp, np, cube_f, kx, ky, kz)

      do ii = 1, nn
        EE = ii*step
        do ith = 0, Ntheta
          theta = ith * Dtheta
          do iph = 0, Nphi
            phi = iph * Dphi
            KK(1) = sqrt(2*EE)*sin(phi)*cos(theta) 
            KK(2) = sqrt(2*EE)*sin(phi)*sin(theta)
            KK(3) = sqrt(2*EE)*cos(phi)
            pes(ii) = pes(ii) + qshep_interpolate(interp, cube_f, KK(1:3))
            pes(ii) = pes(ii) * sin(phi)
          end do 
        end do
!        pes(ii) = pes(ii) * 2*EE    
        pes(ii) = pes(ii) * sqrt(2*EE)    
      end do

      pes = pes * Dtheta * Dphi

      call kill_qshep(interp)
   
      SAFE_DEALLOCATE_A(cube_f)    
      SAFE_DEALLOCATE_A(kx)    
      SAFE_DEALLOCATE_A(ky)    
      SAFE_DEALLOCATE_A(kz)    

    end select

  end if 


  do idim=1, mask%mesh%sb%dim
    pes = pes *( mask%spacing(idim)/sqrt(M_TWO*M_PI))**2
  end do
  
  if (mask%interpolate_out) then 
    call  PES_mask_write_power_total(file, step, pes)
  else 
    call  PES_mask_write_power_total(file, step, pes, npoints)
  end if
  

  SAFE_DEALLOCATE_A(pes)
  SAFE_DEALLOCATE_A(npoints)

  POP_SUB(PES_mask_dump_power_total)

  call profiling_out(prof)

  end subroutine PES_mask_dump_power_total

! ---------------------------------------------------------
subroutine PES_mask_write_power_total(file, step, pes, npoints)
  character(len=*),     intent(in) :: file
  FLOAT,                intent(in) :: step
  FLOAT,                intent(in) :: pes(:)
  FLOAT, optional,      intent(in) :: npoints(:)

  integer :: nn, iunit, ii

  PUSH_SUB(PES_mask_write_power_total)


  nn = size(pes,1)

  iunit = io_open(file, action='write')

  !!Header
  write(iunit, '(a)') '##################################################'
  write(iunit, '(a1,a18,a18)') '#', str_center("E", 18), str_center("P(E)", 18)
  write(iunit, '(a1,a18,a18)') &
          '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18), &
          str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 18)
  write(iunit, '(a)') '##################################################'
 

   do ii = 1, nn
     if(present(npoints)) then
       if(npoints(ii) > 0) then
         write(iunit, '(es18.12,2x,es18.12,2x,es18.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), pes(ii), npoints(ii)
       end if
     else
       write(iunit, '(es18.12,2x,es18.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), pes(ii)
     end if
   end do

  call io_close(iunit)



  POP_SUB(PES_mask_write_power_total)
end subroutine PES_mask_write_power_total

! ---------------------------------------------------------
subroutine PES_mask_dump_ARPES(mask, st, file, wfAk)
  type(PES_mask_t),     intent(in) :: mask
  type(states_t),       intent(in) :: st
  character(len=*),     intent(in) :: file
  CMPLX, optional,      intent(in) :: wfAk(:,:,:,:,:,:)

  integer :: ist, ik, ii, ix, iy, iz, iunit,idim
  FLOAT ::  KK(MAX_DIM),vec

  integer :: nn
  FLOAT  :: step,DE
  FLOAT, allocatable :: npoints(:), arpes(:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:), kx(:),ky(:),kz(:)
  FLOAT :: Dtheta, Dphi, theta, phi,EE
  type(qshep_t) :: interp
  integer :: np, Ntheta, Nphi, ith, iph


  type(profile_t), save :: prof
  call profiling_in(prof, "ARPES_dump")

  PUSH_SUB(PES_mask_dump_ARPES)


  step=mask%energyStep
  nn  = 90

  DE= (mask%Lk(2)-mask%Lk(1))**2/M_TWO
  
  np = mask%ll(1)*mask%ll(2)*mask%ll(3)  

  Ntheta = 36
  Dtheta = M_TWO*M_PI/Ntheta

  Nphi = 18
  Dphi = M_PI/Nphi

  SAFE_ALLOCATE(arpes(1:nn))
  arpes = M_ZERO

  SAFE_ALLOCATE(npoints(1:nn))
  npoints = M_ZERO


    do ix = 1, mask%ll(1)
      KK(1) = mask%Lk(ix)
      do iy = 1, mask%ll(2)
        KK(2) = mask%Lk(iy)
        do iz = 1, mask%ll(3)
          KK(3) = mask%Lk(iz)

          if(KK(3)==0 .and. (KK(1).ne.0 .or. KK(2).ne.0)) then
            vec = atan2(KK(2), KK(1))
            ii  = int(abs(vec) * (nn - 1)/M_PI) + 1

            if(ii <= nn) then ! should always be true
              do ik = 1, st%d%kpt%nglobal
                do ist = 1, st%nst

                  if(present(wfAk))then
                    arpes(ii) = arpes(ii) + st%occ(ist, ik) *& 
                             sum(abs(mask%k(ix, iy, iz, :, ist, ik) + wfAk(ix,iy,iz,:, ist, ik)  )**2)
                  else 
                    arpes(ii) = arpes(ii) + st%occ(ist, ik) * sum(abs(mask%k(ix, iy, iz, :, ist, ik)  )**2)
                  end if 
                end do
              end do
            npoints(ii) = npoints(ii) + M_ONE
            end if
           end if

         end do 
       end do
     end do 
  


  do idim=1, mask%mesh%sb%dim
    arpes = arpes *( mask%spacing(idim)/sqrt(M_TWO*M_PI))**2
  end do


  iunit = io_open(file, action='write')

   do ii = 1, nn
     if(npoints(ii) > 0) then
       write(iunit, '(es18.12,2x,es18.12,2x,es18.12)')  (ii-1)*CNST(180.0)/real(nn-1, REAL_PRECISION), arpes(ii), npoints(ii)
      end if
   end do

  call io_close(iunit)

  SAFE_DEALLOCATE_A(arpes)
  SAFE_DEALLOCATE_A(npoints)

  POP_SUB(PES_mask_dump_ARPES)

  call profiling_out(prof)

  end subroutine PES_mask_dump_ARPES

  
! ---------------------------------------------------------
!
! This routine is the main routine dedicated to the output 
! of PES data
!
! ---------------------------------------------------------
subroutine PES_mask_output(mask, mesh, st,outp, file,gr, geo,iter)
  type(PES_mask_t),     intent(in)    :: mask
  type(mesh_t),         intent(in)    :: mesh
  type(states_t),       intent(in)    :: st
  character(len=*),     intent(in)    :: file
  type(output_t), intent(in)    :: outp
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  integer,              intent(in)    :: iter

  CMPLX, allocatable :: wf(:,:,:),wfAk(:,:,:,:,:,:) 
  FLOAT :: PESK(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3))
  integer :: ist, ik, ii,  iunit, idim, ierr
  character(len=100) :: fn
  character(len=256) :: dir


  PUSH_SUB(PES_mask_output)

  !FIXME: this could be dumped at the beginning of the simlulation 
  !Dump info for easy post-process
  call PES_mask_write_info(mask, tmpdir)
 

  !Photoelectron wavefunction and density in real space
  if(iand(outp%what, C_OUTPUT_PES_WFS) .ne. 0  .or.  iand(outp%what, C_OUTPUT_PES_DENSITY) .ne. 0 ) then
    write(dir, '(a,i7.7)') "td.", iter  ! name of directory
    call  PES_mask_output_states(st, gr, geo, dir, outp, mask)
  end if

  !Dump the output in the td.00iter directories
  dir = file 
  if(iand(outp%what, C_OUTPUT_PES) .ne. 0 ) then
    write(dir, '(a,i7.7,a)') "td.", iter,"/PES"  ! name of directory
  end if


  if(mask%add_psia) then 
    !The contribution of \Psi_A(x,t2) to the PES 
    SAFE_ALLOCATE(wf(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
    SAFE_ALLOCATE(wfAk(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3),1:st%d%dim,1:st%nst,1:st%d%nik))

    do ik = 1,st%d%nik
      do ist =  1, st%nst
        do idim = 1, st%d%dim
          call PES_mask_mesh_to_square(mask,mesh,st%zpsi(:, idim, ist, ik),wf,MaskHow = 3,Const = M_ONE)
          wf = (M_ONE-mask%M**10)*wf ! mask^10 is practically a box function
          call PES_mask_X_to_K(mask,mesh,wf,wfAk(:,:,:,idim, ist, ik))
        end do
      end do
    end do
  end if 

  !Create the full momentum-resolved PES matrix
  if(mask%add_psia) then 
    call PES_mask_create_full_map(mask,st,PESK,wfAk)
  else 
    call PES_mask_create_full_map(mask,st,PESK)
  end if
  
  ! Dump the full matrix in binary format for subsequent post-processing 
  write(fn, '(a,a)') trim(dir), '_map.obf'
  call io_binary_write(fn,mask%ll(1)*mask%ll(2)*mask%ll(3),PESK, ierr)
  

  ! Dump the k resolved PES on plane kz=0
  write(fn, '(a,a)') trim(dir), '_map.z=0'
!  call PES_mask_dump_full_map(mask, st, outp, fn, dir = 3)
  call PES_mask_dump_full_mapM(PESK, fn, mask%Lk, mask%mesh%sb%dim, dir = 3)


  ! Total power spectrum 
  write(fn, '(a,a)') trim(dir), '_power.sum'
!  if(mask%add_psia) then 
!    call PES_mask_dump_power_total(mask, st, fn, wfAk)
!  else 
!    call PES_mask_dump_power_total(mask, st, fn)
!  end if

   call PES_mask_dump_power_totalM(PESK,fn, mask%Lk, mask%mesh%sb%dim, mask%energyMax, mask%energyStep, mask%interpolate_out)

  if(mask%add_psia) then 
   SAFE_DEALLOCATE_A(wf)
   SAFE_DEALLOCATE_A(wfAk)
  end if

  POP_SUB(PES_mask_output)
end subroutine PES_mask_output

! ---------------------------------------------------------
! Read PES info. This is needed by the photoelectron-spectrum 
! utility. 
! ---------------------------------------------------------
subroutine PES_mask_read_info(dir, dim, Emax, Estep, ll, Lk)
  character(len=*), intent(in) :: dir
  integer,          intent(out):: dim  
  integer,          intent(out):: ll
  FLOAT,            intent(out):: Emax
  FLOAT,            intent(out):: Estep
  FLOAT, pointer,   intent(out):: Lk(:)


  character(len=256) :: filename,dummy

  integer :: iunit,ii


  PUSH_SUB(PES_mask_read_info)


  filename=trim(dir)//'td/pes'
  iunit = io_open(filename, action='read', status='old', die=.false.)

  if (iunit < 0) then
    write(message(1),'(a,a)') "Couldn't find file ",filename
    call messages_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *) dummy, dim
  read(iunit, *) dummy, Emax
  read(iunit, *) dummy, Estep
  read(iunit, *) 

  read(iunit, *) dummy, ll
   
  SAFE_DEALLOCATE_P(lk)
  SAFE_ALLOCATE(lk(1:ll))

  do ii=1, ll
    read(iunit, '(es19.12)') Lk(ii)
  end do 


  call io_close(iunit)       



  POP_SUB(PES_mask_read_info)
end subroutine PES_mask_read_info

! ---------------------------------------------------------
! Dump PES info needed for post processing with photoelecron-spectrum
! utility.
! ---------------------------------------------------------
subroutine PES_mask_write_info(mask, dir)
  type(PES_mask_t), intent(in) :: mask
  character(len=*), intent(in) :: dir

  character(len=256) :: filename

  integer :: iunit,ii
  integer :: ll(MAX_DIM)


  PUSH_SUB(PES_mask_write_info)


  filename=trim(dir)//'td/pes'

  iunit = io_open(filename, action='write')

  write(iunit, '(a,2x,i2)') 'dim', mask%mesh%sb%dim
  write(iunit, '(a,2x,es18.12)') 'Emax', mask%energyMax
  write(iunit, '(a,2x,es18.12)') 'Estep', mask%energyStep
  write(iunit, '(a)') '-------'

  write(iunit, '(a,2x,i18)') 'nK', mask%ll(1)
  do ii=1, mask%ll(1)
    write(iunit, '(es19.12)')  mask%Lk(ii)
  end do 


  call io_close(iunit)       



  POP_SUB(PES_mask_write_info)
end subroutine PES_mask_write_info


! ---------------------------------------------------------
!
! ---------------------------------------------------------
subroutine PES_mask_restart_write(mask, mesh, st)
  type(PES_mask_t), intent(in) :: mask
  type(mesh_t),     intent(in) :: mesh
  type(states_t),   intent(in) :: st

  character(len=80) :: filename, dir ,path

  integer :: itot, ik, ist, idim , np, ierr, i
  integer :: ll(MAX_DIM)


  PUSH_SUB(PES_mask_restart_write)


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
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


  ll(1:MAX_DIM) = mask%ll(1:MAX_DIM)
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

