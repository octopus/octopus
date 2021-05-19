!! Copyright (C) 2011 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!



! This needed in order to flip the sign of each Kpoit and 
! still preserve an array ordering on the kpoint mesh 
! such that the lowest index touple is associated with the smaller 
! (negative) kpoint value.
subroutine flip_sign_Lkpt_idx( dim, nk, idx)
   integer, intent(out) :: idx(:,:)
   integer, intent(in)  :: dim, nk(:)

   integer :: idx_tmp(1:maxval(nk(1:3)), 1:3)
   integer :: idir, ii

   PUSH_SUB(flip_sign_Lkpt_idx)
   
   
   idx(:,:) = 1
   do idir = 1, dim
     ! fill it with the range 1:nk
     do ii = 1, nk(idir)
       idx_tmp(ii,idir) = ii
     end do
    
     
     do ii = 1, nk(idir)
       idx(ii,idir) = nk(idir) - idx_tmp(ii,idir) + 1 
     end do
     
!          do ii = 1, nk(idir)
!            print *,idir, ii, "idx = ", idx(ii,idir),"idx_tmp =", idx_tmp(ii,idir),mod(nk(idir),2)
!          end do
   end do
   
   
   POP_SUB(flip_sign_Lkpt_idx)      
end subroutine flip_sign_Lkpt_idx

!< Generate the momentum-space mesh (p) and the arrays mapping the 
!< the mask and the kpoint meshes in p.
subroutine pes_mask_pmesh(namespace, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp)
  type(namespace_t), intent(in)    :: namespace
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             !< ll(1:dim): the dimensions of the mask-mesh
  FLOAT,             intent(in)    :: LG(:,:)           !< LG(1:maxval(ll),1:dim): the mask-mesh points  
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    !< pmesh(i1,i2,i3,1:dim): contains the positions of point
                                                        !< in the final mesh in momentum space "p" combining the 
  integer,           intent(out) :: idxZero(:)          !< The triplet identifying the zero of the coordinates           

  integer,           intent(in)  :: krng(:)             !< The range identifying the zero-weight path 
                                                        !< mask-mesh with kpoints. 
  integer,  dimension(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3), intent(inout) :: Lp  
                                                        !< maps a mask-mesh triplet of indices together with a kpoint 
                                                        !< index into a triplet on the combined momentum space mesh.

   


  integer :: ik, j1, j2, j3, nk(1:3), ip1, ip2, ip3, idir, err
  FLOAT :: kpt(1:3),GG(1:3)

  integer, allocatable :: Lkpt(:,:), idx(:,:), idx_inv(:,:), ikidx(:,:)
  FLOAT, allocatable   :: LG_(:,:)

  integer :: nkpt, kpth_dir
  FLOAT :: zero_thr


  PUSH_SUB(pes_mask_pmesh)
        
        
  nkpt = krng(2) - krng(1) + 1

  SAFE_ALLOCATE(Lkpt(krng(1):krng(2),1:3))
      
  nk(:) = 1  
  nk(1:dim) = kpoints%nik_axis(1:dim)

  Lkpt(:,:) = 1
  kpt(:) = M_ZERO
      
  zero_thr = M_EPSILON    
      
  if ( kpoints%have_zero_weight_path()) then 
    ! supporting paths only along the kx and ky directions in 
    ! reciprocal space
    kpth_dir = -1 
    if (size(pmesh, 1) > ll(1)) kpth_dir = 1
    if (size(pmesh, 2) > ll(2)) kpth_dir = 2
    ASSERT (kpth_dir /= -1 )
    
    nk(:) = 1  
    nk(kpth_dir) = nkpt
    do ik = 1 , nkpt
      Lkpt(krng(1)+ik-1,kpth_dir) = ik
      kpt(1:dim) = kpoints%get_point(krng(1) + ik -1) 
    end do
        
  else  
    
    call kpoints_grid_generate(dim, kpoints%nik_axis(1:dim), kpoints%full%nshifts, &
           kpoints%full%shifts(1:dim,1:kpoints%full%nshifts), kpoints%full%red_point,  &
           Lkpt(:,1:dim))

!       do ik = 1, kpoints_number(kpoints)
!         kpt(1:sb%dim) = kpoints_get_point(kpoints, ik, absolute_coordinates = .true.)
!         print *, ik, "Lkpt(ik)= ", Lkpt(ik,:), "-- kpt= ",kpt
!       end do


  end if

  SAFE_ALLOCATE(ikidx(maxval(nk(1:3)),1:3))
  call flip_sign_Lkpt_idx(dim, nk(:), ikidx(:,:))
  
  if (debug%info) then
    print *,"reordered"
    do ik = krng(1),krng(2)
      kpt(1:dim) = kpoints%get_point(ik, absolute_coordinates = .false.)
      print *, ik, "Lkpt(ik)= [", ikidx(Lkpt(ik,1),1), ikidx(Lkpt(ik,2),2), ikidx(Lkpt(ik,3),3),"] -- kpt= ",kpt(1)
    end do

    print *,"----"
    print *,"ll(:)", ll(:)
    print *,"----"
  end if
  
  ! We want the results to be sorted on a cube i,j,k
  ! with the first triplet associated with the smallest positions 
  SAFE_ALLOCATE(idx(1:maxval(ll(:)), 1:3))
  SAFE_ALLOCATE(idx_inv(1:maxval(ll(:)), 1:3))
  SAFE_ALLOCATE(LG_(1:maxval(ll(:)), 1:3))
  idx(:,:)=1
  idx_inv(:,:)=1
  do idir = 1, dim
    LG_(:,idir) = LG(:,idir)
    call sort(LG_(1:ll(idir), idir), idx(1:ll(idir), idir)) 
    idx_inv(:,idir) = idx(:,idir)
    call sort(idx(1:ll(idir),idir),idx_inv(1:ll(idir),idir))
  end do  
  
  if(debug%info) then
    do idir = 1, dim
      print *, "*** direction =", idir  
      do j1 = 1, ll(idir)
        print *,j1, "LG = ",LG(j1,idir),"LG_ = ",LG_(j1,idir), "idx = ", idx(j1,idir), "idx_inv = ", idx_inv(j1,idir)
      end do
      print *, "*** *** ***"  
    end do
  end if

  pmesh(:, :, :, :) = M_HUGE      
  pmesh(:, :, :, dim+1) = M_ZERO   
  Lp(:,:,:,:,:) = 0   
  err = -1
  
  ! Generate the p-space mesh and populate Lp.
  ! The grid is filled combining G-points and K-points according to the following sketch: 
  ! 
  !          x        x        x
  !
  !       o  o  o
  !       o  x  o     x        x
  !       o  o  o
  !
  !       (G = x and Kpt = o)
  ! 
  ! The grid represents the final momentum p = G + Kpt. 
  ! The lower left corner correspond to the minimum value of p and the lowest 
  ! index-touple value (ip1,ip2,ip3) = (1,1,1). 
  do ik = krng(1),krng(2)
    kpt(1:dim) = kpoints%get_point(ik) 
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 

          GG(1:3)= (/LG_(j1,1),LG_(j2,2),LG_(j3,3)/)
          
          ip1 = (j1 - 1) * nk(1) + ikidx(Lkpt(ik,1), 1)
          ip2 = (j2 - 1) * nk(2) + ikidx(Lkpt(ik,2), 2)
          ip3 = (j3 - 1) * nk(3) + ikidx(Lkpt(ik,3), 3)
          
          Lp(idx_inv(j1,1),idx_inv(j2,2),idx_inv(j3,3),ik,1:3) =  (/ip1,ip2,ip3/)
          
          ! The final momentum corresponds to p = G - K. 
          ! This is due to the convention for the sign of the Bloch phase associated to 
          ! each kpoint subspace in hamilonian_update wich is exp(-i K*r) instead of the 
          ! most commonly used (in textbooks) exp(i K*r). 
          ! Thus in order to solve the propagation equations only for the periodic part 
          ! of each KS orbital u_K(r) I nees to project on plane waves with <p|=<G-K|. 
          ! Note: in order to preserve the correct ordering of the final mesh in p 
          ! we remap each kpoint grid-index with ikidx(:,1:3).
          pmesh(ip1, ip2, ip3, 1:dim) = GG(1:dim) - kpt(1:dim)
          pmesh(ip1, ip2, ip3, dim+1) = pmesh(ip1, ip2, ip3, dim+1) + 1 
          
!               print *,idx_inv(j1,1),idx_inv(j2,2),idx_inv(j3,3),ik,"  Lp(i1,i2,i3,ik,1:dim) = ",  (/ip1,ip2,ip3/)
!               print *, "pmesh = ",pmesh(ip1, ip2, ip3, :) !,"  GG = ",  GG (1:dim), "  kpt = ", kpt(1:dim)



          ! Sanity checks
          if (sum(pmesh(ip1, ip2, ip3, 1:dim)**2)<=zero_thr) then
            err = err + 1 
            !Find the indices identifying the center of the coordinates 
            idxZero(1:3) = (/ip1,ip2,ip3/)
          end if
          
          if (pmesh(ip1, ip2, ip3, dim+1) > 1 ) then
            err = -2 
          end if


        end do 
      end do 
!           print *,idx_inv(j1,1),ik,"  Lp(i1,ll(2),ll(3),ik,1) = ", ip1, "pmesh = ",pmesh(ip1, ip2, ip3, 1)

    end do 
  end do
  
!       do ip1 = 1, ll(1) * nk(1)
!         print *,ip1, "Pmesh", pmesh(ip1, 1, 1, 1)
!       end do

  if (kpoints%have_zero_weight_path()) then 
  ! With a path we just need to get the correct the zero index on the in-plane direction  
  ! perpendicular to the path since is along this direction that we are going 
  ! to slice with pes_mask_output_full_mapM_cut. Since on this direction we only 
  ! have G points I simply need to look for the zero index of the G-grid.
  ! Note that the G-grid must always include the (0,0,0) point. 
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 

          GG(1:3)= (/LG_(j1,1),LG_(j2,2),LG_(j3,3)/)
          if (sum(GG(1:3)**2)<=M_EPSILON) idxZero(1:3) = (/j1,j2,j3/)
        
        end do
      end do
    end do
    
  else   
    
    if (err == -1) then
      call messages_write('Illformed momentum-space mesh: could not find p = 0 coordinate.')
      call messages_fatal(namespace=namespace)
    end if 

    if (err > 1) then
      call messages_write('Illformed momentum-space mesh: more than one point with p = 0 coordinate.')
      call messages_write('This can happen only if the kpoint mesh does not contain gamma.')
      call messages_warning(namespace=namespace)
    end if 

  end if

  if(debug%info) then
    print * ,"idxZero(1:3)=", idxZero(1:3)
  end if


  if (err == -2) then
    call messages_write('Illformed momentum-space mesh: two or more points with the same p.')
    call messages_fatal(namespace=namespace)
  end if 
  
 
 

  SAFE_DEALLOCATE_A(Lkpt)
  SAFE_DEALLOCATE_A(idx)
  SAFE_DEALLOCATE_A(idx_inv)
  SAFE_DEALLOCATE_A(ikidx)      
  SAFE_DEALLOCATE_A(LG_)
  
  
  POP_SUB(pes_mask_pmesh)
end subroutine pes_mask_pmesh


!< Build the photoemission map form the restart files
subroutine pes_mask_map_from_states(restart, st, ll, pesK, krng, Lp, istin)
  type(restart_t),    intent(in) :: restart
  type(states_elec_t),intent(in) :: st
  integer,            intent(in) :: ll(:)
  FLOAT, target,     intent(out) :: pesK(:,:,:,:)
  integer,           intent(in)  :: krng(:) 
  integer,  dimension(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3), intent(in) :: Lp
  integer, optional, intent(in)  :: istin 
  
  integer :: ik, ist, idim, itot, nkpt, ispin
  integer :: i1, i2, i3, ip(1:3)
  integer :: idone, ntodo
  CMPLX   :: psiG1(1:ll(1),1:ll(2),1:ll(3)), psiG2(1:ll(1),1:ll(2),1:ll(3))
  FLOAT   :: weight 
  integer :: istart, iend, nst

  PUSH_SUB(pes_mask_map_from_states)

  istart = 1
  iend = st%nst
  nst = st%nst
  if(present(istin)) then
    istart = istin
    iend = istin
    nst = 1
  end if

  nkpt =  krng(2)-krng(1)+1
!       ntodo = st%d%kpt%nglobal * st%nst * st%d%dim
  ntodo = nkpt * nst 
  idone = 0 
  call loct_progress_bar(-1, ntodo)
  
  pesK = M_ZERO
  do ik = krng(1), krng(2)
    ispin = st%d%get_spin_index(ik)
    
    do ist = istart, iend

      if (st%d%kweights(ik) < M_EPSILON) then
        ! we have a zero-weight path
        weight = st%occ(ist, ik)!/nkpt
      else
        weight = st%occ(ist, ik) * st%d%kweights(ik)
      end if
      
      if(st%d%ispin /= SPINORS) then 

        do idim = 1, st%d%dim
          itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
          call pes_mask_map_from_state(restart, itot, ll, psiG1)
        
          do i1 = 1, ll(1)
            do i2 = 1, ll(2)
              do i3 = 1, ll(3)
                ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
              
                  pesK(ip(1),ip(2),ip(3), ispin) = pesK(ip(1),ip(2),ip(3), ispin) &
                                                 + abs(psiG1(i1,i2,i3))**2 * weight 
                
              end do
            end do
          end do
                
        end do
      else ! SPINORS
        idim = 1
        itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
        call pes_mask_map_from_state(restart, itot, ll, psiG1)
        idim = 2
        itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
        call pes_mask_map_from_state(restart, itot, ll, psiG2)
            
        do i1 = 1, ll(1)
          do i2 = 1, ll(2)
            do i3 = 1, ll(3)
              ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
            
                pesK(ip(1),ip(2),ip(3), 1) = pesK(ip(1),ip(2),ip(3), 1) &
                                               + abs(psiG1(i1,i2,i3))**2 * weight 

                pesK(ip(1),ip(2),ip(3), 2) = pesK(ip(1),ip(2),ip(3), 2) &
                                               + abs(psiG2(i1,i2,i3))**2 * weight

                pesK(ip(1),ip(2),ip(3), 3) = pesK(ip(1),ip(2),ip(3), 3) &
                                               + TOFLOAT(psiG1(i1,i2,i3)*conjg(psiG2(i1,i2,i3))) * weight
                                               
                pesK(ip(1),ip(2),ip(3), 4) = pesK(ip(1),ip(2),ip(3), 4) &
                                               + aimag(psiG1(i1,i2,i3)*conjg(psiG2(i1,i2,i3))) * weight
            end do
          end do
        end do
          
          
      end if
      
      idone = idone +1 
      call loct_progress_bar(idone, ntodo)
      
    end do
  end do

  write(stdout, '(1x)')

  POP_SUB(pes_mask_map_from_states)
end subroutine pes_mask_map_from_states


subroutine pes_mask_map_from_state(restart, idx, ll, psiG)
  type(restart_t),  intent(in)  :: restart
  integer,          intent(in)  :: idx
  integer,          intent(in)  :: ll(:)
  CMPLX, target,    intent(out) :: psiG(:,:,:)

  character(len=80) :: filename
  integer ::  np, err

  PUSH_SUB(pes_mask_map_from_state)

  psiG = M_Z0
  np = product(ll(:))
  
  write(filename,'(a,i10.10)') 'pes_', idx

  call zrestart_read_binary(restart, filename, np, psiG(:,:,:), err)

  POP_SUB(pes_mask_map_from_state)  
end subroutine pes_mask_map_from_state



! ---------------------------------------------------------
!> Write the photoelectron wavefunctions in real space
! ---------------------------------------------------------
subroutine pes_mask_output_states(namespace, space, st, gr, ions, dir, outp, mask)
  type(namespace_t),     intent(in)    :: namespace
  type(space_t),         intent(in)    :: space
  type(states_elec_t),   intent(in)    :: st
  type(grid_t),          intent(in)    :: gr
  type(ions_t),          intent(in)    :: ions
  character(len=*),      intent(in)    :: dir
  type(output_t),        intent(in)    :: outp
  type(pes_mask_t),      intent(inout) :: mask

  integer :: ik, ist, idim, is, ierr
  character(len=80) :: fname
  type(unit_t) :: fn_unit

  CMPLX, allocatable :: PsiAB(:,:,:,:), psi(:)
  FLOAT,allocatable :: RhoAB(:,:) 
  type(cube_function_t) :: cf
  type(mesh_t):: mesh   
  
  type(wfs_elec_t)        :: psib
  type(density_calc_t) :: dens_calc

  PUSH_SUB(pes_mask_output_states)
  
  mesh= gr%mesh

  SAFE_ALLOCATE(PsiAB(1:mesh%np_part, 1:st%d%dim,1:st%nst,1:st%d%nik))
  SAFE_ALLOCATE(RhoAB(1:mesh%np_part, 1:st%d%nspin))
  SAFE_ALLOCATE(psi(1:mesh%np_part))

  call zcube_function_alloc_RS(mask%cube, cf, force_alloc = .true.)

  RhoAB= M_ZERO
  
  !Calculate the pes density \Psi_A + \Psi_B on the simulation Box  
  call density_calc_init(dens_calc, st, gr, RhoAB)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        cf%zRs=M_z0

        call pes_mask_K_to_X(mask, mask%k(:,:,:, idim, ist, ik),cf%zRs)

        call pes_mask_cube_to_mesh(mask, cf, PsiAB(:, idim, ist, ik))        

        call states_elec_get_state(st, gr%mesh, idim, ist, ik, psi)

        if (mask%mode /= PES_MASK_MODE_PASSIVE) then 
          PsiAB(:, idim, ist, ik) = PsiAB(:, idim, ist, ik) + psi(:)
        end if
        
      end do
    end do
     
    call wfs_elec_init(psib, st%d%dim, st%st_start, st%st_end, PsiAB(:, :, st%st_start:, ik), ik)
    call density_calc_accumulate(dens_calc, psib)
    call psib%end()

  end do
  
  call density_calc_end(dens_calc)

  ! THE OUTPUT 
  if (outp%what(OPTION__OUTPUT__PES_DENSITY)) then
    fn_unit = units_out%length**(-space%dim)
    do is = 1, st%d%nspin
      if(st%d%nspin == 1) then
        write(fname, '(a)') 'pes_den'
      else
        write(fname, '(a,i1)') 'pes_den-sp', is
      end if
      call dio_function_output(outp%how(OPTION__OUTPUT__PES_DENSITY), dir, fname, namespace, space, gr%mesh, &
        RhoAB(:, is), fn_unit, ierr, ions = ions, grp = st%dom_st_kpt_mpi_grp)
    end do
  end if


  if (outp%what(OPTION__OUTPUT__PES_WFS)) then
    fn_unit = sqrt(units_out%length**(-space%dim))
    do ist = st%st_start, st%st_end
!        if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            if(st%d%nik > 1) then
              if(st%d%dim > 1) then
                write(fname, '(a,i3.3,a,i4.4,a,i1)') 'pes_wf-k', ik, '-st', ist, '-sd', idim
              else
                write(fname, '(a,i3.3,a,i4.4)')      'pes_wf-k', ik, '-st', ist
              end if
            else
              if(st%d%dim > 1) then
                write(fname, '(a,i4.4,a,i1)')        'pes_wf-st', ist, '-sd', idim
              else
                write(fname, '(a,i4.4)')             'pes_wf-st', ist
              end if
            end if
              
            call zio_function_output(outp%how(OPTION__OUTPUT__PES_WFS), dir, fname, namespace, space, gr%mesh, &
              PsiAB(1:, idim, ist, ik), fn_unit, ierr, ions = ions)

          end do
        end do
 !       end if
    end do
  end if

  SAFE_DEALLOCATE_A(PsiAB)
  SAFE_DEALLOCATE_A(RhoAB)

  call zcube_function_free_RS(mask%cube, cf)


  POP_SUB(pes_mask_output_states)
end subroutine pes_mask_output_states


! ---------------------------------------------------------
!
!> Calculates the momentum-resolved photoelectron probability
!!\f[
!!            P(k) = \sum_i |\Psi_{B,i}(k)|^2 
!!\f]
! ---------------------------------------------------------
subroutine pes_mask_fullmap(mask, space, st, ik, pesK, wfAk)
  type(pes_mask_t),    intent(in)  :: mask
  type(space_t),       intent(in)  :: space
  type(states_elec_t), intent(in)  :: st
  integer,             intent(in)  :: ik  
  FLOAT, target,       intent(out) :: pesK(:,:,:)
  CMPLX, optional,     intent(in)  :: wfAk(:,:,:,:,:,:)

  integer :: ist, kx, ky, kz, idim
  FLOAT   :: scale
  FLOAT, pointer :: pesKloc(:,:,:)


  PUSH_SUB(pes_mask_fullmap)

  pesK = M_ZERO
  if (mask%cube%parallel_in_domains) then
    SAFE_ALLOCATE(pesKloc(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3)))
    pesKloc = M_ZERO
  else 
    pesKloc => pesK
  end if

  do ist = st%st_start, st%st_end

    do kx = 1, mask%ll(1)
      do ky = 1, mask%ll(2)
        do kz = 1, mask%ll(3)

          if(present(wfAk))then
            pesKloc(kx, ky, kz) = pesKloc(kx, ky, kz) + st%occ(ist, ik) * &
              sum(abs(mask%k(kx, ky, kz, :, ist, ik) + wfAk(kx,ky,kz,:, ist, ik)  )**2)
          else
            pesKloc(kx, ky, kz) = pesKloc(kx, ky, kz) + st%occ(ist, ik) * sum(abs(mask%k(kx, ky, kz, :, ist, ik)  )**2)
          end if
      
        end do
      end do
    end do

  end do


!   if(st%parallel_in_states .or. st%d%kpt%parallel) then
!     call comm_allreduce(st%st_kpt_mpi_grp, pesKloc)
!   end if

  if(st%parallel_in_states) then
    call comm_allreduce(st%mpi_grp, pesKloc)
  end if

  
  if (mask%cube%parallel_in_domains) then

    call dcube_function_allgather(mask%cube, pesK, pesKloc, gatherfs = .true.)
    
!     if(mask%pw_map_how .eq. PW_MAP_PFFT) then
!       call dcube_function_allgather(mask%cube, pesK, pesKloc, & 
!                                     order = (/2,3,1/))
! 
!     else if(mask%pw_map_how .eq. PW_MAP_PNFFT) then
!       
!       
! !       pesK(mask%fs_istart(1):mask%fs_istart(1)+mask%ll(1)-1, &
! !            mask%fs_istart(2):mask%fs_istart(2)+mask%ll(2)-1, &
! !            mask%fs_istart(3):mask%fs_istart(3)+mask%ll(3)-1) = &
! !         pesKloc(1:mask%ll(1),1:mask%ll(2),1:mask%ll(3))
! !       
! !       write(file, '(i3.3)') mpi_world%rank
! !       file = "fullmap_"//trim(file)
! !       print *,"write file: ",file
! !       call pes_mask_dump_full_mapM(pesK, file, mask%Lk)
! 
! 
!       
!       call dcube_function_allgather(mask%cube, pesK, pesKloc, gatherfs = .true.)
! !                                     order = (/1,3,2/))
! !     else
! !       call dcube_function_allgather(mask%cube, pesK, pesKloc)
!       
!     end if
    
!     call dcube_function_allgather(mask%cube, pesK, pesKloc, gatherfs = .true.)

    SAFE_DEALLOCATE_P(pesKloc) 
  end if

  ! This is needed in order to normalize the Fourier integral 
  ! since along the periodi dimensions the discrete Fourier transform is 
  ! exact we need to renormalize only along the non periodic ones
  scale = M_ONE
  do idim = space%periodic_dim + 1, space%dim
    scale = scale *( mask%spacing(idim)/sqrt(M_TWO*M_PI))**2
  end do
  pesK = pesK *scale 

  POP_SUB(pes_mask_fullmap)
end subroutine pes_mask_fullmap

! --------------------------------------------------------
!
!>  Qshep interpolation helper function initialization.
!!  Generates the linearized version of pesK (cube_f) and the associated
!!  qshep interpolator opbject (interp).
!
! ---------------------------------------------------------
subroutine pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp, pmesh)
  type(namespace_t),  intent(in)    :: namespace
  FLOAT,              intent(in)    :: pesK(:,:,:)
  FLOAT,              intent(in)    :: Lk(:,:)
  integer,            intent(in)    :: ll(:)
  integer,            intent(in)    :: dim
  FLOAT, allocatable, intent(out)   :: cube_f(:)
  type(qshep_t),      intent(out)   :: interp
  FLOAT, optional,    intent(in)    :: pmesh(:,:,:,:)
  
  integer :: np, ii, ix, iy, iz
  FLOAT   :: KK(3)
  FLOAT, allocatable ::  kx(:),ky(:),kz(:)

  PUSH_SUB(pes_mask_interpolator_init)


  call messages_write("Initializing Qshep interpolator. Be patient it may take a while... ")
  call messages_info()  
  
  np = ll(1)*ll(2)*ll(3)  

  !check dim
  if (dim  <  2 .or. dim > 3) then
    message(1) = "This interpolator works only for 2 <= dim <= 3." 
    call messages_fatal(1, namespace=namespace)
  end if
  
  SAFE_ALLOCATE(cube_f(1:np))

  SAFE_ALLOCATE(kx(1:np))
  SAFE_ALLOCATE(ky(1:np))
  SAFE_ALLOCATE(kz(1:np))

  cube_f = M_ZERO
  kx = M_ZERO
  ky = M_ZERO
  kz = M_ZERO


  ii=1
  if (present(pmesh)) then
    do ix = 1, ll(1)
      do iy = 1, ll(2)
        do iz = 1, ll(3)

          cube_f(ii) =  pesK(ix,iy,iz)

          kx(ii) = pmesh(ix,iy,iz,1)
          ky(ii) = pmesh(ix,iy,iz,2)
          kz(ii) = pmesh(ix,iy,iz,3)

          ii = ii +1
        end do 
      end do
    end do 
    
  else
    
    do ix = 1, ll(1)
      KK(1) = Lk(ix, 1)
      do iy = 1, ll(2)
        KK(2) = Lk(iy, 2)
        do iz = 1, ll(3)
          KK(3) = Lk(iz, 3)

          cube_f(ii) =  pesK(ix,iy,iz)

          kx(ii) = KK(1)
          ky(ii) = KK(2)
          kz(ii) = KK(3)

          ii = ii +1
        end do 
      end do
    end do 
  end if  


  select case(dim)
    case (2)
      call qshep_init(interp, np, cube_f, kx, ky) 
    case (3)
      call qshep_init(interp, np, cube_f, kx, ky, kz) 
  end select
 
  SAFE_DEALLOCATE_A(kx)    
  SAFE_DEALLOCATE_A(ky)    
  SAFE_DEALLOCATE_A(kz)    
 
  call messages_write("done")
  call messages_new_line()
  call messages_info()
  

  POP_SUB(pes_mask_interpolator_init)
end subroutine pes_mask_interpolator_init

! ---------------------------------------------------------
!>  Destroy the interpolation objects
! ---------------------------------------------------------
subroutine pes_mask_interpolator_end(cube_f, interp)
  FLOAT, allocatable, intent(inout) :: cube_f(:)
  type(qshep_t),      intent(inout) :: interp

  PUSH_SUB(pes_mask_interpolator_end)
  
  call qshep_end(interp)
  
  SAFE_DEALLOCATE_A(cube_f)
  
  POP_SUB(pes_mask_interpolator_end)
end subroutine pes_mask_interpolator_end



! ---------------------------------------------------------
subroutine pes_mask_output_full_mapM(pesK, file, namespace, space, Lk, ll, how, pmesh)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  type(space_t),     intent(in) :: space
  FLOAT,             intent(in) :: Lk(:,:)
  integer,           intent(in) :: ll(:)  
  integer,           intent(in) :: how
  FLOAT, optional,   intent(in) :: pmesh(:,:,:,:)  
  
  integer :: iunit
  integer :: ierr
  character(len=512) :: filename
  type(cube_t) :: cube
  type(cube_function_t) :: cf
  integer :: ii
  FLOAT :: dk(3)  

  PUSH_SUB(pes_mask_output_full_mapM)

  call cube_init(cube, ll, namespace, space)
  call dcube_function_alloc_RS(cube, cf, force_alloc = .true.)
  cf%dRS = pesK
  
  
  if (.not. present(pmesh) ) then
    ! Ignore Lk and use pmesh
    dk(:) = M_ZERO
    dk(1:space%dim) = abs(Lk(2,1:space%dim)-Lk(1,1:space%dim))
    do ii = 1, space%dim
      dk(ii) = units_from_atomic(sqrt(units_out%energy), dk(ii))
    end do
  end if
  
#if defined(HAVE_NETCDF)  
  
  if(bitand(how, OPTION__OUTPUTFORMAT__NETCDF) /= 0) then
    filename = trim(file)//".ncdf"
    write(message(1), '(a)') 'Writing netcdf format file: '
    call messages_info(1)
  
    call dout_cf_netcdf(filename, ierr, cf, cube, space, dk(:), .false., sqrt(units_out%energy)**space%dim)

  end if

#endif
  
  if(bitand(how, OPTION__OUTPUTFORMAT__VTK) /= 0)  then
    filename = trim(file)//".vtk"
    write(message(1), '(a)') 'Writing vtk format file: '
    call messages_info(1)
    
    if (present(pmesh)) then          
      call dvtk_out_cf_structured(filename, namespace, 'PES_mapM', ierr, cf, cube,& 
        sqrt(units_out%energy)**space%dim, pmesh, ascii = .false.)
    else 
      call dvtk_out_cf(filename, namespace, 'PES_mapM', ierr, cf, cube, dk(:),& 
        sqrt(units_out%energy)**space%dim)
    end if        
      
   else
     write(message(1), '(a)') 'Writing ASCII format file: '
     call messages_info(1)
     call out_ascii()
  end if
  
  call cube_end(cube)
  call dcube_function_free_RS(cube, cf)

  POP_SUB(pes_mask_output_full_mapM)
  
contains  

! ---------------------------------------------------------
! Just dump the matrix in ascii form x,y,z, val. As the output
! can be pretty big is probably better to use netcdf instead.
! ---------------------------------------------------------
  subroutine out_ascii()

    integer :: ii, ix, iy, iz
    FLOAT ::  KK(3)
    integer :: ll(3)

    PUSH_SUB(pes_mask_output_full_mapM.out_ascii)

    iunit = io_open(file, namespace, action='write')


    ll = 1
    do ii = 1, 3
      ll(ii) = size(pesK,ii) 
    end do

    do ix = 1, ll(1)
      KK(1) = Lk(ix,1)
      do iy = 1, ll(2)
        KK(2) = Lk(iy,2)
        do iz = 1, ll(3)
          KK(3) = Lk(iz,3)
 
            write(iunit, '(es19.12,2x,es19.12,2x,es19.12,2x,es19.12)') &
                    units_from_atomic(sqrt(units_out%energy), KK(1)),&
                    units_from_atomic(sqrt(units_out%energy), KK(2)),&
                    units_from_atomic(sqrt(units_out%energy), KK(3)),&
                    pesK(ix,iy,iz) 
 
          end do  
        end do      
      end do

      call io_close(iunit)
  
      POP_SUB(pes_mask_output_full_mapM.out_ascii)
  end subroutine out_ascii
  
end subroutine pes_mask_output_full_mapM



! ---------------------------------------------------------
subroutine pes_mask_output_full_mapM_cut(pesK, file, namespace, ll, dim, pol, dir, integrate, pos, Lk, pmesh)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: pol(3)
  integer,           intent(in) :: dir
  integer,           intent(in) :: integrate
  integer, optional, intent(in) :: pos(3)
  FLOAT, optional,   intent(in) :: Lk(:,:)
  FLOAT, optional,   intent(in) :: pmesh(:,:,:,:)

  integer              :: ii, ix, iy, iunit,ldir(2), icut(3)
  FLOAT                :: KK(3),temp
  integer, allocatable :: idx(:,:)
  FLOAT, allocatable   :: Lk_(:,:)
  FLOAT                :: rotation(1:dim,1:dim)
  logical              :: aligned_axis
! integration
  FLOAT                :: K, KKK(3), theta, phi, Dphi, Dk(3)
  integer              :: iph, Nphi
! progress
  integer              :: idone, ntodo

  FLOAT, allocatable :: cube_f(:)
  type(qshep_t) :: interp

  PUSH_SUB(pes_mask_output_full_mapM_cut)
  
  iunit = io_open(file, namespace, action='write')


  
  ASSERT(size(pesK, 1) == ll(1))

  if (.not. present(pmesh)) then
    ASSERT(present(Lk))
    
    SAFE_ALLOCATE(idx(1:maxval(ll(:)), 1:3))
    SAFE_ALLOCATE(Lk_(1:maxval(ll(:)), 1:3))

    Dk(:) = M_ZERO
    Dk(1:dim) = abs(Lk(2,1:dim)-Lk(1,1:dim))
  
    do ii = 1, 3
      Lk_(:,ii) = Lk(:,ii)
      call sort(Lk_(1:ll(ii), ii), idx(1:ll(ii), ii)) !We need to sort the k-vectors in order to dump in gnuplot format
    end do  
  end if
  
  aligned_axis = sum((pol-(/0 ,0 ,1/))**2)  <= M_EPSILON  .or. &
                 sum((pol-(/0 ,1 ,0/))**2)  <= M_EPSILON  .or. &
                 sum((pol-(/1 ,0 ,0/))**2)  <= M_EPSILON  
  
  if (present(pos)) then
    icut(1:3) = pos(1:3) 
  else
    icut(1:3) = ll(1:3)/2 + 1
  end if
  
  if (aligned_axis .and. integrate == INTEGRATE_NONE) then !no need to rotate and interpolate
    
    select case (dir)
      case (1)
        ldir(:) =(/2,3/)
      case (2)
        ldir(:) =(/1,3/)
      case (3)
        ldir(:) =(/1,2/)
        
    end select
    
    
    if (present(pmesh)) then
      do ix = 1, ll(ldir(1))
        do iy = 1, ll(ldir(2))
      
          select case (dir)
            case (1)
              temp = pesK(icut(dir), ix, iy)    
              KK(1) = pmesh(icut(dir), ix, iy, ldir(1))
              KK(2) = pmesh(icut(dir), ix, iy, ldir(2))
            case (2)
              temp = pesK(ix, icut(dir), iy)    
              KK(1) = pmesh(ix, icut(dir), iy, ldir(1))
              KK(2) = pmesh(ix, icut(dir), iy, ldir(2))
            case (3)
              temp = pesK(ix, iy, icut(dir))        
              KK(1) = pmesh(ix, iy, icut(dir), ldir(1))
              KK(2) = pmesh(ix, iy, icut(dir), ldir(2))
          end select
        
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                  units_from_atomic(sqrt(units_out%energy), KK(1)),&
                  units_from_atomic(sqrt(units_out%energy), KK(2)),&
                  temp
        end do
        write(iunit, *)  
      end do

    else 
      do ix = 1, ll(ldir(1))
        KK(1) = Lk_(ix, ldir(1))
        do iy = 1, ll(ldir(2))
          KK(2) = Lk_(iy, ldir(2))
      
          select case (dir)
            case (1)
              temp = pesK(idx(icut(dir), 1), idx(ix, 2), idx(iy, 3))    
            case (2)
              temp = pesK(idx(ix, 1), idx(icut(dir), 2), idx(iy, 3))    
            case (3)
              temp = pesK(idx(ix, 1), idx(iy, 2), idx(icut(dir), 3))        
          end select
        
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                  units_from_atomic(sqrt(units_out%energy), KK(1)),&
                  units_from_atomic(sqrt(units_out%energy), KK(2)),&
                  temp
        end do
        write(iunit, *)  
      end do
    end if
    
  else 
    ! We set the z-axis along the pol vector 
    call generate_rotation_matrix(rotation, (/M_ZERO, M_ZERO, M_ONE/), pol )

    if(debug%info) then
      print *,"Rotate z-axis over the zenith axis"
      print *,rotation(1,:)
      print *,rotation(2,:)
      print *,rotation(3,:)
    end if

    call pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp, pmesh)

    ntodo = product(ll(1:2))
    idone = 0 
    call loct_progress_bar(-1, ntodo)

    do ix = 1, ll(1)
     do iy = 1, ll(2)

       !cut 
       select case (dir)
         case (1)
           KK(1) = M_ZERO 
           if (present(pmesh)) then
             KK(2) = pmesh(ix, 1, 1, 1)
             KK(3) = pmesh(1, iy, 1, 2)
           else       
             KK(2) = Lk_(ix, 1)
             KK(3) = Lk_(iy, 2)
           end if
         case (2)
           KK(2) = M_ZERO
           if(present(pmesh)) then
             KK(1) = pmesh(ix,  1, 1, 1)
             KK(3) = pmesh(1,  iy, 1, 2)
           else
             KK(1) = Lk_(ix, 1)
             KK(3) = Lk_(iy, 2)
           end if
 
         case (3)
           KK(3) = M_ZERO
           if(present(pmesh)) then
             KK(1) = pmesh(ix,  1, 1, 1)
             KK(2) = pmesh(1,  iy, 1, 2)
           else 
             KK(1) = Lk_(ix, 1)
             KK(2) = Lk_(iy, 2)
           end if

       end select

       temp = qshep_interpolate(interp, cube_f, matmul(rotation,KK(1:3)) )

       select case (integrate)
         case (INTEGRATE_PHI)
           temp = M_ZERO
           K = sqrt(KK(1)**2 + KK(2)**2 + KK(3)**2)

           Nphi = 360
           Dphi = M_TWO * M_PI/Nphi

           do iph = 0, Nphi
             phi = iph * Dphi
             theta = atan2(sqrt(KK(1)**2+KK(2)**2),KK(3))

             KKK(1) = K *sin(theta)*cos(phi) 
             KKK(2) = K *sin(theta)*sin(phi)
             KKK(3) = K *cos(theta)
        
             temp = temp + &
                           abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
           end do
           temp = temp * Dphi     
           
         case (INTEGRATE_KX)
           temp = M_ZERO
           do ii =1, ll(1)
              KKK(:) = KK(:) + (/Lk(ii,1), M_ZERO, M_ZERO/)
              temp = temp + &
                            abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
           end do
           temp = temp * Dk(1)  
         
         case (INTEGRATE_KY)
           temp = M_ZERO
           do ii =1, ll(2)
              KKK(:) = KK(:) + (/M_ZERO, Lk(ii,2), M_ZERO/)
              temp = temp + &
                            abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
           end do
           temp = temp * Dk(2)  
         
         case (INTEGRATE_KZ)
           temp = M_ZERO
           do ii =1, ll(3)
              KKK(:) = KK(:) + (/M_ZERO, M_ZERO, Lk(ii,3)/)
              temp = temp + &
                            abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
           end do
           temp = temp * Dk(3)  

       end select
   
       write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
               units_from_atomic(sqrt(units_out%energy), KK(1)),&
               units_from_atomic(sqrt(units_out%energy), KK( 2)),&
               temp

  
       idone = idone +1 
       call loct_progress_bar(idone, ntodo)
       
     end do
     write(iunit, *)  
    end do

    write(stdout, '(1x)')


  end if

  call io_close(iunit)
  
  SAFE_DEALLOCATE_A(idx) 
  SAFE_DEALLOCATE_A(Lk_) 
  
  POP_SUB(pes_mask_output_full_mapM_cut)
end subroutine pes_mask_output_full_mapM_cut

! ---------------------------------------------------------
subroutine pes_mask_output_ar_polar_M(pesK, file, namespace, Lk, ll, dim, dir, Emax, Estep)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  FLOAT,             intent(in) :: Lk(:,:)
  integer,           intent(in) :: ll(:)
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: Emax
  FLOAT,             intent(in) :: Estep
  FLOAT,             intent(in) :: dir(:) 

  FLOAT ::  KK(3)

  integer :: nn, ie
  FLOAT  :: step
  FLOAT, allocatable ::  pesM(:,:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:)
  type(qshep_t) :: interp

  FLOAT :: Dtheta, Dphi, theta, phi, EE
  integer :: Nphi, ith, iph, Nth
  FLOAT :: vref(1:dim), rotation(1:dim,1:dim)
  FLOAT :: eGrid(3), thGrid(3), phiBounds(2)
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(pesMask_ar_polar)))

  PUSH_SUB(pes_mask_output_ar_polar_M)

  ! we do the calculation assuming the polarization along the x-axis 
  vref = M_ZERO
  vref(1) = M_ONE
  !the rotation matrix from the x-axis to the actual polarization direction
  call generate_rotation_matrix(rotation, vref, dir)
  
  step= Estep
  nn  = int(Emax/step)

  Nth = 300 
  Dtheta = M_PI/Nth

  Nphi = 360
  Dphi = M_TWO * M_PI/Nphi

  SAFE_ALLOCATE(pesM(1:Nth,1:nn))
  pesM = M_ZERO


  !in 1D we do not interpolate 
  if (  (dim  ==  1) ) then 
    message(1)="Impossible to obtain angle-dependent quantities in 1D."
    call messages_fatal(1, namespace=namespace)

  else

    call pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp)

    select case(dim)
    case(2)

  


    case(3)

      do ith = 1, Nth
        theta = (ith - 1) * Dtheta 
        do ie = 1, nn
          EE = (ie -1) * step

          do iph = 0, Nphi
            phi = iph * Dphi

            KK(2) = sqrt(M_TWO*EE)*sin(theta)*cos(phi) 
            KK(3) = sqrt(M_TWO*EE)*sin(theta)*sin(phi)
            KK(1) = sqrt(M_TWO*EE)*cos(theta)
            
            !sometimes the interpolator gives negative values therefore the abs()  
            pesM(ith,ie) = pesM(ith,ie) + &
                          abs(qshep_interpolate(interp, cube_f, matmul(rotation,KK(1:3)) ))
          end do
          pesM(ith,ie) = pesM(ith,ie) * sqrt(M_TWO*EE)     

        end do
      end do

      pesM = pesM * Dphi

    end select

    call pes_mask_interpolator_end(cube_f, interp)

  end if


  eGrid(1)= M_ZERO
  eGrid(2)= Emax
  eGrid(3)= step

  thGrid(1)= M_ZERO
  thGrid(2)= M_PI
  thGrid(3)= Dtheta

  phiBounds(1) = M_ZERO
  phiBounds(2) = M_TWO * M_PI

  call  pes_mask_write_2D_map(file, namespace, pesM, 2, thGrid, eGrid, dir, phiBounds)

  SAFE_DEALLOCATE_A(pesM)

  POP_SUB(pes_mask_output_ar_polar_M)

  call profiling_out(prof)

end subroutine pes_mask_output_ar_polar_M


! ---------------------------------------------------------
subroutine pes_mask_output_ar_plane_M(pesK, file, namespace, Lk, ll, dim, dir, Emax, Estep)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  FLOAT,             intent(in) :: Lk(:,:)
  integer,           intent(in) :: ll(:)
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: Emax
  FLOAT,             intent(in) :: Estep
  FLOAT,             intent(in) :: dir(:) 

  integer :: ix, iy
  FLOAT ::  KK(3)

  integer :: nn, nx, ny
  FLOAT  :: step, eGrid(3)
  FLOAT, allocatable ::  pesM(:,:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:)
  type(qshep_t) :: interp

  FLOAT :: Dphi, theta, phi, Ex, Ey, EE, phiBounds(2)
  integer :: Nphi, iph
  FLOAT :: vref(1:dim), rotation(1:dim,1:dim)
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(PESMask_ar_plane)))

  PUSH_SUB(pes_mask_output_ar_plane_M)

  ! we do the calculation assuming the polarization along the x-axis 
  vref = M_ZERO
  vref(1) = M_ONE
  !the rotation matrix from the x-axis to the actual polarization direction
  call generate_rotation_matrix(rotation, vref, dir)
  
  step= Estep
  nn  = int(Emax/step)

  nx = 2*nn
  ny = nn  


  Nphi = 360
  Dphi = M_TWO * M_PI/Nphi

  SAFE_ALLOCATE(pesM(1:2*nn,1:nn))
  pesM = M_ZERO


  !in 1D we do not interpolate 
  if (  (dim  ==  1) ) then 
    message(1)="Impossible to obtain angle-dependent quantities in 1D."
    call messages_fatal(1, namespace=namespace)

  else

    call pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp)

    select case(dim)
    case(2)

      do ix = 1, nx
        Ex = (ix -nx/2 -1)*step 
        do iy = 1, ny
          Ey = (iy-1)*step
          EE =sqrt(Ex**2 + Ey**2)
          theta = atan2(Ey,Ex)

          KK(1) = sqrt(2*EE)*cos(theta) 
          KK(2) = sqrt(2*EE)*sin(theta)

          pesM(ix,iy) = pesM(ix,iy) + qshep_interpolate(interp, cube_f, KK(1:2))

        end do
      end do


    case(3)

      do ix = 1, nx
        Ex = (ix -nx/2 -1)*step 
        do iy = 1, ny
          Ey = (iy-1)*step
          EE =sqrt(Ex**2 + Ey**2)

          do iph = 0, Nphi
            phi = iph * Dphi
            theta = atan2(Ey,Ex)

            KK(2) = sqrt(M_TWO*EE)*sin(theta)*cos(phi) 
            KK(3) = sqrt(M_TWO*EE)*sin(theta)*sin(phi)
            KK(1) = sqrt(M_TWO*EE)*cos(theta)
            
            !sometimes the interpolator gives negative values therefore the abs()  
            pesM(ix,iy) = pesM(ix,iy) + &
                          abs(qshep_interpolate(interp, cube_f, matmul(rotation,KK(1:3)) ))
          end do
          pesM(ix,iy) = pesM(ix,iy) * sqrt(M_TWO*EE)     

        end do
      end do

      pesM = pesM * Dphi

    end select

    call pes_mask_interpolator_end(cube_f, interp)

  end if


  eGrid(1)= M_ZERO
  eGrid(2)= Emax
  eGrid(3)= step
  
  phiBounds(1) = M_ZERO
  phiBounds(2) = M_TWO * M_PI
  
  call pes_mask_write_2D_map(file, namespace, pesM, 1, eGrid, eGrid, dir, phiBounds)

  SAFE_DEALLOCATE_A(pesM)

  POP_SUB(pes_mask_output_ar_plane_M)

  call profiling_out(prof)

end subroutine pes_mask_output_ar_plane_M

! ---------------------------------------------------------
subroutine pes_mask_output_ar_spherical_cut_M(pesK, file, namespace, Lk, ll, dim, dir, Emin, Emax, Estep)
  FLOAT,             intent(in) :: pesK(:,:,:)
  type(namespace_t), intent(in) :: namespace
  character(len=*),  intent(in) :: file
  FLOAT,             intent(in) :: Lk(:,:)
  integer,           intent(in) :: ll(:)
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: Emin
  FLOAT,             intent(in) :: Emax
  FLOAT,             intent(in) :: Estep
  FLOAT,             intent(in) :: dir(:) 

  FLOAT ::  KK(3)

  integer :: nn, ie
  FLOAT  :: step
  FLOAT, allocatable ::  pesM(:,:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:)
  type(qshep_t) :: interp

  FLOAT :: Dtheta, Dphi, theta, phi, EE
  integer :: Nphi, ith, iph, Nth
  FLOAT :: vref(1:dim), rotation(1:dim,1:dim)
  FLOAT :: phGrid(3), thGrid(3), eBounds(2)
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(PESMask_ar_spherical_cut)))

  PUSH_SUB(pes_mask_output_ar_spherical_cut_M)

  ! by default the polarization is along the z-axis 
  vref = M_ZERO
  vref(3) = M_ONE
  !the rotation matrix from the x-axis to the actual polarization direction
  call generate_rotation_matrix(rotation, vref, dir)
  
  step= Estep
  nn  = int(abs(Emax-Emin)/step)

  Nth = 300 
  Dtheta = M_PI/Nth

  Nphi = 360
  Dphi = M_TWO * M_PI/Nphi

  SAFE_ALLOCATE(pesM(1:Nphi,1:Nth))
  pesM = M_ZERO


  !in 1D we do not interpolate 
  if (  (dim  ==  1) ) then 
    message(1)="Impossible to obtain angle-dependent quantities in 1D."
    call messages_fatal(1, namespace=namespace)

  else

    call pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp)

    select case(dim)
    case(2)
  


    case(3)

      do iph = 1, Nphi
        phi = (iph - 1) * Dphi

        do ith = 1, Nth
          theta = (ith - 1) * Dtheta 

          do ie = 0, nn
            EE = Emin + ie * step

            KK(1) = sqrt(M_TWO*EE)*sin(theta)*cos(phi) 
            KK(2) = sqrt(M_TWO*EE)*sin(theta)*sin(phi)
            KK(3) = sqrt(M_TWO*EE)*cos(theta)
            
            !sometimes the interpolator gives negative values therefore the abs()  
            pesM(iph,ith) = pesM(iph,ith) + &
                          abs(qshep_interpolate(interp, cube_f, matmul(rotation,KK(1:3)) ))
          end do
          pesM(iph,ith) = pesM(iph,ith) * sqrt(M_TWO*EE)     

        end do
      end do

      pesM = pesM * Dphi

    end select

    call pes_mask_interpolator_end(cube_f, interp)

  end if


  phGrid(1)= M_ZERO
  phGrid(2)= M_TWO * M_PI
  phGrid(3)= Dphi

  thGrid(1)= M_ZERO
  thGrid(2)= M_PI
  thGrid(3)= Dtheta

  eBounds(1) = Emin
  eBounds(2) = Emax

  call  pes_mask_write_2D_map(file, namespace, pesM, 4, phGrid, thGrid, dir, eBounds)

  SAFE_DEALLOCATE_A(pesM)

  POP_SUB(pes_mask_output_ar_spherical_cut_M)

  call profiling_out(prof)

end subroutine pes_mask_output_ar_spherical_cut_M


! ========================================================================
!>  Common interface to write 2D maps in gnuplot with header files for 
!!  different objects. The modes are:
!!
!!  - 1 Angle- and energy-resolved on cartesian coordinates
!!  - 2 Angle- and energy-resolved in polar coordinates
!!  - 3 Velocity map on a plane
!!  - 4 Spherical cut
!
! ========================================================================
subroutine pes_mask_write_2D_map(file, namespace, pesM, mode, xGrid, yGrid, vv, intSpan)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  FLOAT,             intent(in) :: pesM(:,:)
  integer,           intent(in) :: mode
  FLOAT,             intent(in) :: xGrid(:)   !< max min and step for the x axis
  FLOAT,             intent(in) :: yGrid(:)   !< max min and step for the y axis
  FLOAT,             intent(in) :: vv(:)      !< for mode=1,2 indicate the Zenith axis for mode 3 the cutting plane
  FLOAT, optional,   intent(in) :: intSpan(:) !< for integrated quantities indicate the integral region    

  integer :: nx,ny, iunit, ix,iy

  PUSH_SUB(pes_mask_write_2D_map)

  nx = size(pesM,1)
  ny = size(pesM,2)

  iunit = io_open(file, namespace, action='write')

  select case (mode)
    case(1)
      !!Angle Energy Cartesian
        write(iunit, '(a)') '##################################################'

        write(iunit, '(a)') '#'        
        write(iunit, '(a1,a20,a1,f10.2,a2,f10.2,a2,f10.2,a1)')&
                                              "#"," Zenith axis: ","(",&
                                              vv(1),", ",vv(2),", ",vv(3),")"
        write(iunit, '(a1,a20,a1,f10.2,a2,f10.2,a1)')&
                                              "#"," Integrated in phi: ","(",&
                                              intSpan(1),", ",intSpan(2),")"
        write(iunit, '(a)') '#'        


        write(iunit, '(a1,a19,2x,a19,2x,a19)') '#', str_center("E*cos(th)", 19),&
                                              str_center("E*sin(th)", 19),&
                                              str_center("P(E)", 19)
        write(iunit, '(a1,a19,2x,a19,2x,a19)') &
          '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19), &
               str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19), & 
               str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 19)
        write(iunit, '(a)') '##################################################'
    case (2)
    !!Angle Energy polar
      write(iunit, '(a)') '##################################################'

      write(iunit, '(a)') '#'        
      write(iunit, '(a1,a20,a1,f10.2,a2,f10.2,a2,f10.2,a1)')&
                                            "#"," Zenith axis: ","(",&
                                            vv(1),", ",vv(2),", ",vv(3),")"
      write(iunit, '(a1,a20,a1,f10.2,a2,f10.2,a1)')&
                                            "#"," Integrated in phi: ","(",&
                                            intSpan(1),", ",intSpan(2),")"
      write(iunit, '(a)') '#'        
      
      
      write(iunit, '(a1,a19,2x,a19,2x,a19)') '#', str_center("Theta", 19),&
                                            str_center("E", 19),&
                                            str_center("P(E)", 19)
      write(iunit, '(a1,a19,2x,a19,2x,a19)') &
        '#', str_center('[rad]', 19), &
             str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19), & 
             str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 19)
      write(iunit, '(a)') '##################################################'

    case (3)
    !!Velocity map
      write(iunit, '(a)') '##################################################'
      write(iunit, '(a1,a19,2x,a19,2x,a19)') '#', str_center("p1", 19),&
                                            str_center("p2", 19),&
                                            str_center("P(p)", 19)
      write(iunit, '(a1,a19,2x,a19,2x,a19)') &
        '#', str_center('[sqrt('//trim(units_abbrev(units_out%energy)) // ')]', 19), & 
             str_center('[sqrt('//trim(units_abbrev(units_out%energy)) // ')]', 19), & 
             str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 19)
      write(iunit, '(a)') '##################################################'

      case (4)
      !!Spherical Cut
        write(iunit, '(a)') '##################################################'

        write(iunit, '(a)') '#'        
        write(iunit, '(a1,a20,a1,f10.2,a2,f10.2,a2,f10.2,a1)')&
                                              "#"," Zenith axis: ","(",&
                                              vv(1),", ",vv(2),", ",vv(3),")"
        write(iunit, '(a1,a20,a1,es19.12,a2,es19.12,a1,a19)')&
                                              "#"," Integrated in E: ","(",&
                                              intSpan(1),", ",intSpan(2),")",&
                      str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19)
        write(iunit, '(a)') '#'        
      
      
        write(iunit, '(a1,a19,2x,a19,2x,a19)') '#', str_center("Phi", 19),&
                                              str_center("Theta", 19),&
                                              str_center("P(Phi,Theta)", 19)
        write(iunit, '(a1,a19,2x,a19,2x,a19)') &
          '#', str_center('[rad]', 19), &
               str_center('[rad]', 19), & 
               str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 19)
        write(iunit, '(a)') '##################################################'

    
  end select

  do ix = 1, nx
    do iy = 1, ny

      select case (mode)
        case (1)
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                  units_from_atomic(units_out%energy, (ix -nx/2 - 1) * xGrid(3) + xGrid(1) ),&
                  units_from_atomic(units_out%energy, (iy - 1) * yGrid(3) + yGrid(1) ), pesM(ix,iy)
        case(2)
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &   
                  (ix - 1) * xGrid(3) + xGrid(1), &
                  units_from_atomic(units_out%energy, (iy - 1) * yGrid(3) + yGrid(1) ), pesM(ix,iy)
        case(3)
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &   
                  (ix - 1) * xGrid(3) + xGrid(1), &
                  (iy - 1) * yGrid(3) + yGrid(1), pesM(ix,iy)

        case(4)
          write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &   
                  (ix - 1) * xGrid(3) + xGrid(1), &
                  (iy - 1) * yGrid(3) + yGrid(1), pesM(ix,iy)
    
      end select

    end do
    write(iunit,*) 
  end do
  
  call io_close(iunit)


  POP_SUB(pes_mask_write_2D_map)
end subroutine pes_mask_write_2D_map




! ---------------------------------------------------------
subroutine pes_mask_output_power_totalM(pesK, file, namespace, Lk, ll, dim, Emax, Estep, interpolate)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  FLOAT,             intent(in) :: Lk(:,:)
  integer,           intent(in) :: ll(:)
  integer,           intent(in) :: dim
  FLOAT,             intent(in) :: Emax
  FLOAT,             intent(in) :: Estep
  logical,           intent(in) :: interpolate

  integer :: ii, ix, iy, iz
  FLOAT ::  KK(3),vec

  integer :: nn
  FLOAT  :: step
  FLOAT, allocatable :: npoints(:), pes(:)

  ! needed for interpolation in 2D and 3D 
  FLOAT, allocatable :: cube_f(:)
  type(qshep_t) :: interp

  FLOAT :: Dtheta, Dphi, theta, phi,EE
  integer :: Ntheta, Nphi, ith, iph

  PUSH_SUB(pes_mask_output_power_totalM)


  step= Estep
  nn  = int(Emax/step)


  Ntheta = 360
  Dtheta = M_TWO*M_PI/Ntheta

  Nphi = 180
  Dphi = M_PI/Nphi

  SAFE_ALLOCATE(pes(1:nn))
  pes = M_ZERO

  SAFE_ALLOCATE(npoints(1:nn))
  npoints = M_ZERO

  !in 1D we do not interpolate 
  if ( (.not. interpolate) .or.  (dim  ==  1) ) then 

    do ix = 1,ll(1)
      KK(1) = Lk(ix, 1)
      do iy = 1, ll(2)
        KK(2) = Lk(iy, 2)
        do iz = 1, ll(3)
          KK(3) = Lk(iz, 3)

          if(KK(1) /= 0 .or. KK(2) /= 0 .or. KK(3) /= 0) then
            ! the power spectrum
            vec = sum(KK(1:dim)**2) / M_TWO
            ii = int(vec / step) + 1

            if(ii <= nn) then

              pes(ii) = pes(ii)+pesK(ix,iy,iz)
              npoints(ii) = npoints(ii) + M_ONE

            end if
          end if

        end do
      end do
    end do

    do ii = 2, nn
      ! npoints==0.0 when pes==0.0
      if(pes(ii)/= M_ZERO)then
        EE = (ii-1)*step
        !Multiply for the correct Jacobian factor
        pes(ii) = pes(ii)*sqrt(M_TWO*EE)**(dim - 2) 
        pes(ii) = pes(ii) / npoints(ii)
      end if
    end do


    ! Interpolate the output
  else

    call pes_mask_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp)

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

            KK(2) = sqrt(M_TWO*EE)*sin(phi)*cos(theta) 
            KK(3) = sqrt(M_TWO*EE)*sin(phi)*sin(theta)
            KK(1) = sqrt(M_TWO*EE)*cos(phi)

            pes(ii) = pes(ii) + qshep_interpolate(interp, cube_f, KK(1:3))
          end do
        end do
        pes(ii) = pes(ii) * sqrt(M_TWO*EE)    
      end do

      pes = pes * Dtheta * Dphi

    end select

    call pes_mask_interpolator_end(cube_f, interp)

  end if


  if (interpolate) then 
    call pes_mask_write_power_total(file, namespace, step, pes)
  else 
    call pes_mask_write_power_total(file, namespace, step, pes, npoints)
  end if

  SAFE_DEALLOCATE_A(pes)
  SAFE_DEALLOCATE_A(npoints)

  POP_SUB(pes_mask_output_power_totalM)


end subroutine pes_mask_output_power_totalM


! ---------------------------------------------------------
subroutine pes_mask_write_power_total(file, namespace, step, pes, npoints)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  FLOAT,             intent(in) :: step
  FLOAT,             intent(in) :: pes(:)
  FLOAT, optional,   intent(in) :: npoints(:)

  integer :: nn, iunit, ii

  PUSH_SUB(pes_mask_write_power_total)

  nn = size(pes,1)

  iunit = io_open(file, namespace, action='write')

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
        write(iunit, '(es19.12,2x,es19.12,2x,es19.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), pes(ii), npoints(ii)
      end if
    else
      write(iunit, '(es19.12,2x,es19.12)')  units_from_atomic(units_out%energy, (ii - 1) * step), pes(ii)
    end if
  end do
  
  call io_close(iunit)


  POP_SUB(pes_mask_write_power_total)
end subroutine pes_mask_write_power_total


  
! ---------------------------------------------------------
!
!> This routine is the main routine dedicated to the output 
!! of PES data
!
! ---------------------------------------------------------
subroutine pes_mask_output(mask, mesh, st, outp, namespace, space, file, gr, ions, iter)
  type(pes_mask_t),    intent(inout)    :: mask
  type(mesh_t),        intent(in)       :: mesh
  type(states_elec_t), intent(in)       :: st
  type(output_t),      intent(in)       :: outp
  type(namespace_t),   intent(in)       :: namespace
  type(space_t),       intent(in)       :: space
  character(len=*),    intent(in)       :: file
  type(grid_t),        intent(in)       :: gr
  type(ions_t),        intent(in)       :: ions
  integer,             intent(in)       :: iter

  CMPLX, allocatable :: wfAk(:,:,:,:,:,:), psi(:)
  FLOAT :: pesK(1:mask%fs_n_global(1),1:mask%fs_n_global(2),1:mask%fs_n_global(3)),pol(3)
  integer :: ist, ik, idim, ierr, st1, st2, k1, k2
  character(len=100) :: fn
  character(len=256) :: dir
  type(cube_function_t) :: cf1

  type(profile_t), save :: prof

  PUSH_SUB(pes_mask_output)
  call profiling_in(prof, TOSTRING(X(PESMASK_out)))
  
  !Output info for easy post-process
  if(mpi_grp_is_root(mpi_world)) call pes_mask_write_info(mask, "td.general", namespace)
 

  !Photoelectron wavefunction and density in real space
  if (outp%what(OPTION__OUTPUT__PES_WFS) .or. outp%what(OPTION__OUTPUT__PES_DENSITY)) then
    write(dir, '(a,i7.7)') "td.", iter  ! name of directory
    call  pes_mask_output_states(namespace, space, st, gr, ions, dir, outp, mask)
  end if
  
  if (space%is_periodic()) then
    ! For periodic systems the results must be obtained using
    ! the oct-photoelectron-spectrum routine
    call profiling_out(prof)
    POP_SUB(pes_mask_output)
    return
  end if

  !Write the output in the td.00iter directories
  dir = file 
  if (outp%what(OPTION__OUTPUT__PES)) then
    write(dir, '(a,i7.7,a)') "td.", iter,"/PESM"  ! name of directory
  end if

 !The contribution of \Psi_A(x,t2) to the PES 
  if(mask%add_psia) then 
    st1 = st%st_start
    st2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end
    SAFE_ALLOCATE(wfAk(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3),1:st%d%dim,st1:st2,k1:k2))
    wfAk = M_z0
  
    call zcube_function_alloc_RS(mask%cube, cf1, force_alloc = .true.) 
    call cube_function_alloc_FS(mask%cube, cf1, force_alloc = .true.) 

    SAFE_ALLOCATE(psi(1:mesh%np_part))

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist =  st%st_start, st%st_end
        do idim = 1, st%d%dim
          call states_elec_get_state(st, mesh, idim, ist, ik, psi)
          call pes_mask_mesh_to_cube(mask, psi, cf1)
          cf1%zRs = (M_ONE - mask%cM%zRs**10) * cf1%zRs ! mask^10 is practically a box function
          call pes_mask_X_to_K(mask, cf1%zRs, cf1%Fs)
          wfAk(:,:,:,idim, ist, ik) = cf1%Fs
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi)

    call zcube_function_free_RS(mask%cube, cf1)
    call  cube_function_free_FS(mask%cube, cf1)
  end if 

  !Create the full momentum-resolved PES matrix
  pesK = M_ZERO
  ! This loop over kpoints should not be reached since for periodic systems 
  ! this part is skipped. I keep this here for the moment just for the debug purposes.
  do ik = st%d%kpt%start, st%d%kpt%end

    if(mask%add_psia) then 
      call pes_mask_fullmap(mask, space, st, ik, pesK, wfAk)
    else 
      call pes_mask_fullmap(mask, space, st, ik, pesK)
    end if
    
    ! only the root node of the domain and state parallelization group writes the output
    if(mpi_grp_is_root(st%dom_st_mpi_grp)) then 
      ! Output the full matrix in binary format for subsequent post-processing 
      if(st%d%nik == 1) then
      write(fn, '(a,a)') trim(dir), '_map.obf'
      call io_binary_write(io_workpath(fn, namespace), &
        mask%fs_n_global(1)*mask%fs_n_global(2)*mask%fs_n_global(3), pesK, ierr)
                           
         
      ! Total power spectrum 
      write(fn, '(a,a)') trim(dir), '_power.sum'
      call pes_mask_output_power_totalM(pesK,fn, namespace, mask%Lk, mask%ll, space%dim, & 
                                       mask%energyMax, mask%energyStep, .false.)
      end if

      ! Output the p resolved PES on plane pz=0
      if(st%d%nik > 1) then
        write(fn, '(a,a,i3.3,a)') trim(dir), '_map-k',ik ,'.pz=0'
      else
        write(fn, '(a,a)') trim(dir), '_map.pz=0'
      end if
      pol = (/M_ZERO, M_ZERO, M_ONE/)
      call pes_mask_output_full_mapM_cut(pesK, fn, namespace, mask%ll, space%dim, &
        pol = pol, dir = 3, integrate = INTEGRATE_NONE, Lk = mask%Lk)
                                     
    end if
  end do

  if(mask%add_psia) then 
    SAFE_DEALLOCATE_A(wfAk)
  end if

  call profiling_out(prof)
  
  POP_SUB(pes_mask_output)
end subroutine pes_mask_output

! ---------------------------------------------------------
!> Read pes info.
! ---------------------------------------------------------
subroutine pes_mask_read_info(dir, namespace, dim, Emax, Estep, ll, Lk,RR)
  character(len=*),   intent(in)  :: dir
  type(namespace_t),  intent(in)  :: namespace
  integer,            intent(out) :: dim  
  FLOAT,              intent(out) :: Emax
  FLOAT,              intent(out) :: Estep
  integer,            intent(out) :: ll(:)
  FLOAT, allocatable, intent(out) :: Lk(:,:)
  FLOAT, allocatable, intent(out) :: RR(:)


  character(len=256) :: filename, dummy
  integer :: iunit, ii, idim, ierr


  PUSH_SUB(pes_mask_read_info)


  filename = trim(dir)//'pes'
  iunit = io_open(filename, namespace, action='read', status='old')

  SAFE_ALLOCATE(RR(1:2))

  rewind(iunit)
  read(iunit, *) dummy, dim
  read(iunit,'(a10,2x,es19.12)') dummy, RR(1)
  read(iunit,'(a10,2x,es19.12)') dummy, RR(2)
  read(iunit, *) dummy, Emax
  read(iunit, *) dummy, Estep
  read(iunit, *) 

  read(iunit, '(4x,i18)',  advance='no', iostat = ierr) ll(1)
  idim=2
  do while( (ierr == 0) .and. (idim <= dim) )
    read(iunit, '(2x,i18)',  advance='no', iostat = ierr) ll(idim)
    idim=idim+1
  end do
  read(iunit, *) 
   
  SAFE_ALLOCATE(Lk(1:maxval(ll(:)), 1:3))
  Lk = M_ZERO

  do ii = 1, maxval(ll(:))
    do idim = 1, dim
      read(iunit, '(2x,es19.12)', advance='no') Lk(ii, idim)
    end do
    read(iunit, *)
  end do 
  

  call io_close(iunit)       

  POP_SUB(pes_mask_read_info)
end subroutine pes_mask_read_info


! ---------------------------------------------------------
!> Output pes info
! ---------------------------------------------------------
subroutine pes_mask_write_info(mask, dir, namespace)
  type(pes_mask_t),  intent(in) :: mask
  character(len=*),  intent(in) :: dir
  type(namespace_t), intent(in) :: namespace

  character(len=256) :: filename

  integer :: iunit,ii,idim

  PUSH_SUB(pes_mask_write_info)

  filename = trim(dir)//'/pes'

  iunit = io_open(filename, namespace, action='write')

  write(iunit, '(a10,2x,i2)') 'dim', mask%mesh%sb%dim
  write(iunit, '(a10,2x,es19.12)') 'Mask R1', mask%mask_R(1)
  write(iunit, '(a10,2x,es19.12)') 'Mask R2', mask%mask_R(2)
  write(iunit, '(a10,2x,es19.12)') 'Emax', mask%energyMax
  write(iunit, '(a10,2x,es19.12)') 'Estep', mask%energyStep
  write(iunit, '(a)') '-------'

  write(iunit, '(a)', advance='no') 'nK'
  do idim = 1, mask%mesh%sb%dim
    write(iunit, '(2x,i18)', advance='no') mask%ll(idim)
  end do 
  write(iunit, '(1x)')
     
  do ii = 1, maxval(mask%ll(:))
    do idim = 1, mask%mesh%sb%dim
      write(iunit, '(2x,es19.12)', advance='no')  mask%Lk(ii, idim)
    end do
    write(iunit, '(1x)')
  end do 


  call io_close(iunit)       

  POP_SUB(pes_mask_write_info)
end subroutine pes_mask_write_info


! ---------------------------------------------------------
!
! ---------------------------------------------------------
subroutine pes_mask_dump(mask, namespace, restart, st, ierr)
  type(pes_mask_t), target, intent(in)  :: mask
  type(namespace_t),        intent(in)  :: namespace
  type(restart_t),          intent(in)  :: restart
  type(states_elec_t),      intent(in)  :: st
  integer,                  intent(out) :: ierr

  character(len=80) :: filename, path, lines(2)
  integer :: itot, ik, ist, idim, ll(3), np, iunit, err, err2
  CMPLX, pointer :: gwf(:,:,:)
  type(profile_t), save :: prof
   
  PUSH_SUB(pes_mask_dump)

  ierr = 0

  if (debug%info) then
    message(1) = "Debug: Writing PES mask restart."
    call messages_info(1)
  end if

  call profiling_in(prof, TOSTRING(X(PESMASK_dump)))

  iunit = restart_open(restart, 'pes_mask')
  write(lines(1), '(a10,2x,es19.12)') 'Mask R1', mask%mask_r(1)
  write(lines(2), '(a10,2x,es19.12)') 'Mask R2', mask%mask_r(2)
  call restart_write(restart, iunit, lines, 2, err)
  if (err /= 0) ierr = ierr + 1
  call restart_close(restart, iunit)

  ll(1:3) = mask%fs_n_global(1:3)
  np = ll(1)*ll(2)*ll(3) 

  err2 = 0
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim

        itot = ist + (ik - 1)*st%nst +  (idim - 1)*st%nst*st%d%kpt%nglobal

        write(filename,'(i10.10)') itot

        !FIXME: the following should not use io_binary_write directly. Instead, the
        !task of writing the function to the file should be done by the restart module.

        path = trim(restart_dir(restart))//'/pes_'//trim(filename)//'.obf'

        if (mask%cube%parallel_in_domains) then
          SAFE_ALLOCATE(gwf(1:ll(1),1:ll(2),1:ll(3))) 
          !call zcube_function_allgather(mask%cube, gwf, wf, transpose = .true.)
          call zcube_function_allgather(mask%cube, gwf, mask%k(:,:,:, idim, ist, ik), gatherfs = .true.)
          
          !if(mask%pw_map_how .eq. PW_MAP_PFFT) then
          !  call zcube_function_allgather(mask%cube, gwf, wf, order = (/2,3,1/))
          ! 
          !else if(mask%pw_map_how .eq. PW_MAP_PNFFT) then
          !  call zcube_function_allgather(mask%cube, gwf, wf, order = (/2,3,1/))
          !                                      
          !end if
          
        else
          gwf => mask%k(:,:,:, idim, ist, ik)
        end if
        
        !we need to check both cube and mesh mpi_grp to make sure only the root node writes the
        !output and preserve state parallelization 
        if (mpi_grp_is_root(mask%cube%mpi_grp) .and. mpi_grp_is_root(mask%mesh%mpi_grp)) then 
          call io_binary_write(path, np, gwf(:,:,:), err)
          if (err /= 0) then
            err2 = err2 + 1
            message(1) = "Unable to write PES mask restart data to '"//trim(path)//"'."
            call messages_warning(1, namespace=namespace)
          end if
        end if

        if (mask%cube%parallel_in_domains) then
          SAFE_DEALLOCATE_P(gwf)
        end if

      end do
    end do
  end do
  if (err2 /= 0) ierr = ierr + 2

  if (debug%info) then
    message(1) = "Debug: Writing PES mask restart done."
    call messages_info(1)
  end if

  call profiling_out(prof)

  POP_SUB(pes_mask_dump)
end subroutine pes_mask_dump

! ---------------------------------------------------------
subroutine pes_mask_load(mask, namespace, restart, st, ierr)
  type(pes_mask_t),    intent(inout) :: mask
  type(namespace_t),   intent(in)    :: namespace
  type(restart_t),     intent(in)    :: restart
  type(states_elec_t), intent(inout) :: st
  integer,             intent(out)   :: ierr

  character(len=80) :: filename
  integer :: itot, ik, ist, idim , np, err, err2, iunit, ll(3)
  character(len=128) :: lines(2)
  character(len=7) :: dummy
  FLOAT, allocatable :: rr(:)

  PUSH_SUB(pes_mask_load)

  ierr = 0

  if (restart_skip(restart)) then
    ierr = -1
    POP_SUB(pes_mask_load)
    return
  end if

  if (debug%info) then
    message(1) = "Debug: Reading PES mask restart."
    call messages_info(1)
  end if


  SAFE_ALLOCATE(rr(1:2))
  iunit = restart_open(restart, 'pes_mask')
  call restart_read(restart, iunit, lines, 2, err)
  if (err /= 0) then    
    ierr = ierr + 1
  else
    read(lines(1),'(a10,2x,es19.12)') dummy, rr(1)
    read(lines(2),'(a10,2x,es19.12)') dummy, rr(2)
  end if
  call restart_close(restart, iunit)


  ll(1:3) = mask%ll(1:3)
  np =ll(1)*ll(2)*ll(3) 

  err2 = 0
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim
        itot = ist + (ik-1) * st%nst+  (idim-1) * st%nst*st%d%kpt%nglobal

        write(filename,'(a,i10.10)') 'pes_', itot

        call zrestart_read_binary(restart, filename, np, mask%k(:,:,:, idim, ist, ik), err)
        if (err /= 0) err2 = err2 + 1

      end do
    end do
  end do
  if (err2 /= 0) ierr = ierr + 2

  if (ierr == 0) then
    if (rr(1) /= mask%mask_r(1) .or. rr(2) /= mask%mask_r(2)) then
      message(1) = "PhotoElectronSpectrum = pes_mask : The mask parameters have changed."
      message(2) = "I will restart mapping from the previous context."
      call messages_warning(2, namespace=namespace)
      call pes_mask_restart_map(mask, namespace, st, rr)
    end if
  end if

  SAFE_DEALLOCATE_A(rr)

  if (debug%info) then
    message(1) = "Debug: Reading PES mask restart done."
    call messages_info(1)
  end if

  POP_SUB(pes_mask_load)
end subroutine pes_mask_load


! ---------------------------------------------------------
subroutine pes_mask_restart_map(mask, namespace, st, RR)
  type(pes_mask_t),    intent(inout) :: mask
  type(namespace_t),   intent(in)    :: namespace
  type(states_elec_t), intent(inout) :: st
  FLOAT,               intent(in)    :: RR(2)

  integer :: ik, ist, idim
  CMPLX, allocatable :: psi(:)
  FLOAT, allocatable :: M_old(:,:,:)
  type(cube_function_t):: cf1,cf2

  PUSH_SUB(pes_mask_restart_map)

  SAFE_ALLOCATE(M_old(1:mask%ll(1), 1:mask%ll(2), 1:mask%ll(3)))
  call zcube_function_alloc_RS(mask%cube, cf1, force_alloc = .true.) 
  call  cube_function_alloc_FS(mask%cube, cf1, force_alloc = .true.) 
  call zcube_function_alloc_RS(mask%cube, cf2, force_alloc = .true.)

  SAFE_ALLOCATE(psi(1:mask%mesh%np))

  call pes_mask_generate_mask_function(mask, namespace, mask%mesh, mask%shape, RR, M_old)
  
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      do idim = 1, st%d%dim
        cf1%zRs = M_z0
        call pes_mask_K_to_X(mask, mask%k(:,:,:, idim, ist, ik), cf1%zRs)
        call states_elec_get_state(st, mask%mesh, idim, ist, ik, psi)
        call pes_mask_mesh_to_cube(mask, psi, cf2)
        cf2%zRs = cf1%zRs + cf2%zRs ! the whole pes orbital in real space 
        cf1%zRs = cf2%zRs* mask%cM%zRs !modify the orbital in A
        call pes_mask_cube_to_mesh(mask, cf1, psi)
        call states_elec_set_state(st, mask%mesh, idim, ist, ik, psi)
        cf2%zRs = cf2%zRs * (mask%cM%zRs-M_old) ! modify the k-orbital in B 
        call pes_mask_X_to_K(mask, cf2%zRs, cf1%Fs)
        mask%k(:,:,:, idim, ist, ik) = mask%k(:,:,:, idim, ist, ik) - cf1%Fs
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(M_old)
  SAFE_DEALLOCATE_A(psi)

  call zcube_function_free_RS(mask%cube, cf1)
  call  cube_function_free_FS(mask%cube, cf1)
  call zcube_function_free_RS(mask%cube, cf2)

  POP_SUB(pes_mask_restart_map)
end subroutine pes_mask_restart_map

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
