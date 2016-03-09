!! Copyright (C) 2015 P. Wopperer and U. De Giovannini
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
     
!          do ii=1, nk(idir)
!            print *,idir, ii, "idx = ", idx(ii,idir),"idx_tmp =", idx_tmp(ii,idir),mod(nk(idir),2)
!          end do
   end do
   
   
   POP_SUB(flip_sign_Lkpt_idx)      
end subroutine flip_sign_Lkpt_idx


!< Generate the momentum-space mesh (p) and the arrays mapping the 
!< the mask and the kpoint meshes in p.
subroutine pes_flux_pmesh(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp)
  type(pes_flux_t),  intent(in)    :: this
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

  integer :: nkpt, kpth_dir, ig
  FLOAT :: zero_thr


  PUSH_SUB(pes_flux_pmesh)
        
        
  nkpt = krng(2) - krng(1) + 1

  SAFE_ALLOCATE(Lkpt(krng(1):krng(2),1:3))
      
  nk(:) = 1  
  nk(1:dim) = kpoints%nik_axis(1:dim)

  Lkpt(:,:) = 1
  kpt(:) = M_ZERO
      
  zero_thr = M_EPSILON    
      
  if ( kpoints_have_zero_weight_path(kpoints)) then 
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
      kpt(1:dim) = kpoints_get_point(kpoints, krng(1) + ik -1) 
    end do
        
  else  
    
    call kpoints_grid_generate(dim, kpoints%nik_axis(1:dim), kpoints%shifts(1:dim), &
                               kpoints%full%red_point,  Lkpt(:,1:dim))


  end if

  SAFE_ALLOCATE(ikidx(maxval(nk(1:3)),1:3))
  call flip_sign_Lkpt_idx(dim, nk(:), ikidx(:,:))
  
  if (debug%info) then
    print *,"reordered"
    do ik = krng(1),krng(2)
      kpt(1:dim) = kpoints_get_point(kpoints, ik, absolute_coordinates = .false.)
      print *, ik, "Lkpt(ik)= [", ikidx(Lkpt(ik,1),1), ikidx(Lkpt(ik,2),2), ikidx(Lkpt(ik,3),3),"] -- kpt= ",kpt(1)
    end do

    print *,"----"
    print *,"ll(:)", ll(:)
    print *,"----"
  end if
  
!   ! We want the results to be sorted on a cube i,j,k
!   ! with the first triplet associated with the smallest positions
!   cSAFE_ALLOCATE(idx(1:maxval(ll(:)), 1:3))
!   cSAFE_ALLOCATE(idx_inv(1:maxval(ll(:)), 1:3))
!   cSAFE_ALLOCATE(LG_(1:maxval(ll(:)), 1:3))
!   idx(:,:)=1
!   idx_inv(:,:)=1
!   do idir = 1, dim
!     LG_(:,idir) = LG(:,idir)
!     call sort(LG_(1:ll(idir), idir), idx(1:ll(idir), idir))
!     idx_inv(:,idir) = idx(:,idir)
!     call sort(idx(1:ll(idir),idir),idx_inv(1:ll(idir),idir))
!   end do
!
!   if(debug%info) then
!     do idir=1, dim
!       print *, "*** direction =", idir
!       do j1 = 1, ll(idir)
!         print *,j1, "LG = ",LG(j1,idir),"LG_ = ",LG_(j1,idir), "idx = ", idx(j1,idir), "idx_inv = ", idx_inv(j1,idir)
!       end do
!       print *, "*** *** ***"
!     end do
!   end if

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
    kpt(1:dim) = kpoints_get_point(kpoints, ik) 
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 
          
          GG(:) = M_ZERO 
          ig = (j3-1)*ll(1)*ll(2) + (j2-1)*ll(1) + j1
!           ig = (j1-1)*ll(3)*ll(2) + (j2-1)*ll(3) + j3
          GG(1:dim) = this%kcoords_cub(1:dim, ig, ik)
!           print *, ik, j1, j2, j3, "GG(:) = ", GG(:) , ig
!           GG(1:3)= (/LG_(j1,1),LG_(j2,2),LG_(j3,3)/)
          
          ip1 = (j1 - 1) * nk(1) + ikidx(Lkpt(ik,1), 1)
          ip2 = (j2 - 1) * nk(2) + ikidx(Lkpt(ik,2), 2)
          ip3 = (j3 - 1) * nk(3) + ikidx(Lkpt(ik,3), 3)
          
!           Lp(idx_inv(j1,1),idx_inv(j2,2),idx_inv(j3,3),ik,1:3) =  (/ip1,ip2,ip3/)
          Lp(j1,j2,j3,ik,1:3) =  (/ip1,ip2,ip3/)
          
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
          if (sum(pmesh(ip1, ip2, ip3, 1:dim-1)**2)<=zero_thr) then
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

  if ( kpoints_have_zero_weight_path(kpoints)) then 
  ! With a path we just need to get the correct the zero index on the in-plane direction  
  ! perpendicular to the path since is along this direction that we are going 
  ! to slice with pes_flux_output_full_mapM_cut. Since on this direction we only 
  ! have G points I simply need to look for the zero index of the G-grid.
  ! Note that the G-grid must always include the (0,0,0) point. 
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 

          ig = (j3-1)*ll(1)*ll(2) + (j2-1)*ll(1) + j1
          GG(1:dim) = this%kcoords_cub(1:dim, ig, ik)
          if (sum(GG(1:dim-1)**2)<=M_EPSILON) idxZero(1:3) = (/j1,j2,j3/)
        
        end do
      end do
    end do
    
  else   
    
    if (err == -1) then
      call messages_write('Illformed momentum-space mesh: could not find p = 0 coordinate.')
      call messages_fatal()
    end if 

    if (err > 1) then
      call messages_write('Illformed momentum-space mesh: more than one point with p = 0 coordinate.')
      call messages_write('This can happen only if the kpoint mesh does not contain gamma.')
      call messages_warning()
    end if 

  end if

  if(debug%info) then
    print * ,"idxZero(1:3)=", idxZero(1:3)
  end if


  if (err == -2) then
    call messages_write('Illformed momentum-space mesh: two or more points with the same p.')
    call messages_fatal()
  end if 
  
 
 

  SAFE_DEALLOCATE_A(Lkpt)
  SAFE_DEALLOCATE_A(ikidx)      
!   cSAFE_DEALLOCATE_A(LG_)
!   cSAFE_DEALLOCATE_A(idx)
!   cSAFE_DEALLOCATE_A(idx_inv)
  
  
  POP_SUB(pes_flux_pmesh)
end subroutine pes_flux_pmesh


!< Build the photoemission map form the restart files
subroutine pes_flux_map_from_states(this, restart, st, ll, pesP, krng, Lp, istin)
  type(pes_flux_t),   intent(in) :: this
  type(restart_t),    intent(in) :: restart
  type(states_t),     intent(in) :: st
  integer,            intent(in) :: ll(:)
  FLOAT, target,     intent(out) :: pesP(:,:,:,:)
  integer,           intent(in)  :: krng(:) 
  integer,  dimension(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3), intent(in) :: Lp
  integer, optional, intent(in)  :: istin 
  
  integer :: ik, ist, idim, itot, nkpt, ispin
  integer :: i1, i2, i3, ip(1:3)
  integer :: idone, ntodo
  CMPLX   :: psiG1(this%nkpnts), psiG2(this%nkpnts)
  FLOAT   :: weight 
  integer :: istart, iend, nst, ig

  PUSH_SUB(pes_flux_map_from_states)

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
  
  pesP = M_ZERO
  do ik = krng(1), krng(2)
    ispin = states_dim_get_spin_index(st%d, ik)
    
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
          call pes_flux_map_from_state(restart, itot, this%nkpnts, psiG1)
        
          do i1=1, ll(1)
            do i2=1, ll(2)
              do i3=1, ll(3)
                ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
                ig = (i3-1)*ll(1)*ll(2) + (i2-1)*ll(1) + i1
              
                  pesP(ip(1),ip(2),ip(3), ispin) = pesP(ip(1),ip(2),ip(3), ispin) &
                                                 + abs(psiG1(ig))**2 * weight 
!                   if (all(ip(1:2)==(/155,17/))) then
!                     print *, itot, ig, ip(1:2), "psiG1(ig) =", psiG1(ig), & 
!                              "pesP(ip(1),ip(2),ip(3), ispin) =", pesP(ip(1),ip(2),ip(3), ispin) ,&
!                               "abs(psiG1(ig))**2 * weight =", abs(psiG1(ig))**2 * weight
!                   end if
              end do
            end do
          end do
                
        end do
      else ! SPINORS
        idim = 1
        itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
        call pes_flux_map_from_state(restart, itot, this%nkpnts, psiG1)
        idim = 2
        itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
        call pes_flux_map_from_state(restart, itot, this%nkpnts, psiG2)
            
        do i1=1, ll(1)
          do i2=1, ll(2)
            do i3=1, ll(3)
              ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
              ig = (i3-1)*ll(1)*ll(2) + (i2-1)*ll(1) + i1
            
                pesP(ip(1),ip(2),ip(3), 1) = pesP(ip(1),ip(2),ip(3), 1) &
                                               + abs(psiG1(ig))**2 * weight 

                pesP(ip(1),ip(2),ip(3), 2) = pesP(ip(1),ip(2),ip(3), 2) &
                                               + abs(psiG2(ig))**2 * weight

                pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                               + real(psiG1(ig)*conjg(psiG2(ig)), REAL_PRECISION) * weight
                                               
                pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                               + aimag(psiG1(ig)*conjg(psiG2(ig))) * weight
            end do
          end do
        end do
          
          
      end if
      
      idone = idone +1 
      call loct_progress_bar(idone, ntodo)
      
    end do
  end do

  write(stdout, '(1x)')

  POP_SUB(pes_flux_map_from_states)
end subroutine pes_flux_map_from_states


subroutine pes_flux_map_from_state(restart, idx, np, psiG)
  type(restart_t),  intent(in)  :: restart
  integer,          intent(in)  :: idx
  integer,          intent(in)  :: np
  CMPLX, target,    intent(out) :: psiG(:)

  character(len=80) :: filename, path
  integer ::  err, iunit 
  character(len=128) :: lines(2)

  PUSH_SUB(pes_flux_map_from_state)

  psiG = M_Z0
  
  write(filename,'(i10.10)') idx

  path = trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf"
  
  
  call io_binary_read(path, np, psiG(:), err)
  if (err /= 0) then
    message(1) = "Unable to read PES mask restart data from '"//trim(path)//"'."
    call messages_warning(1)
  end if

  POP_SUB(pes_flux_map_from_state)  
end subroutine pes_flux_map_from_state



! ---------------------------------------------------------
subroutine pes_flux_output(this, mesh, sb, st, dt)
  type(pes_flux_t), intent(inout)    :: this
  type(mesh_t),        intent(in)    :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(states_t),      intent(in)    :: st
  FLOAT,               intent(in)    :: dt

  integer            :: stst, stend, kptst, kptend, sdim, mdim
  integer            :: ist, ik, isdim
  integer            :: ikp, iomk, ikp_save, iomk_save
  integer            :: ikk, ith, iph, iphi
  FLOAT              :: phik, thetak, kact

  integer            :: iunitone, iunittwo
  FLOAT, allocatable :: spctrout_cub(:), spctrout_sph(:,:)
  FLOAT              :: weight, spctrsum

  PUSH_SUB(pes_flux_output)

  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = mesh%sb%dim

  select case (this%shape)

  case (M_SPHERICAL)
    SAFE_ALLOCATE(spctrout_sph(1:this%nk, 1:this%nstepsomegak))
    spctrout_sph = M_ZERO

  case (M_CUBIC)
    SAFE_ALLOCATE(spctrout_cub(1:this%nkpnts))
    spctrout_cub = M_ZERO

  case (M_PLANES)
    POP_SUB(pes_flux_output)
    return 
        
  end select

  ! calculate the total spectrum
  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim
        if(this%shape == M_SPHERICAL) then
          spctrout_sph(1:this%nk, 1:this%nstepsomegak) = spctrout_sph(1:this%nk, 1:this%nstepsomegak) + &
            abs(this%spctramp_sph(ist, isdim, ik, 1:this%nk, 1:this%nstepsomegak))**M_TWO * (dt * this%tdstepsinterval)**M_TWO &
            * st%occ(ist, ik)
        else
          spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
            abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO * (dt * this%tdstepsinterval)**M_TWO &
            * st%occ(ist, ik)
        end if
      end do
    end do
  end do

  if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
    if(this%shape == M_SPHERICAL) then
      call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_sph)
    else
      call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_cub)
    end if
#endif
  end if

  if(mpi_grp_is_root(mpi_world)) then
    iunittwo = io_open('td.general/PES_flux.distribution.out', action='write', position='rewind')
    iunitone = io_open('td.general/'//'PES_flux.power.sum', action='write', position='rewind')
    write(iunitone, '(a19)') '# E, total spectrum'

    if(this%shape == M_SPHERICAL) then
      write(iunittwo, '(a29)') '# k, theta, phi, distribution'
      do ikk = 1, this%nk 
        kact = ikk * this%dk
        iomk = 0
        spctrsum = M_ZERO

        do ith = 0, this%nstepsthetak
          thetak = ith * M_PI / this%nstepsthetak 

          if(ith == 0 .or. ith == this%nstepsthetak) then
            weight = (M_ONE - cos(M_PI / this%nstepsthetak / M_TWO)) * M_TWO * M_PI
          else
            weight = abs(cos(thetak - M_PI / this%nstepsthetak / M_TWO) - cos(thetak + M_PI / this%nstepsthetak / M_TWO)) &
              * M_TWO * M_PI / this%nstepsphik
          end if

          do iph = 0, this%nstepsphik - 1
            iomk = iomk + 1
            spctrsum = spctrsum + spctrout_sph(ikk, iomk) * weight 
            phik = iph * M_TWO * M_PI / this%nstepsphik
            if(iph == 0) iomk_save = iomk
            write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)

            ! just repeat the result for output
            if(this%nstepsphik > 1 .and. iph == (this%nstepsphik - 1)) &
              write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_sph(ikk, iomk_save)

            ! just repeat the result for output and exit
            if(ith == 0 .or. ith == this%nstepsthetak) then
              if(this%nstepsphik > 1) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)
                end do
              end if
              exit
            end if
          end do

          if(this%nstepsphik > 1 .or. ith == this%nstepsthetak) write(iunittwo, '(1x)', advance='yes')
        end do
        write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
      end do

    else ! this%shape == M_CUBIC
      select case(mdim)
      case(1)
        write(iunittwo, '(a17)') '# k, distribution'
        do ikp = 1, this%nkpnts
          write(iunittwo, '(5(1x,e18.10E3))') this%kcoords_cub(1, ikp, 1), spctrout_cub(ikp)
        end do

        do ikp = this%nk + 1, this%nkpnts
          kact = this%kcoords_cub(1, ikp, 1)
          write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, &
            (spctrout_cub(ikp) + spctrout_cub(this%nkpnts + 1 - ikp)) / M_TWO * kact
        end do

      case(2)
        write(iunittwo, '(a29)') '# k, phi, distribution'
        ikp = 0
        do ikk = 1, this%nk
          kact = ikk * this%dk
          
          spctrsum = M_ZERO
          do iph = 0, this%nstepsphik - 1
            ikp = ikp + 1
            if(iph == 0) ikp_save = ikp
            spctrsum = spctrsum + spctrout_cub(ikp) * M_TWO * M_PI / this%nstepsphik
            phik = iph * M_TWO * M_PI / this%nstepsphik
            write(iunittwo,'(5(1x,e18.10E3))') kact, phik, spctrout_cub(ikp)
          end do
          ! just repeat the result for output
          write(iunittwo,'(5(1x,e18.10E3))') kact, M_TWO * M_PI, spctrout_cub(ikp_save)
          write(iunittwo,'(1x)', advance = 'yes')
          write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
        end do

      case(3)
        write(iunittwo, '(a29)') '# k, theta, phi, distribution'
        ikp    = 0
        do ikk = 1, this%nk
          kact = ikk * this%dk
          spctrsum = M_ZERO

          do ith = 0, this%nstepsthetak
            thetak = ith * M_PI / this%nstepsthetak 

            if(ith == 0 .or. ith == this%nstepsthetak) then
              weight = (M_ONE - cos(M_PI / this%nstepsthetak / M_TWO)) * M_TWO * M_PI
            else
              weight = abs(cos(thetak - M_PI / this%nstepsthetak / M_TWO) - cos(thetak + M_PI / this%nstepsthetak / M_TWO)) &
                * M_TWO * M_PI / this%nstepsphik
            end if

            do iph = 0, this%nstepsphik - 1
              ikp = ikp + 1
              spctrsum = spctrsum + spctrout_cub(ikp) * weight

              phik = iph * M_TWO * M_PI / this%nstepsphik
              if(iph == 0) ikp_save = ikp
              write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)

              ! just repeat the result for output
              if(iph == (this%nstepsphik - 1)) &
                write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_cub(ikp_save)

              ! just repeat the result for output and exit
              if(ith == 0 .or. ith == this%nstepsthetak) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunittwo,'(5(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)
                end do
                exit
              end if
            end do

            write(iunittwo, '(1x)', advance='yes')
          end do
          write(iunitone, '(2(1x,e18.10E3))') kact**M_TWO / M_TWO, spctrsum * kact
        end do
      end select
    end if

    call io_close(iunittwo)
    call io_close(iunitone)
  end if

  SAFE_DEALLOCATE_A(spctrout_cub)
  SAFE_DEALLOCATE_A(spctrout_sph)

  POP_SUB(pes_flux_output)
end subroutine pes_flux_output

! ---------------------------------------------------------
subroutine pes_flux_dump(restart, this, mesh, st, ierr)
  type(restart_t),  intent(in)  :: restart
  type(pes_flux_t), intent(in)  :: this
  type(mesh_t),     intent(in)  :: mesh
  type(states_t),   intent(in)  :: st
  integer,          intent(out) :: ierr

  integer          :: stst, stend, kptst, kptend, sdim, mdim
  integer          :: ist, ik, isdim, itot
  integer          :: err
  character(len=128) :: filename

  PUSH_SUB(pes_flux_dump)

  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = mesh%sb%dim

  if(restart_skip(restart)) then
    POP_SUB(pes_flux_dump)
    return
  end if

  if(debug%info) then
    message(1) = "Debug: Writing pes_flux restart."
    call messages_info(1)
  end if

  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim
!           write(filename, '(i2.2, a, i2.2, a, i2.2)') ik, '.', ist, '.', isdim
        itot = ist + (ik-1) * st%nst+  (isdim-1) * st%nst*st%d%kpt%nglobal
        write(filename,'(i10.10)') itot

        if(mpi_grp_is_root(mesh%mpi_grp)) then

          if(this%shape == M_SPHERICAL) then
            call io_binary_write(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
              this%nk * this%nstepsomegak, this%spctramp_sph(ist, isdim, ik, :, :), err)
          else
            call io_binary_write(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
              this%nkpnts, this%spctramp_cub(ist, isdim, ik, :), err)
          end if

          call io_binary_write(trim(restart_dir(restart))//"/pesflux2."//trim(filename)//".obf", &
            this%nsrfcpnts * this%tdsteps, this%wf(ist, isdim, ik, :, :), err)
          call io_binary_write(trim(restart_dir(restart))//"/pesflux3."//trim(filename)//".obf", &
            this%nsrfcpnts * this%tdsteps * mdim, this%gwf(ist, isdim, ik, :, :, :), err)

        end if
        
      end do
    end do
  end do

  if(this%shape == M_SPHERICAL) then
    call zrestart_write_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
  else
    call zrestart_write_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub, err)
  end if

  call drestart_write_binary(restart, 'pesflux5', this%tdsteps * mdim, this%veca, err) 

  if(err /= 0) ierr = ierr + 1

  if(debug%info) then
    message(1) = "Debug: Writing pes_flux restart done."
    call messages_info(1)
  end if

  POP_SUB(pes_flux_dump)
end subroutine pes_flux_dump

! ---------------------------------------------------------
subroutine pes_flux_load(restart, this, mesh, st, ierr)
  type(restart_t),     intent(in)    :: restart
  type(pes_flux_t),    intent(inout) :: this
  type(mesh_t),        intent(in)    :: mesh
  type(states_t),      intent(in)    :: st
  integer,             intent(out)   :: ierr

  integer          :: stst, stend, kptst, kptend, sdim, mdim
  integer          :: ist, ik, isdim, itot
  integer          :: err
  character(len=128) :: filename

  PUSH_SUB(pes_flux_load)

  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = mesh%sb%dim

  if(restart_skip(restart)) then
    ierr = -1
    POP_SUB(pes_flux_load)
    return
  end if

  if(debug%info) then
    message(1) = "Debug: Reading pes_flux restart."
    call messages_info(1)
  end if

  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim
!           write(filename, '(i2.2, a, i2.2, a, i2.2)') ik, '.', ist, '.', isdim
        itot = ist + (ik-1) * st%nst+  (isdim-1) * st%nst*st%d%kpt%nglobal
        write(filename,'(i10.10)') itot

        if(this%shape == M_SPHERICAL) then
          call io_binary_read(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
            this%nk * this%nstepsomegak, this%spctramp_sph(ist, isdim, ik, :, :), err)
        else
          call io_binary_read(trim(restart_dir(restart))//"/pesflux1."//trim(filename)//".obf", &
            this%nkpnts, this%spctramp_cub(ist, isdim, ik, :), err)
        end if

        call io_binary_read(trim(restart_dir(restart))//"/pesflux2."//trim(filename)//".obf", &
          this%nsrfcpnts * this%tdsteps, this%wf(ist, isdim, ik, :, :), err)
        call io_binary_read(trim(restart_dir(restart))//"/pesflux3."//trim(filename)//".obf", &
          this%nsrfcpnts * this%tdsteps * mdim, this%gwf(ist, isdim, ik, :, :, :), err)
      end do
    end do
  end do

  if(this%shape == M_SPHERICAL) then
    call zrestart_read_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
  else
    call zrestart_read_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub, err)
  end if

  call drestart_read_binary(restart, 'pesflux5', this%tdsteps * mdim, this%veca, err) 

  if(err /= 0) ierr = ierr + 1
 
  if(debug%info) then
    message(1) = "Debug: Reading pes_flux restart done."
    call messages_info(1)
  end if

  POP_SUB(pes_flux_load)
end subroutine pes_flux_load
