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



! Wrapper function
subroutine pes_flux_pmesh(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp, Ekin)
  type(pes_flux_t),  intent(in)    :: this
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)            
  FLOAT,             intent(in)    :: LG(:,:)           
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    
  integer,           intent(out)   :: idxZero(:)                
  integer,           intent(in)    :: krng(:)             
  integer,  pointer, intent(inout) :: Lp(:,:,:,:,:) 
  FLOAT,  optional,  intent(out)   :: Ekin(:,:,:)  
  

  PUSH_SUB(pes_flux_pmesh)


  select case (this%shape)
  
  case (M_SPHERICAL)
    call pes_flux_pmesh_sph(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp)
  
  case (M_CUBIC)
  ! not implemented

  case (M_PLANES)
    call pes_flux_pmesh_pln(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp, Ekin)
  
  end select
  


  POP_SUB(pes_flux_pmesh)  
end subroutine pes_flux_pmesh


! Wrapper function
subroutine pes_flux_map_from_states(this, restart, st, ll, pesP, krng, Lp, istin)
  type(pes_flux_t),    intent(in) :: this
  type(restart_t),     intent(in) :: restart
  type(states_elec_t), intent(in) :: st
  integer,             intent(in) :: ll(:)
  FLOAT, target,       intent(out) :: pesP(:,:,:,:)
  integer,             intent(in)  :: krng(:) 
  integer,  pointer,   intent(in) :: Lp(:,:,:,:,:)  
  integer, optional,   intent(in)  :: istin 

  PUSH_SUB(pes_flux_map_from_states)

  select case (this%shape)
  
  case (M_SPHERICAL)
    call pes_flux_map_from_states_elec_sph(this, restart, st, ll, pesP, krng, Lp, istin)
  
  case (M_CUBIC)

  case (M_PLANES)
    call pes_flux_map_from_states_elec_pln(this, restart, st, ll, pesP, krng, Lp, istin)
  
  end select

  POP_SUB(pes_flux_map_from_states)

end subroutine pes_flux_map_from_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PLANES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This needed in order to flip the sign of each Kpoit and 
! still preserve an array ordering on the kpoint mesh 
! such that the lowest index touple is associated with the smaller 
! (negative) kpoint value.
subroutine flip_sign_Lkpt_idx( dim, nk, idx, do_nothing)
   integer, intent(out) :: idx(:,:)
   integer, intent(in)  :: dim, nk(:)
   logical, intent(in)  :: do_nothing

   integer :: idx_tmp(1:maxval(nk(1:3)), 1:3)
   integer :: idir, ii

   PUSH_SUB(flip_sign_Lkpt_idx)
   
   idx(:,:) = 1
   do idir = 1, dim
     ! fill it with the range 1:nk
     do ii = 1, nk(idir)
       idx_tmp(ii,idir) = ii
     end do
     
     if (do_nothing) then
       idx(1:nk(idir),idir)  = idx_tmp(1:nk(idir),idir)
     else 
       do ii = 1, nk(idir)
         idx(ii,idir) = nk(idir) - idx_tmp(ii,idir) + 1 
       end do
     end if
     
   end do
   
   POP_SUB(flip_sign_Lkpt_idx)      
end subroutine flip_sign_Lkpt_idx

! ---------------------------------------------------------
integer pure function flatten_indices(i1,i2,i3,igpt, ll, ngpt) result(ii)
  integer, intent(in) :: i1
  integer, intent(in) :: i2
  integer, intent(in) :: i3
  integer, intent(in) :: igpt
  integer, intent(in) :: ll(:)
  integer, intent(in) :: ngpt
  
  ii = (i3-1)*ll(1)*ll(2)*ngpt + (i2-1)*ll(1)*ngpt + (i1-1)*ngpt + igpt 
  
end function flatten_indices



!< Generate the momentum-space mesh (p) and the arrays mapping the 
!< the mask and the kpoint meshes in p.
subroutine pes_flux_pmesh_pln(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp, Ekin)
  type(pes_flux_t),  intent(in)    :: this
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             !< ll(1:dim): the dimensions of the gpoint-mesh
  FLOAT,             intent(in)    :: LG(:,:)           !< LG(1:maxval(ll),1:dim): the  gpoints  
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    !< pmesh(i1,i2,i3,1:dim): contains the positions of point
                                                        !< in the final mesh in momentum space "p" combining the 
  integer,           intent(out) :: idxZero(:)          !< The triplet identifying the zero of the coordinates           

  integer,           intent(in)  :: krng(:)             !< The range identifying the zero-weight path 
                                                        !< mask-mesh with kpoints. 
  integer, pointer,  intent(out) :: Lp(:,:,:,:,:)       !< Allocated inside this subroutine
                                                        !< maps a mask-mesh triplet of indices together with a kpoint 
                                                        !< index into a triplet on the combined momentum space mesh.

  FLOAT,  optional,  intent(out) :: Ekin(:,:,:)         !< The total kinetic energy associated with the momentum p
                                                        !< this is needed when using a kpoint path 
   


  integer :: ik, j1, j2, j3, nk(1:3), ip1, ip2, ip3, idir, err
  FLOAT :: kpt(1:3),GG(1:3)

  integer, allocatable :: Lkpt(:,:), idx(:,:), idx_inv(:,:), ikidx(:,:)
  FLOAT, allocatable   :: LG_(:,:)

  integer :: nkpt, kpth_dir, ig, igpt
  FLOAT :: zero_thr
  
  logical :: use_zweight_path


  PUSH_SUB(pes_flux_pmesh_pln)
        
        
  SAFE_ALLOCATE(Lp(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3))
        
        
  nkpt = krng(2) - krng(1) + 1

  SAFE_ALLOCATE(Lkpt(krng(1):krng(2),1:3))
      
  nk(:) = 1  
  nk(1:dim) = kpoints%nik_axis(1:dim)

  Lkpt(:,:) = 1
  kpt(:) = M_ZERO
      
  zero_thr = M_EPSILON    

  use_zweight_path = (krng(1)>1) ! a dirty way to check whether we want to use the path or not 
      
  if ( kpoints_have_zero_weight_path(kpoints) .and. use_zweight_path) then     
    ! supporting paths only along the kx and ky directions in 
    ! reciprocal space
    
    kpth_dir = 1
    
    nk(:) = 1  
    nk(kpth_dir) = nkpt
    do ik = 1 , nkpt
      Lkpt(krng(1)+ik-1,kpth_dir) = ik
      kpt(1:dim) = kpoints_get_point(kpoints, krng(1) + ik -1) 
    end do
        
        
  else  
    
    call kpoints_grid_generate(dim, kpoints%nik_axis(1:dim), kpoints%full%nshifts, &
            kpoints%full%shifts(1:dim,1:kpoints%full%nshifts),  kpoints%full%red_point,  Lkpt(:,1:dim))


  end if
  SAFE_ALLOCATE(ikidx(maxval(nk(1:3)),1:3))
  call flip_sign_Lkpt_idx(dim, nk(:), ikidx(:,:), do_nothing = .true.)
  
  
  if (debug%info) then
    print *,"reordered"
    do ik = krng(1),krng(2)
      kpt(1:dim) = kpoints_get_point(kpoints, ik, absolute_coordinates = .true.)
      print *, ik, "Lkpt(ik)= [", ikidx(Lkpt(ik,1),1), ikidx(Lkpt(ik,2),2), ikidx(Lkpt(ik,3),3),&
                "] -- kpt= ",kpt(1:dim)
    end do

    print *,"----"
    print *,"ll(:)", ll(:)
    print *,"----"
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
  ! The grid represents the final momentum p =  Kpt -G. 
  ! The lower left corner correspond to the minimum value of p and the lowest 
  ! index-touple value (ip1,ip2,ip3) = (1,1,1). 
  do ik = krng(1),krng(2)
    kpt(1:dim) = kpoints_get_point(kpoints, ik) 
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 
          do igpt = 1, this%ngpt
          
            GG(:) = M_ZERO 
            ig = flatten_indices(j1,j2,j3, igpt, ll, this%ngpt) 

            GG(1:dim) = this%kcoords_cub(1:dim, ig, ik)
          
            ip1 = (j1 - 1) * nk(1) + ikidx(Lkpt(ik,1), 1)
            ip2 = (j2 - 1) * nk(2) + ikidx(Lkpt(ik,2), 2)
            ip3 = (j3 - 1) * nk(3) + ikidx(Lkpt(ik,3), 3)
          
            Lp(j1,j2,j3,ik,1:3) =  (/ip1,ip2,ip3/)
          
            ! The final momentum corresponds to p = K-G. 
            ! I have to subtract G only along the periodic dimensions
            pmesh(ip1, ip2, ip3, 1:this%pdim) = kpt(1:this%pdim) - GG(1:this%pdim)
            pmesh(ip1, ip2, ip3, dim)         = GG(dim)
            pmesh(ip1, ip2, ip3, dim+1)       = pmesh(ip1, ip2, ip3, dim+1) + 1 
          
          
            if (present(Ekin)) Ekin(ip1, ip2, ip3) = sign(M_ONE,pmesh(ip1,ip2,ip3,dim)) &
                                                     * sum(pmesh(ip1,ip2,ip3,1:dim)**2)/M_TWO

            if (debug%info) then
              print *,j1,j2,j3,ik, "pmesh = ",pmesh(ip1, ip2, ip3, :), "Ekin=", Ekin(ip1, ip2, ip3)
            end if

            ! Sanity checks
            if (sum(pmesh(ip1, ip2, ip3, 1:dim-1)**2)<=zero_thr) then
              err = err + 1 
              !Find the indices identifying the center of the coordinates 
              idxZero(1:3) = (/ip1,ip2,ip3/)
            end if
          
            if (pmesh(ip1, ip2, ip3, dim+1) > 1 ) then
              write(message(1),'(a)')'This condition should never happen something bad is going on.'
              call messages_fatal(1)
              err = -2 
            end if
      
          end do

        end do 
      end do 
    end do 
    
  end do
  
  if ( kpoints_have_zero_weight_path(kpoints)) then 
  ! With a path we just need to get the correct the zero index on the in-plane direction  
  ! perpendicular to the path since is along this direction that we are going 
  ! to slice with pes_flux_output_full_mapM_cut. Since on this direction we only 
  ! have G points I simply need to look for the zero index of the G-grid.
  ! Note that the G-grid must always include the (0,0,0) point. 
    do ik = krng(1),krng(2)
      do j1 = 1, ll(1) 
        do j2 = 1, ll(2) 
          do j3 = 1, ll(3) 
            do igpt = 1, this%ngpt
              
              ig = flatten_indices(j1,j2,j3, igpt, ll, this%ngpt) 
              GG(1:dim) = this%kcoords_cub(1:dim, ig, ik)
              if (sum(GG(1:dim-1)**2)<=M_EPSILON) idxZero(1:3) = (/j1,j2,j3/)
            end do
          end do
        end do
      end do
    end do
    
  else   
    
    if (err == -1) then
      call messages_write('Illformed momentum-space mesh: could not find p = 0 coordinate.')
      call messages_fatal()
    end if 

    if (err > 1) then
      call messages_write('More than one point with p = 0 coordinate.')
      call messages_new_line()
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

  
  POP_SUB(pes_flux_pmesh_pln)
end subroutine pes_flux_pmesh_pln




!< Build the photoemission map form the restart files
subroutine pes_flux_map_from_states_elec_pln(this, restart, st, ll, pesP, krng, Lp, istin)
  type(pes_flux_t),   intent(in) :: this
  type(restart_t),    intent(in) :: restart
  type(states_elec_t),intent(in) :: st
  integer,            intent(in) :: ll(:)
  FLOAT, target,     intent(out) :: pesP(:,:,:,:)
  integer,           intent(in)  :: krng(:) 
  integer,  dimension(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3), intent(in) :: Lp
  integer, optional, intent(in)  :: istin 
  
  integer :: ik, ist, idim, itot, nkpt, ispin
  integer :: i1, i2, i3, ip(1:3)
  integer :: idone, ntodo
  CMPLX   :: psiG1(1:this%nkpnts), psiG2(1:this%nkpnts)
  FLOAT   :: weight 
  integer :: istart, iend, nst, ig, igpt

  PUSH_SUB(pes_flux_map_from_states_elec_pln)

  istart = 1
  iend = st%nst
  nst = st%nst
  if(present(istin)) then
    istart = istin
    iend = istin
    nst = 1
  end if

  nkpt =  krng(2)-krng(1)+1
  ntodo = nkpt * nst 
  idone = 0 
  call loct_progress_bar(-1, ntodo)
  
  pesP = M_ZERO
  do ik = krng(1), krng(2)
    ispin = states_elec_dim_get_spin_index(st%d, ik)
    
    do ist = istart, iend

      if (st%d%kweights(ik) < M_EPSILON) then
        ! we have a zero-weight path
        ! the st%occ(ist, ik) factor is already in psiG[1,2]
        weight = M_ONE !/nkpt
      else
        weight = M_ONE * st%d%kweights(ik)
      end if
      
      if(st%d%ispin /= SPINORS) then 

        do idim = 1, st%d%dim
          itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
          call pes_flux_map_from_state_1(restart, itot, this%nkpnts, psiG1)
        
        
          do i1=1, ll(1)
            do i2=1, ll(2)
              do i3=1, ll(3)
                do igpt = 1, this%ngpt
                  ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
                  ig = flatten_indices(i1,i2,i3, igpt, ll, this%ngpt) 
              
                    pesP(ip(1),ip(2),ip(3), ispin) = pesP(ip(1),ip(2),ip(3), ispin) &
                                                   + abs(psiG1(ig))**2 * weight 
                
                end do
              end do
            end do
          end do
                
        end do
      else ! SPINORS
        idim = 1
        itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
        call pes_flux_map_from_state_1(restart, itot, this%nkpnts, psiG1)
        idim = 2
        itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
        call pes_flux_map_from_state_1(restart, itot, this%nkpnts, psiG2)
            
        do i1=1, ll(1)
          do i2=1, ll(2)
            do i3=1, ll(3)
              do igpt = 1, this%ngpt
                ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
                ig = flatten_indices(i1,i2,i3, igpt, ll, this%ngpt) 
            
                  pesP(ip(1),ip(2),ip(3), 1) = pesP(ip(1),ip(2),ip(3), 1) &
                                                 + abs(psiG1(ig))**2 * weight 

                  pesP(ip(1),ip(2),ip(3), 2) = pesP(ip(1),ip(2),ip(3), 2) &
                                                 + abs(psiG2(ig))**2 * weight

                  pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                                 + real(psiG1(ig)*conjg(psiG2(ig)), REAL_PRECISION) * weight
                                               
                  pesP(ip(1),ip(2),ip(3), 4) = pesP(ip(1),ip(2),ip(3), 4) &
                                                 + aimag(psiG1(ig)*conjg(psiG2(ig))) * weight
              end do
            end do
          end do
        end do
          
          
      end if
      
      idone = idone +1 
      call loct_progress_bar(idone, ntodo)
      
    end do
  end do

  write(stdout, '(1x)')

  POP_SUB(pes_flux_map_from_states_elec_pln)
end subroutine pes_flux_map_from_states_elec_pln



subroutine pes_flux_map_from_state_1(restart, idx, np, psiG)
  type(restart_t),  intent(in)  :: restart
  integer,          intent(in)  :: idx
  integer,          intent(in)  :: np
  CMPLX, target,    intent(out) :: psiG(:)

  character(len=80) :: filename, path
  integer ::  err, iunit 
  character(len=128) :: lines(2)

  PUSH_SUB(pes_flux_map_from_state_1)

  psiG = M_Z0
  
  write(filename,'(i10.10)') idx

  path = "pesflux1."//trim(filename)
  
  call zrestart_read_binary(restart, path, np, psiG(:), err)

  POP_SUB(pes_flux_map_from_state_1)  
end subroutine pes_flux_map_from_state_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPHERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine pes_flux_pmesh_sph(this, dim, kpoints, ll, LG, pmesh, idxZero, krng, Lp)
  type(pes_flux_t),  intent(in)    :: this
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             
  FLOAT,             intent(in)    :: LG(:,:)           
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    
  integer,           intent(out) :: idxZero(:)             
  integer,           intent(in)  :: krng(:)                                                                     
  integer, pointer,  intent(out) :: Lp(:,:,:,:,:)       
                                                       
                                                        

  integer            :: iomk
  integer            :: ikk, ith, iph, iphi
  FLOAT              :: phik, thetak, kact, kvec(1:3)
  
  integer            :: ip1, ip2, ip3
  

  PUSH_SUB(pes_flux_pmesh_sph)
                                   
  SAFE_ALLOCATE(Lp(1:this%nk, 1:this%nstepsomegak, 1, krng(1):krng(2), 1:3))                                 

  idxZero(1:3) =(/0,0,0/) 
  
  
  do ikk = 1, this%nk 
    kact = ikk * this%dk
    iomk = 0

    do ith = 0, this%nstepsthetak
      thetak = ith * M_PI / this%nstepsthetak 
      
      do iph = 0, this%nstepsphik - 1
        iomk = iomk + 1
        
        phik = iph * M_TWO * M_PI / this%nstepsphik

        if(ith == 0 .or. ith == this%nstepsthetak) then 
          ! Mark singular points on the sphere with -1 index
          Lp(ikk, iomk, 1, :, 1) = -1           
          exit
        end if


        ip1 = ikk
        ip2 = iph+1
        ip3 = ith 

        Lp(ikk, iomk, 1, :, 1) = ip1  
        Lp(ikk, iomk, 1, :, 2) = ip2
        Lp(ikk, iomk, 1, :, 3) = ip3
        
        kvec(1) = cos(phik) * sin(thetak)
        kvec(2) = sin(phik) * sin(thetak)
        kvec(3) = cos(thetak)
        
        kvec(:) = kvec(:) * kact
         
        
        pmesh(ip1, ip2, ip3, 1:3) = kvec(1:3)
        
        
      end do

    end do
  end do
  
                                                      
  POP_SUB(pes_flux_pmesh_sph)
  
end subroutine pes_flux_pmesh_sph



subroutine pes_flux_map_from_states_elec_sph(this, restart, st, ll, pesP, krng, Lp, istin)
  type(pes_flux_t),   intent(in) :: this
  type(restart_t),    intent(in) :: restart
  type(states_elec_t),intent(in) :: st
  integer,            intent(in) :: ll(:)
  FLOAT, target,     intent(out) :: pesP(:,:,:,:)
  integer,           intent(in)  :: krng(:) 
  integer,  dimension(1:this%nk,1:this%nstepsomegak,1,krng(1):krng(2),1:3), intent(in) :: Lp
  integer, optional, intent(in)  :: istin 

  integer :: ik, ist, idim, itot, nkpt, ispin
  integer :: i1, i2, i3, ip(1:3)
  integer :: idone, ntodo
  CMPLX   :: psiG1(1:this%nk, 1:this%nstepsomegak)
  CMPLX   :: psiG2(1:this%nk, 1:this%nstepsomegak)
  FLOAT   :: weight 
  integer :: istart, iend, nst, ig
  
  PUSH_SUB(pes_flux_map_from_states_elec_sph)

  istart = 1
  iend = st%nst
  nst = st%nst
  if(present(istin)) then
    istart = istin
    iend = istin
    nst = 1
  end if

  nkpt =  krng(2)-krng(1)+1
  ntodo = nkpt * nst 
  idone = 0 
  call loct_progress_bar(-1, ntodo)
  
  pesP = M_ZERO
  do ik = krng(1), krng(2)
    ispin = states_elec_dim_get_spin_index(st%d, ik)
    
    do ist = istart, iend

      if (st%d%kweights(ik) < M_EPSILON) then
        ! we have a zero-weight path
        weight = M_ONE!/nkpt
      else
        weight = M_ONE * st%d%kweights(ik)
      end if
      
      if(st%d%ispin /= SPINORS) then 

        do idim = 1, st%d%dim
          itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
          call pes_flux_map_from_state_2(restart, itot, this%nkpnts, psiG1)
        
          do i1=1, this%nk
            do i2=1, this%nstepsomegak
              
                ip(1:3) = Lp(i1, i2, 1, ik, 1:3) 
                if (ip(1) < 0) cycle
              
                pesP(ip(1),ip(2),ip(3), ispin) = pesP(ip(1),ip(2),ip(3), ispin) &
                                               + abs(psiG1(i1,i2))**2 * weight 
            end do
          end do
                
        end do
      else ! SPINORS
        idim = 1
        itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
        call pes_flux_map_from_state_2(restart, itot, this%nkpnts, psiG1)
        idim = 2
        itot = idim + (ist-1)*st%d%dim + (ik-1)*st%d%dim* st%nst
        call pes_flux_map_from_state_2(restart, itot, this%nkpnts, psiG2)
            
        do i1=1, this%nk
          do i2=1, this%nstepsomegak
            
              ip(1:3) = Lp(i1, i2, 1, ik, 1:3) 
              if (ip(1) < 0) cycle
            
              pesP(ip(1),ip(2),ip(3), 1) = pesP(ip(1),ip(2),ip(3), 1) &
                                             + abs(psiG1(i1,i2))**2 * weight 

              pesP(ip(1),ip(2),ip(3), 2) = pesP(ip(1),ip(2),ip(3), 2) &
                                             + abs(psiG2(i1,i2))**2 * weight

              pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                             + real(psiG1(i1,i2)*conjg(psiG2(i1,i2)), REAL_PRECISION) * weight
                                             
              pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                               + aimag(psiG1(i1,i2)*conjg(psiG2(i1,i2))) * weight
          end do
        end do
          
          
      end if
      
      idone = idone +1 
      call loct_progress_bar(idone, ntodo)
      
    end do
  end do

  write(stdout, '(1x)')


  POP_SUB(pes_flux_map_from_states_elec_sph)
  
end subroutine pes_flux_map_from_states_elec_sph



subroutine pes_flux_map_from_state_2(restart, idx, np, psiG)
  type(restart_t),  intent(in)  :: restart
  integer,          intent(in)  :: idx
  integer,          intent(in)  :: np
  CMPLX, target,    intent(out) :: psiG(:,:)

  character(len=80) :: filename, path
  integer ::  err, iunit 
  character(len=128) :: lines(2)

  PUSH_SUB(pes_flux_map_from_state_2)

  psiG = M_Z0
  
  write(filename,'(i10.10)') idx

  path = "pesflux1."//trim(filename)
  
  
  call zrestart_read_binary(restart, path, np, psiG(:,:), err)

  POP_SUB(pes_flux_map_from_state_2)  
end subroutine pes_flux_map_from_state_2



! ---------------------------------------------------------
subroutine pes_flux_out_energy(this, pesK, file, namespace, ll, pmesh, Ekin)
  type(pes_flux_t),  intent(in) :: this
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  FLOAT,             intent(in) :: Ekin(:,:,:)


  PUSH_SUB(pes_flux_out_energy)

  select case (this%shape)
  
  case (M_SPHERICAL)
    call messages_not_implemented("Energy-resolved PES for the flux method with spherical surfaces")
  ! not needed
  
  case (M_CUBIC)
    call messages_not_implemented("Energy-resolved PES for the flux method with cubic surfaces") 

  case (M_PLANES)
    call pes_flux_out_energy_pln(pesK, file, namespace, ll, pmesh, Ekin)
  
  end select


  POP_SUB(pes_flux_out_energy)
end subroutine pes_flux_out_energy


! ---------------------------------------------------------
subroutine pes_flux_out_energy_pln(arpes, file,namespace, ll, pmesh, Ekin)
  FLOAT,             intent(in) :: arpes(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  FLOAT,             intent(in) :: Ekin(:,:,:)
  
  integer :: iunit, ip, ie
  FLOAT   :: dp, length, pp(1:2), pdiff(1:2)
  
  
  PUSH_SUB(pes_flux_out_energy_pln)

  
  iunit = io_open(file, namespace, action='write')
  write(iunit, '(a)') '##################################################'
  write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                    str_center("E", 18),  str_center("P[E]", 18)

  write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                    str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18), &
                                    str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 18)
  write(iunit, '(a)') '##################################################'
  

  do ie = 1, ll(3) 
    write(iunit, '(es19.12,2x,es19.12,2x,es19.12)')   &
                                    units_from_atomic(units_out%energy, Ekin(1,1,ie)), &
                                    sum(arpes(:,:,ie))
  end do      
  
  
  call io_close(iunit)


  POP_SUB(pes_flux_out_energy_pln)
end subroutine pes_flux_out_energy_pln


subroutine pes_flux_out_energy_sph(pes, file, ll, pmesh, Ekin)
  FLOAT,             intent(in) :: pes(:,:,:)
  character(len=*),  intent(in) :: file
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  FLOAT,             intent(in) :: Ekin(:,:,:)
  
  PUSH_SUB(pes_flux_out_energy_sph)
  
  
  
  POP_SUB(pes_flux_out_energy_sph)
  
  
end subroutine pes_flux_out_energy_sph






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OUTPUT ON THE RUN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ---------------------------------------------------------
subroutine pes_flux_output(this, mesh, sb, st, namespace, dt)
  type(pes_flux_t), intent(inout)    :: this
  type(mesh_t),        intent(in)    :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(states_elec_t), intent(in)    :: st
  type(namespace_t),   intent(in)    :: namespace
  FLOAT,               intent(in)    :: dt

  integer            :: stst, stend, kptst, kptend, sdim, mdim
  integer            :: ist, ik, isdim
  integer            :: ikp, iomk, ikp_save, iomk_save
  integer            :: ikk, ith, iph, iphi
  FLOAT              :: phik, thetak, kact

  integer            :: iunitone, iunittwo
  FLOAT, allocatable :: spctrout_cub(:), spctrout_sph(:,:)
  FLOAT, allocatable :: spctrsum(:,:,:,:)
  FLOAT              :: weight
  
  ! M_PLANES debug
  integer            :: itot
  character(len=80)  :: filename

  PUSH_SUB(pes_flux_output)

  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = mesh%sb%dim

  SAFE_ALLOCATE(spctrsum(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nk))
  spctrsum = M_ZERO

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

  ! calculate spectra & total distribution
  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim

        ! orbital spectra
        select case (this%shape)
        case (M_SPHERICAL)

          do ikk = 1, this%nk 
            iomk = 0

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
                spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                  abs(this%spctramp_sph(ist, isdim, ik, ikk, iomk))**M_TWO * dt**M_TWO * weight

                if(ith == 0 .or. ith == this%nstepsthetak) exit
              end do
            end do
          end do
          ! distribution
          spctrout_sph(1:this%nk, 1:this%nstepsomegak) = spctrout_sph(1:this%nk, 1:this%nstepsomegak) + &
            abs(this%spctramp_sph(ist, isdim, ik, 1:this%nk, 1:this%nstepsomegak))**M_TWO * dt**M_TWO

        case (M_CUBIC)

          select case(mdim)
          case(1)
            weight = M_HALF

            ikk = 0
            do ikp = this%nk + 1, this%nkpnts
              ikk = ikk + 1
              spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                (abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO + &
                 abs(this%spctramp_cub(ist, isdim, ik, this%nkpnts + 1 - ikp))**M_TWO) * dt**M_TWO * weight
            end do
     
          case(2)
            weight = M_TWO * M_PI / this%nstepsphik

            ikp = 0
            do ikk = 1, this%nk
              do iph = 0, this%nstepsphik - 1
                ikp = ikp + 1
                spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + dt**M_TWO * &
                  abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO * dt**M_TWO * weight
              end do
            end do
     
          case(3)
            ikp  = 0
            do ikk = 1, this%nk
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
                  spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                    abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO * dt**M_TWO * weight
     
                  if(ith == 0 .or. ith == this%nstepsthetak) exit
                end do
              end do
            end do
          end select
          ! distribution
          spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
            abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO * dt**M_TWO

        end select
      end do
    end do
  end do

  if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
    ! total spectrum = sum over all states
    if(this%shape == M_SPHERICAL) then
      call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_sph)
    else
      call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout_cub)
    end if

    ! orbital spectra
    call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrsum)
#endif
  end if

  ! -----------------------------------------------------------------
  ! OUTPUT 
  ! -----------------------------------------------------------------
  if(mpi_grp_is_root(mpi_world)) then
    if (this%shape /= M_PLANES) then
      iunittwo = io_open('td.general/PES_flux.distribution.out', namespace, action='write', position='rewind')
      iunitone = io_open('td.general/'//'PES_flux.power.sum', namespace, action='write', position='rewind')
      write(iunitone, '(a19)') '# E, total spectrum'
    end if
    
    select case (this%shape)
    case (M_SPHERICAL)

      write(iunittwo, '(a29)') '# k, theta, phi, distribution'
      do ikk = 1, this%nk 
        kact = ikk * this%dk
        iomk = 0

        do ith = 0, this%nstepsthetak
          thetak = ith * M_PI / this%nstepsthetak 

          do iph = 0, this%nstepsphik - 1
            iomk = iomk + 1
            phik = iph * M_TWO * M_PI / this%nstepsphik
            if(iph == 0) iomk_save = iomk
            write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)

            ! just repeat the result for output
            if(this%nstepsphik > 1 .and. iph == (this%nstepsphik - 1)) &
              write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_sph(ikk, iomk_save)

            ! just repeat the result for output and exit
            if(ith == 0 .or. ith == this%nstepsthetak) then
              if(this%nstepsphik > 1) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, phik, spctrout_sph(ikk, iomk)
                end do
              end if
              exit
            end if
          end do

          if(this%nstepsphik > 1 .or. ith == this%nstepsthetak) write(iunittwo, '(1x)', advance='yes')
        end do
        write(iunitone, '(2(1x,e18.10E3))', advance='no') &
          kact**M_TWO / M_TWO, sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) * kact
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do isdim = 1, st%d%dim
              write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) * kact
            end do
          end do
        end do
        write(iunitone, '(1x)', advance='yes')
      end do

    case (M_CUBIC)

      select case(mdim)
      case(1)
        write(iunittwo, '(a17)') '# k, distribution'
        do ikp = 1, this%nkpnts
          write(iunittwo, '(2(1x,e18.10E3))') this%kcoords_cub(1, ikp, 1), spctrout_cub(ikp)
        end do

        do ikk = 1, this%nk
          kact = this%kcoords_cub(1, this%nk + ikk, 1)
          write(iunitone, '(2(1x,e18.10E3))', advance='no') &
            kact**M_TWO / M_TWO, sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) * kact
          
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) * kact
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')
        end do

      case(2)
        write(iunittwo, '(a22)') '# k, phi, distribution'
        ikp = 0
        do ikk = 1, this%nk
          kact = ikk * this%dk
          
          do iph = 0, this%nstepsphik - 1
            ikp = ikp + 1
            if(iph == 0) ikp_save = ikp
            phik = iph * M_TWO * M_PI / this%nstepsphik
            write(iunittwo,'(3(1x,e18.10E3))') kact, phik, spctrout_cub(ikp)
          end do
          ! just repeat the result for output
          write(iunittwo,'(3(1x,e18.10E3))') kact, M_TWO * M_PI, spctrout_cub(ikp_save)
          write(iunittwo, '(1x)', advance='yes')

          write(iunitone, '(2(1x,e18.10E3))', advance='no') &
            kact**M_TWO / M_TWO, sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) * kact
          
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) * kact
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')
        end do

      case(3)
        write(iunittwo, '(a29)') '# k, theta, phi, distribution'
        ikp    = 0
        do ikk = 1, this%nk
          kact = ikk * this%dk
          spctrsum = M_ZERO

          do ith = 0, this%nstepsthetak
            thetak = ith * M_PI / this%nstepsthetak 

            do iph = 0, this%nstepsphik - 1
              ikp = ikp + 1

              phik = iph * M_TWO * M_PI / this%nstepsphik
              if(iph == 0) ikp_save = ikp
              write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)

              ! just repeat the result for output
              if(iph == (this%nstepsphik - 1)) &
                write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, M_TWO * M_PI, spctrout_cub(ikp_save)

              ! just repeat the result for output and exit
              if(ith == 0 .or. ith == this%nstepsthetak) then
                do iphi = 1, this%nstepsphik
                  phik = iphi * M_TWO * M_PI / this%nstepsphik
                  write(iunittwo,'(4(1x,e18.10E3))') kact, thetak, phik, spctrout_cub(ikp)
                end do
                exit
              end if
            end do

            write(iunittwo, '(1x)', advance='yes')
          end do
          write(iunitone, '(2(1x,e18.10E3))', advance='no') &
            kact**M_TWO / M_TWO, sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) * kact
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) * kact
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')
        end do
      end select
    
    case (M_PLANES)
    
      do ik = kptst, kptend
        do ist = stst, stend
          do isdim = 1, sdim
            itot = ist + (ik-1) * st%nst +  (isdim-1) * st%nst * st%d%kpt%nglobal
            write(filename,'(i10.10)') itot
            
            iunitone = io_open('td.general/'//'PES_flux.distribution_'//trim(filename)//'.out', &
              namespace, action='write', position='rewind')
            write(iunitone, '(a29)') '# gx, gy, gz distribution'
            
            do ikp = 1, this%nkpnts
              
              select case(mdim)
              case (2)
                write(iunitone,'(3(1x,e18.10E3))') this%kcoords_cub(1, ikp, ik), &
                                                   this%kcoords_cub(2, ikp, ik), &
                                 abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO
              case (3)
                write(iunitone,'(4(1x,e18.10E3))') this%kcoords_cub(1, ikp, ik), &
                                                   this%kcoords_cub(2, ikp, ik), &
                                                   this%kcoords_cub(3, ikp, ik), &
                                 abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO
              end select
             
            end do
            
            call io_close(iunitone)
          end do      
        end do
      end do
      
    end select
    
    if (this%shape /= M_PLANES) then
      call io_close(iunittwo)
      call io_close(iunitone)
    end if
    
  end if

  SAFE_DEALLOCATE_A(spctrsum)
  SAFE_DEALLOCATE_A(spctrout_cub)
  SAFE_DEALLOCATE_A(spctrout_sph)

  POP_SUB(pes_flux_output)
end subroutine pes_flux_output





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESTART
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! ---------------------------------------------------------
subroutine pes_flux_dump(restart, this, mesh, st, ierr)
  type(restart_t),     intent(in)  :: restart
  type(pes_flux_t),    intent(in)  :: this
  type(mesh_t),        intent(in)  :: mesh
  type(states_elec_t), intent(in)  :: st
  integer,             intent(out) :: ierr

  integer          :: ist, ik, idim, itot
  integer          :: err
  integer          :: root(1:P_STRATEGY_MAX)
  character(len=128) :: filename

  CMPLX, pointer    :: psi1(:), psi2(:,:)
  
  
  PUSH_SUB(pes_flux_dump)


  if(debug%info) then
    message(1) = "Debug: Writing pes_flux restart."
    call messages_info(1)
  end if

  ierr = 0
  root = -1
  itot = 1
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim
        root(P_STRATEGY_DOMAINS) = mod(itot - 1, mesh%mpi_grp%size)
        write(filename,'(a,i10.10)') "pesflux1.", itot

        if (st%st_start <= ist .and. ist <= st%st_end .and. st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
          if(this%shape == M_SPHERICAL) then
            SAFE_ALLOCATE(psi2(1:this%nk, 1:this%nstepsomegak))
            psi2(:, :) = this%spctramp_sph(ist, idim, ik, :, :)
            call zrestart_write_binary(restart, filename, this%nk * this%nstepsomegak, psi2(:,:), err, root = root)
            SAFE_DEALLOCATE_P(psi2)
            
          else
            SAFE_ALLOCATE(psi1(1:this%nkpnts))
            psi1(:) = this%spctramp_cub(ist, idim, ik, :)
            call zrestart_write_binary(restart, filename, this%nkpnts, psi1(:), err, root = root)
            SAFE_DEALLOCATE_P(psi1)

          end if
        else 
          err = 0  
        end if

        if (err /= 0) ierr = ierr + 1
        itot = itot + 1
      end do
    end do
  end do

  if(this%shape == M_PLANES) then
    root(P_STRATEGY_MAX) = 0
    root(P_STRATEGY_KPOINTS) = -1
    do ik = 1, st%d%nik
      root(P_STRATEGY_DOMAINS) = mod(ik - 1, mesh%mpi_grp%size)
      write(filename,'(a,i5.5)') "pesflux4-kpt", ik

      if (st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
        SAFE_ALLOCATE(psi1(this%nkpnts))
        psi1(:)=this%conjgphase_prev_cub(:,ik)
        call zrestart_write_binary(restart, filename, this%nkpnts, psi1(:), err, root = root)
        SAFE_DEALLOCATE_P(psi1)
      else
        err = 0
      end if
      if (err /= 0) ierr = ierr + 1
    end do
  
  end if



  if(this%shape == M_SPHERICAL) then
    call zrestart_write_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
  else 
    if (this%shape /= M_PLANES) &
      call zrestart_write_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub(:,:), err)
  end if
  if(err /= 0) ierr = ierr + 2

  if(debug%info) then
    message(1) = "Debug: Writing pes_flux restart done."
    call messages_info(1)
  end if

  #ifdef HAVE_MPI
      call MPI_Barrier(mpi_world%comm, mpi_err)
      if(mpi_err /= 0) ierr = ierr + 3
  #endif


  POP_SUB(pes_flux_dump)
end subroutine pes_flux_dump

! ---------------------------------------------------------
subroutine pes_flux_load(restart, this, mesh, st, ierr)
  type(restart_t),     intent(in)    :: restart
  type(pes_flux_t),    intent(inout) :: this
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  integer,             intent(out)   :: ierr

  integer          :: stst, stend, kptst, kptend, sdim, mdim
  integer          :: ist, ik, idim, itot
  integer          :: err
  character(len=128) :: filename

  CMPLX, pointer    :: psi1(:), psi2(:,:)

  PUSH_SUB(pes_flux_load)

  if(restart_skip(restart)) then
    ierr = -1
    POP_SUB(pes_flux_load)
    return
  end if

  if(debug%info) then
    message(1) = "Debug: Reading pes_flux restart."
    call messages_info(1)
  end if

  
  ierr = 0
  itot = 1
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim
        write(filename,'(a,i10.10)') "pesflux1.", itot

        if (st%st_start <= ist .and. ist <= st%st_end .and. st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
          if(this%shape == M_SPHERICAL) then
            SAFE_ALLOCATE(psi2(1:this%nk, 1:this%nstepsomegak))
            call zrestart_read_binary(restart, filename, this%nk * this%nstepsomegak, psi2(:,:), err)
            this%spctramp_sph(ist, idim, ik, :, :) = psi2(:, :)
            SAFE_DEALLOCATE_P(psi2)
          
          else
            SAFE_ALLOCATE(psi1(1:this%nkpnts))
            call zrestart_read_binary(restart, filename, this%nkpnts, psi1(:), err)
            this%spctramp_cub(ist, idim, ik, :) =  psi1(:)
            SAFE_DEALLOCATE_P(psi1)

          end if
        else
          err = 0  
        end if
        
        if (err /= 0) ierr = ierr + 1        
        itot = itot + 1
      end do
    end do
  end do


  if(this%shape == M_PLANES) then
    do ik = 1, st%d%nik
      write(filename,'(a,i5.5)') "pesflux4-kpt", ik      
    
      if (st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
        SAFE_ALLOCATE(psi1(this%nkpnts))
        call zrestart_read_binary(restart, filename, this%nkpnts, psi1(:), err)
        this%conjgphase_prev_cub(:,ik)=psi1(:)
        SAFE_DEALLOCATE_P(psi1)
      else
        err = 0
      end if
      if (err /= 0) ierr = ierr + 1
    end do
  end if



  if(this%shape == M_SPHERICAL) then
    call zrestart_read_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev_sph, err)
  else
    if (this%shape /= M_PLANES) &
      call zrestart_read_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev_cub(:,:), err)
  end if
  if(err /= 0) ierr = ierr + 2
 
  if(debug%info) then
    message(1) = "Debug: Reading pes_flux restart done."
    call messages_info(1)
  end if
  
  #ifdef HAVE_MPI
      call MPI_Barrier(mpi_world%comm, mpi_err)
      if(mpi_err /= 0) ierr = ierr + 3
  #endif

  POP_SUB(pes_flux_load)
end subroutine pes_flux_load
