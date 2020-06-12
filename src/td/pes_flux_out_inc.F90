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
subroutine pes_flux_pmesh(this, namespace, dim, kpoints, ll, pmesh, idxZero, krng, Lp, Ekin)
  type(pes_flux_t),  intent(in)    :: this
  type(namespace_t), intent(in)    :: namespace
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)            
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    
  integer,           intent(out)   :: idxZero(:)                
  integer,           intent(in)    :: krng(:)             
  integer,  pointer, intent(inout) :: Lp(:,:,:,:,:) 
  FLOAT,  optional,  intent(out)   :: Ekin(:,:,:)  
  

  PUSH_SUB(pes_flux_pmesh)


  select case (this%kgrid)
  
  case (M_POLAR)
    call pes_flux_pmesh_sph(this, dim, kpoints, ll, pmesh, idxZero, krng, Lp)
  
  case (M_CARTESIAN)
    if (kpoints_have_zero_weight_path(kpoints)) then
      call pes_flux_pmesh_pln(this, namespace, dim, kpoints, ll, pmesh, idxZero, krng, Lp, Ekin)
    else
      call pes_flux_pmesh_cub(this, namespace, dim, kpoints, ll, pmesh, idxZero, krng, Lp, Ekin)
    end if
  
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
  integer,  pointer,   intent(in)  :: Lp(:,:,:,:,:)  
  integer, optional,   intent(in)  :: istin 

  PUSH_SUB(pes_flux_map_from_states)

  select case (this%kgrid)
  
  case (M_POLAR)
    call pes_flux_map_from_states_elec_sph(this, restart, st, ll, pesP, krng, Lp, istin)
  
  case (M_CARTESIAN)
    call pes_flux_map_from_states_elec_pln(this, restart, st, ll, pesP, krng, Lp, istin)
  
  end select

  POP_SUB(pes_flux_map_from_states)

end subroutine pes_flux_map_from_states




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PLANES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This needed in order to flip the sign of each Kpoint and 
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
integer pure function flatten_indices(i1,i2,i3,ll) result(ii)
  integer, intent(in) :: i1
  integer, intent(in) :: i2
  integer, intent(in) :: i3
  integer, intent(in) :: ll(:)
  
  ii = (i3-1)*ll(1)*ll(2) + (i2-1)*ll(1) + (i1-1) + 1
  
end function flatten_indices



!< Generate the momentum-space mesh (p) and the arrays mapping the 
!< the mask and the kpoint meshes in p.
subroutine pes_flux_pmesh_pln(this, namespace, dim, kpoints, ll, pmesh, idxZero, krng, Lp, Ekin)
  type(pes_flux_t),  intent(in)    :: this
  type(namespace_t), intent(in)    :: namespace
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             !< ll(1:dim): the dimensions of the gpoint-mesh
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
  FLOAT :: kpt(1:3)

  integer, allocatable :: Lkpt(:,:), idx(:,:), idx_inv(:,:), ikidx(:,:)

  integer :: nkpt, kpth_dir, ikp
  FLOAT :: zero_thr
  

  PUSH_SUB(pes_flux_pmesh_pln)
        
        
  SAFE_ALLOCATE(Lp(1:ll(1),1:ll(2),1:ll(3),krng(1):krng(2),1:3))
        
        
  nkpt = krng(2) - krng(1) + 1

  SAFE_ALLOCATE(Lkpt(krng(1):krng(2),1:3))
      
  nk(:) = 1  
  nk(1:dim) = kpoints%nik_axis(1:dim)

  Lkpt(:,:) = 1
  kpt(:) = M_ZERO
      
  zero_thr = M_EPSILON    

  
  kpth_dir = 1
  
  nk(:) = 1  
  nk(kpth_dir) = nkpt
  do ik = 1 , nkpt
    Lkpt(krng(1)+ik-1,kpth_dir) = ik
    kpt(1:dim) = kpoints_get_point(kpoints, krng(1) + ik -1) 
  end do
        
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
  ! The lower left corner correspond to the minimum value of p and the lowest 
  ! index-touple value (ip1,ip2,ip3) = (1,1,1). 
  do ik = krng(1),krng(2)
    do j1 = 1, ll(1) 
      do j2 = 1, ll(2) 
        do j3 = 1, ll(3) 
          
          kpt(:) = M_ZERO 
          ikp = flatten_indices(j1,j2,j3,ll)

          kpt(1:dim) = this%kcoords_cub(1:dim, ikp, ik)
        
          ip1 = (j1 - 1) * nk(1) + ikidx(Lkpt(ik,1), 1)
          ip2 = (j2 - 1) * nk(2) + ikidx(Lkpt(ik,2), 2)
          ip3 = (j3 - 1) * nk(3) + ikidx(Lkpt(ik,3), 3)
        
          Lp(j1,j2,j3,ik,1:3) =  (/ip1,ip2,ip3/)
        
          ! The final momentum corresponds to p = K. 
          pmesh(ip1, ip2, ip3, 1:dim)         = kpt(1:dim)

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
            call messages_fatal(1, namespace=namespace)
            err = -2 
          end if
      

        end do 
      end do 
    end do 
    
  end do
  

  if(debug%info) then
    print * ,"idxZero(1:3)=", idxZero(1:3)
  end if

  if (err == -2) then
    call messages_write('Malformed momentum-space mesh: two or more points with the same p.')
    call messages_fatal(namespace=namespace)
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
  integer :: istart, iend, nst, iflt

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
        ! we have a zero-weifltht path
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
                ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
                iflt = flatten_indices(i1,i2,i3, ll) 
            
                pesP(ip(1),ip(2),ip(3), ispin) = pesP(ip(1),ip(2),ip(3), ispin) &
                                               + abs(psiG1(iflt))**2 * weight 
                
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
              ip(1:3) = Lp(i1, i2, i3, ik, 1:3) 
              iflt = flatten_indices(i1,i2,i3, ll) 
          
              pesP(ip(1),ip(2),ip(3), 1) = pesP(ip(1),ip(2),ip(3), 1) &
                                             + abs(psiG1(iflt))**2 * weight 

              pesP(ip(1),ip(2),ip(3), 2) = pesP(ip(1),ip(2),ip(3), 2) &
                                             + abs(psiG2(iflt))**2 * weight

              pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                             + real(psiG1(iflt)*conjg(psiG2(iflt)), REAL_PRECISION) * weight
                                           
              pesP(ip(1),ip(2),ip(3), 4) = pesP(ip(1),ip(2),ip(3), 4) &
                                             + aimag(psiG1(iflt)*conjg(psiG2(iflt))) * weight
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
  integer ::  err 

  PUSH_SUB(pes_flux_map_from_state_1)

  psiG = M_Z0
  
  write(filename,'(i10.10)') idx

  path = "pesflux1."//trim(filename)
  
  call zrestart_read_binary(restart, path, np, psiG(:), err)

  POP_SUB(pes_flux_map_from_state_1)  
end subroutine pes_flux_map_from_state_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUBE (parallelepiped in spirit)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pes_flux_pmesh_cub(this, namespace, dim, kpoints, ll, pmesh, idxZero, krng, Lp, Ekin)
  type(pes_flux_t),  intent(in)    :: this
  type(namespace_t), intent(in)    :: namespace
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             !< ll(1:dim): the dimensions of the gpoint-mesh
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    !< pmesh(i1,i2,i3,1:dim): contains the positions of point
                                                        !< in the final mesh in momentum space "p" combining the 
  integer,           intent(out) :: idxZero(:)          !< The triplet identifying the zero of the coordinates           

  integer,           intent(in)  :: krng(:)             !< The range identifying the zero-weight path 
                                                        !< mask-mesh with kpoints. 
  integer, pointer,  intent(out) :: Lp(:,:,:,:,:)       !< Allocated inside this subroutine
                                                        !< maps a mask-mesh triplet of indices together with a kpoint 
                                                        !< index into a triplet on the combined momentum space mesh.

  FLOAT,  optional,  intent(out) :: Ekin(:,:,:)         !< The total kinetic energy associated with the momentum p
  
  
  
  integer :: ikpt,ikp, ik1, ik2, ik3
  FLOAT   :: vec(1:3), tmp, min
  
  PUSH_SUB(pes_flux_pmesh_cub)
  
  idxZero(1:3) =(/(this%ll(1)+1)/2,(this%ll(2)+1)/2,(this%ll(3)+1)/2 /) 
  
  SAFE_ALLOCATE(Lp(1:this%ll(1), 1:this%ll(2), this%ll(3), krng(1):krng(2), 1:3))          
  
  min = M_HUGE
  
  ikpt = 1
  ikp = 0
  do ik3 = 1, this%ll(3)
    do ik2 = 1, this%ll(2)
      do ik1 = 1, this%ll(1)
        ikp = ikp + 1 
        

        Lp(ik1, ik2, ik3, :, 1) = ik1  
        Lp(ik1, ik2, ik3, :, 2) = ik2  
        Lp(ik1, ik2, ik3, :, 3) = ik3  
        
        pmesh(ik1, ik2, ik3, 1:dim) = this%kcoords_cub(1:dim, ikp, ikpt)

        Ekin(ik1, ik2, ik3) = sum(this%kcoords_cub(1:dim, ikp, ikpt)**2)*M_HALF
        
!         ! get the origin index
!         tmp=sum(pmesh(ik1, ik2, ik3, 1:dim)**2)
!         if (tmp<min) then
!           min = tmp
!           idxZero(1:3) = (/ik1, ik2, ik3/)
!         end if

      end do
    end do
  end do
  
  
  
  POP_SUB(pes_flux_pmesh_cub)
end subroutine pes_flux_pmesh_cub



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SPHERE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine pes_flux_pmesh_sph(this, dim, kpoints, ll, pmesh, idxZero, krng, Lp)
  type(pes_flux_t),  intent(in)    :: this
  integer,           intent(in)    :: dim
  type(kpoints_t),   intent(inout) :: kpoints 
  integer,           intent(in)    :: ll(:)             
  FLOAT,             intent(out)   :: pmesh(:,:,:,:)    
  integer,           intent(out)   :: idxZero(:)             
  integer,           intent(in)    :: krng(:)                                                                     
  integer, pointer,  intent(out)   :: Lp(:,:,:,:,:)       
                                                       
                                                        

  integer            :: iomk
  integer            :: ikk, ith, iph, iphi
  FLOAT              :: phik, thetak, kact, kvec(1:3), Dthetak, Dphik
  
  integer            :: ip1, ip2, ip3

  PUSH_SUB(pes_flux_pmesh_sph)

  SAFE_ALLOCATE(Lp(1:this%nk, 1:this%nstepsomegak, 1, krng(1):krng(2), 1:3))

  idxZero(1:3) = (/0,0,0/)
  
  Dthetak  = M_ZERO
  if (dim ==3) Dthetak = abs(this%thetak_rng(2) - this%thetak_rng(1))/(this%nstepsthetak)
  Dphik = abs(this%phik_rng(2) - this%phik_rng(1))/(this%nstepsphik)
  
  do ikk = 1, this%nk 
    kact = this%klinear(ikk,1)
    iomk = 0

    do ith = 0, this%nstepsthetak
!       thetak = ith * M_PI / this%nstepsthetak
      thetak = ith * Dthetak + this%thetak_rng(1)
      
      do iph = 0, this%nstepsphik - 1
        iomk = iomk + 1

!         phik = iph * M_TWO * M_PI / this%nstepsphik
        phik = iph * Dphik + this%phik_rng(1)        

!         if(ith == 0 .or. ith == this%nstepsthetak) then 
        if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) then  
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
  FLOAT, target,       intent(out) :: pesP(:,:,:,:)
  integer,             intent(in)  :: krng(:) 
  integer,             intent(in)  :: Lp(1:this%nk,1:this%nstepsomegak,1,krng(1):krng(2),1:3)
  integer, optional,   intent(in)  :: istin 

  integer :: ik, ist, idim, itot, nkpt, ispin
  integer :: i1, i2, ip(1:3)
  integer :: idone, ntodo
  CMPLX   :: psiG1(1:this%nk, 1:this%nstepsomegak)
  CMPLX   :: psiG2(1:this%nk, 1:this%nstepsomegak)
  FLOAT   :: weight 
  integer :: istart, iend, nst
  
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
        
          do i1 = 1, this%nk
            do i2 = 1, this%nstepsomegak
              
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
            
        do i1 = 1, this%nk
          do i2 = 1, this%nstepsomegak
            
              ip(1:3) = Lp(i1, i2, 1, ik, 1:3) 
              if (ip(1) < 0) cycle
            
              pesP(ip(1),ip(2),ip(3), 1) = pesP(ip(1),ip(2),ip(3), 1) &
                                             + abs(psiG1(i1,i2))**2 * weight 

              pesP(ip(1),ip(2),ip(3), 2) = pesP(ip(1),ip(2),ip(3), 2) &
                                             + abs(psiG2(i1,i2))**2 * weight

              pesP(ip(1),ip(2),ip(3), 3) = pesP(ip(1),ip(2),ip(3), 3) &
                                             + TOFLOAT(psiG1(i1,i2)*conjg(psiG2(i1,i2))) * weight
                                             
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
  integer ::  err

  PUSH_SUB(pes_flux_map_from_state_2)

  psiG = M_Z0
  
  write(filename,'(i10.10)') idx

  path = "pesflux1."//trim(filename)
  
  
  call zrestart_read_binary(restart, path, np, psiG(:,:), err)

  POP_SUB(pes_flux_map_from_state_2)  
end subroutine pes_flux_map_from_state_2



! ---------------------------------------------------------
subroutine pes_flux_out_energy(this, pesK, file, namespace, ll, pmesh, Ekin, dim)
  type(pes_flux_t),  intent(in) :: this
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  FLOAT,             intent(in) :: Ekin(:,:,:)
  integer,           intent(in) :: dim


  PUSH_SUB(pes_flux_out_energy)

  select case (this%kgrid)
  
  case (M_POLAR)
    call messages_not_implemented("Energy-resolved PES for the flux method polar momentum grids", namespace=namespace)
!     call  pes_flux_out_polar_ascii(this, st, namespace, dim, efile = file)
  
  case (M_CARTESIAN)
    if (this%surf_shape == M_CUBIC) then
      call pes_flux_out_energy_pln(pesK, file, namespace, ll, pmesh, Ekin, dim)
    else 
      call messages_not_implemented("Energy-resolved PES for the flux method with cartesian momentum grids", namespace=namespace)
    end if
    
  end select


  POP_SUB(pes_flux_out_energy)
end subroutine pes_flux_out_energy


! ---------------------------------------------------------
subroutine pes_flux_out_energy_pln(arpes, file,namespace, ll, pmesh, Ekin, dim)
  FLOAT,             intent(in) :: arpes(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  FLOAT,             intent(in) :: Ekin(:,:,:)
  integer,           intent(in) :: dim
  
  integer :: iunit, ie
  
  PUSH_SUB(pes_flux_out_energy_pln)

  
  iunit = io_open(file, namespace, action='write')
  write(iunit, '(a)') '##################################################'
  write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                    str_center("E", 18),  str_center("P[E]", 18)

  write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                    str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18), &
                                    str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 18)
  write(iunit, '(a)') '##################################################'
  

  select case(dim)
  case (2)
    do ie = 1, ll(2) 
      write(iunit, '(es19.12,2x,es19.12,2x,es19.12)')   &
                                      units_from_atomic(units_out%energy, Ekin(1,ie,1)), &
                                      sum(arpes(:,ie,:))
    end do      
  case (3)
    do ie = 1, ll(3) 
      write(iunit, '(es19.12,2x,es19.12,2x,es19.12)')   &
                                      units_from_atomic(units_out%energy, Ekin(1,1,ie)), &
                                      sum(arpes(:,:,ie))
    end do      
  end select
  
  
  
  call io_close(iunit)


  POP_SUB(pes_flux_out_energy_pln)
end subroutine pes_flux_out_energy_pln


! ---------------------------------------------------------
subroutine pes_flux_out_vmap(this, pesK, file, namespace, ll, pmesh, dim)
  type(pes_flux_t), intent(inout)    :: this
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  integer,           intent(in) :: dim
  
  PUSH_SUB(pes_flux_out_vmap)
  
  select case (this%kgrid)

    case (M_POLAR)

    case (M_CARTESIAN)
      if (this%surf_shape == M_CUBIC) &
        call pes_flux_out_vmap_cub(pesK, file, namespace, ll, pmesh, dim)
      
  end select
  
  
  POP_SUB(pes_flux_out_vmap)
end subroutine pes_flux_out_vmap



! ---------------------------------------------------------
subroutine pes_flux_out_vmap_cub(pesK, file, namespace, ll, pmesh, dim)
  FLOAT,             intent(in) :: pesK(:,:,:)
  character(len=*),  intent(in) :: file
  type(namespace_t), intent(in) :: namespace
  integer,           intent(in) :: ll(:)  
  FLOAT,             intent(in) :: pmesh(:,:,:,:)  
  integer,           intent(in) :: dim
  
  integer :: ik1,ik2,ik3, idir, idim
  integer :: iunit
  
  PUSH_SUB(pes_flux_out_vmap_cub)

  iunit = io_open(file, namespace, action='write', position='rewind')
  write(iunit, '(a)') '##################################################'                                        
  if (dim == 3) then 
    write(iunit, '(a1,a18,2x,a18,2x,a18,2x,a18)') '#', &
                                      str_center("kz", 18), str_center("ky", 18),&
                                      str_center("kx", 18),  str_center("P(kz,ky,kx)", 18)
    write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)    
  end if
  
  if (dim == 2) then 
    write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
                                    str_center("ky", 18), str_center("kx", 18),&
                                    str_center("P(ky,kx)", 18)
    write(iunit, '(a1,a18,2x,a18)') '#', &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)    

  end if  
                                    

  if (dim == 1) then
    write(iunit, '(a1,a18,2x,a18)') '#', &
                                    str_center("kx", 18),  str_center("P(kx)", 18)
    write(iunit, '(a1,a18)') '#', &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)    
                                    
  end if
  write(iunit, '(a)') '##################################################'
  
  

  do ik3 = 1, ll(3)
    do ik2 = 1, ll(2)
      do ik1 = 1, ll(1)
        
        do idir = dim, 1, -1
          write(iunit, '(1x,e18.10E3)', advance='no') & 
            units_from_atomic(unit_one/units_out%length,pmesh(ik1,ik2,ik3, idir))
        end do
        write(iunit, '(1x,e18.10E3)', advance='no') pesK(ik1,ik2,ik3)
        write(iunit, '(1x)', advance='yes')
        
      end do
      write(iunit, '(1x)', advance='yes')
    end do
    write(iunit, '(1x)', advance='yes')
  end do
  
  call io_close(iunit)
  
  

  POP_SUB(pes_flux_out_vmap_cub)
    
end subroutine pes_flux_out_vmap_cub




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OUTPUT ON THE RUN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pes_flux_output(this, mesh, sb, st, namespace)
  type(pes_flux_t), intent(inout)    :: this
  type(mesh_t),        intent(in)    :: mesh
  type(simul_box_t),   intent(in)    :: sb
  type(states_elec_t), intent(in)    :: st
  type(namespace_t),   intent(in)    :: namespace

  PUSH_SUB(pes_flux_output)
  
  if ( .not. this%runtime_output) then
    POP_SUB(pes_flux_output)
    return
  end if
  
  select case (this%kgrid)
  
  case (M_POLAR)
    call  pes_flux_out_polar_ascii(this, st, namespace, sb%dim,&
                                  mfile = io_workpath("td.general/PES_flux.distribution.out", namespace), &
                                  efile = io_workpath("td.general/PES_flux.power.sum", namespace))
  
  case (M_CARTESIAN)
    call pes_flux_out_cartesian_ascii(this, st, namespace, sb%dim, io_workpath("td.general/", namespace))
    
  case default
    !empty      
    
  end select
  
  POP_SUB(pes_flux_output)
end subroutine pes_flux_output


! ---------------------------------------------------------
subroutine pes_flux_out_cartesian_ascii(this, st, namespace, dim, path )
  type(pes_flux_t), intent(inout)    :: this
  type(states_elec_t), intent(in)    :: st
  type(namespace_t),   intent(in)    :: namespace
  integer,             intent(in)    :: dim
  character(len=*),    intent(in)    :: path 
  
  integer            :: stst, stend, kptst, kptend, sdim, mdim, idir
  integer            :: iunit
  integer            :: ik, ist, isdim, ikp, ikpt, ik1, ik2, ik3
  FLOAT, allocatable ::  spctrout(:,:,:), pmesh(:,:,:,:)
  
  PUSH_SUB(pes_flux_out_cartesian_ascii)
  
  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = dim

  SAFE_ALLOCATE(spctrout(1:this%ll(1), 1:this%ll(2), 1:this%ll(3)))
  SAFE_ALLOCATE(pmesh(1:this%ll(1), 1:this%ll(2), 1:this%ll(3),1:dim))
  spctrout = M_ZERO

  ! calculate spectra & total distribution
  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim

        ikp = 0

        do ik3 = 1, this%ll(3)
          do ik2 = 1, this%ll(2)
            do ik1 = 1, this%ll(1)
              ikp = ikp + 1 

              spctrout(ik1,ik2,ik3) = spctrout(ik1,ik2,ik3) + &
                abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO
              pmesh(ik1,ik2,ik3, 1:dim) =  this%kcoords_cub(1:dim, ikp, 1)
                
            end do
          end do
        end do
                  
        
      end do 
    end do
  end do
  
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
      call comm_allreduce(st%st_kpt_mpi_grp%comm, spctrout)
#endif
  end if
  
  
  if(mpi_grp_is_root(mpi_world)) then
    call pes_flux_out_vmap_cub(spctrout, io_workpath("td.general/PES_flux.distribution.out", namespace), &
                               namespace, this%ll(:), pmesh, mdim)

!     iunit = io_open(trim(path)//'PES_flux.distribution.out', namespace, action='write', position='rewind')
!     write(iunit, '(a)') '##################################################'
!     if (mdim == 3) then
!       write(iunit, '(a1,a18,2x,a18,2x,a18,2x,a18)') '#', &
!                                         str_center("kz", 18), str_center("ky", 18),&
!                                         str_center("kx", 18),  str_center("P(kz,ky,kx)", 18)
!       write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)
!     end if
!
!     if (mdim == 2) then
!       write(iunit, '(a1,a18,2x,a18,2x,a18)') '#', &
!                                       str_center("ky", 18), str_center("kx", 18),&
!                                       str_center("P(ky,kx)", 18)
!       write(iunit, '(a1,a18,2x,a18)') '#', &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)
!
!     end if
!
!
!     if (mdim == 1) then
!       write(iunit, '(a1,a18,2x,a18)') '#', &
!                                       str_center("kx", 18),  str_center("P(kx)", 18)
!       write(iunit, '(a1,a18)') '#', &
!                                         str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)
!
!     end if
!     write(iunit, '(a)') '##################################################'
!
!
!     ikpt = 1
!     ikp = 0
!
!     do ik3 = 1, this%ll(3)
!       do ik2 = 1, this%ll(2)
!         do ik1 = 1, this%ll(1)
!           ikp = ikp + 1
!
!           do idir = mdim, 1, -1
!             write(iunit, '(1x,e18.10E3)', advance='no') &
!               units_from_atomic(unit_one/units_out%length,this%kcoords_cub(idir, ikp, ikpt))
!           end do
!           write(iunit, '(1x,e18.10E3)', advance='no') spctrout(ik1,ik2,ik3)
!           write(iunit, '(1x)', advance='yes')
!
!         end do
!         write(iunit, '(1x)', advance='yes')
!       end do
!       write(iunit, '(1x)', advance='yes')
!     end do
!
!     call io_close(iunit)
  end if
  
   SAFE_DEALLOCATE_A(spctrout)
  SAFE_DEALLOCATE_A(pmesh)
  
  
  POP_SUB(pes_flux_out_cartesian_ascii)
end subroutine pes_flux_out_cartesian_ascii


! ---------------------------------------------------------
subroutine pes_flux_out_polar_ascii(this, st, namespace, dim, efile, mfile)
  type(pes_flux_t), intent(inout)    :: this
  type(states_elec_t), intent(in)    :: st
  type(namespace_t),   intent(in)    :: namespace
  integer,             intent(in)    :: dim
  character(len=*), optional,   intent(in)    :: efile
  character(len=*), optional,   intent(in)    :: mfile
    
  integer            :: stst, stend, kptst, kptend, sdim, mdim
  integer            :: ist, ik, isdim
  integer            :: ikp, iomk, ikp_save, iomk_save
  integer            :: ikk, ith, iph, iphi
  FLOAT              :: phik, thetak, kact,kmin, Dthetak, Dphik, Lphik

  integer            :: iunitone, iunittwo
  FLOAT, allocatable :: spctrout_cub(:), spctrout_sph(:,:)
  FLOAT, allocatable :: spctrsum(:,:,:,:)
  FLOAT              :: weight
  logical           :: energy_resolved, momentum_resolved

  
  integer            :: itot

  PUSH_SUB(pes_flux_out_polar_ascii)

  stst   = st%st_start
  stend  = st%st_end
  kptst  = st%d%kpt%start
  kptend = st%d%kpt%end
  sdim   = st%d%dim
  mdim   = dim

  energy_resolved   = .false.
  momentum_resolved = .false.
  if (present(efile)) energy_resolved   = .true.
  if (present(mfile)) momentum_resolved = .true.
   
  if (.not. energy_resolved .and. .not. momentum_resolved) then
    POP_SUB(pes_flux_out_polar_ascii)
    return    
  end if

  SAFE_ALLOCATE(spctrsum(1:st%nst, 1:sdim, 1:st%d%nik, 1:this%nk))
  spctrsum = M_ZERO

  if (this%surf_shape == M_SPHERICAL) then
    SAFE_ALLOCATE(spctrout_sph(1:this%nk, 1:this%nstepsomegak))
    spctrout_sph = M_ZERO
  else
    SAFE_ALLOCATE(spctrout_cub(1:this%nkpnts))
    spctrout_cub = M_ZERO
  end if


  Dthetak  = M_ZERO
  if (mdim ==3)  Dthetak = abs(this%thetak_rng(2) - this%thetak_rng(1))/(this%nstepsthetak)
  Dphik = abs(this%phik_rng(2) - this%phik_rng(1))/(this%nstepsphik)
  Lphik = abs(this%phik_rng(2) - this%phik_rng(1))

  ! calculate spectra & total distribution
  do ik = kptst, kptend
    do ist = stst, stend
      do isdim = 1, sdim

        ! orbital spectra
        select case (this%surf_shape)
        case (M_SPHERICAL)

          do ikk = 1, this%nk 
            iomk = 0

            do ith = 0, this%nstepsthetak
              thetak = ith * Dthetak + this%thetak_rng(1)

              if(ith == 0 .or. ith == this%nstepsthetak) then
                weight = (M_ONE - cos(Dthetak / M_TWO)) * Lphik
              else
                weight = abs(cos(thetak - Dthetak / M_TWO) - cos(thetak + Dthetak / M_TWO)) * Dphik
              end if

              do iph = 0, this%nstepsphik - 1
                iomk = iomk + 1
                spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                  abs(this%spctramp_sph(ist, isdim, ik, ikk, iomk))**M_TWO * weight

                if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) exit
              end do
            end do
          end do
          ! distribution
          spctrout_sph(1:this%nk, 1:this%nstepsomegak) = spctrout_sph(1:this%nk, 1:this%nstepsomegak) + &
            abs(this%spctramp_sph(ist, isdim, ik, 1:this%nk, 1:this%nstepsomegak))**M_TWO

        case default !planes or cub

          select case(mdim)
          case(1)
            weight = M_HALF

            ikk = 0
            do ikp = this%nk + 1, this%nkpnts
              ikk = ikk + 1
              spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                (abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO + &
                 abs(this%spctramp_cub(ist, isdim, ik, this%nkpnts + 1 - ikp))**M_TWO) * weight
            end do
     
          case(2)
            weight = Lphik / this%nstepsphik

            ikp = 0
            do ikk = 1, this%nk
              do iph = 0, this%nstepsphik - 1
                ikp = ikp + 1
                spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + & 
                  abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO * weight
              end do
            end do
     
          case(3)
            ikp  = 0
            do ikk = 1, this%nk
              do ith = 0, this%nstepsthetak
                thetak = ith * Dthetak + this%thetak_rng(1)
     
                if(ith == 0 .or. ith == this%nstepsthetak) then
                  weight = (M_ONE - cos(Dthetak / M_TWO)) * Lphik
                else
                  weight = abs(cos(thetak - Dthetak / M_TWO) - cos(thetak + Dthetak / M_TWO)) * Dphik  
                end if
     
                do iph = 0, this%nstepsphik - 1
                  ikp = ikp + 1
                  spctrsum(ist, isdim, ik, ikk) = spctrsum(ist, isdim, ik, ikk) + &
                    abs(this%spctramp_cub(ist, isdim, ik, ikp))**M_TWO * weight
                  if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) exit
                end do
              end do
            end do
          end select
          
          ! distribution
          spctrout_cub(1:this%nkpnts) = spctrout_cub(1:this%nkpnts) + &
            abs(this%spctramp_cub(ist, isdim, ik, 1:this%nkpnts))**M_TWO

        end select
      end do
    end do
  end do

  if(st%parallel_in_states .or. st%d%kpt%parallel) then
#if defined(HAVE_MPI)
    ! total spectrum = sum over all states
    if(this%surf_shape == M_SPHERICAL) then
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
    
    if(energy_resolved) then 
!       iunitone = io_open(trim(path)//'/PES_flux.power.sum', namespace, action='write', position='rewind')
      iunitone = io_open(trim(efile), namespace, action='write', position='rewind')
      write(iunitone, '(a)') '##################################################'                                        
      write(iunitone, '(a1,a18,2x,a18)') '#', str_center("E", 18), str_center("P(E)", 18)
      write(iunitone, '(a1,a18)') '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18)    
      write(iunitone, '(a)') '##################################################'
    end if

    if(momentum_resolved) then
!       iunittwo = io_open(trim(path)//'/PES_flux.distribution.out', namespace, action='write', position='rewind')
      iunittwo = io_open(trim(mfile), namespace, action='write', position='rewind')
      write(iunittwo, '(a)') '##################################################'                               
    end if         

    if (this%surf_shape==M_SPHERICAL) then
      if (momentum_resolved) then 
        write(iunittwo, '(a1,a18,2x,a18,2x,a18,2x,a18)') '#', &
                                          str_center("p", 18), str_center("theta", 18),&
                                          str_center("phi", 18),  str_center("P(p,theta,phi)", 18)
        write(iunittwo, '(a1,a18,2x,a18,2x,a18)') '#', &
                                          str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                          str_center('[None]', 18), &
                                          str_center('[None]', 18)    
        write(iunittwo, '(a)') '##################################################'                                        
      end if
                                        
      do ikk = 1, this%nk 
        kact = this%klinear(ikk,1)
        iomk = 0
        
        if (momentum_resolved) then
          do ith = 0, this%nstepsthetak
            thetak = ith * Dthetak + this%thetak_rng(1)

            do iph = 0, this%nstepsphik - 1
              iomk = iomk + 1
              phik = iph * Dphik + this%phik_rng(1)
              if(iph == 0) iomk_save = iomk
              write(iunittwo,'(4(1x,e18.10E3))') & 
                units_from_atomic(unit_one/units_out%length,kact), thetak, phik, spctrout_sph(ikk, iomk)

              ! just repeat the result for output
              if(this%nstepsphik > 1 .and. iph == (this%nstepsphik - 1)) &
                write(iunittwo,'(4(1x,e18.10E3))') &
                  units_from_atomic(unit_one/units_out%length,kact), thetak, Lphik, spctrout_sph(ikk, iomk_save)

              ! just repeat the result for output and exit
              if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) then
                if(this%nstepsphik > 1) then
                  do iphi = 1, this%nstepsphik
                    phik = iph * Dphik + this%phik_rng(1)
                    write(iunittwo,'(4(1x,e18.10E3))') &
                      units_from_atomic(unit_one/units_out%length,kact), thetak, phik, spctrout_sph(ikk, iomk)
                  end do
                end if
                exit
              end if
            end do

            if(this%nstepsphik > 1 .or. ith == this%nstepsthetak) write(iunittwo, '(1x)', advance='yes')
          end do
        end if        
        
        if (energy_resolved) then
          write(iunitone, '(2(1x,e18.10E3))', advance='no') &
            units_from_atomic(units_out%energy,kact**M_TWO / M_TWO), sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) &
            * units_from_atomic(unit_one/units_out%length,kact)
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do isdim = 1, st%d%dim
                write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) &
                  * units_from_atomic(unit_one/units_out%length,kact)
              end do
            end do
          end do
          write(iunitone, '(1x)', advance='yes')
        end if
        
      end do

    else
      
      
      select case(mdim)
      case(1)
        if (momentum_resolved) then 
          write(iunittwo, '(a1,a18,2x,a18)') '#', &
                                            str_center("p", 18), str_center("P(p)", 18)
          write(iunittwo, '(a1,a18)') '#', str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18)
          write(iunittwo, '(a)') '##################################################'                                        
                                          

          do ikp = 1, this%nkpnts
            write(iunittwo, '(2(1x,e18.10E3))') &
              units_from_atomic(unit_one/units_out%length,this%kcoords_cub(1, ikp, 1)), spctrout_cub(ikp)
          end do
        end if

        if (energy_resolved) then
          do ikk = 1, this%nk
            kact = this%kcoords_cub(1, this%nk + ikk, 1)
            write(iunitone, '(2(1x,e18.10E3))', advance='no') &
              units_from_atomic(units_out%energy,kact**M_TWO / M_TWO), sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) &
              * units_from_atomic(unit_one/units_out%length,kact)
          
            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do isdim = 1, st%d%dim
                  write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) &
                    * units_from_atomic(unit_one/units_out%length,kact)
                end do
              end do
            end do
            write(iunitone, '(1x)', advance='yes')
          end do
        end if

      case(2)
        if (momentum_resolved) then
          write(iunittwo, '(a1,a18,2x,a18,2x,a18)') '#', &
                                            str_center("p", 18), str_center("theta", 18), str_center("P(p,theta)", 18)
          write(iunittwo, '(a1,a18,2x,a18)') '#', &
                                            str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                            str_center('[None]', 18)                                            
          write(iunittwo, '(a)') '##################################################'                                        
        end if

        ikp = 0
        do ikk = 1, this%nk
          kact = this%klinear(ikk,1)
          
          if (momentum_resolved) then
            do iph = 0, this%nstepsphik - 1
              ikp = ikp + 1
              if(iph == 0) ikp_save = ikp
              phik = iph * M_TWO * M_PI / this%nstepsphik
              write(iunittwo,'(3(1x,e18.10E3))') units_from_atomic(unit_one/units_out%length,kact), phik, spctrout_cub(ikp)
            end do
            ! just repeat the result for output
            write(iunittwo,'(3(1x,e18.10E3))') &
              units_from_atomic(unit_one/units_out%length,kact), M_TWO * M_PI, spctrout_cub(ikp_save)
            write(iunittwo, '(1x)', advance='yes')
          end if

          if (energy_resolved) then
            write(iunitone, '(2(1x,e18.10E3))', advance='no') &
              units_from_atomic(units_out%energy,kact**M_TWO / M_TWO), sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) &
              * units_from_atomic(unit_one/units_out%length,kact)
          
            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do isdim = 1, st%d%dim
                  write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) &
                    * units_from_atomic(unit_one/units_out%length,kact)
                end do
              end do
            end do
            write(iunitone, '(1x)', advance='yes')
          end if
        end do

      case(3)
        if (momentum_resolved) then 
          write(iunittwo, '(a1,a18,2x,a18,2x,a18,2x,a18)') '#', &
                                            str_center("p", 18), str_center("theta", 18),&
                                            str_center("phi", 18),  str_center("P(p,theta,phi)", 18)
          write(iunittwo, '(a1,a18,2x,a18,2x,a18)') '#', &
                                            str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                            str_center('[None]', 18), &
                                            str_center('[None]', 18)                                            
          write(iunittwo, '(a)') '##################################################'                                        
        end if
        
        ikp    = 0
        do ikk = 1, this%nk
          kact = this%klinear(ikk,1)

          if(momentum_resolved) then
            do ith = 0, this%nstepsthetak
              thetak = ith * Dthetak + this%thetak_rng(1) 

              do iph = 0, this%nstepsphik - 1
                ikp = ikp + 1

                phik = iph * Dphik + this%phik_rng(1)
                if(iph == 0) ikp_save = ikp
                write(iunittwo,'(4(1x,e18.10E3))') units_from_atomic(unit_one/units_out%length,kact), thetak, phik, spctrout_cub(ikp)

                ! just repeat the result for output
                if(iph == (this%nstepsphik - 1)) &
                  write(iunittwo,'(4(1x,e18.10E3))') &
                    units_from_atomic(unit_one/units_out%length,kact), thetak, Lphik, spctrout_cub(ikp_save)

                ! just repeat the result for output and exit
                if(thetak < M_EPSILON .or. abs(thetak-M_PI) < M_EPSILON) then  
                  do iphi = 1, this%nstepsphik
                    phik = iphi * M_TWO * M_PI / this%nstepsphik
                    write(iunittwo,'(4(1x,e18.10E3))') &
                      units_from_atomic(unit_one/units_out%length, kact), thetak, phik, spctrout_cub(ikp)
                  end do
                  exit
                end if
              end do
              write(iunittwo, '(1x)', advance='yes')
            end do
          end if
          
          if (energy_resolved) then
            write(iunitone, '(2(1x,e18.10E3))', advance='no') &
              units_from_atomic(units_out%energy, kact**M_TWO / M_TWO), sum(sum(sum(spctrsum(:,:,:,ikk),1),1),1) &
              * units_from_atomic(unit_one/units_out%length,kact)
            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do isdim = 1, st%d%dim
                  write(iunitone, '(1x,e18.10E3)', advance='no') spctrsum(ist, isdim, ik, ikk) &
                  * units_from_atomic(unit_one/units_out%length,kact)
                end do
              end do
            end do
            write(iunitone, '(1x)', advance='yes')
          end if
          
        end do

      end select

    end if
    

    call io_close(iunittwo)
    call io_close(iunitone)
    
  end if

  SAFE_DEALLOCATE_A(spctrsum)
  SAFE_DEALLOCATE_A(spctrout_cub)
  SAFE_DEALLOCATE_A(spctrout_sph)

  POP_SUB(pes_flux_out_polar_ascii)
end subroutine pes_flux_out_polar_ascii





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
          if(this%surf_shape == M_SPHERICAL) then
            SAFE_ALLOCATE(psi2(1:this%nk, 1:this%nstepsomegak))
            psi2(:, :) = this%spctramp_sph(ist, idim, ik, :, :)
            call zrestart_write_binary(restart, filename, this%nk * this%nstepsomegak, psi2(:,:), err, root = root)
            SAFE_DEALLOCATE_P(psi2)
            
          else
            SAFE_ALLOCATE(psi1(1:this%nkpnts))
            psi1(:) = this%spctramp_cub(ist, idim, ik, :)
            call zrestart_write_binary(restart, filename, this%nkpnts, psi1(:), err, root = root)
            SAFE_DEALLOCATE_P(psi1)

! print *,mpi_world%rank, filename, err
          end if
        else 
          err = 0  
        end if

        if (err /= 0) ierr = ierr + 1
        itot = itot + 1
      end do
    end do
  end do

  if(this%surf_shape == M_PLANE) then
    root(P_STRATEGY_MAX) = 0
    root(P_STRATEGY_KPOINTS) = -1
    do ik = 1, st%d%nik
      root(P_STRATEGY_DOMAINS) = mod(ik - 1, mesh%mpi_grp%size)
      write(filename,'(a,i5.5)') "pesflux4-kpt", ik

      if (st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
        SAFE_ALLOCATE(psi1(this%nkpnts))
        psi1(:)=this%conjgphase_prev(:,ik)
        call zrestart_write_binary(restart, filename, this%nkpnts, psi1(:), err, root = root)
! print *,mpi_world%rank, filename, err
        SAFE_DEALLOCATE_P(psi1)
      else
        err = 0
      end if
      if (err /= 0) ierr = ierr + 1
    end do
  
  end if

! print *,mpi_world%rank,"suuca"

  if(this%surf_shape == M_SPHERICAL) then
    call zrestart_write_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev, err)
  else 
    if (this%surf_shape /= M_PLANE) &
      call zrestart_write_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev(:,:), err)
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
subroutine pes_flux_load(restart, this, st, ierr)
  type(restart_t),     intent(in)    :: restart
  type(pes_flux_t),    intent(inout) :: this
  type(states_elec_t), intent(in)    :: st
  integer,             intent(out)   :: ierr

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
          if(this%surf_shape == M_SPHERICAL) then
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


  if(this%surf_shape == M_PLANE) then
    do ik = 1, st%d%nik
      write(filename,'(a,i5.5)') "pesflux4-kpt", ik      
    
      if (st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
        SAFE_ALLOCATE(psi1(this%nkpnts))
        call zrestart_read_binary(restart, filename, this%nkpnts, psi1(:), err)
        this%conjgphase_prev(:,ik)=psi1(:)
        SAFE_DEALLOCATE_P(psi1)
      else
        err = 0
      end if
      if (err /= 0) ierr = ierr + 1
    end do
  end if



  if(this%surf_shape == M_SPHERICAL) then
    call zrestart_read_binary(restart, 'pesflux4', this%nk * this%nstepsomegak, this%conjgphase_prev, err)
  else
    if (this%surf_shape /= M_PLANE) &
      call zrestart_read_binary(restart, 'pesflux4', this%nkpnts, this%conjgphase_prev(:,:), err)
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
