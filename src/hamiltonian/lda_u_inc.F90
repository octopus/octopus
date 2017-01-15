!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$


subroutine X(lda_u_apply)(this, mesh, d, ik, psib, hpsib, has_phase)
  type(lda_u_t),      intent(in) :: this
  type(mesh_t),       intent(in) :: mesh
  integer,            intent(in) :: ik
  type(batch_t),      intent(in) :: psib
  type(batch_t),   intent(inout) :: hpsib
  type(states_dim_t), intent(in) :: d
  logical,            intent(in) :: has_phase !True if the wavefunction has an associated phase

  integer :: ibatch, idim, ios, imp, im, ispin
  integer :: bind, is
  R_TYPE  :: weight, reduced
  R_TYPE, allocatable :: psi(:,:), hpsi(:,:)
  R_TYPE, allocatable :: dot(:), epsi(:,:)
  type(orbital_set_t), pointer  :: os

  PUSH_SUB(lda_u_apply)

  SAFE_ALLOCATE(psi(1:mesh%np, 1:d%dim))
  SAFE_ALLOCATE(epsi(1:this%max_np, 1:d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np, 1:d%dim))
  SAFE_ALLOCATE(dot(1:this%maxnorbs))

  ispin = states_dim_get_spin_index(d, ik)

  do ibatch = 1, psib%nst
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib,ibatch, mesh%np, hpsi)

    do ios = 1, this%norbsets 
      ! We have to compute 
      ! hpsi> += sum_m |phi m> sum_m' Vmm' <phi m' | psi >
      !
      ! We first compute <phi m | psi> for all orbitals of the atom
      !
      os => this%orbsets(ios)
      !If we need to add the phase, we explicitly do the operation using the sphere
      !This does not change anything if the sphere occupies the full mesh or not
      if(has_phase) then
        do idim = 1, d%dim
          do is = 1, os%sphere%np
            epsi(is,idim) = psi(os%sphere%map(is), idim)*os%phase(is, ik)
          end do
        end do
      end if

      do im = 1, os%norbs
        !If we need to add the phase, we explicitly do the operation using the sphere
        !This does not change anything if the sphere occupies the full mesh or not
        if(has_phase) then
          dot(im) = X(mf_dotp)(os%sphere%mesh, os%orbitals(im)%X(orb),&
                               epsi(1:os%sphere%np,1), reduce = .false., np = os%sphere%np)
        else
          dot(im) = submesh_to_mesh_dotp(os%sphere, 1, os%orbitals(im)%X(orb),&
                               psi(1:mesh%np,1:d%dim))
        end if
      end do
      !
      do im = 1, os%norbs
        ! sum_m' Vmm' <phi m' | psi >
        reduced = M_ZERO
        do imp = 1, os%norbs
          reduced = reduced + this%X(V)(im,imp,ispin,ios)*dot(imp)
        end do
      
        !In case of phase, we have to apply the conjugate of the phase here
        if(has_phase) then
          do is = 1, os%sphere%np
            epsi(is,1) = os%orbitals(im)%X(orb)(is)*conjg(os%phase(is, ik))
          end do
          do idim = 1, d%dim
           call submesh_add_to_mesh(os%sphere, epsi(1:os%sphere%np,1),&
                                    hpsi(:, idim), reduced)
          end do
        else
          do idim = 1, d%dim
            call submesh_add_to_mesh(os%sphere, os%orbitals(im)%X(orb), &
                                    hpsi(:, idim), reduced)
          end do !idim
        end if
      end do !im
    end do !ios
    call batch_set_state(hpsib, ibatch, mesh%np, hpsi)
  end do !ibatch

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(epsi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(dot)

  POP_SUB(lda_u_apply)

end subroutine X(lda_u_apply)

! ---------------------------------------------------------
!> This routine compute the values of the occupation matrices
! ---------------------------------------------------------
subroutine X(update_occ_matrices)(this, mesh, st, lda_u_energy, phase)
  type(lda_u_t), intent(inout)         :: this
  type(mesh_t),     intent(in)         :: mesh
  type(states_t),  intent(in)          :: st
  FLOAT, intent(inout)                 :: lda_u_energy
  CMPLX, pointer, optional             :: phase(:,:) 

  integer :: ios, im, ik, ist, ispin, norbs, ip, idim, is
  R_TYPE, allocatable :: psi(:,:), epsi(:)
  R_TYPE, allocatable :: dot(:,:)
  FLOAT   :: weight, renorm_weight
  type(orbital_set_t), pointer :: os
  
  PUSH_SUB(update_occ_matrices)

  this%X(n)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(epsi(1:this%max_np))
  SAFE_ALLOCATE(dot(1:this%maxnorbs,1:this%norbsets))

  if(this%useACBN0) &
    this%X(n_alt)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)

  !TODO: use symmetries of the occupation matrices
  do ik = st%d%kpt%start, st%d%kpt%end
    ispin =  states_dim_get_spin_index(st%d,ik)
    do ist = st%st_start, st%st_end
      
      if(this%useACBN0) then
        this%renorm_occ(:,:,:,ist,ik) = M_ZERO
      else
        this%renorm_occ(:,:,:,ist,ik) = M_ONE !st%occ(ist, ik)
      end if
      weight = st%d%kweights(ik)*st%occ(ist, ik) 

      call states_get_state(st, mesh, ist, ik, psi )  
      if(present(phase)) then 
        ! Apply the phase that contains both the k-point and vector-potential terms.
        do idim = 1, st%d%dim
          do ip = 1, mesh%np
            psi(ip, idim) = phase(ip, ik)*psi(ip, idim)
          end do
        end do
      end if
      
      do ios = 1, this%norbsets
        os => this%orbsets(ios)
        norbs = os%norbs
        if(present(phase)) then
          do is = 1, os%sphere%np
            epsi(is) = psi(os%sphere%map(is),1)*os%phase(is, ik)
          end do
        end if
        !We first compute the matrix elemets <\psi | orb_m>
        !taking into account phase correction if needed 
        do im = 1, norbs
          if(present(phase)) then
            dot(im,ios) = X(mf_dotp)(os%sphere%mesh, epsi(1:os%sphere%np),&
                                 os%orbitals(im)%X(orb), reduce = .false., np = os%sphere%np)
          else
            dot(im,ios) = submesh_to_mesh_dotp(os%sphere, 1, os%orbitals(im)%X(orb), &
                                               psi(1:mesh%np,1:st%d%dim))
          end if 
          !We compute the on-site occupation of the site, if needed
          if(this%useACBN0) then
            this%renorm_occ(species_index(os%spec),os%nn,os%ll,ist,ik) = &
                this%renorm_occ(species_index(os%spec),os%nn, os%ll,ist,ik) + abs(dot(im,ios))**2
          end if
        end do
      end do 
     
      !We can compute the (renormalized) occupation matrices
      do ios = 1, this%norbsets
        os => this%orbsets(ios)
        norbs = this%orbsets(ios)%norbs
        do im = 1, norbs
            this%X(n)(1:norbs,im,ispin,ios) = this%X(n)(1:norbs,im,ispin,ios) &
                                         + weight*dot(1:norbs,ios)*R_CONJ(dot(im,ios))
            !We compute the renomalized occupation matrices
            if(this%useACBN0) then
              renorm_weight = this%renorm_occ(species_index(os%spec),os%nn,os%ll,ist,ik)*weight
              this%X(n_alt)(1:norbs,im,ispin,ios) = this%X(n_alt)(1:norbs,im,ispin,ios) &
                                         + renorm_weight*dot(1:norbs,ios)*R_CONJ(dot(im,ios))
            end if 
        end do
       !  call lalg_her( norbs, weight, this%X(n)(1:norbs,1:norbs,ispin,ios), dot)
      end do
    end do
  end do

  
  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(epsi)
  SAFE_DEALLOCATE_A(psi)

#if defined(HAVE_MPI)        
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp%comm, this%X(n))
    if(this%useACBN0) &
      call comm_allreduce(st%st_kpt_mpi_grp%comm, this%X(n_alt))
  end if
#endif      

  if(this%useACBN0 .and. .not.this%freeze_u) & 
    call X(compute_ACBNO_U)(this, st)

  call X(compute_dudarev_energy)(this, lda_u_energy)
  call X(lda_u_update_potential)(this)

  POP_SUB(update_occ_matrices)
end subroutine X(update_occ_matrices)

! ---------------------------------------------------------
!> This routine compute the value of the double counting term in the LDA+U energy
! ---------------------------------------------------------
subroutine X(compute_dudarev_energy)(this, lda_u_energy)
  type(lda_u_t), intent(inout)    :: this
  FLOAT, intent(inout)            :: lda_u_energy
 
  integer :: ios, imp, im, ispin

  PUSH_SUB(compute_dudarev_energy)

  lda_u_energy = M_ZERO

  do ios = 1, this%norbsets
    do ispin = 1, this%nspins
      !TODO: These are matrix operations, that could be optimized
      do im = 1, this%orbsets(ios)%norbs
        do imp = 1, this%orbsets(ios)%norbs
          lda_u_energy = lda_u_energy - CNST(0.5)*this%orbsets(ios)%Ueff*abs(this%X(n)(im,imp,ispin,ios))**2
        end do
        lda_u_energy = lda_u_energy + CNST(0.5)*this%orbsets(ios)%Ueff*this%X(n)(im,im,ispin,ios)
      end do
    end do
  end do

  POP_SUB(compute_dudarev_energy)
end subroutine X(compute_dudarev_energy)


! ---------------------------------------------------------
!> This routine compute the potential that, once multiplied
!> by the projector Pmm' and summed over m and m' for all the atoms
!> gives the full Hubbard potential
! ---------------------------------------------------------
subroutine X(lda_u_update_potential)(this)
  type(lda_u_t), intent(inout)    :: this

  integer :: ios, im, ispin, norbs

  PUSH_SUB(lda_u_update_potential)

  this%X(V) = M_ZERO

  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    do ispin = 1, this%nspins
      do im = 1, norbs
        this%X(V)(1:norbs,im,ispin,ios) = - this%orbsets(ios)%Ueff*this%X(n)(1:norbs,im,ispin,ios)
        this%X(V)(im,im,ispin,ios) = this%X(V)(im,im,ispin,ios) + CNST(0.5)*this%orbsets(ios)%Ueff
      end do
    end do
  end do

  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp%comm, this%X(V))
  end if

  POP_SUB(lda_u_update_potential)
end subroutine X(lda_u_update_potential)


! ---------------------------------------------------------
!> This routine compute the effective U following the expression 
!> given in Agapito et al., Phys. Rev. X 5, 011006 (2015)
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U)(this, st)
  type(lda_u_t), intent(inout)    :: this
  type(states_t),  intent(in)     :: st
  
  integer :: ios, im, imp, impp, imppp, ispin1, ispin2, norbs
  FLOAT   :: numU, numJ, denomU, denomJ, tmpU, tmpJ

  PUSH_SUB(compute_ACBNO_U)

  do ios = 1, this%norbsets
    norbs = this%orbsets(ios)%norbs
    numU = M_ZERO
    numJ = M_ZERO
    denomU = M_ZERO
    denomJ = M_ZERO

    if(norbs > 1) then

      do im = 1, norbs
      do imp = 1,norbs
        do impp = 1, norbs
        do imppp = 1, norbs
          ! We first compute the terms
          ! sum_{alpha,beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
          ! sum_{alpha} P^alpha_{mmp}P^alpha_{mpp,mppp}
          tmpU = M_ZERO
          tmpJ = M_ZERO
          do ispin1 = 1, st%d%nspin
            do ispin2 = 1, st%d%nspin
              tmpU = tmpU + this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin2,ios)
            end do
            tmpJ = tmpJ + this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin1,ios)
          end do
          ! These are the numerator of the ACBN0 U and J
          numU = numU + tmpU*this%coulomb(im,imp,impp,imppp,ios)
          numJ = numJ + tmpJ*this%coulomb(im,imppp,impp,imp,ios) 
        end do
        end do

        ! We compute the term
        ! sum_{alpha} sum_{m,mp/=m} N^alpha_{m}N^alpha_{mp}
        tmpJ = M_ZERO
        if(imp/=im) then
          do ispin1 = 1, st%d%nspin
            tmpJ = tmpJ + this%X(n)(im,im,ispin1,ios)*this%X(n)(imp,imp,ispin1,ios)
          end do
        end if
        denomJ = denomJ + tmpJ
        denomU = denomU + tmpJ

        ! We compute the term
        ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
        tmpU = M_ZERO
        do ispin1 = 1, st%d%nspin
          do ispin2 = 1, st%d%nspin
            if(ispin1 /= ispin2) then
              tmpU = tmpU + this%X(n)(im,im,ispin1,ios)*this%X(n)(imp,imp,ispin2,ios)
            end if
          end do
        end do
        denomU = denomU + tmpU
 
      end do
      end do

      this%orbsets(ios)%Ueff = numU/denomU - numJ/denomJ
 
  else !In the case of s orbitals, the expression is different
    if(st%d%nspin > 1) then
      ! sum_{alpha/=beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
      tmpU = M_ZERO
      do ispin1 = 1, st%d%nspin
        do ispin2 = 1, st%d%nspin
          if(ispin1 /= ispin2) then
            tmpU = tmpU + this%X(n_alt)(1,1,ispin1,ios)*this%X(n_alt)(1,1,ispin2,ios)
          end if
        end do
      end do
      ! These are the numerator of the ACBN0 U
      numU = tmpU*this%coulomb(1,1,1,1,ios)

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      denomU = M_ZERO
      do ispin1 = 1, st%d%nspin
        do ispin2 = 1, st%d%nspin
          if(ispin1 /= ispin2) then
            denomU = denomU + this%X(n)(1,1,ispin1,ios)*this%X(n)(1,1,ispin2,ios)
          end if
        end do
      end do

     this%orbsets(ios)%Ueff = numU/denomU
   else
     this%orbsets(ios)%Ueff = M_ZERO
   end if
  end if
  end do

  POP_SUB(compute_ACBNO_U)  
end subroutine X(compute_ACBNO_U)

! ---------------------------------------------------------
! TODO: Merge this with the two_body routine in system/output_me_inc.F90
subroutine X(compute_coulomb_integrals) (this, mesh, st)
  type(lda_u_t), intent(inout)  :: this
  type(mesh_t),     intent(in)  :: mesh
  type(states_t),   intent(in)  :: st

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: norbs, np_sphere, ios, ip
  FLOAT, allocatable :: tmp(:), vv(:), nn(:)
  type(orbital_t), pointer :: orbi, orbj, orbk, orbl

  PUSH_SUB(X(compute_coulomb_integrals))

  ASSERT(.not. st%parallel_in_states)
  
  SAFE_ALLOCATE(nn(1:this%max_np))
  SAFE_ALLOCATE(vv(1:this%max_np))
  SAFE_ALLOCATE(tmp(1:this%max_np))

  SAFE_ALLOCATE(this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%norbsets))
  this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%norbsets) = M_ZERO

  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    np_sphere = this%orbsets(ios)%sphere%np  
 
    ijst=0
    do ist = 1, norbs
      orbi => this%orbsets(ios)%orbitals(ist) 
      if( st%d%nspin == 1 .and. np_sphere == 1 ) cycle
      
      do jst = 1, norbs
        if(jst > ist) cycle
        ijst=ijst+1
        orbj => this%orbsets(ios)%orbitals(jst)

        nn(1:np_sphere)  = real(orbi%X(orb)(1:np_sphere))*real(orbj%X(orb)(1:np_sphere))
        !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
        call dpoisson_solve_sm(psolver, this%orbsets(ios)%sphere, vv(1:np_sphere), nn(1:np_sphere))

        klst=0
        do kst = 1, norbs
          orbk => this%orbsets(ios)%orbitals(kst)
          do lst = 1, norbs
            if(lst > kst) cycle
            klst=klst+1
           if(klst > ijst) cycle

            orbl => this%orbsets(ios)%orbitals(lst)

            do ip=1,np_sphere
             tmp(ip) = vv(ip)*real(orbl%X(orb)(ip))*real(orbk%X(orb)(ip))
            end do
            this%coulomb(ist,jst,kst,lst,ios) = dsm_integrate(mesh, this%orbsets(ios)%sphere, tmp(1:np_sphere))

            if(abs(this%coulomb(ist,jst,kst,lst,ios))<CNST(1.0e-12)) then
              this%coulomb(ist,jst,kst,lst,ios) = M_ZERO
            else
              this%coulomb(kst,lst,ist,jst,ios) = this%coulomb(ist,jst,kst,lst,ios)              
              this%coulomb(jst,ist,lst,kst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(lst,kst,jst,ist,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(jst,ist,kst,lst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(lst,kst,ist,jst,ios) = this%coulomb(ist,jst,kst,lst,ios)              
              this%coulomb(ist,jst,lst,kst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(kst,lst,jst,ist,ios) = this%coulomb(ist,jst,kst,lst,ios)              
            end if
          end do !lst
        end do !kst
      end do !jst
    end do !ist
  end do !ia

  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp%comm, this%coulomb)
  end if
 
  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(compute_coulomb_integrals))
end subroutine X(compute_coulomb_integrals)

 ! ---------------------------------------------------------
 !> This routine computes [r,V_lda+u].
 ! ---------------------------------------------------------
 subroutine X(lda_u_commute_r)( ) !this, mesh, ik, psi, gpsi)
 !   type(scissor_t), intent(in)    :: this
 !   type(mesh_t),    intent(in)    :: mesh 
 !   R_TYPE,          intent(in)    :: psi(:,:)
 !   integer,         intent(in)    :: ik
 !   R_TYPE,          intent(inout) :: gpsi(:, :, :)
 !
 !   integer :: ist, idim, idir
 !   R_TYPE  :: dot
 !   R_TYPE, allocatable :: gspsi(:,:), tmpstate(:,:), psi_r(:)
 !
    PUSH_SUB(lda_u_commute_r)

    POP_SUB(lda_u_commute_r)
 end subroutine X(lda_u_commute_r)

! ---------------------------------------------------------
!> This routine is an interface for constructing the orbital basis.
! ---------------------------------------------------------
subroutine X(construct_orbital_basis)(this, geo, mesh, st)
  type(lda_u_t),             intent(inout)    :: this
  type(geometry_t), target,  intent(in)       :: geo
  type(mesh_t),              intent(in)       :: mesh
  type(states_t),            intent(in)       :: st 

  integer :: ia, iorb, norb, ispin, offset
  integer ::  hubbardl, ii, nn, ll, mm, work, work2, iorbset
  FLOAT   :: norm
  logical :: hasSorbitals
  type(orbital_set_t), pointer :: os

  PUSH_SUB(X(construct_orbital_basis))

  offset = 0
  if(this%skipSOrbitals) offset = 1

  write(message(1),'(a)')    'Building the LDA+U localized orbital basis.'
  call messages_info(1)

  !We first count the number of orbital sets we have to treat
  norb = 0
  if( .not. this%useAllOrbitals ) then
    do ia = 1, geo%natoms
      hubbardl = species_hubbard_l(geo%atom(ia)%species)
      if( hubbardl .eq. 0 ) cycle
      norb = norb+1
    end do
  else
    do ia = 1, geo%natoms
      work = 0
      hasSorbitals = .false.
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm) 
        if(ll == 0) hasSorbitals = .true.
        work = max(work, ii)
      end do
      if(this%skipSOrbitals .and. hasSorbitals ) work = work-1
      norb = norb + work
    end do
  end if

  write(message(1),'(a, i3, a)')    'Found ', norb, ' orbital sets.'
  call messages_info(1)

  this%norbsets = norb
  SAFE_ALLOCATE(this%orbsets(1:norb))

  iorbset = 0
  do ia = 1, geo%natoms
    if(.not. this%useAllOrbitals) then
      hubbardl = species_hubbard_l(geo%atom(ia)%species)
      if( hubbardl .eq. 0 ) cycle
      iorbset = iorbset + 1
      os => this%orbsets(iorbset)
      norb = 0
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
        call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
        if(ll .eq. hubbardl ) then 
          norb = norb + 1
          os%ll = hubbardl
          os%nn = nn
        end if
      end do
      os%norbs = norb
      SAFE_ALLOCATE(os%orbitals(1:norb))
      os%Ueff = species_hubbard_u(geo%atom(ia)%species)
      os%spec => geo%atom(ia)%species
      call submesh_null(os%sphere)
      norb = 0
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
        if(ll .eq. hubbardl ) then
          norb  = norb + 1
         ! We obtain the orbital
          call X(get_atomic_orbital)(geo, mesh, os%sphere, ia, iorb, 1, os%orbitals(norb), &
                        this%truncation, this%orbitals_threshold)
          ! We have to normalize the orbitals, 
          ! in case the orbitals that comes out of the pseudo are not properly normalised
          norm = X(sm_nrm2)(os%sphere, os%orbitals(norb)%X(orb)(1:os%sphere%np))
          os%orbitals(norb)%X(orb)(1:os%sphere%np) =  &
                os%orbitals(norb)%X(orb)(1:os%sphere%np) /sqrt(norm)
        endif
      end do
    else !useAllOrbitals
      work = 0
      hasSorbitals = .false.
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
        if(ll == 0) hasSorbitals = .true.
        work = max(work, ii)          
      end do
      if(this%skipSOrbitals .and. hasSorbitals ) work = work-1

      !We loop over the orbital sets of the atom ia
      do norb = 1, work
        os => this%orbsets(iorbset+norb)
        !We count the number of orbital for this orbital set
        work2 = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
          if(ii == norb+offset) then
            work2 = work2 + 1
            os%ll = ll
            os%nn = nn
          end if
        end do
        os%norbs = work2
        SAFE_ALLOCATE(os%orbitals(1:os%norbs))
        os%Ueff = species_hubbard_u(geo%atom(ia)%species)
        os%spec => geo%atom(ia)%species
        call submesh_null(os%sphere)        
 
        work2 = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          if(ii == norb+offset) then
            work2  = work2 + 1
            ! We obtain the orbital
            call X(get_atomic_orbital)(geo, mesh, os%sphere, ia, iorb, 1, os%orbitals(work2),&
                         this%truncation, this%orbitals_threshold)
            ! We have to normalize the orbitals, 
            ! in case the orbitals that comes out of the pseudo are not properly normalised
            norm = X(sm_nrm2)(os%sphere, os%orbitals(work2)%X(orb)(1:os%sphere%np))
            os%orbitals(work2)%X(orb)(1:os%sphere%np) =  &
                 os%orbitals(work2)%X(orb)(1:os%sphere%np) /sqrt(norm)
          endif
        end do !iorb
      end do !norb
      iorbset = iorbset + work
    end if
  end do

  this%maxnorbs = 0
  this%max_np = 0
  do iorbset = 1, this%norbsets
    os => this%orbsets(iorbset)
    if( os%norbs > this%maxnorbs ) this%maxnorbs = os%norbs

    ! In case of complex wavefunction, we allocate the array for the phase correction
  #ifdef R_TCOMPLEX
    SAFE_ALLOCATE(os%phase(1:os%sphere%np, st%d%kpt%start:st%d%kpt%end))
    os%phase(:,:) = M_ZERO
  #endif
    ! We need to know the maximum number of points in order to allocate a temporary array
    ! to apply the phase in lda_u_apply
    if(os%sphere%np > this%max_np) this%max_np = os%sphere%np
  end do  

  do iorbset = 1, this%norbsets
    write(message(1),'(a,i2,a,f8.5,a)')    'Orbital set ', iorbset, ' has a value of U of ',&
                         this%orbsets(iorbset)%Ueff   , ' Ha'
    write(message(2),'(a,i2,a)')    ' and', this%orbsets(iorbset)%norbs, ' orbitals.'
     call messages_info(2)
  end do 
 
  POP_SUB(X(construct_orbital_basis))

end subroutine X(construct_orbital_basis)


! ---------------------------------------------------------
!> This routine returns the atomic orbital basis -- provided
!! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(get_atomic_orbital) (geo, mesh, sm, iatom, iorb, ispin, orb, truncation, threshold)
  type(mesh_t),             intent(in)    :: mesh
  type(geometry_t), target, intent(in)    :: geo
  type(submesh_t),          intent(inout) :: sm
  integer,                  intent(in)    :: iatom
  integer,                  intent(in)    :: iorb
  integer,                  intent(in)    :: ispin
  type(orbital_t),          intent(inout) :: orb
  integer,                  intent(in)    :: truncation
  FLOAT,                    intent(in)    :: threshold

  type(species_t), pointer :: spec
  integer :: ii, ll, mm, ispin_
  FLOAT :: radius
  logical :: complex_ylms

  PUSH_SUB(X(get_atomic_orbital))

  ispin_ = ispin 

  spec => geo%atom(iatom)%species
  ASSERT(iorb <= species_niwfs(spec))

  nullify(orb%dorb)
  nullify(orb%zorb)

  call species_iwf_ilm(spec, iorb, ispin_, ii, ll, mm)

  if(sm%np == -1) then
    if(truncation == OPTION__ORBITALSTRUNCATIONMETHOD__FULL) then
      radius = species_get_iwf_radius(spec, ii, ispin_, threshold) 
    else
      radius = species_get_iwf_radius(spec, ii, ispin_)
   
      if(truncation == OPTION__ORBITALSTRUNCATIONMETHOD__BOX) then
        ! if the orbital is larger than the size of the box, we restrict it to this size, 
        ! otherwise the orbital will overlap more than one time with the simulation box.
        ! This would induces phase problem if the complete mesh is used instead of the sphere
        radius = min(radius, minval(mesh%sb%lsize(1:mesh%sb%dim)-mesh%spacing(1:mesh%sb%dim)*CNST(1.01)))
      else
        !If asked, we truncate the orbital to the radius on the projector spheres 
        !of the NL part of the pseudopotential.
        !This is a way to garanty no overlap between orbitals of different atoms.
        if(species_is_ps(spec)) &
          radius = min(radius,species_get_ps_radius(spec))
      end if
    end if
    ! make sure that if the spacing is too large, the orbitals fit in a few points at least
    radius = max(radius, CNST(2.0)*maxval(mesh%spacing(1:mesh%sb%dim)))
 
    !We initialise the submesh corresponding to the orbital 
    call submesh_init(sm, mesh%sb, mesh, geo%atom(iatom)%x, radius)

    if(radius >= minval(mesh%sb%lsize(1:mesh%sb%dim)-mesh%spacing(1:mesh%sb%dim))) then
      write(message(1),'(a,i5,a)')  'This is an extended orbital, with ', sm%np, ' grid points.'
      write(message(2),'(a,f8.5,a)') 'The radius is ', radius, ' Bohr.'
      call messages_info(2)
    end if
  end if

  !We allocate both the orbital on the submesh and on the complete mesh
  SAFE_ALLOCATE(orb%X(orb)(1:sm%np))
  orb%X(orb) = M_ZERO

  !This is a bit dirty.
  complex_ylms = .false.

  !We get the orbital from the pseudopotential
  #ifdef R_TCOMPLEX
  if(.not. complex_ylms) then
    !In this case we want to get a real orbital and to store it in complex array
    SAFE_ALLOCATE(orb%dorb(1:sm%np))
    call dspecies_get_orbital_submesh(spec, sm, ii, ll, mm, ispin_, geo%atom(iatom)%x, &
                                            orb%dorb)
    orb%X(orb)(1:sm%np) = orb%dorb(1:sm%np)
    SAFE_DEALLOCATE_P(orb%dorb)
  else
  #endif
    call X(species_get_orbital_submesh)(spec, sm, ii, ll, mm, ispin_, geo%atom(iatom)%x,&
                                         orb%X(orb))
  #ifdef R_TCOMPLEX
  end if
  #endif

  POP_SUB(X(get_atomic_orbital))

end subroutine X(get_atomic_orbital)

 ! ---------------------------------------------------------
  subroutine X(lda_u_set_occupations)(this, occ)
    type(lda_u_t),  intent(inout) :: this
    R_TYPE,         intent(in)    :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs

    PUSH_SUB(X(lda_u_set_occupations))

    ind = 0
    do ios = 1, this%norbsets
      norbs = this%orbsets(ios)%norbs
      do ispin = 1, this%nspins
        do im = 1, norbs
          do imp = 1,norbs
            ind = ind + 1
            this%X(n)(im,imp,ispin,ios) = occ(ind)
          end do
        end do
      end do
    end do 

    POP_SUB(X(lda_u_set_occupations))
  end subroutine X(lda_u_set_occupations)

   ! ---------------------------------------------------------
  subroutine X(lda_u_get_occupations)(this, occ)
    type(lda_u_t),  intent(in) :: this
    R_TYPE,      intent(inout) :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs

    PUSH_SUB(X(lda_u_get_occupations))

    ind = 0
    do ios = 1, this%norbsets
      norbs = this%orbsets(ios)%norbs
      do ispin = 1, this%nspins
        do im = 1, norbs
          do imp = 1,norbs
            ind = ind + 1
            occ(ind) = this%X(n)(im,imp,ispin,ios)
          end do
        end do
      end do
    end do

    POP_SUB(X(lda_u_get_occupations))
  end subroutine X(lda_u_get_occupations)

